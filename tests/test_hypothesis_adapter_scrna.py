"""S6 guards for the scRNA evidence adapter (ADR-057).

The adapter flattens AUDITED scRNA findings into EvidenceSignal items across
pseudobulk DE (genes), differential abundance (cell types), and LIANA comms
(ligand/receptor genes + source/target cell types), inheriting the unaddressed
technical confounds (ambient/doublets/batch) and emitting a layer only when its
ledger node ran.
"""

from __future__ import annotations

from aria.agents.narrative.hypothesis import scrna_evidence


def _agent_results() -> dict:
    return {
        "scrna_agent": {
            "findings": {
                # QC without ambient correction or doublet handling -> confounds.
                "qc": {"n_cells": 5000},
                "pseudobulk_de": {
                    "per_group": {
                        "CD8 T": {
                            "per_comparison": {
                                "old_vs_young": {
                                    "status": "success",
                                    "top_genes": [
                                        {"gene": "GZMB", "log2fc": 1.9},
                                        {"gene": "IL7R", "log2fc": -1.4},
                                    ],
                                }
                            }
                        }
                    }
                },
                "differential_abundance": {
                    "per_comparison": {
                        "old_vs_young": {
                            "status": "success",
                            "per_cell_type": [
                                {"name": "CD8 T", "direction": "up",
                                 "significant": True},
                                {"name": "B cells", "direction": "down",
                                 "significant": False},  # skipped
                            ],
                        }
                    }
                },
                "cell_communication": {
                    "status": "success",
                    "top_interactions": [
                        {"source": "Monocyte", "target": "CD8 T",
                         "ligand": "IL15", "receptor": "IL2RB"},
                    ],
                },
            }
        }
    }


def _ledger(de="ran", da="ran", ccc="ran") -> dict:
    return {
        "entries": [
            {"node_id": "ledger://scRNA/pseudobulk_de", "status": de},
            {"node_id": "ledger://scRNA/differential_abundance", "status": da},
            {"node_id": "ledger://scRNA/cell_communication", "status": ccc},
        ]
    }


def test_adapter_empty_without_scrna_findings():
    assert scrna_evidence({}) == []


def test_pseudobulk_genes_with_direction_and_node():
    signals = scrna_evidence(_agent_results(), _ledger())
    genes = {s.entity: s for s in signals
             if s.audited_node_ref == "ledger://scRNA/pseudobulk_de"}
    assert set(genes) == {"GZMB", "IL7R"}
    assert genes["GZMB"].direction == "up"
    assert genes["IL7R"].direction == "down"
    assert genes["GZMB"].entity_kind == "gene"


def test_only_significant_abundance_celltypes():
    signals = scrna_evidence(_agent_results(), _ledger())
    cts = {s.entity for s in signals
           if s.audited_node_ref == "ledger://scRNA/differential_abundance"}
    assert cts == {"CD8 T"}  # B cells not significant
    ct = next(s for s in signals
              if s.audited_node_ref == "ledger://scRNA/differential_abundance")
    assert ct.entity_kind == "celltype"
    assert ct.direction == "up"


def test_liana_emits_ligand_receptor_and_celltypes():
    signals = scrna_evidence(_agent_results(), _ledger())
    ccc = [s for s in signals
           if s.audited_node_ref == "ledger://scRNA/cell_communication"]
    genes = {s.entity for s in ccc if s.entity_kind == "gene"}
    cts = {s.entity for s in ccc if s.entity_kind == "celltype"}
    assert genes == {"IL15", "IL2RB"}
    assert "Monocyte" in cts and "CD8 T" in cts


def test_global_confounds_inherited():
    signals = scrna_evidence(_agent_results(), _ledger())
    gzmb = next(s for s in signals if s.entity == "GZMB")
    blob = " ".join(gzmb.caveats_inherited).lower()
    assert "ambient" in blob
    assert "doublet" in blob


def test_layer_skipped_when_ledger_node_not_run():
    signals = scrna_evidence(_agent_results(), _ledger(de="skipped"))
    assert not any(
        s.audited_node_ref == "ledger://scRNA/pseudobulk_de" for s in signals
    )
    assert any(
        s.audited_node_ref == "ledger://scRNA/differential_abundance"
        for s in signals
    )


def test_celltype_appears_once_per_node():
    # CD8 T is both a DA cell type and a LIANA target; the per-(entity,node)
    # dedup keeps one per node, not a single global collapse.
    signals = scrna_evidence(_agent_results(), _ledger())
    cd8_nodes = [
        s.audited_node_ref for s in signals if s.entity == "CD8 T"
    ]
    assert "ledger://scRNA/differential_abundance" in cd8_nodes
    assert "ledger://scRNA/cell_communication" in cd8_nodes
