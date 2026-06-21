"""S8 guards for the scATAC evidence adapter (ADR-057).

The adapter flattens audited scATAC findings into EvidenceSignal items across
per-cluster TF-motif enrichment (tf_motif), regulatory gene-activity scores
(gene, moderate caveat), and peak-to-gene links (gene, associative caveat),
emitting a layer only when its ledger node ran.
"""

from __future__ import annotations

from aria.agents.narrative.hypothesis import scatac_evidence


def _agent_results() -> dict:
    return {
        "chromatin_agent": {
            "findings": {
                "motifs": {
                    "ran": True,
                    "padj_max": 0.05,
                    "per_group": {
                        "cluster_3": {
                            "top_motifs": [
                                {"name": "EOMES", "log2_enrichment": 2.7,
                                 "padj": 1e-4},
                                {"name": "NOISE", "log2_enrichment": 0.1,
                                 "padj": 0.4},  # not enriched
                            ]
                        }
                    },
                },
                "regulatory": {
                    "ran": True,
                    "gene_scores": {
                        "ran": True,
                        "top_genes": [
                            {"gene": "CD8A", "mean_gene_activity_score": 4.1},
                            {"gene": "GZMK", "mean_gene_activity_score": 3.3},
                        ],
                    },
                    "peak_to_gene": {
                        "ran": True,
                        "top_links": [
                            {"gene": "IFNG", "correlation": 0.62,
                             "direction": "positive"},
                            {"gene": "LEF1", "correlation": -0.55,
                             "direction": "negative"},
                        ],
                    },
                },
            }
        }
    }


def _ledger(motif="ran", reg="ran") -> dict:
    return {
        "entries": [
            {"node_id": "ledger://chromatin/motif_enrichment", "status": motif},
            {"node_id": "ledger://chromatin/regulatory_layers", "status": reg},
        ]
    }


def test_adapter_empty_without_findings():
    assert scatac_evidence({}) == []


def test_per_cluster_tf_motifs():
    signals = scatac_evidence(_agent_results(), _ledger())
    tfs = {s.entity for s in signals if s.entity_kind == "tf_motif"}
    assert tfs == {"EOMES"}  # NOISE excluded (padj >= padj_max)
    eomes = next(s for s in signals if s.entity == "EOMES")
    assert "cluster_3" in eomes.measure
    assert eomes.audited_node_ref == "ledger://chromatin/motif_enrichment"


def test_gene_activity_signals_with_moderate_caveat():
    signals = scatac_evidence(_agent_results(), _ledger())
    ga = {s.entity: s for s in signals if s.measure == "gene_activity"}
    assert set(ga) == {"CD8A", "GZMK"}
    assert ga["CD8A"].entity_kind == "gene"
    assert ga["CD8A"].audited_node_ref == "ledger://chromatin/regulatory_layers"
    assert "moderate" in " ".join(ga["CD8A"].caveats_inherited).lower()


def test_peak2gene_signals_with_direction_and_caveat():
    signals = scatac_evidence(_agent_results(), _ledger())
    p2g = {s.entity: s for s in signals if s.measure == "peak2gene"}
    assert set(p2g) == {"IFNG", "LEF1"}
    assert p2g["IFNG"].direction == "up"
    assert p2g["LEF1"].direction == "down"
    assert "associative" in " ".join(p2g["IFNG"].caveats_inherited).lower()


def test_regulatory_skipped_when_node_not_run():
    signals = scatac_evidence(_agent_results(), _ledger(reg="skipped"))
    assert not any(
        s.audited_node_ref == "ledger://chromatin/regulatory_layers"
        for s in signals
    )
    assert any(s.entity_kind == "tf_motif" for s in signals)


def test_motifs_skipped_when_node_not_run():
    signals = scatac_evidence(_agent_results(), _ledger(motif="skipped"))
    assert not any(s.entity_kind == "tf_motif" for s in signals)
    assert any(s.measure == "gene_activity" for s in signals)
