"""S5 guards for the bulk RNA evidence adapter (ADR-057).

The adapter is pure and deterministic: it flattens AUDITED bulk RNA findings into
EvidenceSignal items (real gene symbols + pathway terms, each tagged with the
run-ledger node that produced it and the confounds that node flags), and emits a
layer only when its ledger node actually ran.
"""

from __future__ import annotations

from aria.agents.narrative.hypothesis import bulk_rna_evidence


def _agent_results() -> dict:
    return {
        "bulk_rna_agent": {
            "findings": {
                "padj_threshold": 0.05,
                "low_power_warning": True,
                "contrasts": [
                    {
                        "name": "old_vs_young",
                        "n_significant": 2,
                        "top_genes": [
                            {"symbol": "GATA1", "log2fc": 2.3, "padj": 1e-4},
                            {"gene": "MYC", "log2fc": -1.8, "padj": 1e-3},
                            {"symbol": "", "log2fc": 1.0},  # skipped: no symbol
                            {"symbol": "NODATA"},           # skipped: no lfc
                        ],
                        "pathways": {
                            "ora": [
                                {"term": "Cell cycle"},
                                {"name": "Apoptosis"},
                            ]
                        },
                    }
                ],
            }
        }
    }


def _ledger(de="ran", ora="ran") -> dict:
    return {
        "entries": [
            {"node_id": "ledger://bulk/differential_expression", "status": de},
            {"node_id": "ledger://bulk/pathway_enrichment", "status": ora},
        ]
    }


def test_adapter_returns_empty_without_bulk_findings():
    assert bulk_rna_evidence({}) == []
    assert bulk_rna_evidence({"bulk_rna_agent": {"findings": {}}}) == []


def test_adapter_flattens_genes_with_direction_and_node():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    genes = {s.entity: s for s in signals if s.entity_kind == "gene"}
    assert set(genes) == {"GATA1", "MYC"}
    assert genes["GATA1"].direction == "up"
    assert genes["GATA1"].value == 2.3
    assert genes["MYC"].direction == "down"
    assert all(
        g.audited_node_ref == "ledger://bulk/differential_expression"
        for g in genes.values()
    )


def test_adapter_inherits_low_replication_caveat():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    gata1 = next(s for s in signals if s.entity == "GATA1")
    assert "low replication" in gata1.caveats_inherited


def test_adapter_emits_pathway_terms():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    terms = {s.entity for s in signals if s.entity_kind == "pathway"}
    assert terms == {"Cell cycle", "Apoptosis"}
    assert all(
        s.audited_node_ref == "ledger://bulk/pathway_enrichment"
        for s in signals
        if s.entity_kind == "pathway"
    )


def test_adapter_skips_layer_whose_ledger_node_did_not_run():
    # DE skipped -> no gene signals; ORA ran -> pathway signals only.
    signals = bulk_rna_evidence(_agent_results(), _ledger(de="skipped"))
    assert not any(s.entity_kind == "gene" for s in signals)
    assert any(s.entity_kind == "pathway" for s in signals)


def test_adapter_without_ledger_emits_audited_inputs():
    signals = bulk_rna_evidence(_agent_results(), run_ledger=None)
    assert any(s.entity == "GATA1" for s in signals)


def test_adapter_flags_unmodeled_batch():
    results = _agent_results()
    exp_ctx = {"design": {"has_batch": True, "covariates": ["age"]}}
    signals = bulk_rna_evidence(results, _ledger(), exp_ctx)
    gata1 = next(s for s in signals if s.entity == "GATA1")
    assert any("batch" in c for c in gata1.caveats_inherited)
