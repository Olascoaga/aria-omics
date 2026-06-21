"""S7 guards for the bulk ATAC evidence adapter (ADR-057).

Bulk ATAC's inline, named, hypothesis-relevant entities are the TFs whose motifs
are enriched in a condition's DA peaks. The adapter flattens that audited
enrichment into tf_motif EvidenceSignal items (node ledger://chromatin/
motif_enrichment), keeps only motifs enriched at the run's threshold, carries the
peak group in `measure`, and inherits the permanent bulk ATAC caveats (motif =
associative, not binding; low replication when DA flagged it).
"""

from __future__ import annotations

from aria.agents.narrative.hypothesis import bulk_atac_evidence


def _agent_results(low_power=True) -> dict:
    return {
        "chromatin_agent": {
            "findings": {
                "motifs": {
                    "ran": True,
                    "padj_max": 0.05,
                    "per_group": {
                        "K562_up": {
                            "top_motifs": [
                                {"name": "KLF1", "log2_enrichment": 3.2,
                                 "padj": 1e-5},
                                {"name": "GATA1", "log2_enrichment": 2.8,
                                 "padj": 1e-4},
                                {"name": "NOISE", "log2_enrichment": 0.1,
                                 "padj": 0.4},  # not enriched -> skipped
                            ]
                        },
                        "GM12878_up": {
                            "top_motifs": [
                                {"name": "IRF1", "log2_enrichment": 2.5,
                                 "padj": 1e-3},
                            ]
                        },
                    },
                },
                "differential_accessibility": {
                    "comparisons": [
                        {"low_power_warning": low_power},
                    ]
                },
            }
        }
    }


def _ledger(motif="ran") -> dict:
    return {
        "entries": [
            {"node_id": "ledger://chromatin/motif_enrichment", "status": motif},
        ]
    }


def test_adapter_empty_without_chromatin_findings():
    assert bulk_atac_evidence({}) == []


def test_adapter_empty_when_motifs_not_run():
    results = _agent_results()
    results["chromatin_agent"]["findings"]["motifs"]["ran"] = False
    assert bulk_atac_evidence(results, _ledger()) == []


def test_flattens_enriched_tfs_only():
    signals = bulk_atac_evidence(_agent_results(), _ledger())
    tfs = {s.entity for s in signals}
    assert tfs == {"KLF1", "GATA1", "IRF1"}  # NOISE excluded (padj >= padj_max)
    assert all(s.entity_kind == "tf_motif" for s in signals)
    assert all(
        s.audited_node_ref == "ledger://chromatin/motif_enrichment"
        for s in signals
    )


def test_measure_carries_peak_group_and_value():
    signals = bulk_atac_evidence(_agent_results(), _ledger())
    klf1 = next(s for s in signals if s.entity == "KLF1")
    assert "K562_up" in klf1.measure
    assert klf1.value == 3.2
    irf1 = next(s for s in signals if s.entity == "IRF1")
    assert "GM12878_up" in irf1.measure


def test_inherits_associative_and_low_replication_caveats():
    signals = bulk_atac_evidence(_agent_results(low_power=True), _ledger())
    blob = " ".join(signals[0].caveats_inherited).lower()
    assert "associative" in blob
    assert "low replication" in blob


def test_low_replication_caveat_absent_when_not_flagged():
    signals = bulk_atac_evidence(_agent_results(low_power=False), _ledger())
    blob = " ".join(signals[0].caveats_inherited).lower()
    assert "associative" in blob
    assert "low replication" not in blob


def test_layer_skipped_when_ledger_node_not_run():
    assert bulk_atac_evidence(_agent_results(), _ledger(motif="skipped")) == []
