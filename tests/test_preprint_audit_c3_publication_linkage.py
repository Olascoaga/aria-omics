"""C3: public result claims fail closed without typed evidence/ledger nodes."""

from __future__ import annotations

import json

import pytest

from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
from aria.agents.narrative.run_ledger import (
    LedgerLinkageError,
    build_run_ledger,
    ensure_report_ledger_nodes,
    verify_blocks_against_ledger,
)
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock


def _result_block(*, analysis: str = "peak_annotation") -> NarrativeBlock:
    return NarrativeBlock(
        id=f"chromatin.{analysis}",
        modality="chromatin",
        analysis=analysis,
        block_type="result",
        title="Chromatin result",
        status="success",
        confidence="medium",
        claim="Genomic annotation classified 10 peaks.",
        evidence=[EvidenceItem(
            label="Peaks annotated",
            value=10,
            source="chromatin_peak_annotation",
        )],
    )


def test_public_compiler_withholds_supported_claim_without_ledger_node():
    block = _result_block()

    compiled = compile_public_claims([block], {}, run_ledger={"entries": []})

    assert compiled.blocks == []
    assert compiled.claims == []
    assert compiled.summary()["n_withheld"] == 1
    assert block.claim not in str(compiled.withheld)


def test_public_claim_resolves_to_typed_ledger_and_evidence_nodes():
    block = _result_block()
    results = {"chromatin_agent": {"findings": {
        "peak_annotation": {"ran": True, "status": "success"},
    }}}
    ledger = build_run_ledger({}, results)

    compiled = compile_public_claims([block], {}, run_ledger=ledger)

    assert len(compiled.claims) == 1
    claim = compiled.claims[0]
    assert claim["ledger_node_id"] == "ledger://chromatin/peak_annotation"
    assert claim["ledger_node_type"] == "analysis_run"
    assert claim["ledger_provenance"]
    assert claim["evidence_card_id"] == "chromatin.peak_annotation#evidence"
    card = claim["verification"]["evidence_card"]
    assert card["node_type"] == "evidence_card"
    assert card["provenance"]["block_id"] == block.id
    assert card["provenance"]["sources"] == ["chromatin_peak_annotation"]


@pytest.mark.parametrize("missing_field", ["node_type", "provenance"])
def test_public_compiler_withholds_untyped_or_unprovenanced_ledger_node(
    missing_field,
):
    node = {
        "node_id": "ledger://chromatin/peak_annotation",
        "node_type": "analysis_run",
        "provenance": {"kind": "structured_agent_result"},
        "status": "ran",
    }
    node.pop(missing_field)

    compiled = compile_public_claims(
        [_result_block()], {}, run_ledger={"entries": [node]}
    )

    assert compiled.claims == []
    assert compiled.summary()["n_withheld"] == 1


def test_active_verifier_rejects_a_result_claim_with_no_ledger_node():
    block = _result_block(analysis="unknown_active_family")
    block.metadata["claim"] = {"tier": "descriptive"}
    block.metadata["claim_verification"] = {
        "status": "supported",
        "evidence_card": {
            "evidence_card_id": f"{block.id}#evidence",
            "node_type": "evidence_card",
            "provenance": {"block_id": block.id, "sources": ["test"]},
        },
    }

    with pytest.raises(LedgerLinkageError, match="no ledger node"):
        verify_blocks_against_ledger([block], {"entries": []}, strict=True)


def _active_chromatin_findings() -> dict:
    return {
        "qc": {
            "status": "success", "data_type": "bulk_ATAC",
            "n_samples": 6, "mito_fraction": 0.04,
        },
        "peaks": {
            "status": "success", "data_type": "bulk_ATAC",
            "n_peaks": 100, "genome": "test_assembly",
        },
        "peak_counts": {
            "status": "success", "data_type": "bulk_ATAC",
            "counting_method": "interval overlap", "n_peaks": 100,
            "n_samples": 6,
        },
        "lsi": {
            "status": "success", "n_clusters": 3, "n_cells_used": 120,
            "n_peaks": 100, "n_components_used": 10,
            "dropped_components": [0],
        },
        "differential_accessibility": {
            "status": "success", "ran": True, "data_type": "bulk_ATAC",
            "method": "DESeq2", "n_comparisons_success": 1,
            "n_replicate_samples": 6, "n_peaks_tested": 100,
            "n_sig_total": 4, "padj_max": 0.05, "lfc_min": 1.0,
        },
        "motifs": {
            "status": "success", "ran": True, "method": "Fisher",
            "motif_source": {"collection": "versioned motifs", "n_motifs": 20},
            "per_group": {"group_a": {"n_enriched": 2}},
        },
        "peak_annotation": {
            "status": "success", "ran": True, "method": "nearest TSS",
            "gtf": "annotation.gtf", "promoter_upstream": 2000,
            "promoter_downstream": 500,
            "feature_distribution_overall": {
                "Promoter": 2, "Exonic": 1, "Intronic": 3,
                "Distal Intergenic": 4,
            },
        },
        "peak_ora": {
            "status": "success", "ran": True, "method": "nearest-TSS ORA",
            "organism": "test organism", "gtf": "annotation.gtf",
            "comparisons": [{
                "ora_method": "local hypergeometric",
                "pathways_summary": {"versioned_db": 2},
            }],
        },
        "regulatory": {
            "status": "success",
            "peak_to_gene": {
                "ran": True, "n_links": 8, "validation_level": "beta",
            },
        },
        "footprinting": {
            "status": "success", "ran": True,
            "group_a": "condition_a", "group_b": "condition_b",
            "differential_summary": {
                "parsed": True, "n_motifs_tested": 20,
                "n_significant": 2,
                "inference": {
                    "status": "success", "test": "welch_t",
                    "replicates_per_condition": {
                        "condition_a": 3, "condition_b": 3,
                    },
                    "null_label_controls": {"n_permutations": 20},
                },
            },
        },
    }


def test_every_active_chromatin_result_claim_has_typed_provenance_nodes():
    findings = _active_chromatin_findings()
    agent_result = {"findings": findings}
    blocks = ChromatinNarrator().collect("chromatin_agent", agent_result)
    result_blocks = [block for block in blocks if block.claim]
    ledger = build_run_ledger(
        {}, {"chromatin_agent": agent_result}
    )
    ensure_report_ledger_nodes(ledger, result_blocks)

    compiled = compile_public_claims(result_blocks, {}, run_ledger=ledger)

    assert compiled.summary()["n_withheld"] == 0
    assert len(compiled.claims) == len(result_blocks)
    assert {
        claim["analysis"] for claim in compiled.claims
    } >= {
        "peak_annotation", "peak_ora", "differential_tf_footprinting",
    }
    for claim in compiled.claims:
        assert claim["ledger_node_id"]
        assert claim["ledger_node_type"]
        assert claim["ledger_provenance"]
        assert claim["evidence_card_id"]
        assert claim["verification"]["evidence_card"]["node_type"]
        assert claim["verification"]["evidence_card"]["provenance"]


def test_rendered_results_resolve_to_existing_typed_nodes(tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    report = agent._render_html_report(
        experiment_id="c3_chromatin_e2e",
        exp_ctx={"organism": "test organism", "genome": "test_assembly"},
        intent={"summary": "C3 linkage E2E"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={
            "high": [], "medium": [], "low": [], "insufficient": [],
        },
        methods="methods",
        decisions=[],
        agent_results={
            "chromatin_agent": {"findings": _active_chromatin_findings()},
        },
        report_dir=tmp_path / "report",
    )

    methodology = json.loads(
        (report.parent / "methodology.json").read_text(encoding="utf-8")
    )
    assert methodology["claim_compilation"]["n_withheld"] == 0
    assert methodology["run_ledger"]["claim_ledger_verification"] == {
        "status": "ok", "n_violations": 0, "violations": [],
    }
    analyses = {claim["analysis"] for claim in methodology["claims"]}
    assert {
        "peak_annotation", "peak_ora", "differential_tf_footprinting",
    } <= analyses
    for claim in methodology["claims"]:
        assert claim["ledger_node_id"]
        assert claim["ledger_node_type"]
        assert claim["ledger_provenance"]
        assert claim["evidence_card_id"]
        card = claim["verification"]["evidence_card"]
        assert card["node_type"] == "evidence_card"
        assert card["provenance"]
