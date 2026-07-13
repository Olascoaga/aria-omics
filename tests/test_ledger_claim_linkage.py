"""W-LEDGER (Tier-A, pre-4.6 polish): close the loop so every report claim is
linked both to an evidence card (W-CLAIM) AND to a node of the run ledger, with
active verification.

The run ledger already reconciles planned-vs-executed analyses. W-LEDGER gives
each ledger entry a stable ``node_id`` and links every claim to its node, then
actively verifies that a claim of associative-or-stronger tier never cites a
ledger node the run marked not-run/skipped/error (the thin-report contradiction:
a downstream DE claim that survived while DE was silently skipped). Pure
bookkeeping; no LLM, no biology.
"""

import pytest

from aria.agents.narrative.types import NarrativeBlock, EvidenceItem
from aria.agents.narrative.evidence_verifier import build_evidence_card
from aria.agents.narrative.run_ledger import (
    build_run_ledger,
    node_id_for,
    link_claims_to_ledger,
    verify_blocks_against_ledger,
    LedgerLinkageError,
)


def _ran_results():
    return {"scrna_agent": {"findings": {
        "qc": {"status": "success", "n_cells": 100},
        "pseudobulk_de": {"per_group": {"A": {}}},
    }}}


def _skipped_de_results():
    # The thin-report gap: DE was planned + skipped, yet a DE claim survives.
    return {"scrna_agent": {"findings": {
        "qc": {"status": "success", "n_cells": 100},
        "pseudobulk_de": {"status": "skipped", "reason": "no_confirmed_contrast"},
    }}}


_PB_PLAN = {"design_intelligence": {
    "recommended": ["Donor-level pseudobulk DESeq2 between conditions."],
    "optional": []}}


# ── node ids ─────────────────────────────────────────────────────────────────

def test_ledger_entries_carry_stable_node_ids():
    ledger = build_run_ledger(_PB_PLAN, _ran_results())
    nids = {e["node_id"] for e in ledger["entries"]}
    assert "ledger://scRNA/pseudobulk_de" in nids
    assert "ledger://scRNA/qc" in nids
    # Deterministic and consistent with the public helper.
    assert node_id_for("scRNA-seq", "pseudobulk_de") == "ledger://scRNA/pseudobulk_de"
    assert node_id_for("scRNA-seq", "sample_qc") == "ledger://scRNA/qc"


# ── claim ↔ node linkage ─────────────────────────────────────────────────────

def test_claim_links_to_ran_node():
    ledger = build_run_ledger(_PB_PLAN, _ran_results())
    claims = [{"claim_id": "scrna.pb.A", "modality": "scRNA-seq",
               "analysis": "pseudobulk_de", "tier": "associative",
               "evidence_card_id": "scrna.pb.A#evidence"}]
    summary = link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_node_id"] == "ledger://scRNA/pseudobulk_de"
    assert claims[0]["ledger_status"] == "ran"
    assert claims[0]["ledger_linked"] is True
    assert summary["n_violations"] == 0
    assert summary["linked"] == 1


def test_claim_citing_not_run_analysis_is_a_violation():
    ledger = build_run_ledger(_PB_PLAN, _skipped_de_results())
    claims = [{"claim_id": "scrna.pb.A", "modality": "scRNA-seq",
               "analysis": "pseudobulk_de", "tier": "associative",
               "evidence_card_id": "scrna.pb.A#evidence"}]
    summary = link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_status"] == "skipped"
    assert summary["n_violations"] == 1
    assert summary["violations"][0]["node_id"] == "ledger://scRNA/pseudobulk_de"
    assert summary["violations"][0]["ledger_status"] == "skipped"


def test_no_ledger_node_is_a_publication_violation():
    # C3: an absent node is recorded honestly and fails public eligibility.
    ledger = build_run_ledger(
        {}, {"scrna_agent": {"findings": {"qc": {"status": "success"}}}})
    claims = [{"claim_id": "scrna.markers", "modality": "scRNA-seq",
               "analysis": "differential_expression", "tier": "associative",
               "evidence_card_id": "scrna.markers#evidence"}]
    summary = link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_linked"] is False
    assert claims[0]["ledger_status"] == "no_ledger_node"
    assert claims[0]["ledger_node_id"] is None
    assert summary["n_violations"] == 1


def test_descriptive_result_claim_to_not_run_node_is_a_violation():
    ledger = build_run_ledger(_PB_PLAN, _skipped_de_results())
    claims = [{"claim_id": "scrna.pb.A", "modality": "scRNA-seq",
               "analysis": "pseudobulk_de", "tier": "descriptive",
               "evidence_card_id": "scrna.pb.A#evidence"}]
    summary = link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_status"] == "skipped"
    assert summary["n_violations"] == 1


# ── bulk coverage ────────────────────────────────────────────────────────────

def test_bulk_claims_link_to_bulk_nodes():
    exp_ctx = {"design_intelligence": {
        "recommended": ["DESeq2 differential expression between conditions."],
        "optional": []}}
    agent_results = {"bulk_rna_agent": {"findings": {
        "sample_qc": {"n_samples": 6},
        "contrasts": [{"name": "t_vs_c", "status": "success",
                       "pathways": {"GO_BP": [{"term": "x"}]}}],
    }}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    assert "bulk" in ledger["modalities"]
    nids = {e["node_id"] for e in ledger["entries"]}
    assert "ledger://bulk/differential_expression" in nids
    claims = [{"claim_id": "bulk.de.t_vs_c", "modality": "bulk RNA-seq",
               "analysis": "differential_expression", "tier": "associative"}]
    link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_node_id"] == "ledger://bulk/differential_expression"
    assert claims[0]["ledger_status"] == "ran"


def test_bulk_pathway_node_ran_for_gsea_only_running_sums():
    # Regression: a bulk contrast with GSEA running-sum plots but NO gsea_table
    # and NO ORA pathways still builds a GSEA block in the narrator, so the
    # pathway_enrichment ledger node must read "ran" (not a false not-run that
    # would flag the associative GSEA claim as a violation).
    agent_results = {"bulk_rna_agent": {"findings": {
        "sample_qc": {"n_samples": 6},
        "contrasts": [{"name": "t_vs_c", "status": "success",
                       "plots": {"gsea_running_sums": ["/p1.png", "/p2.png"]}}],
    }}}
    ledger = build_run_ledger({}, agent_results)
    node = next(e for e in ledger["entries"]
                if e["node_id"] == "ledger://bulk/pathway_enrichment")
    assert node["status"] == "ran"
    # And a GSEA-only claim links to that ran node without a violation.
    claims = [{"claim_id": "bulk.gsea.t_vs_c", "modality": "bulk RNA-seq",
               "analysis": "gsea_preranked", "tier": "associative",
               "evidence_card_id": "bulk.gsea.t_vs_c#evidence"}]
    summary = link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_status"] == "ran"
    assert summary["n_violations"] == 0


def test_bulk_pathway_claim_links_when_pathways_present():
    exp_ctx = {"design_intelligence": {
        "recommended": ["Pathway/ORA enrichment."], "optional": []}}
    agent_results = {"bulk_rna_agent": {"findings": {
        "sample_qc": {"n_samples": 6},
        "contrasts": [{"name": "t_vs_c", "status": "success",
                       "pathways": {"GO_BP": [{"term": "x"}]}}],
    }}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    claims = [{"claim_id": "bulk.pathway.t_vs_c", "modality": "bulk RNA-seq",
               "analysis": "pathway_enrichment", "tier": "associative",
               "evidence_card_id": "bulk.pathway.t_vs_c#evidence"}]
    summary = link_claims_to_ledger(claims, ledger)
    assert claims[0]["ledger_node_id"] == "ledger://bulk/pathway_enrichment"
    assert claims[0]["ledger_status"] == "ran"
    assert summary["n_violations"] == 0


# ── active verification over rendered blocks ─────────────────────────────────

def _pb_block(tier="associative", analysis="pseudobulk_de"):
    b = NarrativeBlock(
        id="scrna.pseudobulk.A", modality="scRNA-seq", analysis=analysis,
        block_type="result", title="DE", status="success", confidence="medium",
        claim="A shows differential expression between treat and ctrl.",
        evidence=[EvidenceItem(label="n_sig", value=42)],
    )
    b.metadata["claim"] = {"tier": tier}
    b.metadata["claim_verification"] = {
        "status": "supported",
        "evidence_card": build_evidence_card(b).as_dict(),
    }
    return b


def test_active_verification_raises_on_block_claiming_not_run():
    ledger = build_run_ledger(_PB_PLAN, _skipped_de_results())
    with pytest.raises(LedgerLinkageError):
        verify_blocks_against_ledger([_pb_block()], ledger, strict=True)


def test_active_verification_passes_when_node_ran():
    ledger = build_run_ledger(_PB_PLAN, _ran_results())
    manifest = verify_blocks_against_ledger([_pb_block()], ledger, strict=True)
    assert manifest["status"] == "ok"
    assert manifest["n_violations"] == 0


def test_active_verification_non_strict_collects_without_raising():
    ledger = build_run_ledger(_PB_PLAN, _skipped_de_results())
    manifest = verify_blocks_against_ledger([_pb_block()], ledger, strict=False)
    assert manifest["status"] == "violated"
    assert manifest["n_violations"] == 1


def test_descriptive_result_block_raises_for_not_run_node():
    ledger = build_run_ledger(_PB_PLAN, _skipped_de_results())
    block = _pb_block(tier="descriptive")
    with pytest.raises(LedgerLinkageError):
        verify_blocks_against_ledger([block], ledger, strict=True)
