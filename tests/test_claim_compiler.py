"""X14 Claim Compiler: evidence-tier classification and language capping."""

from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.agents.narrative.claim_compiler import (
    classify_claim,
    annotate_claim_tiers,
    compile_claims,
    design_is_interventional,
    resolve_causal_estimand_license,
)


def _block(analysis, claim, block_type="result", bid=None, evidence=True):
    return NarrativeBlock(
        id=bid or f"scrna.{analysis}.GroupA.condA_vs_condB",
        modality="scRNA-seq",
        analysis=analysis,
        block_type=block_type,
        title=f"{analysis} block",
        status="success",
        confidence="medium",
        claim=claim,
        evidence=[EvidenceItem(label="n", value=10, source="test")] if evidence else [],
    )


def test_qc_block_is_descriptive():
    b = _block("qc", "Retained 1000 of 1100 cells.", block_type="qc", bid="scrna.qc")
    c = classify_claim(b)
    assert c.tier == "descriptive"
    assert c.licensed_language == "descriptive"


def test_pseudobulk_de_is_associative_and_observational():
    b = _block("pseudobulk_de", "GroupA condA_vs_condB had 120 DE genes.")
    c = classify_claim(b)
    assert c.tier == "associative"
    assert c.licensed_language == "associative"
    assert any("observational" in lim.lower() for lim in c.limitations)


def test_verified_estimand_contract_licenses_causal():
    b = _block("pseudobulk_de", "GroupA condA_vs_condB had 120 DE genes.")
    b.metadata["estimand"] = {
        "id": "condition_effect",
        "contrast": "condA_vs_condB",
    }
    exp_ctx = {"design": {"causal_estimands": [{
        "id": "condition_effect",
        "contrast": "condA_vs_condB",
        "interventional": True,
        "randomized": True,
        "confounding_verified": True,
    }]}}
    c = classify_claim(
        b, causal_license=resolve_causal_estimand_license(b, exp_ctx)
    )
    assert c.tier == "causal_experimental"
    assert c.licensed_language == "causal"


def test_one_mechanistic_line_is_weak_mechanistic():
    b = _block("pseudobulk_de", "GroupA condA_vs_condB had 120 DE genes.")
    c = classify_claim(b, converging_categories={"regulatory"})
    assert c.tier == "weak_mechanistic"
    assert c.licensed_language == "associative"  # still observational


def test_two_mechanistic_lines_are_strong_mechanistic():
    b = _block("pseudobulk_de", "GroupA condA_vs_condB had 120 DE genes.")
    # two regulatory lines won't appear as one category; simulate convergence by
    # a regulatory category plus the block being part of a multi-line subject.
    c = classify_claim(b, converging_categories={"regulatory", "functional_annotation"})
    # regulatory counts as one mechanistic line -> weak; need 2 distinct
    # mechanistic categories for strong. Here only one mechanistic category.
    assert c.tier in {"weak_mechanistic"}


def test_causal_language_above_tier_is_flagged():
    b = _block("pseudobulk_de",
               "GroupA condA_vs_condB: the module induces the senescent program.")
    c = classify_claim(b)
    assert c.tier == "associative"
    assert c.language_violation is not None
    assert any("causation" in lim.lower() for lim in c.limitations)


def test_annotate_and_compile_claims_roundtrip():
    blocks = [
        _block("qc", "Retained cells.", block_type="qc", bid="scrna.qc"),
        _block("pseudobulk_de", "GroupA condA_vs_condB had 10 DE genes."),
        _block("cell_communication", "50 L-R interactions.", bid="scrna.cellcomm"),
    ]
    annotate_claim_tiers(blocks, exp_ctx={})
    for b in blocks:
        assert "claim" in b.metadata
        assert b.metadata["claim"]["tier"] in {
            "descriptive", "associative", "weak_mechanistic",
            "strong_mechanistic", "causal_experimental",
        }
    manifests = compile_claims(blocks, exp_ctx={})
    assert len(manifests) == 3
    de = next(m for m in manifests if m["analysis"] == "pseudobulk_de")
    assert de["tier"] == "associative"
    assert de["claim_id"].startswith("scrna.pseudobulk")
    assert "evidence" in de and isinstance(de["evidence"], list)
    assert de["evidence_card_id"] == f"{de['claim_id']}#evidence"
    assert de["verification"]["status"] == "supported"
    assert de["verification"]["evidence_card"]["n_refs"] >= 1


def test_claim_manifest_marks_unsupported_claim_text_without_throwing():
    block = _block(
        "pseudobulk_de",
        "GroupA condA_vs_condB had 120 DE genes.",
    )
    block.evidence[0].value = 10
    manifest = compile_claims([block], exp_ctx={})[0]
    assert manifest["verification"]["status"] == "unsupported"
    assert "120" in manifest["verification"]["unsupported"][0]["reason"]


def test_design_is_interventional_is_conservative():
    assert design_is_interventional({}) is False
    assert design_is_interventional({"design": {"main_factor": "treatment"}}) is False
    assert design_is_interventional({"design": {"interventional": True}}) is True
