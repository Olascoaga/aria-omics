"""C2 (audit 2026-06-11): the W-CLAIM numeric gate must decide correctly WHICH
number to verify, in both directions.

(b) False NEGATIVE — effect sizes below |2| were blanket-exempt, so a fabricated
    correlation / log2FC / NES / odds ratio sailed through. They must now be
    verified against the evidence card.
(a) False POSITIVE — an integer that is part of a group LABEL ("cluster 3") was
    required on the card, which aborted the report under strict verification.
    Group-label and inequality-threshold numbers must stay exempt.

These are failing-first guards: before the fix, the effect-size case does NOT
raise (the bug). The label/threshold cases are the regression fence.
"""

import pytest

from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.agents.narrative.evidence_verifier import (
    EvidenceVerificationError,
    verify_block_claim_support,
)


def _block(claim, evidence):
    return NarrativeBlock(
        id="x.synth", modality="bulk RNA-seq", analysis="contrast",
        block_type="result", title="Synthetic block", status="success",
        confidence="medium", claim=claim, evidence=evidence,
    )


# ── (b) effect sizes must be verified, even below |2| ────────────────────────

def test_unsupported_correlation_below_two_is_flagged():
    """A correlation of 0.98 absent from the evidence card must raise."""
    block = _block(
        "The two contrasts show a correlation of 0.98 across genes.",
        [EvidenceItem(label="DE genes", value=42, source="contrast")],
    )
    with pytest.raises(EvidenceVerificationError):
        verify_block_claim_support(block, strict=True)


def test_unsupported_log2fc_below_two_is_flagged():
    """A reported log2FC of 1.7 not on the card must raise."""
    block = _block(
        "The marker shifted with a log2 fold change of 1.7.",
        [EvidenceItem(label="DE genes", value=12, source="contrast")],
    )
    with pytest.raises(EvidenceVerificationError):
        verify_block_claim_support(block, strict=True)


def test_supported_correlation_below_two_passes():
    """The same 0.98 IS supported when it is on the evidence card."""
    block = _block(
        "The two contrasts show a correlation of 0.98 across genes.",
        [
            EvidenceItem(label="DE genes", value=42, source="contrast"),
            EvidenceItem(label="contrast correlation", value=0.98,
                         source="contrast"),
        ],
    )
    manifest = verify_block_claim_support(block, strict=True)
    assert manifest["status"] == "supported"


# ── (a) group-label and threshold numbers stay exempt (regression fence) ─────

def test_group_label_integer_is_not_required_on_card():
    """'cluster 3' is a label, not a measurement: it must not raise even though
    3 is absent from the evidence card (the production crash this fixes).

    analysis is 'pathway_enrichment' so the pathway claim matches the block
    family — the only number under test is the group-label '3'."""
    block = NarrativeBlock(
        id="bulk.pathway.x", modality="bulk RNA-seq",
        analysis="pathway_enrichment", block_type="result",
        title="Pathway enrichment cluster 3", status="success",
        confidence="medium",
        claim="Pathway enrichment found 2 term(s) for cluster 3.",
        evidence=[EvidenceItem(label="enriched terms", value=2,
                               source="bulk_pathways")],
        metrics={"n_terms": 2},
    )
    manifest = verify_block_claim_support(block, strict=True)
    assert manifest["status"] == "supported"


def test_inequality_threshold_value_is_not_required_on_card():
    """A threshold like 'FDR < 0.25' is a decision rule, not a measured claim."""
    block = _block(
        "Signals were reported at FDR < 0.25 for this contrast.",
        [EvidenceItem(label="DE genes", value=7, source="contrast")],
    )
    manifest = verify_block_claim_support(block, strict=True)
    assert manifest["status"] == "supported"


# ── T11 (tri-audit 2026-06-14): mixed-case (mouse/rat) gene symbols must be
# verified too. The all-caps-only entity pattern silently skipped Sox2/Olig2, so a
# fabricated mixed-case gene in prose was never checked against the evidence card.

def test_unsupported_mixed_case_mouse_gene_is_flagged():
    block = _block(
        "The population was marked by Sox2 and Olig2 expression.",
        [EvidenceItem(label="DE genes", value=12, source="contrast")],
    )
    with pytest.raises(EvidenceVerificationError):
        verify_block_claim_support(block, strict=True)


def test_supported_mixed_case_mouse_gene_passes():
    block = _block(
        "The population was marked by Sox2 expression.",
        [
            EvidenceItem(label="DE genes", value=12, source="contrast"),
            EvidenceItem(label="top gene Sox2", value=2.1, source="contrast"),
        ],
    )
    manifest = verify_block_claim_support(block, strict=True)
    assert manifest["status"] == "supported"


def test_deseq2_tool_token_is_not_treated_as_gene():
    """A digit-bearing tool token (DESeq2) must not be required on the card as if
    it were a gene — it is a stopword."""
    from aria.agents.narrative.evidence_verifier import _claim_entities

    ents = _claim_entities("Counts were modeled with DESeq2 and Sox2 was tested.")
    assert "Sox2" in ents
    assert "DESeq2" not in ents
