"""Structured caveat codes for the SPECULATIVE tier (ADR-057 rail #5; round-3 H16).

Codex blocker A1: the adversarial gate (``devils_advocate.py``) only recognised
five TECHNICAL confound categories (batch / composition / ambient / doublet / low
replication). The per-modality adapters also inherit SCIENTIFIC caveats — an
enriched motif is not TF binding, a gene-activity score is a moderate chromatin
proxy, a peak-to-gene link is associative — but those were free prose the gate did
not match, so a hypothesis could lean on them without ever owning them.

This module makes every caveat a canonical CODE. ``EvidenceSignal.caveats_inherited``
now carries codes (the structured source of truth the gate enforces), and the
human-readable GLOSS is rendered at the boundaries the human/LLM see (the proposer
prompt and the report). A hypothesis must acknowledge EVERY code on the evidence it
uses, scientific caveats included.

Pure data: no biology, no LLM. These are fixed methodological/epistemic categories
(ADR-011's gene-name policy does not apply), so the list is a closed, auditable set.
"""

from __future__ import annotations

# Technical confounds (QC-level).
BATCH = "batch"
COMPOSITION = "composition"
AMBIENT = "ambient"
DOUBLET = "doublet"
LOW_REPLICATION = "low_replication"
# Scientific / interpretive caveats (assay-level associativity limits).
MOTIF_NOT_BINDING = "motif_not_binding"
GENE_ACTIVITY_PROXY = "gene_activity_proxy"
PEAK2GENE_ASSOCIATIVE = "peak2gene_associative"

# code -> human-readable gloss (rendered in the prompt + the report).
CAVEAT_GLOSSES: dict[str, str] = {
    BATCH: "residual batch effect — condition is partially confounded with batch",
    COMPOSITION: "cell-type composition shift may drive the apparent signal",
    AMBIENT: "ambient RNA contamination was not corrected",
    DOUBLET: "doublets were not ruled out",
    LOW_REPLICATION: (
        "low replication / underpowered — effects are directional, not "
        "FDR-calibrated"
    ),
    MOTIF_NOT_BINDING: (
        "an enriched motif is an associative database match, not evidence of TF "
        "binding or activity"
    ),
    GENE_ACTIVITY_PROXY: (
        "gene-activity score is a moderate chromatin-accessibility proxy "
        "(~0.51 external concordance), not RNA expression"
    ),
    PEAK2GENE_ASSOCIATIVE: (
        "a peak-to-gene link is an associative cross-cell correlation, not an "
        "established regulatory mechanism"
    ),
}

CAVEAT_CODES = frozenset(CAVEAT_GLOSSES)


def caveat_gloss(code: str) -> str:
    """Human-readable gloss for a caveat code (the code itself if unknown)."""
    return CAVEAT_GLOSSES.get(str(code), str(code))


def is_caveat_code(value) -> bool:
    """True iff ``value`` is a known structured caveat code."""
    return isinstance(value, str) and value in CAVEAT_CODES
