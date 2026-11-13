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
# Differential TF footprinting (ADR-078): Tn5-protected motif occupancy.
FOOTPRINT_ASSOCIATIVE = "footprint_associative"
# Paired cross-modal integration caveats (ADR-072 / V48 Slice I9).
PARALOG_MOTIF_AMBIGUITY = "paralog_motif_ambiguity"
DEPTH_CONFOUNDED = "depth_confounded"
MOFA_LATENT_ASSOCIATIVE = "mofa_latent_associative"
WNN_GUARDRAIL_DEGRADED = "wnn_guardrail_degraded"
# Non-same-cell cross-modal concordance caveat (ADR-080 / ADR-028 Slice 2).
CROSS_MODAL_UNPAIRED_COHORT = "cross_modal_unpaired_cohort"
# 3D genome focal chromatin loop caveat (ADR-082).
HIC_LOOP_ASSOCIATIVE = "hic_loop_associative"

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
    FOOTPRINT_ASSOCIATIVE: (
        "differential TF footprinting measures Tn5-protected motif occupancy, "
        "which is associative (paralog TFs share motifs) and not causal gene "
        "regulation"
    ),
    PARALOG_MOTIF_AMBIGUITY: (
        "candidate eRegulon shares motif sequence similarity and cis peaks with "
        "co-expressed paralog TFs — regulator identity is ambiguous without "
        "direct binding data"
    ),
    DEPTH_CONFOUNDED: (
        "TF-target correlation attenuates substantially after conditioning on "
        "sequencing depth and may reflect library-size covariation"
    ),
    MOFA_LATENT_ASSOCIATIVE: (
        "MOFA+ factor loading reflects undirected latent co-variation with "
        "arbitrary factor sign, not a directional perturbation effect"
    ),
    WNN_GUARDRAIL_DEGRADED: (
        "joint WNN embedding triggered anti-degradation guardrail warnings "
        "relative to single-modality structure (lens, not replacement)"
    ),
    CROSS_MODAL_UNPAIRED_COHORT: (
        "cross-modal concordance comes from non-same-cell assays (separate "
        "libraries or cohorts) and reflects condition- or donor-level "
        "co-variation, not single-cell cis regulation"
    ),
    HIC_LOOP_ASSOCIATIVE: (
        "3D chromatin loop contact enrichment indicates spatial proximity above "
        "local distance background, not proof of causal enhancer-promoter "
        "gene regulation"
    ),
}

# Recognition tokens for caveat codes whose free-text acknowledgement the
# adversarial gate must parse but that post-date its built-in token table. The
# hypothesis package registers these into ``devils_advocate._CONFOUND_TOKENS``
# at import time (additive, idempotent), so a hypothesis can acknowledge the
# caveat by code, by its gloss, or by its own phrasing.
EXTRA_CONFOUND_TOKENS: dict[str, tuple[str, ...]] = {
    CROSS_MODAL_UNPAIRED_COHORT: (
        "cross_modal_unpaired_cohort",
        "unpaired cohort",
        "unpaired assay",
        "non-same-cell",
        "not same-cell",
        "not from the same cells",
        "separate libraries",
        "separate cohorts",
        "condition-level co-variation",
        "donor-level co-variation",
        "condition- or donor-level",
        "not single-cell cis",
    ),
    HIC_LOOP_ASSOCIATIVE: (
        "hic_loop_associative",
        "spatial proximity",
        "physical proximity",
        "contact frequency",
        "3d loop",
        "chromatin loop",
        "loop contact",
        "not causal regulation",
        "proximity alone",
    ),
}

CAVEAT_CODES = frozenset(CAVEAT_GLOSSES)


def caveat_gloss(code: str) -> str:
    """Human-readable gloss for a caveat code (the code itself if unknown)."""
    return CAVEAT_GLOSSES.get(str(code), str(code))


def is_caveat_code(value) -> bool:
    """True iff ``value`` is a known structured caveat code."""
    return isinstance(value, str) and value in CAVEAT_CODES
