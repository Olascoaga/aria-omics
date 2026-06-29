"""Adversarial gate for the SPECULATIVE tier (ADR-057 rail #5).

The hypothesis-layer analog of ``aria/agents/narrative/devils_advocate.py``:
even a speculation must declare its alternatives. That module runs a
deterministic, single-shot challenge over report blocks using a fixed set of
TECHNICAL confounders (batch / composition / ambient RNA / doublets / low
replication — methodological vocabulary, not biological content; ADR-011 does
not apply to fixed QC/technical terms). This module reuses the same philosophy
one level up, over a ``Hypothesis`` and the audited ``EvidenceSignal`` it cites.

Each hypothesis must (a) offer a simpler/competing explanation and (b) own every
confound the AUDITED evidence it cites already flags (carried in
``EvidenceSignal.caveats_inherited``, populated by the per-modality adapters in
S5-S8). A hypothesis that hides a known confound, or offers no competing
explanation, has not earned publication and is REJECTED — the single-clean-
narrative failure mode ARIA exists to prevent. The formal parsimony ranking is
S9; this gate is the binary adversarial check.
"""

from __future__ import annotations

from typing import Mapping

from .gates import GateResult
from .types import EvidenceSignal, Hypothesis

# Technical confounder categories — mirror the alternatives enumerated in
# narrative/devils_advocate.py. Methodological vocabulary only.
_CONFOUND_TOKENS: dict[str, tuple[str, ...]] = {
    "batch": ("batch",),
    "composition": (
        "composition",
        "compositional",
        "proportion",
        "cell-type shift",
        "celltype shift",
        "cell type shift",
    ),
    "ambient": ("ambient", "soupx", "decontx", "contamination"),
    "doublet": ("doublet", "scrublet"),
    "low_replication": (
        "low replicate",
        "low replication",
        "low power",
        "low-power",
        "n=2",
        "n = 2",
        "single-sample",
        "single sample",
        "underpowered",
    ),
}


def _categories(text: str) -> set[str]:
    t = str(text or "").lower()
    cats = {
        cat
        for cat, toks in _CONFOUND_TOKENS.items()
        if any(tok in t for tok in toks)
    }
    # Also accept a bare category name (e.g. a declared confound "batch").
    norm = t.strip()
    if norm in _CONFOUND_TOKENS:
        cats.add(norm)
    return cats


def visible_confounds(
    hyp: Hypothesis,
    evidence_index: Mapping[str, EvidenceSignal | list[EvidenceSignal]],
) -> set[str]:
    """Confound categories the audited evidence the hypothesis cites already flags.

    Accepts either an entity->signal index or an entity->[signals] index (H4): a
    confound flagged on the entity in ANY context it was measured in must be
    owned, so all of an entity's context-distinct signals are unioned.
    """
    cats: set[str] = set()
    for ent in hyp.entities or []:
        value = evidence_index.get(str(ent).strip().lower())
        if value is None:
            continue
        sigs = value if isinstance(value, list) else [value]
        for sig in sigs:
            for caveat in getattr(sig, "caveats_inherited", None) or []:
                cats |= _categories(caveat)
    return cats


def declared_confounds(hyp: Hypothesis) -> set[str]:
    """Confound categories the hypothesis explicitly acknowledges."""
    da = hyp.devils_advocate or {}
    cats: set[str] = set()
    for item in da.get("confounds") or []:
        cats |= _categories(item)
    return cats


def check_devils_advocate(
    hyp: Hypothesis, evidence_index: Mapping[str, EvidenceSignal]
) -> GateResult:
    """Reject a hypothesis that offers no alternative or hides a known confound."""
    da = hyp.devils_advocate or {}
    simpler = str(da.get("simpler_explanation") or "").strip()
    if not simpler:
        return GateResult(
            "devils_advocate",
            False,
            "no simpler/competing explanation declared; a single clean "
            "narrative is not publishable in this tier",
        )
    unacknowledged = sorted(
        visible_confounds(hyp, evidence_index) - declared_confounds(hyp)
    )
    if unacknowledged:
        return GateResult(
            "devils_advocate",
            False,
            "confounds flagged by the cited audited evidence are not "
            f"acknowledged: {unacknowledged}",
        )
    return GateResult("devils_advocate", True)
