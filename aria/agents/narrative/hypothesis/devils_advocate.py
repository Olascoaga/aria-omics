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

from .caveats import CAVEAT_CODES, is_caveat_code
from .gates import GateResult
from .grounding import _canon
from .types import EvidenceSignal, Hypothesis

# Recognition tokens per caveat CODE (H16). The inherited caveats on the evidence
# are now structured codes (caveats.py), so visible_confounds reads them directly.
# These tokens exist only to parse the model's free-text DECLARED confounds back
# to a code — it may echo a code, the code's gloss, or its own phrasing. Matching
# is intentionally PERMISSIVE substring (unlike the word-boundary lints in
# gates.py): the goal is to never MISS an acknowledged confound. Now covers the
# SCIENTIFIC caveats Codex flagged (motif != binding, gene-activity proxy,
# peak-to-gene associative), not only the five technical ones.
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
        "low_replication",
        "low replicate",
        "low replication",
        "low power",
        "low-power",
        "n=2",
        "n = 2",
        "single-sample",
        "single sample",
        "underpowered",
        "not fdr",
        "directional, not",
        "directional not",
    ),
    "motif_not_binding": (
        "motif_not_binding",
        "not binding",
        "not evidence of tf binding",
        "not evidence of binding",
        "associative database match",
        "not direct binding",
        "does not imply binding",
        "binding or activity",
        "not tf binding",
        "motif is associative",
        "not actual binding",
    ),
    "gene_activity_proxy": (
        "gene_activity_proxy",
        "gene activity",
        "gene-activity",
        "accessibility proxy",
        "chromatin proxy",
        "not rna expression",
        "moderate",
    ),
    "peak2gene_associative": (
        "peak2gene_associative",
        "peak-to-gene",
        "peak to gene",
        "peak2gene",
        "associative correlation",
        "not regulatory mechanism",
        "not establish regulatory",
        "associative cross-cell",
    ),
    "footprint_associative": (
        "footprint_associative",
        "footprint",
        "occupancy",
        "tn5",
        "protected motif",
        "not causal regulation",
        "associative occupancy",
        "paralog",
    ),
    "paralog_motif_ambiguity": (
        "paralog_motif_ambiguity",
        "paralog",
        "motif family",
        "motif-family",
        "shared motif",
        "paralogous",
        "regulator identity is ambiguous",
    ),
    "depth_confounded": (
        "depth_confounded",
        "depth confound",
        "depth-confounded",
        "sequencing depth",
        "library depth",
        "library-size",
        "library size",
        "total counts",
    ),
    "mofa_latent_associative": (
        "mofa_latent_associative",
        "mofa",
        "latent factor",
        "latent co-variation",
        "latent covariation",
        "factor sign",
        "undirected latent",
    ),
    "wnn_guardrail_degraded": (
        "wnn_guardrail_degraded",
        "wnn guardrail",
        "anti-degradation",
        "guardrail warning",
        "guardrail degraded",
        "joint embedding warning",
        "lens, not replacement",
    ),
    "cross_modal_unpaired_cohort": (
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
    "hic_loop_associative": (
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


def _categories(text: str) -> set[str]:
    """Map a free-text declared confound (or a bare code) to caveat codes."""
    t = str(text or "").lower()
    cats = {
        cat
        for cat, toks in _CONFOUND_TOKENS.items()
        if any(tok in t for tok in toks)
    }
    # Also accept a bare code (e.g. the model declares "batch" or "motif_not_binding").
    norm = t.strip()
    if norm in _CONFOUND_TOKENS:
        cats.add(norm)
    return cats


def _extract_caveat_codes(sig: Any) -> set[str]:
    cats: set[str] = set()
    for caveat in getattr(sig, "caveats_inherited", None) or []:
        if is_caveat_code(caveat):
            cats.add(caveat)
        else:
            cats |= _categories(caveat)
    return cats


def visible_confounds(
    hyp: Hypothesis,
    evidence_index: Mapping[str, EvidenceSignal | list[EvidenceSignal]],
    signals_by_id: Mapping[str, EvidenceSignal] | None = None,
) -> set[str]:
    """Caveat codes the audited evidence the hypothesis cites already flags.

    Accepts either an entity->signal index or an entity->[signals] index (H4): a
    caveat flagged on the entity in ANY context it was measured in must be owned,
    so all of an entity's context-distinct signals are unioned. Keys and queried
    entities are canonicalised via ``_canon`` so formatting/Greek variants
    (``NKX2_1`` <-> ``NKX2-1``, ``REV-ERBa`` <-> ``REV-ERBα``) cannot bypass the
    gate. Caveats from signals cited in ``hyp.observed_claims`` (by ``signal_id``)
    and node-wide caveats on ``hyp.observation_refs`` are also unioned. H16: the
    inherited caveats are structured CODES, read directly; a legacy free-text
    caveat is still mapped through ``_categories`` so external callers do not
    break.
    """
    canon_index: dict[str, list[EvidenceSignal]] = {}
    id_index: dict[str, EvidenceSignal] = {}
    node_signals: dict[str, list[EvidenceSignal]] = {}

    def _register(sig: Any, key_hint: str | None = None) -> None:
        if sig is None:
            return
        if key_hint:
            ckey = _canon(key_hint) or str(key_hint).strip().lower()
            if ckey:
                canon_index.setdefault(ckey, []).append(sig)
        ent = getattr(sig, "entity", None)
        if ent:
            ec = _canon(ent) or str(ent).strip().lower()
            if ec and ec != (_canon(key_hint) if key_hint else ""):
                canon_index.setdefault(ec, []).append(sig)
        sid = getattr(sig, "signal_id", None)
        if sid:
            id_index.setdefault(str(sid), sig)
        nref = getattr(sig, "audited_node_ref", None)
        if nref:
            node_signals.setdefault(str(nref), []).append(sig)

    for key, value in (evidence_index or {}).items():
        if value is None:
            continue
        sigs = value if isinstance(value, list) else [value]
        for sig in sigs:
            _register(sig, key_hint=str(key))

    for sid, sig in (signals_by_id or {}).items():
        if sig is not None:
            if sid:
                id_index.setdefault(str(sid), sig)
            _register(sig)

    cats: set[str] = set()
    for ent in hyp.entities or []:
        ckey = _canon(ent) or str(ent).strip().lower()
        for sig in canon_index.get(ckey) or []:
            cats |= _extract_caveat_codes(sig)

    for claim in hyp.observed_claims or []:
        if not isinstance(claim, dict):
            continue
        sid = str(claim.get("signal_id") or "").strip()
        if sid and sid in id_index:
            cats |= _extract_caveat_codes(id_index[sid])

    for ref in hyp.observation_refs or []:
        sigs = node_signals.get(str(ref))
        if not sigs:
            continue
        # Union node-wide caveats shared by every audited signal from this node
        # so citing a node via an alias entity still inherits its mandatory
        # modality/node caveats without leaking unrelated entity-specific flags.
        common = set.intersection(*(_extract_caveat_codes(s) for s in sigs))
        cats |= common

    return cats


def inherited_caveat_codes(
    hyp: Hypothesis,
    evidence_index: Mapping[str, EvidenceSignal | list[EvidenceSignal]],
    signals_by_id: Mapping[str, EvidenceSignal] | None = None,
) -> list[str]:
    """Sorted list of the caveat codes inherited by a hypothesis (for rendering)."""
    return sorted(
        c
        for c in visible_confounds(hyp, evidence_index, signals_by_id=signals_by_id)
        if c in CAVEAT_CODES
    )


def declared_confounds(hyp: Hypothesis) -> set[str]:
    """Confound categories the hypothesis explicitly acknowledges."""
    da = hyp.devils_advocate or {}
    cats: set[str] = set()
    for item in da.get("confounds") or []:
        cats |= _categories(item)
    return cats


def check_devils_advocate(
    hyp: Hypothesis,
    evidence_index: Mapping[str, EvidenceSignal | list[EvidenceSignal]],
    signals_by_id: Mapping[str, EvidenceSignal] | None = None,
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
        visible_confounds(hyp, evidence_index, signals_by_id=signals_by_id)
        - declared_confounds(hyp)
    )
    if unacknowledged:
        return GateResult(
            "devils_advocate",
            False,
            "confounds flagged by the cited audited evidence are not "
            f"acknowledged: {unacknowledged}",
        )
    return GateResult("devils_advocate", True)
