"""Competitive ranking for the SPECULATIVE tier (ADR-057 rail #4).

A ranked set of rival hypotheses, NOT one clean narrative. The ranking is
deliberately evidence-driven and penalises confidence: a hypothesis ranks higher
when more INDEPENDENT audited observations converge on it and the cited effects
are stronger, and lower when the evidence it leans on carries more inherited
confounds. Narrative elegance is never a factor.

``Hypothesis.rank_evidence`` is a HEURISTIC ORDERING KEY, not a scientific
confidence score. It says which of two rival speculations rests on more
independent audited lines — it never asserts how likely either is to be true.

H5 hardening (audit Codex issues 5 + 6):
  - Effects are NORMALISED per measure type before they are combined: raw
    log2FC, correlation, enrichment and gene-activity live on different scales
    and must not be averaged directly. Each is mapped to a comparable [0, 1]
    magnitude.
  - The set is made genuinely competitive: exact-duplicate hypotheses collapse,
    and every hypothesis records the rivals it ``competing_with`` (others that
    explain overlapping evidence differently).

H18 hardening (audit Codex issue 5/M1):
  - Independence, effect and confound load are scored over the exact audited
    signals the hypothesis CITES in ``observed_claims`` (resolved by
    ``signal_id``) — NOT over every signal the named entities were ever measured
    in. Naming an entity broadly can no longer inherit independence the
    hypothesis never read: a hypothesis that cites one contrast counts as one
    line even if the gene was measured in three. The cited signals are already
    grounded (H15 rejects unknown/misattributed/contradicting citations before a
    hypothesis reaches the ranking), so this cannot be padded either.

The per-hypothesis basis is recorded in ``Hypothesis.rank_evidence`` so the
report can show why one outranks another.
"""

from __future__ import annotations

import math

from .types import EvidenceSignal, Hypothesis


def _norm(entity: str) -> str:
    return str(entity or "").strip().lower()


def _num(value) -> float | None:
    if isinstance(value, bool):
        return None
    return float(value) if isinstance(value, (int, float)) else None


def _normalized_effect(signal) -> float | None:
    """Map a signal's effect to a comparable [0, 1] magnitude, or None.

    Effects arrive on incompatible scales: log2FC / motif log2-enrichment are
    log-scale (unbounded), correlations and gene-activity scores already sit near
    [0, 1]. They must not be averaged raw. Log-scale effects are squashed with a
    parameter-free ``tanh`` (|effect|=1 -> 0.76, 2 -> 0.96); already-bounded
    measures are clamped. Non-numeric signals (ORA terms, abundance, L-R pairs)
    have no magnitude and contribute no effect.
    """
    value = _num(getattr(signal, "value", None))
    if value is None:
        return None
    measure = str(getattr(signal, "measure", "") or "")
    if measure.startswith("peak2gene") or measure == "gene_activity":
        return min(abs(value), 1.0)  # correlation / bounded score
    return math.tanh(abs(value))  # log2FC, motif enrichment, generic log-scale


def _signals_by_id(signals: list) -> dict[str, EvidenceSignal]:
    """Index audited signals by ``signal_id`` (first-wins, deterministic)."""
    out: dict[str, EvidenceSignal] = {}
    for sig in signals or []:
        if isinstance(sig, EvidenceSignal) and sig.signal_id:
            out.setdefault(sig.signal_id, sig)
    return out


def _rank_basis(hyp: Hypothesis, by_id: dict) -> dict:
    """Heuristic evidence basis for one hypothesis, scored ONLY over cited signals.

    H18 (Codex M1): independence, effect and confound load come from the exact
    audited signals the hypothesis cites in ``observed_claims`` (resolved by
    ``signal_id``), never from every signal its named entities were measured in.
    An independent line = a distinct ``(audited node, context)`` among the CITED
    signals: the same entity cited in two contrasts is two lines; cited in one is
    one, even if it was measured in three. This is an ordering key, not a
    confidence score.
    """
    lines: set[tuple[str, str]] = set()
    effects: list[float] = []
    caveats: set[str] = set()
    for claim in hyp.observed_claims or []:
        sid = str((claim or {}).get("signal_id", ""))
        sig = by_id.get(sid)
        if sig is None:
            continue
        lines.add((str(sig.audited_node_ref or ""), str(sig.context or "")))
        effect = _normalized_effect(sig)
        if effect is not None:
            effects.append(effect)
        for caveat in sig.caveats_inherited or []:
            caveats.add(str(caveat).strip().lower())
    mean_effect = round(sum(effects) / len(effects), 4) if effects else 0.0
    return {
        "n_independent_lines": len(lines),
        "n_entities": len({_norm(e) for e in (hyp.entities or [])}),
        "mean_effect_norm": mean_effect,
        "evidence_caveat_load": len(caveats),
        "basis": "cited_signals",  # heuristic ordering key, not a confidence score
    }


def _signature(hyp: Hypothesis) -> tuple:
    """Exact-duplicate signature: same entities, mechanism and predicted direction."""
    entities = frozenset(_norm(e) for e in (hyp.entities or []))
    mechanism = str(getattr(hyp, "mechanism", "") or "").strip().lower()
    exp = getattr(hyp, "experiment", None)
    direction = str(getattr(exp, "predicted_direction", "") or "").strip().lower()
    return (entities, mechanism, direction)


def _dedupe(ranked: list[Hypothesis]) -> list[Hypothesis]:
    """Drop exact-duplicate hypotheses, keeping the higher-ranked one.

    Conservative on purpose: only collapses hypotheses with the SAME entities,
    mechanism and predicted direction. Genuine alternatives (same entities, a
    different proposed connection) are rivals, not duplicates, and are kept.
    """
    seen: set[tuple] = set()
    out: list[Hypothesis] = []
    for hyp in ranked:
        sig = _signature(hyp)
        if sig in seen:
            continue
        seen.add(sig)
        out.append(hyp)
    return out


def _annotate_competing(ranked: list[Hypothesis]) -> None:
    """Mark each hypothesis with the rivals that explain overlapping evidence.

    Two surviving hypotheses that share at least one named entity are competing
    explanations; recording them makes the set a genuine rival set instead of one
    clean narrative (rail #4). Deterministic, sorted by id.
    """
    ent_sets = [
        (hyp, {_norm(e) for e in (hyp.entities or [])}) for hyp in ranked
    ]
    for hyp, ents in ent_sets:
        rivals = sorted(
            str(other.id)
            for other, other_ents in ent_sets
            if other is not hyp and ents & other_ents
        )
        hyp.competing_with = rivals


def rank_hypotheses(
    hypotheses: list[Hypothesis], signals: list
) -> list[Hypothesis]:
    """Populate ``rank_evidence`` / ``competing_with`` and return ranked order.

    Order key (descending): independent converging observations, then mean
    NORMALISED effect, then FEWER inherited confounds. Independence and effects
    are derived from the audited signals the hypothesis CITES in
    ``observed_claims`` (H18) — never from the LLM's ``observation_refs`` nor from
    every signal its entities were measured in. Exact duplicates are collapsed
    and the surviving rivals are cross-linked. Deterministic and stable.
    """
    by_id = _signals_by_id(signals)
    decorated: list[tuple[tuple, int, Hypothesis]] = []
    for i, hyp in enumerate(hypotheses or []):
        basis = _rank_basis(hyp, by_id)
        hyp.rank_evidence = basis
        key = (
            basis["n_independent_lines"],
            basis["mean_effect_norm"],
            -basis["evidence_caveat_load"],
        )
        decorated.append((key, i, hyp))
    # Sort by score desc, breaking ties by original order (stable, deterministic).
    decorated.sort(key=lambda d: (tuple(-x for x in d[0]), d[1]))
    ranked = _dedupe([hyp for _, _, hyp in decorated])
    _annotate_competing(ranked)
    return ranked
