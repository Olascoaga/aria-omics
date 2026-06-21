"""Competitive ranking for the SPECULATIVE tier (ADR-057 rail #4).

A ranked set of rival hypotheses, NOT one clean narrative. The ranking is
deliberately evidence-driven and penalises confidence: a hypothesis ranks higher
when more INDEPENDENT audited lines converge on it and the cited effects are
stronger, and lower when the evidence it leans on carries more inherited
confounds. Narrative elegance is never a factor. The per-hypothesis basis is
recorded in ``Hypothesis.rank_evidence`` so the report can show why one outranks
another.
"""

from __future__ import annotations

from .grounding import build_evidence_index
from .types import Hypothesis


def _norm(entity: str) -> str:
    return str(entity or "").strip().lower()


def _rank_basis(hyp: Hypothesis, index: dict) -> dict:
    nodes: set[str] = set(hyp.observation_refs or [])
    effects: list[float] = []
    caveat_load = 0
    for ent in hyp.entities or []:
        sig = index.get(_norm(ent))
        if sig is None:
            continue
        if sig.audited_node_ref:
            nodes.add(sig.audited_node_ref)
        if isinstance(sig.value, (int, float)):
            effects.append(abs(float(sig.value)))
        caveat_load += len(sig.caveats_inherited or [])
    mean_effect = round(sum(effects) / len(effects), 4) if effects else 0.0
    return {
        "n_independent_lines": len(nodes),
        "n_entities": len({_norm(e) for e in (hyp.entities or [])}),
        "mean_abs_effect": mean_effect,
        "evidence_caveat_load": caveat_load,
    }


def rank_hypotheses(
    hypotheses: list[Hypothesis], signals: list
) -> list[Hypothesis]:
    """Populate ``rank_evidence`` and return the hypotheses in ranked order.

    Order key (descending): independent converging lines, then mean absolute
    effect, then FEWER inherited confounds on the cited evidence. Deterministic
    and stable; never reorders by mechanism prose.
    """
    index = build_evidence_index(signals)
    decorated: list[tuple[tuple, int, Hypothesis]] = []
    for i, hyp in enumerate(hypotheses or []):
        basis = _rank_basis(hyp, index)
        hyp.rank_evidence = basis
        key = (
            basis["n_independent_lines"],
            basis["mean_abs_effect"],
            -basis["evidence_caveat_load"],
        )
        decorated.append((key, i, hyp))
    # Sort by score desc, breaking ties by original order (stable, deterministic).
    decorated.sort(key=lambda d: (tuple(-x for x in d[0]), d[1]))
    return [hyp for _, _, hyp in decorated]
