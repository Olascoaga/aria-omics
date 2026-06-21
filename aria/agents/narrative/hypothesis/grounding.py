"""Grounding verifier (ADR-057 rail #7): zero invention over the facts.

Every entity a hypothesis names must resolve to a real audited
``EvidenceSignal``, and every observation it arises from must cite a run-ledger
node that actually ran. A hypothesis that fails is REJECTED, never caveated.
This is the mechanical wall that lets the LLM be free over the *connection*
while the *facts* stay real — "LLM proposes, code guarantees" on the most
dangerous layer.

Reuses the existing W-LEDGER machinery (``run_ledger._node_index`` +
``_NOT_RUN_STATUSES``) so a hypothesis cannot arise from an analysis the run
marked not-run/skipped/error — the same contradiction W-LEDGER catches for
claims. S1 does not modify the ledger; it only reads it.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from aria.agents.narrative.run_ledger import _NOT_RUN_STATUSES, _node_index

from .types import EvidenceSignal, Hypothesis


def _norm(entity: str) -> str:
    return str(entity or "").strip().lower()


def build_evidence_index(
    signals: list[EvidenceSignal],
) -> dict[str, EvidenceSignal]:
    """Index audited evidence by normalized entity — the grounding universe.

    Only real ``EvidenceSignal`` items with a non-empty entity enter the index;
    anything else is ignored (no fabrication of a grounding target).
    """
    index: dict[str, EvidenceSignal] = {}
    for sig in signals or []:
        if isinstance(sig, EvidenceSignal) and _norm(sig.entity):
            index[_norm(sig.entity)] = sig
    return index


@dataclass
class GroundingResult:
    """Outcome of grounding one hypothesis against audited evidence."""

    grounded: bool
    missing_entities: list[str] = field(default_factory=list)
    not_run_refs: list[dict] = field(default_factory=list)
    reason: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def verify_hypothesis_grounding(
    hypothesis: Hypothesis,
    signals: list[EvidenceSignal],
    run_ledger: dict | None = None,
) -> GroundingResult:
    """Reject a hypothesis that invents facts.

    Two mechanical checks:

    1. Every entity in ``hypothesis.entities`` must exist in the audited
       evidence index (built from ``signals``). An entity not measured by any
       audited signal is invented — the hypothesis is rejected.
    2. Every ``hypothesis.observation_refs`` ledger node must exist AND have run
       (status not in not-run/skipped/error). Reusing W-LEDGER, a hypothesis
       cannot arise from an analysis that did not actually produce results.

    The check on (2) only runs when a ``run_ledger`` is supplied; entity
    grounding (1) is always enforced.
    """
    evidence = build_evidence_index(signals)
    missing = [
        ent
        for ent in (hypothesis.entities or [])
        if _norm(ent) not in evidence
    ]

    not_run: list[dict] = []
    if run_ledger is not None:
        index = _node_index(run_ledger)
        for ref in hypothesis.observation_refs or []:
            node = index.get(ref)
            if node is None:
                not_run.append({"node_id": ref, "status": "no_ledger_node"})
            elif node.get("status") in _NOT_RUN_STATUSES:
                not_run.append(
                    {
                        "node_id": ref,
                        "status": node.get("status"),
                        "reason": node.get("reason"),
                    }
                )

    grounded = not missing and not not_run
    reason = None
    if not grounded:
        parts: list[str] = []
        if missing:
            parts.append(f"ungrounded entities: {sorted(set(missing))}")
        if not_run:
            parts.append(
                "observations not run: "
                f"{[n['node_id'] for n in not_run]}"
            )
        reason = "; ".join(parts)
    return GroundingResult(
        grounded=grounded,
        missing_entities=missing,
        not_run_refs=not_run,
        reason=reason,
    )
