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

The wall also guards the rendered ``mechanism`` prose, not only the structured
``entities`` field: an LLM can list grounded entities yet still smuggle an
un-measured entity into the free-text mechanism the reader actually sees. This
reuses W-CLAIM's own named-entity check (``evidence_verifier._claim_entities``)
so a gene-like token named in the mechanism but absent from the audited evidence
is rejected, the same way W-CLAIM rejects it for an audited claim.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from aria.agents.narrative.evidence_verifier import _claim_entities
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
    ungrounded_mechanism_entities: list[str] = field(default_factory=list)
    reason: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def verify_hypothesis_grounding(
    hypothesis: Hypothesis,
    signals: list[EvidenceSignal],
    run_ledger: dict | None = None,
) -> GroundingResult:
    """Reject a hypothesis that invents facts.

    Three mechanical checks:

    1. Every entity in ``hypothesis.entities`` must exist in the audited
       evidence index (built from ``signals``). An entity not measured by any
       audited signal is invented — the hypothesis is rejected.
    2. Every gene-like entity named in the ``hypothesis.mechanism`` PROSE must
       also resolve to the audited evidence. The structured ``entities`` field
       is not the only surface the reader sees; an entity smuggled only into the
       free-text mechanism would otherwise evade check (1). Reuses W-CLAIM's
       ``_claim_entities`` extractor (which already drops non-gene acronyms via
       ``_GENE_STOPWORDS``); an undeclared prose entity absent from evidence is
       rejected.
    3. Every ``hypothesis.observation_refs`` ledger node must exist AND have run
       (status not in not-run/skipped/error). Reusing W-LEDGER, a hypothesis
       cannot arise from an analysis that did not actually produce results.

    The check on (3) only runs when a ``run_ledger`` is supplied; entity
    grounding (1) and mechanism-prose grounding (2) are always enforced.
    """
    evidence = build_evidence_index(signals)
    declared = {_norm(ent) for ent in (hypothesis.entities or [])}
    missing = [
        ent
        for ent in (hypothesis.entities or [])
        if _norm(ent) not in evidence
    ]

    # (2) Prose grounding: gene-like tokens named in the mechanism that are
    # neither in the evidence universe nor among the declared entities are the
    # evasion path the structured-only check (1) misses. Declared-but-missing
    # entities are already reported by (1), so exclude them here to keep the
    # signals distinct.
    ungrounded_prose = sorted(
        {
            token
            for token in _claim_entities(str(hypothesis.mechanism or ""))
            if _norm(token) not in evidence and _norm(token) not in declared
        }
    )

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

    grounded = not missing and not not_run and not ungrounded_prose
    reason = None
    if not grounded:
        parts: list[str] = []
        if missing:
            parts.append(f"ungrounded entities: {sorted(set(missing))}")
        if ungrounded_prose:
            parts.append(f"ungrounded mechanism entities: {ungrounded_prose}")
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
        ungrounded_mechanism_entities=ungrounded_prose,
        reason=reason,
    )
