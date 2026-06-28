"""HypothesisAgent (ADR-057): the SPECULATIVE epistemic tier.

Runs LAST, downstream of W-CLAIM + W-LEDGER passing, over audited artifacts
ONLY. Creativity in the reasoning, zero invention in the facts. The agent is
READ-ONLY: it triggers no compute, changes no thresholds, never reopens DE/DA,
and never writes the audited claim manifest. Its output is quarantined to the
SPECULATIVE tier and (from S4) to a ``hypothesis://`` ledger node that nothing
downstream can promote to a claim.

S1 scope: the agent core + the modality-agnostic Evidence Interface
(``EvidenceSignal``) + the grounding verifier. The proposer (where the LLM lives)
is injected so the core stays deterministic and testable. The default proposer is
honest-null: until a per-modality evidence adapter exists (S5-S8) the agent
proposes nothing rather than fabricating. The falsifiability gate (S2),
devils_advocate (S3), and ledger quarantine (S4) land in later slices.

See ``memory/audit/ARIA_PLAN_HypothesisAgent_2026-06-21.md``.
"""

from __future__ import annotations

from typing import Callable

from aria.agents.narrative.hypothesis.devils_advocate import (
    check_devils_advocate,
)
from aria.agents.narrative.hypothesis.gates import (
    check_falsifiability,
    check_language,
)
from aria.agents.narrative.hypothesis.grounding import (
    build_evidence_index,
    verify_hypothesis_grounding,
)
from aria.agents.narrative.hypothesis.quarantine import quarantine_hypotheses
from aria.agents.narrative.hypothesis.ranking import rank_hypotheses
from aria.agents.narrative.hypothesis.types import EvidenceSignal, Hypothesis

# A proposer turns grounded audited evidence into candidate hypotheses. The LLM
# lives HERE (injected at call time), so the agent core has no LLM dependency and
# can be tested with a deterministic fake. The default proposer is honest-null.
Proposer = Callable[[list[EvidenceSignal], dict], list[Hypothesis]]


def _null_proposer(
    signals: list[EvidenceSignal], exp_ctx: dict
) -> list[Hypothesis]:
    """Honest-null: no evidence adapter yet (S5-S8) -> no hypotheses."""
    return []


def _gate_failure_counts(rejected: list[dict]) -> dict[str, int]:
    """Aggregate why nothing survived, per gate (for honest-null reporting)."""
    counts: dict[str, int] = {}
    for record in rejected:
        for failure in record.get("failures", []) or []:
            gate = failure.get("gate")
            if gate:
                counts[gate] = counts.get(gate, 0) + 1
    return counts


class HypothesisAgent:
    """Generates SPECULATIVE, grounded, machine hypotheses from audited evidence."""

    TIER = "SPECULATIVE"

    def __init__(self, proposer: Proposer | None = None) -> None:
        self._proposer = proposer or _null_proposer

    def generate(
        self,
        signals: list[EvidenceSignal] | None,
        run_ledger: dict | None = None,
        exp_ctx: dict | None = None,
        *,
        w_claim_passed: bool = True,
        w_ledger_passed: bool = True,
    ) -> dict:
        """Produce a grounded SPECULATIVE hypothesis set (read-only).

        Gate #1 (causal gate, ADR-057 rail #1): hypotheses are only generated
        downstream of W-CLAIM + W-LEDGER passing; if verification aborted, the
        agent does not speculate. Every candidate must then clear all publication
        gates — grounding (rail #7), falsifiability (rail #8), language lint
        (rail #11) and the adversarial devils_advocate gate (rail #5). A
        candidate failing any gate is rejected with the failing reasons, never
        caveated into the output. When every candidate is rejected the result is
        honest-null with an aggregated per-gate summary.
        """
        if not (w_claim_passed and w_ledger_passed):
            return {
                "tier": self.TIER,
                "ran": False,
                "requires_ack": True,
                "reason": "verification_gate_not_passed",
                "hypotheses": [],
                "rejected": [],
                "honest_null": True,
            }

        signal_list = list(signals or [])
        evidence_index = build_evidence_index(signal_list)
        candidates = self._proposer(signal_list, exp_ctx or {}) or []
        # ADR-057 rail #10 (C): surface the proposer's parse diagnostic so an
        # honest-null can never silently hide a model that DID answer but whose
        # response we failed to parse (e.g. JSON truncated by the token budget).
        proposer_diag = getattr(self._proposer, "last_diagnostics", None)

        accepted: list[Hypothesis] = []
        rejected: list[dict] = []
        for hyp in candidates:
            failures: list[dict] = []

            grounding = verify_hypothesis_grounding(
                hyp, signal_list, run_ledger
            )
            if not grounding.grounded:
                failures.append(
                    {
                        "gate": "grounding",
                        "passed": False,
                        "reason": grounding.reason,
                        "missing_entities": grounding.missing_entities,
                        "ungrounded_mechanism_entities": (
                            grounding.ungrounded_mechanism_entities
                        ),
                        "not_run_refs": grounding.not_run_refs,
                    }
                )

            for gate in (
                check_falsifiability(hyp),
                check_language(hyp),
                check_devils_advocate(hyp, evidence_index),
            ):
                if not gate.passed:
                    failures.append(gate.to_dict())

            if failures:
                rejected.append(
                    {
                        "hypothesis_id": getattr(hyp, "id", None),
                        "failures": failures,
                    }
                )
            else:
                accepted.append(hyp)

        # Competitive ranking (ADR-057 rail #4): order the surviving set by
        # independent converging evidence / effect strength, penalising confound
        # load — never by narrative elegance. Stamps each rank_evidence.
        accepted = rank_hypotheses(accepted, signal_list)

        # Quarantine: assign each accepted hypothesis its hypothesis:// node and
        # build the non-promotable manifest (ADR-057 rail #6). Done before
        # serialization so the emitted hypotheses carry their ledger_node.
        quarantine = quarantine_hypotheses(accepted)

        result = {
            "tier": self.TIER,
            "ran": True,
            "requires_ack": True,
            "n_evidence": len(signal_list),
            "n_candidates": len(candidates),
            "hypotheses": [h.to_dict() for h in accepted],
            "quarantine": quarantine,
            "rejected": rejected,
            "honest_null": not accepted,
            "proposer_diagnostics": proposer_diag,
        }
        if not accepted and candidates:
            # Honest-null with reasons: don't fail mute when the proposer offered
            # candidates but none earned publication.
            result["null_summary"] = _gate_failure_counts(rejected)
        elif not candidates and isinstance(proposer_diag, dict) and (
            proposer_diag.get("status") not in (None, "ok", "no_candidates")
        ):
            # The proposer produced no candidates for a non-benign reason (e.g.
            # a truncated/invalid response). Name it so the honest-null is honest.
            result["null_reason"] = proposer_diag.get("status")
        return result
