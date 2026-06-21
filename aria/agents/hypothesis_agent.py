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

from aria.agents.narrative.hypothesis.grounding import (
    verify_hypothesis_grounding,
)
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
        agent does not speculate. Every candidate from the proposer must pass the
        grounding verifier (rail #7) — ungrounded candidates are rejected with a
        reason, never caveated into the output.
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
        candidates = self._proposer(signal_list, exp_ctx or {}) or []

        accepted: list[Hypothesis] = []
        rejected: list[dict] = []
        for hyp in candidates:
            grounding = verify_hypothesis_grounding(
                hyp, signal_list, run_ledger
            )
            if grounding.grounded:
                accepted.append(hyp)
            else:
                rejected.append(
                    {
                        "hypothesis_id": getattr(hyp, "id", None),
                        "grounding": grounding.to_dict(),
                    }
                )

        return {
            "tier": self.TIER,
            "ran": True,
            "requires_ack": True,
            "n_evidence": len(signal_list),
            "n_candidates": len(candidates),
            "hypotheses": [h.to_dict() for h in accepted],
            "rejected": rejected,
            "honest_null": not accepted,
        }
