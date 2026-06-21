"""Publication gates for the SPECULATIVE tier (ADR-057 rails #8 + #11).

A hypothesis must EARN the freedom to be published. Beyond grounding (rail #7,
grounding.py), each candidate must clear two deterministic, code-side gates (no
LLM); a failing hypothesis is REJECTED, never caveated:

- Falsifiability (#8): a complete, concrete discriminating experiment. A missing
  or vacuous (c) is rejected. This is what separates a HypothesisAgent from a
  SpeculationAgent — the LLM only earns the freedom to speculate if it also says
  how it could be wrong.
- Language lint (#11): the speculative register. No assertive/causal finding
  verbs, and an explicit hedge marker, so the section can never read as a
  finding. The single place mechanism is allowed lives here, always stamped as
  speculation.
"""

from __future__ import annotations

from dataclasses import dataclass

from .types import Hypothesis


@dataclass
class GateResult:
    gate: str
    passed: bool
    reason: str | None = None

    def to_dict(self) -> dict:
        return {"gate": self.gate, "passed": self.passed, "reason": self.reason}


# A (c) experiment that says nothing — an LLM's favorite non-falsifiable filler.
_VACUOUS_EXPERIMENT_PHRASES = (
    "further study",
    "further studies",
    "further experiment",
    "more experiments",
    "more research",
    "more studies",
    "additional experiments",
    "additional studies",
    "needs validation",
    "should be validated",
    "remains to be",
    "future work",
    "to be determined",
    "tbd",
    "n/a",
)

# A real predicted direction must name a measurable change, not restate a vibe.
_VALID_DIRECTION_TOKENS = (
    "up",
    "down",
    "increase",
    "decrease",
    "higher",
    "lower",
    "gain",
    "loss",
    "lost",
    "positive",
    "negative",
    "no change",
    "unchanged",
    "reduced",
    "elevated",
    "depleted",
    "enriched",
    "abolish",
    "rescue",
)


def check_falsifiability(hyp: Hypothesis) -> GateResult:
    """Reject a hypothesis without a complete, concrete discriminating experiment."""
    exp = getattr(hyp, "experiment", None)
    if exp is None or not exp.is_complete():
        return GateResult(
            "falsifiability",
            False,
            "discriminating experiment incomplete: perturbation, readout, "
            "predicted_direction and refuting_outcome are all required",
        )
    blob = " ".join(
        [
            exp.perturbation,
            exp.readout,
            exp.predicted_direction,
            exp.refuting_outcome,
        ]
    ).lower()
    for phrase in _VACUOUS_EXPERIMENT_PHRASES:
        if phrase in blob:
            return GateResult(
                "falsifiability",
                False,
                f"discriminating experiment is non-specific ({phrase!r})",
            )
    direction = str(exp.predicted_direction or "").lower()
    if not any(tok in direction for tok in _VALID_DIRECTION_TOKENS):
        return GateResult(
            "falsifiability",
            False,
            "predicted_direction does not state a concrete, measurable direction",
        )
    return GateResult("falsifiability", True)


# Assertive/causal verbs banned in the SPECULATIVE tier — they would make a
# speculation read as a finding.
_FORBIDDEN_LANGUAGE = (
    "we find",
    "we found",
    "we show",
    "we demonstrate",
    "demonstrates",
    "we prove",
    "proves",
    "proven",
    "we establish",
    "establishes",
    "we confirm",
    "confirms",
    "confirmed",
    "clearly shows",
    "this shows that",
    "the data show that",
    "causes",
    "is caused by",
    "drives",
    "results in",
    "leads to",
    "is responsible for",
    "determines",
)

# At least one explicit hedge must be present so the register stays speculative.
_SPECULATIVE_MARKERS = (
    "may ",
    "might",
    "could",
    "would",
    "suggest",
    "consistent with",
    "hypothesiz",
    "propose",
    "possibl",
    "potentially",
    "predict",
    "candidate",
    "plausibl",
    "raises the possibility",
)


def check_language(hyp: Hypothesis) -> GateResult:
    """Reject mechanism prose that is assertive/causal or not hedged."""
    text = str(getattr(hyp, "mechanism", "") or "").lower()
    if not text.strip():
        return GateResult("language", False, "empty mechanism")
    hit = next((t for t in _FORBIDDEN_LANGUAGE if t in text), None)
    if hit is not None:
        return GateResult(
            "language",
            False,
            f"mechanism uses assertive/causal language ({hit!r}); the "
            "SPECULATIVE tier requires a hedged register",
        )
    if not any(m in text for m in _SPECULATIVE_MARKERS):
        return GateResult(
            "language",
            False,
            "mechanism lacks an explicit speculative marker "
            "(e.g. may/could/suggests/we hypothesize)",
        )
    return GateResult("language", True)
