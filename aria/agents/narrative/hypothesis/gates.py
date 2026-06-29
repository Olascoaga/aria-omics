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

import re
from dataclasses import dataclass

from .types import Hypothesis

# Best-effort, deterministic lints (no LLM). They are word-boundary matched, not
# raw substrings: a substring lint both MISSED hedges (the old ``"may "`` needed a
# literal trailing space, so ``"...may."`` failed to count as hedged) and FIRED
# falsely (a direction token ``"up"`` matched inside ``"upstream"``). Word
# boundaries fix both; STEM lints still match inflections on purpose
# (``"suggest"`` -> suggests/suggested). These remain a heuristic register check,
# not a semantic judge — a determined paraphrase can still evade them.


def _word_patterns(terms: tuple[str, ...]) -> tuple[re.Pattern, ...]:
    """Match each term as a whole word/phrase (``\\bterm\\b``)."""
    return tuple(re.compile(r"\b" + re.escape(t) + r"\b") for t in terms)


def _stem_patterns(terms: tuple[str, ...]) -> tuple[re.Pattern, ...]:
    """Match each term at a word start, allowing inflected suffixes (``\\bstem``)."""
    return tuple(re.compile(r"\b" + re.escape(t)) for t in terms)


def _first_match(text: str, patterns: tuple[re.Pattern, ...]) -> str | None:
    for pattern in patterns:
        m = pattern.search(text)
        if m:
            return m.group(0)
    return None


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

_VACUOUS_PATTERNS = _word_patterns(_VACUOUS_EXPERIMENT_PHRASES)
_DIRECTION_PATTERNS = _word_patterns(_VALID_DIRECTION_TOKENS)


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
    vacuous = _first_match(blob, _VACUOUS_PATTERNS)
    if vacuous is not None:
        return GateResult(
            "falsifiability",
            False,
            f"discriminating experiment is non-specific ({vacuous!r})",
        )
    direction = str(exp.predicted_direction or "").lower()
    if _first_match(direction, _DIRECTION_PATTERNS) is None:
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
_FORBIDDEN_PATTERNS = _word_patterns(_FORBIDDEN_LANGUAGE)

# At least one explicit hedge must be present so the register stays speculative.
# Whole-word markers are ambiguous as substrings (``may`` is inside "mayonnaise",
# so it needs a boundary); stem markers must still catch inflections
# (``suggest`` -> suggests/suggested), so they are word-START anchored only.
_SPECULATIVE_MARKER_WORDS = (
    "may",
    "might",
    "could",
    "would",
    "candidate",
    "consistent with",
    "raises the possibility",
)
_SPECULATIVE_MARKER_STEMS = (
    "suggest",
    "hypothesiz",
    "propos",
    "possibl",
    "potential",
    "predict",
    "plausibl",
)
_MARKER_PATTERNS = (
    _word_patterns(_SPECULATIVE_MARKER_WORDS)
    + _stem_patterns(_SPECULATIVE_MARKER_STEMS)
)


def check_language(hyp: Hypothesis) -> GateResult:
    """Reject mechanism prose that is assertive/causal or not hedged."""
    text = str(getattr(hyp, "mechanism", "") or "").lower()
    if not text.strip():
        return GateResult("language", False, "empty mechanism")
    hit = _first_match(text, _FORBIDDEN_PATTERNS)
    if hit is not None:
        return GateResult(
            "language",
            False,
            f"mechanism uses assertive/causal language ({hit!r}); the "
            "SPECULATIVE tier requires a hedged register",
        )
    if _first_match(text, _MARKER_PATTERNS) is None:
        return GateResult(
            "language",
            False,
            "mechanism lacks an explicit speculative marker "
            "(e.g. may/could/suggests/we hypothesize)",
        )
    return GateResult("language", True)
