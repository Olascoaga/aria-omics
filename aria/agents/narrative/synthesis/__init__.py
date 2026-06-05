"""Evidence-governed biological synthesis (BiologicalSynthesisAgent, Slice 1).

Deterministic cross-analysis pattern detection + a discussion composer that emits
``NarrativeBlock``s, so the integrated discussion inherits ARIA's existing
governance (claim tiering, strict evidence verification, devil's advocate, and
run-ledger linkage) instead of a parallel validator. No LLM decides patterns;
the math is pure set/sign operations over the structured DE + pathway results.
"""

from aria.agents.narrative.synthesis.pattern_detector import (
    detect_bulk_patterns,
    CrossContrastPattern,
    WithinContrastPattern,
)

__all__ = [
    "detect_bulk_patterns",
    "CrossContrastPattern",
    "WithinContrastPattern",
]
