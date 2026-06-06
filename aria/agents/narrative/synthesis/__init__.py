"""Evidence-governed biological synthesis (BiologicalSynthesisAgent, Slice 1).

Deterministic cross-analysis pattern detection + a discussion composer that emits
``NarrativeBlock``s, so the integrated discussion inherits ARIA's existing
governance (claim tiering, strict evidence verification, devil's advocate, and
run-ledger linkage) instead of a parallel validator. No LLM decides patterns;
the math is pure set/sign operations over the structured DE + pathway results.
"""

from aria.agents.narrative.synthesis.pattern_detector import (
    detect_bulk_patterns,
    detect_scrna_patterns,
    CrossContrastPattern,
    WithinContrastPattern,
)
from aria.agents.narrative.synthesis.discussion_composer import (
    compose_discussion_blocks,
    compose_scrna_discussion_blocks,
)

__all__ = [
    "detect_bulk_patterns",
    "detect_scrna_patterns",
    "compose_discussion_blocks",
    "compose_scrna_discussion_blocks",
    "CrossContrastPattern",
    "WithinContrastPattern",
]
