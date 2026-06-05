"""BiologicalSynthesisAgent — evidence-governed scientific synthesis (Slice 1).

Integrates the structured results ARIA already produced and emits an Integrated
Biological Discussion as governed ``NarrativeBlock``s. It does NOT read raw data,
search external literature, propose unmeasured mechanisms, or add general
knowledge: it organizes evidence, detects cross-analysis patterns deterministically,
and writes only what the counts support. ``data_only=True`` is the default and the
only mode in Slice 1 (a ``literature_augmented`` mode is intentionally not built).
"""

from __future__ import annotations

import logging

from aria.agents.narrative.synthesis.pattern_detector import detect_bulk_patterns
from aria.agents.narrative.synthesis.discussion_composer import (
    compose_discussion_blocks,
)

log = logging.getLogger("aria.synthesis")


class BiologicalSynthesisAgent:
    """Compose the integrated biological discussion from structured evidence."""

    def __init__(self, data_only: bool = True):
        # Slice 1 is data-only by contract; the flag exists so a future
        # literature-augmented mode is an explicit, gated opt-in, never a default.
        self.data_only = data_only

    def synthesize(self, agent_results: dict,
                   exp_ctx: dict | None = None) -> list:
        """Return the integration NarrativeBlocks (empty when unsupported)."""
        bulk = (agent_results or {}).get("bulk_rna_agent", {}) or {}
        findings = bulk.get("findings", bulk) if isinstance(bulk, dict) else {}
        contrasts = (findings or {}).get("contrasts", []) or []
        if not contrasts:
            return []
        patterns = detect_bulk_patterns(contrasts)
        return compose_discussion_blocks(patterns)
