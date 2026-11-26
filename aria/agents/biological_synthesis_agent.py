"""BiologicalSynthesisAgent — evidence-governed scientific synthesis.

Slice 1a (bulk RNA) and Slice 1b (scRNA) integrate single-modality results;
Slice 2 (ADR-080 / ADR-028) adds the non-same-cell RNA x ATAC cross-modal
synthesis (``integration.cross_modal_*`` blocks), backed by the
``cross_modal_concordance`` engine when it ran, else an in-memory deterministic
fallback. Paired same-cell Multiome triangulation stays with
``IntegrationNarrator`` and is never duplicated here.

Integrates the structured results ARIA already produced and emits an Integrated
Biological Discussion as governed ``NarrativeBlock``s. It does NOT read raw data,
search external literature, propose unmeasured mechanisms, or add general
knowledge: it organizes evidence, detects cross-analysis patterns deterministically,
and writes only what the counts support. ``data_only=True`` is the default and the
only mode in Slice 1 (a ``literature_augmented`` mode is intentionally not built).
"""

from __future__ import annotations

import logging

from aria.agents.narrative.synthesis.pattern_detector import (
    detect_bulk_patterns,
    detect_scrna_patterns,
)
from aria.agents.narrative.synthesis.discussion_composer import (
    compose_discussion_blocks,
    compose_scrna_discussion_blocks,
)
from aria.agents.narrative.synthesis.cross_modal_synthesis import (
    compose_cross_modal_discussion_blocks,
    detect_cross_modal_patterns,
)
from aria.agents.narrative.types import Caveat, NarrativeBlock

log = logging.getLogger("aria.synthesis")


def _synthesis_error_block(slice_key: str, label: str, exc: Exception) -> NarrativeBlock:
    err_msg = f"{type(exc).__name__}: {exc}"
    return NarrativeBlock(
        id=f"integration.{slice_key}_error",
        modality="integrated",
        analysis="integration_limitations",
        block_type="error",
        title=f"Integrated synthesis error ({label})",
        status="error",
        confidence="insufficient",
        claim="",
        error=err_msg,
        caveats=[
            Caveat(
                text=f"Biological synthesis ({label}) failed: {err_msg}",
                severity="blocking",
            )
        ],
        metadata={
            "degradation_signals": [
                {
                    "scope": f"synthesis.{slice_key}",
                    "field": "exception",
                    "severity": "blocking",
                    "message": err_msg,
                }
            ]
        },
    )


class BiologicalSynthesisAgent:
    """Compose the integrated biological discussion from structured evidence."""

    def __init__(self, data_only: bool = True):
        # Slice 1 is data-only by contract; the flag exists so a future
        # literature-augmented mode is an explicit, gated opt-in, never a default.
        self.data_only = data_only

    def synthesize(self, agent_results: dict,
                   exp_ctx: dict | None = None) -> list:
        """Return the integration NarrativeBlocks (empty when unsupported)."""
        blocks = []
        # Slice 1a — bulk RNA.
        try:
            bulk = (agent_results or {}).get("bulk_rna_agent", {}) or {}
            findings = bulk.get("findings", bulk) if isinstance(bulk, dict) else {}
            contrasts = (findings or {}).get("contrasts", []) or []
            if contrasts:
                patterns = detect_bulk_patterns(contrasts)
                blocks.extend(compose_discussion_blocks(patterns))
        except Exception as exc:  # aria: allow-broad-except isolate synthesis slice and emit error block
            log.warning("Slice 1a (bulk) synthesis failed: %s", exc, exc_info=True)
            blocks.append(_synthesis_error_block("bulk_synthesis", "bulk RNA", exc))

        # Slice 1b — scRNA.
        try:
            scrna = (
                (agent_results or {}).get("scrna_agent")
                or (agent_results or {}).get("rna_agent")
                or {}
            )
            if isinstance(scrna, dict):
                patterns = detect_scrna_patterns(scrna)
                blocks.extend(compose_scrna_discussion_blocks(patterns))
        except Exception as exc:  # aria: allow-broad-except isolate synthesis slice and emit error block
            log.warning("Slice 1b (scRNA) synthesis failed: %s", exc, exc_info=True)
            blocks.append(_synthesis_error_block("scrna_synthesis", "scRNA", exc))

        # Slice 2 — non-same-cell cross-modal RNA x ATAC (ADR-080).
        try:
            patterns = detect_cross_modal_patterns(agent_results or {}, exp_ctx)
            blocks.extend(
                compose_cross_modal_discussion_blocks(patterns, agent_results)
            )
        except Exception as exc:  # aria: allow-broad-except isolate synthesis slice and emit error block
            log.warning("Slice 2 (cross-modal) synthesis failed: %s", exc,
                        exc_info=True)
            blocks.append(_synthesis_error_block("cross_modal_synthesis", "cross-modal", exc))
        return blocks
