"""Registry for modality-specific narrative plugins."""

from __future__ import annotations

import logging
from typing import Iterable

from aria.agents.narrative.protocols import ModalityNarrator
from aria.agents.narrative.types import Caveat, NarrativeBlock
from aria.agents.narrative.validators import validate_blocks
from aria.agents.narrative.warning_router import (
    _AGENT_META,
    _safe_id,
    route_degradation_signals,
)

log = logging.getLogger("aria.narrative")


def _narrator_error_block(agent_name: str, exc: Exception) -> NarrativeBlock:
    """Build an isolated error block when a modality narrator raises."""
    prefix, modality = _AGENT_META.get(
        agent_name, (_safe_id(agent_name), str(agent_name))
    )
    err_msg = f"{type(exc).__name__}: {exc}"
    return NarrativeBlock(
        id=f"{prefix}.narrator_error",
        modality=modality,
        analysis="narrator_error",
        block_type="error",
        title=f"{modality} narrative collection error",
        status="error",
        confidence="insufficient",
        claim="",
        error=err_msg,
        caveats=[
            Caveat(
                text=f"Narrative collection failed for {agent_name}: {err_msg}",
                severity="blocking",
            )
        ],
        metadata={
            "degradation_signals": [
                {
                    "scope": "__narrator__",
                    "field": "exception",
                    "severity": "blocking",
                    "message": err_msg,
                }
            ]
        },
    )


class NarrativeRegistry:
    def __init__(self) -> None:
        self._narrators: list[ModalityNarrator] = []

    def register(self, narrator: ModalityNarrator) -> None:
        self._narrators.append(narrator)

    @property
    def narrators(self) -> tuple[ModalityNarrator, ...]:
        return tuple(self._narrators)

    def collect_blocks(self, agent_results: dict,
                       context: dict | None = None) -> list[NarrativeBlock]:
        blocks: list[NarrativeBlock] = []
        ctx = {**(context or {}), "agent_results": agent_results or {}}
        for agent_name, result in (agent_results or {}).items():
            accepted = False
            for narrator in self._narrators:
                try:  # isolate narrator failures per modality (D5-001)
                    if not narrator.accepts(agent_name, result or {}):
                        continue
                    accepted = True
                    agent_blocks = narrator.collect(agent_name, result or {}, ctx)
                    routed = route_degradation_signals(
                        agent_name, result or {}, agent_blocks
                    )
                    blocks.extend(validate_blocks(routed))
                except Exception as exc:  # aria: allow-broad-except isolate narrator failure per modality
                    accepted = True
                    log.warning(
                        "Narrator %s failed for %s: %s",
                        getattr(narrator, "name", type(narrator).__name__),
                        agent_name,
                        exc,
                        exc_info=True,
                    )
                    blocks.extend(validate_blocks([_narrator_error_block(agent_name, exc)]))
                break
            if (
                not accepted
                and isinstance(result, dict)
                and str(result.get("status") or "").strip().lower() in {"error", "failed"}
            ):
                routed = route_degradation_signals(agent_name, result, [])
                if routed:
                    blocks.extend(validate_blocks(routed))
        return blocks


def registry_with(narrators: Iterable[ModalityNarrator]) -> NarrativeRegistry:
    registry = NarrativeRegistry()
    for narrator in narrators:
        registry.register(narrator)
    return registry
