"""Registry for modality-specific narrative plugins."""

from __future__ import annotations

from typing import Iterable

from aria.agents.narrative.protocols import ModalityNarrator
from aria.agents.narrative.types import NarrativeBlock
from aria.agents.narrative.validators import validate_blocks


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
        for agent_name, result in (agent_results or {}).items():
            for narrator in self._narrators:
                if narrator.accepts(agent_name, result or {}):
                    blocks.extend(
                        narrator.collect(agent_name, result or {}, context or {})
                    )
                    break
        return validate_blocks(blocks)


def registry_with(narrators: Iterable[ModalityNarrator]) -> NarrativeRegistry:
    registry = NarrativeRegistry()
    for narrator in narrators:
        registry.register(narrator)
    return registry
