"""Protocols for modality-specific narrative plugins."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

from aria.agents.narrative.types import NarrativeBlock


class ModalityNarrator(Protocol):
    name: str

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        """Return True when this narrator can handle the agent result."""

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        """Convert an agent result into structured narrative blocks."""

    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        """Return copy-pasteable methods text for this modality."""

    def figures(self, agent_name: str, agent_result: dict,
                report_dir: Path | None = None) -> list[dict]:
        """Return figure descriptors, optionally staging files first."""

    def tables(self, agent_name: str, agent_result: dict,
               report_dir: Path | None = None) -> list[dict]:
        """Return table descriptors, optionally staging files first."""
