"""
ARIA BaseAgent
--------------
All ARIA agents inherit from this class.
Handles: message bus integration, caveman mode, memory access,
         LLM calls (Claude API), and structured finding output.
"""

from __future__ import annotations
import json
import uuid
from abc import ABC, abstractmethod
from typing import Optional

from aria.bus.message_bus import (
    Message, MessageType, Confidence, CavemanMode, bus
)
from aria.memory.memory import ARIAMemory
from aria.llm.provider import LLMProvider, TaskTier


# Caveman prompt injected for internal calls
CAVEMAN_PROMPT = """
Respond terse. Caveman mode ON.
No filler. No pleasantries. No hedging. No articles if avoidable.
Technical terms exact. Code blocks unchanged.
Fragments OK. Short synonyms (big not extensive, fix not implement).
Pattern: [finding] [evidence] [confidence]. [next step].
""".strip()


class BaseAgent(ABC):
    """
    Base class for all ARIA agents.

    Uses LLMProvider for all LLM calls — completely provider-agnostic.
    Caveman compression injected automatically for internal (non-user) calls.
    """

    name: str = "base_agent"
    description: str = ""

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider = None,
                 api_key: str = None):
        self.memory = memory
        # Accept injected provider or build default
        self.llm = llm or LLMProvider.from_config()
        bus.register(self.name, self)

    # ── LLM INTERFACE ────────────────────────────────────────────────────

    def think(self, prompt: str, system: str = "",
              caveman: CavemanMode = CavemanMode.FULL,
              tier: TaskTier = TaskTier.MEDIUM,
              max_tokens: int = 1024) -> str:
        """
        Core LLM call. Provider-agnostic via LLMProvider.
        Caveman compression injected for internal calls.
        """
        system_parts = []
        if system:
            system_parts.append(system)

        if caveman != CavemanMode.OFF:
            system_parts.append(CAVEMAN_PROMPT)
            if caveman == CavemanMode.ULTRA:
                system_parts.append(
                    "ULTRA: arrows→causality. abbrev: DB/fn/impl/req. "
                    "drop conjunctions. one word > two."
                )

        full_system = "\n\n".join(system_parts)
        return self.llm.complete(
            prompt=prompt,
            system=full_system,
            tier=tier,
            max_tokens=max_tokens,
        )

    def think_structured(self, prompt: str, system: str = "",
                         schema_hint: str = "",
                         tier: TaskTier = TaskTier.MEDIUM) -> dict:
        """LLM call returning structured JSON. Always ULTRA caveman."""
        json_instruction = (
            f"Respond ONLY with valid JSON. No preamble. No markdown. "
            f"No backticks. {schema_hint}"
        )
        raw = self.think(
            f"{prompt}\n\n{json_instruction}",
            system=system,
            caveman=CavemanMode.ULTRA,
            tier=tier,
            max_tokens=2048,
        )
        raw = raw.strip()
        if raw.startswith("```"):
            raw = raw.split("```")[1]
            if raw.startswith("json"):
                raw = raw[4:]
        raw = raw.strip("` \n")
        return json.loads(raw)

    # ── MESSAGE BUS ──────────────────────────────────────────────────────

    def publish_finding(self, experiment_id: str, content: dict,
                        confidence: Confidence,
                        receiver: str = "orchestrator"):
        msg = Message(
            id=str(uuid.uuid4())[:8],
            sender=self.name,
            receiver=receiver,
            type=MessageType.FINDING,
            confidence=confidence,
            payload=content,
            caveman_mode=CavemanMode.FULL,
            experiment_id=experiment_id
        )
        bus.publish(msg)
        return msg

    def publish_escalation(self, experiment_id: str, checkpoint: int | float | str,
                           question: str, options: list[str],
                           context: dict = None):
        """Trigger a user checkpoint — pause and wait for decision."""
        msg = Message(
            id=str(uuid.uuid4())[:8],
            sender=self.name,
            receiver="orchestrator",
            type=MessageType.ESCALATION,
            confidence=Confidence.HIGH,
            payload={
                "checkpoint":  checkpoint,
                "question":    question,
                "options":     options,
                "context":     context or {},
                "resolved":    False
            },
            caveman_mode=CavemanMode.OFF,  # User-facing: normal language
            checkpoint=checkpoint,
            experiment_id=experiment_id
        )
        bus.publish(msg)
        return msg

    def publish_status(self, experiment_id: str, message: str,
                       progress: float = None):
        msg = Message(
            sender=self.name,
            receiver="orchestrator",
            type=MessageType.STATUS,
            payload={"message": message, "progress": progress},
            experiment_id=experiment_id
        )
        bus.publish(msg)

    def receive(self, message: Message):
        """Called by MessageBus when this agent receives a message."""
        pass  # Override in subclasses if needed

    # ── ABSTRACT INTERFACE ───────────────────────────────────────────────

    @abstractmethod
    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Main execution method for the agent.
        Must return a dict with at minimum:
          { "status": "done|failed|pending", "findings": [...] }
        """
        ...
