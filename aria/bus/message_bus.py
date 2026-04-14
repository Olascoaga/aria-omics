"""
ARIA MessageBus
---------------
Inter-agent communication layer.
All internal messages use CavemanMode compression to minimize token usage.
Only NarrativeAgent outputs are decompressed for the user.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any, Optional
import uuid


class MessageType(Enum):
    FINDING    = "finding"
    REQUEST    = "request"
    ESCALATION = "escalation"
    STATUS     = "status"
    ERROR      = "error"


class Confidence(Enum):
    HIGH         = "high"
    MEDIUM       = "medium"
    LOW          = "low"
    INSUFFICIENT = "insufficient"


class CavemanMode(Enum):
    OFF   = "off"    # Normal prose — user-facing only
    LITE  = "lite"   # Drop filler, keep grammar
    FULL  = "full"   # Fragments OK, short synonyms — DEFAULT internal
    ULTRA = "ultra"  # Max compression, arrows, abbreviations


@dataclass
class Message:
    """
    Universal message unit between ARIA agents.
    Every finding MUST carry a confidence level.
    Every escalation MUST carry a checkpoint number.
    """
    id:            str          = field(default_factory=lambda: str(uuid.uuid4())[:8])
    timestamp:     datetime     = field(default_factory=datetime.now)
    sender:        str          = ""
    receiver:      str          = ""
    type:          MessageType  = MessageType.STATUS
    confidence:    Confidence   = Confidence.MEDIUM
    payload:       dict         = field(default_factory=dict)
    caveman_mode:  CavemanMode  = CavemanMode.FULL
    checkpoint:    Optional[int] = None
    experiment_id: str          = ""

    def to_dict(self) -> dict:
        return {
            "id":            self.id,
            "timestamp":     self.timestamp.isoformat(),
            "sender":        self.sender,
            "receiver":      self.receiver,
            "type":          self.type.value,
            "confidence":    self.confidence.value,
            "payload":       self.payload,
            "caveman_mode":  self.caveman_mode.value,
            "checkpoint":    self.checkpoint,
            "experiment_id": self.experiment_id,
        }


class MessageBus:
    """
    Central message broker for all ARIA agents.
    Agents register themselves and subscribe to message types.
    All internal messages are compressed (CavemanMode.FULL).
    The bus logs every message for full traceability.
    """

    def __init__(self):
        self._subscribers: dict[str, list] = {}
        self._log:         list[Message]   = []
        self._agents:      dict[str, Any]  = {}

    def register(self, agent_name: str, agent_instance: Any):
        self._agents[agent_name] = agent_instance
        self._subscribers.setdefault(agent_name, [])

    def publish(self, message: Message) -> None:
        self._log.append(message)
        receiver = message.receiver
        if receiver == "all":
            for name, agent in self._agents.items():
                if name != message.sender and hasattr(agent, "receive"):
                    agent.receive(message)
        elif receiver in self._agents:
            agent = self._agents[receiver]
            if hasattr(agent, "receive"):
                agent.receive(message)

    def get_log(self, experiment_id: str = None) -> list[Message]:
        if experiment_id:
            return [m for m in self._log if m.experiment_id == experiment_id]
        return self._log

    def get_findings(self, experiment_id: str) -> list[Message]:
        return [
            m for m in self._log
            if m.type == MessageType.FINDING
            and m.experiment_id == experiment_id
        ]

    def get_pending_checkpoints(self) -> list[Message]:
        return [
            m for m in self._log
            if m.type == MessageType.ESCALATION
            and m.payload.get("resolved") is not True
        ]

    def resolve_checkpoint(self, message_id: str, user_decision: dict) -> None:
        for m in self._log:
            if m.id == message_id and m.type == MessageType.ESCALATION:
                m.payload["resolved"]        = True
                m.payload["user_decision"]   = user_decision
                break


# Global bus instance — imported by all agents
bus = MessageBus()
