"""
ARIA MessageBus
---------------
Inter-agent communication layer.
All internal messages use CavemanMode compression to minimize token usage.
Only NarrativeAgent outputs are decompressed for the user.
"""

from __future__ import annotations
import threading
import logging
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any, Optional
import uuid

log = logging.getLogger("aria.bus")


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
    Every escalation MUST carry a checkpoint identifier.
    """
    id:            str          = field(default_factory=lambda: str(uuid.uuid4())[:8])
    timestamp:     datetime     = field(default_factory=datetime.now)
    sender:        str          = ""
    receiver:      str          = ""
    type:          MessageType  = MessageType.STATUS
    confidence:    Confidence   = Confidence.MEDIUM
    payload:       dict         = field(default_factory=dict)
    caveman_mode:  CavemanMode  = CavemanMode.FULL
    checkpoint:    Optional[int | float | str] = None
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

    # Maximum messages kept in memory. Old messages are evicted (FIFO).
    # 100k covers very long sessions; each Message is small (~1-2KB).
    MAX_LOG_SIZE = 100_000

    def __init__(self, max_log_size: int = None):
        self._subscribers: dict[str, list] = {}
        self._log: deque[Message] = deque(
            maxlen=max_log_size or self.MAX_LOG_SIZE
        )
        self._agents: dict[str, Any] = {}
        self._lock = threading.RLock()
        self._checkpoint_condition = threading.Condition(self._lock)

    def register(self, agent_name: str, agent_instance: Any):
        with self._lock:
            self._agents[agent_name] = agent_instance
            self._subscribers.setdefault(agent_name, [])

    def publish(self, message: Message) -> None:
        # Snapshot receivers under the lock, then dispatch outside it so a
        # slow agent.receive() cannot block other publishers.
        with self._lock:
            self._log.append(message)
            receiver = message.receiver
            if receiver == "all":
                targets = [(n, a) for n, a in self._agents.items()
                           if n != message.sender]
            elif receiver in self._agents:
                targets = [(receiver, self._agents[receiver])]
            else:
                targets = []

        for name, agent in targets:
            if hasattr(agent, "receive"):
                try:
                    agent.receive(message)
                except Exception:
                    log.exception(
                        "MessageBus receiver failed; continuing fan-out "
                        "message_id=%s receiver=%s sender=%s type=%s",
                        message.id,
                        name,
                        message.sender,
                        message.type.value,
                    )

    def get_log(self, experiment_id: str = None) -> list[Message]:
        with self._lock:
            snapshot = list(self._log)
        if experiment_id:
            return [m for m in snapshot if m.experiment_id == experiment_id]
        return snapshot

    def get_findings(self, experiment_id: str) -> list[Message]:
        with self._lock:
            snapshot = list(self._log)
        return [
            m for m in snapshot
            if m.type == MessageType.FINDING
            and m.experiment_id == experiment_id
        ]

    def get_pending_checkpoints(self) -> list[Message]:
        with self._lock:
            snapshot = list(self._log)
        return [
            m for m in snapshot
            if m.type == MessageType.ESCALATION
            and m.payload.get("resolved") is not True
        ]

    def resolve_checkpoint(self, message_id: str, user_decision: dict) -> None:
        with self._lock:
            for m in self._log:
                if m.id == message_id and m.type == MessageType.ESCALATION:
                    m.payload["resolved"]      = True
                    m.payload["user_decision"] = user_decision
                    self._checkpoint_condition.notify_all()
                    break

    def wait_for_checkpoint_resolution(
        self,
        message_id: str,
        timeout: float | None = None,
    ) -> dict | None:
        """Block until an escalation has a user_decision, or timeout expires."""
        with self._checkpoint_condition:
            def _resolved_message() -> Message | None:
                for m in self._log:
                    if m.id == message_id and m.type == MessageType.ESCALATION:
                        if m.payload.get("resolved") is True:
                            return m
                        return None
                return None

            msg = _resolved_message()
            if msg is None:
                self._checkpoint_condition.wait_for(
                    lambda: _resolved_message() is not None,
                    timeout=timeout,
                )
                msg = _resolved_message()
            if msg is None:
                return None
            return msg.payload.get("user_decision") or {}


# Global bus instance — imported by all agents
bus = MessageBus()
