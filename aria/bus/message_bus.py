"""
ARIA MessageBus
---------------
Inter-agent communication layer.
All internal messages use CavemanMode compression to minimize token usage.
Only NarrativeAgent outputs are decompressed for the user.
"""

from __future__ import annotations
import json
import threading
import logging
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
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

    def __init__(self, max_log_size: int = None, persist_path: str | None = None):
        self._subscribers: dict[str, list] = {}
        self._log: deque[Message] = deque(
            maxlen=max_log_size or self.MAX_LOG_SIZE
        )
        self._agents: dict[str, Any] = {}
        self._lock = threading.RLock()
        self._checkpoint_condition = threading.Condition(self._lock)
        # R6: incremental indices so the 0.5s TUI/headless polls do not scan the
        # full (up to 100k) deque on every tick. Kept eviction-consistent with
        # the FIFO deque.
        self._findings_by_exp: dict[str, list[Message]] = {}
        self._escalations: dict[str, Message] = {}
        # R6: optional append-only durability so a mid-run crash does not lose
        # the findings/decisions of completed stages. Per-run path keeps
        # concurrent runs from sharing a log file.
        self._persist_path: Path | None = (
            Path(persist_path) if persist_path else None
        )
        if self._persist_path is not None:
            try:
                self._persist_path.parent.mkdir(parents=True, exist_ok=True)
            except OSError:
                self._persist_path = None

    def enable_persistence(self, persist_path: str) -> None:
        """Turn on append-only durability for this execution's broker."""
        with self._lock:
            path = Path(persist_path)
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
            except OSError:
                return
            self._persist_path = path

    # ── Index maintenance (R6) ───────────────────────────────────────────
    def _index(self, m: Message) -> None:
        if m.type == MessageType.FINDING:
            self._findings_by_exp.setdefault(m.experiment_id, []).append(m)
        elif m.type == MessageType.ESCALATION:
            self._escalations[m.id] = m

    def _deindex(self, m: Message) -> None:
        if m.type == MessageType.FINDING:
            lst = self._findings_by_exp.get(m.experiment_id)
            if lst:
                try:
                    lst.remove(m)
                except ValueError:
                    pass
                if not lst:
                    self._findings_by_exp.pop(m.experiment_id, None)
        elif m.type == MessageType.ESCALATION:
            self._escalations.pop(m.id, None)

    def _persist(self, m: Message) -> None:
        if self._persist_path is None:
            return
        try:
            with self._persist_path.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(m.to_dict(), default=str) + "\n")
        except Exception:
            log.debug("bus persistence write failed", exc_info=True)

    def register(self, agent_name: str, agent_instance: Any):
        with self._lock:
            self._agents[agent_name] = agent_instance
            self._subscribers.setdefault(agent_name, [])

    def publish(self, message: Message) -> None:
        # Snapshot receivers under the lock, then dispatch outside it so a
        # slow agent.receive() cannot block other publishers.
        with self._lock:
            # Keep the indices eviction-consistent: if the deque is full, the
            # leftmost message is about to be dropped by append().
            if (self._log.maxlen is not None
                    and len(self._log) == self._log.maxlen):
                self._deindex(self._log[0])
            self._log.append(message)
            self._index(message)
            self._persist(message)
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
        # R6: served from the per-experiment index, not a full deque scan.
        with self._lock:
            return list(self._findings_by_exp.get(experiment_id, []))

    def get_pending_checkpoints(
        self, experiment_id: str | None = None
    ) -> list[Message]:
        # R6: scan only the (small) escalation index, and optionally scope to a
        # single experiment so concurrent runs sharing the global bus do not
        # read each other's pending checkpoints.
        with self._lock:
            return [
                m for m in self._escalations.values()
                if m.payload.get("resolved") is not True
                and (experiment_id is None
                     or m.experiment_id == experiment_id)
            ]

    def resolve_checkpoint(self, message_id: str, user_decision: dict) -> None:
        with self._lock:
            m = self._escalations.get(message_id)
            if m is not None and m.type == MessageType.ESCALATION:
                m.payload["resolved"]      = True
                m.payload["user_decision"] = user_decision
                self._checkpoint_condition.notify_all()

    def wait_for_checkpoint_resolution(
        self,
        message_id: str,
        timeout: float | None = None,
    ) -> dict | None:
        """Block until an escalation has a user_decision, or timeout expires."""
        with self._checkpoint_condition:
            def _resolved_message() -> Message | None:
                m = self._escalations.get(message_id)
                if m is not None and m.payload.get("resolved") is True:
                    return m
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

    @staticmethod
    def replay(persist_path: str, experiment_id: str | None = None) -> list[dict]:
        """
        Read back a persisted bus log (R6) for crash postmortem / resume. Returns
        the message dicts in order, optionally scoped to one experiment. Missing
        or unreadable files yield an empty list rather than raising.
        """
        path = Path(persist_path)
        if not path.exists():
            return []
        out: list[dict] = []
        try:
            with path.open("r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if (experiment_id is None
                            or rec.get("experiment_id") == experiment_id):
                        out.append(rec)
        except OSError:
            return out
        return out


class _LazyMessageBus:
    """Compatibility router over isolated per-experiment MessageBus instances.

    Legacy callers (TUI/headless/BaseAgent) import one ``bus`` accessor. A3 keeps
    that API but routes scoped operations to the owning run instead of storing all
    executions in one broker.
    """

    def __init__(self):
        self._instance: MessageBus | None = None
        self._lock = threading.Lock()
        self._routing_lock = threading.RLock()
        self._experiment_buses: dict[str, MessageBus] = {}

    def _get(self) -> MessageBus:
        if self._instance is None:
            with self._lock:
                if self._instance is None:
                    self._instance = MessageBus()
        return self._instance

    def bind_experiment(
        self, experiment_id: str, message_bus: MessageBus
    ) -> MessageBus:
        """Bind one execution to its broker, preserving early legacy messages."""
        if not experiment_id:
            raise ValueError("experiment_id is required to bind a MessageBus")
        with self._routing_lock:
            current = self._experiment_buses.get(experiment_id)
            if current is message_bus:
                return message_bus
            early = self._get().get_log(experiment_id)
            self._experiment_buses[experiment_id] = message_bus
        existing = {m.id for m in message_bus.get_log(experiment_id)}
        for message in early:
            if message.id not in existing:
                message_bus.publish(message)
        return message_bus

    def unbind_experiment(self, experiment_id: str) -> None:
        with self._routing_lock:
            self._experiment_buses.pop(experiment_id, None)

    def for_experiment(self, experiment_id: str) -> MessageBus:
        with self._routing_lock:
            return self._experiment_buses.get(experiment_id) or self._get()

    def publish(self, message: Message) -> None:
        self.for_experiment(message.experiment_id).publish(message)

    def get_log(self, experiment_id: str = None) -> list[Message]:
        if experiment_id:
            return self.for_experiment(experiment_id).get_log(experiment_id)
        with self._routing_lock:
            brokers = [self._get(), *self._experiment_buses.values()]
        messages: dict[str, Message] = {}
        for broker in brokers:
            for message in broker.get_log():
                messages.setdefault(message.id, message)
        return sorted(messages.values(), key=lambda message: message.timestamp)

    def get_findings(self, experiment_id: str) -> list[Message]:
        return self.for_experiment(experiment_id).get_findings(experiment_id)

    def get_pending_checkpoints(
        self, experiment_id: str | None = None
    ) -> list[Message]:
        if experiment_id:
            return self.for_experiment(experiment_id).get_pending_checkpoints(
                experiment_id
            )
        with self._routing_lock:
            brokers = [self._get(), *self._experiment_buses.values()]
        pending: dict[str, Message] = {}
        for broker in brokers:
            for message in broker.get_pending_checkpoints():
                pending.setdefault(message.id, message)
        return sorted(pending.values(), key=lambda message: message.timestamp)

    def _broker_for_message(self, message_id: str) -> MessageBus:
        with self._routing_lock:
            brokers = [*self._experiment_buses.values(), self._get()]
        for broker in brokers:
            if any(message.id == message_id for message in broker.get_log()):
                return broker
        return self._get()

    def resolve_checkpoint(self, message_id: str, user_decision: dict) -> None:
        self._broker_for_message(message_id).resolve_checkpoint(
            message_id, user_decision
        )

    def wait_for_checkpoint_resolution(
        self, message_id: str, timeout: float | None = None
    ) -> dict | None:
        return self._broker_for_message(message_id).wait_for_checkpoint_resolution(
            message_id, timeout=timeout
        )

    def __getattr__(self, name: str):
        return getattr(self._get(), name)


# Global accessor imported by agents; construction is lazy to avoid import-time
# broker state creation in modules that only need enum/dataclass definitions.
bus = _LazyMessageBus()
