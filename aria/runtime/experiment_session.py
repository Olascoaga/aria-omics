"""Per-experiment runtime state.

This object is intentionally small and mutable. It keeps orchestration state
scoped by experiment_id so concurrent runs do not share design-agent or pending
dispatch state through one process-global slot.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from aria.bus.message_bus import MessageBus, bus
from aria.utils.privacy import EgressPolicy


@dataclass
class ExperimentSession:
    experiment_id: str
    context: dict[str, Any] = field(default_factory=dict)
    intent: dict[str, Any] = field(default_factory=dict)
    exp_context: dict[str, Any] = field(default_factory=dict)
    plan: dict[str, Any] = field(default_factory=dict)
    status: str = "new"
    findings: list[Any] = field(default_factory=list)
    agent_results: dict[str, Any] = field(default_factory=dict)
    design_agent: Any = None
    pending_dispatch: tuple[dict, dict] | None = None
    log_handler: Any = None
    message_bus: MessageBus = field(default_factory=MessageBus)
    egress_policy: EgressPolicy = field(default_factory=EgressPolicy.from_environment)
    llm_provider: Any = None
    usage_log: Any = None
    run_ledger: dict[str, Any] = field(default_factory=dict)
    cache: dict[str, Any] = field(default_factory=dict)
    locks: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        bus.bind_experiment(self.experiment_id, self.message_bus)

    def as_plan_record(self) -> dict[str, Any]:
        return {
            "user_question": self.context.get("user_question", ""),
            "intent": self.intent,
            "context": self.context,
            "status": self.status,
            "findings": self.findings,
            "exp_context": self.exp_context,
            "plan": self.plan,
        }
