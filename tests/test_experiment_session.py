from __future__ import annotations

import sys
import types

litellm_stub = types.ModuleType("litellm")
litellm_stub.completion = lambda *args, **kwargs: None
sys.modules.setdefault("litellm", litellm_stub)

from aria.agents.orchestrator_agent import OrchestratorAgent
from aria.runtime.experiment_session import ExperimentSession


class DummyDesignAgent:
    def __init__(self, label: str):
        self.label = label
        self.calls = []

    def handle_user_response(self, experiment_id, checkpoint_num, choice):
        self.calls.append((experiment_id, checkpoint_num, choice))
        return {"status": "next", "step": self.label}


def _orchestrator() -> OrchestratorAgent:
    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch._sessions = {}
    orch._experiment_plans = {
        "exp-A": {"intent": {}, "exp_context": {}},
        "exp-B": {"intent": {}, "exp_context": {}},
    }
    orch._pending_dispatch = {}
    orch._active_design_agent = None
    return orch


def test_experiment_session_scopes_design_agent_per_experiment():
    orch = _orchestrator()
    agent_a = DummyDesignAgent("A")
    agent_b = DummyDesignAgent("B")
    orch._sessions["exp-A"] = ExperimentSession(
        experiment_id="exp-A", design_agent=agent_a)
    orch._sessions["exp-B"] = ExperimentSession(
        experiment_id="exp-B", design_agent=agent_b)

    result = orch._handle_design_checkpoint(
        experiment_id="exp-A",
        checkpoint=2.2,
        user_decision="continue",
    )

    assert result == {"status": "design_in_progress", "next_step": "A"}
    assert agent_a.calls == [("exp-A", 2.2, "continue")]
    assert agent_b.calls == []


def test_experiment_session_scopes_pending_dispatch_by_experiment(monkeypatch):
    orch = _orchestrator()
    orch._sessions["exp-A"] = ExperimentSession(
        experiment_id="exp-A",
        pending_dispatch=({"steps": ["A"]}, {"modalities": {"bulk_RNA": ["a"]}}),
    )
    orch._sessions["exp-B"] = ExperimentSession(
        experiment_id="exp-B",
        pending_dispatch=({"steps": ["B"]}, {"modalities": {"bulk_RNA": ["b"]}}),
    )
    started = []

    class ImmediateThread:
        def __init__(self, target, args, daemon):
            self.target = target
            self.args = args
            self.daemon = daemon

        def start(self):
            started.append(self.args)

    monkeypatch.setattr("aria.agents.orchestrator_agent.threading.Thread",
                        ImmediateThread)
    orch.publish_status = lambda *args, **kwargs: None
    result = orch._after_audit_checkpoint("exp-A", "Proceed anyway", None)

    assert result["status"] == "analysis_running"
    assert started == [("exp-A", {"steps": ["A"]}, {"modalities": {"bulk_RNA": ["a"]}})]
    assert orch._sessions["exp-A"].pending_dispatch is None
    assert orch._sessions["exp-B"].pending_dispatch[0] == {"steps": ["B"]}
