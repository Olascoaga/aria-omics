"""Preprint-readiness audit A3: per-execution runtime isolation.

Two interleaved runs in one process must not share bus state/persistence, provider
mutable state, LLM usage provenance, or egress policy.
"""
from __future__ import annotations

import json
import threading
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

import pytest

from aria.bus.message_bus import Message, MessageBus, MessageType, bus
from aria.agents.orchestrator_agent import OrchestratorAgent
from aria.llm.provider import LLMProvider, ModelConfig, TaskTier
from aria.runtime.experiment_session import ExperimentSession
from aria.utils.privacy import (
    EgressPolicy,
    egress_allowed,
    execution_environment,
    use_egress_policy,
)
from aria.utils.provenance import collect_llm_usage


def _message(experiment_id: str, label: str) -> Message:
    return Message(
        sender="probe",
        receiver="orchestrator",
        type=MessageType.STATUS,
        payload={"label": label},
        experiment_id=experiment_id,
    )


def _response(text: str):
    usage = SimpleNamespace(
        prompt_tokens=2, completion_tokens=3, total_tokens=5
    )
    choice = SimpleNamespace(
        message=SimpleNamespace(content=text), finish_reason="stop"
    )
    return SimpleNamespace(choices=[choice], usage=usage)


def test_interleaved_sessions_have_distinct_bus_state_and_persistence(tmp_path):
    session_a = ExperimentSession(
        "exp-A", message_bus=MessageBus(persist_path=tmp_path / "A.jsonl")
    )
    session_b = ExperimentSession(
        "exp-B", message_bus=MessageBus(persist_path=tmp_path / "B.jsonl")
    )
    try:
        assert session_a.message_bus is not session_b.message_bus

        # Publish in interleaved order through the compatibility router used by
        # TUI/headless/BaseAgent. The underlying stores remain per execution.
        barrier = threading.Barrier(2)

        def publish_run(experiment_id: str, labels: tuple[str, ...]) -> None:
            for label in labels:
                barrier.wait(timeout=2)
                bus.publish(_message(experiment_id, label))

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [
                pool.submit(publish_run, "exp-A", ("A1", "A2")),
                pool.submit(publish_run, "exp-B", ("B1", "B2")),
            ]
            for future in futures:
                future.result(timeout=5)

        assert [m.payload["label"] for m in session_a.message_bus.get_log()] == [
            "A1", "A2"
        ]
        assert [m.payload["label"] for m in session_b.message_bus.get_log()] == [
            "B1", "B2"
        ]
        assert [m.payload["label"] for m in bus.get_log("exp-A")] == ["A1", "A2"]
        assert [m.payload["label"] for m in bus.get_log("exp-B")] == ["B1", "B2"]

        persisted_a = [json.loads(line) for line in (tmp_path / "A.jsonl").read_text().splitlines()]
        persisted_b = [json.loads(line) for line in (tmp_path / "B.jsonl").read_text().splitlines()]
        assert {item["experiment_id"] for item in persisted_a} == {"exp-A"}
        assert {item["experiment_id"] for item in persisted_b} == {"exp-B"}
    finally:
        bus.unbind_experiment("exp-A")
        bus.unbind_experiment("exp-B")


def test_orchestrator_owns_provider_bus_policy_and_usage_log_per_execution(
    monkeypatch,
):
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")
    base = LLMProvider(models={
        TaskTier.MEDIUM: [
            ModelConfig(
                "ollama", "local-model", 8000, is_local=True,
                api_base="http://localhost:11434",
            )
        ]
    })
    orchestrator = OrchestratorAgent.__new__(OrchestratorAgent)
    orchestrator.llm = base
    orchestrator._sessions = {}

    session_a = orchestrator._get_session("owned-A")
    session_b = orchestrator._get_session("owned-B")
    try:
        assert session_a.message_bus is not session_b.message_bus
        assert session_a.egress_policy is not session_b.egress_policy
        assert session_a.llm_provider is not session_b.llm_provider
        assert session_a.llm_provider.experiment_id == "owned-A"
        assert session_b.llm_provider.experiment_id == "owned-B"
        assert session_a.usage_log != session_b.usage_log
        assert session_a.llm_provider.usage_log == session_a.usage_log
        assert session_b.llm_provider.usage_log == session_b.usage_log
    finally:
        bus.unbind_experiment("owned-A")
        bus.unbind_experiment("owned-B")


def test_execution_egress_policies_do_not_mutate_process_environment(monkeypatch):
    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)
    open_policy = EgressPolicy(air_gapped=False)
    sealed_policy = EgressPolicy(air_gapped=True, reason="sensitive")

    with use_egress_policy(open_policy):
        assert egress_allowed() is True
        assert execution_environment().get("ARIA_AIR_GAPPED") is None
        with use_egress_policy(sealed_policy):
            assert egress_allowed() is False
            assert execution_environment()["ARIA_AIR_GAPPED"] == "1"
        assert egress_allowed() is True

    barrier = threading.Barrier(2)

    def observe(policy: EgressPolicy) -> tuple[bool, str | None]:
        with use_egress_policy(policy):
            barrier.wait(timeout=2)
            return egress_allowed(), execution_environment().get(
                "ARIA_AIR_GAPPED"
            )

    with ThreadPoolExecutor(max_workers=2) as pool:
        open_result = pool.submit(observe, open_policy)
        sealed_result = pool.submit(observe, sealed_policy)
        assert open_result.result(timeout=5) == (True, None)
        assert sealed_result.result(timeout=5) == (False, "1")

    assert "ARIA_AIR_GAPPED" not in __import__("os").environ


def test_cp1_air_gap_choice_seals_only_the_selected_orchestrator_session(
    monkeypatch,
):
    from aria.agents import design_agent as design_module
    from aria.utils.sensitivity import AIR_GAPPED_OPTION

    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)

    class Memory:
        def store_decision(self, **kwargs):
            return None

    class FailingDesignAgent:
        name = "design_agent"

        def __init__(self, memory, llm):
            self.llm = llm

        def start_design(self, **kwargs):
            return {"status": "failed", "reason": "probe complete"}

    orchestrator = OrchestratorAgent.__new__(OrchestratorAgent)
    orchestrator.memory = Memory()
    orchestrator._sessions = {
        "cp1-A": ExperimentSession("cp1-A"),
        "cp1-B": ExperimentSession("cp1-B"),
    }
    orchestrator._experiment_plans = {
        "cp1-A": {"intent": {}},
        "cp1-B": {"intent": {}},
    }
    orchestrator.publish_finding = lambda *args, **kwargs: None
    monkeypatch.setattr(design_module, "DesignAgent", FailingDesignAgent)
    message = Message(
        sender="data_audit_agent",
        receiver="orchestrator",
        type=MessageType.ESCALATION,
        checkpoint=1,
        experiment_id="cp1-A",
        payload={
            "question": "Confirm audit",
            "context": {"exp_context": {}, "sensitivity": {"level": "high"}},
        },
    )

    try:
        result = orchestrator._after_checkpoint_1(
            "cp1-A", AIR_GAPPED_OPTION, message
        )
        assert result == {"status": "cancelled", "reason": "probe complete"}
        assert orchestrator._sessions["cp1-A"].egress_policy.air_gapped is True
        assert orchestrator._sessions["cp1-B"].egress_policy.air_gapped is False
        assert "ARIA_AIR_GAPPED" not in __import__("os").environ
    finally:
        bus.unbind_experiment("cp1-A")
        bus.unbind_experiment("cp1-B")


def test_environment_manager_injects_only_the_active_execution_policy(
    tmp_path, monkeypatch
):
    from aria.utils import environment_manager as environment_module

    seen = {}
    script = tmp_path / "a3_noop.py"
    script.write_text("# subprocess body is replaced by FakeProcess\n")

    class FakeProcess:
        def __init__(self, cmd, **kwargs):
            seen["env"] = kwargs["env"]
            self.output_file = __import__("pathlib").Path(cmd[-1])
            self.returncode = 0

        def communicate(self, timeout):
            self.output_file.write_text('{"status": "success"}')
            return "", ""

    manager = environment_module.EnvironmentManager.__new__(
        environment_module.EnvironmentManager
    )
    manager.workspace = tmp_path / "workspace"
    manager.workspace.mkdir()
    (manager.workspace / "failed").mkdir()
    manager._resolve_env = lambda stack: "aria-env"
    monkeypatch.setattr(environment_module.subprocess, "Popen", FakeProcess)

    with use_egress_policy(EgressPolicy(True, "sensitive")):
        result = manager.run_in_stack(
            "rna", str(script), {"probe": "A3"}, timeout=120
        )

    assert result["status"] == "success"
    assert seen["env"]["ARIA_AIR_GAPPED"] == "1"


def test_per_execution_providers_and_usage_logs_do_not_cross_contaminate(
    tmp_path, monkeypatch
):
    from aria.llm import provider as provider_module

    monkeypatch.setenv("ARIA_LLM_CACHE", "0")
    barrier = threading.Barrier(2)

    def complete_concurrently(**kwargs):
        barrier.wait(timeout=2)
        return _response(str(kwargs["model"]))

    monkeypatch.setattr(provider_module, "completion", complete_concurrently)
    monkeypatch.setattr(
        provider_module.litellm,
        "completion_cost",
        lambda **kwargs: 0.001,
        raising=False,
    )
    models = {
        TaskTier.MEDIUM: [
            ModelConfig("openai", "cloud-model", 8000, is_local=False),
            ModelConfig(
                "ollama", "local-model", 8000, is_local=True,
                api_base="http://localhost:11434",
            ),
        ]
    }
    base = LLMProvider(models=models)
    provider_a = base.for_execution(
        "exp-A", tmp_path / "A-usage.jsonl", EgressPolicy(False)
    )
    provider_b = base.for_execution(
        "exp-B", tmp_path / "B-usage.jsonl", EgressPolicy(True, "sensitive")
    )

    assert provider_a is not provider_b
    with ThreadPoolExecutor(max_workers=2) as pool:
        result_a = pool.submit(provider_a.complete, "A", tier=TaskTier.MEDIUM)
        result_b = pool.submit(provider_b.complete, "B", tier=TaskTier.MEDIUM)
        assert result_a.result(timeout=5) == "cloud-model"
        assert result_b.result(timeout=5) == "local-model"
    assert provider_a._last_completion["model"] == "cloud-model"
    assert provider_b._last_completion["model"] == "local-model"

    usage_a = collect_llm_usage(
        experiment_id="exp-A", usage_log=tmp_path / "A-usage.jsonl"
    )
    usage_b = collect_llm_usage(
        experiment_id="exp-B", usage_log=tmp_path / "B-usage.jsonl"
    )
    assert usage_a["calls"] == 1 and usage_a["models"] == ["cloud-model"]
    assert usage_b["calls"] == 1 and usage_b["models"] == ["local-model"]

    records_a = [json.loads(line) for line in (tmp_path / "A-usage.jsonl").read_text().splitlines()]
    records_b = [json.loads(line) for line in (tmp_path / "B-usage.jsonl").read_text().splitlines()]
    assert {item["experiment_id"] for item in records_a} == {"exp-A"}
    assert {item["experiment_id"] for item in records_b} == {"exp-B"}
