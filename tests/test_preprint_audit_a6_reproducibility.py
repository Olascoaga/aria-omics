"""Preprint-readiness audit A6: immutable runs and honest LLM determinism.

Report bundles are execution artifacts: a later execution must never reuse,
delete, or overwrite an earlier execution's directory.  LLM provenance must
describe controls the answering backend actually received; temperature zero is
not evidence that a backend was seed-deterministic when no seed was sent.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest


def test_reproducible_report_directories_are_immutable_per_run(tmp_path):
    from aria.agents.narrative.report_builder import ReportBuilderMixin
    from aria.runtime.experiment_view import _resolve_report_dir

    class _Host(ReportBuilderMixin):
        def __init__(self, reports_dir):
            self.reports_dir = Path(reports_dir)

        def _build_slug(self, intent, exp_ctx):
            return "slug"

    host = _Host(tmp_path)
    exp_ctx = {
        "reproducible_mode": True,
        "input_files": [{"sha256": "abc123def4567890"}],
    }

    first = host._build_report_dir("run-alpha", {}, exp_ctx)
    first_report = first / "report.html"
    first_manifest = first / "methodology.json"
    first_report.write_text("first output", encoding="utf-8")
    first_manifest.write_text('{"run":"alpha"}', encoding="utf-8")

    with pytest.raises(FileExistsError):
        host._build_report_dir("run-alpha", {}, exp_ctx)

    second = host._build_report_dir("run-beta", {}, exp_ctx)
    (second / "report.html").write_text("second output", encoding="utf-8")
    (second / "methodology.json").write_text(
        '{"run":"beta"}', encoding="utf-8"
    )

    assert second != first
    assert first_report.read_text(encoding="utf-8") == "first output"
    assert first_manifest.read_text(encoding="utf-8") == '{"run":"alpha"}'
    assert (second / "report.html").read_text(encoding="utf-8") == "second output"
    assert not (second / "figures" / "residue-from-alpha.png").exists()
    assert _resolve_report_dir("run-alpha", tmp_path) == first


def test_backend_without_seed_is_not_reported_seed_deterministic(
    monkeypatch, tmp_path
):
    from aria.llm import provider as provider_mod
    from aria.llm.provider import LLMProvider, ModelConfig, TaskTier
    from aria.utils.provenance import collect_llm_usage

    calls = []

    def fake_completion(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            choices=[
                SimpleNamespace(
                    message=SimpleNamespace(content="answer"),
                    finish_reason="stop",
                )
            ],
            usage=SimpleNamespace(
                prompt_tokens=4,
                completion_tokens=2,
                total_tokens=6,
            ),
        )

    monkeypatch.setattr(provider_mod, "completion", fake_completion)
    monkeypatch.setattr(
        provider_mod.litellm, "completion_cost", lambda **_kw: 0.0, raising=False
    )
    monkeypatch.setenv("ARIA_LLM_CACHE", "1")
    usage_log = tmp_path / "usage.jsonl"
    llm = LLMProvider(
        models={
            TaskTier.MEDIUM: [
                ModelConfig("anthropic", "test-anthropic", 128_000),
            ]
        },
        cache_dir=str(tmp_path / "cache"),
        experiment_id="run-anthropic",
        usage_log=usage_log,
    )

    assert llm.complete("prompt", tier=TaskTier.MEDIUM) == "answer"
    assert llm.complete("prompt", tier=TaskTier.MEDIUM) == "answer"
    assert len(calls) == 1
    assert "seed" not in calls[0]

    usage = collect_llm_usage(
        experiment_id="run-anthropic", usage_log=usage_log
    )
    assert usage["calls"] == 2
    assert usage["seed"] is None
    assert usage["seed_applied"] is False
    assert usage["seed_deterministic"] is False
    assert usage["deterministic"] is False
    assert usage["by_model"]["test-anthropic"]["seed_applied"] is False
    assert usage["by_model"]["test-anthropic"]["seed_deterministic"] is False


def test_seed_capable_backend_records_applied_control(monkeypatch, tmp_path):
    from aria.llm import provider as provider_mod
    from aria.llm.provider import LLMProvider, ModelConfig, TaskTier
    from aria.utils.provenance import collect_llm_usage

    calls = []

    def fake_completion(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            choices=[
                SimpleNamespace(
                    message=SimpleNamespace(content="answer"),
                    finish_reason="stop",
                )
            ],
            usage=SimpleNamespace(
                prompt_tokens=4,
                completion_tokens=2,
                total_tokens=6,
            ),
        )

    monkeypatch.setattr(provider_mod, "completion", fake_completion)
    monkeypatch.setattr(
        provider_mod.litellm, "completion_cost", lambda **_kw: 0.0, raising=False
    )
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")
    usage_log = tmp_path / "usage.jsonl"
    llm = LLMProvider(
        models={
            TaskTier.MEDIUM: [ModelConfig("openai", "test-openai", 128_000)]
        },
        experiment_id="run-openai",
        usage_log=usage_log,
    )

    assert llm.complete("prompt", tier=TaskTier.MEDIUM) == "answer"
    assert calls[0]["seed"] == 0

    usage = collect_llm_usage(experiment_id="run-openai", usage_log=usage_log)
    assert usage["seed"] == 0
    assert usage["seed_applied"] is True
    assert usage["seed_deterministic"] is True
    assert usage["deterministic"] is True


def test_report_discloses_backend_specific_determinism():
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    html = agent._build_provenance_section(
        provenance={},
        input_files=[],
        agent_results={},
        llm_usage={
            "calls": 1,
            "temperature": 0.0,
            "temperature_controlled": True,
            "seed": None,
            "seed_applied": False,
            "seed_deterministic": False,
            "deterministic": False,
            "models": ["test-anthropic"],
            "tiers": ["medium"],
            "by_model": {
                "test-anthropic": {
                    "provider": "anthropic",
                    "temperature": 0.0,
                    "temperature_controlled": True,
                    "seed": None,
                    "seed_applied": False,
                    "seed_deterministic": False,
                    "deterministic": False,
                }
            },
        },
    )

    assert "seed_applied" in html
    assert "seed_deterministic" in html
    assert "test-anthropic" in html
    assert "anthropic" in html
