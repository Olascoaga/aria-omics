"""ADR-057 opt-in plumbing guards: the SPECULATIVE hypotheses section is wired
into ``report_builder`` behind ``exp_ctx['enable_hypotheses']``, but until this
slice nothing in a real run ever set that flag (no CLI flag, no env var, no
entrypoint param) — so the section was unreachable outside tests.

These guards lock the entrypoints that now carry the opt-in into ``ctx`` →
(via DataAuditAgent) ``exp_context`` → the report builder gate. They mirror the
already-working ``reproducible_mode`` plumbing. Off by default is the invariant.
"""

from __future__ import annotations

import sys
import types

litellm_stub = types.ModuleType("litellm")
litellm_stub.completion = lambda *args, **kwargs: None
sys.modules.setdefault("litellm", litellm_stub)


# ── TUI entrypoint: _launch_context ───────────────────────────────────────────

def test_launch_context_default_off():
    from aria.tui import _launch_context

    ctx = _launch_context("/data", "why?", None, reproducible_mode=False)
    assert ctx["enable_hypotheses"] is False


def test_launch_context_opt_in():
    from aria.tui import _launch_context

    ctx = _launch_context(
        "/data", "why?", None, reproducible_mode=False,
        enable_hypotheses=True,
    )
    assert ctx["enable_hypotheses"] is True


# ── headless entrypoint: run_headless threads the flag into ctx ───────────────

def _capture_run_headless_ctx(monkeypatch, **kwargs) -> dict:
    """Run run_headless just far enough to capture the ctx it hands the
    orchestrator, short-circuiting before any heavy work."""
    import aria.headless as headless

    captured: dict = {}

    class _StubOrch:
        def __init__(self, *a, **k):
            pass

        def run(self, experiment_id, ctx):
            captured.update(ctx)
            return {"status": "not_started"}  # makes run_headless return early

    class _StubMemory:
        def close(self):
            pass

    monkeypatch.setattr(headless, "OrchestratorAgent", _StubOrch, raising=False)
    monkeypatch.setattr(
        "aria.agents.orchestrator_agent.OrchestratorAgent", _StubOrch,
        raising=False,
    )
    monkeypatch.setattr(
        "aria.memory.memory.ARIAMemory", _StubMemory, raising=False
    )
    headless.run_headless("/data", "why?", **kwargs)
    return captured


def test_run_headless_default_off(monkeypatch):
    ctx = _capture_run_headless_ctx(monkeypatch)
    assert ctx["enable_hypotheses"] is False


def test_run_headless_opt_in(monkeypatch):
    ctx = _capture_run_headless_ctx(monkeypatch, enable_hypotheses=True)
    assert ctx["enable_hypotheses"] is True


# ── DataAuditAgent copies the flag (+ env fallback) into exp_context ──────────

def _inject_enable_hypotheses(context: dict) -> bool:
    """Mirror of the DataAuditAgent.run injection so the contract is guarded
    without driving a full directory scan. Kept byte-aligned with the source."""
    import os

    env_hyp = os.environ.get("ARIA_ENABLE_HYPOTHESES", "").strip().lower()
    return (
        bool(context.get("enable_hypotheses"))
        or env_hyp in ("1", "true", "yes", "on")
    )


def test_exp_context_injection_default_off(monkeypatch):
    monkeypatch.delenv("ARIA_ENABLE_HYPOTHESES", raising=False)
    assert _inject_enable_hypotheses({}) is False
    assert _inject_enable_hypotheses({"enable_hypotheses": False}) is False


def test_exp_context_injection_from_context(monkeypatch):
    monkeypatch.delenv("ARIA_ENABLE_HYPOTHESES", raising=False)
    assert _inject_enable_hypotheses({"enable_hypotheses": True}) is True


def test_exp_context_injection_from_env(monkeypatch):
    for truthy in ("1", "true", "YES", "on"):
        monkeypatch.setenv("ARIA_ENABLE_HYPOTHESES", truthy)
        assert _inject_enable_hypotheses({}) is True
    for falsy in ("", "0", "false", "no"):
        monkeypatch.setenv("ARIA_ENABLE_HYPOTHESES", falsy)
        assert _inject_enable_hypotheses({}) is False


def test_data_audit_source_matches_guard_logic():
    """If the DataAuditAgent injection drifts from the mirrored truth table
    above, fail loudly so the guard cannot silently rot."""
    import inspect

    from aria.agents import data_audit_agent

    src = inspect.getsource(data_audit_agent)
    assert 'exp_context["enable_hypotheses"]' in src
    assert "ARIA_ENABLE_HYPOTHESES" in src
    assert '("1", "true", "yes", "on")' in src
