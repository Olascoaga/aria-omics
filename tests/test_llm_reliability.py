"""Stage 4 reliability hardening (audit 2026-05-29):

- R3: every LLM call carries a wall-clock timeout so a hung provider cannot
  block the dispatch thread forever.
- R4: tier fallbacks are recorded as model degradation in provenance, and the
  default model IDs are current (no claude-opus-4-7 / gemini-1.5-pro).
- R5: a subprocess timeout kills the whole process group, not just `conda run`.
"""

import os
import signal
import subprocess
import time
from datetime import datetime, timezone

import pytest


# ── R3: LLM call timeout ─────────────────────────────────────────────────────

def _fake_response(text="ok"):
    msg = type("M", (), {"content": text})()
    choice = type("C", (), {"message": msg})()
    return type("R", (), {"choices": [choice], "usage": None})()


def test_completion_call_gets_a_timeout(monkeypatch):
    from aria.llm import provider as pv
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")
    monkeypatch.setattr(pv, "record_llm_usage", lambda e: None)

    seen = {}

    def fake_completion(**kwargs):
        seen.update(kwargs)
        return _fake_response()

    monkeypatch.setattr(pv, "completion", fake_completion)
    prov = pv.LLMProvider(api_keys={"anthropic": "x"})
    prov.complete("hi", tier=pv.TaskTier.LIGHT, max_tokens=8)

    assert "timeout" in seen, "LLM call must pass a timeout to litellm (R3)"
    assert seen["timeout"] == prov._timeout_s


def test_timeout_env_override():
    from aria.llm.provider import LLMProvider
    old = os.environ.get("ARIA_LLM_TIMEOUT")
    os.environ["ARIA_LLM_TIMEOUT"] = "37"
    try:
        assert LLMProvider(api_keys={"anthropic": "x"})._timeout_s == 37.0
    finally:
        if old is None:
            os.environ.pop("ARIA_LLM_TIMEOUT", None)
        else:
            os.environ["ARIA_LLM_TIMEOUT"] = old


# ── R4: model degradation provenance + current model IDs ─────────────────────

def test_default_model_ids_are_current():
    from aria.llm.provider import DEFAULT_MODELS, TaskTier
    heavy = [m.model for m in DEFAULT_MODELS[TaskTier.HEAVY]]
    assert "claude-opus-4-8" in heavy
    assert "claude-opus-4-7" not in heavy
    assert all("gemini-1.5" not in m for m in heavy)


def test_complete_records_fallback_degradation(monkeypatch):
    from aria.llm import provider as pv
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")
    captured = []
    monkeypatch.setattr(pv, "record_llm_usage", lambda e: captured.append(e))

    calls = {"n": 0}

    def fake_completion(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("primary down (429 rate limit)")
        return _fake_response("fallback-answer")

    monkeypatch.setattr(pv, "completion", fake_completion)
    prov = pv.LLMProvider(api_keys={"anthropic": "x", "openai": "x",
                                    "gemini": "x"})
    out = prov.complete("hi", tier=pv.TaskTier.HEAVY, max_tokens=16)

    assert out == "fallback-answer"
    live = [e for e in captured if not e.get("cache_hit")]
    assert live, "a usage event must be recorded for the answering call"
    ev = live[-1]
    assert ev["is_fallback"] is True
    assert ev["fallback_attempt"] == 1
    assert ev["fallback_from"]                       # the skipped primary model
    assert "429" in (ev["fallback_reason"] or "")


def test_collect_llm_usage_flags_degradation(tmp_path, monkeypatch):
    from aria.utils import provenance as prov
    monkeypatch.setattr(prov, "USAGE_LOG", tmp_path / "usage.jsonl")
    ts = datetime.now(timezone.utc).isoformat()
    prov.record_llm_usage({
        "timestamp_utc": ts, "model": "primary", "tier": "heavy",
        "deterministic": True, "temperature": 0.0, "seed": 0,
        "is_fallback": False,
    })
    prov.record_llm_usage({
        "timestamp_utc": ts, "model": "fallback-m", "tier": "heavy",
        "deterministic": True, "temperature": 0.0, "seed": 0,
        "is_fallback": True, "fallback_from": "primary",
        "fallback_reason": "429 rate limit",
    })
    summ = prov.collect_llm_usage()
    assert summ["degraded"] is True
    assert summ["fallback_calls"] == 1
    assert summ["fallbacks"][0]["fallback_from"] == "primary"
    assert summ["fallbacks"][0]["model"] == "fallback-m"


# ── R5: subprocess timeout kills the whole process group ─────────────────────

def test_terminate_process_tree_kills_grandchildren(tmp_path):
    from aria.utils.environment_manager import EnvironmentManager

    pidfile = tmp_path / "child.pid"
    # A group leader (bash) that spawns a grandchild `sleep` and waits. Both
    # belong to the new session, so killpg must reap both.
    proc = subprocess.Popen(
        ["bash", "-c", f"sleep 60 & echo $! > {pidfile}; wait"],
        start_new_session=True,
    )
    for _ in range(50):
        if pidfile.exists() and pidfile.read_text().strip():
            break
        time.sleep(0.1)
    child_pid = int(pidfile.read_text().strip())
    os.kill(child_pid, 0)  # grandchild is alive (raises if not)

    EnvironmentManager._terminate_process_tree(proc)
    time.sleep(0.5)

    assert proc.poll() is not None, "group leader must be reaped"
    with pytest.raises(ProcessLookupError):
        os.kill(child_pid, 0)  # grandchild must be gone, not orphaned
