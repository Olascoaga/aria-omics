"""Real-run bug (2026-06-04): a user whose only API key is GEMINI hit
`RuntimeError: All models failed for tier medium ... OllamaException - Connection
refused`. Two gaps: (1) the MEDIUM/LIGHT default fallback chains never included
Gemini (only HEAVY did), so a Gemini-only user fell through to the dead local
Ollama; (2) the failure surfaced as a raw traceback instead of an actionable hint.
"""

import importlib

import pytest

pytest.importorskip("litellm")

from aria.llm.provider import DEFAULT_MODELS, TaskTier, diagnose_llm_failure


def test_every_tier_default_has_gemini_before_local_ollama():
    # A user with only a GEMINI key must be able to reach a cloud model on every
    # tier, and reach it BEFORE the local Ollama fallback (which is off by default).
    for tier in (TaskTier.HEAVY, TaskTier.MEDIUM, TaskTier.LIGHT):
        providers = [m.provider for m in DEFAULT_MODELS[tier]]
        assert "gemini" in providers, f"{tier.value} lacks a gemini fallback: {providers}"
        assert providers.index("gemini") < providers.index("ollama"), (
            f"{tier.value}: gemini must precede ollama, got {providers}")


def test_diagnose_points_to_the_api_key_the_user_actually_has(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("GEMINI_API_KEY", "x")
    hint = diagnose_llm_failure(
        "All models failed for tier medium. Last error: OllamaException - "
        "[Errno 111] Connection refused")
    low = hint.lower()
    assert "gemini" in low                      # names the available provider
    assert "aria doctor --llm" in low           # points at the diagnostic
    assert "traceback" not in low               # it's a hint, not a dump


def test_seed_is_not_sent_to_providers_that_reject_it(monkeypatch):
    # Real-run bug: forcing seed=0 on every call made gemini/anthropic raise
    # litellm.UnsupportedParamsError (only openai/ollama support `seed`), so the
    # whole tier fell through to the dead Ollama fallback.
    import aria.llm.provider as prov
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")

    captured = {}

    def fake_completion(**kwargs):
        captured.clear()
        captured.update(kwargs)

        class _Msg: content = "ok"
        class _Choice: message = _Msg()
        class _Usage: prompt_tokens = completion_tokens = total_tokens = 1
        class _Resp:
            choices = [_Choice()]
            usage = _Usage()
        return _Resp()

    monkeypatch.setattr(prov, "completion", fake_completion)

    gemini = prov.LLMProvider(
        models={prov.TaskTier.MEDIUM: [
            prov.ModelConfig("gemini", "gemini/gemini-2.5-flash", 1_000_000)]},
        api_keys={"gemini": "x"})
    gemini.complete_medium("hi", max_tokens=50)
    assert "seed" not in captured            # gemini rejects seed -> not sent
    assert captured["temperature"] == 0.0    # determinism via temperature stays

    openai = prov.LLMProvider(
        models={prov.TaskTier.MEDIUM: [
            prov.ModelConfig("openai", "gpt-4o-mini", 128_000)]},
        api_keys={"openai": "x"})
    openai.complete_medium("hi", max_tokens=50)
    assert captured.get("seed") == 0         # openai supports seed -> still sent


def test_diagnose_when_no_key_is_present(monkeypatch):
    for k in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GEMINI_API_KEY",
              "GOOGLE_API_KEY"):
        monkeypatch.delenv(k, raising=False)
    hint = diagnose_llm_failure("All models failed for tier medium.")
    low = hint.lower()
    assert "api key" in low
    assert "ollama" in low or "start" in low    # mentions the local fallback option
    assert "aria doctor --llm" in low
