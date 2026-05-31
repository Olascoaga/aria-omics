"""P0-2 regression: LLMProvider.complete() must survive a partial model config.

`complete()` resolved candidates with
`self.models.get(tier, self.models[TaskTier.MEDIUM])`. Python evaluates the
default argument eagerly, so when the configured `models` dict has no MEDIUM tier,
`self.models[TaskTier.MEDIUM]` raises KeyError on EVERY call — even for a tier that
IS configured (e.g. HEAVY-only setups). Resolution must be lazy: requested tier
first, then a MEDIUM fallback, then an explicit RuntimeError.
"""

import pytest


def _fake_response(text="ok"):
    msg = type("M", (), {"content": text})()
    choice = type("C", (), {"message": msg})()
    return type("R", (), {"choices": [choice], "usage": None})()


def _provider(monkeypatch, models):
    from aria.llm import provider as pv
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")
    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)
    monkeypatch.setattr(pv, "record_llm_usage", lambda e: None)
    monkeypatch.setattr(pv, "completion", lambda **kw: _fake_response())
    return pv, pv.LLMProvider(models=models, api_keys={"anthropic": "x"})


def test_heavy_only_config_does_not_keyerror_on_present_tier(monkeypatch):
    from aria.llm.provider import TaskTier, ModelConfig
    pv, prov = _provider(monkeypatch, {
        TaskTier.HEAVY: [ModelConfig("anthropic", "claude-opus-4-8", 200_000)],
    })
    # HEAVY is configured; the absent MEDIUM tier must not crash the call.
    assert prov.complete("hi", tier=TaskTier.HEAVY, max_tokens=8) == "ok"


def test_absent_tier_falls_back_to_medium(monkeypatch):
    from aria.llm.provider import TaskTier, ModelConfig
    pv, prov = _provider(monkeypatch, {
        TaskTier.HEAVY:  [ModelConfig("anthropic", "claude-opus-4-8", 200_000)],
        TaskTier.MEDIUM: [ModelConfig("anthropic", "claude-sonnet-4-6", 200_000)],
    })
    # LIGHT is not configured -> resolve via the MEDIUM fallback.
    assert prov.complete("hi", tier=TaskTier.LIGHT, max_tokens=8) == "ok"


def test_no_tier_and_no_medium_raises_explicit_runtimeerror(monkeypatch):
    from aria.llm.provider import TaskTier, ModelConfig
    pv, prov = _provider(monkeypatch, {
        TaskTier.LIGHT: [ModelConfig("anthropic", "claude-haiku-4-5-20251001",
                                     200_000)],
    })
    # HEAVY absent and no MEDIUM fallback -> explicit RuntimeError, not KeyError.
    with pytest.raises(RuntimeError):
        prov.complete("hi", tier=TaskTier.HEAVY, max_tokens=8)
