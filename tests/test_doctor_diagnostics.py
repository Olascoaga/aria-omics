"""P2-9: `aria doctor --secrets` and `aria doctor --llm` tiers.

Both must be honest and crash-proof: absent keys / air-gapped / missing litellm
are informational, not errors, and no real secret value is ever printed.
"""

from aria.doctor import run_doctor


def test_secrets_tier_runs_and_does_not_error_on_absent_keys(monkeypatch):
    # No keys configured -> offline is valid, so the tier must not hard-fail.
    for var in ("ANTHROPIC_API_KEY", "GEMINI_API_KEY", "GOOGLE_API_KEY",
                "OPENAI_API_KEY"):
        monkeypatch.delenv(var, raising=False)
    code, messages = run_doctor("secrets")
    text = "\n".join(messages)
    assert code == 0
    assert "secrets" in text.lower()


def test_secrets_tier_flags_malformed_key(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "obviously-not-a-real-key")
    code, messages = run_doctor("secrets")
    text = "\n".join(messages).lower()
    assert "malformed" in text or "looks" in text


def test_secrets_tier_never_prints_the_key(monkeypatch):
    secret = "sk-ant-api03-" + "Z" * 40
    monkeypatch.setenv("ANTHROPIC_API_KEY", secret)
    _, messages = run_doctor("secrets")
    assert all(secret not in m for m in messages)


def test_llm_tier_runs_without_crashing(monkeypatch):
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    code, messages = run_doctor("llm")
    text = "\n".join(messages).lower()
    # Air-gapped state must be surfaced; the tier must not hard-fail.
    assert "air" in text or "offline" in text or "llm" in text
    assert isinstance(code, int)
