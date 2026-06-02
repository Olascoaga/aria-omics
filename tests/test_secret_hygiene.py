"""P2-9: secret-hygiene detection (pure, no network).

Generic credential-FORMAT detection (an ADR-011 technical-detection exception,
like sensitivity.py) + key classification + a bounded project-file scan. No
hardcoded real secrets; the fixtures are synthetic look-alikes.
"""

from pathlib import Path

from aria.utils.secret_hygiene import (
    classify_key,
    detect_key_patterns,
    mask_secret,
    scan_paths_for_secrets,
)

# Synthetic look-alikes (NOT real credentials).
FAKE_ANTHROPIC = "sk-ant-api03-" + "A" * 40
FAKE_GOOGLE = "AIza" + "B" * 35
FAKE_OPENAI = "sk-" + "C" * 40


def test_detect_anthropic_pattern():
    kinds = detect_key_patterns(f"export ANTHROPIC_API_KEY={FAKE_ANTHROPIC}")
    assert "anthropic" in kinds


def test_detect_google_pattern():
    assert "google" in detect_key_patterns(f"key = '{FAKE_GOOGLE}'")


def test_no_false_positive_on_plain_text():
    assert detect_key_patterns("the quick brown fox sk-ant- jumps") == []
    assert detect_key_patterns("AIza is a prefix but too short") == []


def test_classify_key_states():
    assert classify_key("anthropic", FAKE_ANTHROPIC) == "ok"
    assert classify_key("anthropic", "not-a-key") == "malformed"
    assert classify_key("google", FAKE_GOOGLE) == "ok"
    assert classify_key("google", "AIzaShort") == "malformed"
    assert classify_key("anthropic", "") == "absent"
    assert classify_key("anthropic", None) == "absent"


def test_mask_secret_hides_body():
    masked = mask_secret(FAKE_ANTHROPIC)
    assert FAKE_ANTHROPIC not in masked
    assert masked.startswith("sk-ant")
    assert "…" in masked or "*" in masked


def test_scan_finds_planted_secret_and_ignores_clean(tmp_path: Path):
    bad = tmp_path / "leaky.py"
    bad.write_text(f'ANTHROPIC_API_KEY = "{FAKE_ANTHROPIC}"\n')
    good = tmp_path / "clean.py"
    good.write_text("x = 1  # no secrets here\n")
    hits = scan_paths_for_secrets([bad, good])
    hit_paths = {Path(h["path"]).name for h in hits}
    assert "leaky.py" in hit_paths
    assert "clean.py" not in hit_paths
    # The report must NOT contain the raw secret (masked only).
    for h in hits:
        assert FAKE_ANTHROPIC not in h.get("match", "")


def test_scan_skips_binary_and_missing(tmp_path: Path):
    binf = tmp_path / "blob.bin"
    binf.write_bytes(b"\x00\x01\x02" + FAKE_OPENAI.encode() + b"\xff")
    missing = tmp_path / "nope.py"
    # Should not raise on binary or missing files.
    hits = scan_paths_for_secrets([binf, missing])
    assert isinstance(hits, list)
