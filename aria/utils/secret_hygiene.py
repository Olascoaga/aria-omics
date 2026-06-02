"""Secret-hygiene detection for `aria doctor --secrets` (P2-9).

Detects API-credential FORMATS in text/files and classifies configured keys.
The patterns match credential SHAPES (provider key prefixes), an ADR-011
technical-detection exception like `aria/utils/sensitivity.py` — no biological
content, no hardcoded real secrets. Reporting is always MASKED: a real key value
is never returned or printed.

This complements the CI `gitleaks` job (P2-2) with a local runtime check a user
can run before sharing a repo, and reconciles configured keys with the
providers ARIA actually uses.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

# Credential-shape patterns (provider -> compiled regex). Conservative lengths
# avoid matching the bare prefixes that appear in docs/strings.
_PATTERNS: dict[str, re.Pattern[str]] = {
    "anthropic": re.compile(r"sk-ant-[A-Za-z0-9_\-]{24,}"),
    "google": re.compile(r"AIza[A-Za-z0-9_\-]{35}"),
    # OpenAI keys are `sk-` + a long body; exclude the anthropic `sk-ant-` form
    # by requiring the char after `sk-` not to start the `ant-` prefix.
    "openai": re.compile(r"sk-(?!ant-)[A-Za-z0-9]{20,}"),
}

# Env var each provider key is read from (mirrors LLMProvider._inject_api_keys).
PROVIDER_ENV: dict[str, str] = {
    "anthropic": "ANTHROPIC_API_KEY",
    "google": "GOOGLE_API_KEY",
    "gemini": "GEMINI_API_KEY",
    "openai": "OPENAI_API_KEY",
}

_TEXT_SUFFIXES = {
    ".py", ".sh", ".yml", ".yaml", ".toml", ".cfg", ".ini", ".env", ".txt",
    ".md", ".json", ".lock", ".conf", "",
}
_MAX_SCAN_BYTES = 1_000_000  # skip files larger than ~1MB


def detect_key_patterns(text: str) -> list[str]:
    """Return the provider kinds whose credential shape appears in `text`."""
    if not text:
        return []
    found = []
    for kind, pat in _PATTERNS.items():
        if pat.search(text):
            found.append(kind)
    return found


def mask_secret(value: str, show: int = 6) -> str:
    """Mask a secret for display: keep a short prefix, hide the rest."""
    if not value:
        return ""
    value = str(value)
    if len(value) <= show:
        return "*" * len(value)
    return f"{value[:show]}…{'*' * 4} ({len(value)} chars)"


def classify_key(provider: str, value: str | None) -> str:
    """Classify a configured key as 'ok' | 'malformed' | 'absent'.

    Format-only (no network). 'google'/'gemini' share the AIza shape.
    """
    if value is None or not str(value).strip():
        return "absent"
    value = str(value).strip()
    provider = (provider or "").lower()
    pat = _PATTERNS.get("google" if provider == "gemini" else provider)
    if pat is None:
        # Unknown provider: accept any sufficiently long opaque token.
        return "ok" if len(value) >= 20 else "malformed"
    # Anchor the shape to the whole value (a key field should BE the key).
    return "ok" if pat.fullmatch(value) else "malformed"


def _looks_text(path: Path) -> bool:
    if path.suffix.lower() not in _TEXT_SUFFIXES:
        return False
    try:
        if path.stat().st_size > _MAX_SCAN_BYTES:
            return False
    except OSError:
        return False
    return True


def scan_paths_for_secrets(paths: list[Path | str]) -> list[dict[str, Any]]:
    """Scan files for committed credential shapes.

    Returns a list of {path, kind, match} dicts where `match` is MASKED (the raw
    secret is never returned). Missing/binary/oversized files are skipped
    silently — this is a hygiene aid, not a parser.
    """
    hits: list[dict[str, Any]] = []
    for raw in paths:
        path = Path(raw)
        if not path.is_file() or not _looks_text(path):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue  # binary or unreadable
        for kind, pat in _PATTERNS.items():
            m = pat.search(text)
            if m:
                hits.append({
                    "path": str(path),
                    "kind": kind,
                    "match": mask_secret(m.group(0)),
                })
    return hits
