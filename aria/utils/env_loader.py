"""
ARIA environment loader.

Loads private runtime settings from ~/.aria/.env without requiring the user's
interactive shell to source that file manually. The parser intentionally
supports only the simple shell-assignment format ARIA writes:

    ANTHROPIC_API_KEY=...
    export OPENAI_API_KEY="..."

Existing environment variables win over file values, so one-off terminal
overrides remain possible.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
from pathlib import Path
from typing import Optional

log = logging.getLogger("aria.env_loader")

_LOADED_PATHS: set[Path] = set()
_ENV_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def load_aria_env(path: Optional[str | Path] = None,
                  override: bool = False) -> dict[str, str]:
    """
    Load KEY=VALUE pairs from ~/.aria/.env into os.environ.

    Args:
        path: Optional env-file path. Defaults to ARIA_ENV_FILE or
              ~/.aria/.env.
        override: When False, existing os.environ values are preserved.

    Returns:
        Dict of variables loaded or overwritten during this call.
    """
    env_path = Path(
        path or os.environ.get("ARIA_ENV_FILE") or Path.home() / ".aria" / ".env"
    ).expanduser()
    if not env_path.exists():
        return {}
    if not override and env_path in _LOADED_PATHS:
        return {}

    loaded: dict[str, str] = {}
    try:
        lines = env_path.read_text(encoding="utf-8").splitlines()
    except OSError as e:
        log.warning(f"Could not read ARIA env file {env_path}: {e}")
        return {}

    for raw in lines:
        parsed = _parse_env_line(raw)
        if not parsed:
            continue
        key, value = parsed
        if not override and key in os.environ:
            continue
        os.environ[key] = value
        loaded[key] = value

    _LOADED_PATHS.add(env_path)
    return loaded


def _parse_env_line(raw: str) -> Optional[tuple[str, str]]:
    """Parse one .env line, accepting optional leading 'export'."""
    line = raw.strip()
    if not line or line.startswith("#"):
        return None

    try:
        parts = shlex.split(line, comments=True, posix=True)
    except ValueError:
        return None
    if not parts:
        return None
    if parts[0] == "export":
        parts = parts[1:]
    if len(parts) != 1 or "=" not in parts[0]:
        return None

    key, value = parts[0].split("=", 1)
    key = key.strip()
    if not _ENV_NAME.match(key):
        return None
    return key, value
