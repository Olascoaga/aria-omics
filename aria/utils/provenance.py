"""Provenance and hashing helpers for ARIA reports."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _git(args: list[str]) -> str:
    try:
        result = subprocess.run(
            ["git", *args],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return "unknown"


def collect_provenance() -> dict[str, Any]:
    """Collect runtime provenance for a report."""
    try:
        from aria import __version__ as aria_version
    except Exception:
        aria_version = "unknown"
    porcelain = _git(["status", "--porcelain"])
    return {
        "aria_version": aria_version,
        "git_sha": _git(["rev-parse", "HEAD"]),
        "git_dirty": bool(porcelain and porcelain != "unknown"),
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "conda_env": os.environ.get("CONDA_PREFIX", ""),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def hash_file(path: str | Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Streaming SHA-256 for large input files."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_params(params: dict) -> str:
    """Order-invariant SHA-256 for JSON-serializable parameter dicts."""
    payload = json.dumps(
        params,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()

