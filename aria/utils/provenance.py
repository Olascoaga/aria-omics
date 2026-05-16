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

USAGE_LOG = Path.home() / ".aria" / "llm_usage.jsonl"


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


def record_llm_usage(event: dict) -> None:
    """Append one LLM usage event for later report provenance."""
    try:
        USAGE_LOG.parent.mkdir(parents=True, exist_ok=True)
        with USAGE_LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(event, sort_keys=True, default=str) + "\n")
    except Exception:
        pass


def collect_llm_usage(since_utc: str | None = None) -> dict[str, Any]:
    """Summarize LLM token/cost events since an ISO UTC timestamp."""
    totals = {
        "calls": 0,
        "cache_hits": 0,
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "total_tokens": 0,
        "estimated_cost_usd": 0.0,
        "by_model": {},
    }
    since = None
    if since_utc:
        try:
            since = datetime.fromisoformat(str(since_utc).replace("Z", "+00:00"))
        except Exception:
            return totals
    if not USAGE_LOG.exists():
        return totals
    try:
        with USAGE_LOG.open("r", encoding="utf-8") as fh:
            for line in fh:
                try:
                    event = json.loads(line)
                except Exception:
                    continue
                if since is not None:
                    try:
                        ts = datetime.fromisoformat(
                            str(event.get("timestamp_utc", "")).replace("Z", "+00:00")
                        )
                    except Exception:
                        continue
                    if ts < since:
                        continue
                model = str(event.get("model") or "unknown")
                totals["calls"] += 1
                if event.get("cache_hit"):
                    totals["cache_hits"] += 1
                for key in ("prompt_tokens", "completion_tokens", "total_tokens"):
                    totals[key] += int(event.get(key) or 0)
                totals["estimated_cost_usd"] += float(
                    event.get("estimated_cost_usd") or 0.0
                )
                by_model = totals["by_model"].setdefault(
                    model,
                    {
                        "calls": 0,
                        "prompt_tokens": 0,
                        "completion_tokens": 0,
                        "total_tokens": 0,
                        "estimated_cost_usd": 0.0,
                    },
                )
                by_model["calls"] += 1
                for key in ("prompt_tokens", "completion_tokens", "total_tokens"):
                    by_model[key] += int(event.get(key) or 0)
                by_model["estimated_cost_usd"] += float(
                    event.get("estimated_cost_usd") or 0.0
                )
    except Exception:
        return totals
    totals["estimated_cost_usd"] = round(totals["estimated_cost_usd"], 6)
    for model_totals in totals["by_model"].values():
        model_totals["estimated_cost_usd"] = round(
            model_totals["estimated_cost_usd"], 6
        )
    return totals
