"""Provenance and hashing helpers for ARIA reports."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any

USAGE_LOG = Path.home() / ".aria" / "llm_usage.jsonl"


def collect_provenance() -> dict[str, Any]:
    """Collect runtime provenance for a report."""
    try:
        from aria.version import collect_version_metadata
        version_metadata = collect_version_metadata()
    except Exception:
        version_metadata = {
            "aria_version": "unknown",
            "version": "unknown",
            "version_source": "aria.version.__version__",
            "git_sha": "unknown",
            "git_commit": "unknown",
            "git_dirty": False,
            "git_tree_sha": "unknown",
            "git_describe": "unknown",
            "workflow_hash": "unknown",
            "workflow_hash_algorithm": "unknown",
            "image": {
                "containerized": False,
                "kind": None,
                "digest": None,
                "reference": None,
                "revision": None,
                "env_lock_sha256": None,
                "validation": None,
            },
            "environment": {
                "conda_prefix": None,
                "env_name": None,
                "env_lock_file": None,
                "env_lock_sha256": None,
                "source": "no_conda_env",
            },
        }
    return {
        **version_metadata,
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


def record_llm_usage(
    event: dict,
    *,
    experiment_id: str | None = None,
    usage_log: str | Path | None = None,
) -> None:
    """Append one LLM usage event to its execution-owned provenance log."""
    record = dict(event or {})
    embedded_usage_log = record.pop("_usage_log", None)
    path = Path(usage_log or embedded_usage_log or USAGE_LOG)
    if experiment_id:
        record["experiment_id"] = experiment_id
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(record, sort_keys=True, default=str) + "\n")
    except Exception:
        pass


def collect_llm_usage(
    since_utc: str | None = None,
    grace_seconds: int = 5,
    *,
    experiment_id: str | None = None,
    usage_log: str | Path | None = None,
) -> dict[str, Any]:
    """Summarize one execution's LLM events since an ISO UTC timestamp."""
    totals = {
        "calls": 0,
        "cache_hits": 0,
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "total_tokens": 0,
        "estimated_cost_usd": 0.0,
        "by_model": {},
        "models": [],
        "tiers": [],
        "temperature": None,
        "seed": None,
        "deterministic": True,
        # R4: model-degradation provenance. `degraded` is True if any section
        # ran on a tier fallback instead of the primary model.
        "fallback_calls": 0,
        "degraded": False,
        "fallbacks": [],
    }
    models_seen = set()
    tiers_seen = set()
    temperatures_seen = set()
    seeds_seen = set()
    since = None
    if since_utc:
        try:
            since = datetime.fromisoformat(str(since_utc).replace("Z", "+00:00"))
            if grace_seconds:
                since = since - timedelta(seconds=max(0, int(grace_seconds)))
        except Exception:
            return totals
    path = Path(usage_log) if usage_log is not None else USAGE_LOG
    if not path.exists():
        return totals
    try:
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                try:
                    event = json.loads(line)
                except Exception:
                    continue
                if (experiment_id is not None
                        and event.get("experiment_id") != experiment_id):
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
                tier = str(event.get("tier") or "unknown")
                totals["calls"] += 1
                models_seen.add(model)
                tiers_seen.add(tier)
                if "temperature" in event:
                    try:
                        temperatures_seen.add(float(event.get("temperature")))
                    except Exception:
                        totals["deterministic"] = False
                else:
                    totals["deterministic"] = False
                if "seed" in event:
                    try:
                        seeds_seen.add(int(event.get("seed")))
                    except Exception:
                        totals["deterministic"] = False
                else:
                    totals["deterministic"] = False
                if event.get("deterministic") is not True:
                    totals["deterministic"] = False
                if event.get("cache_hit"):
                    totals["cache_hits"] += 1
                if event.get("is_fallback"):
                    totals["fallback_calls"] += 1
                    totals["degraded"] = True
                    totals["fallbacks"].append({
                        "model": model,
                        "tier": tier,
                        "fallback_from": event.get("fallback_from"),
                        "fallback_reason": event.get("fallback_reason"),
                    })
                for key in ("prompt_tokens", "completion_tokens", "total_tokens"):
                    totals[key] += int(event.get(key) or 0)
                totals["estimated_cost_usd"] += float(
                    event.get("estimated_cost_usd") or 0.0
                )
                by_model = totals["by_model"].setdefault(
                    model,
                    {
                        "calls": 0,
                        "tiers": [],
                        "temperature": event.get("temperature"),
                        "seed": event.get("seed"),
                        "deterministic": event.get("deterministic") is True,
                        "prompt_tokens": 0,
                        "completion_tokens": 0,
                        "total_tokens": 0,
                        "estimated_cost_usd": 0.0,
                    },
                )
                by_model["calls"] += 1
                if tier not in by_model["tiers"]:
                    by_model["tiers"].append(tier)
                    by_model["tiers"].sort()
                if by_model.get("temperature") != event.get("temperature"):
                    by_model["temperature"] = "mixed"
                if by_model.get("seed") != event.get("seed"):
                    by_model["seed"] = "mixed"
                if event.get("deterministic") is not True:
                    by_model["deterministic"] = False
                for key in ("prompt_tokens", "completion_tokens", "total_tokens"):
                    by_model[key] += int(event.get(key) or 0)
                by_model["estimated_cost_usd"] += float(
                    event.get("estimated_cost_usd") or 0.0
                )
    except Exception:
        return totals
    totals["estimated_cost_usd"] = round(totals["estimated_cost_usd"], 6)
    totals["models"] = sorted(models_seen)
    totals["tiers"] = sorted(tiers_seen)
    if len(temperatures_seen) == 1:
        totals["temperature"] = next(iter(temperatures_seen))
    elif len(temperatures_seen) > 1:
        totals["temperature"] = "mixed"
        totals["deterministic"] = False
    if len(seeds_seen) == 1:
        totals["seed"] = next(iter(seeds_seen))
    elif len(seeds_seen) > 1:
        totals["seed"] = "mixed"
        totals["deterministic"] = False
    if temperatures_seen and temperatures_seen != {0.0}:
        totals["deterministic"] = False
    for model_totals in totals["by_model"].values():
        model_totals["estimated_cost_usd"] = round(
            model_totals["estimated_cost_usd"], 6
        )
    return totals
