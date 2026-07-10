"""Privacy helpers for local IPC archives and LLM/cache controls."""

from __future__ import annotations

import os
import re
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_PATH_KEYS = {
    "path",
    "data_path",
    "input_path",
    "output_path",
    "output_dir",
    "fastq_dir",
    "genome_dir",
    "genome_fasta",
    "gtf_file",
    "bam_files",
    "samples",
}
_SECRET_RE = re.compile(r"(api[_-]?key|token|secret|password|credential)", re.I)


@dataclass
class EgressPolicy:
    """Mutable network policy owned by exactly one experiment execution."""

    air_gapped: bool = False
    reason: str | None = None

    @classmethod
    def from_environment(cls) -> "EgressPolicy":
        enabled = os.environ.get("ARIA_AIR_GAPPED", "0").strip().lower() in {
            "1", "true", "yes", "on",
        }
        return cls(air_gapped=enabled, reason="environment" if enabled else None)

    def enable(self, reason: str = "user_decision") -> None:
        self.air_gapped = True
        self.reason = reason


_ACTIVE_EGRESS_POLICY: ContextVar[EgressPolicy | None] = ContextVar(
    "aria_active_egress_policy", default=None
)


@contextmanager
def use_egress_policy(policy: EgressPolicy):
    """Bind one execution's egress policy for the current thread/context."""
    token = _ACTIVE_EGRESS_POLICY.set(policy)
    try:
        yield policy
    finally:
        _ACTIVE_EGRESS_POLICY.reset(token)


def active_egress_policy() -> EgressPolicy | None:
    return _ACTIVE_EGRESS_POLICY.get()


def air_gapped_enabled() -> bool:
    """Return True when ARIA must avoid cloud/network LLM calls."""
    policy = active_egress_policy()
    if policy is not None:
        return bool(policy.air_gapped)
    return os.environ.get("ARIA_AIR_GAPPED", "0").strip().lower() in {
        "1", "true", "yes", "on",
    }


# Tracks whether air-gapped mode was turned on at runtime (vs. preset in the
# environment), purely for honest report/methodology disclosure.
_runtime_air_gapped_reason: str | None = None


def enable_air_gapped_runtime(
    reason: str = "user_decision",
    policy: EgressPolicy | None = None,
) -> None:
    """Turn on air-gapped mode for one run (P1-8a/A3).

    An explicit or context-bound policy stays execution-local and is injected
    only into that run's child-process environment. Calls without either retain
    the legacy process-wide behavior for backward compatibility.
    """
    target = policy or active_egress_policy()
    if target is not None:
        target.enable(reason)
        return
    global _runtime_air_gapped_reason
    os.environ["ARIA_AIR_GAPPED"] = "1"
    _runtime_air_gapped_reason = reason


def air_gapped_runtime_reason(
    policy: EgressPolicy | None = None,
) -> str | None:
    """Reason air-gapped was enabled at runtime, or None if preset/off."""
    target = policy or active_egress_policy()
    if target is not None:
        return target.reason
    return _runtime_air_gapped_reason


def execution_environment(base: dict[str, str] | None = None) -> dict[str, str]:
    """Return a subprocess environment carrying only the active run's policy."""
    env = dict(os.environ if base is None else base)
    policy = active_egress_policy()
    if policy is None:
        return env
    if policy.air_gapped:
        env["ARIA_AIR_GAPPED"] = "1"
    else:
        env.pop("ARIA_AIR_GAPPED", None)
    return env


class EgressBlocked(RuntimeError):
    """Raised when a network egress is attempted under air-gapped mode."""


def egress_allowed() -> bool:
    """W-PRIV: network egress (Enrichr ORA, GEO/SRA fetch, …) is allowed only
    when ARIA is not air-gapped. Air-gapped governs ALL egress, not just the LLM
    layer, so a sensitive run never ships gene lists or fetches remote data."""
    return not air_gapped_enabled()


def assert_egress_allowed(channel: str) -> None:
    """Refuse a network egress when air-gapped. `channel` names the destination
    (e.g. 'enrichr', 'GEO/SRA') for the error message."""
    if not egress_allowed():
        raise EgressBlocked(
            f"ARIA_AIR_GAPPED is enabled; refusing network egress to '{channel}'. "
            f"Disable air-gapped mode to allow it, or use a local alternative."
        )


def _redact_scalar(value: Any) -> Any:
    if isinstance(value, (str, os.PathLike)):
        text = str(value)
        if "/" in text or "\\" in text:
            return f"<path:{Path(text).name}>"
    return value


def redact_sensitive_params(value: Any, *, key: str = "") -> Any:
    """Recursively redact paths/secrets from diagnostic JSON.

    Runtime IPC still needs the real paths. This helper is for persisted
    postmortem artifacts, logs, and metadata that do not need path contents.
    """
    key_l = key.lower()
    if _SECRET_RE.search(key_l):
        return "<redacted>"
    if key_l in _PATH_KEYS or key_l.endswith("_path") or key_l.endswith("_dir"):
        if isinstance(value, list):
            return [
                redact_sensitive_params(v, key="")
                if isinstance(v, (dict, list)) else _redact_scalar(v)
                for v in value
            ]
        if isinstance(value, dict):
            return redact_sensitive_params(value)
        return _redact_scalar(value)
    if isinstance(value, dict):
        return {
            str(k): redact_sensitive_params(v, key=str(k))
            for k, v in value.items()
        }
    if isinstance(value, list):
        return [redact_sensitive_params(v, key=key) for v in value]
    return value
