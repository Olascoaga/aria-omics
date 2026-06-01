"""Privacy helpers for local IPC archives and LLM/cache controls."""

from __future__ import annotations

import os
import re
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


def air_gapped_enabled() -> bool:
    """Return True when ARIA must avoid cloud/network LLM calls."""
    return os.environ.get("ARIA_AIR_GAPPED", "0").strip().lower() in {
        "1", "true", "yes", "on",
    }


# Tracks whether air-gapped mode was turned on at runtime (vs. preset in the
# environment), purely for honest report/methodology disclosure.
_runtime_air_gapped_reason: str | None = None


def enable_air_gapped_runtime(reason: str = "user_decision") -> None:
    """Turn on air-gapped mode for the rest of the run (P1-8a).

    Sets ``ARIA_AIR_GAPPED`` in the process environment so BOTH in-process egress
    (LLM) and dispatched subprocesses (which inherit the env via ``conda run``)
    refuse network egress. Used when the user opts into air-gapped at the
    sensitivity checkpoint; ARIA never flips this on without the user's choice.
    """
    global _runtime_air_gapped_reason
    os.environ["ARIA_AIR_GAPPED"] = "1"
    _runtime_air_gapped_reason = reason


def air_gapped_runtime_reason() -> str | None:
    """Reason air-gapped was enabled at runtime, or None if preset/off."""
    return _runtime_air_gapped_reason


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
