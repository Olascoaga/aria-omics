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
