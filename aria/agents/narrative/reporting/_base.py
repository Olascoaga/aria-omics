"""Shared imports and pure helpers for report-building mixins."""

from __future__ import annotations

import json
import logging
import html as _html
import hashlib
import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from aria import __version__ as ARIA_VERSION
from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock
from aria.agents.narrative.validators import (
    collect_named_entities,
    find_causal_language,
)
from aria.utils.provenance import collect_llm_usage, collect_provenance

log = logging.getLogger("aria.narrative")

_EXEC_SUMMARY_RANGE_LABEL_RE = re.compile(r"\d+(?:-\d+)+")
_EXEC_SUMMARY_THOUSANDS_RE = re.compile(r"(?<=\d),(?=\d{3}(?:\D|$))")
_EXEC_SUMMARY_NUMBER_RE = re.compile(
    r"(?<![A-Za-z])[-+]?\d+(?:,\d{3})*(?:\.\d+)?%?(?![A-Za-z])"
)


def _normalize_exec_summary_number(raw: str) -> str:
    raw = str(raw).strip()
    pct = raw.endswith("%")
    if pct:
        raw = raw[:-1]
    try:
        value = float(raw)
    except ValueError:
        return raw
    if value.is_integer():
        out = str(int(value))
    else:
        out = f"{value:.6g}"
    return out + ("%" if pct else "")


def _executive_summary_numbers(text: str) -> set[str]:
    cleaned = _EXEC_SUMMARY_RANGE_LABEL_RE.sub(" ", str(text or ""))
    cleaned = _EXEC_SUMMARY_THOUSANDS_RE.sub("", cleaned)
    return {
        _normalize_exec_summary_number(match.group(0))
        for match in _EXEC_SUMMARY_NUMBER_RE.finditer(cleaned)
    }

__all__ = [
    "ARIA_VERSION",
    "Caveat",
    "EvidenceItem",
    "NarrativeBlock",
    "Optional",
    "Path",
    "_executive_summary_numbers",
    "_html",
    "_normalize_exec_summary_number",
    "collect_llm_usage",
    "collect_named_entities",
    "collect_provenance",
    "datetime",
    "find_causal_language",
    "hashlib",
    "json",
    "log",
    "logging",
    "re",
]

