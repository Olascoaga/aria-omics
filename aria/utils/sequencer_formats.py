"""Supported raw-ingestion boundary and sequencer-native format rejection (E1).

ARIA's supported raw-ingestion entry boundary is **FASTQ** (plus already-supported
downstream forms: 10x MEX matrices and canonical ``.h5ad``). ARIA does NOT perform
demultiplexing (Illumina BCL → FASTQ) or basecalling (Oxford Nanopore POD5/FAST5 →
FASTQ); those are separate upstream stages. A sequencer-native input must be
rejected with a clear boundary message instead of silently falling through to a
generic "no supported inputs" skip.

This module is a pure, dependency-light detector: it classifies a path as an
unsupported sequencer-native format by extension or by a well-known run-folder
marker, and names the upstream stage the user must run first.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

SUPPORTED_ENTRY_BOUNDARY = "FASTQ"

# Each unsupported format: file extensions, run-folder marker basenames, path
# segment markers, the platform, and the upstream stage that produces FASTQ.
_UNSUPPORTED_FORMATS: dict[str, dict[str, Any]] = {
    "BCL": {
        "suffixes": (".bcl", ".bcl.gz", ".cbcl"),
        "markers": ("runinfo.xml", "rtacomplete.txt", "runparameters.xml"),
        "path_markers": ("data/intensities/basecalls",),
        "platform": "Illumina",
        "upstream_stage": "demultiplexing (bcl2fastq / BCL Convert)",
    },
    "POD5": {
        "suffixes": (".pod5",),
        "markers": (),
        "path_markers": (),
        "platform": "Oxford Nanopore",
        "upstream_stage": "basecalling (Dorado)",
    },
    "FAST5": {
        "suffixes": (".fast5",),
        "markers": (),
        "path_markers": (),
        "platform": "Oxford Nanopore",
        "upstream_stage": "basecalling (Dorado / Guppy)",
    },
}


def classify_unsupported(path: str | Path) -> dict[str, Any] | None:
    """Classify one path as an unsupported sequencer-native format, or ``None``.

    Recognises by file extension, by a run-folder marker basename (e.g.
    ``RunInfo.xml``), or by a path segment (e.g. ``Data/Intensities/BaseCalls``).
    """
    raw = str(path)
    lower = raw.lower()
    normalized = lower.replace("\\", "/")
    base = os.path.basename(normalized.rstrip("/"))
    for fmt, spec in _UNSUPPORTED_FORMATS.items():
        if lower.endswith(spec["suffixes"]):
            return _hit(raw, fmt, spec, "extension")
        if base in spec["markers"]:
            return _hit(raw, fmt, spec, "run_folder_marker")
        if any(marker in normalized for marker in spec["path_markers"]):
            return _hit(raw, fmt, spec, "run_folder_path")
    return None


def _hit(path: str, fmt: str, spec: dict, matched_by: str) -> dict[str, Any]:
    return {
        "path": path,
        "format": fmt,
        "platform": spec["platform"],
        "upstream_stage": spec["upstream_stage"],
        "matched_by": matched_by,
    }


def detect_unsupported_inputs(paths) -> list[dict[str, Any]]:
    """Return one detection record per unsupported sequencer-native input path.

    Order-preserving and de-duplicated by path; supported/neutral paths (FASTQ,
    ``.h5ad``, 10x MEX, arbitrary text) yield no records.
    """
    hits: list[dict[str, Any]] = []
    seen: set[str] = set()
    for path in (paths or []):
        record = classify_unsupported(path)
        if record and record["path"] not in seen:
            seen.add(record["path"])
            hits.append(record)
    return hits
