"""Pure-stdlib bulk RNA sample QC contract helpers shared by scripts and agents.

Lives under ``aria.utils`` so agents and narrators never import compute scripts
under ``aria.scripts`` (which may carry heavy scientific stack dependencies).
"""

from __future__ import annotations

import math


def sample_qc_failed(sqc: dict | None) -> bool:
    """Return True when bulk RNA sample QC failed or left outlier screening incomplete."""
    if not isinstance(sqc, dict) or not sqc:
        return False
    if sqc.get("status") == "error" or bool(sqc.get("error")):
        return True
    if "outliers" in sqc and sqc.get("outliers") is None:
        return True
    if "candidate_outliers" in sqc and sqc.get("candidate_outliers") is None:
        return True
    return False


def sample_qc_size_ratio(sqc: dict | None) -> float | None:
    """Return finite numeric library-size ratio from sample_qc, or None if unavailable."""
    if not isinstance(sqc, dict):
        return None
    val = sqc.get("size_ratio")
    if val is None or isinstance(val, bool):
        return None
    try:
        num = float(val)
    except (TypeError, ValueError):
        return None
    return num if math.isfinite(num) else None


__all__ = ["sample_qc_failed", "sample_qc_size_ratio"]
