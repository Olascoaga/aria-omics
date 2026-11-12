"""Build a donor-aware scATAC cohort with interval-overlap consensus peaks."""

from __future__ import annotations

import os
from pathlib import Path
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.utils.scatac_cohort_union import merge_scatac_matrices


def chromatin_merge_matrices(params: dict) -> dict:
    output = params.get("output_path") or str(
        Path(params.get("output_dir") or "aria_chromatin") / "merged_scatac_cohort.h5ad")
    result = merge_scatac_matrices(
        params.get("files"), str(output), sample_ids=params.get("sample_ids"),
        sample_metadata=params.get("sample_metadata"),
    )
    if result.get("status") == "success":
        result["validation_level"] = "production"
    return result


if __name__ == "__main__":
    run_script(chromatin_merge_matrices)
