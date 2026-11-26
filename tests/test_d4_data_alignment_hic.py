"""Regression tests for D4 silent data misalignment fixes in Hi-C."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _make_cooler(path: Path, seed: int):
    import cooler

    rng = np.random.default_rng(seed)
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 20,
        "start": [i * 100_000 for i in range(20)],
        "end": [(i + 1) * 100_000 for i in range(20)],
    })
    rows = [
        (i, j, int(max(1, round(200 / (1 + j - i) + rng.integers(0, 10)))))
        for i in range(20)
        for j in range(i, 20)
    ]
    cooler.create_cooler(
        str(path), bins, pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"]),
        assembly="hg38",
    )


def test_d4_009_hic_qc_and_topology_fail_closed_on_missing_file_or_sample_id_mismatch(
    tmp_path: Path,
):
    """D4-009: Missing input files or sample_ids length mismatch must fail closed
    in both hic_qc_and_balance and hic_topology._compute_differential_3d."""
    pytest.importorskip("cooler")
    from aria.scripts.hic_qc_and_balance import hic_qc_and_balance
    from aria.scripts.hic_topology import _compute_differential_3d

    c1, c2 = tmp_path / "sample_A.cool", tmp_path / "sample_B.cool"
    _make_cooler(c1, 1)
    _make_cooler(c2, 2)
    missing = tmp_path / "missing.cool"

    # 1. hic_qc_and_balance: missing file must fail closed (not drop + discard sample_ids).
    res_qc_missing = hic_qc_and_balance({
        "files": [str(missing), str(c1), str(c2)],
        "sample_ids": ["WT_rep1", "WT_rep2", "KO_rep1"],
        "resolution": 100_000,
        "output_dir": str(tmp_path / "out_qc_missing"),
        "resume": False,
    })
    assert res_qc_missing["status"] == "error", (
        f"Missing file silently dropped: {res_qc_missing.get('balanced_sample_ids')}"
    )
    assert "missing.cool" in res_qc_missing.get("details", "")

    # 2. hic_qc_and_balance: mismatched sample_ids length must fail closed.
    res_qc_len = hic_qc_and_balance({
        "files": [str(c1), str(c2)],
        "sample_ids": ["WT_rep1"],
        "resolution": 100_000,
        "output_dir": str(tmp_path / "out_qc_len"),
        "resume": False,
    })
    assert res_qc_len["status"] == "error" and res_qc_len.get("error_type") == "SampleIdLengthMismatch"

    # 3. hic_topology._compute_differential_3d: missing file must fail closed.
    res_d3_missing = _compute_differential_3d(
        {
            "files": [str(missing), str(c1), str(c2)],
            "sample_ids": ["WT_rep1", "WT_rep2", "KO_rep1"],
            "sample_metadata": {
                "WT_rep1": {"condition": "WT"},
                "WT_rep2": {"condition": "WT"},
                "KO_rep1": {"condition": "KO"},
            },
            "comparisons": [{"numerator": "KO", "denominator": "WT"}],
            "resolution": 100_000,
        },
        str(tmp_path / "out_d3_missing"),
    )
    assert res_d3_missing["status"] == "error"
    assert "missing.cool" in res_d3_missing.get("details", "")

    # 4. hic_topology._compute_differential_3d: mismatched sample_ids length must fail closed.
    res_d3_len = _compute_differential_3d(
        {
            "files": [str(c1), str(c2)],
            "sample_ids": ["WT_rep1"],
            "sample_metadata": {"WT_rep1": {"condition": "WT"}},
            "comparisons": [{"numerator": "KO", "denominator": "WT"}],
            "resolution": 100_000,
        },
        str(tmp_path / "out_d3_len"),
    )
    assert res_d3_len["status"] == "error" and res_d3_len.get("error_type") == "SampleIdLengthMismatch"
