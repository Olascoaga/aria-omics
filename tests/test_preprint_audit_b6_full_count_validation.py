"""Preprint-readiness audit B6: full vectorized raw-count validation.

`classify_matrix` decides `is_raw_counts` from a <=200-row random sample, and
`rna_bulk_de._load_counts` then rounds the WHOLE matrix to int. A fractional,
negative, NaN, or inf value in a row OUTSIDE the sampled window was therefore
accepted as raw counts and silently coerced (or crashed `.astype(int)` on NaN).

After B6 a dedicated full-matrix validator (`validate_raw_count_matrix`) inspects
every value; any non-integer / negative / NaN / inf ANYWHERE blocks the raw path.
The classifier stays a fast sampled heuristic; the full validator is the
deterministic gate that runs before rounding for DESeq2.
"""
from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")


def _validate():
    from aria.utils.count_classifier import validate_raw_count_matrix
    return validate_raw_count_matrix


# ── The full vectorized validator ─────────────────────────────────────────────

def test_all_integer_nonnegative_is_valid():
    mat = np.arange(0, 600, dtype=float).reshape(200, 3)
    result = _validate()(mat)
    assert result["valid"] is True
    assert result["n_noninteger"] == 0
    assert result["n_negative"] == 0
    assert result["n_nonfinite"] == 0


def test_fractional_value_outside_sampled_window_is_caught():
    # 500 rows, all integer except a single fractional value deep in the matrix
    # (row 400) — well outside a 200-row sample window. classify_matrix can miss
    # it; the full validator must not.
    mat = np.ones((500, 4), dtype=float)
    mat[400, 2] = 3.5
    result = _validate()(mat)
    assert result["valid"] is False
    assert result["n_noninteger"] >= 1


def test_negative_value_anywhere_is_invalid():
    mat = np.ones((300, 3), dtype=float)
    mat[250, 0] = -1.0
    result = _validate()(mat)
    assert result["valid"] is False
    assert result["n_negative"] >= 1


def test_nan_is_invalid():
    mat = np.ones((300, 3), dtype=float)
    mat[275, 1] = np.nan
    result = _validate()(mat)
    assert result["valid"] is False
    assert result["n_nonfinite"] >= 1


def test_inf_is_invalid():
    mat = np.ones((10, 3), dtype=float)
    mat[5, 1] = np.inf
    result = _validate()(mat)
    assert result["valid"] is False
    assert result["n_nonfinite"] >= 1


def test_sparse_fractional_in_stored_data_is_caught():
    sparse = pytest.importorskip("scipy.sparse")
    dense = np.zeros((400, 5), dtype=float)
    dense[399, 4] = 2.7   # fractional, in a late row
    mat = sparse.csr_matrix(dense)
    result = _validate()(mat)
    assert result["valid"] is False
    assert result["n_noninteger"] >= 1


def test_pandas_frame_supported():
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame({"s1": [1.0, 2.0, 3.0], "s2": [4.0, 5.0, 6.5]})
    result = _validate()(df)
    assert result["valid"] is False
    assert result["n_noninteger"] >= 1


def test_empty_matrix_reports_no_violations():
    result = _validate()(np.empty((0, 0)))
    # Nothing to reject; the caller handles empty separately.
    assert result["n_noninteger"] == 0
    assert result["n_negative"] == 0
    assert result["n_nonfinite"] == 0


# ── End-to-end: bulk_rna_de refuses instead of silently rounding ──────────────

def _unsampled_row(n_rows: int) -> int:
    """A 1-based gene row index that the seeded classifier sample does NOT visit,
    so the pre-fix code would round it silently and only the B6 full validation
    catches it."""
    from aria.utils.count_classifier import sample_row_indices

    sampled = set(int(i) for i in sample_row_indices(n_rows))
    for zero_based in range(n_rows - 1, -1, -1):
        if zero_based not in sampled:
            return zero_based + 1  # 1-based gene id in the TSV
    raise AssertionError("every row is sampled; enlarge n_rows")


def _write_counts_with_hidden_fractional(path, target_row_1based, n_rows):
    lines = ["gene\ts1\ts2\ts3\ts4"]
    for i in range(1, n_rows + 1):
        vals = [str(10 + i), str(12 + i), str(11 + i), str(9 + i)]
        if i == target_row_1based:
            vals[2] = "11.5"   # fractional value in a deterministically-unsampled row
        lines.append(f"g{i}\t" + "\t".join(vals))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_bulk_rna_de_refuses_late_fractional_instead_of_rounding(tmp_path):
    from aria.scripts.rna_bulk_de import bulk_rna_de

    n_rows = 500
    hidden = _unsampled_row(n_rows)
    counts_path = tmp_path / "counts.tsv"
    _write_counts_with_hidden_fractional(counts_path, hidden, n_rows)
    meta_path = tmp_path / "meta.tsv"
    meta_path.write_text(
        "sample\tcondition\n"
        "s1\tctrl\ns2\tctrl\ns3\ttreat\ns4\ttreat\n"
    )
    result = bulk_rna_de({
        "files": [str(counts_path)],
        "metadata_file": str(meta_path),
        "design_factor": "condition",
        "comparison": {"numerator": "treat", "denominator": "ctrl"},
        "output_dir": str(tmp_path / "out"),
        "run_pathways": False,
    })
    # Must NOT silently round a fractional value into a raw-count DE run.
    assert result["status"] == "error"
    assert result["error_type"] in {"NonRawCounts", "InvalidRawCounts"}
