"""Tests for the shared raw-count classifier (audit 2026-05-29, P-RAWCLASS).

Covers B10 (bulk silently rounds non-raw matrices) and R7 (the probe must not
depend on row order). Runs with numpy only — no scientific stack required.
"""

import numpy as np

from aria.utils.count_classifier import (
    classify_matrix,
    sample_row_indices,
    sample_rows,
)


def _raw_counts(n=500, m=8, seed=1):
    rng = np.random.default_rng(seed)
    return rng.poisson(40, size=(n, m)).astype(float)  # max well above 50


def test_raw_integer_counts_classified_as_raw():
    info = classify_matrix(_raw_counts())
    assert info["is_raw_counts"] is True
    assert info["kind"] == "raw"


def test_integer_dtype_small_max_is_raw():
    mat = np.array([[0, 1, 2], [3, 4, 5]], dtype=int)
    info = classify_matrix(mat)
    assert info["is_raw_counts"] is True  # integer dtype is trusted


def test_lognormalized_matrix_is_not_raw():
    rng = np.random.default_rng(2)
    mat = rng.uniform(0, 8, size=(400, 6))  # non-integer, bounded
    info = classify_matrix(mat)
    assert info["is_raw_counts"] is False
    assert info["kind"] == "lognorm"


def test_tpm_like_continuous_matrix_is_not_raw():
    rng = np.random.default_rng(3)
    mat = rng.uniform(0, 5000, size=(400, 6))  # non-integer, large
    info = classify_matrix(mat)
    assert info["is_raw_counts"] is False
    assert info["kind"] == "continuous"


def test_scaled_matrix_with_negatives_is_not_raw():
    rng = np.random.default_rng(4)
    mat = rng.normal(0, 1, size=(400, 6))  # z-scored: has negatives
    info = classify_matrix(mat)
    assert info["is_raw_counts"] is False
    assert info["kind"] == "scaled"


def test_sampling_is_deterministic_given_seed():
    a = sample_row_indices(10_000, n_rows=200, seed=0)
    b = sample_row_indices(10_000, n_rows=200, seed=0)
    assert np.array_equal(a, b)


def test_classification_is_order_independent_R7():
    """R7: a matrix whose first rows are all-zero (e.g. ordered by cell type)
    must still classify on its real content, not the head."""
    raw = _raw_counts(n=1000)
    blanked = raw.copy()
    blanked[:200] = 0  # the old first-200 probe would see only zeros here
    info = classify_matrix(blanked)
    assert info["is_raw_counts"] is True


def test_sample_rows_densifies_only_the_slice():
    mat = _raw_counts(n=1000)
    block = sample_rows(mat, n_rows=50)
    assert block.shape == (50, mat.shape[1])
