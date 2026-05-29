"""Shared raw-count classifier for the DE entry points.

Audit 2026-05-29, P-RAWCLASS (closes B10 + R7). Both bulk and pseudobulk DE
must decide whether an input matrix is genuine raw counts (safe for DESeq2) or a
transformed matrix (TPM/CPM/FPKM/log-normalized/scaled) that would silently be
coerced into pseudo-counts. This is the single detector plus a deterministic
*random* row sampler.

R7: the previous pseudobulk probes read the first 200 rows (``mat[:200]``), which
is biased when the matrix is ordered by cell type or condition. Sampling is now
random but seeded, so it stays reproducible (``--reproducible`` friendly) while
no longer depending on row order.
"""

from __future__ import annotations

import numpy as np

DEFAULT_PROBE_ROWS = 200
DEFAULT_SEED = 0


def sample_row_indices(n_total: int,
                       n_rows: int = DEFAULT_PROBE_ROWS,
                       seed: int = DEFAULT_SEED) -> np.ndarray:
    """Return up to ``n_rows`` sorted, randomly chosen row indices in [0, n_total).

    Deterministic given ``seed``. Returns an empty array when ``n_total`` is 0.
    Sorting keeps sparse fancy-indexing efficient and the slice stable.
    """
    if n_total <= 0:
        return np.empty(0, dtype=int)
    k = int(min(n_rows, n_total))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_total, size=k, replace=False))


def sample_rows(mat,
                n_rows: int = DEFAULT_PROBE_ROWS,
                seed: int = DEFAULT_SEED) -> np.ndarray:
    """Densify and return up to ``n_rows`` randomly chosen rows of ``mat``.

    Only the sampled slice is densified; the full matrix is never materialized.
    ``mat`` may be a scipy sparse matrix or a dense array-like (pass
    ``DataFrame.values`` for pandas, so indexing selects rows not columns).
    """
    idx = sample_row_indices(mat.shape[0], n_rows=n_rows, seed=seed)
    if idx.size == 0:
        empty = mat[:0]
        return empty.toarray() if hasattr(empty, "toarray") else np.asarray(empty)
    block = mat[idx]
    return block.toarray() if hasattr(block, "toarray") else np.asarray(block)


def _is_integer_like(sample: np.ndarray) -> bool:
    """True when every finite sampled value is (numerically) an integer.

    Mirrors the historical pseudobulk ``np.allclose(sample, round(sample))``
    check so genuine raw counts stored as float still classify as raw.
    """
    if np.issubdtype(sample.dtype, np.integer):
        return True
    finite = sample[np.isfinite(sample)]
    if finite.size == 0:
        return False
    return bool(np.allclose(finite, np.round(finite)))


def classify_matrix(mat, seed: int = DEFAULT_SEED) -> dict:
    """Classify a counts-like matrix from a random sampled slice.

    Returns a dict with:
      - ``kind``: one of ``raw`` (non-negative integers with a large max),
        ``integer_small`` (non-negative integers but small max — ambiguous),
        ``scaled`` (contains negatives — z-scored/standardized),
        ``lognorm`` (non-integer, bounded — log1p-normalized / logCPM),
        ``continuous`` (non-integer, large max — TPM/CPM/FPKM), or ``empty``;
      - ``is_raw_counts``: bool — safe to hand to DESeq2 without coercion;
      - ``max`` / ``min``: sampled extremes (diagnostic);
      - ``integer_like``: bool.

    The ``is_raw_counts`` decision preserves the pre-existing pseudobulk
    semantics: integer dtype, or (non-negative AND integer-valued AND max > 50).
    """
    sample = sample_rows(mat, seed=seed)
    if sample.size == 0:
        return {"kind": "empty", "max": 0.0, "min": 0.0,
                "is_raw_counts": False, "integer_like": False}

    max_val = float(np.nanmax(sample))
    min_val = float(np.nanmin(sample))
    integer_like = _is_integer_like(sample)
    non_negative = min_val >= 0
    is_raw = bool(integer_like and non_negative and
                  (np.issubdtype(sample.dtype, np.integer) or max_val > 50))

    if is_raw:
        kind = "raw"
    elif not non_negative:
        kind = "scaled"
    elif integer_like:
        kind = "integer_small"
    elif max_val <= 50:
        kind = "lognorm"
    else:
        kind = "continuous"

    return {"kind": kind, "max": max_val, "min": min_val,
            "is_raw_counts": is_raw, "integer_like": integer_like}
