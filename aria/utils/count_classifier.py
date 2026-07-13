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
RAW_SCORE_THRESHOLD = 0.75

RAW_TOOL_HINTS = (
    "counts", "count", "featurecounts", "htseq", "raw", "readcount",
)
NONRAW_TOOL_HINTS = (
    "tpm", "fpkm", "cpm", "rpkm", "expected_count", "expected.count",
    "expected-count", "est_counts", "abundance", "quant.sf", "normalized",
    "normalised", "lognorm", "log_norm", "log-normalized", "log-normalised",
)


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


def _finite_values(sample: np.ndarray) -> np.ndarray:
    return np.asarray(sample)[np.isfinite(sample)]


def _integer_fraction(values: np.ndarray) -> float:
    if values.size == 0:
        return 0.0
    if np.issubdtype(values.dtype, np.integer):
        return 1.0
    return float((np.abs(values - np.round(values)) <= 1e-6).mean())


def _is_integer_like(sample: np.ndarray) -> bool:
    """True when every finite sampled value is (numerically) an integer.

    Mirrors the historical pseudobulk ``np.allclose(sample, round(sample))``
    check so genuine raw counts stored as float still classify as raw.
    """
    if np.issubdtype(sample.dtype, np.integer):
        return True
    finite = _finite_values(sample)
    if finite.size == 0:
        return False
    return bool(np.allclose(finite, np.round(finite)))


def _nonnegative_fraction(values: np.ndarray) -> float:
    if values.size == 0:
        return 0.0
    return float((values >= 0).mean())


def _library_size_score(sample: np.ndarray) -> tuple[float, dict]:
    arr = np.asarray(sample, dtype=float)
    if arr.size == 0:
        return 0.0, {"reason": "empty"}
    if arr.ndim == 1:
        arr = arr[:, None]
    finite_rows = np.isfinite(arr).all(axis=1)
    arr = arr[finite_rows]
    if arr.size == 0:
        return 0.0, {"reason": "no_finite_rows"}
    lib_sizes = arr.sum(axis=0)
    finite = lib_sizes[np.isfinite(lib_sizes)]
    positive = finite[finite > 0]
    if positive.size < 2:
        return 0.5 if positive.size == 1 else 0.0, {
            "n_positive_libraries": int(positive.size)
        }
    cv = float(np.std(positive) / max(np.mean(positive), 1e-12))
    # Raw count libraries should have positive depth. Very high CV may be real
    # biology or library-size imbalance, so it reduces confidence but does not
    # veto rawness by itself.
    score = 1.0 if cv <= 1.5 else max(0.25, 1.0 - (cv - 1.5) / 3.0)
    return float(min(1.0, max(0.0, score))), {
        "n_positive_libraries": int(positive.size),
        "library_size_cv": round(cv, 4),
        "min_library_size": float(np.min(positive)),
        "max_library_size": float(np.max(positive)),
    }


def _gene_id_type_score(gene_ids) -> tuple[float, dict]:
    if gene_ids is None:
        return 0.5, {"gene_id_type": "unknown"}
    ids = [str(g) for g in list(gene_ids)[:500] if str(g)]
    if not ids:
        return 0.5, {"gene_id_type": "unknown"}
    upper = [g.upper() for g in ids]
    ensembl = sum(g.startswith(("ENSG", "ENSMUSG", "ENS", "FBGN", "WBG")) for g in upper)
    symbols = sum(g.replace("-", "").replace(".", "").isalnum() for g in ids)
    numeric = sum(g.isdigit() for g in ids)
    if ensembl / len(ids) >= 0.2:
        return 1.0, {"gene_id_type": "ensembl_like", "n_checked": len(ids)}
    if symbols / len(ids) >= 0.8 and numeric / len(ids) < 0.5:
        return 0.85, {"gene_id_type": "symbol_like", "n_checked": len(ids)}
    return 0.4, {"gene_id_type": "weak_or_numeric", "n_checked": len(ids)}


def _tool_signature_score(source_hint: str | None) -> tuple[float, dict]:
    text = str(source_hint or "").lower()
    if not text:
        return 0.5, {"tool_signature": "unknown"}
    if any(h in text for h in NONRAW_TOOL_HINTS):
        return 0.0, {"tool_signature": "normalized_or_expected_count"}
    if any(h in text for h in RAW_TOOL_HINTS):
        return 1.0, {"tool_signature": "raw_count_hint"}
    return 0.5, {"tool_signature": "unknown"}


def _confidence(score: float, sub_scores: dict[str, float]) -> str:
    if score >= 0.9 and min(sub_scores.values()) >= 0.75:
        return "high"
    if score >= RAW_SCORE_THRESHOLD:
        return "medium"
    if score >= 0.5:
        return "low"
    return "insufficient"


def _matrix_values(mat) -> np.ndarray:
    """Flat float view of every value to inspect, without densifying a sparse
    input. Sparse matrices are validated on their stored entries (implicit zeros
    are valid non-negative integers); pandas/numpy densify as they already sit in
    memory."""
    # scipy sparse: stored nonzeros live in ``.data``; avoid a dense copy.
    if hasattr(mat, "tocoo") and hasattr(mat, "data") and not hasattr(mat, "iloc"):
        return np.asarray(mat.data, dtype=float).ravel()
    if hasattr(mat, "to_numpy"):          # pandas DataFrame/Series
        return np.asarray(mat.to_numpy(), dtype=float).ravel()
    return np.asarray(mat, dtype=float).ravel()


def validate_raw_count_matrix(mat, max_examples: int = 5) -> dict:
    """Full vectorized check that EVERY value is a finite, non-negative integer.

    Unlike :func:`classify_matrix` — a fast heuristic that scores a <=200-row
    random sample — this scans the entire matrix, so a single fractional,
    negative, NaN, or inf value ANYWHERE fails the check. It is the deterministic
    gate the DE entry points run before rounding a matrix for DESeq2, so a matrix
    cannot be accepted as raw counts on the strength of a sampled slice.

    Returns ``{"valid", "reason", "n_nonfinite", "n_negative", "n_noninteger",
    "n_offending", "examples"}``. An empty matrix has no violations (``valid`` is
    True); callers handle emptiness separately.
    """
    values = _matrix_values(mat)
    if values.size == 0:
        return {"valid": True, "reason": None, "n_nonfinite": 0,
                "n_negative": 0, "n_noninteger": 0, "n_offending": 0,
                "examples": []}

    finite = np.isfinite(values)
    n_nonfinite = int((~finite).sum())
    fin = values[finite]
    negative_mask = fin < 0
    noninteger_mask = fin != np.round(fin)
    n_negative = int(negative_mask.sum())
    n_noninteger = int(noninteger_mask.sum())
    n_offending = n_nonfinite + n_negative + n_noninteger
    valid = n_offending == 0

    reasons = []
    examples: list[float] = []
    if n_nonfinite:
        reasons.append(f"{n_nonfinite} non-finite (NaN/inf)")
    if n_negative:
        reasons.append(f"{n_negative} negative")
        examples.extend(float(v) for v in fin[negative_mask][:max_examples])
    if n_noninteger:
        reasons.append(f"{n_noninteger} non-integer")
        examples.extend(float(v) for v in fin[noninteger_mask][:max_examples])

    return {
        "valid": valid,
        "reason": None if valid else "; ".join(reasons),
        "n_nonfinite": n_nonfinite,
        "n_negative": n_negative,
        "n_noninteger": n_noninteger,
        "n_offending": n_offending,
        "examples": examples[:max_examples],
    }


def classify_matrix(mat, seed: int = DEFAULT_SEED, gene_ids=None,
                    source_hint: str | None = None) -> dict:
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
      - ``raw_count_score`` and ``sub_scores``: deterministic evidence scores
        for integer-ness, non-negativity, library-size plausibility,
        decimal-fraction, gene-ID type, and tool/signature hints.

    The boolean decision is intentionally conservative: raw requires a high
    score, integer-like sampled values, and non-negative values. This removes the
    old hard dependency on ``max > 50`` so low-depth raw matrices can pass, while
    expected-count / TPM-like decimal matrices remain non-raw.
    """
    sample = sample_rows(mat, seed=seed)
    if sample.size == 0:
        return {"kind": "empty", "max": 0.0, "min": 0.0,
                "is_raw_counts": False, "integer_like": False,
                "raw_count_score": 0.0, "confidence": "insufficient",
                "sub_scores": {
                    "integer": 0.0,
                    "nonnegative": 0.0,
                    "library_size": 0.0,
                    "decimal_fraction": 0.0,
                    "gene_id_type": 0.5 if gene_ids is not None else 0.0,
                    "tool_signature": 0.5 if source_hint else 0.0,
                },
                "score_basis": {"empty": True}}

    max_val = float(np.nanmax(sample))
    min_val = float(np.nanmin(sample))
    finite = _finite_values(sample)
    integer_like = _is_integer_like(sample)
    integer_fraction = _integer_fraction(finite)
    nonnegative_fraction = _nonnegative_fraction(finite)
    non_negative = nonnegative_fraction >= 0.999
    decimal_fraction_score = 1.0 - min(1.0, max(0.0, 1.0 - integer_fraction))
    lib_score, lib_basis = _library_size_score(sample)
    gene_score, gene_basis = _gene_id_type_score(gene_ids)
    tool_score, tool_basis = _tool_signature_score(source_hint)
    sub_scores = {
        "integer": round(integer_fraction, 3),
        "nonnegative": round(nonnegative_fraction, 3),
        "library_size": round(lib_score, 3),
        "decimal_fraction": round(decimal_fraction_score, 3),
        "gene_id_type": round(gene_score, 3),
        "tool_signature": round(tool_score, 3),
    }
    weights = {
        "integer": 0.34,
        "nonnegative": 0.22,
        "library_size": 0.12,
        "decimal_fraction": 0.16,
        "gene_id_type": 0.08,
        "tool_signature": 0.08,
    }
    raw_count_score = round(
        sum(sub_scores[k] * weights[k] for k in weights), 3
    )
    source_is_nonraw = tool_score == 0.0
    is_raw = bool(
        raw_count_score >= RAW_SCORE_THRESHOLD
        and integer_like
        and non_negative
        and not source_is_nonraw
    )

    if is_raw:
        kind = "raw"
    elif not non_negative:
        kind = "scaled"
    elif integer_like:
        kind = "integer_small"
    elif source_is_nonraw and "expected" in str(source_hint or "").lower():
        kind = "expected_count"
    elif max_val <= 50:
        kind = "lognorm"
    else:
        kind = "continuous"

    return {"kind": kind, "max": max_val, "min": min_val,
            "is_raw_counts": is_raw, "integer_like": integer_like,
            "raw_count_score": raw_count_score,
            "confidence": _confidence(raw_count_score, sub_scores),
            "sub_scores": sub_scores,
            "score_basis": {
                "threshold": RAW_SCORE_THRESHOLD,
                "integer_fraction": round(integer_fraction, 4),
                "nonnegative_fraction": round(nonnegative_fraction, 4),
                "sample_dtype": str(sample.dtype),
                "source_hint": source_hint,
                **lib_basis,
                **gene_basis,
                **tool_basis,
            }}
