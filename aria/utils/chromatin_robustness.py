"""scATAC robustness helpers for doublets, batch signals, and peak provenance.

These helpers are deliberately conservative and dependency-light where possible.
They warn or filter only from measured matrix/metadata signals; they never assign
cell identity or infer biology from peak names.
"""

from __future__ import annotations

import math
import re
from typing import Any


def _as_1d(values) -> list[float]:
    try:
        import numpy as np

        arr = np.asarray(values, dtype=float).ravel()
        return [float(x) for x in arr]
    except Exception:
        return [float(x) for x in values]


def _row_sums_and_features(counts) -> tuple[list[float], list[float]]:
    try:
        import numpy as np
        from scipy import sparse

        if sparse.issparse(counts):
            x = counts.tocsr()
            totals = np.asarray(x.sum(axis=1)).ravel()
            features = np.diff(x.indptr)
            return _as_1d(totals), _as_1d(features)
        arr = np.asarray(counts)
        return _as_1d(arr.sum(axis=1)), _as_1d((arr > 0).sum(axis=1))
    except Exception:
        totals = []
        features = []
        for row in counts:
            row_vals = list(row)
            totals.append(float(sum(row_vals)))
            features.append(float(sum(1 for v in row_vals if v > 0)))
        return totals, features


def _robust_z(values: list[float]) -> list[float]:
    import statistics

    if not values:
        return []
    med = statistics.median(values)
    abs_dev = [abs(v - med) for v in values]
    mad = statistics.median(abs_dev)
    scale = 1.4826 * mad
    if scale <= 0:
        try:
            scale = statistics.pstdev(values)
        except statistics.StatisticsError:
            scale = 0.0
    if scale <= 0:
        return [0.0 for _ in values]
    return [(v - med) / scale for v in values]


def detect_atac_doublets(
    counts,
    *,
    min_cells: int = 50,
    depth_z: float = 2.5,
    feature_z: float = 2.5,
    max_rate: float = 0.12,
) -> dict[str, Any]:
    """Conservative scATAC doublet detector from depth + accessible features.

    This is an ATAC-specific outlier screen, not a cell identity model. A cell is
    flagged only when both total fragments and the number of accessible peaks are
    high robust outliers. The final call is capped by ``max_rate`` to avoid
    aggressive filtering on skewed datasets; the highest-scoring cells are kept
    when too many exceed thresholds.
    """
    totals, features = _row_sums_and_features(counts)
    n = len(totals)
    if n < int(min_cells):
        return {
            "ran": False,
            "reason": f"only {n} cells; doublet detection needs >= {min_cells}",
            "method": "robust_depth_feature_outlier",
            "n_cells": n,
            "n_doublets": 0,
            "doublet_rate": 0.0,
            "_doublet_mask": [False] * n,
            "_doublet_score": [0.0] * n,
        }

    log_depth = [math.log1p(max(v, 0.0)) for v in totals]
    log_features = [math.log1p(max(v, 0.0)) for v in features]
    z_depth = _robust_z(log_depth)
    z_features = _robust_z(log_features)
    scores = [
        max(0.0, zd) + max(0.0, zf)
        for zd, zf in zip(z_depth, z_features)
    ]
    candidates = [
        zd >= depth_z and zf >= feature_z
        for zd, zf in zip(z_depth, z_features)
    ]

    max_n = max(1, int(math.ceil(max(0.0, min(float(max_rate), 1.0)) * n)))
    ranked = sorted(range(n), key=lambda i: scores[i], reverse=True)
    allowed = set(ranked[:max_n])
    mask = [bool(candidates[i] and i in allowed) for i in range(n)]
    n_doublets = int(sum(mask))

    finite_scores = [s for s in scores if math.isfinite(s)]
    score_max = max(finite_scores) if finite_scores else 0.0
    return {
        "ran": True,
        "reason": None,
        "method": "robust_depth_feature_outlier",
        "n_cells": n,
        "n_doublets": n_doublets,
        "doublet_rate": round(n_doublets / n, 4),
        "removed": n_doublets,
        "thresholds": {
            "depth_z": float(depth_z),
            "feature_z": float(feature_z),
            "max_rate": float(max_rate),
        },
        "metrics": {
            "median_fragments": _median(totals),
            "median_accessible_peaks": _median(features),
            "max_doublet_score": round(float(score_max), 4),
        },
        "_doublet_mask": mask,
        "_doublet_score": [round(float(s), 6) for s in scores],
    }


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    import statistics

    return round(float(statistics.median(values)), 4)


def public_doublet_summary(result: dict[str, Any], *, removed: int | None = None) -> dict[str, Any]:
    """Drop private per-cell arrays before serializing script output."""
    out = {k: v for k, v in result.items() if not str(k).startswith("_")}
    if removed is not None:
        out["removed"] = int(removed)
    return out


def assess_atac_batch_embedding(
    *,
    obs_columns: list[str] | None,
    embedding=None,
    cluster_labels=None,
    batch_labels=None,
    condition_col: str | None = None,
    replicate_col: str | None = None,
    declared_batch: str | None = None,
    batch_high: float = 0.10,
    cluster_low: float = 0.0,
) -> dict[str, Any]:
    """Warn about hidden batch columns and batch-dominated LSI embeddings."""
    from aria.utils.batch_qc import assess_hidden_batch
    from aria.utils.integration_qc import assess_integration_quality

    hidden = assess_hidden_batch(
        obs_columns,
        condition_col=condition_col,
        replicate_col=replicate_col,
        declared_batch=declared_batch,
        integration_ran=False,
    )
    issues = list(hidden.get("issues") or [])
    metrics: dict[str, Any] = {
        "candidate_batch_columns": hidden.get("candidate_batch_columns", []),
        "batch_col_evaluated": declared_batch,
        "batch_silhouette": None,
        "cluster_silhouette": None,
    }

    batch_sil = _safe_silhouette(embedding, batch_labels)
    cluster_sil = _safe_silhouette(embedding, cluster_labels)
    metrics["batch_silhouette"] = batch_sil
    metrics["cluster_silhouette"] = cluster_sil

    iqc = assess_integration_quality(
        None,
        batch_sil,
        cluster_sil,
        batch_high=batch_high,
        cluster_low=cluster_low,
    )
    issues.extend(iqc.get("issues") or [])

    return {
        "status": "warnings" if issues else "clean",
        "issues": issues,
        "metrics": metrics,
    }


def _safe_silhouette(embedding, labels) -> float | None:
    if embedding is None or labels is None:
        return None
    try:
        import numpy as np
        from sklearn.metrics import silhouette_score

        x = np.asarray(embedding)
        y = np.asarray(labels).astype(str)
        valid = np.array([v not in {"", "nan", "None"} for v in y])
        x = x[valid]
        y = y[valid]
        if x.shape[0] < 3 or len(set(y.tolist())) < 2:
            return None
        if len(set(y.tolist())) >= x.shape[0]:
            return None
        return round(float(silhouette_score(x, y)), 4)
    except Exception:
        return None


_COORD_RE = re.compile(r"^(chr)?[A-Za-z0-9_.]+[:-]\d+[-:]\d+$")


def assess_consensus_peak_provenance(
    peak_names,
    *,
    metadata: dict[str, Any] | None = None,
    input_kind: str | None = None,
) -> dict[str, Any]:
    """Assess whether the peak universe carries a reproducible consensus story."""
    if peak_names is None:
        names = []
    else:
        names = [str(p) for p in list(peak_names)]
    n = len(names)
    coord_like = sum(1 for p in names if _COORD_RE.match(p))
    dup = n - len(set(names))
    coord_fraction = round(coord_like / n, 4) if n else 0.0
    metadata = metadata or {}

    method = str(metadata.get("method") or metadata.get("source") or "").lower()
    n_samples = _to_int(metadata.get("n_samples") or metadata.get("samples"))
    reproducibility = (
        metadata.get("idr")
        or metadata.get("overlap_fraction")
        or metadata.get("reproducibility")
    )
    rare_policy = metadata.get("rare_peak_policy")

    verified = (
        method in {"consensus", "overlap_unified", "idr", "reproducible_consensus"}
        and (n_samples or 0) >= 2
        and reproducibility is not None
    )
    if verified:
        status = "verified"
        reason = None
    elif metadata:
        status = "partial"
        reason = (
            "Peak provenance metadata is present but lacks a full reproducibility "
            "statement (method, n_samples>=2, and overlap/IDR metric)."
        )
    else:
        status = "unverified"
        reason = (
            "No per-sample/condition consensus-peak provenance was provided; "
            "this is a preprocessed peak matrix, so ARIA cannot verify peak "
            "reproducibility or rare peak preservation from the matrix alone."
        )

    issues = []
    if dup:
        issues.append({
            "severity": "warning",
            "check": "duplicate_peak_names",
            "message": f"{dup} duplicate peak name(s) are present.",
            "recommendation": "Use a unique coordinate peak universe before DA.",
        })
    if coord_fraction < 0.95:
        issues.append({
            "severity": "warning",
            "check": "non_coordinate_peak_names",
            "message": (
                f"Only {coord_fraction:.1%} of peaks look like genomic "
                "coordinates."
            ),
            "recommendation": (
                "Use coordinate-style peak identifiers so downstream overlap, "
                "motif, and provenance checks are auditable."
            ),
        })
    if status != "verified":
        issues.append({
            "severity": "warning",
            "check": "consensus_peak_provenance_unverified",
            "message": reason or "Consensus peak provenance is incomplete.",
            "recommendation": (
                "Provide per-sample peak-calling/consensus metadata, including "
                "method, samples, reproducibility metric, and rare-peak policy."
            ),
        })

    return {
        "status": status,
        "reason": reason,
        "n_peaks": n,
        "coordinate_peak_fraction": coord_fraction,
        "n_duplicate_peak_names": dup,
        "input_kind": input_kind,
        "metadata": {
            "method": metadata.get("method") or metadata.get("source"),
            "n_samples": n_samples,
            "reproducibility": reproducibility,
            "rare_peak_policy": rare_policy,
        },
        "issues": issues,
    }


def _to_int(value) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
