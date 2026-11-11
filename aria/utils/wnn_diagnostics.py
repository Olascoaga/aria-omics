"""
ARIA WNN Anti-Degradation Guardrail ("Lens, Not Replacement")
--------------------------------------------------------------
Deterministic, pure NumPy/SciPy diagnostics for paired same-cell WNN joint
embeddings (V48 Slice I3, ADR-064 Principle #2, ADR-066).

Validates five quantitative axes without any dataset-specific assumptions or
hardcoded biological labels (ADR-005 / ADR-011) and without fabricating missing
measurements (ADR-002):

1. Per-modality neighborhood preservation (chance-adjusted kNN recall against
   exact Euclidean neighborhoods in unit-L2 RNA and ATAC latent spaces,
   distinguishing downweighting of uninformative assays from collapse of
   informative modality structure).
2. Rare-population retention (size-normalized local kNN cohesion adjusted for
   hypergeometric random-mixing expectation across rare single-modality clusters
   and optional external audit labels).
3. Spurious intermediate / artificial bridge-state screen (orphan WNN edge
   rate outside native 2k neighborhoods and cross-subpopulation native bridge
   density on bimodally impure joint clusters).
4. Per-cluster modality weight audit (mean, median, IQR, dominance fractions,
   and cluster-level modality collapse flags).
5. Zero-leakage deterministic held-out check (stratified K-fold cross-modal
   prediction against a seeded permutation null, plus variance-balanced
   orthogonal dimension-split validation of both out-of-fold WNN and the
   supplied joint WNN graph).
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree


CRITICAL_ISSUE_CODES: frozenset[str] = frozenset(
    {
        "rna_structure_degraded",
        "atac_structure_degraded",
        "rare_population_degraded",
        "reference_rare_population_degraded",
        "spurious_bridge_candidate",
        "held_out_joint_degraded",
    }
)


def _l2_normalize_rows(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Row-wise unit L2 normalization with explicit zero-norm mask."""
    arr = np.asarray(mat, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D embedding array, got shape {arr.shape}.")
    norms = np.linalg.norm(arr, axis=1, keepdims=True)
    zero_mask = norms.ravel() <= 1e-12
    safe_norms = np.where(norms > 1e-12, norms, 1.0)
    return arr / safe_norms, zero_mask


def _extract_offdiag_neighbors(
    graph: sp.spmatrix | np.ndarray,
    *,
    is_connectivity: bool = False,
    max_k: int | None = None,
) -> list[np.ndarray]:
    """Extract deterministic off-diagonal neighbor index arrays per cell."""
    csr = graph.tocsr() if sp.issparse(graph) else sp.csr_matrix(np.asarray(graph))
    n_cells = csr.shape[0]
    if csr.shape != (n_cells, n_cells):
        raise ValueError(
            f"WNN graph must be square (n_cells, n_cells), got {csr.shape}."
        )

    neighbors: list[np.ndarray] = []
    for i in range(n_cells):
        start, end = csr.indptr[i], csr.indptr[i + 1]
        cols = csr.indices[start:end]
        vals = csr.data[start:end]
        mask = (cols != i) & (vals > 0 if is_connectivity else np.isfinite(vals))
        cols = cols[mask]
        vals = vals[mask]
        if cols.size > 1:
            if is_connectivity:
                order = np.lexsort((cols, -vals))
            else:
                order = np.lexsort((cols, vals))
            cols = cols[order]
        if max_k is not None and cols.size > max_k:
            cols = cols[:max_k]
        neighbors.append(cols.astype(np.int64, copy=False))
    return neighbors


def _exact_knn_indices(emb_l2: np.ndarray, k_neighbors: int) -> np.ndarray:
    """Return exact off-diagonal Euclidean kNN indices of shape (n, k_eff)."""
    n_cells = emb_l2.shape[0]
    k_eff = min(max(1, int(k_neighbors)), max(1, n_cells - 1))
    tree = cKDTree(emb_l2)
    dists, idxs = tree.query(emb_l2, k=min(n_cells, k_eff + 1))
    if idxs.ndim == 1:
        idxs = idxs[:, None]
        dists = dists[:, None]

    out = np.empty((n_cells, k_eff), dtype=np.int64)
    for i in range(n_cells):
        row_idx = idxs[i]
        row_dist = dists[i]
        keep = row_idx != i
        cand_idx = row_idx[keep]
        cand_dist = row_dist[keep]
        if cand_idx.size > 1:
            order = np.lexsort((cand_idx, cand_dist))
            cand_idx = cand_idx[order]
        out[i, :] = cand_idx[:k_eff]
    return out


def _chance_adjust(observed: float, chance: float) -> float | None:
    """Compute (observed - chance) / (1 - chance), or None if degenerate."""
    denom = 1.0 - float(chance)
    if denom <= 1e-9:
        return None
    return float((float(observed) - float(chance)) / denom)


def _stratified_kfold_indices(
    labels: Sequence[str],
    n_folds: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Deterministic cluster-stratified K-fold cell splits."""
    n_cells = len(labels)
    n_folds_eff = min(max(2, int(n_folds)), max(2, n_cells // 4))
    rng = np.random.default_rng(int(seed))

    by_label: dict[str, list[int]] = {}
    for idx, lbl in enumerate(labels):
        by_label.setdefault(str(lbl), []).append(idx)

    fold_buckets: list[list[int]] = [[] for _ in range(n_folds_eff)]
    cursor = 0
    for lbl in sorted(by_label.keys()):
        members = np.asarray(by_label[lbl], dtype=np.int64)
        perm = rng.permutation(members)
        for j, cell_idx in enumerate(perm):
            fold_buckets[(cursor + j) % n_folds_eff].append(int(cell_idx))
        cursor = (cursor + len(perm)) % n_folds_eff

    all_idx = np.arange(n_cells, dtype=np.int64)
    folds: list[tuple[np.ndarray, np.ndarray]] = []
    for f in range(n_folds_eff):
        test_idx = np.sort(np.asarray(fold_buckets[f], dtype=np.int64))
        if test_idx.size == 0 or test_idx.size >= n_cells:
            continue
        train_mask = np.ones(n_cells, dtype=bool)
        train_mask[test_idx] = False
        train_idx = all_idx[train_mask]
        folds.append((train_idx, test_idx))
    return folds


def _centered_r2(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    train_mean: np.ndarray,
) -> float | None:
    """Multivariate R^2 centered on the training-fold mean."""
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - train_mean[None, :]) ** 2))
    if ss_tot <= 1e-12:
        return None
    return float(1.0 - ss_res / ss_tot)


def _predict_from_neighbors(
    target_train: np.ndarray,
    nbr_local_idx: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (raw neighbor mean, L2-normalized neighbor mean) for query cells."""
    # target_train: (n_train, d), nbr_local_idx: (n_query, k)
    gathered = target_train[nbr_local_idx]  # (n_query, k, d)
    raw_mean = np.mean(gathered, axis=1)
    l2_mean, _ = _l2_normalize_rows(raw_mean)
    return raw_mean, l2_mean


def _variance_balanced_dim_split(
    n_dims: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Split ordered latent components into matched-spectrum build/eval halves."""
    rng = np.random.default_rng(int(seed))
    build_dims: list[int] = []
    eval_dims: list[int] = []
    for p in range(0, n_dims - 1, 2):
        if int(rng.integers(0, 2)) == 0:
            build_dims.append(p)
            eval_dims.append(p + 1)
        else:
            build_dims.append(p + 1)
            eval_dims.append(p)
    if n_dims % 2 == 1:
        build_dims.append(n_dims - 1)
    return np.asarray(build_dims, dtype=np.int64), np.asarray(eval_dims, dtype=np.int64)


def _oof_wnn_query_neighbors(
    rna_train: np.ndarray,
    atac_train: np.ndarray,
    rna_test: np.ndarray,
    atac_test: np.ndarray,
    k_oof: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Out-of-fold Hao et al. (2021) WNN neighbor selection strictly in train."""
    n_train = rna_train.shape[0]
    n_test = rna_test.shape[0]
    k_eff = min(max(1, int(k_oof)), n_train)
    k_cand = min(n_train, max(k_eff, 2 * k_eff))

    rna_tree = cKDTree(rna_train)
    atac_tree = cKDTree(atac_train)

    d_rna_cand, idx_rna_cand = rna_tree.query(rna_test, k=k_cand)
    d_atac_cand, idx_atac_cand = atac_tree.query(atac_test, k=k_cand)
    if idx_rna_cand.ndim == 1:
        idx_rna_cand = idx_rna_cand[:, None]
        d_rna_cand = d_rna_cand[:, None]
        idx_atac_cand = idx_atac_cand[:, None]
        d_atac_cand = d_atac_cand[:, None]

    nbrs_r = idx_rna_cand[:, :k_eff]
    nbrs_a = idx_atac_cand[:, :k_eff]

    _, r_hat_from_r = _predict_from_neighbors(rna_train, nbrs_r)
    _, r_hat_from_a = _predict_from_neighbors(rna_train, nbrs_a)
    _, a_hat_from_a = _predict_from_neighbors(atac_train, nbrs_a)
    _, a_hat_from_r = _predict_from_neighbors(atac_train, nbrs_r)

    sigma_r = np.maximum(np.mean(d_rna_cand, axis=1), 1e-6)
    sigma_a = np.maximum(np.mean(d_atac_cand, axis=1), 1e-6)

    dist_r_self = np.linalg.norm(rna_test - r_hat_from_r, axis=1)
    dist_r_cross = np.linalg.norm(rna_test - r_hat_from_a, axis=1)
    dist_a_self = np.linalg.norm(atac_test - a_hat_from_a, axis=1)
    dist_a_cross = np.linalg.norm(atac_test - a_hat_from_r, axis=1)

    # Affinity ratio in log space (Hao et al. 2021 Eq. 3-4)
    score_r = np.maximum(0.0, (dist_r_cross - dist_r_self) / sigma_r)
    score_a = np.maximum(0.0, (dist_a_cross - dist_a_self) / sigma_a)
    max_s = np.maximum(score_r, score_a)
    exp_r = np.exp(np.clip(score_r - max_s, -50.0, 0.0))
    exp_a = np.exp(np.clip(score_a - max_s, -50.0, 0.0))
    w_rna = exp_r / (exp_r + exp_a)
    w_atac = 1.0 - w_rna

    joint_nbrs = np.empty((n_test, k_eff), dtype=np.int64)
    for i in range(n_test):
        cand = np.unique(np.concatenate([idx_rna_cand[i], idx_atac_cand[i]]))
        dr = np.linalg.norm(rna_train[cand] - rna_test[i][None, :], axis=1)
        da = np.linalg.norm(atac_train[cand] - atac_test[i][None, :], axis=1)
        dr0 = float(d_rna_cand[i, 0])
        da0 = float(d_atac_cand[i, 0])
        sim = w_rna[i] * np.exp(-np.maximum(0.0, dr - dr0) / sigma_r[i]) + w_atac[
            i
        ] * np.exp(-np.maximum(0.0, da - da0) / sigma_a[i])
        order = np.lexsort((cand, -sim))
        joint_nbrs[i, :] = cand[order[:k_eff]]

    return joint_nbrs, w_rna


def _compute_held_out_check(
    rna_emb: np.ndarray,
    atac_emb: np.ndarray,
    rna_l2: np.ndarray,
    atac_l2: np.ndarray,
    wnn_neighbors: list[np.ndarray],
    wnn_clusters: Sequence[str],
    rna_weights: np.ndarray | None,
    *,
    k_eval: int,
    n_folds: int,
    seed: int,
) -> dict[str, Any]:
    """Axis 5: Stratified K-fold cross-modal & orthogonal split-dimension check."""
    n_cells = rna_l2.shape[0]
    if n_cells < 12:
        return {
            "evaluable": False,
            "reason": "fewer_than_12_cells",
            "n_cells": int(n_cells),
        }

    folds = _stratified_kfold_indices(wnn_clusters, n_folds=n_folds, seed=seed)
    if len(folds) < 2:
        return {
            "evaluable": False,
            "reason": "insufficient_valid_folds",
            "n_cells": int(n_cells),
        }

    # 5A: Zero-overlap cross-modal out-of-fold prediction
    r_from_a_cos_all: list[np.ndarray] = []
    r_from_a_null_all: list[np.ndarray] = []
    r_from_a_r2_folds: list[float] = []

    a_from_r_cos_all: list[np.ndarray] = []
    a_from_r_null_all: list[np.ndarray] = []
    a_from_r_r2_folds: list[float] = []

    for f_idx, (train_idx, test_idx) in enumerate(folds):
        k_oof = min(max(1, k_eval), max(1, len(train_idx) - 1))
        rng_null = np.random.default_rng(int(seed) + 10_007 + f_idx)
        perms_train = [rng_null.permutation(len(train_idx)) for _ in range(5)]

        atac_tree = cKDTree(atac_l2[train_idx])
        _, nbr_a = atac_tree.query(atac_l2[test_idx], k=k_oof)
        if nbr_a.ndim == 1:
            nbr_a = nbr_a[:, None]

        pred_r_raw, pred_r_l2 = _predict_from_neighbors(rna_l2[train_idx], nbr_a)
        null_r_draws = [
            np.sum(
                rna_l2[test_idx]
                * _predict_from_neighbors(rna_l2[train_idx][p], nbr_a)[1],
                axis=1,
            )
            for p in perms_train
        ]
        r_from_a_cos_all.append(np.sum(rna_l2[test_idx] * pred_r_l2, axis=1))
        r_from_a_null_all.append(np.mean(null_r_draws, axis=0))
        r2_ra = _centered_r2(
            rna_l2[test_idx], pred_r_raw, np.mean(rna_l2[train_idx], axis=0)
        )
        if r2_ra is not None:
            r_from_a_r2_folds.append(r2_ra)

        rna_tree = cKDTree(rna_l2[train_idx])
        _, nbr_r = rna_tree.query(rna_l2[test_idx], k=k_oof)
        if nbr_r.ndim == 1:
            nbr_r = nbr_r[:, None]

        pred_a_raw, pred_a_l2 = _predict_from_neighbors(atac_l2[train_idx], nbr_r)
        null_a_draws = [
            np.sum(
                atac_l2[test_idx]
                * _predict_from_neighbors(atac_l2[train_idx][p], nbr_r)[1],
                axis=1,
            )
            for p in perms_train
        ]
        a_from_r_cos_all.append(np.sum(atac_l2[test_idx] * pred_a_l2, axis=1))
        a_from_r_null_all.append(np.mean(null_a_draws, axis=0))
        r2_ar = _centered_r2(
            atac_l2[test_idx], pred_a_raw, np.mean(atac_l2[train_idx], axis=0)
        )
        if r2_ar is not None:
            a_from_r_r2_folds.append(r2_ar)

    r_from_a_cos = float(np.mean(np.concatenate(r_from_a_cos_all)))
    r_from_a_null = float(np.mean(np.concatenate(r_from_a_null_all)))
    r_from_a_excess = float(r_from_a_cos - r_from_a_null)

    a_from_r_cos = float(np.mean(np.concatenate(a_from_r_cos_all)))
    a_from_r_null = float(np.mean(np.concatenate(a_from_r_null_all)))
    a_from_r_excess = float(a_from_r_cos - a_from_r_null)

    cross_modal = {
        "evaluable": True,
        "n_folds": len(folds),
        "rna_from_atac_cosine": round(r_from_a_cos, 4),
        "rna_from_atac_null_cosine": round(r_from_a_null, 4),
        "rna_from_atac_excess_cosine": round(r_from_a_excess, 4),
        "rna_from_atac_r2": (
            round(float(np.mean(r_from_a_r2_folds)), 4)
            if r_from_a_r2_folds
            else None
        ),
        "atac_from_rna_cosine": round(a_from_r_cos, 4),
        "atac_from_rna_null_cosine": round(a_from_r_null, 4),
        "atac_from_rna_excess_cosine": round(a_from_r_excess, 4),
        "atac_from_rna_r2": (
            round(float(np.mean(a_from_r_r2_folds)), 4)
            if a_from_r_r2_folds
            else None
        ),
    }

    # 5B: Orthogonal variance-balanced dimension split (zero self-prediction leakage)
    d_rna = rna_emb.shape[1]
    d_atac = atac_emb.shape[1]
    if d_rna < 4 or d_atac < 4:
        split_dim: dict[str, Any] = {
            "evaluable": False,
            "reason": "requires_at_least_4_dims_per_modality",
        }
    else:
        r_build_idx, r_eval_idx = _variance_balanced_dim_split(d_rna, seed + 20_011)
        a_build_idx, a_eval_idx = _variance_balanced_dim_split(d_atac, seed + 30_029)

        r_build, _ = _l2_normalize_rows(rna_emb[:, r_build_idx])
        r_eval, _ = _l2_normalize_rows(rna_emb[:, r_eval_idx])
        a_build, _ = _l2_normalize_rows(atac_emb[:, a_build_idx])
        a_eval, _ = _l2_normalize_rows(atac_emb[:, a_eval_idx])

        oof_rna_r_cos: list[np.ndarray] = []
        oof_rna_a_cos: list[np.ndarray] = []
        oof_rna_null_r: list[np.ndarray] = []
        oof_atac_r_cos: list[np.ndarray] = []
        oof_atac_a_cos: list[np.ndarray] = []
        oof_atac_null_a: list[np.ndarray] = []
        oof_joint_r_cos: list[np.ndarray] = []
        oof_joint_a_cos: list[np.ndarray] = []
        oof_weights = np.full(n_cells, np.nan, dtype=np.float64)

        for f_idx, (train_idx, test_idx) in enumerate(folds):
            k_oof = min(max(1, k_eval), max(1, len(train_idx) - 1))
            rng_null = np.random.default_rng(int(seed) + 40_009 + f_idx)
            perms_train = [rng_null.permutation(len(train_idx)) for _ in range(5)]

            r_tree = cKDTree(r_build[train_idx])
            _, nbr_r = r_tree.query(r_build[test_idx], k=k_oof)
            if nbr_r.ndim == 1:
                nbr_r = nbr_r[:, None]

            a_tree = cKDTree(a_build[train_idx])
            _, nbr_a = a_tree.query(a_build[test_idx], k=k_oof)
            if nbr_a.ndim == 1:
                nbr_a = nbr_a[:, None]

            nbr_j, w_r_fold = _oof_wnn_query_neighbors(
                r_build[train_idx],
                a_build[train_idx],
                r_build[test_idx],
                a_build[test_idx],
                k_oof=k_oof,
            )
            oof_weights[test_idx] = w_r_fold

            _, pr_r = _predict_from_neighbors(r_eval[train_idx], nbr_r)
            _, pa_r = _predict_from_neighbors(a_eval[train_idx], nbr_r)
            pr_null_draws = [
                np.sum(
                    r_eval[test_idx]
                    * _predict_from_neighbors(r_eval[train_idx][p], nbr_r)[1],
                    axis=1,
                )
                for p in perms_train
            ]

            _, pr_a = _predict_from_neighbors(r_eval[train_idx], nbr_a)
            _, pa_a = _predict_from_neighbors(a_eval[train_idx], nbr_a)
            pa_null_draws = [
                np.sum(
                    a_eval[test_idx]
                    * _predict_from_neighbors(a_eval[train_idx][p], nbr_a)[1],
                    axis=1,
                )
                for p in perms_train
            ]

            _, pr_j = _predict_from_neighbors(r_eval[train_idx], nbr_j)
            _, pa_j = _predict_from_neighbors(a_eval[train_idx], nbr_j)

            oof_rna_r_cos.append(np.sum(r_eval[test_idx] * pr_r, axis=1))
            oof_rna_a_cos.append(np.sum(a_eval[test_idx] * pa_r, axis=1))
            oof_rna_null_r.append(np.mean(pr_null_draws, axis=0))

            oof_atac_r_cos.append(np.sum(r_eval[test_idx] * pr_a, axis=1))
            oof_atac_a_cos.append(np.sum(a_eval[test_idx] * pa_a, axis=1))
            oof_atac_null_a.append(np.mean(pa_null_draws, axis=0))

            oof_joint_r_cos.append(np.sum(r_eval[test_idx] * pr_j, axis=1))
            oof_joint_a_cos.append(np.sum(a_eval[test_idx] * pa_j, axis=1))

        rna_only_r = float(np.mean(np.concatenate(oof_rna_r_cos)))
        rna_only_a = float(np.mean(np.concatenate(oof_rna_a_cos)))
        rna_null_r = float(np.mean(np.concatenate(oof_rna_null_r)))

        atac_only_r = float(np.mean(np.concatenate(oof_atac_r_cos)))
        atac_only_a = float(np.mean(np.concatenate(oof_atac_a_cos)))
        atac_null_a = float(np.mean(np.concatenate(oof_atac_null_a)))

        joint_r = float(np.mean(np.concatenate(oof_joint_r_cos)))
        joint_a = float(np.mean(np.concatenate(oof_joint_a_cos)))

        # Evaluate actual supplied wnn_distances graph on held-out D_eval dimensions
        graph_r_cos_list: list[float] = []
        graph_a_cos_list: list[float] = []
        for i, nbrs in enumerate(wnn_neighbors):
            if nbrs.size == 0:
                continue
            pred_r = np.mean(r_eval[nbrs], axis=0, keepdims=True)
            pred_r_l2, _ = _l2_normalize_rows(pred_r)
            pred_a = np.mean(a_eval[nbrs], axis=0, keepdims=True)
            pred_a_l2, _ = _l2_normalize_rows(pred_a)
            graph_r_cos_list.append(float(np.dot(r_eval[i], pred_r_l2[0])))
            graph_a_cos_list.append(float(np.dot(a_eval[i], pred_a_l2[0])))

        graph_r_cos = float(np.mean(graph_r_cos_list)) if graph_r_cos_list else 0.0
        graph_a_cos = float(np.mean(graph_a_cos_list)) if graph_a_cos_list else 0.0
        graph_bimodal_cos = 0.5 * (graph_r_cos + graph_a_cos)

        rna_bimodal_cos = 0.5 * (rna_only_r + rna_only_a)
        atac_bimodal_cos = 0.5 * (atac_only_r + atac_only_a)
        joint_bimodal_cos = 0.5 * (joint_r + joint_a)

        weight_oof_mae: float | None = None
        if rna_weights is not None and np.all(np.isfinite(oof_weights)):
            weight_oof_mae = round(
                float(np.mean(np.abs(rna_weights - oof_weights))), 4
            )

        split_dim = {
            "evaluable": True,
            "n_rna_build_dims": int(r_build_idx.size),
            "n_rna_eval_dims": int(r_eval_idx.size),
            "n_atac_build_dims": int(a_build_idx.size),
            "n_atac_eval_dims": int(a_eval_idx.size),
            "rna_within_eval_cosine": round(rna_only_r, 4),
            "rna_within_null_cosine": round(rna_null_r, 4),
            "rna_within_excess_cosine": round(rna_only_r - rna_null_r, 4),
            "atac_within_eval_cosine": round(atac_only_a, 4),
            "atac_within_null_cosine": round(atac_null_a, 4),
            "atac_within_excess_cosine": round(atac_only_a - atac_null_a, 4),
            "rna_only_bimodal_eval_cosine": round(rna_bimodal_cos, 4),
            "atac_only_bimodal_eval_cosine": round(atac_bimodal_cos, 4),
            "oof_joint_rna_eval_cosine": round(joint_r, 4),
            "oof_joint_atac_eval_cosine": round(joint_a, 4),
            "oof_joint_bimodal_eval_cosine": round(joint_bimodal_cos, 4),
            "supplied_graph_rna_eval_cosine": round(graph_r_cos, 4),
            "supplied_graph_atac_eval_cosine": round(graph_a_cos, 4),
            "supplied_graph_bimodal_eval_cosine": round(graph_bimodal_cos, 4),
            "oof_vs_insample_rna_weight_mae": weight_oof_mae,
        }

    # Determine whether each modality carries genuine manifold structure
    rna_within_exc = float(
        split_dim.get("rna_within_excess_cosine", 0.0)
        if split_dim.get("evaluable")
        else 0.0
    )
    atac_within_exc = float(
        split_dim.get("atac_within_excess_cosine", 0.0)
        if split_dim.get("evaluable")
        else 0.0
    )
    rna_informative = bool(
        r_from_a_excess > 0.10 or a_from_r_excess > 0.10 or rna_within_exc > 0.12
    )
    atac_informative = bool(
        r_from_a_excess > 0.10 or a_from_r_excess > 0.10 or atac_within_exc > 0.12
    )

    joint_degraded = False
    if split_dim.get("evaluable"):
        best_single = max(
            float(split_dim["rna_only_bimodal_eval_cosine"]),
            float(split_dim["atac_only_bimodal_eval_cosine"]),
        )
        supplied_bi = float(split_dim["supplied_graph_bimodal_eval_cosine"])
        oof_bi = float(split_dim["oof_joint_bimodal_eval_cosine"])
        if best_single > 0.12 and (
            supplied_bi < best_single - 0.10 or oof_bi < best_single - 0.10
        ):
            joint_degraded = True

    return {
        "evaluable": True,
        "seed": int(seed),
        "rna_manifold_informative": rna_informative,
        "atac_manifold_informative": atac_informative,
        "joint_bimodal_degraded": joint_degraded,
        "cross_modal_oof": cross_modal,
        "split_dimension_oof": split_dim,
    }


def _evaluate_rare_groups(
    group_labels: Sequence[str],
    wnn_clusters: Sequence[str],
    wnn_neighbors: list[np.ndarray],
    native_neighbors_primary: np.ndarray,
    native_neighbors_secondary: np.ndarray | None,
    *,
    source_name: str,
    max_rare_count: int,
    min_rare_cells: int,
    min_evaluable_rare_cells: int,
) -> dict[str, Any]:
    """Axis 2 helper: audit local kNN cohesion retention of rare groups."""
    n_cells = len(group_labels)
    lbl_arr = np.asarray([str(x) for x in group_labels], dtype=object)
    wnn_arr = np.asarray([str(x) for x in wnn_clusters], dtype=object)

    unique_groups = sorted(set(lbl_arr.tolist()))
    groups_report: list[dict[str, Any]] = []
    n_degraded = 0
    n_diluted = 0
    n_retained = 0
    n_not_evaluable = 0

    for grp in unique_groups:
        idx = np.flatnonzero(lbl_arr == grp)
        n_grp = int(idx.size)
        if n_grp < min_rare_cells or n_grp > max_rare_count:
            continue

        grp_set = set(idx.tolist())
        wnn_coh_cells: list[float] = []
        wnn_null_cells: list[float] = []
        nat1_coh_cells: list[float] = []
        nat1_null_cells: list[float] = []
        nat2_coh_cells: list[float] = []
        nat2_null_cells: list[float] = []

        for cell_i in idx:
            nbr_w = wnn_neighbors[int(cell_i)]
            k_w = int(nbr_w.size)
            if k_w > 0:
                denom_w = min(k_w, n_grp - 1)
                hits_w = sum(1 for j in nbr_w if int(j) in grp_set and int(j) != int(cell_i))
                wnn_coh_cells.append(hits_w / denom_w)
                wnn_null_cells.append((k_w * (n_grp - 1) / max(1, n_cells - 1)) / denom_w)

            nbr_1 = native_neighbors_primary[int(cell_i)]
            k_1 = int(nbr_1.size)
            denom_1 = min(k_1, n_grp - 1)
            hits_1 = sum(1 for j in nbr_1 if int(j) in grp_set and int(j) != int(cell_i))
            nat1_coh_cells.append(hits_1 / denom_1)
            nat1_null_cells.append((k_1 * (n_grp - 1) / max(1, n_cells - 1)) / denom_1)

            if native_neighbors_secondary is not None:
                nbr_2 = native_neighbors_secondary[int(cell_i)]
                k_2 = int(nbr_2.size)
                denom_2 = min(k_2, n_grp - 1)
                hits_2 = sum(
                    1 for j in nbr_2 if int(j) in grp_set and int(j) != int(cell_i)
                )
                nat2_coh_cells.append(hits_2 / denom_2)
                nat2_null_cells.append(
                    (k_2 * (n_grp - 1) / max(1, n_cells - 1)) / denom_2
                )

        wnn_raw = float(np.mean(wnn_coh_cells)) if wnn_coh_cells else 0.0
        wnn_null = float(np.mean(wnn_null_cells)) if wnn_null_cells else 0.0
        wnn_adj = _chance_adjust(wnn_raw, wnn_null)

        nat1_raw = float(np.mean(nat1_coh_cells))
        nat1_null = float(np.mean(nat1_null_cells))
        nat1_adj = _chance_adjust(nat1_raw, nat1_null)

        if native_neighbors_secondary is not None and nat2_coh_cells:
            nat2_raw = float(np.mean(nat2_coh_cells))
            nat2_null = float(np.mean(nat2_null_cells))
            nat2_adj = _chance_adjust(nat2_raw, nat2_null)
            if (nat2_adj or -1.0) > (nat1_adj or -1.0):
                nat_raw, nat_null, nat_adj = nat2_raw, nat2_null, nat2_adj
            else:
                nat_raw, nat_null, nat_adj = nat1_raw, nat1_null, nat1_adj
        else:
            nat_raw, nat_null, nat_adj = nat1_raw, nat1_null, nat1_adj

        # Joint cluster concentration
        wnn_sub = wnn_arr[idx]
        vals, counts = np.unique(wnn_sub, return_counts=True)
        best_pos = int(np.argmax(counts))
        dom_cluster = str(vals[best_pos])
        dom_count = int(counts[best_pos])
        dom_share = float(dom_count / n_grp)
        dom_cluster_total = int(np.sum(wnn_arr == dom_cluster))
        purity_in_dom = float(dom_count / max(1, dom_cluster_total))

        entry: dict[str, Any] = {
            "source": source_name,
            "group_id": str(grp),
            "n_cells": n_grp,
            "fraction_of_cells": round(n_grp / n_cells, 5),
            "native_cohesion_raw": round(nat_raw, 4),
            "native_cohesion_null": round(nat_null, 4),
            "native_cohesion_adjusted": (
                round(nat_adj, 4) if nat_adj is not None else None
            ),
            "wnn_cohesion_raw": round(wnn_raw, 4),
            "wnn_cohesion_null": round(wnn_null, 4),
            "wnn_cohesion_adjusted": (
                round(wnn_adj, 4) if wnn_adj is not None else None
            ),
            "dominant_wnn_cluster": dom_cluster,
            "dominant_wnn_cluster_share": round(dom_share, 4),
            "group_purity_in_dominant_wnn_cluster": round(purity_in_dom, 4),
        }

        if n_grp < min_evaluable_rare_cells:
            entry["evaluable"] = False
            entry["status"] = "not_evaluable"
            entry["reason"] = "below_min_evaluable_rare_cells"
            entry["retention_ratio_adjusted"] = None
            n_not_evaluable += 1
        elif nat_adj is None or wnn_adj is None or nat_adj < 0.20:
            entry["evaluable"] = False
            entry["status"] = "not_evaluable"
            entry["reason"] = "weak_native_cohesion"
            entry["retention_ratio_adjusted"] = None
            n_not_evaluable += 1
        else:
            retention_ratio = float(wnn_adj / nat_adj)
            entry["evaluable"] = True
            entry["retention_ratio_adjusted"] = round(retention_ratio, 4)
            entry["merged_at_global_resolution"] = bool(
                purity_in_dom < 0.50 and dom_share >= 0.70 and retention_ratio >= 0.50
            )
            if nat_adj >= 0.30 and (wnn_adj < 0.20 or retention_ratio < 0.40):
                entry["status"] = "degraded"
                n_degraded += 1
            elif retention_ratio < 0.65:
                entry["status"] = "warn"
                n_diluted += 1
            else:
                entry["status"] = "retained"
                n_retained += 1

        groups_report.append(entry)

    return {
        "source": source_name,
        "n_rare_groups_detected": len(groups_report),
        "n_evaluable": len(groups_report) - n_not_evaluable,
        "n_retained": n_retained,
        "n_warn": n_diluted,
        "n_degraded": n_degraded,
        "n_not_evaluable": n_not_evaluable,
        "groups": groups_report,
    }


def _cross_subpop_bridge_density(
    idx_a: np.ndarray,
    idx_b: np.ndarray,
    native_knn: np.ndarray,
) -> float:
    """Symmetric fraction of native kNN edges crossing between two sub-groups."""
    if idx_a.size == 0 or idx_b.size == 0:
        return 0.0
    set_a = set(idx_a.tolist())
    set_b = set(idx_b.tolist())

    fracs_ab = [
        sum(1 for j in native_knn[int(i)] if int(j) in set_b)
        / max(1, native_knn[int(i)].size)
        for i in idx_a
    ]
    fracs_ba = [
        sum(1 for j in native_knn[int(i)] if int(j) in set_a)
        / max(1, native_knn[int(i)].size)
        for i in idx_b
    ]
    return float(0.5 * (np.mean(fracs_ab) + np.mean(fracs_ba)))


def compute_wnn_anti_degradation_diagnostics(
    *,
    rna_embedding: np.ndarray,
    atac_embedding: np.ndarray,
    wnn_distances: sp.spmatrix | np.ndarray,
    wnn_clusters: Sequence[str],
    rna_clusters: Sequence[str],
    atac_clusters: Sequence[str],
    rna_weights: np.ndarray | Sequence[float] | None = None,
    atac_weights: np.ndarray | Sequence[float] | None = None,
    audit_labels: Sequence[str] | None = None,
    doublet_scores: np.ndarray | Sequence[float] | None = None,
    graph_is_connectivity: bool = False,
    k: int = 20,
    seed: int = 0,
    rare_fraction_threshold: float = 0.05,
    min_rare_cells: int = 3,
    min_evaluable_rare_cells: int = 5,
    n_folds: int = 5,
) -> dict[str, Any]:
    """Compute the 5-axis WNN anti-degradation guardrail report."""
    rna_emb = np.asarray(rna_embedding, dtype=np.float64)
    atac_emb = np.asarray(atac_embedding, dtype=np.float64)
    n_cells = rna_emb.shape[0]
    if atac_emb.shape[0] != n_cells:
        raise ValueError(
            f"RNA ({n_cells}) and ATAC ({atac_emb.shape[0]}) cell counts differ."
        )
    if (
        len(wnn_clusters) != n_cells
        or len(rna_clusters) != n_cells
        or len(atac_clusters) != n_cells
    ):
        raise ValueError("Cluster label vectors must all have length n_cells.")
    if n_cells < 6:
        raise ValueError(f"At least 6 shared cells required, got {n_cells}.")

    rna_l2, rna_zero = _l2_normalize_rows(rna_emb)
    atac_l2, atac_zero = _l2_normalize_rows(atac_emb)

    wnn_nbrs = _extract_offdiag_neighbors(
        wnn_distances, is_connectivity=graph_is_connectivity
    )
    degrees = np.asarray([nbr.size for nbr in wnn_nbrs], dtype=np.int64)
    n_isolated = int(np.sum(degrees == 0))
    nonzero_deg = degrees[degrees > 0]
    k_med = int(np.median(nonzero_deg)) if nonzero_deg.size > 0 else max(1, min(k, n_cells - 1))
    k_max = int(np.max(nonzero_deg)) if nonzero_deg.size > 0 else k_med

    rna_w: np.ndarray | None = None
    atac_w: np.ndarray | None = None
    if rna_weights is not None and atac_weights is not None:
        rna_w = np.asarray(rna_weights, dtype=np.float64)
        atac_w = np.asarray(atac_weights, dtype=np.float64)
        if rna_w.shape != (n_cells,) or atac_w.shape != (n_cells,):
            raise ValueError("Modality weight arrays must have shape (n_cells,).")

    dbl_arr: np.ndarray | None = None
    if doublet_scores is not None:
        dbl_arr = np.asarray(doublet_scores, dtype=np.float64)
        if dbl_arr.shape != (n_cells,):
            dbl_arr = None

    issues: list[dict[str, Any]] = []
    if n_isolated > 0:
        issues.append(
            {
                "code": "isolated_wnn_cells",
                "severity": "warn",
                "detail": f"{n_isolated} cell(s) have zero off-diagonal WNN neighbors.",
            }
        )

    # ── Axis 5 first: independent manifold & held-out cross-prediction ────
    held_out = _compute_held_out_check(
        rna_emb=rna_emb,
        atac_emb=atac_emb,
        rna_l2=rna_l2,
        atac_l2=atac_l2,
        wnn_neighbors=wnn_nbrs,
        wnn_clusters=wnn_clusters,
        rna_weights=rna_w,
        k_eval=k_med,
        n_folds=n_folds,
        seed=seed,
    )
    if held_out.get("evaluable") and held_out.get("joint_bimodal_degraded"):
        issues.append(
            {
                "code": "held_out_joint_degraded",
                "severity": "degraded",
                "detail": (
                    "Joint WNN graph underperforms the best single modality on "
                    "held-out orthogonal bimodal latent reconstruction."
                ),
            }
        )

    # ── Axis 1: Per-modality neighborhood preservation ────────────────────
    rna_knn_max = _exact_knn_indices(rna_l2, k_max)
    atac_knn_max = _exact_knn_indices(atac_l2, k_max)

    r_recalls: list[float] = []
    a_recalls: list[float] = []
    u_recalls: list[float] = []
    p0_list: list[float] = []
    pu0_list: list[float] = []
    r_dom_recalls: list[float] = []
    a_dom_recalls: list[float] = []

    for i, w_nbr in enumerate(wnn_nbrs):
        ki = int(w_nbr.size)
        if ki == 0:
            continue
        w_set = set(w_nbr.tolist())
        r_set = set(rna_knn_max[i, :ki].tolist())
        a_set = set(atac_knn_max[i, :ki].tolist())

        rr = len(w_set & r_set) / ki
        ar = len(w_set & a_set) / ki
        ur = len(w_set & (r_set | a_set)) / ki
        p0 = ki / max(1, n_cells - 1)
        pu0 = 1.0 - (1.0 - p0) ** 2

        r_recalls.append(rr)
        a_recalls.append(ar)
        u_recalls.append(ur)
        p0_list.append(p0)
        pu0_list.append(pu0)

        if rna_w is not None and atac_w is not None:
            if rna_w[i] > 0.55:
                r_dom_recalls.append(rr)
            if atac_w[i] > 0.55:
                a_dom_recalls.append(ar)

    mean_p0 = float(np.mean(p0_list)) if p0_list else 1.0
    mean_pu0 = float(np.mean(pu0_list)) if pu0_list else 1.0
    mean_r_raw = float(np.mean(r_recalls)) if r_recalls else 0.0
    mean_a_raw = float(np.mean(a_recalls)) if a_recalls else 0.0
    mean_u_raw = float(np.mean(u_recalls)) if u_recalls else 0.0

    mean_r_adj = _chance_adjust(mean_r_raw, mean_p0)
    mean_a_adj = _chance_adjust(mean_a_raw, mean_p0)
    mean_u_adj = _chance_adjust(mean_u_raw, mean_pu0)

    rna_inf = bool(held_out.get("rna_manifold_informative", True))
    atac_inf = bool(held_out.get("atac_manifold_informative", True))
    mean_rw = float(np.mean(rna_w)) if rna_w is not None else None
    mean_aw = float(np.mean(atac_w)) if atac_w is not None else None

    if mean_r_adj is not None and mean_r_adj < 0.10:
        if rna_inf:
            issues.append(
                {
                    "code": "rna_structure_degraded",
                    "severity": "degraded",
                    "detail": (
                        f"RNA carries manifold signal but chance-adjusted WNN "
                        f"kNN recall collapsed to {mean_r_adj:.3f}."
                    ),
                }
            )
        else:
            issues.append(
                {
                    "code": "rna_uninformative_downweighted",
                    "severity": "warn",
                    "detail": (
                        "RNA shows near-null manifold predictability and low "
                        "neighborhood retention in joint WNN."
                    ),
                }
            )
    elif not rna_inf and mean_rw is not None and mean_rw < 0.25:
        issues.append(
            {
                "code": "rna_uninformative_downweighted",
                "severity": "warn",
                "detail": (
                    f"RNA lacks independent manifold signal and was downweighted "
                    f"(mean_rna_weight={mean_rw:.3f})."
                ),
            }
        )

    if mean_a_adj is not None and mean_a_adj < 0.10:
        if atac_inf:
            issues.append(
                {
                    "code": "atac_structure_degraded",
                    "severity": "degraded",
                    "detail": (
                        f"ATAC carries manifold signal but chance-adjusted WNN "
                        f"kNN recall collapsed to {mean_a_adj:.3f}."
                    ),
                }
            )
        else:
            issues.append(
                {
                    "code": "atac_uninformative_downweighted",
                    "severity": "warn",
                    "detail": (
                        "ATAC shows near-null manifold predictability and low "
                        "neighborhood retention in joint WNN."
                    ),
                }
            )
    elif not atac_inf and mean_aw is not None and mean_aw < 0.25:
        issues.append(
            {
                "code": "atac_uninformative_downweighted",
                "severity": "warn",
                "detail": (
                    f"ATAC lacks independent manifold signal and was downweighted "
                    f"(mean_atac_weight={mean_aw:.3f})."
                ),
            }
        )

    modality_preservation = {
        "evaluable": bool(mean_r_adj is not None and mean_p0 < 0.80),
        "median_wnn_degree": k_med,
        "n_isolated_cells": n_isolated,
        "n_zero_norm_rna_cells": int(np.sum(rna_zero)),
        "n_zero_norm_atac_cells": int(np.sum(atac_zero)),
        "chance_expected_single_modality_recall": round(mean_p0, 4),
        "chance_expected_union_recall": round(mean_pu0, 4),
        "mean_rna_recall_raw": round(mean_r_raw, 4),
        "mean_rna_recall_adjusted": (
            round(mean_r_adj, 4) if mean_r_adj is not None else None
        ),
        "mean_atac_recall_raw": round(mean_a_raw, 4),
        "mean_atac_recall_adjusted": (
            round(mean_a_adj, 4) if mean_a_adj is not None else None
        ),
        "mean_union_recall_raw": round(mean_u_raw, 4),
        "mean_union_recall_adjusted": (
            round(mean_u_adj, 4) if mean_u_adj is not None else None
        ),
        "rna_dominant_cells_rna_recall": (
            round(float(np.mean(r_dom_recalls)), 4) if r_dom_recalls else None
        ),
        "atac_dominant_cells_atac_recall": (
            round(float(np.mean(a_dom_recalls)), 4) if a_dom_recalls else None
        ),
    }

    # ── Axis 2: Rare-population retention ─────────────────────────────────
    base_max_rare = int(np.floor(float(rare_fraction_threshold) * n_cells))
    if n_cells < 200 and rare_fraction_threshold <= 0.10:
        max_rare_count = max(base_max_rare, min(15, int(np.floor(0.15 * n_cells))))
    else:
        max_rare_count = max(min_evaluable_rare_cells, base_max_rare)

    rna_knn_med = _exact_knn_indices(rna_l2, k_med)
    atac_knn_med = _exact_knn_indices(atac_l2, k_med)

    rare_rna = _evaluate_rare_groups(
        rna_clusters,
        wnn_clusters,
        wnn_nbrs,
        rna_knn_med,
        None,
        source_name="leiden_rna",
        max_rare_count=max_rare_count,
        min_rare_cells=min_rare_cells,
        min_evaluable_rare_cells=min_evaluable_rare_cells,
    )
    rare_atac = _evaluate_rare_groups(
        atac_clusters,
        wnn_clusters,
        wnn_nbrs,
        atac_knn_med,
        None,
        source_name="leiden_atac",
        max_rare_count=max_rare_count,
        min_rare_cells=min_rare_cells,
        min_evaluable_rare_cells=min_evaluable_rare_cells,
    )

    if rare_rna["n_degraded"] > 0 or rare_atac["n_degraded"] > 0:
        degraded_ids = [
            f"{g['source']}:{g['group_id']}"
            for g in [*rare_rna["groups"], *rare_atac["groups"]]
            if g.get("status") == "degraded"
        ]
        issues.append(
            {
                "code": "rare_population_degraded",
                "severity": "degraded",
                "detail": (
                    f"{len(degraded_ids)} rare per-modality population(s) lost "
                    f"local kNN cohesion in joint WNN: {', '.join(degraded_ids)}."
                ),
            }
        )
    elif rare_rna["n_warn"] > 0 or rare_atac["n_warn"] > 0:
        warn_ids = [
            f"{g['source']}:{g['group_id']}"
            for g in [*rare_rna["groups"], *rare_atac["groups"]]
            if g.get("status") == "warn"
        ]
        issues.append(
            {
                "code": "rare_population_partially_diluted",
                "severity": "warn",
                "detail": (
                    f"{len(warn_ids)} rare per-modality population(s) showed "
                    f"partial kNN cohesion dilution: {', '.join(warn_ids)}."
                ),
            }
        )

    rare_reference: dict[str, Any] | None = None
    if audit_labels is not None:
        if len(audit_labels) != n_cells:
            raise ValueError("audit_labels must have length n_cells.")
        rare_reference = _evaluate_rare_groups(
            audit_labels,
            wnn_clusters,
            wnn_nbrs,
            rna_knn_med,
            atac_knn_med,
            source_name="audit_labels",
            max_rare_count=max_rare_count,
            min_rare_cells=min_rare_cells,
            min_evaluable_rare_cells=min_evaluable_rare_cells,
        )
        if rare_reference["n_degraded"] > 0:
            ref_ids = [
                str(g["group_id"])
                for g in rare_reference["groups"]
                if g.get("status") == "degraded"
            ]
            issues.append(
                {
                    "code": "reference_rare_population_degraded",
                    "severity": "degraded",
                    "detail": (
                        f"{len(ref_ids)} audited rare reference group(s) lost "
                        f"local kNN cohesion in joint WNN: {', '.join(ref_ids)}."
                    ),
                }
            )

    rare_population_retention: dict[str, Any] = {
        "max_rare_group_cells": int(max_rare_count),
        "min_rare_cells": int(min_rare_cells),
        "min_evaluable_rare_cells": int(min_evaluable_rare_cells),
        "rna_rare_clusters": rare_rna,
        "atac_rare_clusters": rare_atac,
    }
    if rare_reference is not None:
        rare_population_retention["reference_rare_groups"] = rare_reference

    # ── Axis 3: Spurious intermediate / artificial bridge screen ──────────
    k_anchor = min(2 * k_med, max(k_med, max(1, (n_cells - 1) // 5)))
    k_anchor = min(k_anchor, max(1, n_cells - 1))
    p_anchor_0 = 1.0 - (1.0 - (k_anchor / max(1, n_cells - 1))) ** 2
    anchor_evaluable = bool(p_anchor_0 < 0.75)

    rna_anchor_knn = _exact_knn_indices(rna_l2, k_anchor)
    atac_anchor_knn = _exact_knn_indices(atac_l2, k_anchor)

    orphan_frac_per_cell = np.zeros(n_cells, dtype=np.float64)
    for i, w_nbr in enumerate(wnn_nbrs):
        if w_nbr.size == 0:
            orphan_frac_per_cell[i] = 1.0
            continue
        native_union = set(rna_anchor_knn[i].tolist()) | set(atac_anchor_knn[i].tolist())
        orphans = sum(1 for j in w_nbr if int(j) not in native_union)
        orphan_frac_per_cell[i] = orphans / w_nbr.size

    wnn_arr = np.asarray([str(x) for x in wnn_clusters], dtype=object)
    rna_arr = np.asarray([str(x) for x in rna_clusters], dtype=object)
    atac_arr = np.asarray([str(x) for x in atac_clusters], dtype=object)

    bridge_clusters: list[dict[str, Any]] = []
    bridge_candidates: list[str] = []

    for cl in sorted(set(wnn_arr.tolist())):
        c_idx = np.flatnonzero(wnn_arr == cl)
        n_cl = int(c_idx.size)

        r_vals, r_cnts = np.unique(rna_arr[c_idx], return_counts=True)
        r_top_pos = int(np.argmax(r_cnts))
        r_modal = str(r_vals[r_top_pos])
        r_purity = float(r_cnts[r_top_pos] / max(1, n_cl))

        a_vals, a_cnts = np.unique(atac_arr[c_idx], return_counts=True)
        a_top_pos = int(np.argmax(a_cnts))
        a_modal = str(a_vals[a_top_pos])
        a_purity = float(a_cnts[a_top_pos] / max(1, n_cl))

        mean_orphan = float(np.mean(orphan_frac_per_cell[c_idx]))
        mean_anchor = 1.0 - mean_orphan
        anchor_adj = (
            _chance_adjust(mean_anchor, p_anchor_0) if anchor_evaluable else None
        )

        cl_entry: dict[str, Any] = {
            "joint_cluster": str(cl),
            "n_cells": n_cl,
            "modal_rna_cluster": r_modal,
            "modal_rna_purity": round(r_purity, 4),
            "modal_atac_cluster": a_modal,
            "modal_atac_purity": round(a_purity, 4),
            "mean_orphan_edge_fraction": round(mean_orphan, 4),
            "mean_union_anchor_adjusted": (
                round(anchor_adj, 4) if anchor_adj is not None else None
            ),
        }
        if dbl_arr is not None:
            cl_entry["mean_doublet_score"] = round(float(np.mean(dbl_arr[c_idx])), 4)

        if n_cl < 5:
            cl_entry["evaluable"] = False
            cl_entry["flagged_spurious_bridge"] = False
            cl_entry["reasons"] = []
            bridge_clusters.append(cl_entry)
            continue

        reasons: list[str] = []
        if anchor_evaluable and mean_orphan > 0.50 and (anchor_adj or 0.0) < 0.25:
            reasons.append("high_orphan_wnn_edges")

        if r_purity < 0.65 and a_purity < 0.65:
            sub_r1 = c_idx[rna_arr[c_idx] == r_modal]
            sub_r2 = c_idx[rna_arr[c_idx] != r_modal]
            rna_bridge_r = _cross_subpop_bridge_density(sub_r1, sub_r2, rna_knn_med)
            atac_bridge_r = _cross_subpop_bridge_density(sub_r1, sub_r2, atac_knn_med)

            sub_a1 = c_idx[atac_arr[c_idx] == a_modal]
            sub_a2 = c_idx[atac_arr[c_idx] != a_modal]
            rna_bridge_a = _cross_subpop_bridge_density(sub_a1, sub_a2, rna_knn_med)
            atac_bridge_a = _cross_subpop_bridge_density(sub_a1, sub_a2, atac_knn_med)

            cl_entry["rna_split_native_bridge_max"] = round(
                max(rna_bridge_r, atac_bridge_r), 4
            )
            cl_entry["atac_split_native_bridge_max"] = round(
                max(rna_bridge_a, atac_bridge_a), 4
            )

            if (
                max(rna_bridge_r, atac_bridge_r) < 0.05
                and max(rna_bridge_a, atac_bridge_a) < 0.05
            ):
                reasons.append("bimodal_mixture_without_native_knn_bridge")

        flagged = len(reasons) > 0
        cl_entry["evaluable"] = True
        cl_entry["flagged_spurious_bridge"] = flagged
        cl_entry["reasons"] = reasons
        if flagged:
            bridge_candidates.append(str(cl))
        bridge_clusters.append(cl_entry)

    if bridge_candidates:
        issues.append(
            {
                "code": "spurious_bridge_candidate",
                "severity": "degraded",
                "detail": (
                    f"{len(bridge_candidates)} joint cluster(s) lack native "
                    f"single-modality neighborhood support and may represent "
                    f"fabricated intermediate/bridge states: "
                    f"{', '.join(bridge_candidates)}."
                ),
            }
        )

    spurious_bridge_screen = {
        "anchor_k": int(k_anchor),
        "anchor_chance_expected": round(p_anchor_0, 4),
        "anchor_evaluable": anchor_evaluable,
        "global_mean_orphan_edge_fraction": round(
            float(np.mean(orphan_frac_per_cell)), 4
        ),
        "n_flagged_clusters": len(bridge_candidates),
        "flagged_clusters": bridge_candidates,
        "clusters": bridge_clusters,
    }

    # ── Axis 4: Per-cluster modality weight audit ─────────────────────────
    if rna_w is None or atac_w is None:
        issues.append(
            {
                "code": "missing_modality_weights",
                "severity": "warn",
                "detail": "Per-cell modality weights were not provided.",
            }
        )
        weight_audit: dict[str, Any] = {
            "evaluable": False,
            "reason": "missing_modality_weights",
            "clusters": {},
        }
    else:
        g_mean_rna = float(np.mean(rna_w))
        g_med_rna = float(np.median(rna_w))
        g_iqr_rna = float(np.percentile(rna_w, 75) - np.percentile(rna_w, 25))
        if g_mean_rna > 0.85 or g_mean_rna < 0.15:
            issues.append(
                {
                    "code": "global_modality_imbalance",
                    "severity": "warn",
                    "detail": (
                        f"Global mean RNA weight={g_mean_rna:.3f} indicates "
                        f"strong dominance of a single modality."
                    ),
                }
            )

        cluster_weights: dict[str, dict[str, Any]] = {}
        collapsed_clusters: list[str] = []
        for cl in sorted(set(wnn_arr.tolist())):
            c_idx = np.flatnonzero(wnn_arr == cl)
            n_cl = int(c_idx.size)
            rw_c = rna_w[c_idx]
            aw_c = atac_w[c_idx]
            m_rw = float(np.mean(rw_c))
            m_aw = float(np.mean(aw_c))
            med_rw = float(np.median(rw_c))
            iqr_rw = float(np.percentile(rw_c, 75) - np.percentile(rw_c, 25))
            frac_r_dom = float(np.mean(rw_c > 0.55))
            frac_a_dom = float(np.mean(aw_c > 0.55))
            frac_bal = float(np.mean(np.abs(rw_c - aw_c) <= 0.10))
            collapsed = bool(n_cl >= 5 and (m_rw > 0.90 or m_rw < 0.10))
            if collapsed:
                collapsed_clusters.append(str(cl))

            dom_label = (
                "rna"
                if m_rw > m_aw + 0.05
                else ("atac" if m_aw > m_rw + 0.05 else "balanced")
            )
            cluster_weights[str(cl)] = {
                "n_cells": n_cl,
                "mean_rna_weight": round(m_rw, 4),
                "median_rna_weight": round(med_rw, 4),
                "iqr_rna_weight": round(iqr_rw, 4),
                "mean_atac_weight": round(m_aw, 4),
                "fraction_rna_dominant": round(frac_r_dom, 4),
                "fraction_atac_dominant": round(frac_a_dom, 4),
                "fraction_balanced": round(frac_bal, 4),
                "dominant_modality": dom_label,
                "modality_collapsed": collapsed,
            }

        if collapsed_clusters:
            issues.append(
                {
                    "code": "cluster_modality_collapse",
                    "severity": "warn",
                    "detail": (
                        f"{len(collapsed_clusters)} joint cluster(s) have extreme "
                        f"single-modality weight (>0.90 or <0.10): "
                        f"{', '.join(collapsed_clusters)}."
                    ),
                }
            )

        weight_audit = {
            "evaluable": True,
            "global_mean_rna_weight": round(g_mean_rna, 4),
            "global_median_rna_weight": round(g_med_rna, 4),
            "global_iqr_rna_weight": round(g_iqr_rna, 4),
            "global_mean_atac_weight": round(float(np.mean(atac_w)), 4),
            "n_collapsed_clusters": len(collapsed_clusters),
            "collapsed_clusters": collapsed_clusters,
            "clusters": cluster_weights,
        }

    # ── Final Guardrail Status ────────────────────────────────────────────
    issue_codes = [str(it["code"]) for it in issues]
    has_critical = any(
        it.get("severity") == "degraded" or it.get("code") in CRITICAL_ISSUE_CODES
        for it in issues
    )
    if has_critical:
        guardrail_status = "degraded"
    elif issues:
        guardrail_status = "warn"
    else:
        guardrail_status = "pass"

    return {
        "guardrail_status": guardrail_status,
        "principle": "lens_not_replacement",
        "n_cells": int(n_cells),
        "seed": int(seed),
        "issue_codes": issue_codes,
        "issues": issues,
        "modality_preservation": modality_preservation,
        "rare_population_retention": rare_population_retention,
        "spurious_bridge_screen": spurious_bridge_screen,
        "per_cluster_weight_audit": weight_audit,
        "held_out_check": held_out,
    }


def sanitize_for_hdf5_uns(obj: Any) -> Any:
    """Convert diagnostic nested structures into AnnData/MuData HDF5-safe types.

    Drops keys whose value is ``None`` (since ``anndata>=0.11`` raises
    ``IORegistryError`` on ``None`` inside ``.uns``), replaces ``/`` in HDF5
    group keys, converts lists of dicts or heterogeneous sequences into
    string-indexed dicts (since ``h5py`` rejects ``object`` arrays of dicts),
    and casts NumPy scalars/arrays to native Python primitives.
    """
    if isinstance(obj, dict):
        cleaned: dict[str, Any] = {}
        for k, v in obj.items():
            if v is None:
                continue
            safe_key = str(k).replace("/", "_") or "_"
            cleaned[safe_key] = sanitize_for_hdf5_uns(v)
        return cleaned
    if isinstance(obj, np.ndarray):
        return sanitize_for_hdf5_uns(obj.tolist())
    if isinstance(obj, (list, tuple)):
        items = [sanitize_for_hdf5_uns(v) for v in obj if v is not None]
        if not items:
            return []
        if all(isinstance(v, bool) for v in items):
            return [bool(v) for v in items]
        if all(isinstance(v, int) and not isinstance(v, bool) for v in items):
            return [int(v) for v in items]
        if all(
            isinstance(v, (int, float)) and not isinstance(v, bool) for v in items
        ):
            return [float(v) for v in items]
        if all(isinstance(v, str) for v in items):
            return [str(v) for v in items]
        return {str(idx): val for idx, val in enumerate(items)}
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj
