"""V48 Slice I2 — External concordance scoring for Weighted Nearest Neighbors.

Pure, dependency-light (NumPy only) scoring layer comparing ARIA's paired WNN
joint-state outputs (`aria/scripts/integration_wnn.py`) against an external
multimodal reference implementation (e.g., Seurat v5 `FindMultiModalNeighbors`).

Metrics evaluated on deterministically aligned shared barcodes:
  1. Joint cluster concordance: Adjusted Rand Index (ARI) and Normalized Mutual
     Information (NMI).
  2. Per-cell modality weight concordance: Pearson r, Spearman rho, mean/median
     absolute difference, and dominant-modality (`w_rna >= 0.5`) agreement.
  3. Weighted kNN neighborhood overlap: per-cell Jaccard and recall over
     non-self neighbor sets.
  4. Seed determinism / stability: pairwise ARI and maximum weight deviation
     across repeated ARIA runs.
  5. Optional reference-partition alignment (neutral label recovery, ADR-011):
     ARI and NMI of ARIA WNN, external reference WNN, ARIA RNA-only, and ARIA
     ATAC-only partitions against an external cell-level partition column when
     present, without interpreting biological label strings.

Honesty contract (ADR-002): returns `status="not_run"` with an explicit reason
when external reference outputs are absent; never fabricates concordance values.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np

from aria.benchmarks.concordance_atac import (
    adjusted_rand_index,
    normalized_mutual_info,
)


def _rankdata_average(values: np.ndarray) -> np.ndarray:
    """Average-tie rank transform in pure NumPy (1-based ranks)."""
    arr = np.asarray(values, dtype=np.float64).ravel()
    n = arr.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    order = np.argsort(arr, kind="mergesort")
    sorted_vals = arr[order]
    ranks = np.empty(n, dtype=np.float64)

    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_vals[j] == sorted_vals[i]:
            j += 1
        avg_rank = 0.5 * (i + (j - 1)) + 1.0
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def pearson_correlation(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Pearson correlation over finite paired values; None when undefined."""
    xa = np.asarray(list(a), dtype=np.float64).ravel()
    xb = np.asarray(list(b), dtype=np.float64).ravel()
    if xa.size != xb.size or xa.size < 2:
        return None
    mask = np.isfinite(xa) & np.isfinite(xb)
    if int(mask.sum()) < 2:
        return None
    xa = xa[mask]
    xb = xb[mask]
    sa = float(np.std(xa))
    sb = float(np.std(xb))
    if sa == 0.0 or sb == 0.0:
        return None
    r_val = float(np.corrcoef(xa, xb)[0, 1])
    return r_val if np.isfinite(r_val) else None


def spearman_correlation(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Spearman rank correlation in pure NumPy; None when constant or <2 items."""
    xa = np.asarray(list(a), dtype=np.float64).ravel()
    xb = np.asarray(list(b), dtype=np.float64).ravel()
    if xa.size != xb.size or xa.size < 2:
        return None
    mask = np.isfinite(xa) & np.isfinite(xb)
    if int(mask.sum()) < 2:
        return None
    ra = _rankdata_average(xa[mask])
    rb = _rankdata_average(xb[mask])
    return pearson_correlation(ra, rb)


def modality_weight_concordance(
    aria_rna_weights: Sequence[float],
    ref_rna_weights: Sequence[float],
) -> dict[str, Any]:
    """Compare per-cell RNA modality weights on shared cells."""
    wa = np.asarray(list(aria_rna_weights), dtype=np.float64).ravel()
    wr = np.asarray(list(ref_rna_weights), dtype=np.float64).ravel()
    if wa.size != wr.size:
        raise ValueError(
            f"Weight vectors must match in length: {wa.size} != {wr.size}"
        )
    mask = np.isfinite(wa) & np.isfinite(wr)
    n_valid = int(mask.sum())
    if n_valid == 0:
        return {
            "n_cells_compared": 0,
            "pearson_r": None,
            "spearman_rho": None,
            "mae": None,
            "median_abs_diff": None,
            "dominant_modality_agreement": None,
            "aria_mean_rna_weight": None,
            "ref_mean_rna_weight": None,
        }

    wa_v = wa[mask]
    wr_v = wr[mask]
    abs_diff = np.abs(wa_v - wr_v)
    dom_agree = float(np.mean((wa_v >= 0.5) == (wr_v >= 0.5)))
    p_r = pearson_correlation(wa_v, wr_v)
    s_rho = spearman_correlation(wa_v, wr_v)

    return {
        "n_cells_compared": n_valid,
        "pearson_r": round(float(p_r), 4) if p_r is not None else None,
        "spearman_rho": round(float(s_rho), 4) if s_rho is not None else None,
        "mae": round(float(np.mean(abs_diff)), 4),
        "median_abs_diff": round(float(np.median(abs_diff)), 4),
        "dominant_modality_agreement": round(dom_agree, 4),
        "aria_mean_rna_weight": round(float(np.mean(wa_v)), 4),
        "ref_mean_rna_weight": round(float(np.mean(wr_v)), 4),
    }


def knn_jaccard_overlap(
    aria_neighbors: Sequence[Sequence[Any]],
    ref_neighbors: Sequence[Sequence[Any]],
    *,
    cell_ids: Sequence[Any] | None = None,
    match_effective_k: bool = True,
) -> dict[str, Any]:
    """Compute per-cell weighted kNN Jaccard and recall on non-self neighbors.

    When `match_effective_k=True`, if one implementation returns `k` non-self
    neighbors and the other returns `k - 1` (because self occupied slot 1),
    both distance-ordered neighbor lists are truncated to `min(len(a), len(b))`
    per cell so Jaccard is not biased by self-inclusion conventions.
    """
    if len(aria_neighbors) != len(ref_neighbors):
        raise ValueError(
            f"Neighbor lists must have equal length: "
            f"{len(aria_neighbors)} != {len(ref_neighbors)}"
        )
    n = len(aria_neighbors)
    if n == 0:
        return {
            "n_cells_compared": 0,
            "mean_jaccard": None,
            "median_jaccard": None,
            "mean_recall": None,
        }

    jaccards: list[float] = []
    recalls: list[float] = []
    for idx in range(n):
        self_id = cell_ids[idx] if cell_ids is not None else idx
        na_list = [x for x in aria_neighbors[idx] if x != self_id]
        nr_list = [x for x in ref_neighbors[idx] if x != self_id]
        if match_effective_k and na_list and nr_list:
            k_eval = min(len(na_list), len(nr_list))
            na_list = na_list[:k_eval]
            nr_list = nr_list[:k_eval]
        sa = set(na_list)
        sr = set(nr_list)
        if not sa and not sr:
            continue
        inter = len(sa & sr)
        union = len(sa | sr)
        jaccards.append(float(inter) / float(union) if union else 0.0)
        if sr:
            recalls.append(float(inter) / float(len(sr)))

    if not jaccards:
        return {
            "n_cells_compared": 0,
            "mean_jaccard": None,
            "median_jaccard": None,
            "mean_recall": None,
        }

    return {
        "n_cells_compared": int(len(jaccards)),
        "mean_jaccard": round(float(np.mean(jaccards)), 4),
        "median_jaccard": round(float(np.median(jaccards)), 4),
        "mean_recall": (
            round(float(np.mean(recalls)), 4) if recalls else None
        ),
    }


def score_seed_stability(
    cluster_runs: Sequence[Sequence[Any]],
    weight_runs: Sequence[Sequence[float]] | None = None,
) -> dict[str, Any]:
    """Score determinism / stability across repeated ARIA WNN executions."""
    n_runs = len(cluster_runs)
    if n_runs < 2:
        return {
            "status": "not_run",
            "n_runs": n_runs,
            "mean_pairwise_ari": None,
            "min_pairwise_ari": None,
            "max_weight_abs_diff": None,
        }

    aris: list[float] = []
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            aris.append(adjusted_rand_index(cluster_runs[i], cluster_runs[j]))

    max_w_diff: float | None = None
    if weight_runs is not None and len(weight_runs) >= 2:
        diffs: list[float] = []
        for i in range(len(weight_runs)):
            wi = np.asarray(list(weight_runs[i]), dtype=np.float64)
            for j in range(i + 1, len(weight_runs)):
                wj = np.asarray(list(weight_runs[j]), dtype=np.float64)
                if wi.size == wj.size and wi.size > 0:
                    diffs.append(float(np.max(np.abs(wi - wj))))
        if diffs:
            max_w_diff = round(float(max(diffs)), 8)

    return {
        "status": "success",
        "n_runs": int(n_runs),
        "mean_pairwise_ari": round(float(np.mean(aris)), 6),
        "min_pairwise_ari": round(float(np.min(aris)), 6),
        "max_weight_abs_diff": max_w_diff,
    }


def score_wnn_concordance(
    *,
    aria_barcodes: Sequence[str],
    aria_clusters: Sequence[Any],
    aria_rna_weights: Sequence[float] | None = None,
    ref_barcodes: Sequence[str] | None = None,
    ref_clusters: Sequence[Any] | None = None,
    ref_rna_weights: Sequence[float] | None = None,
    aria_knn: Mapping[str, Sequence[str]] | None = None,
    ref_knn: Mapping[str, Sequence[str]] | None = None,
    aria_rna_only_clusters: Sequence[Any] | None = None,
    aria_atac_only_clusters: Sequence[Any] | None = None,
    reference_labels: Sequence[Any] | None = None,
    seed_cluster_runs: Sequence[Sequence[Any]] | None = None,
    seed_weight_runs: Sequence[Sequence[float]] | None = None,
) -> dict[str, Any]:
    """Score ARIA WNN outputs against an external multimodal WNN reference.

    Returns `status="not_run"` when external reference clusters/barcodes are
    absent (ADR-002). Neutral labels only (ADR-011).
    """
    if not ref_barcodes or ref_clusters is None or len(ref_clusters) == 0:
        return {
            "status": "not_run",
            "reason": "External reference WNN outputs were not provided.",
            "seed_stability": (
                score_seed_stability(seed_cluster_runs, seed_weight_runs)
                if seed_cluster_runs is not None
                else {"status": "not_run", "n_runs": 0}
            ),
        }

    aria_bc = [str(x) for x in aria_barcodes]
    ref_bc = [str(x) for x in ref_barcodes]
    if len(aria_bc) != len(aria_clusters):
        raise ValueError("aria_barcodes and aria_clusters lengths do not match.")
    if len(ref_bc) != len(ref_clusters):
        raise ValueError("ref_barcodes and ref_clusters lengths do not match.")

    aria_idx = {bc: i for i, bc in enumerate(aria_bc)}
    ref_idx = {bc: i for i, bc in enumerate(ref_bc)}
    shared_bc = [bc for bc in aria_bc if bc in ref_idx]

    if not shared_bc:
        return {
            "status": "not_run",
            "reason": "Zero shared cell barcodes between ARIA and reference WNN.",
            "n_aria_cells": len(aria_bc),
            "n_ref_cells": len(ref_bc),
            "n_shared_cells": 0,
        }

    a_take = [aria_idx[bc] for bc in shared_bc]
    r_take = [ref_idx[bc] for bc in shared_bc]

    a_cl = [str(aria_clusters[i]) for i in a_take]
    r_cl = [str(ref_clusters[i]) for i in r_take]

    cluster_ari = adjusted_rand_index(a_cl, r_cl)
    cluster_nmi = normalized_mutual_info(a_cl, r_cl)

    weights_block: dict[str, Any] = {"status": "not_run"}
    if aria_rna_weights is not None and ref_rna_weights is not None:
        a_w = [float(aria_rna_weights[i]) for i in a_take]
        r_w = [float(ref_rna_weights[i]) for i in r_take]
        weights_block = {
            "status": "success",
            **modality_weight_concordance(a_w, r_w),
        }

    knn_block: dict[str, Any] = {"status": "not_run"}
    if aria_knn is not None and ref_knn is not None:
        shared_set = set(shared_bc)
        a_lists = [
            [str(nb) for nb in aria_knn.get(bc, ()) if str(nb) in shared_set]
            for bc in shared_bc
        ]
        r_lists = [
            [str(nb) for nb in ref_knn.get(bc, ()) if str(nb) in shared_set]
            for bc in shared_bc
        ]
        knn_block = {
            "status": "success",
            **knn_jaccard_overlap(
                a_lists, r_lists, cell_ids=shared_bc, match_effective_k=True
            ),
        }

    label_block: dict[str, Any] = {"status": "not_run"}
    if reference_labels is not None and len(reference_labels) == len(aria_bc):
        raw_lbls = [reference_labels[i] for i in a_take]
        valid_pos = [
            idx for idx, val in enumerate(raw_lbls)
            if val is not None and str(val).strip() not in ("", "nan", "None")
        ]
        if len(valid_pos) >= 2:
            gt = [str(raw_lbls[idx]) for idx in valid_pos]
            a_sub = [a_cl[idx] for idx in valid_pos]
            r_sub = [r_cl[idx] for idx in valid_pos]
            entry: dict[str, Any] = {
                "status": "success",
                "n_labeled_cells": int(len(valid_pos)),
                "n_unique_labels": int(len(set(gt))),
                "aria_wnn_ari": round(adjusted_rand_index(a_sub, gt), 4),
                "aria_wnn_nmi": round(normalized_mutual_info(a_sub, gt), 4),
                "ref_wnn_ari": round(adjusted_rand_index(r_sub, gt), 4),
                "ref_wnn_nmi": round(normalized_mutual_info(r_sub, gt), 4),
            }
            if (
                aria_rna_only_clusters is not None
                and len(aria_rna_only_clusters) == len(aria_bc)
            ):
                rna_all = [str(aria_rna_only_clusters[i]) for i in a_take]
                rna_sub = [rna_all[idx] for idx in valid_pos]
                entry["aria_rna_only_ari"] = round(
                    adjusted_rand_index(rna_sub, gt), 4
                )
                entry["aria_rna_only_nmi"] = round(
                    normalized_mutual_info(rna_sub, gt), 4
                )
            if (
                aria_atac_only_clusters is not None
                and len(aria_atac_only_clusters) == len(aria_bc)
            ):
                atac_all = [str(aria_atac_only_clusters[i]) for i in a_take]
                atac_sub = [atac_all[idx] for idx in valid_pos]
                entry["aria_atac_only_ari"] = round(
                    adjusted_rand_index(atac_sub, gt), 4
                )
                entry["aria_atac_only_nmi"] = round(
                    normalized_mutual_info(atac_sub, gt), 4
                )
            label_block = entry

    seed_block = (
        score_seed_stability(seed_cluster_runs, seed_weight_runs)
        if seed_cluster_runs is not None
        else {"status": "not_run", "n_runs": 0}
    )

    return {
        "status": "success",
        "n_aria_cells": int(len(aria_bc)),
        "n_ref_cells": int(len(ref_bc)),
        "n_shared_cells": int(len(shared_bc)),
        "n_aria_clusters": int(len(set(a_cl))),
        "n_ref_clusters": int(len(set(r_cl))),
        "cluster_ari": round(float(cluster_ari), 4),
        "cluster_nmi": round(float(cluster_nmi), 4),
        "modality_weights": weights_block,
        "knn_overlap": knn_block,
        "label_recovery": label_block,
        "seed_stability": seed_block,
    }
