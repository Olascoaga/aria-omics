"""
ARIA MOFA+ Integration Script (V48 Slice I4 — Tier 2a Interpretable Factors)
----------------------------------------------------------------------------
Deterministic Multi-Omics Factor Analysis Plus (MOFA+, Argelaguet et al. 2020,
Genome Biology; implemented via pure-Python variational Bayes `mofapy2`) over
same-cell paired multiome data (`.h5mu` preferred, or aligned multi-modality
matrices). Executed inside `aria-integration-env` by `EnvironmentManager`.

Pipeline & Governance Guarantees (ADR-002, ADR-005, ADR-011, ADR-064):
  1. Verified raw-count ingestion via `read_h5mu_paired` / `_select_raw_counts_matrix`
     and canonical barcode ordering (never `set` iteration order, never double
     log-normalization).
  2. Deterministic stratified subsampling (`max_cells`, default 5000) with a
     rare-population floor (`min_stratum_cells`, default 5), inverse-probability
     inclusion weights, and exact closed-form Gaussian variational posterior
     projection (`Z = (sum_v Y_v T_v W_v)(I + sum_v W_v^T T_v W_v)^{-1}`) for
     out-of-subset cells (validated at `|r| >= 0.99` on training cells so
     `obsm["X_mofa"]` contains zero `NaN`s).
  3. Modality-appropriate preprocessing fitted strictly on the training subset:
     - RNA: `normalize_total(1e4)` -> `log1p` -> top `n_top_genes` HVGs.
     - ATAC: peak prevalence filter -> Signac log-TF-log-IDF -> top `n_top_peaks`
       variable peaks.
     - View variance scaling (`scale_views=True`) so differing feature counts do
       not bias the variational ELBO.
  4. Canonical factor ordering (descending positive R^2 sum across views),
     deterministic loading-anchored sign orientation, and hybrid absolute +
     view-relative factor classification (`shared`, `modality_specific`,
     `low_variance`).
  5. Measured-only confounder detection on the training subset (no embedded
     biological gene panels):
     - Sequencing depth (`|r| >= depth_corr_cutoff` against pre-filter
       `log10(1 + total_counts)` per modality).
     - Cell cycle (`|r| >= cell_cycle_corr_cutoff` against existing `.obs`
       `S_score`/`G2M_score` columns or caller-supplied `cell_cycle_genes`;
       honestly listed in `metrics_not_computed` when absent).
"""

from __future__ import annotations

import contextlib
import io
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.benchmarks.concordance_mofa import (
    is_factor_active_in_view,
    pearson_correlation,
    spearman_correlation,
)
from aria.scripts._base import mocks_allowed, run_script
from aria.utils.stage_manifest import (
    build_stage_manifest,
    stage_is_current,
    write_stage_manifest,
)
from aria.utils.wnn_diagnostics import sanitize_for_hdf5_uns

_STRATA_COL_PREFERENCE = (
    "leiden_rna",
    "leiden",
    "leiden_wnn",
    "leiden_atac",
)
_S_SCORE_COLS = ("S_score", "s_score")
_G2M_SCORE_COLS = ("G2M_score", "g2m_score")
_AUTO_BATCH_CANDIDATES = ("batch", "donor_id", "sample_id")
_AUTO_CONDITION_CANDIDATES = (
    "condition",
    "disease",
    "group",
    "treatment",
    "age_group",
    "phenotype",
)
_STAGE_VERSION = "aria.integration_mofa.v2"


def _expand_manifest_inputs(role_prefix: str, raw_path: str | None) -> list[tuple[str, str]]:
    """Expand a file or directory (e.g. 10x MEX dir) into deterministic file roles."""
    if not raw_path:
        return []
    p = Path(str(raw_path))
    if not p.exists():
        return []
    if p.is_file():
        return [(role_prefix, str(p.resolve()))]
    if p.is_dir():
        out: list[tuple[str, str]] = []
        for child in sorted(p.rglob("*")):
            if child.is_file():
                rel = child.relative_to(p).as_posix()
                out.append((f"{role_prefix}:{rel}", str(child.resolve())))
        return out
    return []


def _resolve_categorical_obs(
    obs_dfs: list[Any],
    explicit_col: str | None,
    fallback_candidates: tuple[str, ...],
    *,
    exclude_col: str | None = None,
) -> tuple[str | None, list[str] | None, str]:
    """Resolve a categorical column with >= 2 levels across aligned `.obs` frames."""
    if explicit_col is not None and str(explicit_col).strip():
        col = str(explicit_col).strip()
        for df in obs_dfs:
            if df is not None and col in df.columns:
                vals = [str(v) for v in df[col].tolist()]
                if len(set(vals)) >= 2 and len(set(vals)) < len(vals):
                    return col, vals, "explicit"
                return col, None, "explicit_insufficient_levels"
        return col, None, "explicit_missing"

    for col in fallback_candidates:
        if exclude_col and col == exclude_col:
            continue
        for df in obs_dfs:
            if df is not None and col in df.columns:
                vals = [str(v) for v in df[col].tolist()]
                uniq = set(vals)
                if len(uniq) >= 2 and len(uniq) < len(vals):
                    return col, vals, f"auto:{col}"
    return None, None, "absent"


def _compute_categorical_eta2(
    values,
    group_labels: list[str],
    *,
    weights=None,
) -> tuple[float | None, float | None]:
    """Compute one-way ANOVA eta^2 (weighted if weights given) and unweighted F-stat."""
    import numpy as np

    z = np.asarray(values, dtype=np.float64)
    n = int(z.size)
    if n < 3 or len(group_labels) != n:
        return None, None
    uniq = sorted(set(group_labels))
    b = len(uniq)
    if b < 2 or b >= n:
        return None, None
    if float(np.std(z)) <= 1e-12:
        return 0.0, 0.0

    w = (
        np.asarray(weights, dtype=np.float64)
        if weights is not None
        else np.ones(n, dtype=np.float64)
    )
    w_sum = float(np.sum(w))
    if w_sum <= 0.0:
        return None, None

    grand_mean = float(np.sum(w * z) / w_sum)
    ss_tot = float(np.sum(w * ((z - grand_mean) ** 2)))
    if ss_tot <= 1e-15:
        return 0.0, 0.0

    ss_between = 0.0
    # Unweighted ANOVA F-statistic alongside weighted eta^2
    mean_unw = float(np.mean(z))
    ss_tot_unw = float(np.sum((z - mean_unw) ** 2))
    ss_bet_unw = 0.0
    for g in uniq:
        mask = np.asarray([lbl == g for lbl in group_labels], dtype=bool)
        if not np.any(mask):
            continue
        w_g = float(np.sum(w[mask]))
        if w_g > 0.0:
            mean_g = float(np.sum(w[mask] * z[mask]) / w_g)
            ss_between += w_g * ((mean_g - grand_mean) ** 2)
        n_g = int(np.sum(mask))
        mean_g_unw = float(np.mean(z[mask]))
        ss_bet_unw += float(n_g) * ((mean_g_unw - mean_unw) ** 2)

    eta2 = max(0.0, min(1.0, ss_between / ss_tot))
    ss_within_unw = max(0.0, ss_tot_unw - ss_bet_unw)
    df_bet = b - 1
    df_with = n - b
    if df_bet > 0 and df_with > 0 and ss_within_unw > 1e-15:
        f_stat = (ss_bet_unw / float(df_bet)) / (ss_within_unw / float(df_with))
    else:
        f_stat = 0.0
    return round(float(eta2), 4), round(float(f_stat), 4)


def _classify_modality_kind(name: str) -> str:
    low = str(name).strip().lower()
    if any(k in low for k in ("rna", "gex", "expr")):
        return "rna"
    if any(k in low for k in ("atac", "peak", "chrom", "acc")):
        return "atac"
    return "other"


def _deterministic_subsample(
    barcodes: list[str],
    obs_dfs: list[Any],
    *,
    max_cells: int,
    min_stratum_cells: int,
    explicit_strata_col: str | None,
    seed: int,
) -> dict[str, Any]:
    """Select training cell indices deterministically, preserving rare strata."""
    import numpy as np

    n_total = len(barcodes)
    if max_cells <= 0 or n_total <= max_cells:
        return {
            "subsampled": False,
            "method": "none",
            "strata_col": None,
            "train_indices": np.arange(n_total, dtype=int),
            "inclusion_weights": np.ones(n_total, dtype=np.float64),
            "exceeded_max_cells_for_rare_floor": False,
        }

    chosen_col: str | None = None
    strata_values: list[str] | None = None

    candidates = (
        (explicit_strata_col,)
        if explicit_strata_col
        else _STRATA_COL_PREFERENCE
    )
    for col in candidates:
        if not col:
            continue
        for df in obs_dfs:
            if df is not None and col in df.columns:
                vals = [str(v) for v in df[col].values]
                if len(set(vals)) >= 2:
                    chosen_col = str(col)
                    strata_values = vals
                    break
        if chosen_col is not None:
            break

    rng = np.random.default_rng(int(seed))

    if strata_values is None:
        # Sort by barcode before uniform sampling so choice is row-order invariant.
        lex_order = np.argsort(np.asarray(barcodes, dtype=object), kind="mergesort")
        picked_lex = rng.choice(n_total, size=int(max_cells), replace=False)
        train_idx = np.sort(lex_order[picked_lex]).astype(int)
        weights = np.ones(n_total, dtype=np.float64)
        return {
            "subsampled": True,
            "method": "uniform",
            "strata_col": None,
            "train_indices": train_idx,
            "inclusion_weights": weights,
            "exceeded_max_cells_for_rare_floor": False,
        }

    buckets: dict[str, list[int]] = {}
    for i, lbl in enumerate(strata_values):
        buckets.setdefault(lbl, []).append(i)

    sorted_labels = sorted(buckets.keys())
    floors: dict[str, int] = {}
    quotas: dict[str, int] = {}
    for lbl in sorted_labels:
        n_s = len(buckets[lbl])
        fl = min(n_s, max(1, int(min_stratum_cells)))
        floors[lbl] = fl
        prop = int(round(float(max_cells) * float(n_s) / float(n_total)))
        quotas[lbl] = min(n_s, max(fl, prop))

    floor_sum = sum(floors.values())
    exceeded_floor = False
    if floor_sum >= int(max_cells):
        quotas = dict(floors)
        exceeded_floor = floor_sum > int(max_cells)
    else:
        # Rebalance quotas deterministically so sum(quotas) == max_cells.
        cur_sum = sum(quotas.values())
        while cur_sum > int(max_cells):
            reducible = [lbl for lbl in sorted_labels if quotas[lbl] > floors[lbl]]
            if not reducible:
                break
            reducible.sort(key=lambda l: (-(quotas[l] - floors[l]), -quotas[l], l))
            quotas[reducible[0]] -= 1
            cur_sum -= 1
        while cur_sum < int(max_cells):
            expandable = [lbl for lbl in sorted_labels if quotas[lbl] < len(buckets[lbl])]
            if not expandable:
                break
            expandable.sort(key=lambda l: (-(len(buckets[l]) - quotas[l]), l))
            quotas[expandable[0]] += 1
            cur_sum += 1

    picked_all: list[int] = []
    n_used = sum(quotas.values())
    raw_weights = np.ones(n_total, dtype=np.float64)

    for lbl in sorted_labels:
        idxs = buckets[lbl]
        idxs_sorted = sorted(idxs, key=lambda idx: barcodes[idx])
        q_s = quotas[lbl]
        if q_s >= len(idxs_sorted):
            chosen = idxs_sorted
        else:
            pos = rng.choice(len(idxs_sorted), size=q_s, replace=False)
            chosen = [idxs_sorted[int(p)] for p in pos]
        picked_all.extend(chosen)

        # Inverse inclusion probability weight (Horvitz-Thompson style), normalized below.
        pi_s = float(q_s) / float(len(idxs_sorted))
        w_s = (1.0 / pi_s) if pi_s > 0.0 else 1.0
        for idx in idxs_sorted:
            raw_weights[idx] = w_s

    train_idx = np.asarray(sorted(picked_all), dtype=int)
    mean_w_train = float(np.mean(raw_weights[train_idx])) if train_idx.size else 1.0
    if mean_w_train > 0.0:
        raw_weights = raw_weights / mean_w_train

    return {
        "subsampled": True,
        "method": "stratified",
        "strata_col": chosen_col,
        "train_indices": train_idx,
        "inclusion_weights": raw_weights,
        "exceeded_max_cells_for_rare_floor": bool(exceeded_floor),
        "n_strata": len(sorted_labels),
        "n_used": int(n_used),
    }


def _prepare_view_matrices(
    *,
    mod_name: str,
    mod_kind: str,
    matrix_all: Any,
    var_names: list[str],
    is_raw_counts: bool,
    train_indices: Any,
    n_top_genes: int,
    n_top_peaks: int,
    scale_views: bool,
) -> dict[str, Any]:
    """Fit feature selection & scaling on `train_indices` and transform all cells."""
    import anndata as ad
    import numpy as np
    import scanpy as sc
    import scipy.sparse as sp

    is_sparse = sp.issparse(matrix_all)
    x_csr = sp.csr_matrix(matrix_all) if is_sparse else None
    n_all, n_vars = matrix_all.shape
    n_train = int(len(train_indices))

    if is_sparse:
        cell_totals = np.asarray(x_csr.sum(axis=1)).ravel().astype(np.float64)
    else:
        arr_all = np.asarray(matrix_all, dtype=np.float64)
        cell_totals = np.sum(np.maximum(arr_all, 0.0), axis=1)

    if mod_kind == "rna" and is_raw_counts:
        # Fit HVG selection strictly on training cells using log1p(CP10k).
        x_train_raw = (
            x_csr[train_indices, :].copy().astype(np.float64)
            if is_sparse
            else sp.csr_matrix(np.asarray(matrix_all[train_indices, :], dtype=np.float64))
        )
        inv_tot_tr = np.zeros(n_train, dtype=np.float64)
        tot_tr = cell_totals[train_indices]
        nz_tr = tot_tr > 0.0
        inv_tot_tr[nz_tr] = 1e4 / tot_tr[nz_tr]
        x_train_norm = (sp.diags(inv_tot_tr) @ x_train_raw).tocsr()
        x_train_norm.data = np.log1p(x_train_norm.data)

        n_pick = min(max(2, int(n_top_genes)), n_vars)
        if n_pick < n_vars:
            tmp_ad = ad.AnnData(X=x_train_norm.copy())
            tmp_ad.var_names = [str(v) for v in var_names]
            try:
                sc.pp.highly_variable_genes(
                    tmp_ad,
                    n_top_genes=n_pick,
                    flavor="seurat",
                    subset=False,
                )
                hvg_mask = np.asarray(tmp_ad.var["highly_variable"].values, dtype=bool)
                sel_idx = np.flatnonzero(hvg_mask)
                if sel_idx.size == 0:
                    raise RuntimeError("No HVGs selected")
            except Exception:
                # Deterministic variance fallback for small synthetic matrices.
                mean_f = np.asarray(x_train_norm.mean(axis=0)).ravel()
                sq_f = x_train_norm.copy()
                sq_f.data = sq_f.data ** 2
                var_f = np.maximum(0.0, np.asarray(sq_f.mean(axis=0)).ravel() - mean_f ** 2)
                sel_idx = np.sort(np.argsort(var_f)[::-1][:n_pick])
        else:
            sel_idx = np.arange(n_vars, dtype=int)

        y_train = x_train_norm[:, sel_idx].toarray().astype(np.float64)

        if n_train == n_all and np.array_equal(train_indices, np.arange(n_all)):
            y_all = y_train
        else:
            x_sub_all = (
                x_csr[:, sel_idx].copy().astype(np.float64)
                if is_sparse
                else sp.csr_matrix(np.asarray(matrix_all[:, sel_idx], dtype=np.float64))
            )
            inv_tot_all = np.zeros(n_all, dtype=np.float64)
            nz_all = cell_totals > 0.0
            inv_tot_all[nz_all] = 1e4 / cell_totals[nz_all]
            x_norm_all = (sp.diags(inv_tot_all) @ x_sub_all).tocsr()
            x_norm_all.data = np.log1p(x_norm_all.data)
            y_all = x_norm_all.toarray().astype(np.float64)

    elif mod_kind == "atac" and is_raw_counts:
        # Signac log-TF-log-IDF fitted strictly on training cells over full cell totals.
        x_train_raw = (
            x_csr[train_indices, :].copy().astype(np.float64)
            if is_sparse
            else sp.csr_matrix(np.asarray(matrix_all[train_indices, :], dtype=np.float64))
        )
        # Filter ultra-sparse and near-ubiquitous peaks on training subset.
        peak_nnz = np.asarray((x_train_raw > 0).sum(axis=0)).ravel()
        min_cells_pk = max(3, int(math.ceil(0.005 * n_train)))
        max_cells_pk = max(min_cells_pk + 1, int(math.floor(0.99 * n_train)))
        cand_idx = np.flatnonzero((peak_nnz >= min_cells_pk) & (peak_nnz <= max_cells_pk))
        if cand_idx.size < min(20, n_vars):
            cand_idx = np.flatnonzero(peak_nnz >= 1)
        if cand_idx.size == 0:
            cand_idx = np.arange(min(n_vars, max(2, int(n_top_peaks))), dtype=int)

        peak_sums_tr = np.asarray(x_train_raw[:, cand_idx].sum(axis=0)).ravel()
        idf_cand = np.log1p(float(n_train) / np.maximum(peak_sums_tr, 1.0))

        tot_tr = cell_totals[train_indices]
        inv_tot_tr = np.zeros(n_train, dtype=np.float64)
        nz_tr = tot_tr > 0.0
        inv_tot_tr[nz_tr] = 1.0 / tot_tr[nz_tr]

        tfidf_cand = (
            sp.diags(inv_tot_tr) @ x_train_raw[:, cand_idx] @ sp.diags(idf_cand)
        ).tocsr()
        tfidf_cand.data = np.log1p(tfidf_cand.data * 1e4)

        n_pick = min(max(2, int(n_top_peaks)), int(cand_idx.size))
        if n_pick < cand_idx.size:
            mean_p = np.asarray(tfidf_cand.mean(axis=0)).ravel()
            sq_p = tfidf_cand.copy()
            sq_p.data = sq_p.data ** 2
            var_p = np.maximum(0.0, np.asarray(sq_p.mean(axis=0)).ravel() - mean_p ** 2)
            local_pick = np.sort(np.argsort(var_p)[::-1][:n_pick])
        else:
            local_pick = np.arange(cand_idx.size, dtype=int)

        sel_idx = cand_idx[local_pick]
        idf_sel = idf_cand[local_pick]
        y_train = tfidf_cand[:, local_pick].toarray().astype(np.float64)

        if n_train == n_all and np.array_equal(train_indices, np.arange(n_all)):
            y_all = y_train
        else:
            x_sub_all = (
                x_csr[:, sel_idx].copy().astype(np.float64)
                if is_sparse
                else sp.csr_matrix(np.asarray(matrix_all[:, sel_idx], dtype=np.float64))
            )
            inv_tot_all = np.zeros(n_all, dtype=np.float64)
            nz_all = cell_totals > 0.0
            inv_tot_all[nz_all] = 1.0 / cell_totals[nz_all]
            tfidf_all = (sp.diags(inv_tot_all) @ x_sub_all @ sp.diags(idf_sel)).tocsr()
            tfidf_all.data = np.log1p(tfidf_all.data * 1e4)
            y_all = tfidf_all.toarray().astype(np.float64)

    else:
        # Generic continuous view or pre-normalized matrix: select top variable features.
        x_tr = (
            x_csr[train_indices, :].toarray().astype(np.float64)
            if is_sparse
            else np.asarray(matrix_all[train_indices, :], dtype=np.float64)
        )
        n_pick = min(max(2, int(n_top_genes)), n_vars)
        if n_pick < n_vars:
            var_f = np.var(x_tr, axis=0)
            sel_idx = np.sort(np.argsort(var_f)[::-1][:n_pick])
        else:
            sel_idx = np.arange(n_vars, dtype=int)
        y_train = x_tr[:, sel_idx]
        if n_train == n_all and np.array_equal(train_indices, np.arange(n_all)):
            y_all = y_train
        else:
            y_all = (
                x_csr[:, sel_idx].toarray().astype(np.float64)
                if is_sparse
                else np.asarray(matrix_all[:, sel_idx], dtype=np.float64)
            )

    mu_train = np.mean(y_train, axis=0, keepdims=True)
    y_train_centered = y_train - mu_train
    sigma_train = float(np.std(y_train_centered)) if scale_views else 1.0
    if not np.isfinite(sigma_train) or sigma_train <= 1e-12:
        sigma_train = 1.0

    selected_features = [str(var_names[int(i)]) for i in sel_idx]
    return {
        "mod_name": mod_name,
        "mod_kind": mod_kind,
        "selected_indices": sel_idx,
        "selected_features": selected_features,
        "n_features_total": int(n_vars),
        "n_features_selected": int(len(selected_features)),
        "cell_totals": cell_totals,
        "y_train": y_train,
        "y_all": y_all,
        "mu_train": mu_train,
        "sigma_train": sigma_train,
    }


def _extract_or_compute_cell_cycle_scores(
    *,
    obs_dfs: list[Any],
    rna_raw_train: Any | None,
    rna_var_names: list[str] | None,
    cell_cycle_genes: dict[str, Any] | None,
    seed: int,
) -> dict[str, Any]:
    """Extract existing `.obs` S/G2M scores or score caller-supplied gene lists."""
    import numpy as np

    # 1. Prefer pre-computed S/G2M scores already present in `.obs`.
    for df in obs_dfs:
        if df is None:
            continue
        s_col = next((c for c in _S_SCORE_COLS if c in df.columns), None)
        g2m_col = next((c for c in _G2M_SCORE_COLS if c in df.columns), None)
        if s_col is not None and g2m_col is not None:
            try:
                s_vals = np.asarray(df[s_col].values, dtype=np.float64)
                g_vals = np.asarray(df[g2m_col].values, dtype=np.float64)
                if np.all(np.isfinite(s_vals)) and np.all(np.isfinite(g_vals)):
                    return {
                        "available": True,
                        "source": f"obs:{s_col},{g2m_col}",
                        "scope": "all_cells",
                        "s_scores": s_vals,
                        "g2m_scores": g_vals,
                    }
            except Exception:
                pass

    # 2. Score caller-provided `cell_cycle_genes={"s_genes": [...], "g2m_genes": [...]}`.
    if (
        isinstance(cell_cycle_genes, dict)
        and rna_raw_train is not None
        and rna_var_names is not None
    ):
        s_in = [str(g).strip() for g in (cell_cycle_genes.get("s_genes") or []) if str(g).strip()]
        g_in = [
            str(g).strip()
            for g in (cell_cycle_genes.get("g2m_genes") or [])
            if str(g).strip()
        ]
        if s_in and g_in:
            upper_map: dict[str, str] = {}
            for vn in rna_var_names:
                upper_map.setdefault(str(vn).upper(), str(vn))
            s_matched = [upper_map[g.upper()] for g in s_in if g.upper() in upper_map]
            g_matched = [upper_map[g.upper()] for g in g_in if g.upper() in upper_map]
            if len(s_matched) >= 3 and len(g_matched) >= 3:
                try:
                    import anndata as ad
                    import scanpy as sc
                    import scipy.sparse as sp

                    tmp = ad.AnnData(
                        X=(
                            sp.csr_matrix(rna_raw_train).copy().astype(np.float64)
                            if sp.issparse(rna_raw_train)
                            else np.asarray(rna_raw_train, dtype=np.float64).copy()
                        )
                    )
                    tmp.var_names = [str(v) for v in rna_var_names]
                    sc.pp.normalize_total(tmp, target_sum=1e4)
                    sc.pp.log1p(tmp)
                    sc.tl.score_genes_cell_cycle(
                        tmp,
                        s_genes=s_matched,
                        g2m_genes=g_matched,
                        random_state=int(seed),
                    )
                    return {
                        "available": True,
                        "source": "caller_cell_cycle_genes",
                        "scope": "train_cells",
                        "s_scores": np.asarray(tmp.obs["S_score"].values, dtype=np.float64),
                        "g2m_scores": np.asarray(tmp.obs["G2M_score"].values, dtype=np.float64),
                        "n_s_genes_matched": int(len(s_matched)),
                        "n_g2m_genes_matched": int(len(g_matched)),
                    }
                except Exception as exc:
                    return {
                        "available": False,
                        "reason": f"score_genes_cell_cycle_failed:{type(exc).__name__}",
                    }
            return {
                "available": False,
                "reason": (
                    f"insufficient_matched_cell_cycle_genes:"
                    f"s={len(s_matched)},g2m={len(g_matched)}"
                ),
            }

    return {
        "available": False,
        "reason": "no_obs_cell_cycle_scores_or_caller_gene_sets",
    }


def integration_mofa(params: dict) -> dict:
    import numpy as np

    h5mu_param = params.get("h5mu_path") or params.get("data_path")
    modalities = dict(params.get("modalities") or {})
    genome = str(params.get("genome") or "hg38")
    organism = str(params.get("organism") or "Homo sapiens")
    n_factors_req = max(1, int(params.get("n_factors", 10)))
    max_cells = int(params.get("max_cells", 5000))
    min_stratum_cells = max(1, int(params.get("min_stratum_cells", 5)))
    subsample_strata_col = params.get("subsample_strata_col")
    n_top_genes = max(10, int(params.get("n_top_genes", 3000)))
    n_top_peaks = max(10, int(params.get("n_top_peaks", 5000)))
    scale_views = bool(params.get("scale_views", True))
    max_iter = max(10, int(params.get("max_iter", 1000)))
    convergence_mode = str(params.get("convergence_mode") or "fast")
    drop_r2_param = params.get("drop_r2")
    drop_r2 = float(drop_r2_param) if drop_r2_param is not None else None
    active_r2_min = float(params.get("active_r2_min", 0.01))
    rel_r2_min = float(params.get("rel_r2_min", 0.05))
    min_floor_r2 = float(params.get("min_floor_r2", 0.0008))
    depth_corr_cutoff = float(params.get("depth_corr_cutoff", 0.70))
    cell_cycle_corr_cutoff = float(params.get("cell_cycle_corr_cutoff", 0.50))
    cell_cycle_genes = params.get("cell_cycle_genes")
    batch_key_param = params.get("batch_key")
    condition_col_param = params.get("condition_col")
    batch_eta2_cutoff = float(params.get("batch_eta2_cutoff", 0.65))
    exclude_batch_factors = bool(params.get("exclude_batch_factors", True))
    top_n_loadings = max(5, int(params.get("top_n_loadings", 25)))
    min_shared_cells = max(10, int(params.get("min_cells", 100)))
    seed = int(params.get("seed", 0))
    output_dir = str(params.get("output_dir") or "/tmp/aria_mofa")
    write_h5mu = bool(params.get("write_h5mu", bool(h5mu_param)))
    save_model_data = bool(params.get("save_model_data", False))
    resume = bool(params.get("resume", True))
    allow_mock = mocks_allowed(params)
    warnings: list[str] = []
    metrics_not_computed: list[str] = []

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Resolve paired .h5mu vs multi-modality file map.
    resolved_h5mu: str | None = None
    if h5mu_param is not None and str(h5mu_param).strip():
        cand = Path(str(h5mu_param))
        if not cand.exists():
            return {
                "status": "error",
                "error_type": "MissingPairedH5mu",
                "details": f"Paired .h5mu file does not exist: {h5mu_param}",
            }
        resolved_h5mu = str(cand)
    else:
        if len(modalities) < 2:
            return {
                "status": "error",
                "error_type": "InsufficientModalities",
                "details": (
                    f"MOFA+ requires at least 2 modalities (or paired h5mu_path). "
                    f"Got: {list(modalities.keys())}"
                ),
            }
        # If all modality entries point to the same .h5mu container, route via .h5mu.
        h5mu_candidates: list[Path] = []
        all_are_h5mu = True
        for files in modalities.values():
            valid = [Path(str(f)) for f in (files or []) if Path(str(f)).exists()]
            if valid and valid[0].suffix.lower() == ".h5mu":
                h5mu_candidates.append(valid[0].resolve())
            else:
                all_are_h5mu = False
        if all_are_h5mu and len(h5mu_candidates) >= 2:
            if len(set(h5mu_candidates)) > 1:
                return {
                    "status": "error",
                    "error_type": "MismatchedPairedH5mu",
                    "details": (
                        "Modalities point to different .h5mu containers: "
                        f"{[str(p) for p in h5mu_candidates]}"
                    ),
                }
            resolved_h5mu = str(h5mu_candidates[0])

    stage_inputs: list[tuple[str, str | None]] = []
    if resolved_h5mu is not None:
        stage_inputs.extend(_expand_manifest_inputs("paired_input", resolved_h5mu))
    else:
        for mod_name in sorted(modalities.keys()):
            valid_m = [str(f) for f in (modalities[mod_name] or []) if Path(str(f)).exists()]
            for idx, vf in enumerate(valid_m):
                stage_inputs.extend(_expand_manifest_inputs(f"mod:{mod_name}:{idx}", vf))

    manifest_params = {
        "genome": genome,
        "organism": organism,
        "n_factors": n_factors_req,
        "max_cells": max_cells,
        "min_stratum_cells": min_stratum_cells,
        "subsample_strata_col": (
            str(subsample_strata_col) if subsample_strata_col is not None else None
        ),
        "n_top_genes": n_top_genes,
        "n_top_peaks": n_top_peaks,
        "scale_views": scale_views,
        "max_iter": max_iter,
        "convergence_mode": convergence_mode,
        "drop_r2": drop_r2,
        "active_r2_min": active_r2_min,
        "rel_r2_min": rel_r2_min,
        "min_floor_r2": min_floor_r2,
        "depth_corr_cutoff": depth_corr_cutoff,
        "cell_cycle_corr_cutoff": cell_cycle_corr_cutoff,
        "cell_cycle_genes": cell_cycle_genes if isinstance(cell_cycle_genes, dict) else None,
        "batch_key": (
            str(batch_key_param).strip()
            if batch_key_param is not None and str(batch_key_param).strip()
            else None
        ),
        "condition_col": (
            str(condition_col_param).strip()
            if condition_col_param is not None and str(condition_col_param).strip()
            else None
        ),
        "batch_eta2_cutoff": batch_eta2_cutoff,
        "exclude_batch_factors": exclude_batch_factors,
        "top_n_loadings": top_n_loadings,
        "min_cells": min_shared_cells,
        "seed": seed,
        "write_h5mu": write_h5mu,
        "save_model_data": save_model_data,
    }
    out_dir_candidate = Path(output_dir)
    manifest_path = out_dir_candidate / "integration_mofa.stage_manifest.json"
    cached_result_path = out_dir_candidate / ".integration_mofa.result.json"
    req_outputs = [
        out_dir_candidate / "mofa_model.hdf5",
        out_dir_candidate / "mofa_factor_scores.csv",
        out_dir_candidate / "mofa_loadings.json",
    ]
    if write_h5mu and resolved_h5mu is not None:
        req_outputs.append(out_dir_candidate / "mofa_integrated.h5mu")

    if (
        resume
        and stage_inputs
        and cached_result_path.is_file()
        and all(p.is_file() and p.stat().st_size > 0 for p in req_outputs)
    ):
        is_curr, _why = stage_is_current(
            manifest_path,
            inputs=stage_inputs,
            params=manifest_params,
            tool_version=_STAGE_VERSION,
        )
        if is_curr:
            try:
                cached_payload = json.loads(
                    cached_result_path.read_text(encoding="utf-8")
                )
                if cached_payload.get("status") == "success":
                    cached_payload["resumed"] = True
                    cached_payload["warnings"] = list(
                        cached_payload.get("warnings") or []
                    ) + [
                        "[resume] Reused current MOFA+ stage output verified by SHA-256 stage manifest."
                    ]
                    return cached_payload
            except Exception:
                pass

    try:
        import anndata as ad
        import mofapy2  # noqa: F401
        from mofapy2.run.entry_point import entry_point
        import pandas as pd
        import scanpy as sc  # noqa: F401
        import scipy.sparse as sp

        from aria.utils.mudata_io import (
            _select_raw_counts_matrix,
            read_h5mu_paired,
        )

        raw_views: list[dict[str, Any]] = []
        obs_dfs_aligned: list[Any] = []
        src_mdata = None
        rna_mod_key_in_mdata: str | None = None
        atac_mod_key_in_mdata: str | None = None

        if resolved_h5mu is not None:
            paired = read_h5mu_paired(resolved_h5mu, require_raw_counts=True)
            if paired.get("status") != "success":
                return {
                    "status": "error",
                    "error_type": paired.get("error_type", "H5muReadFailed"),
                    "details": paired.get("details", "Could not read paired .h5mu."),
                    "warnings": warnings,
                }
            common_cells = list(paired["shared_barcodes"])
            if len(common_cells) < min_shared_cells:
                return {
                    "status": "error",
                    "error_type": "InsufficientSharedCells",
                    "details": (
                        f"Only {len(common_cells)} shared cells in paired .h5mu; "
                        f"MOFA+ requires at least {min_shared_cells}."
                    ),
                }

            src_mdata = paired["mdata"]
            src_rna = paired["rna_adata"]
            src_atac = paired["atac_adata"]
            rna_mod_key_in_mdata = str(paired["rna_modality"])
            atac_mod_key_in_mdata = str(paired["atac_modality"])

            r_map = {str(bc): i for i, bc in enumerate(src_rna.obs_names)}
            a_map = {str(bc): i for i, bc in enumerate(src_atac.obs_names)}
            r_take = [r_map[bc] for bc in common_cells]
            a_take = [a_map[bc] for bc in common_cells]

            if src_mdata.obs is not None and len(src_mdata.obs.columns) > 0:
                m_map = {str(bc): i for i, bc in enumerate(src_mdata.obs_names)}
                if all(bc in m_map for bc in common_cells):
                    m_take = [m_map[bc] for bc in common_cells]
                    obs_dfs_aligned.append(src_mdata.obs.iloc[m_take].copy())

            obs_dfs_aligned.append(src_rna.obs.iloc[r_take].copy())
            obs_dfs_aligned.append(src_atac.obs.iloc[a_take].copy())

            raw_views.append(
                {
                    "name": rna_mod_key_in_mdata,
                    "kind": "rna",
                    "matrix": sp.csr_matrix(paired["rna_counts_matrix"][r_take, :]),
                    "var_names": [str(v) for v in src_rna.var_names],
                    "is_raw_counts": True,
                    "counts_source": str(paired["rna_counts_source"]),
                }
            )
            raw_views.append(
                {
                    "name": atac_mod_key_in_mdata,
                    "kind": "atac",
                    "matrix": sp.csr_matrix(paired["atac_counts_matrix"][a_take, :]),
                    "var_names": [str(v) for v in src_atac.var_names],
                    "is_raw_counts": True,
                    "counts_source": str(paired["atac_counts_source"]),
                }
            )
        else:
            loaded_adatas: dict[str, Any] = {}
            for mod_name, files in modalities.items():
                valid_files = [str(f) for f in (files or []) if Path(str(f)).exists()]
                if not valid_files:
                    warnings.append(f"No valid files for {mod_name} — skipping.")
                    continue
                adata_loaded = _load_modality(mod_name, valid_files[0])
                if adata_loaded is not None:
                    loaded_adatas[str(mod_name)] = adata_loaded

            if len(loaded_adatas) < 2:
                if allow_mock:
                    return _mock_mofa(n_factors_req, reason="Could not load 2+ modalities")
                return {
                    "status": "error",
                    "error_type": "InsufficientModalities",
                    "details": "Could not load at least two valid modalities for MOFA+.",
                    "warnings": warnings,
                }

            # Canonical deterministic barcode intersection preserving first modality order.
            mod_keys = list(loaded_adatas.keys())
            first_obs = [str(bc) for bc in loaded_adatas[mod_keys[0]].obs_names]
            shared_set = set(first_obs)
            for mk in mod_keys[1:]:
                shared_set &= {str(bc) for bc in loaded_adatas[mk].obs_names}
            seen_bc: set[str] = set()
            common_cells = []
            for bc in first_obs:
                if bc in shared_set and bc not in seen_bc:
                    common_cells.append(bc)
                    seen_bc.add(bc)

            if len(common_cells) < min_shared_cells:
                return {
                    "status": "error",
                    "error_type": "InsufficientSharedCells",
                    "details": (
                        f"Only {len(common_cells)} cells shared across modalities. "
                        f"MOFA+ requires at least {min_shared_cells} shared cells."
                    ),
                }

            for mk in mod_keys:
                ad_obj = loaded_adatas[mk]
                bc_map = {str(bc): i for i, bc in enumerate(ad_obj.obs_names)}
                take = [bc_map[bc] for bc in common_cells]
                ad_sub = ad_obj[take].copy()
                obs_dfs_aligned.append(ad_sub.obs.copy())

                kind = _classify_modality_kind(mk)
                raw_chk = _select_raw_counts_matrix(
                    ad_sub, modality_label=mk, require_raw_counts=False
                )
                is_raw = bool(raw_chk.get("is_raw_counts", False))
                mat_sel = raw_chk["matrix"] if raw_chk.get("matrix") is not None else ad_sub.X
                if not is_raw and kind in ("rna", "atac"):
                    warnings.append(
                        f"{mk}: non-integer input matrix detected; using values as "
                        "provided without re-applying count normalization."
                    )
                raw_views.append(
                    {
                        "name": mk,
                        "kind": kind,
                        "matrix": mat_sel,
                        "var_names": [str(v) for v in ad_sub.var_names],
                        "is_raw_counts": is_raw,
                        "counts_source": str(raw_chk.get("source") or "X"),
                    }
                )

        # ── Step 1: Deterministic subsampling (preserving rare strata) ────
        sub_info = _deterministic_subsample(
            common_cells,
            obs_dfs_aligned,
            max_cells=max_cells,
            min_stratum_cells=min_stratum_cells,
            explicit_strata_col=(
                str(subsample_strata_col) if subsample_strata_col else None
            ),
            seed=seed,
        )
        train_indices = sub_info["train_indices"]
        inclusion_weights_all = sub_info["inclusion_weights"]
        train_weights = inclusion_weights_all[train_indices]
        n_cells_total = int(len(common_cells))
        n_cells_trained = int(len(train_indices))
        train_barcodes = [common_cells[int(i)] for i in train_indices]

        if sub_info["subsampled"]:
            warnings.append(
                f"Subsampled {n_cells_trained}/{n_cells_total} cells deterministically "
                f"(method={sub_info['method']}, strata={sub_info['strata_col']}, seed={seed}) "
                "for MOFA+ training; remaining cells projected via exact posterior mean."
            )

        # ── Step 2: Fit view preprocessing strictly on training subset ────
        prepared_views: list[dict[str, Any]] = []
        for rv in raw_views:
            prep = _prepare_view_matrices(
                mod_name=rv["name"],
                mod_kind=rv["kind"],
                matrix_all=rv["matrix"],
                var_names=rv["var_names"],
                is_raw_counts=rv["is_raw_counts"],
                train_indices=train_indices,
                n_top_genes=n_top_genes,
                n_top_peaks=n_top_peaks,
                scale_views=scale_views,
            )
            prepared_views.append(prep)

        min_dim = min(
            [n_cells_trained - 1]
            + [max(1, pv["n_features_selected"] - 1) for pv in prepared_views]
        )
        n_factors_fit = max(1, min(n_factors_req, int(min_dim)))
        if n_factors_fit < n_factors_req:
            warnings.append(
                f"Clamped n_factors from {n_factors_req} to {n_factors_fit} "
                f"based on training dimensions (n_cells={n_cells_trained})."
            )

        # ── Step 3: Train variational Bayes MOFA+ model ───────────────────
        view_names = [pv["mod_name"] for pv in prepared_views]
        with contextlib.redirect_stdout(io.StringIO()):
            ent = entry_point()
            ent.set_data_options(
                scale_groups=False,
                scale_views=bool(scale_views),
                center_groups=True,
            )
            ent.set_data_matrix(
                data=[[pv["y_train"].copy()] for pv in prepared_views],
                views_names=view_names,
                groups_names=["group1"],
                samples_names=[train_barcodes],
                features_names=[pv["selected_features"] for pv in prepared_views],
            )
            ent.set_model_options(
                factors=int(n_factors_fit),
                spikeslab_weights=True,
                spikeslab_factors=False,
                ard_factors=False,
                ard_weights=True,
            )
            ent.set_train_options(
                iter=int(max_iter),
                convergence_mode=convergence_mode,
                startELBO=1,
                freqELBO=5,
                dropR2=drop_r2,
                startDrop=50 if drop_r2 is not None else 1,
                verbose=False,
                quiet=True,
                seed=int(seed),
            )
            ent.build()
            ent.run()

        model = ent.model
        z_train_raw = np.asarray(model.nodes["Z"].getExpectations()["E"], dtype=np.float64)
        w_expectations = model.nodes["W"].getExpectations()
        tau_expectations = model.nodes["Tau"].getExpectations()
        w_raw_list = [
            np.asarray(w_expectations[v]["E"], dtype=np.float64)
            for v in range(len(view_names))
        ]

        r2_raw_matrix = np.asarray(
            model.calculate_variance_explained()[0], dtype=np.float64
        )
        if r2_raw_matrix.ndim == 1:
            r2_raw_matrix = r2_raw_matrix.reshape(len(view_names), -1)
        r2_tot_vec = np.asarray(
            model.calculate_variance_explained(total=True)[0], dtype=np.float64
        ).ravel()

        k_actual = int(z_train_raw.shape[1])

        # ── Step 4: Exact posterior mean projection onto all cells ────────
        tau_vecs: list[np.ndarray] = []
        for v in range(len(view_names)):
            tau_v = np.asarray(tau_expectations[v]["E"], dtype=np.float64)
            tau_vec = tau_v[0] if tau_v.ndim == 2 else tau_v.ravel()
            tau_vecs.append(tau_vec)

        prec = np.eye(k_actual, dtype=np.float64)
        rhs_all = np.zeros((n_cells_total, k_actual), dtype=np.float64)
        for v, pv in enumerate(prepared_views):
            yv_scaled = (pv["y_all"] - pv["mu_train"]) / pv["sigma_train"]
            wv = w_raw_list[v]
            tv = tau_vecs[v]
            prec += (wv.T * tv) @ wv
            rhs_all += (yv_scaled * tv) @ wv

        z_proj_all = np.linalg.solve(prec, rhs_all.T).T
        z_proj_train = z_proj_all[train_indices, :]

        proj_check_rs: list[float] = []
        for k in range(k_actual):
            if float(np.std(z_train_raw[:, k])) > 1e-10:
                pr = pearson_correlation(z_train_raw[:, k], z_proj_train[:, k])
                if pr is not None:
                    proj_check_rs.append(abs(float(pr)))
        min_proj_r = float(min(proj_check_rs)) if proj_check_rs else 1.0

        if sub_info["subsampled"]:
            z_all_raw = z_proj_all.copy()
            # Preserve exact variational posterior expectations on training subset.
            z_all_raw[train_indices, :] = z_train_raw
            projection_method = "exact_variational_posterior_mean"
            if min_proj_r < 0.99:
                warnings.append(
                    f"Out-of-subset posterior projection min |r|={min_proj_r:.4f} (< 0.99)."
                )
        else:
            z_all_raw = z_train_raw.copy()
            projection_method = "none_all_cells_trained"

        # ── Step 5: Canonical factor ordering & deterministic sign ────────
        r2_pos_matrix = np.maximum(0.0, r2_raw_matrix)
        factor_r2_sums = np.sum(r2_pos_matrix, axis=0)
        order = np.argsort(-factor_r2_sums, kind="mergesort")

        z_train = z_train_raw[:, order].copy()
        z_all = z_all_raw[:, order].copy()
        w_list = [w[:, order].copy() for w in w_raw_list]
        r2_raw = r2_raw_matrix[:, order].copy()
        r2_pos = r2_pos_matrix[:, order].copy()

        rna_view_idx = next(
            (i for i, pv in enumerate(prepared_views) if pv["mod_kind"] == "rna"),
            0,
        )
        for k in range(k_actual):
            # Anchor sign on highest-R^2 view for factor k (tie-breaking to RNA).
            best_v = int(np.argmax(r2_pos[:, k]))
            if r2_pos[best_v, k] <= 0.0:
                best_v = rna_view_idx
            w_anchor = w_list[best_v][:, k]
            if w_anchor.size > 0:
                top_f = int(np.argmax(np.abs(w_anchor)))
                if float(w_anchor[top_f]) < 0.0:
                    z_train[:, k] *= -1.0
                    z_all[:, k] *= -1.0
                    for v in range(len(view_names)):
                        w_list[v][:, k] *= -1.0

        # ── Step 6: Variance explained & factor classification ────────────
        variance_explained: dict[str, float] = {}
        variance_explained_additive: dict[str, float] = {}
        view_pos_sums: dict[str, float] = {}
        for v, v_name in enumerate(view_names):
            tot_v = max(0.0, float(r2_tot_vec[v]))
            add_v = float(np.sum(r2_pos[v, :]))
            variance_explained[v_name] = round(tot_v, 4)
            variance_explained_additive[v_name] = round(add_v, 4)
            view_pos_sums[v_name] = max(tot_v, add_v)
            if tot_v < 0.05:
                warnings.append(
                    f"{v_name} explains {tot_v:.1%} total variance across MOFA+ factors "
                    "(note: sparse chromatin views naturally exhibit lower global R^2 "
                    "than transcriptomic HVG views)."
                )

        # ── Step 7: Quantitative confounder screening (depth & cell cycle) ─
        depth_log10_train: dict[str, np.ndarray] = {
            pv["mod_name"]: np.log10(1.0 + pv["cell_totals"][train_indices])
            for pv in prepared_views
        }

        rna_rv = next((rv for rv in raw_views if rv["kind"] == "rna"), None)
        rna_raw_train = (
            rna_rv["matrix"][train_indices, :]
            if (rna_rv is not None and rna_rv["is_raw_counts"])
            else None
        )
        obs_dfs_train = [df.iloc[train_indices].copy() for df in obs_dfs_aligned if df is not None]
        cc_info = _extract_or_compute_cell_cycle_scores(
            obs_dfs=obs_dfs_train,
            rna_raw_train=rna_raw_train,
            rna_var_names=rna_rv["var_names"] if rna_rv is not None else None,
            cell_cycle_genes=cell_cycle_genes if isinstance(cell_cycle_genes, dict) else None,
            seed=seed,
        )
        if not cc_info.get("available", False):
            metrics_not_computed.append("cell_cycle_factor_correlation")

        explicit_batch_col = (
            str(batch_key_param).strip()
            if batch_key_param is not None and str(batch_key_param).strip()
            else None
        )
        explicit_cond_col = (
            str(condition_col_param).strip()
            if condition_col_param is not None and str(condition_col_param).strip()
            else None
        )
        resolved_batch_col, batch_labels_train, batch_source = _resolve_categorical_obs(
            obs_dfs_train,
            explicit_batch_col,
            _AUTO_BATCH_CANDIDATES,
        )
        if explicit_batch_col and batch_labels_train is None:
            warnings.append(
                f"Requested batch_key='{explicit_batch_col}' was missing or had <2 levels "
                "on MOFA+ training cells; batch factor screening skipped."
            )
        if batch_labels_train is None:
            metrics_not_computed.append("batch_factor_association")

        resolved_cond_col, cond_labels_train, _cond_source = _resolve_categorical_obs(
            obs_dfs_train,
            explicit_cond_col,
            _AUTO_CONDITION_CANDIDATES,
            exclude_col=resolved_batch_col,
        )
        batch_condition_confounded = False
        if batch_labels_train is not None and cond_labels_train is not None:
            cond_to_batches: dict[str, set[str]] = {}
            for b_lbl, c_lbl in zip(batch_labels_train, cond_labels_train):
                cond_to_batches.setdefault(c_lbl, set()).add(b_lbl)
            if cond_to_batches and max(len(bs) for bs in cond_to_batches.values()) <= 1:
                batch_condition_confounded = True
                warnings.append(
                    f"Batch/donor column '{resolved_batch_col}' is 1:1 confounded with "
                    f"condition column '{resolved_cond_col}' (1 donor/batch per condition); "
                    "donor-specific technical effects cannot be disentangled from condition."
                )

        top_factors: list[dict[str, Any]] = []
        shared_factors: list[dict[str, Any]] = []
        modality_specific_factors: dict[str, list[dict[str, Any]]] = {
            v_name: [] for v_name in view_names
        }
        excluded_factors: list[dict[str, Any]] = []
        loadings_payload: dict[str, Any] = {}

        best_tech_fid: int | None = None
        best_tech_abs_r = -1.0
        cell_cycle_factor_ids: list[int] = []
        batch_factor_ids: list[int] = []
        n_active_factors = 0

        # Precompute depth and batch associations per factor so sequencing-depth or
        # extreme batch-outlier artifacts do not inflate the biological view
        # denominator `view_bio_sums` and mask genuine cross-modal factors.
        depth_info_per_k: list[tuple[dict[str, dict[str, float | None]], float, str | None, bool]] = []
        batch_info_per_k: list[tuple[float | None, float | None, float | None, bool]] = []
        for k in range(k_actual):
            depth_corrs_k: dict[str, dict[str, float | None]] = {}
            max_depth_r = 0.0
            max_depth_view: str | None = None
            for v_name, d_vec in depth_log10_train.items():
                r_unw = pearson_correlation(z_train[:, k], d_vec)
                r_w = pearson_correlation(z_train[:, k], d_vec, weights=train_weights)
                rho = spearman_correlation(z_train[:, k], d_vec)
                depth_corrs_k[v_name] = {
                    "pearson_r": round(float(r_unw), 4) if r_unw is not None else None,
                    "weighted_pearson_r": (
                        round(float(r_w), 4) if r_w is not None else None
                    ),
                    "spearman_rho": round(float(rho), 4) if rho is not None else None,
                }
                eff_r = max(
                    abs(float(r_w)) if r_w is not None else 0.0,
                    abs(float(r_unw)) if r_unw is not None else 0.0,
                )
                if eff_r > max_depth_r:
                    max_depth_r = eff_r
                    max_depth_view = v_name
            tentative_depth = bool(
                float(np.std(z_train[:, k])) > 1e-10 and max_depth_r >= depth_corr_cutoff
            )
            depth_info_per_k.append(
                (depth_corrs_k, max_depth_r, max_depth_view, tentative_depth)
            )

            b_eta2, b_fstat = (
                _compute_categorical_eta2(
                    z_train[:, k], batch_labels_train, weights=train_weights
                )
                if batch_labels_train is not None
                else (None, None)
            )
            c_eta2, _ = (
                _compute_categorical_eta2(
                    z_train[:, k], cond_labels_train, weights=train_weights
                )
                if cond_labels_train is not None
                else (None, None)
            )
            tentative_batch = bool(
                b_eta2 is not None
                and b_eta2 >= max(0.80, batch_eta2_cutoff)
                and (
                    c_eta2 is None
                    or batch_condition_confounded
                    or (b_eta2 - c_eta2) >= 0.25
                )
            )
            batch_info_per_k.append((b_eta2, b_fstat, c_eta2, tentative_batch))

        non_depth_mask = np.asarray(
            [
                not (depth_info_per_k[idx][3] or batch_info_per_k[idx][3])
                for idx in range(k_actual)
            ],
            dtype=bool,
        )
        view_bio_sums: dict[str, float] = {}
        for v, v_name in enumerate(view_names):
            non_depth_sum = (
                float(np.sum(r2_pos[v, non_depth_mask]))
                if np.any(non_depth_mask)
                else 0.0
            )
            view_bio_sums[v_name] = (
                non_depth_sum if non_depth_sum > 0.0 else view_pos_sums[v_name]
            )

        for k in range(k_actual):
            fid = k + 1
            depth_corrs_k, max_depth_r, max_depth_view, tentative_depth = depth_info_per_k[k]
            b_eta2, b_fstat, c_eta2, tentative_batch = batch_info_per_k[k]
            active_views_k: list[str] = []
            r2_by_view_k: dict[str, float] = {}
            r2_raw_by_view_k: dict[str, float] = {}

            for v, v_name in enumerate(view_names):
                val_pos = float(r2_pos[v, k])
                val_raw = float(r2_raw[v, k])
                r2_by_view_k[v_name] = round(val_pos, 4)
                r2_raw_by_view_k[v_name] = round(val_raw, 6)
                denom_v = (
                    view_pos_sums[v_name]
                    if (tentative_depth or tentative_batch)
                    else view_bio_sums[v_name]
                )
                if is_factor_active_in_view(
                    val_pos,
                    denom_v,
                    active_r2_min=active_r2_min,
                    rel_r2_min=rel_r2_min,
                    min_floor_r2=min_floor_r2,
                ):
                    active_views_k.append(v_name)

            is_active = len(active_views_k) > 0 and float(np.std(z_train[:, k])) > 1e-10
            if is_active:
                n_active_factors += 1

            if len(active_views_k) >= 2:
                category = "shared"
                act_vals = [r2_by_view_k[vn] for vn in active_views_k]
                mx = max(act_vals)
                sharing_balance = round(min(act_vals) / mx, 4) if mx > 0 else 0.0
                dominant_mod = view_names[int(np.argmax(r2_pos[:, k]))]
            elif len(active_views_k) == 1:
                category = "modality_specific"
                sharing_balance = 0.0
                dominant_mod = active_views_k[0]
            else:
                category = "low_variance"
                sharing_balance = 0.0
                dominant_mod = None

            is_depth_flagged = bool(is_active and max_depth_r >= depth_corr_cutoff)
            if is_depth_flagged and max_depth_r > best_tech_abs_r:
                best_tech_abs_r = max_depth_r
                best_tech_fid = fid

            # Cell-cycle correlation on training subset when scores are available.
            cc_corrs_k: dict[str, float | None] | None = None
            max_cc_r = 0.0
            if cc_info.get("available", False):
                s_sc = cc_info["s_scores"]
                g_sc = cc_info["g2m_scores"]
                r_s = pearson_correlation(z_train[:, k], s_sc, weights=train_weights)
                r_g = pearson_correlation(z_train[:, k], g_sc, weights=train_weights)
                rho_s = spearman_correlation(z_train[:, k], s_sc)
                rho_g = spearman_correlation(z_train[:, k], g_sc)
                cc_corrs_k = {
                    "s_score_weighted_r": round(float(r_s), 4) if r_s is not None else None,
                    "g2m_score_weighted_r": (
                        round(float(r_g), 4) if r_g is not None else None
                    ),
                    "s_score_spearman_rho": (
                        round(float(rho_s), 4) if rho_s is not None else None
                    ),
                    "g2m_score_spearman_rho": (
                        round(float(rho_g), 4) if rho_g is not None else None
                    ),
                }
                max_cc_r = max(
                    abs(float(r_s)) if r_s is not None else 0.0,
                    abs(float(r_g)) if r_g is not None else 0.0,
                )

            is_cc_flagged = bool(
                is_active
                and cc_info.get("available", False)
                and max_cc_r >= cell_cycle_corr_cutoff
            )
            if is_cc_flagged:
                cell_cycle_factor_ids.append(fid)

            is_batch_flagged = bool(
                is_active
                and batch_labels_train is not None
                and b_eta2 is not None
                and b_eta2 >= batch_eta2_cutoff
                and (
                    c_eta2 is None
                    or batch_condition_confounded
                    or (b_eta2 - c_eta2) >= 0.25
                )
            )
            if is_batch_flagged:
                batch_factor_ids.append(fid)

            batch_assoc_k: dict[str, Any] | None = None
            if batch_labels_train is not None:
                batch_assoc_k = {
                    "batch_key": resolved_batch_col,
                    "eta2": b_eta2,
                    "f_stat": b_fstat,
                    "condition_col": resolved_cond_col,
                    "condition_eta2": c_eta2,
                    "batch_condition_confounded": bool(batch_condition_confounded),
                }

            exclusion_reasons: list[str] = []
            if is_depth_flagged:
                exclusion_reasons.append(
                    f"sequencing_depth_correlated:{max_depth_view}(|r|={max_depth_r:.3f})"
                )
            if is_cc_flagged:
                exclusion_reasons.append(
                    f"cell_cycle_correlated(|r|={max_cc_r:.3f})"
                )
            if is_batch_flagged and exclude_batch_factors:
                exclusion_reasons.append(
                    f"batch_or_donor_confounded:{resolved_batch_col}(eta2={b_eta2:.3f})"
                )
            excluded = len(exclusion_reasons) > 0

            factor_entry: dict[str, Any] = {
                "factor_id": fid,
                "active": bool(is_active),
                "category": category,
                "active_views": active_views_k,
                "dominant_modality": dominant_mod,
                "sharing_balance": sharing_balance,
                "r2_by_view": r2_by_view_k,
                "r2_raw_by_view": r2_raw_by_view_k,
                "depth_correlations": depth_corrs_k,
                "max_depth_abs_r": round(max_depth_r, 4),
                "technical_depth_flag": is_depth_flagged,
                "cell_cycle_correlations": cc_corrs_k,
                "cell_cycle_flag": is_cc_flagged,
                "batch_association": batch_assoc_k,
                "batch_eta2": b_eta2,
                "batch_flag": is_batch_flagged,
                "excluded_from_biological_interpretation": excluded,
                "exclusion_reasons": exclusion_reasons,
            }
            for v_name, r2_val in r2_by_view_k.items():
                factor_entry[f"r2_{v_name}"] = r2_val

            top_factors.append(factor_entry)
            if category == "shared":
                shared_factors.append(factor_entry)
            elif category == "modality_specific" and dominant_mod in modality_specific_factors:
                modality_specific_factors[dominant_mod].append(factor_entry)
            if excluded:
                excluded_factors.append(factor_entry)

            # Signed top loadings per view for factor k.
            factor_loadings_views: dict[str, list[dict[str, Any]]] = {}
            for v, pv in enumerate(prepared_views):
                wv_col = w_list[v][:, k]
                feats_v = pv["selected_features"]
                top_idx = np.argsort(np.abs(wv_col))[::-1][: min(top_n_loadings, len(feats_v))]
                factor_loadings_views[pv["mod_name"]] = [
                    {
                        "feature": feats_v[int(idx)],
                        "weight": round(float(wv_col[int(idx)]), 6),
                        "abs_weight": round(abs(float(wv_col[int(idx)])), 6),
                        "direction": "positive" if float(wv_col[int(idx)]) >= 0.0 else "negative",
                    }
                    for idx in top_idx
                ]
            loadings_payload[f"Factor{fid}"] = factor_loadings_views

        technical_factor_flag = best_tech_fid is not None
        technical_factor_id = best_tech_fid
        cell_cycle_factor_flag = len(cell_cycle_factor_ids) > 0
        batch_factor_flag = len(batch_factor_ids) > 0

        if technical_factor_flag:
            warnings.append(
                f"Factor {technical_factor_id} strongly correlates with library sequencing "
                f"depth (|r|={best_tech_abs_r:.3f} >= {depth_corr_cutoff:.2f}); marked "
                "excluded_from_biological_interpretation."
            )
        if cell_cycle_factor_flag:
            warnings.append(
                f"Factor(s) {cell_cycle_factor_ids} correlate with S/G2M phase scores "
                f"(|r| >= {cell_cycle_corr_cutoff:.2f}); marked "
                "excluded_from_biological_interpretation."
            )
        if batch_factor_flag:
            warnings.append(
                f"Factor(s) {batch_factor_ids} strongly associate with batch/donor "
                f"'{resolved_batch_col}' (ANOVA eta^2 >= {batch_eta2_cutoff:.2f})"
                + (
                    " and were marked excluded_from_biological_interpretation."
                    if exclude_batch_factors
                    else "."
                )
            )

        rna_or_first_name = prepared_views[rna_view_idx]["mod_name"]
        factor1_top = [
            item["feature"]
            for item in loadings_payload.get("Factor1", {}).get(rna_or_first_name, [])[:10]
        ]

        # ── Step 8: Atomic output persistence ─────────────────────────────
        out_dir_path = Path(output_dir)
        model_path = out_dir_path / "mofa_model.hdf5"
        scores_path = out_dir_path / "mofa_factor_scores.csv"
        loadings_path = out_dir_path / "mofa_loadings.json"

        tmp_model = out_dir_path / ".mofa_model.hdf5.tmp"
        if tmp_model.exists():
            tmp_model.unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            ent.save(
                str(tmp_model),
                save_data=bool(save_model_data),
                save_parameters=False,
            )
        os.replace(tmp_model, model_path)

        tmp_scores = out_dir_path / ".mofa_factor_scores.csv.tmp"
        pd.DataFrame(
            z_all,
            index=common_cells,
            columns=[f"Factor{i + 1}" for i in range(k_actual)],
        ).to_csv(tmp_scores)
        os.replace(tmp_scores, scores_path)

        tmp_loadings = out_dir_path / ".mofa_loadings.json.tmp"
        tmp_loadings.write_text(
            json.dumps(
                {
                    "seed": int(seed),
                    "n_factors": int(k_actual),
                    "views": view_names,
                    "canonical_factor_order_from_raw": [int(i) for i in order.tolist()],
                    "factors": loadings_payload,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        os.replace(tmp_loadings, loadings_path)

        h5mu_output_path: str | None = None
        if write_h5mu and src_mdata is not None:
            import mudata as md

            h5mu_dest = out_dir_path / "mofa_integrated.h5mu"
            tmp_h5mu = out_dir_path / ".mofa_integrated.h5mu.tmp"
            if tmp_h5mu.exists():
                tmp_h5mu.unlink()

            # Build clean aligned MuData preserving modalities and writing X_mofa + LFs.
            mod_Copies: dict[str, Any] = {}
            for v, (rv, pv) in enumerate(zip(raw_views, prepared_views)):
                m_key = rv["name"]
                orig_ad = src_mdata.mod[m_key]
                bc_map = {str(bc): i for i, bc in enumerate(orig_ad.obs_names)}
                take = [bc_map[bc] for bc in common_cells]
                ad_copy = orig_ad[take].copy()
                lfs_full = np.zeros((ad_copy.n_vars, k_actual), dtype=np.float32)
                lfs_full[pv["selected_indices"], :] = w_list[v].astype(np.float32)
                ad_copy.varm["LFs"] = lfs_full
                var_used = np.zeros(ad_copy.n_vars, dtype=bool)
                var_used[pv["selected_indices"]] = True
                ad_copy.var["mofa_feature_used"] = var_used
                mod_Copies[m_key] = ad_copy

            out_mdata = md.MuData(mod_Copies)
            if src_mdata.obs is not None and len(src_mdata.obs.columns) > 0:
                m_map = {str(bc): i for i, bc in enumerate(src_mdata.obs_names)}
                if all(bc in m_map for bc in common_cells):
                    m_take = [m_map[bc] for bc in common_cells]
                    for col in src_mdata.obs.columns:
                        out_mdata.obs[col] = src_mdata.obs[col].iloc[m_take].values
            for obsm_k, obsm_v in src_mdata.obsm.items():
                if getattr(obsm_v, "shape", (0,))[0] == len(common_cells):
                    out_mdata.obsm[obsm_k] = obsm_v.copy()

            train_mask = np.zeros(n_cells_total, dtype=bool)
            train_mask[train_indices] = True
            out_mdata.obs["mofa_training_subset"] = train_mask
            out_mdata.obs["mofa_inclusion_weight"] = inclusion_weights_all.astype(np.float32)
            out_mdata.obsm["X_mofa"] = z_all.astype(np.float32)
            out_mdata.uns["aria_mofa"] = sanitize_for_hdf5_uns(
                {
                    "seed": int(seed),
                    "n_factors": int(k_actual),
                    "n_active_factors": int(n_active_factors),
                    "scale_views": bool(scale_views),
                    "subsampled": bool(sub_info["subsampled"]),
                    "subsample_method": sub_info["method"],
                    "subsample_strata_col": sub_info["strata_col"],
                    "projection_method": projection_method,
                    "projection_min_train_abs_r": round(min_proj_r, 6),
                    "variance_explained": variance_explained,
                    "variance_explained_additive": variance_explained_additive,
                    "shared_factor_ids": [int(f["factor_id"]) for f in shared_factors],
                    "excluded_factor_ids": [int(f["factor_id"]) for f in excluded_factors],
                    "technical_factor_flag": bool(technical_factor_flag),
                    "technical_factor_id": technical_factor_id,
                    "cell_cycle_factor_flag": bool(cell_cycle_factor_flag),
                    "cell_cycle_factor_ids": cell_cycle_factor_ids,
                    "batch_evaluated": bool(batch_labels_train is not None),
                    "batch_key": resolved_batch_col if batch_labels_train is not None else None,
                    "batch_condition_confounded": bool(batch_condition_confounded),
                    "batch_factor_flag": bool(batch_factor_flag),
                    "batch_factor_ids": batch_factor_ids,
                }
            )
            out_mdata.write_h5mu(str(tmp_h5mu))
            os.replace(tmp_h5mu, h5mu_dest)
            h5mu_output_path = str(h5mu_dest)

        subsampling_summary = {
            "subsampled": bool(sub_info["subsampled"]),
            "method": str(sub_info["method"]),
            "strata_col": sub_info["strata_col"],
            "n_cells_total": int(n_cells_total),
            "n_cells_trained": int(n_cells_trained),
            "max_cells": int(max_cells),
            "exceeded_max_cells_for_rare_floor": bool(
                sub_info.get("exceeded_max_cells_for_rare_floor", False)
            ),
            "projection_method": projection_method,
            "projection_min_train_abs_r": round(min_proj_r, 6),
        }

        result_payload = {
            "status": "success",
            "genome": genome,
            "organism": organism,
            "seed": int(seed),
            "n_factors_requested": int(n_factors_req),
            "n_factors": int(k_actual),
            "n_active_factors": int(n_active_factors),
            "n_cells": int(n_cells_total),
            "n_cells_trained": int(n_cells_trained),
            "subsampled": bool(sub_info["subsampled"]),
            "subsampling": subsampling_summary,
            "scale_views": bool(scale_views),
            "views": {
                pv["mod_name"]: {
                    "kind": pv["mod_kind"],
                    "n_features_total": pv["n_features_total"],
                    "n_features_selected": pv["n_features_selected"],
                }
                for pv in prepared_views
            },
            "top_factors": top_factors,
            "shared_factors": shared_factors,
            "modality_specific_factors": modality_specific_factors,
            "excluded_factors": excluded_factors,
            "variance_explained": variance_explained,
            "variance_explained_additive": variance_explained_additive,
            "factor1_top_features": factor1_top,
            "technical_factor_flag": bool(technical_factor_flag),
            "technical_factor_id": technical_factor_id,
            "cell_cycle_evaluated": bool(cc_info.get("available", False)),
            "cell_cycle_source": cc_info.get("source"),
            "cell_cycle_factor_flag": bool(cell_cycle_factor_flag),
            "cell_cycle_factor_ids": cell_cycle_factor_ids,
            "batch_evaluated": bool(batch_labels_train is not None),
            "batch_key": resolved_batch_col if batch_labels_train is not None else None,
            "batch_source": batch_source,
            "condition_col": resolved_cond_col if cond_labels_train is not None else None,
            "batch_condition_confounded": bool(batch_condition_confounded),
            "batch_factor_flag": bool(batch_factor_flag),
            "batch_factor_ids": batch_factor_ids,
            "factor_scores": str(scores_path),
            "factor_loadings_path": str(loadings_path),
            "output_path": str(model_path),
            "h5mu_output_path": h5mu_output_path,
            "metrics_not_computed": metrics_not_computed,
            "resumed": False,
            "warnings": warnings,
        }

        if stage_inputs:
            try:
                tmp_res = cached_result_path.with_name(
                    f".integration_mofa.result.tmp.{os.getpid()}.json"
                )
                tmp_res.write_text(
                    json.dumps(result_payload, indent=2), encoding="utf-8"
                )
                os.replace(tmp_res, cached_result_path)
                write_stage_manifest(
                    manifest_path,
                    build_stage_manifest(
                        stage="integration_mofa",
                        inputs=stage_inputs,
                        params=manifest_params,
                        tool_version=_STAGE_VERSION,
                        extra={"output_path": str(model_path)},
                    ),
                )
            except Exception as manifest_exc:
                warnings.append(
                    f"Could not write MOFA+ stage manifest: {manifest_exc}"
                )

        return result_payload

    except ImportError as e:
        if allow_mock:
            return _mock_mofa(n_factors_req, reason=str(e))
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"MOFA+ dependencies are required: {e}",
            "warnings": warnings,
        }
    except Exception as e:
        return {
            "status": "error",
            "error_type": "MOFAFailed",
            "details": str(e)[:500],
            "warnings": warnings,
        }


def _load_modality(mod_name: str, path: str):
    """Load a single-modality matrix from `.h5ad`, `.h5`, or 10x `.mtx` dir."""
    try:
        import scanpy as sc

        p = Path(path)
        if p.suffix == ".h5ad":
            return sc.read_h5ad(str(p))
        if p.suffix == ".h5":
            return sc.read_10x_h5(str(p))
        if p.is_dir():
            return sc.read_10x_mtx(str(p), var_names="gene_symbols", cache=True)
        return None
    except Exception:
        return None


def _mock_mofa(n_factors: int, reason: str) -> dict:
    """Mock MOFA+ result when mofapy2 is unavailable and mocks are explicitly allowed."""
    mock_factors = [
        {
            "factor_id": i + 1,
            "r2_scRNA": round(0.15 - i * 0.01, 3),
            "r2_scATAC": round(0.08 - i * 0.005, 3),
        }
        for i in range(min(n_factors, 5))
    ]
    return {
        "status": "success",
        "n_factors": int(n_factors),
        "n_cells": 5000,
        "top_factors": mock_factors,
        "variance_explained": {"scRNA": 0.42, "scATAC": 0.23},
        "factor1_top_features": [f"feature_{i}" for i in range(1, 9)],
        "technical_factor_flag": False,
        "technical_factor_id": None,
        "factor_scores": None,
        "output_path": None,
        "warnings": [
            f"Mock MOFA+ — install aria-integration-env. ({reason})",
        ],
        "note": f"Mock MOFA+ — {reason}",
    }


if __name__ == "__main__":
    run_script(integration_mofa)
