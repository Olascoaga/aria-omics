"""Pure NumPy/SciPy Leave-One-Biological-Unit-Out (LOBO) Jackknife Engine.

Zero imports from ``aria.agents.*``, ``litellm``, or ``pydantic`` so this module
executes hermetically inside every scientific Conda stack (``aria-rna-env``,
``aria-chromatin-env``, ``aria-hic-env``, ``aria-env``).

Key statistical guarantees:
1. **Biological-unit granularity strictly enforced:** Rejects duplicate unit IDs
   (unaggregated technical replicates) and refuses ``min(n_num, n_den) < 3``
   (or ``n_pairs < 3`` for paired-donor designs) with
   ``status="not_testable_insufficient_replicates"``.
2. **Smyth (2004) empirical-Bayes variance moderation on genome-wide matrices
   (``M >= 20`` features):** Borrows prior residual variance ``s0^2 = median(s_g^2)``
   with ``d0 = 4.0`` prior degrees of freedom so ``3 vs 3`` leave-one-out folds
   (``2 vs 3``, ``df_obs = 3 -> df_post = 7``) retain full genome-wide BH-FDR
   calibration without falsely flagging stable ``3v3`` discoveries, while still
   collapsing ``100%`` when a single outlier donor/sample is held out. Small
   matrices (``M < 20``) use exact unmoderated Welch / paired t-tests.
3. **Full-universe BH-FDR with bounded IPC serialization:** Multiplicity
   adjustment runs across ALL ``M`` tested features on every jackknife fold;
   per-feature records are emitted for ``target_features`` (capped at
   ``max_reported_features``) so JSON IPC remains compact.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Deterministic Benjamini-Hochberg adjusted p-values preserving NaNs."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    m = int(valid.sum())
    if m == 0:
        return out
    pv = np.clip(p[valid], 0.0, 1.0)
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    adj = ranked * float(m) / np.arange(1, m + 1, dtype=float)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    out[valid] = adj[inv]
    return out


def _unpaired_fold_stats(
    mat: np.ndarray,
    n_idx: list[int],
    d_idx: list[int],
    moderate: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(lfc, pvals, padj_bh)`` for one unpaired fold ``n_idx`` vs ``d_idx``."""
    from scipy import stats

    x_num = mat[:, n_idx]
    x_den = mat[:, d_idx]
    n1 = x_num.shape[1]
    n2 = x_den.shape[1]
    lfc = np.mean(x_num, axis=1) - np.mean(x_den, axis=1)

    if not moderate:
        with np.errstate(divide="ignore", invalid="ignore"):
            _, pvals = stats.ttest_ind(x_num, x_den, axis=1, equal_var=False)
        pvals = np.where(np.isfinite(pvals), pvals, 1.0)
        return lfc, pvals, bh_adjust(pvals)

    v1 = np.var(x_num, axis=1, ddof=1)
    v2 = np.var(x_den, axis=1, ddof=1)
    df_obs = float(n1 + n2 - 2)
    s2_pooled = ((n1 - 1) * v1 + (n2 - 1) * v2) / max(1.0, df_obs)
    pos = s2_pooled[np.isfinite(s2_pooled) & (s2_pooled > 1e-12)]
    s0_sq = float(np.median(pos)) if pos.size > 0 else 1e-4
    d0 = 4.0
    s2_post = (d0 * s0_sq + df_obs * s2_pooled) / (d0 + df_obs)
    se = np.sqrt(np.maximum(s2_post * (1.0 / n1 + 1.0 / n2), 1e-15))
    t_stat = lfc / se
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df_obs + d0)
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    return lfc, pvals, bh_adjust(pvals)


def _paired_fold_stats(
    diff_mat: np.ndarray,
    pair_indices: list[int],
    moderate: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(lfc, pvals, padj_bh)`` for one paired donor fold."""
    from scipy import stats

    sub = diff_mat[:, pair_indices]
    k = sub.shape[1]
    lfc = np.mean(sub, axis=1)
    df_obs = float(k - 1)
    v_diff = np.var(sub, axis=1, ddof=1)

    if not moderate:
        se = np.sqrt(np.maximum(v_diff / k, 1e-15))
        with np.errstate(divide="ignore", invalid="ignore"):
            t_stat = np.where(v_diff > 0, lfc / se, np.where(np.abs(lfc) > 0, np.inf, 0.0))
        pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df_obs)
        pvals = np.where(np.isfinite(pvals), pvals, 1.0)
        return lfc, pvals, bh_adjust(pvals)

    pos = v_diff[np.isfinite(v_diff) & (v_diff > 1e-12)]
    s0_sq = float(np.median(pos)) if pos.size > 0 else 1e-4
    d0 = 4.0
    s2_post = (d0 * s0_sq + df_obs * v_diff) / (d0 + df_obs)
    se = np.sqrt(np.maximum(s2_post / k, 1e-15))
    t_stat = lfc / se
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df_obs + d0)
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    return lfc, pvals, bh_adjust(pvals)


def compute_lobo_jackknife(
    feature_ids: list[str],
    log_norm_matrix: np.ndarray,
    unit_ids: list[str],
    unit_to_group: dict[str, str],
    *,
    numerator: str,
    denominator: str,
    alpha: float = 0.05,
    target_features: list[str] | None = None,
    primary_significant_features: set[str] | list[str] | None = None,
    unit_to_donor: dict[str, str] | None = None,
    max_reported_features: int | None = None,
) -> dict[str, Any]:
    """Compute Leave-One-Biological-Unit-Out (LOBO) jackknife stability."""
    mat = np.asarray(log_norm_matrix, dtype=float)
    n_features_universe = len(feature_ids)
    if mat.ndim != 2 or mat.shape[0] != n_features_universe or mat.shape[1] != len(unit_ids):
        raise ValueError("log_norm_matrix shape must match (len(feature_ids), len(unit_ids))")
    if len(set(unit_ids)) != len(unit_ids):
        raise ValueError("unit_ids must represent unique biological units (aggregate technical replicates first)")

    moderate = bool(n_features_universe >= 20)

    # Check paired-donor mode first if unit_to_donor was provided
    paired_mode = False
    donor_pairs: list[tuple[str, int, int]] = []
    if unit_to_donor:
        num_by_donor: dict[str, int] = {}
        den_by_donor: dict[str, int] = {}
        for col_i, uid in enumerate(unit_ids):
            grp = unit_to_group.get(uid)
            donor = str(unit_to_donor.get(uid) or "")
            if not donor:
                continue
            if grp == numerator and donor not in num_by_donor:
                num_by_donor[donor] = col_i
            elif grp == denominator and donor not in den_by_donor:
                den_by_donor[donor] = col_i
        shared_donors = sorted(set(num_by_donor) & set(den_by_donor))
        if shared_donors:
            paired_mode = True
            donor_pairs = [(d, num_by_donor[d], den_by_donor[d]) for d in shared_donors]

    if paired_mode:
        n_pairs = len(donor_pairs)
        if n_pairs < 3:
            return {
                "status": "not_testable_insufficient_replicates",
                "design_mode": "paired_donor",
                "biological_unit_enforced": True,
                "numerator": numerator,
                "denominator": denominator,
                "n_units_numerator": n_pairs,
                "n_units_denominator": n_pairs,
                "minimum_required_per_arm": 3,
                "n_features_universe": n_features_universe,
                "per_feature": {},
            }
        diff_mat = np.stack(
            [mat[:, num_c] - mat[:, den_c] for _, num_c, den_c in donor_pairs],
            axis=1,
        )
        all_pair_idx = list(range(n_pairs))
        full_lfc, _, full_padj = _paired_fold_stats(diff_mat, all_pair_idx, moderate)
        folds = []
        for drop_p in range(n_pairs):
            sub_idx = [p for p in all_pair_idx if p != drop_p]
            f_lfc, _, f_padj = _paired_fold_stats(diff_mat, sub_idx, moderate)
            folds.append((donor_pairs[drop_p][0], f_lfc, f_padj))
        n_num_out = n_pairs
        n_den_out = n_pairs
        design_mode = "paired_donor"
    else:
        num_cols = [i for i, u in enumerate(unit_ids) if unit_to_group.get(u) == numerator]
        den_cols = [i for i, u in enumerate(unit_ids) if unit_to_group.get(u) == denominator]
        n_num_out = len(num_cols)
        n_den_out = len(den_cols)
        if min(n_num_out, n_den_out) < 3:
            return {
                "status": "not_testable_insufficient_replicates",
                "design_mode": "unpaired_biological_replicate",
                "biological_unit_enforced": True,
                "numerator": numerator,
                "denominator": denominator,
                "n_units_numerator": n_num_out,
                "n_units_denominator": n_den_out,
                "minimum_required_per_arm": 3,
                "n_features_universe": n_features_universe,
                "per_feature": {},
            }
        full_lfc, _, full_padj = _unpaired_fold_stats(mat, num_cols, den_cols, moderate)
        folds = []
        for drop_i in num_cols:
            sub_num = [c for c in num_cols if c != drop_i]
            f_lfc, _, f_padj = _unpaired_fold_stats(mat, sub_num, den_cols, moderate)
            folds.append((unit_ids[drop_i], f_lfc, f_padj))
        for drop_i in den_cols:
            sub_den = [c for c in den_cols if c != drop_i]
            f_lfc, _, f_padj = _unpaired_fold_stats(mat, num_cols, sub_den, moderate)
            folds.append((unit_ids[drop_i], f_lfc, f_padj))
        design_mode = "unpaired_biological_replicate"

    fold_lfcs = np.stack([f[1] for f in folds], axis=1)
    fold_padjs = np.stack([f[2] for f in folds], axis=1)
    fold_units = [f[0] for f in folds]

    id_to_idx = {str(fid): idx for idx, fid in enumerate(feature_ids)}
    if target_features is not None:
        eval_ids = [str(f) for f in target_features if str(f) in id_to_idx]
    else:
        eval_ids = [str(f) for f in feature_ids]

    primary_sig_set = (
        {str(x) for x in primary_significant_features}
        if primary_significant_features is not None
        else None
    )

    n_stable = 0
    n_single_unit_driven = 0
    per_feature: dict[str, dict[str, Any]] = {}

    for rank_i, fid in enumerate(eval_ids):
        row_idx = id_to_idx[fid]
        f_lfc_val = float(full_lfc[row_idx])
        f_padj_val = float(full_padj[row_idx])
        ref_sign = float(np.sign(f_lfc_val))
        row_lfcs = fold_lfcs[row_idx, :]
        row_padjs = fold_padjs[row_idx, :]

        sign_frac = (
            float(np.mean(np.sign(row_lfcs) == ref_sign))
            if ref_sign != 0.0
            else 0.0
        )
        q75, q25 = np.percentile(row_lfcs, [75, 25])
        worst_idx = int(np.argmax(row_padjs))
        max_padj = float(row_padjs[worst_idx])
        min_abs_lfc = float(np.min(np.abs(row_lfcs)))
        abs_full = abs(f_lfc_val)
        severe_effect_collapse = bool(abs_full > 0.25 and min_abs_lfc < 0.35 * abs_full)

        if primary_sig_set is not None:
            is_sig = fid in primary_sig_set
            # For DESeq2-called hits, require 100% sign retention, <65% single-unit
            # LFC collapse, AND (when full log-linear padj < alpha) max_padj < max(alpha, 0.10).
            padj_collapse = bool((f_padj_val < alpha) and (max_padj >= max(alpha * 2.0, 0.10)))
            single_unit_driven = bool(
                is_sig and (sign_frac < 1.0 or severe_effect_collapse or padj_collapse)
            ) or (not is_sig and sign_frac < 1.0)
        else:
            was_sig = bool(f_padj_val < alpha)
            single_unit_driven = bool(
                (was_sig and max_padj >= alpha) or sign_frac < 1.0 or severe_effect_collapse
            )

        if single_unit_driven:
            n_single_unit_driven += 1
        else:
            n_stable += 1

        if max_reported_features is None or rank_i < int(max_reported_features):
            per_feature[fid] = {
                "full_lfc": round(f_lfc_val, 6),
                "full_padj": float(f_padj_val),
                "n_jackknife_folds": len(folds),
                "sign_consistent_fraction": round(sign_frac, 4),
                "lfc_min": round(float(np.min(row_lfcs)), 6),
                "lfc_max": round(float(np.max(row_lfcs)), 6),
                "lfc_iqr": round(float(q75 - q25), 6),
                "max_jackknife_padj": float(max_padj),
                "most_influential_unit": fold_units[worst_idx],
                "single_unit_driven": single_unit_driven,
            }

    n_eval = len(eval_ids)
    return {
        "status": "measured",
        "design_mode": design_mode,
        "variance_moderation": "empirical_bayes_smyth2004" if moderate else "exact_welch_or_paired_t",
        "biological_unit_enforced": True,
        "numerator": numerator,
        "denominator": denominator,
        "n_units_numerator": n_num_out,
        "n_units_denominator": n_den_out,
        "n_jackknife_folds": len(folds),
        "n_features_universe": n_features_universe,
        "n_features_evaluated": n_eval,
        "n_features_lobo_stable": n_stable,
        "n_single_unit_driven": n_single_unit_driven,
        "lobo_stable_fraction": round(n_stable / n_eval, 4) if n_eval > 0 else None,
        "per_feature_truncated": bool(
            max_reported_features is not None and n_eval > int(max_reported_features)
        ),
        "per_feature": per_feature,
    }


def compute_counts_lobo_jackknife(
    counts_df,
    meta_df,
    condition_col: str,
    numerator: str,
    denominator: str,
    *,
    alpha: float = 0.05,
    target_features: list[str] | None = None,
    primary_significant_features: set[str] | list[str] | None = None,
    size_factors: dict[str, float] | Any | None = None,
    paired_unit_col: str | None = None,
    max_reported_features: int = 200,
) -> dict[str, Any]:
    """Normalize a ``features x biological_units`` integer count DataFrame and run LOBO jackknife."""
    samples = [str(c) for c in counts_df.columns]
    meta_aligned = meta_df.loc[counts_df.columns]
    unit_to_group = {
        str(idx): str(meta_aligned.at[idx, condition_col])
        for idx in counts_df.columns
    }

    raw = counts_df.to_numpy(dtype=float)
    normalization = "log2_cpm_plus_1"
    if size_factors is not None:
        if isinstance(size_factors, dict):
            sf_vec = np.array([float(size_factors.get(s, 1.0)) for s in samples], dtype=float)
        else:
            sf_vec = np.asarray(size_factors, dtype=float).ravel()
        if sf_vec.shape[0] == len(samples) and np.all(np.isfinite(sf_vec)) and np.all(sf_vec > 0):
            scaled = raw / sf_vec[None, :]
            log_mat = np.log2(np.maximum(scaled, 0.0) + 1.0)
            normalization = "spikein_or_external_size_factor_log2"
        else:
            lib_sizes = np.maximum(raw.sum(axis=0), 1.0)
            log_mat = np.log2((raw / lib_sizes[None, :]) * 1e6 + 1.0)
    else:
        lib_sizes = np.maximum(raw.sum(axis=0), 1.0)
        log_mat = np.log2((raw / lib_sizes[None, :]) * 1e6 + 1.0)

    unit_to_donor = None
    if paired_unit_col and paired_unit_col in meta_aligned.columns:
        unit_to_donor = {
            str(idx): str(meta_aligned.at[idx, paired_unit_col])
            for idx in counts_df.columns
        }

    out = compute_lobo_jackknife(
        feature_ids=[str(i) for i in counts_df.index],
        log_norm_matrix=log_mat,
        unit_ids=samples,
        unit_to_group=unit_to_group,
        numerator=str(numerator),
        denominator=str(denominator),
        alpha=float(alpha),
        target_features=target_features,
        primary_significant_features=primary_significant_features,
        unit_to_donor=unit_to_donor,
        max_reported_features=max_reported_features,
    )
    out["normalization"] = normalization
    return out
