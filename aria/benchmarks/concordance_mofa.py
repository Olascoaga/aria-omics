"""V48 Slice I4 — External concordance scoring for paired MOFA+ decompositions.

Pure NumPy/SciPy scoring layer comparing ARIA's paired MOFA+ decomposition
(`aria/scripts/integration_mofa.py`) against an external reference execution
(e.g., `muon.tl.mofa` / standalone `mofapy2`).

Latent factor models are identifiable only up to factor permutation, sign flip,
and (when eigenvalues/variances are near-degenerate) orthogonal rotation within
a degenerate subspace. Accordingly, this module evaluates:
  1. Rectangular Hungarian 1-to-1 matching (`linear_sum_assignment` on `-|r|`)
     restricted to genuinely active factors under a hybrid absolute + view-relative
     variance threshold, weighted strictly by the reference model's R^2.
  2. Rotation-invariant subspace concordance via canonical correlations
     (principal cosines between orthonormal QR bases of active factor scores).
  3. Per-view feature loading concordance (Pearson r, Spearman rho, sign
     consistency `sign(r_Z) == sign(r_W)`, and top-k feature Jaccard) on shared
     feature universes.
  4. Total variance explained (`R^2_total`) agreement per modality.
  5. Seed determinism (`max |delta Z| == 0` on repeated fixed-seed runs) and
     cross-seed variational convergence ceiling (`seed_a` vs `seed_b`).

Honesty contract (ADR-002): returns `status="not_run"` with an explicit reason
when external reference outputs are absent; never fabricates concordance values.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np


def _rankdata_average(values: np.ndarray) -> np.ndarray:
    """Average-tie 1-based rank transform in pure NumPy."""
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
        ranks[order[i:j]] = 0.5 * (i + (j - 1)) + 1.0
        i = j
    return ranks


def pearson_correlation(
    a: Sequence[float] | np.ndarray,
    b: Sequence[float] | np.ndarray,
    weights: Sequence[float] | np.ndarray | None = None,
) -> float | None:
    """Standard or inclusion-weighted Pearson r over finite pairs; None if undefined."""
    xa = np.asarray(a, dtype=np.float64).ravel()
    xb = np.asarray(b, dtype=np.float64).ravel()
    if xa.size != xb.size or xa.size < 2:
        return None

    if weights is not None:
        w = np.asarray(weights, dtype=np.float64).ravel()
        if w.size != xa.size:
            raise ValueError(
                f"weights length ({w.size}) must match input length ({xa.size})"
            )
        mask = np.isfinite(xa) & np.isfinite(xb) & np.isfinite(w) & (w > 0.0)
        if int(mask.sum()) < 2:
            return None
        xa = xa[mask]
        xb = xb[mask]
        w = w[mask]
        w_sum = float(np.sum(w))
        if w_sum <= 0.0:
            return None
        w_norm = w / w_sum
        ma = float(np.sum(w_norm * xa))
        mb = float(np.sum(w_norm * xb))
        da = xa - ma
        db = xb - mb
        va = float(np.sum(w_norm * da * da))
        vb = float(np.sum(w_norm * db * db))
        if va <= 1e-24 or vb <= 1e-24:
            return None
        cov = float(np.sum(w_norm * da * db))
        r_val = cov / float(np.sqrt(va * vb))
    else:
        mask = np.isfinite(xa) & np.isfinite(xb)
        if int(mask.sum()) < 2:
            return None
        xa = xa[mask]
        xb = xb[mask]
        sa = float(np.std(xa))
        sb = float(np.std(xb))
        if sa <= 1e-15 or sb <= 1e-15:
            return None
        r_val = float(np.corrcoef(xa, xb)[0, 1])

    if not np.isfinite(r_val):
        return None
    return float(np.clip(r_val, -1.0, 1.0))


def spearman_correlation(
    a: Sequence[float] | np.ndarray,
    b: Sequence[float] | np.ndarray,
) -> float | None:
    """Spearman rank correlation in pure NumPy; None when constant or <2 items."""
    xa = np.asarray(a, dtype=np.float64).ravel()
    xb = np.asarray(b, dtype=np.float64).ravel()
    if xa.size != xb.size or xa.size < 2:
        return None
    mask = np.isfinite(xa) & np.isfinite(xb)
    if int(mask.sum()) < 2:
        return None
    ra = _rankdata_average(xa[mask])
    rb = _rankdata_average(xb[mask])
    return pearson_correlation(ra, rb)


def is_factor_active_in_view(
    r2_vk: float,
    r2_view_sum: float,
    *,
    active_r2_min: float = 0.01,
    rel_r2_min: float = 0.05,
    min_floor_r2: float = 0.0008,
) -> bool:
    """Hybrid absolute + view-relative variance criterion.

    Prevents sparse views (such as scATAC TF-IDF peaks, where total R^2 across
    all factors is typically 0.3%–5%, an order of magnitude lower than scRNA
    HVGs) from having genuine shared or view-specific factors misclassified as
    inactive while still enforcing a hard noise floor (`min_floor_r2`).
    """
    val = max(0.0, float(r2_vk))
    v_tot = max(0.0, float(r2_view_sum))
    if val >= float(active_r2_min):
        return True
    if (
        v_tot >= max(float(min_floor_r2), 0.002)
        and val >= float(min_floor_r2)
        and (val / v_tot) >= float(rel_r2_min)
    ):
        return True
    return False


def _active_factor_indices(
    scores: np.ndarray,
    r2_per_view: Mapping[str, Sequence[float]] | None,
    *,
    active_r2_min: float = 0.01,
    rel_r2_min: float = 0.05,
    min_floor_r2: float = 0.0008,
) -> tuple[list[int], np.ndarray]:
    """Return 0-based active factor indices and per-factor positive R^2 sums."""
    mat = np.asarray(scores, dtype=np.float64)
    if mat.ndim != 2:
        raise ValueError(f"Expected 2D factor score matrix, got shape {mat.shape}")
    k_tot = mat.shape[1]
    stds = np.std(mat, axis=0) if mat.shape[0] >= 2 else np.zeros(k_tot, dtype=np.float64)

    if not r2_per_view:
        active = [k for k in range(k_tot) if np.isfinite(stds[k]) and stds[k] > 1e-12]
        return active, np.ones(k_tot, dtype=np.float64)

    r2_sum_per_k = np.zeros(k_tot, dtype=np.float64)
    view_arrays: dict[str, np.ndarray] = {}
    view_totals: dict[str, float] = {}
    for v_name, vals in r2_per_view.items():
        arr = np.maximum(0.0, np.asarray(list(vals), dtype=np.float64).ravel())
        if arr.size != k_tot:
            raise ValueError(
                f"r2_per_view[{v_name!r}] length {arr.size} != n_factors {k_tot}"
            )
        view_arrays[str(v_name)] = arr
        view_totals[str(v_name)] = float(np.sum(arr))
        r2_sum_per_k += arr

    active = []
    for k in range(k_tot):
        if not np.isfinite(stds[k]) or stds[k] <= 1e-12:
            continue
        if any(
            is_factor_active_in_view(
                float(view_arrays[v][k]),
                view_totals[v],
                active_r2_min=active_r2_min,
                rel_r2_min=rel_r2_min,
                min_floor_r2=min_floor_r2,
            )
            for v in view_arrays
        ):
            active.append(k)
    return active, r2_sum_per_k


def _linear_assignment_maximize(score_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Solve rectangular maximum-weight assignment via SciPy (greedy fallback)."""
    try:
        from scipy.optimize import linear_sum_assignment

        return linear_sum_assignment(-np.asarray(score_matrix, dtype=np.float64))
    except ImportError:
        mat = np.asarray(score_matrix, dtype=np.float64).copy()
        n_rows, n_cols = mat.shape
        n_pick = min(n_rows, n_cols)
        used_r: set[int] = set()
        used_c: set[int] = set()
        rows: list[int] = []
        cols: list[int] = []
        for _ in range(n_pick):
            best_val = -1.0
            best_rc = (-1, -1)
            for r in range(n_rows):
                if r in used_r:
                    continue
                for c in range(n_cols):
                    if c in used_c:
                        continue
                    if mat[r, c] > best_val:
                        best_val = float(mat[r, c])
                        best_rc = (r, c)
            if best_rc[0] < 0:
                break
            used_r.add(best_rc[0])
            used_c.add(best_rc[1])
            rows.append(best_rc[0])
            cols.append(best_rc[1])
        order = np.argsort(rows)
        return np.asarray(rows, dtype=int)[order], np.asarray(cols, dtype=int)[order]


def _max_offdiag_abs_corr(mat: np.ndarray) -> float | None:
    """Maximum |Pearson r| between distinct columns of `mat`."""
    arr = np.asarray(mat, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] < 2:
        return None
    k = arr.shape[1]
    best = 0.0
    found = False
    for i in range(k):
        for j in range(i + 1, k):
            r = pearson_correlation(arr[:, i], arr[:, j])
            if r is not None:
                found = True
                if abs(r) > best:
                    best = abs(r)
    return round(float(best), 4) if found else None


def subspace_canonical_concordance(
    aria_scores: np.ndarray,
    ref_scores: np.ndarray,
) -> dict[str, Any]:
    """Rotation-invariant principal cosines between active factor subspaces."""
    za = np.asarray(aria_scores, dtype=np.float64)
    zr = np.asarray(ref_scores, dtype=np.float64)
    if za.ndim != 2 or zr.ndim != 2 or za.shape[0] != zr.shape[0]:
        raise ValueError(
            f"Incompatible score shapes for subspace concordance: {za.shape} vs {zr.shape}"
        )
    if za.shape[0] < 2 or za.shape[1] == 0 or zr.shape[1] == 0:
        return {
            "n_dims_aria": 0,
            "n_dims_ref": 0,
            "canonical_correlations": [],
            "mean_canonical_corr": None,
            "min_top3_canonical_corr": None,
        }

    za_c = za - np.mean(za, axis=0, keepdims=True)
    zr_c = zr - np.mean(zr, axis=0, keepdims=True)
    keep_a = np.std(za_c, axis=0) > 1e-12
    keep_r = np.std(zr_c, axis=0) > 1e-12
    za_c = za_c[:, keep_a]
    zr_c = zr_c[:, keep_r]

    if za_c.shape[1] == 0 or zr_c.shape[1] == 0:
        return {
            "n_dims_aria": int(za_c.shape[1]),
            "n_dims_ref": int(zr_c.shape[1]),
            "canonical_correlations": [],
            "mean_canonical_corr": None,
            "min_top3_canonical_corr": None,
        }

    qa, _ = np.linalg.qr(za_c, mode="reduced")
    qr, _ = np.linalg.qr(zr_c, mode="reduced")
    # Guard against rank-deficient columns after QR.
    sv = np.linalg.svd(qa.T @ qr, compute_uv=False)
    sv = np.clip(sv, 0.0, 1.0)
    top3 = sv[: min(3, sv.size)]

    return {
        "n_dims_aria": int(qa.shape[1]),
        "n_dims_ref": int(qr.shape[1]),
        "canonical_correlations": [round(float(x), 4) for x in sv],
        "mean_canonical_corr": round(float(np.mean(sv)), 4) if sv.size else None,
        "min_top3_canonical_corr": round(float(np.min(top3)), 4) if top3.size else None,
    }


def match_factors_hungarian(
    aria_scores: np.ndarray,
    ref_scores: np.ndarray,
    *,
    aria_r2_per_view: Mapping[str, Sequence[float]] | None = None,
    ref_r2_per_view: Mapping[str, Sequence[float]] | None = None,
    cell_weights: Sequence[float] | np.ndarray | None = None,
    active_r2_min: float = 0.01,
    rel_r2_min: float = 0.05,
    min_floor_r2: float = 0.0008,
) -> dict[str, Any]:
    """Optimal 1-to-1 matching of active factors invariant to permutation & sign."""
    za = np.asarray(aria_scores, dtype=np.float64)
    zr = np.asarray(ref_scores, dtype=np.float64)
    if za.ndim != 2 or zr.ndim != 2 or za.shape[0] != zr.shape[0]:
        raise ValueError(
            f"Score matrices must be 2D with matching n_cells: {za.shape} vs {zr.shape}"
        )

    act_a, _ = _active_factor_indices(
        za,
        aria_r2_per_view,
        active_r2_min=active_r2_min,
        rel_r2_min=rel_r2_min,
        min_floor_r2=min_floor_r2,
    )
    act_r, ref_r2_sums = _active_factor_indices(
        zr,
        ref_r2_per_view,
        active_r2_min=active_r2_min,
        rel_r2_min=rel_r2_min,
        min_floor_r2=min_floor_r2,
    )

    ka = len(act_a)
    kr = len(act_r)
    if ka == 0 or kr == 0:
        return {
            "n_total_aria": int(za.shape[1]),
            "n_total_ref": int(zr.shape[1]),
            "n_active_aria": ka,
            "n_active_ref": kr,
            "n_matched": 0,
            "n_unmatched_aria": ka,
            "n_unmatched_ref": kr,
            "matches": [],
            "mean_matched_abs_r": None,
            "median_matched_abs_r": None,
            "min_matched_abs_r": None,
            "ref_r2_weighted_abs_r": None,
            "max_offdiag_abs_r_aria": _max_offdiag_abs_corr(za[:, act_a]) if ka >= 2 else None,
            "max_offdiag_abs_r_ref": _max_offdiag_abs_corr(zr[:, act_r]) if kr >= 2 else None,
        }

    r_mat = np.zeros((ka, kr), dtype=np.float64)
    for ia, col_a in enumerate(act_a):
        for ir, col_r in enumerate(act_r):
            r_val = pearson_correlation(za[:, col_a], zr[:, col_r], weights=cell_weights)
            r_mat[ia, ir] = float(r_val) if r_val is not None else 0.0

    row_ind, col_ind = _linear_assignment_maximize(np.abs(r_mat))
    matches: list[dict[str, Any]] = []
    for ia, ir in zip(row_ind.tolist(), col_ind.tolist()):
        col_a = act_a[ia]
        col_r = act_r[ir]
        p_r = float(r_mat[ia, ir])
        s_rho = spearman_correlation(za[:, col_a], zr[:, col_r])
        sign = 1 if p_r >= 0.0 else -1
        matches.append(
            {
                "aria_factor_idx": int(col_a),
                "aria_factor_id": int(col_a + 1),
                "ref_factor_idx": int(col_r),
                "ref_factor_id": int(col_r + 1),
                "pearson_r": round(p_r, 4),
                "abs_pearson_r": round(abs(p_r), 4),
                "spearman_rho": round(float(s_rho), 4) if s_rho is not None else None,
                "abs_spearman_rho": (
                    round(abs(float(s_rho)), 4) if s_rho is not None else None
                ),
                "sign": int(sign),
                "ref_r2_sum": round(float(ref_r2_sums[col_r]), 4),
            }
        )

    matches.sort(key=lambda m: (-m["ref_r2_sum"], -m["abs_pearson_r"], m["ref_factor_id"]))
    abs_rs = np.asarray([m["abs_pearson_r"] for m in matches], dtype=np.float64)
    ref_w = np.asarray([max(0.0, m["ref_r2_sum"]) for m in matches], dtype=np.float64)
    if float(np.sum(ref_w)) > 1e-12:
        weighted_abs_r = float(np.sum(ref_w * abs_rs) / np.sum(ref_w))
    else:
        weighted_abs_r = float(np.mean(abs_rs))

    return {
        "n_total_aria": int(za.shape[1]),
        "n_total_ref": int(zr.shape[1]),
        "n_active_aria": int(ka),
        "n_active_ref": int(kr),
        "n_matched": int(len(matches)),
        "n_unmatched_aria": int(max(0, ka - len(matches))),
        "n_unmatched_ref": int(max(0, kr - len(matches))),
        "matches": matches,
        "mean_matched_abs_r": round(float(np.mean(abs_rs)), 4),
        "median_matched_abs_r": round(float(np.median(abs_rs)), 4),
        "min_matched_abs_r": round(float(np.min(abs_rs)), 4),
        "ref_r2_weighted_abs_r": round(weighted_abs_r, 4),
        "max_offdiag_abs_r_aria": _max_offdiag_abs_corr(za[:, act_a]),
        "max_offdiag_abs_r_ref": _max_offdiag_abs_corr(zr[:, act_r]),
    }


def _canonical_modality_family(name: str) -> str:
    low = str(name).strip().lower()
    if any(k in low for k in ("rna", "gex", "expr")):
        return "rna"
    if any(k in low for k in ("atac", "peak", "chrom", "acc")):
        return "atac"
    return low


def _align_view_keys(
    aria_map: Mapping[str, Any],
    ref_map: Mapping[str, Any],
) -> list[tuple[str, str, str]]:
    """Return `(canonical_name, aria_key, ref_key)` pairs across two view maps."""
    pairs: list[tuple[str, str, str]] = []
    used_ref: set[str] = set()
    for ak in aria_map:
        if ak in ref_map:
            pairs.append((str(ak), str(ak), str(ak)))
            used_ref.add(str(ak))
    for ak in aria_map:
        if any(p[1] == str(ak) for p in pairs):
            continue
        fam_a = _canonical_modality_family(ak)
        for rk in ref_map:
            if str(rk) in used_ref:
                continue
            if _canonical_modality_family(rk) == fam_a:
                pairs.append((fam_a, str(ak), str(rk)))
                used_ref.add(str(rk))
                break
    return pairs


def score_loading_concordance(
    aria_loadings: Mapping[str, tuple[Sequence[str], np.ndarray]],
    ref_loadings: Mapping[str, tuple[Sequence[str], np.ndarray]],
    matches: Sequence[Mapping[str, Any]],
    *,
    min_shared_features: int = 50,
    top_k: int = 50,
) -> dict[str, Any]:
    """Evaluate sign-aligned feature weight concordance per view."""
    views_out: dict[str, Any] = {}
    not_computed: list[str] = []

    for canon_view, ak, rk in _align_view_keys(aria_loadings, ref_loadings):
        a_feats, a_w = aria_loadings[ak]
        r_feats, r_w = ref_loadings[rk]
        a_w_arr = np.asarray(a_w, dtype=np.float64)
        r_w_arr = np.asarray(r_w, dtype=np.float64)

        a_idx_map = {str(f): i for i, f in enumerate(a_feats)}
        r_idx_map = {str(f): i for i, f in enumerate(r_feats)}
        shared = [str(f) for f in a_feats if str(f) in r_idx_map]
        n_shared = len(shared)

        if n_shared < int(min_shared_features):
            views_out[canon_view] = {
                "status": "insufficient_feature_overlap",
                "n_shared_features": int(n_shared),
                "min_required": int(min_shared_features),
            }
            not_computed.append(f"loading_concordance:{canon_view}")
            continue

        a_take = [a_idx_map[f] for f in shared]
        r_take = [r_idx_map[f] for f in shared]
        wa_sub = a_w_arr[a_take, :]
        wr_sub = r_w_arr[r_take, :]
        k_eff = max(1, min(int(top_k), n_shared // 2))

        per_factor: list[dict[str, Any]] = []
        prs: list[float] = []
        srhos: list[float] = []
        jaccs: list[float] = []
        sign_ok: list[bool] = []

        for m in matches:
            ia = int(m["aria_factor_idx"])
            ir = int(m["ref_factor_idx"])
            sign_z = float(m.get("sign", 1))
            va = wa_sub[:, ia]
            vr = sign_z * wr_sub[:, ir]

            p_r = pearson_correlation(va, vr)
            s_rho = spearman_correlation(va, vr)

            top_a = set(np.argsort(np.abs(va))[::-1][:k_eff].tolist())
            top_r = set(np.argsort(np.abs(vr))[::-1][:k_eff].tolist())
            union = len(top_a | top_r)
            jacc = float(len(top_a & top_r)) / float(union) if union else 0.0
            is_sign_ok = bool(p_r is not None and p_r > 0.0)

            if p_r is not None:
                prs.append(float(p_r))
                sign_ok.append(is_sign_ok)
            if s_rho is not None:
                srhos.append(float(s_rho))
            jaccs.append(jacc)

            per_factor.append(
                {
                    "aria_factor_id": int(m["aria_factor_id"]),
                    "ref_factor_id": int(m["ref_factor_id"]),
                    "pearson_r": round(float(p_r), 4) if p_r is not None else None,
                    "spearman_rho": round(float(s_rho), 4) if s_rho is not None else None,
                    "top_k_jaccard": round(jacc, 4),
                    "sign_consistent": is_sign_ok,
                }
            )

        views_out[canon_view] = {
            "status": "success",
            "n_shared_features": int(n_shared),
            "top_k_evaluated": int(k_eff),
            "mean_pearson_r": round(float(np.mean(prs)), 4) if prs else None,
            "mean_spearman_rho": round(float(np.mean(srhos)), 4) if srhos else None,
            "mean_top_k_jaccard": round(float(np.mean(jaccs)), 4) if jaccs else None,
            "sign_consistency_rate": (
                round(float(np.mean(sign_ok)), 4) if sign_ok else None
            ),
            "per_factor": per_factor,
        }

    return {
        "views": views_out,
        "metrics_not_computed": not_computed,
    }


def score_mofa_seed_stability(
    score_runs: Sequence[np.ndarray],
    r2_runs: Sequence[Mapping[str, float]] | None = None,
) -> dict[str, Any]:
    """Score exact fixed-seed determinism or cross-seed variational stability."""
    n_runs = len(score_runs)
    if n_runs < 2:
        return {
            "status": "not_run",
            "n_runs": int(n_runs),
            "max_elementwise_abs_diff": None,
            "mean_pairwise_matched_abs_r": None,
            "min_pairwise_matched_abs_r": None,
            "mean_pairwise_canonical_corr": None,
            "max_r2_total_abs_diff": None,
        }

    elem_diffs: list[float] = []
    matched_rs: list[float] = []
    canon_rs: list[float] = []

    for i in range(n_runs):
        zi = np.asarray(score_runs[i], dtype=np.float64)
        for j in range(i + 1, n_runs):
            zj = np.asarray(score_runs[j], dtype=np.float64)
            if zi.shape == zj.shape and zi.size > 0:
                elem_diffs.append(float(np.max(np.abs(zi - zj))))
            if zi.shape[0] == zj.shape[0] and zi.shape[0] >= 2:
                mf = match_factors_hungarian(zi, zj)
                if mf["mean_matched_abs_r"] is not None:
                    matched_rs.append(float(mf["mean_matched_abs_r"]))
                sc = subspace_canonical_concordance(zi, zj)
                if sc["mean_canonical_corr"] is not None:
                    canon_rs.append(float(sc["mean_canonical_corr"]))

    max_r2_diff: float | None = None
    if r2_runs is not None and len(r2_runs) >= 2:
        r2_diffs: list[float] = []
        for i in range(len(r2_runs)):
            for j in range(i + 1, len(r2_runs)):
                for _, ak, rk in _align_view_keys(r2_runs[i], r2_runs[j]):
                    r2_diffs.append(abs(float(r2_runs[i][ak]) - float(r2_runs[j][rk])))
        if r2_diffs:
            max_r2_diff = round(float(max(r2_diffs)), 8)

    return {
        "status": "success",
        "n_runs": int(n_runs),
        "max_elementwise_abs_diff": (
            round(float(max(elem_diffs)), 8) if elem_diffs else None
        ),
        "mean_pairwise_matched_abs_r": (
            round(float(np.mean(matched_rs)), 6) if matched_rs else None
        ),
        "min_pairwise_matched_abs_r": (
            round(float(np.min(matched_rs)), 6) if matched_rs else None
        ),
        "mean_pairwise_canonical_corr": (
            round(float(np.mean(canon_rs)), 6) if canon_rs else None
        ),
        "max_r2_total_abs_diff": max_r2_diff,
    }


def score_mofa_concordance(
    *,
    aria_barcodes: Sequence[str],
    aria_scores: np.ndarray,
    aria_r2_total: Mapping[str, float],
    aria_r2_per_view: Mapping[str, Sequence[float]] | None = None,
    aria_loadings: Mapping[str, tuple[Sequence[str], np.ndarray]] | None = None,
    ref_barcodes: Sequence[str] | None = None,
    ref_scores: np.ndarray | None = None,
    ref_r2_total: Mapping[str, float] | None = None,
    ref_r2_per_view: Mapping[str, Sequence[float]] | None = None,
    ref_loadings: Mapping[str, tuple[Sequence[str], np.ndarray]] | None = None,
    cell_weights: Sequence[float] | np.ndarray | None = None,
    min_shared_features: int = 50,
    top_k: int = 50,
    active_r2_min: float = 0.01,
    rel_r2_min: float = 0.05,
    min_floor_r2: float = 0.0008,
    reference_name: str = "muon.tl.mofa",
) -> dict[str, Any]:
    """Compare ARIA MOFA+ outputs against an external reference decomposition."""
    if ref_barcodes is None or ref_scores is None:
        return {
            "status": "not_run",
            "reason": "External MOFA+ reference outputs not provided.",
            "reference_name": str(reference_name),
        }

    a_bcs = [str(x) for x in aria_barcodes]
    r_bcs = [str(x) for x in ref_barcodes]
    if len(set(a_bcs)) != len(a_bcs):
        raise ValueError("Duplicate barcodes detected in aria_barcodes.")
    if len(set(r_bcs)) != len(r_bcs):
        raise ValueError("Duplicate barcodes detected in ref_barcodes.")

    za_full = np.asarray(aria_scores, dtype=np.float64)
    zr_full = np.asarray(ref_scores, dtype=np.float64)
    if za_full.shape[0] != len(a_bcs):
        raise ValueError(
            f"aria_scores rows ({za_full.shape[0]}) != len(aria_barcodes) ({len(a_bcs)})"
        )
    if zr_full.shape[0] != len(r_bcs):
        raise ValueError(
            f"ref_scores rows ({zr_full.shape[0]}) != len(ref_barcodes) ({len(r_bcs)})"
        )

    a_idx = {bc: i for i, bc in enumerate(a_bcs)}
    r_idx = {bc: i for i, bc in enumerate(r_bcs)}
    shared_bcs = [bc for bc in a_bcs if bc in r_idx]
    if len(shared_bcs) < 10:
        return {
            "status": "error",
            "error_type": "InsufficientOverlap",
            "details": (
                f"Need >= 10 shared cells to score MOFA+ concordance; "
                f"found {len(shared_bcs)}."
            ),
            "n_shared_cells": int(len(shared_bcs)),
            "reference_name": str(reference_name),
        }

    a_take = [a_idx[bc] for bc in shared_bcs]
    r_take = [r_idx[bc] for bc in shared_bcs]
    za = za_full[a_take, :]
    zr = zr_full[r_take, :]
    w_aligned = None
    if cell_weights is not None:
        w_arr = np.asarray(cell_weights, dtype=np.float64).ravel()
        if w_arr.size != len(a_bcs):
            raise ValueError("cell_weights length must match aria_barcodes.")
        w_aligned = w_arr[a_take]

    matching = match_factors_hungarian(
        za,
        zr,
        aria_r2_per_view=aria_r2_per_view,
        ref_r2_per_view=ref_r2_per_view,
        cell_weights=w_aligned,
        active_r2_min=active_r2_min,
        rel_r2_min=rel_r2_min,
        min_floor_r2=min_floor_r2,
    )

    act_a, _ = _active_factor_indices(
        za,
        aria_r2_per_view,
        active_r2_min=active_r2_min,
        rel_r2_min=rel_r2_min,
        min_floor_r2=min_floor_r2,
    )
    act_r, _ = _active_factor_indices(
        zr,
        ref_r2_per_view,
        active_r2_min=active_r2_min,
        rel_r2_min=rel_r2_min,
        min_floor_r2=min_floor_r2,
    )
    subspace_active = subspace_canonical_concordance(
        za[:, act_a] if act_a else za[:, :0],
        zr[:, act_r] if act_r else zr[:, :0],
    )

    not_computed: list[str] = []
    var_conc: dict[str, Any] = {"per_view_abs_diff": {}, "max_abs_diff": None}
    if ref_r2_total:
        diffs: dict[str, float] = {}
        for canon_view, ak, rk in _align_view_keys(aria_r2_total, ref_r2_total):
            diffs[canon_view] = round(
                abs(float(aria_r2_total[ak]) - float(ref_r2_total[rk])), 4
            )
        var_conc = {
            "aria_r2_total": {str(k): round(float(v), 4) for k, v in aria_r2_total.items()},
            "ref_r2_total": {str(k): round(float(v), 4) for k, v in ref_r2_total.items()},
            "per_view_abs_diff": diffs,
            "max_abs_diff": round(float(max(diffs.values())), 4) if diffs else None,
        }
    else:
        not_computed.append("variance_explained_concordance")

    loading_conc: dict[str, Any] | None = None
    if aria_loadings is not None and ref_loadings is not None:
        loading_conc = score_loading_concordance(
            aria_loadings,
            ref_loadings,
            matching["matches"],
            min_shared_features=min_shared_features,
            top_k=top_k,
        )
        not_computed.extend(loading_conc.get("metrics_not_computed", []))
    else:
        not_computed.append("loading_concordance")

    return {
        "status": "success",
        "reference_name": str(reference_name),
        "n_aria_cells": int(len(a_bcs)),
        "n_ref_cells": int(len(r_bcs)),
        "n_shared_cells": int(len(shared_bcs)),
        "factor_matching": matching,
        "subspace_concordance": subspace_active,
        "variance_explained_concordance": var_conc,
        "loading_concordance": loading_conc,
        "metrics_not_computed": not_computed,
    }
