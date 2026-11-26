"""ARIA Hi-C / Micro-C QC, ICE Balancing, P(s) Decay & Replicate SCC (ADR-082).

validation_level: production (stack="hic").

All operations are strictly out-of-core (chunked pixel iteration and
chromosome-scoped ``cooltools`` expected curves). Never loads a full genome-wide
dense matrix into RAM and never fabricates fallback QC numbers unless explicit
test mock mode (``mocks_allowed``) is requested.

Outputs:
  * Exact ``n_valid_pairs`` (sum of ``count`` across pixels, not ``nnz``),
    ``cis_reads``, ``trans_reads``, ``cis_trans_ratio``,
    ``cis_short_lt20kb``, ``cis_long_ge20kb``, ``long_range_cis``.
  * Writable ICE-balanced ``.cool`` / ``.mcool`` paths in ``output_dir``
    (never mutates read-only user input files in place).
  * Distance-decay ``P(s)`` curve TSV (``cooltools.api.expected.expected_cis``)
    and log-log polymer contact decay slope in [20 kb, 2 Mb].
  * Multi-sample HiCRep stratum-adjusted correlation coefficient (``SCC``,
    Yang et al. 2017) matrix and ``mean_replicate_scc`` when >= 2 samples share
    a bin grid.
"""

from __future__ import annotations

import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts._base import mocks_allowed, run_script
from aria.scripts.hic_preprocess import _check_hic_stage_cache, _save_hic_stage_cache
from aria.utils.hic_constants import list_mcool_resolutions

STAGE_VERSION = 1
PIXEL_CHUNK_SIZE = 5_000_000
EXCLUDED_CHROMS = frozenset({"chrM", "MT", "chrMT", "M", "chrEBV", "EBV"})


def _nuclear_view_df(clr):
    import pandas as pd

    rows = []
    for ch in clr.chromnames:
        s = str(ch)
        if s in EXCLUDED_CHROMS or "random" in s.lower() or s.lower().startswith("chrun"):
            continue
        length = int(clr.chromsizes[ch])
        if length > 0:
            rows.append({"chrom": s, "start": 0, "end": length, "name": s})
    if not rows:
        for ch in clr.chromnames:
            length = int(clr.chromsizes[ch])
            if length > 0:
                rows.append({"chrom": str(ch), "start": 0, "end": length, "name": str(ch)})
    return pd.DataFrame(rows)


def _resolve_cooler_uri(path: str, resolution: int) -> tuple[str, int]:
    import cooler

    if path.lower().endswith(".mcool"):
        res_map = list_mcool_resolutions(path)
        if not res_map:
            raise ValueError(f"No resolutions inside .mcool: {path}")
        avail = sorted(res_map.keys())
        chosen = int(resolution) if int(resolution) in res_map else min(
            avail, key=lambda r: abs(r - int(resolution))
        )
        return f"{path}::{res_map[chosen]}", chosen
    clr = cooler.Cooler(path)
    return path, int(clr.binsize)


def _resolve_scc_uri(
    path: str,
    base_resolution: int,
    explicit_scc_resolution: int | None = None,
) -> tuple[str, int]:
    """Select the resolution used for HiCRep SCC (defaults .mcool to >= 40 kb to prevent Micro-C OOM)."""
    if path.lower().endswith(".mcool"):
        res_map = list_mcool_resolutions(path)
        if not res_map:
            raise ValueError(f"No resolutions inside .mcool: {path}")
        avail = sorted(res_map.keys())
        if explicit_scc_resolution is not None:
            target = int(explicit_scc_resolution)
            chosen = target if target in res_map else min(avail, key=lambda r: abs(r - target))
        elif int(base_resolution) < 40_000:
            ge_40k = [r for r in avail if r >= 40_000]
            chosen = ge_40k[0] if ge_40k else avail[-1]
        else:
            chosen = int(base_resolution) if int(base_resolution) in res_map else min(
                avail, key=lambda r: abs(r - int(base_resolution))
            )
        return f"{path}::{res_map[chosen]}", chosen
    return _resolve_cooler_uri(path, base_resolution)


def _has_valid_weights(clr) -> bool:
    import numpy as np

    bins_cols = list(clr.bins().columns)
    if "weight" not in bins_cols:
        return False
    w = np.asarray(clr.bins()[:]["weight"], dtype=float)
    finite = w[np.isfinite(w)]
    return bool(finite.size > 0 and float(np.max(finite)) > 0.0)


def _measure_cooler_pixels(clr) -> dict[str, Any]:
    """Exact out-of-core chunked contact statistics over ``clr.pixels()``."""
    import numpy as np

    bins_df = clr.bins()[:][["chrom", "start"]]
    bin_chrom = bins_df["chrom"].astype(str).to_numpy()
    bin_start = bins_df["start"].to_numpy(dtype=np.int64)

    nnz = int(clr.info.get("nnz", len(clr.pixels())))
    total_pairs = 0.0
    cis_pairs = 0.0
    trans_pairs = 0.0
    cis_short = 0.0
    cis_long = 0.0

    px_selector = clr.pixels()
    for lo in range(0, nnz, PIXEL_CHUNK_SIZE):
        hi = min(nnz, lo + PIXEL_CHUNK_SIZE)
        chunk = px_selector[lo:hi]
        b1 = chunk["bin1_id"].to_numpy(dtype=np.int64)
        b2 = chunk["bin2_id"].to_numpy(dtype=np.int64)
        cnt_raw = chunk["count"].to_numpy(dtype=np.float64)
        cnt = np.where(np.isfinite(cnt_raw), cnt_raw, 0.0)

        total_pairs += float(np.sum(cnt, dtype=np.float64))
        same_chrom = bin_chrom[b1] == bin_chrom[b2]
        cis_cnt = cnt[same_chrom]
        cis_pairs += float(np.sum(cis_cnt, dtype=np.float64))
        trans_pairs += float(np.sum(cnt[~same_chrom], dtype=np.float64))

        if cis_cnt.size:
            dist = np.abs(bin_start[b2[same_chrom]] - bin_start[b1[same_chrom]])
            long_mask = dist >= 20_000
            cis_long += float(np.sum(cis_cnt[long_mask], dtype=np.float64))
            cis_short += float(np.sum(cis_cnt[~long_mask], dtype=np.float64))

    denom = total_pairs if total_pairs > 0.0 else 1.0
    return {
        "n_valid_pairs": int(round(total_pairs)),
        "n_nonzero_pixels": int(nnz),
        "n_bins": int(clr.info.get("nbins", len(bins_df))),
        "cis_reads": int(round(cis_pairs)),
        "trans_reads": int(round(trans_pairs)),
        "cis_trans_ratio": round(float(cis_pairs / denom), 6),
        "cis_short_lt20kb": int(round(cis_short)),
        "cis_long_ge20kb": int(round(cis_long)),
        "long_range_cis": round(float(cis_long / denom), 6),
    }


def _stage_and_balance_cooler(
    input_path: str,
    sample_id: str,
    resolution: int,
    out_dir: Path,
    force_rebalance: bool,
    warnings: list[str],
    extra_resolutions: list[int] | None = None,
) -> tuple[str | None, int, dict[str, Any]]:
    import cooler
    import numpy as np

    is_mcool = input_path.lower().endswith(".mcool")
    res_map = list_mcool_resolutions(input_path) if is_mcool else {}
    uri, chosen_res = _resolve_cooler_uri(input_path, resolution)
    clr = cooler.Cooler(uri)

    # Check if all needed resolutions already carry ICE weights
    resolutions_to_check = [chosen_res]
    if is_mcool:
        avail = sorted(res_map.keys())
        for r in avail:
            if r >= 100_000 and r not in resolutions_to_check:
                resolutions_to_check.append(r)
        for er in extra_resolutions or []:
            if int(er) in res_map and int(er) not in resolutions_to_check:
                resolutions_to_check.append(int(er))

    if not force_rebalance:
        all_weighted = True
        for r in resolutions_to_check:
            r_uri = f"{input_path}::{res_map[r]}" if is_mcool else input_path
            if not _has_valid_weights(cooler.Cooler(r_uri)):
                all_weighted = False
                break
        if all_weighted:
            w = np.asarray(clr.bins()[:]["weight"], dtype=float)
            n_valid = int(np.isfinite(w).sum())
            return input_path, chosen_res, {
                "already_balanced": True,
                "n_balanced_bins": n_valid,
                "frac_valid_bins": round(n_valid / max(len(w), 1), 4),
            }

    ext = ".mcool" if is_mcool else ".cool"
    staged_path = out_dir / f"{sample_id}.balanced{ext}"
    if Path(input_path).resolve() != staged_path.resolve():
        shutil.copy2(input_path, staged_path)
        staged_path.chmod(0o644)

    balance_meta: dict[str, Any] = {"already_balanced": False}
    converged_by_res: dict[str, bool] = {}
    all_converged = True
    for r in resolutions_to_check:
        r_uri = f"{staged_path}::{res_map[r]}" if is_mcool else str(staged_path)
        r_clr = cooler.Cooler(r_uri)
        if not force_rebalance and _has_valid_weights(r_clr):
            converged_by_res[str(r)] = True
            continue
        nbins = int(r_clr.info.get("nbins", 0))
        min_nnz = 10 if nbins >= 500 else max(1, min(3, nbins // 10))
        ignore_diags = 2 if nbins >= 20 else 1
        try:
            weights, stats = cooler.balance_cooler(
                r_clr,
                cis_only=False,
                ignore_diags=ignore_diags,
                min_nnz=min_nnz,
                mad_max=5,
                max_iters=300,
                store=True,
            )
            r_conv = bool(stats.get("converged", True))
            converged_by_res[str(r)] = r_conv
            if not r_conv:
                all_converged = False
                warnings.append(
                    f"ICE balancing did not converge for {sample_id} at {r:,}bp."
                )
            if r == chosen_res:
                n_valid = int(np.isfinite(weights).sum())
                balance_meta.update({
                    "n_balanced_bins": n_valid,
                    "frac_valid_bins": round(n_valid / max(len(weights), 1), 4),
                })
        except Exception as exc:
            converged_by_res[str(r)] = False
            all_converged = False
            warnings.append(
                f"ICE balancing failed for {sample_id} at {r:,}bp: {str(exc)[:160]}"
            )
            if r == chosen_res:
                return None, chosen_res, {
                    "already_balanced": False,
                    "converged": False,
                    "converged_by_resolution": converged_by_res,
                }

    balance_meta["converged"] = all_converged
    balance_meta["converged_by_resolution"] = converged_by_res
    final_uri = f"{staged_path}::{res_map[chosen_res]}" if is_mcool else str(staged_path)
    final_clr = cooler.Cooler(final_uri)
    if _has_valid_weights(final_clr):
        w = np.asarray(final_clr.bins()[:]["weight"], dtype=float)
        n_valid = int(np.isfinite(w).sum())
        balance_meta.setdefault("n_balanced_bins", n_valid)
        balance_meta.setdefault("frac_valid_bins", round(n_valid / max(len(w), 1), 4))
        return str(staged_path), chosen_res, balance_meta
    warnings.append(
        f"ICE balancing produced zero valid finite weights for {sample_id} at {chosen_res:,}bp."
    )
    balance_meta["converged"] = False
    balance_meta["n_balanced_bins"] = 0
    balance_meta["frac_valid_bins"] = 0.0
    return None, chosen_res, balance_meta


def _compute_ps_decay(clr, out_tsv: Path) -> dict[str, Any]:
    import numpy as np
    import pandas as pd
    from cooltools.api import expected

    view_df = _nuclear_view_df(clr)
    if view_df.empty:
        return {"ps_curve_tsv": None, "ps_decay_slope": None}

    weight_name = "weight" if _has_valid_weights(clr) else None
    ignore_diags = 2 if int(clr.info.get("nbins", 0)) >= 20 else 1
    exp_df = expected.expected_cis(
        clr,
        view_df=view_df,
        smooth=False,
        aggregate_smoothed=False,
        clr_weight_name=weight_name,
        ignore_diags=ignore_diags,
        nproc=1,
    )
    use_bal_col = bool(weight_name and "balanced.avg" in exp_df.columns)
    sum_col = "balanced.sum" if use_bal_col else "count.sum"
    val_col = "balanced.avg" if use_bal_col else "count.avg"
    denom_col = "n_valid" if "n_valid" in exp_df.columns else "n_total"

    if sum_col in exp_df.columns and denom_col in exp_df.columns:
        work = exp_df.copy()
        s_val = pd.to_numeric(work[sum_col], errors="coerce").to_numpy(dtype=float)
        d_val = pd.to_numeric(work[denom_col], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(s_val) & np.isfinite(d_val) & (d_val > 0.0)
        work["_sum_ok"] = np.where(ok, s_val, 0.0)
        work["_den_ok"] = np.where(ok, d_val, 0.0)
        agg = (
            work.groupby("dist", as_index=False)
            .agg(
                dist_bp=("dist_bp", "first"),
                n_valid=("_den_ok", "sum"),
                _sum_val=("_sum_ok", "sum"),
            )
            .sort_values("dist")
        )
        denom = agg["n_valid"].to_numpy(dtype=float)
        num = agg["_sum_val"].to_numpy(dtype=float)
        agg["contact_prob"] = np.divide(
            num, denom, out=np.full_like(num, np.nan, dtype=float), where=denom > 0.0
        )
        agg = agg.drop(columns=["_sum_val"])
    else:
        # Fallback: n_valid-weighted average across chromosomes.
        work = exp_df.copy()
        w = pd.to_numeric(work.get(denom_col, 1.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
        v = pd.to_numeric(work[val_col], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(w) & (w > 0.0) & np.isfinite(v)
        work["_w"] = np.where(ok, w, 0.0)
        work["_wv"] = np.where(ok, w * v, 0.0)
        agg = (
            work.groupby("dist", as_index=False)
            .agg(
                dist_bp=("dist_bp", "first"),
                n_valid=("_w", "sum"),
                _wv=("_wv", "sum"),
            )
            .sort_values("dist")
        )
        denom = agg["n_valid"].to_numpy(dtype=float)
        wv_arr = agg["_wv"].to_numpy(dtype=float)
        agg["contact_prob"] = np.divide(
            wv_arr, denom, out=np.full_like(wv_arr, np.nan, dtype=float), where=denom > 0.0
        )
        agg = agg.drop(columns=["_wv"])
    agg.to_csv(out_tsv, sep="\t", index=False)

    mask = (
        (agg["dist_bp"] >= max(2 * int(clr.binsize), 20_000))
        & (agg["dist_bp"] <= 2_000_000)
        & np.isfinite(agg["contact_prob"])
        & (agg["contact_prob"] > 0)
    )
    sub = agg.loc[mask]
    if len(sub) < 3:
        mask_fb = (
            (agg["dist"] >= 2)
            & np.isfinite(agg["contact_prob"])
            & (agg["contact_prob"] > 0)
        )
        sub = agg.loc[mask_fb]

    slope = None
    if len(sub) >= 3:
        x = np.log10(sub["dist_bp"].to_numpy(dtype=float))
        y = np.log10(sub["contact_prob"].to_numpy(dtype=float))
        fin = np.isfinite(x) & np.isfinite(y)
        x, y = x[fin], y[fin]
        if x.size >= 3 and np.ptp(x) > 0:
            s_raw = float(np.polyfit(x, y, 1)[0])
            if math.isfinite(s_raw):
                slope = round(s_raw, 4)

    return {
        "ps_curve_tsv": str(out_tsv),
        "ps_decay_slope": slope,
        "n_distance_strata": int(len(agg)),
    }


def _build_box_kernel(n_dim: int, h_win: int):
    """Row-stochastic 1D box-smoothing operator K so K @ M @ K.T is symmetric and preserves constants."""
    import numpy as np
    import scipy.sparse as sp

    diags = [np.ones(n_dim - abs(d), dtype=float) for d in range(-h_win, h_win + 1)]
    k_mat = sp.diags(
        diags,
        offsets=list(range(-h_win, h_win + 1)),
        shape=(n_dim, n_dim),
        format="csr",
    )
    row_sums = np.asarray(k_mat.sum(axis=1), dtype=float).ravel()
    inv = np.where(row_sums > 0.0, 1.0 / row_sums, 0.0)
    return (sp.diags(inv) @ k_mat).tocsr()


def _box_smooth_2d(m_sym, n_dim: int, h_win: int, kernel_cache: dict[tuple[int, int], Any] | None = None):
    """Symmetric 2D mean filter ``K @ M @ K.T`` preserving symmetry and constant matrices."""
    if h_win <= 0:
        return m_sym.tocoo()
    key = (n_dim, h_win)
    if kernel_cache is not None and key in kernel_cache:
        k_op = kernel_cache[key]
    else:
        k_op = _build_box_kernel(n_dim, h_win)
        if kernel_cache is not None:
            kernel_cache[key] = k_op
    return (k_op @ m_sym @ k_op.T).tocoo()


def _compute_hicrep_scc(
    sample_uris: list[tuple[str, str]],
    max_distance_bp: int = 2_000_000,
    max_diag_cap: int = 200,
) -> tuple[dict[str, dict[str, float | None]], float | None]:
    """HiCRep stratum-adjusted correlation coefficient (Yang et al. 2017) across
    coolers sharing a bin coordinate grid, using symmetric sparse banded 2D mean
    smoothing (``K @ M @ K.T``), per-chromosome streaming (zero genome-wide OOM),
    and rank-based vstran stratum weights."""
    import cooler
    import numpy as np
    import scipy.sparse as sp
    from scipy.stats import rankdata

    if len(sample_uris) < 2:
        return {}, None

    coolers = [(sid, cooler.Cooler(uri)) for sid, uri in sample_uris]
    base_bins = coolers[0][1].bins()[:][["chrom", "start", "end"]].reset_index(drop=True)
    binsize = int(coolers[0][1].binsize)
    for _, clr in coolers[1:]:
        b_other = clr.bins()[:][["chrom", "start", "end"]].reset_index(drop=True)
        if int(clr.binsize) != binsize or len(b_other) != len(base_bins):
            return {}, None
        if not (
            np.array_equal(
                b_other["chrom"].astype(str).to_numpy(),
                base_bins["chrom"].astype(str).to_numpy(),
            )
            and np.array_equal(
                b_other["start"].to_numpy(dtype=np.int64),
                base_bins["start"].to_numpy(dtype=np.int64),
            )
        ):
            return {}, None

    max_diag = max(3, int(max_distance_bp // max(binsize, 1)))
    nbins_total = int(coolers[0][1].info.get("nbins", len(base_bins)))
    start_diag = 2 if nbins_total >= 20 else 1
    if max_diag_cap > 0:
        max_diag = min(max_diag, int(max_diag_cap))
    if max_diag < start_diag:
        max_diag = start_diag
    # Yang et al. (2017) resolution-scaled 2D mean-filter half-window h.
    if binsize >= 1_000_000 or nbins_total < 20:
        h_base = 0
    elif binsize >= 100_000:
        h_base = 1
    else:
        h_base = 2

    view_df = _nuclear_view_df(coolers[0][1])
    chroms = list(view_df["chrom"]) if not view_df.empty else list(coolers[0][1].chromnames)

    # Cache 1D sparse box-smoothing operator K by (n, h) across samples.
    kernel_cache: dict[tuple[int, int], Any] = {}
    sids = [sid for sid, _ in coolers]
    use_bal_by_sid = {sid: _has_valid_weights(clr) for sid, clr in coolers}
    pair_num: dict[tuple[int, int], float] = {
        (i, j): 0.0 for i in range(len(sids)) for j in range(i + 1, len(sids))
    }
    pair_den: dict[tuple[int, int], float] = {
        (i, j): 0.0 for i in range(len(sids)) for j in range(i + 1, len(sids))
    }

    # Stream chromosome by chromosome so dense diagonals are freed after each chromosome.
    max_elems_per_chrom_sample = 10_000_000
    for ch in chroms:
        n = int(len(coolers[0][1].bins().fetch(ch)))
        if n <= start_diag:
            continue
        k_budget = max(start_diag, max_elems_per_chrom_sample // max(n, 1))
        k_hi = min(max_diag, n - 1, k_budget)
        if k_hi < start_diag:
            continue
        h_ch = min(h_base, max(0, (n - start_diag - 1) // 4))

        ch_strata: dict[str, dict[int, np.ndarray]] = {}
        for sid, clr in coolers:
            use_bal = use_bal_by_sid[sid]
            bins_ch = clr.bins().fetch(ch)
            bad_bin = (
                ~np.isfinite(np.asarray(bins_ch["weight"], dtype=float))
                if (use_bal and "weight" in bins_ch.columns)
                else np.zeros(n, dtype=bool)
            )
            coo = clr.matrix(balance=use_bal, sparse=True).fetch(ch).tocoo()
            rows = coo.row
            cols = coo.col
            vals = np.asarray(coo.data, dtype=float)
            offs = cols - rows
            # Exclude k < start_diag before smoothing so uncalibrated ICE diagonals never bleed in.
            band = (offs >= start_diag) & (offs <= k_hi + 2 * h_ch) & np.isfinite(vals)
            rows_b = rows[band]
            cols_b = cols[band]
            vals_b = vals[band]
            if rows_b.size > 0 and np.any(bad_bin):
                ok_in = ~(bad_bin[rows_b] | bad_bin[cols_b])
                rows_b = rows_b[ok_in]
                cols_b = cols_b[ok_in]
                vals_b = vals_b[ok_in]

            if h_ch > 0 and rows_b.size > 0:
                m_up = sp.coo_matrix((vals_b, (rows_b, cols_b)), shape=(n, n)).tocsr()
                m_sym = m_up + m_up.T
                sm = _box_smooth_2d(m_sym, n, h_ch, kernel_cache=kernel_cache)
                sm_offs = sm.col - sm.row
                sm_mask = (sm_offs >= start_diag) & (sm_offs <= k_hi) & np.isfinite(sm.data)
                rows_b = sm.row[sm_mask]
                offs_b = sm_offs[sm_mask]
                vals_b = np.asarray(sm.data[sm_mask], dtype=float)
            else:
                offs_b = cols_b - rows_b

            diag_arrays: dict[int, np.ndarray] = {}
            for k in range(start_diag, k_hi + 1):
                arr = np.zeros(n - k, dtype=float)
                if np.any(bad_bin):
                    arr[bad_bin[: n - k] | bad_bin[k:]] = np.nan
                diag_arrays[k] = arr
            if rows_b.size > 0:
                for k in range(start_diag, k_hi + 1):
                    km = offs_b == k
                    if np.any(km):
                        rk = rows_b[km]
                        vk = vals_b[km]
                        if np.any(bad_bin):
                            ok_b = ~(bad_bin[rk] | bad_bin[rk + k])
                            rk = rk[ok_b]
                            vk = vk[ok_b]
                        if rk.size > 0:
                            diag_arrays[k][rk] = vk
            ch_strata[sid] = diag_arrays

        for i in range(len(sids)):
            s1 = sids[i]
            for j in range(i + 1, len(sids)):
                s2 = sids[j]
                for k in range(start_diag, k_hi + 1):
                    x = ch_strata[s1][k]
                    y = ch_strata[s2][k]
                    if x.size < 3:
                        continue
                    # Yang et al. (2017) Section 2.2: filter both-zero pairs per stratum
                    # before variance-stabilizing rank weighting and Pearson correlation.
                    valid = np.isfinite(x) & np.isfinite(y) & ((x > 0.0) | (y > 0.0))
                    n_k = int(valid.sum())
                    if n_k < 3:
                        continue
                    xv = x[valid]
                    yv = y[valid]
                    if float(np.std(xv)) <= 1e-12 or float(np.std(yv)) <= 1e-12:
                        continue
                    r_k = float(np.corrcoef(xv, yv)[0, 1])
                    if not math.isfinite(r_k):
                        continue
                    rx = rankdata(xv, method="average") / float(n_k)
                    ry = rankdata(yv, method="average") / float(n_k)
                    sx = float(np.std(rx))
                    sy = float(np.std(ry))
                    if sx <= 1e-12 or sy <= 1e-12:
                        continue
                    w_k = float(n_k) * sx * sy
                    pair_num[(i, j)] += w_k * r_k
                    pair_den[(i, j)] += w_k
        kernel_cache.clear()

    matrix: dict[str, dict[str, float | None]] = {s: {} for s in sids}
    off_diag_vals: list[float] = []
    for i, s1 in enumerate(sids):
        matrix[s1][s1] = 1.0
        for j in range(i + 1, len(sids)):
            s2 = sids[j]
            num = pair_num[(i, j)]
            den = pair_den[(i, j)]
            if den > 0.0 and math.isfinite(num / den):
                scc: float | None = round(float(num / den), 4)
                off_diag_vals.append(scc)
            else:
                scc = None
            matrix[s1][s2] = scc
            matrix[s2][s1] = scc

    mean_scc = (
        round(float(sum(off_diag_vals) / len(off_diag_vals)), 4)
        if off_diag_vals else None
    )
    return matrix, mean_scc


def hic_qc_and_balance(params: dict) -> dict:
    import numpy as np

    files = [str(f) for f in (params.get("files") or []) if f]
    genome = str(params.get("genome") or "hg38")
    resolution = int(params.get("resolution") or 40_000)
    force_rebalance = bool(params.get("force_rebalance", False))
    allow_mock = mocks_allowed(params)
    warnings: list[str] = []

    if not files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": "No valid Hi-C files found.",
        }
    missing_files = [f for f in files if not Path(f).exists()]
    if missing_files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": f"Declared Hi-C input file(s) not found: {missing_files}",
        }
    valid_files = list(files)

    out_dir = Path(
        str(params.get("output_dir") or (Path(valid_files[0]).parent / "aria_hic_qc"))
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    declared_sids = params.get("sample_ids")
    if declared_sids is not None:
        if not isinstance(declared_sids, list) or len(declared_sids) != len(valid_files):
            return {
                "status": "error",
                "error_type": "SampleIdLengthMismatch",
                "details": (
                    f"Declared sample_ids length "
                    f"({len(declared_sids) if isinstance(declared_sids, list) else type(declared_sids).__name__}) "
                    f"does not match files length ({len(valid_files)})."
                ),
            }
        raw_sids = [str(s) for s in declared_sids]
    else:
        raw_sids = [Path(f).stem.replace(".mcool", "").replace(".cool", "") for f in valid_files]
    sid_counts: dict[str, int] = {}
    for s in raw_sids:
        sid_counts[s] = sid_counts.get(s, 0) + 1
    seen_sids: set[str] = set()
    sample_ids: list[str] = []
    for idx_s, (s, fpath_s) in enumerate(zip(raw_sids, valid_files)):
        cand = s
        if sid_counts.get(s, 0) > 1:
            parent_tag = Path(fpath_s).parent.name.strip()
            cand = f"{parent_tag}_{s}" if parent_tag else f"{s}_{idx_s + 1}"
            if cand in seen_sids:
                cand = f"{s}_{idx_s + 1}"
        while cand in seen_sids:
            cand = f"{cand}_{idx_s + 1}"
        seen_sids.add(cand)
        sample_ids.append(cand)

    min_cis_trans = float(params.get("min_cis_trans_ratio", 0.60))
    min_valid_pairs = int(params.get("min_valid_pairs", 1_000_000))
    chrom_sizes_param = params.get("chrom_sizes") or params.get("chromsizes")
    genome_fasta_param = params.get("genome_fasta") or params.get("reference_fasta")
    genome_fai_param = params.get("genome_fai")

    raw_input_files = list(valid_files)
    for ref_val in (chrom_sizes_param, genome_fasta_param, genome_fai_param):
        if isinstance(ref_val, (str, Path)) and str(ref_val) and Path(str(ref_val)).is_file():
            raw_input_files.append(str(ref_val))
    explicit_scc_res = (
        int(params["scc_resolution"])
        if params.get("scc_resolution") is not None
        else None
    )
    manifest_params = {
        "genome": genome,
        "resolution": resolution,
        "scc_resolution": explicit_scc_res,
        "force_rebalance": force_rebalance,
        "sample_ids": sample_ids,
        "min_cis_trans_ratio": min_cis_trans,
        "min_valid_pairs": min_valid_pairs,
        "chrom_sizes": (
            chrom_sizes_param
            if isinstance(chrom_sizes_param, dict)
            else str(chrom_sizes_param or "")
        ),
        "genome_fasta": str(genome_fasta_param or ""),
        "genome_fai": str(genome_fai_param or ""),
    }
    if bool(params.get("resume", True)):
        cached = _check_hic_stage_cache(
            out_dir,
            stage="hic_qc_and_balance",
            version=STAGE_VERSION,
            params=manifest_params,
            inputs=raw_input_files,
        )
        if cached is not None and isinstance(cached.get("result"), dict):
            resumed = dict(cached["result"])
            resumed["resumed"] = True
            return resumed

    try:
        import cooler  # noqa: F401
    except ImportError:
        if allow_mock:
            return _mock_hic_qc_bundle(valid_files, resolution)
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": "cooler and cooltools are required in aria-hic-env.",
        }

    # If raw .hic or .pairs files arrived directly without hic_preprocess, convert first
    needs_prep = any(
        Path(f).name.lower().endswith((".hic", ".pairs", ".pairs.gz", ".bam", ".sam"))
        for f in valid_files
    )
    if needs_prep:
        from aria.scripts.hic_preprocess import hic_preprocess

        prep_base_res = (
            resolution if resolution <= 25_000 else math.gcd(resolution, 10_000)
        )
        if prep_base_res <= 0:
            prep_base_res = min(resolution, 10_000)
        prep_res = hic_preprocess({
            "files": valid_files,
            "sample_ids": sample_ids,
            "genome": genome,
            "resolution": prep_base_res,
            "zoom_resolutions": sorted({prep_base_res, resolution, 100_000}),
            "chrom_sizes": chrom_sizes_param,
            "genome_fasta": genome_fasta_param,
            "genome_fai": genome_fai_param,
            "output_dir": str(out_dir / "preprocessed"),
        })
        if prep_res.get("status") != "success":
            return prep_res
        valid_files = list(prep_res["cooler_files"])

    import cooler

    results_per_file: list[dict[str, Any]] = []
    balanced_files: list[str] = []
    balanced_sample_ids: list[str] = []
    sample_uris_for_scc: list[tuple[str, str]] = []
    scc_res_used: int | None = None

    for sid, fpath in zip(sample_ids, valid_files, strict=True):
        try:
            extra_res: list[int] = []
            if fpath.lower().endswith(".mcool"):
                _, target_scc_res = _resolve_scc_uri(fpath, resolution, explicit_scc_res)
                extra_res.append(target_scc_res)
            balanced_path, chosen_res, bal_meta = _stage_and_balance_cooler(
                fpath,
                sid,
                resolution,
                out_dir,
                force_rebalance,
                warnings,
                extra_resolutions=extra_res,
            )
            active_path = balanced_path or fpath
            uri, actual_res = _resolve_cooler_uri(active_path, chosen_res)
            clr = cooler.Cooler(uri)
            px_stats = _measure_cooler_pixels(clr)
            ps_info = _compute_ps_decay(clr, out_dir / f"{sid}.expected_cis.tsv")

            rec = {
                "sample_id": sid,
                "file": fpath,
                "format": "mcool" if fpath.lower().endswith(".mcool") else "cool",
                "resolution": actual_res,
                **px_stats,
                "balanced_path": balanced_path,
                "balance_info": bal_meta,
                **ps_info,
                "status": "success",
            }
            results_per_file.append(rec)
            if balanced_path:
                balanced_files.append(balanced_path)
                balanced_sample_ids.append(sid)
                scc_uri, scc_res_cur = _resolve_scc_uri(
                    balanced_path, actual_res, explicit_scc_res
                )
                if scc_res_used is None:
                    scc_res_used = scc_res_cur
                sample_uris_for_scc.append((sid, scc_uri))
        except (OSError, ValueError, RuntimeError, KeyError) as exc:
            results_per_file.append({
                "sample_id": sid,
                "file": fpath,
                "status": "error",
                "error": str(exc)[:300],
            })

    failed_records = [r for r in results_per_file if r.get("status") != "success"]
    if failed_records or not results_per_file:
        return {
            "status": "error",
            "error_type": "ProcessingFailed",
            "details": (
                f"Could not QC/balance {len(failed_records)} of {len(valid_files)} "
                f"Hi-C file(s): {failed_records or results_per_file}"
            ),
            "per_file": results_per_file,
            "warnings": warnings,
        }
    ok_records = results_per_file

    all_valid_pairs = int(sum(r.get("n_valid_pairs", 0) for r in ok_records))
    avg_cis_trans = float(np.mean([r["cis_trans_ratio"] for r in ok_records]))
    avg_long_cis = float(np.mean([r.get("long_range_cis", 0.0) for r in ok_records]))
    slopes = [r["ps_decay_slope"] for r in ok_records if r.get("ps_decay_slope") is not None]
    mean_ps_slope = round(float(np.mean(slopes)), 4) if slopes else None

    scc_matrix, mean_scc = _compute_hicrep_scc(sample_uris_for_scc)
    scc_tsv: str | None = None
    if scc_matrix:
        import pandas as pd

        scc_df = pd.DataFrame(scc_matrix)
        scc_path = out_dir / "hicrep_scc_matrix.tsv"
        scc_df.to_csv(scc_path, sep="\t")
        scc_tsv = str(scc_path)

    pass_qc = True
    if len(balanced_files) < len(ok_records):
        warnings.append(
            f"ICE balancing failed or produced zero valid weights for "
            f"{len(ok_records) - len(balanced_files)} of {len(ok_records)} sample(s)."
        )
        pass_qc = False
    if any(
        isinstance(r.get("balance_info"), dict) and r["balance_info"].get("converged") is False
        for r in ok_records
    ):
        pass_qc = False
    if avg_cis_trans < min_cis_trans:
        warnings.append(
            f"Cis/trans ratio {avg_cis_trans:.2%} < {min_cis_trans:.0%}. "
            "Low cis fraction suggests elevated random ligation or background."
        )
        pass_qc = False
    if all_valid_pairs < min_valid_pairs:
        warnings.append(
            f"Total valid contact pairs ({all_valid_pairs:,}) is below the "
            f"configured floor ({min_valid_pairs:,}) for {resolution:,}bp analysis."
        )
        pass_qc = False

    primary_ps = next(
        (r.get("ps_curve_tsv") for r in ok_records if r.get("ps_curve_tsv")),
        None,
    )
    result = {
        "status": "success",
        "ran": True,
        "data_type": "HiC",
        "analysis": "hic_qc_and_balance",
        "validation_level": "production",
        "genome": genome,
        "resolution": int(ok_records[0]["resolution"]),
        "n_samples": len(ok_records),
        "sample_ids": [r["sample_id"] for r in ok_records],
        "n_valid_pairs": all_valid_pairs,
        "cis_trans_ratio": round(avg_cis_trans, 4),
        "long_range_cis": round(avg_long_cis, 4),
        "ps_decay_slope": mean_ps_slope,
        "ps_curve_tsv": primary_ps,
        "scc_matrix": scc_matrix,
        "mean_replicate_scc": mean_scc,
        "scc_resolution": scc_res_used,
        "scc_tsv": scc_tsv,
        "balanced_files": balanced_files,
        "balanced_sample_ids": balanced_sample_ids,
        "per_file": results_per_file,
        "thresholds": {
            "min_cis_trans_ratio": min_cis_trans,
            "min_valid_pairs": min_valid_pairs,
        },
        "pass_qc": bool(pass_qc),
        "warnings": warnings,
        "output_dir": str(out_dir),
    }

    all_ps_tsvs = [
        r.get("ps_curve_tsv") for r in ok_records if r.get("ps_curve_tsv")
    ]
    manifest_outputs = list(
        dict.fromkeys(p for p in [*balanced_files, *all_ps_tsvs, scc_tsv] if p)
    )
    _save_hic_stage_cache(
        out_dir,
        stage="hic_qc_and_balance",
        version=STAGE_VERSION,
        params=manifest_params,
        inputs=raw_input_files,
        outputs=manifest_outputs,
        extra={"result": result},
    )
    return result


def _mock_hic_qc_bundle(files: list[str], resolution: int) -> dict:
    return {
        "status": "success",
        "ran": True,
        "data_type": "HiC",
        "analysis": "hic_qc_and_balance",
        "validation_level": "scaffold",
        "resolution": resolution,
        "n_valid_pairs": 500_000_000,
        "cis_trans_ratio": 0.72,
        "long_range_cis": 0.45,
        "balanced_files": [],
        "balanced_sample_ids": [],
        "per_file": [
            {
                "file": f,
                "format": "mock",
                "resolution": resolution,
                "n_valid_pairs": 500_000_000,
                "cis_trans_ratio": 0.72,
                "status": "success",
            }
            for f in files
        ],
        "pass_qc": True,
        "warnings": ["Mock QC — aria-hic-env not active"],
    }


if __name__ == "__main__":
    run_script(hic_qc_and_balance)
