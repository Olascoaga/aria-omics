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
        avail = sorted(
            int(p.split("/")[-1])
            for p in cooler.fileops.list_coolers(path)
            if p.split("/")[-1].isdigit()
        )
        if not avail:
            raise ValueError(f"No resolutions inside .mcool: {path}")
        chosen = int(resolution) if int(resolution) in avail else min(
            avail, key=lambda r: abs(r - int(resolution))
        )
        return f"{path}::resolutions/{chosen}", chosen
    clr = cooler.Cooler(path)
    return path, int(clr.binsize)


def _has_valid_weights(clr) -> bool:
    import numpy as np

    bins_cols = list(clr.bins().columns)
    if "weight" not in bins_cols:
        return False
    w = np.asarray(clr.bins()[:]["weight"], dtype=float)
    return bool(np.isfinite(w).any() and np.nanmax(w) > 0)


def _measure_cooler_pixels(clr) -> dict[str, Any]:
    """Exact out-of-core chunked contact statistics over ``clr.pixels()``."""
    import numpy as np

    bins_df = clr.bins()[:][["chrom", "start"]]
    bin_chrom = bins_df["chrom"].astype(str).to_numpy()
    bin_start = bins_df["start"].to_numpy(dtype=np.int64)

    nnz = int(clr.info.get("nnz", len(clr.pixels())))
    total_pairs = 0
    cis_pairs = 0
    trans_pairs = 0
    cis_short = 0
    cis_long = 0

    px_selector = clr.pixels()
    for lo in range(0, nnz, PIXEL_CHUNK_SIZE):
        hi = min(nnz, lo + PIXEL_CHUNK_SIZE)
        chunk = px_selector[lo:hi]
        b1 = chunk["bin1_id"].to_numpy(dtype=np.int64)
        b2 = chunk["bin2_id"].to_numpy(dtype=np.int64)
        cnt = chunk["count"].to_numpy(dtype=np.int64)

        total_pairs += int(cnt.sum())
        same_chrom = bin_chrom[b1] == bin_chrom[b2]
        cis_cnt = cnt[same_chrom]
        cis_pairs += int(cis_cnt.sum())
        trans_pairs += int(cnt[~same_chrom].sum())

        if cis_cnt.size:
            dist = np.abs(bin_start[b2[same_chrom]] - bin_start[b1[same_chrom]])
            long_mask = dist >= 20_000
            cis_long += int(cis_cnt[long_mask].sum())
            cis_short += int(cis_cnt[~long_mask].sum())

    denom = max(total_pairs, 1)
    return {
        "n_valid_pairs": int(total_pairs),
        "n_nonzero_pixels": int(nnz),
        "n_bins": int(clr.info.get("nbins", len(bins_df))),
        "cis_reads": int(cis_pairs),
        "trans_reads": int(trans_pairs),
        "cis_trans_ratio": round(float(cis_pairs / denom), 6),
        "cis_short_lt20kb": int(cis_short),
        "cis_long_ge20kb": int(cis_long),
        "long_range_cis": round(float(cis_long / denom), 6),
    }


def _stage_and_balance_cooler(
    input_path: str,
    sample_id: str,
    resolution: int,
    out_dir: Path,
    force_rebalance: bool,
    warnings: list[str],
) -> tuple[str | None, int, dict[str, Any]]:
    import cooler
    import numpy as np

    uri, chosen_res = _resolve_cooler_uri(input_path, resolution)
    clr = cooler.Cooler(uri)

    # Check if all needed resolutions already carry ICE weights
    resolutions_to_check = [chosen_res]
    if input_path.lower().endswith(".mcool"):
        avail = sorted(
            int(p.split("/")[-1])
            for p in cooler.fileops.list_coolers(input_path)
            if p.split("/")[-1].isdigit()
        )
        for r in avail:
            if r >= 100_000 and r not in resolutions_to_check:
                resolutions_to_check.append(r)

    if not force_rebalance:
        all_weighted = True
        for r in resolutions_to_check:
            r_uri = f"{input_path}::resolutions/{r}" if input_path.lower().endswith(".mcool") else input_path
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

    ext = ".mcool" if input_path.lower().endswith(".mcool") else ".cool"
    staged_path = out_dir / f"{sample_id}.balanced{ext}"
    if Path(input_path).resolve() != staged_path.resolve():
        shutil.copy2(input_path, staged_path)
        staged_path.chmod(0o644)

    balance_meta: dict[str, Any] = {"already_balanced": False}
    for r in resolutions_to_check:
        r_uri = f"{staged_path}::resolutions/{r}" if ext == ".mcool" else str(staged_path)
        r_clr = cooler.Cooler(r_uri)
        if not force_rebalance and _has_valid_weights(r_clr):
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
            if r == chosen_res:
                n_valid = int(np.isfinite(weights).sum())
                balance_meta.update({
                    "converged": bool(stats.get("converged", True)),
                    "n_balanced_bins": n_valid,
                    "frac_valid_bins": round(n_valid / max(len(weights), 1), 4),
                })
        except Exception as exc:
            warnings.append(
                f"ICE balancing failed for {sample_id} at {r:,}bp: {str(exc)[:160]}"
            )
            if r == chosen_res:
                return None, chosen_res, {"already_balanced": False, "converged": False}

    final_uri = f"{staged_path}::resolutions/{chosen_res}" if ext == ".mcool" else str(staged_path)
    final_clr = cooler.Cooler(final_uri)
    if _has_valid_weights(final_clr):
        w = np.asarray(final_clr.bins()[:]["weight"], dtype=float)
        n_valid = int(np.isfinite(w).sum())
        balance_meta.setdefault("n_balanced_bins", n_valid)
        balance_meta.setdefault("frac_valid_bins", round(n_valid / max(len(w), 1), 4))
        return str(staged_path), chosen_res, balance_meta
    return None, chosen_res, balance_meta


def _compute_ps_decay(clr, out_tsv: Path) -> dict[str, Any]:
    import numpy as np
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
    val_col = "balanced.avg" if weight_name and "balanced.avg" in exp_df.columns else "count.avg"
    agg = (
        exp_df.groupby("dist", as_index=False)
        .agg(
            dist_bp=("dist_bp", "first"),
            n_valid=("n_valid", "sum"),
            contact_prob=(val_col, "mean"),
        )
        .sort_values("dist")
    )
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
        if np.ptp(x) > 0:
            slope = round(float(np.polyfit(x, y, 1)[0]), 4)

    return {
        "ps_curve_tsv": str(out_tsv),
        "ps_decay_slope": slope,
        "n_distance_strata": int(len(agg)),
    }


def _compute_hicrep_scc(
    sample_uris: list[tuple[str, str]],
    max_distance_bp: int = 2_000_000,
) -> tuple[dict[str, dict[str, float]], float | None]:
    """HiCRep stratum-adjusted correlation coefficient (Yang et al. 2017) across
    balanced coolers sharing a bin coordinate grid."""
    import cooler
    import numpy as np

    if len(sample_uris) < 2:
        return {}, None

    coolers = [(sid, cooler.Cooler(uri)) for sid, uri in sample_uris]
    base_bins = coolers[0][1].bins()[:][["chrom", "start", "end"]].reset_index(drop=True)
    binsize = int(coolers[0][1].binsize)
    for _, clr in coolers[1:]:
        b_other = clr.bins()[:][["chrom", "start", "end"]].reset_index(drop=True)
        if int(clr.binsize) != binsize or len(b_other) != len(base_bins):
            return {}, None
        if not (b_other["chrom"].astype(str).equals(base_bins["chrom"].astype(str))
                and np.array_equal(b_other["start"].to_numpy(), base_bins["start"].to_numpy())):
            return {}, None

    max_diag = max(3, int(max_distance_bp // max(binsize, 1)))
    view_df = _nuclear_view_df(coolers[0][1])
    chroms = list(view_df["chrom"]) if not view_df.empty else list(coolers[0][1].chromnames)

    # Extract cis strata vectors per sample: dict[diag_k -> list of 1D arrays]
    strata_by_sample: dict[str, dict[int, np.ndarray]] = {}
    for sid, clr in coolers:
        use_bal = _has_valid_weights(clr)
        per_diag: dict[int, list[np.ndarray]] = {k: [] for k in range(1, max_diag + 1)}
        for ch in chroms:
            mat = clr.matrix(balance=use_bal).fetch(ch)
            n = mat.shape[0]
            if n < 3:
                continue
            for k in range(1, min(max_diag + 1, n)):
                d = np.asarray(np.diag(mat, k=k), dtype=float)
                per_diag[k].append(d)
        strata_by_sample[sid] = {
            k: np.concatenate(arrs) if arrs else np.empty(0, dtype=float)
            for k, arrs in per_diag.items()
        }

    sids = [sid for sid, _ in coolers]
    matrix: dict[str, dict[str, float]] = {s: {} for s in sids}
    off_diag_vals: list[float] = []

    for i, s1 in enumerate(sids):
        matrix[s1][s1] = 1.0
        for j in range(i + 1, len(sids)):
            s2 = sids[j]
            num = 0.0
            den = 0.0
            for k in range(1, max_diag + 1):
                x = strata_by_sample[s1][k]
                y = strata_by_sample[s2][k]
                if x.size < 3:
                    continue
                valid = np.isfinite(x) & np.isfinite(y)
                n_k = int(valid.sum())
                if n_k < 3:
                    continue
                xv = x[valid]
                yv = y[valid]
                # Variance-stabilizing weight N_k * std(r_x) * std(r_y) / N_total_approx
                sx = float(np.std(xv))
                sy = float(np.std(yv))
                if sx <= 1e-12 or sy <= 1e-12:
                    continue
                r_k = float(np.corrcoef(xv, yv)[0, 1])
                if not math.isfinite(r_k):
                    continue
                w_k = float(n_k) * sx * sy
                num += w_k * r_k
                den += w_k
            scc = round(float(num / den), 4) if den > 0 else 0.0
            matrix[s1][s2] = scc
            matrix[s2][s1] = scc
            off_diag_vals.append(scc)

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

    valid_files = [f for f in files if Path(f).exists()]
    if not valid_files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": "No valid Hi-C files found.",
        }

    out_dir = Path(
        str(params.get("output_dir") or (Path(valid_files[0]).parent / "aria_hic_qc"))
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    declared_sids = params.get("sample_ids")
    if isinstance(declared_sids, list) and len(declared_sids) == len(valid_files):
        sample_ids = [str(s) for s in declared_sids]
    else:
        sample_ids = [Path(f).stem.replace(".mcool", "").replace(".cool", "") for f in valid_files]

    manifest_params = {
        "genome": genome,
        "resolution": resolution,
        "force_rebalance": force_rebalance,
        "sample_ids": sample_ids,
    }
    if bool(params.get("resume", True)):
        cached = _check_hic_stage_cache(
            out_dir,
            stage="hic_qc_and_balance",
            version=STAGE_VERSION,
            params=manifest_params,
            inputs=valid_files,
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

        prep_res = hic_preprocess({
            "files": valid_files,
            "sample_ids": sample_ids,
            "genome": genome,
            "resolution": min(resolution, 10_000),
            "zoom_resolutions": sorted({min(resolution, 10_000), resolution, 100_000}),
            "chrom_sizes": params.get("chrom_sizes"),
            "genome_fasta": params.get("genome_fasta"),
            "genome_fai": params.get("genome_fai"),
            "output_dir": str(out_dir / "preprocessed"),
        })
        if prep_res.get("status") != "success":
            return prep_res
        valid_files = list(prep_res["cooler_files"])

    import cooler

    results_per_file: list[dict[str, Any]] = []
    balanced_files: list[str] = []
    sample_uris_for_scc: list[tuple[str, str]] = []

    for sid, fpath in zip(sample_ids, valid_files):
        try:
            balanced_path, chosen_res, bal_meta = _stage_and_balance_cooler(
                fpath, sid, resolution, out_dir, force_rebalance, warnings
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
                sample_uris_for_scc.append((sid, uri))
        except Exception as exc:
            results_per_file.append({
                "sample_id": sid,
                "file": fpath,
                "status": "error",
                "error": str(exc)[:300],
            })

    ok_records = [r for r in results_per_file if r.get("status") == "success"]
    if not ok_records:
        return {
            "status": "error",
            "error_type": "ProcessingFailed",
            "details": f"Could not QC/balance any Hi-C file: {results_per_file}",
        }

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

    min_cis_trans = float(params.get("min_cis_trans_ratio", 0.60))
    min_valid_pairs = int(params.get("min_valid_pairs", 1_000_000))

    pass_qc = True
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
        "scc_tsv": scc_tsv,
        "balanced_files": balanced_files,
        "per_file": results_per_file,
        "thresholds": {
            "min_cis_trans_ratio": min_cis_trans,
            "min_valid_pairs": min_valid_pairs,
        },
        "pass_qc": bool(pass_qc),
        "warnings": warnings,
        "output_dir": str(out_dir),
    }

    manifest_outputs = [p for p in [*balanced_files, primary_ps, scc_tsv] if p]
    _save_hic_stage_cache(
        out_dir,
        stage="hic_qc_and_balance",
        version=STAGE_VERSION,
        params=manifest_params,
        inputs=valid_files,
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
