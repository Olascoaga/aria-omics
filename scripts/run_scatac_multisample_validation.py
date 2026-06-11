#!/usr/bin/env python
"""
ARIA scATAC multi-sample pseudobulk DA validation harness (v4.6, ADR-041 follow-up)
-----------------------------------------------------------------------------------
Exercises ARIA's REAL chromatin pseudobulk differential-accessibility lane
(`aria.scripts.chromatin_diffacc`, shared `rna_bulk_de._run_deseq2` core) on a
genuine multi-donor, multi-condition scATAC peak matrix — the run HC11
(single-sample) could only honestly skip and that ADR-041 only validated against
synthetic truth.

Input: Samael's Erosion CONSENSUS peak universe (`results/02_consensus/`), which
already solves the cross-donor peak-alignment problem. The raw per-sample
`.h5mu` files cannot be concatenated (Cell Ranger ARC calls peaks per sample →
disjoint peak sets), but the consensus groups share ONE peak set
(`peaks.tsv` identical md5 across groups) and the barcodes retain donor identity.

For one cell type (default Oligo) the consensus has 8 groups = 4 age_group x 2
sex, each a common 142,228-peak matrix, 10 distinct donors per age_group.

Pipeline (all REAL ARIA code). It MUST run in an env that has BOTH the scanpy
LSI/Wilcoxon stack AND pydeseq2 — that is `aria-rna-env` (scanpy 1.11 + pydeseq2
0.5.4), NOT `aria-chromatin-env` (which lacks pydeseq2, so its pseudobulk lane
can only error out — the v1 run used that env and would have wasted the whole
Wilcoxon before failing the DESeq2 step it set out to validate).
  1. build a combined cells x peaks AnnData from the 8 consensus MEX groups
     (obs: age_group=condition, donor=biological replicate, sex, cell_type),
     then SUBSAMPLE to a tractable validation size: a stratified ~8k-cell draw
     by (age_group, sex, donor) with a per-donor floor so the 10-vs-10 donor
     structure of the contrast survives, plus a peak filter (drop near-empty
     peaks, cap to the most-accessible). A validation exercises the code path,
     not the full 96k x 142k matrix;
  2. chromatin_lsi_clustering at a low leiden resolution (default 0.4 -> tens of
     clusters, not ~143) -> structural groupby (LSI assigns no biology, ADR-011);
  3. chromatin_diffacc with condition_col=age_group, replicate_col=donor and an
     explicit contrast (default 80-100 vs 20-39): per-cluster (descriptive) +
     pseudobulk DA (inferential, ~10 vs ~10 donor replicates) via the shared
     validated DESeq2 core, parallelized over n_cpus.

This is a dataset-gated VALIDATION harness, not a product feature and not a
promotion of scATAC beyond alpha. It writes only to a workspace dir.

Usage (RNA env — has scanpy AND pydeseq2):
    conda run -n aria-rna-env python scripts/run_scatac_multisample_validation.py

Env overrides:
    ARIA_SCATAC_CONSENSUS_DIR   (default /home/medusa/Samael/Erosion/results/02_consensus)
    ARIA_SCATAC_VAL_OUT         (default ~/.aria/workspace/scatac_multisample_validation)
    ARIA_SCATAC_CELL_TYPE       (default Oligo)
    ARIA_SCATAC_TEST_AGE        (default 80-100)
    ARIA_SCATAC_REF_AGE         (default 20-39)
    ARIA_SCATAC_MAX_CELLS       (default 8000)  stratified subsample target
    ARIA_SCATAC_FLOOR_PER_DONOR (default 25)    min cells kept per donor stratum
    ARIA_SCATAC_MIN_CELLS_PEAK  (default 10)    drop peaks open in fewer cells
    ARIA_SCATAC_MAX_PEAKS       (default 25000) cap to most-accessible peaks
    ARIA_SCATAC_LEIDEN_RES      (default 0.4)   leiden resolution
    ARIA_SCATAC_N_CPUS          (default os.cpu_count())  DESeq2 + BLAS threads
    ARIA_SCATAC_SEED            (default 0)      subsample RNG seed
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

# Multi-core: cap every thread pool BEFORE numpy/anndata/scanpy import so the
# LSI SVD, scanpy ops and pydeseq2's numpy land use all requested cores. Set as
# defaults (setdefault) so an explicit outer env still wins.
_N_CPUS = int(os.environ.get("ARIA_SCATAC_N_CPUS") or (os.cpu_count() or 1))
for _tv in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_tv, str(_N_CPUS))

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

CONSENSUS_DIR = Path(
    os.environ.get("ARIA_SCATAC_CONSENSUS_DIR",
                   "/home/medusa/Samael/Erosion/results/02_consensus")
)
OUT_DIR = Path(
    os.environ.get("ARIA_SCATAC_VAL_OUT",
                   str(Path.home() / ".aria/workspace/scatac_multisample_validation"))
)
CELL_TYPE = os.environ.get("ARIA_SCATAC_CELL_TYPE", "Oligo")
AGES = ["20-39", "40-59", "60-79", "80-100"]
SEXES = ["F", "M"]
TEST_AGE = os.environ.get("ARIA_SCATAC_TEST_AGE", "80-100")
REF_AGE = os.environ.get("ARIA_SCATAC_REF_AGE", "20-39")
# Coarse Leiden resolution on purpose: this is a structural groupby for the
# descriptive per-cluster lane only (the validation target is the pseudobulk
# DA lane). A high resolution over-clusters the multi-donor pool into 100+ tiny
# clusters, which makes the per-cluster pct_in/pct_out CSR column slicing
# pathologically slow without adding meaning. 0.4 -> tens of clusters.
LEIDEN_RES = float(os.environ.get("ARIA_SCATAC_LEIDEN_RES", "0.4"))

# Validation-scale knobs. The full pool is ~96k cells x 142k peaks; a validation
# only needs to exercise the real DA code path, so we draw a stratified subsample
# and filter near-empty peaks (see _subsample_and_filter). Neutral filters only.
MAX_CELLS = int(os.environ.get("ARIA_SCATAC_MAX_CELLS", "8000"))
FLOOR_PER_DONOR = int(os.environ.get("ARIA_SCATAC_FLOOR_PER_DONOR", "25"))
MIN_CELLS_PEAK = int(os.environ.get("ARIA_SCATAC_MIN_CELLS_PEAK", "10"))
MAX_PEAKS = int(os.environ.get("ARIA_SCATAC_MAX_PEAKS", "25000"))
SEED = int(os.environ.get("ARIA_SCATAC_SEED", "0"))

# Barcode = "<donor>_<16bp cell barcode>-<digit>"; donor may itself contain "_"
# (e.g. "hc29_deep"), so strip the trailing _<ACGTN...>-<n> instead of split("_").
_BC_RE = re.compile(r"_[ACGTN]+-\d+$")


def _donor_of(barcode: str) -> str:
    return _BC_RE.sub("", barcode)


def _read_lines(p: Path) -> list[str]:
    with open(p) as f:
        return [ln.strip() for ln in f if ln.strip()]


def build_combined_adata():
    import anndata as ad
    import pandas as pd
    from scipy.io import mmread

    groups = []
    for age in AGES:
        for sex in SEXES:
            gdir = CONSENSUS_DIR / f"{CELL_TYPE}_{age}_{sex}"
            if not gdir.is_dir():
                print(f"  [skip] missing group dir: {gdir}", flush=True)
                continue
            mtx = gdir / "matrix.mtx"
            barcodes = _read_lines(gdir / "barcodes.tsv")
            peaks = _read_lines(gdir / "peaks.tsv")
            # MEX is peaks x cells -> transpose to AnnData cells x peaks.
            X = mmread(str(mtx)).T.tocsr()
            assert X.shape == (len(barcodes), len(peaks)), (
                f"{gdir.name}: matrix {X.shape} != "
                f"(cells {len(barcodes)}, peaks {len(peaks)})")
            donors = [_donor_of(bc) for bc in barcodes]
            obs = pd.DataFrame({
                "barcode": barcodes,
                "donor": donors,
                "age_group": age,
                "sex": sex,
                "cell_type": CELL_TYPE,
            })
            obs.index = [f"{age}_{sex}_{bc}" for bc in barcodes]
            a = ad.AnnData(X=X, obs=obs)
            a.var_names = peaks
            groups.append(a)
            print(f"  [load] {gdir.name}: {a.n_obs} cells x {a.n_vars} peaks, "
                  f"{len(set(donors))} donors", flush=True)

    if not groups:
        raise SystemExit(f"No consensus groups found under {CONSENSUS_DIR}")

    # All groups share the identical peak set -> inner join is loss-free.
    combined = ad.concat(groups, join="inner", index_unique=None)
    combined.var_names = groups[0].var_names  # concat preserves; be explicit
    combined.obs["age_group"] = combined.obs["age_group"].astype("category")
    combined.obs["donor"] = combined.obs["donor"].astype("category")
    return combined


def _subsample_and_filter(adata, report):
    """Stratified cell subsample + peak filter -> tractable validation matrix.

    Cells: a ~MAX_CELLS draw stratified by (age_group, sex, donor), each stratum
    floored at FLOOR_PER_DONOR so every donor keeps enough cells for a stable
    pseudobulk aggregate (the 10-vs-10 contrast must survive). Peaks: drop those
    open in < MIN_CELLS_PEAK cells, then cap to the MAX_PEAKS most-accessible.
    Both are biologically NEUTRAL filters (composition / accessibility only, no
    effect-size peeking), so they don't bias the DA they feed.
    """
    import numpy as np
    import scipy.sparse as sp

    n0_cells, n0_peaks = int(adata.n_obs), int(adata.n_vars)

    if adata.n_obs > MAX_CELLS:
        rng = np.random.default_rng(SEED)
        obs = adata.obs
        strata = (obs["age_group"].astype(str) + "|" + obs["sex"].astype(str)
                  + "|" + obs["donor"].astype(str))
        total = adata.n_obs
        keep = []
        for _, pos in obs.groupby(strata, observed=True).indices.items():
            n_avail = len(pos)
            quota = max(FLOOR_PER_DONOR, round(MAX_CELLS * n_avail / total))
            keep.append(rng.choice(pos, size=min(n_avail, quota), replace=False))
        keep = np.sort(np.concatenate(keep))
        adata = adata[keep, :].copy()

    # Peak filter on the subsampled cells (a peak is "open" in a cell when > 0).
    Xc = adata.X.tocsc() if sp.issparse(adata.X) else sp.csc_matrix(adata.X)
    cells_per_peak = np.asarray((Xc > 0).sum(axis=0)).ravel()
    keep_mask = cells_per_peak >= MIN_CELLS_PEAK
    if MAX_PEAKS and int(keep_mask.sum()) > MAX_PEAKS:
        idx = np.where(keep_mask)[0]
        top = idx[np.argsort(cells_per_peak[idx])[::-1][:MAX_PEAKS]]
        keep_mask = np.zeros(adata.n_vars, dtype=bool)
        keep_mask[top] = True
    adata = adata[:, keep_mask].copy()

    # Drop unused donor/age_group/sex categories so downstream design matrices
    # and groupbys see only levels actually present after subsetting.
    for col in ("age_group", "donor", "sex"):
        if col in adata.obs and str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()

    report["subsample"] = {
        "seed": SEED,
        "cells_before": n0_cells, "cells_after": int(adata.n_obs),
        "max_cells_target": MAX_CELLS, "floor_per_donor": FLOOR_PER_DONOR,
        "peaks_before": n0_peaks, "peaks_after": int(adata.n_vars),
        "min_cells_per_peak": MIN_CELLS_PEAK, "max_peaks_cap": MAX_PEAKS,
    }
    print(f"      subsampled: {adata.n_obs} cells x {adata.n_vars} peaks "
          f"(from {n0_cells} x {n0_peaks})", flush=True)
    return adata


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    report: dict = {"cell_type": CELL_TYPE, "contrast": f"{TEST_AGE}_vs_{REF_AGE}"}

    # ── Step 1: build (or reuse) full combined AnnData, then subsample ──────
    import anndata as ad
    combined_full_path = OUT_DIR / f"{CELL_TYPE.lower()}_consensus_combined.h5ad"
    if combined_full_path.exists():
        print(f"[1/3] reuse full combined AnnData: {combined_full_path}", flush=True)
        combined = ad.read_h5ad(str(combined_full_path))
    else:
        print("[1/3] building combined multi-donor AnnData ...", flush=True)
        combined = build_combined_adata()
        combined.write_h5ad(str(combined_full_path))

    # Subsample + peak-filter to the validation matrix, then persist it; LSI and
    # DA both run on this so they share one tractable, reproducible input.
    combined = _subsample_and_filter(combined, report)
    combined_path = OUT_DIR / f"{CELL_TYPE.lower()}_consensus_sub.h5ad"
    combined.write_h5ad(str(combined_path))

    by_age = combined.obs.groupby("age_group", observed=True)["donor"].nunique()
    report["combined"] = {
        "n_cells": int(combined.n_obs),
        "n_peaks": int(combined.n_vars),
        "donors_per_age_group": {k: int(v) for k, v in by_age.items()},
        "n_donors_total": int(combined.obs["donor"].nunique()),
        "path": str(combined_path),
    }
    print(f"      combined: {combined.n_obs} cells x {combined.n_vars} peaks; "
          f"donors/age = {report['combined']['donors_per_age_group']}", flush=True)

    # ── Step 2: REAL LSI/clustering for a structural groupby ────────────────
    print("[2/3] chromatin_lsi_clustering (real ARIA script) ...", flush=True)
    from aria.scripts.chromatin_lsi_clustering import chromatin_lsi_clustering
    lsi = chromatin_lsi_clustering({
        "data_path": str(combined_path),
        "output_dir": str(OUT_DIR),
        "resolution": LEIDEN_RES,
    })
    report["lsi"] = {k: lsi.get(k) for k in (
        "status", "n_cells_used", "n_peaks", "n_components_used",
        "dropped_components", "n_clusters", "output_path", "warnings")}
    if lsi.get("status") != "success":
        report["lsi"]["details"] = lsi.get("details")
        _save(report)
        raise SystemExit(f"LSI failed: {lsi.get('details')}")
    clustered_path = lsi["output_path"]
    print(f"      LSI ok: {lsi['n_clusters']} clusters, dropped "
          f"{lsi['dropped_components']}", flush=True)

    # ── Step 3: REAL differential accessibility (per-cluster + pseudobulk) ──
    print("[3/3] chromatin_diffacc (real ARIA script, pseudobulk lane) ...",
          flush=True)
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    da = chromatin_diffacc({
        "data_path": clustered_path,
        "groupby": "leiden",
        "condition_col": "age_group",
        "replicate_col": "donor",
        "comparisons": [{"test": TEST_AGE, "reference": REF_AGE}],
        "output_dir": str(OUT_DIR),
        "n_cpus": _N_CPUS,
    })
    report["diffacc"] = {
        "status": da.get("status"),
        "n_clusters": da.get("n_clusters"),
        "padj_max": da.get("padj_max"),
        "lfc_min": da.get("lfc_min"),
        "per_cluster": {k: da.get("per_cluster", {}).get(k) for k in (
            "ran", "reason", "n_da_total", "output_csv")},
        "pseudobulk": da.get("pseudobulk"),
        "warnings": (da.get("warnings") or [])[:10],
    }
    pb = da.get("pseudobulk", {})
    print(f"      per-cluster ran={da.get('per_cluster', {}).get('ran')} "
          f"n_da={da.get('per_cluster', {}).get('n_da_total')}", flush=True)
    print(f"      pseudobulk ran={pb.get('ran')} "
          f"n_pseudosamples={pb.get('n_pseudosamples')} "
          f"n_peaks_tested={pb.get('n_peaks_tested')}", flush=True)
    for c in pb.get("comparisons", []):
        print(f"        {c.get('test')} vs {c.get('reference')}: "
              f"status={c.get('status')} n_sig={c.get('n_sig')} "
              f"n_up={c.get('n_up')} n_down={c.get('n_down')} "
              f"n_rep={c.get('n_replicates')} "
              f"low_power={c.get('low_power_warning')}", flush=True)

    _save(report)
    print(f"\nValidation report -> {OUT_DIR / 'validation_report.json'}",
          flush=True)


def _save(report: dict):
    out = OUT_DIR / "validation_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)


if __name__ == "__main__":
    main()
