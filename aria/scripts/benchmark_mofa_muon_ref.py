"""V48 Slice I4 — External reference comparator (`muon.tl.mofa`).

Executes `muon`'s native `muon.tl.mofa` workflow (`scanpy` HVG log-normalization +
`muon.atac.pp.tfidf` + `muon.tl.mofa`) on a paired `.h5mu` file and exports
aligned cell barcodes, factor scores `X_mofa`, per-view variance explained, and
feature loadings `LFs` for external concordance scoring against ARIA's
`aria/scripts/integration_mofa.py`.
"""

from __future__ import annotations

import contextlib
import io
import math
import os
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def benchmark_mofa_muon_ref(params: dict[str, Any]) -> dict[str, Any]:
    import anndata as ad
    import h5py
    import mudata as md
    import muon as mu
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    from aria.utils.mudata_io import read_h5mu_paired

    h5mu_path = str(params["h5mu_path"])
    output_dir = Path(str(params["output_dir"]))
    n_factors = int(params.get("n_factors", 10))
    n_top_genes = int(params.get("n_top_genes", 3000))
    n_top_peaks = int(params.get("n_top_peaks", 5000))
    scale_views = bool(params.get("scale_views", True))
    seed = int(params.get("seed", 0))
    output_dir.mkdir(parents=True, exist_ok=True)

    paired = read_h5mu_paired(h5mu_path, require_raw_counts=True)
    if paired.get("status") != "success":
        return {
            "status": "error",
            "error_type": paired.get("error_type", "H5muReadFailed"),
            "details": paired.get("details", "Failed to load paired .h5mu"),
        }

    barcodes = list(paired["shared_barcodes"])
    src_rna = paired["rna_adata"]
    src_atac = paired["atac_adata"]
    rna_mod = str(paired["rna_modality"])
    atac_mod = str(paired["atac_modality"])

    r_map = {str(bc): i for i, bc in enumerate(src_rna.obs_names)}
    a_map = {str(bc): i for i, bc in enumerate(src_atac.obs_names)}
    r_take = [r_map[bc] for bc in barcodes]
    a_take = [a_map[bc] for bc in barcodes]

    rna_ad = ad.AnnData(
        X=sp.csr_matrix(paired["rna_counts_matrix"][r_take, :]).astype(np.float64),
        obs=src_rna.obs.iloc[r_take].copy(),
        var=src_rna.var.copy(),
    )
    atac_ad = ad.AnnData(
        X=sp.csr_matrix(paired["atac_counts_matrix"][a_take, :]).astype(np.float64),
        obs=src_atac.obs.iloc[a_take].copy(),
        var=src_atac.var.copy(),
    )

    # Standard Scanpy + Muon preprocessing before muon.tl.mofa
    sc.pp.normalize_total(rna_ad, target_sum=1e4)
    sc.pp.log1p(rna_ad)
    n_g_pick = min(n_top_genes, rna_ad.n_vars)
    if n_g_pick < rna_ad.n_vars:
        sc.pp.highly_variable_genes(
            rna_ad, n_top_genes=n_g_pick, flavor="seurat", subset=True
        )
    rna_ad.X = (
        rna_ad.X.toarray() if sp.issparse(rna_ad.X) else np.asarray(rna_ad.X)
    ).astype(np.float64)

    n_cells = atac_ad.n_obs
    peak_nnz = np.asarray((atac_ad.X > 0).sum(axis=0)).ravel()
    min_pk = max(3, int(math.ceil(0.005 * n_cells)))
    max_pk = max(min_pk + 1, int(math.floor(0.99 * n_cells)))
    keep_pk = np.flatnonzero((peak_nnz >= min_pk) & (peak_nnz <= max_pk))
    if keep_pk.size < min(20, atac_ad.n_vars):
        keep_pk = np.flatnonzero(peak_nnz >= 1)
    atac_ad = atac_ad[:, keep_pk].copy()

    mu.atac.pp.tfidf(atac_ad, scale_factor=1e4)
    X_atac_dense = (
        atac_ad.X.toarray() if sp.issparse(atac_ad.X) else np.asarray(atac_ad.X)
    ).astype(np.float64)
    n_p_pick = min(n_top_peaks, atac_ad.n_vars)
    if n_p_pick < atac_ad.n_vars:
        p_vars = np.var(X_atac_dense, axis=0)
        top_p = np.sort(np.argsort(p_vars)[::-1][:n_p_pick])
        atac_ad = atac_ad[:, top_p].copy()
        X_atac_dense = X_atac_dense[:, top_p]
    atac_ad.X = X_atac_dense

    mdata = md.MuData({rna_mod: rna_ad, atac_mod: atac_ad})
    model_hdf5 = output_dir / "muon_ref_model.hdf5"

    # Note: muon.tl.mofa prints its banner inside entry_point().__init__()
    with contextlib.redirect_stdout(io.StringIO()):
        mu.tl.mofa(
            mdata,
            use_var=None,
            n_factors=n_factors,
            scale_views=scale_views,
            scale_groups=False,
            center_groups=True,
            ard_factors=False,
            ard_weights=True,
            spikeslab_weights=True,
            spikeslab_factors=False,
            convergence_mode="fast",
            seed=seed,
            outfile=str(model_hdf5),
            save_data=False,
            save_parameters=False,
            quiet=True,
            verbose=False,
        )

    scores = np.asarray(mdata.obsm["X_mofa"], dtype=np.float64)
    scores_csv = output_dir / "muon_ref_scores.csv"
    pd.DataFrame(
        scores,
        index=barcodes,
        columns=[f"Factor{i + 1}" for i in range(scores.shape[1])],
    ).to_csv(scores_csv)

    r2_per_view: dict[str, list[float]] = {}
    r2_total: dict[str, float] = {}
    loadings_npz = output_dir / "muon_ref_loadings.npz"

    with h5py.File(model_hdf5, "r") as hf:
        views = [
            v.decode("utf-8") if isinstance(v, bytes) else str(v)
            for v in hf["views"]["views"][:]
        ]
        groups = [
            g.decode("utf-8") if isinstance(g, bytes) else str(g)
            for g in hf["groups"]["groups"][:]
        ]
        g0 = groups[0]
        r2_pf_mat = np.asarray(
            hf["variance_explained"]["r2_per_factor"][g0][:], dtype=np.float64
        )
        # Note: mofapy2 HDF5 stores r2_per_factor and r2_total in percentage (0..100)!
        # Verify scale via r2_total vs 1.0: mofapy2's save() multiplies R2 by 100 in HDF5!
        r2_tot_arr = np.asarray(
            hf["variance_explained"]["r2_total"][g0][:], dtype=np.float64
        )
        scale_div = 100.0 if (np.max(r2_tot_arr) > 1.0 or np.max(r2_pf_mat) > 1.0) else 1.0
        for vi, vname in enumerate(views):
            r2_per_view[vname] = [
                round(float(x / scale_div), 6) for x in r2_pf_mat[vi, :]
            ]
            r2_total[vname] = round(float(r2_tot_arr[vi] / scale_div), 6)

        w_rna = np.asarray(hf["expectations"]["W"][rna_mod][:], dtype=np.float64).T
        w_atac = np.asarray(hf["expectations"]["W"][atac_mod][:], dtype=np.float64).T

    np.savez_compressed(
        loadings_npz,
        rna_features=np.asarray([str(v) for v in rna_ad.var_names], dtype=object),
        rna_weights=w_rna,
        atac_features=np.asarray([str(v) for v in atac_ad.var_names], dtype=object),
        atac_weights=w_atac,
    )

    return {
        "status": "success",
        "n_cells": int(len(barcodes)),
        "n_factors": int(scores.shape[1]),
        "rna_modality": rna_mod,
        "atac_modality": atac_mod,
        "r2_total": r2_total,
        "r2_per_view": r2_per_view,
        "scores_csv": str(scores_csv),
        "loadings_npz": str(loadings_npz),
        "model_hdf5": str(model_hdf5),
        "r2_hdf5_scale_divisor": scale_div,
    }


if __name__ == "__main__":
    run_script(benchmark_mofa_muon_ref)
