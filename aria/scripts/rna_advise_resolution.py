"""
ARIA Leiden Resolution Advisor Script
--------------------------------------
Evaluates N Leiden resolutions on a preprocessed .h5ad and returns objective
metrics so the agent can pick the best resolution before user approval at CP3.

Executed inside aria-rna-env by EnvironmentManager. The agent no longer
imports scanpy in-process for this — it just orchestrates.

Input params:
    data_path:     str          — path to clustered or preprocessed .h5ad
                                   (must have neighbors; PCA computed on it)
    resolutions:   list[float]  — resolutions to evaluate (default: [0.2,0.5,0.8,1.2])
    rep:           str  (opt)   — embedding key for silhouette (default: auto)
                                   prefers X_pca_harmony, falls back to X_pca
    max_cells:     int  (opt)   — sketch size for objective scoring

Output:
    {
      "status":      "success" | "error",
      "rep_used":    str,
      "candidates":  [
          {
            "resolution":           float,
            "n_clusters":           int,
            "silhouette":           float,
            "n_singleton_clusters": int,
            "min_cluster_size":     int,
            "max_cluster_size":     int,
          },
          ...
      ],
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def rna_advise_resolution(params: dict) -> dict:
    import scanpy as sc
    import numpy as np
    import json
    import hashlib
    from pathlib import Path
    from sklearn.metrics import silhouette_score

    data_path   = params["data_path"]
    resolutions = params.get("resolutions") or [0.2, 0.5, 0.8, 1.2]
    rep_pref    = params.get("rep")
    seed        = int(params.get("seed", 0))
    n_hvg       = int(params.get("n_hvg", 3000))
    n_pcs       = int(params.get("n_pcs", 40))
    max_cells   = int(params.get("max_cells", 100_000))

    data_file = Path(data_path)
    cache_key = hashlib.sha1(json.dumps({
        "data_path": str(data_file.resolve()),
        "resolutions": [float(r) for r in resolutions],
        "rep": rep_pref,
        "seed": seed,
        "n_hvg": n_hvg,
        "n_pcs": n_pcs,
        "max_cells": max_cells,
    }, sort_keys=True).encode()).hexdigest()[:12]
    cache_path = data_file.parent / f"leiden_advice_{cache_key}.json"
    if (cache_path.exists() and data_file.exists() and
            cache_path.stat().st_mtime >= data_file.stat().st_mtime):
        try:
            with open(cache_path) as f:
                cached = json.load(f)
            cached["resumed"] = True
            cached["warnings"] = (
                cached.get("warnings", []) +
                ["[resume] Leiden resolution advice loaded from cache"]
            )
            return cached
        except Exception:
            pass

    adata = sc.read_h5ad(data_path)
    n_cells_total = int(adata.n_obs)
    sketch_used = False
    if adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        if "batch" in adata.obs.columns:
            batches = adata.obs["batch"].astype(str)
            keep = []
            per_batch = max(1, max_cells // max(1, int(batches.nunique())))
            for batch in sorted(batches.unique()):
                idx = np.flatnonzero((batches == batch).to_numpy())
                if len(idx) > per_batch:
                    idx = rng.choice(idx, size=per_batch, replace=False)
                keep.extend(idx.tolist())
            if len(keep) > max_cells:
                keep = rng.choice(np.array(keep), size=max_cells,
                                  replace=False).tolist()
            keep = np.array(sorted(keep), dtype=int)
        else:
            keep = np.sort(rng.choice(adata.n_obs, size=max_cells,
                                      replace=False))
        adata = adata[keep, :].copy()
        sketch_used = True

    # If PCA is missing (e.g. caller is at qc_filtered.h5ad stage),
    # do the preprocessing in-memory. Mirrors rna_clustering.py so the
    # downstream clustering step will skip these stages (idempotency).
    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        if adata.raw is None:
            adata.raw = adata.copy()
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs,
                  random_state=seed)

    # Pick embedding: harmony-corrected if available, else PCA.
    if rep_pref and rep_pref in adata.obsm:
        rep = rep_pref
    elif "X_pca_harmony" in adata.obsm:
        rep = "X_pca_harmony"
    else:
        rep = "X_pca"

    # Ensure neighbor graph exists so leiden can run on every candidate.
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, n_pcs=n_pcs,
                        random_state=seed)

    pca_emb = adata.obsm[rep][:, :20] if adata.obsm[rep].shape[1] >= 20 \
              else adata.obsm[rep]

    candidates = []
    for res in resolutions:
        try:
            res_val = float(res)
        except (TypeError, ValueError):
            continue
        # Use a scratch key so we don't overwrite an existing "leiden" column.
        sc.tl.leiden(
            adata, resolution=res_val,
            key_added="_advise_leiden",
            flavor="igraph", n_iterations=2, directed=False,
            random_state=0,
        )
        labels       = adata.obs["_advise_leiden"].astype(str).values
        unique       = list(set(labels))
        n_clusters   = len(unique)
        cluster_sizes = [int((labels == u).sum()) for u in unique]
        n_singletons = int(sum(1 for s in cluster_sizes if s < 10))

        sil = 0.0
        try:
            if 1 < n_clusters < len(labels):
                sil = float(silhouette_score(
                    pca_emb, labels, sample_size=min(2000, len(labels)),
                    random_state=0,
                ))
        except Exception:
            pass

        candidates.append({
            "resolution":           res_val,
            "n_clusters":           n_clusters,
            "silhouette":           round(sil, 4),
            "n_singleton_clusters": n_singletons,
            "min_cluster_size":     int(min(cluster_sizes)),
            "max_cluster_size":     int(max(cluster_sizes)),
        })

    # Don't leave the scratch column on the AnnData on disk; we never write
    # it back, but other consumers reading from the same path shouldn't see
    # half-evaluated state.
    if "_advise_leiden" in adata.obs.columns:
        adata.obs.drop(columns="_advise_leiden", inplace=True)

    result = {
        "status":        "success",
        "rep_used":      rep,
        "candidates":    candidates,
        "sketch_used":   sketch_used,
        "n_cells_total": n_cells_total,
        "n_cells_used":  int(adata.n_obs),
    }
    try:
        with open(cache_path, "w") as f:
            json.dump(result, f, indent=2)
    except Exception:
        pass
    return result


if __name__ == "__main__":
    run_script(rna_advise_resolution)
