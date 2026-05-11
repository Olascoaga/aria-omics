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
    from sklearn.metrics import silhouette_score

    data_path   = params["data_path"]
    resolutions = params.get("resolutions") or [0.2, 0.5, 0.8, 1.2]
    rep_pref    = params.get("rep")

    adata = sc.read_h5ad(data_path)

    # Pick embedding: harmony-corrected if available, else PCA, else fail.
    if rep_pref and rep_pref in adata.obsm:
        rep = rep_pref
    elif "X_pca_harmony" in adata.obsm:
        rep = "X_pca_harmony"
    elif "X_pca" in adata.obsm:
        rep = "X_pca"
    else:
        return {
            "status":     "error",
            "error_type": "NoEmbedding",
            "details":    "Neither X_pca_harmony nor X_pca found in obsm. "
                          "Run clustering or integration first.",
        }

    # Ensure neighbor graph exists so leiden can run on every candidate.
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, n_pcs=30)

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

    return {
        "status":     "success",
        "rep_used":   rep,
        "candidates": candidates,
    }


if __name__ == "__main__":
    run_script(rna_advise_resolution)
