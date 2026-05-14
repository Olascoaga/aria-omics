"""
ARIA RNA Clustering Script
--------------------------
Normalization, dimensionality reduction, and Leiden clustering on a
QC-filtered (or integrated) AnnData object.

Idempotent: if the input already has X_pca / X_pca_harmony / neighbors /
X_umap, those stages are skipped and the existing embeddings are reused.
This lets the agent run integration first, then clustering on the
corrected embedding, without re-doing PCA.

Reproducible: random_state pinned everywhere it matters.

Preserves raw counts in .raw before HVG subsetting so downstream DE,
pathway, and cell-comm scripts can access all genes.

Executed inside aria-rna-env by EnvironmentManager.

Input params:
    data_path:    str   — path to QC-filtered or integrated .h5ad
    resolution:   float — Leiden resolution (from ParameterAdvisor / CP3)
    n_hvg:        int   (optional) — highly variable genes (default: 3000)
    n_pcs:        int   (optional) — PCA components (default: 40)
    n_neighbors:  int   (optional) — kNN neighbors (default: 15)
    max_cells:    int   (optional) — sketch size for very large datasets
    seed:         int   (optional) — random_state (default: 0)

Output:
    {
      "status":        "success",
      "n_clusters":    int,
      "cluster_sizes": {cluster_id: n_cells},
      "top_markers":   {cluster_id: [gene1, gene2, ...]},
      "rep_used":      str   — "X_pca_harmony" if present, else "X_pca",
      "output_path":   str   — path to clustered .h5ad,
      "resolution":    float,
      "preprocessing_skipped": bool,
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _cache_params(params: dict) -> dict:
    return {
        "resolution": float(params["resolution"]),
        "n_hvg": int(params.get("n_hvg", 3000)),
        "n_pcs": int(params.get("n_pcs", 40)),
        "n_neighbors": int(params.get("n_neighbors", 15)),
        "max_cells": int(params.get("max_cells", 100_000)),
        "seed": int(params.get("seed", 0)),
    }


def _cache_matches(cached: dict, expected: dict) -> bool:
    return cached.get("cache_version") == 3 and cached.get("cache_params") == expected


def rna_clustering(params: dict) -> dict:
    import scanpy as sc
    import numpy as np
    import json
    from pathlib import Path

    data_path   = params["data_path"]
    resolution  = float(params["resolution"])
    n_hvg       = int(params.get("n_hvg", 3000))
    n_pcs       = int(params.get("n_pcs", 40))
    n_neighbors = int(params.get("n_neighbors", 15))
    max_cells   = int(params.get("max_cells", 100_000))
    seed        = int(params.get("seed", 0))
    cache_params = _cache_params(params)

    input_path = Path(data_path)
    output_name = "clustered_sketch.h5ad"
    full_output = input_path.parent / "clustered.h5ad"
    sketch_output = input_path.parent / output_name
    for existing in (sketch_output, full_output):
        summary_path = existing.with_suffix(".summary.json")
        if existing.exists() and existing.stat().st_mtime >= input_path.stat().st_mtime:
            try:
                if summary_path.exists():
                    with open(summary_path) as f:
                        cached = json.load(f)
                    if not _cache_matches(cached, cache_params):
                        continue
                    cached["resumed"] = True
                    cached["warnings"] = (
                        cached.get("warnings", []) +
                        [f"[resume] clustering skipped "
                         f"(valid {existing.name} exists)"]
                    )
                    return cached
                # Legacy cached h5ads lack the exact parameter signature; rerun
                # them instead of risking a stale resolution/sketch result.
                continue
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

    preprocessing_skipped = "X_pca" in adata.obsm

    if not preprocessing_skipped:
        # Preserve all genes in .raw BEFORE HVG subsetting so downstream
        # DE / pathway / cell-comm can still address non-HVG genes (housekeeping
        # ligand-receptor pairs, MT, ribo, etc.).
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        if adata.raw is None:
            adata.raw = adata.copy()
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs,
                  random_state=seed)

    # Prefer harmony-corrected embedding if integration ran first.
    rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"

    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                        use_rep=rep, random_state=seed)

    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata, random_state=seed)

    sc.tl.leiden(
        adata, resolution=resolution,
        flavor="igraph", n_iterations=2, directed=False,
        random_state=seed,
    )

    # Marker genes per cluster — use raw counts if available so Wilcoxon
    # operates on the full feature set, not just HVGs.
    use_raw = adata.raw is not None
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon", use_raw=use_raw,
    )

    n_clusters    = int(adata.obs["leiden"].nunique())
    cluster_sizes = {
        str(c): int((adata.obs["leiden"] == c).sum())
        for c in adata.obs["leiden"].unique()
    }
    top_markers = {}
    for cluster in adata.obs["leiden"].unique():
        try:
            genes = list(adata.uns["rank_genes_groups"]["names"][str(cluster)][:10])
            top_markers[str(cluster)] = [g for g in genes if g and str(g) != "nan"]
        except Exception:
            top_markers[str(cluster)] = []

    output_name = "clustered_sketch.h5ad" if sketch_used else "clustered.h5ad"
    output_path = Path(data_path).parent / output_name
    adata.uns["aria_n_cells_total"] = n_cells_total
    adata.write_h5ad(str(output_path))

    result = {
        "status":                "success",
        "n_clusters":            n_clusters,
        "cluster_sizes":         cluster_sizes,
        "top_markers":           top_markers,
        "rep_used":              rep,
        "output_path":           str(output_path),
        "resolution":            resolution,
        "preprocessing_skipped": preprocessing_skipped,
        "sketch_used":           sketch_used,
        "n_cells_total":         n_cells_total,
        "n_cells_used":          int(adata.n_obs),
        "cache_version":         3,
        "cache_params":          cache_params,
    }
    try:
        with open(output_path.with_suffix(".summary.json"), "w") as f:
            json.dump(result, f, indent=2)
    except Exception:
        pass
    return result


if __name__ == "__main__":
    run_script(rna_clustering)
