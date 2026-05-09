"""
ARIA RNA Clustering Script
--------------------------
Runs normalization, dimensionality reduction, and Leiden clustering
on a QC-filtered AnnData object.
Executed inside aria-rna-env by EnvironmentManager.

Input params:
    data_path:    str   — path to QC-filtered .h5ad
    resolution:   float — Leiden resolution (from ParameterAdvisor)
    n_hvg:        int   (optional) — highly variable genes (default: 3000)
    n_pcs:        int   (optional) — PCA components (default: 40)
    n_neighbors:  int   (optional) — kNN neighbors (default: 15)

Output:
    {
      "status":       "success",
      "n_clusters":   int,
      "cluster_sizes": {cluster_id: n_cells},
      "top_markers":  {cluster_id: [gene1, gene2, ...]},
      "output_path":  str — path to clustered .h5ad
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def rna_clustering(params: dict) -> dict:
    import scanpy as sc
    from pathlib import Path

    data_path   = params["data_path"]
    resolution  = float(params["resolution"])
    n_hvg       = int(params.get("n_hvg", 3000))
    n_pcs       = int(params.get("n_pcs", 40))
    n_neighbors = int(params.get("n_neighbors", 15))

    adata = sc.read_h5ad(data_path)

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
    sc.pp.scale(adata, max_value=10)

    # Dimensionality reduction
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(
        adata, resolution=resolution,
        flavor="igraph", n_iterations=2, directed=False,
    )

    # Marker genes per cluster
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

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

    output_path = str(Path(data_path).parent / "clustered.h5ad")
    adata.write_h5ad(output_path)

    return {
        "status":        "success",
        "n_clusters":    n_clusters,
        "cluster_sizes": cluster_sizes,
        "top_markers":   top_markers,
        "output_path":   output_path,
        "resolution":    resolution,
    }


if __name__ == "__main__":
    run_script(rna_clustering)
