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


def rna_clustering(params: dict) -> dict:
    import scanpy as sc
    from pathlib import Path

    data_path   = params["data_path"]
    resolution  = float(params["resolution"])
    n_hvg       = int(params.get("n_hvg", 3000))
    n_pcs       = int(params.get("n_pcs", 40))
    n_neighbors = int(params.get("n_neighbors", 15))
    seed        = int(params.get("seed", 0))

    adata = sc.read_h5ad(data_path)

    preprocessing_skipped = "X_pca" in adata.obsm

    if not preprocessing_skipped:
        # Preserve all genes in .raw BEFORE HVG subsetting so downstream
        # DE / pathway / cell-comm can still address non-HVG genes (housekeeping
        # ligand-receptor pairs, MT, ribo, etc.).
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        if adata.raw is None:
            adata.raw = adata
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

    output_path = str(Path(data_path).parent / "clustered.h5ad")
    adata.write_h5ad(output_path)

    return {
        "status":                "success",
        "n_clusters":            n_clusters,
        "cluster_sizes":         cluster_sizes,
        "top_markers":           top_markers,
        "rep_used":              rep,
        "output_path":           output_path,
        "resolution":            resolution,
        "preprocessing_skipped": preprocessing_skipped,
    }


if __name__ == "__main__":
    run_script(rna_clustering)
