"""
ARIA RNA Integration Script
----------------------------
Harmony batch correction for multi-sample / multi-batch scRNA-seq.
Executed inside aria-rna-env by EnvironmentManager (standalone entry point).

Input params:
    data_path:    str   — path to preprocessed .h5ad (PCA already computed)
    batch_col:    str   — obs column with batch/sample labels (default: "batch")
    output_dir:   str   (optional) — where to write integrated.h5ad

Output:
    {
      "status":              "success" | "skipped" | "error",
      "method":              "harmony",
      "n_batches":           int,
      "silhouette_before":   float,   # batch separation before (lower = better)
      "silhouette_after":    float,   # batch separation after
      "batch_correction_delta": float,
      "rep_used":            str,     # embedding key in obsm
      "output_path":         str,
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def rna_integration(params: dict) -> dict:
    import scanpy as sc
    import numpy as np
    from pathlib import Path

    data_path  = params["data_path"]
    batch_col  = params.get("batch_col", "batch")
    output_dir = params.get("output_dir", str(Path(data_path).parent))

    adata = sc.read_h5ad(data_path)

    if batch_col not in adata.obs.columns:
        return {"status": "skipped",
                "reason":  f"batch column '{batch_col}' not found in obs"}

    n_batches = int(adata.obs[batch_col].nunique())
    if n_batches < 2:
        return {"status": "skipped",
                "reason":  "only one batch — no integration needed"}

    # Ensure PCA is present
    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    # Silhouette score — measures batch separation (lower is better after correction)
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    batch_labels = le.fit_transform(adata.obs[batch_col].astype(str))
    pca20 = adata.obsm["X_pca"][:, :20]
    sil_before = float(silhouette_score(pca20, batch_labels))

    # Harmony integration
    try:
        sc.external.pp.harmony_integrate(
            adata, batch_col,
            basis="X_pca",
            adjusted_basis="X_pca_harmony",
        )
        rep = "X_pca_harmony"
    except Exception:
        # Direct harmonypy fallback
        import harmonypy as hm
        ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_col)
        adata.obsm["X_pca_harmony"] = ho.Z_corr.T
        rep = "X_pca_harmony"

    sil_after = float(silhouette_score(adata.obsm[rep][:, :20], batch_labels))

    # Recompute graph + UMAP on corrected embedding
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)

    output_path = str(Path(output_dir) / "integrated.h5ad")
    adata.write_h5ad(output_path)

    return {
        "status":                 "success",
        "method":                 "harmony",
        "n_batches":              n_batches,
        "batch_col":              batch_col,
        "silhouette_before":      round(sil_before, 4),
        "silhouette_after":       round(sil_after, 4),
        "batch_correction_delta": round(sil_before - sil_after, 4),
        "rep_used":               rep,
        "output_path":            output_path,
    }


if __name__ == "__main__":
    run_script(rna_integration)
