"""
ARIA RNA Per-Cluster Differential Expression Script
----------------------------------------------------
Wilcoxon rank-sum DE per Leiden cluster (cluster vs rest) on a clustered
.h5ad. Respects user-confirmed thresholds (padj, |log2FC|) from CP3.

Executed inside aria-rna-env by EnvironmentManager.

Input params:
    data_path:    str    — path to clustered .h5ad (must contain "leiden" obs)
    groupby:      str    (optional) — obs column to group by (default: "leiden")
    padj_max:     float  (optional) — adjusted p-value threshold (default: 0.05)
    lfc_min:      float  (optional) — minimum |log2FC| (default: 0.5)
    min_pct_in:   float  (optional) — minimum fraction of cells in-cluster
                                      detecting the marker (default: 0.10)
    top_n:        int    (optional) — top hits per cluster to return (default: 20)
    output_dir:   str    (optional)

Output:
    {
      "status":               "success" | "error",
      "groupby":              str,
      "marker_source":        str,
      "padj_max":             float,
      "lfc_min":              float,
      "n_significant_total":  int,
      "n_clusters":           int,
      "de_genes_by_cluster":  {cluster_id: [{gene, log2fc, padj, pct_in, pct_out}, ...]},
      "n_sig_by_cluster":     {cluster_id: int},
      "output_csv":           str — path to per-cluster DE table
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def rna_de_per_cluster(params: dict) -> dict:
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from aria.utils.safe_h5ad import read_h5ad

    data_path  = params["data_path"]
    groupby    = params.get("groupby", "leiden")
    padj_max   = float(params.get("padj_max", 0.05))
    lfc_min    = float(params.get("lfc_min", 0.5))
    min_pct_in = float(params.get("min_pct_in", 0.10))
    top_n      = int(params.get("top_n", 20))
    output_dir = params.get("output_dir", str(Path(data_path).parent))

    adata = read_h5ad(data_path)

    if groupby not in adata.obs.columns:
        return {
            "status":     "error",
            "error_type": "NoGroupBy",
            "details":    f"Column '{groupby}' not found in obs. "
                          f"Available: {list(adata.obs.columns)[:15]}",
        }

    n_clusters = int(adata.obs[groupby].nunique())
    if n_clusters < 2:
        return {
            "status": "skipped",
            "reason": f"Only {n_clusters} cluster(s) — DE requires ≥2.",
        }

    cluster_ids = {str(c) for c in adata.obs[groupby].unique()}

    def _compatible_rank_genes_groups(candidate: dict | None) -> bool:
        if not candidate:
            return False
        required = ("names", "logfoldchanges", "pvals_adj")
        if any(k not in candidate for k in required):
            return False
        field_names = getattr(getattr(candidate["names"], "dtype", None), "names", None)
        return field_names is not None and {str(c) for c in field_names} == cluster_ids

    # Wilcoxon marker ranking — cluster vs rest. Reuse clustering's compatible
    # rank_genes_groups table when present; otherwise compute it here. The .raw
    # matrix produced by rna_clustering is log-normalized full-feature data, not
    # raw counts, and lets marker ranking address non-HVG genes.
    use_raw = adata.raw is not None
    if _compatible_rank_genes_groups(adata.uns.get("rank_genes_groups")):
        rgg = adata.uns["rank_genes_groups"]
        marker_source = "rank_genes_groups"
    else:
        sc.tl.rank_genes_groups(
            adata, groupby=groupby, method="wilcoxon",
            use_raw=use_raw, key_added="aria_de",
        )
        rgg = adata.uns["aria_de"]
        marker_source = "computed_aria_de"

    # Compute pct_in / pct_out on the same matrix backing marker lookup. This
    # guards against marker calls that look strong by LFC alone but are detected
    # in <10% of cluster cells (likely noise).
    X = adata.raw.X if use_raw else adata.X
    var_names = list(adata.raw.var_names) if use_raw else list(adata.var_names)
    var_idx   = {g: i for i, g in enumerate(var_names)}

    cluster_membership = {
        str(c): np.asarray(adata.obs[groupby].astype(str) == str(c))
        for c in adata.obs[groupby].unique()
    }

    def _pct(gene_idx: int, mask: np.ndarray) -> float:
        col = X[mask, gene_idx]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        return float((np.asarray(col) > 0).mean())

    de_genes_by_cluster: dict = {}
    rows_for_csv = []
    n_sig_by_cluster: dict = {}
    n_significant_total = 0

    for cluster_id in rgg["names"].dtype.names:
        cl_str = str(cluster_id)
        mask_in  = cluster_membership.get(cl_str)
        mask_out = ~mask_in if mask_in is not None else None

        names = rgg["names"][cluster_id]
        lfc   = rgg["logfoldchanges"][cluster_id]
        padj  = rgg["pvals_adj"][cluster_id]

        sig = []
        for gene, l, p in zip(names, lfc, padj):
            gene = str(gene)
            try:
                p_val = float(p)
                l_val = float(l)
            except (TypeError, ValueError):
                continue
            if (np.isnan(p_val) or np.isnan(l_val) or
                    gene in ("nan", "") or gene not in var_idx):
                continue
            if p_val >= padj_max or abs(l_val) < lfc_min:
                continue
            gi = var_idx[gene]
            pct_in  = _pct(gi, mask_in)  if mask_in  is not None else None
            pct_out = _pct(gi, mask_out) if mask_out is not None else None
            if pct_in is not None and pct_in < min_pct_in:
                continue
            record = {
                "cluster": cl_str,
                "gene":    gene,
                "log2fc":  round(l_val, 3),
                "padj":    round(p_val, 6),
                "pct_in":  round(pct_in, 3)  if pct_in  is not None else None,
                "pct_out": round(pct_out, 3) if pct_out is not None else None,
            }
            sig.append(record)
            rows_for_csv.append(record)

        n_sig_by_cluster[cl_str]  = len(sig)
        n_significant_total      += len(sig)
        de_genes_by_cluster[cl_str] = sig[:top_n]

    # Persist full table as CSV for the report
    csv_path = str(Path(output_dir) / "de_per_cluster.csv")
    try:
        pd.DataFrame(rows_for_csv).to_csv(csv_path, index=False)
    except Exception:
        csv_path = None

    return {
        "status":                "success",
        "groupby":               groupby,
        "marker_source":         marker_source,
        "padj_max":              padj_max,
        "lfc_min":               lfc_min,
        "min_pct_in":            min_pct_in,
        "n_clusters":            n_clusters,
        "n_significant_total":   n_significant_total,
        "n_sig_by_cluster":      n_sig_by_cluster,
        "de_genes_by_cluster":   de_genes_by_cluster,
        "output_csv":            csv_path,
    }


if __name__ == "__main__":
    run_script(rna_de_per_cluster)
