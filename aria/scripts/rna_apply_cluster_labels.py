"""
Apply per-cluster biological labels to an AnnData obs column.

Used by scRNAAgent when database-backed annotation is unavailable but a
conservative marker-based label map exists.
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def rna_apply_cluster_labels(params: dict) -> dict:
    import os
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/aria_numba_cache")
    import scanpy as sc
    from pathlib import Path

    data_path = params["data_path"]
    labels = params.get("labels") or {}
    cluster_col = params.get("cluster_col", "leiden")
    label_col = params.get("label_col", "cell_type_marker")
    output_dir = params.get("output_dir", str(Path(data_path).parent))
    out_dir = Path(output_dir)

    adata = sc.read_h5ad(data_path)
    if cluster_col not in adata.obs.columns:
        return {
            "status": "error",
            "error_type": "MissingClusterColumn",
            "details": f"Column '{cluster_col}' not found in obs.",
        }

    def label_for(cluster):
        info = labels.get(str(cluster), {})
        if isinstance(info, dict):
            return str(info.get("cell_type") or info.get("label") or cluster)
        return str(info or cluster)

    adata.obs[label_col] = [
        label_for(cl) for cl in adata.obs[cluster_col].astype(str)
    ]

    out_dir.mkdir(parents=True, exist_ok=True)
    out = str(out_dir / "annotated_marker.h5ad")
    adata.write_h5ad(out)
    return {
        "status": "success",
        "output_path": out,
        "label_col": label_col,
        "n_labels": int(adata.obs[label_col].nunique()),
    }


if __name__ == "__main__":
    run_script(rna_apply_cluster_labels)
