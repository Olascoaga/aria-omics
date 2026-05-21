"""
ARIA RNA CellTypist Annotation Script
--------------------------------------
Database-backed cell type annotation. Implements the "code guarantees"
half of the LLM-proposes-code-guarantees pattern: CellTypist gives a
defensible label per cell from a pretrained model, the agent's LLM only
reinterprets / refines on top of that label.

Executed inside aria-rna-env by EnvironmentManager (celltypist 1.7.1+).

Input params:
    data_path:      str  — path to clustered .h5ad (log1p(CPM), sum=1e4)
    organism:       str  — "Homo sapiens" | "Mus musculus" (default: Homo sapiens)
    model:          str  (optional) — explicit CellTypist model name (.pkl)
                                       If omitted, picked from tissue_hint.
    tissue_hint:    str  (optional) — "pbmc" | "blood" | "brain" | "fetal"
                                       | "kidney" | "lung" | "intestine"
                                       Picks a default model if `model` is empty.
    cluster_col:    str  (optional) — obs column for cluster IDs (default: "leiden")
    majority_voting: bool (optional) — collapse labels per cluster (default: True)
    over_clustering: str (optional) — column for majority voting groups
                                       (default: cluster_col)
    output_dir:     str  (optional)

Output:
    {
      "status":              "success" | "skipped" | "error",
      "model_used":          str,
      "predictions_path":    str   — CSV with per-cell labels
      "per_cluster": {
          cluster_id: {
              "label":        str,
              "frequency":    float,    # fraction of cluster with this label
              "n_cells":      int,
              "alt_labels":   [{label, frequency}, ...]   # runner-ups
          },
          ...
      },
      "output_path":         str   — annotated .h5ad (with cell_type_celltypist obs)
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


# Tissue hint → default CellTypist model.
# These models ship with celltypist's hub; the script downloads them on
# first use (~30MB each, cached under ~/.celltypist/data/models).
_DEFAULT_MODELS = {
    "pbmc":      "Immune_All_Low.pkl",     # 98 fine-grained immune types
    "blood":     "Immune_All_Low.pkl",
    "immune":    "Immune_All_Low.pkl",
    "broad":     "Immune_All_High.pkl",    # 32 broad immune types
    "fetal":     "Pan_Fetal_Human.pkl",
    "brain":     "Adult_Human_PrefrontalCortex.pkl",
    "cortex":    "Adult_Human_PrefrontalCortex.pkl",
    "kidney":    "Adult_Human_Kidney.pkl",
    "lung":      "Human_Lung_Atlas.pkl",
    "intestine": "Cells_Intestinal_Tract.pkl",
    "gut":       "Cells_Intestinal_Tract.pkl",
    "skin":      "Adult_Human_Skin.pkl",
}

_MOUSE_MODELS = {
    "brain":  "Mouse_Whole_Brain.pkl",
    "immune": "Mouse_Isocortex_Hippocampus.pkl",
}


def rna_celltypist(params: dict) -> dict:
    import os
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/aria_numba_cache")
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from pathlib import Path
    from aria.utils.safe_h5ad import read_h5ad

    data_path       = params["data_path"]
    organism        = params.get("organism", "Homo sapiens").lower()
    model_explicit  = params.get("model")
    tissue_hint     = (params.get("tissue_hint") or "").lower().strip()
    cluster_col     = params.get("cluster_col", "leiden")
    majority_voting = bool(params.get("majority_voting", True))
    output_dir      = params.get("output_dir", str(Path(data_path).parent))
    out_dir         = Path(output_dir)

    # Resolve model
    if model_explicit:
        model_name = model_explicit
    else:
        catalog = _MOUSE_MODELS if "musculus" in organism else _DEFAULT_MODELS
        model_name = catalog.get(tissue_hint, "Immune_All_Low.pkl")

    try:
        import celltypist
        from celltypist import models
    except ImportError:
        return {
            "status":     "error",
            "error_type": "CellTypistMissing",
            "details":    "celltypist not installed in aria-rna-env. "
                          "Run: pip install celltypist",
        }

    # Download model if not cached. download_models with model_name skips
    # download when present, so this is cheap on subsequent runs.
    try:
        models.download_models(model=model_name, force_update=False)
    except Exception as e:
        return {
            "status":     "error",
            "error_type": "ModelDownloadFailed",
            "details":    f"Could not download CellTypist model "
                          f"'{model_name}': {e}",
        }

    adata = read_h5ad(data_path)

    def _sample_matrix(a):
        n = min(int(a.n_obs), 2000)
        if n <= 0:
            return a.X
        idx = np.linspace(0, int(a.n_obs) - 1, n, dtype=int)
        return a.X[idx, :]

    def _looks_lognorm_10k(a) -> bool:
        try:
            x = _sample_matrix(a)
            if hasattr(x, "toarray"):
                x = x.toarray()
            if np.nanmin(x) < -1e-8 or np.nanmax(x) > 30:
                return False
            sums = np.expm1(x).sum(axis=1)
            sums = sums[np.isfinite(sums) & (sums > 0)]
            if len(sums) == 0:
                return False
            med = float(np.median(sums))
            return 7000 <= med <= 13000
        except Exception:
            return False

    # CellTypist validates both .X and .raw.X. Clustering leaves .X scaled
    # after PCA, so use a clean log1p-normalized AnnData copy and drop .raw
    # before annotation to avoid validating stale scaled/raw matrices.
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        raw_adata.obs = adata.obs.copy()
        if _looks_lognorm_10k(raw_adata):
            adata = raw_adata

    if not _looks_lognorm_10k(adata):
        if adata.X.min() < 0:
            return {
                "status": "error",
                "error_type": "InvalidExpressionMatrix",
                "details": (
                    "Expression matrix has negative/scaled values and no "
                    "valid log1p-normalized raw matrix for CellTypist."
                ),
            }
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    adata.raw = None

    # Predict
    pred = celltypist.annotate(
        adata,
        model=model_name,
        majority_voting=majority_voting,
        over_clustering=(cluster_col if cluster_col in adata.obs.columns else None),
    )
    pred_df = pred.predicted_labels.copy()  # cells × {predicted_labels, ...}

    # Write per-cell labels back into adata.obs
    prediction_label_col = (
        "majority_voting" if "majority_voting" in pred_df.columns
        else "predicted_labels"
    )
    adata.obs["cell_type_celltypist"] = pred_df[prediction_label_col].astype(str).values

    # Aggregate to per-cluster summary
    per_cluster: dict = {}
    if cluster_col in adata.obs.columns:
        for cl in sorted(map(str, adata.obs[cluster_col].unique())):
            mask = adata.obs[cluster_col].astype(str) == cl
            labels = adata.obs.loc[mask, "cell_type_celltypist"]
            n = int(mask.sum())
            if n == 0:
                continue
            value_counts = labels.value_counts()
            top_label    = str(value_counts.index[0])
            top_freq     = round(float(value_counts.iloc[0]) / n, 3)
            alt = [
                {"label": str(label), "frequency": round(float(c) / n, 3)}
                for label, c in value_counts.iloc[1:4].items()
            ]
            per_cluster[cl] = {
                "label":      top_label,
                "frequency":  top_freq,
                "n_cells":    n,
                "alt_labels": alt,
            }

    # Persist outputs
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_csv = str(out_dir / "celltypist_predictions.csv")
    pred_df.to_csv(pred_csv, index=True)

    annotated_path = str(out_dir / "annotated.h5ad")
    adata.write_h5ad(annotated_path)

    return {
        "status":           "success",
        "model_used":       model_name,
        "label_col":        "cell_type_celltypist",
        "prediction_label_col": prediction_label_col,
        "n_cells":          int(adata.n_obs),
        "n_unique_labels":  int(adata.obs["cell_type_celltypist"].nunique()),
        "per_cluster":      per_cluster,
        "predictions_path": pred_csv,
        "output_path":      annotated_path,
    }


if __name__ == "__main__":
    run_script(rna_celltypist)
