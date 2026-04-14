"""
ARIA RNA QC Script
------------------
Runs quality control on scRNA-seq or bulk RNA-seq data.
Executed inside aria-rna-env by EnvironmentManager.

Input params:
    data_path:       str  — path to .h5ad, .h5 (10x), or MEX directory
    organism:        str  — "Homo sapiens", "Mus musculus", etc.
    mt_threshold:    float (optional) — max mitochondrial % (default: adaptive MAD)
    min_genes:       int  (optional) — min genes per cell (default: 200)
    min_cells:       int  (optional) — min cells per gene (default: 3)
    biological_context: dict (optional) — from OrchestratorAgent intent
                                          used to adjust thresholds for
                                          stress/senescence phenotypes

Output:
    {
      "status":           "success",
      "n_cells_before":   int,
      "n_cells_after":    int,
      "n_genes_before":   int,
      "n_genes_after":    int,
      "pct_removed":      float,
      "mt_threshold_used": float,
      "min_genes_used":   int,
      "warnings":         [str],
      "output_path":      str   — path to filtered .h5ad
    }
"""

from __future__ import annotations

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from aria.scripts._base import run_script


def rna_qc(params: dict) -> dict:
    import numpy as np
    import scanpy as sc
    from pathlib import Path

    data_path = params["data_path"]
    organism  = params.get("organism", "Homo sapiens")
    bio_ctx   = params.get("biological_context", {})

    # ── Load data ─────────────────────────────────────────────────────────
    path = Path(data_path)
    
 # Validate path exists before h5py tries to open it
    if not path.exists():
        return {
            "status":     "error",
            "error_type": "FileNotFound",
            "details":    f"Path does not exist: {data_path}",
        }
    
    if path.suffix == ".h5ad":
        adata = sc.read_h5ad(str(path))
    elif path.suffix == ".h5":
        adata = sc.read_10x_h5(str(path))
    elif path.is_dir():
        adata = sc.read_10x_mtx(str(path), var_names="gene_symbols", cache=True)
    else:
        return {
            "status":     "error",
            "error_type": "UnsupportedFormat",
            "details":    f"Cannot load: {data_path}. "
                          f"Supported: .h5ad, .h5, MEX directory.",
        }

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    warnings_list  = []

    # ── Mitochondrial gene prefix ─────────────────────────────────────────
    # Human and mouse use MT- / mt- respectively
    is_human = "sapiens" in organism.lower()
    is_mouse = "musculus" in organism.lower()
    mt_prefix = "MT-" if is_human else "mt-" if is_mouse else "MT-"

    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"],
        percent_top=None, log1p=False, inplace=True,
    )

    # ── Adaptive MAD thresholds ───────────────────────────────────────────
    def mad_bounds(values: np.ndarray, n_mad: float = 3.0):
        median = np.median(values)
        mad    = np.median(np.abs(values - median))
        return median - n_mad * mad, median + n_mad * mad

    counts_low, counts_high = mad_bounds(adata.obs["total_counts"].values)
    genes_low,  _           = mad_bounds(adata.obs["n_genes_by_counts"].values)
    _,          mt_high     = mad_bounds(adata.obs["pct_counts_mt"].values)

    # ── Biological context adjustments ───────────────────────────────────
    # Senescent, stressed, or apoptotic cells have elevated MT% — do not
    # discard them if user explicitly mentioned these phenotypes
    STRESS_KEYWORDS = [
        "senescen", "stress", "hypox", "apoptot", "dying",
        "activat", "inflam", "exhausted",
    ]
    question = bio_ctx.get("user_question", "").lower()
    is_stress_context = any(kw in question for kw in STRESS_KEYWORDS)

    if is_stress_context:
        mt_ceiling = 35.0   # relaxed for stress phenotypes
        warnings_list.append(
            "MT% threshold relaxed to 35% due to stress/activation context. "
            "Review MT% distribution before publishing."
        )
    else:
        mt_ceiling = 25.0   # standard

    # Apply user overrides if provided
    mt_threshold = float(params.get("mt_threshold", min(mt_high, mt_ceiling)))
    min_genes    = int(params.get("min_genes", max(int(genes_low), 200)))
    min_cells    = int(params.get("min_cells", 3))

    # ── Filter ────────────────────────────────────────────────────────────
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] <= mt_threshold].copy()

    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    pct_removed   = round(
        (n_cells_before - n_cells_after) / n_cells_before * 100, 2
    )

    # ── Quality checks ────────────────────────────────────────────────────
    if n_cells_after < 100:
        warnings_list.append(
            f"Only {n_cells_after} cells passed QC. "
            f"Analysis results will be unreliable. "
            f"Consider relaxing thresholds."
        )

    lib_sizes = adata.obs["total_counts"]
    size_ratio = lib_sizes.max() / max(lib_sizes.min(), 1)
    if size_ratio > 10:
        warnings_list.append(
            f"Library size range: {size_ratio:.1f}x. "
            f"Verify normalization strategy."
        )

    if pct_removed > 40:
        warnings_list.append(
            f"{pct_removed}% of cells removed. "
            f"Review raw data quality or consider relaxing QC thresholds."
        )

    # ── Save filtered data ────────────────────────────────────────────────
    output_path = str(Path(data_path).parent / "qc_filtered.h5ad")
    adata.write_h5ad(output_path)

    return {
        "status":            "success",
        "n_cells_before":    int(n_cells_before),
        "n_cells_after":     int(n_cells_after),
        "n_genes_before":    int(n_genes_before),
        "n_genes_after":     int(n_genes_after),
        "pct_removed":       float(pct_removed),
        "mt_threshold_used": float(mt_threshold),
        "min_genes_used":    int(min_genes),
        "stress_context":    bool(is_stress_context),
        "warnings":          warnings_list,
        "output_path":       output_path,
        "mt_stats": {
            "mean": round(float(adata.obs["pct_counts_mt"].mean()), 3),
            "max":  round(float(adata.obs["pct_counts_mt"].max()), 3),
        },
    }


if __name__ == "__main__":
    run_script(rna_qc)
