"""
ARIA RNA QC Script
------------------
Runs quality control on scRNA-seq or bulk RNA-seq data.
Executed inside aria-rna-env by EnvironmentManager.

QC pipeline:
  1. MT% / total_counts / n_genes_by_counts filtering (MAD-adaptive,
     stress-context-aware via biological_context).
  2. Doublet detection via Scrublet (single-cell only, opt-out with
     `run_scrublet=False`). Operates on raw counts BEFORE normalization.
  3. Persist filtered .h5ad with doublet_score / predicted_doublet in obs.

Input params:
    data_path:       str  — path to .h5ad, .h5 (10x), or MEX directory
    organism:        str  — "Homo sapiens", "Mus musculus", etc.
    mt_threshold:    float (optional) — max mitochondrial % (default: adaptive MAD)
    min_genes:       int  (optional) — min genes per cell (default: 200)
    min_cells:       int  (optional) — min cells per gene (default: 3)
    run_scrublet:    bool (optional) — disable doublet detection (default: True)
    expected_doublet_rate: float (optional) — 10x rate-of-thumb (default: 0.06)
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
      "scrublet": {
          "ran":               bool,
          "n_doublets":        int,
          "doublet_rate":      float,   # observed
          "threshold_used":    float,
          "expected_rate":     float,
      },
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

    data_path             = params["data_path"]
    organism              = params.get("organism", "Homo sapiens")
    bio_ctx               = params.get("biological_context", {})
    run_scrublet          = bool(params.get("run_scrublet", True))
    expected_doublet_rate = float(params.get("expected_doublet_rate", 0.06))

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

    # CellRanger 10x matrices commonly contain duplicate gene symbols (multiple
    # Ensembl IDs collapsing to the same HGNC name). Without dedup, downstream
    # HVG/PCA/Scanpy ops emit warnings and can produce non-deterministic var
    # ordering after filter_genes.
    adata.var_names_make_unique()

    warnings_list  = []
    n_barcodes_raw = adata.n_obs

    # ── Drop empty droplets BEFORE estimating MAD thresholds ─────────────
    # When data_path points to a raw_feature_bc_matrix (CellRanger raw output),
    # ~95% of barcodes are empty drops. Their dominance in the count
    # distribution collapses median→0 and MAD→0, producing mt_high=NaN and
    # nuking every real cell. Apply a coarse initial filter so the adaptive
    # MAD thresholds are estimated on plausible cells only.
    initial_min_genes = int(params.get("initial_min_genes", 200))
    initial_min_cells = int(params.get("initial_min_cells", 3))
    sc.pp.filter_cells(adata, min_genes=initial_min_genes)
    sc.pp.filter_genes(adata, min_cells=initial_min_cells)

    if adata.n_obs == 0:
        return {
            "status":     "error",
            "error_type": "NoCellsAfterInitialFilter",
            "details":    (f"No barcodes survived initial filter "
                           f"(min_genes={initial_min_genes}, "
                           f"min_cells={initial_min_cells}) from "
                           f"{n_barcodes_raw} raw barcodes."),
        }

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    if n_barcodes_raw > n_cells_before * 1.5:
        warnings_list.append(
            f"Input looks like a raw_feature_bc_matrix: "
            f"{n_barcodes_raw} barcodes → {n_cells_before} after empty-drop "
            f"removal (min_genes={initial_min_genes})."
        )

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
    # Use nanmedian: pct_counts_mt is NaN when a cell has total_counts=0,
    # which can still slip through if initial_min_genes is lowered.
    def mad_bounds(values: np.ndarray, n_mad: float = 3.0):
        median = float(np.nanmedian(values))
        mad    = float(np.nanmedian(np.abs(values - median)))
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

    # ── Doublet detection (Scrublet) ─────────────────────────────────────
    # Run BEFORE normalization on raw counts. Skip if scrublet not
    # installed or caller explicitly disabled it.
    scrublet_report = {
        "ran":            False,
        "n_doublets":     0,
        "doublet_rate":   0.0,
        "threshold_used": None,
        "expected_rate":  expected_doublet_rate,
    }
    if run_scrublet and adata.n_obs >= 50:
        try:
            import scrublet as scr
            counts = adata.X
            if hasattr(counts, "toarray"):
                # Scrublet wants a sparse or dense numpy matrix; CSR is fine.
                pass
            scrub = scr.Scrublet(
                counts,
                expected_doublet_rate=expected_doublet_rate,
                random_state=0,
            )
            doublet_scores, predicted = scrub.scrub_doublets(
                min_counts=2, min_cells=3,
                min_gene_variability_pctl=85,
                n_prin_comps=30, verbose=False,
            )
            # If Scrublet failed to converge on a threshold it returns None.
            # Fall back to a conservative 0.5 cutoff so we don't drop cells
            # we can't justify dropping.
            if predicted is None:
                predicted = doublet_scores > 0.5
                thr_used  = 0.5
                warnings_list.append(
                    "Scrublet could not auto-derive a doublet threshold; "
                    "fell back to score>0.5. Inspect doublet_score in obs."
                )
            else:
                thr_used = float(getattr(scrub, "threshold_", 0.0))

            adata.obs["doublet_score"]     = doublet_scores
            adata.obs["predicted_doublet"] = predicted.astype(bool)

            n_doublets = int(predicted.sum())
            obs_rate   = round(float(n_doublets) / adata.n_obs, 4)

            # Drop predicted doublets from the filtered AnnData.
            adata = adata[~adata.obs["predicted_doublet"]].copy()

            scrublet_report.update({
                "ran":            True,
                "n_doublets":     n_doublets,
                "doublet_rate":   obs_rate,
                "threshold_used": round(thr_used, 4),
            })

            # Warn if the observed rate is wildly off from the 10x rule of
            # thumb — possible cell-type composition issue or batch artefact.
            if obs_rate > expected_doublet_rate * 3:
                warnings_list.append(
                    f"Scrublet flagged {obs_rate*100:.1f}% doublets "
                    f"(expected ~{expected_doublet_rate*100:.1f}%). "
                    f"Review raw barcode QC and loading concentration."
                )
        except ImportError:
            warnings_list.append(
                "scrublet not available — doublet detection skipped. "
                "Install scrublet in aria-rna-env to enable."
            )
        except Exception as e:
            warnings_list.append(
                f"Scrublet failed ({type(e).__name__}: {str(e)[:120]}); "
                f"continuing without doublet filtering."
            )

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
        "scrublet":          scrublet_report,
        "warnings":          warnings_list,
        "output_path":       output_path,
        "mt_stats": {
            "mean": round(float(adata.obs["pct_counts_mt"].mean()), 3),
            "max":  round(float(adata.obs["pct_counts_mt"].max()), 3),
        },
    }


if __name__ == "__main__":
    run_script(rna_qc)
