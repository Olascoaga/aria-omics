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
    run_ambient:     bool (optional) — opt-in ambient-RNA decontamination
                                        (SoupX/decontX). Default False; honest
                                        no-op (records reason) when no backend.
    ambient_method:  str  (optional) — "auto"|"decontx"|"soupx" (default: auto)
    expected_doublet_rate: float (optional) — 10x rate-of-thumb (default: 0.06)
    biological_context: dict (optional) — from OrchestratorAgent intent
                                          used to adjust thresholds for
                                          stress/senescence phenotypes
    sample_id:       str (optional) — per-sample label injected into
                                        obs["sample_id"]/obs["batch"]; also
                                        used to name the output as
                                        qc_filtered_{sample_id}.h5ad
    output_dir:      str (optional) — write the filtered .h5ad here instead
                                        of next to the input

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


def _run_ambient_decontamination(adata, params: dict) -> dict:
    """Optional, OPT-IN ambient-RNA decontamination (SoupX/decontX).

    P1-4: this is the documented optional decontamination step. It is honest by
    design (ADR-002): if no supported backend is importable, or no validated
    ARIA wrapper exists for an importable one, it reports ``ran=False`` with the
    reason and NEVER fabricates a corrected matrix. ARIA does not silently alter
    counts; building/validating a real SoupX/decontX wrapper is deferred
    modality work, so until then the step records why it did not run rather than
    inventing a result.
    """
    import importlib

    method = str(params.get("ambient_method", "auto")).lower()
    # Supported backends are heavy/optional and not bundled by default.
    candidates = []
    if method in ("auto", "decontx"):
        candidates.append(("decontx", "decontx"))
    if method in ("auto", "soupx"):
        candidates.append(("soupx", "SoupX"))

    for mod, label in candidates:
        try:
            importlib.import_module(mod)
        except Exception:
            continue
        # A backend is present, but we will not fabricate a correction without a
        # validated wrapper — report honestly that it was not applied.
        return {
            "ran": False,
            "method": None,
            "available_backend": mod,
            "reason": (f"{label} is importable but ARIA has no validated "
                       f"decontamination wrapper yet; counts were left "
                       f"uncorrected (no fabricated correction)."),
        }

    return {
        "ran": False,
        "method": None,
        "available_backend": None,
        "reason": ("no supported ambient-correction backend installed "
                   "(decontx/SoupX); install one and rerun with "
                   "run_ambient=True to enable decontamination."),
    }


def _cache_params(params: dict) -> dict:
    bio_ctx = params.get("biological_context") or {}
    return {
        "organism": str(params.get("organism", "Homo sapiens")),
        "mt_threshold": (
            None if "mt_threshold" not in params
            else float(params.get("mt_threshold"))
        ),
        "min_genes": int(params.get("min_genes", 200)),
        "min_cells": int(params.get("min_cells", 3)),
        "initial_min_genes": int(params.get("initial_min_genes", 200)),
        "initial_min_cells": int(params.get("initial_min_cells", 3)),
        "run_scrublet": bool(params.get("run_scrublet", True)),
        "run_ambient": bool(params.get("run_ambient", False)),
        "ambient_method": str(params.get("ambient_method", "auto")),
        "expected_doublet_rate": float(params.get("expected_doublet_rate", 0.06)),
        "sample_id": None if params.get("sample_id") is None else str(params.get("sample_id")),
        "user_question": str(bio_ctx.get("user_question", "")),
    }


def _cache_matches(cached: dict, expected: dict) -> bool:
    return cached.get("cache_version") == 3 and cached.get("cache_params") == expected


def rna_qc(params: dict) -> dict:
    import numpy as np
    import scanpy as sc
    from pathlib import Path
    from aria.utils.safe_h5ad import h5ad_is_readable, read_h5ad

    data_path             = params["data_path"]
    organism              = params.get("organism", "Homo sapiens")
    bio_ctx               = params.get("biological_context", {})
    run_scrublet          = bool(params.get("run_scrublet", True))
    expected_doublet_rate = float(params.get("expected_doublet_rate", 0.06))
    # Multi-sample workflow: caller injects a per-sample label so the output
    # filename and obs columns are sample-aware (avoiding qc_filtered.h5ad
    # collisions when multiple inputs share an output directory).
    sample_id             = params.get("sample_id")
    output_dir            = params.get("output_dir")
    cache_params          = _cache_params(params)

    # ── Load data ─────────────────────────────────────────────────────────
    path = Path(data_path)
    
 # Validate path exists before h5py tries to open it
    if not path.exists():
        return {
            "status":     "error",
            "error_type": "FileNotFound",
            "details":    f"Path does not exist: {data_path}",
        }

    out_dir  = Path(output_dir) if output_dir else Path(data_path).parent
    out_name = (f"qc_filtered_{sample_id}.h5ad"
                if sample_id else "qc_filtered.h5ad")
    output_path = out_dir / out_name
    summary_path = output_path.with_suffix(".summary.json")

    if output_path.exists() and output_path.stat().st_mtime >= path.stat().st_mtime:
        try:
            if summary_path.exists():
                import json
                with open(summary_path) as f:
                    cached = json.load(f)
                if not _cache_matches(cached, cache_params):
                    raise ValueError("stale QC cache")
                if not h5ad_is_readable(output_path):
                    raise ValueError("unreadable QC cache")
                cached["resumed"] = True
                cached["warnings"] = (
                    cached.get("warnings", []) +
                    [f"[resume] scRNA QC skipped for {sample_id or path.stem} "
                     f"(valid filtered h5ad exists)"]
                )
                return cached

            # Legacy filtered h5ads without signed summaries are rerun so
            # threshold changes cannot silently reuse stale QC decisions.
            pass
        except Exception:
            pass
    
    if path.suffix == ".h5ad":
        adata = read_h5ad(str(path))
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

    if adata.n_obs == 0 or adata.n_vars == 0:
        return {
            "status":     "error",
            "error_type": "EmptyAnnData",
            "details":    (
                f"Input contains no cells or genes: {adata.n_obs} cells x "
                f"{adata.n_vars} genes. This is often a stale ARIA "
                "intermediate such as qc_filtered.h5ad from a failed run."
            ),
            "n_cells_before": int(adata.n_obs),
            "n_cells_after":  0,
            "n_genes_before": int(adata.n_vars),
            "n_genes_after":  0,
        }

    # CellRanger 10x matrices commonly contain duplicate gene symbols (multiple
    # Ensembl IDs collapsing to the same HGNC name). Without dedup, downstream
    # HVG/PCA/Scanpy ops emit warnings and can produce non-deterministic var
    # ordering after filter_genes.
    adata.var_names_make_unique()

    # Annotate per-sample provenance for downstream concat / batch correction.
    # We always stamp these columns when sample_id is provided so rna_concat
    # can stack without further annotation work.
    if sample_id:
        adata.obs["sample_id"] = str(sample_id)
        adata.obs["batch"]     = str(sample_id)

    warnings_list  = []

    # ── Optional ambient-RNA decontamination (P1-4, opt-in) ───────────────
    # Off by default. Runs on raw counts before downstream steps. Honest: if no
    # supported backend is available it records why it did not run and leaves
    # counts untouched (never fabricates a correction).
    ambient_report = {"ran": False, "method": None, "reason": "not_requested"}
    if bool(params.get("run_ambient", False)):
        ambient_report = _run_ambient_decontamination(adata, params)
        if not ambient_report.get("ran"):
            warnings_list.append(
                f"Ambient-RNA decontamination requested but not applied: "
                f"{ambient_report.get('reason', 'unavailable')}"
            )

    n_barcodes_raw = adata.n_obs
    existing_gene_col = next(
        (c for c in ("nFeature_RNA", "n_genes_by_counts", "n_genes")
         if c in adata.obs.columns),
        None,
    )
    existing_count_col = next(
        (c for c in ("nCount_RNA", "total_counts", "n_counts")
         if c in adata.obs.columns),
        None,
    )
    existing_mt_col = next(
        (c for c in ("percent.mt", "pct_counts_mt", "percent_mito")
         if c in adata.obs.columns),
        None,
    )
    use_existing_h5ad_qc = bool(
        path.suffix == ".h5ad"
        and existing_gene_col
        and existing_count_col
        and existing_mt_col
    )

    if use_existing_h5ad_qc:
        warnings_list.append(
            "Using existing h5ad obs QC metrics "
            f"({existing_gene_col}, {existing_count_col}, {existing_mt_col}); "
            "skipping X-based empty-droplet/gene filtering."
        )
        adata.obs["n_genes_by_counts"] = adata.obs[existing_gene_col].astype(float)
        adata.obs["total_counts"] = adata.obs[existing_count_col].astype(float)
        adata.obs["pct_counts_mt"] = adata.obs[existing_mt_col].astype(float)
        n_cells_before = adata.n_obs
        n_genes_before = adata.n_vars

        # Preprocessed Seurat/Scanpy h5ads often store scaled/log-normalized
        # values in X. Scrublet requires raw counts and would be invalid here.
        if run_scrublet:
            run_scrublet = False
            warnings_list.append(
                "Scrublet skipped because h5ad provides precomputed QC metrics; "
                "doublet status from obs is preserved when present."
            )
    else:
        # ── Drop empty droplets BEFORE estimating MAD thresholds ─────────
        # When data_path points to a raw_feature_bc_matrix (CellRanger raw
        # output), ~95% of barcodes are empty drops. Their dominance in the
        # count distribution collapses median→0 and MAD→0, producing
        # mt_high=NaN and nuking every real cell. Apply a coarse initial filter
        # so adaptive MAD thresholds are estimated on plausible cells only.
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

    if not use_existing_h5ad_qc:
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
    if use_existing_h5ad_qc:
        qc_mask = (
            np.isfinite(adata.obs["n_genes_by_counts"].astype(float))
            & np.isfinite(adata.obs["pct_counts_mt"].astype(float))
            & (adata.obs["n_genes_by_counts"].astype(float) >= min_genes)
            & (adata.obs["pct_counts_mt"].astype(float) <= mt_threshold)
        )
        adata = adata[qc_mask].copy()
    else:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata = adata[adata.obs["pct_counts_mt"] <= mt_threshold].copy()

    if adata.n_obs == 0:
        return {
            "status":     "error",
            "error_type": "NoCellsAfterQC",
            "details":    (f"No cells survived QC "
                           f"(min_genes={min_genes}, "
                           f"mt_threshold={mt_threshold:.3f})."),
            "n_cells_before": int(n_cells_before),
            "n_cells_after":  0,
            "min_genes_used": int(min_genes),
            "mt_threshold_used": float(mt_threshold),
            "warnings": warnings_list,
        }

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
    if n_cells_after == 0:
        return {
            "status":     "error",
            "error_type": "NoCellsAfterQC",
            "details":    "No cells remained after doublet/QC filtering.",
            "n_cells_before": int(n_cells_before),
            "n_cells_after":  0,
            "min_genes_used": int(min_genes),
            "mt_threshold_used": float(mt_threshold),
            "warnings": warnings_list,
        }
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
    out_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(output_path))

    result = {
        "status":            "success",
        "sample_id":         sample_id,
        "n_cells_before":    int(n_cells_before),
        "n_cells_after":     int(n_cells_after),
        "n_genes_before":    int(n_genes_before),
        "n_genes_after":     int(n_genes_after),
        "pct_removed":       float(pct_removed),
        "mt_threshold_used": float(mt_threshold),
        "min_genes_used":    int(min_genes),
        "stress_context":    bool(is_stress_context),
        "ambient_correction": ambient_report,
        "scrublet":          scrublet_report,
        "warnings":          warnings_list,
        "output_path":       str(output_path),
        "mt_stats": {
            "mean": round(float(adata.obs["pct_counts_mt"].mean()), 3),
            "max":  round(float(adata.obs["pct_counts_mt"].max()), 3),
        },
        "cache_version": 3,
        "cache_params": cache_params,
    }
    try:
        import json
        with open(summary_path, "w") as f:
            json.dump(result, f, indent=2)
    except Exception:
        pass
    return result


if __name__ == "__main__":
    run_script(rna_qc)
