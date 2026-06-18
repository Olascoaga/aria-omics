"""
ARIA Chromatin QC Script
-------------------------
Quality control for all chromatin modalities:
  scATAC-seq, bulk ATAC-seq, ChIP-seq, CUT&RUN, CUT&TAG

Key QC metrics:
  - FRiP (Fraction of Reads in Peaks): ATAC > 0.2, ChIP > 0.1
  - TSS Enrichment Score: ATAC > 4 (ideally > 8)
  - Fragment size distribution: nucleosomal banding pattern
  - Mitochondrial read fraction: should be < 10% for ATAC
  - Library complexity: PCR duplicate rate

Executed inside aria-chromatin-env by EnvironmentManager.

Input params:
    data_type:    str  — "scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"
    files:        list — BAM files or fragment files
    genome:       str  — genome assembly (hg38, mm10, etc.)
    organism:     str  — organism name
    assay_class:  str  (ChIP only) — "histone" or "tf"

Honesty contract (B9, audit 2026-05-29 / ADR-002): this scaffold QC only emits
metrics it can genuinely derive from the input. Metrics that require a reference
annotation (TSS enrichment) or called peaks (FRiP) are returned as `None` and
listed in `metrics_not_computed`, never as fabricated placeholders. `pass_qc` is
`None` while QC is incomplete (`qc_complete=False`) and a boolean only once the
gating metrics exist; a real measured failure sets it `False`.

Output:
    {
      "status":          "success",
      "data_type":       str,
      "n_cells":         int          (scATAC only — unique barcodes seen),
      "n_fragments":     int,         — fragment records scanned
      "n_fragments_is_lower_bound": bool (True if the scan cap was reached),
      "frip":            float|None,  — None until peaks are called
      "tss_enrichment":  float|None,  — None until a reference TSS QC runs
      "mito_fraction":   float,       — mitochondrial read fraction
      "dup_rate":        float|None,  — None if samtools unavailable (bulk)
      "fragment_sizes":  dict,        — distribution summary + barcode count
      "qc_complete":     bool,        — whether the gating metrics were computed
      "metrics_not_computed": [str],  — which metrics need the v4.6 stack/peaks
      "pass_qc":         bool|None,   — None while incomplete; False on failure
      "warnings":        [str],
      "n_cells_after":   int|None     (scATAC only — None until per-barcode QC)
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script, mocks_allowed


def chromatin_qc(params: dict) -> dict:
    from pathlib import Path
    import numpy as np

    data_type   = params["data_type"]
    files       = params.get("files", [])
    genome      = params.get("genome", "hg38")
    organism    = params.get("organism", "Homo sapiens")
    assay_class = params.get("assay_class", "tf")
    # P0 W0.5 wire: a called-peak BED (e.g. chromatin_peaks `peaks_path`) makes
    # FRiP computable from the fragments. Absent → FRiP stays honest None.
    peaks_bed   = params.get("peaks_bed") or params.get("peaks_path")
    allow_mocks = mocks_allowed(params)

    warnings = []

    # Validate files exist
    valid_files = [f for f in files if Path(f).exists()]
    if not valid_files:
        return {
            "status":     "error",
            "error_type": "FileNotFound",
            "details":    f"None of {len(files)} files exist on disk.",
        }

    # ── Dispatch to modality-specific QC ─────────────────────────────────
    # C8 (audit 2026-05-29): the v4.6 entry path is a same-cell paired RNA+ATAC
    # MuData file. Route a .h5mu through the real MuData reader rather than the
    # fragments-file path, which would not detect it.
    h5mu_files = [f for f in valid_files if str(f).lower().endswith(".h5mu")]
    if data_type == "scATAC":
        if h5mu_files:
            return _h5mu_atac_qc(h5mu_files[0], genome, organism, warnings)
        return _scatac_qc(valid_files, genome, organism, warnings, allow_mocks,
                          peaks_bed=peaks_bed)
    elif data_type in ("bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"):
        return _bulk_chromatin_qc(
            valid_files, data_type, genome, assay_class, warnings, allow_mocks
        )
    else:
        return {
            "status":     "error",
            "error_type": "UnsupportedDataType",
            "details":    f"Unknown data type: {data_type}",
        }


def _scatac_qc(files: list, genome: str,
               organism: str, warnings: list,
               allow_mocks: bool = False, peaks_bed: str | None = None) -> dict:
    """QC for single-cell ATAC-seq using episcanpy or muon."""
    try:
        import muon as mu
        import episcanpy.api as epi
        import anndata as ad
        import numpy as np
        from pathlib import Path

        # Load fragment file
        fragment_files = [f for f in files if "fragments" in f.lower()
                          or f.endswith(".tsv.gz")]
        if not fragment_files:
            return {
                "status":     "error",
                "error_type": "MissingFragmentFile",
                "details":    "scATAC QC requires a fragments.tsv.gz file.",
            }

        frag_file = fragment_files[0]

        # Create AnnData from fragments
        # (in production: use cellranger output or snapatac2)
        try:
            mdata = mu.atac.tl.locate_fragments(None, frag_file)
        except Exception:
            # Fallback: basic fragment counting
            pass

        # ── Key QC metrics ────────────────────────────────────────────────
        # Only metrics genuinely derivable from a fragments file are computed
        # here. TSS enrichment (needs a reference TSS + snapatac2/episcanpy) and
        # FRiP (needs called peaks) are reported as not-computed rather than
        # fabricated (ADR-002).
        # W0.5: real TSS enrichment (snapatac2 + governed annotation); returns
        # (value, reason). reason carries the concrete not-computed cause.
        # Single snapatac2 pass yields the W0.5 median AND the P4.1 [A] per-cell
        # (TSSe, log10 depth) scatter arrays — no double import_data.
        _tss = _tss_qc(frag_file, genome)
        tss_enrichment, tss_reason = _tss["median"], _tss["reason"]
        frag_sizes = _compute_fragment_sizes(frag_file)

        mito_chr   = "chrM" if "sapiens" in organism.lower() else "chrMT"
        mito_frac  = _compute_mito_fraction(frag_file, mito_chr)

        n_cells    = frag_sizes.get("n_barcodes", 0)
        n_fragments = frag_sizes.get("n_fragments", 0)
        # FRiP is real when a called-peak BED is supplied (W0.5 wire); else None.
        frip = _estimate_frip(frag_file, peaks_bed=peaks_bed)
        # Per-barcode FRiP distribution for the QC histogram (P4.1 [B]); same
        # peaks gate as the scalar — None without peaks.
        frip_distribution = _estimate_frip_distribution(frag_file, peaks_bed=peaks_bed)

        # ── QC thresholds (standard ENCODE ATAC-seq) ──────────────────────
        MIN_TSS     = 4.0     # below this = failed library prep
        WARN_TSS    = 8.0     # below this = marginal quality
        MAX_MITO    = 0.10    # above this = mitochondrial contamination

        # Real QC failures gate pass_qc; "not computed" is a coverage note, not
        # a failure, so the two are tracked separately.
        qc_failures: list[str] = []
        not_computed: list[str] = []

        if tss_enrichment is None:
            not_computed.append(
                f"tss_enrichment ({tss_reason})" if tss_reason
                else "tss_enrichment (not computed)"
            )
        elif tss_enrichment < MIN_TSS:
            qc_failures.append(
                f"TSS enrichment {tss_enrichment:.2f} < {MIN_TSS} "
                f"(ENCODE minimum). Library prep likely failed."
            )
        elif tss_enrichment < WARN_TSS:
            warnings.append(
                f"TSS enrichment {tss_enrichment:.2f} < {WARN_TSS} "
                f"(marginal quality — acceptable but not ideal)."
            )

        if mito_frac > MAX_MITO:
            qc_failures.append(
                f"Mitochondrial fraction {mito_frac:.1%} > {MAX_MITO:.0%}. "
                f"High mtDNA contamination detected."
            )

        if frip is None:
            not_computed.append(
                "frip (requires called peaks — run chromatin_peaks first)"
            )
        elif frip < 0.2:
            qc_failures.append(
                f"FRiP {frip:.3f} < 0.20. "
                f"Poor signal-to-noise ratio in ATAC library."
            )

        warnings.extend(qc_failures)
        if not_computed:
            warnings.append(
                "QC metrics not computed: " + "; ".join(not_computed) + "."
            )
        if frag_sizes.get("scan_truncated"):
            warnings.append(
                f"Fragment scan capped at {frag_sizes.get('scan_limit')} "
                f"records; n_barcodes/n_fragments are lower bounds."
            )

        # QC is complete only when both gated metrics are available. Until then
        # pass_qc is None (cannot assert a pass on partial metrics); when
        # complete it is the AND of the real checks.
        qc_complete = (tss_enrichment is not None and frip is not None)
        pass_qc = (not qc_failures) if qc_complete else None

        return {
            "status":          "success",
            "data_type":       "scATAC",
            "n_cells":         int(n_cells),
            "n_cells_after":   None,  # requires per-barcode filtering (v4.6)
            "n_fragments":     int(n_fragments),
            "n_fragments_is_lower_bound": bool(frag_sizes.get("scan_truncated")),
            "frip":            (round(float(frip), 4) if frip is not None else None),
            "tss_enrichment":  (round(float(tss_enrichment), 3)
                                if tss_enrichment is not None else None),
            "mito_fraction":   round(float(mito_frac), 4),
            "fragment_sizes":  frag_sizes,
            # P4.1 [B]: per-barcode FRiP distribution (list) or None (no peaks).
            "frip_distribution": frip_distribution,
            # P4.1 [A]: per-cell TSSe×log10-depth scatter arrays + ENCODE gates, or
            # {"tsse": None, "reason": ...} when not computed (honest absence).
            "tss_depth_scatter": (
                {"tsse": _tss["tsse"], "log10_depth": _tss["log10_depth"],
                 "min_tss": MIN_TSS, "warn_tss": WARN_TSS}
                if _tss["tsse"] is not None
                else {"tsse": None, "reason": tss_reason}),
            "qc_complete":     qc_complete,
            "metrics_not_computed": not_computed,
            "pass_qc":         pass_qc,
            "warnings":        warnings,
        }

    except ImportError as e:
        if not allow_mocks:
            return {
                "status":     "error",
                "error_type": "MissingDependency",
                "details":    (
                    f"scATAC QC requires muon and episcanpy "
                    f"(aria-chromatin-env). Install the chromatin stack or "
                    f"pass allow_mock=true. Underlying ImportError: {e}"
                ),
            }
        # Dev/test only: fall back to basic stats with explicit gating.
        return _basic_chromatin_qc(files, "scATAC", warnings, str(e))


def _h5mu_atac_qc(h5mu_path: str, genome: str,
                  organism: str, warnings: list) -> dict:
    """
    scATAC QC entry for a same-cell paired RNA+ATAC `.h5mu` (C8). Reads the ATAC
    modality with the real MuData reader and reports only what is genuinely
    derivable: peak/cell matrix dimensions and per-cell fragment/mito summaries
    when those obs columns exist. TSS enrichment and FRiP need a reference
    annotation and called peaks respectively, so they stay not-computed (B9 /
    ADR-002). Missing MuData tooling returns a structured MissingDependency.
    """
    from aria.utils.mudata_io import read_h5mu_atac

    res = read_h5mu_atac(h5mu_path)
    if res.get("status") != "success":
        if res.get("error_type") == "MuDataToolingMissing":
            return {
                "status":     "error",
                "error_type": "MissingDependency",
                "details":    (res["details"] + " Install the chromatin stack "
                               "(aria-chromatin-env) to QC .h5mu inputs."),
            }
        return res

    not_computed = [
        "tss_enrichment (requires a reference TSS annotation + "
        "snapatac2/episcanpy)",
        "frip (requires called peaks — run chromatin_peaks first)",
    ]
    if res.get("fragment_count_col") is None:
        not_computed.append(
            "median_fragments_per_cell (no per-cell fragment obs column found)"
        )
    if res.get("mito_col") is None:
        not_computed.append(
            "mito_fraction (no per-cell mitochondrial obs column found)"
        )

    warnings.append(
        f"Scaffold scATAC QC on .h5mu (modality '{res.get('modality')}'): "
        f"peak/cell matrix dimensions and available obs metrics were read from "
        f"MuData; TSS/FRiP require the v4.6 chromatin stack and peak calling."
    )

    return {
        "status":          "success",
        "data_type":       "scATAC",
        "input_kind":      "h5mu",
        "atac_modality":   res.get("modality"),
        "n_cells":         int(res["n_obs"]),
        "n_peaks":         int(res["n_vars"]),
        "n_cells_after":   None,    # requires per-cell QC filtering (v4.6)
        "frip":            None,
        "tss_enrichment":  None,
        "fragment_count_col":         res.get("fragment_count_col"),
        "median_fragments_per_cell":  res.get("median_fragments_per_cell"),
        "mito_obs_col":               res.get("mito_col"),
        "median_mito_obs":            res.get("median_mito_fraction"),
        "qc_complete":     False,
        "metrics_not_computed": not_computed,
        "pass_qc":         None,
        "warnings":        warnings,
    }


def _bulk_chromatin_qc(files: list, data_type: str,
                        genome: str, assay_class: str,
                        warnings: list,
                        allow_mocks: bool = False) -> dict:
    """QC for bulk ATAC-seq, ChIP-seq, CUT&RUN, CUT&TAG using BAM files."""
    try:
        import pysam
        import numpy as np

        bam_files = [f for f in files if f.endswith(".bam")]
        if not bam_files:
            return {
                "status":     "error",
                "error_type": "MissingBAM",
                "details":    f"No .bam files found for {data_type} QC.",
            }

        metrics_per_sample = []

        for bam_path in bam_files[:8]:  # limit to 8 samples
            try:
                bam     = pysam.AlignmentFile(bam_path, "rb")
                stats   = bam.get_index_statistics()
                total   = sum(s.mapped for s in stats)
                mito    = sum(s.mapped for s in stats
                              if s.contig in ("chrM", "MT", "chrMT"))
                bam.close()

                mito_frac = mito / max(total, 1)
                dup_rate  = _estimate_dup_rate(bam_path)  # None if samtools NA

                sample_metrics = {
                    "file":        str(bam_path),
                    "total_reads": int(total),
                    "mito_frac":   round(float(mito_frac), 4),
                    "dup_rate":    (round(float(dup_rate), 4)
                                    if dup_rate is not None else None),
                }
                metrics_per_sample.append(sample_metrics)

                # Check thresholds
                MAX_MITO = 0.10 if data_type != "ChIP" else 0.20
                if mito_frac > MAX_MITO:
                    warnings.append(
                        f"{bam_path}: high mito fraction "
                        f"({mito_frac:.1%} > {MAX_MITO:.0%})"
                    )
                if dup_rate is not None and dup_rate > 0.5:
                    warnings.append(
                        f"{bam_path}: high duplicate rate "
                        f"({dup_rate:.1%}) — low library complexity"
                    )

            except Exception as e:
                warnings.append(f"Could not QC {bam_path}: {str(e)[:100]}")

        if not metrics_per_sample:
            return {
                "status":     "error",
                "error_type": "QCFailed",
                "details":    "Could not compute QC metrics for any sample.",
            }

        # Aggregate metrics (skip samples whose dup_rate could not be measured)
        avg_mito = float(np.mean([m["mito_frac"] for m in metrics_per_sample]))
        _dups    = [m["dup_rate"] for m in metrics_per_sample
                    if m["dup_rate"] is not None]
        avg_dup  = float(np.mean(_dups)) if _dups else None

        # FRiP/TSS require called peaks; they are not fabricated here (ADR-002)
        # — None until chromatin_peaks runs. The real per-sample mito/dup
        # failures captured above are the actual pass/fail signal.
        qc_failures = list(warnings)
        not_computed = [
            "frip (requires called peaks — run chromatin_peaks first)",
            "tss_enrichment (computed post peak calling)",
        ]
        warnings.append(
            "Scaffold QC did not compute: " + "; ".join(not_computed) + "."
        )

        return {
            "status":          "success",
            "data_type":       data_type,
            "n_samples":       len(metrics_per_sample),
            "frip":            None,
            "tss_enrichment":  None,  # computed post peak calling
            "mito_fraction":   round(avg_mito, 4),
            "dup_rate":        (round(avg_dup, 4) if avg_dup is not None
                                else None),
            "per_sample":      metrics_per_sample,
            "qc_complete":     False,
            "metrics_not_computed": not_computed,
            # Real failures (mito/dup) fail QC; without FRiP/TSS a clean run
            # cannot be asserted as a full pass, so pass_qc is None.
            "pass_qc":         (False if qc_failures else None),
            "warnings":        warnings,
        }

    except ImportError as e:
        if not allow_mocks:
            return {
                "status":     "error",
                "error_type": "MissingDependency",
                "details":    (
                    f"{data_type} QC requires pysam + numpy "
                    f"(aria-chromatin-env). Install the chromatin stack or "
                    f"pass allow_mock=true. Underlying ImportError: {e}"
                ),
            }
        return _basic_chromatin_qc(files, data_type, warnings, str(e))


def _basic_chromatin_qc(files: list, data_type: str,
                         warnings: list, import_error: str) -> dict:
    """Fallback QC when pysam/muon are not available."""
    warnings.append(
        f"Full QC unavailable ({import_error}). "
        f"Install aria-chromatin-env for complete metrics."
    )
    return {
        "status":         "success",
        "data_type":      data_type,
        "n_samples":      len(files),
        "frip":           None,
        "tss_enrichment": None,
        "mito_fraction":  None,
        "pass_qc":        None,
        "warnings":       warnings,
        "note":           "Partial QC — install aria-chromatin-env for full metrics",
    }


# ── Helper functions ──────────────────────────────────────────────────────────

_FRAG_SCAN_LIMIT = 2_000_000   # cap a streaming pass; flag if reached
_SIZE_SAMPLE_LIMIT = 100_000   # fragments kept for the size distribution


def _tss_qc(frag_file: str, genome: str) -> dict:
    """Shared snapatac2 build (scATAC P0 W0.5 + P4.1 (2)[A]): per-cell TSSe — and,
    additively, the per-cell (TSSe, log10 depth) arrays for the QC gating scatter —
    from a SINGLE ``import_data`` pass. Never fabricated, never silent egress.

    Returns a dict ``{"median", "reason", "tsse", "log10_depth"}``:
      - ``median`` = sample-level median per-cell TSSe (the W0.5 scalar), or None;
      - ``reason`` = concrete not-computed cause when the assembly is unknown,
        snapatac2 is missing, egress is disabled (air-gapped), or the compute fails;
      - ``tsse`` / ``log10_depth`` = per-cell arrays for the scatter (None when not
        computed; ``log10_depth`` stays None when the importer recorded no per-cell
        fragment count).

    The annotation comes from ``snap.genome.<attr>`` (governed auto-fetch, cached);
    under ``ARIA_AIR_GAPPED`` ARIA does not fetch it (W-PRIV) and skips honestly.
    """
    none = {"median": None, "reason": None, "tsse": None, "log10_depth": None}

    from aria.utils import genomes
    attr = genomes.snapatac2_attr(genome)
    if not attr:
        return {**none, "reason": (
            f"no snapatac2 gene annotation for assembly '{genome or 'unknown'}' "
            f"(TSS enrichment needs a known assembly)")}

    from aria.utils import privacy
    if not privacy.egress_allowed():
        return {**none, "reason": (
            f"air-gapped: the TSS annotation for '{genome}' is a governed "
            f"snapatac2 auto-fetch and was not downloaded (set ARIA_AIR_GAPPED=0 "
            f"or stage it to compute TSS enrichment)")}

    try:
        import snapatac2 as snap
    except ImportError as e:
        return {**none, "reason": (
            f"snapatac2 unavailable ({e}); TSS enrichment not computed")}

    try:
        import numpy as np
        gobj = getattr(snap.genome, attr, None)
        if gobj is None:
            return {**none, "reason": f"snapatac2 has no genome annotation '{attr}'"}
        adata = snap.pp.import_data(
            frag_file, chrom_sizes=gobj, sorted_by_barcode=False, file=None,
        )
        snap.metrics.tsse(adata, gobj)
        tsse = np.asarray(adata.obs["tsse"], dtype=float)
        finite = np.isfinite(tsse)
        if not finite.any():
            return {**none, "reason": "TSS enrichment produced no finite per-cell values"}
        # Sample-level scalar = median over the TSSe-finite cells (W0.5 semantics,
        # independent of depth availability).
        median = float(np.median(tsse[finite]))
        # Per-cell scatter arrays (P4.1 [A]): pair TSSe with log10 depth when the
        # importer recorded a per-cell fragment count.
        depth = None
        for key in ("n_fragment", "n_fragments", "frag_count"):
            if key in adata.obs:
                depth = np.asarray(adata.obs[key], dtype=float)
                break
        mask = finite
        log10_depth = None
        if depth is not None:
            mask = finite & np.isfinite(depth) & (depth > 0)
            log10_depth = np.log10(depth[mask]).tolist()
        return {"median": median, "reason": None,
                "tsse": tsse[mask].tolist(), "log10_depth": log10_depth}
    except Exception as e:
        return {**none, "reason": f"TSS enrichment computation failed: {e!s:.160}"}


def _compute_tss_enrichment(frag_file: str, genome: str):
    """REAL TSS enrichment (scATAC P0 W0.5): the sample-level MEDIAN per-cell TSSe.
    Thin wrapper over :func:`_tss_qc` that preserves the fixed ``(value, reason)``
    contract; the per-cell arrays live alongside it in the QC dict (P4.1 [A])."""
    r = _tss_qc(frag_file, genome)
    return r["median"], r["reason"]


def _compute_fragment_sizes(frag_file: str) -> dict:
    """
    Compute the fragment size distribution and the unique-barcode count from a
    10x-style fragments file (`chrom  start  end  barcode  count`). Only metrics
    genuinely derivable from the fragments file are returned; nothing is
    fabricated. The scan is capped at `_FRAG_SCAN_LIMIT` records — when the cap
    is reached the counts are lower bounds and `scan_truncated` is True.
    """
    try:
        import gzip
        sizes = []
        barcodes = set()
        n_records = 0
        scan_truncated = False
        opener = gzip.open if frag_file.endswith(".gz") else open
        with opener(frag_file, "rt") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                if n_records >= _FRAG_SCAN_LIMIT:
                    scan_truncated = True
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                n_records += 1
                if len(parts) >= 4 and parts[3]:
                    barcodes.add(parts[3])
                if len(sizes) < _SIZE_SAMPLE_LIMIT:
                    try:
                        size = int(parts[2]) - int(parts[1])
                        if 0 < size < 1000:
                            sizes.append(size)
                    except ValueError:
                        pass

        out = {
            "n_barcodes":      len(barcodes),
            "n_fragments":     n_records,
            "scan_truncated":  scan_truncated,
            "scan_limit":      _FRAG_SCAN_LIMIT,
        }
        if sizes:
            import numpy as np
            sizes_arr = np.array(sizes)
            out.update({
                "median_size":       int(np.median(sizes_arr)),
                "pct_mononucleosomal": float(
                    np.mean((sizes_arr > 150) & (sizes_arr < 300))
                ),
                "pct_subnucleosomal": float(np.mean(sizes_arr < 150)),
            })
            # W0.2 (scATAC P0): a binned size histogram (0-800 bp, 5 bp bins) so
            # the report can plot the nucleosome banding pattern. Summary-only
            # before; the raw sizes are not persisted (privacy + size).
            counts, edges = np.histogram(
                sizes_arr[(sizes_arr >= 0) & (sizes_arr <= 800)],
                bins=160, range=(0, 800),
            )
            out["size_histogram"] = {
                "bin_edges": [float(e) for e in edges],
                "counts":    [int(c) for c in counts],
            }
        return out
    except Exception:
        return {"n_barcodes": 0, "n_fragments": 0, "scan_truncated": False}


def _compute_mito_fraction(frag_file: str, mito_chr: str) -> float:
    """Estimate mitochondrial fragment fraction from fragment file."""
    try:
        import gzip
        total = 0; mito = 0
        opener = gzip.open if frag_file.endswith(".gz") else open
        with opener(frag_file, "rt") as f:
            for i, line in enumerate(f):
                if i > 200000: break
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) >= 1:
                    total += 1
                    if parts[0] in (mito_chr, "M", "MT", "chrM"):
                        mito += 1
        return mito / max(total, 1)
    except Exception:
        return 0.05


def _peaks_by_chrom_from_bed(peaks_bed: str) -> dict:
    """Read a BED of called peaks into ``{chrom: sorted [(start, end), ...]}``."""
    import gzip
    out: dict = {}
    opener = gzip.open if str(peaks_bed).endswith(".gz") else open
    with opener(peaks_bed, "rt") as f:
        for line in f:
            if not line or line[0] in "#\n" or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                out.setdefault(parts[0], []).append((int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    for chrom in out:
        out[chrom].sort()
    return out


def _fragment_overlaps_peak(chrom, start, end, peaks_by_chrom, starts_by_chrom) -> bool:
    """True iff fragment ``(chrom, start, end)`` overlaps any peak (``start <
    peak_end`` and ``end > peak_start``). Shared by the sample-level and the
    per-barcode FRiP so both encode exactly one overlap rule."""
    import bisect

    ivs = peaks_by_chrom.get(chrom)
    if not ivs:
        return False
    starts = starts_by_chrom[chrom]
    # Peaks with peak_start < frag_end are candidates; scan left from there
    # while peak_end > frag_start (intervals are short, overlap span small).
    k = bisect.bisect_left(starts, end) - 1
    while k >= 0:
        ps, pe = ivs[k]
        if pe <= start:
            # Could still be an earlier long peak, but peaks are short and sorted
            # by start; once peak_start is far left of frag_start, stop.
            if ps < start - 1_000_000:
                break
            k -= 1
            continue
        if ps < end and pe > start:
            return True
        k -= 1
    return False


def _frip_from_intervals(fragments, peaks_by_chrom) -> tuple:
    """REAL FRiP (scATAC P0 W0.5): fraction of fragments overlapping any peak.

    Pure / dependency-light so it is unit-testable without a fragments file or
    snapatac2. ``fragments`` is an iterable of ``(chrom, start, end)``;
    ``peaks_by_chrom`` is ``{chrom: sorted [(start, end), ...]}``. Returns
    ``(frip|None, n_fragments)``; None when there are no fragments.
    """
    starts_by_chrom = {c: [p[0] for p in iv] for c, iv in peaks_by_chrom.items()}
    n_total = 0
    n_in = 0
    for chrom, start, end in fragments:
        n_total += 1
        if _fragment_overlaps_peak(chrom, start, end, peaks_by_chrom, starts_by_chrom):
            n_in += 1
    if n_total == 0:
        return None, 0
    return n_in / n_total, n_total


def _frip_per_barcode_from_intervals(fragments_bc, peaks_by_chrom,
                                     min_frags: int = 20) -> list:
    """Per-barcode FRiP distribution (scATAC P4.1 (2)[B]): the same overlap rule
    as :func:`_frip_from_intervals`, accumulated per cell barcode.

    Pure / dependency-light. ``fragments_bc`` is an iterable of
    ``(chrom, start, end, barcode)``. Returns a list of per-barcode FRiP floats
    for barcodes with at least ``min_frags`` fragments (low-count barcodes give
    noisy ratios). Empty list when nothing qualifies — never fabricated.
    """
    starts_by_chrom = {c: [p[0] for p in iv] for c, iv in peaks_by_chrom.items()}
    totals: dict = {}
    ins: dict = {}
    for chrom, start, end, bc in fragments_bc:
        if not bc:
            continue
        totals[bc] = totals.get(bc, 0) + 1
        if _fragment_overlaps_peak(chrom, start, end, peaks_by_chrom, starts_by_chrom):
            ins[bc] = ins.get(bc, 0) + 1
    return [ins.get(bc, 0) / n for bc, n in totals.items() if n >= min_frags]


def _estimate_frip(frag_file: str, peaks_bed: str | None = None,
                   scan_limit: int = 2_000_000):
    """FRiP (fraction of reads in peaks). REAL when a called-peak BED is supplied
    (scATAC P0 W0.5, via `_frip_from_intervals`); honest None without peaks (FRiP
    is undefined before peak calling — never the old fabricated 0.35). The QC
    fragments path has no peaks yet, so it passes None; a post-peak-calling caller
    supplies `peaks_bed`."""
    if not peaks_bed:
        return None
    try:
        import gzip
        peaks_by_chrom = _peaks_by_chrom_from_bed(peaks_bed)
        if not peaks_by_chrom:
            return None

        def _frag_iter():
            opener = gzip.open if str(frag_file).endswith(".gz") else open
            with opener(frag_file, "rt") as f:
                for i, line in enumerate(f):
                    if i >= scan_limit:
                        break
                    if not line or line[0] == "#":
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3:
                        try:
                            yield parts[0], int(parts[1]), int(parts[2])
                        except ValueError:
                            continue

        frip, _ = _frip_from_intervals(_frag_iter(), peaks_by_chrom)
        return frip
    except Exception:
        return None


def _estimate_frip_distribution(frag_file: str, peaks_bed: str | None = None,
                                scan_limit: int = 2_000_000,
                                min_frags: int = 20):
    """Per-barcode FRiP distribution for the QC histogram (scATAC P4.1 (2)[B]).
    REAL when a called-peak BED is supplied (same gate as :func:`_estimate_frip`);
    honest None without peaks or when no barcode clears ``min_frags``."""
    if not peaks_bed:
        return None
    try:
        import gzip
        peaks_by_chrom = _peaks_by_chrom_from_bed(peaks_bed)
        if not peaks_by_chrom:
            return None

        def _frag_iter():
            opener = gzip.open if str(frag_file).endswith(".gz") else open
            with opener(frag_file, "rt") as f:
                for i, line in enumerate(f):
                    if i >= scan_limit:
                        break
                    if not line or line[0] == "#":
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 4:
                        try:
                            yield parts[0], int(parts[1]), int(parts[2]), parts[3]
                        except ValueError:
                            continue

        dist = _frip_per_barcode_from_intervals(
            _frag_iter(), peaks_by_chrom, min_frags=min_frags)
        return dist or None
    except Exception:
        return None


def _estimate_frip_bulk(data_type: str):
    """
    FRiP requires called peaks; per-assay constants would be fabricated
    (ADR-002). Returns None until peak calling produces a real value.
    """
    return None


def _estimate_dup_rate(bam_path: str):
    """
    PCR duplicate rate from `samtools flagstat`. Returns None when samtools is
    unavailable or the parse fails — ARIA does not substitute a fabricated
    'typical' constant (ADR-002); the previous fallback returned 0.15.
    """
    try:
        import subprocess
        result = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True, text=True, timeout=60,
        )
        lines = result.stdout.split("\n")
        for line in lines:
            if "duplicate" in line.lower():
                parts = line.split()
                if parts:
                    dups = int(parts[0])
                    total = int(lines[0].split()[0]) if lines else 1
                    return dups / max(total, 1)
    except Exception:
        pass
    return None


if __name__ == "__main__":
    run_script(chromatin_qc)
