"""
ARIA Chromatin Peak Calling Script
------------------------------------
Calls peaks using MACS3 for all chromatin modalities.
Executed inside aria-chromatin-env by EnvironmentManager.

Assay-specific parameters (pre-set by ChromatinAgent):
  ATAC (bulk/single):   --nomodel --extsize 200 --keep-dup all
  ChIP (TF):            narrow peaks, with input control if available
  ChIP (histone):       --broad, H3K4me1/H3K27ac enhancers
  CUT&RUN / CUT&TAG:   --nomodel --nolambda (very low background)

Input params:
    data_type:    str   — modality
    files:        list  — BAM or fragment files
    genome:       str   — assembly
    macs3_params: dict  — assay-specific MACS3 flags
    output_dir:   str   — where to write peaks
    control_files: list (optional) — input/IgG for ChIP

Output:
    {
      "status":               "success",
      "n_peaks":              int,
      "peaks_path":           str,   — path to .narrowPeak or .broadPeak
      "consensus_peaks_path": str,   — merged peaks across replicates
      "frip":                 float|None, — actual FRiP after peak calling
      "warnings":             [str]
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from pathlib import Path
from aria.scripts.chromatin_peak_counts import (
    file_record, load_valid_stage_result, stage_manifest_path, write_stage_manifest,
    _save_stage_result, _open_alignment, _interval_index, _overlaps, _primary_mapped,
)


def chromatin_peaks(params: dict) -> dict:
    from pathlib import Path

    data_type     = params["data_type"]
    files         = params.get("files", [])
    genome        = params.get("genome", "hg38")
    macs3_params  = dict(params.get("macs3_params", {}))
    for key in ("qvalue", "broad", "format"):
        if key in params:
            macs3_params[key] = params[key]
    if params.get("reference_fasta"):
        macs3_params["reference_fasta"] = params["reference_fasta"]
    output_dir    = params.get("output_dir", "/tmp/aria_peaks")
    control_files = params.get("control_files", [])
    assay_class   = str(params.get("assay_class") or ("histone" if macs3_params.get("broad") else "tf"))
    requested_caller = str(params.get("peak_caller") or "macs3").strip().lower()
    if requested_caller == "auto":
        if assay_class == "broad_domain":
            peak_caller = "broad_domain"
        elif data_type in {"CUT_AND_RUN", "CUT_AND_TAG"}:
            peak_caller = "seacr"
        else:
            peak_caller = "macs3"
    else:
        peak_caller = requested_caller
    if peak_caller not in {"macs3", "seacr", "broad_domain"}:
        return {
            "status": "error",
            "error_type": "InvalidPeakCaller",
            "details": f"Unsupported peak_caller={peak_caller!r}; expected macs3, seacr, broad_domain, or auto.",
        }
    seacr_mode = str(params.get("seacr_mode") or ("relaxed" if assay_class in {"histone", "broad"} else "stringent"))
    seacr_norm = str(params.get("seacr_norm") or "norm")
    seacr_top_quantile = float(params.get("seacr_top_quantile", 0.01))
    target_scale_factors = params.get("target_scale_factors")
    control_scale_factors = params.get("control_scale_factors")
    run_idr = bool(
        params.get(
            "run_idr",
            data_type in {"ChIP", "CUT_AND_RUN", "CUT_AND_TAG"}
            and assay_class == "tf"
            and not macs3_params.get("broad", False),
        )
    )
    idr_threshold = float(params.get("idr_threshold", 0.05))
    write_bigwig = bool(params.get("write_bigwig", False))
    run_replicate_peak_calling = params.get(
        "run_replicate_peak_calling",
        data_type in {"bulk_ATAC", "scATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"},
    )
    min_replicate_support = int(params.get("min_replicates_for_consensus",
                                           params.get("min_replicate_peak_support", 2)))

    warnings = []
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Validate input files
    valid_files = [f for f in files if Path(f).exists()]
    if not valid_files:
        return {
            "status":     "error",
            "error_type": "FileNotFound",
            "details":    f"No valid input files for peak calling.",
        }

    if len(valid_files) != len(files):
        return {"status": "error", "error_type": "FileNotFound",
                "details": "Some peak-calling inputs are missing."}
    blacklist = params.get("blacklist_bed")
    reference = params.get("reference_fasta")
    inputs = [file_record(f) for f in valid_files]
    for extra in [blacklist, reference, *control_files]:
        if extra and Path(extra).is_file():
            inputs.append(file_record(extra))
    parameters = {"bam_files": list(files), "genome": genome,
                  "qvalue": macs3_params.get("qvalue", 0.05),
                  "broad": macs3_params.get("broad", False),
                  "format": macs3_params.get("format", "BAMPE"),
                  "blacklist_bed": blacklist,
                  "min_replicates_for_consensus": min_replicate_support,
                  "macs3_params": macs3_params, "control_files": control_files,
                  "run_replicate_peak_calling": bool(run_replicate_peak_calling),
                  "data_type": data_type}
    if peak_caller != "macs3":
        parameters["peak_caller"] = peak_caller
        parameters["seacr_mode"] = seacr_mode
        parameters["seacr_norm"] = seacr_norm
        parameters["seacr_top_quantile"] = seacr_top_quantile
    if run_idr:
        parameters["run_idr"] = True
        parameters["idr_threshold"] = idr_threshold
    if target_scale_factors:
        parameters["target_scale_factors"] = list(target_scale_factors)
    if write_bigwig:
        parameters["write_bigwig"] = True
    manifest = stage_manifest_path(output_dir, "chromatin_peaks")
    cached = load_valid_stage_result(
        manifest, expected_inputs=inputs, expected_parameters=parameters)
    if (cached and Path(cached.get("peaks_path", "")).is_file()
            and (not cached.get("consensus_peaks_path")
                 or Path(cached["consensus_peaks_path"]).is_file())):
        return {**cached, "reused_from_stage_manifest": True}
    Path(manifest).unlink(missing_ok=True)

    genome_size = _get_genome_size(genome)
    sample_name = Path(valid_files[0]).stem.replace(".bam", "")
    out_dir = Path(output_dir)

    valid_ctrl = []
    if control_files:
        valid_ctrl = [f for f in control_files if Path(f).exists()]
        if not valid_ctrl:
            warnings.append(
                "Control files specified but not found. "
                "Running without input control."
            )

    cmd: list[str] = []
    effective_format = macs3_params.get("format", "BAMPE")
    seacr_diagnostics = None

    if peak_caller == "seacr":
        from aria.utils.seacr_idr import call_seacr_peaks

        peaks_ext = f".{seacr_mode}.bed"
        peaks_file = out_dir / f"{sample_name}_peaks{peaks_ext}"
        peaks_file.unlink(missing_ok=True)
        seacr_res = call_seacr_peaks(
            valid_files,
            peaks_file,
            control_inputs=valid_ctrl or None,
            mode=seacr_mode,
            norm=seacr_norm,
            top_quantile=seacr_top_quantile,
            target_scale_factors=target_scale_factors,
            control_scale_factors=control_scale_factors,
            keep_duplicates=bool(params.get("keep_duplicates", True)),
            reference_fasta=reference,
            work_dir=out_dir / "_seacr_pooled",
        )
        if seacr_res.get("status") == "error":
            return seacr_res
        warnings.extend(seacr_res.get("warnings") or [])
        seacr_diagnostics = seacr_res.get("diagnostics")
        effective_format = "SEACR_BEDGRAPH"
    elif peak_caller == "broad_domain":
        from aria.utils.seacr_idr import call_broad_domains

        peaks_ext = ".broadDomain.bed"
        peaks_file = out_dir / f"{sample_name}_peaks{peaks_ext}"
        peaks_file.unlink(missing_ok=True)
        dom_res = call_broad_domains(
            valid_files,
            peaks_file,
            control_bams=valid_ctrl or None,
            bin_size=int(params.get("domain_bin_size", 10_000)),
            min_log2_enrichment=float(params.get("domain_min_log2fc", 0.585)),
            max_gap_bp=int(params.get("domain_max_gap_bp", 30_000)),
            min_domain_bp=int(params.get("domain_min_bp", 20_000)),
            ip_scale_factors=target_scale_factors,
            control_scale_factors=control_scale_factors,
            reference_fasta=reference,
        )
        if dom_res.get("status") == "error":
            return dom_res
        effective_format = "BROAD_DOMAIN_BINS"
    else:
        effective_params, effective_format = _effective_macs3_params(
            valid_files, macs3_params, data_type=data_type
        )
        if effective_format == "BAM_SE_TN5_SHIFT":
            warnings.append(
                "Single-end alignments detected: MACS3 uses -f BAM --nomodel "
                "--shift -100 --extsize 200 (Tn5-centered SE ATAC). "
                "Pooled inputs containing SE reads use this treatment for all reads.")
        cmd = _build_macs3_cmd(
            valid_files, valid_ctrl, genome_size, sample_name, output_dir,
            effective_params, data_type=data_type,
        )
        broad = macs3_params.get("broad", False)
        peaks_ext = ".broadPeak" if broad else ".narrowPeak"
        peaks_file = out_dir / f"{sample_name}_peaks{peaks_ext}"
        peaks_file.unlink(missing_ok=True)

        # A failed rerun must never reuse stale peak files after invalidation.
        peak_run = _run_macs3(cmd, timeout=7200)
        if peak_run["status"] == "error" and data_type == "ChIP" and not effective_params.get("nomodel", False):
            # Sparse/single-contig ChIP BAMs can fail MACS3 paired-peak model building;
            # fall back transparently to --nomodel --extsize 200 with an explicit warning.
            fallback_params = {**effective_params, "nomodel": True, "extsize": int(effective_params.get("extsize", 200))}
            warnings.append(
                "MACS3 paired-peak model building failed on input BAM(s); "
                f"retried with --nomodel --extsize {fallback_params['extsize']}."
            )
            cmd = _build_macs3_cmd(
                valid_files, valid_ctrl, genome_size, sample_name, output_dir,
                fallback_params, data_type=data_type,
            )
            peak_run = _run_macs3(cmd, timeout=7200)
        if peak_run["status"] == "error":
            return peak_run
        if peak_run.get("warning"):
            warnings.append(str(peak_run["warning"]))

    if not peaks_file.is_file():
        return {"status": "error", "error_type": "NoPeaksFile",
                "details": f"Expected {peaks_file} not found after {peak_caller}."}

    # ── Count peaks and compute basic stats ──────────────────────────────
    n_blacklist_filtered = _filter_blacklist(str(peaks_file), blacklist)
    n_peaks = _count_peaks(str(peaks_file))

    if n_peaks < 1000:
        warnings.append(
            f"Only {n_peaks} peaks called. "
            f"Expected >10,000 for typical ATAC-seq. "
            f"Check library quality and alignment rate."
        )
    elif n_peaks > 500000:
        warnings.append(
            f"{n_peaks:,} peaks is unusually high. "
            f"Consider increasing q-value threshold."
        )

    # ── Replicate reproducibility policy (C5 / ADR-081) ──────────────────
    peak_repro = _replicate_reproducibility(
        valid_files=valid_files,
        control_files=valid_ctrl,
        genome_size=genome_size,
        macs3_params=macs3_params,
        output_dir=output_dir,
        peaks_ext=peaks_ext,
        pooled_peaks_path=str(peaks_file),
        min_support=min_replicate_support,
        run_replicate_peak_calling=bool(run_replicate_peak_calling),
        blacklist_bed=blacklist,
        peak_caller=peak_caller,
        data_type=data_type,
        seacr_mode=seacr_mode,
        seacr_norm=seacr_norm,
        seacr_top_quantile=seacr_top_quantile,
        target_scale_factors=target_scale_factors,
        control_scale_factors=control_scale_factors,
        run_idr=run_idr,
        idr_threshold=idr_threshold,
        reference_fasta=reference,
        extra_params=params,
    )
    warnings.extend(peak_repro.get("warnings") or [])
    consensus_path = peak_repro.get("consensus_peaks_path")

    n_blacklist_filtered += peak_repro.get("n_blacklist_filtered", 0)
    if consensus_path and Path(consensus_path).is_file():
        n_blacklist_filtered += _filter_blacklist(consensus_path, blacklist)
        overlap = peak_repro.get("overlap") or {}
        if overlap.get("ran"):
            overlap["n_reproducible_regions"] = _count_peaks(consensus_path)
            candidates = overlap.get("n_candidate_regions", 0)
            overlap["fraction_reproducible"] = (
                round(overlap["n_reproducible_regions"] / candidates, 4) if candidates else None)
        if _count_peaks(consensus_path) == 0:
            peak_repro["status"] = "unverified"
            peak_repro["reason"] = "no_reproducible_regions_after_blacklist"
    frip_per_sample = {}
    for path in valid_files:
        value = _compute_frip(path, str(peaks_file), reference_fasta=reference)
        if value is not None:
            frip_per_sample[Path(path).stem] = value
    frip = sum(frip_per_sample.values()) / len(frip_per_sample) if frip_per_sample else None
    if frip is None and any(str(p).lower().endswith((".bam", ".cram")) for p in valid_files):
        warnings.append("FRiP unavailable: no mapped nonduplicate primary reads or unreadable alignments.")

    bigwig_tracks: dict[str, str] = {}
    if write_bigwig:
        from aria.utils.seacr_idr import write_normalized_bigwig

        bw_dir = out_dir / "bigwig"
        for idx, path in enumerate(valid_files):
            if not str(path).lower().endswith((".bam", ".cram")):
                continue
            sid = Path(path).stem
            sf = (
                float(target_scale_factors[idx])
                if (target_scale_factors and idx < len(target_scale_factors))
                else None
            )
            try:
                bw_info = write_normalized_bigwig(
                    path,
                    bw_dir / f"{sid}.bw",
                    normalization=str(params.get("bigwig_normalization") or "CPM"),
                    spikein_scale_factor=sf,
                    reference_fasta=reference,
                )
                bigwig_tracks[sid] = bw_info["bigwig_path"]
            except Exception as exc:
                warnings.append(f"BigWig export skipped for {sid}: {str(exc)[:140]}")

    result = {
        "status": "success", "data_type": data_type,
        "assay_class": assay_class,
        "peak_caller": peak_caller,
        "validation_level": "production",
        "n_peaks": int(n_peaks), "peaks_path": str(peaks_file),
        "consensus_peaks_path": consensus_path,
        "n_consensus_peaks": _count_peaks(consensus_path) if consensus_path else None,
        "n_blacklist_filtered": n_blacklist_filtered,
        "frip": frip, "frip_per_sample": frip_per_sample,
        "genome": genome,
        "macs3_cmd": " ".join(cmd) if cmd else None,
        "seacr_diagnostics": seacr_diagnostics,
        "bigwig_tracks": bigwig_tracks or None,
        "effective_format": effective_format,
        "peak_calling_strategy": peak_repro.get("strategy"),
        "peak_reproducibility": peak_repro, "warnings": warnings,
        "reused_from_stage_manifest": False,
    }
    _save_stage_result(manifest, "chromatin_peaks", inputs, parameters, result)
    return result


# ── Helper functions ──────────────────────────────────────────────────────────

def _get_genome_size(genome: str) -> str:
    """Return MACS3 genome size shorthand."""
    SIZE_MAP = {
        "hg38":    "hs", "hg19":    "hs", "GRCh38": "hs",
        "mm10":    "mm", "mm39":    "mm", "GRCm38": "mm",
        "dm6":     "dm",
        "ce11":    "ce",
        "danRer11": "2e9",
        "sacCer3":  "1.2e7",
    }
    return SIZE_MAP.get(genome, "hs")


def _build_macs3_cmd(input_files: list[str],
                     control_files: list[str],
                     genome_size: str,
                     sample_name: str,
                     output_dir: str,
                     macs3_params: dict,
                     data_type: str = "bulk_ATAC") -> list[str]:
    """Build the exact MACS3 command used for a pooled or per-replicate call."""
    macs3_params, _ = _effective_macs3_params(input_files, macs3_params, data_type=data_type)
    cmd = ["macs3", "callpeak", "-t"] + list(input_files)
    if control_files:
        cmd += ["-c"] + list(control_files)
    cmd += ["-f", macs3_params.get("format", "BAMPE")]
    cmd += ["-g", genome_size, "-q", str(macs3_params.get("qvalue", 0.05))]
    cmd += ["-n", sample_name, "--outdir", output_dir]
    if macs3_params.get("nomodel", False):
        cmd.append("--nomodel")
        cmd += ["--extsize", str(macs3_params.get("extsize", 200))]
    if "shift" in macs3_params:
        cmd += ["--shift", str(macs3_params["shift"])]
    if macs3_params.get("nolambda", False):
        cmd.append("--nolambda")
    if macs3_params.get("broad", False):
        cmd.append("--broad")
    cmd += ["--keep-dup", str(macs3_params.get("keep_dup", "1"))]
    cmd += ["-B", "--SPMR"]
    return cmd


def _run_macs3(cmd: list[str], timeout: int) -> dict:
    import subprocess

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {
            "status":     "error",
            "error_type": "Timeout",
            "details":    "MACS3 exceeded time limit.",
        }
    except FileNotFoundError:
        return {
            "status":     "error",
            "error_type": "MACS3NotFound",
            "details":    "MACS3 not installed.",
        }

    if result.returncode != 0:
        stderr = result.stderr or ""
        if "Traceback" in stderr or "Error" in stderr:
            return {
                "status":     "error",
                "error_type": "MACS3Failed",
                "details":    stderr[-1000:],
            }
        return {
            "status": "error",
            "error_type": "MACS3NonZero",
            "details": stderr[-1000:],
        }
    return {"status": "success"}


def _sample_label(path: str, index: int | None = None) -> str:
    from pathlib import Path
    import re

    name = Path(path).name
    for suffix in (".bam", ".cram", ".sam", ".fragments.tsv.gz",
                   ".tsv.gz", ".bed.gz", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-") or "sample"
    return f"rep{index + 1}_{label}" if index is not None else label


def _find_peak_file(output_dir: str, sample_name: str, peaks_ext: str,
                    allow_fallback: bool = True) -> str | None:
    from pathlib import Path

    out_dir = Path(output_dir)
    expected = out_dir / f"{sample_name}_peaks{peaks_ext}"
    if expected.exists():
        return str(expected)
    if not allow_fallback:
        return None
    candidates = sorted(out_dir.glob(f"*{peaks_ext}"))
    return str(candidates[0]) if candidates else None


def _replicate_reproducibility(valid_files: list[str],
                               genome_size: str,
                               macs3_params: dict,
                               output_dir: str,
                               peaks_ext: str,
                               pooled_peaks_path: str,
                               min_support: int,
                               run_replicate_peak_calling: bool,
                               blacklist_bed: str | None = None,
                               control_files: list[str] | None = None,
                               peak_caller: str = "macs3",
                               data_type: str = "bulk_ATAC",
                               seacr_mode: str = "stringent",
                               seacr_norm: str = "norm",
                               seacr_top_quantile: float = 0.01,
                               target_scale_factors: list[float] | None = None,
                               control_scale_factors: list[float] | None = None,
                               run_idr: bool = False,
                               idr_threshold: float = 0.05,
                               reference_fasta: str | None = None,
                               extra_params: dict | None = None) -> dict:
    from pathlib import Path

    control_files = list(control_files or [])
    extra_params = extra_params or {}
    strategy = f"pooled_{peak_caller}"
    if len(valid_files) < 2:
        return {
            "status": "not_applicable",
            "strategy": strategy,
            "reason": "single_input_file",
            "pooled_peak_calling": True,
            "pooled_peaks_path": pooled_peaks_path,
            "replicate_peak_calling": {"ran": False, "reason": "single_input_file"},
            "overlap": {"ran": False, "reason": "single_input_file"},
            "idr": {"ran": False, "reason": "requires_replicate_peak_sets"},
            "warnings": [],
        }

    repro = {
        "status": "unverified",
        "strategy": (
            f"pooled_{peak_caller}_with_copula_idr"
            if run_idr
            else f"pooled_{peak_caller}_with_replicate_overlap_qc"
        ),
        "pooled_peak_calling": True,
        "pooled_peaks_path": pooled_peaks_path,
        "replicate_peak_calling": {
            "ran": False,
            "n_inputs": len(valid_files),
            "peak_files": [],
            "failed": [],
        },
        "overlap": {"ran": False, "reason": "replicate_peak_sets_missing"},
        "idr": {
            "ran": False,
            "reason": "idr_not_run_overlap_reproducibility_policy",
        },
        "warnings": [],
    }
    if not run_replicate_peak_calling:
        repro["reason"] = "replicate_peak_calling_disabled"
        repro["warnings"].append(
            f"Peak calling used a pooled {peak_caller.upper()} call across multiple inputs; "
            "per-replicate overlap/IDR reproducibility was not run."
        )
        return repro

    rep_dir = Path(output_dir) / "replicate_peaks"
    rep_dir.mkdir(parents=True, exist_ok=True)
    peak_files: list[str] = []
    failed: list[dict] = []
    for i, input_file in enumerate(valid_files):
        label = _sample_label(input_file, i)
        rep_out = rep_dir / f"{label}_peaks{peaks_ext}"
        rep_out.unlink(missing_ok=True)
        sf_i = (
            [float(target_scale_factors[i])]
            if (target_scale_factors and i < len(target_scale_factors))
            else None
        )
        if peak_caller == "seacr":
            from aria.utils.seacr_idr import call_seacr_peaks

            run = call_seacr_peaks(
                [input_file],
                rep_out,
                control_inputs=control_files or None,
                mode=seacr_mode,
                norm=seacr_norm,
                top_quantile=seacr_top_quantile,
                target_scale_factors=sf_i,
                control_scale_factors=control_scale_factors,
                keep_duplicates=bool(extra_params.get("keep_duplicates", True)),
                reference_fasta=reference_fasta,
                work_dir=rep_dir / f"_seacr_{label}",
            )
            peak_file = str(rep_out) if run.get("status") == "success" and rep_out.is_file() else None
        elif peak_caller == "broad_domain":
            from aria.utils.seacr_idr import call_broad_domains

            run = call_broad_domains(
                [input_file],
                rep_out,
                control_bams=control_files or None,
                bin_size=int(extra_params.get("domain_bin_size", 10_000)),
                min_log2_enrichment=float(extra_params.get("domain_min_log2fc", 0.585)),
                max_gap_bp=int(extra_params.get("domain_max_gap_bp", 30_000)),
                min_domain_bp=int(extra_params.get("domain_min_bp", 20_000)),
                ip_scale_factors=sf_i,
                control_scale_factors=control_scale_factors,
                reference_fasta=reference_fasta,
            )
            peak_file = str(rep_out) if run.get("status") == "success" and rep_out.is_file() else None
        else:
            # For IDR, relax MACS3 q-value slightly (0.10) on individual replicates if strict
            rep_macs3 = dict(macs3_params)
            if run_idr and float(rep_macs3.get("qvalue", 0.05)) < 0.05:
                rep_macs3["qvalue"] = 0.05
            cmd = _build_macs3_cmd(
                [input_file], control_files, genome_size, label, str(rep_dir),
                rep_macs3, data_type=data_type,
            )
            run = _run_macs3(cmd, timeout=7200)
            if run.get("status") != "success" and data_type == "ChIP" and not rep_macs3.get("nomodel", False):
                fb = {**rep_macs3, "nomodel": True, "extsize": int(rep_macs3.get("extsize", 200))}
                cmd = _build_macs3_cmd(
                    [input_file], control_files, genome_size, label, str(rep_dir),
                    fb, data_type=data_type,
                )
                run = _run_macs3(cmd, timeout=7200)
            peak_file = (
                _find_peak_file(str(rep_dir), label, peaks_ext, allow_fallback=False)
                if run.get("status") == "success" else None
            )

        if run.get("status") != "success":
            failed.append({
                "input_file": input_file,
                "error_type": run.get("error_type", "PeakCallerFailed"),
                "details": run.get("details"),
            })
            continue
        if peak_file:
            removed = _filter_blacklist(peak_file, blacklist_bed)
            repro["n_blacklist_filtered"] = repro.get("n_blacklist_filtered", 0) + removed
            peak_files.append(peak_file)
        else:
            failed.append({
                "input_file": input_file,
                "error_type": "NoReplicatePeaksFile",
                "details": f"Expected {label}_peaks{peaks_ext}",
            })

    repro["replicate_peak_calling"] = {
        "ran": bool(peak_files),
        "n_inputs": len(valid_files),
        "n_peak_files": len(peak_files),
        "peak_files": peak_files,
        "failed": failed,
    }
    if failed:
        repro["warnings"].append(
            f"Some per-replicate {peak_caller.upper()} peak calls failed; reproducibility "
            "metrics are partial."
        )
    if len(peak_files) < 2:
        repro["reason"] = "fewer_than_two_replicate_peak_sets"
        repro["warnings"].append(
            "Per-replicate peak reproducibility could not be assessed because "
            "fewer than two replicate peak files were available."
        )
        return repro

    consensus = str(Path(output_dir) / f"reproducible_consensus{peaks_ext}")
    overlap = _write_overlap_consensus(
        peak_files, consensus, min_support=max(2, min_support),
    )
    repro["overlap"] = overlap

    if run_idr:
        from aria.utils.seacr_idr import compute_idr_reproducible_peaks

        idr_out = str(Path(output_dir) / f"idr_reproducible_consensus{peaks_ext}")
        idr_res = compute_idr_reproducible_peaks(
            peak_files, idr_out, idr_threshold=idr_threshold
        )
        repro["idr"] = idr_res
        if idr_res.get("ran") and idr_res.get("n_reproducible_peaks", 0) > 0:
            repro["consensus_peaks_path"] = idr_out
            repro["status"] = "verified" if not failed else "partial"
            return repro
        repro["warnings"].append(
            f"Copula IDR did not yield reproducible peaks ({idr_res.get('reason') or '0 passing'}); "
            "falling back to replicate overlap consensus."
        )

    repro["consensus_peaks_path"] = consensus if overlap.get("ran") else None
    if overlap.get("ran") and overlap.get("n_reproducible_regions", 0) > 0:
        repro["status"] = "verified" if not failed else "partial"
    else:
        repro["status"] = "unverified"
        repro["reason"] = "no_reproducible_overlap_regions"
    if not run_idr:
        repro["warnings"].append(
            f"IDR was not run; ARIA used an overlap-support reproducible peak "
            f"policy across per-replicate {peak_caller.upper()} calls."
        )
    return repro


def _read_peak_intervals(peaks_file: str) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []
    try:
        with open(peaks_file, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start = int(float(parts[1]))
                    end = int(float(parts[2]))
                except ValueError:
                    continue
                if end > start:
                    intervals.append((parts[0], start, end))
    except OSError:
        return []
    return intervals


def _merge_intervals_with_support(
    peak_files: list[str],
) -> list[tuple[str, int, int, set[int]]]:
    tagged: list[tuple[str, int, int, int]] = []
    for idx, path in enumerate(peak_files):
        tagged.extend((chrom, start, end, idx)
                      for chrom, start, end in _read_peak_intervals(path))
    tagged.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    merged: list[tuple[str, int, int, set[int]]] = []
    for chrom, start, end, idx in tagged:
        if not merged or merged[-1][0] != chrom or start > merged[-1][2]:
            merged.append((chrom, start, end, {idx}))
            continue
        prev_chrom, prev_start, prev_end, support = merged[-1]
        support.add(idx)
        merged[-1] = (prev_chrom, prev_start, max(prev_end, end), support)
    return merged


def _write_overlap_consensus(peak_files: list[str],
                             output_path: str,
                             min_support: int = 2) -> dict:
    from pathlib import Path

    merged = _merge_intervals_with_support(peak_files)
    kept = [m for m in merged if len(m[3]) >= min_support]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as out:
        for chrom, start, end, support in kept:
            out.write(
                f"{chrom}\t{start}\t{end}\t"
                f"support_{len(support)}\t{len(support)}\n"
            )
    n_candidates = len(merged)
    return {
        "ran": True,
        "method": "overlap_support_consensus",
        "support_threshold": int(min_support),
        "n_replicate_peak_files": len(peak_files),
        "n_candidate_regions": n_candidates,
        "n_reproducible_regions": len(kept),
        "fraction_reproducible": (
            round(len(kept) / n_candidates, 4) if n_candidates else None
        ),
        "output_path": output_path,
    }


def _count_peaks(peaks_file: str) -> int:
    """Count number of peaks in a narrowPeak/broadPeak file."""
    try:
        with open(peaks_file) as f:
            return sum(1 for line in f
                       if line.strip() and not line.startswith("#"))
    except Exception:
        return 0


def _effective_macs3_params(files, params, data_type: str = "bulk_ATAC"):
    """Inspect at most 1,000 mapped primary reads per BAM before selecting BAMPE."""
    result = dict(params)
    fmt = result.get("format", "BAMPE")
    if fmt != "BAMPE" and data_type != "ChIP":
        return result, ("BAM_SE_TN5_SHIFT" if fmt == "BAM" and result.get("shift") == -100 else fmt)
    for path in files:
        try:
            mapped = paired = 0
            with _open_alignment(path, params.get("reference_fasta")) as bam:
                for read in bam.fetch(until_eof=True):
                    if not _primary_mapped(read):
                        continue
                    mapped += 1
                    paired += int(read.is_paired)
                    if mapped >= 1000:
                        break
            if mapped and paired and fmt == "BAM":
                result["format"] = "BAMPE"
                return result, "BAMPE"
            if mapped and not paired:
                if data_type == "ChIP":
                    result["format"] = "BAM"
                    return result, "BAM_SE_CHIP"
                result.update(format="BAM", nomodel=True, shift=-100, extsize=200)
                return result, "BAM_SE_TN5_SHIFT"
        except (ImportError, OSError, ValueError):
            continue
    return result, result.get("format", fmt)


def _filter_blacklist(peaks_path, blacklist_bed):
    """Exclude any overlap, retaining every nonblacklisted peak column."""
    if not blacklist_bed or not Path(blacklist_bed).is_file():
        return 0
    index = _interval_index(_read_peak_intervals(str(blacklist_bed)))
    kept = []
    removed = 0
    with open(peaks_path) as handle:
        for line in handle:
            cols = line.rstrip().split("\t")
            if line.startswith(("#", "track", "browser")) or len(cols) < 3:
                kept.append(line)
                continue
            if _overlaps(index, cols[0], int(cols[1]), int(cols[2])):
                removed += 1
            else:
                kept.append(line)
    if removed:
        Path(peaks_path).write_text("".join(kept))
    return removed


def _compute_frip(bam_file: str, peaks_file: str,
                  reference_fasta: str | None = None) -> float | None:
    """FRiP among mapped nonduplicate primary reads; no samtools dependency.

    Stream a filtered BAM into bedtools so its BED count uses exactly the same
    population as the denominator. CRAM and missing/failed bedtools use the
    interval-bisect result accumulated during the pysam pass.
    """
    import shutil
    import subprocess
    import tempfile
    try:
        import pysam
        intervals = _read_peak_intervals(peaks_file)
        if not intervals:
            return None
        index = _interval_index(intervals)
        bedtools = shutil.which("bedtools")
        use_bedtools = bedtools and not str(bam_file).lower().endswith(".cram")
        total = overlaps = 0
        with tempfile.TemporaryDirectory(prefix="aria_frip_") as tmp:
            filtered = str(Path(tmp) / "primary_nondup.bam")
            with _open_alignment(bam_file, reference_fasta) as bam:
                out = pysam.AlignmentFile(filtered, "wb", template=bam) if use_bedtools else None
                try:
                    for read in bam.fetch(until_eof=True):
                        if not _primary_mapped(read) or read.is_duplicate:
                            continue
                        total += 1
                        if read.reference_end is not None and _overlaps(
                                index, read.reference_name, read.reference_start, read.reference_end):
                            overlaps += 1
                        if out is not None:
                            out.write(read)
                finally:
                    if out is not None:
                        out.close()
            if not total:
                return None
            if use_bedtools:
                try:
                    with tempfile.TemporaryFile(mode="w+t") as output:
                        proc = subprocess.run(
                            [bedtools, "intersect", "-u", "-abam", filtered,
                             "-b", str(peaks_file), "-bed"],
                            stdout=output, stderr=subprocess.PIPE, text=True, timeout=7200)
                        if proc.returncode == 0:
                            output.seek(0)
                            overlaps = sum(1 for line in output if line.strip())
                except (OSError, subprocess.TimeoutExpired):
                    pass  # The exact same interval population was counted above.
        return max(0.0, min(1.0, overlaps / total))
    except (ImportError, OSError, ValueError):
        return None


if __name__ == "__main__":
    run_script(chromatin_peaks)
