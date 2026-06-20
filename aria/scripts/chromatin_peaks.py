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


def chromatin_peaks(params: dict) -> dict:
    from pathlib import Path

    data_type     = params["data_type"]
    files         = params.get("files", [])
    genome        = params.get("genome", "hg38")
    macs3_params  = params.get("macs3_params", {})
    output_dir    = params.get("output_dir", "/tmp/aria_peaks")
    control_files = params.get("control_files", [])
    run_replicate_peak_calling = params.get(
        "run_replicate_peak_calling",
        data_type in {"bulk_ATAC", "scATAC"},
    )
    min_replicate_support = int(params.get("min_replicate_peak_support", 2))

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

    # ── Build MACS3 command ───────────────────────────────────────────────
    genome_size = _get_genome_size(genome)
    sample_name = Path(valid_files[0]).stem.replace(".bam", "")

    valid_ctrl = []
    if control_files:
        valid_ctrl = [f for f in control_files if Path(f).exists()]
        if not valid_ctrl:
            warnings.append(
                "Control files specified but not found. "
                "Running without input control."
            )

    cmd = _build_macs3_cmd(
        valid_files, valid_ctrl, genome_size, sample_name, output_dir,
        macs3_params,
    )

    # ── Run MACS3 ────────────────────────────────────────────────────────
    peak_run = _run_macs3(cmd, timeout=7200)
    if peak_run["status"] == "error":
        if peak_run.get("error_type") == "MACS3NonZero":
            warnings.append(f"MACS3 non-zero exit: {peak_run.get('details', '')[-200:]}")
        else:
            return peak_run
    elif peak_run.get("warning"):
        warnings.append(str(peak_run["warning"]))

    if peak_run.get("error_type") == "MACS3NotFound":
        return {
            "status":     "error",
            "error_type": "MACS3NotFound",
            "details":    "MACS3 not installed. Activate aria-chromatin-env.",
        }

    # ── Find output peaks file ────────────────────────────────────────────
    out_dir    = Path(output_dir)
    broad      = macs3_params.get("broad", False)
    peaks_ext  = ".broadPeak" if broad else ".narrowPeak"
    peaks_file = out_dir / f"{sample_name}_peaks{peaks_ext}"

    if not peaks_file.exists():
        # Try alternate naming
        candidates = list(out_dir.glob(f"*{peaks_ext}"))
        if candidates:
            peaks_file = candidates[0]
        else:
            return {
                "status":     "error",
                "error_type": "NoPeaksFile",
                "details":    f"Expected {peaks_file} not found after MACS3.",
            }

    # ── Count peaks and compute basic stats ──────────────────────────────
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

    # ── Replicate reproducibility policy (C5) ────────────────────────────
    peak_repro = _replicate_reproducibility(
        valid_files=valid_files,
        genome_size=genome_size,
        macs3_params=macs3_params,
        output_dir=output_dir,
        peaks_ext=peaks_ext,
        pooled_peaks_path=str(peaks_file),
        min_support=min_replicate_support,
        run_replicate_peak_calling=bool(run_replicate_peak_calling),
    )
    warnings.extend(peak_repro.get("warnings") or [])
    consensus_path = peak_repro.get("consensus_peaks_path")

    # ── Compute actual FRiP ───────────────────────────────────────────────
    frip = _compute_frip(valid_files[0], str(peaks_file))
    if frip is None:
        warnings.append(
            "FRiP was not computed because samtools/bedtools counting failed "
            "or was unavailable."
        )

    return {
        "status":               "success",
        "data_type":            data_type,
        "n_peaks":              int(n_peaks),
        "peaks_path":           str(peaks_file),
        "consensus_peaks_path": consensus_path,
        "frip":                 (round(float(frip), 4)
                                  if frip is not None else None),
        "genome":               genome,
        "macs3_cmd":            " ".join(cmd),  # for reproducibility
        "peak_calling_strategy": peak_repro.get("strategy"),
        "peak_reproducibility":  peak_repro,
        "warnings":             warnings,
    }


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
                     macs3_params: dict) -> list[str]:
    """Build the exact MACS3 command used for a pooled or per-replicate call."""
    cmd = ["macs3", "callpeak", "-t"] + list(input_files)
    if control_files:
        cmd += ["-c"] + list(control_files)
    cmd += ["-f", macs3_params.get("format", "BAMPE")]
    cmd += ["-g", genome_size]
    cmd += ["-n", sample_name, "--outdir", output_dir]
    if macs3_params.get("nomodel", False):
        cmd.append("--nomodel")
        cmd += ["--extsize", str(macs3_params.get("extsize", 200))]
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
                               run_replicate_peak_calling: bool) -> dict:
    from pathlib import Path

    strategy = "pooled_macs3"
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
        "strategy": "pooled_macs3_with_replicate_overlap_qc",
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
            "Peak calling used a pooled MACS3 call across multiple inputs; "
            "per-replicate overlap/IDR reproducibility was not run."
        )
        return repro

    rep_dir = Path(output_dir) / "replicate_peaks"
    rep_dir.mkdir(parents=True, exist_ok=True)
    peak_files: list[str] = []
    failed: list[dict] = []
    for i, input_file in enumerate(valid_files):
        label = _sample_label(input_file, i)
        cmd = _build_macs3_cmd(
            [input_file], [], genome_size, label, str(rep_dir), macs3_params,
        )
        run = _run_macs3(cmd, timeout=7200)
        if run.get("status") != "success":
            failed.append({
                "input_file": input_file,
                "error_type": run.get("error_type", "MACS3Failed"),
                "details": run.get("details"),
            })
            continue
        peak_file = _find_peak_file(
            str(rep_dir), label, peaks_ext, allow_fallback=False)
        if peak_file:
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
            "Some per-replicate MACS3 peak calls failed; reproducibility "
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
    repro["consensus_peaks_path"] = consensus if overlap.get("ran") else None
    if overlap.get("ran") and overlap.get("n_reproducible_regions", 0) > 0:
        repro["status"] = "verified" if not failed else "partial"
    else:
        repro["status"] = "unverified"
        repro["reason"] = "no_reproducible_overlap_regions"
    repro["warnings"].append(
        "IDR was not run; ARIA used an overlap-support reproducible peak "
        "policy across per-replicate MACS3 calls."
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


def _compute_frip(bam_file: str, peaks_file: str) -> float | None:
    """
    Compute FRiP (Fraction of Reads in Peaks).
    Requires samtools and bedtools in PATH.
    """
    try:
        import subprocess

        # Count reads in peaks
        cmd = [
            "bedtools", "intersect",
            "-a", bam_file,
            "-b", peaks_file,
            "-u", "-f", "0.5",
        ]
        reads_in_peaks = subprocess.run(
            ["samtools", "view", "-c", "-F", "4"],
            input=subprocess.run(cmd, capture_output=True).stdout,
            capture_output=True, text=True, timeout=120,
        )

        # Total mapped reads
        total_reads = subprocess.run(
            ["samtools", "view", "-c", "-F", "4", bam_file],
            capture_output=True, text=True, timeout=120,
        )

        in_peaks = int(reads_in_peaks.stdout.strip() or 0)
        total    = int(total_reads.stdout.strip() or 1)
        return in_peaks / max(total, 1)

    except Exception:
        return None


if __name__ == "__main__":
    run_script(chromatin_peaks)
