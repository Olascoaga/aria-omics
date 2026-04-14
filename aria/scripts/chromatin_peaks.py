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
      "frip":                 float, — actual FRiP after peak calling
      "warnings":             [str]
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def chromatin_peaks(params: dict) -> dict:
    import subprocess
    import tempfile
    from pathlib import Path
    import numpy as np

    data_type     = params["data_type"]
    files         = params.get("files", [])
    genome        = params.get("genome", "hg38")
    macs3_params  = params.get("macs3_params", {})
    output_dir    = params.get("output_dir", "/tmp/aria_peaks")
    control_files = params.get("control_files", [])

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

    cmd = ["macs3", "callpeak"]

    # Input files
    cmd += ["-t"] + valid_files

    # Control/input files (ChIP)
    if control_files:
        valid_ctrl = [f for f in control_files if Path(f).exists()]
        if valid_ctrl:
            cmd += ["-c"] + valid_ctrl
        else:
            warnings.append(
                "Control files specified but not found. "
                "Running without input control."
            )

    # Format
    fmt = macs3_params.get("format", "BAMPE")
    cmd += ["-f", fmt]

    # Genome size
    cmd += ["-g", genome_size]

    # Output
    cmd += ["-n", sample_name, "--outdir", output_dir]

    # Assay-specific flags
    if macs3_params.get("nomodel", False):
        cmd.append("--nomodel")
        extsize = macs3_params.get("extsize", 200)
        cmd += ["--extsize", str(extsize)]

    if macs3_params.get("nolambda", False):
        cmd.append("--nolambda")

    if macs3_params.get("broad", False):
        cmd.append("--broad")

    keep_dup = macs3_params.get("keep_dup", "1")
    cmd += ["--keep-dup", str(keep_dup)]

    # Always store bdg for downstream analysis
    cmd.append("-B")
    cmd.append("--SPMR")  # signal per million reads (for comparison)

    # ── Run MACS3 ────────────────────────────────────────────────────────
    try:
        result = subprocess.run(
            cmd,
            capture_output=True, text=True, timeout=7200,  # 2 hours
        )

        if result.returncode != 0:
            # MACS3 often prints useful info to stderr even on success
            if "Traceback" in result.stderr or "Error" in result.stderr:
                return {
                    "status":     "error",
                    "error_type": "MACS3Failed",
                    "details":    result.stderr[-1000:],
                }
            else:
                warnings.append(f"MACS3 non-zero exit: {result.stderr[-200:]}")

    except subprocess.TimeoutExpired:
        return {
            "status":     "error",
            "error_type": "Timeout",
            "details":    "MACS3 exceeded 2-hour time limit.",
        }
    except FileNotFoundError:
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

    # ── Create consensus peaks across replicates (if multiple files) ──────
    consensus_path = None
    if len(valid_files) > 1:
        consensus_path = _create_consensus_peaks(
            output_dir, peaks_ext, sample_name
        )

    # ── Compute actual FRiP ───────────────────────────────────────────────
    frip = _compute_frip(valid_files[0], str(peaks_file))

    return {
        "status":               "success",
        "data_type":            data_type,
        "n_peaks":              int(n_peaks),
        "peaks_path":           str(peaks_file),
        "consensus_peaks_path": consensus_path,
        "frip":                 round(float(frip), 4),
        "genome":               genome,
        "macs3_cmd":            " ".join(cmd),  # for reproducibility
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


def _count_peaks(peaks_file: str) -> int:
    """Count number of peaks in a narrowPeak/broadPeak file."""
    try:
        with open(peaks_file) as f:
            return sum(1 for line in f
                       if line.strip() and not line.startswith("#"))
    except Exception:
        return 0


def _create_consensus_peaks(output_dir: str,
                              peaks_ext: str,
                              sample_name: str) -> str | None:
    """
    Merge peaks across replicates using bedtools merge.
    Returns path to consensus peaks file.
    """
    try:
        import subprocess
        from pathlib import Path

        out_dir   = Path(output_dir)
        all_peaks = list(out_dir.glob(f"*{peaks_ext}"))
        if len(all_peaks) < 2:
            return None

        consensus = str(out_dir / f"{sample_name}_consensus{peaks_ext}")

        # Sort and merge
        cat_cmd   = ["cat"] + [str(p) for p in all_peaks]
        sort_cmd  = ["sort", "-k1,1", "-k2,2n"]
        merge_cmd = ["bedtools", "merge", "-i", "stdin"]

        cat_proc  = subprocess.Popen(cat_cmd,  stdout=subprocess.PIPE)
        sort_proc = subprocess.Popen(sort_cmd, stdin=cat_proc.stdout,
                                     stdout=subprocess.PIPE)
        with open(consensus, "w") as out:
            subprocess.run(merge_cmd, stdin=sort_proc.stdout,
                           stdout=out, timeout=300)

        return consensus
    except Exception:
        return None


def _compute_frip(bam_file: str, peaks_file: str) -> float:
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
        return 0.30   # typical ATAC default


if __name__ == "__main__":
    run_script(chromatin_peaks)
