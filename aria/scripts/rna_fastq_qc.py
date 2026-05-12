"""
ARIA RNA-seq FASTQ QC Script
------------------------------
Runs fastp (trimming + QC) on all samples, then MultiQC for summary.
Executed inside aria-rnaseq-env via EnvironmentManager.

Detects paired-end vs single-end automatically from file naming.
Naming patterns supported:
  sample_1.fq.gz / sample_2.fq.gz
  sample_R1.fastq.gz / sample_R2.fastq.gz
  sample_R1_001.fastq.gz / sample_R2_001.fastq.gz

Input params:
  fastq_dir:      str   — directory containing FASTQ files
  output_dir:     str   — where to write trimmed FASTQs + QC reports
  threads:        int   — threads per sample (default: 4)
  min_length:     int   — minimum read length after trimming (default: 36)
  quality:        int   — phred quality threshold (default: 20)

Output:
  {
    "status":        "success",
    "samples":       [{"name", "r1_trimmed", "r2_trimmed", "r1_raw", "r2_raw",
                       "paired", "n_reads_raw", "n_reads_trimmed", "pct_passed"}],
    "multiqc_report": str,   — path to MultiQC HTML
    "output_dir":    str,
    "warnings":      [str]
  }
"""

from __future__ import annotations
import sys, os, subprocess, json, re
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script


def rna_fastq_qc(params: dict) -> dict:
    fastq_dir  = Path(params["fastq_dir"])
    output_dir = Path(params.get("output_dir",
                                  str(fastq_dir.parent / "aria_qc")))
    threads    = int(params.get("threads", 4))
    min_len    = int(params.get("min_length", 36))
    quality    = int(params.get("quality", 20))
    allow_mock = mocks_allowed(params)
    warnings   = []

    trimmed_dir = output_dir / "trimmed"
    qc_dir      = output_dir / "fastqc"
    fastp_dir   = output_dir / "fastp_reports"
    for d in [trimmed_dir, qc_dir, fastp_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ── Detect samples ────────────────────────────────────────────────────
    samples = _detect_samples(fastq_dir, warnings)
    if not samples:
        return {
            "status":     "error",
            "error_type": "NoFASTQFound",
            "details":    f"No FASTQ files found in {fastq_dir}",
        }

    # ── Run fastp per sample ──────────────────────────────────────────────
    processed = []
    for sample in samples:
        result = _run_fastp(
            sample=sample,
            trimmed_dir=trimmed_dir,
            fastp_dir=fastp_dir,
            threads=threads,
            min_len=min_len,
            quality=quality,
            warnings=warnings,
            allow_mock=allow_mock,
        )
        processed.append(result)

    # ── MultiQC summary ───────────────────────────────────────────────────
    multiqc_report = _run_multiqc(fastp_dir, output_dir, warnings)

    return {
        "status":         "success",
        "samples":        processed,
        "n_samples":      len(processed),
        "multiqc_report": str(multiqc_report) if multiqc_report else None,
        "output_dir":     str(output_dir),
        "trimmed_dir":    str(trimmed_dir),
        "warnings":       warnings,
    }


def _detect_samples(fastq_dir: Path, warnings: list) -> list[dict]:
    """
    Detect samples and pair R1/R2 files.
    Returns list of dicts: {name, r1, r2, paired}
    """
    # Find all FASTQ files
    fq_files = sorted([
        f for f in fastq_dir.iterdir()
        if f.suffix in (".gz", ".fastq", ".fq")
        or f.name.endswith((".fastq.gz", ".fq.gz"))
    ])

    if not fq_files:
        return []

    # Pair R1/R2 using common naming patterns
    r1_patterns = [r"_R1_\d+", r"_R1\b", r"_1\b"]
    r2_patterns = [r"_R2_\d+", r"_R2\b", r"_2\b"]

    r1_files = {}  # sample_name → path
    r2_files = {}

    for f in fq_files:
        stem = f.name
        # Remove all FASTQ extensions
        for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            stem = stem.replace(ext, "")

        # Try to identify R1
        for pat in r1_patterns:
            if re.search(pat, stem):
                sample = re.sub(pat, "", stem)
                r1_files[sample] = f
                break
        else:
            # Try R2
            for pat in r2_patterns:
                if re.search(pat, stem):
                    sample = re.sub(pat, "", stem)
                    r2_files[sample] = f
                    break

    samples = []
    all_sample_names = set(r1_files) | set(r2_files)

    for name in sorted(all_sample_names):
        r1 = r1_files.get(name)
        r2 = r2_files.get(name)
        paired = r1 is not None and r2 is not None

        if r1 is None and r2 is not None:
            # Only R2 found — treat as single-end with warning
            warnings.append(
                f"Sample '{name}': only R2 found, no R1. Skipping."
            )
            continue

        if r1 is None:
            continue

        samples.append({
            "name":   name,
            "r1":     str(r1),
            "r2":     str(r2) if r2 else None,
            "paired": paired,
        })

    if not samples:
        warnings.append(
            "Could not pair R1/R2 files. "
            "Expected naming: sample_1.fq.gz / sample_2.fq.gz "
            "or sample_R1.fastq.gz / sample_R2.fastq.gz"
        )

    return samples


def _fastp_outputs_valid(r1_out: str, r2_out: str, json_out: str,
                          paired: bool) -> bool:
    """
    Check if fastp outputs are present and valid for skipping.
    Requires:
      - R1 trimmed file > 1 KB (non-empty gzip)
      - R2 trimmed file > 1 KB if paired
      - JSON report > 100 bytes (non-empty report)
    Note: a valid gzip file has a minimum empty-payload size of ~20 bytes;
    we use 1 KB as a conservative real-data floor.
    """
    try:
        r1 = Path(r1_out)
        if not r1.exists() or r1.stat().st_size < 1024:
            return False
        if paired:
            if not r2_out:
                return False
            r2 = Path(r2_out)
            if not r2.exists() or r2.stat().st_size < 1024:
                return False
        j = Path(json_out)
        if not j.exists() or j.stat().st_size < 100:
            return False
        return True
    except Exception:
        return False


def _run_fastp(sample: dict, trimmed_dir: Path, fastp_dir: Path,
               threads: int, min_len: int, quality: int,
               warnings: list,
               allow_mock: bool = False) -> dict:
    """Run fastp on one sample. Idempotent: skips if outputs already valid."""
    name   = sample["name"]
    r1_in  = sample["r1"]
    r2_in  = sample.get("r2")
    paired = sample["paired"]

    r1_out = str(trimmed_dir / f"{name}_R1_trimmed.fq.gz")
    r2_out = str(trimmed_dir / f"{name}_R2_trimmed.fq.gz") if paired else None
    json_out = str(fastp_dir / f"{name}_fastp.json")
    html_out = str(fastp_dir / f"{name}_fastp.html")

    result = {
        "name":       name,
        "r1_raw":     r1_in,
        "r2_raw":     r2_in,
        "r1_trimmed": r1_out,
        "r2_trimmed": r2_out,
        "paired":     paired,
    }

    # ── Resume check: do all outputs exist and parse cleanly? ────────────
    if _fastp_outputs_valid(r1_out, r2_out, json_out, paired):
        try:
            with open(json_out) as f:
                fastp_data = json.load(f)
            summary = fastp_data.get("summary", {})
            before  = summary.get("before_filtering", {})
            after   = summary.get("after_filtering", {})
            n_raw     = int(before.get("total_reads", 0))
            n_trimmed = int(after.get("total_reads",  0))
            pct_pass  = round(n_trimmed / max(n_raw, 1) * 100, 1)
            result.update({
                "status":          "success",
                "n_reads_raw":     n_raw,
                "n_reads_trimmed": n_trimmed,
                "pct_passed":      pct_pass,
                "q30_rate":        round(float(after.get("q30_rate", 0))
                                          * 100, 1),
                "resumed":         True,
            })
            warnings.append(
                f"[resume] fastp skipped for {name} "
                f"(valid outputs exist: {pct_pass}% passed previously)"
            )
            return result
        except Exception:
            # JSON corrupted — fall through to re-run fastp
            pass

    cmd = [
        "fastp",
        "--in1",  r1_in,
        "--out1", r1_out,
        "--json", json_out,
        "--html", html_out,
        "--thread",      str(threads),
        "--length_required", str(min_len),
        "--qualified_quality_phred", str(quality),
        "--detect_adapter_for_pe" if paired else "--disable_adapter_trimming",
        "--correction" if paired else "--disable_quality_filtering",
        "--overrepresentation_analysis",
    ]
    if paired and r2_in and r2_out:
        cmd += ["--in2", r2_in, "--out2", r2_out]
    cmd = [c for c in cmd if c]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )
        if proc.returncode != 0:
            warnings.append(
                f"fastp failed for {name}: {proc.stderr[-200:]}"
            )
            result["status"] = "failed"
            return result

        # Parse fastp JSON for QC metrics
        if Path(json_out).exists():
            with open(json_out) as f:
                fastp_data = json.load(f)
            summary = fastp_data.get("summary", {})
            before  = summary.get("before_filtering", {})
            after   = summary.get("after_filtering", {})

            n_raw     = int(before.get("total_reads", 0))
            n_trimmed = int(after.get("total_reads",  0))
            pct_pass  = round(n_trimmed / max(n_raw, 1) * 100, 1)

            result.update({
                "status":          "success",
                "n_reads_raw":     n_raw,
                "n_reads_trimmed": n_trimmed,
                "pct_passed":      pct_pass,
                "q30_rate":        round(float(
                    after.get("q30_rate", 0)) * 100, 1),
            })

            if pct_pass < 70:
                warnings.append(
                    f"Sample {name}: only {pct_pass}% reads passed QC. "
                    f"Check library quality."
                )
        else:
            result["status"] = "success"

    except subprocess.TimeoutExpired:
        warnings.append(f"fastp timed out for {name}")
        result["status"] = "timeout"
    except FileNotFoundError:
        if allow_mock:
            warnings.append(
                "fastp not found — using mock QC result (explicit mock mode)"
            )
            result["status"] = "missing_tool"
            result.update(_mock_fastp_result(name, r1_in, r2_in,
                                              r1_out, r2_out, paired))
        else:
            warnings.append("fastp not found — FASTQ QC failed.")
            result["status"] = "failed"
            result["error_type"] = "MissingDependency"
            result["details"] = "fastp is required for raw FASTQ preprocessing."

    return result


def _run_multiqc(fastp_dir: Path, output_dir: Path,
                  warnings: list) -> Path | None:
    """Run MultiQC over fastp reports."""
    try:
        proc = subprocess.run(
            ["multiqc", str(fastp_dir),
             "--outdir", str(output_dir),
             "--filename", "multiqc_report",
             "--force", "--quiet"],
            capture_output=True, text=True, timeout=300
        )
        report = output_dir / "multiqc_report.html"
        if report.exists():
            return report
        warnings.append(f"MultiQC ran but report not found: {proc.stderr[-100:]}")
    except FileNotFoundError:
        warnings.append("multiqc not found — skipping summary report")
    except Exception as e:
        warnings.append(f"MultiQC failed: {e}")
    return None


def _mock_fastp_result(name, r1_in, r2_in, r1_out, r2_out, paired) -> dict:
    """Mock result when fastp not installed."""
    import random
    random.seed(hash(name) % 1000)
    n_raw = random.randint(20_000_000, 50_000_000)
    pct   = random.uniform(92, 98)
    return {
        "status":          "success",
        "n_reads_raw":     n_raw,
        "n_reads_trimmed": int(n_raw * pct / 100),
        "pct_passed":      round(pct, 1),
        "q30_rate":        round(random.uniform(88, 95), 1),
        "note":            "mock — fastp not installed",
    }


if __name__ == "__main__":
    run_script(rna_fastq_qc)
