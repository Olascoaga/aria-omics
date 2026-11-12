"""
ARIA RNA-seq STAR Alignment Script
------------------------------------
Aligns trimmed FASTQs to reference genome using STAR.
Executed inside aria-rnaseq-env via EnvironmentManager.

Handles:
  - Genome index building if not present
  - Paired-end and single-end alignment
  - Coordinate-sorted BAM output
  - Alignment QC metrics (% mapped, multimappers)

Input params:
  samples:        list  — sample/lane/read-layout manifest from rna_fastq_qc
  genome_dir:     str   — STAR genome index directory
  genome_fasta:   str   — (optional) FASTA for building index if missing
  gtf_file:       str   — GTF annotation file
  output_dir:     str
  threads:        int   (default: 8)
  two_pass:       bool  (default: True) — STAR 2-pass for novel junctions

Output:
  {
    "status":   "success",
    "bam_files": [{"name", "bam", "read_layout", "pct_unique", "n_reads"}],
    "output_dir": str,
    "warnings":  [str]
  }
"""

from __future__ import annotations
import math, sys, os, subprocess, re
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script
from aria.utils.stage_manifest import (
    build_stage_manifest, write_stage_manifest, stage_is_current,
)

_INDEX_MANIFEST = "aria_index.stage.json"


def _star_version() -> str:
    """Best-effort STAR version string for the A4 stage manifests."""
    try:
        proc = subprocess.run(["STAR", "--version"],
                              capture_output=True, text=True, timeout=30)
        return (proc.stdout or proc.stderr or "").strip()
    except Exception:
        return ""


def _sa_index_nbases(fasta_path: str) -> int:
    """Compute STAR --genomeSAindexNbases per STAR manual:
    min(14, floor(log2(GenomeLength)/2 - 1)), bounded below by 4."""
    if not fasta_path:
        return 14
    try:
        p = Path(fasta_path)
        if not p.is_file():
            return 14
        genome_len = 0
        with open(p, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith(">"):
                    genome_len += len(line.strip())
        if genome_len <= 0:
            return 14
        return min(14, max(4, int(math.floor(math.log2(genome_len) / 2.0 - 1.0))))
    except Exception:
        return 14


def rna_align(params: dict) -> dict:
    samples    = params["samples"]           # from rna_fastq_qc
    genome_dir = Path(params["genome_dir"])
    gtf_file   = params.get("gtf_file", "")
    output_dir = Path(params.get("output_dir",
                                  str(Path(samples[0]["r1_trimmed"]).parent.parent
                                      / "aria_aligned")))
    threads    = int(params.get("threads", 8))
    library_protocol = str(params.get("library_protocol") or "full_length").strip().lower()
    default_two_pass = library_protocol == "full_length"
    two_pass   = bool(params["two_pass"]) if "two_pass" in params else default_two_pass
    genome_fasta = params.get("genome_fasta", "")
    allow_mock = mocks_allowed(params)
    warnings   = []

    output_dir.mkdir(parents=True, exist_ok=True)

    if not samples:
        return {
            "status": "error",
            "error_type": "NoSamples",
            "details": "No QC sample manifest was provided for alignment.",
            "warnings": warnings,
        }

    upstream_failed = [
        sample.get("name", "unknown")
        for sample in samples
        if sample.get("status") not in (None, "success")
    ]
    if upstream_failed:
        return {
            "status": "error",
            "error_type": "UpstreamFastqQCFailure",
            "details": "STAR alignment requires every QC sample to succeed.",
            "failed_samples": upstream_failed,
            "warnings": warnings,
        }

    # ── Build genome index if needed ──────────────────────────────────────
    sa_nbases = _sa_index_nbases(genome_fasta)
    index_manifest = str(genome_dir / _INDEX_MANIFEST)
    index_inputs = [("fasta", genome_fasta or None), ("gtf", gtf_file or None)]
    index_params = {"sa_index_nbases": sa_nbases}
    star_version = _star_version()
    need_build = not _index_exists(genome_dir)
    # A4: when we have the source (FASTA + GTF), a changed reference/GTF/STAR
    # version must rebuild the index instead of silently reusing a stale one. An
    # external index provided WITHOUT its source keeps the existence check (no
    # source to rebuild or hash the reference from).
    if not need_build and genome_fasta and gtf_file:
        index_current, index_reason = stage_is_current(
            index_manifest, inputs=index_inputs, params=index_params,
            tool_version=star_version,
        )
        if not index_current:
            warnings.append(f"[A4] rebuilding STAR index (stale: {index_reason})")
            need_build = True
    if need_build:
        if genome_fasta and gtf_file:
            build_result = _build_star_index(
                genome_dir, genome_fasta, gtf_file, threads, warnings,
                sa_index_nbases=sa_nbases,
            )
            if build_result != "ok":
                return {
                    "status":     "error",
                    "error_type": "IndexBuildFailed",
                    "details":    build_result,
                }
            try:
                write_stage_manifest(index_manifest, build_stage_manifest(
                    stage="star_index", inputs=index_inputs,
                    params=index_params, tool_version=star_version,
                ))
            except Exception as exc:
                warnings.append(f"[A4] could not write STAR index manifest: {exc}")
        else:
            return {
                "status":     "error",
                "error_type": "MissingGenomeIndex",
                "details": (
                    f"STAR index not found at {genome_dir}. "
                    f"Provide genome_fasta + gtf_file to build it, "
                    f"or set genome_dir to an existing index. "
                    f"For human/mouse, run: aria download-genome hg38|mm39"
                ),
            }

    # ── Align each sample ─────────────────────────────────────────────────
    bam_files = []
    for sample in samples:
        result = _align_sample(
            sample=sample,
            genome_dir=genome_dir,
            output_dir=output_dir,
            threads=threads,
            two_pass=two_pass,
            warnings=warnings,
            allow_mock=allow_mock,
            index_ref=(index_manifest if os.path.exists(index_manifest)
                       else str(genome_dir / "SAindex")),
            star_version=star_version,
            library_protocol=library_protocol,
        )
        bam_files.append(result)

    n_ok = sum(1 for b in bam_files if b.get("status") == "success")
    failed_samples = [
        result.get("name", "unknown")
        for result in bam_files
        if result.get("status") != "success"
    ]
    if failed_samples:
        return {
            "status":     "error",
            "error_type": "PartialAlignmentFailure",
            "details":    (
                "STAR alignment is all-or-fail; one or more samples did not "
                "complete successfully."
            ),
            "failed_samples": failed_samples,
            "bam_files": bam_files,
            "warnings":   warnings,
        }

    return {
        "status":    "success",
        "bam_files": bam_files,
        "n_aligned": n_ok,
        "output_dir": str(output_dir),
        "warnings":  warnings,
    }


def _index_exists(genome_dir: Path) -> bool:
    """Check if STAR index exists by looking for key files."""
    required = ["SA", "SAindex", "Genome"]
    return genome_dir.exists() and all(
        (genome_dir / f).exists() for f in required
    )


def _build_star_index(genome_dir: Path, fasta: str, gtf: str,
                       threads: int, warnings: list,
                       *,
                       sa_index_nbases: int = 14) -> str:
    """Build STAR genome index."""
    genome_dir.mkdir(parents=True, exist_ok=True)
    warnings.append(f"Building STAR index at {genome_dir} (SAindexNbases={sa_index_nbases})")

    cmd = [
        "STAR", "--runMode", "genomeGenerate",
        "--genomeDir",       str(genome_dir),
        "--genomeFastaFiles", fasta,
        "--sjdbGTFfile",     gtf,
        "--runThreadN",      str(threads),
        "--genomeSAindexNbases", str(sa_index_nbases),
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=7200
        )
        if proc.returncode != 0:
            return (proc.stderr or proc.stdout or "STAR genomeGenerate failed")[-300:]
        return "ok"
    except FileNotFoundError:
        return "STAR not found — is aria-rnaseq-env active?"
    except subprocess.TimeoutExpired:
        return "STAR index build timed out (>2h)"


def _star_output_valid(bam_path: Path, log_file: Path) -> bool:
    """
    Check if STAR output is present and complete.
    Requires:
      - Non-empty BAM file exists (> 128 bytes)
      - Log.final.out present (STAR writes this only after finishing)
      - samtools quickcheck passes (BAM EOF marker intact, not truncated)
    """
    try:
        if not bam_path.exists() or bam_path.stat().st_size < 128:
            return False
        if not log_file.exists():
            return False
        # samtools quickcheck exit 0 = complete, non-zero = truncated/corrupt
        proc = subprocess.run(
            ["samtools", "quickcheck", str(bam_path)],
            capture_output=True, timeout=60
        )
        return proc.returncode == 0
    except Exception:
        return False


def _align_sample(sample: dict, genome_dir: Path, output_dir: Path,
                  threads: int, two_pass: bool,
                  warnings: list,
                  allow_mock: bool = False,
                  index_ref: str | None = None,
                  star_version: str = "",
                  *,
                  library_protocol: str = "full_length") -> dict:
    """Align one sample with STAR. Idempotent: skips if BAM already valid."""
    name   = sample["name"]
    r1_files = sample.get("r1_trimmed_files") or [
        sample.get("r1_trimmed") or sample.get("r1")
    ]
    r2_files = sample.get("r2_trimmed_files") or [
        path for path in [sample.get("r2_trimmed") or sample.get("r2")]
        if path
    ]
    read_layout = sample.get(
        "read_layout",
        "paired-end" if sample.get("paired", bool(r2_files)) else "single-end",
    )
    paired = read_layout == "paired-end"
    r1 = ",".join(r1_files)
    r2 = ",".join(r2_files) if r2_files else None

    sample_dir = output_dir / name
    sample_dir.mkdir(exist_ok=True)
    prefix = str(sample_dir / f"{name}_")

    result = {
        "name":        name,
        "status":      "pending",
        "read_layout": read_layout,
        "paired":      paired,
        "library_protocol": library_protocol,
    }

    # A4: content-addressed resume manifest for this alignment stage. Inputs are
    # the trimmed reads plus the genome-index reference (its own manifest when
    # present, so a rebuilt index invalidates alignment through the DAG).
    manifest_path = str(sample_dir / f"{name}.align.stage.json")
    stage_inputs = [
        *[(f"r1_lane_{index}", path)
          for index, path in enumerate(r1_files, start=1)],
        *[(f"r2_lane_{index}", path)
          for index, path in enumerate(r2_files, start=1)],
        ("index", index_ref),
    ]
    stage_params: dict = {"two_pass": two_pass, "read_layout": read_layout}
    if library_protocol != "full_length":
        stage_params["library_protocol"] = library_protocol

    # ── Resume check: valid BAM AND the stage manifest still matches? ─────
    bam_path = Path(f"{prefix}Aligned.sortedByCoord.out.bam")
    log_file = Path(f"{prefix}Log.final.out")
    manifest_current, manifest_reason = stage_is_current(
        manifest_path, inputs=stage_inputs, params=stage_params,
        tool_version=star_version,
    )
    if _star_output_valid(bam_path, log_file) and manifest_current:
        stats = _parse_star_log(log_file)
        warnings.append(
            f"[resume] STAR skipped for {name} "
            f"(valid BAM exists: {stats.get('pct_unique',0)}% unique)"
        )
        result.update({
            "status":      "success",
            "bam":         str(bam_path),
            "log":         str(log_file),
            "pct_unique":  stats.get("pct_unique", 0),
            "pct_multi":   stats.get("pct_multi",  0),
            "pct_unmapped":stats.get("pct_unmapped", 0),
            "n_input":     stats.get("n_input", 0),
            "n_reads":     stats.get("n_reads", 0),
            "resumed":     True,
        })
        return result

    is_gzipped = str(r1_files[0] or "").lower().endswith(".gz")
    cmd = [
        "STAR",
        "--runThreadN",          str(threads),
        "--genomeDir",           str(genome_dir),
        "--readFilesIn",         r1,
    ]
    if paired and r2:
        cmd.append(r2)
    if is_gzipped:
        cmd += ["--readFilesCommand", "zcat"]
    cmd += [
        "--outSAMtype",          "BAM", "SortedByCoordinate",
        "--outSAMattributes",    "NH", "HI", "AS", "NM", "MD",
        "--outFileNamePrefix",   prefix,
        "--outSAMstrandField",   "intronMotif",
        "--outSAMunmapped",      "Within",
        "--quantMode",           "GeneCounts",
    ]
    # ADR-075: 3'-tag, 5'-tag, and tag-UMI protocols align primarily within
    # terminal exons/UTRs/TSSs with local soft-clipping; omit RemoveNoncanonical
    # so short-overhang or unspliced terminal tags are retained.
    if library_protocol == "full_length":
        cmd += ["--outFilterIntronMotifs", "RemoveNoncanonical"]

    if two_pass:
        cmd += ["--twopassMode", "Basic"]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=7200
        )

        if proc.returncode != 0:
            warnings.append(
                f"STAR failed for {name}: {proc.stderr[-200:]}"
            )
            result["status"] = "failed"
            return result

        # Find output BAM
        bam = Path(f"{prefix}Aligned.sortedByCoord.out.bam")
        if not bam.exists():
            warnings.append(f"BAM not found for {name}: {prefix}")
            result["status"] = "failed"
            return result

        # Index BAM
        subprocess.run(
            ["samtools", "index", str(bam)],
            capture_output=True, timeout=300
        )

        # Parse alignment stats from Log.final.out
        log_file = Path(f"{prefix}Log.final.out")
        stats    = _parse_star_log(log_file)

        if stats.get("pct_unique", 0) < 50:
            warnings.append(
                f"Sample {name}: low unique mapping rate "
                f"({stats.get('pct_unique', 0):.1f}%). "
                f"Check library quality or genome compatibility."
            )

        result.update({
            "status":      "success",
            "bam":         str(bam),
            "counts_tab":  str(Path(f"{prefix}ReadsPerGene.out.tab")),
            **stats,
        })
        # A4: record the content-addressed manifest for this genuine run so a
        # later resume can tell whether reads/index/params/version still match.
        try:
            write_stage_manifest(manifest_path, build_stage_manifest(
                stage="star_align", inputs=stage_inputs, params=stage_params,
                tool_version=star_version,
            ))
        except Exception as exc:
            warnings.append(f"[A4] could not write STAR align manifest for "
                            f"{name}: {exc}")

    except FileNotFoundError:
        if allow_mock:
            warnings.append(
                f"STAR not found for {name} — using mock alignment "
                "(explicit mock mode)"
            )
            result.update(_mock_alignment(name, prefix))
        else:
            warnings.append(f"STAR not found for {name}; alignment failed.")
            result["status"] = "failed"
            result["error_type"] = "MissingDependency"
            result["details"] = "STAR is required for RNA-seq alignment."

    except subprocess.TimeoutExpired:
        warnings.append(f"STAR timed out for {name} (>2h)")
        result["status"] = "timeout"

    return result


def _parse_star_log(log_file: Path) -> dict:
    """Parse STAR Log.final.out for alignment metrics."""
    if not log_file.exists():
        return {}

    stats = {}
    patterns = {
        "n_reads":       r"Number of input reads \|\s+([\d,]+)",
        "pct_unique":    r"Uniquely mapped reads % \|\s+([\d.]+)%",
        "pct_multi":     r"% of reads mapped to multiple loci \|\s+([\d.]+)%",
        "pct_unmapped":  r"% of reads unmapped: too many mismatches \|\s+([\d.]+)%",
    }

    text = log_file.read_text()
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            val = m.group(1).replace(",", "")
            stats[key] = float(val) if "." in val else int(val)

    return stats


def _mock_alignment(name: str, prefix: str) -> dict:
    """Mock alignment result when STAR not installed."""
    import random
    import zlib
    # A2: zlib.crc32 is stable across processes; builtin hash() is randomized by
    # PYTHONHASHSEED, so the same sample produced different mock numbers per run.
    random.seed(zlib.crc32(str(name).encode()))
    return {
        "status":      "success",
        "bam":         f"{prefix}Aligned.sortedByCoord.out.bam",
        "counts_tab":  f"{prefix}ReadsPerGene.out.tab",
        "n_reads":     random.randint(20_000_000, 50_000_000),
        "pct_unique":  round(random.uniform(80, 92), 1),
        "pct_multi":   round(random.uniform(5, 12), 1),
        "note":        "mock — STAR not installed",
    }


if __name__ == "__main__":
    run_script(rna_align)
