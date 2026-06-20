"""
ARIA bulk ATAC-seq FASTQ -> BAM Alignment Script
-------------------------------------------------
Aligns trimmed ATAC-seq FASTQs to a reference genome with bwa-mem2 and produces
the coordinate-sorted, duplicate-marked, ATAC-filtered BAM that the bulk ATAC
lane (QC -> MACS3 peaks -> DA -> ...) consumes.

Executed inside aria-atacseq-env via EnvironmentManager.

This is the ATAC analogue of ``rna_align.py`` (STAR): it consumes the output of
``rna_fastq_qc.py`` (assay-agnostic fastp trimming) and emits one filtered BAM
per sample. It does NOT call peaks or quantify; downstream chromatin scripts do.

ATAC-specific post-alignment filtering (ENCODE-style):
  - mark + remove PCR/optical duplicates (samtools fixmate -m + markdup)
  - keep properly paired alignments (-f 2), drop unmapped/mate-unmapped/
    secondary/duplicate/QC-fail reads (-F 1804, the canonical ENCODE value;
    note it does NOT exclude supplementary 0x800)
  - MAPQ threshold (default 30)
  - drop mitochondrial reads (chrM / MT)

Input params:
  samples:        list  — from rna_fastq_qc: [{name, r1_trimmed|r1,
                          r2_trimmed|r2, paired}]
  genome_fasta:   str   — reference FASTA (used to build the bwa-mem2 index)
  bwa_index:      str   — (optional) existing bwa-mem2 index prefix; defaults to
                          genome_fasta (bwa-mem2 indexes alongside the FASTA)
  output_dir:     str
  threads:        int   (default: 8)
  mapq:           int   (default: 30)
  remove_mito:    bool  (default: True) — drop chrM/MT reads
  mito_contigs:   list  (default: ["chrM", "MT", "chrMT", "M"])

Output:
  {
    "status":   "success",
    "bam_files": [{"name", "bam", "n_reads", "pct_mapped", "pct_dup",
                   "pct_mito", "n_final"}],
    "output_dir": str,
    "aligner":   "bwa-mem2",
    "filter":    {"mapq": int, "samtools_F": 1804, "samtools_f": 2,
                  "remove_mito": bool},
    "warnings":  [str]
  }
"""

from __future__ import annotations
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script

# Canonical ENCODE ATAC read filter: drop unmapped (0x4), mate-unmapped (0x8),
# secondary (0x100), QC-fail (0x200), duplicate (0x400). NOTE: 1804 does NOT
# include supplementary (0x800) — this matches the ENCODE ATAC-seq spec, which
# retains supplementary alignments at this step (they are a tiny fraction and
# are gated by the -f 2 proper-pair requirement).
SAMTOOLS_EXCLUDE_FLAGS = 1804
SAMTOOLS_REQUIRE_FLAGS = 2  # properly paired
DEFAULT_MAPQ = 30
DEFAULT_MITO_CONTIGS = ("chrM", "MT", "chrMT", "M", "chrMt")


def atac_align(params: dict) -> dict:
    samples = params.get("samples") or []
    genome_fasta = params.get("genome_fasta", "")
    bwa_index = params.get("bwa_index") or genome_fasta
    output_dir = Path(params.get("output_dir",
                                 str(Path("aria_atac_aligned"))))
    threads = int(params.get("threads", 8))
    mapq = int(params.get("mapq", DEFAULT_MAPQ))
    remove_mito = bool(params.get("remove_mito", True))
    mito_contigs = tuple(params.get("mito_contigs") or DEFAULT_MITO_CONTIGS)
    allow_mock = mocks_allowed(params)
    warnings: list[str] = []

    if not samples:
        return {
            "status": "error",
            "error_type": "NoSamples",
            "details": "atac_align requires a non-empty 'samples' list "
                       "(from rna_fastq_qc).",
        }
    if not genome_fasta and not params.get("bwa_index"):
        return {
            "status": "error",
            "error_type": "MissingGenomeFasta",
            "details": "atac_align requires 'genome_fasta' (to build/locate the "
                       "bwa-mem2 index) or an explicit 'bwa_index' prefix.",
        }

    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Build bwa-mem2 index if missing ───────────────────────────────────
    if not _bwa_index_exists(bwa_index):
        if not genome_fasta or not Path(genome_fasta).exists():
            return {
                "status": "error",
                "error_type": "MissingGenomeIndex",
                "details": (
                    f"bwa-mem2 index not found at '{bwa_index}' and no readable "
                    f"genome_fasta to build it. Provide genome_fasta or a built "
                    f"index prefix."
                ),
            }
        build = _build_bwa_index(genome_fasta, threads, warnings)
        if build != "ok":
            if build == "missing_tool" and allow_mock:
                warnings.append("bwa-mem2 not found — emitting mock BAM records "
                                "(explicit mock mode).")
                return _mock_result(samples, output_dir, mapq, remove_mito,
                                    warnings)
            return {
                "status": "error",
                "error_type": ("MissingDependency"
                               if build == "missing_tool"
                               else "IndexBuildFailed"),
                "details": (
                    "bwa-mem2 is required for ATAC alignment "
                    "(install aria-atacseq-env)." if build == "missing_tool"
                    else build),
                "warnings": warnings,
            }
        bwa_index = genome_fasta

    # ── Align + filter each sample ────────────────────────────────────────
    bam_files = []
    for sample in samples:
        if sample.get("status") == "failed":
            warnings.append(f"Skipping {sample.get('name')} — QC failed upstream")
            continue
        result = _align_sample(
            sample=sample,
            bwa_index=bwa_index,
            output_dir=output_dir,
            threads=threads,
            mapq=mapq,
            remove_mito=remove_mito,
            mito_contigs=mito_contigs,
            warnings=warnings,
            allow_mock=allow_mock,
        )
        bam_files.append(result)

    n_ok = sum(1 for b in bam_files if b.get("status") == "success")
    if n_ok == 0:
        return {
            "status": "error",
            "error_type": "AllAlignmentsFailed",
            "details": "No ATAC samples aligned successfully.",
            "bam_files": bam_files,
            "warnings": warnings,
        }

    return {
        "status": "success",
        "bam_files": bam_files,
        "n_aligned": n_ok,
        "output_dir": str(output_dir),
        "aligner": "bwa-mem2",
        "filter": {
            "mapq": mapq,
            "samtools_F": SAMTOOLS_EXCLUDE_FLAGS,
            "samtools_f": SAMTOOLS_REQUIRE_FLAGS,
            "remove_mito": remove_mito,
            "mito_contigs": list(mito_contigs) if remove_mito else [],
        },
        "warnings": warnings,
    }


def _bwa_index_exists(prefix: str) -> bool:
    """bwa-mem2 writes <prefix>.bwt.2bit.64 + .pac/.ann/.amb/.0123 next to it."""
    if not prefix:
        return False
    p = Path(prefix)
    return (p.with_suffix(p.suffix + ".bwt.2bit.64").exists()
            or Path(str(prefix) + ".bwt.2bit.64").exists())


def _build_bwa_index(fasta: str, threads: int, warnings: list) -> str:
    warnings.append(f"Building bwa-mem2 index for {fasta} — this can take a while")
    try:
        proc = subprocess.run(
            ["bwa-mem2", "index", fasta],
            capture_output=True, text=True, timeout=14400,
        )
        if proc.returncode != 0:
            return proc.stderr[-300:]
        return "ok"
    except FileNotFoundError:
        return "missing_tool"
    except subprocess.TimeoutExpired:
        return "bwa-mem2 index build timed out (>4h)"


def _align_sample(sample: dict, bwa_index: str, output_dir: Path,
                  threads: int, mapq: int, remove_mito: bool,
                  mito_contigs: tuple, warnings: list,
                  allow_mock: bool = False) -> dict:
    """Align one ATAC sample with bwa-mem2 and apply ATAC filtering.

    Idempotent: skips if the final filtered BAM is already valid.
    """
    name = sample["name"]
    r1 = sample.get("r1_trimmed") or sample.get("r1")
    r2 = sample.get("r2_trimmed") or sample.get("r2")
    paired = sample.get("paired", r2 is not None)

    sample_dir = output_dir / name
    sample_dir.mkdir(exist_ok=True)
    final_bam = sample_dir / f"{name}.filtered.bam"

    result = {"name": name, "status": "pending", "paired": paired}

    if not r1:
        warnings.append(f"{name}: no R1 FASTQ — skipping")
        result["status"] = "failed"
        return result

    # ── Resume check ──────────────────────────────────────────────────────
    if _bam_valid(final_bam):
        warnings.append(f"[resume] atac_align skipped for {name} "
                        f"(valid filtered BAM exists)")
        result.update({"status": "success", "bam": str(final_bam),
                       "resumed": True})
        result.update(_bam_metrics(final_bam, mito_contigs, warnings))
        return result

    namesort = sample_dir / f"{name}.namesort.bam"
    fixmate = sample_dir / f"{name}.fixmate.bam"
    coordsort = sample_dir / f"{name}.coordsort.bam"
    markdup = sample_dir / f"{name}.markdup.bam"

    try:
        # 1. bwa-mem2 mem -> name-sorted BAM (name sort required by fixmate)
        align_cmd = ["bwa-mem2", "mem", "-t", str(threads), bwa_index, r1]
        if paired and r2:
            align_cmd.append(r2)
        sort_n = ["samtools", "sort", "-n", "-@", str(threads),
                  "-o", str(namesort), "-"]
        if not _pipe(align_cmd, sort_n, warnings, name, "align"):
            result["status"] = "failed"
            return result

        # 2. fixmate -m (adds mate score tags for markdup)
        if not _run(["samtools", "fixmate", "-m", "-@", str(threads),
                     str(namesort), str(fixmate)], warnings, name, "fixmate"):
            result["status"] = "failed"
            return result

        # 3. coordinate sort
        if not _run(["samtools", "sort", "-@", str(threads),
                     "-o", str(coordsort), str(fixmate)],
                    warnings, name, "coordsort"):
            result["status"] = "failed"
            return result

        # 4. mark + remove duplicates
        if not _run(["samtools", "markdup", "-r", "-@", str(threads),
                     str(coordsort), str(markdup)],
                    warnings, name, "markdup"):
            result["status"] = "failed"
            return result
        _run(["samtools", "index", str(markdup)], warnings, name, "index_markdup")

        # 5. ATAC filter: MAPQ + proper-pair + flag exclude (+ drop mito)
        view_cmd = ["samtools", "view", "-b", "-@", str(threads),
                    "-q", str(mapq),
                    "-f", str(SAMTOOLS_REQUIRE_FLAGS),
                    "-F", str(SAMTOOLS_EXCLUDE_FLAGS),
                    str(markdup)]
        if remove_mito:
            keep = _non_mito_contigs(markdup, mito_contigs, warnings)
            if keep is not None:
                view_cmd += keep
        view_cmd += ["-o", str(final_bam)]
        if not _run(view_cmd, warnings, name, "filter"):
            result["status"] = "failed"
            return result
        _run(["samtools", "index", str(final_bam)], warnings, name, "index_final")

        if not _bam_valid(final_bam):
            warnings.append(f"{name}: filtered BAM missing/invalid after run")
            result["status"] = "failed"
            return result

        result.update({"status": "success", "bam": str(final_bam)})
        result.update(_bam_metrics(markdup, mito_contigs, warnings,
                                   final_bam=final_bam))
        # Tidy intermediates (keep final + markdup for transparency)
        for tmp in (namesort, fixmate, coordsort):
            try:
                tmp.unlink()
            except OSError:
                pass

    except FileNotFoundError:
        if allow_mock:
            warnings.append("samtools/bwa-mem2 not found — mock BAM "
                            "(explicit mock mode)")
            result.update({"status": "success", "bam": str(final_bam),
                           "note": "mock — aligner not installed"})
        else:
            warnings.append("bwa-mem2/samtools not found — ATAC alignment failed")
            result["status"] = "failed"
            result["error_type"] = "MissingDependency"
            result["details"] = ("bwa-mem2 + samtools are required "
                                 "(install aria-atacseq-env).")
    return result


def _non_mito_contigs(bam: Path, mito_contigs: tuple,
                      warnings: list) -> list[str] | None:
    """Return the list of non-mito reference contigs (region args) for samtools.

    Uses idxstats so the filter only keeps reads on real (non-mito) contigs.
    Returns None when contigs cannot be listed (caller skips mito removal).
    """
    try:
        proc = subprocess.run(["samtools", "idxstats", str(bam)],
                              capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            warnings.append(f"idxstats failed; mito reads not removed: "
                            f"{proc.stderr[-120:]}")
            return None
        contigs = []
        mito_lower = {m.lower() for m in mito_contigs}
        for line in proc.stdout.splitlines():
            fields = line.split("\t")
            if not fields or fields[0] in ("*", ""):
                continue
            if fields[0].lower() in mito_lower:
                continue
            contigs.append(fields[0])
        return contigs or None
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        warnings.append(f"idxstats unavailable; mito reads not removed: {e}")
        return None


def _bam_metrics(pre_filter_bam: Path, mito_contigs: tuple, warnings: list,
                 final_bam: Path | None = None) -> dict:
    """Collect alignment QC: total/mapped/dup from flagstat, mito% from idxstats."""
    metrics: dict = {}
    try:
        fs = subprocess.run(["samtools", "flagstat", str(pre_filter_bam)],
                            capture_output=True, text=True, timeout=300)
        if fs.returncode == 0:
            total = mapped = dup = 0
            for line in fs.stdout.splitlines():
                if " in total" in line:
                    total = int(line.split()[0])
                elif " mapped (" in line and "primary" not in line:
                    mapped = int(line.split()[0])
                elif " duplicates" in line:
                    dup = int(line.split()[0])
            metrics["n_reads"] = total
            metrics["pct_mapped"] = round(mapped / max(total, 1) * 100, 2)
            metrics["pct_dup"] = round(dup / max(total, 1) * 100, 2)
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError) as e:
        warnings.append(f"flagstat metrics unavailable: {e}")

    try:
        idx = subprocess.run(["samtools", "idxstats", str(pre_filter_bam)],
                             capture_output=True, text=True, timeout=300)
        if idx.returncode == 0:
            mito_lower = {m.lower() for m in mito_contigs}
            mito = nuc = 0
            for line in idx.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 3 or f[0] in ("*", ""):
                    continue
                n = int(f[2])
                if f[0].lower() in mito_lower:
                    mito += n
                else:
                    nuc += n
            metrics["pct_mito"] = round(mito / max(mito + nuc, 1) * 100, 2)
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError) as e:
        warnings.append(f"idxstats metrics unavailable: {e}")

    if final_bam is not None:
        try:
            fs = subprocess.run(["samtools", "view", "-c", str(final_bam)],
                                capture_output=True, text=True, timeout=300)
            if fs.returncode == 0:
                metrics["n_final"] = int(fs.stdout.strip())
        except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
            pass
    return metrics


def _bam_valid(bam: Path) -> bool:
    try:
        if not bam.exists() or bam.stat().st_size < 1024:
            return False
        proc = subprocess.run(["samtools", "quickcheck", str(bam)],
                              capture_output=True, timeout=60)
        return proc.returncode == 0
    except Exception:
        return False


def _run(cmd: list, warnings: list, name: str, step: str) -> bool:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
        if proc.returncode != 0:
            warnings.append(f"{name}: {step} failed: {proc.stderr[-200:]}")
            return False
        return True
    except subprocess.TimeoutExpired:
        warnings.append(f"{name}: {step} timed out")
        return False


def _pipe(cmd1: list, cmd2: list, warnings: list, name: str, step: str) -> bool:
    """Run `cmd1 | cmd2`, surfacing either side's failure."""
    try:
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        p2 = subprocess.run(cmd2, stdin=p1.stdout, capture_output=True,
                            text=True, timeout=14400)
        p1.stdout.close()
        err1 = p1.stderr.read().decode(errors="replace")
        p1.wait()
        if p1.returncode != 0:
            warnings.append(f"{name}: {step} (bwa-mem2) failed: {err1[-200:]}")
            return False
        if p2.returncode != 0:
            warnings.append(f"{name}: {step} (sort) failed: {p2.stderr[-200:]}")
            return False
        return True
    except subprocess.TimeoutExpired:
        warnings.append(f"{name}: {step} timed out")
        return False


def _mock_result(samples: list, output_dir: Path, mapq: int,
                 remove_mito: bool, warnings: list) -> dict:
    bam_files = []
    for s in samples:
        bam_files.append({
            "name": s.get("name"),
            "status": "success",
            "bam": str(output_dir / s.get("name", "sample") /
                       f"{s.get('name', 'sample')}.filtered.bam"),
            "note": "mock — bwa-mem2 not installed",
        })
    return {
        "status": "success",
        "bam_files": bam_files,
        "n_aligned": len(bam_files),
        "output_dir": str(output_dir),
        "aligner": "bwa-mem2",
        "filter": {"mapq": mapq, "remove_mito": remove_mito},
        "warnings": warnings,
        "note": "mock — aligner not installed",
    }


if __name__ == "__main__":
    run_script(atac_align)
