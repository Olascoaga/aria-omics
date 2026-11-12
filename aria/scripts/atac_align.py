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
  - for PE only, keep properly paired alignments (-f 2); drop unmapped/mate-unmapped/
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
  mito_contigs:   list  (default: ["chrM", "MT", "M"])

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
try:
    from aria.utils.stage_manifest import (
        file_record, load_valid_stage_result, stage_manifest_path, write_stage_manifest,
    )
except ImportError:
    from aria.scripts.chromatin_peak_counts import (
        file_record, load_valid_stage_result, stage_manifest_path, write_stage_manifest,
    )
from aria.scripts.chromatin_peak_counts import _save_stage_result

# Canonical ENCODE ATAC read filter: drop unmapped (0x4), mate-unmapped (0x8),
# secondary (0x100), QC-fail (0x200), duplicate (0x400). NOTE: 1804 does NOT
# include supplementary (0x800) — this matches the ENCODE ATAC-seq spec, which
# retains supplementary alignments at this step (they are a tiny fraction and
# are gated by the -f 2 proper-pair requirement).
SAMTOOLS_EXCLUDE_FLAGS = 1804
SAMTOOLS_REQUIRE_FLAGS = 2  # properly paired
DEFAULT_MIN_MAPQ = 30
DEFAULT_MAPQ = DEFAULT_MIN_MAPQ
DEFAULT_MITO = frozenset({"chrM", "MT", "M"})
DEFAULT_MITO_CONTIGS = tuple(sorted(DEFAULT_MITO))


def atac_align(params: dict) -> dict:
    samples = params.get("samples") or []
    genome_fasta = params.get("genome_fasta", "")
    bwa_index = params.get("bwa_index") or genome_fasta
    output_dir = Path(params.get("output_dir",
                                 str(Path("aria_atac_aligned"))))
    threads = int(params.get("threads", 8))
    mapq = int(params.get("min_mapq", params.get("mapq", DEFAULT_MIN_MAPQ)))
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

    index_built = False

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
        index_built = True

    # ── Align + filter each sample ────────────────────────────────────────
    bam_files = []
    for sample in samples:
        if sample.get("status") == "failed":
            warnings.append(f"Skipping {sample.get('name')} — QC failed upstream")
            continue
        sid = str(sample.get("sample_id") or sample.get("name") or "sample")
        r1 = sample.get("r1_trimmed") or sample.get("r1")
        r2 = sample.get("r2_trimmed") or sample.get("r2")
        is_paired = bool(sample.get("paired", bool(r2)) and r2)
        if not r1 or not Path(r1).is_file() or (is_paired and not Path(r2).is_file()):
            bam_files.append({"name": sid, "sample_id": sid, "status": "failed",
                              "bam": None, "paired": is_paired, "error_type": "FileNotFound",
                              "reused_from_stage_manifest": False,
                              "mapped_reads_post_filter": None})
            warnings.append(f"{sid}: input FASTQ missing on disk")
            continue
        inputs = [file_record(r1)]
        if is_paired:
            inputs.append(file_record(r2))
        if genome_fasta:
            inputs.append(file_record(genome_fasta))
        mito = set(mito_contigs) if remove_mito else set()
        parameters = {"sample_id": sid, "paired": is_paired, "min_mapq": mapq,
                      "mito_contigs": sorted(mito), "bwa_index": bwa_index}
        manifest = stage_manifest_path(output_dir, f"atac_align_{sid}")
        cached = load_valid_stage_result(
            manifest, expected_inputs=inputs, expected_parameters=parameters)
        if (cached and Path(cached.get("bam", "")).is_file()
                and Path(cached["bam"] + ".bai").is_file()):
            result = {**cached, "reused_from_stage_manifest": True}
        else:
            Path(manifest).unlink(missing_ok=True)
            result = _align_sample(
                sample={**sample, "name": sid}, bwa_index=bwa_index,
                output_dir=output_dir, threads=threads, mapq=mapq,
                remove_mito=remove_mito, mito_contigs=tuple(mito),
                warnings=warnings, allow_mock=allow_mock,
            )
            result.update({"sample_id": sid, "paired": is_paired,
                           "mapped_reads_post_filter": result.get("n_final"),
                           "reused_from_stage_manifest": False})
            if (result.get("status") == "success" and result.get("bam")
                    and Path(result["bam"]).is_file()
                    and Path(result["bam"] + ".bai").is_file()):
                _save_stage_result(manifest, "atac_align", inputs, parameters, result)
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

    layouts = {b["paired"] for b in bam_files if b.get("status") == "success"}
    all_paired = layouts == {True}
    return {
        "status": "success",
        "validation_level": "production",
        "per_sample": bam_files,
        "read_layout": "mixed" if len(layouts) > 1 else "paired" if all_paired else "single",
        "index_built": index_built,
        "filters": {"min_mapq": mapq,
                    "samtools_f": SAMTOOLS_REQUIRE_FLAGS if all_paired else 0,
                    "samtools_F": SAMTOOLS_EXCLUDE_FLAGS,
                    "markdup_removed": True,
                    "mito_stripped": sorted(mito_contigs) if remove_mito else []},
        "bam_files": bam_files,
        "n_aligned": n_ok,
        "output_dir": str(output_dir),
        "aligner": "bwa-mem2",
        "filter": {
            "mapq": mapq,
            "samtools_F": SAMTOOLS_EXCLUDE_FLAGS,
            "samtools_f": SAMTOOLS_REQUIRE_FLAGS if all_paired else 0,
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

    Resume is governed by the caller's content-addressed A4 manifest.
    """
    name = sample["name"]
    r1 = sample.get("r1_trimmed") or sample.get("r1")
    r2 = sample.get("r2_trimmed") or sample.get("r2")
    paired = bool(sample.get("paired", bool(r2)) and r2)

    sample_dir = output_dir / name
    sample_dir.mkdir(exist_ok=True)
    final_bam = sample_dir / f"{name}.filtered.bam"

    result = {"name": name, "status": "pending", "paired": paired, "bam": None}

    if not r1:
        warnings.append(f"{name}: no R1 FASTQ — skipping")
        result["status"] = "failed"
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
        sort_cmd = ["samtools", "sort", "-@", str(threads)]
        if paired:
            sort_cmd += ["-n"]
        sort_cmd += ["-o", str(namesort if paired else coordsort), "-"]
        if not _pipe(align_cmd, sort_cmd, warnings, name, "align"):
            result["status"] = "failed"
            return result

        if paired:
            if not _run(["samtools", "fixmate", "-m", "-@", str(threads),
                         str(namesort), str(fixmate)], warnings, name, "fixmate"):
                result["status"] = "failed"
                return result
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
        # Filter through SAM text so mitochondrial stripping never depends on an index.
        view_cmd = ["samtools", "view", "-h", "-q", str(mapq)]
        if paired:
            view_cmd += ["-f", str(SAMTOOLS_REQUIRE_FLAGS)]
        view_cmd += ["-F", str(SAMTOOLS_EXCLUDE_FLAGS), str(markdup)]
        if not _filter_bam(view_cmd, final_bam, mito_contigs if remove_mito else (),
                           warnings, name):
            result["status"] = "failed"
            return result
        if not _run(["samtools", "index", str(final_bam)], warnings, name, "index_final"):
            result["status"] = "failed"
            return result

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


def _strip_mito_sam(lines, mito=DEFAULT_MITO):
    """Yield SAM headers and alignments outside the explicit mitochondrial set."""
    mito = set(mito)
    for line in lines:
        fields = line.split("\t", 3)
        if line.startswith("@") or (len(fields) >= 3 and fields[2] not in mito):
            yield line


def _filter_bam(view_cmd, out_bam, mito, warnings, name):
    """Stream SAM through the mitochondrial filter, checking both processes."""
    import tempfile
    with tempfile.TemporaryFile(mode="w+t") as err1, tempfile.TemporaryFile(mode="w+t") as err2:
        with subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=err1, text=True) as src:
            with subprocess.Popen(["samtools", "view", "-b", "-o", str(out_bam), "-"],
                                  stdin=subprocess.PIPE, stderr=err2, text=True) as dst:
                try:
                    for line in _strip_mito_sam(src.stdout, mito):
                        dst.stdin.write(line)
                except BrokenPipeError:
                    pass
                finally:
                    dst.stdin.close()
                    src.stdout.close()
                dst.wait()
            src.wait()
        if src.returncode or dst.returncode:
            err1.seek(0)
            err2.seek(0)
            warnings.append(f"{name}: filter failed: {(err1.read() + err2.read())[-300:]}")
            return False
    return True


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
        if not bam.is_file() or bam.stat().st_size == 0:
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
    """Run cmd1 | cmd2 without blocking on a full aligner stderr pipe."""
    import tempfile
    with tempfile.TemporaryFile() as stderr:
        with subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=stderr) as src:
            try:
                proc = subprocess.run(cmd2, stdin=src.stdout, capture_output=True,
                                      text=True, timeout=14400)
            except subprocess.TimeoutExpired:
                src.kill()
                src.wait()
                warnings.append(f"{name}: {step} timed out")
                return False
            finally:
                src.stdout.close()
            src.wait()
        if src.returncode or proc.returncode:
            stderr.seek(0)
            detail = stderr.read().decode(errors="replace") + (proc.stderr or "")
            warnings.append(f"{name}: {step} failed: {detail[-300:]}")
            return False
    return True


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
