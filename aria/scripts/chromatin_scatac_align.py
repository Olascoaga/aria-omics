"""
ARIA scATAC FASTQ -> fragments alignment (chromap)
---------------------------------------------------
Aligns raw single-cell ATAC FASTQ (10x-style: genomic R1/R3 + a barcode read)
to a reference with chromap and emits the barcoded, position-sorted
`fragments.tsv.gz` (+ tabix index) that the fragments->matrix bridge
(chromatin_fragments_to_matrix.py) consumes. Together they take scATAC from raw
sequencer FASTQ all the way to the validated matrix pipeline.

Executed inside aria-atacseq-env (chromap + samtools) via EnvironmentManager.

chromap is the open-source single-cell ATAC aligner (it handles the cell barcode
+ whitelist correction + Tn5 fragment generation in one pass); this is the
single-cell analogue of atac_align.py (bulk, bwa-mem2).

Input params:
  r1_fastq:        str  — genomic read 1 (10x ATAC R1)
  r3_fastq:        str  — genomic read 2 (10x ATAC R3; "r2_fastq" also accepted
                          for the genomic mate when no separate barcode key)
  barcode_fastq:   str  — cell-barcode read (10x ATAC R2)
  genome_fasta:    str  — reference FASTA
  chromap_index:   str  — (optional) prebuilt chromap index; default <fasta>.index
  barcode_whitelist: str — 10x barcode whitelist (plain or .gz)
  output_dir:      str
  threads:         int  (default 8)

Output:
  {
    "status": "success",
    "fragments_file": str,    — bgzipped, tabix-indexed fragments.tsv.gz
    "aligner": "chromap",
    "n_fragments": int,
    "warnings": [str]
  }
  Missing chromap/samtools or required inputs -> structured not-run (never a
  fabricated fragments file).
"""

from __future__ import annotations
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def chromatin_scatac_align(params: dict) -> dict:
    r1 = params.get("r1_fastq")
    r3 = params.get("r3_fastq") or params.get("r2_fastq")
    barcode = params.get("barcode_fastq")
    genome_fasta = params.get("genome_fasta", "")
    chromap_index = params.get("chromap_index") or (
        f"{genome_fasta}.index" if genome_fasta else "")
    whitelist = params.get("barcode_whitelist")
    output_dir = Path(params.get("output_dir") or "aria_scatac_aligned")
    threads = int(params.get("threads", 8))
    warnings: list[str] = []

    def not_run(reason: str, **extra) -> dict:
        return {"status": "skipped", "ran": False,
                "analysis": "scatac_fastq_to_fragments", "reason": reason,
                "warnings": warnings, **extra}

    # Required inputs (honest, before touching any tool).
    for key, val in (("r1_fastq", r1), ("genomic mate (r3/r2)", r3),
                     ("barcode_fastq", barcode)):
        if not val or not Path(val).exists():
            return not_run("missing_input",
                           message=f"required scATAC FASTQ not found: {key}={val}")
    if not genome_fasta or not Path(genome_fasta).exists():
        return not_run("missing_genome_fasta",
                       message=f"genome_fasta not found: {genome_fasta}")
    if not whitelist or not Path(whitelist).exists():
        return not_run("missing_barcode_whitelist",
                       message=("a 10x barcode whitelist is required for "
                                "single-cell barcode correction"))

    output_dir.mkdir(parents=True, exist_ok=True)
    fragments_raw = output_dir / "fragments.tsv"
    fragments_gz = output_dir / "fragments.tsv.gz"

    if _fragments_valid(fragments_gz):
        warnings.append("[resume] reused existing fragments.tsv.gz")
        return {"status": "success", "fragments_file": str(fragments_gz),
                "aligner": "chromap", "resumed": True, "warnings": warnings}

    # 1. Build the chromap index if missing.
    if not Path(chromap_index).exists():
        idx = _run(["chromap", "-i", "-r", genome_fasta, "-o", chromap_index],
                   warnings, "index")
        if idx == "missing_tool":
            return not_run("chromap_unavailable",
                           message="chromap not found (install aria-atacseq-env)")
        if idx != "ok":
            return not_run("index_build_failed", message=idx)

    # 2. Align: chromap --preset atac -> fragments TSV (barcode-aware).
    cmd = [
        "chromap", "--preset", "atac",
        "-x", chromap_index, "-r", genome_fasta,
        "-1", r1, "-2", r3, "-b", barcode,
        "--barcode-whitelist", str(whitelist),
        "-o", str(fragments_raw),
        "-t", str(threads),
    ]
    res = _run(cmd, warnings, "align")
    if res == "missing_tool":
        return not_run("chromap_unavailable",
                       message="chromap not found (install aria-atacseq-env)")
    if res != "ok":
        return not_run("alignment_failed", message=res)
    if not fragments_raw.exists() or fragments_raw.stat().st_size == 0:
        return not_run("no_fragments_produced",
                       message="chromap produced no fragments")

    # 3. bgzip + tabix so the fragments are random-access for the bridge.
    bg = _bgzip_tabix(fragments_raw, fragments_gz, warnings)
    if bg != "ok":
        # The plain fragments still exist; the bridge can read uncompressed.
        warnings.append(f"bgzip/tabix unavailable ({bg}); leaving plain TSV")
        n = _count_lines(fragments_raw)
        return {"status": "success", "fragments_file": str(fragments_raw),
                "aligner": "chromap", "n_fragments": n,
                "compressed": False, "warnings": warnings}

    n = _count_gz_lines(fragments_gz)
    return {
        "status": "success",
        "fragments_file": str(fragments_gz),
        "aligner": "chromap",
        "n_fragments": n,
        "compressed": True,
        "tabix_indexed": Path(str(fragments_gz) + ".tbi").exists(),
        "genome_fasta": genome_fasta,
        "validation_level": "beta",
        "warnings": warnings,
    }


def _bgzip_tabix(plain: Path, gz: Path, warnings: list) -> str:
    try:
        with open(gz, "wb") as out:
            p = subprocess.run(["bgzip", "-c", str(plain)], stdout=out,
                               stderr=subprocess.PIPE, timeout=3600)
        if p.returncode != 0:
            return p.stderr.decode(errors="replace")[-200:]
        t = subprocess.run(["tabix", "-p", "bed", str(gz)],
                           capture_output=True, text=True, timeout=1800)
        if t.returncode != 0:
            return t.stderr[-200:]
        try:
            plain.unlink()
        except OSError:
            pass
        return "ok"
    except FileNotFoundError:
        return "missing_tool"
    except subprocess.TimeoutExpired:
        return "bgzip/tabix timed out"


def _run(cmd: list, warnings: list, step: str) -> str:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
        if p.returncode != 0:
            warnings.append(f"{step} failed: {p.stderr[-200:]}")
            return p.stderr[-300:]
        return "ok"
    except FileNotFoundError:
        return "missing_tool"
    except subprocess.TimeoutExpired:
        warnings.append(f"{step} timed out")
        return f"{step} timed out"


def _fragments_valid(gz: Path) -> bool:
    try:
        return gz.exists() and gz.stat().st_size > 1024
    except OSError:
        return False


def _count_lines(path: Path) -> int:
    try:
        with open(path) as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


def _count_gz_lines(path: Path) -> int:
    import gzip
    try:
        with gzip.open(path, "rt") as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


if __name__ == "__main__":
    run_script(chromatin_scatac_align)
