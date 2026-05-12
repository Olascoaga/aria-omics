"""
ARIA RNA-seq Quantification Script
-------------------------------------
Runs featureCounts on aligned BAMs → counts matrix ready for DESeq2.
Executed inside aria-rnaseq-env via EnvironmentManager.

Input params:
  bam_files:   list  — from rna_align output
  gtf_file:    str   — GTF annotation
  output_dir:  str
  strand:      int   — 0=unstranded, 1=stranded, 2=reverse-stranded
                       (default: 0, auto-detect if "auto")
  threads:     int   (default: 8)
  feature:     str   (default: "gene") — exon/gene
  paired:      bool  (default: True)

Output:
  {
    "status":        "success",
    "counts_matrix": str,   — path to TSV (genes × samples)
    "summary":       str,   — featureCounts summary
    "n_genes":       int,
    "n_samples":     int,
    "strand_used":   int,
    "warnings":      [str]
  }
"""

from __future__ import annotations
import sys, os, subprocess, re
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script


def rna_quantify(params: dict) -> dict:
    bam_files  = params["bam_files"]         # from rna_align
    gtf_file   = params["gtf_file"]
    output_dir = Path(params.get("output_dir",
                                  str(Path(bam_files[0]["bam"]).parent.parent
                                      / "aria_counts")))
    threads    = int(params.get("threads", 8))
    strand     = params.get("strand", 0)     # 0/1/2 or "auto"
    paired     = bool(params.get("paired", True))
    # Default to "exon" — standard for RNA-seq.
    # Reads overlapping any exon of a gene contribute to that gene's count.
    # Using -t gene would also include intronic reads (DNA contamination noise).
    feature    = params.get("feature", "exon")
    allow_mock = mocks_allowed(params)
    warnings   = []

    output_dir.mkdir(parents=True, exist_ok=True)

    # Filter successful BAMs
    valid_bams = [b for b in bam_files
                  if b.get("status") == "success" and b.get("bam")
                  and Path(b["bam"]).exists()]

    if not valid_bams:
        # Check for mock BAMs (bam path doesn't exist but has "note": "mock")
        mock_bams = [b for b in bam_files if "mock" in str(b.get("note", ""))]
        if mock_bams and allow_mock:
            return _mock_counts_matrix(mock_bams, output_dir, warnings)
        return {
            "status":     "error",
            "error_type": "NoBamsFound",
            "details":    "No valid BAM files found from alignment step.",
        }

    # ── Auto-detect strandedness if requested ─────────────────────────────
    if strand == "auto":
        strand, strand_warn = _detect_strandedness(
            valid_bams[0]["bam"], gtf_file, warnings
        )
        warnings.extend(strand_warn)

    # ── Run featureCounts ─────────────────────────────────────────────────
    bam_paths    = [b["bam"] for b in valid_bams]
    sample_names = [b["name"] for b in valid_bams]
    counts_file  = str(output_dir / "counts.txt")
    matrix_file  = output_dir / "counts_matrix.tsv"
    summary_file = counts_file + ".summary"

    # ── Resume check: clean matrix + summary already present? ───────────
    # We check for the cleaned-up matrix (counts_matrix.tsv) because that
    # is the final deliverable. If it's present AND has all samples as
    # columns AND the BAMs are older than the matrix → skip featureCounts.
    if _counts_outputs_valid(matrix_file, sample_names, bam_paths):
        try:
            import pandas as _pd
            df = _pd.read_csv(matrix_file, sep="\t", index_col=0, nrows=5)
            n_genes_resumed = sum(1 for _ in open(matrix_file)) - 1
            warnings.append(
                f"[resume] featureCounts skipped "
                f"(valid counts_matrix.tsv exists: "
                f"{n_genes_resumed} genes × {len(sample_names)} samples)"
            )
            return {
                "status":        "success",
                "counts_matrix": str(matrix_file),
                "summary":       summary_file if Path(summary_file).exists() else None,
                "n_genes":       n_genes_resumed,
                "n_samples":     len(valid_bams),
                "strand_used":   strand,
                "sample_names":  sample_names,
                "warnings":      warnings,
                "resumed":       True,
            }
        except Exception:
            pass  # fall through to full run

    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", gtf_file,
        "-o", counts_file,
        "-g", "gene_id",
        "-t", feature,             # "exon" by default, "gene" if user forces
        "-s", str(strand),
    ]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += bam_paths

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )
        if proc.returncode != 0:
            # Try without --countReadPairs (older subread)
            if "--countReadPairs" in cmd:
                cmd.remove("--countReadPairs")
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=3600
                )

        if proc.returncode != 0:
            return {
                "status":     "error",
                "error_type": "featureCountsFailed",
                "details":    proc.stderr[-400:],
                "warnings":   warnings,
            }

        # ── Clean up counts matrix + annotate symbols ─────────────────
        matrix_path, n_genes = _clean_counts_matrix(
            counts_file, sample_names, output_dir, warnings,
            gtf_file=gtf_file,
        )

        symbols_path = output_dir / "counts_with_symbols.tsv"
        return {
            "status":        "success",
            "counts_matrix": str(matrix_path),
            "counts_with_symbols": (str(symbols_path)
                                      if symbols_path.exists() else None),
            "summary":       summary_file if Path(summary_file).exists() else None,
            "n_genes":       n_genes,
            "n_samples":     len(valid_bams),
            "strand_used":   strand,
            "sample_names":  sample_names,
            "warnings":      warnings,
        }

    except FileNotFoundError:
        if allow_mock:
            warnings.append(
                "featureCounts not found — using mock counts matrix "
                "(explicit mock mode)"
            )
            return _mock_counts_matrix(valid_bams, output_dir, warnings)
        return {
            "status":     "error",
            "error_type": "MissingDependency",
            "details":    "featureCounts is required for RNA-seq quantification.",
            "warnings":   warnings,
        }

    except subprocess.TimeoutExpired:
        return {
            "status":     "error",
            "error_type": "featureCountsTimeout",
            "details":    "featureCounts timed out (>1h)",
            "warnings":   warnings,
        }


def _counts_outputs_valid(matrix_file: Path, expected_samples: list,
                           bam_paths: list) -> bool:
    """
    Check if counts_matrix.tsv is present, complete, and newer than BAMs.
    """
    try:
        if not matrix_file.exists() or matrix_file.stat().st_size < 1024:
            return False
        with open(matrix_file) as f:
            header = f.readline().rstrip("\n")
        cols = header.split("\t")
        cols_set = set(cols)
        missing = [s for s in expected_samples if s not in cols_set]
        if missing:
            return False
        matrix_mtime = matrix_file.stat().st_mtime
        for b in bam_paths:
            bp = Path(b)
            if bp.exists() and bp.stat().st_mtime > matrix_mtime:
                return False
        return True
    except Exception:
        return False


def _clean_counts_matrix(counts_file: str, sample_names: list,
                           output_dir: Path, warnings: list,
                           gtf_file: str = None) -> tuple:
    """
    featureCounts output has extra columns (Chr, Start, End, Strand, Length).
    Clean to genes × samples TSV that DESeq2 expects.

    If gtf_file is provided, also annotate Ensembl IDs with gene symbols
    and write counts_with_symbols.tsv (used by pathway enrichment which
    needs HGNC symbols, not ENSG IDs).
    """
    import pandas as pd

    df = pd.read_csv(counts_file, sep="\t", comment="#")

    # First 6 columns: Geneid, Chr, Start, End, Strand, Length
    count_cols = df.columns[6:].tolist()

    counts = df[["Geneid"] + count_cols].copy()
    counts.columns = ["gene_id"] + sample_names
    counts = counts.set_index("gene_id")

    # Remove zero-count genes
    n_before = len(counts)
    counts   = counts[counts.sum(axis=1) > 0]
    n_after  = len(counts)

    if n_before - n_after > 0:
        warnings.append(
            f"{n_before - n_after} genes with zero counts across all "
            f"samples removed."
        )

    # Standard output: counts_matrix.tsv (matches resume check filename)
    out_path = output_dir / "counts_matrix.tsv"
    counts.to_csv(str(out_path), sep="\t")

    # ── Gene symbol annotation (for pathway enrichment) ──────────────────
    # Enrichr/gseapy require HGNC symbols, not Ensembl IDs.
    if gtf_file and Path(gtf_file).exists():
        symbol_map = _build_ensembl_to_symbol_map(gtf_file, warnings)
        if symbol_map:
            n_mapped = sum(1 for gid in counts.index if gid in symbol_map)
            symbols  = pd.Series(
                [symbol_map.get(gid, gid) for gid in counts.index],
                index=counts.index,
                name="gene_symbol",
            )
            counts_sym = counts.copy()
            counts_sym.insert(0, "gene_symbol", symbols)

            sym_path = output_dir / "counts_with_symbols.tsv"
            counts_sym.to_csv(str(sym_path), sep="\t")

            warnings.append(
                f"Annotated {n_mapped:,}/{len(counts):,} genes with HGNC "
                f"symbols from GTF."
            )

    return out_path, len(counts)


def _build_ensembl_to_symbol_map(gtf_file: str, warnings: list) -> dict:
    """
    Parse a GTF to build {ensembl_id: gene_symbol}.

    Handles both gzipped and plain GTF. Reads only `gene` feature lines
    (one per gene, fast even for full Ensembl GTF on hg38).

    Returns empty dict on failure (graceful — pathway enrichment will
    fall back to using Ensembl IDs and report no matches).
    """
    import gzip as _gzip
    import re as _re

    mapping = {}
    opener = _gzip.open if str(gtf_file).endswith(".gz") else open

    try:
        with opener(gtf_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attrs = fields[8]
                gid_match = _re.search(r'gene_id "([^"]+)"', attrs)
                sym_match = _re.search(r'gene_name "([^"]+)"', attrs)
                if gid_match and sym_match:
                    # Strip Ensembl version suffix (ENSG00000123.4 → ENSG00000123)
                    gid = gid_match.group(1).split(".")[0]
                    mapping[gid] = sym_match.group(1)
    except Exception as e:
        warnings.append(f"GTF symbol mapping failed: {e}")
        return {}

    return mapping


def _detect_strandedness(bam: str, gtf: str,
                          warnings: list) -> tuple[int, list]:
    """
    Detect library strandedness by sampling reads and comparing to gene strand.

    Algorithm:
      1. Load gene coordinates from GTF (first ~5000 protein-coding genes is enough).
      2. Sample reads from the BAM (up to 200000 reads from chr1).
      3. For each read overlapping a gene, count whether it aligns in the
         "sense" direction (same as gene strand) or "antisense".
      4. Compute the sense fraction:
           ~0.50  → unstranded     (return 0)
           >0.80  → forward (FR)   (return 1) — Ligation kits like NEB Ultra II Directional
           <0.20  → reverse (RF)   (return 2) — TruSeq Stranded mRNA (most common)

    Falls back to 0 (unstranded) if pysam not available or no overlap found.
    """
    try:
        import pysam
    except ImportError:
        warnings.append(
            "pysam not available — cannot auto-detect strandedness. "
            "Defaulting to unstranded. Install pysam in aria-rnaseq-env."
        )
        return 0, []

    if not Path(bam).exists():
        warnings.append(f"BAM not found for strand detection: {bam}")
        return 0, []
    if not Path(gtf).exists():
        warnings.append(f"GTF not found for strand detection: {gtf}")
        return 0, []

    # ── 1. Load gene intervals from GTF ────────────────────────────────────
    # We only need ~5000 genes worth of intervals to get a stable signal.
    # Index by chromosome → list of (start, end, strand) intervals.
    import gzip as _gzip
    chrom_genes: dict[str, list] = {}
    n_genes = 0
    GENE_TARGET = 5000

    opener = _gzip.open if str(gtf).endswith(".gz") else open
    try:
        with opener(gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if n_genes >= GENE_TARGET:
                    break
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                chrom  = fields[0]
                start  = int(fields[3]) - 1   # GTF is 1-based, BAM is 0-based
                end    = int(fields[4])
                strand = fields[6]            # "+" or "-"
                if strand not in ("+", "-"):
                    continue
                chrom_genes.setdefault(chrom, []).append((start, end, strand))
                n_genes += 1
    except Exception as e:
        warnings.append(f"GTF parse failed for strand detection: {e}")
        return 0, []

    if not chrom_genes:
        warnings.append("No genes parsed from GTF — using unstranded default.")
        return 0, []

    # Sort intervals per chromosome (binary search later)
    for chrom in chrom_genes:
        chrom_genes[chrom].sort()

    # ── 2. Sample reads from the BAM ───────────────────────────────────────
    # We focus on chr1 (or first available chromosome with genes) to keep IO low.
    # For paired-end, count only read1 to avoid double-counting.
    READ_TARGET     = 200_000
    sense, anti     = 0, 0
    reads_examined  = 0

    try:
        bam_obj = pysam.AlignmentFile(bam, "rb")
    except Exception as e:
        warnings.append(f"Cannot open BAM for strand detection: {e}")
        return 0, []

    # Find a chromosome that exists in both BAM and GTF
    bam_chroms = set(bam_obj.references)
    target_chrom = None
    for candidate in ["chr1", "1", "chrI"]:
        if candidate in bam_chroms and candidate in chrom_genes:
            target_chrom = candidate
            break
    if target_chrom is None:
        # Try any chromosome with genes in the BAM
        for c in sorted(bam_chroms):
            if c in chrom_genes:
                target_chrom = c
                break

    if target_chrom is None:
        warnings.append(
            "No matching chromosome between BAM and GTF for strand detection. "
            "Defaulting to unstranded."
        )
        bam_obj.close()
        return 0, []

    intervals = chrom_genes[target_chrom]

    # Extract just starts for binary search
    starts = [iv[0] for iv in intervals]
    import bisect

    try:
        for read in bam_obj.fetch(target_chrom):
            if reads_examined >= READ_TARGET:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_paired and read.is_read2:
                continue   # count only R1 for paired-end
            if read.mapping_quality < 30:
                continue   # high-quality alignments only

            r_start = read.reference_start
            r_end   = read.reference_end if read.reference_end else r_start + 50

            # Binary search: find genes that could overlap this read
            idx = bisect.bisect_right(starts, r_end) - 1
            if idx < 0:
                continue

            # Walk back briefly looking for an overlapping gene
            best_strand = None
            for j in range(max(0, idx - 5), idx + 1):
                g_start, g_end, g_strand = intervals[j]
                if g_start <= r_end and g_end >= r_start:
                    best_strand = g_strand
                    break
            if best_strand is None:
                continue

            reads_examined += 1
            # Sense = read aligns in same direction as gene transcription
            read_on_plus  = not read.is_reverse
            gene_on_plus  = (best_strand == "+")
            if read_on_plus == gene_on_plus:
                sense += 1
            else:
                anti += 1
    except Exception as e:
        warnings.append(f"Error reading BAM for strand detection: {e}")
    finally:
        bam_obj.close()

    total = sense + anti
    if total < 1000:
        warnings.append(
            f"Strand detection: only {total} informative reads found "
            f"on {target_chrom}. Defaulting to unstranded."
        )
        return 0, []

    sense_frac = sense / total

    # ── 3. Decide ──────────────────────────────────────────────────────────
    # Thresholds based on standard library prep behavior.
    # Real datasets very rarely fall in the ambiguous middle band.
    if sense_frac >= 0.80:
        warnings.append(
            f"Strandedness auto-detected: forward stranded "
            f"({sense_frac*100:.1f}% sense reads, n={total}). "
            f"Using featureCounts -s 1."
        )
        return 1, []
    elif sense_frac <= 0.20:
        warnings.append(
            f"Strandedness auto-detected: reverse stranded "
            f"({sense_frac*100:.1f}% sense reads, n={total}). "
            f"Using featureCounts -s 2 (typical for TruSeq Stranded mRNA)."
        )
        return 2, []
    elif 0.40 <= sense_frac <= 0.60:
        warnings.append(
            f"Strandedness auto-detected: unstranded "
            f"({sense_frac*100:.1f}% sense reads, n={total})."
        )
        return 0, []
    else:
        # Ambiguous — flag for user attention but pick the closer side
        if sense_frac > 0.5:
            inferred = 1
            label    = "forward"
        else:
            inferred = 2
            label    = "reverse"
        warnings.append(
            f"Strandedness ambiguous ({sense_frac*100:.1f}% sense reads, "
            f"n={total}). Tentatively using -s {inferred} ({label}). "
            f"VERIFY this matches your library prep — wrong strand = "
            f"~50% read loss."
        )
        return inferred, []


def _mock_counts_matrix(bam_files: list, output_dir: Path,
                         warnings: list) -> dict:
    """Generate mock counts matrix when featureCounts not available."""
    import numpy as np
    import pandas as pd

    rng      = np.random.default_rng(42)
    n_genes  = 20000
    samples  = [b["name"] for b in bam_files]
    genes    = [f"ENSMUSG{i:011d}" for i in range(n_genes)]

    # Simulate realistic NB counts with ~500 DE genes
    base     = rng.negative_binomial(20, 0.3, (n_genes, len(samples)))
    de_idx   = list(range(500))

    # Samples from groups other than first get 2-8x expression for DE genes
    n_grp    = len(samples) // 2
    for i in de_idx[:250]:
        base[i, n_grp:] = base[i, n_grp:] * rng.integers(2, 8)
    for i in de_idx[250:]:
        base[i, :n_grp] = base[i, :n_grp] * rng.integers(2, 8)

    df = pd.DataFrame(base, index=genes, columns=samples)
    out_path = output_dir / "counts_clean.tsv"
    df.to_csv(str(out_path), sep="\t")

    warnings.append(
        "MOCK counts matrix generated — featureCounts not available. "
        "Install aria-rnaseq-env for real quantification."
    )

    return {
        "status":        "success",
        "counts_matrix": str(out_path),
        "n_genes":       n_genes,
        "n_samples":     len(samples),
        "strand_used":   0,
        "sample_names":  samples,
        "warnings":      warnings,
        "note":          "mock",
    }


if __name__ == "__main__":
    run_script(rna_quantify)
