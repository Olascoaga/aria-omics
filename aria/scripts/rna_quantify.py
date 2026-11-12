"""
ARIA RNA-seq Quantification Script (ADR-075)
--------------------------------------------
Runs featureCounts on aligned BAMs -> counts matrix ready for DESeq2.
Executed inside ``aria-rnaseq-env`` via ``EnvironmentManager``.

Features:
  - Content-addressed A4 stage resume via ``counts_matrix.stage.json``
  - Fail-closed multi-contig ``pysam`` library strandedness auto-detection
    (never silently falls back to ``-s 0`` when ``strand="auto"`` fails or
    hits an ambiguous 20-40% / 60-80% sense band)
  - Structured ``assignment_summary`` parsed from ``counts.txt.summary``
  - Ensembl -> HGNC gene symbol sidecar table (``counts_with_symbols.tsv``)
"""

from __future__ import annotations

import bisect
import gzip
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script
from aria.utils.stage_manifest import (
    build_stage_manifest,
    stage_is_current,
    write_stage_manifest,
)

_COUNTS_MANIFEST = "counts_matrix.stage.json"
_COUNTS_METRICS = "counts_matrix.metrics.json"


def _featurecounts_version() -> str:
    """Best-effort featureCounts version string for A4 stage manifests."""
    try:
        proc = subprocess.run(
            ["featureCounts", "-v"],
            capture_output=True,
            text=True,
            timeout=15,
        )
        raw = (proc.stdout or proc.stderr or "").strip()
        for line in raw.splitlines():
            s = line.strip()
            if s:
                return s
        return ""
    except Exception:
        return ""


def rna_quantify(params: dict) -> dict:
    bam_files = params["bam_files"]
    gtf_file = params["gtf_file"]
    output_dir = Path(
        params.get(
            "output_dir",
            str(Path(bam_files[0]["bam"]).parent.parent / "aria_counts"),
        )
    )
    threads = int(params.get("threads", 8))
    strand_in = params.get("strand", 0)  # 0/1/2 or "auto"
    paired_param = params.get("paired")
    min_strand_reads = int(params.get("min_strand_reads", 50))
    # Default to "exon" — standard for RNA-seq (avoids intronic gDNA noise).
    feature = str(params.get("feature", "exon"))
    allow_mock = mocks_allowed(params)
    warnings: list[str] = []

    output_dir.mkdir(parents=True, exist_ok=True)

    failed_bams = [
        b.get("name", "unknown")
        for b in bam_files
        if b.get("status") != "success"
    ]
    if failed_bams:
        return {
            "status": "error",
            "error_type": "PartialAlignmentInput",
            "details": (
                "Quantification is all-or-fail and received unsuccessful "
                "alignment records."
            ),
            "failed_samples": failed_bams,
        }

    valid_bams = [
        b
        for b in bam_files
        if b.get("status") == "success"
        and b.get("bam")
        and Path(b["bam"]).exists()
    ]

    missing_bams = [
        b.get("name", "unknown")
        for b in bam_files
        if b.get("status") == "success"
        and not (b.get("bam") and Path(b["bam"]).exists())
        and "mock" not in str(b.get("note", ""))
    ]
    if missing_bams:
        return {
            "status": "error",
            "error_type": "MissingBamInputs",
            "details": (
                "Quantification is all-or-fail and aligned BAMs are missing "
                f"for: {', '.join(missing_bams)}."
            ),
            "failed_samples": missing_bams,
        }

    if not valid_bams:
        mock_bams = [b for b in bam_files if "mock" in str(b.get("note", ""))]
        if mock_bams and allow_mock:
            mock_layout = _read_layout_for_bams(mock_bams, paired_param)
            if mock_layout == "mixed":
                return {
                    "status": "error",
                    "error_type": "MixedReadLayouts",
                    "details": "Mock BAM records contain mixed read layouts.",
                }
            return _mock_counts_matrix(
                mock_bams, output_dir, warnings, read_layout=mock_layout
            )
        return {
            "status": "error",
            "error_type": "NoBamsFound",
            "details": "No valid BAM files found from alignment step.",
        }

    layouts = {
        b.get(
            "read_layout",
            "paired-end" if b.get("paired", paired_param is not False) else "single-end",
        )
        for b in valid_bams
    }
    if not layouts <= {"paired-end", "single-end"} or len(layouts) != 1:
        return {
            "status": "error",
            "error_type": "MixedReadLayouts",
            "details": (
                "featureCounts requires one explicit read layout per run; "
                f"received {sorted(layouts)}."
            ),
        }
    read_layout = next(iter(layouts))
    paired = read_layout == "paired-end"

    bam_paths = [str(b["bam"]) for b in valid_bams]
    sample_names = [str(b["name"]) for b in valid_bams]

    # ── Resolve or auto-detect library strandedness (fail-closed) ─────────
    if isinstance(strand_in, str) and strand_in.strip().lower() == "auto":
        strand, strand_err = _detect_strandedness(
            bam_paths,
            gtf_file,
            warnings,
            min_informative_reads=min_strand_reads,
        )
        if strand_err is not None or strand is None:
            err_payload = dict(strand_err or {
                "status": "error",
                "error_type": "StrandDetectionFailed",
                "details": "Library strandedness auto-detection failed.",
            })
            err_payload["warnings"] = warnings
            return err_payload
    else:
        try:
            strand = int(strand_in)
        except (TypeError, ValueError):
            strand = -1
        if strand not in (0, 1, 2):
            return {
                "status": "error",
                "error_type": "InvalidStrandParameter",
                "details": (
                    f"strand must be 0 (unstranded), 1 (forward), 2 (reverse), "
                    f"or 'auto'; got {strand_in!r}."
                ),
                "warnings": warnings,
            }

    # ── Content-addressed A4 resume check ─────────────────────────────────
    counts_file = str(output_dir / "counts.txt")
    matrix_file = output_dir / "counts_matrix.tsv"
    symbols_file = output_dir / "counts_with_symbols.tsv"
    summary_file = counts_file + ".summary"
    manifest_path = output_dir / _COUNTS_MANIFEST
    metrics_path = output_dir / _COUNTS_METRICS

    stage_inputs: list[tuple[str, str | None]] = [
        ("gtf", str(gtf_file) if gtf_file else None),
        *[(f"bam_{idx}", p) for idx, p in enumerate(bam_paths, start=1)],
    ]
    stage_params = {
        "strand": int(strand),
        "feature": feature,
        "read_layout": read_layout,
        "sample_names": sample_names,
    }

    if manifest_path.exists() and _counts_outputs_valid(matrix_file, sample_names):
        fc_ver = _featurecounts_version()
        current, _ = stage_is_current(
            str(manifest_path),
            inputs=stage_inputs,
            params=stage_params,
            tool_version=fc_ver,
        )
        if current:
            n_genes_resumed = max(0, sum(1 for _ in open(matrix_file, encoding="utf-8")) - 1)
            cached_summary = None
            if metrics_path.exists():
                try:
                    cached = json.loads(metrics_path.read_text(encoding="utf-8"))
                    cached_summary = cached.get("assignment_summary")
                except Exception:
                    cached_summary = None
            warnings.append(
                f"[resume] featureCounts skipped "
                f"(valid counts_matrix.tsv + A4 manifest: "
                f"{n_genes_resumed} genes × {len(sample_names)} samples)"
            )
            return {
                "status": "success",
                "counts_matrix": str(matrix_file),
                "counts_with_symbols": (
                    str(symbols_file) if symbols_file.exists() else None
                ),
                "summary": summary_file if Path(summary_file).exists() else None,
                "assignment_summary": cached_summary,
                "n_genes": n_genes_resumed,
                "n_samples": len(valid_bams),
                "strand_used": int(strand),
                "read_layout": read_layout,
                "sample_names": sample_names,
                "warnings": warnings,
                "resumed": True,
            }

    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", str(gtf_file),
        "-o", counts_file,
        "-g", "gene_id",
        "-t", feature,
        "-s", str(strand),
    ]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += bam_paths

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )
        if proc.returncode != 0 and "--countReadPairs" in cmd:
            # Older subread fallback without --countReadPairs
            cmd.remove("--countReadPairs")
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600
            )

        if proc.returncode != 0:
            return {
                "status": "error",
                "error_type": "featureCountsFailed",
                "details": (proc.stderr or proc.stdout or "featureCounts failed")[-400:],
                "warnings": warnings,
            }

        matrix_path, n_genes = _clean_counts_matrix(
            counts_file,
            sample_names,
            output_dir,
            warnings,
            gtf_file=gtf_file,
        )
        assignment_summary = _parse_featurecounts_summary(
            summary_file, sample_names, warnings
        )

        if Path(matrix_path).exists() and Path(gtf_file).exists():
            fc_ver = _featurecounts_version()
            metrics_path.write_text(
                json.dumps(
                    {
                        "n_genes": n_genes,
                        "n_samples": len(valid_bams),
                        "strand_used": int(strand),
                        "read_layout": read_layout,
                        "assignment_summary": assignment_summary,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            write_stage_manifest(
                str(manifest_path),
                build_stage_manifest(
                    stage="featurecounts",
                    inputs=stage_inputs,
                    params=stage_params,
                    tool_version=fc_ver,
                ),
            )

        return {
            "status": "success",
            "counts_matrix": str(matrix_path),
            "counts_with_symbols": (
                str(symbols_file) if symbols_file.exists() else None
            ),
            "summary": summary_file if Path(summary_file).exists() else None,
            "assignment_summary": assignment_summary,
            "n_genes": n_genes,
            "n_samples": len(valid_bams),
            "strand_used": int(strand),
            "read_layout": read_layout,
            "sample_names": sample_names,
            "warnings": warnings,
        }

    except FileNotFoundError:
        if allow_mock:
            warnings.append(
                "featureCounts not found — using mock counts matrix "
                "(explicit mock mode)"
            )
            return _mock_counts_matrix(
                valid_bams, output_dir, warnings, read_layout=read_layout
            )
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": "featureCounts is required for RNA-seq quantification.",
            "warnings": warnings,
        }

    except subprocess.TimeoutExpired:
        return {
            "status": "error",
            "error_type": "featureCountsTimeout",
            "details": "featureCounts timed out (>1h)",
            "warnings": warnings,
        }


def _counts_outputs_valid(matrix_file: Path, expected_samples: list[str]) -> bool:
    """Verify counts_matrix.tsv exists, is non-empty, and carries all sample columns."""
    try:
        if not matrix_file.exists() or matrix_file.stat().st_size < 16:
            return False
        with open(matrix_file, encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n")
        cols_set = set(header.split("\t"))
        return all(s in cols_set for s in expected_samples)
    except Exception:
        return False


def _parse_featurecounts_summary(
    summary_file: str,
    sample_names: list[str],
    warnings: list[str],
) -> dict[str, Any] | None:
    """Parse featureCounts ``.summary`` table into structured assignment metrics."""
    p = Path(summary_file)
    if not p.is_file():
        return None
    try:
        rows: dict[str, list[int]] = {}
        with open(p, encoding="utf-8") as fh:
            header = fh.readline()
            if not header:
                return None
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                status = parts[0].strip()
                vals = [int(x) for x in parts[1:]]
                rows[status] = vals
        if not rows or "Assigned" not in rows:
            return None
        assigned_vec = rows["Assigned"]
        n_cols = len(assigned_vec)
        totals_vec = [
            sum(vals[i] for vals in rows.values() if i < len(vals))
            for i in range(n_cols)
        ]
        total_assigned = int(sum(assigned_vec))
        total_all = int(sum(totals_vec))
        pct_assigned = (
            round(total_assigned / total_all * 100.0, 2) if total_all > 0 else 0.0
        )
        per_sample_pct: dict[str, float] = {}
        for i in range(min(len(sample_names), n_cols)):
            denom = totals_vec[i]
            per_sample_pct[sample_names[i]] = (
                round(assigned_vec[i] / denom * 100.0, 2) if denom > 0 else 0.0
            )
        if total_all > 0 and pct_assigned < 20.0:
            warnings.append(
                f"Low featureCounts assignment rate ({pct_assigned:.1f}% assigned). "
                "Verify library strandedness and GTF annotation compatibility."
            )
        return {
            "n_assigned": total_assigned,
            "n_total": total_all,
            "pct_assigned": pct_assigned,
            "per_sample_pct_assigned": per_sample_pct,
        }
    except Exception:
        return None


def _clean_counts_matrix(
    counts_file: str,
    sample_names: list,
    output_dir: Path,
    warnings: list,
    gtf_file: str | None = None,
) -> tuple[Path, int]:
    """Clean featureCounts output into a ``genes × samples`` TSV."""
    import pandas as pd

    df = pd.read_csv(counts_file, sep="\t", comment="#")
    count_cols = df.columns[6:].tolist()

    counts = df[["Geneid"] + count_cols].copy()
    counts.columns = ["gene_id"] + sample_names
    counts = counts.set_index("gene_id")

    n_before = len(counts)
    counts = counts[counts.sum(axis=1) > 0]
    n_after = len(counts)

    if n_before - n_after > 0:
        warnings.append(
            f"{n_before - n_after} genes with zero counts across all "
            "samples removed."
        )

    out_path = output_dir / "counts_matrix.tsv"
    counts.to_csv(str(out_path), sep="\t")

    if gtf_file and Path(gtf_file).exists():
        symbol_map = _build_ensembl_to_symbol_map(gtf_file, warnings)
        if symbol_map:
            n_mapped = sum(1 for gid in counts.index if str(gid).split(".")[0] in symbol_map)
            symbols = pd.Series(
                [symbol_map.get(str(gid).split(".")[0], gid) for gid in counts.index],
                index=counts.index,
                name="gene_symbol",
            )
            counts_sym = counts.copy()
            counts_sym.insert(0, "gene_symbol", symbols)

            sym_path = output_dir / "counts_with_symbols.tsv"
            counts_sym.to_csv(str(sym_path), sep="\t")
            warnings.append(
                f"Annotated {n_mapped:,}/{len(counts):,} genes with HGNC "
                "symbols from GTF."
            )

    return out_path, len(counts)


def _build_ensembl_to_symbol_map(gtf_file: str, warnings: list) -> dict:
    """Parse a GTF to build ``{ensembl_id_without_version: gene_symbol}``."""
    mapping = {}
    opener = gzip.open if str(gtf_file).endswith(".gz") else open

    try:
        with opener(gtf_file, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] not in ("gene", "exon"):
                    continue
                if fields[2] == "exon" and mapping:
                    # Prefer gene lines when present, but allow exon-only mini GTFs
                    pass
                attrs = fields[8]
                gid_match = re.search(r'gene_id "([^"]+)"', attrs)
                sym_match = re.search(r'gene_name "([^"]+)"', attrs)
                if gid_match and sym_match:
                    gid = gid_match.group(1).split(".")[0]
                    mapping.setdefault(gid, sym_match.group(1))
    except Exception as exc:
        warnings.append(f"GTF symbol mapping failed: {exc}")
        return {}

    return mapping


def _overlaps_sorted_intervals(
    intervals: list[tuple[int, int]],
    starts: list[int],
    r_start: int,
    r_end: int,
) -> bool:
    """Return True if ``[r_start, r_end)`` overlaps any merged ``(start, end)`` interval."""
    if not intervals:
        return False
    idx = bisect.bisect_right(starts, r_end - 1) - 1
    if idx < 0:
        return False
    iv_start, iv_end = intervals[idx]
    return iv_end > r_start and iv_start < r_end


def _merge_intervals(raw: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], list[int]]:
    if not raw:
        return [], []
    ordered = sorted(raw)
    merged: list[tuple[int, int]] = [ordered[0]]
    for s, e in ordered[1:]:
        prev_s, prev_e = merged[-1]
        if s <= prev_e:
            merged[-1] = (prev_s, max(prev_e, e))
        else:
            merged.append((s, e))
    return merged, [iv[0] for iv in merged]


def _load_strand_intervals(gtf: str) -> tuple[dict[str, dict[str, Any]], str | None]:
    """Load merged ``+`` and ``-`` exon (or gene) intervals across all GTF contigs."""
    opener = gzip.open if str(gtf).endswith(".gz") else open
    exon_plus: dict[str, list[tuple[int, int]]] = {}
    exon_minus: dict[str, list[tuple[int, int]]] = {}
    gene_plus: dict[str, list[tuple[int, int]]] = {}
    gene_minus: dict[str, list[tuple[int, int]]] = {}

    MAX_EXONS = 100_000
    MAX_GENES = 20_000
    n_exons = 0
    n_genes = 0

    try:
        with opener(gtf, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feat = fields[2]
                if feat not in ("exon", "gene"):
                    continue
                strand = fields[6]
                if strand not in ("+", "-"):
                    continue
                chrom = fields[0]
                start = max(0, int(fields[3]) - 1)
                end = int(fields[4])
                if end <= start:
                    continue
                if feat == "exon" and n_exons < MAX_EXONS:
                    target = exon_plus if strand == "+" else exon_minus
                    target.setdefault(chrom, []).append((start, end))
                    n_exons += 1
                elif feat == "gene" and n_genes < MAX_GENES:
                    target = gene_plus if strand == "+" else gene_minus
                    target.setdefault(chrom, []).append((start, end))
                    n_genes += 1
                if n_exons >= MAX_EXONS and n_genes >= MAX_GENES:
                    break
    except Exception as exc:
        return {}, f"GTF parse failed for strand detection: {exc}"

    src_plus = exon_plus if (exon_plus or exon_minus) else gene_plus
    src_minus = exon_minus if (exon_plus or exon_minus) else gene_minus
    all_chroms = sorted(set(src_plus) | set(src_minus))
    if not all_chroms:
        return {}, "No stranded exon/gene intervals found in GTF."

    compiled: dict[str, dict[str, Any]] = {}
    for chrom in all_chroms:
        p_ivs, p_starts = _merge_intervals(src_plus.get(chrom, []))
        m_ivs, m_starts = _merge_intervals(src_minus.get(chrom, []))
        compiled[chrom] = {
            "plus": p_ivs,
            "plus_starts": p_starts,
            "minus": m_ivs,
            "minus_starts": m_starts,
        }
    return compiled, None


def _detect_strandedness(
    bam: str | list[str],
    gtf: str,
    warnings: list,
    *,
    min_informative_reads: int = 50,
) -> tuple[int | None, dict[str, Any] | None]:
    """Fail-closed library strandedness auto-detection across shared BAM/GTF contigs.

    Returns ``(strand_int, None)`` on unambiguous inference (`0`, `1`, or `2`), or
    ``(None, error_dict)`` when ``pysam`` is missing, inputs are unreadable,
    fewer than ``min_informative_reads`` strand-unambiguous alignments exist, or
    the sense fraction falls inside the ambiguous 20-40% / 60-80% bands.
    """
    try:
        import pysam
    except ImportError:
        return None, {
            "status": "error",
            "error_type": "MissingDependency",
            "details": (
                "pysam is required in aria-rnaseq-env when strand='auto'. "
                "Install pysam or supply an explicit strand parameter (0, 1, or 2)."
            ),
        }

    bam_list = [bam] if isinstance(bam, (str, Path)) else [str(p) for p in bam]
    existing_bams = [str(p) for p in bam_list if Path(str(p)).exists()]
    if not existing_bams:
        return None, {
            "status": "error",
            "error_type": "MissingStrandReference",
            "details": f"BAM file(s) not found for strand auto-detection: {bam_list}",
        }
    if not gtf or not Path(gtf).exists():
        return None, {
            "status": "error",
            "error_type": "MissingStrandReference",
            "details": f"GTF file not found for strand auto-detection: {gtf}",
        }

    chrom_intervals, gtf_err = _load_strand_intervals(str(gtf))
    if gtf_err or not chrom_intervals:
        return None, {
            "status": "error",
            "error_type": "EmptyGTFForStrandDetection",
            "details": gtf_err or "No stranded intervals found in GTF.",
        }

    READ_TARGET = 200_000
    sense = 0
    anti = 0
    shared_contigs_seen: set[str] = set()

    def _score_read(read, iv_map: dict[str, Any]) -> int:
        """Return +1 for sense, -1 for antisense, 0 if uninformative/ambiguous."""
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            return 0
        if read.is_paired and read.is_read2:
            return 0
        if read.mapping_quality < 30:
            return 0
        r_start = int(read.reference_start)
        r_end = int(read.reference_end) if read.reference_end else (r_start + max(1, read.query_length or 50))
        if r_end <= r_start:
            return 0
        hits_plus = _overlaps_sorted_intervals(
            iv_map["plus"], iv_map["plus_starts"], r_start, r_end
        )
        hits_minus = _overlaps_sorted_intervals(
            iv_map["minus"], iv_map["minus_starts"], r_start, r_end
        )
        if hits_plus == hits_minus:
            # Either intergenic (neither) or overlapping sense+antisense genes (both)
            return 0
        read_on_plus = not read.is_reverse
        gene_on_plus = hits_plus
        return 1 if (read_on_plus == gene_on_plus) else -1

    for bam_path in existing_bams:
        if (sense + anti) >= READ_TARGET:
            break
        try:
            bam_obj = pysam.AlignmentFile(bam_path, "rb")
        except Exception as exc:
            return None, {
                "status": "error",
                "error_type": "UnreadableBamForStrandDetection",
                "details": f"Cannot open BAM {bam_path} for strand detection: {exc}",
            }

        try:
            bam_refs = list(bam_obj.references)
            shared = [c for c in bam_refs if c in chrom_intervals]
            if not shared:
                continue
            shared_contigs_seen.update(shared)

            # Prioritize canonical large chromosomes when indexed, then remaining contigs
            priority = [c for c in ("chr1", "1", "chrI") if c in chrom_intervals and c in shared]
            ordered_contigs = priority + [c for c in shared if c not in priority]

            has_idx = False
            try:
                has_idx = bool(bam_obj.has_index())
            except Exception:
                has_idx = False

            if has_idx:
                for contig in ordered_contigs:
                    if (sense + anti) >= READ_TARGET:
                        break
                    iv_map = chrom_intervals[contig]
                    for read in bam_obj.fetch(contig):
                        verdict = _score_read(read, iv_map)
                        if verdict == 1:
                            sense += 1
                        elif verdict == -1:
                            anti += 1
                        if (sense + anti) >= READ_TARGET:
                            break
            else:
                for read in bam_obj.fetch(until_eof=True):
                    if read.is_unmapped:
                        continue
                    contig = read.reference_name
                    if not contig or contig not in chrom_intervals:
                        continue
                    verdict = _score_read(read, chrom_intervals[contig])
                    if verdict == 1:
                        sense += 1
                    elif verdict == -1:
                        anti += 1
                    if (sense + anti) >= READ_TARGET:
                        break
        finally:
            bam_obj.close()

    if not shared_contigs_seen:
        return None, {
            "status": "error",
            "error_type": "ContigMismatchForStrandDetection",
            "details": (
                "No chromosome/contig names match between aligned BAM(s) and GTF "
                "(e.g. 'chr1' vs '1' mismatch). Verify FASTA and GTF build match."
            ),
        }

    total = sense + anti
    if total < min_informative_reads:
        return None, {
            "status": "error",
            "error_type": "InsufficientReadsForStrandDetection",
            "details": (
                f"Strand auto-detection found only {total} informative primary "
                f"alignments across {len(shared_contigs_seen)} shared contig(s) "
                f"(minimum {min_informative_reads} required). Supply an explicit "
                "strand parameter (0=unstranded, 1=forward, 2=reverse)."
            ),
            "informative_reads": total,
        }

    sense_frac = sense / total
    if sense_frac >= 0.80:
        warnings.append(
            f"Strandedness auto-detected: forward stranded "
            f"({sense_frac * 100:.1f}% sense reads, n={total}). "
            "Using featureCounts -s 1."
        )
        return 1, None
    if sense_frac <= 0.20:
        warnings.append(
            f"Strandedness auto-detected: reverse stranded "
            f"({sense_frac * 100:.1f}% sense reads, n={total}). "
            "Using featureCounts -s 2."
        )
        return 2, None
    if 0.40 <= sense_frac <= 0.60:
        warnings.append(
            f"Strandedness auto-detected: unstranded "
            f"({sense_frac * 100:.1f}% sense reads, n={total}). "
            "Using featureCounts -s 0."
        )
        return 0, None

    return None, {
        "status": "error",
        "error_type": "AmbiguousStrandedness",
        "details": (
            f"Library strandedness is ambiguous ({sense_frac * 100:.1f}% sense "
            f"reads across n={total} informative alignments; expected >=80% for "
            "forward -s 1, <=20% for reverse -s 2, or 40-60% for unstranded -s 0). "
            "Inspect genomic DNA contamination or specify explicit strand (0, 1, or 2)."
        ),
        "sense_fraction": round(sense_frac, 4),
        "informative_reads": total,
    }


def _read_layout_for_bams(bam_files: list, paired_param=None) -> str:
    """Resolve one homogeneous read layout, preserving the legacy default."""
    layouts = {
        b.get(
            "read_layout",
            "paired-end" if b.get("paired", paired_param is not False) else "single-end",
        )
        for b in bam_files
    }
    return next(iter(layouts)) if len(layouts) == 1 else "mixed"


def _mock_counts_matrix(
    bam_files: list,
    output_dir: Path,
    warnings: list,
    read_layout: str = "paired-end",
) -> dict:
    """Generate mock counts matrix when featureCounts is unavailable in mock mode."""
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(42)
    n_genes = 20000
    samples = [b["name"] for b in bam_files]
    genes = [f"ENSMUSG{i:011d}" for i in range(n_genes)]

    base = rng.negative_binomial(20, 0.3, (n_genes, len(samples)))
    de_idx = list(range(500))

    n_grp = len(samples) // 2
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
        "status": "success",
        "counts_matrix": str(out_path),
        "n_genes": n_genes,
        "n_samples": len(samples),
        "strand_used": 0,
        "read_layout": read_layout,
        "sample_names": samples,
        "warnings": warnings,
        "note": "mock",
    }


if __name__ == "__main__":
    run_script(rna_quantify)
