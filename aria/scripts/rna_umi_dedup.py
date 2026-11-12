"""
ARIA Bulk RNA-seq Inline-UMI Deduplication Script (ADR-075)
-----------------------------------------------------------
Deduplicates coordinate-sorted STAR BAMs whose read names carry fastp-extracted
inline UMIs (``<qname>:<UMI>``), using the directional UMI adjacency criterion
(Smith et al., Genome Research 2017 — ``counts[u] >= 2 * counts[v] - 1`` within
Hamming distance 1) per 5'-anchor genomic locus and strand.

Executed inside ``aria-rnaseq-env`` via ``EnvironmentManager``.
Zero mock fabrication in production; content-addressed A4 resume via
``stage_manifest.py``.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.utils.stage_manifest import (
    build_stage_manifest,
    stage_is_current,
    write_stage_manifest,
)

_UMI_RE = re.compile(r"^[ACGTN]{4,32}$")


def _pysam_version() -> str:
    try:
        import pysam

        return f"pysam={pysam.__version__}"
    except Exception:
        return ""


def _extract_umi(query_name: str, separator: str = ":") -> str | None:
    if not query_name or separator not in query_name:
        return None
    token = query_name.rsplit(separator, 1)[-1].strip().upper()
    if "/" in token:
        token = token.split("/", 1)[0]
    if _UMI_RE.fullmatch(token):
        return token
    return None


def _hamming_le_1(a: str, b: str) -> bool:
    if len(a) != len(b):
        return False
    mismatches = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            mismatches += 1
            if mismatches > 1:
                return False
    return True


def _directional_winning_umis(umi_counts: Counter[str]) -> set[str]:
    """Return canonical molecule UMIs under the UMI-tools directional graph rule."""
    if len(umi_counts) <= 1:
        return set(umi_counts.keys())

    # Sort descending by count, then lexicographically for strict determinism.
    ordered = [
        umi for umi, _ in sorted(umi_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    ]
    adj: dict[str, list[str]] = defaultdict(list)
    n = len(ordered)
    for i in range(n):
        u = ordered[i]
        cu = umi_counts[u]
        for j in range(i + 1, n):
            v = ordered[j]
            cv = umi_counts[v]
            if cu >= 2 * cv - 1 and _hamming_le_1(u, v):
                adj[u].append(v)

    visited: set[str] = set()
    winners: set[str] = set()
    for root in ordered:
        if root in visited:
            continue
        winners.add(root)
        stack = [root]
        visited.add(root)
        while stack:
            curr = stack.pop()
            for nxt in adj.get(curr, ()):
                if nxt not in visited:
                    visited.add(nxt)
                    stack.append(nxt)
    return winners


def _five_prime_coord(read) -> int:
    if read.is_reverse:
        end = read.reference_end
        return (end - 1) if end is not None else int(read.reference_start)
    return int(read.reference_start)


def _dedup_single_bam(
    sample_rec: dict[str, Any],
    output_dir: Path,
    separator: str,
    warnings: list[str],
    tool_version: str,
) -> dict[str, Any]:
    import pysam

    name = str(sample_rec.get("name", "sample"))
    in_bam = Path(str(sample_rec["bam"]))
    read_layout = str(
        sample_rec.get(
            "read_layout",
            "paired-end" if sample_rec.get("paired") else "single-end",
        )
    )
    paired = read_layout == "paired-end"

    out_bam = output_dir / f"{name}.umi_dedup.sorted.bam"
    manifest_path = output_dir / f"{name}.umi_dedup.stage.json"
    metrics_path = output_dir / f"{name}.umi_dedup.metrics.json"
    stage_inputs = [("input_bam", str(in_bam))]
    stage_params = {
        "method": "directional_1mm",
        "separator": separator,
        "read_layout": read_layout,
    }

    current, _ = stage_is_current(
        str(manifest_path),
        inputs=stage_inputs,
        params=stage_params,
        tool_version=tool_version,
    )
    if current and out_bam.exists() and out_bam.stat().st_size > 128 and metrics_path.exists():
        import json

        saved = json.loads(metrics_path.read_text(encoding="utf-8"))
        warnings.append(
            f"[resume] UMI deduplication skipped for {name} "
            f"({saved.get('n_molecules_out', 0)} molecules retained)"
        )
        resumed_rec = dict(sample_rec)
        resumed_rec.update(saved)
        resumed_rec["bam"] = str(out_bam)
        resumed_rec["resumed"] = True
        return resumed_rec

    # Pass 1: scan primary alignments (R1 for paired-end, all primary for SE) to
    # determine winning query_names per (chrom, strand, 5'-anchor[, tlen]) + UMI.
    locus_umi_counts: dict[tuple, Counter[str]] = defaultdict(Counter)
    # Best representative per (locus, umi): (score_tuple, query_name)
    best_qname_for_locus_umi: dict[tuple[tuple, str], tuple[tuple[int, int, str], str]] = {}

    n_primary_checked = 0
    n_with_umi = 0

    with pysam.AlignmentFile(str(in_bam), "rb") as reader:
        for read in reader.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if paired and read.is_paired and read.is_read2:
                continue
            n_primary_checked += 1
            qname = read.query_name or ""
            umi = _extract_umi(qname, separator=separator)
            if umi is None:
                continue
            n_with_umi += 1

            anchor = _five_prime_coord(read)
            if paired and read.is_paired:
                locus = (
                    int(read.reference_id),
                    bool(read.is_reverse),
                    anchor,
                    int(read.template_length),
                )
            else:
                locus = (
                    int(read.reference_id),
                    bool(read.is_reverse),
                    anchor,
                )

            locus_umi_counts[locus][umi] += 1
            mean_q = int(sum(read.query_qualities or (0,)))
            # Negate qname tie-breaker? Keep deterministic: higher MAPQ, higher sum(baseQ),
            # lexicographically smallest qname wins.
            cand_key = (-int(read.mapping_quality), -mean_q, qname)
            existing = best_qname_for_locus_umi.get((locus, umi))
            if existing is None or cand_key < existing[0]:
                best_qname_for_locus_umi[(locus, umi)] = (cand_key, qname)

    if n_primary_checked == 0:
        return {
            "name": name,
            "status": "failed",
            "error_type": "EmptyAlignment",
            "details": f"No primary mapped reads found in {in_bam} for UMI deduplication.",
        }

    umi_frac = n_with_umi / max(n_primary_checked, 1)
    if umi_frac < 0.90:
        return {
            "name": name,
            "status": "failed",
            "error_type": "MissingReadHeaderUMIs",
            "details": (
                f"Sample '{name}' has valid header UMIs on only "
                f"{umi_frac * 100:.1f}% of mapped reads ({n_with_umi}/{n_primary_checked}). "
                "Verify fastp UMI extraction parameters (umi_loc / umi_len)."
            ),
        }

    keep_qnames: set[str] = set()
    for locus, counts in locus_umi_counts.items():
        winning_umis = _directional_winning_umis(counts)
        for umi in winning_umis:
            keep_qnames.add(best_qname_for_locus_umi[(locus, umi)][1])

    # Pass 2: stream coordinate-sorted input BAM to coordinate-sorted output BAM,
    # keeping primary alignments whose query_name won directional UMI collapse.
    unsorted_or_sorted = str(out_bam)
    with pysam.AlignmentFile(str(in_bam), "rb") as reader, pysam.AlignmentFile(
        unsorted_or_sorted, "wb", template=reader
    ) as writer:
        for read in reader.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if (read.query_name or "") in keep_qnames:
                writer.write(read)

    subprocess.run(
        ["samtools", "index", str(out_bam)],
        capture_output=True,
        timeout=300,
        check=False,
    )

    n_molecules = len(keep_qnames)
    dedup_rate = round((1.0 - (n_molecules / max(n_with_umi, 1))) * 100.0, 2)
    metrics = {
        "status": "success",
        "n_reads_pre_dedup": n_with_umi,
        "n_molecules_out": n_molecules,
        "umi_dedup_rate_pct": dedup_rate,
        "umi_deduplicated": True,
    }

    import json

    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    write_stage_manifest(
        str(manifest_path),
        build_stage_manifest(
            stage="umi_dedup",
            inputs=stage_inputs,
            params=stage_params,
            tool_version=tool_version,
        ),
    )

    updated = dict(sample_rec)
    updated.update(metrics)
    updated["bam"] = str(out_bam)
    return updated


def rna_umi_dedup(params: dict[str, Any]) -> dict[str, Any]:
    try:
        import pysam  # noqa: F401
    except ImportError:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": (
                "pysam is required in aria-rnaseq-env for inline-UMI "
                "deduplication (protocol='umi_bulk')."
            ),
        }

    bam_files = params.get("bam_files") or []
    if not bam_files:
        return {
            "status": "error",
            "error_type": "NoBamsFound",
            "details": "No BAM records were provided to rna_umi_dedup.",
        }

    first_bam = Path(str(bam_files[0].get("bam", "/tmp/sample.bam")))
    output_dir = Path(
        params.get("output_dir") or str(first_bam.parent.parent / "umi_dedup")
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    separator = str(params.get("umi_separator") or ":")
    warnings: list[str] = []
    tool_version = _pysam_version()

    dedup_bams: list[dict[str, Any]] = []
    for rec in bam_files:
        if rec.get("status") != "success" or not rec.get("bam") or not Path(rec["bam"]).exists():
            return {
                "status": "error",
                "error_type": "InvalidAlignmentInput",
                "details": f"Cannot UMI-deduplicate invalid BAM record: {rec.get('name')}",
            }
        out_rec = _dedup_single_bam(
            rec,
            output_dir=output_dir,
            separator=separator,
            warnings=warnings,
            tool_version=tool_version,
        )
        dedup_bams.append(out_rec)

    failed = [b.get("name", "unknown") for b in dedup_bams if b.get("status") != "success"]
    if failed:
        first_err = next(b for b in dedup_bams if b.get("status") != "success")
        return {
            "status": "error",
            "error_type": first_err.get("error_type", "UmiDedupFailed"),
            "details": first_err.get("details", "UMI deduplication failed."),
            "failed_samples": failed,
            "bam_files": dedup_bams,
            "warnings": warnings,
        }

    mean_dedup = round(
        sum(float(b.get("umi_dedup_rate_pct", 0.0)) for b in dedup_bams)
        / max(len(dedup_bams), 1),
        2,
    )
    return {
        "status": "success",
        "bam_files": dedup_bams,
        "n_deduplicated": len(dedup_bams),
        "mean_umi_dedup_rate_pct": mean_dedup,
        "output_dir": str(output_dir),
        "warnings": warnings,
    }


if __name__ == "__main__":
    run_script(rna_umi_dedup)
