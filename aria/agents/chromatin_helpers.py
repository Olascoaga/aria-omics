"""Pure helper functions shared by chromatin agent lanes."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional


def bulk_da_motif_regions(comparisons, *, max_per_group: int = 5000):
    """Build direction-split motif regions from bulk ATAC DA comparisons.

    Returns ``(regions, background, warnings)``. Empty ``regions`` is an honest
    "nothing to interpret" signal for callers.
    """
    regions: dict[str, list[str]] = {}
    background: list[str] = []
    seen_bg: set[str] = set()
    warnings: list[str] = []

    for comp in comparisons or []:
        if not isinstance(comp, dict) or comp.get("status") != "success":
            continue
        csv_path = comp.get("full_results_csv")
        if not csv_path or not Path(str(csv_path)).is_file():
            warnings.append(
                f"comparison '{comp.get('test')}_vs_{comp.get('reference')}': "
                f"DA results CSV not found for motif enrichment.")
            continue
        test = str(comp.get("test", "test"))
        ref = str(comp.get("reference", "reference"))
        comp_key = f"{test}_vs_{ref}"
        up_test: list[tuple[str, float, float]] = []
        up_ref: list[tuple[str, float, float]] = []
        try:
            with open(str(csv_path), newline="", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                cols = reader.fieldnames or []
                if "peak" not in cols or "log2FoldChange" not in cols:
                    warnings.append(
                        f"comparison '{comp_key}': DA CSV missing peak/"
                        f"log2FoldChange columns; skipped for motifs.")
                    continue
                for row in reader:
                    peak = str(row.get("peak", "")).strip()
                    if not peak:
                        continue
                    if peak not in seen_bg:
                        seen_bg.add(peak)
                        background.append(peak)
                    if str(row.get("significant", "")).strip().lower() != "true":
                        continue
                    try:
                        lfc = float(row.get("log2FoldChange"))
                    except (TypeError, ValueError):
                        continue
                    try:
                        padj = float(row.get("padj"))
                    except (TypeError, ValueError):
                        padj = float("inf")
                    ranked = (peak, padj, abs(lfc))
                    if lfc > 0:
                        up_test.append(ranked)
                    elif lfc < 0:
                        up_ref.append(ranked)
        except OSError as exc:
            warnings.append(
                f"comparison '{comp_key}': could not read DA CSV ({exc}).")
            continue
        if up_test:
            regions[f"{comp_key}::up_in_{test}"] = rank_bulk_da_motif_peaks(
                up_test, max_per_group, f"{comp_key}::up_in_{test}", warnings)
        if up_ref:
            regions[f"{comp_key}::up_in_{ref}"] = rank_bulk_da_motif_peaks(
                up_ref, max_per_group, f"{comp_key}::up_in_{ref}", warnings)

    return regions, background, warnings


def rank_bulk_da_motif_peaks(rows, max_per_group: int, group: str,
                             warnings: list[str]) -> list[str]:
    ranked = sorted(rows, key=lambda item: (item[1], -item[2], item[0]))
    if len(ranked) > max_per_group:
        warnings.append(
            f"group '{group}': capped {len(ranked)} -> {max_per_group} peaks "
            "for motif scanning after ranking by padj then |log2FoldChange|.")
    return [peak for peak, _padj, _abs_lfc in ranked[:max_per_group]]


def positive_int(value, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def is_fastq(files: list) -> bool:
    """True when the input list is raw FASTQ."""
    if not files:
        return False
    return str(files[0]).lower().endswith(FASTQ_SUFFIXES)


def pick_read(fastqs: list, tokens: tuple) -> Optional[str]:
    """Pick the FASTQ whose name contains one of the read tokens."""
    for path in fastqs:
        low = str(path).lower()
        if any(token.lower() in low for token in tokens):
            return path
    return None
