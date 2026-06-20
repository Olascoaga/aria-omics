"""F10 (preprint audit 2026-06-19): reproducible peak-to-gene marker mapping for
bulk ATAC DA sanity checks.

The V47 DA validation artifact reported a per-marker `mean_log2fc` (GATA1/HBB/... →
K562, MS4A1/PAX5/... → GM12878) as a directional sanity check. That mapping was
hand-curated (e.g. CD19's "narrow 2-peak window"), so the numbers were not
reproducible from declared parameters. This module makes the *method* reproducible and
deterministic: given DA rows (peak + log2FoldChange), marker gene TSS coordinates, and
an explicit window, it selects the peaks within the window of each gene's TSS and
returns the mean log2FoldChange.

It is a directional sanity helper, NOT a calibration/sensitivity claim: FDR and recall
calibration are validated separately by the synthetic recovery gate
(`aria/benchmarks/synthetic_atac_da.py`, `tests/test_benchmark_scatac_da.py`).
"""
from __future__ import annotations

from typing import Iterable, Optional


def parse_peak_interval(peak: str) -> Optional[tuple[str, int, int]]:
    """Parse a ``chrom:start-end`` peak id into ``(chrom, start, end)``.

    Returns ``None`` for anything that does not match (never fabricates coordinates).
    """
    if not isinstance(peak, str) or ":" not in peak or "-" not in peak:
        return None
    chrom, _, span = peak.partition(":")
    start_s, _, end_s = span.partition("-")
    try:
        start, end = int(start_s), int(end_s)
    except ValueError:
        return None
    if not chrom or end < start:
        return None
    return chrom, start, end


def map_peaks_to_gene_markers(
    da_rows: Iterable[dict],
    genes: Iterable[dict],
    window_bp: int,
    *,
    peak_key: str = "peak",
    lfc_key: str = "log2FoldChange",
) -> list:
    """Map DA peaks to marker genes by a declared TSS window (pure, deterministic).

    ``da_rows``: rows with a peak id (``peak_key``, ``chrom:start-end``) and a
    log2 fold change (``lfc_key``). ``genes``: dicts with ``gene``, ``chrom``, ``tss``
    (int) and optional pass-through fields (e.g. ``expected``, ``lineage``).
    A peak matches a gene when it is on the same chromosome and its midpoint is within
    ``window_bp`` of the gene's TSS. Returns one result dict per gene with the selected
    peaks (sorted by genomic distance, then peak id), the count, and the mean
    log2FoldChange (``None`` when no peak falls in the window — never imputed).
    """
    if window_bp < 0:
        raise ValueError("window_bp must be >= 0")

    parsed: list[tuple[str, int, float]] = []  # (chrom, midpoint, lfc)
    for row in da_rows:
        iv = parse_peak_interval(str(row.get(peak_key, "")))
        if iv is None:
            continue
        try:
            lfc = float(row.get(lfc_key))
        except (TypeError, ValueError):
            continue
        chrom, start, end = iv
        parsed.append((chrom, (start + end) // 2, lfc))

    out: list = []
    for g in genes:
        chrom = str(g.get("chrom", ""))
        try:
            tss = int(g.get("tss"))
        except (TypeError, ValueError):
            tss = None
        selected = []
        if chrom and tss is not None:
            for pc, mid, lfc in parsed:
                if pc == chrom and abs(mid - tss) <= window_bp:
                    selected.append((abs(mid - tss), mid, lfc))
        selected.sort(key=lambda x: (x[0], x[1]))
        lfcs = [lfc for _d, _m, lfc in selected]
        mean_lfc = round(sum(lfcs) / len(lfcs), 4) if lfcs else None
        entry = {k: v for k, v in g.items()}
        entry.update({
            "n_peaks_in_window": len(lfcs),
            "mean_log2fc": mean_lfc,
            "window_bp": int(window_bp),
        })
        if mean_lfc is not None:
            entry["called_direction"] = "positive" if mean_lfc > 0 else (
                "negative" if mean_lfc < 0 else "zero")
        out.append(entry)
    return out
