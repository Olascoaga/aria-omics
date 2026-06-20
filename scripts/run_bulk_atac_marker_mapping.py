#!/usr/bin/env python3
"""Reproducible bulk ATAC peak-to-gene marker mapping (F10).

Thin CLI over ``aria.benchmarks.bulk_atac_marker_mapping.map_peaks_to_gene_markers``.
Given a DA results CSV (peak + log2FoldChange columns) and a marker-gene TSS TSV
(gene, chrom, tss[, expected, lineage]) plus an explicit ``--window-bp``, it emits a
deterministic per-marker mean log2FoldChange — the directional sanity check the V47
artifact reported, now regenerable from declared parameters rather than hand-curated.

This is a directional sanity helper, NOT a calibration/sensitivity claim.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.bulk_atac_marker_mapping import map_peaks_to_gene_markers  # noqa: E402


def _read_rows(path: str) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--da-csv", required=True,
                        help="DA results CSV with peak + log2FoldChange columns.")
    parser.add_argument("--genes-tsv", required=True,
                        help="TSV with gene, chrom, tss columns (optional extras).")
    parser.add_argument("--window-bp", type=int, default=50000,
                        help="TSS +/- window for peak selection (declared param).")
    parser.add_argument("--peak-key", default="peak")
    parser.add_argument("--lfc-key", default="log2FoldChange")
    parser.add_argument("--output", default=None, help="Optional JSON output path.")
    args = parser.parse_args(argv)

    da_rows = _read_rows(args.da_csv)
    genes = []
    with open(args.genes_tsv, newline="", encoding="utf-8") as fh:
        for g in csv.DictReader(fh, delimiter="\t"):
            genes.append(g)

    mapping = map_peaks_to_gene_markers(
        da_rows, genes, args.window_bp,
        peak_key=args.peak_key, lfc_key=args.lfc_key)
    manifest = {
        "analysis": "bulk_atac_marker_peak_to_gene_mapping",
        "window_bp": args.window_bp,
        "n_da_rows": len(da_rows),
        "n_genes": len(genes),
        "markers": mapping,
    }
    payload = json.dumps(manifest, indent=2)
    if args.output:
        Path(args.output).write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
