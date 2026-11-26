"""Regression tests for D4 genomic coordinate and contig normalization fixes."""
from __future__ import annotations

from pathlib import Path

import pytest


def test_d4_014_peak_annotation_converts_gtf_1based_to_bed_0based_once(tmp_path: Path):
    """D4-014: parse_gtf_gene_models must convert 1-based closed GTF [start, end]
    to 0-based half-open BED [start-1, end) once at parse time so BED peak boundaries
    overlap exons accurately and TSS coordinates are 0-based."""
    from aria.scripts.chromatin_peak_annotation import (
        DEFAULT_PROMOTER_DOWNSTREAM,
        DEFAULT_PROMOTER_UPSTREAM,
        annotate_one_peak,
        parse_gtf_gene_models,
    )

    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tHAVANA\tgene\t1000\t10000\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n'
        'chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n'
        'chr1\tHAVANA\texon\t5000\t5200\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n',
        encoding="utf-8",
    )
    models = parse_gtf_gene_models(gtf)
    tss_index = {c: [g["tss"] for g in glist] for c, glist in models.items()}

    # BED [4800, 5000) covers 0-based base 4999 == 1-based base 5000 (first base of exon 5000..5200).
    ann_overlap = annotate_one_peak(
        "chr1", 4800, 5000, models, tss_index,
        DEFAULT_PROMOTER_UPSTREAM, DEFAULT_PROMOTER_DOWNSTREAM,
    )
    assert ann_overlap["feature"] == "Exonic", f"Expected Exonic at 1-based base 5000, got {ann_overlap}"

    # BED [4800, 4999) ends immediately BEFORE 0-based base 4999 -> Intronic.
    ann_adjacent = annotate_one_peak(
        "chr1", 4800, 4999, models, tss_index,
        DEFAULT_PROMOTER_UPSTREAM, DEFAULT_PROMOTER_DOWNSTREAM,
    )
    assert ann_adjacent["feature"] == "Intronic", f"Expected Intronic for [4800, 4999), got {ann_adjacent}"
