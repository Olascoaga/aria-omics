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


def test_d4_016_count_cram_normalizes_chr_vs_ensembl_contigs(tmp_path: Path):
    """D4-016: chromatin_peak_counts._count_cram must normalize UCSC vs Ensembl contigs
    via canonical genomic_coords helper instead of returning silent zero counts."""
    pysam = pytest.importorskip("pysam")
    from aria.scripts.chromatin_peak_counts import _count_cram

    bam = tmp_path / "ensembl.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": 50_000}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i, pos in enumerate((1000, 1020, 1050)):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag = f"r{i}", "A" * 100, 0
            a.reference_id, a.reference_start, a.mapping_quality = 0, pos, 60
            a.cigar = [(0, 100)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(a)
    pysam.index(str(bam))

    bed = tmp_path / "peaks_ucsc.bed"
    bed.write_text("chr1\t950\t1200\tpeak_1\n", encoding="utf-8")
    assert _count_cram(bed, bam) == [3]


def test_d4_017_peak2gene_gc_uses_0based_bed_start_without_1bp_shift(tmp_path: Path):
    """D4-017: compute_peak_gc_from_indexed_fasta must treat peak coords as 0-based
    half-open BED [start, end) without subtracting 1 from start."""
    pysam = pytest.importorskip("pysam")
    from aria.utils.peak2gene_core import compute_peak_gc_from_indexed_fasta

    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + ("G" * 10) + ("A" * 10) + ("T" * 20) + "\n", encoding="utf-8")
    pysam.faidx(str(fasta))

    gc, info = compute_peak_gc_from_indexed_fasta(str(fasta), [("chr1", 10, 20)])
    assert info["available"] is True
    assert gc[0] == pytest.approx(0.0, abs=1e-9), f"Expected GC=0.0 for [10,20), got {gc[0]:.6f}"

    # Verify [0, 10) is pure G (GC=1.0) and [9, 11) is 'GA' (GC=0.5).
    gc_edge, _ = compute_peak_gc_from_indexed_fasta(str(fasta), [("chr1", 0, 10), ("chr1", 9, 11)])
    assert gc_edge[0] == pytest.approx(1.0, abs=1e-9)
    assert gc_edge[1] == pytest.approx(0.5, abs=1e-9)
