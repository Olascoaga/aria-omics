"""Permanent regression tests for D1 silent-failure fixes in chromatin compute kernels."""
from __future__ import annotations

import gzip

import pytest


def test_d1_008_fragment_sizes_corrupt_gzip_must_fail_closed(tmp_path):
    """D1-008: Corrupt fragments.tsv.gz must raise instead of returning 0 barcodes / 0 fragments."""
    from aria.scripts.chromatin_qc import _compute_fragment_sizes

    bad_gz = tmp_path / "fragments.tsv.gz"
    bad_gz.write_bytes(b"NOT_A_VALID_GZIP_STREAM")
    with pytest.raises((OSError, gzip.BadGzipFile)):
        _compute_fragment_sizes(str(bad_gz))


def test_d1_015_frip_stats_corrupt_fragments_must_fail_closed(tmp_path):
    """D1-015: Corrupt fragments file with supplied peaks BED must fail closed in _estimate_frip_stats."""
    from aria.scripts.chromatin_qc import _estimate_frip_stats

    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"NOT_GZIP")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t100\t300\tp1\n")
    with pytest.raises((OSError, gzip.BadGzipFile)):
        _estimate_frip_stats(str(frag), peaks_bed=str(peaks))


def test_d1_015_frip_distribution_corrupt_fragments_must_fail_closed(tmp_path):
    """D1-015 sibling: Corrupt fragments file must fail closed in _estimate_frip_distribution."""
    from aria.scripts.chromatin_qc import _estimate_frip_distribution

    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"NOT_GZIP")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t100\t300\tp1\n")
    with pytest.raises((OSError, gzip.BadGzipFile)):
        _estimate_frip_distribution(str(frag), peaks_bed=str(peaks))


def test_d1_021_parse_gtf_gene_models_corrupt_must_fail_closed(tmp_path):
    """D1-021: Corrupt GTF must raise instead of silently returning {}."""
    from aria.scripts.chromatin_peak_annotation import parse_gtf_gene_models

    bad = tmp_path / "corrupt.gtf.gz"
    bad.write_bytes(b"NOT_GZIP_BYTES")
    with pytest.raises((OSError, gzip.BadGzipFile)):
        parse_gtf_gene_models(str(bad))


def test_d1_022_count_peaks_unreadable_must_raise(tmp_path):
    """D1-022: Unreadable/missing peak file must raise OSError instead of returning 0."""
    from aria.scripts.chromatin_peaks import _count_peaks

    with pytest.raises(OSError):
        _count_peaks(str(tmp_path / "missing.narrowPeak"))


def test_d1_018_bulk_chromatin_qc_must_fail_closed_on_corrupt_bam(tmp_path):
    """D1-018: Corrupt BAM in a bulk chromatin cohort must fail closed instead of averaging a subset."""
    pysam = pytest.importorskip("pysam")
    from aria.scripts.chromatin_qc import _bulk_chromatin_qc

    good_bam = tmp_path / "rep1.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 100000, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(good_bam), "wb", header=header) as outf:
        for i in range(10):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag = f"r{i}", "A" * 50, 0
            a.reference_id, a.reference_start, a.mapping_quality = 0, 1000 + i * 100, 60
            a.cigar, a.query_qualities = [(0, 50)], pysam.qualitystring_to_array("I" * 50)
            outf.write(a)
    bad_bam = tmp_path / "rep2_corrupt.bam"
    bad_bam.write_bytes(b"NOT_A_BAM_FILE")
    res = _bulk_chromatin_qc(
        [str(good_bam), str(bad_bam)], "bulk_ATAC", "hg38", "open_chromatin", [], params={}
    )
    assert res["status"] == "error"


def test_d1_016_h5ad_atac_qc_frag_parse_error_must_not_skip_mito(tmp_path):
    """D1-016: Non-numeric fragment obs column must not abort mito obs coercion or omit metrics_not_computed."""
    ad = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd
    from aria.scripts.chromatin_qc import _h5ad_atac_qc

    obs = pd.DataFrame(
        {"n_fragments": ["1,000", "2,000", "3,000"], "pct_mito": [0.04, 0.06, 0.08]},
        index=["bc1", "bc2", "bc3"],
    )
    p = tmp_path / "atac.h5ad"
    ad.AnnData(X=np.ones((3, 4), dtype=np.float32), obs=obs).write_h5ad(p)
    res = _h5ad_atac_qc(str(p), "hg38", "Homo sapiens", [])
    assert res["median_mito_obs"] == pytest.approx(0.06)
    assert any("median_fragments_per_cell" in m for m in res["metrics_not_computed"])
