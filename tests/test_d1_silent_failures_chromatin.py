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
