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


def test_r2_lobo_index_error_does_not_abort_chromatin_bulk_diffacc(tmp_path, monkeypatch):
    """R2: Unexpected LOBO IndexError must record lobo_jackknife error without aborting bulk chromatin DA."""
    pytest.importorskip("pydeseq2")
    import numpy as np
    import pandas as pd
    import aria.scripts.chromatin_bulk_diffacc as cbd

    rng = np.random.default_rng(14)
    peaks = [f"chr1:{i * 200}-{i * 200 + 100}" for i in range(40)]
    units = ["ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"]
    mat = rng.poisson(60, size=(40, 6))
    mat[:8, 3:] *= 5
    counts_p = tmp_path / "peaks.tsv"
    meta_p = tmp_path / "meta.tsv"
    pd.DataFrame(mat, index=peaks, columns=units).to_csv(counts_p, sep="\t", index_label="peak_id")
    pd.DataFrame(
        {"condition": ["ctrl"] * 3 + ["treat"] * 3, "replicate": [f"r{i}" for i in range(6)]},
        index=units,
    ).to_csv(meta_p, sep="\t", index_label="sample_id")
    monkeypatch.setattr(
        cbd, "compute_counts_lobo_jackknife",
        lambda *a, **kw: (_ for _ in ()).throw(IndexError("lobo fold index out of bounds")),
    )
    res = cbd.chromatin_bulk_diffacc({
        "counts_matrix_path": str(counts_p),
        "sample_metadata_path": str(meta_p),
        "comparisons": [["treat", "ctrl"]],
        "lfc_min": 0.0,
        "n_cpus": 1,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success"
    assert res["comparisons"][0]["status"] == "success"
    assert res["comparisons"][0]["lobo_jackknife"]["status"] == "error"
    assert "index out of bounds" in res["comparisons"][0]["lobo_jackknife"]["details"]


def test_r2_lobo_index_error_does_not_abort_chromatin_diffacc(tmp_path, monkeypatch):
    """R2: Unexpected LOBO IndexError must record lobo_jackknife error without aborting scATAC pseudobulk DA."""
    pytest.importorskip("pydeseq2")
    ad = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd
    from scipy import sparse
    import aria.scripts.chromatin_diffacc as cd

    rng = np.random.default_rng(15)
    peaks = [f"chr1:{i * 200}-{i * 200 + 100}" for i in range(40)]
    donors = ["d1", "d2", "d3", "d4", "d5", "d6"]
    conds = ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]
    X_rows, obs_rows = [], []
    for d, c in zip(donors, conds):
        base = rng.poisson(30, size=(4, 40))
        if c == "treat":
            base[:, :6] = rng.poisson(140, size=(4, 6))
        X_rows.append(base)
        obs_rows.extend({"leiden": "0", "condition": c, "donor": d} for _ in range(4))
    h5 = tmp_path / "scatac.h5ad"
    ad.AnnData(
        X=sparse.csr_matrix(np.vstack(X_rows).astype(np.int32)),
        obs=pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(24)]),
        var=pd.DataFrame(index=peaks),
    ).write_h5ad(h5)
    monkeypatch.setattr(
        cd, "compute_counts_lobo_jackknife",
        lambda *a, **kw: (_ for _ in ()).throw(IndexError("lobo fold index out of bounds")),
    )
    res = cd.chromatin_diffacc({
        "data_path": str(h5),
        "groupby": "leiden",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["treat", "ctrl"]],
        "skip_per_cluster": True,
        "min_cells_per_pseudosample": 2,
        "min_replicates_per_condition": 2,
        "lfc_min": 0.0,
        "n_cpus": 1,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success"
    comp = res["pseudobulk"]["comparisons"][0]
    assert comp["status"] == "success"
    assert comp["lobo_jackknife"]["status"] == "error"
    assert "index out of bounds" in comp["lobo_jackknife"]["details"]
