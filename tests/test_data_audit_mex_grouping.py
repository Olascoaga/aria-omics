"""E2E-1 (production verification 2026-06-12): a raw 10x MEX directory must be
treated as ONE scRNA sample, not split into its component files. Otherwise the
per-sample QC path tries to load barcodes.tsv as a sample and the whole pipeline
aborts with PerSampleQCFailed."""

from aria.agents.data_audit_agent import _collapse_mex_directories


def test_single_mex_triplet_collapses_to_directory():
    paths = [
        "/data/pbmc/matrix.mtx",
        "/data/pbmc/barcodes.tsv",
        "/data/pbmc/genes.tsv",
    ]
    out = _collapse_mex_directories(paths)
    assert out == ["/data/pbmc"]


def test_gzipped_v3_features_triplet_collapses():
    paths = [
        "/data/s1/matrix.mtx.gz",
        "/data/s1/barcodes.tsv.gz",
        "/data/s1/features.tsv.gz",
    ]
    out = _collapse_mex_directories(paths)
    assert out == ["/data/s1"]


def test_two_mex_subdirs_collapse_to_two_directories():
    paths = [
        "/data/s1/matrix.mtx", "/data/s1/barcodes.tsv", "/data/s1/genes.tsv",
        "/data/s2/matrix.mtx", "/data/s2/barcodes.tsv", "/data/s2/genes.tsv",
    ]
    out = _collapse_mex_directories(paths)
    assert out == ["/data/s1", "/data/s2"]


def test_h5ad_and_h5_are_untouched():
    paths = ["/data/a.h5ad", "/data/b.h5"]
    assert _collapse_mex_directories(paths) == paths


def test_mex_dir_plus_h5ad_keeps_both():
    paths = [
        "/data/mex/matrix.mtx", "/data/mex/barcodes.tsv", "/data/mex/genes.tsv",
        "/data/other.h5ad",
    ]
    out = _collapse_mex_directories(paths)
    assert "/data/mex" in out
    assert "/data/other.h5ad" in out
    assert not any(p.endswith("barcodes.tsv") for p in out)


def test_loose_tsv_without_matrix_is_not_collapsed():
    # No matrix.mtx in the dir -> not a MEX sample, leave files as-is.
    paths = ["/data/loose/barcodes.tsv"]
    assert _collapse_mex_directories(paths) == paths
