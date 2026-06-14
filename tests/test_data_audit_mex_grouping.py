"""E2E-1 (production verification 2026-06-12): a raw 10x MEX directory must be
treated as ONE scRNA sample, not split into its component files. Otherwise the
per-sample QC path tries to load barcodes.tsv as a sample and the whole pipeline
aborts with PerSampleQCFailed."""

from aria.agents.data_audit_agent import (
    _collapse_mex_directories,
    _incomplete_mex_warnings,
)


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


# ── T5 (tri-auditoría 2026-06-14): collapsing on matrix.mtx ALONE promotes an
# INCOMPLETE MEX directory (missing barcodes/features) to a sample, which then
# fails late in rna_qc (sc.read_10x_mtx) instead of at CP1. Only complete triplets
# may collapse; an incomplete MEX must be left untouched AND reported at CP1 with
# the missing component named.

def test_incomplete_mex_matrix_only_is_not_collapsed():
    paths = ["/data/partial/matrix.mtx"]
    # Not a complete triplet -> must NOT collapse to the directory.
    assert _collapse_mex_directories(paths) == paths


def test_incomplete_mex_matrix_plus_barcodes_without_features_is_not_collapsed():
    paths = ["/data/partial/matrix.mtx", "/data/partial/barcodes.tsv"]
    out = _collapse_mex_directories(paths)
    # Missing features/genes -> left as component files, not the directory.
    assert "/data/partial" not in out
    assert set(out) == set(paths)


def test_complete_triplet_collapses_and_emits_no_warning():
    paths = [
        "/data/ok/matrix.mtx", "/data/ok/barcodes.tsv", "/data/ok/genes.tsv",
    ]
    assert _collapse_mex_directories(paths) == ["/data/ok"]
    assert _incomplete_mex_warnings(paths) == []


def test_incomplete_mex_warning_names_missing_components():
    # matrix.mtx only -> barcodes AND features/genes missing.
    warns = _incomplete_mex_warnings(["/data/partial/matrix.mtx"])
    assert len(warns) == 1
    text = warns[0].lower()
    assert "/data/partial" in warns[0]
    assert "barcodes" in text
    assert "features" in text or "genes" in text


def test_incomplete_mex_warning_names_only_missing_features():
    warns = _incomplete_mex_warnings(
        ["/data/p2/matrix.mtx", "/data/p2/barcodes.tsv"])
    assert len(warns) == 1
    text = warns[0].lower()
    assert "features" in text or "genes" in text
    # barcodes is present, so it must not be reported missing.
    assert "barcodes" not in text


def test_gzipped_components_count_toward_completeness():
    paths = [
        "/data/g/matrix.mtx.gz", "/data/g/barcodes.tsv.gz",
        "/data/g/features.tsv.gz",
    ]
    assert _collapse_mex_directories(paths) == ["/data/g"]
    assert _incomplete_mex_warnings(paths) == []
