"""Preprint-readiness audit B5: total metadata correspondence.

Partial metadata must never silently reduce the analysis.  When an explicit
metadata TSV covers only a subset of the count-matrix columns (>=2 matching, so
the loader accepts it), the unmatched columns used to be dropped without a trace:
``_load_or_infer_metadata`` returned ``meta.loc[common]`` and ``_run_deseq2`` then
subset ``counts[meta_sub.index]``, so DE ran on a reduced sample set with no error
and no disclosure.

After B5 every count column must have an aligned metadata row (total
correspondence), OR be named in an explicit, audited ``excluded_samples`` list.
An unmatched column that is not explicitly excluded fails closed with the exact
list; only genuine orphan columns may be excluded (excluding a sample that has
metadata, or a name that is not a count column, is a misuse and also fails
closed).
"""
from __future__ import annotations

import pytest

pd = pytest.importorskip("pandas")


def _counts(cols):
    return pd.DataFrame({c: [1, 2, 3] for c in cols}, index=["g1", "g2", "g3"])


def _meta(samples, factor="condition"):
    return pd.DataFrame(
        {factor: ["ctrl" if i % 2 == 0 else "treat" for i in range(len(samples))]},
        index=list(samples),
    )


def _enforce():
    from aria.scripts.rna_bulk_de import _enforce_metadata_correspondence
    return _enforce_metadata_correspondence


def test_full_correspondence_passes():
    fn = _enforce()
    counts = _counts(["s1", "s2", "s3", "s4"])
    meta = _meta(["s1", "s2", "s3", "s4"])
    kept, disclosure = fn(counts, meta, None)
    assert list(kept.columns) == ["s1", "s2", "s3", "s4"]
    assert disclosure["n_excluded"] == 0
    assert disclosure["n_analyzed"] == 4


def test_partial_metadata_errors_with_exact_list():
    fn = _enforce()
    counts = _counts(["s1", "s2", "s3", "s4"])
    meta = _meta(["s1", "s2", "s3"])  # s4 has no metadata row
    with pytest.raises(ValueError) as exc:
        fn(counts, meta, None)
    # The unmatched column must be named explicitly; never a silent reduction.
    assert "s4" in str(exc.value)


def test_explicit_exclusion_drops_only_named_orphan():
    fn = _enforce()
    counts = _counts(["s1", "s2", "s3", "s4"])
    meta = _meta(["s1", "s2", "s3"])
    kept, disclosure = fn(counts, meta, ["s4"])
    assert list(kept.columns) == ["s1", "s2", "s3"]
    assert disclosure["excluded_samples"] == ["s4"]
    assert disclosure["n_excluded"] == 1
    assert disclosure["n_analyzed"] == 3


def test_excluding_a_sample_with_metadata_is_misuse():
    fn = _enforce()
    counts = _counts(["s1", "s2", "s3", "s4"])
    meta = _meta(["s1", "s2", "s3", "s4"])
    # s2 HAS metadata; excluded_samples is only for orphan columns, so this is a
    # misuse and must fail closed rather than silently drop a described sample.
    with pytest.raises(ValueError) as exc:
        fn(counts, meta, ["s2"])
    assert "s2" in str(exc.value)


def test_excluding_a_nonexistent_column_is_misuse():
    fn = _enforce()
    counts = _counts(["s1", "s2", "s3"])
    meta = _meta(["s1", "s2"])
    # "s9" is not a count column at all; only genuine orphans (here s3) qualify.
    with pytest.raises(ValueError) as exc:
        fn(counts, meta, ["s9"])
    assert "s9" in str(exc.value)


def test_partial_exclusion_still_errors_on_remaining_orphans():
    fn = _enforce()
    counts = _counts(["s1", "s2", "s3", "s4", "s5"])
    meta = _meta(["s1", "s2", "s3"])  # s4 AND s5 are orphans
    # Excluding only s4 leaves s5 unauthorized → still fails closed on s5.
    with pytest.raises(ValueError) as exc:
        fn(counts, meta, ["s4"])
    assert "s5" in str(exc.value)


def test_run_reports_metadata_correspondence_error(tmp_path):
    """End-to-end: bulk_rna_de returns a MetadataCorrespondenceError (before any
    DE) when the metadata TSV covers only a subset of the count columns."""
    from aria.scripts.rna_bulk_de import bulk_rna_de

    counts_path = tmp_path / "counts.tsv"
    # More genes than samples so the orientation heuristic does not transpose.
    rows = "\n".join(
        f"g{i}\t{10 + i}\t{12 + i}\t{11 + i}\t{9 + i}" for i in range(1, 9)
    )
    counts_path.write_text("gene\ts1\ts2\ts3\ts4\n" + rows + "\n")
    meta_path = tmp_path / "meta.tsv"
    # s4 deliberately omitted.
    meta_path.write_text(
        "sample\tcondition\n"
        "s1\tctrl\n"
        "s2\tctrl\n"
        "s3\ttreat\n"
    )
    result = bulk_rna_de({
        "files": [str(counts_path)],
        "metadata_file": str(meta_path),
        "design_factor": "condition",
        "comparison": {"numerator": "treat", "denominator": "ctrl"},
        "output_dir": str(tmp_path / "out"),
        "run_pathways": False,
    })
    assert result["status"] == "error"
    assert result["error_type"] == "MetadataCorrespondenceError"
    assert "s4" in result["details"]
