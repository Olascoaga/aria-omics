"""F9 (preprint audit 2026-06-19): public artifact CSV writes must disclose
failures, not swallow them into a silent ``path=None`` downstream skip.

The DA tables (`chromatin_diffacc`) and the motif enrichment table (`chromatin_motifs`)
route their writes through `persist_csv`; this guard pins its contract so a failed
write surfaces a structured `CsvWriteFailed` warning while still returning None (honest
skip, never a fabricated path).
"""
import pandas as pd

from aria.utils.csv_artifacts import persist_csv


def test_persist_csv_success_returns_path_without_warning(tmp_path):
    warnings = []
    path = str(tmp_path / "ok.csv")
    out = persist_csv(
        path, "ok.csv",
        lambda: pd.DataFrame([{"a": 1, "b": 2}]).to_csv(path, index=False),
        warnings)
    assert out == path
    assert warnings == []
    assert pd.read_csv(path)["a"].tolist() == [1]


def test_persist_csv_failure_records_warning_and_returns_none(tmp_path):
    # parent directory does not exist -> the real pandas write raises, mirroring an
    # unwritable output dir at the chromatin_diffacc call sites.
    warnings = []
    bad = str(tmp_path / "missing_dir" / "x.csv")
    out = persist_csv(
        bad, "chromatin_da_pseudobulk.csv",
        lambda: pd.DataFrame([{"a": 1}]).to_csv(bad, index=False),
        warnings)
    assert out is None
    assert len(warnings) == 1
    assert warnings[0].startswith("CsvWriteFailed")
    assert "chromatin_da_pseudobulk.csv" in warnings[0]


def test_persist_csv_surfaces_writer_exception_type_and_message():
    # mirrors the DictWriter motif site: any write_fn exception is disclosed.
    warnings = []

    def _boom():
        raise OSError("disk full")

    out = persist_csv("/n/a.csv", "chromatin_motif_enrichment.csv", _boom, warnings)
    assert out is None
    assert any(
        "CsvWriteFailed" in w and "OSError" in w and "disk full" in w
        for w in warnings)
