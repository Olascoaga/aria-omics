"""V47 bulk ATAC peak-count matrix guards."""

from pathlib import Path


def test_peak_counts_contract_registered():
    from aria.utils.script_contracts import contract_for_script

    contract = contract_for_script("aria/scripts/chromatin_peak_counts.py")
    assert contract is not None
    assert contract.validation_level == "scaffold"
    assert {f.name for f in contract.success_outputs} >= {
        "counts_matrix_path", "sample_metadata_path", "n_peaks", "n_samples",
    }


def test_peak_count_helpers_parse_bed_and_coverage(tmp_path):
    from aria.scripts.chromatin_peak_counts import (
        _parse_bedtools_counts,
        _read_peak_ids,
        _resolve_sample_ids,
    )

    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\tp1\nchr1\t20\t30\tp2\n")

    assert _read_peak_ids(peaks) == ["chr1:10-20", "chr1:20-30"]
    assert _parse_bedtools_counts(
        "chr1\t10\t20\tp1\t7\nchr1\t20\t30\tp2\t0\n",
        ["chr1:10-20", "chr1:20-30"],
    ) == [7, 0]
    assert _resolve_sample_ids(
        ["/data/A.bam", "/data/B.fragments.tsv.gz"]
    ) == ["A", "B"]


def test_peak_counts_missing_bedtools_is_honest_error(tmp_path, monkeypatch):
    from aria.scripts.chromatin_peak_counts import chromatin_peak_counts

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\n")
    monkeypatch.setattr("shutil.which", lambda name: None)

    res = chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam)],
        "peaks_path": str(peaks),
    })

    assert res["status"] == "error"
    assert res["error_type"] == "MissingDependency"


def test_peak_counts_writes_matrix_and_metadata(tmp_path, monkeypatch):
    import subprocess

    from aria.scripts.chromatin_peak_counts import chromatin_peak_counts

    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\nchr1\t20\t30\n")
    bam_a = tmp_path / "A.bam"
    bam_b = tmp_path / "B.bam"
    bam_a.write_bytes(b"")
    bam_b.write_bytes(b"")

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/bedtools")

    def fake_run(cmd, capture_output, text, check):
        assert cmd[:4] == ["/usr/bin/bedtools", "coverage", "-counts", "-a"]
        stdout = (
            "chr1\t10\t20\t5\nchr1\t20\t30\t1\n"
            if cmd[-1] == str(bam_a)
            else "chr1\t10\t20\t2\nchr1\t20\t30\t9\n"
        )
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    res = chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam_a), str(bam_b)],
        "sample_ids": ["ctrl_1", "treated_1"],
        "sample_metadata": {
            "ctrl_1": {"condition": "ctrl", "replicate": "r1"},
            "treated_1": {"condition": "treated", "replicate": "r1"},
        },
        "peaks_path": str(peaks),
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert res["validation_level"] == "scaffold"
    assert res["n_peaks"] == 2
    assert res["n_samples"] == 2
    matrix = Path(res["counts_matrix_path"]).read_text().splitlines()
    assert matrix == [
        "peak_id\tctrl_1\ttreated_1",
        "chr1:10-20\t5\t2",
        "chr1:20-30\t1\t9",
    ]
    metadata = Path(res["sample_metadata_path"]).read_text()
    assert "condition" in metadata
    assert "replicate" in metadata
    assert res["warnings"] == []


def test_run_ledger_has_peak_count_matrix_node():
    from aria.agents.narrative.run_ledger import build_run_ledger

    ledger = build_run_ledger(
        {"design_intelligence": {"recommended": ["build peak count matrix"]}},
        {"chromatin_agent": {"status": "done", "findings": {
            "bulk_ATAC": {"status": "done", "findings": {
                "peak_counts": {"status": "success"},
            }}
        }}},
    )

    nodes = {entry["node_id"]: entry for entry in ledger["entries"]}
    assert nodes["ledger://chromatin/peak_count_matrix"]["planned"] is True
    assert nodes["ledger://chromatin/peak_count_matrix"]["status"] == "ran"
