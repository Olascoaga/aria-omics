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
    # Empty stub BAMs have no readable header -> honest fallback warning.
    assert any("genome" in w for w in res["warnings"])


def test_peak_counts_uses_sorted_low_memory_coverage(tmp_path, monkeypatch):
    """Large bulk ATAC BAMs OOM `bedtools coverage` unless the low-memory
    sweeping algorithm (`-sorted -g`) is used. When a genome can be read from
    the alignment header, counting must engage `-sorted` and emit peaks in
    genome order. Counts are identical to the non-sorted call."""
    import subprocess

    from aria.scripts import chromatin_peak_counts as mod

    peaks = tmp_path / "peaks.bed"
    # Deliberately NOT in genome order, to prove we sort before -sorted counting.
    peaks.write_text("chr1\t20\t30\nchr1\t10\t20\n")
    bam = tmp_path / "A.bam"
    bam.write_bytes(b"")

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/bedtools")
    monkeypatch.setattr(mod, "_genome_sizes_from_bam", lambda f: [("chr1", 1000)])

    seen = {}

    def fake_run(cmd, capture_output, text, check):
        seen["cmd"] = cmd
        # bedtools emits one row per -a feature, in the (sorted) -a order.
        return subprocess.CompletedProcess(
            cmd, 0, stdout="chr1\t10\t20\t3\nchr1\t20\t30\t8\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    res = mod.chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam)],
        "sample_ids": ["s1"],
        "peaks_path": str(peaks),
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert "-sorted" in seen["cmd"] and "-g" in seen["cmd"]
    assert "sorted" in res["counting_method"]
    matrix = Path(res["counts_matrix_path"]).read_text().splitlines()
    assert matrix == ["peak_id\ts1", "chr1:10-20\t3", "chr1:20-30\t8"]


def test_peak_counts_falls_back_to_unsorted_without_genome(tmp_path, monkeypatch):
    """When no genome can be read from the header, counting falls back to the
    plain (non-sorted) call rather than failing."""
    import subprocess

    from aria.scripts import chromatin_peak_counts as mod

    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\n")
    bam = tmp_path / "A.bam"
    bam.write_bytes(b"")

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/bedtools")
    monkeypatch.setattr(mod, "_genome_sizes_from_bam", lambda f: None)

    def fake_run(cmd, capture_output, text, check):
        assert cmd[:4] == ["/usr/bin/bedtools", "coverage", "-counts", "-a"]
        return subprocess.CompletedProcess(cmd, 0, stdout="chr1\t10\t20\t4\n",
                                           stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    res = mod.chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam)],
        "sample_ids": ["s1"],
        "peaks_path": str(peaks),
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert res["counting_method"] == "bedtools coverage -counts"


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
