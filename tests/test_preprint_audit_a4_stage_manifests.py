"""Preprint-readiness audit A4: per-stage content-addressed resume manifests.

The FASTQ QC (fastp) and alignment (STAR index + align) resume gates validated
only output existence/integrity (`_fastp_outputs_valid`, `_star_output_valid`,
`_index_exists`). A changed input FASTQ, GTF/reference, tool version, or parameter
was therefore never noticed: the stale stage output was silently reused.

After A4 each stage writes a content-addressed manifest (inputs signature + params
hash + tool version); the resume gate reuses a stage only when that manifest still
matches. Changing an input / param / reference / version selectively invalidates
the affected stage. The signature is hybrid: content sha256, with an (size,
mtime) fast path so an unchanged large file is not re-hashed on every resume.
"""
from __future__ import annotations

import pytest


def _mod():
    return pytest.importorskip("aria.utils.stage_manifest")


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_current_manifest_is_reused(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "read-data-v1")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="fastp", inputs=[("r1", fq)],
            params={"quality": 20, "min_len": 36}, tool_version="fastp 0.23.4",
        ),
    )
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq)],
        params={"quality": 20, "min_len": 36}, tool_version="fastp 0.23.4",
    )
    assert current is True, reason


def test_changed_input_content_invalidates(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "read-data-v1")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="fastp", inputs=[("r1", fq)],
            params={"quality": 20}, tool_version="fastp 0.23.4"),
    )
    # Same path, DIFFERENT content and a different size → must invalidate.
    _write(tmp_path / "r1.fq", "read-data-v2-longer")
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq)],
        params={"quality": 20}, tool_version="fastp 0.23.4")
    assert current is False
    assert "input" in reason


def test_same_size_different_content_invalidates_via_rehash(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "AAAA")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="fastp", inputs=[("r1", fq)], params={},
            tool_version="v1"),
    )
    # Rewrite with the SAME byte length but different content and bump mtime so
    # the fast path falls through to a real content re-hash.
    import os
    import time
    _write(tmp_path / "r1.fq", "TTTT")
    st = os.stat(fq)
    os.utime(fq, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000))
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq)], params={}, tool_version="v1")
    assert current is False


def test_changed_param_invalidates(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "data")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="fastp", inputs=[("r1", fq)],
            params={"quality": 20}, tool_version="v1"),
    )
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq)],
        params={"quality": 30}, tool_version="v1")  # quality changed
    assert current is False
    assert "param" in reason


def test_changed_tool_version_invalidates(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "data")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="star", inputs=[("r1", fq)], params={}, tool_version="STAR 2.7.10a"),
    )
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq)], params={}, tool_version="STAR 2.7.11b")
    assert current is False
    assert "version" in reason


def test_changed_reference_input_invalidates(tmp_path):
    # The GTF / reference case: a genome-index manifest over FASTA + GTF must
    # invalidate when the GTF content changes.
    sm = _mod()
    fasta = _write(tmp_path / "genome.fa", ">chr1\nACGT")
    gtf = _write(tmp_path / "genes.gtf", "gene1\texon")
    manifest = tmp_path / "index.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="star_index", inputs=[("fasta", fasta), ("gtf", gtf)],
            params={"sjdbOverhang": 100}, tool_version="STAR 2.7.10a"),
    )
    _write(tmp_path / "genes.gtf", "gene1\texon\ngene2\texon")  # GTF edited
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("fasta", fasta), ("gtf", gtf)],
        params={"sjdbOverhang": 100}, tool_version="STAR 2.7.10a")
    assert current is False
    assert "input" in reason


def test_missing_or_corrupt_manifest_is_not_current(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "data")
    missing = tmp_path / "nope.stage.json"
    current, reason = sm.stage_is_current(
        str(missing), inputs=[("r1", fq)], params={}, tool_version="v1")
    assert current is False
    corrupt = tmp_path / "bad.stage.json"
    corrupt.write_text("{not json")
    current2, reason2 = sm.stage_is_current(
        str(corrupt), inputs=[("r1", fq)], params={}, tool_version="v1")
    assert current2 is False


def test_missing_input_file_is_not_current(tmp_path):
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "data")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="fastp", inputs=[("r1", fq)], params={}, tool_version="v1"),
    )
    import os
    os.remove(fq)
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq)], params={}, tool_version="v1")
    assert current is False


def test_fastp_resume_gate_invalidates_on_param_change(tmp_path):
    """Integration: the fastp resume gate reuses outputs only when the stage
    manifest matches; changing the quality param must stop the skip."""
    import os
    import json as _json
    from aria.scripts import rna_fastq_qc as qc

    trimmed_dir = tmp_path / "trimmed"
    fastp_dir = tmp_path / "fastp"
    trimmed_dir.mkdir()
    fastp_dir.mkdir()

    r1_in = tmp_path / "sampleA_R1.fastq.gz"
    r1_in.write_bytes(os.urandom(4096))

    name = "sampleA"
    # _fastp_outputs_valid only checks file size (> 1 KB), not gzip validity.
    r1_out = trimmed_dir / f"{name}_R1_trimmed.fq.gz"
    r1_out.write_bytes(os.urandom(4096))   # > 1 KB
    json_out = fastp_dir / f"{name}_fastp.json"
    json_out.write_text(_json.dumps({"summary": {
        "before_filtering": {"total_reads": 1000},
        "after_filtering": {"total_reads": 950, "q30_rate": 0.95},
    }}))

    sample = {"name": name, "r1": str(r1_in), "r2": None, "paired": False}
    warnings: list = []

    # First: manifest matches the (min_len=36, quality=20) params → resume skip.
    qc.write_stage_manifest(
        str(fastp_dir / f"{name}.fastp.stage.json"),
        qc.build_stage_manifest(
            stage="fastp", inputs=[("r1", str(r1_in)), ("r2", None)],
            params={"min_len": 36, "quality": 20, "paired": False},
            tool_version=qc._fastp_version(),
        ),
    )
    resumed = qc._run_fastp(sample, trimmed_dir, fastp_dir, 4, 36, 20, warnings)
    assert resumed.get("resumed") is True

    # Now change the quality param → manifest no longer matches → must NOT skip.
    # fastp is not installed in aria-env, so a non-resume path fails/mocks out;
    # the key assertion is that it did not blindly reuse the stale output.
    warnings2: list = []
    rerun = qc._run_fastp(sample, trimmed_dir, fastp_dir, 4, 36, 30, warnings2,
                          allow_mock=True)
    assert not rerun.get("resumed")


def test_none_paths_are_skipped(tmp_path):
    # Single-end: r2 is None and must be ignored, not crash.
    sm = _mod()
    fq = _write(tmp_path / "r1.fq", "data")
    manifest = tmp_path / "s.stage.json"
    sm.write_stage_manifest(
        str(manifest),
        sm.build_stage_manifest(
            stage="fastp", inputs=[("r1", fq), ("r2", None)],
            params={}, tool_version="v1"),
    )
    current, reason = sm.stage_is_current(
        str(manifest), inputs=[("r1", fq), ("r2", None)],
        params={}, tool_version="v1")
    assert current is True, reason
