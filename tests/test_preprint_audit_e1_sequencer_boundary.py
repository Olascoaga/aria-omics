"""Preprint-readiness audit E1: FASTQ is the supported ingestion boundary.

ARIA starts at most from FASTQ — it does no demultiplexing (BCL) or basecalling
(POD5/FAST5). Before E1 those sequencer-native inputs matched no scanner and fell
through to a generic ``status="skipped"`` / "No supported raw-ingestion inputs
detected", indistinguishable from an empty directory and offering no guidance.

After E1 a pure detector recognises unsupported sequencer-native formats and
`RawIngestionAgent` fails closed with `UnsupportedSequencerFormat`, naming the
required upstream stage and pointing at FASTQ as the supported entry boundary.
"""
from __future__ import annotations

import pytest


# ── The pure detector ─────────────────────────────────────────────────────────

def _detect():
    from aria.utils.sequencer_formats import detect_unsupported_inputs
    return detect_unsupported_inputs


@pytest.mark.parametrize("path,fmt", [
    ("/data/lane1/L001_C1.1.bcl", "BCL"),
    ("/data/run/0001.cbcl", "BCL"),
    ("/data/run/RunInfo.xml", "BCL"),
    ("/data/run/RTAComplete.txt", "BCL"),
    ("/data/run/Data/Intensities/BaseCalls/", "BCL"),
    ("/data/reads/signal.pod5", "POD5"),
    ("/data/reads/batch0.fast5", "FAST5"),
])
def test_unsupported_formats_are_detected(path, fmt):
    hits = _detect()([path])
    assert len(hits) == 1
    assert hits[0]["format"] == fmt
    assert hits[0]["path"] == path
    assert hits[0].get("upstream_stage")  # names the required upstream tool


@pytest.mark.parametrize("path", [
    "/data/sampleA_R1.fastq.gz",
    "/data/sampleA_R2.fq",
    "/data/matrix.h5ad",
    "/data/filtered_feature_bc_matrix/matrix.mtx.gz",
    "/data/notes.txt",
])
def test_supported_and_neutral_paths_not_flagged(path):
    assert _detect()([path]) == []


def test_detector_dedupes_and_handles_empty():
    assert _detect()([]) == []
    hits = _detect()(["/x/a.pod5", "/x/a.pod5"])
    assert len(hits) == 1


# ── The RawIngestionAgent boundary gate ───────────────────────────────────────

def _agent():
    pytest.importorskip("litellm")
    from aria.agents.raw_ingestion_agent import RawIngestionAgent
    return RawIngestionAgent


def test_agent_rejects_pod5_with_fastq_boundary_message(tmp_path):
    RawIngestionAgent = _agent()
    exp_ctx = {"modalities": {"scRNA": [str(tmp_path / "run.pod5")]}}
    rejection = RawIngestionAgent._unsupported_input_rejection(exp_ctx, str(tmp_path))
    assert rejection is not None
    assert rejection["error_type"] == "UnsupportedSequencerFormat"
    assert rejection["supported_boundary"] == "FASTQ"
    assert any(u["format"] == "POD5" for u in rejection["unsupported_inputs"])
    assert "FASTQ" in rejection["details"]


def test_agent_rejects_illumina_run_folder(tmp_path):
    RawIngestionAgent = _agent()
    # A BCL run folder discovered under data_dir (not a declared modality file).
    basecalls = tmp_path / "Data" / "Intensities" / "BaseCalls"
    basecalls.mkdir(parents=True)
    (tmp_path / "RunInfo.xml").write_text("<RunInfo/>")
    rejection = RawIngestionAgent._unsupported_input_rejection(
        {"modalities": {}}, str(tmp_path))
    assert rejection is not None
    assert rejection["error_type"] == "UnsupportedSequencerFormat"
    assert any(u["format"] == "BCL" for u in rejection["unsupported_inputs"])


def test_agent_does_not_reject_fastq_inputs(tmp_path):
    RawIngestionAgent = _agent()
    exp_ctx = {"modalities": {"scRNA": [str(tmp_path / "s_R1.fastq.gz")]}}
    rejection = RawIngestionAgent._unsupported_input_rejection(exp_ctx, str(tmp_path))
    assert rejection is None
