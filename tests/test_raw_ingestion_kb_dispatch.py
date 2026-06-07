from pathlib import Path
import json


def test_rna_kb_count_maps_incomplete_inputs_to_structured_error(tmp_path):
    from aria.scripts.rna_kb_count import rna_kb_count

    result = rna_kb_count({
        "fastq_files": [],
        "output_dir": str(tmp_path / "kb"),
    })

    assert result["status"] == "error"
    assert result["error_type"] == "KbInputBlocked"
    assert "FASTQ files" in result["details"]
    assert result["mode"] == "fastq_kb_count"


def test_raw_ingestion_agent_dispatches_kb_count_through_rna_stack(
    tmp_path, monkeypatch
):
    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    output_h5ad = tmp_path / "kb_count.h5ad"
    output_h5ad.write_bytes(b"fake-h5ad")
    calls = []

    class FakeEnv:
        def run_in_stack(self, *, stack, script_path, params, timeout=None):
            calls.append({
                "stack": stack,
                "script_path": script_path,
                "params": params,
                "timeout": timeout,
            })
            return {
                "status": "success",
                "mode": "fastq_kb_count",
                "output_h5ad": str(output_h5ad),
                "output_sha256": "abc123",
                "command": ["kb", "count"],
                "reference": {"chemistry": "10xv3"},
            }

    monkeypatch.setattr(raw_agent, "discover_10x_mtx_triplets", lambda root: [])
    monkeypatch.setattr(raw_agent, "scan_fastq_plan", lambda root: {
        "status": "blocked",
        "mode": "fastq_kb_plan",
        "fastq_count": 2,
        "samples": [{"sample_id": "sample1"}],
        "blockers": ["FASTQ ingestion requires explicit chemistry."],
    })
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = FakeEnv()
    agent.publish_status = lambda *args, **kwargs: None

    result = agent.run(
        "exp123",
        {
            "exp_context": {
                "data_dir": str(tmp_path),
                "modalities": {},
                "input_files": [],
                "raw_ingestion_kb": {
                    "execute": True,
                    "fastq_files": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
                    "index_path": str(tmp_path / "idx"),
                    "index_sha256": "idx-sha",
                    "t2g_path": str(tmp_path / "t2g.txt"),
                    "t2g_sha256": "t2g-sha",
                    "chemistry": "10xv3",
                    "output_dir": str(tmp_path / "kb"),
                    "timeout_seconds": 123,
                },
            }
        },
    )

    assert result["status"] == "done"
    assert calls == [{
        "stack": "ingestion",
        "script_path": "aria/scripts/rna_kb_count.py",
        "params": {
            "execute": True,
            "fastq_files": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
            "index_path": str(tmp_path / "idx"),
            "index_sha256": "idx-sha",
            "t2g_path": str(tmp_path / "t2g.txt"),
            "t2g_sha256": "t2g-sha",
            "chemistry": "10xv3",
            "output_dir": str(tmp_path / "kb"),
            "timeout_seconds": 123,
        },
        "timeout": 123,
    }]
    updates = result["exp_context_updates"]
    assert updates["modalities"]["scRNA"] == [str(output_h5ad)]
    assert updates["input_files"][0]["source_mode"] == "fastq_kb_count"
    assert updates["input_files"][0]["path"] == str(output_h5ad)


def test_raw_ingestion_fastq_checkpoint_skips_without_explicit_kb_params(
    tmp_path, monkeypatch
):
    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    class FailEnv:
        def run_in_stack(self, **kwargs):
            raise AssertionError("kb count should not run when checkpoint skips")

    monkeypatch.setattr(raw_agent, "discover_10x_mtx_triplets", lambda root: [])
    monkeypatch.setattr(raw_agent, "scan_fastq_plan", lambda root: {
        "status": "blocked",
        "mode": "fastq_kb_plan",
        "fastq_count": 2,
        "samples": [{"sample_id": "sample1", "files": {"R1": [], "R2": []}}],
        "blockers": ["FASTQ ingestion requires explicit chemistry."],
    })
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = FailEnv()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_blocking_escalation = lambda **kwargs: (
        None,
        {"choice": "Skip FASTQ quantification for now (Recommended)"},
    )

    result = agent.run(
        "exp123",
        {"exp_context": {"data_dir": str(tmp_path), "modalities": {}}},
    )

    assert result["status"] == "done"
    updates = result["exp_context_updates"]
    assert updates["modalities"] == {}
    checkpoint = updates["raw_ingestion"][1]
    assert checkpoint["mode"] == "fastq_kb_checkpoint"
    assert checkpoint["status"] == "skipped"


def test_raw_ingestion_fastq_checkpoint_json_enables_kb_dispatch(
    tmp_path, monkeypatch
):
    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    output_h5ad = tmp_path / "kb_count.h5ad"
    output_h5ad.write_bytes(b"fake-h5ad")
    calls = []

    class FakeEnv:
        def run_in_stack(self, *, stack, script_path, params, timeout=None):
            calls.append((stack, script_path, params, timeout))
            return {
                "status": "success",
                "mode": "fastq_kb_count",
                "output_h5ad": str(output_h5ad),
                "output_sha256": "abc123",
                "command": ["kb", "count"],
                "reference": {"chemistry": params["chemistry"]},
            }

    kb_json = {
        "fastq_files": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
        "index_path": str(tmp_path / "transcriptome.idx"),
        "index_sha256": "idx-sha",
        "t2g_path": str(tmp_path / "t2g.txt"),
        "t2g_sha256": "t2g-sha",
        "chemistry": "10xv3",
        "output_dir": str(tmp_path / "kb"),
    }

    monkeypatch.setattr(raw_agent, "discover_10x_mtx_triplets", lambda root: [])
    monkeypatch.setattr(raw_agent, "scan_fastq_plan", lambda root: {
        "status": "blocked",
        "mode": "fastq_kb_plan",
        "fastq_count": 2,
        "samples": [{"sample_id": "sample1", "files": {"R1": [], "R2": []}}],
        "blockers": ["FASTQ ingestion requires explicit chemistry."],
    })
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = FakeEnv()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_blocking_escalation = lambda **kwargs: (
        None,
        {"choice": json.dumps(kb_json)},
    )

    result = agent.run(
        "exp123",
        {"exp_context": {"data_dir": str(tmp_path), "modalities": {}}},
    )

    assert result["status"] == "done"
    assert calls[0][0] == "ingestion"
    assert calls[0][1] == "aria/scripts/rna_kb_count.py"
    assert calls[0][2]["execute"] is True
    assert calls[0][2]["chemistry"] == "10xv3"
    assert result["exp_context_updates"]["modalities"]["scRNA"] == [str(output_h5ad)]


def test_raw_ingestion_blocks_scrna_fastq_only_when_kb_is_skipped(
    tmp_path, monkeypatch
):
    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    class FailEnv:
        def run_in_stack(self, **kwargs):
            raise AssertionError("kb count should not run when checkpoint skips")

    monkeypatch.setattr(raw_agent, "discover_10x_mtx_triplets", lambda root: [])
    monkeypatch.setattr(raw_agent, "scan_fastq_plan", lambda root: {
        "status": "blocked",
        "mode": "fastq_kb_plan",
        "fastq_count": 2,
        "samples": [{"sample_id": "sample1", "files": {"R1": [], "R2": []}}],
        "blockers": ["FASTQ ingestion requires explicit chemistry."],
    })
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = FailEnv()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_blocking_escalation = lambda **kwargs: (
        None,
        {"choice": "Skip FASTQ quantification for now (Recommended)"},
    )

    result = agent.run(
        "exp123",
        {
            "exp_context": {
                "data_dir": str(tmp_path),
                "modalities": {
                    "scRNA": [
                        str(tmp_path / "sample_R1.fastq.gz"),
                        str(tmp_path / "sample_R2.fastq.gz"),
                    ]
                },
            }
        },
    )

    assert result["status"] == "error"
    assert result["error_type"] == "RawIngestionFailed"
    assert "no canonical .h5ad was generated" in result["details"]
    assert "will not dispatch scRNAAgent on raw FASTQ files" in result["details"]
    assert result["records"][1]["status"] == "skipped"
    assert result["errors"][-1]["error_type"] == "CanonicalH5adMissing"


def test_raw_ingestion_blocks_scrna_fastq_only_when_kb_fails(
    tmp_path, monkeypatch
):
    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    class ErrorEnv:
        def run_in_stack(self, *, stack, script_path, params, timeout=None):
            return {
                "status": "error",
                "mode": "fastq_kb_count",
                "error_type": "KbCountFailed",
                "details": "kb count exited with return code 1.",
            }

    monkeypatch.setattr(raw_agent, "discover_10x_mtx_triplets", lambda root: [])
    monkeypatch.setattr(raw_agent, "scan_fastq_plan", lambda root: {
        "status": "blocked",
        "mode": "fastq_kb_plan",
        "fastq_count": 2,
        "samples": [{"sample_id": "sample1"}],
        "blockers": [],
    })
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = ErrorEnv()
    agent.publish_status = lambda *args, **kwargs: None

    result = agent.run(
        "exp123",
        {
            "exp_context": {
                "data_dir": str(tmp_path),
                "modalities": {
                    "scRNA": [
                        str(tmp_path / "sample_R1.fastq.gz"),
                        str(tmp_path / "sample_R2.fastq.gz"),
                    ]
                },
                "raw_ingestion_kb": {
                    "execute": True,
                    "fastq_files": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
                    "index_path": str(tmp_path / "idx"),
                    "index_sha256": "idx-sha",
                    "t2g_path": str(tmp_path / "t2g.txt"),
                    "t2g_sha256": "t2g-sha",
                    "chemistry": "10xv3",
                    "output_dir": str(tmp_path / "kb"),
                },
            }
        },
    )

    assert result["status"] == "error"
    assert result["errors"][0]["error_type"] == "KbCountFailed"
    assert "kb count exited with return code 1" in result["details"]
    assert result["errors"][-1]["error_type"] == "CanonicalH5adMissing"
