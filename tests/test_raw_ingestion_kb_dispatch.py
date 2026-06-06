from pathlib import Path


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
        "stack": "rna",
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
