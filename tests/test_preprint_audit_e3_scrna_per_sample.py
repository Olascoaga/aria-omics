from pathlib import Path


def _sample_plan(tmp_path: Path) -> dict:
    def paths(sample: str, role: str) -> list[str]:
        return [
            str(tmp_path / f"{sample}_S1_L001_{role}_001.fastq.gz"),
            str(tmp_path / f"{sample}_S1_L002_{role}_001.fastq.gz"),
        ]

    samples = []
    for sample in ("sample_a", "sample_b"):
        samples.append({
            "sample_id": sample,
            "files": {
                "I1": paths(sample, "I1"),
                "R1": paths(sample, "R1"),
                "R2": paths(sample, "R2"),
            },
            "lanes": ["L001", "L002"],
        })
    return {
        "status": "blocked",
        "mode": "fastq_kb_plan",
        "fastq_count": 12,
        "samples": samples,
        "blockers": ["FASTQ ingestion requires explicit chemistry."],
    }


def test_e3_quantifies_each_sample_then_unions_with_explicit_identity(
    tmp_path, monkeypatch
):
    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    plan = _sample_plan(tmp_path)
    selected_fastqs = []
    for sample in plan["samples"]:
        selected_fastqs.extend(sample["files"]["R1"])
        selected_fastqs.extend(sample["files"]["R2"])

    per_sample_outputs = {
        sample: tmp_path / f"{sample}.h5ad"
        for sample in ("sample_a", "sample_b")
    }
    for output in per_sample_outputs.values():
        output.write_bytes(b"per-sample-h5ad")
    union_output = tmp_path / "union" / "concatenated.h5ad"
    union_output.parent.mkdir()
    union_output.write_bytes(b"union-h5ad")

    calls = []

    class FakeEnv:
        def run_in_stack(self, *, stack, script_path, params, timeout=None):
            calls.append({
                "stack": stack,
                "script_path": script_path,
                "params": params,
                "timeout": timeout,
            })
            if script_path == "aria/scripts/rna_kb_count.py":
                sample_id = params["sample_id"]
                return {
                    "status": "success",
                    "mode": "fastq_kb_count",
                    "sample_id": sample_id,
                    "library": sample_id,
                    "output_h5ad": str(per_sample_outputs[sample_id]),
                    "output_sha256": f"sha-{sample_id}",
                    "command": ["kb", "count", *params["fastq_files"]],
                    "reference": {"chemistry": params["chemistry"]},
                }
            assert script_path == "aria/scripts/rna_concat.py"
            return {
                "status": "success",
                "output_path": str(union_output),
                "n_samples": 2,
                "n_cells_total": 5,
                "n_genes_shared": 2,
                "per_sample": [
                    {"sample_id": "sample_a", "n_cells": 2, "n_genes": 2},
                    {"sample_id": "sample_b", "n_cells": 3, "n_genes": 2},
                ],
                "batch_col": "batch",
                "join": "inner",
            }

    monkeypatch.setattr(raw_agent, "discover_10x_mtx_triplets", lambda root: [])
    monkeypatch.setattr(raw_agent, "scan_fastq_plan", lambda root: plan)
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = FakeEnv()
    agent.publish_status = lambda *args, **kwargs: None

    result = agent.run(
        "exp-e3",
        {
            "exp_context": {
                "data_dir": str(tmp_path),
                "modalities": {},
                "input_files": [],
                "raw_ingestion_kb": {
                    "execute": True,
                    "fastq_files": selected_fastqs,
                    "index_path": str(tmp_path / "transcriptome.idx"),
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
    kb_calls = [
        call for call in calls
        if call["script_path"] == "aria/scripts/rna_kb_count.py"
    ]
    assert [call["params"]["sample_id"] for call in kb_calls] == [
        "sample_a", "sample_b"
    ]
    for call in kb_calls:
        sample_id = call["params"]["sample_id"]
        fastqs = call["params"]["fastq_files"]
        assert len(fastqs) == 4
        assert all(sample_id in path for path in fastqs)
        assert all("_I1_" not in path for path in fastqs)
        assert [Path(path).name for path in fastqs] == [
            f"{sample_id}_S1_L001_R1_001.fastq.gz",
            f"{sample_id}_S1_L001_R2_001.fastq.gz",
            f"{sample_id}_S1_L002_R1_001.fastq.gz",
            f"{sample_id}_S1_L002_R2_001.fastq.gz",
        ]
        assert call["params"]["output_dir"].endswith(sample_id)

    union_calls = [
        call for call in calls
        if call["script_path"] == "aria/scripts/rna_concat.py"
    ]
    assert len(union_calls) == 1
    assert union_calls[0]["stack"] == "ingestion"
    assert union_calls[0]["params"]["samples"] == [
        {
            "path": str(per_sample_outputs["sample_a"]),
            "sample_id": "sample_a",
            "library": "sample_a",
        },
        {
            "path": str(per_sample_outputs["sample_b"]),
            "sample_id": "sample_b",
            "library": "sample_b",
        },
    ]
    assert result["output_h5ads"] == [str(union_output)]
    updates = result["exp_context_updates"]
    assert updates["modalities"]["scRNA"] == [str(union_output)]
    union_record = next(
        record for record in updates["raw_ingestion"]
        if record.get("mode") == "fastq_kb_union"
    )
    assert union_record["sample_manifest"] == union_calls[0]["params"]["samples"]
    assert union_record["n_samples"] == 2


def test_e3_explicit_union_preserves_per_sample_counts(tmp_path):
    import anndata as ad
    import numpy as np

    from aria.scripts.rna_concat import rna_concat

    sample_a = ad.AnnData(
        X=np.array([[1, 2], [3, 4]], dtype=np.int64),
    )
    sample_a.obs_names = ["cell1", "cell2"]
    sample_a.var_names = ["gene1", "gene2"]
    sample_b = ad.AnnData(
        X=np.array([[5, 6], [7, 8], [9, 10]], dtype=np.int64),
    )
    sample_b.obs_names = ["cell1", "cell2", "cell3"]
    sample_b.var_names = ["gene1", "gene2"]

    path_a = tmp_path / "sample_a.h5ad"
    path_b = tmp_path / "sample_b.h5ad"
    sample_a.write_h5ad(path_a)
    sample_b.write_h5ad(path_b)

    result = rna_concat({
        "samples": [
            {"path": str(path_a), "sample_id": "sample_a", "library": "lib_a"},
            {"path": str(path_b), "sample_id": "sample_b", "library": "lib_b"},
        ],
        "output_dir": str(tmp_path / "union"),
        "join": "inner",
    })

    assert result["status"] == "success"
    combined = ad.read_h5ad(result["output_path"])
    assert combined.obs["sample_id"].tolist() == [
        "sample_a", "sample_a", "sample_b", "sample_b", "sample_b"
    ]
    assert combined.obs["library"].tolist() == [
        "lib_a", "lib_a", "lib_b", "lib_b", "lib_b"
    ]
    values = np.asarray(combined.X)
    assert int(values[:2].sum()) == int(np.asarray(sample_a.X).sum())
    assert int(values[2:].sum()) == int(np.asarray(sample_b.X).sum())


def test_e3_refuses_partial_lane_selection(tmp_path):
    from aria.agents.raw_ingestion_agent import RawIngestionAgent

    plan = _sample_plan(tmp_path)
    selected_fastqs = []
    for sample in plan["samples"]:
        selected_fastqs.extend(sample["files"]["R1"])
        selected_fastqs.extend(sample["files"]["R2"])
    selected_fastqs.remove(plan["samples"][1]["files"]["R2"][1])

    params, error = RawIngestionAgent._prepare_sample_kb_params(
        plan,
        {
            "fastq_files": selected_fastqs,
            "output_dir": str(tmp_path / "kb"),
        },
    )

    assert params == []
    assert error["status"] == "error"
    assert error["error_type"] == "KbSamplePlanInvalid"
    assert error["reason"] == "unbalanced_fastq_sample_selection"
    assert error["sample_id"] == "sample_b"
    assert error["n_r1"] == 2
    assert error["n_r2"] == 1
