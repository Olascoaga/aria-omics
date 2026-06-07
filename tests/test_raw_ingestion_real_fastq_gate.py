import json
import os
from pathlib import Path

import pytest


REQUIRED_KB_FIELDS = (
    "fastq_files",
    "index_path",
    "index_sha256",
    "t2g_path",
    "t2g_sha256",
    "chemistry",
    "output_dir",
)


def _expanded(path_value):
    return str(Path(str(path_value)).expanduser())


def _load_validation_config():
    config_path = os.environ.get("ARIA_SCRNA_FASTQ_KB_VALIDATION_JSON")
    if not config_path:
        pytest.skip(
            "Set ARIA_SCRNA_FASTQ_KB_VALIDATION_JSON to run the gated real "
            "scRNA FASTQ -> kb -> h5ad -> QC/clustering validation."
        )

    path = Path(config_path).expanduser()
    if not path.exists():
        pytest.fail(
            f"ARIA_SCRNA_FASTQ_KB_VALIDATION_JSON does not exist: {path}"
        )
    with path.open(encoding="utf-8") as fh:
        config = json.load(fh)
    if not isinstance(config, dict):
        pytest.fail("FASTQ validation config must be a JSON object.")

    kb_params = dict(config.get("raw_ingestion_kb") or config)
    missing = [
        field for field in REQUIRED_KB_FIELDS
        if kb_params.get(field) in (None, "", [])
    ]
    if missing:
        pytest.fail(
            "FASTQ validation config is missing raw_ingestion_kb fields: "
            + ", ".join(missing)
        )

    kb_params["execute"] = True
    kb_params["fastq_files"] = [
        _expanded(path) for path in kb_params.get("fastq_files", [])
    ]
    for field in ("index_path", "t2g_path", "output_dir"):
        kb_params[field] = _expanded(kb_params[field])

    for fastq in kb_params["fastq_files"]:
        if not Path(fastq).exists():
            pytest.fail(f"FASTQ validation file is absent: {fastq}")
    for field in ("index_path", "t2g_path"):
        if not Path(kb_params[field]).exists():
            pytest.fail(f"FASTQ validation {field} is absent: {kb_params[field]}")

    return config, kb_params


def _validation_data_dir(config: dict, kb_params: dict) -> Path:
    if config.get("data_dir"):
        return Path(_expanded(config["data_dir"]))
    fastq_parents = [str(Path(path).resolve().parent)
                     for path in kb_params["fastq_files"]]
    return Path(os.path.commonpath(fastq_parents))


def test_real_scrna_fastq_kb_validation_hands_off_to_scrna_qc_and_clustering(
    tmp_path, monkeypatch
):
    """
    Gated real-data validation for ADR-029 Slice 4.

    Enable with a JSON file:
      {
        "data_dir": "/path/containing/fastqs",
        "organism": "Homo sapiens",
        "raw_ingestion_kb": {
          "fastq_files": ["/path/sample_R1.fastq.gz", "/path/sample_R2.fastq.gz"],
          "index_path": "/path/transcriptome.idx",
          "index_sha256": "...",
          "t2g_path": "/path/t2g.txt",
          "t2g_sha256": "...",
          "chemistry": "10xv3",
          "output_dir": "/tmp/aria_scrna_fastq_kb_validation/kb",
          "threads": 4
        }
      }

    The test intentionally does not infer chemistry/reference assets and skips
    unless the explicit local validation bundle is supplied.
    """
    config, kb_params = _load_validation_config()

    import aria.agents.raw_ingestion_agent as raw_agent
    from aria.agents.raw_ingestion_agent import RawIngestionAgent
    from aria.utils.environment_manager import EnvironmentManager

    data_dir = _validation_data_dir(config, kb_params)
    if not data_dir.exists():
        pytest.fail(f"FASTQ validation data_dir is absent: {data_dir}")

    env = EnvironmentManager(workspace_dir=str(tmp_path / "env_workspace"))
    monkeypatch.setattr(raw_agent.Path, "home", lambda: tmp_path)

    agent = RawIngestionAgent.__new__(RawIngestionAgent)
    agent.env = env
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_blocking_escalation = lambda **kwargs: pytest.fail(
        "real FASTQ validation must use complete explicit raw_ingestion_kb "
        "params and must not publish a metadata checkpoint"
    )

    result = agent.run(
        "real_scrna_fastq_kb_validation",
        {
            "exp_context": {
                "data_dir": str(data_dir),
                "modalities": {},
                "input_files": [],
                "raw_ingestion_kb": kb_params,
            }
        },
    )
    if result.get("status") != "done":
        pytest.fail(json.dumps(result, indent=2, sort_keys=True))

    records = result.get("records", [])
    kb_records = [
        record for record in records
        if record.get("mode") == "fastq_kb_count"
    ]
    assert len(kb_records) == 1
    assert kb_records[0]["status"] == "success"

    updates = result["exp_context_updates"]
    output_h5ads = updates["modalities"]["scRNA"]
    assert len(output_h5ads) == 1
    output_h5ad = Path(output_h5ads[0])
    assert output_h5ad.exists()
    assert updates["input_files"][0]["modality"] == "scRNA_ingested_h5ad"
    assert updates["input_files"][0]["source_mode"] == "fastq_kb_count"

    qc_dir = tmp_path / "scrna_qc"
    qc_params = {
        "data_path": str(output_h5ad),
        "organism": config.get("organism", "Homo sapiens"),
        "sample_id": config.get("sample_id", "kb_validation"),
        "output_dir": str(qc_dir),
        "initial_min_genes": 1,
        "initial_min_cells": 1,
        "min_genes": 1,
        "min_cells": 1,
        "run_scrublet": False,
    }
    qc_params.update(config.get("qc_params") or {})
    qc_result = env.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_qc.py",
        params=qc_params,
        timeout=int(config.get("qc_timeout_seconds") or 1800),
    )
    if qc_result.get("status") != "success":
        pytest.fail(json.dumps(qc_result, indent=2, sort_keys=True))
    qc_output = Path(qc_result["output_path"])
    assert qc_output.exists()

    n_cells = int(qc_result.get("n_cells_after") or 0)
    n_genes = int(qc_result.get("n_genes_after") or 0)
    if n_cells < 3 or n_genes < 3:
        pytest.fail(
            "FASTQ validation fixture produced too few post-QC cells/genes "
            f"for clustering handoff: {n_cells} cells x {n_genes} genes."
        )
    n_pcs = max(2, min(10, n_cells - 1, n_genes - 1))
    cluster_params = {
        "data_path": str(qc_output),
        "resolution": 0.2,
        "n_hvg": min(2000, n_genes),
        "n_pcs": n_pcs,
        "n_neighbors": max(2, min(5, n_cells - 1)),
        "output_dir": str(tmp_path / "scrna_clustering"),
        "seed": 0,
    }
    cluster_params.update(config.get("clustering_params") or {})
    cluster_result = env.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_clustering.py",
        params=cluster_params,
        timeout=int(config.get("clustering_timeout_seconds") or 1800),
    )
    if cluster_result.get("status") != "success":
        pytest.fail(json.dumps(cluster_result, indent=2, sort_keys=True))
    assert Path(cluster_result["output_path"]).exists()
