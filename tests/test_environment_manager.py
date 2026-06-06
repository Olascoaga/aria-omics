"""ARIA EnvironmentManager — native pytest (P1-11 follow-up).

Converted from the legacy script-style harness (top-level ok()/fail()/sys.exit,
which crashed pytest collection and was never run in CI). Validates the IPC
architecture without real conda envs: JSON contract, structured errors, fallback
logic, and the script base contract. Heavy/live cases are scanpy + dataset gated.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

# EnvironmentManager -> aria.utils.script_contracts imports pydantic (a core dep,
# present in the PR lane via `pip install -e .`). Skip cleanly where it is absent
# (e.g. the slim aria-rna-env heavy lane) instead of erroring.
pytest.importorskip("pydantic")


# ── EnvironmentManager surface + status report ───────────────────────────────

def test_env_manager_has_required_methods():
    from aria.utils.environment_manager import env_manager
    for m in ("run_in_stack", "check_environments", "get_status_report"):
        assert hasattr(env_manager, m)


def test_status_report_keys():
    from aria.utils.environment_manager import env_manager
    report = env_manager.get_status_report()
    for key in ("conda_available", "environments", "ready_stacks", "missing_stacks"):
        assert key in report


def test_invalid_stack_returns_structured_error():
    from aria.utils.environment_manager import env_manager
    result = env_manager.run_in_stack(
        stack="invalid_stack", script_path="nonexistent.py", params={},
    )
    assert result["status"] == "error"
    assert result["error_type"] == "UnknownStack"


def test_check_environments_shape():
    from aria.utils.environment_manager import env_manager
    envs = env_manager.check_environments()
    assert isinstance(envs, dict)
    assert {"rna", "ingestion", "chromatin", "hic", "integration"} <= set(envs)
    assert all(isinstance(v, bool) for v in envs.values())


def test_ingestion_stack_has_dedicated_environment():
    from aria.utils.environment_manager import EnvironmentManager
    assert EnvironmentManager.STACKS["ingestion"] == "aria-ingestion-env"
    assert "ingestion" in EnvironmentManager.TIMEOUTS


# ── Script base IPC contract ─────────────────────────────────────────────────

def test_base_json_serializer_handles_numpy():
    np = pytest.importorskip("numpy")
    from aria.scripts._base import _json_serializer
    assert _json_serializer(np.int64(42)) == 42
    assert _json_serializer(np.float64(3.14)) == 3.14
    assert _json_serializer(np.array([1, 2, 3])) == [1, 2, 3]


def test_ipc_cycle_roundtrip():
    from aria.scripts._base import run_script  # noqa: F401 (import contract)
    with tempfile.TemporaryDirectory() as tmp:
        input_file = Path(tmp) / "input.json"
        output_file = Path(tmp) / "output.json"
        input_file.write_text(json.dumps({"test_key": "test_value", "number": 42}))
        loaded = json.loads(input_file.read_text())
        output_file.write_text(json.dumps(
            {"status": "success", "echo": loaded["test_key"],
             "doubled": loaded["number"] * 2}
        ))
        final = json.loads(output_file.read_text())
    assert final == {"status": "success", "echo": "test_value", "doubled": 84}


# ── Script imports (lazy scanpy — import is light) ───────────────────────────

def test_rna_qc_importable():
    from aria.scripts.rna_qc import rna_qc  # noqa: F401


def test_rna_clustering_importable():
    from aria.scripts.rna_clustering import rna_clustering  # noqa: F401


def test_rna_qc_invalid_path_is_structured():
    pytest.importorskip("scanpy")
    from aria.scripts.rna_qc import rna_qc
    result = rna_qc({"data_path": "/nonexistent/path.h5ad"})
    assert "status" in result          # structured error, not an exception


# ── Workspace ─────────────────────────────────────────────────────────────────

def test_workspace_exists():
    from aria.utils.environment_manager import env_manager
    assert env_manager.workspace.exists()
    assert env_manager.workspace.is_dir()


# ── live PBMC 3k QC (scanpy + local dataset) ─────────────────────────────────

def test_rna_qc_on_pbmc3k():
    sc = pytest.importorskip("scanpy")
    pbmc = None
    for c in (Path.home() / "aria-data" / "pbmc3k_test",
              Path.home() / "aria-data" / "pbmc3k_test" / "hg19"):
        if c.exists() and list(c.rglob("*.mtx*")):
            pbmc = c
            break
    if pbmc is None:
        pytest.skip("PBMC 3k dataset not present (run install.sh to download)")
    adata = sc.read_10x_mtx(str(pbmc), var_names="gene_symbols", cache=True)
    assert adata.n_obs > 0
    from aria.scripts.rna_qc import rna_qc
    result = rna_qc({"data_path": str(pbmc), "organism": "Homo sapiens"})
    assert result["status"] == "success"
    assert result["n_cells_after"] > 2000
    qc_out = Path(result.get("output_path", ""))
    if qc_out.exists():
        qc_out.unlink()
