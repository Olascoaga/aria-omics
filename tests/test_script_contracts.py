import json
from pathlib import Path


def _manager(tmp_path):
    from aria.utils.environment_manager import EnvironmentManager

    mgr = EnvironmentManager.__new__(EnvironmentManager)
    mgr.workspace = tmp_path / "workspace"
    mgr.workspace.mkdir()
    (mgr.workspace / "failed").mkdir()
    mgr._conda_ok = True
    mgr._env_cache = {
        "rna": True,
        "rnaseq": False,
        "chromatin": False,
        "hic": False,
        "integration": False,
        "spatial": True,
    }
    mgr._env_cache_key = 0.0
    mgr.MAX_FAILED_RUNS = 20
    return mgr


def test_script_contract_rejects_missing_required_input(tmp_path):
    mgr = _manager(tmp_path)

    result = mgr.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_qc.py",
        params={},
    )

    assert result["status"] == "error"
    assert result["error_type"] == "InvalidScriptParams"
    assert result["contract_stage"] == "input"
    assert result["contract_issues"][0]["field"] == "data_path"


def test_script_contract_rejects_version_mismatch(tmp_path):
    mgr = _manager(tmp_path)
    data_path = tmp_path / "input.h5ad"
    data_path.write_text("placeholder")

    result = mgr.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_qc.py",
        params={
            "data_path": str(data_path),
            "_aria_contract_version": "999.0",
        },
    )

    assert result["status"] == "error"
    assert result["error_type"] == "IncompatibleScriptContract"
    assert "expected 1.0" in result["details"]


def test_environment_manager_validates_script_output_contract(
    tmp_path, monkeypatch
):
    from aria.utils import environment_manager as em

    mgr = _manager(tmp_path)
    data_path = tmp_path / "input.h5ad"
    data_path.write_text("placeholder")

    class _Result:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(cmd, *args, **kwargs):
        output_file = Path(cmd[-1])
        output_file.write_text(json.dumps({
            "status": "success",
            "n_cells_before": 10,
        }))
        return _Result()

    monkeypatch.setattr(em.subprocess, "run", fake_run)

    result = mgr.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_qc.py",
        params={"data_path": str(data_path)},
    )

    assert result["status"] == "error"
    assert result["error_type"] == "IncompatibleScriptContract"
    assert result["contract_stage"] == "output"
    assert any(
        issue["field"] == "n_cells_after"
        for issue in result["contract_issues"]
    )


def test_environment_manager_attaches_contract_metadata_on_success(
    tmp_path, monkeypatch
):
    from aria.utils import environment_manager as em

    mgr = _manager(tmp_path)
    data_path = tmp_path / "input.h5ad"
    data_path.write_text("placeholder")

    class _Result:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(cmd, *args, **kwargs):
        output_file = Path(cmd[-1])
        output_file.write_text(json.dumps({
            "status": "success",
            "n_cells_before": 10,
            "n_cells_after": 9,
        }))
        return _Result()

    monkeypatch.setattr(em.subprocess, "run", fake_run)

    result = mgr.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_qc.py",
        params={"data_path": str(data_path)},
    )

    assert result["status"] == "success"
    assert result["ipc_contract"]["script_path"] == "aria/scripts/rna_qc.py"
    assert result["ipc_contract"]["contract_version"] == "1.0"
