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


def test_bulk_rna_production_contract_requires_metadata_and_contrasts(tmp_path):
    from aria.utils.script_contracts import contract_for_script

    counts_path = tmp_path / "counts.tsv"
    counts_path.write_text("gene\tA_1\tB_1\nG1\t10\t20\n", encoding="utf-8")
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "sample\tcondition\nA_1\tA\nB_1\tB\n",
        encoding="utf-8",
    )

    contract = contract_for_script("aria/scripts/rna_bulk_de.py")
    assert contract is not None

    issues = contract.validate_params({
        "files": [str(counts_path)],
        "design_factor": "condition",
    })
    assert {issue.field for issue in issues} == {"metadata_file", "contrasts"}

    issues = contract.validate_params({
        "files": [str(counts_path)],
        "metadata_file": str(metadata_path),
        "design_factor": "condition",
        "contrasts": [],
    })
    assert any(issue.field == "contrasts" for issue in issues)

    issues = contract.validate_params({
        "files": [str(counts_path)],
        "metadata_file": str(metadata_path),
        "design_factor": "condition",
        "contrasts": [{"numerator": "B", "denominator": "A"}],
    })
    assert issues == []


# T10 (tri-audit 2026-06-14): the script subprocess is launched with
# subprocess.Popen + start_new_session (ADR-020/R5 process-group reaping), not
# subprocess.run. These tests mocked `subprocess.run`, so the mock no longer
# intercepted, the real `conda run` executed, and they failed with SubprocessFailed
# — silently dropping coverage of the output-contract path. The fake now patches
# Popen with the same cmd[-1] → output_file convention the manager uses.

def _fake_popen_factory(payload: dict):
    class _FakePopen:
        def __init__(self, cmd, *args, **kwargs):
            self._cmd = list(cmd)
            self.returncode = 0
            # run_in_stack also enumerates envs via `Popen(["conda","env","list",
            # "--json"])`; ONLY the script execution (`conda run ... <output.json>`)
            # writes the result. Gating on "run" stops the enumeration call from
            # creating a stray ./--json (cmd[-1] == "--json") in the repo root.
            if "run" in self._cmd:
                Path(self._cmd[-1]).write_text(json.dumps(payload))

        def communicate(self, timeout=None):
            return ("", "")

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    return _FakePopen


def test_environment_manager_validates_script_output_contract(
    tmp_path, monkeypatch
):
    from aria.utils import environment_manager as em

    mgr = _manager(tmp_path)
    monkeypatch.setattr(mgr, "_resolve_env", lambda stack: "aria-rna-env")
    data_path = tmp_path / "input.h5ad"
    data_path.write_text("placeholder")

    monkeypatch.setattr(em.subprocess, "Popen", _fake_popen_factory({
        "status": "success",
        "n_cells_before": 10,
    }))

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
    monkeypatch.setattr(mgr, "_resolve_env", lambda stack: "aria-rna-env")
    data_path = tmp_path / "input.h5ad"
    data_path.write_text("placeholder")

    monkeypatch.setattr(em.subprocess, "Popen", _fake_popen_factory({
        "status": "success",
        "n_cells_before": 10,
        "n_cells_after": 9,
    }))

    result = mgr.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_qc.py",
        params={"data_path": str(data_path)},
    )

    assert result["status"] == "success"
    assert result["ipc_contract"]["script_path"] == "aria/scripts/rna_qc.py"
    assert result["ipc_contract"]["contract_version"] == "1.0"
