"""A5: one portable, fail-closed scientific environment contract."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from aria.utils.environment_specs import (
    ENVIRONMENT_SPECS,
    lock_environment_names,
    portable_pip_lock_lines,
    setup_environment_definitions,
)
from aria.utils.environment_audit import explicit_artifacts, pip_pins, required_binaries


ROOT = Path(__file__).resolve().parents[1]
ENVS = ROOT / "envs"


def test_one_registry_drives_every_runtime_stack_and_timeout():
    from aria.utils.environment_manager import EnvironmentManager

    assert set(EnvironmentManager.STACKS) == set(ENVIRONMENT_SPECS)
    assert EnvironmentManager.STACKS == {
        stack: spec.env_name for stack, spec in ENVIRONMENT_SPECS.items()
    }
    assert EnvironmentManager.TIMEOUTS == {
        stack: spec.timeout_s for stack, spec in ENVIRONMENT_SPECS.items()
    }


def test_setup_reads_versioned_yamls_instead_of_inline_copies():
    from aria.agents.setup_agent import ARIA_ENVS

    assert ARIA_ENVS == setup_environment_definitions()
    for env_name, definition in ARIA_ENVS.items():
        yml_path = Path(definition["yml_path"])
        assert yml_path.exists()
        assert yml_path.read_text(encoding="utf-8").startswith(f"name: {env_name}\n")
        assert "yml" not in definition


def test_active_scientific_and_raw_stacks_have_committed_locks():
    required = set(lock_environment_names())
    assert {
        "aria-rna-env",
        "aria-ingestion-env",
        "aria-rnaseq-env",
        "aria-atacseq-env",
        "aria-chromatin-env",
        "aria-tobias-env",
        "aria-bench-env",
    } <= required
    for env_name in required:
        lock = ENVS / f"{env_name}.linux-64.lock"
        assert lock.exists(), f"missing explicit lock for {env_name}"
        text = lock.read_text(encoding="utf-8")
        assert "@EXPLICIT" in text and "https://" in text


def test_all_pip_locks_are_portable_exact_version_pins():
    for lock in ENVS.glob("*.pip.lock"):
        lines = [
            line.strip() for line in lock.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert lines, f"empty pip lock should not be committed: {lock.name}"
        assert all("==" in line or (" @ git+https://" in line and ".git@" in line) for line in lines), lock.name
        assert not any(" @ file:" in line or "file:///" in line for line in lines)


def test_portable_pip_lock_uses_only_pypi_records_and_canonical_names():
    records = [
        {"name": "NumPy", "version": "2.4.3", "channel": "conda-forge"},
        {"name": "My_Package", "version": "1.2.0", "channel": "pypi"},
        {"name": "other.package", "version": "3.0", "channel": "pypi"},
    ]
    assert portable_pip_lock_lines(records) == [
        "my-package==1.2.0",
        "other-package==3.0",
    ]
    assert portable_pip_lock_lines(records, {
        "my-package": {
            "url": "https://example.test/my_package.git",
            "vcs_info": {"vcs": "git", "commit_id": "abc123"},
        }
    }) == [
        "my-package @ git+https://example.test/my_package.git@abc123",
        "other-package==3.0",
    ]
    with pytest.raises(ValueError, match="non-portable direct install"):
        portable_pip_lock_lines(records, {
            "my-package": {"url": "file:///tmp/build/my_package"},
        })


def test_exact_auditor_normalizes_lock_metadata_and_owns_binary_contracts():
    assert explicit_artifacts("# platform: linux-64\n@EXPLICIT\nhttps://x/pkg.conda\n") == [
        "https://x/pkg.conda"
    ]
    assert pip_pins("# generated\nMy_Package==1.2\n") == ["my-package==1.2"]
    assert {"fastp", "STAR", "featureCounts", "samtools"} <= set(
        required_binaries("aria-rnaseq-env")
    )
    assert required_binaries("aria-hic-env") == []


def test_missing_scientific_environment_fails_closed(monkeypatch):
    from aria.utils.environment_manager import (
        EnvironmentManager,
        MissingEnvironment,
    )

    manager = EnvironmentManager.__new__(EnvironmentManager)
    monkeypatch.setattr(
        manager,
        "check_environments",
        lambda: {stack: False for stack in EnvironmentManager.STACKS},
    )

    with pytest.raises(MissingEnvironment, match="aria-rna-env"):
        manager._resolve_env("rna")
    assert not hasattr(EnvironmentManager, "FALLBACK_ENV")


def test_run_in_stack_reports_missing_environment_without_launching(tmp_path, monkeypatch):
    from aria.utils import environment_manager as module

    script = tmp_path / "noop.py"
    script.write_text("# never executed\n", encoding="utf-8")
    manager = module.EnvironmentManager.__new__(module.EnvironmentManager)
    manager.workspace = tmp_path / "workspace"
    manager.workspace.mkdir()
    (manager.workspace / "failed").mkdir()
    monkeypatch.setattr(
        manager,
        "_resolve_env",
        lambda stack: (_ for _ in ()).throw(
            module.MissingEnvironment(stack, module.EnvironmentManager.STACKS[stack])
        ),
    )
    monkeypatch.setattr(
        module.subprocess,
        "Popen",
        lambda *args, **kwargs: pytest.fail("missing env must not launch subprocess"),
    )

    result = manager.run_in_stack("rna", str(script), {}, timeout=120)

    assert result["status"] == "error"
    assert result["error_type"] == "MissingEnvironment"
    assert result["stack"] == "rna"
    assert result["environment"] == "aria-rna-env"


def test_lock_generator_reads_registry_not_a_second_whitelist():
    source = (ROOT / "scripts" / "generate_locks.sh").read_text(encoding="utf-8")
    assert "lock-envs" in source
    assert "WHITELIST=" not in source
    assert "pip freeze" not in source
    assert "conda list --name" in source and "--json" in source


def test_resolver_source_contains_no_fallback_or_alias_path():
    from aria.utils.environment_manager import EnvironmentManager

    source = inspect.getsource(EnvironmentManager._resolve_env)
    assert "fallback" not in source.lower()
    assert "env_aliases" not in source
