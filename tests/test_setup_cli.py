"""Tests for the lock-first `aria setup` CLI and shared environment installer."""

from __future__ import annotations

from pathlib import Path

import pytest

from aria.setup_cli import (
    SETUP_PROFILES,
    build_install_plan,
    install_environment,
    main,
    resolve_setup_targets,
)


def test_default_setup_target_is_single_rna_env():
    specs = resolve_setup_targets([])
    assert [s.env_name for s in specs] == ["aria-rna-env"]


def test_setup_profiles_resolve_to_valid_unique_envs():
    for profile in SETUP_PROFILES:
        specs = resolve_setup_targets([profile])
        names = [s.env_name for s in specs]
        assert len(names) == len(set(names))
        assert len(names) >= 1

    all_names = [s.env_name for s in resolve_setup_targets(["all"])]
    assert "aria-rna-env" in all_names
    assert "aria-hic-env" in all_names
    assert "aria-bench-env" not in all_names

    all_with_bench = [
        s.env_name for s in resolve_setup_targets(["all"], include_bench=True)
    ]
    assert "aria-bench-env" in all_with_bench


def test_unknown_setup_target_raises_clear_error():
    with pytest.raises(ValueError, match="Unknown setup target"):
        resolve_setup_targets(["nonexistent_profile"])


def test_build_install_plan_prefers_explicit_lock_on_linux_x64():
    plan = build_install_plan("aria-rna-env", "conda", use_locks=True, is_linux_x64=True)
    assert plan["method"] == "explicit_lock"
    assert len(plan["commands"]) == 2
    assert plan["commands"][0][:4] == ["conda", "create", "--name", "aria-rna-env"]
    assert plan["commands"][0][4] == "--file"
    assert plan["commands"][0][5].endswith("aria-rna-env.linux-64.lock")
    assert plan["commands"][1][:6] == [
        "conda", "run", "--name", "aria-rna-env", "python", "-m",
    ]
    assert "--no-deps" in plan["commands"][1]
    assert plan["commands"][1][-1].endswith("aria-rna-env.pip.lock")


def test_build_install_plan_skips_empty_comment_only_pip_lock(tmp_path):
    (tmp_path / "aria-rnaseq-env.linux-64.lock").write_text("@EXPLICIT\n", encoding="utf-8")
    (tmp_path / "aria-rnaseq-env.pip.lock").write_text("# no pypi deps\n", encoding="utf-8")

    plan = build_install_plan(
        "aria-rnaseq-env", "conda", envs_dir=tmp_path, use_locks=True, is_linux_x64=True
    )
    assert plan["method"] == "explicit_lock"
    assert len(plan["commands"]) == 1
    assert plan["commands"][0][-2].endswith("aria-rnaseq-env.linux-64.lock")


def test_build_install_plan_falls_back_to_yml_when_requested_or_non_linux():
    plan_yml = build_install_plan(
        "aria-rna-env", "mamba", use_locks=False, is_linux_x64=True
    )
    assert plan_yml["method"] == "environment_yml"
    assert plan_yml["commands"] == [[
        "mamba", "env", "create", "-f", plan_yml["source_files"][0],
        "--name", "aria-rna-env", "-y",
    ]]

    plan_mac = build_install_plan(
        "aria-rna-env", "conda", use_locks=True, is_linux_x64=False
    )
    assert plan_mac["method"] == "environment_yml"


def test_install_environment_respects_air_gapped_mode():
    res = install_environment(
        "aria-rna-env",
        "conda",
        air_gapped=True,
        installed_envs=set(),
    )
    assert res["status"] == "blocked_air_gapped"
    assert "ARIA_AIR_GAPPED" in res["error"]


def test_install_environment_skips_existing_unless_forced():
    res = install_environment(
        "aria-rna-env",
        "conda",
        installed_envs={"aria-rna-env"},
    )
    assert res["status"] == "already_installed"
    assert res["commands"] == []

    res_force = install_environment(
        "aria-rna-env",
        "conda",
        force=True,
        dry_run=True,
        installed_envs={"aria-rna-env"},
    )
    assert res_force["status"] == "dry_run"
    assert res_force["commands"][0] == ["conda", "env", "remove", "-n", "aria-rna-env", "-y"]


def test_cli_list_and_dry_run(capsys, monkeypatch):
    monkeypatch.setattr("aria.setup_cli.find_conda_executable", lambda: "conda")
    monkeypatch.setattr("aria.setup_cli.list_installed_envs", lambda _c: set())

    rc_list = main(["--list"])
    assert rc_list == 0
    out_list = capsys.readouterr().out
    assert "aria setup rna" in out_list
    assert "aria-rna-env" in out_list

    rc_dry = main(["rna", "--dry-run"])
    assert rc_dry == 0
    out_dry = capsys.readouterr().out
    assert "[DRY-RUN] aria-rna-env" in out_dry
    assert "aria-rna-env.linux-64.lock" in out_dry
