"""One-command lock-first installer for ARIA's isolated scientific environments.

Exposes ``aria setup`` so users never have to type per-environment
``conda create --file ...`` + ``pip install --no-deps -r ...`` commands by hand.

Design principles:
  * Single source of truth: every environment resolves from
    ``aria.utils.environment_specs.ENVIRONMENT_SPECS``.
  * Minimal default footprint: ``aria setup`` (or ``aria setup rna``) installs
    ONLY ``aria-rna-env``, which covers Bulk RNA count matrices, scRNA-seq,
    pseudobulk DE, ORA/GSEA, LIANA, and trajectory. Additional modality stacks
    (``raw-rna``, ``chromatin``, ``tobias``, ``multiome``, ``hic``, ``all``)
    are installed only when requested or on demand by ``SetupAgent``.
  * Lock-first reproducibility on Linux ``x86_64``: uses committed ``@EXPLICIT``
    ``envs/<env>.linux-64.lock`` + ``envs/<env>.pip.lock`` (fast, zero SAT
    solving, bit-reproducible) and falls back cleanly to ``envs/<env>.yml``
    when requested or on non-Linux architectures.
  * Air-gap honest: refuses network installation fail-closed when
    ``ARIA_AIR_GAPPED=1`` (or ``air_gapped=True``).
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from aria.utils.environment_specs import (
    ENVS_DIR,
    ENVIRONMENT_SPECS,
    EnvironmentSpec,
)

# User-facing modality profiles -> ordered tuple of stack keys in ENVIRONMENT_SPECS.
SETUP_PROFILES: dict[str, dict[str, Any]] = {
    "rna": {
        "stacks": ("rna",),
        "description": "Bulk RNA matrices, scRNA-seq, pseudobulk DE, ORA/GSEA, LIANA, trajectory (1 env — recommended default)",
    },
    "raw-rna": {
        "stacks": ("rna", "rnaseq", "ingestion"),
        "description": "Bulk RNA & scRNA raw FASTQ pipelines (STAR, featureCounts, kb-python) + RNA analysis",
    },
    "chromatin": {
        "stacks": ("rna", "chromatin", "atacseq"),
        "description": "scATAC-seq, Bulk ATAC-seq, ChIP-seq, CUT&RUN, CUT&TAG (bwa-mem2, chromap, MACS3, SEACR, SnapATAC2)",
    },
    "tobias": {
        "stacks": ("tobias",),
        "description": "Optional TOBIAS Tn5-bias-corrected differential TF footprinting",
    },
    "multiome": {
        "stacks": ("rna", "chromatin", "integration"),
        "description": "Paired 10x Multiome (WNN, MOFA+, peak2gene, eRegulons) + Cross-Modal concordance",
    },
    "hic": {
        "stacks": ("hic", "atacseq"),
        "description": "3D genome architecture (Hi-C / Micro-C: cooler, cooltools, pairtools, HiCRep)",
    },
    "all": {
        "stacks": ("rna", "rnaseq", "ingestion", "atacseq", "chromatin", "tobias", "integration", "hic"),
        "description": "All 8 user-facing scientific stacks (excludes developer-only R benchmark comparator env)",
    },
}


def find_conda_executable() -> str | None:
    """Locate ``mamba`` or ``conda`` in PATH or standard Miniforge/Anaconda locations."""
    for cmd in ("mamba", "conda"):
        found = shutil.which(cmd)
        if found:
            return found
    for candidate in (
        Path.home() / "miniforge3" / "bin" / "mamba",
        Path.home() / "miniforge3" / "bin" / "conda",
        Path.home() / "mambaforge" / "bin" / "mamba",
        Path.home() / "mambaforge" / "bin" / "conda",
        Path.home() / "anaconda3" / "bin" / "conda",
        Path.home() / "miniconda3" / "bin" / "conda",
    ):
        if candidate.is_file():
            return str(candidate)
    return None


def list_installed_envs(conda_cmd: str | None) -> set[str]:
    """Return the set of installed Conda environment basenames."""
    if not conda_cmd:
        return set()
    try:
        proc = subprocess.run(
            [conda_cmd, "env", "list", "--json"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if proc.returncode != 0:
            return set()
        data = json.loads(proc.stdout)
        return {Path(p).name for p in data.get("envs", [])}
    except Exception:
        return set()


def _spec_for_env_name(env_name: str) -> EnvironmentSpec | None:
    for spec in ENVIRONMENT_SPECS.values():
        if spec.env_name == env_name:
            return spec
    return None


def _has_noncomment_requirements(path: Path) -> bool:
    if not path.is_file():
        return False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return True
    return False


def _is_linux_x86_64() -> bool:
    return (
        platform.system().lower() == "linux"
        and platform.machine().lower() in {"x86_64", "amd64"}
    )


def resolve_setup_targets(
    targets: list[str] | None,
    *,
    include_bench: bool = False,
) -> list[EnvironmentSpec]:
    """Resolve profile names, stack keys, or env names into unique EnvironmentSpecs."""
    raw = list(targets) if targets else ["rna"]
    ordered: list[EnvironmentSpec] = []
    seen_envs: set[str] = set()

    def _add_spec(spec: EnvironmentSpec) -> None:
        if spec.env_name not in seen_envs:
            seen_envs.add(spec.env_name)
            ordered.append(spec)

    for token in raw:
        key = token.strip().lower()
        if not key:
            continue
        if key in SETUP_PROFILES:
            for stack in SETUP_PROFILES[key]["stacks"]:
                _add_spec(ENVIRONMENT_SPECS[stack])
            if key == "all" and include_bench:
                _add_spec(ENVIRONMENT_SPECS["benchmark"])
            continue
        if key in ENVIRONMENT_SPECS:
            _add_spec(ENVIRONMENT_SPECS[key])
            continue
        matched = _spec_for_env_name(token.strip())
        if matched is not None:
            _add_spec(matched)
            continue
        valid = ", ".join(sorted(SETUP_PROFILES))
        raise ValueError(
            f"Unknown setup target {token!r}. Valid profiles: {valid} "
            "(or any stack/environment name from `aria setup --list`)."
        )

    if include_bench:
        _add_spec(ENVIRONMENT_SPECS["benchmark"])
    return ordered


def build_install_plan(
    env_name: str,
    conda_cmd: str,
    *,
    envs_dir: Path | None = None,
    use_locks: bool = True,
    force: bool = False,
    is_linux_x64: bool | None = None,
) -> dict[str, Any]:
    """Build the ordered command plan to install one ARIA scientific environment.

    On Linux ``x86_64``, prefers the committed ``@EXPLICIT`` ``.linux-64.lock``
    plus ``.pip.lock`` (fast, deterministic, no SAT solver). Falls back cleanly
    to ``envs/<env>.yml`` when ``use_locks=False``, on non-Linux hosts, or when
    no lockfile is present.
    """
    spec = _spec_for_env_name(env_name)
    if spec is None:
        raise ValueError(f"Unknown ARIA environment: {env_name!r}")

    base_dir = Path(envs_dir) if envs_dir is not None else ENVS_DIR
    conda_lock = base_dir / f"{env_name}.linux-64.lock"
    pip_lock = base_dir / f"{env_name}.pip.lock"
    yml_path = base_dir / spec.yml_name

    linux_ok = _is_linux_x86_64() if is_linux_x64 is None else bool(is_linux_x64)
    commands: list[list[str]] = []

    if force:
        commands.append([conda_cmd, "env", "remove", "-n", env_name, "-y"])

    if use_locks and linux_ok and conda_lock.is_file():
        method = "explicit_lock"
        commands.append(
            [conda_cmd, "create", "--name", env_name, "--file", str(conda_lock), "-y"]
        )
        if _has_noncomment_requirements(pip_lock):
            commands.append([
                conda_cmd, "run", "--name", env_name,
                "python", "-m", "pip", "install", "--no-deps",
                "--requirement", str(pip_lock),
            ])
        source_files = [str(conda_lock)] + (
            [str(pip_lock)] if _has_noncomment_requirements(pip_lock) else []
        )
    elif yml_path.is_file():
        method = "environment_yml"
        commands.append(
            [conda_cmd, "env", "create", "-f", str(yml_path), "--name", env_name, "-y"]
        )
        source_files = [str(yml_path)]
    else:
        raise FileNotFoundError(
            f"No lockfile or YAML spec found for {env_name} under {base_dir}"
        )

    return {
        "env_name": env_name,
        "stack": spec.stack,
        "description": spec.description,
        "method": method,
        "source_files": source_files,
        "commands": commands,
    }


def install_environment(
    env_name: str,
    conda_cmd: str,
    *,
    envs_dir: Path | None = None,
    use_locks: bool = True,
    force: bool = False,
    dry_run: bool = False,
    air_gapped: bool | None = None,
    timeout: int = 1800,
    installed_envs: set[str] | None = None,
) -> dict[str, Any]:
    """Install (or verify already installed) a single ARIA Conda environment."""
    is_air_gapped = (
        os.environ.get("ARIA_AIR_GAPPED", "").strip().lower() in {"1", "true", "yes"}
        if air_gapped is None
        else bool(air_gapped)
    )
    installed = (
        list_installed_envs(conda_cmd)
        if installed_envs is None
        else set(installed_envs)
    )
    if env_name in installed and not force:
        return {
            "env_name": env_name,
            "status": "already_installed",
            "method": "existing",
            "commands": [],
        }

    plan = build_install_plan(
        env_name,
        conda_cmd,
        envs_dir=envs_dir,
        use_locks=use_locks,
        force=force,
    )

    if dry_run:
        return {
            "env_name": env_name,
            "status": "dry_run",
            "method": plan["method"],
            "commands": plan["commands"],
            "source_files": plan["source_files"],
        }

    if is_air_gapped:
        return {
            "env_name": env_name,
            "status": "blocked_air_gapped",
            "method": plan["method"],
            "commands": plan["commands"],
            "error": (
                f"ARIA_AIR_GAPPED is enabled; refusing network installation of "
                f"'{env_name}'. Provision it offline from its lockfile before launch."
            ),
        }

    for cmd in plan["commands"]:
        # Ignore removal error if --force was passed on an absent environment.
        is_remove = len(cmd) >= 3 and cmd[1] == "env" and cmd[2] == "remove"
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return {
                "env_name": env_name,
                "status": "error",
                "method": plan["method"],
                "commands": plan["commands"],
                "error": f"Timeout installing {env_name} (>{timeout}s)",
            }
        except Exception as exc:
            return {
                "env_name": env_name,
                "status": "error",
                "method": plan["method"],
                "commands": plan["commands"],
                "error": str(exc),
            }
        if proc.returncode != 0 and not is_remove:
            stderr_tail = (proc.stderr or proc.stdout or "").strip()[-400:]
            return {
                "env_name": env_name,
                "status": "error",
                "method": plan["method"],
                "commands": plan["commands"],
                "error": f"Command failed (rc={proc.returncode}): {stderr_tail}",
            }

    return {
        "env_name": env_name,
        "status": "installed",
        "method": plan["method"],
        "commands": plan["commands"],
        "source_files": plan["source_files"],
    }


def _format_status_table(installed: set[str]) -> str:
    lines = [
        "ARIA Scientific Environment Profiles (`aria setup <profile>`):",
        "",
        "  Profiles:",
    ]
    for name, info in SETUP_PROFILES.items():
        env_names = [ENVIRONMENT_SPECS[s].env_name for s in info["stacks"]]
        ready = sum(1 for e in env_names if e in installed)
        lines.append(
            f"    {name:<11} [{ready}/{len(env_names)} installed]  {info['description']}"
        )
    lines.extend(["", "  Individual environments:"])
    seen: set[str] = set()
    for spec in ENVIRONMENT_SPECS.values():
        if spec.env_name in seen:
            continue
        seen.add(spec.env_name)
        state = "installed" if spec.env_name in installed else "not installed"
        lock_ok = (ENVS_DIR / f"{spec.env_name}.linux-64.lock").is_file()
        lock_tag = "lock" if lock_ok else "yml"
        lines.append(
            f"    {spec.env_name:<22} ({spec.stack:<11}) [{state:<13} | {lock_tag}]  {spec.description}"
        )
    lines.extend([
        "",
        "  Quick start (covers Bulk RNA matrices + scRNA-seq with 1 environment):",
        "    aria setup rna",
        "",
        "  Install additional modalities only when you need them:",
        "    aria setup chromatin | multiome | hic | raw-rna | all",
    ])
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="aria setup",
        description=(
            "Install ARIA's lock-pinned scientific Conda environments on demand. "
            "Defaults to `rna` (1 environment: Bulk RNA matrices + scRNA-seq)."
        ),
    )
    parser.add_argument(
        "targets",
        nargs="*",
        help=(
            "Profile(s) or environment(s) to install: "
            "rna (default), raw-rna, chromatin, tobias, multiome, hic, all, "
            "or any specific env name (e.g. aria-rna-env)."
        ),
    )
    parser.add_argument(
        "--list",
        "--status",
        dest="list_only",
        action="store_true",
        help="Show available profiles, individual environments, and installation status.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the exact lockfile/YAML commands that would run without executing them.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recreate target environments even if they already exist.",
    )
    parser.add_argument(
        "--from-yml",
        action="store_true",
        help="Use envs/*.yml instead of @EXPLICIT linux-64 lockfiles.",
    )
    parser.add_argument(
        "--include-bench",
        action="store_true",
        help="Also include the developer R comparator environment (aria-bench-env).",
    )
    parser.add_argument(
        "--audit",
        action="store_true",
        help="Run `aria.utils.environment_audit` on each target after installation.",
    )
    args = parser.parse_args(argv)

    conda_cmd = find_conda_executable()
    installed = list_installed_envs(conda_cmd) if conda_cmd else set()

    if args.list_only:
        print(_format_status_table(installed))
        return 0

    try:
        specs = resolve_setup_targets(args.targets, include_bench=args.include_bench)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    effective_conda = conda_cmd or ("conda" if args.dry_run else None)
    if effective_conda is None:
        print(
            "Error: Neither `mamba` nor `conda` was found on PATH.\n"
            "Install Miniforge first: https://github.com/conda-forge/miniforge/releases/latest",
            file=sys.stderr,
        )
        return 1

    any_failed = False
    for spec in specs:
        res = install_environment(
            spec.env_name,
            effective_conda,
            use_locks=(not args.from_yml),
            force=args.force,
            dry_run=args.dry_run,
            installed_envs=installed,
        )
        status = res["status"]
        if status == "already_installed":
            print(f"[OK] {spec.env_name} ({spec.stack}): already installed.")
        elif status == "dry_run":
            print(f"[DRY-RUN] {spec.env_name} ({spec.stack}) via {res['method']}:")
            for cmd in res["commands"]:
                print("  " + " ".join(cmd))
        elif status == "installed":
            print(f"[INSTALLED] {spec.env_name} ({spec.stack}) via {res['method']}.")
            installed.add(spec.env_name)
        else:
            any_failed = True
            print(
                f"[FAILED] {spec.env_name} ({spec.stack}): {res.get('error', status)}",
                file=sys.stderr,
            )

    if args.audit and not args.dry_run and not any_failed:
        for spec in specs:
            proc = subprocess.run(
                [sys.executable, "-m", "aria.utils.environment_audit", "--env", spec.env_name],
                check=False,
            )
            if proc.returncode != 0:
                any_failed = True

    return 1 if any_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
