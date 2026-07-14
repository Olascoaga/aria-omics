"""Audit ARIA scientific environments against their committed exact locks.

The auditor supports both installed named environments and clean prefixes created
directly from the explicit Conda locks.  It compares the complete Conda artifact
set, the complete portable PyPI pin set, and the required command-line binaries.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable

from aria.utils.environment_specs import (
    ENVIRONMENT_SPECS,
    ENVS_DIR,
    lock_environment_names,
    portable_pip_lock_lines,
)


def explicit_artifacts(text: str) -> list[str]:
    """Return the canonical artifact URLs from a Conda explicit lock/output."""
    return sorted(
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith(("#", "@"))
    )


def pip_pins(text: str) -> list[str]:
    """Return normalized non-comment pins from a portable pip lock."""
    pins = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        separator = " @ " if " @ " in line else "=="
        name, value = line.split(separator, 1)
        name = re.sub(r"[-_.]+", "-", name.strip()).lower()
        pins.append(f"{name}{separator}{value.strip()}")
    return sorted(pins, key=str.lower)


def required_binaries(env_name: str) -> list[str]:
    """Union binary requirements for every stack routed to ``env_name``."""
    return sorted({
        binary
        for spec in ENVIRONMENT_SPECS.values()
        if spec.env_name == env_name and spec.lock_required
        for binary in spec.required_binaries
    })


def _conda(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["conda", *command],
        check=False,
        capture_output=True,
        text=True,
    )


def _direct_urls(selector: list[str]) -> dict[str, dict[str, Any]]:
    script = (
        "import importlib.metadata as m,json,re; out={}; "
        "[(out.__setitem__(re.sub(r'[-_.]+','-',d.metadata['Name']).lower(),json.loads(x))) "
        "for d in m.distributions() if d.metadata.get('Name') "
        "and (x:=d.read_text('direct_url.json'))]; print(json.dumps(out))"
    )
    probe = _conda(["run", *selector, "python", "-c", script])
    return json.loads(probe.stdout) if probe.returncode == 0 else {}


def _selector(*, env_name: str | None = None, prefix: Path | None = None) -> list[str]:
    if (env_name is None) == (prefix is None):
        raise ValueError("provide exactly one of env_name or prefix")
    return ["--name", env_name] if env_name is not None else ["--prefix", str(prefix)]


def audit_environment(
    lock_env_name: str,
    *,
    env_name: str | None = None,
    prefix: Path | None = None,
) -> dict[str, Any]:
    """Compare one live environment with its lock and required binaries."""
    selector = _selector(env_name=env_name, prefix=prefix)
    conda_lock = ENVS_DIR / f"{lock_env_name}.linux-64.lock"
    pip_lock = ENVS_DIR / f"{lock_env_name}.pip.lock"
    result: dict[str, Any] = {
        "environment": lock_env_name,
        "target": env_name if env_name is not None else str(prefix),
        "conda_lock": str(conda_lock),
        "pip_lock": str(pip_lock) if pip_lock.exists() else None,
    }

    explicit = _conda(["list", *selector, "--explicit"])
    records = _conda(["list", *selector, "--json"])
    if explicit.returncode or records.returncode:
        result.update({
            "status": "error",
            "reason": "environment_unavailable",
            "details": (explicit.stderr or records.stderr).strip(),
        })
        return result

    expected_conda = explicit_artifacts(conda_lock.read_text(encoding="utf-8"))
    actual_conda = explicit_artifacts(explicit.stdout)
    expected_pip = pip_pins(pip_lock.read_text(encoding="utf-8")) if pip_lock.exists() else []
    actual_pip = portable_pip_lock_lines(json.loads(records.stdout), _direct_urls(selector))
    expected_pip_normalized = expected_pip

    binaries: dict[str, bool] = {}
    for binary in required_binaries(lock_env_name):
        probe = _conda(["run", *selector, "which", binary])
        binaries[binary] = probe.returncode == 0 and bool(probe.stdout.strip())

    conda_match = expected_conda == actual_conda
    pip_match = expected_pip_normalized == actual_pip
    binaries_match = all(binaries.values())
    result.update({
        "status": "pass" if conda_match and pip_match and binaries_match else "fail",
        "conda_exact": conda_match,
        "conda_expected": len(expected_conda),
        "conda_actual": len(actual_conda),
        "conda_missing": sorted(set(expected_conda) - set(actual_conda)),
        "conda_extra": sorted(set(actual_conda) - set(expected_conda)),
        "pip_exact": pip_match,
        "pip_expected": len(expected_pip_normalized),
        "pip_actual": len(actual_pip),
        "pip_missing": sorted(set(expected_pip_normalized) - set(actual_pip)),
        "pip_extra": sorted(set(actual_pip) - set(expected_pip_normalized)),
        "binaries": binaries,
    })
    return result


def create_clean_environment(env_name: str, root: Path) -> tuple[Path | None, str | None]:
    """Create an isolated prefix from committed locks without re-solving."""
    prefix = root / env_name
    if prefix.exists():
        return None, f"clean prefix already exists: {prefix}"
    lock = ENVS_DIR / f"{env_name}.linux-64.lock"
    created = _conda(["create", "--yes", "--prefix", str(prefix), "--file", str(lock)])
    if created.returncode:
        return None, created.stderr.strip() or created.stdout.strip()

    pip_lock = ENVS_DIR / f"{env_name}.pip.lock"
    if pip_lock.exists():
        installed = _conda([
            "run", "--prefix", str(prefix), "python", "-m", "pip", "install",
            "--no-deps", "--requirement", str(pip_lock),
        ])
        if installed.returncode:
            return None, installed.stderr.strip() or installed.stdout.strip()
    return prefix, None


def audit_many(env_names: Iterable[str], clean_root: Path | None = None) -> list[dict[str, Any]]:
    reports = []
    for env_name in env_names:
        if clean_root is None:
            reports.append(audit_environment(env_name, env_name=env_name))
            continue
        prefix, error = create_clean_environment(env_name, clean_root)
        if error is not None:
            reports.append({
                "environment": env_name,
                "target": str(clean_root / env_name),
                "status": "error",
                "reason": "clean_install_failed",
                "details": error,
            })
            continue
        reports.append(audit_environment(env_name, prefix=prefix))
    return reports


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit exact Conda/PyPI versions and required ARIA binaries.",
    )
    parser.add_argument("--env", action="append", dest="env_names")
    parser.add_argument(
        "--clean-root", type=Path,
        help="Create fresh prefixes below this absent/empty directory before auditing.",
    )
    args = parser.parse_args(argv)
    if shutil.which("conda") is None:
        parser.error("conda is required")
    env_names = args.env_names or lock_environment_names()
    unknown = sorted(set(env_names) - set(lock_environment_names()))
    if unknown:
        parser.error(f"not active lock targets: {', '.join(unknown)}")
    if args.clean_root is not None:
        args.clean_root.mkdir(parents=True, exist_ok=True)
    reports = audit_many(env_names, clean_root=args.clean_root)
    print(json.dumps({"environments": reports}, indent=2, sort_keys=True))
    return 0 if all(report["status"] == "pass" for report in reports) else 1


if __name__ == "__main__":
    raise SystemExit(_main())
