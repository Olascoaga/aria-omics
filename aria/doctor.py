"""ARIA installation and integrity diagnostics."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import stat
import sys

from aria import __version__
from aria.utils.registry_integrity import IntegrityIssue, check_registry_integrity


def run_doctor(tier: str = "smoke") -> tuple[int, list[str]]:
    messages = [f"ARIA doctor {tier} (v{__version__})"]
    issues: list[IntegrityIssue] = []

    issues.extend(check_registry_integrity())
    issues.extend(_check_env_file_permissions())

    if tier in {"synthetic", "benchmark"}:
        issues.extend(_check_synthetic_assets())
    if tier == "benchmark":
        issues.append(IntegrityIssue(
            "benchmark_not_executed",
            "Benchmark tier requires an explicit dataset run; use tests/test_pbmc_e2e.py or the PBMC headless validation harness.",
            severity="warning",
        ))

    errors = [issue for issue in issues if issue.severity == "error"]
    warnings = [issue for issue in issues if issue.severity != "error"]

    for issue in errors:
        messages.append(f"ERROR {issue.code}: {issue.message}")
    for issue in warnings:
        messages.append(f"WARN {issue.code}: {issue.message}")

    if errors:
        messages.append(f"Result: failed ({len(errors)} errors, {len(warnings)} warnings)")
        return 1, messages
    messages.append(f"Result: passed ({len(warnings)} warnings)")
    return 0, messages


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="aria doctor")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--smoke", action="store_true", help="Fast import/registry/security checks.")
    group.add_argument("--synthetic", action="store_true", help="Smoke checks plus synthetic-test asset checks.")
    group.add_argument("--benchmark", action="store_true", help="Smoke/synthetic checks plus benchmark readiness warnings.")
    args = parser.parse_args(argv)

    tier = "smoke"
    if args.synthetic:
        tier = "synthetic"
    elif args.benchmark:
        tier = "benchmark"

    code, messages = run_doctor(tier)
    print("\n".join(messages))
    return code


def console_main() -> None:
    raise SystemExit(main())


def _check_env_file_permissions() -> list[IntegrityIssue]:
    env_path = Path(os.environ.get("ARIA_ENV_FILE") or Path.home() / ".aria" / ".env")
    if not env_path.exists():
        return []

    mode = stat.S_IMODE(env_path.stat().st_mode)
    if mode & (stat.S_IRWXG | stat.S_IRWXO):
        return [IntegrityIssue(
            "env_file_too_permissive",
            f"{env_path} should be readable only by the owner (chmod 600).",
        )]
    return []


def _check_synthetic_assets() -> list[IntegrityIssue]:
    root = Path(__file__).resolve().parents[1]
    required = (
        root / "tests" / "test_headless_design_e2e.py",
        root / "tests" / "test_pseudobulk_gate.py",
    )
    return [
        IntegrityIssue(
            "synthetic_asset_missing",
            f"Required synthetic validation asset missing: {path}",
        )
        for path in required
        if not path.exists()
    ]


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
