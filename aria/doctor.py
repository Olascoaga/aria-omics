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

    if tier == "secrets":
        # P2-9: focused secret-hygiene diagnostics (no registry/benchmark).
        issues.extend(_check_secrets())
        issues.extend(_check_env_file_permissions())
    elif tier == "llm":
        # P2-9: focused LLM provider/offline/fallback diagnostics.
        issues.extend(_check_llm())
    else:
        issues.extend(check_registry_integrity())
        issues.extend(_check_env_file_permissions())
        if tier in {"synthetic", "benchmark"}:
            issues.extend(_check_synthetic_assets())
        if tier == "benchmark":
            issues.extend(_run_synthetic_de_benchmark())

    errors = [i for i in issues if i.severity == "error"]
    warnings = [i for i in issues if i.severity == "warning"]
    infos = [i for i in issues if i.severity == "info"]

    for issue in infos:
        messages.append(f"INFO {issue.code}: {issue.message}")
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
    group.add_argument("--secrets", action="store_true", help="Secret hygiene: missing/malformed keys + leaked credentials in project files.")
    group.add_argument("--llm", action="store_true", help="LLM provider/model/offline/fallback diagnostics (optional gated latency probe).")
    args = parser.parse_args(argv)

    tier = "smoke"
    if args.synthetic:
        tier = "synthetic"
    elif args.benchmark:
        tier = "benchmark"
    elif args.secrets:
        tier = "secrets"
    elif args.llm:
        tier = "llm"

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


def _run_synthetic_de_benchmark() -> list[IntegrityIssue]:
    """X6: run the synthetic ground-truth DE benchmark when the DE stack is
    available. Skips gracefully (warning) when pydeseq2/anndata are absent
    (e.g. the base aria-env); the gate then runs in aria-rna-env / CI."""
    try:
        import anndata  # noqa: F401
        import pydeseq2  # noqa: F401
    except Exception:
        return [IntegrityIssue(
            "benchmark_skipped",
            "Synthetic DE benchmark skipped: pydeseq2/anndata not in this env. "
            "Run in aria-rna-env (or CI) to exercise numerical-accuracy recovery.",
            severity="warning",
        )]
    issues: list[IntegrityIssue] = []
    # W-CALIB single source: recovery (recall + empirical FDR) AND the
    # label-permutation negative control (false-positive rate under the null)
    # for both DE paths. `quick` uses small/fast configs; the pytest gate uses
    # the full, powered configs.
    try:
        from aria.benchmarks.synthetic_de import run_calibration_suite
        manifest = run_calibration_suite(seed=11, quick=True)
    except Exception as exc:  # never let doctor crash on the benchmark
        return [IntegrityIssue(
            "benchmark_error",
            f"Synthetic DE calibration suite could not run: {exc}",
            severity="warning",
        )]
    for path, blocks in manifest.get("paths", {}).items():
        recovery = blocks.get("recovery", {})
        if recovery.get("status") != "pass":
            msgs = recovery.get("messages") or [recovery]
            issues.append(IntegrityIssue(
                "benchmark_failed",
                f"Synthetic {path} DE recovery out of tolerance: {msgs[0]}",
                severity="error",
            ))
        neg = blocks.get("negative_control", {})
        if neg.get("status") != "pass":
            msgs = neg.get("messages") or [neg]
            issues.append(IntegrityIssue(
                "calibration_failed",
                f"Synthetic {path} DE negative control out of tolerance "
                f"(over-calls under permuted null): {msgs[0]}",
                severity="error",
            ))
    # Pass: clean (no issue). The pytest gate records the positive metrics.
    return issues


def _check_secrets() -> list[IntegrityIssue]:
    """P2-9: API key presence/format + leaked-credential scan of project files.

    Absent keys are informational (offline is valid); a malformed key is a
    warning; a credential committed to a project file is an error. No key value
    is ever printed (masked only).
    """
    from aria.utils.secret_hygiene import (
        PROVIDER_ENV, classify_key, mask_secret, scan_paths_for_secrets,
    )

    issues: list[IntegrityIssue] = []
    for provider, env_var in PROVIDER_ENV.items():
        value = os.environ.get(env_var)
        state = classify_key(provider, value)
        if state == "absent":
            issues.append(IntegrityIssue(
                "key_absent", f"{env_var} not set ({provider}).",
                severity="info"))
        elif state == "ok":
            issues.append(IntegrityIssue(
                "key_present",
                f"{env_var} present ({mask_secret(value)}).",
                severity="info"))
        else:
            issues.append(IntegrityIssue(
                "key_malformed",
                f"{env_var} is set but looks malformed for {provider} "
                f"({mask_secret(value)}); verify it was pasted correctly.",
                severity="warning"))

    # Leaked-credential scan over the project tree (bounded). The legitimate key
    # store ~/.aria/.env is OUTSIDE the repo, so any hit here is a real concern.
    repo_root = Path(__file__).resolve().parents[1]
    _SKIP_DIRS = {".git", "memory", "graphify-out", "build", "dist",
                  ".eggs", "__pycache__", ".pytest_cache", "node_modules"}
    candidates: list[Path] = []
    for path in repo_root.rglob("*"):
        if not path.is_file():
            continue
        if any(part in _SKIP_DIRS for part in path.relative_to(repo_root).parts):
            continue
        candidates.append(path)
        if len(candidates) >= 5000:  # bound the walk
            break
    for hit in scan_paths_for_secrets(candidates):
        rel = Path(hit["path"]).relative_to(repo_root)
        issues.append(IntegrityIssue(
            "credential_in_project_file",
            f"Possible {hit['kind']} credential committed in {rel} "
            f"(matched {hit['match']}). Remove it and rotate the key.",
            severity="error"))
    return issues


def _check_llm() -> list[IntegrityIssue]:
    """P2-9: provider/model/offline/fallback diagnostics.

    litellm-free: reads ~/.aria/config.yaml + env directly so it runs in any
    env. A latency probe is OFF unless egress is allowed and a key is present;
    it is a plain TCP reachability timing (no LLM call, no cost).
    """
    from aria.utils.secret_hygiene import PROVIDER_ENV, classify_key

    issues: list[IntegrityIssue] = []
    try:
        from aria.utils.privacy import air_gapped_enabled, egress_allowed
    except Exception:
        air_gapped_enabled = lambda: bool(os.environ.get("ARIA_AIR_GAPPED"))  # noqa: E731
        egress_allowed = lambda *a, **k: not air_gapped_enabled()  # noqa: E731

    air_gapped = bool(air_gapped_enabled())
    issues.append(IntegrityIssue(
        "llm_mode",
        f"air-gapped mode is {'ON (cloud egress refused)' if air_gapped else 'OFF'}.",
        severity="info"))

    # Configured providers (key presence only, masked-free; never the value).
    configured = []
    for provider, env_var in PROVIDER_ENV.items():
        if classify_key(provider, os.environ.get(env_var)) == "ok":
            configured.append(provider)
    issues.append(IntegrityIssue(
        "llm_keys",
        f"providers with a valid-looking key: "
        f"{', '.join(sorted(set(configured))) or 'none (offline only)'}.",
        severity="info"))

    # Tier -> model map from config.yaml when present (no provider import).
    tier_models = _read_llm_config_models()
    if tier_models:
        for tier, entries in tier_models.items():
            chain = " -> ".join(entries)
            issues.append(IntegrityIssue(
                f"llm_tier_{tier}",
                f"{tier}: {chain}", severity="info"))
    else:
        issues.append(IntegrityIssue(
            "llm_config",
            "no ~/.aria/config.yaml llm section; built-in defaults are used "
            "(heavy=frontier, medium/light fall back to local when air-gapped).",
            severity="info"))

    # Offline readiness: air-gapped needs a local model somewhere.
    if air_gapped and tier_models and not any(
        "local" in e or "ollama" in e for entries in tier_models.values()
        for e in entries
    ):
        issues.append(IntegrityIssue(
            "llm_offline_unready",
            "air-gapped is ON but no local/ollama model is configured; "
            "LLM calls will be refused.", severity="warning"))

    # Optional, gated latency probe: plain TCP, no LLM call, no cost.
    if os.environ.get("ARIA_DOCTOR_LLM_PROBE") == "1":
        if air_gapped or not egress_allowed("llm"):
            issues.append(IntegrityIssue(
                "llm_latency",
                "latency probe skipped (egress not allowed).", severity="info"))
        else:
            for provider, host in (("anthropic", "api.anthropic.com"),
                                   ("openai", "api.openai.com"),
                                   ("google", "generativelanguage.googleapis.com")):
                if provider not in configured:
                    continue
                issues.append(_probe_latency(provider, host))
    else:
        issues.append(IntegrityIssue(
            "llm_latency",
            "latency probe off (set ARIA_DOCTOR_LLM_PROBE=1 to TCP-ping "
            "configured providers; never makes a billed call).",
            severity="info"))
    return issues


def _read_llm_config_models() -> dict[str, list[str]]:
    """Read tier -> ['provider/model', ...] from ~/.aria/config.yaml if present."""
    cfg = Path(os.environ.get("ARIA_CONFIG_FILE")
               or Path.home() / ".aria" / "config.yaml")
    if not cfg.exists():
        return {}
    try:
        import yaml
        data = yaml.safe_load(cfg.read_text()) or {}
    except Exception:
        return {}
    llm = (data.get("llm") or {}) if isinstance(data, dict) else {}
    out: dict[str, list[str]] = {}
    for tier in ("heavy", "medium", "light"):
        spec = llm.get(tier)
        if not isinstance(spec, dict):
            continue
        provider = str(spec.get("provider", "?"))
        model = str(spec.get("model", "?"))
        tag = "local" if spec.get("api_base") else "cloud"
        out[tier] = [f"{provider}/{model} ({tag})"]
    return out


def _probe_latency(provider: str, host: str) -> IntegrityIssue:
    import socket
    import time
    start = time.monotonic()
    try:
        with socket.create_connection((host, 443), timeout=3.0):
            ms = (time.monotonic() - start) * 1000.0
        return IntegrityIssue(
            "llm_latency",
            f"{provider} ({host}) reachable in {ms:.0f} ms (TCP).",
            severity="info")
    except Exception as exc:
        return IntegrityIssue(
            "llm_latency",
            f"{provider} ({host}) not reachable: {type(exc).__name__}.",
            severity="warning")


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
