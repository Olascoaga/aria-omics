"""Single source of truth for ARIA's isolated scientific environments.

The registry owns stack-to-environment routing, timeouts, setup YAML identity,
lock policy, and required command-line tools. Dependency declarations remain in
the versioned ``envs/*.yml`` files; SetupAgent reads those files directly rather
than maintaining inline copies.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
ENVS_DIR = REPOSITORY_ROOT / "envs"


@dataclass(frozen=True)
class EnvironmentSpec:
    stack: str
    env_name: str
    yml_name: str
    description: str
    timeout_s: int
    required_binaries: tuple[str, ...] = ("python",)
    setup_managed: bool = True
    lock_required: bool = True
    lock_policy: str = "active_scientific_runtime"

    @property
    def yml_path(self) -> Path:
        return ENVS_DIR / self.yml_name


ENVIRONMENT_SPECS: dict[str, EnvironmentSpec] = {
    "rna": EnvironmentSpec(
        "rna", "aria-rna-env", "aria-rna-env.yml",
        "scRNA-seq, bulk DE, and shared DESeq2-backed DA", 3600,
    ),
    "ingestion": EnvironmentSpec(
        "ingestion", "aria-ingestion-env", "aria-ingestion-env.yml",
        "Raw scRNA ingestion with kallisto-bustools", 10800,
        ("python", "kb", "kallisto", "bustools"),
    ),
    "rnaseq": EnvironmentSpec(
        "rnaseq", "aria-rnaseq-env", "aria-rnaseq-env.yml",
        "Raw bulk RNA FASTQ alignment and quantification", 10800,
        ("python", "fastp", "STAR", "featureCounts", "samtools", "multiqc", "fastqc"),
    ),
    "atacseq": EnvironmentSpec(
        "atacseq", "aria-atacseq-env", "aria-atacseq-env.yml",
        "Raw bulk/scATAC FASTQ alignment", 14400,
        ("python", "bwa-mem2", "samtools", "chromap"),
    ),
    "chromatin": EnvironmentSpec(
        "chromatin", "aria-chromatin-env", "aria-chromatin-env.yml",
        "Chromatin QC, peak calling, clustering, and regulatory layers", 7200,
        ("python", "bedtools", "macs3"),
    ),
    "tobias": EnvironmentSpec(
        "tobias", "aria-tobias-env", "aria-tobias-env.yml",
        "Tn5-bias-corrected footprinting", 14400,
        ("python", "TOBIAS", "samtools"),
    ),
    "hic": EnvironmentSpec(
        "hic", "aria-hic-env", "aria-hic-env.yml",
        "Experimental Hi-C analysis", 14400,
        ("python", "cooler"), lock_required=False,
        lock_policy="scaffold_dispatch_disabled",
    ),
    "integration": EnvironmentSpec(
        "integration", "aria-integration-env", "aria-integration-env.yml",
        "Scaffold multi-omics integration", 7200,
        ("python", "R"), lock_required=False,
        lock_policy="scaffold_not_dispatchable",
    ),
    "benchmark": EnvironmentSpec(
        "benchmark", "aria-bench-env", "aria-bench-env.yml",
        "External RNA reference comparators", 14400,
        ("python", "Rscript"), setup_managed=False,
    ),
    "spatial": EnvironmentSpec(
        "spatial", "aria-rna-env", "aria-rna-env.yml",
        "Spatial transcriptomics on the validated RNA stack", 3600,
        setup_managed=False, lock_required=False,
        lock_policy="alias_of_rna",
    ),
}


def environment_names_by_stack() -> dict[str, str]:
    return {stack: spec.env_name for stack, spec in ENVIRONMENT_SPECS.items()}


def timeouts_by_stack() -> dict[str, int]:
    return {stack: spec.timeout_s for stack, spec in ENVIRONMENT_SPECS.items()}


def setup_environment_definitions() -> dict[str, dict[str, Any]]:
    """Return unique SetupAgent definitions backed by versioned YAML files."""
    definitions: dict[str, dict[str, Any]] = {}
    for spec in ENVIRONMENT_SPECS.values():
        if not spec.setup_managed:
            continue
        current = definitions.get(spec.env_name)
        definition = {
            "description": spec.description,
            "stack": spec.stack,
            "yml_path": str(spec.yml_path),
        }
        if current is not None and current["yml_path"] != definition["yml_path"]:
            raise ValueError(f"Conflicting YAMLs for {spec.env_name}")
        definitions.setdefault(spec.env_name, definition)
    return definitions


def lock_environment_names() -> list[str]:
    """Unique active environment names whose exact locks are release artifacts."""
    return sorted({
        spec.env_name for spec in ENVIRONMENT_SPECS.values()
        if spec.lock_required
    })


_CANONICAL_NAME_RE = re.compile(r"[-_.]+")


def _installed_direct_urls() -> dict[str, dict[str, Any]]:
    direct_urls = {}
    for distribution in importlib.metadata.distributions():
        raw = distribution.read_text("direct_url.json")
        name = distribution.metadata.get("Name")
        if not raw or not name:
            continue
        try:
            direct_urls[_CANONICAL_NAME_RE.sub("-", name).lower()] = json.loads(raw)
        except json.JSONDecodeError:
            continue
    return direct_urls


def portable_pip_lock_lines(
    records: list[dict[str, Any]],
    direct_urls: dict[str, dict[str, Any]] | None = None,
) -> list[str]:
    """Convert conda's PyPI records to portable, immutable requirement pins."""
    direct_urls = direct_urls or {}
    pins = set()
    for record in records:
        if str(record.get("channel", "")).lower() != "pypi":
            continue
        name = _CANONICAL_NAME_RE.sub("-", str(record.get("name", "")).strip()).lower()
        version = str(record.get("version", "")).strip()
        if not name or not version:
            continue
        direct_url = direct_urls.get(name, {})
        vcs_info = direct_url.get("vcs_info", {})
        commit = str(vcs_info.get("commit_id", "")).strip()
        vcs = str(vcs_info.get("vcs", "")).strip()
        url = str(direct_url.get("url", "")).strip()
        if vcs and url and commit:
            pins.add(f"{name} @ {vcs}+{url}@{commit}")
        elif url.startswith("file:"):
            raise ValueError(f"non-portable direct install for {name}: {url}")
        elif url:
            archive_info = direct_url.get("archive_info", {})
            hashes = archive_info.get("hashes", {})
            sha256 = str(hashes.get("sha256", "")).strip()
            if not sha256:
                legacy_hash = str(archive_info.get("hash", "")).strip()
                sha256 = legacy_hash.removeprefix("sha256=")
            if not url.startswith("https://") or not sha256:
                raise ValueError(f"unverifiable direct install for {name}: {url}")
            pins.add(f"{name} @ {url}#sha256={sha256}")
        else:
            pins.add(f"{name}=={version}")
    return sorted(pins, key=str.lower)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("lock-envs", "pip-lock"))
    args = parser.parse_args(argv)
    if args.action == "lock-envs":
        print("\n".join(lock_environment_names()))
        return 0
    records = json.load(sys.stdin)
    lines = portable_pip_lock_lines(records, _installed_direct_urls())
    if lines:
        sys.stdout.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
