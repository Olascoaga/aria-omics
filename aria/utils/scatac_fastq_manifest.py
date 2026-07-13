"""Typed 10x scATAC FASTQ manifest loading and validation.

The manifest is the authoritative mapping from biological libraries to the
three 10x ATAC read roles. Filenames remain non-binding hints (ADR-001/E2).
"""

from __future__ import annotations

import json
import re
from copy import deepcopy
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1"
REQUIRED_ROLES = ("R1", "R2", "R3")
DISCOVERY_NAMES = (
    "scatac_fastq_manifest.json",
    "scatac_manifest.json",
)
_SAFE_LIBRARY_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def resolve_scatac_fastq_manifest(
    context: dict[str, Any] | None,
    *,
    require_paths: bool = False,
) -> dict[str, Any]:
    """Resolve and validate an inline, referenced, or discovered manifest."""
    context = context or {}
    source = context.get("scatac_fastq_manifest")
    source_path = context.get("scatac_fastq_manifest_path")
    data_dir = Path(context.get("data_dir") or ".")

    if source is None and source_path:
        source = source_path
    if source is None:
        for filename in DISCOVERY_NAMES:
            candidate = data_dir / filename
            if candidate.is_file():
                source = str(candidate)
                break
    if source is None:
        return {"status": "missing", "manifest": None, "errors": []}

    base_dir = data_dir
    resolved_source = "inline"
    if isinstance(source, (str, Path)):
        path = Path(source).expanduser()
        if not path.is_absolute():
            path = data_dir / path
        resolved_source = str(path)
        base_dir = path.parent
        try:
            source = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            return {
                "status": "invalid",
                "manifest": None,
                "source": resolved_source,
                "errors": [f"cannot read manifest JSON: {exc}"],
            }

    manifest, errors = _canonicalize_manifest(
        source,
        base_dir=base_dir,
        require_paths=require_paths,
    )
    return {
        "status": "valid" if not errors else "invalid",
        "manifest": manifest if not errors else None,
        "source": resolved_source,
        "errors": errors,
    }


def _canonicalize_manifest(
    raw: Any,
    *,
    base_dir: Path,
    require_paths: bool,
) -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    if not isinstance(raw, dict):
        return {}, ["manifest must be a JSON object"]
    version = str(raw.get("schema_version") or "")
    if version != SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {SCHEMA_VERSION!r}; got {version!r}"
        )
    libraries = raw.get("libraries")
    if not isinstance(libraries, list) or not libraries:
        return {}, errors + ["libraries must be a non-empty list"]

    canonical = []
    library_ids: set[str] = set()
    assigned_fastqs: dict[str, str] = {}
    for index, row in enumerate(libraries, start=1):
        label = f"libraries[{index - 1}]"
        if not isinstance(row, dict):
            errors.append(f"{label} must be an object")
            continue
        library_id = str(row.get("library_id") or "").strip()
        sample_id = str(row.get("sample_id") or "").strip()
        metadata = deepcopy(row.get("metadata") or {})
        if not isinstance(metadata, dict):
            errors.append(f"{label}.metadata must be an object")
            metadata = {}
        donor_id = str(row.get("donor_id") or metadata.get("donor_id") or "").strip()
        if not library_id:
            errors.append(f"{label}.library_id is required")
        elif not _SAFE_LIBRARY_ID.fullmatch(library_id):
            errors.append(
                f"{label}.library_id must use only letters, numbers, '.', "
                "'_' or '-'"
            )
        elif library_id in library_ids:
            errors.append(f"duplicate library_id: {library_id}")
        else:
            library_ids.add(library_id)
        if not sample_id:
            errors.append(f"{label}.sample_id is required")
        if not donor_id:
            errors.append(f"{label}.donor_id is required")
        metadata["donor_id"] = donor_id

        fastqs = row.get("fastqs")
        if not isinstance(fastqs, dict):
            errors.append(f"{label}.fastqs must be an object with R1/R2/R3")
            fastqs = {}
        normalized_fastqs = {}
        for role in REQUIRED_ROLES:
            value = fastqs.get(role) or fastqs.get(role.lower())
            if isinstance(value, list):
                if len(value) != 1:
                    errors.append(
                        f"{label}.fastqs.{role} must name exactly one FASTQ"
                    )
                    continue
                value = value[0]
            path = _resolve_path(value, base_dir)
            if not path:
                errors.append(f"{label}.fastqs.{role} is required")
                continue
            if path in assigned_fastqs:
                errors.append(
                    f"FASTQ assigned more than once: {path} "
                    f"({assigned_fastqs[path]} and {library_id}/{role})"
                )
            else:
                assigned_fastqs[path] = f"{library_id}/{role}"
            if require_paths and not Path(path).is_file():
                errors.append(f"{label}.fastqs.{role} not found: {path}")
            normalized_fastqs[role] = path

        whitelist = _resolve_path(row.get("barcode_whitelist"), base_dir)
        if not whitelist:
            errors.append(f"{label}.barcode_whitelist is required")
        elif require_paths and not Path(whitelist).is_file():
            errors.append(f"{label}.barcode_whitelist not found: {whitelist}")

        canonical.append({
            "library_id": library_id,
            "sample_id": sample_id,
            "donor_id": donor_id,
            "fastqs": normalized_fastqs,
            "barcode_whitelist": whitelist,
            "metadata": metadata,
        })

    return {
        "schema_version": SCHEMA_VERSION,
        "libraries": canonical,
    }, errors


def _resolve_path(value: Any, base_dir: Path) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return str(path)


def manifest_library_types(manifest: dict[str, Any]) -> dict[str, str]:
    """Return authoritative per-file declarations for DataAudit classification."""
    declarations = {}
    for row in (manifest or {}).get("libraries") or []:
        for path in (row.get("fastqs") or {}).values():
            if path:
                declarations[str(path)] = "scATAC"
    return declarations


def manifest_fastq_files(manifest: dict[str, Any]) -> list[str]:
    """Return FASTQs in stable library then R1/R2/R3 order."""
    paths = []
    for row in (manifest or {}).get("libraries") or []:
        roles = row.get("fastqs") or {}
        paths.extend(str(roles[role]) for role in REQUIRED_ROLES if roles.get(role))
    return paths
