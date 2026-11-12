"""Typed, multi-platform scATAC FASTQ manifest loading and validation.

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
PLATFORM_ROLES = {
    "10x_atac": REQUIRED_ROLES,
    "10x_multiome_atac": REQUIRED_ROLES,
    "share_seq": REQUIRED_ROLES,
    "sci_atac": ("R1", "R2"),
    "surecell": ("R1", "R2"),
    "custom": ("R1",),
}
ORIENTATION_ALIASES = {
    "auto": "auto", "direct": "direct", "forward": "direct", "fwd": "direct",
    "reverse_complement": "reverse_complement", "rc": "reverse_complement",
    "revcomp": "reverse_complement",
}
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
        _structural_keys = {
            "library_id", "sample_id", "donor_id", "fastqs", "fastq_lanes",
            "platform", "barcode_orientation", "read_format",
            "barcode_whitelist", "metadata", "barcode_prefix",
            "fragments_file", "status", "alignment",
        }
        for extra_k, extra_v in row.items():
            if (
                extra_k not in _structural_keys
                and isinstance(extra_v, (str, int, float, bool))
                and extra_v is not None
            ):
                metadata.setdefault(extra_k, extra_v)

        platform = row.get("platform", raw.get("platform", "10x_atac"))
        if not isinstance(platform, str) or platform not in PLATFORM_ROLES:
            errors.append(f"{label}.platform is unknown: {platform!r}")
            platform = "10x_atac"
        orientation = row.get(
            "barcode_orientation", raw.get("barcode_orientation", "auto")
        )
        orientation = ORIENTATION_ALIASES.get(str(orientation).lower().strip())
        if orientation is None:
            errors.append(f"{label}.barcode_orientation is invalid")
            orientation = "auto"
        read_format = row.get("read_format", raw.get("read_format"))
        if read_format is not None:
            if not isinstance(read_format, str) or not read_format.strip():
                errors.append(f"{label}.read_format must be a non-empty string")
            else:
                read_format = read_format.strip()
        if platform == "custom" and not read_format:
            errors.append(f"{label}.read_format is required for custom platform")

        fastqs = row.get("fastqs")
        if not isinstance(fastqs, dict):
            errors.append(f"{label}.fastqs must be an object with R1/R2/R3")
            fastqs = {}
        normalized_fastqs = {}
        fastq_lanes = {}
        for role in REQUIRED_ROLES:
            value = fastqs.get(role, fastqs.get(role.lower()))
            if value is None and role not in PLATFORM_ROLES[platform]:
                continue
            lanes = value if isinstance(value, list) else [value]
            if not lanes:
                errors.append(f"{label}.fastqs.{role} must be a non-empty list")
            paths = []
            for lane in lanes:
                if not isinstance(lane, str) or not lane.strip():
                    errors.append(f"{label}.fastqs.{role} requires non-empty path strings")
                    continue
                path = _resolve_path(lane, base_dir)
                if path in assigned_fastqs:
                    errors.append(
                        f"FASTQ assigned more than once: {path} "
                        f"({assigned_fastqs[path]} and {library_id}/{role})"
                    )
                else:
                    assigned_fastqs[path] = f"{library_id}/{role}"
                if require_paths and not Path(path).is_file():
                    errors.append(f"{label}.fastqs.{role} not found: {path}")
                paths.append(path)
            if paths:
                fastq_lanes[role] = paths
                normalized_fastqs[role] = paths[0] if len(paths) == 1 else paths

        whitelist = _resolve_path(row.get("barcode_whitelist"), base_dir)
        if not whitelist and platform in {"10x_atac", "10x_multiome_atac"}:
            errors.append(f"{label}.barcode_whitelist is required")
        elif whitelist and require_paths and not Path(whitelist).is_file():
            errors.append(f"{label}.barcode_whitelist not found: {whitelist}")

        canonical.append({
            "library_id": library_id,
            "sample_id": sample_id,
            "donor_id": donor_id,
            "fastqs": normalized_fastqs,
            "fastq_lanes": fastq_lanes,
            "platform": platform,
            "barcode_orientation": orientation,
            **({"read_format": read_format} if read_format else {}),
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
    return str(path.resolve())


def manifest_library_types(manifest: dict[str, Any]) -> dict[str, str]:
    """Return authoritative per-file declarations for DataAudit classification."""
    return {path: "scATAC" for path in manifest_fastq_files(manifest)}


def manifest_fastq_files(manifest: dict[str, Any]) -> list[str]:
    """Return FASTQs in stable library then R1/R2/R3 order."""
    paths = []
    for row in (manifest or {}).get("libraries") or []:
        roles = row.get("fastq_lanes") or row.get("fastqs") or {}
        for role in REQUIRED_ROLES:
            value = roles.get(role)
            if value:
                paths.extend(value if isinstance(value, list) else [value])
    return paths


def first_manifest_fastq(manifest: dict[str, Any]) -> str | None:
    """Return the first declared lane in stable library/read/lane order."""
    return next(iter(manifest_fastq_files(manifest)), None)
