"""Reference-resource checksum integrity checks.

ARIA stores large reusable references (GMT gene sets, MEME motif collections,
genome FASTAs) outside the repository. When a resource has a manifest with an
expected SHA-256, the runtime must verify it before using the file.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


_TERMINAL_BAD = {"file_missing", "manifest_unreadable", "checksum_mismatch"}


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest for ``path`` using bounded memory."""
    digest = hashlib.sha256()
    with open(Path(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_manifest(path: Path) -> tuple[dict[str, Any] | None, str | None]:
    if not path.is_file():
        return None, None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"
    return data if isinstance(data, dict) else {}, None


def _sha_from_entry(entry: Any) -> str | None:
    if isinstance(entry, str) and entry.strip():
        return entry.strip()
    if isinstance(entry, dict):
        value = entry.get("sha256") or entry.get("sha256sum")
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def expected_sha256(manifest: dict[str, Any], resource_path: str | Path) -> str | None:
    """Find the expected digest for ``resource_path`` in a flexible manifest."""
    direct = _sha_from_entry(manifest.get("sha256") or manifest.get("sha256sum"))
    if direct:
        return direct

    path = Path(resource_path)
    keys = {path.name, path.as_posix(), str(path)}
    files = manifest.get("files")
    if isinstance(files, dict):
        for key, entry in files.items():
            key_text = str(key)
            if key_text in keys or Path(key_text).name == path.name:
                digest = _sha_from_entry(entry)
                if digest:
                    return digest
    elif isinstance(files, list):
        for entry in files:
            if not isinstance(entry, dict):
                continue
            entry_path = entry.get("path") or entry.get("name") or entry.get("file")
            if entry_path and (str(entry_path) in keys or Path(str(entry_path)).name == path.name):
                digest = _sha_from_entry(entry)
                if digest:
                    return digest

    for key in (f"{path.stem}_sha256", f"{path.name}_sha256"):
        digest = _sha_from_entry(manifest.get(key))
        if digest:
            return digest
    return None


def verify_reference_file(
    resource_path: str | Path,
    *,
    manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """Verify ``resource_path`` against an adjacent or explicit manifest.

    Status values:
    - ``ok``: expected SHA-256 exists and matches.
    - ``checksum_mismatch``: expected SHA-256 exists and differs.
    - ``checksum_not_declared``: manifest exists but declares no SHA-256.
    - ``manifest_missing``: no manifest exists.
    - ``manifest_unreadable`` / ``file_missing``: resource cannot be verified.
    """
    path = Path(resource_path).expanduser()
    manifest = Path(manifest_path).expanduser() if manifest_path else path.parent / "manifest.json"
    result: dict[str, Any] = {
        "path": str(path),
        "manifest_path": str(manifest),
        "status": "unknown",
        "expected_sha256": None,
        "observed_sha256": None,
    }
    if not path.is_file():
        result["status"] = "file_missing"
        return result

    manifest_data, error = _load_manifest(manifest)
    if error:
        result["status"] = "manifest_unreadable"
        result["error"] = error
        return result
    if manifest_data is None:
        result["status"] = "manifest_missing"
        return result

    expected = expected_sha256(manifest_data, path)
    result["manifest"] = manifest_data
    if not expected:
        result["status"] = "checksum_not_declared"
        return result

    observed = sha256_file(path)
    result["expected_sha256"] = expected.lower()
    result["observed_sha256"] = observed.lower()
    result["status"] = "ok" if observed.lower() == expected.lower() else "checksum_mismatch"
    return result


def public_integrity_result(result: dict[str, Any] | None) -> dict[str, Any] | None:
    """Return a compact, report-safe integrity payload."""
    if not result:
        return None
    return {
        key: result.get(key)
        for key in (
            "path",
            "manifest_path",
            "status",
            "expected_sha256",
            "observed_sha256",
            "error",
        )
        if result.get(key) is not None
    }


def reference_is_usable(result: dict[str, Any] | None) -> bool:
    """True when the integrity result does not prove the resource is bad."""
    if not result:
        return True
    return result.get("status") not in _TERMINAL_BAD


def _finding(severity: str, check: str, message: str, recommendation: str) -> dict[str, Any]:
    return {
        "severity": severity,
        "check": check,
        "message": message,
        "recommendation": recommendation,
        "modality": "reference",
    }


def assess_reference_integrity(exp_context: dict[str, Any]) -> dict[str, Any]:
    """Data-light preflight over automatically resolved local references."""
    checks: list[dict[str, Any]] = []
    findings: list[dict[str, Any]] = []

    def add(kind: str, name: str, result: dict[str, Any]) -> None:
        public = public_integrity_result(result) or {}
        public.update({"kind": kind, "name": name})
        checks.append(public)
        status = result.get("status")
        if status in {"checksum_mismatch", "manifest_unreadable"}:
            findings.append(_finding(
                "blocking",
                f"reference_{status}",
                f"{kind} reference '{name}' failed checksum integrity validation ({status}).",
                "Re-stage the reference from its governed fetch script or fix the manifest.",
            ))
        elif status == "checksum_not_declared":
            findings.append(_finding(
                "warning",
                "reference_checksum_not_declared",
                f"{kind} reference '{name}' could not be fully checksum-verified.",
                "Use a manifest with a SHA-256 generated by ARIA's reference bootstrap scripts.",
            ))

    try:
        from aria.utils import ora
        base = ora.genesets_dir()
        if base.exists():
            for gmt in sorted(base.glob("*/*.gmt")):
                add("geneset", gmt.parent.name, verify_reference_file(gmt))
    except Exception as exc:
        findings.append(_finding(
            "warning",
            "reference_integrity_genesets_uninspectable",
            f"Could not inspect local GMT references: {type(exc).__name__}: {exc}",
            "Check ARIA_GMT_DIR permissions and manifest format.",
        ))

    modalities = set((exp_context.get("modalities") or {}).keys())
    chromatin = bool(modalities & {"scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"})
    if chromatin:
        try:
            from aria.utils import motifs
            collection = exp_context.get("motif_collection") or motifs.DEFAULT_COLLECTION
            meme_path, _manifest_path = motifs.collection_paths(str(collection))
            if meme_path.is_file():
                add("motif", str(collection), verify_reference_file(meme_path))
        except Exception as exc:
            findings.append(_finding(
                "warning",
                "reference_integrity_motifs_uninspectable",
                f"Could not inspect local motif references: {type(exc).__name__}: {exc}",
                "Check ARIA_MOTIF_DIR permissions and manifest format.",
            ))

        assembly = exp_context.get("genome") or exp_context.get("assembly")
        if assembly:
            try:
                from aria.utils import genomes
                resolved = getattr(genomes, "resolve_local_genome_fasta_with_integrity", None)
                if resolved:
                    fasta, source, integrity = resolved(str(assembly))
                    if integrity:
                        name = f"{assembly} ({source or 'unresolved'})"
                        add("genome", name, integrity)
                    elif fasta:
                        add("genome", str(assembly), verify_reference_file(fasta))
            except Exception as exc:
                findings.append(_finding(
                    "warning",
                    "reference_integrity_genome_uninspectable",
                    f"Could not inspect local genome reference: {type(exc).__name__}: {exc}",
                    "Check ARIA_GENOME_DIR permissions and manifest format.",
                ))

    if not checks and not findings:
        return {"status": "not_applicable", "checks": [], "findings": []}

    if any((check.get("status") in {"checksum_mismatch", "manifest_unreadable"}) for check in checks):
        status = "red"
    elif findings:
        status = "yellow"
    else:
        status = "green"
    return {"status": status, "checks": checks, "findings": findings}
