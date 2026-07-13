"""Atomic download, archive extraction, and cache-manifest primitives.

These helpers are intentionally assay-agnostic.  A connector may stage many
files, but no staged byte becomes a public cache entry until every payload has
passed content validation and the retrieval manifest has been written.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import stat
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable
from urllib.parse import urlparse
from uuid import uuid4


MANIFEST_NAME = "retrieval_manifest.json"
SCHEMA_VERSION = "1"
DEFAULT_MAX_ARCHIVE_MEMBERS = 100_000
DEFAULT_MAX_EXTRACTED_BYTES = 500 * 1024**3


class RetrievalError(RuntimeError):
    """A public-data retrieval could not be validated or published."""


def open_url_with_retry(
    url: str,
    *,
    timeout: float = 120.0,
    max_attempts: int = 4,
):
    """Open a public URL with bounded retry for rate limits/transient failures."""
    if urlparse(url).scheme.lower() not in {"", "file"}:
        from aria.utils.privacy import assert_egress_allowed

        assert_egress_allowed("public data retrieval")
    last_error = None
    for attempt in range(max_attempts):
        try:
            return urllib.request.urlopen(url, timeout=timeout)
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code != 429 and not 500 <= exc.code < 600:
                raise
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            delay = float(retry_after) if retry_after and retry_after.isdigit() else 1.0
        except urllib.error.URLError as exc:
            last_error = exc
            delay = 1.0
        if attempt + 1 < max_attempts:
            time.sleep(delay * (attempt + 1))
    assert last_error is not None
    raise last_error


def download_atomic(
    url: str,
    destination: str | Path,
    *,
    expected_md5: str | None = None,
    expected_size: int | None = None,
    timeout: float = 120.0,
) -> dict:
    """Stream a URL to a private part file, verify it, then atomically replace."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.parent / f".{destination.name}.part-{uuid4().hex}"
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    size = 0
    try:
        with open_url_with_retry(url, timeout=timeout) as response, open(part, "wb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                md5.update(chunk)
                sha256.update(chunk)
                size += len(chunk)
            out.flush()
            os.fsync(out.fileno())
        if size <= 0:
            raise RetrievalError(f"downloaded payload is empty: {url}")
        if expected_size is not None and size != int(expected_size):
            raise RetrievalError(
                f"size mismatch for {url}: expected {expected_size}, got {size}"
            )
        observed_md5 = md5.hexdigest()
        if expected_md5 and observed_md5.lower() != expected_md5.strip().lower():
            raise RetrievalError(
                f"MD5 mismatch for {url}: expected {expected_md5}, got {observed_md5}"
            )
        os.replace(part, destination)
        return {
            "path": str(destination),
            "source_url": url,
            "size": size,
            "md5": observed_md5,
            "sha256": sha256.hexdigest(),
            "expected_md5": expected_md5,
            "expected_size": expected_size,
        }
    except RetrievalError:
        raise
    except Exception as exc:
        raise RetrievalError(f"download failed for {url}: {exc}") from exc
    finally:
        try:
            part.unlink()
        except FileNotFoundError:
            pass


def decompress_gzip_atomic(
    source: str | Path,
    destination: str | Path,
) -> dict:
    """Decompress one gzip payload privately, then atomically publish it."""
    source = Path(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.parent / f".{destination.name}.part-{uuid4().hex}"
    digest = hashlib.sha256()
    size = 0
    try:
        with gzip.open(source, "rb") as handle, open(part, "wb") as out:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                digest.update(chunk)
                size += len(chunk)
            out.flush()
            os.fsync(out.fileno())
        if size <= 0:
            raise RetrievalError(f"decompressed payload is empty: {source}")
        os.replace(part, destination)
        return {
            "path": str(destination),
            "compressed_source": str(source),
            "size": size,
            "sha256": digest.hexdigest(),
        }
    except RetrievalError:
        raise
    except (OSError, EOFError) as exc:
        raise RetrievalError(f"gzip decompression failed for {source}: {exc}") from exc
    finally:
        try:
            part.unlink()
        except FileNotFoundError:
            pass


def safe_extract_archive(
    archive_path: str | Path,
    destination: str | Path,
    *,
    max_members: int = DEFAULT_MAX_ARCHIVE_MEMBERS,
    max_unpacked_bytes: int = DEFAULT_MAX_EXTRACTED_BYTES,
) -> list[Path]:
    """Extract tar/zip without traversal, links, devices, or partial publication."""
    archive_path = Path(archive_path)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.extract-", dir=str(destination.parent)
    ))
    extracted: list[Path] = []
    try:
        if tarfile.is_tarfile(archive_path):
            with tarfile.open(archive_path, "r:*") as archive:
                members = archive.getmembers()
                _check_archive_limits(
                    len(members), sum(max(0, member.size) for member in members),
                    max_members, max_unpacked_bytes,
                )
                for member in members:
                    target = _safe_member_target(temp_root, member.name)
                    if member.isdir():
                        target.mkdir(parents=True, exist_ok=True)
                        continue
                    if not member.isfile():
                        raise RetrievalError(
                            f"unsafe archive member type: {member.name}"
                        )
                    target.parent.mkdir(parents=True, exist_ok=True)
                    source = archive.extractfile(member)
                    if source is None:
                        raise RetrievalError(f"cannot read archive member: {member.name}")
                    with source, open(target, "wb") as out:
                        shutil.copyfileobj(source, out, length=1024 * 1024)
                    extracted.append(target)
        elif zipfile.is_zipfile(archive_path):
            with zipfile.ZipFile(archive_path) as archive:
                members = archive.infolist()
                _check_archive_limits(
                    len(members), sum(max(0, member.file_size) for member in members),
                    max_members, max_unpacked_bytes,
                )
                for member in members:
                    target = _safe_member_target(temp_root, member.filename)
                    mode = member.external_attr >> 16
                    if stat.S_ISLNK(mode):
                        raise RetrievalError(
                            f"unsafe archive member type: {member.filename}"
                        )
                    if member.is_dir():
                        target.mkdir(parents=True, exist_ok=True)
                        continue
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(member, "r") as source, open(target, "wb") as out:
                        shutil.copyfileobj(source, out, length=1024 * 1024)
                    extracted.append(target)
        else:
            raise RetrievalError(f"unsupported or corrupt archive: {archive_path}")

        if destination.exists():
            raise RetrievalError(f"archive destination already exists: {destination}")
        relative = [path.relative_to(temp_root) for path in extracted]
        os.replace(temp_root, destination)
        return [destination / path for path in relative]
    except RetrievalError:
        raise
    except Exception as exc:
        raise RetrievalError(f"archive extraction failed for {archive_path}: {exc}") from exc
    finally:
        if temp_root.exists():
            shutil.rmtree(temp_root, ignore_errors=True)


def validate_payload(path: str | Path, kind: str | None = None) -> dict:
    """Perform a bounded format check; manifest hashing supplies full-byte integrity."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size <= 0:
        raise RetrievalError(f"payload is missing or empty: {path}")

    kind = kind or "unknown"
    suffix = path.suffix.lower()
    # BAM is BGZF (concatenated gzip members), so its magic appears only after
    # decompression. CRAM is uncompressed and shares the connector's ``bam``
    # bucket for downstream routing.
    is_bam = kind == "bam" and suffix == ".bam"
    is_cram = kind in {"bam", "cram"} and suffix == ".cram"
    is_gzip = path.name.lower().endswith(".gz") or is_bam
    try:
        if is_gzip:
            # Read the whole stream so gzip CRC/truncation errors cannot be cached.
            with gzip.open(path, "rb") as handle:
                first = handle.read(8192)
                while handle.read(1024 * 1024):
                    pass
        else:
            with open(path, "rb") as handle:
                first = handle.read(8192)
    except (OSError, EOFError) as exc:
        raise RetrievalError(f"gzip validation failed for {path}: {exc}") from exc

    if kind in {"h5", "h5ad"} and not first.startswith(b"\x89HDF\r\n\x1a\n"):
        raise RetrievalError(f"invalid HDF5 signature: {path}")
    if kind == "mtx" and not first.lstrip().startswith(b"%%MatrixMarket"):
        raise RetrievalError(f"invalid MatrixMarket header: {path}")
    if (kind == "cram" or is_cram) and not first.startswith(b"CRAM"):
        raise RetrievalError(f"invalid CRAM signature: {path}")
    if kind == "bam" and not is_cram and not first.startswith(b"BAM\x01"):
        raise RetrievalError(f"invalid BAM signature: {path}")
    if kind == "fastq":
        _validate_fastq(path)
    if not first:
        raise RetrievalError(f"payload has no readable content: {path}")
    return {"path": str(path), "kind": kind, "size": path.stat().st_size}


def write_retrieval_manifest(
    root: str | Path,
    *,
    accession: str,
    payloads: Iterable[str | Path],
    source_accessions: Iterable[str] | None = None,
    sources: dict[str, dict] | None = None,
) -> Path:
    """Hash every published payload and atomically write the retrieval manifest."""
    root = Path(root).resolve()
    records = []
    seen: set[str] = set()
    for raw_path in payloads:
        path = Path(raw_path).resolve()
        try:
            relative = path.relative_to(root).as_posix()
        except ValueError as exc:
            raise RetrievalError(f"manifest payload escapes retrieval root: {path}") from exc
        if relative == MANIFEST_NAME or relative in seen:
            continue
        seen.add(relative)
        if not path.is_file():
            raise RetrievalError(f"manifest payload missing: {path}")
        record = {
            "path": relative,
            "size": path.stat().st_size,
            "sha256": _sha256_file(path),
        }
        if sources and relative in sources:
            record["source"] = sources[relative]
        records.append(record)
    if not records:
        raise RetrievalError("retrieval manifest requires at least one payload")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "accession": accession,
        "source_accessions": list(source_accessions or [accession]),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "files": sorted(records, key=lambda row: row["path"]),
    }
    path = root / MANIFEST_NAME
    temp = root / f".{MANIFEST_NAME}.tmp-{uuid4().hex}"
    temp.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temp, path)
    return path


def validate_retrieval_manifest(
    root: str | Path,
    *,
    expected_accession: str | None = None,
) -> dict:
    """Validate schema, containment, size, and SHA-256 for a published cache entry."""
    root = Path(root).resolve()
    manifest_path = root / MANIFEST_NAME
    errors: list[str] = []
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"status": "invalid", "errors": [f"cannot read manifest: {exc}"]}
    if str(manifest.get("schema_version")) != SCHEMA_VERSION:
        errors.append("unsupported retrieval manifest schema_version")
    if (
        expected_accession is not None
        and str(manifest.get("accession") or "").upper()
        != str(expected_accession).upper()
    ):
        errors.append(
            "accession mismatch: expected "
            f"{expected_accession}, got {manifest.get('accession')!r}"
        )
    files = manifest.get("files")
    if not isinstance(files, list) or not files:
        errors.append("retrieval manifest has no files")
        files = []
    manifested_paths: set[str] = set()
    for record in files:
        relative = str(record.get("path") or "")
        if relative in manifested_paths:
            errors.append(f"duplicate manifest payload: {relative}")
            continue
        manifested_paths.add(relative)
        try:
            path = _safe_member_target(root, relative)
        except RetrievalError as exc:
            errors.append(str(exc))
            continue
        if not path.is_file():
            errors.append(f"missing payload: {relative}")
            continue
        observed_size = path.stat().st_size
        if observed_size != record.get("size"):
            errors.append(
                f"size mismatch for {relative}: expected {record.get('size')}, "
                f"got {observed_size}"
            )
            continue
        observed_sha = _sha256_file(path)
        if observed_sha != record.get("sha256"):
            errors.append(
                f"sha256 mismatch for {relative}: expected {record.get('sha256')}, "
                f"got {observed_sha}"
            )
    actual_paths = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.name != MANIFEST_NAME
    }
    for relative in sorted(actual_paths - manifested_paths):
        errors.append(f"unmanifested payload: {relative}")
    return {
        "status": "valid" if not errors else "invalid",
        "errors": errors,
        "manifest": manifest,
        "manifest_path": str(manifest_path),
    }


def publish_directory_atomic(staging: str | Path, target: str | Path) -> Path:
    """Publish a fully validated staging tree, rolling back an invalid old cache."""
    staging = Path(staging)
    target = Path(target)
    if validate_retrieval_manifest(staging).get("status") != "valid":
        raise RetrievalError(f"staging manifest is invalid: {staging}")
    if not target.exists():
        os.replace(staging, target)
        return target
    if validate_retrieval_manifest(target).get("status") == "valid":
        shutil.rmtree(staging, ignore_errors=True)
        return target

    backup = target.parent / f".{target.name}.invalid-{uuid4().hex}"
    os.replace(target, backup)
    try:
        os.replace(staging, target)
    except Exception:
        os.replace(backup, target)
        raise
    shutil.rmtree(backup, ignore_errors=True)
    return target


def _validate_fastq(path: Path, max_records: int = 1000) -> None:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    records = 0
    try:
        with opener(path, "rt", encoding="ascii", errors="strict") as handle:
            while records < max_records:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip("\r\n")
                plus = handle.readline()
                quality = handle.readline().rstrip("\r\n")
                if (
                    not header.startswith("@")
                    or not plus.startswith("+")
                    or not sequence
                    or len(sequence) != len(quality)
                ):
                    raise RetrievalError(f"invalid FASTQ record in {path}")
                records += 1
    except (OSError, UnicodeError) as exc:
        raise RetrievalError(f"FASTQ validation failed for {path}: {exc}") from exc
    if records == 0:
        raise RetrievalError(f"FASTQ contains no records: {path}")


def _check_archive_limits(
    n_members: int,
    unpacked_bytes: int,
    max_members: int,
    max_unpacked_bytes: int,
) -> None:
    if n_members > max_members:
        raise RetrievalError(
            f"archive has {n_members} members; limit is {max_members}"
        )
    if unpacked_bytes > max_unpacked_bytes:
        raise RetrievalError(
            f"archive expands to {unpacked_bytes} bytes; limit is {max_unpacked_bytes}"
        )


def _safe_member_target(root: Path, name: str) -> Path:
    member = PurePosixPath(str(name).replace("\\", "/"))
    if member.is_absolute() or not member.parts or any(part in {"", ".", ".."} for part in member.parts):
        raise RetrievalError(f"unsafe archive member: {name}")
    target = root.joinpath(*member.parts).resolve()
    try:
        target.relative_to(root.resolve())
    except ValueError as exc:
        raise RetrievalError(f"unsafe archive member: {name}") from exc
    return target


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()
