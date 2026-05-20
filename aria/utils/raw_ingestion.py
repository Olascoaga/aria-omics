"""Deterministic raw scRNA ingestion helpers."""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from aria.utils.provenance import hash_file, hash_params


TENX_BARCODE = "barcodes.tsv.gz"
TENX_MATRIX = "matrix.mtx.gz"
TENX_FEATURES = ("features.tsv.gz", "genes.tsv.gz")


@dataclass(frozen=True)
class TenXTriplet:
    directory: Path
    barcodes: Path
    features: Path
    matrix: Path
    sample_id: str
    n_features: int
    n_barcodes: int
    n_nonzero: int

    def source_paths(self) -> list[Path]:
        return [self.barcodes, self.features, self.matrix]


def discover_10x_mtx_triplets(root: str | Path) -> list[TenXTriplet]:
    """Find valid Cell Ranger-style 10X matrix triplets under *root*."""
    root = Path(root)
    if root.is_file():
        root = root.parent
    candidates: set[Path] = set()
    for matrix in root.rglob(TENX_MATRIX):
        candidates.add(matrix.parent)
    if (root / TENX_MATRIX).exists():
        candidates.add(root)

    triplets: list[TenXTriplet] = []
    seen: set[Path] = set()
    for directory in sorted(candidates):
        try:
            triplet = validate_10x_mtx_triplet(directory)
        except ValueError:
            continue
        if triplet.directory not in seen:
            triplets.append(triplet)
            seen.add(triplet.directory)
    return triplets


def validate_10x_mtx_triplet(directory: str | Path) -> TenXTriplet:
    """Validate one 10X matrix directory and return its metadata."""
    directory = Path(directory)
    barcodes = directory / TENX_BARCODE
    matrix = directory / TENX_MATRIX
    features = next((directory / name for name in TENX_FEATURES
                     if (directory / name).exists()), None)

    missing = []
    if not barcodes.exists():
        missing.append(TENX_BARCODE)
    if features is None:
        missing.append("features.tsv.gz or genes.tsv.gz")
    if not matrix.exists():
        missing.append(TENX_MATRIX)
    if missing:
        raise ValueError(
            f"Incomplete 10X matrix triplet in {directory}: missing "
            + ", ".join(missing)
        )

    assert features is not None
    for path in (barcodes, features, matrix):
        _validate_gzip(path)

    n_barcodes = _count_gzip_lines(barcodes)
    n_features = _count_gzip_lines(features)
    mtx_rows, mtx_cols, n_nonzero = _read_matrix_market_shape(matrix)
    if mtx_rows != n_features or mtx_cols != n_barcodes:
        raise ValueError(
            "10X matrix dimensions do not match feature/barcode counts: "
            f"matrix={mtx_rows}x{mtx_cols}, "
            f"features={n_features}, barcodes={n_barcodes}"
        )

    return TenXTriplet(
        directory=directory,
        barcodes=barcodes,
        features=features,
        matrix=matrix,
        sample_id=_sample_id_for_10x_dir(directory),
        n_features=n_features,
        n_barcodes=n_barcodes,
        n_nonzero=n_nonzero,
    )


def ingest_10x_mtx_triplet(
    directory: str | Path,
    output_dir: str | Path,
    sample_id: str | None = None,
) -> dict[str, Any]:
    """Convert a validated 10X matrix directory into a canonical h5ad."""
    triplet = validate_10x_mtx_triplet(directory)
    sample = _safe_sample_id(sample_id or triplet.sample_id)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / f"{sample}.h5ad"

    params = {
        "mode": "10x_mtx",
        "directory": str(triplet.directory),
        "sample_id": sample,
        "reader": "scanpy.read_10x_mtx",
        "var_names": "gene_symbols",
        "make_unique": True,
    }

    os.environ.setdefault(
        "NUMBA_CACHE_DIR",
        "/tmp/aria_numba_cache",
    )
    Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

    import anndata as ad  # noqa: F401 - version captured by report tooling
    import scanpy as sc

    adata = sc.read_10x_mtx(
        str(triplet.directory),
        var_names="gene_symbols",
        make_unique=True,
        cache=False,
    )
    adata.uns.setdefault("aria_ingestion", {})
    adata.uns["aria_ingestion"] = {
        "mode": "10x_mtx",
        "source_directory": str(triplet.directory),
        "sample_id": sample,
        "source_files_json": json.dumps([
            {
                "path": str(path),
                "sha256": hash_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in triplet.source_paths()
        ], sort_keys=True),
        "params_sha256": hash_params(params),
    }
    adata.write_h5ad(output_path)

    source_files = [
        {
            "path": str(path),
            "sha256": hash_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in triplet.source_paths()
    ]
    return {
        "status": "success",
        "mode": "10x_mtx",
        "sample_id": sample,
        "source_directory": str(triplet.directory),
        "source_files": source_files,
        "output_h5ad": str(output_path),
        "output_sha256": hash_file(output_path),
        "n_features": triplet.n_features,
        "n_barcodes": triplet.n_barcodes,
        "n_nonzero": triplet.n_nonzero,
        "params": params,
        "params_sha256": hash_params(params),
    }


def scan_fastq_plan(root: str | Path) -> dict[str, Any]:
    """
    Detect FASTQ groups and return a non-executing ingestion plan.

    This deliberately does not infer chemistry, references, or indices.
    Missing requirements are reported as blockers for a future explicit
    checkpoint.
    """
    root = Path(root)
    fastqs = sorted(
        p for p in root.rglob("*")
        if p.is_file() and _is_fastq_name(p.name)
    )
    groups: dict[str, dict[str, Any]] = {}
    warnings: list[str] = []

    for path in fastqs:
        try:
            if path.name.endswith(".gz"):
                _validate_gzip(path)
        except ValueError as exc:
            warnings.append(str(exc))
            continue
        parsed = _parse_fastq_name(path.name)
        if not parsed:
            warnings.append(f"Could not parse FASTQ read role: {path}")
            continue
        sample = parsed["sample"]
        role = parsed["role"]
        lane = parsed.get("lane")
        rec = groups.setdefault(sample, {"sample_id": sample, "files": {}})
        rec["files"].setdefault(role, []).append(str(path))
        if lane:
            rec.setdefault("lanes", set()).add(lane)

    normalized_groups = []
    blockers = []
    for sample, rec in sorted(groups.items()):
        files = rec["files"]
        if not files.get("R1") or not files.get("R2"):
            blockers.append(
                f"FASTQ sample '{sample}' is missing an R1/R2 pair."
            )
        normalized_groups.append({
            "sample_id": sample,
            "files": {k: sorted(v) for k, v in sorted(files.items())},
            "lanes": sorted(rec.get("lanes", set())),
        })

    if fastqs:
        blockers.extend([
            "FASTQ ingestion requires explicit chemistry.",
            "FASTQ ingestion requires organism/genome or transcriptome.",
            "FASTQ ingestion requires a local kb index path and SHA-256.",
            "FASTQ ingestion requires kb/kallisto/bustools availability.",
        ])

    return {
        "status": "blocked" if blockers else "no_fastq_found",
        "mode": "fastq_kb_plan",
        "fastq_count": len(fastqs),
        "samples": normalized_groups,
        "warnings": warnings,
        "blockers": blockers,
        "tool_versions": _kb_tool_versions(),
        "params_sha256": hash_params({
            "mode": "fastq_kb_plan",
            "root": str(root),
            "fastq_count": len(fastqs),
            "samples": normalized_groups,
        }),
    }


def build_kb_count_command(
    fastq_files: list[str],
    index_path: str,
    output_dir: str,
    chemistry: str,
    threads: int = 8,
    technology: str | None = None,
) -> list[str]:
    """Build an explicit kb count command without executing it."""
    tech = technology or chemistry
    if not fastq_files:
        raise ValueError("kb count requires FASTQ files")
    if not index_path:
        raise ValueError("kb count requires a local index path")
    if not chemistry:
        raise ValueError("kb count requires explicit chemistry")
    return [
        "kb", "count",
        "-i", str(index_path),
        "-x", str(tech),
        "-o", str(output_dir),
        "-t", str(int(threads)),
        *map(str, fastq_files),
    ]


def _validate_gzip(path: Path) -> None:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
    except Exception as exc:
        raise ValueError(f"Invalid gzip file {path}: {exc}") from exc


def _count_gzip_lines(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return sum(1 for _ in fh)


def _read_matrix_market_shape(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().strip()
        if not header.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Invalid Matrix Market header in {path}")
        for line in fh:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid Matrix Market shape line in {path}")
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"Missing Matrix Market shape line in {path}")


def _sample_id_for_10x_dir(directory: Path) -> str:
    if directory.name in {"filtered_feature_bc_matrix", "raw_feature_bc_matrix"}:
        return _safe_sample_id(directory.parent.name)
    return _safe_sample_id(directory.name)


def _safe_sample_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return cleaned.strip("._-") or "sample"


def _is_fastq_name(name: str) -> bool:
    lower = name.lower()
    return lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def _parse_fastq_name(name: str) -> dict[str, str] | None:
    base = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name, flags=re.IGNORECASE)
    match = re.match(
        r"(?P<sample>.+?)(?:_S\d+)?(?:_(?P<lane>L\d{3}))?"
        r"_(?P<role>R[12]|I[12])(?:_\d{3})?$",
        base,
        flags=re.IGNORECASE,
    )
    if not match:
        match = re.match(
            r"(?P<sample>.+?)_(?P<role>[12])$",
            base,
            flags=re.IGNORECASE,
        )
    if not match:
        return None
    role = match.group("role").upper()
    if role == "1":
        role = "R1"
    elif role == "2":
        role = "R2"
    return {
        "sample": _safe_sample_id(match.group("sample")),
        "lane": match.groupdict().get("lane") or "",
        "role": role,
    }


def _kb_tool_versions() -> dict[str, str]:
    versions = {}
    for tool in ("kb", "kallisto", "bustools"):
        executable = shutil.which(tool)
        if not executable:
            versions[tool] = "not installed"
            continue
        try:
            result = subprocess.run(
                [tool, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
            text = (result.stdout or result.stderr or "").strip().splitlines()
            versions[tool] = text[0] if text else "unknown"
        except Exception as exc:
            versions[tool] = f"unavailable: {exc}"
    return versions
