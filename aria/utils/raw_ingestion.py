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
from typing import Any, Callable

from aria.utils.provenance import hash_file, hash_params


TENX_BARCODE = "barcodes.tsv.gz"
TENX_MATRIX = "matrix.mtx.gz"
TENX_FEATURES = ("features.tsv.gz", "genes.tsv.gz")
ProgressCallback = Callable[[str, dict[str, Any]], None]


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


def discover_10x_mtx_triplets(
    root: str | Path,
    progress_callback: ProgressCallback | None = None,
) -> list[TenXTriplet]:
    """Find valid Cell Ranger-style 10X matrix triplets under *root*."""
    root = Path(root)
    if root.is_file():
        root = root.parent
    candidates: set[Path] = set()
    files_seen = 0
    dirs_seen = 0
    last_report = 0
    for path in _walk_files(root):
        files_seen += 1
        if path.name == TENX_MATRIX:
            candidates.add(path.parent)
        parent_count = len(candidates)
        if progress_callback and files_seen - last_report >= 1000:
            last_report = files_seen
            progress_callback(
                "scan_10x_mtx",
                {
                    "files_seen": files_seen,
                    "candidate_dirs": parent_count,
                },
            )
    if (root / TENX_MATRIX).exists():
        candidates.add(root)
    candidate_count = len(candidates)
    if progress_callback:
        progress_callback(
            "scan_10x_mtx",
            {"files_seen": files_seen, "candidate_dirs": candidate_count},
        )

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
            if progress_callback:
                progress_callback(
                    "validate_10x_mtx",
                    {
                        "validated_triplets": len(triplets),
                        "candidate_dirs": candidate_count,
                    },
                )
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
    barcode_summary = _inspect_10x_barcodes(barcodes)
    if not _is_likely_10x_mtx_dir(directory, barcode_summary):
        raise ValueError(
            "10X matrix triplet lacks 10X-like cell barcode evidence: "
            f"{directory}"
        )
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


def scan_fastq_plan(
    root: str | Path,
    progress_callback: ProgressCallback | None = None,
) -> dict[str, Any]:
    """
    Detect FASTQ groups and return a non-executing ingestion plan.

    This deliberately does not infer chemistry, references, or indices.
    Missing requirements are reported as blockers for a future explicit
    checkpoint.
    """
    root = Path(root)
    fastqs = []
    files_seen = 0
    last_report = 0
    for path in _walk_files(root):
        files_seen += 1
        if _is_fastq_name(path.name):
            fastqs.append(path)
        if progress_callback and files_seen - last_report >= 1000:
            last_report = files_seen
            progress_callback(
                "scan_fastq",
                {"files_seen": files_seen, "fastq_count": len(fastqs)},
            )
    fastqs = sorted(fastqs)
    if progress_callback:
        progress_callback(
            "scan_fastq",
            {"files_seen": files_seen, "fastq_count": len(fastqs)},
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


def _walk_files(root: Path):
    if not root.exists():
        return
    if root.is_file():
        yield root
        return
    stack = [root]
    while stack:
        current = stack.pop()
        try:
            with os.scandir(current) as entries:
                for entry in entries:
                    path = Path(entry.path)
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            stack.append(path)
                        elif entry.is_file(follow_symlinks=False):
                            yield path
                    except OSError:
                        continue
        except OSError:
            continue


def build_kb_count_command(
    fastq_files: list[str],
    index_path: str,
    t2g_path: str,
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
    if not t2g_path:
        raise ValueError("kb count requires a local transcript-to-gene path")
    if not chemistry:
        raise ValueError("kb count requires explicit chemistry")
    return [
        "kb", "count",
        "-i", str(index_path),
        "-g", str(t2g_path),
        "-x", str(tech),
        "-o", str(output_dir),
        "-t", str(int(threads)),
        *map(str, fastq_files),
    ]


def execute_kb_count(params: dict[str, Any]) -> dict[str, Any]:
    """
    Run kb count only when every deterministic input is explicit.

    Required params:
      fastq_files, index_path, index_sha256, t2g_path, t2g_sha256,
      chemistry, output_dir
    """
    blockers = []
    fastq_files = [str(p) for p in params.get("fastq_files", []) or []]
    index_path = str(params.get("index_path") or "")
    t2g_path = str(params.get("t2g_path") or "")
    output_dir = str(params.get("output_dir") or "")
    chemistry = str(params.get("chemistry") or "")
    threads = int(params.get("threads") or 8)
    expected_index_hash = str(params.get("index_sha256") or "")
    expected_t2g_hash = str(params.get("t2g_sha256") or "")

    if not fastq_files:
        blockers.append("kb execution requires FASTQ files.")
    for fq in fastq_files:
        if not Path(fq).exists():
            blockers.append(f"FASTQ not found: {fq}")
    if not index_path or not Path(index_path).exists():
        blockers.append("kb execution requires an existing local index_path.")
    if not expected_index_hash:
        blockers.append("kb execution requires index_sha256.")
    if not t2g_path or not Path(t2g_path).exists():
        blockers.append("kb execution requires an existing local t2g_path.")
    if not expected_t2g_hash:
        blockers.append("kb execution requires t2g_sha256.")
    if not chemistry:
        blockers.append("kb execution requires explicit chemistry.")
    if not output_dir:
        blockers.append("kb execution requires output_dir.")
    if shutil.which("kb") is None:
        blockers.append("kb executable is not installed.")

    if not blockers:
        try:
            if hash_file(index_path) != expected_index_hash:
                blockers.append("index_sha256 does not match index_path.")
        except Exception as exc:
            blockers.append(f"Could not hash index_path: {exc}")
    if not blockers:
        try:
            if hash_file(t2g_path) != expected_t2g_hash:
                blockers.append("t2g_sha256 does not match t2g_path.")
        except Exception as exc:
            blockers.append(f"Could not hash t2g_path: {exc}")

    command = []
    if not blockers:
        command = build_kb_count_command(
            fastq_files=fastq_files,
            index_path=index_path,
            t2g_path=t2g_path,
            output_dir=output_dir,
            chemistry=chemistry,
            threads=threads,
            technology=params.get("technology"),
        )
    if blockers:
        return {
            "status": "blocked",
            "mode": "fastq_kb_count",
            "blockers": blockers,
            "tool_versions": _kb_tool_versions(),
            "params_sha256": hash_params(params),
        }

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = out_dir / "kb_count.stdout.txt"
    stderr_path = out_dir / "kb_count.stderr.txt"
    with stdout_path.open("w", encoding="utf-8") as stdout, \
            stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(
            command,
            stdout=stdout,
            stderr=stderr,
            text=True,
            timeout=int(params.get("timeout_seconds") or 86400),
            check=False,
        )
    if result.returncode != 0:
        return {
            "status": "failed",
            "mode": "fastq_kb_count",
            "command": command,
            "returncode": result.returncode,
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
            "tool_versions": _kb_tool_versions(),
            "params_sha256": hash_params(params),
        }

    matrix_dir = _find_kb_matrix_dir(out_dir)
    if matrix_dir is None:
        return {
            "status": "failed",
            "mode": "fastq_kb_count",
            "command": command,
            "returncode": result.returncode,
            "details": "kb completed but no 10X-style matrix directory was found.",
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
            "tool_versions": _kb_tool_versions(),
            "params_sha256": hash_params(params),
        }
    sample_id = _safe_sample_id(params.get("sample_id") or "kb_count")
    h5ad_result = ingest_10x_mtx_triplet(
        matrix_dir,
        out_dir / "h5ad",
        sample_id=sample_id,
    )
    return {
        "status": "success",
        "mode": "fastq_kb_count",
        "command": command,
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "matrix_dir": str(matrix_dir),
        "output_h5ad": h5ad_result["output_h5ad"],
        "output_sha256": h5ad_result["output_sha256"],
        "source_files": [
            {"path": fq, "sha256": hash_file(fq), "size_bytes": Path(fq).stat().st_size}
            for fq in fastq_files
        ],
        "reference": {
            "index_path": index_path,
            "index_sha256": expected_index_hash,
            "t2g_path": t2g_path,
            "t2g_sha256": expected_t2g_hash,
            "chemistry": chemistry,
        },
        "tool_versions": _kb_tool_versions(),
        "params_sha256": hash_params(params),
    }


def _validate_gzip(path: Path) -> None:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
    except Exception as exc:
        raise ValueError(f"Invalid gzip file {path}: {exc}") from exc


def _inspect_10x_barcodes(path: Path) -> dict[str, Any]:
    checked = 0
    tenx_like = 0
    examples: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            value = line.strip().split("\t", 1)[0]
            if not value:
                continue
            checked += 1
            if len(examples) < 3:
                examples.append(value)
            if _looks_like_10x_cell_barcode(value):
                tenx_like += 1
            if checked >= 64:
                break
    return {
        "barcode_rows_checked": checked,
        "tenx_like_barcode_fraction": (
            round(tenx_like / checked, 3) if checked else 0.0
        ),
        "barcode_examples": examples,
    }


def _is_likely_10x_mtx_dir(directory: Path, barcode_summary: dict[str, Any]) -> bool:
    if directory.name.lower() in {"filtered_feature_bc_matrix", "raw_feature_bc_matrix"}:
        return True
    checked = int(barcode_summary.get("barcode_rows_checked") or 0)
    fraction = float(barcode_summary.get("tenx_like_barcode_fraction") or 0.0)
    return checked > 0 and fraction >= 0.80


def _looks_like_10x_cell_barcode(value: str) -> bool:
    return re.fullmatch(r"[ACGTNacgtn]{8,32}(?:-\d+)?", value.strip()) is not None


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


def _find_kb_matrix_dir(output_dir: Path) -> Path | None:
    candidates = [
        output_dir / "counts_filtered",
        output_dir / "counts_unfiltered",
        output_dir,
    ]
    candidates.extend(p for p in output_dir.rglob("*") if p.is_dir())
    for candidate in candidates:
        if (candidate / TENX_MATRIX).exists() and (candidate / TENX_BARCODE).exists():
            if any((candidate / name).exists() for name in TENX_FEATURES):
                return candidate
    return None


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
