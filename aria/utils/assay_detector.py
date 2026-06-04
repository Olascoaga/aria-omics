"""Content-based assay detection for DataAudit.

The detector is conservative: it only returns a modality when the file's
internal structure is recognizable. Filename signatures remain a fallback in
DataAudit for formats whose content alone cannot identify the assay.
"""

from __future__ import annotations

import csv
import gzip
import math
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"
TEXT_SNIFF_BYTES = 65536


@dataclass(frozen=True)
class AssayDetection:
    modality: str
    confidence: str
    reason: str
    evidence: dict[str, object] = field(default_factory=dict)
    possible_alternatives: tuple[str, ...] = ()
    blocking_issues: tuple[str, ...] = ()


class AssayDetector:
    """Detect supported assays from lightweight file content inspection."""

    def detect_file(self, path: str | Path) -> AssayDetection | None:
        p = Path(path)
        if not p.is_file():
            return None

        try:
            prefix = _read_prefix(p, len(HDF5_MAGIC))
        except OSError:
            return None

        if prefix == HDF5_MAGIC:
            return self._detect_hdf5(p)
        alignment = self._detect_alignment(p, prefix)
        if alignment is not None:
            return alignment

        text = _read_text_prefix(p)
        if not text:
            return None
        matrix = self._detect_matrix_market(p, text)
        if matrix is not None:
            return matrix
        return self._detect_count_table(text)

    def is_supported_file(self, path: str | Path) -> bool:
        return self.detect_file(path) is not None

    def _detect_hdf5(self, path: Path) -> AssayDetection | None:
        try:
            import h5py
        except ImportError:
            return None

        try:
            with h5py.File(path, "r") as h5:
                keys = {str(k).lower() for k in h5.keys()}
                if {"chroms", "bins", "pixels"}.issubset(keys):
                    return AssayDetection(
                        modality="HiC",
                        confidence="high",
                        reason="HDF5 contains cooler chroms/bins/pixels groups.",
                        evidence={"format": "cooler_hdf5"},
                    )
                if "mod" in keys:
                    mod = h5["mod"]
                    if hasattr(mod, "keys"):
                        mod_names = {str(k).lower() for k in mod.keys()}
                        has_atac = _contains_any(
                            mod_names, ("atac", "peaks", "chromatin", "accessibility")
                        )
                        has_rna = _contains_any(
                            mod_names, ("rna", "gex", "expression")
                        )
                        if has_atac:
                            return AssayDetection(
                                modality="scATAC",
                                confidence="high" if has_rna else "medium",
                                reason="HDF5 MuData mod group contains an ATAC-like modality.",
                                evidence={
                                    "format": "h5mu",
                                    "modalities": sorted(mod_names),
                                    "paired_rna": has_rna,
                                },
                                possible_alternatives=(("scRNA",) if has_rna else ()),
                            )
                if {"obs", "var"}.issubset(keys) and (
                    "x" in keys or "layers" in keys or "raw" in keys
                ):
                    return AssayDetection(
                        modality="scRNA",
                        confidence="high",
                        reason="HDF5 contains AnnData obs/var/X-like structure.",
                        evidence={"format": "h5ad"},
                    )
                if "matrix" in keys:
                    matrix = h5["matrix"]
                    if hasattr(matrix, "keys"):
                        matrix_keys = {str(k).lower() for k in matrix.keys()}
                        if {"barcodes", "data"}.issubset(matrix_keys) and (
                            "features" in matrix_keys or "genes" in matrix_keys
                        ):
                            tenx = _inspect_10x_h5_matrix(matrix)
                            if tenx["has_peaks"]:
                                return AssayDetection(
                                    modality="scATAC",
                                    confidence="high",
                                    reason=(
                                        "HDF5 contains a 10X matrix with ATAC peak "
                                        "features."
                                    ),
                                    evidence={
                                        "format": (
                                            "10x_multiome_h5"
                                            if tenx["has_gene_expression"]
                                            else "10x_atac_h5"
                                        ),
                                        **tenx,
                                    },
                                    possible_alternatives=(
                                        ("scRNA",) if tenx["has_gene_expression"]
                                        else ()
                                    ),
                                )
                            return AssayDetection(
                                modality="scRNA",
                                confidence="high",
                                reason="HDF5 contains Cell Ranger 10X matrix group.",
                                evidence={"format": "10x_h5", **tenx},
                            )
        except OSError:
            return None
        return None

    def _detect_matrix_market(
        self, path: Path, text: str
    ) -> AssayDetection | None:
        first = _first_nonempty_line(text)
        if first is None or not first.startswith("%%MatrixMarket"):
            return None
        try:
            siblings = {p.name.lower(): p for p in path.parent.iterdir()}
        except OSError:
            return None
        has_barcodes = bool({"barcodes.tsv", "barcodes.tsv.gz"} & siblings.keys())
        feature_name = next(
            (
                name for name in (
                    "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"
                )
                if name in siblings
            ),
            None,
        )
        has_features = feature_name is not None
        if has_barcodes and has_features:
            feature_summary = _inspect_10x_feature_table(
                siblings[feature_name] if feature_name else None
            )
            if feature_summary["has_peaks"]:
                return AssayDetection(
                    modality="scATAC",
                    confidence="high",
                    reason="MatrixMarket file is part of a 10X ATAC/Multiome MEX triplet.",
                    evidence={"format": "10x_mtx", **feature_summary},
                    possible_alternatives=(
                        ("scRNA",) if feature_summary["has_gene_expression"] else ()
                    ),
                )
            return AssayDetection(
                modality="scRNA",
                confidence="high",
                reason="MatrixMarket file is part of a 10X MEX triplet.",
                evidence={"format": "10x_mtx", **feature_summary},
            )
        return None

    def _detect_count_table(self, text: str) -> AssayDetection | None:
        lines = [
            line for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ][:8]
        if len(lines) < 3:
            return None

        delimiter = _pick_delimiter(lines[0])
        if delimiter is None:
            return None
        try:
            rows = list(csv.reader(lines, delimiter=delimiter))
        except csv.Error:
            return None
        if len(rows) < 3 or len(rows[0]) < 3:
            return None

        header = [c.strip().lower() for c in rows[0]]
        if _is_salmon_quant_header(header):
            return AssayDetection(
                modality="bulk_RNA",
                confidence="high",
                reason="Table header matches Salmon quant.sf transcript quantification.",
                evidence={"format": "salmon_quant", "columns": header},
            )
        if _is_kallisto_abundance_header(header):
            return AssayDetection(
                modality="bulk_RNA",
                confidence="high",
                reason="Table header matches Kallisto abundance.tsv quantification.",
                evidence={"format": "kallisto_abundance", "columns": header},
            )

        first_col = header[0]
        if first_col not in {
            "gene", "gene_id", "geneid", "genes", "symbol", "feature",
            "feature_id", "target_id",
        }:
            return None

        numeric_rows = 0
        integer_values = 0
        nonnegative_values = 0
        zero_values = 0
        total_values = 0
        for row in rows[1:]:
            if len(row) < 3:
                continue
            values = [_as_float(v) for v in row[1:]]
            present = [v for v in values if v is not None and math.isfinite(v)]
            if len(present) >= max(2, len(values) - 1):
                numeric_rows += 1
            for value in present:
                total_values += 1
                if value >= 0:
                    nonnegative_values += 1
                if float(value).is_integer():
                    integer_values += 1
                if value == 0:
                    zero_values += 1

        if numeric_rows >= 2:
            integer_fraction = (
                integer_values / total_values if total_values else 0.0
            )
            return AssayDetection(
                modality="bulk_RNA",
                confidence="high" if integer_fraction >= 0.95 else "medium",
                reason="Delimited table has gene identifiers and numeric sample columns.",
                evidence={
                    "format": "count_table",
                    "delimiter": delimiter,
                    "numeric_rows_checked": numeric_rows,
                    "integer_fraction": round(integer_fraction, 3),
                    "nonnegative_fraction": (
                        round(nonnegative_values / total_values, 3)
                        if total_values else 0.0
                    ),
                    "sparsity": (
                        round(zero_values / total_values, 3)
                        if total_values else 0.0
                    ),
                },
            )
        return None

    def _detect_alignment(
        self, path: Path, prefix: bytes
    ) -> AssayDetection | None:
        header = ""
        if prefix.startswith(b"BAM\x01"):
            header = _read_bam_header_from_stream(path)
        elif prefix.startswith(b"\x1f\x8b"):
            decompressed = _read_gzip_prefix(path, 32768)
            # Only a gzipped BAM is alignment content. A gzipped FASTQ or count
            # table is NOT: scanning its decompressed bytes for assay keywords is
            # unsafe because read sequences contain "atac"/"star" 4-mers by chance
            # (this misclassified bulk RNA FASTQs as bulk_ATAC). Leave the header
            # empty so detection falls through to the filename/path signal.
            if decompressed.startswith("BAM\x01"):
                header = _parse_bam_header_text(
                    decompressed.encode("latin1", "ignore"))
        else:
            text = _read_text_prefix(path)
            if text.startswith("@HD") or text.startswith("@SQ"):
                header = text[:32768]

        if not header:
            return None
        lower = header.lower()
        if "star" in lower or "hisat" in lower or "rna-seq" in lower:
            return AssayDetection(
                modality="bulk_RNA",
                confidence="medium",
                reason="Alignment header indicates RNA-seq aligner/program metadata.",
                evidence={"format": "bam_or_sam", "program_hint": _program_hint(lower)},
                possible_alternatives=("bulk_ATAC", "ChIP"),
            )
        if "atac" in lower:
            return AssayDetection(
                modality="bulk_ATAC",
                confidence="medium",
                reason="Alignment header contains ATAC assay metadata.",
                evidence={"format": "bam_or_sam", "program_hint": _program_hint(lower)},
                possible_alternatives=("ChIP", "CUT_AND_RUN", "CUT_AND_TAG"),
            )
        return AssayDetection(
            modality="bulk_ATAC",
            confidence="low",
            reason=(
                "Alignment content was detected, but assay-specific header "
                "metadata was absent."
            ),
            evidence={"format": "bam_or_sam", "program_hint": _program_hint(lower)},
            possible_alternatives=("bulk_RNA", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"),
            blocking_issues=("assay_specific_alignment_metadata_missing",),
        )


def _read_prefix(path: Path, n: int) -> bytes:
    with path.open("rb") as fh:
        return fh.read(n)


def _read_text_prefix(path: Path) -> str:
    try:
        opener = gzip.open if "".join(path.suffixes).lower().endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read(TEXT_SNIFF_BYTES)
    except (OSError, UnicodeError):
        return ""


def _read_gzip_prefix(path: Path, n: int) -> str:
    try:
        with gzip.open(path, "rb") as fh:
            return fh.read(n).decode("latin1", "ignore")
    except OSError:
        return ""


def _read_bam_header_from_stream(path: Path) -> str:
    try:
        with path.open("rb") as fh:
            return _parse_bam_header_text(fh.read(32768))
    except OSError:
        return ""


def _parse_bam_header_text(raw: bytes) -> str:
    if not raw.startswith(b"BAM\x01") or len(raw) < 8:
        return ""
    try:
        header_len = struct.unpack("<i", raw[4:8])[0]
    except struct.error:
        return ""
    if header_len <= 0:
        return ""
    start = 8
    end = min(len(raw), start + header_len)
    return raw[start:end].decode("latin1", "ignore")


def _contains_any(values: Iterable[str], needles: Iterable[str]) -> bool:
    return any(any(needle in value for needle in needles) for value in values)


def _first_nonempty_line(text: str) -> str | None:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return None


def _pick_delimiter(header: str) -> str | None:
    counts = {delim: header.count(delim) for delim in ("\t", ",", ";")}
    delim, count = max(counts.items(), key=lambda kv: kv[1])
    return delim if count >= 2 else None


def _as_float(value: str) -> float | None:
    try:
        return float(value.strip())
    except ValueError:
        return None


def _decode_h5_values(dataset, limit: int = 200) -> list[str]:
    try:
        values = dataset[:limit]
    except Exception:
        return []
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "ignore"))
        else:
            decoded.append(str(value))
    return decoded


def _inspect_10x_h5_matrix(matrix) -> dict[str, object]:
    feature_types: list[str] = []
    genomes: list[str] = []
    try:
        features = matrix["features"]
        if hasattr(features, "keys"):
            if "feature_type" in features:
                feature_types = _decode_h5_values(features["feature_type"])
            if "genome" in features:
                genomes = _decode_h5_values(features["genome"])
    except Exception:
        pass
    return _feature_type_summary(feature_types, genomes)


def _inspect_10x_feature_table(path: Path | None) -> dict[str, object]:
    if path is None:
        return _feature_type_summary([], [])
    text = _read_text_prefix(path)
    feature_types = []
    genomes = []
    for line in text.splitlines()[:500]:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            feature_types.append(parts[2])
        if len(parts) >= 4:
            genomes.append(parts[3])
    return _feature_type_summary(feature_types, genomes)


def _feature_type_summary(
    feature_types: list[str], genomes: list[str]
) -> dict[str, object]:
    lowered = [ft.lower() for ft in feature_types]
    has_peaks = any("peak" in ft or "atac" in ft for ft in lowered)
    has_gene_expression = any(
        "gene expression" in ft or ft in {"gene", "genes", "expression"}
        for ft in lowered
    )
    return {
        "feature_types": sorted({ft for ft in feature_types if ft})[:12],
        "genomes": sorted({g for g in genomes if g})[:12],
        "has_peaks": has_peaks,
        "has_gene_expression": has_gene_expression,
    }


def _is_salmon_quant_header(header: list[str]) -> bool:
    return header[:5] == ["name", "length", "effectivelength", "tpm", "numreads"]


def _is_kallisto_abundance_header(header: list[str]) -> bool:
    return header[:5] == ["target_id", "length", "eff_length", "est_counts", "tpm"]


def _program_hint(header_lower: str) -> str | None:
    for hint in ("star", "hisat", "bowtie", "bwa", "minimap", "chromap"):
        if hint in header_lower:
            return hint
    return None
