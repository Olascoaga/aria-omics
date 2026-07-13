"""Pre-dispatch assay contract validation.

These checks are deliberately data-light. They validate structural assumptions
that would otherwise fail late inside expensive modality scripts: assembly,
sample metadata, barcode namespace, replicate structure, and coordinate-feature
shape. Optional downstream resources such as motif collections or GTFs remain
script-gated honest skips.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any


_CHROMATIN_MODALITIES = {"scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"}
_FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
_COORD_SUFFIXES = (
    ".bed",
    ".bed.gz",
    ".narrowpeak",
    ".narrowpeak.gz",
    ".broadpeak",
    ".broadpeak.gz",
    ".fragments.tsv",
    ".fragments.tsv.gz",
    "fragments.tsv",
    "fragments.tsv.gz",
)
_UNKNOWN = {"", "unknown", "none", "null", "na", "n/a"}


def validate_assay_contract(
    exp_context: dict[str, Any],
    modality: str,
    files: list[str] | tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Validate a modality's pre-dispatch assay contract.

    Returns a readiness-card-shaped update:
    ``{"status": "green|yellow|red", "checks": {"assay_contract": ...},
    "findings": [...]}``.
    """
    files = list(files or (exp_context.get("modalities") or {}).get(modality) or [])
    findings: list[dict[str, Any]] = []
    checks: dict[str, Any] = {
        "modality": modality,
        "genome": _genome_check(exp_context, modality),
        "sample_metadata": _sample_metadata_check(exp_context, files),
        "barcode_namespace": _barcode_check(exp_context, modality, files),
        "replicate_structure": _replicate_check(exp_context, modality, files),
        "feature_coordinates": _coordinate_check(files),
    }

    if modality in _CHROMATIN_MODALITIES and not checks["genome"]["valid"]:
        findings.append(_finding(
            "blocking",
            "assay_contract_genome_missing",
            (
                f"{modality} dispatch needs a confirmed genome/assembly before "
                "reference-dependent chromatin analysis."
            ),
            (
                "Confirm the genome assembly at the design checkpoint or provide "
                "`genome`/`assembly` in the run context before dispatch."
            ),
            modality,
        ))

    barcode = checks["barcode_namespace"]
    if barcode["status"] == "blocked":
        errors = barcode.get("errors") or []
        findings.append(_finding(
            "blocking",
            "assay_contract_scatac_fastq_barcode_missing",
            (
                "scATAC FASTQ input needs a valid typed library manifest with "
                "explicit R1/R2/R3 roles, barcode whitelist, and donor/sample "
                "identity. " + ("; ".join(errors) if errors else
                "The barcode read cannot be inferred from genomic reads.")
            ),
            (
                "Provide/edit `scatac_fastq_manifest` (schema_version 1), or "
                "for one legacy library provide `barcode_fastq` and "
                "`barcode_whitelist`; alternatively start from fragments or a "
                "peak matrix."
            ),
            modality,
        ))

    coords = checks["feature_coordinates"]
    if coords["status"] == "blocked":
        findings.append(_finding(
            "blocking",
            "assay_contract_coordinate_features_invalid",
            (
                f"{modality} coordinate features could not be parsed as "
                "chrom/start/end BED-like intervals."
            ),
            (
                "Use BED/narrowPeak/fragments-style 0-based half-open "
                "coordinates with numeric start/end columns, and ensure contig "
                "names match the selected genome."
            ),
            modality,
        ))

    metadata = checks["sample_metadata"]
    if metadata["status"] == "warning":
        findings.append(_finding(
            "warning",
            "assay_contract_sample_metadata_incomplete",
            (
                f"{modality} has incomplete sample metadata for explicit "
                "condition/replicate-aware analysis."
            ),
            (
                "Provide sample metadata keyed by sample id with condition and "
                "biological replicate columns. Without it, inferential DA/DE "
                "layers will skip honestly."
            ),
            modality,
        ))

    reps = checks["replicate_structure"]
    if reps["status"] == "warning":
        findings.append(_finding(
            "warning",
            "assay_contract_replicates_low_or_missing",
            (
                f"{modality} replicate structure is incomplete or low-powered: "
                f"{reps.get('replicates_per_condition') or 'not available'}."
            ),
            (
                "Use at least two biological replicates per condition for a "
                "supported inferential run; three or more is preferred for "
                "publication-grade designs."
            ),
            modality,
        ))

    status = (
        "red" if any(f["severity"] == "blocking" for f in findings)
        else "yellow" if any(f["severity"] == "warning" for f in findings)
        else "green"
    )
    checks["status"] = status
    return {
        "status": status,
        "checks": {"assay_contract": checks},
        "findings": findings,
    }


def _finding(
    severity: str,
    check: str,
    message: str,
    recommendation: str,
    modality: str,
) -> dict[str, Any]:
    return {
        "severity": severity,
        "check": check,
        "message": message,
        "recommendation": recommendation,
        "modality": modality,
    }


def _genome_check(exp_context: dict[str, Any], modality: str) -> dict[str, Any]:
    genome = (
        exp_context.get("genome")
        or exp_context.get("assembly")
        or (exp_context.get("design") or {}).get("genome")
        or (exp_context.get("inferred_design") or {}).get("genome")
    )
    organism = (
        exp_context.get("organism")
        or (exp_context.get("design") or {}).get("organism")
        or (exp_context.get("inferred_design") or {}).get("organism")
    )
    genome_text = str(genome or "").strip()
    return {
        "genome": genome_text or None,
        "organism": str(organism).strip() if organism else None,
        "valid": genome_text.lower() not in _UNKNOWN,
        "required": modality in _CHROMATIN_MODALITIES,
    }


def _barcode_check(
    exp_context: dict[str, Any],
    modality: str,
    files: list[str],
) -> dict[str, Any]:
    fastq = any(_is_fastq(p) for p in files)
    if modality != "scATAC" or not fastq:
        return {"status": "pass", "required": False, "fastq_input": fastq}
    from aria.utils.scatac_fastq_manifest import resolve_scatac_fastq_manifest

    manifest_result = exp_context.get("scatac_fastq_manifest_validation")
    if not isinstance(manifest_result, dict) or manifest_result.get("status") == "missing":
        manifest_result = resolve_scatac_fastq_manifest(
            exp_context, require_paths=False
        )
    if manifest_result.get("status") == "valid":
        manifest = manifest_result["manifest"]
        return {
            "status": "pass",
            "required": True,
            "fastq_input": True,
            "contract": "scatac_fastq_manifest/v1",
            "n_libraries": len(manifest["libraries"]),
            "libraries": [
                {
                    "library_id": row["library_id"],
                    "sample_id": row["sample_id"],
                    "donor_id": row["donor_id"],
                    "roles": list(row["fastqs"]),
                }
                for row in manifest["libraries"]
            ],
        }
    if manifest_result.get("status") == "invalid":
        return {
            "status": "blocked",
            "required": True,
            "fastq_input": True,
            "contract": "scatac_fastq_manifest/v1",
            "errors": list(manifest_result.get("errors") or []),
        }

    barcode = exp_context.get("barcode_fastq") or exp_context.get("cell_barcode_fastq")
    whitelist = (
        exp_context.get("barcode_whitelist")
        or exp_context.get("cell_barcode_whitelist")
    )
    if barcode and whitelist:
        return {
            "status": "pass",
            "required": True,
            "fastq_input": True,
            "barcode_fastq": str(barcode),
            "barcode_whitelist": str(whitelist),
        }
    return {
        "status": "blocked",
        "required": True,
        "fastq_input": True,
        "missing": [
            name for name, value in (
                ("barcode_fastq", barcode),
                ("barcode_whitelist", whitelist),
            )
            if not value
        ],
    }


def _sample_metadata_check(exp_context: dict[str, Any], files: list[str]) -> dict[str, Any]:
    comparisons = exp_context.get("comparisons")
    if not comparisons:
        return {"status": "pass", "required_for_comparison": False}

    metadata = exp_context.get("sample_metadata") or {}
    sample_ids = _sample_ids(exp_context, files)
    condition_col = exp_context.get("condition_col", "condition")
    replicate_col = exp_context.get("replicate_col") or exp_context.get(
        "replicate_column", "replicate"
    )
    if not isinstance(metadata, dict) or not metadata:
        return {
            "status": "warning",
            "required_for_comparison": True,
            "sample_ids": sample_ids,
            "reason": "missing_sample_metadata",
        }

    missing = []
    for sample_id in sample_ids:
        row = metadata.get(sample_id) or {}
        if not isinstance(row, dict):
            missing.append(sample_id)
            continue
        if not row.get(condition_col) or not row.get(replicate_col):
            missing.append(sample_id)

    return {
        "status": "warning" if missing else "pass",
        "required_for_comparison": True,
        "sample_ids": sample_ids,
        "condition_col": condition_col,
        "replicate_col": replicate_col,
        "missing_or_incomplete_samples": missing,
    }


def _replicate_check(
    exp_context: dict[str, Any],
    modality: str,
    files: list[str],
) -> dict[str, Any]:
    comparisons = exp_context.get("comparisons")
    if not comparisons or modality not in {"bulk_ATAC", "scATAC", "bulk_RNA", "scRNA"}:
        return {"status": "pass", "required_for_comparison": bool(comparisons)}

    metadata = exp_context.get("sample_metadata") or {}
    condition_col = exp_context.get("condition_col", "condition")
    replicate_col = exp_context.get("replicate_col") or exp_context.get(
        "replicate_column", "replicate"
    )
    counts: dict[str, set[str]] = {}
    if isinstance(metadata, dict) and metadata:
        for sample_id in _sample_ids(exp_context, files):
            row = metadata.get(sample_id) or {}
            if not isinstance(row, dict):
                continue
            condition = row.get(condition_col) or row.get("condition")
            replicate = row.get(replicate_col) or row.get("replicate")
            if condition and replicate:
                counts.setdefault(str(condition), set()).add(str(replicate))
    else:
        design = exp_context.get("design") or exp_context.get("inferred_design") or {}
        groups = design.get("analysis_groups") or design.get("groups") or {}
        for condition, reps in groups.items():
            if isinstance(reps, (list, tuple, set)):
                counts[str(condition)] = {str(r) for r in reps if str(r).strip()}
            elif reps:
                counts[str(condition)] = {str(reps)}

    per_condition = {condition: len(reps) for condition, reps in counts.items()}
    min_reps = min(per_condition.values()) if per_condition else 0
    status = "warning" if not per_condition or min_reps < 2 else "pass"
    return {
        "status": status,
        "required_for_comparison": True,
        "replicates_per_condition": per_condition,
        "minimum_supported": 2,
        "publication_ready_min": 3,
    }


def _coordinate_check(files: list[str]) -> dict[str, Any]:
    coord_files = [p for p in files if _is_coordinate_like(p)]
    if not coord_files:
        return {"status": "pass", "coordinate_files": []}

    inspected = []
    invalid = []
    contig_styles = set()
    for path in coord_files:
        parsed = _inspect_coordinate_file(path)
        inspected.append(parsed)
        if parsed["status"] == "invalid":
            invalid.append(path)
        style = parsed.get("contig_style")
        if style:
            contig_styles.add(style)

    status = "blocked" if invalid else "pass"
    return {
        "status": status,
        "coordinate_files": coord_files,
        "inspected": inspected,
        "contig_style": (
            next(iter(contig_styles)) if len(contig_styles) == 1
            else "mixed" if contig_styles else "unchecked"
        ),
    }


def _inspect_coordinate_file(path: str) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        return {"path": path, "status": "unchecked_missing_path"}
    try:
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or stripped.startswith("track"):
                    continue
                parts = stripped.split("\t")
                if len(parts) < 3:
                    return {"path": path, "status": "invalid", "reason": "too_few_columns"}
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                if start < 0 or end <= start:
                    return {
                        "path": path,
                        "status": "invalid",
                        "reason": "invalid_interval",
                    }
                return {
                    "path": path,
                    "status": "valid",
                    "coordinate_convention": "bed_0_based_half_open",
                    "contig_style": _contig_style(chrom),
                }
        return {"path": path, "status": "invalid", "reason": "no_interval_rows"}
    except Exception as exc:
        return {"path": path, "status": "invalid", "reason": type(exc).__name__}


def _sample_ids(exp_context: dict[str, Any], files: list[str]) -> list[str]:
    sample_ids = exp_context.get("sample_ids")
    if isinstance(sample_ids, (list, tuple)) and sample_ids:
        return [str(s) for s in sample_ids]
    return [Path(f).stem.replace(".bam", "") for f in files]


def _is_fastq(path: str) -> bool:
    lower = str(path).lower()
    return lower.endswith(_FASTQ_SUFFIXES)


def _is_coordinate_like(path: str) -> bool:
    lower = str(path).lower()
    return any(lower.endswith(suffix) for suffix in _COORD_SUFFIXES)


def _contig_style(chrom: str) -> str:
    if chrom.startswith("chr"):
        return "ucsc_chr"
    return "ensembl_no_chr"
