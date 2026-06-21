"""Pre-integration multiome object contracts.

S10 validates each assay lane independently. This module validates the
cross-modal object that V48 integration will need: whether RNA and ATAC are
paired in the same cells, whether their barcode namespace is explicit, and
whether sample identities agree. It is intentionally data-light and does not run
integration algorithms.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


_RNA_MODALITIES = {"scRNA", "bulk_RNA", "bulk_RNA_raw"}
_ATAC_MODALITIES = {"scATAC", "bulk_ATAC"}
_H5MU_FORMATS = {"h5mu", "10x_multiome_h5"}


def infer_multiome_contract(exp_context: dict[str, Any]) -> dict[str, Any]:
    """Infer a conservative multiome contract from audited modalities.

    The result is an initial contract, not proof of integration validity. A
    paired `.h5mu` with RNA+ATAC evidence is treated as same-cell; split scRNA
    and scATAC inputs are deliberately marked as unconfirmed unless the caller
    supplies an explicit contract.
    """
    explicit = exp_context.get("multiome_contract")
    if isinstance(explicit, dict) and explicit:
        return dict(explicit)

    modalities = exp_context.get("modalities") or {}
    rna = sorted(m for m in modalities if m in _RNA_MODALITIES)
    atac = sorted(m for m in modalities if m in _ATAC_MODALITIES)
    paired_h5mu = _paired_h5mu_records(exp_context)

    if paired_h5mu:
        paths = [str(r.get("path")) for r in paired_h5mu if r.get("path")]
        return {
            "status": "inferred",
            "source": "assay_detector",
            "object_type": "paired_mudata",
            "same_cell": True,
            "rna_modality": "scRNA",
            "atac_modality": "scATAC",
            "cell_namespace": "shared_mudata_obs",
            "barcode_namespace": "shared_mudata_obs",
            "evidence_paths": paths,
            "degradation": None,
        }

    if "scRNA" in modalities and "scATAC" in modalities:
        return {
            "status": "inferred",
            "source": "modalities",
            "object_type": "split_single_cell_modalities",
            "same_cell": None,
            "rna_modality": "scRNA",
            "atac_modality": "scATAC",
            "rna_files": list(modalities.get("scRNA") or []),
            "atac_files": list(modalities.get("scATAC") or []),
            "cell_namespace": None,
            "barcode_namespace": None,
            "degradation": "pairing_unconfirmed",
        }

    if rna or atac:
        return {
            "status": "not_applicable",
            "source": "modalities",
            "object_type": "single_assay_or_unpaired",
            "same_cell": False,
            "rna_modalities": rna,
            "atac_modalities": atac,
            "degradation": "missing_rna_or_atac_partner",
        }

    return {
        "status": "not_applicable",
        "source": "modalities",
        "object_type": "no_rna_atac_modalities",
        "same_cell": False,
        "degradation": "no_multiome_candidate",
    }


def validate_multiome_contract(exp_context: dict[str, Any]) -> dict[str, Any]:
    """Validate the RNA+ATAC object contract before integration work.

    Returns a structure suitable for `capability_matrix["contracts"]["multiome"]`.
    Findings intentionally target integration readiness, not the independent RNA
    or ATAC dispatch lanes.
    """
    contract = infer_multiome_contract(exp_context)
    findings: list[dict[str, Any]] = []
    checks = {
        "object_type": contract.get("object_type"),
        "same_cell": contract.get("same_cell"),
        "modalities": _modality_presence(exp_context),
        "cell_namespace": _cell_namespace_check(contract),
        "sample_alignment": _sample_alignment_check(contract),
        "genome_feature_space": _genome_feature_space_check(exp_context, contract),
    }

    if contract.get("status") == "not_applicable":
        return {
            "status": "not_applicable",
            "contract": contract,
            "checks": checks,
            "findings": [],
        }

    if checks["cell_namespace"]["status"] == "blocked":
        findings.append(_finding(
            "blocking",
            "multiome_contract_barcode_namespace_missing",
            (
                "The multiome contract claims same-cell RNA+ATAC, but no shared "
                "cell-barcode namespace or barcode map is declared."
            ),
            (
                "Provide a paired `.h5mu`, set `cell_namespace`/"
                "`barcode_namespace`, or provide `barcode_map`/`cell_barcode_map` "
                "before enabling same-cell integration."
            ),
        ))

    if checks["sample_alignment"]["status"] == "blocked":
        findings.append(_finding(
            "blocking",
            "multiome_contract_sample_mismatch",
            (
                "RNA and ATAC sample identifiers in the multiome contract do "
                "not match."
            ),
            (
                "Correct the multiome contract so RNA and ATAC sample/donor/"
                "condition identifiers refer to the same biological units, or "
                "degrade the run to unpaired cross-modal analysis."
            ),
        ))
    elif checks["sample_alignment"]["status"] == "warning":
        findings.append(_finding(
            "warning",
            "multiome_contract_sample_alignment_unconfirmed",
            (
                "RNA and ATAC sample alignment is not fully declared for this "
                "cross-modal run."
            ),
            (
                "Provide sample-level RNA/ATAC pairing metadata before running "
                "integration that assumes matched samples or donors."
            ),
        ))

    if checks["cell_namespace"]["status"] == "warning":
        findings.append(_finding(
            "warning",
            "multiome_contract_pairing_unconfirmed",
            (
                "scRNA and scATAC inputs are present, but ARIA cannot prove they "
                "come from the same cells."
            ),
            (
                "Use a paired `.h5mu` or provide an explicit multiome contract "
                "with shared barcodes / a barcode map. Until then, WNN must stay "
                "disabled and ARIA should treat the modalities as parallel lanes."
            ),
        ))

    if checks["genome_feature_space"]["status"] == "warning":
        findings.append(_finding(
            "warning",
            "multiome_contract_genome_feature_space_unconfirmed",
            (
                "Genome or feature-space compatibility is not fully declared for "
                "the RNA+ATAC object."
            ),
            (
                "Confirm the genome/assembly and ATAC feature space before "
                "running peak-to-gene or same-cell integration."
            ),
        ))

    status = (
        "red" if any(f["severity"] == "blocking" for f in findings)
        else "yellow" if any(f["severity"] == "warning" for f in findings)
        else "green"
    )
    return {
        "status": status,
        "contract": {**contract, "validation_status": status},
        "checks": checks,
        "findings": findings,
    }


def _paired_h5mu_records(exp_context: dict[str, Any]) -> list[dict[str, Any]]:
    records = []
    for rec in exp_context.get("assay_detections") or []:
        evidence = rec.get("evidence") or {}
        fmt = str(evidence.get("format") or "").lower()
        if fmt in _H5MU_FORMATS and evidence.get("paired_rna"):
            records.append(rec)
    if records:
        return records

    # Filename fallback only marks a candidate. It does not inspect the file; a
    # later script still reads the `.h5mu` and can return structured blockers.
    for path in (exp_context.get("modalities") or {}).get("scATAC", []) or []:
        if str(path).lower().endswith(".h5mu"):
            records.append({
                "path": str(path),
                "evidence": {"format": "h5mu", "paired_rna": True},
                "confidence": "medium",
            })
    return records


def _modality_presence(exp_context: dict[str, Any]) -> dict[str, Any]:
    modalities = exp_context.get("modalities") or {}
    rna = sorted(m for m in modalities if m in _RNA_MODALITIES)
    atac = sorted(m for m in modalities if m in _ATAC_MODALITIES)
    return {
        "rna_modalities": rna,
        "atac_modalities": atac,
        "has_single_cell_rna": "scRNA" in modalities,
        "has_single_cell_atac": "scATAC" in modalities,
        "has_rna_and_atac": bool(rna and atac),
    }


def _cell_namespace_check(contract: dict[str, Any]) -> dict[str, Any]:
    same_cell = contract.get("same_cell")
    namespace = (
        contract.get("cell_namespace")
        or contract.get("barcode_namespace")
        or contract.get("shared_barcode_namespace")
    )
    barcode_map = contract.get("barcode_map") or contract.get("cell_barcode_map")
    paired_obs = bool(contract.get("paired_obs_names"))

    if same_cell is True:
        ok = bool(namespace or barcode_map or paired_obs)
        return {
            "status": "pass" if ok else "blocked",
            "required_for_same_cell": True,
            "cell_namespace": namespace,
            "has_barcode_map": bool(barcode_map),
            "paired_obs_names": paired_obs,
        }

    if same_cell is None:
        return {
            "status": "warning",
            "required_for_same_cell": True,
            "reason": "same_cell_pairing_unconfirmed",
        }

    return {
        "status": "pass",
        "required_for_same_cell": False,
        "reason": "not_same_cell_integration",
    }


def _sample_alignment_check(contract: dict[str, Any]) -> dict[str, Any]:
    rna_ids = _ids(contract, "rna_sample_ids")
    atac_ids = _ids(contract, "atac_sample_ids")
    shared_ids = _ids(contract, "sample_ids")
    same_cell = contract.get("same_cell")

    if shared_ids and not rna_ids and not atac_ids:
        return {
            "status": "pass",
            "sample_ids": shared_ids,
            "shared_sample_ids": True,
        }
    if rna_ids and atac_ids:
        match = set(rna_ids) == set(atac_ids)
        return {
            "status": "pass" if match else ("blocked" if same_cell else "warning"),
            "rna_sample_ids": rna_ids,
            "atac_sample_ids": atac_ids,
            "missing_in_atac": sorted(set(rna_ids) - set(atac_ids)),
            "missing_in_rna": sorted(set(atac_ids) - set(rna_ids)),
        }
    if same_cell is True:
        return {
            "status": "pass",
            "reason": "same_cell_barcode_namespace_supplies_alignment",
        }
    if same_cell is None:
        return {
            "status": "warning",
            "reason": "sample_alignment_unconfirmed",
        }
    return {"status": "pass", "reason": "not_required"}


def _genome_feature_space_check(
    exp_context: dict[str, Any],
    contract: dict[str, Any],
) -> dict[str, Any]:
    genome = (
        exp_context.get("genome")
        or exp_context.get("assembly")
        or contract.get("genome")
        or contract.get("assembly")
    )
    feature_space = (
        contract.get("atac_feature_space")
        or contract.get("feature_space")
        or (
            "peak_matrix"
            if contract.get("object_type") == "paired_mudata"
            else None
        )
    )
    if contract.get("same_cell") is False:
        return {
            "status": "pass",
            "genome": genome,
            "feature_space": feature_space,
            "required": False,
        }
    return {
        "status": "pass" if genome and feature_space else "warning",
        "genome": str(genome) if genome else None,
        "feature_space": feature_space,
        "required": True,
    }


def _ids(contract: dict[str, Any], key: str) -> list[str]:
    value = contract.get(key)
    if isinstance(value, (list, tuple, set)):
        return sorted({str(v) for v in value if str(v).strip()})
    if isinstance(value, str) and value.strip():
        return [value.strip()]
    return []


def _finding(
    severity: str,
    check: str,
    message: str,
    recommendation: str,
) -> dict[str, Any]:
    return {
        "severity": severity,
        "check": check,
        "message": message,
        "recommendation": recommendation,
        "modality": "multiome",
    }


def multiome_evidence_paths(contract: dict[str, Any]) -> list[str]:
    """Return normalized local evidence paths for tests/reporting helpers."""
    return [
        str(Path(p))
        for p in contract.get("evidence_paths") or []
        if str(p).strip()
    ]
