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
    inspected_h5mu_paths: set[str] = set()
    for rec in exp_context.get("assay_detections") or []:
        evidence = rec.get("evidence") or {}
        fmt = str(evidence.get("format") or "").lower()
        if fmt in _H5MU_FORMATS:
            if rec.get("path"):
                inspected_h5mu_paths.add(str(rec.get("path")))
            if evidence.get("paired_rna"):
                records.append(rec)
    if records:
        return records

    # Filename fallback only marks a candidate when AssayDetector has not
    # already inspected the `.h5mu` path and determined `paired_rna` is false.
    seen_paths: set[str] = set()
    mods = exp_context.get("modalities") or {}
    for mod_key in ("scATAC", "scRNA"):
        for path in mods.get(mod_key, []) or []:
            p_str = str(path)
            if (
                p_str.lower().endswith(".h5mu")
                and p_str not in seen_paths
                and p_str not in inspected_h5mu_paths
            ):
                seen_paths.add(p_str)
                records.append({
                    "path": p_str,
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


_SAME_CELL_BLOCKING_CHECKS = frozenset({
    "multiome_contract_barcode_namespace_missing",
    "multiome_contract_sample_mismatch",
    "multiome_contract_pairing_unconfirmed",
    "multiome_contract_sample_alignment_unconfirmed",
})


def enforce_same_cell_integration_gate(
    exp_context: dict[str, Any],
    *,
    require_genome: bool = False,
) -> dict[str, Any]:
    """Hard pre-flight gate for same-cell vertical multiome integration (V48 I1).

    `validate_multiome_contract` intentionally emits `yellow` warnings when split
    scRNA + scATAC inputs lack confirmed cell pairing so independent single-assay
    lanes remain dispatchable. Same-cell integration (paired `.h5mu` joint QC,
    WNN, same-cell MOFA+, cell-level peak-to-gene) is stricter: it MUST fail
    closed unless `same_cell is True`, cell-barcode identity is explicit, and no
    blocking or unconfirmed-pairing issues exist.
    """
    validation = validate_multiome_contract(exp_context)
    status = str(validation.get("status") or "not_applicable")
    contract = dict(validation.get("contract") or {})
    checks = dict(validation.get("checks") or {})
    raw_findings = list(validation.get("findings") or [])

    blockers: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []

    if status == "not_applicable" or contract.get("same_cell") is not True:
        # Promote any unconfirmed-pairing findings to explicit gate blockers, and
        # always emit a primary same-cell requirement blocker when absent.
        for finding in raw_findings:
            check = str(finding.get("check") or "")
            if (
                finding.get("severity") == "blocking"
                or check in _SAME_CELL_BLOCKING_CHECKS
            ):
                blockers.append({**finding, "severity": "blocking"})
            else:
                warnings.append(dict(finding))
        if not blockers:
            blockers.append(_finding(
                "blocking",
                "multiome_same_cell_required",
                (
                    "Same-cell multi-omics integration requires confirmed "
                    "same-cell RNA+ATAC pairing (`same_cell=True`), which is "
                    f"not satisfied (object_type={contract.get('object_type')!r}, "
                    f"status={status!r})."
                ),
                (
                    "Provide a paired `.h5mu` containing RNA and ATAC modalities "
                    "over a shared cell-barcode namespace, or supply an explicit "
                    "`multiome_contract` with `same_cell=True` and verified "
                    "shared barcodes."
                ),
            ))
    else:
        for finding in raw_findings:
            check = str(finding.get("check") or "")
            if (
                finding.get("severity") == "blocking"
                or check in _SAME_CELL_BLOCKING_CHECKS
                or (
                    require_genome
                    and check == "multiome_contract_genome_feature_space_unconfirmed"
                )
            ):
                blockers.append({**finding, "severity": "blocking"})
            else:
                warnings.append(dict(finding))

        evidence_paths = multiome_evidence_paths(contract)
        has_explicit_map = bool(
            contract.get("barcode_map")
            or contract.get("cell_barcode_map")
            or contract.get("paired_obs_names")
        )
        if (
            contract.get("object_type") == "paired_mudata"
            and not evidence_paths
            and not has_explicit_map
        ):
            blockers.append(_finding(
                "blocking",
                "multiome_paired_path_missing",
                (
                    "The multiome contract declares `paired_mudata`, but no "
                    "backing `.h5mu` file path is recorded."
                ),
                (
                    "Include at least one readable `.h5mu` path in "
                    "`evidence_paths` or `modalities` before running paired "
                    "integration."
                ),
            ))
        elif (
            contract.get("object_type") == "paired_mudata"
            and len(evidence_paths) > 1
            and not has_explicit_map
        ):
            blockers.append(_finding(
                "blocking",
                "multiome_multiple_mudata_unsupported",
                (
                    f"Multiple paired `.h5mu` containers ({len(evidence_paths)}) "
                    "were detected without an explicit merged cell-barcode "
                    "namespace."
                ),
                (
                    "Merge shared-donor `.h5mu` containers with collision-safe "
                    "barcodes first, or select a single paired `.h5mu` container "
                    "before running same-cell integration."
                ),
            ))

    allowed = len(blockers) == 0
    return {
        "allowed": allowed,
        "status": "allowed" if allowed else "blocked",
        "reason": None if allowed else "multiome_same_cell_gate_blocked",
        "validation_status": status,
        "contract": contract,
        "checks": checks,
        "blockers": blockers,
        "warnings": warnings,
        "evidence_paths": multiome_evidence_paths(contract),
    }
