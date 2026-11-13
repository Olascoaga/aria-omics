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
# ADR-079: verified paired same-cell containers (MuData, Cell Ranger ARC
# filtered `.h5`, and Cell Ranger ARC MEX directories).
_H5MU_FORMATS = {"h5mu", "10x_multiome_h5", "10x_multiome_mtx"}

# Object types a multiome contract may declare. `paired_mudata` with
# `cohort_merge_required=True` is the ADR-079 multi-donor paired cohort form
# (kept as `paired_mudata` for caller compatibility); `paired_cohort_mudata` is
# accepted as an explicit synonym.
_ALLOWED_OBJECT_TYPES = frozenset({
    "paired_mudata",
    "paired_cohort_mudata",
    "raw_multiome_fastq",
    "explicit_same_cell",
    "explicit_pairing",
    "split_single_cell_modalities",
    "single_assay_or_unpaired",
    "no_rna_atac_modalities",
})

# Cell-barcode namespaces that establish same-cell identity across RNA and ATAC.
_ALLOWED_CELL_NAMESPACES = frozenset({
    "shared_mudata_obs",
    "shared_barcodes",
    "shared_10x_barcodes",
    "10x_multiome_barcodes",
    "10x_arc_whitelist_pair",
    "barcode_map_tsv",
})

# Platforms whose chemistry reads RNA and ATAC from the same nucleus, so a raw
# FASTQ manifest is same-cell by construction (barcode translation still has to
# be declared or derivable from the platform whitelist pair).
_PAIRED_MULTIOME_PLATFORMS = frozenset({
    "10x_multiome_arc_v1",
    "10x_multiome",
    "share_seq",
    "paired_tag",
    "issaac_seq",
    "dogma_seq",
})
_10X_ARC_PLATFORMS = frozenset({"10x_multiome_arc_v1", "10x_multiome"})

# Count-container suffixes that represent an independent RNA or ATAC matrix. A
# paired container mixed with one of these (without an explicit contract) is
# ambiguous; ATAC sidecars (fragments, peaks, indexes, metrics) are not.
_COUNT_CONTAINER_SUFFIXES = (
    ".h5ad", ".h5mu", ".h5", ".loom", ".rds", ".mtx", ".mtx.gz", ".zarr",
)


def infer_multiome_contract(exp_context: dict[str, Any]) -> dict[str, Any]:
    """Infer a conservative multiome contract from audited modalities.

    The result is an initial contract, not proof of integration validity. A
    paired `.h5mu` / Cell Ranger ARC container with RNA+ATAC evidence is treated
    as same-cell; several verified paired containers form a multi-donor paired
    cohort (ADR-079, `cohort_merge_required=True`); a complete raw Multiome FASTQ
    manifest is same-cell by chemistry. Split scRNA and scATAC inputs are
    deliberately marked as unconfirmed unless the caller supplies an explicit
    contract.
    """
    explicit = exp_context.get("multiome_contract")
    if isinstance(explicit, dict) and explicit:
        if str(explicit.get("object_type") or "") == "raw_multiome_fastq":
            # Explicit raw-FASTQ contracts are re-validated against their own
            # FASTQ lists; declaring `same_cell=True` is not sufficient.
            return _merge_explicit_raw_contract(explicit)
        return dict(explicit)

    manifest = exp_context.get("multiome_fastq_manifest")
    if isinstance(manifest, dict) and manifest:
        return _raw_fastq_contract(manifest)

    modalities = exp_context.get("modalities") or {}
    rna = sorted(m for m in modalities if m in _RNA_MODALITIES)
    atac = sorted(m for m in modalities if m in _ATAC_MODALITIES)
    paired_h5mu = _paired_h5mu_records(exp_context)

    if paired_h5mu:
        paired_h5mu = _drop_arc_raw_siblings(paired_h5mu)
        paths = [str(r.get("path")) for r in paired_h5mu if r.get("path")]
        unpaired = _unpaired_container_files(exp_context, paths)
        contract: dict[str, Any] = {
            "status": "inferred",
            "source": "assay_detector",
            "object_type": "paired_mudata",
            "same_cell": True,
            "rna_modality": "scRNA",
            "atac_modality": "scATAC",
            "cell_namespace": "shared_mudata_obs",
            "barcode_namespace": "shared_mudata_obs",
            "evidence_paths": paths,
            "container_formats": sorted({
                str((r.get("evidence") or {}).get("format") or "h5mu").lower()
                for r in paired_h5mu
            }),
            "n_containers": len(paths),
            "cohort_merge_required": len(paths) >= 2,
            "degradation": None,
        }
        if unpaired:
            # A paired container mixed with independent RNA/ATAC count files is
            # ambiguous without an explicit contract: fail closed at the
            # same-cell gate (single-assay lanes are unaffected).
            contract["same_cell"] = False
            contract["cohort_merge_required"] = False
            contract["unpaired_files"] = unpaired
            contract["degradation"] = "mixed_paired_and_unpaired_containers"
        return contract

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
        "object_type_recognized": contract.get("object_type") in _ALLOWED_OBJECT_TYPES,
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

    if contract.get("degradation") == "mixed_paired_and_unpaired_containers":
        findings.append(_finding(
            "blocking",
            "multiome_mixed_paired_and_unpaired_containers",
            (
                "Verified paired RNA+ATAC container(s) were detected together "
                "with separate unpaired scRNA/scATAC count files "
                f"({', '.join(contract.get('unpaired_files') or [])}); same-cell "
                "identity of the combined input cannot be established."
            ),
            (
                "Run same-cell integration on the paired container(s) only, or "
                "supply an explicit `multiome_contract` declaring how the "
                "separate files pair with the paired cells. The separate files "
                "remain available to their single-assay lanes."
            ),
        ))

    raw_blocker = contract.get("raw_fastq_blocker")
    if raw_blocker == "multiome_raw_fastq_incomplete":
        findings.append(_finding(
            "blocking",
            "multiome_raw_fastq_incomplete",
            (
                "The raw Multiome FASTQ manifest does not provide both RNA and "
                "ATAC FASTQs for every sample (incomplete: "
                f"{', '.join(contract.get('incomplete_samples') or [])})."
            ),
            (
                "List `rna_fastqs` and `atac_fastqs` for every sample in "
                "`multiome_fastq_manifest` before same-cell assembly."
            ),
        ))
    elif raw_blocker == "multiome_not_same_cell":
        findings.append(_finding(
            "blocking",
            "multiome_not_same_cell",
            (
                "The raw Multiome FASTQ manifest does not establish same-cell "
                "RNA+ATAC chemistry (same_cell=False, or no known paired "
                f"platform / barcode translation; platform={contract.get('platform')!r})."
            ),
            (
                "Declare a paired-nucleus platform (e.g. `10x_multiome_arc_v1`), "
                "or provide `barcode_map_path` or the "
                "`gex_whitelist_path`+`atac_whitelist_path` pair. Unpaired "
                "RNA/ATAC FASTQs remain separate single-assay lanes."
            ),
        ))

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


def _as_path_list(value: Any) -> list[str]:
    if isinstance(value, (list, tuple)):
        out: list[str] = []
        for item in value:
            if isinstance(item, dict):
                out.extend(
                    str(v) for k, v in item.items()
                    if k in {"r1", "r2", "r3", "i1", "i2", "path", "fastq"} and v
                )
            elif item is not None and str(item).strip():
                out.append(str(item))
        return out
    if isinstance(value, dict):
        return [str(v) for v in value.values() if v is not None and str(v).strip()]
    if isinstance(value, str) and value.strip():
        return [value.strip()]
    return []


def _merge_explicit_raw_contract(explicit: dict[str, Any]) -> dict[str, Any]:
    """Re-derive same-cell status for an explicit `raw_multiome_fastq` contract."""
    derived = _raw_fastq_contract(explicit)
    merged = {**explicit}
    for key in (
        "same_cell", "cell_namespace", "barcode_namespace", "sample_ids",
        "n_samples", "cohort_merge_required", "evidence_paths",
        "incomplete_samples", "raw_fastq_blocker", "platform",
    ):
        if key in derived:
            merged[key] = derived[key]
    merged.setdefault("atac_feature_space", "fragments")
    merged.setdefault("rna_modality", "scRNA")
    merged.setdefault("atac_modality", "scATAC")
    if derived.get("raw_fastq_blocker"):
        merged["degradation"] = derived["raw_fastq_blocker"]
    return merged


def _raw_fastq_contract(manifest: dict[str, Any]) -> dict[str, Any]:
    """Infer a raw Multiome FASTQ contract (ADR-079) from a FASTQ manifest.

    Same-cell identity is only asserted when (a) every sample carries BOTH RNA
    and ATAC FASTQs and (b) the platform is a known paired-nucleus chemistry or
    an explicit barcode translation (map or ARC whitelist pair) is supplied.
    Otherwise the contract is emitted with `same_cell=False` and a blocker code
    so the same-cell gate fails closed.
    """
    platform = str(manifest.get("platform") or "").strip().lower() or None
    samples_raw = manifest.get("samples")
    samples: list[dict[str, Any]] = []
    if isinstance(samples_raw, list) and samples_raw:
        for idx, entry in enumerate(samples_raw):
            if isinstance(entry, dict):
                samples.append({
                    "sample_id": str(entry.get("sample_id") or f"sample_{idx}"),
                    "rna_fastqs": _as_path_list(entry.get("rna_fastqs")),
                    "atac_fastqs": _as_path_list(entry.get("atac_fastqs")),
                })
    else:
        samples.append({
            "sample_id": str(manifest.get("sample_id") or "sample_0"),
            "rna_fastqs": _as_path_list(manifest.get("rna_fastqs")),
            "atac_fastqs": _as_path_list(manifest.get("atac_fastqs")),
        })

    barcode_map = manifest.get("barcode_map_path") or manifest.get("barcode_map")
    gex_wl = manifest.get("gex_whitelist_path")
    atac_wl = manifest.get("atac_whitelist_path")
    declared_ns = manifest.get("cell_namespace") or manifest.get("barcode_namespace")
    if barcode_map:
        namespace = "barcode_map_tsv"
    elif (gex_wl and atac_wl) or platform in _10X_ARC_PLATFORMS:
        namespace = "10x_arc_whitelist_pair"
    elif platform in _PAIRED_MULTIOME_PLATFORMS:
        namespace = "shared_barcodes"
    elif declared_ns in _ALLOWED_CELL_NAMESPACES:
        # Explicit contract with a recognized shared-barcode namespace.
        namespace = str(declared_ns)
    else:
        namespace = None

    explicit_same_cell = manifest.get("same_cell")
    if explicit_same_cell is False:
        same_cell: bool | None = False
    elif explicit_same_cell is True:
        same_cell = True if namespace else None
    else:
        same_cell = True if (platform in _PAIRED_MULTIOME_PLATFORMS or barcode_map
                             or (gex_wl and atac_wl)) else None

    incomplete = [
        s["sample_id"] for s in samples
        if not s["rna_fastqs"] or not s["atac_fastqs"]
    ] if samples else ["sample_0"]

    blocker: str | None = None
    if incomplete:
        blocker = "multiome_raw_fastq_incomplete"
        same_cell = False
    elif same_cell is not True:
        blocker = "multiome_not_same_cell"
        same_cell = False

    contract: dict[str, Any] = {
        "status": "inferred",
        "source": "multiome_fastq_manifest",
        "object_type": "raw_multiome_fastq",
        "same_cell": same_cell,
        "rna_modality": "scRNA",
        "atac_modality": "scATAC",
        "platform": platform,
        "cell_namespace": namespace if same_cell else None,
        "barcode_namespace": namespace if same_cell else None,
        "atac_feature_space": "fragments",
        "sample_ids": [s["sample_id"] for s in samples],
        "n_samples": len(samples),
        "cohort_merge_required": len(samples) >= 2,
        "evidence_paths": [
            p for s in samples for p in (*s["rna_fastqs"], *s["atac_fastqs"])
        ],
        "degradation": blocker,
    }
    if barcode_map:
        contract["barcode_map"] = str(barcode_map)
    if gex_wl:
        contract["gex_whitelist_path"] = str(gex_wl)
    if atac_wl:
        contract["atac_whitelist_path"] = str(atac_wl)
    if incomplete:
        contract["incomplete_samples"] = incomplete
    if blocker:
        contract["raw_fastq_blocker"] = blocker
    return contract


def _is_arc_matrix(path: str, kind: str) -> bool:
    return Path(path.rstrip("/")).name.lower().startswith(f"{kind}_feature_bc_matrix")


def _arc_parent(path: str) -> str:
    return str(Path(path.rstrip("/")).parent)


def _drop_arc_raw_siblings(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Cell Ranger ARC `outs/` ships `raw_` and `filtered_` matrices of ONE run.

    The unfiltered `raw_feature_bc_matrix` next to a `filtered_feature_bc_matrix`
    in the same directory is the same library (not a second donor), so it must
    not turn a single run into a multi-container cohort.
    """
    filtered_dirs = {
        _arc_parent(str(r.get("path")))
        for r in records
        if r.get("path") and _is_arc_matrix(str(r.get("path")), "filtered")
    }
    return [
        r for r in records
        if not (
            r.get("path")
            and _is_arc_matrix(str(r.get("path")), "raw")
            and _arc_parent(str(r.get("path"))) in filtered_dirs
        )
    ]


def _unpaired_container_files(
    exp_context: dict[str, Any],
    paired_paths: list[str],
) -> list[str]:
    """Independent RNA/ATAC count containers present next to paired containers.

    Files inside a paired container directory (Cell Ranger ARC MEX) and ATAC
    sidecars (fragments, peaks, indexes, metrics) are not counted.
    """
    paired = {str(p) for p in paired_paths}
    paired_dirs = tuple(p.rstrip("/") + "/" for p in paired)
    out: list[str] = []
    mods = exp_context.get("modalities") or {}
    for mod_key in ("scRNA", "scATAC"):
        for path in mods.get(mod_key, []) or []:
            p_str = str(path)
            low = p_str.lower()
            if p_str in paired or p_str.startswith(paired_dirs):
                continue
            if _is_arc_matrix(p_str, "raw") and any(
                _arc_parent(p_str) == _arc_parent(pp) for pp in paired
            ):
                # Unfiltered sibling of a Cell Ranger ARC run (same library).
                continue
            if any(low.endswith(sfx) for sfx in _COUNT_CONTAINER_SUFFIXES):
                if p_str not in out:
                    out.append(p_str)
    return sorted(out)


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
            # Informational: whether the namespace is one ARIA knows how to
            # construct (ADR-079); custom namespaces are not blocked.
            "namespace_recognized": (
                namespace in _ALLOWED_CELL_NAMESPACES if namespace else None
            ),
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
            if contract.get("object_type") in _PAIRED_CONTAINER_OBJECT_TYPES
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


_PAIRED_CONTAINER_OBJECT_TYPES = frozenset({"paired_mudata", "paired_cohort_mudata"})

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
        paired_object = (
            contract.get("object_type") in _PAIRED_CONTAINER_OBJECT_TYPES
        )
        verified_cohort = (
            contract.get("cohort_merge_required") is True
            and not contract.get("unpaired_files")
        )
        if (
            paired_object
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
            paired_object
            and len(evidence_paths) > 1
            and not has_explicit_map
            and not verified_cohort
        ):
            # ADR-079: several verified paired containers form a multi-donor
            # cohort only when the contract declares `cohort_merge_required`
            # (always set by inference when every container is verified
            # paired). A hand-written contract listing several containers
            # without that declaration stays fail-closed.
            blockers.append(_finding(
                "blocking",
                "multiome_multiple_mudata_unsupported",
                (
                    f"Multiple paired containers ({len(evidence_paths)}) were "
                    "declared without `cohort_merge_required=True` or an "
                    "explicit merged cell-barcode namespace."
                ),
                (
                    "Declare the containers as a paired cohort "
                    "(`cohort_merge_required=True`; merged onto consensus peaks "
                    "with sample-prefixed barcodes), or select a single paired "
                    "container before running same-cell integration."
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
