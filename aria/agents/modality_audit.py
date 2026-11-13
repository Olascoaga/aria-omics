"""Per-modality readiness cards and capability matrix.

The checks here are deliberately pre-dispatch and data-light. They decide
whether a modality can be dispatched automatically, requires explicit user
acknowledgement, or must be blocked before expensive analysis starts.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from aria.utils.assay_contracts import validate_assay_contract
from aria.utils.design_power import assess_design_power
from aria.utils.multiome_contracts import (
    enforce_same_cell_integration_gate,
    validate_multiome_contract,
)
from aria.utils.reference_integrity import assess_reference_integrity


_STATUS_RANK = {"green": 0, "yellow": 1, "red": 2}
_POLICY_FOR_STATUS = {
    "green": "auto",
    "yellow": "requires_ack",
    "red": "blocked",
}

_CHROMATIN_MODALITIES = {
    "scATAC",
    "bulk_ATAC",
    "ChIP",
    "CUT_AND_RUN",
    "CUT_AND_TAG",
}

# Technical (engineering/science) blockers to ATAC production. These were the
# C1-C6 audit items; all closed. BUT closing them is NOT a tier promotion — see
# _ATAC_GOVERNANCE_BLOCKERS + ADR-056.
_ATAC_PRODUCTION_BLOCKERS = {
    "scATAC": (),
    "bulk_ATAC": (),
}

# ADR-056: governance criteria that must ALL be met by an explicit promotion ADR
# before ATAC leaves beta. Technical blockers being clear does not promote the
# tier — production is never inferred from the audit.
_ATAC_GOVERNANCE_BLOCKERS = (
    "end_to_end_chromatin_ci_green",       # real env solve + chr20 smoke in CI
    "reproducible_preprint_freeze",        # artifacts regenerated from clean checkout
    "multi_annotator_claim_review",        # independent expert review of conclusions
    "explicit_production_promotion_adr",   # a deliberate promotion decision
)

# ADR-056 + ADR-079: paired same-cell RNA+ATAC integration is promoted to
# `production` by ADR-079, so no governance blocker remains for the production
# tier. The tuple below is reported ONLY when a caller explicitly configures the
# lane below production (e.g. a `beta` override): those are the criteria a
# promotion ADR has to close, and production is never inferred from the audit.
_INTEGRATION_GOVERNANCE_BLOCKERS: tuple[str, ...] = ()
_INTEGRATION_SUBPRODUCTION_GOVERNANCE_BLOCKERS = (
    "external_multi_donor_multiome_replication",
    "end_to_end_integration_ci_green",
    "reproducible_preprint_freeze",
    "explicit_production_promotion_adr",
)

# Associative-only sub-layers of the production integration lane (ADR-070 /
# ADR-079): reported with zero causal verbs; HypothesisAgent S10 keeps any
# mechanistic extrapolation quarantined behind SPECULATIVE.
_INTEGRATION_CAVEATED_SUBLAYERS = {
    "peak_to_gene": "production_associative",
    "eregulons": "production_associative",
    "mofa_factors": "production_associative",
}

# WNN bandwidth selection and the cis peak-to-gene background both need a
# minimum number of shared cells for stable estimates.
_INTEGRATION_MIN_SHARED_CELLS = 200
_INTEGRATION_MIN_SHARED_BARCODE_FRACTION = 0.50

_DOUBLET_HINTS = ("doublet", "scrublet", "demuxlet")
_MITO_HINTS = ("pct_counts_mt", "pct_mito", "percent_mito", "mito_frac")


def _finding(
    severity: str,
    check: str,
    message: str,
    recommendation: str,
    *,
    modality: str,
) -> dict[str, Any]:
    return {
        "severity": severity,
        "check": check,
        "message": message,
        "recommendation": recommendation,
        "modality": modality,
    }


def _base_card(
    modality: str,
    *,
    agent: str | None = None,
    validation_level: str = "unknown",
    status: str = "green",
    reason: str = "",
) -> dict[str, Any]:
    return {
        "modality": modality,
        "agent": agent,
        "validation_level": validation_level,
        "status": status,
        "dispatch_policy": _POLICY_FOR_STATUS[status],
        "reason": reason,
        "checks": {},
        "findings": [],
    }


def _set_status(card: dict[str, Any], status: str) -> None:
    if _STATUS_RANK[status] > _STATUS_RANK[card["status"]]:
        card["status"] = status
        card["dispatch_policy"] = _POLICY_FOR_STATUS[status]


def _merge_cards(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    card = deepcopy(base)
    _set_status(card, update.get("status", "green"))
    card["checks"].update(update.get("checks", {}))
    card["findings"].extend(update.get("findings", []))
    if update.get("agent"):
        card["agent"] = update["agent"]
    if update.get("validation_level"):
        card["validation_level"] = update["validation_level"]
    if update.get("reason"):
        card["reason"] = update["reason"]
    return card


def _combined_design(exp_context: dict[str, Any]) -> dict[str, Any]:
    inferred = deepcopy(exp_context.get("inferred_design") or {})
    design = deepcopy(exp_context.get("design") or {})
    if inferred and design:
        merged = inferred
        merged.update(design)
        pb = {}
        pb.update(inferred.get("pseudobulk") or {})
        pb.update(design.get("pseudobulk") or {})
        if pb:
            merged["pseudobulk"] = pb
        return merged
    return design or inferred


def _obs_column_names(design: dict[str, Any]) -> list[str]:
    cols = design.get("obs_columns") or []
    if isinstance(cols, dict):
        names = list(cols.keys())
    else:
        names = list(cols or [])
    return [str(c) for c in names]


def _has_column_hint(columns: list[str], hints: tuple[str, ...]) -> bool:
    lowered = [c.lower() for c in columns]
    return any(any(hint in col for hint in hints) for col in lowered)


_SCATAC_CAVEATED_SUBLAYERS = {
    "per_cluster_marker_peaks": "exploratory",
    "gene_scores": "exploratory",
    "motif_activity": "scaffold",
    "label_transfer": "exploratory",
    "differential_tf_footprinting": "production_associative",
}

_BULK_ATAC_CAVEATED_SUBLAYERS = {
    "differential_tf_footprinting": "production_associative",
}

_EPIGENOMIC_CAVEATED_SUBLAYERS = {
    "motif_enrichment": "production_associative",
    "peak_annotation": "production_associative",
    "peak_ora": "production_associative",
}

_HIC_CAVEATED_SUBLAYERS = {
    "compartments": "production_associative",
    "tads": "production_associative",
    "loops": "production_associative",
    "differential_3d": "production_associative",
}


def _atac_production_readiness(
    modality: str,
    *,
    tier_override: str | None = None,
) -> dict[str, Any]:
    """Chromatin readiness, ADR-056 / ADR-076 / ADR-077 / ADR-078 / ADR-081-honest.

    `scATAC` is promoted to `production` by `ADR-076`, `bulk_ATAC` by `ADR-077`,
    `TOBIAS` differential TF footprinting (`production_associative`) by `ADR-078`,
    and `ChIP`, `CUT_AND_RUN`, and `CUT_AND_TAG` by `ADR-081`.
    """
    technical_blockers = list(_ATAC_PRODUCTION_BLOCKERS.get(modality, ()))
    effective_tier = tier_override or (
        "production"
        if modality in {"scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"}
        else "beta"
    )
    if modality == "scATAC" and effective_tier == "production":
        return {
            "tier": "production",
            "requires_ack": False,
            "production_promoted": True,
            "technical_blockers_cleared": not technical_blockers,
            "open_blockers": technical_blockers,
            "governance_blockers": [],
            "promotion_policy_adr": "ADR-056",
            "promotion_adr": "ADR-076",
            "caveated_sub_layers": dict(_SCATAC_CAVEATED_SUBLAYERS),
        }
    if modality == "bulk_ATAC" and effective_tier == "production":
        return {
            "tier": "production",
            "requires_ack": False,
            "production_promoted": True,
            "technical_blockers_cleared": not technical_blockers,
            "open_blockers": technical_blockers,
            "governance_blockers": [],
            "promotion_policy_adr": "ADR-056",
            "promotion_adr": "ADR-077",
            "caveated_sub_layers": dict(_BULK_ATAC_CAVEATED_SUBLAYERS),
        }
    if modality in {"ChIP", "CUT_AND_RUN", "CUT_AND_TAG"} and effective_tier == "production":
        return {
            "tier": "production",
            "requires_ack": False,
            "production_promoted": True,
            "technical_blockers_cleared": True,
            "open_blockers": [],
            "governance_blockers": [],
            "promotion_policy_adr": "ADR-056",
            "promotion_adr": "ADR-081",
            "caveated_sub_layers": dict(_EPIGENOMIC_CAVEATED_SUBLAYERS),
        }
    return {
        "tier": effective_tier if effective_tier in {"beta", "scaffold"} else "beta",
        "requires_ack": effective_tier != "scaffold",
        # Production is NEVER inferred from the audit (ADR-056). Flip only via a
        # deliberate promotion ADR that clears every governance blocker.
        "production_promoted": False,
        "technical_blockers_cleared": not technical_blockers and effective_tier != "scaffold",
        "open_blockers": technical_blockers,
        "governance_blockers": list(_ATAC_GOVERNANCE_BLOCKERS),
        "promotion_policy_adr": "ADR-056",
    }


def _hic_production_readiness(*, tier_override: str | None = None) -> dict[str, Any]:
    """3D genome (`HiC` / Micro-C) readiness, ADR-056 / ADR-082-honest."""
    effective_tier = tier_override or "production"
    if effective_tier == "production":
        return {
            "tier": "production",
            "requires_ack": False,
            "production_promoted": True,
            "technical_blockers_cleared": True,
            "open_blockers": [],
            "governance_blockers": [],
            "promotion_policy_adr": "ADR-056",
            "promotion_adr": "ADR-082",
            "caveated_sub_layers": dict(_HIC_CAVEATED_SUBLAYERS),
        }
    return {
        "tier": effective_tier,
        "requires_ack": effective_tier != "scaffold",
        "production_promoted": False,
        "technical_blockers_cleared": effective_tier != "scaffold",
        "open_blockers": [],
        "governance_blockers": list(_ATAC_GOVERNANCE_BLOCKERS),
        "promotion_policy_adr": "ADR-056",
        "caveated_sub_layers": dict(_HIC_CAVEATED_SUBLAYERS),
    }


def _replicate_counts(groups: dict[str, Any]) -> dict[str, int]:
    counts = {}
    for condition, reps in (groups or {}).items():
        if reps is None:
            counts[str(condition)] = 0
        elif isinstance(reps, (list, tuple, set)):
            counts[str(condition)] = len([r for r in reps if str(r).strip()])
        else:
            counts[str(condition)] = 1
    return counts


def _batch_condition_confounded(
    design: dict[str, Any],
    exp_context: dict[str, Any] | None = None,
) -> bool:
    groups = design.get("analysis_groups") or design.get("groups") or {}
    batch_map = design.get("batch_map") or {}
    if groups and batch_map:
        condition_batches = {}
        for condition, samples in groups.items():
            batches = set()
            for sample in samples or []:
                sample_text = str(sample)
                for stem, batch in batch_map.items():
                    if str(stem) in sample_text:
                        batches.add(str(batch))
                        break
            if batches:
                condition_batches[str(condition)] = batches

        if len(condition_batches) >= 2:
            seen: dict[tuple[str, ...], list[str]] = {}
            for condition, batches in condition_batches.items():
                key = tuple(sorted(batches))
                seen.setdefault(key, []).append(condition)
            if all(len(conditions) == 1 for conditions in seen.values()):
                return True

    ctx = exp_context if isinstance(exp_context, dict) else {}
    rows: list[tuple[str, dict[str, Any]]] = []
    raw_meta = ctx.get("sample_metadata") or design.get("sample_metadata")
    if isinstance(raw_meta, dict):
        for sid, row in raw_meta.items():
            if isinstance(row, dict):
                rows.append((str(sid), row))
    elif isinstance(raw_meta, list):
        for idx, row in enumerate(raw_meta):
            if isinstance(row, dict):
                sid = str(row.get("sample_id") or row.get("sample") or idx)
                rows.append((sid, row))

    manifest_samples = (ctx.get("bulk_fastq_manifest") or {}).get("samples")
    if isinstance(manifest_samples, list):
        for idx, entry in enumerate(manifest_samples):
            if isinstance(entry, dict):
                sid = str(entry.get("sample_id") or idx)
                meta_row = entry.get("metadata")
                rows.append((sid, meta_row if isinstance(meta_row, dict) else entry))

    if not rows:
        return False

    sample_to_cond: dict[str, str] = {}
    if isinstance(groups, dict):
        for condition, samples in groups.items():
            for sample in samples or []:
                sample_to_cond[str(sample)] = str(condition)

    cond_col = str(
        design.get("condition_col")
        or design.get("main_factor")
        or ctx.get("condition_col")
        or "condition"
    )
    batch_keys: list[str] = []
    for candidate in (
        *(design.get("batch_factors") or []),
        design.get("batch_col"),
        *(ctx.get("batch_factors") or []),
        ctx.get("batch_col"),
        "batch",
    ):
        if candidate and str(candidate) not in batch_keys:
            batch_keys.append(str(candidate))

    for batch_key in batch_keys:
        cond_batches_meta: dict[str, set[str]] = {}
        for sid, row in rows:
            cond = sample_to_cond.get(sid) or row.get(cond_col) or row.get("condition")
            batch_val = row.get(batch_key)
            if cond and batch_val is not None and str(batch_val).strip():
                cond_batches_meta.setdefault(str(cond), set()).add(str(batch_val).strip())
        if len(cond_batches_meta) < 2:
            continue
        batch_to_conditions: dict[str, set[str]] = {}
        for cond, batches in cond_batches_meta.items():
            for b in batches:
                batch_to_conditions.setdefault(b, set()).add(cond)
        if batch_to_conditions and all(
            len(conds) == 1 for conds in batch_to_conditions.values()
        ):
            return True

    return False


class ScRNAAuditAgent:
    """Pre-dispatch readiness checks for scRNA workflows."""

    modality = "scRNA"

    def audit(self, exp_context: dict[str, Any]) -> dict[str, Any]:
        card = _base_card(
            self.modality,
            agent="scrna_agent",
            validation_level="production",
        )
        design = _combined_design(exp_context)
        pb = design.get("pseudobulk") or {}
        groups = design.get("analysis_groups") or design.get("groups") or {}

        condition_col = pb.get("condition_col") or design.get("condition_col")
        replicate_col = pb.get("replicate_col") or design.get("replicate_col")
        groupby_col = pb.get("groupby_col") or design.get("groupby_col")
        counts = _replicate_counts(groups)
        min_reps = min(counts.values()) if counts else 0

        card["checks"]["pseudobulk_design"] = {
            "condition_col": condition_col,
            "replicate_col": replicate_col,
            "groupby_col": groupby_col,
            "has_groups": bool(groups),
        }
        card["checks"]["pseudobulk_replicates"] = {
            "replicates_per_condition": counts,
            "min_replicates": min_reps,
            "publication_ready_min": 3,
            "minimum_supported": 2,
        }

        missing = [
            name
            for name, value in (
                ("condition_col", condition_col),
                ("replicate_col", replicate_col),
                ("groupby_col", groupby_col),
            )
            if not value
        ]
        if missing or not groups:
            _set_status(card, "yellow")
            card["findings"].append(_finding(
                "warning",
                "scrna_pseudobulk_design_incomplete",
                (
                    "scRNA pseudobulk readiness is incomplete "
                    f"(missing: {', '.join(missing) or 'groups'})."
                ),
                (
                    "Confirm condition, biological replicate/donor, cell-group "
                    "column, and explicit contrasts before running inferential "
                    "pseudobulk DE."
                ),
                modality=self.modality,
            ))
        elif min_reps < 2:
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                "scrna_pseudobulk_under_replicated",
                (
                    "scRNA pseudobulk has condition levels with fewer than two "
                    f"biological replicates: {counts}."
                ),
                (
                    "Do not dispatch inferential scRNA pseudobulk DE. Add "
                    "replicates, revise the design, or run only non-inferential "
                    "QC/clustering workflows."
                ),
                modality=self.modality,
            ))
        elif min_reps < 3:
            _set_status(card, "yellow")
            card["findings"].append(_finding(
                "warning",
                "scrna_pseudobulk_low_power",
                (
                    "scRNA pseudobulk has exactly two biological replicates in "
                    f"at least one condition: {counts}."
                ),
                (
                    "Treat DE as low-power. Dispatch requires explicit user "
                    "acknowledgement; three or more replicates per condition is "
                    "the publication-ready default."
                ),
                modality=self.modality,
            ))

        if _batch_condition_confounded(design, exp_context):
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                "scrna_batch_condition_confounding",
                "Each condition maps to a distinct batch in the confirmed design.",
                (
                    "Do not dispatch condition DE from this design. Add a balanced "
                    "batch design or change the biological question."
                ),
                modality=self.modality,
            ))

        obs_columns = _obs_column_names(design)
        card["checks"]["qc_columns"] = {
            "doublet_column_detected": _has_column_hint(obs_columns, _DOUBLET_HINTS),
            "mitochondrial_column_detected": _has_column_hint(obs_columns, _MITO_HINTS),
            "ambient_check": "deferred_to_scrna_qc",
        }
        return card


class ScATACAuditAgent:
    """Pre-dispatch readiness checks for production scATAC workflows (ADR-076).

    Mirrors ``ScRNAAuditAgent`` for multi-condition inferential pseudobulk DA
    while allowing non-contrasted chromatin-characterization runs (`no_de_design`
    or single-condition clustering/peak/motif runs) to dispatch autonomously
    (`green` / `auto`). Optional exploratory/beta sub-layers (per-cluster
    Wilcoxon marker peaks, gene activity scores, TOBIAS footprinting) remain
    explicitly caveated inside the report without gating the core spine.
    """

    modality = "scATAC"

    def audit(self, exp_context: dict[str, Any]) -> dict[str, Any]:
        card = _base_card(
            self.modality,
            agent="chromatin_agent",
            validation_level="production",
        )
        design = _combined_design(exp_context)
        pb = design.get("pseudobulk") or {}
        groups = design.get("analysis_groups") or design.get("groups") or {}

        condition_col = (
            pb.get("condition_col")
            or design.get("condition_col")
            or design.get("main_factor")
            or exp_context.get("condition_col")
        )
        replicate_col = (
            pb.get("replicate_col")
            or design.get("replicate_col")
            or exp_context.get("replicate_col")
            or exp_context.get("replicate_column")
        )
        comparisons = (
            exp_context.get("comparisons")
            or design.get("comparisons")
            or design.get("contrasts")
            or design.get("plan_contrasts")
        )

        counts = _replicate_counts(groups)
        if not counts:
            metadata = exp_context.get("sample_metadata") or {}
            if isinstance(metadata, dict) and metadata:
                by_cond: dict[str, set[str]] = {}
                c_key = str(condition_col or "condition")
                r_key = str(replicate_col or "replicate")
                for sample_id, row in metadata.items():
                    if not isinstance(row, dict):
                        continue
                    cond = row.get(c_key) or row.get("condition")
                    rep = (
                        row.get(r_key)
                        or row.get("donor_id")
                        or row.get("replicate")
                        or sample_id
                    )
                    if cond and rep:
                        by_cond.setdefault(str(cond), set()).add(str(rep))
                counts = {k: len(v) for k, v in by_cond.items()}
        if not counts:
            libraries = (
                (exp_context.get("scatac_fastq_manifest") or {}).get("libraries")
                or []
            )
            if isinstance(libraries, list) and libraries:
                by_cond = {}
                c_key = str(condition_col or "condition")
                for row in libraries:
                    if not isinstance(row, dict):
                        continue
                    meta = row.get("metadata") or {}
                    cond = meta.get(c_key) or meta.get("condition")
                    donor = row.get("donor_id") or meta.get("donor_id") or row.get("sample_id")
                    if cond and donor:
                        by_cond.setdefault(str(cond), set()).add(str(donor))
                counts = {k: len(v) for k, v in by_cond.items()}

        min_reps = min(counts.values()) if counts else 0
        has_inferential_contrast = bool(comparisons) or len(counts) >= 2

        card["checks"]["validation"] = {
            "validated": True,
            "stage": "production",
            "live_orchestrator_run": "done",
            "external_concordance": "done",
            "promotion_adr": "ADR-076",
        }
        card["checks"]["production_readiness"] = _atac_production_readiness(
            self.modality, tier_override="production"
        )
        card["checks"]["pseudobulk_design"] = {
            "condition_col": condition_col,
            "replicate_col": replicate_col,
            "has_groups": bool(counts),
            "inferential_contrast_requested": has_inferential_contrast,
        }
        card["checks"]["pseudobulk_replicates"] = {
            "replicates_per_condition": counts,
            "min_replicates": min_reps,
            "publication_ready_min": 3,
            "minimum_supported": 2,
        }

        if has_inferential_contrast:
            missing = [
                name
                for name, value in (
                    ("condition_col", condition_col),
                    ("replicate_col", replicate_col),
                )
                if not value
            ]
            if not counts or (comparisons and missing):
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    "scatac_pseudobulk_design_incomplete",
                    (
                        "scATAC inferential pseudobulk DA readiness is incomplete "
                        f"(missing: {', '.join(missing) or 'biological replicate groups'})."
                    ),
                    (
                        "Confirm condition, biological replicate/donor column, "
                        "and explicit contrasts before running inferential "
                        "pseudobulk DA."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 2:
                _set_status(card, "red")
                card["findings"].append(_finding(
                    "blocking",
                    "scatac_pseudobulk_under_replicated",
                    (
                        "scATAC pseudobulk contrast has condition levels with "
                        f"fewer than two biological replicates: {counts}."
                    ),
                    (
                        "Do not dispatch inferential scATAC pseudobulk DA with "
                        "n<2 donors/replicates per condition. Add biological "
                        "replicates or remove the cross-condition contrast."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 3:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    "scatac_pseudobulk_low_power",
                    (
                        "scATAC pseudobulk has exactly two biological replicates "
                        f"in at least one condition: {counts}."
                    ),
                    (
                        "Treat pseudobulk DA as low-power. Dispatch requires "
                        "explicit user acknowledgement; three or more biological "
                        "replicates per condition is the publication-ready default."
                    ),
                    modality=self.modality,
                ))

        if _batch_condition_confounded(design, exp_context):
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                "scatac_batch_condition_confounding",
                "Each condition maps to a distinct batch in the confirmed design.",
                (
                    "Do not dispatch condition DA from a batch-confounded design. "
                    "Add a balanced batch design or run non-contrasted chromatin "
                    "characterization."
                ),
                modality=self.modality,
            ))

        return card


class BulkATACAuditAgent:
    """Pre-dispatch readiness checks for production bulk ATAC workflows (ADR-077).

    Allows non-contrasted chromatin profiling runs (QC, MACS3 peak calling,
    consensus peaks, motif enrichment) to dispatch autonomously (`green` /
    `auto`), while enforcing biological-replicate power and batch-confounding
    gates on inferential cross-condition differential accessibility (`n >= 3`
    `green` / `auto`, `n == 2` `yellow` / `requires_ack` low-power gate such as
    ENCODE isogenic duplicates, `n < 2` or batch-confounded `red` / `blocked`).
    Optional `TOBIAS` differential footprinting remains explicitly `beta` in
    `caveated_sub_layers` without gating the core bulk ATAC spine.
    """

    modality = "bulk_ATAC"

    def audit(self, exp_context: dict[str, Any]) -> dict[str, Any]:
        card = _base_card(
            self.modality,
            agent="chromatin_agent",
            validation_level="production",
        )
        design = _combined_design(exp_context)
        groups = design.get("analysis_groups") or design.get("groups") or {}

        condition_col = (
            design.get("condition_col")
            or design.get("main_factor")
            or exp_context.get("condition_col")
        )
        replicate_col = (
            design.get("replicate_col")
            or exp_context.get("replicate_col")
            or exp_context.get("replicate_column")
        )
        comparisons = (
            exp_context.get("comparisons")
            or design.get("comparisons")
            or design.get("contrasts")
            or design.get("plan_contrasts")
        )

        counts = _replicate_counts(groups)
        if not counts:
            metadata = exp_context.get("sample_metadata") or {}
            if isinstance(metadata, dict) and metadata:
                by_cond: dict[str, set[str]] = {}
                c_key = str(condition_col or "condition")
                r_key = str(replicate_col or "replicate")
                for sample_id, row in metadata.items():
                    if not isinstance(row, dict):
                        continue
                    cond = row.get(c_key) or row.get("condition")
                    rep = (
                        row.get(r_key)
                        or row.get("replicate")
                        or row.get("donor_id")
                        or sample_id
                    )
                    if cond and rep:
                        by_cond.setdefault(str(cond), set()).add(str(rep))
                counts = {k: len(v) for k, v in by_cond.items()}
        if not counts:
            manifest_samples = (
                (exp_context.get("bulk_fastq_manifest") or {}).get("samples")
                or []
            )
            if isinstance(manifest_samples, list) and manifest_samples:
                by_cond = {}
                c_key = str(condition_col or "condition")
                r_key = str(replicate_col or "replicate")
                for row in manifest_samples:
                    if not isinstance(row, dict):
                        continue
                    meta = row.get("metadata") or row
                    cond = meta.get(c_key) or meta.get("condition")
                    rep = (
                        meta.get(r_key)
                        or meta.get("replicate")
                        or meta.get("donor_id")
                        or row.get("sample_id")
                    )
                    if cond and rep:
                        by_cond.setdefault(str(cond), set()).add(str(rep))
                counts = {k: len(v) for k, v in by_cond.items()}

        effective_condition_col = condition_col or ("condition" if counts else None)
        effective_replicate_col = replicate_col or ("replicate" if counts else None)
        min_reps = min(counts.values()) if counts else 0
        has_inferential_contrast = bool(comparisons) or len(counts) >= 2

        card["checks"]["validation"] = {
            "validated": True,
            "stage": "production",
            "live_orchestrator_run": "done",
            "external_concordance": "done",
            "promotion_adr": "ADR-077",
        }
        card["checks"]["production_readiness"] = _atac_production_readiness(
            self.modality, tier_override="production"
        )
        card["checks"]["bulk_da_design"] = {
            "condition_col": effective_condition_col,
            "replicate_col": effective_replicate_col,
            "has_groups": bool(counts),
            "inferential_contrast_requested": has_inferential_contrast,
        }
        card["checks"]["bulk_da_replicates"] = {
            "replicates_per_condition": counts,
            "min_replicates": min_reps,
            "publication_ready_min": 3,
            "minimum_supported": 2,
        }

        if has_inferential_contrast:
            if not counts:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    "bulk_atac_design_incomplete",
                    (
                        "Bulk ATAC inferential DA readiness is incomplete "
                        "(missing biological replicate groups / sample metadata)."
                    ),
                    (
                        "Confirm condition, biological replicate mapping, and "
                        "explicit contrasts before running inferential bulk ATAC DA."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 2:
                _set_status(card, "red")
                card["findings"].append(_finding(
                    "blocking",
                    "bulk_atac_under_replicated",
                    (
                        "Bulk ATAC contrast has condition levels with fewer than "
                        f"two biological replicates: {counts}."
                    ),
                    (
                        "Do not dispatch inferential bulk ATAC DA with n<2 "
                        "biological replicates per condition. Add biological "
                        "replicates or run non-contrasted peak/QC profiling."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 3:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    "bulk_atac_low_power",
                    (
                        "Bulk ATAC contrast has exactly two biological replicates "
                        f"in at least one condition: {counts}."
                    ),
                    (
                        "Treat bulk ATAC DA as low-power (e.g., isogenic duplicate "
                        "designs). Dispatch requires explicit user acknowledgement; "
                        "three or more biological replicates per condition is the "
                        "publication-ready default."
                    ),
                    modality=self.modality,
                ))

        if _batch_condition_confounded(design, exp_context):
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                "bulk_atac_batch_condition_confounding",
                "Each condition maps to a distinct batch in the confirmed design.",
                (
                    "Do not dispatch condition DA from a batch-confounded design. "
                    "Add a balanced batch design or run non-contrasted chromatin "
                    "characterization."
                ),
                modality=self.modality,
            ))

        return card


class _BaseEpigenomicAuditAgent:
    """Shared pre-dispatch readiness checks for ADR-081 epigenomic workflows
    (`ChIP`, `CUT_AND_RUN`, `CUT_AND_TAG`).

    Excludes Input/IgG control libraries when counting biological IP replicates
    per condition, dispatches non-contrasted profiling and `n >= 3` contrasts
    automatically (`green` / `auto`), escalates `n == 2` contrasts or explicit
    non-UMI duplicate removal on sparse cleavage assays to `yellow` /
    `requires_ack`, and blocks `n < 2` contrasts or batch-confounded designs
    (`red` / `blocked`).
    """

    modality: str = "ChIP"
    check_prefix: str = "chip"

    @staticmethod
    def _is_control_row(row: dict[str, Any], condition_col: str) -> bool:
        if not isinstance(row, dict):
            return False
        if row.get("is_control") is True:
            return True
        role = str(row.get("role") or row.get("assay_role") or "").strip().lower()
        if role in {"control", "input", "igg", "mock"}:
            return True
        cond = str(row.get(condition_col) or row.get("condition") or "").strip().lower()
        return cond in {"input", "igg"}

    def audit(self, exp_context: dict[str, Any]) -> dict[str, Any]:
        card = _base_card(
            self.modality,
            agent="chromatin_agent",
            validation_level="production",
        )
        design = _combined_design(exp_context)
        raw_groups = design.get("analysis_groups") or design.get("groups") or {}
        groups = (
            {
                k: v for k, v in raw_groups.items()
                if str(k).strip().lower() not in {"input", "igg"}
            }
            if isinstance(raw_groups, dict)
            else {}
        )
        condition_col = (
            design.get("condition_col")
            or design.get("main_factor")
            or exp_context.get("condition_col")
        )
        replicate_col = (
            design.get("replicate_col")
            or exp_context.get("replicate_col")
            or exp_context.get("replicate_column")
        )
        comparisons = (
            exp_context.get("comparisons")
            or design.get("comparisons")
            or design.get("contrasts")
            or design.get("plan_contrasts")
        )

        counts = _replicate_counts(groups)
        has_control = bool(
            exp_context.get("control_files") or exp_context.get("control_bams")
        )
        metadata = (
            exp_context.get("sample_metadata")
            or design.get("sample_metadata")
            or design.get("sample_sheet")
        )
        rows_iter: list[tuple[str, dict[str, Any]]] = []
        if isinstance(metadata, dict) and metadata:
            for sample_id, row in metadata.items():
                if isinstance(row, dict):
                    rows_iter.append((str(sample_id), row))
        elif isinstance(metadata, list) and metadata:
            for idx, row in enumerate(metadata):
                if isinstance(row, dict):
                    sid = str(
                        row.get("sample_id") or row.get("sample") or f"s{idx}"
                    )
                    rows_iter.append((sid, row))

        if rows_iter:
            c_key = str(condition_col or "condition")
            r_key = str(replicate_col or "replicate")
            by_cond: dict[str, set[str]] = {}
            for sample_id, row in rows_iter:
                if self._is_control_row(row, c_key):
                    has_control = True
                    continue
                cond = row.get(c_key) or row.get("condition")
                rep = (
                    row.get(r_key)
                    or row.get("replicate")
                    or row.get("donor_id")
                    or sample_id
                )
                if cond and rep:
                    by_cond.setdefault(str(cond), set()).add(str(rep))
            if by_cond:
                counts = {k: len(v) for k, v in by_cond.items()}

        min_reps = min(counts.values()) if counts else 0
        has_inferential_contrast = bool(comparisons) or len(counts) >= 2

        card["checks"]["validation"] = {
            "validated": True,
            "stage": "production",
            "live_orchestrator_run": "done",
            "external_concordance": "done",
            "promotion_adr": "ADR-081",
        }
        card["checks"]["production_readiness"] = _atac_production_readiness(
            self.modality, tier_override="production"
        )
        card["checks"]["epigenomic_controls"] = {
            "has_input_or_igg_control": has_control,
            "spikein_configured": bool(
                exp_context.get("spikein_counts")
                or exp_context.get("spikein_scale_factors")
                or exp_context.get("spikein_contig_prefix")
                or exp_context.get("spikein_bwa_index")
                or exp_context.get("spikein_fasta")
            ),
        }
        card["checks"]["bulk_da_design"] = {
            "condition_col": condition_col or ("condition" if counts else None),
            "replicate_col": replicate_col or ("replicate" if counts else None),
            "has_groups": bool(counts),
            "inferential_contrast_requested": has_inferential_contrast,
        }
        card["checks"]["bulk_da_replicates"] = {
            "replicates_per_condition": counts,
            "min_replicates": min_reps,
            "publication_ready_min": 3,
            "minimum_supported": 2,
        }

        files = list(
            (exp_context.get("modalities") or {}).get(self.modality) or []
        )
        raw_or_bam = any(
            str(f).lower().endswith(
                (".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz")
            )
            for f in files
        )
        if self.modality == "ChIP" and raw_or_bam and not has_control:
            _set_status(card, "yellow")
            card["findings"].append(_finding(
                "warning",
                "chip_missing_input_control",
                (
                    "ChIP-seq alignment/FASTQ inputs do not include a matched "
                    "Input or IgG background control library."
                ),
                (
                    "Provide `control_files` / `control_bams` (or Input/IgG rows "
                    "in sample metadata) to correct sonication and open-chromatin "
                    "bias, or acknowledge control-free ChIP peak calling."
                ),
                modality=self.modality,
            ))

        if (
            self.modality in {"CUT_AND_RUN", "CUT_AND_TAG"}
            and exp_context.get("remove_duplicates") is True
            and not exp_context.get("has_umi")
        ):
            _set_status(card, "yellow")
            card["findings"].append(_finding(
                "warning",
                f"{self.check_prefix}_duplicate_removal_forced",
                (
                    f"{self.modality} has `remove_duplicates=True` without UMI "
                    "barcodes; enzymatic MNase/Tn5 cleavage concentrates true "
                    "biological fragments at identical 5' ends."
                ),
                (
                    "Leave `remove_duplicates=False` (ARIA default for sparse "
                    "CUT&RUN/CUT&Tag) unless high PCR duplication was confirmed."
                ),
                modality=self.modality,
            ))

        if has_inferential_contrast:
            if not counts:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    f"{self.check_prefix}_design_incomplete",
                    (
                        f"{self.modality} inferential differential binding "
                        "readiness is incomplete (missing biological IP "
                        "replicate groups / sample metadata)."
                    ),
                    (
                        "Confirm condition, biological IP replicate mapping, and "
                        "explicit contrasts before running differential binding."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 2:
                _set_status(card, "red")
                card["findings"].append(_finding(
                    "blocking",
                    f"{self.check_prefix}_under_replicated",
                    (
                        f"{self.modality} contrast has condition levels with "
                        f"fewer than two biological IP replicates: {counts}."
                    ),
                    (
                        f"Do not dispatch inferential {self.modality} differential "
                        "binding with n<2 biological IP replicates per condition. "
                        "Add biological replicates or run non-contrasted profiling."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 3:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    f"{self.check_prefix}_low_power",
                    (
                        f"{self.modality} contrast has exactly two biological IP "
                        f"replicates in at least one condition: {counts}."
                    ),
                    (
                        f"Treat {self.modality} differential binding as low-power "
                        "(e.g., ENCODE isogenic duplicates). Dispatch requires "
                        "explicit acknowledgement; n>=3 is preferred."
                    ),
                    modality=self.modality,
                ))

        if _batch_condition_confounded(design, exp_context):
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                f"{self.check_prefix}_batch_condition_confounding",
                "Each condition maps to a distinct batch in the confirmed design.",
                (
                    "Do not dispatch differential binding from a batch-confounded "
                    "design. Add a balanced batch design or run non-contrasted "
                    "profiling."
                ),
                modality=self.modality,
            ))

        return card


class ChIPAuditAgent(_BaseEpigenomicAuditAgent):
    """Pre-dispatch readiness checks for production ChIP-seq workflows (ADR-081)."""

    modality = "ChIP"
    check_prefix = "chip"


class CutAndRunAuditAgent(_BaseEpigenomicAuditAgent):
    """Pre-dispatch readiness checks for production CUT&RUN workflows (ADR-081)."""

    modality = "CUT_AND_RUN"
    check_prefix = "cut_and_run"


class CutAndTagAuditAgent(_BaseEpigenomicAuditAgent):
    """Pre-dispatch readiness checks for production CUT&TAG workflows (ADR-081)."""

    modality = "CUT_AND_TAG"
    check_prefix = "cut_and_tag"


class HiCAuditAgent:
    """Pre-dispatch readiness checks for production 3D genome (`HiC` / Micro-C)
    workflows (ADR-082)."""

    modality = "HiC"

    def audit(
        self,
        exp_context: dict[str, Any],
        *,
        tier_override: str | None = None,
    ) -> dict[str, Any]:
        effective_tier = tier_override or "production"
        if effective_tier != "production":
            status = "red" if effective_tier == "scaffold" else "yellow"
            card = _base_card(
                self.modality,
                agent="genome_arch_agent",
                validation_level=effective_tier,
                status=status,
                reason=(
                    "Hi-C validation is configured below production."
                    if effective_tier == "scaffold"
                    else f"Hi-C {effective_tier} dispatch requires explicit acknowledgement."
                ),
            )
            card["checks"]["validation"] = {
                "validated": effective_tier != "scaffold",
                "stage": effective_tier,
            }
            card["checks"]["production_readiness"] = _hic_production_readiness(
                tier_override=effective_tier
            )
            if effective_tier == "scaffold":
                card["findings"].append(_finding(
                    "blocking",
                    "hic_readiness_scaffold",
                    "HiC is configured at scaffold tier and is not dispatch-ready.",
                    "Use ADR-082 production configuration before dispatch.",
                    modality=self.modality,
                ))
            else:
                card["findings"].append(_finding(
                    "warning",
                    "hic_beta_ack_required",
                    f"HiC is configured at {effective_tier} tier and requires CP3.5 acknowledgement.",
                    "Proceed only after explicit CP3.5 acknowledgement.",
                    modality=self.modality,
                ))
            return card

        card = _base_card(
            self.modality,
            agent="genome_arch_agent",
            validation_level="production",
        )
        design = _combined_design(exp_context)
        groups = design.get("analysis_groups") or design.get("groups") or {}
        condition_col = (
            design.get("condition_col")
            or design.get("main_factor")
            or exp_context.get("condition_col")
        )
        replicate_col = (
            design.get("replicate_col")
            or exp_context.get("replicate_col")
            or exp_context.get("replicate_column")
        )
        comparisons = (
            exp_context.get("comparisons")
            or design.get("comparisons")
            or design.get("contrasts")
            or design.get("plan_contrasts")
        )
        counts = _replicate_counts(groups)
        metadata = (
            exp_context.get("sample_metadata")
            or design.get("sample_metadata")
            or design.get("sample_sheet")
        )
        rows_iter: list[tuple[str, dict[str, Any]]] = []
        if isinstance(metadata, dict) and metadata:
            for sample_id, row in metadata.items():
                if isinstance(row, dict):
                    rows_iter.append((str(sample_id), row))
        elif isinstance(metadata, list) and metadata:
            for idx, row in enumerate(metadata):
                if isinstance(row, dict):
                    sid = str(
                        row.get("sample_id") or row.get("sample") or f"s{idx}"
                    )
                    rows_iter.append((sid, row))
        if rows_iter:
            c_key = str(condition_col or "condition")
            r_key = str(replicate_col or "replicate")
            by_cond: dict[str, set[str]] = {}
            for sample_id, row in rows_iter:
                cond = row.get(c_key) or row.get("condition")
                rep = (
                    row.get(r_key)
                    or row.get("replicate")
                    or row.get("donor_id")
                    or sample_id
                )
                if cond and rep:
                    by_cond.setdefault(str(cond), set()).add(str(rep))
            if by_cond:
                counts = {k: len(v) for k, v in by_cond.items()}

        files = list((exp_context.get("modalities") or {}).get("HiC") or [])
        raw_like = any(
            str(f).lower().endswith(
                (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".sam", ".pairs", ".pairs.gz")
            )
            for f in files
        )
        has_genome = bool(
            exp_context.get("genome")
            or exp_context.get("assembly")
            or exp_context.get("chromsizes")
            or exp_context.get("genome_fasta")
            or design.get("genome")
        )
        min_reps = min(counts.values()) if counts else 0
        has_inferential_contrast = bool(comparisons)

        card["checks"]["validation"] = {
            "validated": True,
            "stage": "production",
            "live_orchestrator_run": "done",
            "external_concordance": "done",
            "promotion_adr": "ADR-082",
        }
        card["checks"]["production_readiness"] = _hic_production_readiness(
            tier_override="production"
        )
        card["checks"]["hic_design"] = {
            "condition_col": condition_col or ("condition" if counts else None),
            "replicate_col": replicate_col or ("replicate" if counts else None),
            "has_groups": bool(counts),
            "inferential_contrast_requested": has_inferential_contrast,
            "phasing_reference_available": bool(
                exp_context.get("genome_fasta") or exp_context.get("gtf_path")
            ),
        }
        card["checks"]["hic_replicates"] = {
            "replicates_per_condition": counts,
            "min_replicates": min_reps,
            "publication_ready_min": 3,
            "minimum_supported": 2,
        }

        if raw_like and not has_genome:
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                "hic_missing_genome_assembly",
                (
                    "Hi-C raw FASTQ/BAM/pairs intake requires a confirmed genome "
                    "assembly, reference FASTA, or chromsizes file."
                ),
                (
                    "Confirm `genome` / `genome_fasta` / `chromsizes` before "
                    "dispatching raw Hi-C preprocessing."
                ),
                modality=self.modality,
            ))

        req_loop_res = exp_context.get("loop_resolution") or exp_context.get("hic_resolution")
        if req_loop_res is not None:
            try:
                if int(req_loop_res) > 25_000 and exp_context.get("require_loops") is True:
                    _set_status(card, "yellow")
                    card["findings"].append(_finding(
                        "warning",
                        "hic_coarse_resolution_for_loops",
                        (
                            f"Requested Hi-C resolution ({int(req_loop_res)} bp) is "
                            "coarser than the <=25 kb resolution needed for focal "
                            "donut loop calling."
                        ),
                        (
                            "Provide a <=10-25 kb cooler/mcool for chromatin loop "
                            "calling, or acknowledge that loops will skip honestly."
                        ),
                        modality=self.modality,
                    ))
            except (TypeError, ValueError):
                pass

        if has_inferential_contrast:
            if not counts:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    "hic_design_incomplete",
                    (
                        "Hi-C differential 3D contrast was requested without "
                        "biological replicate groups / sample metadata."
                    ),
                    (
                        "Provide condition and biological replicate metadata, or "
                        "acknowledge that differential 3D will skip / run "
                        "descriptively."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 2:
                _set_status(card, "red")
                card["findings"].append(_finding(
                    "blocking",
                    "hic_under_replicated",
                    (
                        "Hi-C differential 3D contrast has condition levels with "
                        f"fewer than two biological replicates: {counts}."
                    ),
                    (
                        "Do not dispatch inferential differential 3D compartment/"
                        "insulation testing with n<2 replicates per condition."
                    ),
                    modality=self.modality,
                ))
            elif min_reps < 3:
                _set_status(card, "yellow")
                card["findings"].append(_finding(
                    "warning",
                    "hic_low_power",
                    (
                        "Hi-C differential 3D contrast has exactly two biological "
                        f"replicates in at least one condition: {counts}."
                    ),
                    (
                        "Treat 2v2 Hi-C differential 3D calls as low-power; "
                        "dispatch requires explicit acknowledgement (n>=3 preferred)."
                    ),
                    modality=self.modality,
                ))

        if _batch_condition_confounded(design, exp_context):
            _set_status(card, "red")
            card["findings"].append(_finding(
                "blocking",
                "hic_batch_condition_confounding",
                "Each condition maps to a distinct batch in the confirmed design.",
                (
                    "Do not dispatch cross-condition differential 3D testing from "
                    "a batch-confounded design."
                ),
                modality=self.modality,
            ))

        return card


class ChromatinAuditAgent:
    """Readiness cards for chromatin modalities (`scATAC`, `bulk_ATAC`, `ChIP`,
    `CUT_AND_RUN`, `CUT_AND_TAG`)."""

    _EPIGENOMIC_AUDITORS = {
        "ChIP": ChIPAuditAgent,
        "CUT_AND_RUN": CutAndRunAuditAgent,
        "CUT_AND_TAG": CutAndTagAuditAgent,
    }

    def audit(
        self,
        exp_context: dict[str, Any],
        modality: str,
        *,
        tier_override: str | None = None,
    ) -> dict[str, Any]:
        if modality == "scATAC":
            effective_tier = tier_override or "production"
            if effective_tier == "production":
                return ScATACAuditAgent().audit(exp_context)
            card = _base_card(
                modality,
                agent="chromatin_agent",
                validation_level=effective_tier,
                status="yellow",
                reason=(
                    f"scATAC {effective_tier} dispatch is available only after "
                    "explicit user acknowledgement."
                ),
            )
            card["checks"]["validation"] = {
                "validated": True,
                "stage": effective_tier,
                "live_orchestrator_run": "done",
                "external_concordance": "done",
            }
            card["checks"]["production_readiness"] = _atac_production_readiness(
                modality, tier_override="beta"
            )
            if effective_tier == "beta":
                card["findings"].append(_finding(
                    "warning",
                    "chromatin_readiness_beta_ack_required",
                    (
                        "scATAC is explicitly configured at beta tier and "
                        "requires CP3.5 acknowledgement."
                    ),
                    (
                        "Proceed only after explicit CP3.5 acknowledgement; use "
                        "the pseudobulk DA lane with biological replicates and "
                        "treat per-cluster wilcoxon markers as exploratory."
                    ),
                    modality=modality,
                ))
            return card

        if modality == "bulk_ATAC":
            effective_tier = tier_override or "production"
            if effective_tier == "production":
                return BulkATACAuditAgent().audit(exp_context)
            card = _base_card(
                modality,
                agent="chromatin_agent",
                validation_level=effective_tier,
                status="yellow",
                reason=(
                    "Bulk ATAC beta dispatch supports measured QC, MACS3 peak "
                    "calling, and replicate-gated DESeq2 DA with explicit "
                    "metadata; it requires explicit acknowledgement."
                ),
            )
            card["checks"]["validation"] = {
                "validated": True,
                "stage": "beta_qc_peak_calling_peak_counts_da",
                "differential_accessibility": "beta_replicate_gated",
            }
            card["checks"]["production_readiness"] = _atac_production_readiness(
                modality, tier_override="beta"
            )
            if effective_tier == "beta":
                card["findings"].append(_finding(
                    "warning",
                    "bulk_atac_beta_ack_required",
                    (
                        "Bulk ATAC is open as a beta lane for QC, peak calling, "
                        "peak-count matrices, and replicate-aware DA."
                    ),
                    (
                        "Proceed only after explicit CP3.5 acknowledgement. DA runs "
                        "only with explicit condition, biological replicate, and "
                        "comparison metadata; under-specified designs skip honestly."
                    ),
                    modality=modality,
                ))
            return card

        if modality in self._EPIGENOMIC_AUDITORS:
            effective_tier = tier_override or "production"
            if effective_tier == "production":
                return self._EPIGENOMIC_AUDITORS[modality]().audit(exp_context)
            if effective_tier == "beta":
                card = _base_card(
                    modality,
                    agent="chromatin_agent",
                    validation_level="beta",
                    status="yellow",
                    reason=(
                        f"{modality} beta dispatch requires explicit CP3.5 "
                        "acknowledgement."
                    ),
                )
                card["checks"]["validation"] = {"validated": True, "stage": "beta"}
                card["checks"]["production_readiness"] = _atac_production_readiness(
                    modality, tier_override="beta"
                )
                card["findings"].append(_finding(
                    "warning",
                    f"{modality.lower()}_beta_ack_required",
                    f"{modality} is configured at beta tier and requires CP3.5 acknowledgement.",
                    "Proceed only after explicit CP3.5 acknowledgement.",
                    modality=modality,
                ))
                return card

        card = _base_card(
            modality,
            agent="chromatin_agent",
            validation_level="scaffold",
            status="red",
            reason="Chromatin modality validation is not closed.",
        )
        card["checks"]["validation"] = {"validated": False, "stage": "scaffold"}
        card["findings"].append(_finding(
            "blocking",
            "chromatin_readiness_scaffold",
            f"{modality} is scaffolded and is not dispatch-ready.",
            "Wait for the v4.6 chromatin validation lane before dispatch.",
            modality=modality,
        ))
        return card


def _integration_production_readiness(
    *,
    tier_override: str | None = None,
    promotion_adr: str | None = None,
) -> dict[str, Any]:
    """Integration readiness, ADR-056 / ADR-079-honest.

    Mirrors `_atac_production_readiness`: paired same-cell RNA+ATAC integration
    (`scope="paired_same_cell_rna_atac"`) is explicitly promoted to `production`
    by `ADR-079` (verified paired `.h5mu`, Cell Ranger ARC `.h5` / MEX, raw
    Multiome FASTQs, and multi-donor paired cohorts with batch-aware WNN).
    Peak-to-gene links, candidate eRegulons, and MOFA+ factors stay
    `production_associative` (zero causal claims). Any explicit sub-production
    tier keeps the legacy acknowledgement gate and governance criteria.
    """
    effective_tier = tier_override or "production"
    if effective_tier == "production":
        return {
            "tier": "production",
            "requires_ack": False,
            "production_promoted": True,
            "technical_blockers_cleared": True,
            "open_blockers": [],
            "governance_blockers": list(_INTEGRATION_GOVERNANCE_BLOCKERS),
            "promotion_policy_adr": "ADR-056",
            "promotion_adr": promotion_adr or "ADR-079",
            "caveated_sub_layers": dict(_INTEGRATION_CAVEATED_SUBLAYERS),
        }
    return {
        "tier": effective_tier,
        "requires_ack": effective_tier != "scaffold",
        # Production is NEVER inferred from the audit (ADR-056).
        "production_promoted": False,
        "technical_blockers_cleared": effective_tier != "scaffold",
        "open_blockers": [],
        "governance_blockers": list(_INTEGRATION_SUBPRODUCTION_GOVERNANCE_BLOCKERS),
        "promotion_policy_adr": "ADR-056",
        "promotion_adr": promotion_adr or "ADR-071",
        "caveated_sub_layers": dict(_INTEGRATION_CAVEATED_SUBLAYERS),
    }


def _as_count(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        n = int(value)
    except (TypeError, ValueError):
        return None
    return n if n >= 0 else None


def _as_fraction(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if 0.0 <= f <= 1.0 else None


class IntegrationAuditAgent:
    """Pre-dispatch readiness for paired same-cell RNA+ATAC integration (ADR-079).

    Mirrors `ScATACAuditAgent` / `BulkATACAuditAgent`: at the `production` tier a
    verified paired multiome input dispatches automatically (`green` / `auto`)
    and is downgraded to `yellow` / `requires_ack` only by a measured data or
    design risk (too few shared cells, partial barcode overlap, donor/batch
    confounding in a multi-container cohort). Inputs that fail the same-cell
    gate (unpaired scRNA+scATAC, bulk cross-modal, ATAC-only `.h5mu`, mixed
    paired+unpaired containers) are `red` / `blocked` for the integration lane
    only — the result is kept out of per-assay `dispatch`, so single-assay lanes
    are never blocked here. An explicit sub-production tier (`beta` override)
    preserves the legacy CP3.5 acknowledgement behaviour.
    """

    modality = "multiome"

    def audit(
        self,
        exp_context: dict[str, Any],
        multiome_contract: dict[str, Any],
        integration_validation: dict[str, Any] | None,
    ) -> dict[str, Any]:
        meta = dict(integration_validation or {})
        tier = str(meta.get("level") or "production")
        default_ack = tier in {"alpha", "beta", "experimental"}
        requires_ack_policy = bool(meta.get("requires_ack", default_ack))
        enabled = bool(meta.get("dispatch_enabled", True))
        scope = str(meta.get("scope") or "paired_same_cell_rna_atac")
        readiness = _integration_production_readiness(
            tier_override=tier,
            promotion_adr=meta.get("promotion_adr"),
        )
        production_tier = tier == "production" and not requires_ack_policy

        def _result(status: str, policy: str, **extra: Any) -> dict[str, Any]:
            out: dict[str, Any] = {
                "status": status,
                "tier": tier,
                "dispatch_enabled": enabled,
                "requires_ack": requires_ack_policy,
                "requires_ack_for_dispatch": False,
                "dispatch_policy": policy,
                "scope": scope,
                "production_promoted": bool(
                    production_tier and readiness["production_promoted"]
                ),
                "promotion_policy_adr": readiness["promotion_policy_adr"],
                "promotion_adr": readiness["promotion_adr"],
                "governance_blockers": (
                    [] if production_tier
                    else list(readiness["governance_blockers"])
                ),
                "caveated_sub_layers": dict(readiness["caveated_sub_layers"]),
                "findings": [],
            }
            out.update(extra)
            return out

        presence = (multiome_contract.get("checks") or {}).get("modalities") or {}
        if (
            multiome_contract.get("status") == "not_applicable"
            and not presence.get("has_rna_and_atac")
        ):
            return _result("not_applicable", "not_applicable")
        # RNA + ATAC present without a same-cell contract (e.g. bulk_RNA +
        # bulk_ATAC) is a cross-modal request the gate must explicitly block.

        gate = enforce_same_cell_integration_gate(exp_context)
        if not enabled:
            return _result(
                "red", "blocked",
                dispatch_enabled=False,
                reason="agent_not_validated",
                same_cell_gate_allowed=bool(gate.get("allowed")),
            )

        if not gate.get("allowed"):
            return _result(
                "red", "blocked",
                reason="multiome_same_cell_gate_blocked",
                same_cell_gate_allowed=False,
                blockers=list(gate.get("blockers") or []),
            )

        if tier == "scaffold":
            return _result(
                "red", "blocked",
                reason="integration_scaffold_tier",
                same_cell_gate_allowed=True,
                findings=[_finding(
                    "blocking",
                    "integration_scaffold_not_dispatchable",
                    "IntegrationAgent is configured at scaffold tier.",
                    "Do not dispatch same-cell integration at scaffold tier.",
                    modality=self.modality,
                )],
            )

        if not production_tier:
            # Legacy beta / explicit-ack behaviour (ADR-071), preserved.
            findings: list[dict[str, Any]] = []
            if integration_validation is not None:
                findings.append(_finding(
                    "warning",
                    "integration_beta_ack_required",
                    (
                        f"Cross-modal IntegrationAgent is {tier} "
                        "(`paired_same_cell_rna_atac`; ADR-056 / "
                        f"{readiness['promotion_adr']}) and requires explicit "
                        "CP3.5 acknowledgement before dispatch."
                    ),
                    (
                        f"Proceed only if you accept paired-multiome {tier} "
                        "limitations (associative peak-to-gene / candidate "
                        "eRegulons; WNN joint topology audited against "
                        "per-modality structure, never causal)."
                    ),
                    modality=self.modality,
                ))
            return _result(
                "yellow", "requires_ack",
                requires_ack_for_dispatch=integration_validation is not None,
                reason=str(meta.get("reason") or ""),
                same_cell_gate_allowed=True,
                findings=findings,
            )

        # ── production (ADR-079): green unless a measured risk is present ──
        contract = dict(gate.get("contract") or {})
        design = _combined_design(exp_context)
        checks = self._risk_checks(exp_context, design, contract)
        findings = []
        if checks["shared_cells"]["status"] == "warning":
            n = checks["shared_cells"]["n_shared_cells"]
            findings.append(_finding(
                "warning",
                "integration_low_cell_count",
                (
                    f"Paired multiome input reports {n} shared cells "
                    f"(< {_INTEGRATION_MIN_SHARED_CELLS}); WNN modality-weight "
                    "bandwidths and the cis peak-to-gene background are "
                    "unstable at this size."
                ),
                (
                    "Proceed only after acknowledging the low cell count, or "
                    "add cells/samples. Treat neighbourhood weights and "
                    "peak-to-gene links as exploratory."
                ),
                modality=self.modality,
            ))
        if checks["barcode_overlap"]["status"] == "warning":
            overlap = checks["barcode_overlap"]
            findings.append(_finding(
                "warning",
                "integration_partial_barcode_overlap",
                (
                    "RNA and ATAC barcodes only partially overlap "
                    f"(shared fraction={overlap['shared_barcode_fraction']}, "
                    f"allow_partial_overlap={overlap['allow_partial_overlap']})."
                ),
                (
                    "Integration will run on the shared-barcode intersection "
                    "only. Confirm barcode translation (ARC whitelist pair / "
                    "barcode map) before acknowledging."
                ),
                modality=self.modality,
            ))
        if checks["donor_design"]["status"] == "warning":
            findings.append(_finding(
                "warning",
                "integration_donor_confounded_design",
                (
                    "Multi-sample paired multiome cohort has condition "
                    "completely confounded with batch/donor (each batch holds "
                    "a single condition)."
                ),
                (
                    "Batch-aware WNN cannot separate condition from batch in "
                    "this design; harmonization may remove the condition "
                    "signal. Proceed only after acknowledgement and report "
                    "cross-condition integrated structure as descriptive."
                ),
                modality=self.modality,
            ))

        notes: list[str] = []
        if contract.get("cohort_merge_required") is True:
            notes.append(
                "Multi-sample paired multiome cohort will be merged onto "
                "sweep-line consensus ATAC peaks prior to joint QC and "
                "batch-aware WNN."
            )
        if contract.get("object_type") == "raw_multiome_fastq":
            notes.append(
                "Raw Multiome FASTQs will be assembled into a paired container "
                "(RNA counts + ATAC fragments over a shared barcode namespace) "
                "before joint QC."
            )

        status = "yellow" if findings else "green"
        return _result(
            status,
            _POLICY_FOR_STATUS[status],
            requires_ack_for_dispatch=bool(findings),
            reason=str(meta.get("reason") or ""),
            same_cell_gate_allowed=True,
            checks=checks,
            notes=notes,
            findings=findings,
        )

    @staticmethod
    def _risk_checks(
        exp_context: dict[str, Any],
        design: dict[str, Any],
        contract: dict[str, Any],
    ) -> dict[str, Any]:
        raw_contract = exp_context.get("multiome_contract")
        sources = [
            exp_context,
            design,
            contract,
            raw_contract if isinstance(raw_contract, dict) else {},
        ]

        # Only reported values are judged; an unreported count is not a risk.
        counts = [
            n
            for src in sources
            for key in ("n_shared_cells", "n_cells")
            if (n := _as_count(src.get(key))) is not None
        ]
        n_shared = min(counts) if counts else None
        if n_shared is None:
            cells_status = "not_reported"
        elif n_shared < _INTEGRATION_MIN_SHARED_CELLS:
            cells_status = "warning"
        else:
            cells_status = "pass"
        shared_cells = {
            "status": cells_status,
            "n_shared_cells": n_shared,
            "minimum": _INTEGRATION_MIN_SHARED_CELLS,
        }

        allow_partial = any(
            src.get("allow_partial_overlap") is True for src in sources
        )
        fractions = [
            f for src in sources
            if (f := _as_fraction(src.get("shared_barcode_fraction"))) is not None
        ]
        fraction = min(fractions) if fractions else None
        partial = allow_partial or (
            fraction is not None
            and fraction < _INTEGRATION_MIN_SHARED_BARCODE_FRACTION
        )
        barcode_overlap = {
            "status": "warning" if partial else "pass",
            "allow_partial_overlap": allow_partial,
            "shared_barcode_fraction": fraction,
            "minimum_fraction": _INTEGRATION_MIN_SHARED_BARCODE_FRACTION,
        }

        n_units = max(
            _as_count(contract.get("n_containers")) or 0,
            _as_count(contract.get("n_samples")) or 0,
        )
        is_cohort = contract.get("cohort_merge_required") is True or n_units >= 2
        design_ctx = dict(exp_context)
        manifest = exp_context.get("multiome_fastq_manifest")
        if (
            not design_ctx.get("sample_metadata")
            and isinstance(manifest, dict)
            and isinstance(manifest.get("samples"), list)
        ):
            # Reuse the shared confounding check on raw-FASTQ manifest rows.
            design_ctx["bulk_fastq_manifest"] = {"samples": manifest["samples"]}
        groups = design.get("analysis_groups") or design.get("groups") or {}
        one_donor_per_condition = bool(
            isinstance(groups, dict)
            and len(groups) >= 2
            and all(isinstance(v, (list, tuple, set)) and len(v) <= 1 for v in groups.values())
        )
        confounded = bool(
            is_cohort
            and (
                one_donor_per_condition
                or _batch_condition_confounded(design, design_ctx)
            )
        )
        donor_design = {
            "status": "warning" if confounded else "pass",
            "multi_sample_cohort": is_cohort,
            "n_samples": n_units or None,
            "condition_batch_confounded": confounded,
        }
        return {
            "shared_cells": shared_cells,
            "barcode_overlap": barcode_overlap,
            "donor_design": donor_design,
        }


def _build_integration_readiness(
    exp_context: dict[str, Any],
    multiome_contract: dict[str, Any],
    integration_validation: dict[str, Any] | None,
) -> dict[str, Any]:
    """Evaluate pre-dispatch readiness for paired same-cell integration (ADR-079).

    Intentionally kept separate from per-assay `dispatch` so unpaired/blocked
    cross-modal contracts never block independent single-modality lanes (`scRNA`,
    `scATAC`, `bulk_RNA`, `bulk_ATAC`).
    """
    return IntegrationAuditAgent().audit(
        exp_context,
        multiome_contract,
        integration_validation,
    )


_DIAGONAL_IMPUTATION_STRATEGIES = frozenset(
    {"diagonal", "bridge", "glue", "multivi", "seurat_bridge", "cell_imputation"}
)
_RNA_CONCORDANCE_MODALITIES = frozenset({"scRNA", "bulk_RNA", "bulk_RNA_raw"})
_ATAC_CONCORDANCE_MODALITIES = frozenset({"scATAC", "bulk_ATAC"})


class CrossModalConcordanceAuditAgent:
    """Pre-dispatch readiness auditor for non-same-cell RNA+ATAC concordance (ADR-080).

    Evaluates deterministic contrast- and sample/donor-level concordance readiness
    across ``bulk_RNA``/``bulk_RNA_raw`` + ``bulk_ATAC`` and unpaired cohort
    ``scRNA`` + ``scATAC`` while keeping synthetic cell-level diagonal/bridge
    imputation (``GLUE``, ``MultiVI``, ``Seurat bridge``) blocked fail-closed.
    """

    name = "cross_modal_concordance_audit_agent"

    def audit(
        self,
        exp_context: dict[str, Any],
        multiome_contract: dict[str, Any],
        cards: dict[str, dict[str, Any]],
        concordance_validation: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        meta = concordance_validation or {}
        tier = str(meta.get("level") or "production")
        enabled = bool(meta.get("dispatch_enabled", True))
        modalities = exp_context.get("modalities") or {}

        strategy = (
            str(
                exp_context.get("cross_modal_strategy")
                or exp_context.get("integration_strategy")
                or ""
            )
            .lower()
            .strip()
        )
        diagonal_requested = bool(
            exp_context.get("diagonal_imputation") is True
            or strategy in _DIAGONAL_IMPUTATION_STRATEGIES
        )
        if diagonal_requested:
            finding = _finding(
                "blocking",
                "diagonal_cell_imputation_blocked",
                (
                    "Unpaired single-cell diagonal/bridge cell imputation "
                    "(e.g. GLUE, MultiVI, Seurat bridge) is not validated in ARIA "
                    "and is blocked fail-closed."
                ),
                (
                    "Remove diagonal_imputation / bridge strategy to run audited "
                    "single-modality pipelines plus ADR-080 contrast- and "
                    "donor-level cross-modal concordance."
                ),
                modality="cross_modal",
            )
            return {
                "tier": tier,
                "evidence_tier": "associative",
                "scope": "contrast_and_sample_matched_rna_atac",
                "promotion_adr": str(meta.get("promotion_adr") or "ADR-080"),
                "diagonal_cell_imputation": "blocked",
                "status": "red",
                "dispatch_policy": "blocked",
                "requires_ack_for_dispatch": False,
                "reason": "diagonal_cell_imputation_not_validated",
                "findings": [finding],
            }

        if not enabled:
            return {
                "tier": tier,
                "evidence_tier": "associative",
                "scope": "contrast_and_sample_matched_rna_atac",
                "promotion_adr": str(meta.get("promotion_adr") or "ADR-080"),
                "diagonal_cell_imputation": "blocked",
                "status": "red",
                "dispatch_policy": "blocked",
                "requires_ack_for_dispatch": False,
                "reason": "cross_modal_concordance_dispatch_disabled",
                "findings": [],
            }

        contract_details = multiome_contract.get("contract") or multiome_contract
        if bool(contract_details.get("same_cell") is True) and (
            multiome_contract.get("status") not in {"red", "blocked"}
        ):
            return {
                "tier": tier,
                "evidence_tier": "associative",
                "scope": "contrast_and_sample_matched_rna_atac",
                "promotion_adr": str(meta.get("promotion_adr") or "ADR-080"),
                "diagonal_cell_imputation": "blocked",
                "status": "not_applicable",
                "dispatch_policy": "superseded_by_paired_same_cell_integration",
                "requires_ack_for_dispatch": False,
                "reason": "paired_same_cell_multiome_uses_integration_agent",
                "findings": [],
            }

        active_rna = sorted(
            m
            for m in (set(modalities) & _RNA_CONCORDANCE_MODALITIES)
            if (cards.get(m) or {}).get("status") != "red"
        )
        active_atac = sorted(
            m
            for m in (set(modalities) & _ATAC_CONCORDANCE_MODALITIES)
            if (cards.get(m) or {}).get("status") != "red"
        )
        if active_rna and active_atac:
            return {
                "tier": tier,
                "evidence_tier": "associative",
                "scope": "contrast_and_sample_matched_rna_atac",
                "promotion_adr": str(meta.get("promotion_adr") or "ADR-080"),
                "diagonal_cell_imputation": "blocked",
                "status": "green",
                "dispatch_policy": "auto",
                "requires_ack_for_dispatch": False,
                "rna_modalities": active_rna,
                "atac_modalities": active_atac,
                "reason": "non_same_cell_rna_atac_contrast_and_sample_concordance_ready",
                "findings": [],
            }

        return {
            "tier": tier,
            "evidence_tier": "associative",
            "scope": "contrast_and_sample_matched_rna_atac",
            "promotion_adr": str(meta.get("promotion_adr") or "ADR-080"),
            "diagonal_cell_imputation": "blocked",
            "status": "not_applicable",
            "dispatch_policy": "not_applicable",
            "requires_ack_for_dispatch": False,
            "reason": "requires_validated_rna_and_atac_modalities",
            "findings": [],
        }


def _build_cross_modal_concordance_readiness(
    exp_context: dict[str, Any],
    multiome_contract: dict[str, Any],
    cards: dict[str, dict[str, Any]],
    concordance_validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return CrossModalConcordanceAuditAgent().audit(
        exp_context,
        multiome_contract,
        cards,
        concordance_validation,
    )


def build_capability_matrix(
    exp_context: dict[str, Any],
    *,
    modality_validation: dict[str, dict[str, Any]] | None = None,
    integration_validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    modalities = exp_context.get("modalities") or {}
    validation = modality_validation or {}
    cards: dict[str, dict[str, Any]] = {}

    for modality in modalities:
        meta = validation.get(modality, {})
        level = meta.get("level", "unknown")
        enabled = bool(meta.get("dispatch_enabled", True))
        if not enabled:
            status = "red"
        elif level in {"alpha", "beta", "experimental"}:
            status = "yellow"
        else:
            status = "green"
        cards[modality] = _base_card(
            modality,
            validation_level=level,
            status=status,
            reason=meta.get("reason", ""),
        )
        if not enabled:
            cards[modality]["findings"].append(_finding(
                "blocking",
                "modality_not_dispatch_validated",
                (
                    f"{modality} is {level or 'unvalidated'} and is not "
                    "validated for dispatch."
                ),
                meta.get(
                    "reason",
                    "Do not dispatch this modality until validation closes.",
                ),
                modality=modality,
            ))
        elif status == "yellow":
            cards[modality]["findings"].append(_finding(
                "warning",
                "modality_requires_ack",
                f"{modality} is {level} and requires explicit acknowledgement.",
                "Proceed only if the user accepts the validation limitations.",
                modality=modality,
            ))

    if "scRNA" in modalities:
        cards["scRNA"] = _merge_cards(
            cards.get("scRNA", _base_card("scRNA")),
            ScRNAAuditAgent().audit(exp_context),
        )

    chromatin_auditor = ChromatinAuditAgent()
    for modality in sorted(set(modalities) & _CHROMATIN_MODALITIES):
        if cards.get(modality, {}).get("status") != "red":
            configured_level = (validation.get(modality) or {}).get("level")
            cards[modality] = _merge_cards(
                cards.get(modality, _base_card(modality)),
                chromatin_auditor.audit(
                    exp_context,
                    modality,
                    tier_override=configured_level,
                ),
            )

    if "HiC" in modalities and cards.get("HiC", {}).get("status") != "red":
        configured_hic_level = (validation.get("HiC") or {}).get("level")
        cards["HiC"] = _merge_cards(
            cards.get("HiC", _base_card("HiC")),
            HiCAuditAgent().audit(exp_context, tier_override=configured_hic_level),
        )

    for modality, files in modalities.items():
        cards[modality] = _merge_cards(
            cards.get(modality, _base_card(modality)),
            validate_assay_contract(exp_context, modality, files),
        )

    contracts = {}
    multiome_contract = validate_multiome_contract(exp_context)
    if multiome_contract.get("status") != "not_applicable":
        contracts["multiome"] = multiome_contract
    integration_readiness = _build_integration_readiness(
        exp_context,
        multiome_contract,
        integration_validation,
    )
    cross_modal_concordance_readiness = _build_cross_modal_concordance_readiness(
        exp_context,
        multiome_contract,
        cards,
    )

    preflight = {}
    design_power = assess_design_power(exp_context)
    if design_power.get("status") != "not_applicable":
        preflight["design_power"] = design_power
    reference_integrity = assess_reference_integrity(exp_context)
    if reference_integrity.get("status") != "not_applicable":
        preflight["reference_integrity"] = reference_integrity

    dispatch = {"allowed": [], "requires_ack": [], "blocked": []}
    findings = []
    for modality, card in cards.items():
        policy = card["dispatch_policy"]
        if policy == "blocked":
            dispatch["blocked"].append(modality)
        elif policy == "requires_ack":
            dispatch["requires_ack"].append(modality)
            dispatch["allowed"].append(modality)
        else:
            dispatch["allowed"].append(modality)
        findings.extend(card.get("findings", []))
    for contract in contracts.values():
        findings.extend(contract.get("findings", []))
    findings.extend(integration_readiness.get("findings", []))
    findings.extend(cross_modal_concordance_readiness.get("findings", []))
    for check in preflight.values():
        findings.extend(check.get("findings", []))

    for key in dispatch:
        dispatch[key] = sorted(dispatch[key])

    return {
        "status": (
            "red" if (
                dispatch["blocked"]
                or cross_modal_concordance_readiness.get("status") == "red"
            )
            else "yellow" if (
                dispatch["requires_ack"]
                or integration_readiness.get("requires_ack_for_dispatch")
            )
            else "green"
        ),
        "cards": cards,
        "contracts": contracts,
        "integration_readiness": integration_readiness,
        "cross_modal_concordance_readiness": cross_modal_concordance_readiness,
        "preflight": preflight,
        "dispatch": dispatch,
        "findings": findings,
    }
