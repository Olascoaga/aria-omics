"""Per-modality readiness cards and capability matrix.

The checks here are deliberately pre-dispatch and data-light. They decide
whether a modality can be dispatched automatically, requires explicit user
acknowledgement, or must be blocked before expensive analysis starts.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from aria.utils.assay_contracts import validate_assay_contract
from aria.utils.multiome_contracts import validate_multiome_contract


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


def _atac_production_readiness(modality: str) -> dict[str, Any]:
    """ATAC readiness, ADR-056-honest.

    The technical C1-C6 blockers are closed, but ATAC stays beta + requires_ack:
    a closed technical blocker list does NOT promote the tier. `production_promoted`
    is therefore False until an explicit promotion ADR clears the governance
    blockers. The old ambiguous `production_ready` field (which meant "no technical
    blockers" and could be misread as "is production") is replaced by the explicit
    `technical_blockers_cleared` + `production_promoted` pair.
    """
    technical_blockers = list(_ATAC_PRODUCTION_BLOCKERS.get(modality, ()))
    return {
        "tier": "beta",
        "requires_ack": True,
        # Production is NEVER inferred from the audit (ADR-056). Flip only via a
        # deliberate promotion ADR that clears every governance blocker.
        "production_promoted": False,
        "technical_blockers_cleared": not technical_blockers,
        "open_blockers": technical_blockers,
        "governance_blockers": list(_ATAC_GOVERNANCE_BLOCKERS),
        "promotion_policy_adr": "ADR-056",
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


def _batch_condition_confounded(design: dict[str, Any]) -> bool:
    groups = design.get("groups") or {}
    batch_map = design.get("batch_map") or {}
    if not groups or not batch_map:
        return False

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

    if len(condition_batches) < 2:
        return False
    seen = {}
    for condition, batches in condition_batches.items():
        key = tuple(sorted(batches))
        seen.setdefault(key, []).append(condition)
    return all(len(conditions) == 1 for conditions in seen.values())


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
        groups = design.get("groups") or {}

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

        if _batch_condition_confounded(design):
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


class ChromatinAuditAgent:
    """Readiness cards for chromatin modalities before full validation closes."""

    def audit(self, exp_context: dict[str, Any], modality: str) -> dict[str, Any]:
        if modality == "scATAC":
            card = _base_card(
                modality,
                agent="chromatin_agent",
                validation_level="beta",
                status="yellow",
                reason=(
                    "scATAC beta dispatch is available only after explicit "
                    "user acknowledgement."
                ),
            )
            card["checks"]["validation"] = {
                "validated": True,
                "stage": "beta",
                # Closed 2026-06-09 (ADR-034): a real headless orchestrator/bus
                # run on the HC11 paired .h5mu drove CP1 -> CP2 -> CP3.5 ack ->
                # live dispatch into aria-chromatin-env (QC -> LSI/clustering ->
                # differential accessibility -> motif enrichment).
                "live_orchestrator_run": "done",
                # De-alpha 2026-06-15 (ADR-048): clustering + motif externally
                # concordant (SnapATAC2/HC11); pseudobulk condition-DA validated
                # vs edgeR/limma/DESeq2 on real biological replicates. scATAC is
                # beta + requires_ack.
                "external_concordance": "done",
            }
            card["checks"]["production_readiness"] = _atac_production_readiness(
                modality
            )
            card["findings"].append(_finding(
                "warning",
                "chromatin_readiness_beta_ack_required",
                (
                    "scATAC is beta: P0-P3 closed, clustering + motif are externally "
                    "concordant, and the pseudobulk condition-DA lane is validated "
                    "against edgeR/limma/DESeq2 on real biological replicates. The "
                    "per-cluster wilcoxon marker path is exploratory (single-sample "
                    "fragile)."
                ),
                (
                    "Proceed only after explicit CP3.5 acknowledgement; use the "
                    "pseudobulk DA lane with biological replicates and treat the "
                    "per-cluster wilcoxon markers as exploratory."
                ),
                modality=modality,
            ))
            return card

        if modality == "bulk_ATAC":
            card = _base_card(
                modality,
                agent="chromatin_agent",
                validation_level="beta",
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
                modality
            )
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


def build_capability_matrix(
    exp_context: dict[str, Any],
    *,
    modality_validation: dict[str, dict[str, Any]] | None = None,
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
            cards[modality] = _merge_cards(
                cards.get(modality, _base_card(modality)),
                chromatin_auditor.audit(exp_context, modality),
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

    for key in dispatch:
        dispatch[key] = sorted(dispatch[key])

    return {
        "status": (
            "red" if dispatch["blocked"]
            else "yellow" if dispatch["requires_ack"]
            else "green"
        ),
        "cards": cards,
        "contracts": contracts,
        "dispatch": dispatch,
        "findings": findings,
    }
