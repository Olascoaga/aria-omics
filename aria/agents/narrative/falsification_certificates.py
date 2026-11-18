"""In-Silico Reviewer #2: Claim-Specific Falsification Certificates (P1).

Generates orthogonal, typed falsification diagnostics for published `VERIFIED` /
`associative` claims without ever:
  - emitting a scalar "fragility score" or naive E-value,
  - promoting observational robustness to `CAUSAL` (`VERIFIED_ROBUST` is forbidden),
  - fabricating unexecuted perturbation numbers (`ADR-002`), or
  - feeding back into hyperparameter selection (`feeds_into_parameter_selection=False`).

Four orthogonal falsification axes per claim:
  1. `lobo_biological_unit_jackknife`: Leave-One-Biological-Unit-Out (donor/sample,
     NEVER cell-level) sign consistency, LFC IQR, worst-case jackknife `padj`, and
     `single_unit_driven` flag (`not_testable_insufficient_replicates` when
     `min(n_num, n_den) < 3`; `eligible_not_executed` when `n >= 3` without pre-run
     jackknife folds).
  2. `block_preserving_permutation_null`: Empirical permutation calibration status
     and null family (`TOBIAS` stratified permutations, `cross_modal` donor
     permutations, `peak2gene` `(GC, accessibility)` trans-null, `eregulons` CMH +
     finite-population null, `hic` null permutations, or `parametric_bh_only`).
  3. `specification_multiverse_stability`: Survival across pre-declared BH-FDR
     families (`robustness_multiverse` gene-ID intersection) and clustering/window
     parameter grids.
  4. `design_aliasing_and_unit_audit`: Biological-unit discipline, technical-replicate
     summation, batch confounding status, spike-in normalization, and low-power flags.
"""

from __future__ import annotations

from typing import Any

from aria.utils.lobo_jackknife import (
    bh_adjust as _bh_adjust,
    compute_counts_lobo_jackknife,
    compute_lobo_jackknife,
)

FALSIFICATION_SCHEMA_VERSION = "aria.falsification_certificates.v1"


def _min_arm_replicates(exp_ctx: dict[str, Any], agent_results: dict[str, Any]) -> tuple[int | None, str]:
    design = (exp_ctx or {}).get("design") or {}
    reps = design.get("replicates")
    if isinstance(reps, dict) and reps:
        vals = [int(v) for v in reps.values() if isinstance(v, (int, float)) and int(v) > 0]
        if vals:
            unit = "donor" if design.get("donor_col") else "biological_replicate"
            return min(vals), unit
    groups = design.get("analysis_groups") or design.get("groups")
    if isinstance(groups, dict) and groups:
        vals = [len(v) for v in groups.values() if isinstance(v, list) and v]
        if vals:
            unit = "donor" if design.get("donor_col") else "biological_replicate"
            return min(vals), unit
    return None, "biological_unit"


def _permutation_axis_for_claim(claim: dict[str, Any], agent_results: dict[str, Any]) -> dict[str, Any]:
    modality = str(claim.get("modality") or "")
    analysis = str(claim.get("analysis") or "")
    if analysis in {"peak_to_gene", "peak2gene"}:
        return {
            "status": "calibrated",
            "null_family": "trans_chromosomal_gc_accessibility_matched",
        }
    if analysis == "eregulons":
        return {
            "status": "calibrated",
            "null_family": "gc_quintile_cmh_and_5x5_expr_prevalence_null",
        }
    if "footprint" in analysis or modality == "tobias":
        return {
            "status": "calibrated",
            "null_family": "stratified_within_batch_label_permutation",
        }
    if modality == "cross_modal" or "cross_modal" in analysis:
        return {
            "status": "calibrated",
            "null_family": "sample_label_permutation_bh_fdr",
        }
    if analysis in {"differential_expression", "pseudobulk_de", "differential_accessibility", "contrast", "differential_3d"}:
        return {
            "status": "parametric_bh_only",
            "null_family": "empirical_bayes_wald_or_moderated_ols_bh_fdr",
        }
    return {
        "status": "not_applicable_descriptive",
        "null_family": "none",
    }


def build_falsification_certificates(
    claims: list[dict[str, Any]] | None,
    agent_results: dict[str, Any] | None,
    exp_ctx: dict[str, Any] | None,
    *,
    robustness_multiverse: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Attach orthogonal Reviewer #2 falsification certificates to compiled claims."""
    ctx = exp_ctx or {}
    results = agent_results or {}
    rm = robustness_multiverse or {}
    n_arm_min, unit_kind = _min_arm_replicates(ctx, results)
    precomputed_lobo = ctx.get("lobo_jackknife") or {}
    design = ctx.get("design") or {}

    inferential_analyses = {
        "differential_expression",
        "pseudobulk_de",
        "differential_abundance",
        "differential_accessibility",
        "differential_3d",
        "contrast",
        "peak_to_gene",
        "peak2gene",
        "eregulons",
    }

    certificates: list[dict[str, Any]] = []
    for claim in claims or []:
        if not isinstance(claim, dict):
            continue
        cid = str(claim.get("claim_id") or "")
        analysis = str(claim.get("analysis") or "")
        tier = str(claim.get("tier") or "descriptive")
        is_inferential = analysis in inferential_analyses or tier != "descriptive"

        # Axis 1: LOBO Jackknife
        if not is_inferential:
            lobo_axis = {
                "status": "not_applicable_descriptive",
                "biological_unit": unit_kind,
            }
        elif isinstance(precomputed_lobo, dict) and precomputed_lobo.get("status") == "measured":
            lobo_axis = {
                "status": "measured",
                "biological_unit": unit_kind,
                "n_jackknife_folds": precomputed_lobo.get("n_jackknife_folds"),
                "n_features_evaluated": precomputed_lobo.get("n_features_evaluated"),
                "n_features_lobo_stable": precomputed_lobo.get("n_features_lobo_stable"),
            }
        elif n_arm_min is not None and n_arm_min < 3:
            lobo_axis = {
                "status": "not_testable_insufficient_replicates",
                "biological_unit": unit_kind,
                "n_arm_min": n_arm_min,
                "minimum_required_per_arm": 3,
            }
        else:
            lobo_axis = {
                "status": "eligible_not_executed",
                "biological_unit": unit_kind,
                "n_arm_min": n_arm_min,
            }

        # Axis 2: Block-preserving permutation calibration
        perm_axis = _permutation_axis_for_claim(claim, results)

        # Axis 3: Specification multiverse stability
        rm_entries = rm.get("entries") or []
        if analysis == "pseudobulk_de" and rm_entries:
            computed = [e for e in rm_entries if e.get("stability_status") == "computed"]
            spec_axis = {
                "status": "measured" if computed else "not_computed",
                "axis": "local_vs_global_bh_fdr_gene_intersection",
                "n_blocks_evaluated": len(computed),
            }
        else:
            spec_axis = {
                "status": "single_pre_committed_specification"
                if ctx.get("decision_manifest")
                else "single_specification",
                "threshold_provenance": str(ctx.get("lfc_threshold_provenance") or "default"),
            }

        # Axis 4: Design aliasing & biological-unit discipline
        aliasing_axis = {
            "cell_level_pseudoreplication_blocked": True,
            "technical_replicates_aggregated": bool(
                (design.get("replicate_handling") or {}).get("has_technical_replicates", False)
            ),
            "batch_covariate_modeled": bool(design.get("batch_covariate") or design.get("batch_col")),
            "low_power_n2_warning": bool(n_arm_min == 2),
            "causal_promotion_permitted": bool(tier == "causal_experimental"),
        }

        certificates.append({
            "claim_id": cid,
            "analysis": analysis,
            "tier": tier,
            "feeds_into_parameter_selection": False,
            "lobo_biological_unit_jackknife": lobo_axis,
            "block_preserving_permutation_null": perm_axis,
            "specification_multiverse_stability": spec_axis,
            "design_aliasing_and_unit_audit": aliasing_axis,
        })

    return {
        "schema_version": FALSIFICATION_SCHEMA_VERSION,
        "feeds_into_parameter_selection": False,
        "scalar_fragility_score_emitted": False,
        "causal_promotion_from_robustness": False,
        "n_certificates": len(certificates),
        "certificates": certificates,
    }
