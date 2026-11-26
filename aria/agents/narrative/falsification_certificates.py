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

import math
import re
from typing import Any

from aria.utils.lobo_jackknife import compute_lobo_jackknife

FALSIFICATION_SCHEMA_VERSION = "aria.falsification_certificates.v1"


def _finite_count(value: Any) -> int:
    """Non-bool finite numeric -> int; anything else (None/NaN/bool/str) -> 0."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return 0
    if not math.isfinite(float(value)):
        return 0
    return int(value)


def _min_arm_replicates(exp_ctx: dict[str, Any], agent_results: dict[str, Any]) -> tuple[int | None, str]:
    design = (exp_ctx or {}).get("design") or {}
    reps = design.get("replicates")
    if isinstance(reps, dict) and reps:
        # bool is an int subclass and NaN/inf cannot be int()-ed: reject both.
        # A 0-replicate arm is a real (untestable) design, so it is kept.
        vals = [
            int(v) for v in reps.values()
            if isinstance(v, (int, float)) and not isinstance(v, bool)
            and math.isfinite(float(v)) and v >= 0
        ]
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


def _stage_ran(stage: Any) -> bool:
    if not isinstance(stage, dict) or not stage:
        return False
    if stage.get("ran") is False:
        return False
    st = str(stage.get("status") or "").lower().strip()
    if st in {"not_run", "skipped", "blocked", "error", "failed"}:
        return False
    return stage.get("ran") is True or st in {"success", "ok", "done"}


def _unwrap_integration_dict(agent_results: dict[str, Any]) -> dict[str, Any]:
    integ = (agent_results or {}).get("integration_agent")
    if not isinstance(integ, dict):
        return {}
    findings = integ.get("findings")
    return findings if isinstance(findings, dict) else integ


def _unwrap_chromatin_dict(agent_results: dict[str, Any]) -> dict[str, Any]:
    chrom = (agent_results or {}).get("chromatin_agent")
    if not isinstance(chrom, dict):
        return {}
    try:
        from aria.agents.narrative.narrators.chromatin import unwrap_chromatin_findings
        return unwrap_chromatin_findings(chrom) or {}
    except Exception:
        findings = chrom.get("findings")
        return findings if isinstance(findings, dict) else chrom


def _permutation_axis_for_claim(claim: dict[str, Any], agent_results: dict[str, Any]) -> dict[str, Any]:
    modality = str(claim.get("modality") or "").lower()
    analysis = str(claim.get("analysis") or "").lower()
    results = agent_results or {}

    if analysis in {"peak_to_gene", "peak2gene"}:
        integ = _unwrap_integration_dict(results)
        chrom = _unwrap_chromatin_dict(results)
        chrom_reg = chrom.get("regulatory") if isinstance(chrom.get("regulatory"), dict) else {}
        candidates = (
            [chrom_reg.get("peak_to_gene"), integ.get("peak_to_gene")]
            if any(tok in modality for tok in ("chromatin", "atac", "chip"))
            else [integ.get("peak_to_gene"), chrom_reg.get("peak_to_gene")]
        )
        p2g = next((c for c in candidates if isinstance(c, dict) and c), None)
        bg = p2g.get("background_method") if isinstance(p2g, dict) else None
        bg_ok = (
            isinstance(bg, str)
            and bg.strip().lower() not in {"", "none", "null", "false"}
        )
        if _stage_ran(p2g) and bg_ok:
            return {
                "status": "calibrated",
                "null_family": "trans_chromosomal_gc_accessibility_matched",
            }
        return {
            "status": "not_calibrated",
            "null_family": "none",
        }

    if analysis == "eregulons":
        integ = _unwrap_integration_dict(results)
        ereg = integ.get("eregulons") if isinstance(integ.get("eregulons"), dict) else None
        if _stage_ran(ereg):
            raw_n_perm = ereg.get("n_null_permutations")
            explicit_zero_perm = raw_n_perm is not None and _finite_count(raw_n_perm) <= 0
            p2g_sum = ereg.get("p2g_summary") if isinstance(ereg.get("p2g_summary"), dict) else {}
            bg_match = str(p2g_sum.get("background_match") or "").strip().lower()
            has_null_cal = (
                _finite_count(raw_n_perm) > 0
                or bool(ereg.get("null_calibration"))
                or bool(ereg.get("cistrome_calibration"))
                or (
                    isinstance(ereg.get("calibration"), dict)
                    and bool(ereg["calibration"])
                )
                or bg_match not in {"", "none", "null", "false"}
                or any(
                    isinstance(row, dict)
                    and any(
                        row.get(k) is not None
                        and not isinstance(row.get(k), bool)
                        and isinstance(row.get(k), (int, float))
                        and math.isfinite(float(row.get(k)))
                        for k in ("cistrome_cmh_padj", "cistrome_cmh_z", "cistrome_log2_fe")
                    )
                    for row in (ereg.get("top_eregulons") or [])
                )
            )
            if has_null_cal and not explicit_zero_perm:
                return {
                    "status": "calibrated",
                    "null_family": "gc_quintile_cmh_and_5x5_expr_prevalence_null",
                }
        return {
            "status": "not_calibrated",
            "null_family": "none",
        }

    if "footprint" in analysis or modality == "tobias":
        chrom = _unwrap_chromatin_dict(results)
        foot = chrom.get("footprinting") if isinstance(chrom.get("footprinting"), dict) else None
        if _stage_ran(foot):
            diff_sum = (
                foot.get("differential_summary")
                if isinstance(foot.get("differential_summary"), dict) else {}
            )
            inference = (
                diff_sum.get("inference")
                if isinstance(diff_sum.get("inference"), dict) else {}
            )
            inf_status = str(inference.get("status") or "").lower().strip()
            fdr_controlled = (
                (diff_sum.get("ranking_basis") or {}).get("fdr_controlled")
                if isinstance(diff_sum.get("ranking_basis"), dict) else None
            )
            null_ctrl = (
                inference.get("null_label_controls")
                or diff_sum.get("null_label_controls")
                or foot.get("null_label_controls")
                or {}
            )
            n_perms = (
                _finite_count(null_ctrl.get("n_permutations"))
                if isinstance(null_ctrl, dict) else 0
            )
            if inf_status == "success" and fdr_controlled is not False:
                if n_perms > 0:
                    return {
                        "status": "calibrated",
                        "null_family": "stratified_within_batch_label_permutation",
                    }
                return {
                    "status": "parametric_bh_only",
                    "null_family": "replicate_welch_t_bh_without_label_permutation",
                }
        return {
            "status": "not_calibrated_descriptive_only",
            "null_family": "none",
        }
    if modality == "cross_modal" or "cross_modal" in analysis:
        # The sample-label permutation null exists ONLY for sample-matched link
        # claims when the sample-matched link layer actually executed; contrast-
        # level DE×DA, TF, and pathway layers use parametric Fisher/Spearman/BH.
        from aria.agents.narrative.run_ledger import cross_modal_findings

        claim_id = str(claim.get("claim_id") or "").lower()
        layer = str(
            claim.get("cross_modal_layer")
            or (claim.get("metadata") or {}).get("cross_modal_layer")
            or ""
        ).lower()
        is_sample_links_claim = (
            "sample_links" in claim_id
            or "sample_matched" in claim_id
            or layer == "sample_matched_links"
        )
        links = cross_modal_findings(agent_results or {}).get("sample_matched_links")
        if is_sample_links_claim and isinstance(links, dict) and links.get("ran") is True:
            return {
                "status": "calibrated",
                "null_family": "sample_label_permutation_bh_fdr",
            }
        return {
            "status": "parametric_bh_only",
            "null_family": "bh_fdr_without_executed_sample_permutation",
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


def _safe_token(value: Any) -> str:
    # Same normalization the narrators use for block ids (``_safe_id``).
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or ""))
    return text.strip("_")


# Block-id prefixes whose remainder names ONE contrast (and group).
_CONTRAST_SCOPED_PREFIXES = ("bulk.contrast.", "scrna.pseudobulk.")
_CLAIM_SCOPE_FIELDS = ("contrast", "comparison", "group", "cluster")


def _claim_has_scope_hint(claim: dict[str, Any]) -> bool:
    if any(str(claim.get(k) or "").strip() for k in _CLAIM_SCOPE_FIELDS):
        return True
    cid = str(claim.get("claim_id") or "")
    return cid.startswith(_CONTRAST_SCOPED_PREFIXES)


def _entry_matches_claim(
    claim: dict[str, Any], group: Any, name: Any,
) -> bool:
    """True when the claim identifies this (group, contrast) entry."""
    hay = "." + str(claim.get("claim_id") or "") + "."
    explicit_names = {
        _safe_token(claim.get(k)) for k in ("contrast", "comparison")
        if str(claim.get(k) or "").strip()
    }
    explicit_groups = {
        _safe_token(claim.get(k)) for k in ("group", "cluster")
        if str(claim.get(k) or "").strip()
    }
    name_tok = _safe_token(name)
    if not name_tok:
        return False
    if explicit_names:
        name_ok = name_tok in explicit_names
    else:
        name_ok = f".{name_tok}." in hay
    if not name_ok:
        return False
    if group is None:
        return True
    group_tok = _safe_token(group)
    if explicit_groups:
        return group_tok in explicit_groups
    return f".{group_tok}." in hay


def _contrast_name(contrast: dict[str, Any], key: Any = None) -> Any:
    for field in ("name", "comparison"):
        if contrast.get(field):
            return contrast.get(field)
    test, ref = contrast.get("test"), contrast.get("reference")
    if test is not None and ref is not None:
        return f"{test}_vs_{ref}"
    num, den = contrast.get("numerator"), contrast.get("denominator")
    if num is not None and den is not None:
        return f"{num}_vs_{den}"
    return key


def _extract_lobo_from_agent_results(
    claim: dict[str, Any], agent_results: dict[str, Any],
) -> list[dict[str, Any]]:
    """Read executed DE/DA summaries from the claim's modality AND contrast.

    When the claim identifies a specific contrast/group (``bulk.contrast.<name>``,
    ``scrna.pseudobulk.<group>.<comparison>`` or explicit ``contrast`` /
    ``group`` fields) only the matching entries are used; a hinted claim that
    matches nothing yields ``[]`` (never another contrast's diagnostics).
    Modality-wide claims fall back to all contrasts of the modality.
    """
    if claim.get("analysis") not in {
        "differential_expression", "pseudobulk_de",
        "differential_accessibility", "contrast",
    }:
        return []
    modality = str(claim.get("modality") or "").lower()
    modality = modality.replace("_", "").replace("-", "").replace(" ", "")
    # (group, contrast_name, contrast_dict)
    entries: list[tuple[Any, Any, Any]] = []
    if modality in {"bulkrna", "bulkrnaseq"}:
        result = agent_results.get("bulk_rna_agent") or {}
        findings = result.get("findings") or result
        for c in findings.get("contrasts") or []:
            name = _contrast_name(c) if isinstance(c, dict) else None
            entries.append((None, name, c))
    elif modality in {"scrna", "scrnaseq"}:
        from aria.agents.narrative.scrna import unwrap_scrna_findings

        result = agent_results.get("scrna_agent") or {}
        findings = unwrap_scrna_findings(result)
        pb = findings.get("pseudobulk_de") or {}
        for gkey, group in (pb.get("per_group") or {}).items():
            if isinstance(group, dict):
                for ckey, c in (group.get("per_comparison") or {}).items():
                    entries.append((gkey, ckey, c))
    elif modality in {
        "chromatin", "scatac", "bulkatac", "atac", "atacseq",
        "chip", "chipseq", "cutandrun", "cutandtag", "cut&run", "cut&tag",
    }:
        from aria.agents.narrative.narrators.chromatin import unwrap_chromatin_findings

        findings = unwrap_chromatin_findings(agent_results.get("chromatin_agent") or {})
        seen = set()

        def _add(group: Any, comparisons: Any) -> None:
            # Aliases can point at the same comparison list in a live result.
            if not comparisons or id(comparisons) in seen:
                return
            seen.add(id(comparisons))
            items = (comparisons.items() if isinstance(comparisons, dict)
                     else enumerate(comparisons))
            for key, c in items:
                if not isinstance(c, dict):
                    continue
                name = _contrast_name(
                    c, key if isinstance(comparisons, dict) else None)
                entries.append((group, name, c))

        for key in ("diffacc", "differential_accessibility", "bulk_diffacc"):
            da = findings.get(key) or {}
            if not isinstance(da, dict):
                continue
            for container in (da.get("pseudobulk") or {}, da):
                if not isinstance(container, dict):
                    continue
                _add(None, container.get("comparisons"))
                # scATAC per-group (cluster / cell-type) DA layout.
                per_group = container.get("per_group")
                if isinstance(per_group, dict):
                    for gkey, ginfo in per_group.items():
                        if not isinstance(ginfo, dict):
                            continue
                        _add(gkey, ginfo.get("comparisons"))
                        _add(gkey, ginfo.get("per_comparison"))
    if _claim_has_scope_hint(claim):
        entries = [
            e for e in entries if _entry_matches_claim(claim, e[0], e[1])
        ]
    return [
        c["lobo_jackknife"] for _, _, c in entries
        if isinstance(c, dict) and isinstance(c.get("lobo_jackknife"), dict)
    ]


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
    precomputed_lobo = ctx.get("lobo_jackknife")
    # Only a real precomputed LOBO summary (a dict with a status) shadows the
    # live agent-result extraction; None / {} must not suppress it.
    has_precomputed_lobo = (
        isinstance(precomputed_lobo, dict)
        and precomputed_lobo.get("status") not in (None, "")
    )
    if not isinstance(precomputed_lobo, dict):
        precomputed_lobo = {}
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

        extracted_lobo = (
            _extract_lobo_from_agent_results(claim, results)
            if not has_precomputed_lobo and is_inferential else []
        )
        measured_lobo = [s for s in extracted_lobo if s.get("status") == "measured"]
        error_lobo = [s for s in extracted_lobo if s.get("status") == "error"]

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
        elif isinstance(precomputed_lobo, dict) and precomputed_lobo.get("status") == "error":
            lobo_axis = {
                "status": "error",
                "biological_unit": unit_kind,
                "error_details": [str(precomputed_lobo.get("details") or "LOBO jackknife failed")],
            }
        elif measured_lobo:
            n_evaluated = sum(_finite_count(s.get("n_features_evaluated")) for s in measured_lobo)
            n_stable = sum(_finite_count(s.get("n_features_lobo_stable")) for s in measured_lobo)
            lobo_axis = {
                "status": "measured",
                "biological_unit": unit_kind,
                "n_contrasts_measured": len(measured_lobo),
                "n_jackknife_folds": max(_finite_count(s.get("n_jackknife_folds")) for s in measured_lobo),
                "n_features_evaluated": n_evaluated,
                "n_features_lobo_stable": n_stable,
                "n_single_unit_driven": sum(_finite_count(s.get("n_single_unit_driven")) for s in measured_lobo),
                "lobo_stable_fraction": round(n_stable / n_evaluated, 4) if n_evaluated else None,
            }
            if error_lobo:
                lobo_axis["n_contrasts_errored"] = len(error_lobo)
                lobo_axis["error_details"] = [
                    str(s.get("details") or "LOBO jackknife failed") for s in error_lobo
                ]
        elif error_lobo:
            lobo_axis = {
                "status": "error",
                "biological_unit": unit_kind,
                "n_contrasts_errored": len(error_lobo),
                "error_details": [
                    str(s.get("details") or "LOBO jackknife failed") for s in error_lobo
                ],
            }
        elif extracted_lobo and all(
            s.get("status") == "not_testable_insufficient_replicates" for s in extracted_lobo
        ):
            lobo_axis = {
                "status": "not_testable_insufficient_replicates",
                "biological_unit": unit_kind,
                "minimum_required_per_arm": 3,
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

        # Axis 4: Design aliasing & biological-unit discipline (D7-008: read
        # canonical DesignAgent fields through resolve_design_contract).
        from aria.utils.design_matrix import resolve_design_contract

        d_contract = resolve_design_contract(design, ctx)
        aliasing_axis = {
            "cell_level_pseudoreplication_blocked": True,
            "technical_replicates_aggregated": bool(
                d_contract["has_technical_replicates"]
            ),
            "batch_covariate_modeled": bool(d_contract["batch_covariate"]),
            "low_power_n2_warning": bool(n_arm_min is not None and n_arm_min <= 2),
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
