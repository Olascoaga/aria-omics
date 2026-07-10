"""Pre-flight design and power advisor for DE/DA requests.

This is an advisory gate, not an inferential engine. It uses the confirmed
sample-level design to flag under-powered, imbalanced, or confounded contrasts
before expensive DE/DA dispatch. Power estimates are approximate and explicitly
state their assumptions so they never masquerade as measured model power.
"""

from __future__ import annotations

import math
from copy import deepcopy
from typing import Any


DEFAULT_TARGET_LOG2FC = 1.0
DEFAULT_ALPHA = 0.05
DEFAULT_POWER_TARGET = 0.80
DEFAULT_ASSUMED_LOG2_SIGMA = 1.0


def assess_design_power(exp_context: dict[str, Any]) -> dict[str, Any]:
    """Assess pre-flight design power from `design`/`inferred_design`.

    Returns a JSON-serializable block:
    `{status, assumptions, contrasts, checks, findings}`.
    """
    design = _combined_design(exp_context)
    groups = _normalize_groups(
        design.get("analysis_groups") or design.get("groups") or {}
    )
    comparisons = _normalize_comparisons(
        design.get("comparisons")
        or (design.get("pseudobulk") or {}).get("comparisons")
        or design.get("plan_contrasts")
        or exp_context.get("comparisons")
        or [],
        groups,
    )
    if not groups or len(groups) < 2:
        return _empty("not_applicable", "no_condition_groups")

    if not comparisons:
        comparisons = _default_comparisons(groups)

    assumptions = {
        "target_log2fc": float(
            exp_context.get("target_log2fc")
            or design.get("target_log2fc")
            or DEFAULT_TARGET_LOG2FC
        ),
        "alpha": float(exp_context.get("global_padj") or DEFAULT_ALPHA),
        "power_target": DEFAULT_POWER_TARGET,
        "assumed_log2_sigma": float(
            exp_context.get("assumed_log2_sigma")
            or design.get("assumed_log2_sigma")
            or DEFAULT_ASSUMED_LOG2_SIGMA
        ),
        "method": "normal_approx_two_sample_log2_scale",
        "advisory_only": True,
    }

    checks = {
        "replicates_per_condition": {
            condition: len(samples)
            for condition, samples in groups.items()
        },
        "balance_ratio": _balance_ratio(groups),
        "batch_condition_confounding": _batch_confounding(design, groups),
    }

    contrast_results = [
        _assess_contrast(test, ref, groups, assumptions)
        for test, ref in comparisons
    ]
    findings = _findings(checks, contrast_results)
    status = (
        "red" if any(f["severity"] == "blocking" for f in findings)
        else "yellow" if findings
        else "green"
    )
    return {
        "status": status,
        "assumptions": assumptions,
        "checks": checks,
        "contrasts": contrast_results,
        "findings": findings,
    }


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


def _normalize_groups(groups: dict[str, Any]) -> dict[str, list[str]]:
    normalized: dict[str, list[str]] = {}
    for condition, samples in (groups or {}).items():
        if samples is None:
            normalized[str(condition)] = []
        elif isinstance(samples, (list, tuple, set)):
            normalized[str(condition)] = [
                str(sample) for sample in samples if str(sample).strip()
            ]
        else:
            normalized[str(condition)] = [str(samples)]
    return normalized


def _normalize_comparisons(
    comparisons: list[Any],
    groups: dict[str, list[str]],
) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for comp in comparisons or []:
        test = ref = None
        if isinstance(comp, dict):
            test = (
                comp.get("test")
                or comp.get("case")
                or comp.get("numerator")
                or comp.get("condition")
            )
            ref = (
                comp.get("reference")
                or comp.get("control")
                or comp.get("denominator")
                or comp.get("baseline")
            )
        elif isinstance(comp, (list, tuple)) and len(comp) >= 2:
            test, ref = comp[0], comp[1]
        if test is None or ref is None:
            continue
        test_s, ref_s = str(test), str(ref)
        if test_s in groups and ref_s in groups and test_s != ref_s:
            out.append((test_s, ref_s))
    return out


def _default_comparisons(groups: dict[str, list[str]]) -> list[tuple[str, str]]:
    names = sorted(groups)
    lower = {name.lower(): name for name in names}
    for key in ("control", "ctrl", "wt", "healthy", "untreated", "young"):
        if key in lower:
            ref = lower[key]
            return [(name, ref) for name in names if name != ref]
    return [
        (b, a)
        for i, a in enumerate(names)
        for b in names[i + 1:]
    ]


def _assess_contrast(
    test: str,
    ref: str,
    groups: dict[str, list[str]],
    assumptions: dict[str, float],
) -> dict[str, Any]:
    n_test = len(groups.get(test, []))
    n_ref = len(groups.get(ref, []))
    se = _standard_error(n_test, n_ref, assumptions["assumed_log2_sigma"])
    if se is None:
        power = 0.0
        mde = None
    else:
        z_alpha = 1.96 if assumptions["alpha"] <= 0.05 else _z_for_alpha(assumptions["alpha"])
        z_power = 0.84
        z_effect = assumptions["target_log2fc"] / se
        power = _normal_cdf(z_effect - z_alpha)
        mde = (z_alpha + z_power) * se

    return {
        "comparison": {"test": test, "reference": ref},
        "n_replicates": {"test": n_test, "reference": n_ref},
        "min_replicates": min(n_test, n_ref),
        "balance_ratio": _pair_balance_ratio(n_test, n_ref),
        "power_estimate_at_target_log2fc": round(float(power), 3),
        "minimum_detectable_log2fc_at_80_power": (
            round(float(mde), 3) if mde is not None else None
        ),
        "assessment": _contrast_assessment(n_test, n_ref, power, mde),
    }


def _standard_error(n_a: int, n_b: int, sigma: float) -> float | None:
    if n_a < 1 or n_b < 1:
        return None
    return float(sigma) * math.sqrt((1.0 / n_a) + (1.0 / n_b))


def _contrast_assessment(
    n_test: int,
    n_ref: int,
    power: float,
    mde: float | None,
) -> str:
    min_n = min(n_test, n_ref)
    if min_n < 2:
        return "unsupported"
    if min_n < 3:
        return "low_power_supported_with_caveat"
    if power < 0.5:
        return "exploratory_low_power"
    if mde is not None and mde > 1.5:
        return "large_effects_only"
    return "adequate_for_target_assumption"


def _findings(
    checks: dict[str, Any],
    contrasts: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    findings = []
    confounding = checks.get("batch_condition_confounding") or {}
    if confounding.get("confounded"):
        findings.append(_finding(
            "blocking",
            "design_power_batch_condition_confounding",
            (
                "Condition is completely confounded with batch in the confirmed "
                f"design: {confounding.get('condition_batches')}."
            ),
            (
                "Do not run condition DE/DA from this design. Add a balanced "
                "batch design, change the contrast, or analyze only descriptive QC."
            ),
        ))

    for contrast in contrasts:
        assessment = contrast.get("assessment")
        comp = contrast.get("comparison") or {}
        label = f"{comp.get('test')} vs {comp.get('reference')}"
        reps = contrast.get("n_replicates") or {}
        if assessment == "unsupported":
            findings.append(_finding(
                "blocking",
                "design_power_unsupported_replicates",
                (
                    f"Contrast {label} has fewer than two biological replicates "
                    f"in at least one condition: {reps}."
                ),
                (
                    "Do not run inferential DE/DA for this contrast. Add "
                    "biological replicates or restrict the analysis to descriptive outputs."
                ),
            ))
        elif assessment == "low_power_supported_with_caveat":
            findings.append(_finding(
                "warning",
                "design_power_n2_low_power",
                (
                    f"Contrast {label} has n=2 in at least one condition: {reps}; "
                    "supported analyses must be treated as low-power."
                ),
                (
                    "Proceed only with an explicit low-power caveat. Three or "
                    "more biological replicates per condition is the preferred "
                    "publication-grade floor."
                ),
            ))
        elif assessment in {"exploratory_low_power", "large_effects_only"}:
            findings.append(_finding(
                "warning",
                "design_power_low_power_or_large_effects_only",
                (
                    f"Contrast {label} has approximate power "
                    f"{contrast.get('power_estimate_at_target_log2fc')} for the "
                    "target log2FC under stated assumptions; detectable effects "
                    f"at 80% power are about "
                    f"{contrast.get('minimum_detectable_log2fc_at_80_power')} log2FC."
                ),
                (
                    "Treat negative DE/DA results as inconclusive for smaller "
                    "effects, or add replicates before running confirmatory analysis."
                ),
            ))

        if (contrast.get("balance_ratio") or 1.0) < 0.5:
            findings.append(_finding(
                "warning",
                "design_power_imbalanced_contrast",
                (
                    f"Contrast {label} is strongly imbalanced: {reps}."
                ),
                (
                    "Prefer balanced biological replicate counts, or disclose "
                    "reduced power and sensitivity to outliers."
                ),
            ))
    return findings


def _batch_confounding(
    design: dict[str, Any],
    groups: dict[str, list[str]],
) -> dict[str, Any]:
    batch_map = design.get("batch_map") or {}
    if not batch_map:
        return {"confounded": False, "reason": "no_batch_map"}

    condition_batches: dict[str, list[str]] = {}
    for condition, samples in groups.items():
        batches = set()
        for sample in samples:
            for stem, batch in batch_map.items():
                if str(stem) in str(sample):
                    batches.add(str(batch))
                    break
        condition_batches[condition] = sorted(batches)

    populated = {
        condition: batches
        for condition, batches in condition_batches.items()
        if batches
    }
    if len(populated) < 2:
        return {
            "confounded": False,
            "condition_batches": condition_batches,
            "reason": "insufficient_batch_labels",
        }

    batch_sets = {tuple(batches) for batches in populated.values()}
    all_single_batch = all(len(batches) == 1 for batches in populated.values())
    distinct_batches = len(batch_sets) == len(populated)
    return {
        "confounded": bool(all_single_batch and distinct_batches),
        "condition_batches": condition_batches,
    }


def _balance_ratio(groups: dict[str, list[str]]) -> float | None:
    counts = [len(samples) for samples in groups.values() if samples]
    if not counts:
        return None
    return round(min(counts) / max(counts), 3)


def _pair_balance_ratio(n_a: int, n_b: int) -> float | None:
    if not n_a or not n_b:
        return None
    return round(min(n_a, n_b) / max(n_a, n_b), 3)


def _normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _z_for_alpha(alpha: float) -> float:
    # A small lookup avoids a scipy dependency. Values are two-sided.
    if alpha <= 0.001:
        return 3.29
    if alpha <= 0.01:
        return 2.58
    if alpha <= 0.05:
        return 1.96
    if alpha <= 0.10:
        return 1.64
    return 1.28


def _finding(
    severity: str,
    check: str,
    message: str,
    recommendation: str,
) -> dict[str, str]:
    return {
        "severity": severity,
        "check": check,
        "message": message,
        "recommendation": recommendation,
        "modality": "design",
    }


def _empty(status: str, reason: str) -> dict[str, Any]:
    return {
        "status": status,
        "reason": reason,
        "assumptions": {},
        "checks": {},
        "contrasts": [],
        "findings": [],
    }
