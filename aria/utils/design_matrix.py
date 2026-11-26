"""Design-matrix sanity checks before DE model fitting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

_CONTRAST_TEST_KEYS = ("test", "numerator", "case", "contrast", "condition")
_CONTRAST_REF_KEYS = ("reference", "ref", "denominator", "control", "baseline")
_CONTRAST_VS_SEPS = ("_vs_", " vs ", "_VS_", " VS ")


def _extract_single_comparison(
    comp: Any, *, keep_name: bool = False
) -> dict[str, str] | None:
    if isinstance(comp, str):
        text = comp.strip()
        for sep in _CONTRAST_VS_SEPS:
            if sep in text:
                left, right = text.split(sep, 1)
                left, right = left.strip(), right.strip()
                if left and right:
                    return {"test": left, "reference": right}
        return None
    if isinstance(comp, dict):
        test = next(
            (
                comp[k]
                for k in _CONTRAST_TEST_KEYS
                if comp.get(k) is not None and str(comp.get(k)).strip() != ""
            ),
            None,
        )
        ref = next(
            (
                comp[k]
                for k in _CONTRAST_REF_KEYS
                if comp.get(k) is not None and str(comp.get(k)).strip() != ""
            ),
            None,
        )
        if test is None or ref is None:
            return None
        rec = {"test": str(test).strip(), "reference": str(ref).strip()}
        if keep_name and comp.get("name") not in (None, ""):
            rec["name"] = str(comp["name"]).strip()
        return rec
    if isinstance(comp, (list, tuple)) and len(comp) >= 2:
        t_val, r_val = comp[0], comp[1]
        if (
            t_val is not None
            and str(t_val).strip() != ""
            and r_val is not None
            and str(r_val).strip() != ""
        ):
            return {"test": str(t_val).strip(), "reference": str(r_val).strip()}
    return None


def normalize_comparisons(
    raw: Any,
    *,
    valid_groups: Any = None,
    keep_name: bool = False,
) -> list[dict[str, str]]:
    """Canonical contrast/comparison normalizer across ARIA design, power, and DE modules.

    Accepts any combination of:
      - dicts using ARIA's contrast vocabulary
        (``test``/``numerator``/``case``/``contrast``/``condition`` vs
        ``reference``/``ref``/``denominator``/``control``/``baseline``),
      - 2-element ``[test, reference]`` / ``(test, reference)`` pairs,
      - ``"A_vs_B"`` / ``"A vs B"`` strings.

    When ``valid_groups`` is provided, only non-self contrasts whose levels are
    both present in ``valid_groups`` are retained.
    """
    if raw in (None, "", []):
        return []
    out: list[dict[str, str]] = []
    if isinstance(raw, (str, dict)):
        rec = _extract_single_comparison(raw, keep_name=keep_name)
        if rec is not None:
            out.append(rec)
    elif isinstance(raw, (list, tuple)):
        for item in raw:
            rec = _extract_single_comparison(item, keep_name=keep_name)
            if rec is not None:
                out.append(rec)
    if valid_groups is not None:
        vg = {str(g) for g in valid_groups}
        out = [
            c
            for c in out
            if c["test"] in vg and c["reference"] in vg and c["test"] != c["reference"]
        ]
    return out


def resolve_design_contract(
    design: dict | None,
    exp_ctx: dict | None = None,
) -> dict[str, Any]:
    """Normalize the confirmed design into the fields the differential lanes need.

    This is the single RNA/ATAC design contract (E7). The confirmed design dict
    from ``DesignAgent._build_design`` (``main_factor`` / ``batch_covariate`` /
    ``pseudobulk.{condition_col,replicate_col}`` / ``covariates``) is the source of
    truth; ``exp_ctx`` supplies explicit run-level overrides. Both bulk and
    single-cell, RNA and ATAC, resolve condition/replicate/covariates through here
    instead of re-deriving them from scattered keys with divergent fallbacks and
    hard-coded literal defaults.

    Returns ``condition_col`` / ``replicate_col`` / ``covariates`` (the confirmed
    batch covariate plus any extra covariates, deduped and order-preserving, never
    the condition or replicate columns) / ``batch_covariate`` / ``comparisons`` /
    ``pseudobulk``. Any field with no confirmed value is ``None`` (or ``[]`` for
    covariates), never a fabricated literal such as ``"condition"``.
    """
    design = design or {}
    exp_ctx = exp_ctx or {}
    pseudobulk = design.get("pseudobulk") or {}

    def _first(*values):
        for value in values:
            if value is not None and str(value).strip():
                return str(value)
        return None

    condition_col = _first(
        exp_ctx.get("condition_col"),
        pseudobulk.get("condition_col"),
        design.get("condition_col"),
        design.get("main_factor"),
    )
    replicate_col = _first(
        exp_ctx.get("replicate_col"),
        exp_ctx.get("replicate_column"),
        pseudobulk.get("replicate_col"),
        design.get("replicate_col"),
    )
    batch_covariate = _first(
        exp_ctx.get("batch_covariate"),
        exp_ctx.get("batch_col"),
        exp_ctx.get("batch_factor"),
        design.get("batch_covariate"),
        design.get("batch_col"),
        design.get("batch_factor"),
    )

    covariates: list[str] = []
    if batch_covariate:
        covariates.append(batch_covariate)
    for source in (design.get("covariates"), exp_ctx.get("covariates")):
        for cov in (source or []):
            if cov and str(cov).strip():
                covariates.append(str(cov))
    # Dedupe, order-preserving, and never fold the condition or replicate columns
    # into the covariate set (they are modelled as the factor / donor block).
    seen: set[str] = set()
    clean_covariates: list[str] = []
    for cov in covariates:
        if cov in (condition_col, replicate_col) or cov in seen:
            continue
        seen.add(cov)
        clean_covariates.append(cov)

    raw_comparisons = exp_ctx.get("comparisons")
    if raw_comparisons is None:
        raw_comparisons = design.get("comparisons")
    if raw_comparisons is None and isinstance(design.get("contrasts"), list):
        norm_comps = normalize_comparisons(design["contrasts"])
        if norm_comps:
            raw_comparisons = norm_comps

    design_formula = _first(
        exp_ctx.get("design_formula"),
        design.get("design_formula"),
        design.get("formula"),
        exp_ctx.get("formula"),
    )

    rep_handling = design.get("replicate_handling")
    if not isinstance(rep_handling, dict):
        rep_handling = {}
    rep_mode = str(rep_handling.get("mode") or "").strip().lower()
    exp_unit = str(design.get("experimental_unit") or "").strip().lower()
    is_paired = bool(
        exp_ctx.get("paired_replicates")
        or design.get("paired")
        or rep_mode == "paired_biological_block"
        or exp_unit == "paired_biological_unit"
    )
    has_technical_replicates = bool(rep_mode == "technical_aggregate")

    # Sample-to-covariate value mappings from DesignAgent (paired_covariate_values,
    # replicate_handling.sample_to_pair, sample_sheet) plus legacy batch_map /
    # covariate_map shapes.
    covariate_values: dict[str, dict[str, str]] = {}

    def _merge_cov_map(cov_name: str | None, mapping: Any) -> None:
        if not cov_name or not isinstance(mapping, dict):
            return
        bucket = covariate_values.setdefault(str(cov_name), {})
        for k, v in mapping.items():
            ks = str(k).strip()
            if ks and v not in (None, ""):
                bucket.setdefault(ks, str(v).strip())

    paired_cov_vals = design.get("paired_covariate_values")
    if isinstance(paired_cov_vals, dict):
        for cov_name, s_map in paired_cov_vals.items():
            _merge_cov_map(str(cov_name).strip(), s_map)

    if rep_handling.get("pair_covariate") and isinstance(
        rep_handling.get("sample_to_pair"), dict
    ):
        _merge_cov_map(
            str(rep_handling["pair_covariate"]).strip(),
            rep_handling["sample_to_pair"],
        )

    cov_map_legacy = design.get("covariate_map")
    if isinstance(cov_map_legacy, dict):
        for cov_name, s_map in cov_map_legacy.items():
            _merge_cov_map(str(cov_name).strip(), s_map)

    legacy_batch_map = design.get("batch_map") or design.get("batches")
    if isinstance(legacy_batch_map, dict) and legacy_batch_map:
        _merge_cov_map(batch_covariate or "batch", legacy_batch_map)

    inferred_design = exp_ctx.get("inferred_design")
    inferred_design = inferred_design if isinstance(inferred_design, dict) else {}
    target_covs = list(dict.fromkeys(
        ([batch_covariate] if batch_covariate else []) + clean_covariates
    ))
    for sheet in (
        design.get("sample_sheet"),
        exp_ctx.get("sample_sheet"),
        inferred_design.get("sample_sheet"),
    ):
        if not isinstance(sheet, list) or not target_covs:
            continue
        for row in sheet:
            if not isinstance(row, dict):
                continue
            sid = str(row.get("sample") or row.get("id") or "").strip()
            if not sid:
                continue
            for cov_name in target_covs:
                val = row.get(cov_name)
                if val not in (None, ""):
                    covariate_values.setdefault(cov_name, {}).setdefault(
                        sid, str(val).strip()
                    )

    confounded_covariates: list[str] = []
    for src in (
        design.get("confounded_covariates"),
        inferred_design.get("confounded_covariates"),
    ):
        for cov in (src or []):
            cs = str(cov).strip()
            if cs and cs not in confounded_covariates:
                confounded_covariates.append(cs)

    return {
        "condition_col": condition_col,
        "replicate_col": replicate_col,
        "covariates": clean_covariates,
        "batch_covariate": batch_covariate,
        "comparisons": raw_comparisons,
        "pseudobulk": {
            "condition_col": condition_col,
            "replicate_col": replicate_col,
        },
        "design_formula": design_formula,
        "is_paired": is_paired,
        "has_technical_replicates": has_technical_replicates,
        "covariate_values": covariate_values,
        "confounded_covariates": confounded_covariates,
    }


@dataclass(frozen=True)
class DesignIssue:
    severity: str
    check: str
    message: str
    recommendation: str

    def as_dict(self) -> dict[str, str]:
        return {
            "severity": self.severity,
            "check": self.check,
            "message": self.message,
            "recommendation": self.recommendation,
        }


def validate_design_matrix(
    metadata,
    condition_col: str,
    covariates: list[str] | None = None,
    *,
    min_replicates_per_condition: int = 2,
    min_replicates_per_design_cell: int = 2,
) -> dict[str, Any]:
    """
    Validate a sample-level DE design before handing it to DESeq2.

    The function is intentionally dependency-light (pandas/numpy only) and
    returns JSON-serializable issues so scripts can surface failures instead of
    letting DESeq2 raise opaque rank/contrast errors.
    """
    import numpy as np
    import pandas as pd

    covariates = [c for c in (covariates or []) if c and c != condition_col]
    issues: list[DesignIssue] = []

    if metadata is None or len(metadata) == 0:
        return _result("blocking", [DesignIssue(
            "blocking",
            "empty_metadata",
            "No sample metadata are available for design validation.",
            "Provide sample-level metadata before running differential expression.",
        )])

    missing = [c for c in [condition_col, *covariates] if c not in metadata.columns]
    if missing:
        return _result("blocking", [DesignIssue(
            "blocking",
            "missing_design_columns",
            f"Design columns missing from metadata: {', '.join(missing)}.",
            "Confirm condition and covariate columns before running DESeq2.",
        )])

    design_cols = [condition_col, *covariates]
    meta = metadata.loc[:, design_cols].copy()
    meta = meta.replace({"": pd.NA, "NA": pd.NA, "nan": pd.NA})
    before = len(meta)
    meta = meta.dropna(axis=0, how="any")
    if len(meta) < before:
        issues.append(DesignIssue(
            "warning",
            "metadata_rows_dropped",
            f"{before - len(meta)} sample(s) have missing design values.",
            "Inspect metadata; samples with missing design values cannot support DE.",
        ))
    if len(meta) == 0:
        issues.append(DesignIssue(
            "blocking",
            "empty_complete_cases",
            "No samples remain after dropping rows with missing design values.",
            "Fill or remove missing design metadata before running DE.",
        ))
        return _result("blocking", issues)

    condition_levels = meta[condition_col].astype(str).nunique(dropna=True)
    if condition_levels < 2:
        issues.append(DesignIssue(
            "blocking",
            "condition_has_one_level",
            f"Condition column '{condition_col}' has fewer than two levels.",
            "Use a condition column with at least two biological groups.",
        ))

    condition_counts = meta[condition_col].astype(str).value_counts().to_dict()
    low_conditions = {
        str(k): int(v)
        for k, v in condition_counts.items()
        if int(v) < min_replicates_per_condition
    }
    if low_conditions:
        severity = "blocking" if min(low_conditions.values()) < 2 else "warning"
        issues.append(DesignIssue(
            severity,
            "insufficient_condition_replicates",
            (
                f"Condition replicate counts below {min_replicates_per_condition}: "
                f"{low_conditions}."
            ),
            "Add biological replicates or lower the replicate floor explicitly with a low-power caveat.",
        ))

    categorical, continuous = _classify_terms(meta, condition_col, covariates, issues)
    categorical_terms = [condition_col] + categorical

    for cov in categorical:
        _check_condition_confounding(meta, condition_col, cov, issues)

    if categorical:
        group_cols = [condition_col, *categorical]
        cell_sizes = meta.groupby(group_cols, dropna=False).size()
        singletons = {
            tuple(map(str, idx if isinstance(idx, tuple) else (idx,))): int(n)
            for idx, n in cell_sizes.items()
            if int(n) < min_replicates_per_design_cell
        }
        if singletons:
            issues.append(DesignIssue(
                "warning",
                "n1_design_cells",
                (
                    f"Some condition x covariate design cells have fewer than "
                    f"{min_replicates_per_design_cell} samples: {singletons}."
                ),
                "Treat affected contrasts as low confidence or simplify the design.",
            ))

    design = _build_design_matrix(meta, categorical_terms, continuous)
    rank = int(np.linalg.matrix_rank(design.to_numpy(dtype=float)))
    n_columns = int(design.shape[1])
    if rank < n_columns:
        issues.append(DesignIssue(
            "blocking",
            "rank_deficient_design",
            f"Design matrix is rank deficient (rank {rank} < {n_columns} columns).",
            "Remove or merge collinear factors before running DESeq2.",
        ))

    status = (
        "blocking" if any(i.severity == "blocking" for i in issues)
        else "warnings" if issues else "clean"
    )
    return {
        "status": status,
        "issues": [i.as_dict() for i in issues],
        "condition_col": condition_col,
        "categorical_factors": categorical_terms,
        "continuous_factors": continuous,
        "design_factors": [condition_col, *covariates],
        "rank": rank,
        "n_columns": n_columns,
        "n_samples": int(len(meta)),
        "residual_degrees_of_freedom": max(0, int(len(meta)) - rank),
        "condition_counts": {str(k): int(v) for k, v in condition_counts.items()},
    }


def _result(status: str, issues: list[DesignIssue]) -> dict[str, Any]:
    return {
        "status": status,
        "issues": [i.as_dict() for i in issues],
        "categorical_factors": [],
        "continuous_factors": [],
        "design_factors": [],
        "rank": 0,
        "n_columns": 0,
        "n_samples": 0,
        "residual_degrees_of_freedom": 0,
        "condition_counts": {},
    }


def _classify_terms(meta, condition_col: str, covariates: list[str],
                    issues: list[DesignIssue]) -> tuple[list[str], list[str]]:
    import pandas as pd

    categorical: list[str] = []
    continuous: list[str] = []

    if pd.api.types.is_numeric_dtype(meta[condition_col]):
        issues.append(DesignIssue(
            "blocking",
            "numeric_condition_factor",
            f"Condition column '{condition_col}' is numeric.",
            "Encode biological groups as categorical labels before DESeq2.",
        ))

    for cov in covariates:
        series = meta[cov]
        n_unique = int(series.nunique(dropna=True))
        if pd.api.types.is_numeric_dtype(series):
            if n_unique > 2:
                continuous.append(cov)
            else:
                categorical.append(cov)
                issues.append(DesignIssue(
                    "warning",
                    "binary_numeric_covariate",
                    f"Covariate '{cov}' is numeric with {n_unique} levels and will be treated as categorical.",
                    "Use explicit labels for categorical batches, or provide a truly continuous covariate.",
                ))
            continue

        numeric = pd.to_numeric(series, errors="coerce")
        numeric_fraction = float(numeric.notna().mean()) if len(series) else 0.0
        if numeric_fraction >= 0.95 and n_unique > 5:
            continuous.append(cov)
            meta[cov] = numeric.astype(float)
            issues.append(DesignIssue(
                "warning",
                "numeric_covariate_stored_as_text",
                f"Covariate '{cov}' is stored as text but appears numeric.",
                "Store continuous covariates as numeric metadata to avoid factor expansion.",
            ))
        else:
            categorical.append(cov)

    return categorical, continuous


def factors_confounded_with_condition(
    sample_sheet: list[dict],
    condition_key: str,
    candidate_keys: list[str],
) -> list[str]:
    """Return candidate factors perfectly aliased with the condition.

    Dependency-light sibling of :func:`_check_condition_confounding` that works on
    a plain list-of-dict sample sheet (no pandas), so the confounding check can run
    at design-assembly time — before any DE metadata frame exists.

    A candidate is *confounded* with the condition when the level mapping is
    bijective: every candidate level appears under exactly one condition level and
    every condition level appears under exactly one candidate level.  Such a term
    cannot be separated from biology and must block, not merely be "ignored".

    Rows missing either key are skipped.  A candidate with fewer than two observed
    levels, or a condition with fewer than two observed levels, is never confounded.
    """
    confounded: list[str] = []
    for candidate in candidate_keys:
        if not candidate or candidate == condition_key:
            continue
        cond_to_cov: dict[str, set] = {}
        cov_to_cond: dict[str, set] = {}
        for row in sample_sheet:
            cond = row.get(condition_key)
            cov = row.get(candidate)
            if cond in (None, "") or cov in (None, ""):
                continue
            cond = str(cond)
            cov = str(cov)
            cond_to_cov.setdefault(cond, set()).add(cov)
            cov_to_cond.setdefault(cov, set()).add(cond)
        if len(cond_to_cov) < 2 or len(cov_to_cond) < 2:
            continue
        cov_to_one_condition = all(len(v) == 1 for v in cov_to_cond.values())
        condition_to_one_cov = all(len(v) == 1 for v in cond_to_cov.values())
        if cov_to_one_condition and condition_to_one_cov:
            confounded.append(candidate)
    return confounded


def classify_pairing(
    sample_sheet: list[dict],
    condition_key: str,
    replicate_key: str,
    *,
    test: str | None = None,
    ref: str | None = None,
) -> dict:
    """Classify the pairing structure of a two-condition contrast.

    Dependency-light (no pandas): works on a plain list-of-dict sample sheet so
    the classification runs before any DE metadata frame exists.  Only the two
    contrast levels are considered when ``test``/``ref`` are given; a replicate
    that pairs across other conditions does not count as paired here.

    A replicate is *paired* when it appears under both contrast conditions.

      * ``complete``    — every replicate present in the contrast is paired.
      * ``partial``     — some replicates are paired, some are nested in one level.
      * ``independent`` — no replicate spans both conditions.
      * ``insufficient``— fewer than two usable samples/conditions.

    Returns paired/unpaired replicate ids plus the sample ids that belong to each,
    so a caller can model the within-subject block on the complete-paired subset
    and disclose the excluded unpaired samples (B3 option B).
    """
    levels = {str(test), str(ref)} if (test is not None and ref is not None) else None

    rep_to_conditions: dict[str, set] = {}
    rep_to_samples: dict[str, list] = {}
    for row in sample_sheet:
        cond = row.get(condition_key)
        rep = row.get(replicate_key)
        sample = row.get("sample", row.get("id"))
        if cond in (None, "") or rep in (None, ""):
            continue
        cond = str(cond)
        rep = str(rep)
        if levels is not None and cond not in levels:
            continue
        rep_to_conditions.setdefault(rep, set()).add(cond)
        rep_to_samples.setdefault(rep, []).append(sample)

    contrast_levels = levels or {
        c for conds in rep_to_conditions.values() for c in conds
    }

    paired_replicates: list[str] = []
    unpaired_replicates: list[str] = []
    for rep in sorted(rep_to_conditions):
        seen = rep_to_conditions[rep] & contrast_levels
        if len(seen) >= 2:
            paired_replicates.append(rep)
        else:
            unpaired_replicates.append(rep)

    paired_samples = [
        s for rep in paired_replicates for s in rep_to_samples[rep]
    ]
    unpaired_samples = [
        s for rep in unpaired_replicates for s in rep_to_samples[rep]
    ]

    n_reps = len(rep_to_conditions)
    n_samples = sum(len(v) for v in rep_to_samples.values())
    if n_reps < 2 or n_samples < 2 or len(contrast_levels) < 2:
        status = "insufficient"
    elif not paired_replicates:
        status = "independent"
    elif not unpaired_replicates:
        status = "complete"
    else:
        status = "partial"

    return {
        "status": status,
        "paired_replicates": paired_replicates,
        "unpaired_replicates": unpaired_replicates,
        "paired_samples": paired_samples,
        "unpaired_samples": unpaired_samples,
    }


def _check_condition_confounding(meta, condition_col: str, covariate: str,
                                 issues: list[DesignIssue]) -> None:
    table = meta.groupby([covariate, condition_col], dropna=False).size().unstack(fill_value=0)
    if table.shape[0] < 2 or table.shape[1] < 2:
        return

    cov_to_one_condition = bool((table > 0).sum(axis=1).eq(1).all())
    condition_to_one_cov = bool((table > 0).sum(axis=0).eq(1).all())
    if cov_to_one_condition and condition_to_one_cov:
        issues.append(DesignIssue(
            "blocking",
            "condition_covariate_confounding",
            (
                f"Condition '{condition_col}' is completely confounded with "
                f"covariate '{covariate}'."
            ),
            "Do not fit both terms together; redesign the comparison or drop the confounded covariate.",
        ))


def _build_design_matrix(meta, categorical_terms: list[str], continuous_terms: list[str]):
    import pandas as pd

    parts = [pd.DataFrame({"Intercept": [1.0] * len(meta)}, index=meta.index)]
    for term in categorical_terms:
        dummies = pd.get_dummies(meta[term].astype(str), prefix=term, drop_first=True)
        if not dummies.empty:
            parts.append(dummies.astype(float))
    for term in continuous_terms:
        vals = pd.to_numeric(meta[term], errors="coerce").astype(float)
        centered = vals - vals.mean()
        parts.append(pd.DataFrame({term: centered}, index=meta.index))
    return pd.concat(parts, axis=1)
