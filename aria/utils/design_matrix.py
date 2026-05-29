"""Design-matrix sanity checks before DE model fitting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


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
