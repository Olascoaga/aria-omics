"""
ARIA Differential Abundance Script
-----------------------------------
Tests whether cell-type proportions differ between conditions. Run BEFORE
pseudobulk DE so that the downstream narrative can:

  1. Report shifts in cell-type composition as a primary observation in
     their own right (e.g. "cell type A increased 1.8x in condition B").
  2. Pass composition_covariate=True to rna_pseudobulk_de.py whenever a
     significant abundance shift exists, so that per-cell-type expression
     contrasts are not confounded by changes in the cell-type mixture.

Why this matters. Pseudobulk DE per cell type aggregates counts within
(donor, cell_type) blocks. If the proportion of a cell type changes
between conditions, the aggregated counts for that type are dominated by
a different sub-population in each group, producing spurious DE that
reflects mixture shift rather than within-type regulation. Standard
practice (Crowell et al. 2020, Cao et al. 2021) is to (a) report
abundance separately and (b) include a composition covariate in the DE
design when abundance is shifting.

Method.
  Primary: donor-level compositional model. ARIA converts each pseudosample's
           cell-type counts to proportions, applies a centered log-ratio (CLR)
           transform with a small pseudocount, and fits per-cell-type OLS with
           HC3 robust standard errors:

               CLR(cell type proportion) ~ condition + covariates (+ donor FE)

           Donor fixed effects are included automatically for paired designs.
           This respects the compositional constraint (cell-type proportions
           sum to one) and tests shifts at the biological replicate level rather
           than treating cells as independent observations (P1-3, audit
           2026-06-02). BH correction is across donor-level CLR tests in the
           comparison. statsmodels.
  Fallback: Fisher's exact test per cell type on a 2x2 table of
           (cells_in_type, cells_in_other_types) x (test, ref). Fisher results
           are cell-level diagnostics and are not included in donor-level FDR.
           Only used when statsmodels is missing or the compositional model
           fails for a cell type.

Input params:
    data_path:        str — path to annotated .h5ad
    groupby:          str — obs column with cell-type/cluster labels
    condition_col:    str — obs column with experimental condition
    replicate_col:    str — obs column with biological replicate ID
    comparisons:      [[test, ref], ...] — list of pairs (matches
                       rna_pseudobulk_de.py's contract)
    covariates:       [str] (optional) — extra obs columns to include
    output_dir:       str (optional)
    significance_alpha: float (default 0.10) — threshold used downstream
                       to decide whether to gate composition_covariate
                       in pseudobulk DE. Reported here for transparency.

Output:
    {
      "status": "success",
      "method": "donor_clr_ols_hc3" | "fisher_exact_fallback",
      "groupby": str,
      "condition_col": str,
      "replicate_col": str,
      "n_replicates_per_group": {cond: int, ...},
      "significance_alpha": float,
      "any_significant": bool,
      "per_comparison": {
        "test_vs_ref": {
          "per_cell_type": [
            {"name": str,
             "n_test": int,         // total cells in this type across test reps
             "n_ref":  int,
             "prop_test": float,    // mean per-replicate proportion in test
             "prop_ref":  float,
             "log2_fold_change": float,  // log2(prop_test / prop_ref)
             "pval": float,
             "clr_log_ratio": float | null,  // condition coefficient on the
                                             // CLR scale; null for fallback
             "dispersion": float,   // retained for backward compatibility;
                                    // NaN for compositional/fallback paths
             "padj": float | null,  // BH within donor-level CLR tests only;
                                    // null for cell-level Fisher diagnostics
             "fdr_included": bool,
             "pval_role": "donor_level_primary"|"cell_level_diagnostic",
             "direction": "up"|"down"|"none",
             "significant": bool},
            ...
          ],
          "n_significant": int,
          "n_replicates":  {"test": int, "ref": int}
        }
      },
      "output_path": str | None,
      "warnings": [str]
    }
"""

from __future__ import annotations

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from aria.scripts._base import run_script
from aria.utils.stats import bh_correct as _bh_correct


def _clr_transform_counts(counts):
    """Return centered log-ratio values for sample x cell-type counts."""
    import numpy as np

    arr = counts.astype(float).to_numpy()
    k = max(arr.shape[1], 1)
    smoothed = arr + 0.5
    props = smoothed / np.clip(smoothed.sum(axis=1, keepdims=True), 1e-12, None)
    log_props = np.log(props)
    clr = log_props - log_props.mean(axis=1, keepdims=True)
    return counts.__class__(clr, index=counts.index, columns=counts.columns)


def _build_clr_design(rep_subset, cond_vec, test_lvl, covariates, rep_to_covs,
                      rep_to_base, paired_design):
    """Build a donor-level design matrix for compositional abundance tests."""
    import pandas as pd
    import statsmodels.api as sm

    X = pd.DataFrame({
        "is_test": (cond_vec.values == test_lvl).astype(float),
    }, index=rep_subset)
    if covariates:
        for cov in covariates:
            X[cov] = [rep_to_covs[r][cov] for r in rep_subset]
    if paired_design:
        donors = pd.Series([rep_to_base.get(r, r) for r in rep_subset],
                           index=rep_subset, name="_aria_donor")
        donor_has_both = (
            pd.DataFrame({"donor": donors, "condition": cond_vec})
              .drop_duplicates()
              .groupby("donor")["condition"]
              .nunique()
        )
        if (donor_has_both > 1).any():
            donor_dummies = pd.get_dummies(
                donors, prefix="_aria_donor", drop_first=True
            )
            X = pd.concat([X, donor_dummies], axis=1)
    if len(X.columns) > 1:
        X = pd.get_dummies(X, drop_first=True)
    X = X.astype(float)
    return sm.add_constant(X, has_constant="add")


def rna_diff_abundance(params: dict) -> dict:
    import math
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from aria.utils.safe_h5ad import read_h5ad

    data_path     = params["data_path"]
    groupby       = params["groupby"]
    condition_col = params["condition_col"]
    replicate_col = params["replicate_col"]
    comparisons   = params["comparisons"]
    covariates    = list(params.get("covariates") or [])
    output_dir    = params.get("output_dir")
    alpha         = float(params.get("significance_alpha", 0.10))

    warnings: list[str] = []

    if not Path(data_path).exists():
        return {"status":     "error",
                "error_type": "FileNotFound",
                "details":    f"data_path does not exist: {data_path}"}

    adata = read_h5ad(data_path)

    for col in [groupby, condition_col, replicate_col, *covariates]:
        if col not in adata.obs.columns:
            return {"status":     "error",
                    "error_type": "MissingObsColumn",
                    "details":    (f"Required obs column '{col}' not found. "
                                   f"Available: {list(adata.obs.columns)}")}

    obs = adata.obs[[groupby, condition_col, replicate_col, *covariates]].copy()
    obs[groupby]       = obs[groupby].astype(str)
    obs[condition_col] = obs[condition_col].astype(str)
    obs[replicate_col] = obs[replicate_col].astype(str)

    rep_condition_counts = obs.groupby(replicate_col)[condition_col].nunique()
    paired_design = bool((rep_condition_counts > 1).any())
    if paired_design:
        obs["_aria_pseudosample_key"] = (
            obs[replicate_col].astype(str)
            + "__"
            + obs[condition_col].astype(str)
        )
    else:
        obs["_aria_pseudosample_key"] = obs[replicate_col].astype(str)

    # (pseudosample, cell_type) -> count
    pivot = (obs.groupby(["_aria_pseudosample_key", groupby])
                .size()
                .unstack(fill_value=0))
    # Pseudosample -> total cells (across all cell types observed in adata)
    total_per_rep = pivot.sum(axis=1)

    # Pseudosample -> condition (constant because paired keys include it)
    rep_to_condition = (obs.drop_duplicates("_aria_pseudosample_key")
                          .set_index("_aria_pseudosample_key")[condition_col]
                          .to_dict())
    rep_to_base = (obs.drop_duplicates("_aria_pseudosample_key")
                     .set_index("_aria_pseudosample_key")[replicate_col]
                     .to_dict())

    # Pseudosample -> covariate values
    rep_to_covs = {}
    if covariates:
        rep_to_covs = (obs.drop_duplicates("_aria_pseudosample_key")
                         .set_index("_aria_pseudosample_key")[covariates]
                         .to_dict(orient="index"))

    method = "donor_clr_ols_hc3"
    try:
        import statsmodels.api as sm
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        sm = None
        multipletests = None
        method = "fisher_exact_fallback"
        warnings.append(
            "statsmodels not available; falling back to per-cell-type "
            "Fisher's exact + manual BH. Install statsmodels for the primary "
            "donor-level CLR OLS model with HC3 robust standard errors."
        )

    per_comparison: dict = {}
    any_significant = False

    for comp in comparisons:
        test_lvl, ref_lvl = comp[0], comp[1]
        comp_key = f"{test_lvl}_vs_{ref_lvl}"

        # Restrict to replicates in this contrast
        rep_subset = [r for r, c in rep_to_condition.items()
                      if c in (test_lvl, ref_lvl)]
        if not rep_subset:
            per_comparison[comp_key] = {
                "status": "skipped",
                "reason": f"no replicates match {test_lvl} or {ref_lvl}",
            }
            continue

        sub_pivot = pivot.loc[rep_subset].copy()
        sub_total = total_per_rep.loc[rep_subset]
        clr_pivot = _clr_transform_counts(sub_pivot)

        cond_vec = pd.Series([rep_to_condition[r] for r in rep_subset],
                             index=rep_subset)
        n_test_reps = int((cond_vec == test_lvl).sum())
        n_ref_reps  = int((cond_vec == ref_lvl).sum())

        if n_test_reps < 2 or n_ref_reps < 2:
            per_comparison[comp_key] = {
                "status": "skipped",
                "reason": (f"need at least 2 replicates per condition "
                           f"({test_lvl}={n_test_reps}, "
                           f"{ref_lvl}={n_ref_reps})"),
            }
            continue

        cell_types = list(sub_pivot.columns)
        rows = []

        for ct in cell_types:
            y = sub_pivot[ct].astype(float).values
            total = sub_total.astype(float).values

            # Skip cell types that are absent everywhere in this contrast
            if y.sum() == 0:
                continue

            # Per-replicate proportions (smoothed by +1 to avoid log(0))
            props_test = (y[cond_vec.values == test_lvl] /
                          np.clip(total[cond_vec.values == test_lvl], 1, None))
            props_ref  = (y[cond_vec.values == ref_lvl] /
                          np.clip(total[cond_vec.values == ref_lvl], 1, None))
            prop_t = float(props_test.mean()) if len(props_test) else 0.0
            prop_r = float(props_ref.mean())  if len(props_ref)  else 0.0
            n_t    = int(y[cond_vec.values == test_lvl].sum())
            n_r    = int(y[cond_vec.values == ref_lvl].sum())
            log2fc = float(math.log2((prop_t + 1e-9) / (prop_r + 1e-9)))

            pval = float("nan")
            dispersion = float("nan")
            clr_log_ratio = None

            if sm is not None:
                try:
                    # P1-3: compositional donor-level model. The response is
                    # the cell type's centered log-ratio abundance in each
                    # pseudosample, so the test respects the sum-to-one
                    # constraint across cell types. HC3 keeps small-n standard
                    # errors conservative; paired designs include donor fixed
                    # effects when donors appear in both conditions.
                    X_design = _build_clr_design(
                        rep_subset, cond_vec, test_lvl, covariates,
                        rep_to_covs, rep_to_base, paired_design,
                    )
                    model = sm.OLS(
                        endog=clr_pivot.loc[rep_subset, ct].astype(float).values,
                        exog=X_design,
                        missing="drop",
                    )
                    fit = model.fit(cov_type="HC3", use_t=True)
                    clr_log_ratio = float(
                        fit.params.get("is_test", float("nan"))
                    )
                    pval = float(fit.pvalues.get("is_test", float("nan")))
                except Exception as exc:
                    # Singular X / zero residual df / convergence failure:
                    # surface via warning and fall through to Fisher for this
                    # cell type as a diagnostic rather than reporting a
                    # fabricated donor-level p-value.
                    warnings.append(
                        f"[{comp_key}/{ct}] donor-level CLR model failed "
                        f"({exc!s:.120}); using Fisher exact as a cell-level "
                        "diagnostic for this cell type."
                    )
                    pval = float("nan")

            if math.isnan(pval):
                # Fisher exact fallback on 2x2
                try:
                    from scipy.stats import fisher_exact
                    other_test = int(total[cond_vec.values == test_lvl].sum() - n_t)
                    other_ref  = int(total[cond_vec.values == ref_lvl].sum()  - n_r)
                    table = [[n_t, other_test], [n_r, other_ref]]
                    _, pval = fisher_exact(table)
                    pval = float(pval)
                except Exception as exc:
                    warnings.append(
                        f"[{comp_key}/{ct}] Fisher exact also failed: "
                        f"{exc!s:.120}"
                    )
                    pval = 1.0

            row_model = (
                "fisher_exact_fallback"
                if clr_log_ratio is None else "donor_clr_ols_hc3"
            )
            fdr_included = bool(
                row_model == "donor_clr_ols_hc3" and math.isfinite(pval)
            )
            rows.append({
                "name":             str(ct),
                "n_test":           n_t,
                "n_ref":            n_r,
                "prop_test":        prop_t,
                "prop_ref":         prop_r,
                "log2_fold_change": log2fc,
                "pval":             pval,
                "dispersion":       dispersion,
                "clr_log_ratio":    clr_log_ratio,
                "model":            row_model,
                "fdr_included":     fdr_included,
                "pval_role":        (
                    "donor_level_primary"
                    if fdr_included else "cell_level_diagnostic"
                ),
            })

        # BH correction across donor-level CLR tests within this comparison.
        # Cell-level Fisher fallbacks are diagnostics only; mixing them into
        # the same FDR family would treat cells as independent for one row while
        # using biological replicates for another.
        if rows:
            primary_rows = [r for r in rows if r["fdr_included"]]
            pvals = [r["pval"] for r in primary_rows]
            if not primary_rows:
                padj_arr = []
            elif multipletests is not None:
                _, padj_arr, _, _ = multipletests(pvals, method="fdr_bh")
            else:
                padj_arr = _bh_correct(pvals)
            for r in rows:
                r["padj"] = None
                r["significant"] = False
                r["direction"] = "none"
            for r, padj in zip(primary_rows, padj_arr):
                r["padj"] = float(padj)
                r["significant"] = bool(r["padj"] < alpha)
                if r["significant"]:
                    any_significant = True
                if r["log2_fold_change"] > 0 and r["significant"]:
                    r["direction"] = "up"
                elif r["log2_fold_change"] < 0 and r["significant"]:
                    r["direction"] = "down"
                else:
                    r["direction"] = "none"

        n_sig = sum(1 for r in rows if r["significant"])
        n_fisher = sum(1 for r in rows if r.get("model") == "fisher_exact_fallback")
        n_clr = sum(1 for r in rows if r.get("fdr_included"))
        degraded = n_fisher > 0
        degradation_reason = None
        caveats = []
        # T9: CLR components are sum-to-zero (compositionally dependent). We fit an
        # independent OLS per cell type and BH-correct across them, which does NOT
        # model that dependence: a genuine increase in one cell type mechanically
        # forces apparent decreases in the others. Disclose the assumption so the
        # per-type directions are read jointly, not as independent shifts. The
        # runtime method is unchanged (a joint compositional model would replace
        # it; out of scope here).
        compositional_dependence_caveat = None
        if n_clr > 0:
            compositional_dependence_caveat = (
                "Compositional dependence: per-cell-type CLR abundance tests are "
                "fit independently and BH-corrected together, but CLR values are "
                "sum-to-zero (compositionally dependent), so an increase in one "
                "cell type mechanically forces apparent decreases in the others. "
                "Interpret the per-type directions jointly as one compositional "
                "shift, not as independent changes; a joint compositional model "
                "(e.g. scCODA or propeller) would model this dependence directly."
            )
            caveats.append(compositional_dependence_caveat)
        if degraded:
            degradation_reason = (
                f"{n_fisher} cell type(s) used Fisher exact as a cell-level "
                "diagnostic and were excluded from donor-level FDR."
            )
            caveats.append(degradation_reason)

        def _row_sort_key(row):
            padj = row.get("padj")
            if padj is None:
                return (float("inf"), str(row.get("name", "")))
            try:
                padj_float = float(padj)
            except (TypeError, ValueError):
                return (float("inf"), str(row.get("name", "")))
            if math.isnan(padj_float):
                return (float("inf"), str(row.get("name", "")))
            return (padj_float, str(row.get("name", "")))

        per_comparison[comp_key] = {
            "status":         "success",
            "model":          method,
            "transform":      "centered_log_ratio",
            "covariance":     "HC3",
            "fdr_family":     "donor_level_clr_only",
            "n_fdr_tests":    sum(1 for r in rows if r.get("fdr_included")),
            "n_fisher_diagnostic": n_fisher,
            "degraded":       degraded,
            "confidence":     "degraded" if degraded else "standard",
            "degradation_reason": degradation_reason,
            "compositional_dependence_caveat": compositional_dependence_caveat,
            "caveats":        caveats,
            "per_cell_type":  sorted(rows, key=_row_sort_key),
            "n_significant":  n_sig,
            "n_replicates":   {"test": n_test_reps, "ref": n_ref_reps},
        }

    # Persist a flat TSV for the report sidecar
    out_path = None
    if output_dir:
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = str(out_dir / "differential_abundance.tsv")
        flat = []
        for comp_key, info in per_comparison.items():
            if info.get("status") != "success":
                continue
            for r in info["per_cell_type"]:
                flat.append({"comparison": comp_key, **r})
        pd.DataFrame(flat).to_csv(out_path, sep="\t", index=False)

    n_reps_per_group = (obs.drop_duplicates("_aria_pseudosample_key")
                          .groupby(condition_col)["_aria_pseudosample_key"]
                          .nunique()
                          .to_dict())

    return {
        "status":                 "success",
        "method":                 method,
        "model": {
            "primary": method,
            "transform": "centered_log_ratio",
            "covariance": "HC3",
            "paired_donor_fixed_effects": paired_design,
            "unit": "biological replicate / donor pseudosample",
        },
        "groupby":                groupby,
        "condition_col":          condition_col,
        "replicate_col":          replicate_col,
        "covariates":             covariates,
        "paired_design":          paired_design,
        "n_replicates_per_group": n_reps_per_group,
        "significance_alpha":     alpha,
        "any_significant":        any_significant,
        "per_comparison":         per_comparison,
        "output_path":            out_path,
        "warnings":               warnings,
    }


if __name__ == "__main__":
    run_script(rna_diff_abundance)
