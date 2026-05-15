"""
ARIA Differential Abundance Script
-----------------------------------
Tests whether cell-type proportions differ between conditions. Run BEFORE
pseudobulk DE so that the downstream narrative can:

  1. Report shifts in cell-type composition as a primary observation in
     their own right (e.g. "microglia increased 1.8x in aged hippocampus").
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
  Primary: Poisson GLM with offset(log total_cells_per_replicate),
           Wald p-value on the condition coefficient, BH correction
           across all cell types tested in the comparison. statsmodels.
  Fallback: Fisher's exact test per cell type on a 2x2 table of
           (cells_in_type, cells_in_other_types) x (test, ref) with
           BH correction. Only used when statsmodels is missing or the
           Poisson GLM fails (e.g. perfect separation).

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
      "method": "poisson_offset_glm" | "fisher_exact_fallback",
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
             "padj": float,         // BH within this comparison
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


def rna_diff_abundance(params: dict) -> dict:
    import math
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from pathlib import Path

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

    adata = sc.read_h5ad(data_path)

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

    # (rep, cell_type) -> count
    pivot = (obs.groupby([replicate_col, groupby])
                .size()
                .unstack(fill_value=0))
    # Replicate -> total cells (across all cell types observed in adata)
    total_per_rep = pivot.sum(axis=1)

    # Replicate -> condition (constant per rep by construction)
    rep_to_condition = (obs.drop_duplicates(replicate_col)
                          .set_index(replicate_col)[condition_col]
                          .to_dict())

    # Replicate -> covariate values (constant per rep)
    rep_to_covs = {}
    if covariates:
        rep_to_covs = (obs.drop_duplicates(replicate_col)
                         .set_index(replicate_col)[covariates]
                         .to_dict(orient="index"))

    method = "poisson_offset_glm"
    try:
        import statsmodels.api as sm
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        sm = None
        multipletests = None
        method = "fisher_exact_fallback"
        warnings.append(
            "statsmodels not available; falling back to per-cell-type "
            "Fisher's exact + manual BH. Install statsmodels for the "
            "primary Poisson-GLM-with-offset method."
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
        pvals = []

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

            if sm is not None:
                try:
                    # Design matrix: intercept + condition dummy (+ covariates)
                    X = pd.DataFrame({
                        "is_test": (cond_vec.values == test_lvl).astype(float),
                    }, index=rep_subset)
                    if covariates:
                        for cov in covariates:
                            X[cov] = [rep_to_covs[r][cov] for r in rep_subset]
                        # One-hot any categorical covariates so statsmodels
                        # does not silently coerce strings to NaN.
                        X = pd.get_dummies(X, drop_first=True).astype(float)
                    X_design = sm.add_constant(X, has_constant="add")
                    model = sm.GLM(
                        endog=y,
                        exog=X_design,
                        exposure=np.clip(total, 1, None),
                        family=sm.families.Poisson(),
                    )
                    fit = model.fit(disp=False, maxiter=200)
                    pval = float(fit.pvalues.get("is_test", float("nan")))
                except Exception as exc:
                    # Perfect separation / singular X / convergence failure:
                    # surface via warning and fall through to Fisher for this
                    # cell type.
                    warnings.append(
                        f"[{comp_key}/{ct}] Poisson GLM failed ({exc!s:.120}); "
                        f"using Fisher exact for this cell type."
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

            rows.append({
                "name":             str(ct),
                "n_test":           n_t,
                "n_ref":            n_r,
                "prop_test":        prop_t,
                "prop_ref":         prop_r,
                "log2_fold_change": log2fc,
                "pval":             pval,
            })
            pvals.append(pval)

        # BH correction across cell types within this comparison
        if rows:
            if multipletests is not None:
                _, padj_arr, _, _ = multipletests(pvals, method="fdr_bh")
            else:
                padj_arr = _bh_correct(pvals)
            for r, padj in zip(rows, padj_arr):
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
        per_comparison[comp_key] = {
            "status":         "success",
            "per_cell_type":  sorted(rows, key=lambda r: r["padj"]),
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

    n_reps_per_group = (obs.drop_duplicates(replicate_col)
                          .groupby(condition_col)[replicate_col]
                          .nunique()
                          .to_dict())

    return {
        "status":                 "success",
        "method":                 method,
        "groupby":                groupby,
        "condition_col":          condition_col,
        "replicate_col":          replicate_col,
        "covariates":             covariates,
        "n_replicates_per_group": n_reps_per_group,
        "significance_alpha":     alpha,
        "any_significant":        any_significant,
        "per_comparison":         per_comparison,
        "output_path":            out_path,
        "warnings":               warnings,
    }


def _bh_correct(pvals):
    """Benjamini-Hochberg correction without statsmodels."""
    import numpy as np
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    # Monotone non-increasing from the largest down
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = adj
    return out


if __name__ == "__main__":
    run_script(rna_diff_abundance)
