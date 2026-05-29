"""
ARIA Pseudobulk Differential Expression
----------------------------------------
Between-condition DE for scRNA-seq, the way the methods sections of real
papers do it: aggregate raw counts per (biological replicate × cell type),
then hand the resulting count matrix to DESeq2.

Why not per-cell Wilcoxon between conditions? Single cells from the same
donor are not independent observations; treating them as replicates
inflates statistical significance by 1-2 orders of magnitude. Pseudobulk
respects the biological replication unit (the donor/sample).

Input params:
    data_path:        str   — preprocessed .h5ad. Must have `.raw` populated
                              with integer counts, OR `X` itself must be raw
                              counts (use_raw=False).
    groupby:          str   — obs column to stratify by (e.g. "subclass",
                              "leiden", "cell_type"). One DESeq2 model is
                              fit per level of this column.
    condition_col:    str   — obs column with the experimental factor
                              (e.g. "age_group", "treatment").
    replicate_col:    str   — obs column identifying the biological
                              replication unit (e.g. "orig.ident", "donor").
                              Aggregation key.
    comparisons:      [[test_level, ref_level], ...]
                              List of pairwise contrasts within
                              condition_col. Each becomes a DESeq2 contrast.
    covariates:       [str] (optional) — extra obs columns to include in
                              the design formula (e.g. ["Gender"]).
    min_cells_per_pseudosample: int (default 10) — drop (replicate × group)
                              combinations smaller than this; their counts
                              are too noisy to be useful.
    min_replicates_per_condition: int (default 3) — DESeq2 floor; groups
                              that don't have ≥ min replicates in BOTH
                              levels of a comparison are reported as skipped.
    use_raw:          bool (default True) — use adata.raw.X for counts.
    padj_max:         float (default 0.05)
    lfc_min:          float (default 0.5)
    output_dir:       str   — CSV destination (default: dirname of input)

Output:
    {
      "status":    "success" | "error",
      "n_groups":  int,
      "per_group": {
          group_id: {
              "n_pseudosamples":   int,
              "per_comparison": {
                  "<test>_vs_<ref>": {
                      "status":       "success" | "skipped",
                      "n_significant": int,
                      "n_up":         int,
                      "n_down":       int,
                      "top_genes":    [ {gene, log2fc, padj}, ... ],
                      "all_sig":      [ {gene, log2fc, padj}, ... ],
                  },
                  ...
              }
          },
          ...
      },
      "output_csv": str,
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


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
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = adj
    return out


def _global_bh(pvals):
    """Return BH-adjusted p-values for one global family of tests."""
    try:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        return padj
    except Exception:
        return _bh_correct(pvals)


def rna_pseudobulk_de(params: dict) -> dict:
    import warnings as _w
    _w.filterwarnings("ignore")

    import numpy as np
    import pandas as pd
    from pathlib import Path
    from scipy import sparse
    import inspect
    from aria.utils.design_matrix import validate_design_matrix
    from aria.utils.power_estimation import pseudobulk_power_estimate
    from aria.utils.safe_h5ad import read_h5ad

    data_path                     = params["data_path"]
    groupby                       = params["groupby"]
    condition_col                 = params["condition_col"]
    replicate_col                 = params["replicate_col"]
    comparisons                   = params["comparisons"]
    covariates                    = params.get("covariates") or []
    min_cells_per_pseudosample    = int(params.get("min_cells_per_pseudosample", 10))
    min_replicates_per_condition  = int(params.get("min_replicates_per_condition", 3))
    use_raw                       = bool(params.get("use_raw", True))
    # T1.1: when True, add log(n_cells_in_group / n_cells_in_replicate) as a
    # continuous covariate in the DESeq2 design. The scrna_agent flips this
    # on whenever rna_diff_abundance flagged any cell type as significantly
    # shifting at the user-configured alpha (default 0.10). Per-block result
    # carries corrected_for_composition: True so the narrative can describe
    # it honestly.
    composition_covariate         = bool(params.get("composition_covariate", False))
    COMPOSITION_COL               = "_aria_composition_log_ratio"
    # Many published h5ads (especially Seurat exports) store log-normalized
    # values in raw.X rather than counts. When True we reverse NormalizeData
    # (counts = round(expm1(x) * lib_size / scale_factor)) using the
    # per-cell library size preserved in obs. Set to False to refuse such
    # inputs explicitly.
    allow_lognorm_recovery        = bool(params.get("allow_lognorm_recovery", True))
    lib_size_col                  = params.get("lib_size_col", "nCount_RNA")
    norm_scale_factor             = float(params.get("norm_scale_factor", 10000.0))
    padj_max                      = float(params.get("padj_max", 0.05))
    lfc_min                       = float(params.get("lfc_min", 0.5))
    top_n                         = int(params.get("top_n", 50))
    # F-SCI-FDR (audit 2026-05-28): which BH family defines "significant".
    # 'per_cluster' = per-block BH (field standard, e.g. muscat/Crowell 2020;
    # decouples unrelated cell types). 'global' = one BH family pooled across
    # gene x cell-type x contrast (whole-experiment FDR control, conservative).
    # Both adjusted p-values are always computed and reported; this selects the
    # primary significance call only. Default per_cluster.
    fdr_strategy                  = str(params.get("fdr_strategy", "per_cluster")).lower()
    if fdr_strategy not in ("per_cluster", "global"):
        fdr_strategy = "per_cluster"
    output_dir                    = params.get("output_dir")
    auto_paired_donor_covariate   = bool(
        params.get("auto_paired_donor_covariate", True)
    )

    if not Path(data_path).exists():
        return {"status": "error",
                "error_type": "FileNotFound",
                "details": f"data_path does not exist: {data_path}"}

    adata = read_h5ad(data_path)

    # ── Resolve count source ──────────────────────────────────────────────
    def _looks_integerlike(mat) -> tuple[bool, float]:
        if hasattr(mat, "toarray"):
            sample = mat[:200].toarray()
        else:
            sample = np.asarray(mat[:200])
        max_val = float(np.nanmax(sample)) if sample.size else 0.0
        if np.issubdtype(sample.dtype, np.integer):
            return True, max_val
        # log-normalized data sits in [0, ~10]; raw counts are non-negative
        # integers spanning [0, thousands]. The combination of "large max"
        # and "values close to integers" is the discriminator.
        rounded_ok = np.allclose(sample, np.round(sample))
        return (max_val > 50 and rounded_ok), max_val

    def _validate_lognorm_recovery(mat, lib_sizes) -> bool:
        """Probe: does (expm1(x) * lib/scale) on a small block produce
        integer-like values? Used to gate per-block recovery without
        materializing the full matrix."""
        if hasattr(mat, "toarray"):
            block = mat[:200].toarray()
        else:
            block = np.asarray(mat[:200])
        if block.size == 0 or block.max() > 50 or block.min() < 0:
            return False
        sample_lib = np.asarray(lib_sizes[:200], dtype=float)
        recovered  = np.expm1(block) * sample_lib[:, None] / norm_scale_factor
        frac_int   = float((np.abs(recovered - np.round(recovered)) < 0.05).mean())
        return frac_int >= 0.85

    if use_raw and adata.raw is not None:
        counts     = adata.raw.X
        gene_names = list(adata.raw.var_names)
    else:
        counts     = adata.X
        gene_names = list(adata.var_names)

    try:
        if sparse.issparse(counts):
            expressed_mask = np.asarray((counts > 0).sum(axis=0)).ravel() > 0
        else:
            expressed_mask = (np.asarray(counts) > 0).sum(axis=0) > 0
        background_genes = [
            str(g) for g, keep in zip(gene_names, expressed_mask) if keep
        ]
    except Exception:
        background_genes = list(map(str, gene_names))

    integerlike, max_val = _looks_integerlike(counts)
    needs_recovery       = False
    if not integerlike:
        if not allow_lognorm_recovery:
            return {"status":     "error",
                    "error_type": "NonIntegerCounts",
                    "details":    (f"counts appear non-integer (max={max_val:.2f}). "
                                   f"Set allow_lognorm_recovery=True with a "
                                   f"valid lib_size_col, or supply an h5ad "
                                   f"with raw counts.")}
        # Tolerate both Seurat-derived (nCount_RNA) and scanpy-derived
        # (total_counts / n_counts) library-size column names — try the
        # requested one first, then the common alternatives.
        candidate_libsize_cols = [lib_size_col,
                                   "nCount_RNA", "total_counts", "n_counts"]
        chosen_libsize_col = None
        for c in candidate_libsize_cols:
            if c and c in adata.obs.columns:
                chosen_libsize_col = c
                break
        if chosen_libsize_col is None:
            return {"status":     "error",
                    "error_type": "MissingLibSizeColumn",
                    "details":    (f"counts look log-normalized (max={max_val:.2f}) "
                                   f"but none of {candidate_libsize_cols} are "
                                   f"in obs. Cannot recover counts.")}
        lib_size_col = chosen_libsize_col
        if not _validate_lognorm_recovery(counts, adata.obs[lib_size_col].values):
            return {"status":     "error",
                    "error_type": "LognormRecoveryFailed",
                    "details":    (f"reversing with lib_size_col="
                                   f"'{lib_size_col}' and scale_factor="
                                   f"{norm_scale_factor} did not produce "
                                   f"integer-like values on a probe slice.")}
        needs_recovery = True
        # Stash for the aggregation loop to use on each per-replicate block.
        # Avoids materializing a 295K x 31K dense recovery matrix.
        lib_sizes_full = np.asarray(adata.obs[lib_size_col].values, dtype=float)

    # Validate obs columns
    for col in [groupby, condition_col, replicate_col, *covariates]:
        if col not in adata.obs.columns:
            return {"status":     "error",
                    "error_type": "MissingObsColumn",
                    "details":    (f"Required obs column '{col}' not found. "
                                   f"Available: {list(adata.obs.columns)}")}

    # ── Aggregate per (replicate × group) ────────────────────────────────
    obs = adata.obs[[groupby, condition_col, replicate_col, *covariates]].copy()
    obs[groupby]       = obs[groupby].astype(str)
    obs[condition_col] = obs[condition_col].astype(str)
    obs[replicate_col] = obs[replicate_col].astype(str)

    is_sparse = sparse.issparse(counts)

    # Paired scRNA designs (same donor under multiple conditions) require a
    # replicate x condition pseudobulk key. Grouping only by donor collapses
    # the two aliquots into one pseudosample and destroys the contrast.
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

    # T1.1: pre-compute total cells per pseudosample ACROSS ALL cell types so
    # each per-group loop iteration can compute the composition log-ratio
    # without re-scanning obs.
    total_cells_per_sample = obs.groupby("_aria_pseudosample_key").size().to_dict()

    groups = sorted(obs[groupby].unique())
    per_group: dict = {}
    csv_rows: list = []
    successful_blocks: list = []

    for group in groups:
        cells_in_group = obs[obs[groupby] == group].index
        if len(cells_in_group) < min_cells_per_pseudosample:
            per_group[group] = {"n_pseudosamples": 0,
                                 "status":          "skipped",
                                 "reason":          f"only {len(cells_in_group)} cells"}
            continue

        # cell indices into the adata
        cell_pos = adata.obs_names.get_indexer(cells_in_group)
        group_counts = counts[cell_pos]
        group_obs    = obs.loc[cells_in_group]

        # Pseudosample key. In paired designs this is donor x condition; in
        # standard cohort designs it remains donor/sample.
        sample_keys = group_obs["_aria_pseudosample_key"].values
        unique_samples = sorted(np.unique(sample_keys))

        # Aggregate sums per pseudosample
        rep_to_sum: dict   = {}
        rep_to_n:   dict   = {}
        rep_to_meta: dict  = {}
        for sample_key in unique_samples:
            mask = sample_keys == sample_key
            n    = int(mask.sum())
            if n < min_cells_per_pseudosample:
                continue
            block = group_counts[mask]
            if needs_recovery:
                # Reverse log1p(normalize_total / scale_factor) on this
                # per-replicate block only. Block size is ~50-500 cells,
                # so memory stays small and we never materialize the full
                # recovered matrix.
                block_dense = block.toarray() if is_sparse else np.asarray(block)
                # cell_pos[mask] maps back into the full adata, indexed
                # against lib_sizes_full
                rep_cell_idx = cell_pos[mask]
                rep_lib      = lib_sizes_full[rep_cell_idx][:, None]
                recovered    = np.expm1(block_dense) * rep_lib / norm_scale_factor
                summed       = np.round(recovered.sum(axis=0)).astype(np.int64)
            elif is_sparse:
                summed = np.asarray(block.sum(axis=0)).ravel().astype(np.int64)
            else:
                summed = np.asarray(block.sum(axis=0)).ravel().astype(np.int64)
            rep_to_sum[sample_key]  = summed
            rep_to_n[sample_key]    = n
            # Use the first cell's metadata for the pseudosample. In paired
            # designs condition is constant because sample_key includes it.
            row              = group_obs[mask].iloc[0]
            meta_entry = {condition_col: row[condition_col],
                          replicate_col: row[replicate_col],
                          **{cv: row[cv] for cv in covariates}}
            if composition_covariate:
                # log( (cells_of_this_type_in_pseudosample + 1) /
                #      (total_cells_in_pseudosample + 1) ).
                # Smoothing keeps the value finite when n is tiny.
                total_for_rep = int(total_cells_per_sample.get(sample_key, n))
                import math as _math
                meta_entry[COMPOSITION_COL] = float(
                    _math.log((n + 1.0) / (total_for_rep + 1.0))
                )
            rep_to_meta[sample_key] = meta_entry

        if len(rep_to_sum) < 2 * min_replicates_per_condition:
            per_group[group] = {"n_pseudosamples": len(rep_to_sum),
                                 "status":          "skipped",
                                 "reason":          (f"only {len(rep_to_sum)} "
                                                     f"pseudosamples survived "
                                                     f"≥{min_cells_per_pseudosample} "
                                                     f"cells filter")}
            continue

        # Build counts (genes × samples) and metadata (samples × covariates)
        sample_order = list(rep_to_sum.keys())
        counts_mat   = np.stack([rep_to_sum[s] for s in sample_order], axis=1)
        meta_df      = pd.DataFrame.from_dict(rep_to_meta, orient="index").loc[sample_order]

        # Drop genes with zero counts across all pseudosamples — DESeq2
        # ignores them anyway but they slow down normalization.
        gene_keep = counts_mat.sum(axis=1) > 0
        counts_df = pd.DataFrame(
            counts_mat[gene_keep],
            index=[g for g, keep in zip(gene_names, gene_keep) if keep],
            columns=sample_order,
        )

        per_group_entry: dict = {
            "n_pseudosamples":     len(sample_order),
            "pseudosample_sizes":  rep_to_n,
            "per_comparison":      {},
        }

        # Run DESeq2 once per comparison.
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds  import DeseqStats

        for comp in comparisons:
            test_lvl, ref_lvl = comp[0], comp[1]
            comp_key  = f"{test_lvl}_vs_{ref_lvl}"
            mask      = meta_df[condition_col].isin([test_lvl, ref_lvl])
            meta_sub  = meta_df[mask].copy()
            counts_sub = counts_df.loc[:, mask].copy()

            n_test = int((meta_sub[condition_col] == test_lvl).sum())
            n_ref  = int((meta_sub[condition_col] == ref_lvl).sum())
            if n_test < min_replicates_per_condition or n_ref < min_replicates_per_condition:
                per_group_entry["per_comparison"][comp_key] = {
                    "status":  "skipped",
                    "reason":  (f"not enough replicates: "
                                f"{test_lvl}={n_test}, {ref_lvl}={n_ref}"),
                }
                continue

            design_factors = [condition_col] + (covariates or [])
            paired_donor_covariate_used = False
            if (
                paired_design
                and auto_paired_donor_covariate
                and replicate_col not in design_factors
            ):
                reps_by_condition = meta_sub.groupby(replicate_col)[
                    condition_col
                ].nunique()
                if len(reps_by_condition) > 0 and bool(
                    (reps_by_condition == 2).all()
                ):
                    design_factors = design_factors + [replicate_col]
                    paired_donor_covariate_used = True
            block_corrected_for_composition = False
            if composition_covariate and COMPOSITION_COL in meta_sub.columns:
                # Only add the covariate if it actually varies across the
                # pseudosamples in this comparison (DESeq2 errors on a
                # constant factor).
                if meta_sub[COMPOSITION_COL].nunique() > 1:
                    design_factors = design_factors + [COMPOSITION_COL]
                    block_corrected_for_composition = True
            design_check = validate_design_matrix(
                meta_sub,
                condition_col=condition_col,
                covariates=[f for f in design_factors if f != condition_col],
                min_replicates_per_condition=min_replicates_per_condition,
            )
            if design_check.get("status") == "blocking":
                per_group_entry["per_comparison"][comp_key] = {
                    "status": "skipped",
                    "reason": "design_matrix_invalid",
                    "design_check": design_check,
                }
                continue
            continuous_factors = design_check.get("continuous_factors", [])
            for col in continuous_factors:
                if col in meta_sub.columns:
                    meta_sub[col] = pd.to_numeric(meta_sub[col], errors="coerce")
            try:
                means = counts_sub.mean(axis=1).astype(float)
                variances = counts_sub.var(axis=1, ddof=1).astype(float)
                disp = ((variances - means) / (means ** 2)).replace(
                    [np.inf, -np.inf], np.nan
                )
                dispersion_estimate = float(
                    np.nanmedian(np.clip(disp.dropna(), 1e-8, None))
                ) if len(disp.dropna()) else 0.1
                mean_expression = float(np.nanmedian(means[means > 0])) \
                    if (means > 0).any() else 1.0
                power_estimate = pseudobulk_power_estimate(
                    n_per_group=(n_test, n_ref),
                    dispersion_estimate=dispersion_estimate,
                    target_log2fc=lfc_min,
                    alpha=padj_max,
                    mean_expression=mean_expression,
                )
                dds_kwargs = {
                    "counts": counts_sub.T,    # samples × genes
                    "metadata": meta_sub,
                    "design_factors": design_factors,
                    "refit_cooks": True,
                }
                if continuous_factors and "continuous_factors" in inspect.signature(DeseqDataSet).parameters:
                    dds_kwargs["continuous_factors"] = continuous_factors
                dds = DeseqDataSet(**dds_kwargs)
                dds.deseq2()
                stat_res = DeseqStats(
                    dds,
                    contrast=[condition_col, test_lvl, ref_lvl],
                )
                stat_res.summary()
                res = stat_res.results_df.dropna(subset=["pvalue"]).copy()
                res["padj_local"] = res.get("padj")
            except Exception as e:
                per_group_entry["per_comparison"][comp_key] = {
                    "status":     "error",
                    "error_type": type(e).__name__,
                    "details":    str(e)[:200],
                }
                continue

            # Low-power flag: even after passing min_replicates_per_condition,
            # n<=2 on either side leaves dispersion estimation noisy. The DE
            # ran, but downstream narrative must caveat the result.
            low_power = (n_test <= 2 or n_ref <= 2)
            low_power_reason = (
                f"n={n_test} vs n={n_ref}: dispersion estimation is unreliable "
                f"with fewer than three replicates per group. DESeq2 produced "
                f"results, but effect-size estimates and FDR are noisy. "
                f"Interpret with caution."
            ) if low_power else None

            per_group_entry["per_comparison"][comp_key] = {
                "status":                    "success",
                "n_significant":             0,
                "n_significant_local":       0,
                "n_significant_global":      0,
                "n_up":                      0,
                "n_up_local":                0,
                "n_up_global":               0,
                "n_down":                    0,
                "n_down_local":              0,
                "n_down_global":             0,
                "n_replicates":              {"test": n_test, "ref": n_ref},
                "low_power_warning":         low_power,
                "low_power_reason":          low_power_reason,
                "dispersion_estimate":        dispersion_estimate,
                "mean_expression_estimate":   mean_expression,
                "power_estimate_at_lfc_min":  power_estimate,
                "corrected_for_composition": block_corrected_for_composition,
                "paired_design":              paired_design,
                "paired_donor_covariate":     paired_donor_covariate_used,
                "design_check":               design_check,
                "top_genes":                 [],
                "all_sig":                   [],
            }
            successful_blocks.append({
                "group":      group,
                "comparison": comp_key,
                "results":    res,
            })

        per_group[group] = per_group_entry

    # ── Global FDR across every gene × group × comparison test ───────────
    global_pvals = []
    for block in successful_blocks:
        res = block["results"]
        pvals = res["pvalue"].astype(float)
        valid = pvals.notna() & np.isfinite(pvals)
        block["valid_mask"] = valid
        global_pvals.extend(pvals[valid].tolist())

    n_tests_global = int(len(global_pvals))
    # F-SCI-POWER (audit 2026-05-28): the empirical BH cutoff is the largest
    # raw p-value whose global-BH adjusted value still clears padj_max. This is
    # the per-test alpha actually applied, which is far stricter than the
    # nominal alpha used for the planning power estimate. Reporting power at
    # this effective alpha reconciles the "power vs decision rule" gap.
    effective_alpha_global = None
    if n_tests_global:
        padj_global = _global_bh(global_pvals)
        _gp = np.asarray(global_pvals, dtype=float)
        _pg = np.asarray(padj_global, dtype=float)
        _passing = _gp[_pg < padj_max]
        effective_alpha_global = float(_passing.max()) if _passing.size else 0.0
        cursor = 0
        for block in successful_blocks:
            res = block["results"]
            valid = block["valid_mask"]
            n_valid = int(valid.sum())
            res["padj_global"] = np.nan
            if n_valid:
                res.loc[valid, "padj_global"] = padj_global[cursor:cursor + n_valid]
                cursor += n_valid

            comp = per_group[block["group"]]["per_comparison"][block["comparison"]]
            sig_local = res[
                (res["padj_local"] < padj_max)
                & (res["log2FoldChange"].abs() > lfc_min)
            ].sort_values("padj_local")
            sig_global = res[
                (res["padj_global"] < padj_max)
                & (res["log2FoldChange"].abs() > lfc_min)
            ].sort_values("padj_global")

            # Primary significance set follows fdr_strategy. 'padj' in each
            # record mirrors the primary adjusted p-value so downstream ORA and
            # narrative use the chosen family; padj_local/padj_global are always
            # carried for audit.
            if fdr_strategy == "global":
                sig_primary = sig_global
                primary_padj_col = "padj_global"
            else:
                sig_primary = sig_local
                primary_padj_col = "padj_local"

            top_records = [
                {"gene":        str(g),
                 "log2fc":      round(float(row["log2FoldChange"]), 3),
                 "pvalue":      float(row["pvalue"]),
                 "padj":        float(row[primary_padj_col]),
                 "padj_local":  float(row["padj_local"]),
                 "padj_global": float(row["padj_global"]),
                 "basemean":    round(float(row["baseMean"]), 1)}
                for g, row in sig_primary.head(top_n).iterrows()
            ]
            all_records = [
                {"gene":        str(g),
                 "log2fc":      round(float(row["log2FoldChange"]), 3),
                 "pvalue":      float(row["pvalue"]),
                 "padj":        float(row[primary_padj_col]),
                 "padj_local":  float(row["padj_local"]),
                 "padj_global": float(row["padj_global"])}
                for g, row in sig_primary.iterrows()
            ]

            comp.update({
                "fdr_strategy":         fdr_strategy,
                "n_significant":        int(len(sig_primary)),
                "n_significant_local":  int(len(sig_local)),
                "n_significant_global": int(len(sig_global)),
                "n_up":                 int((sig_primary["log2FoldChange"] > 0).sum()),
                "n_up_local":           int((sig_local["log2FoldChange"] > 0).sum()),
                "n_up_global":          int((sig_global["log2FoldChange"] > 0).sum()),
                "n_down":               int((sig_primary["log2FoldChange"] < 0).sum()),
                "n_down_local":         int((sig_local["log2FoldChange"] < 0).sum()),
                "n_down_global":        int((sig_global["log2FoldChange"] < 0).sum()),
                "top_genes":            top_records,
                "all_sig":              all_records,
            })

            # Power at the effective (global-BH) threshold, not just nominal.
            if effective_alpha_global and effective_alpha_global > 0:
                nreps = comp.get("n_replicates", {}) or {}
                comp["power_estimate_at_effective_alpha"] = pseudobulk_power_estimate(
                    n_per_group=(nreps.get("test", 0), nreps.get("ref", 0)),
                    dispersion_estimate=comp.get("dispersion_estimate") or 0.1,
                    target_log2fc=lfc_min,
                    alpha=effective_alpha_global,
                    mean_expression=comp.get("mean_expression_estimate") or 100.0,
                )

            for r in all_records:
                csv_rows.append({
                    "group":       block["group"],
                    "comparison":  block["comparison"],
                    "gene":        r["gene"],
                    "log2fc":      r["log2fc"],
                    "pvalue":      r["pvalue"],
                    "padj":        r["padj"],
                    "padj_local":  r["padj_local"],
                    "padj_global": r["padj_global"],
                })

    # ── Write CSV ─────────────────────────────────────────────────────────
    out_dir = Path(output_dir) if output_dir else Path(data_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = str(out_dir / "pseudobulk_de.csv")
    if csv_rows:
        pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
    else:
        # Empty marker file so callers don't trip on a missing path
        pd.DataFrame(columns=[
            "group", "comparison", "gene", "log2fc", "pvalue",
            "padj", "padj_local", "padj_global",
        ]) \
          .to_csv(csv_path, index=False)

    # Count provenance (audit 2026-05-28, F-SCI-LOGNORM): be explicit about
    # whether DESeq2 saw genuine raw counts or counts reverse-engineered from
    # log-normalized values, so the report can disclose it honestly.
    if needs_recovery:
        count_source = "recovered_from_lognorm"
    elif use_raw and adata.raw is not None:
        count_source = "raw_counts"
    else:
        count_source = "X_counts"

    return {
        "status":                          "success",
        "n_groups":                        len(per_group),
        "groupby":                         groupby,
        "condition_col":                   condition_col,
        "replicate_col":                   replicate_col,
        "covariates":                      covariates,
        "count_source":                    count_source,
        "lognorm_recovered":               bool(needs_recovery),
        "norm_scale_factor_used":          (norm_scale_factor if needs_recovery else None),
        "lognorm_lib_size_col":            (lib_size_col if needs_recovery else None),
        "composition_covariate_requested": composition_covariate,
        "paired_design":                   paired_design,
        "auto_paired_donor_covariate":      auto_paired_donor_covariate,
        "background_genes": background_genes,
        "background_size":  len(background_genes),
        "background_source": "dataset_expressed_genes",
        "multiple_testing": {
            "fdr_strategy":   fdr_strategy,
            "primary_family": ("per-cluster BH" if fdr_strategy == "per_cluster"
                               else "global pooled BH"),
            "local_method":   "BH",
            "global_method":  "BH",
            "n_tests_global": n_tests_global,
        },
        # F-SCI-POWER (audit 2026-05-28): power is an approximate NB-Wald value
        # at the NOMINAL per-test alpha. Significance is actually declared with
        # global BH-FDR across all blocks, whose effective per-test threshold is
        # far stricter, so the reported power is an UPPER BOUND on the power at
        # the applied threshold. Disclosed so Methods can state it honestly.
        "power": {
            "alpha_nominal":  padj_max,
            "effective_alpha_global": effective_alpha_global,
            "method":         "approximate NB-Wald (closed form)",
            "applied_threshold": "global BH-FDR",
            "note": (
                "power_estimate_at_lfc_min is computed at the nominal per-test "
                "alpha and is an UPPER BOUND. Significance is declared with "
                "global BH-FDR across all blocks; the empirical per-test cutoff "
                "is effective_alpha_global. power_estimate_at_effective_alpha "
                "reports power at that stricter, actually-applied threshold."
            ),
        },
        "thresholds":     {"padj_max": padj_max, "lfc_min": lfc_min,
                           "min_cells_per_pseudosample": min_cells_per_pseudosample,
                           "min_replicates_per_condition": min_replicates_per_condition},
        "per_group":      per_group,
        "output_csv":     csv_path,
    }


if __name__ == "__main__":
    run_script(rna_pseudobulk_de)
