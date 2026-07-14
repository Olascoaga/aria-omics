"""Bulk RNA-seq DESeq2 core: design formula, covariate resolution, LFC shrink,
fit-warning serialization, the DESeq2 run and the mock result (A7 split).

Re-exported from aria.scripts.rna_bulk_de."""
from __future__ import annotations
import warnings


def _build_design_formula(design_factor: str, covariates: list) -> str:
    """Build a DESeq2 design with covariates first and the factor of interest
    last (DESeq2/pydeseq2 convention). Covariates are deduped and never repeat
    the factor of interest. P0-4: honors the confirmed `~ batch + condition`
    design instead of a hardcoded `~ condition`."""
    terms = [c for c in dict.fromkeys(covariates or []) if c and c != design_factor]
    return "~ " + " + ".join(terms + [design_factor])


def _resolve_covariates(metadata, design_factor: str, covariates: list) -> tuple:
    """Keep only covariates usable in this contrast subset: present as a column,
    distinct from the factor of interest, and varying (>=2 non-null levels).
    Returns (usable, dropped) where dropped is a list of (name, reason). A
    confirmed-but-unusable covariate is DISCLOSED, never silently ignored."""
    usable: list = []
    dropped: list = []
    for cov in dict.fromkeys(covariates or []):
        if not cov or cov == design_factor:
            continue
        if cov not in metadata.columns:
            dropped.append((cov, "not present in the sample metadata"))
            continue
        if metadata[cov].dropna().nunique() < 2:
            dropped.append(
                (cov, "constant within the contrast (no levels to adjust for)")
            )
            continue
        usable.append(cov)
    return usable, dropped


def _shrink_coeff(dds, condition_col: str, test_lvl: str):
    """Find the apeGLM LFC coefficient column for test-vs-ref (P1-1 / ADR-023).

    pydeseq2 names coefficients patsy-style, e.g. ``condition[T.treat]``. With the
    dds reference fixed to the contrast's ref level, the only non-reference level
    is ``test_lvl``, so that column is exactly the test-vs-ref effect.
    """
    try:
        cols = [str(c) for c in dds.varm["LFC"].columns]
    except Exception:
        return None
    exact = f"{condition_col}[T.{test_lvl}]"
    if exact in cols:
        return exact
    for c in cols:
        if c.startswith(str(condition_col)) and str(test_lvl) in c:
            return c
    return None


# F8 (preprint audit 2026-06-19): warning categories that are benign third-party
# API churn (NOT numerical/convergence signal). Everything else emitted during the
# pydeseq2 fit is surfaced into the result's audit trail instead of being silenced.
_BENIGN_WARNING_CATEGORIES = (
    DeprecationWarning, PendingDeprecationWarning, FutureWarning,
    ImportWarning, ResourceWarning,
)


def _serialize_fit_warnings(caught) -> list:
    """Turn captured fit warnings into audit-trail strings, dropping benign ones.

    Non-benign warnings (RuntimeWarning/UserWarning and any pydeseq2/numpy
    convergence or numerical warning) are kept verbatim so a reviewer sees them in
    `result["warnings"]`; benign API-churn categories are dropped to avoid noise.
    """
    out = []
    for wm in caught or []:
        if issubclass(wm.category, _BENIGN_WARNING_CATEGORIES):
            continue
        out.append(f"pydeseq2 fit warning [{wm.category.__name__}]: {str(wm.message)[:300]}")
    return out


def _run_deseq2(counts, metadata, design_factor: str,
                numerator: str, denominator: str,
                padj_thr: float, lfc_thr: float,
                allow_mock: bool = False,
                min_replicates_per_condition: int = 3,
                covariates: list = None,
                lfc_shrink: bool = True,
                n_cpus: int = None,
                expose_convergence: bool = False) -> tuple:
    """
    Run DESeq2 via pydeseq2 with correct design factor.
    Returns (result_dict, warnings_list).

    min_replicates_per_condition: production floor. Default 3. n=2 is
    permitted only when the caller explicitly lowers this (e.g. for
    pilot/legacy datasets); the resulting block is flagged with
    low_power_warning in the return dict.
    """
    import pandas as pd
    warnings = []

    # Filter to comparison groups only
    mask      = metadata[design_factor].isin([numerator, denominator])
    meta_sub  = metadata[mask].copy()
    counts_sub = counts[meta_sub.index].copy()

    group_sizes = meta_sub[design_factor].value_counts().to_dict()
    n_num = int(group_sizes.get(numerator, 0))
    n_den = int(group_sizes.get(denominator, 0))

    if n_num < min_replicates_per_condition or n_den < min_replicates_per_condition:
        return {
            "status":     "error",
            "error_type": "InsufficientReplicates",
            "details":    (
                f"DESeq2 requires at least {min_replicates_per_condition} "
                f"replicates per group. Found {group_sizes}. Lower "
                f"min_replicates_per_condition explicitly (e.g. to 2) to run "
                f"with a low_power_warning instead."
            ),
        }, warnings

    # P0-4: honor confirmed covariates (e.g. batch). Keep only those usable in
    # this subset; disclose any confirmed covariate we cannot adjust for.
    usable_covariates, dropped_covariates = _resolve_covariates(
        meta_sub, design_factor, covariates
    )
    for cov, reason in dropped_covariates:
        warnings.append(
            f"Confirmed covariate '{cov}' was NOT adjusted for: {reason}. "
            f"The fitted model does not control for it."
        )

    from aria.utils.design_matrix import validate_design_matrix
    design_check = validate_design_matrix(
        meta_sub,
        condition_col=design_factor,
        covariates=usable_covariates,
        min_replicates_per_condition=min_replicates_per_condition,
    )
    warnings.extend(
        f"Design matrix {issue['severity']}: {issue['message']}"
        for issue in design_check.get("issues", [])
    )
    if design_check.get("status") == "blocking":
        return {
            "status": "error",
            "error_type": "InvalidDesignMatrix",
            "details": "; ".join(
                issue["message"] for issue in design_check.get("issues", [])
                if issue.get("severity") == "blocking"
            )[:500],
            "design_check": design_check,
        }, warnings

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        from aria.utils.power_estimation import bulk_power_estimate
        import warnings as w

        means = counts_sub.mean(axis=1).astype(float)
        variances = counts_sub.var(axis=1, ddof=1).astype(float)
        disp = ((variances - means) / (means ** 2)).replace(
            [float("inf"), float("-inf")], float("nan")
        )
        disp = disp.dropna()
        dispersion_estimate = float(disp[disp > 0].median()) \
            if (disp > 0).any() else 0.1
        mean_expression = float(means[means > 0].median()) \
            if (means > 0).any() else 1.0
        power_estimate = bulk_power_estimate(
            n_per_group=(n_num, n_den),
            mean_expression=mean_expression,
            dispersion=dispersion_estimate,
            target_log2fc=lfc_thr,
            alpha=padj_thr,
        )

        # P0-4: fit the covariate-adjusted design (covariates first, factor of
        # interest last). The contrast still names `design_factor` explicitly.
        design_formula = _build_design_formula(design_factor, usable_covariates)

        # P1-1/ADR-023: fix the reference to the contrast's denominator so the
        # apeGLM coefficient (design_factor[T.numerator]) is exactly the
        # numerator-vs-denominator effect we shrink. Only when shrinkage is on,
        # so the disabled path is byte-identical to the legacy behavior.
        ref_kwargs = {"ref_level": [design_factor, denominator]} if lfc_shrink else {}

        # Optional multi-core: pydeseq2 parallelizes the per-gene dispersion/LFC
        # fits over n_cpus (DeseqStats reuses dds.inference, so the Wald tests
        # follow). Default None leaves the established RNA path byte-identical;
        # only callers that explicitly request it (e.g. the scATAC pseudobulk
        # validation harness) pay for / benefit from the extra cores.
        par_kwargs = {"n_cpus": int(n_cpus)} if n_cpus else {}

        # pydeseq2 expects samples × genes. The public API changed from
        # design_factors=... to design="~ factor"; support both.
        try:
            dds = DeseqDataSet(
                counts=counts_sub.T,
                metadata=meta_sub,
                design=design_formula,
                refit_cooks=True,
                quiet=True,
                **ref_kwargs,
                **par_kwargs,
            )
        except TypeError:
            dds = DeseqDataSet(
                counts=counts_sub.T,
                metadata=meta_sub,
                design_factors=usable_covariates + [design_factor],
                refit_cooks=True,
                **ref_kwargs,
                **par_kwargs,
            )
        # F8: capture pydeseq2 convergence/numeric warnings during the fit into the
        # audit trail instead of silencing them globally. simplefilter("always")
        # records every warning emitted in this scope; non-benign ones are surfaced
        # in `warnings` below. This changes NO estimation — only warning handling.
        wald_kwargs = {
            "contrast": [design_factor, numerator, denominator],
            "alpha": padj_thr,
            "quiet": True,
        }
        lfc_threshold_test = {
            "requested": float(lfc_thr),
            "applied": False,
            "alt_hypothesis": None,
            "lfc_null": 0.0,
            "significance_rule": "padj < threshold",
        }
        if float(lfc_thr) > 0:
            wald_kwargs.update({
                "lfc_null": float(lfc_thr),
                "alt_hypothesis": "greaterAbs",
            })
            lfc_threshold_test.update({
                "applied": True,
                "alt_hypothesis": "greaterAbs",
                "lfc_null": float(lfc_thr),
            })
        with w.catch_warnings(record=True) as _fit_warns:
            w.simplefilter("always")
            # P1-1b: put the biological effect-size threshold inside the Wald test
            # instead of filtering |log2FC| after BH, matching DESeq2's lfcThreshold.
            dds.deseq2()
            stat_res = DeseqStats(dds, **wald_kwargs)
            stat_res.summary()
        warnings.extend(_serialize_fit_warnings(_fit_warns))

        # P1-1/ADR-023: apeGLM LFC shrinkage. Keep the raw MLE LFC, shrink in
        # place, and fall back to raw if the coefficient is unavailable or
        # shrinkage raises (never break the DE). apeGLM leaves p-values unchanged.
        shrink_applied = False
        shrink_reason = None
        raw_lfc = stat_res.results_df["log2FoldChange"].copy()
        if lfc_shrink:
            coeff = _shrink_coeff(dds, design_factor, numerator)
            if coeff is None:
                shrink_reason = "apeGLM coefficient not found"
            else:
                try:
                    stat_res.lfc_shrink(coeff=coeff)
                    shrink_applied = True
                except Exception as _sx:
                    shrink_reason = f"lfc_shrink failed: {str(_sx)[:120]}"
        else:
            shrink_reason = "disabled"

        results_df = stat_res.results_df.dropna(subset=["padj"])
        results_df["log2FoldChange_raw"] = raw_lfc.reindex(results_df.index)

        # T15 (scATAC P0 W0.4): OPT-IN, additive exposure of pydeseq2's per-gene
        # LFC-fit convergence flag (`dds.var["_LFC_converged"]`). Default OFF, so
        # bulk RNA n_sig/n_up/n_down and every existing consumer are byte-identical;
        # the scATAC DA lane turns it on to gate non-converged peaks (apeGLM LFC≈0
        # with a significant Wald p). This does NOT change DESeq2 estimation.
        if expose_convergence and "_LFC_converged" in dds.var.columns:
            results_df["lfc_converged"] = (
                dds.var["_LFC_converged"].reindex(results_df.index)
            )

        sig = results_df[results_df["padj"] < padj_thr]
        # Sort by padj ascending (most significant first) so that
        # downstream "top N" slicing returns the most reliable hits,
        # not the most extreme log2FC (which can be noisy at low counts).
        sig = sig.sort_values("padj")

        up_genes   = list(sig[sig["log2FoldChange"] > 0].index)
        down_genes = list(sig[sig["log2FoldChange"] < 0].index)

        if len(sig) == 0:
            warnings.append(
                f"No significant genes found at padj < {padj_thr} "
                f"and |log2FC| > {lfc_thr}. "
                f"Consider relaxing thresholds or checking data quality."
            )

        low_power = (n_num <= 2 or n_den <= 2)
        low_power_reason = (
            f"n={n_num} ({numerator}) vs n={n_den} ({denominator}): "
            f"dispersion estimation is unreliable with fewer than three "
            f"replicates per group. DESeq2 produced results, but effect-size "
            f"estimates and FDR are noisy. Interpret with caution."
        ) if low_power else None
        if low_power:
            warnings.append(low_power_reason)

        return {
            "status":            "success",
            "results":           results_df,
            "n_sig":             len(sig),
            "n_up":              len(up_genes),
            "n_down":            len(down_genes),
            "sig_genes":         list(sig.index),
            "up_genes":          up_genes,
            "down_genes":        down_genes,
            "n_replicates":      {"test": n_num, "ref": n_den},
            "low_power_warning": low_power,
            "low_power_reason":  low_power_reason,
            "dispersion_estimate": dispersion_estimate,
            "mean_expression_estimate": mean_expression,
            "power_estimate_at_lfc_min": power_estimate,
            "design_check":      design_check,
            "fitted_design_formula": design_formula,
            "covariates_adjusted":   usable_covariates,
            "covariates_dropped":    [
                {"covariate": c, "reason": r} for c, r in dropped_covariates
            ],
            "lfc_shrinkage": {
                "requested": bool(lfc_shrink),
                "applied":   shrink_applied,
                "method":    "apeGLM" if shrink_applied else None,
                "reason":    shrink_reason,
            },
            "lfc_threshold_test": lfc_threshold_test,
        }, warnings

    except ImportError:
        if allow_mock:
            warnings.append(
                "pydeseq2 not installed — returning mock DE results "
                "(explicit mock mode)."
            )
            return _mock_de_result(padj_thr, lfc_thr), warnings
        return {
            "status":     "error",
            "error_type": "MissingDependency",
            "details":    (
                "pydeseq2 is required for bulk RNA differential expression. "
                "Install it in the RNA environment or rerun with explicit "
                "allow_mock=true for development only."
            ),
        }, warnings

    except Exception as e:
        return {
            "status":     "error",
            "error_type": "DESeq2Failed",
            "details":    str(e)[:500],
        }, warnings + [f"DESeq2 error: {e}"]


def _mock_de_result(padj_thr, lfc_thr) -> dict:
    """Mock DE result for environments without pydeseq2."""
    import pandas as pd
    import numpy as np
    rng   = np.random.default_rng(42)
    genes = [f"GENE_{i:04d}" for i in range(100)]
    mock_df = pd.DataFrame({
        "log2FoldChange": rng.normal(0, 2, 100),
        "padj":           rng.uniform(0, 0.2, 100),
        "pvalue":         rng.uniform(0, 0.1, 100),
        "baseMean":       rng.exponential(500, 100),
    }, index=genes)

    sig  = mock_df[(mock_df["padj"] < padj_thr) &
                   (mock_df["log2FoldChange"].abs() > lfc_thr)]
    up   = list(sig[sig["log2FoldChange"] > 0].index)
    down = list(sig[sig["log2FoldChange"] < 0].index)

    return {
        "status":    "success",
        "results":   mock_df,
        "n_sig":     len(sig),
        "n_up":      len(up),
        "n_down":    len(down),
        "sig_genes": list(sig.index),
        "up_genes":  up,
        "down_genes": down,
        "note":      "mock — pydeseq2 not installed",
    }


# ── Pathway enrichment ────────────────────────────────────────────────────────























# ── Visualizations ────────────────────────────────────────────────────────────







