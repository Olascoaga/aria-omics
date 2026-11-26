"""Bulk RNA-seq DESeq2 core: design formula, covariate resolution, LFC shrink,
fit-warning serialization, the DESeq2 run and the mock result (A7 split).

Re-exported from aria.scripts.rna_bulk_de."""
from __future__ import annotations
import warnings


def _quote_formula_term(name: str) -> str:
    """Quote non-identifier column names for formulaic using backticks."""
    import keyword
    text = str(name)
    if text.isidentifier() and not keyword.iskeyword(text):
        return text
    escaped = text.replace("`", "\\`")
    return f"`{escaped}`"


def _can_wrap_in_formulaic_c(name: str) -> bool:
    """True when formulaic safely parses ``C(<term>, Treatment(...))``."""
    import keyword
    text = str(name)
    return "." not in text and not keyword.iskeyword(text)


def _build_design_formula(
    design_factor: str,
    covariates: list,
    *,
    ref_level: str | None = None,
) -> str:
    """Build a DESeq2 design with covariates first and the factor of interest
    last (DESeq2/pydeseq2 convention). Covariates are deduped and never repeat
    the factor of interest. P0-4: honors the confirmed `~ batch + condition`
    design instead of a hardcoded `~ condition`.

    When ``ref_level`` is supplied, encodes the explicit formulaic treatment
    reference ``C(<factor>, Treatment(reference=<ref_level>))`` so pydeseq2
    0.5+ (where ``ref_level=`` is deprecated and ignored) uses the exact
    contrast denominator as the baseline level.
    """
    terms = [
        _quote_formula_term(c)
        for c in dict.fromkeys(covariates or [])
        if c and c != design_factor
    ]
    q_factor = _quote_formula_term(design_factor)
    if ref_level is not None and _can_wrap_in_formulaic_c(design_factor):
        factor_term = f"C({q_factor}, Treatment(reference={str(ref_level)!r}))"
    else:
        factor_term = q_factor
    return "~ " + " + ".join(terms + [factor_term])


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

    Matches ONLY the exact ``[T.<test_lvl>]`` suffix on ``condition_col``
    (plain, backtick-quoted, or formulaic ``C(<condition_col>, ...)[T.<test_lvl>]``).
    Never falls back to substring containment (which could match ``[T.KO_ctrl]``
    when ``test_lvl == 'KO'`` and invert every LFC sign).
    """
    try:
        cols = [str(c) for c in dds.varm["LFC"].columns]
    except Exception:
        return None
    cond = str(condition_col)
    q_cond = _quote_formula_term(cond)
    suffix = f"[T.{test_lvl}]"
    cond_tokens = {cond, q_cond, f"`{cond}`"}
    for c in cols:
        if not c.endswith(suffix):
            continue
        lhs = c[: -len(suffix)].strip()
        if lhs in cond_tokens:
            return c
    for c in cols:
        if not c.endswith(suffix):
            continue
        lhs = c[: -len(suffix)].strip()
        if lhs.startswith("C(") and lhs.endswith(")"):
            inner = lhs[2:-1].strip()
            for tok in sorted(cond_tokens, key=len, reverse=True):
                if inner.startswith(tok):
                    rest = inner[len(tok):].lstrip()
                    if not rest or rest.startswith(","):
                        return c
    return None


def _detect_shrink_sign_flip(raw_lfc, shrunk_lfc) -> bool:
    """Guard against shrinking the inverted contrast coefficient.

    Returns True when a clear majority of well-separated raw MLEs flipped sign
    after ``lfc_shrink`` (something apeGLM shrinking toward 0 never does on the
    correctly oriented coefficient).
    """
    import numpy as np

    raw = np.asarray(raw_lfc, dtype=float)
    shrunk = np.asarray(shrunk_lfc, dtype=float)
    if raw.shape != shrunk.shape or raw.size == 0:
        return False
    clear = (
        np.isfinite(raw)
        & np.isfinite(shrunk)
        & (np.abs(raw) >= 0.5)
        & (np.abs(shrunk) >= 0.1)
    )
    n_clear = int(clear.sum())
    if n_clear == 0:
        return False
    flipped = np.sign(raw[clear]) != np.sign(shrunk[clear])
    flip_frac = float(flipped.mean())
    if n_clear >= 3 and flip_frac > 0.5:
        return True
    if flip_frac == 1.0 and float(np.abs(raw[clear]).max()) >= 1.0:
        return True
    return False


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
                expose_convergence: bool = False,
                external_size_factors=None) -> tuple:
    """
    Run DESeq2 via pydeseq2 with correct design factor.
    Returns (result_dict, warnings_list).

    min_replicates_per_condition: production floor. Default 3. n=2 is
    permitted only when the caller explicitly lowers this (e.g. for
    pilot/legacy datasets); the resulting block is flagged with
    low_power_warning in the return dict.

    external_size_factors: optional Series indexed by sample ID with positive
    finite heterologous spike-in size factors (ADR-081). When supplied, bypasses
    median-of-ratios ``fit_size_factors()`` while keeping the identical
    formulaic ``Treatment(reference=denominator)`` design, continuous/categorical
    covariate handling, Wald ``lfcThreshold`` test, and apeGLM sign-flip guard.
    """
    import numpy as np
    import pandas as pd
    warnings = []

    num_str = str(numerator)
    den_str = str(denominator)
    if not metadata.index.is_unique or not counts.columns.is_unique:
        dup_m = sorted(set(metadata.index[metadata.index.duplicated()].astype(str).tolist()))
        dup_c = sorted(set(counts.columns[counts.columns.duplicated()].astype(str).tolist()))
        return {
            "status": "error",
            "error_type": "DuplicateSampleIds",
            "details": (
                "DESeq2 requires unique sample identifiers in both metadata "
                f"and counts (duplicate metadata IDs={dup_m}, "
                f"duplicate count columns={dup_c})."
            ),
        }, warnings
    # Filter to comparison groups only
    mask      = metadata[design_factor].astype(str).isin([num_str, den_str])
    meta_sub  = metadata[mask].copy()
    meta_sub[design_factor] = meta_sub[design_factor].astype(str)
    counts_sub = counts[meta_sub.index].copy()

    group_sizes = meta_sub[design_factor].value_counts().to_dict()
    n_num = int(group_sizes.get(num_str, 0))
    n_den = int(group_sizes.get(den_str, 0))

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

    sf_sub = None
    if external_size_factors is not None:
        sf_raw = external_size_factors.reindex(meta_sub.index).astype(float)
        sf_arr = sf_raw.to_numpy(dtype=float)
        if sf_arr.size == 0 or (~np.isfinite(sf_arr)).any() or (sf_arr <= 0).any():
            return {
                "status": "error",
                "error_type": "InvalidSpikeInSizeFactors",
                "details": (
                    "External spike-in size factors must be finite and > 0 "
                    "for all contrast samples."
                ),
            }, warnings
        # Re-center subset size factors to geometric mean = 1.0 (R DESeq2 semantics).
        sf_sub = sf_raw / float(np.exp(np.mean(np.log(sf_arr))))
        if not np.isfinite(sf_sub.to_numpy(dtype=float)).all() or (sf_sub <= 0).any():
            return {
                "status": "error",
                "error_type": "InvalidSpikeInSizeFactors",
                "details": "Re-centred spike-in size factors are not finite and positive.",
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
        # interest last). Keep the readable formula for reporting, and encode
        # the explicit Treatment(reference=denominator) level in the formulaic
        # design passed to pydeseq2 0.5+ (where ref_level= is ignored).
        design_formula = _build_design_formula(design_factor, usable_covariates)
        model_design = _build_design_formula(
            design_factor, usable_covariates, ref_level=den_str
        )

        # Ensure categorical/continuous covariate dtypes match design_check, and
        # order design_factor with denominator first so both formulaic C(...)
        # and plain categorical encodings orient the contrast as num vs den.
        for col in design_check.get("continuous_factors", []) or []:
            if col in meta_sub.columns and col != design_factor:
                meta_sub[col] = pd.to_numeric(meta_sub[col], errors="coerce")
        for col in design_check.get("categorical_factors", []) or []:
            if col in meta_sub.columns and col != design_factor:
                meta_sub[col] = meta_sub[col].astype(str)
        meta_sub[design_factor] = pd.Categorical(
            meta_sub[design_factor].astype(str),
            categories=[den_str, num_str],
        )

        # Legacy pre-0.5 fallback kwargs only (pydeseq2 0.5+ ignores ref_level).
        ref_kwargs = {"ref_level": [design_factor, den_str]} if lfc_shrink else {}

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
                design=model_design,
                refit_cooks=True,
                quiet=True,
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
            "contrast": [design_factor, num_str, den_str],
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
            if sf_sub is not None:
                # Inject external heterologous spike-in size factors directly,
                # bypassing median-of-ratios fit_size_factors() so global
                # unidirectional chromatin shifts are preserved (ADR-081).
                sf_aligned = sf_sub.reindex(dds.obs_names).to_numpy(dtype=float)
                dds.obs["size_factors"] = sf_aligned
                dds.layers["normed_counts"] = dds.X / sf_aligned[:, None]
                dds.var["_normed_means"] = dds.layers["normed_counts"].mean(axis=0)
                dds.fit_genewise_dispersions()
                dds.fit_dispersion_trend()
                dds.fit_dispersion_prior()
                dds.fit_MAP_dispersions()
                dds.fit_LFC()
                dds.calculate_cooks()
                if dds.refit_cooks:
                    dds.refit()
                dds.cooks_outlier()
            else:
                # P1-1b: put the biological effect-size threshold inside the Wald test
                # instead of filtering |log2FC| after BH, matching DESeq2's lfcThreshold.
                dds.deseq2()
            stat_res = DeseqStats(dds, **wald_kwargs)
            stat_res.summary()
        warnings.extend(_serialize_fit_warnings(_fit_warns))

        # P1-1/ADR-023: apeGLM LFC shrinkage. Keep the raw MLE LFC, shrink in
        # place, and fall back to raw if the coefficient is unavailable,
        # shrinkage raises, or the post-shrink sign-flip guard triggers.
        shrink_applied = False
        shrink_reason = None
        raw_lfc = stat_res.results_df["log2FoldChange"].copy()
        raw_lfc_se = (
            stat_res.results_df["lfcSE"].copy()
            if "lfcSE" in stat_res.results_df.columns
            else None
        )
        if lfc_shrink:
            coeff = _shrink_coeff(dds, design_factor, num_str)
            if coeff is None:
                shrink_reason = "apeGLM coefficient not found"
            else:
                try:
                    stat_res.lfc_shrink(coeff=coeff)
                    if _detect_shrink_sign_flip(
                        raw_lfc, stat_res.results_df["log2FoldChange"]
                    ):
                        stat_res.results_df["log2FoldChange"] = raw_lfc
                        if raw_lfc_se is not None:
                            stat_res.results_df["lfcSE"] = raw_lfc_se
                        shrink_reason = (
                            "lfc_shrink sign-flip detected; reverted to raw MLE"
                        )
                        warnings.append(
                            f"apeGLM LFC shrinkage for '{num_str} vs {den_str}' "
                            f"flipped the sign of clear raw MLEs on coeff "
                            f"'{coeff}'; reverted to unshrunken MLE log2FoldChange."
                        )
                    else:
                        shrink_applied = True
                except Exception as _sx:
                    shrink_reason = f"lfc_shrink failed: {str(_sx)[:120]}"
        else:
            shrink_reason = "disabled"

        results_df = stat_res.results_df.dropna(subset=["padj"]).copy()
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

        padj_vals = pd.to_numeric(results_df["padj"], errors="coerce")
        lfc_vals = pd.to_numeric(results_df["log2FoldChange"], errors="coerce")
        sig = results_df[
            np.isfinite(padj_vals)
            & (padj_vals < padj_thr)
            & np.isfinite(lfc_vals)
        ]
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

        out = {
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
        }
        if sf_sub is not None:
            out["size_factors_used"] = {
                str(k): round(float(v), 6) for k, v in dds.obs["size_factors"].items()
            }
            out["normalization_method"] = "heterologous_spikein_size_factors"
        return out, warnings

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







