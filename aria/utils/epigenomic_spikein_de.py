"""External heterologous spike-in size-factor adapter over PyDESeq2 (ADR-081).

Why this module exists:
  In ChIP-Rx, spike-in CUT&RUN (e.g. *S. cerevisiae* ``sacCer3`` / *E. coli*
  carryover DNA), and spike-in CUT&TAG, global unidirectional chromatin loss
  (e.g. EZH2 inhibition wiping out 90% of H3K27me3) violates median-of-ratios
  normalization because median-of-ratios assumes most peaks are invariant.
  Bioconductor DESeq2 handles this via ``sizeFactors(dds) <- s`` where
  ``s_r = N_spike,r / geometric_mean(N_spike)`` (equivalently
  ``(1 / scale_r) / geometric_mean(1 / scale)`` for coverage scale factors
  ``scale_r = C / N_spike,r``).

  This helper preserves ``aria/scripts/rna_bulk/deseq2.py`` at zero diff while
  reproducing DESeq2's exact fixed-size-factor pipeline in ``pydeseq2``:
    1. Inject ``dds.obs["size_factors"] = s`` (normalized to geometric mean 1.0),
       ``dds.layers["normed_counts"] = dds.X / s[:, None]``, and
       ``dds.var["_normed_means"]`` WITHOUT running median-of-ratios
       ``fit_size_factors()``.
    2. Run gene-wise dispersion, trend, prior, MAP dispersion, LFC, and Cook's
       outlier refit (``DeseqDataSet._refit_without_outliers`` preserves
       ``obs["size_factors"]`` automatically).
    3. Run ``DeseqStats`` with ``lfc_null`` / ``greaterAbs`` and apeGLM
       shrinkage identically to ``_run_deseq2``.

Boundary invariant: imported inside ``aria-rna-env``; never imports
``aria.agents.*``, ``litellm``, or ``pydantic``.
"""

from __future__ import annotations

import math
import warnings as w
from typing import Any

import numpy as np
import pandas as pd


def resolve_replicate_spikein_size_factors(
    raw_metadata: pd.DataFrame,
    rep_meta: pd.DataFrame,
    condition_col: str,
    replicate_col: str,
    *,
    spikein_counts: dict[str, float] | None = None,
    spikein_scale_factors: dict[str, float] | None = None,
    spikein_counts_col: str | None = None,
    spikein_scale_col: str | None = None,
) -> tuple[pd.Series | None, dict[str, Any]]:
    """Aggregate sample-level spike-in exposures to biological replicates and
    compute geometric-mean-1 DESeq2 ``size_factors``.

    Aggregation invariant (Codex Round-1 review):
      - Raw peak counts are summed across technical replicates ``i in r``.
      - Therefore spike-in read counts ``N_spike,i`` MUST also be summed across
        technical replicates first: ``N_spike,r = sum_{i in r} N_spike,i``.
      - If inverse coverage scale factors ``a_i = C / N_spike,i`` are supplied
        instead of raw spike-in read counts, individual exposures ``1 / a_i``
        are summed across technical replicates first: ``E_r = sum_{i in r} (1 / a_i)``,
        then ``s_r = E_r / geometric_mean(E)``.
    """
    meta = raw_metadata.copy()
    meta["_replicate_sample"] = (
        meta[condition_col].astype(str) + "::" + meta[replicate_col].astype(str)
    )

    source_kind: str | None = None
    sample_exposures: dict[str, float] = {}

    if isinstance(spikein_counts, dict) and spikein_counts:
        source_kind = "explicit_spikein_counts"
        for k, v in spikein_counts.items():
            if v is not None:
                fv = float(v)
                if fv <= 0 or not math.isfinite(fv):
                    raise ValueError(f"Invalid non-positive spikein_counts for {k}: {v!r}")
                sample_exposures[str(k)] = fv
    elif isinstance(spikein_scale_factors, dict) and spikein_scale_factors:
        source_kind = "explicit_spikein_scale_factors"
        for k, v in spikein_scale_factors.items():
            if v is not None:
                fv = float(v)
                if fv <= 0 or not math.isfinite(fv):
                    raise ValueError(f"Invalid non-positive spikein_scale_factor for {k}: {v!r}")
                sample_exposures[str(k)] = 1.0 / fv
    else:
        count_col_candidates = [
            c for c in (spikein_counts_col, "spikein_reads", "spikein_counts", "spike_in_reads")
            if c and c in meta.columns
        ]
        scale_col_candidates = [
            c for c in (spikein_scale_col, "spikein_scale_factor", "scale_factor")
            if c and c in meta.columns
        ]
        if count_col_candidates:
            col = count_col_candidates[0]
            raw_col = meta[col]
            if raw_col.notna().any():
                vals = pd.to_numeric(raw_col, errors="coerce")
                if not (vals.notna().all() and np.isfinite(vals).all() and (vals > 0).all()):
                    raise ValueError(
                        f"Metadata spike-in column {col!r} contains missing, non-numeric, "
                        "or non-positive values; refusing silent fallback to conventional normalization."
                    )
                source_kind = f"metadata_column:{col}"
                sample_exposures = {str(k): float(v) for k, v in vals.items()}
        elif scale_col_candidates:
            col = scale_col_candidates[0]
            raw_col = meta[col]
            if raw_col.notna().any():
                vals = pd.to_numeric(raw_col, errors="coerce")
                if not (vals.notna().all() and np.isfinite(vals).all() and (vals > 0).all()):
                    raise ValueError(
                        f"Metadata spike-in scale column {col!r} contains missing, non-numeric, "
                        "or non-positive values; refusing silent fallback to conventional normalization."
                    )
                source_kind = f"metadata_scale_column:{col}"
                sample_exposures = {str(k): 1.0 / float(v) for k, v in vals.items()}

    if not sample_exposures:
        return None, {"used": False, "source": None}

    missing = [str(s) for s in meta.index if str(s) not in sample_exposures]
    if missing:
        raise ValueError(
            f"Spike-in normalization requested ({source_kind}) but missing values "
            f"for sample(s): {missing}"
        )

    rep_exposures: dict[str, float] = {}
    for rep_id, rows in meta.groupby("_replicate_sample", sort=False):
        exp_sum = sum(sample_exposures[str(sid)] for sid in rows.index)
        if exp_sum <= 0 or not math.isfinite(exp_sum):
            raise ValueError(f"Non-positive aggregated spike-in exposure for {rep_id}: {exp_sum}")
        rep_exposures[str(rep_id)] = float(exp_sum)

    ordered_reps = [str(r) for r in rep_meta.index]
    exp_vec = np.array([rep_exposures[r] for r in ordered_reps], dtype=float)
    geo_mean = float(np.exp(np.mean(np.log(exp_vec))))
    sf_vec = exp_vec / geo_mean
    sf_series = pd.Series(sf_vec, index=rep_meta.index, dtype=float)

    return sf_series, {
        "used": True,
        "method": "heterologous_spikein_size_factors",
        "source": source_kind,
        "replicate_size_factors": {
            str(k): round(float(v), 6) for k, v in sf_series.items()
        },
    }


def run_deseq2_with_external_size_factors(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    external_size_factors: pd.Series,
    design_factor: str,
    numerator: str,
    denominator: str,
    padj_thr: float,
    lfc_thr: float,
    *,
    min_replicates_per_condition: int = 3,
    covariates: list[str] | None = None,
    lfc_shrink: bool = True,
    n_cpus: int | None = None,
    expose_convergence: bool = True,
) -> tuple[dict[str, Any], list[str]]:
    """Run PyDESeq2 with pre-computed external (spike-in) size factors."""
    from aria.scripts.rna_bulk.deseq2 import (
        _build_design_formula,
        _resolve_covariates,
        _serialize_fit_warnings,
        _shrink_coeff,
    )
    from aria.utils.design_matrix import validate_design_matrix
    from aria.utils.power_estimation import bulk_power_estimate

    warnings_out: list[str] = []

    mask = metadata[design_factor].isin([numerator, denominator])
    meta_sub = metadata[mask].copy()
    counts_sub = counts[meta_sub.index].copy()

    group_sizes = meta_sub[design_factor].value_counts().to_dict()
    n_num = int(group_sizes.get(numerator, 0))
    n_den = int(group_sizes.get(denominator, 0))
    if n_num < min_replicates_per_condition or n_den < min_replicates_per_condition:
        return {
            "status": "error",
            "error_type": "InsufficientReplicates",
            "details": (
                f"DESeq2 requires at least {min_replicates_per_condition} "
                f"replicates per group. Found {group_sizes}."
            ),
        }, warnings_out

    usable_covariates, dropped_covariates = _resolve_covariates(
        meta_sub, design_factor, covariates or []
    )
    for cov, reason in dropped_covariates:
        warnings_out.append(
            f"Confirmed covariate '{cov}' was NOT adjusted for: {reason}. "
            f"The fitted model does not control for it."
        )

    design_check = validate_design_matrix(
        meta_sub,
        condition_col=design_factor,
        covariates=usable_covariates,
        min_replicates_per_condition=min_replicates_per_condition,
    )
    warnings_out.extend(
        f"Design matrix {issue['severity']}: {issue['message']}"
        for issue in design_check.get("issues", [])
    )
    if design_check.get("status") == "blocking":
        return {
            "status": "error",
            "error_type": "InvalidDesignMatrix",
            "details": "; ".join(
                issue["message"]
                for issue in design_check.get("issues", [])
                if issue.get("severity") == "blocking"
            )[:500],
            "design_check": design_check,
        }, warnings_out

    sf_raw = external_size_factors.reindex(meta_sub.index).astype(float)
    if sf_raw.isna().any() or (sf_raw <= 0).any():
        return {
            "status": "error",
            "error_type": "InvalidSpikeInSizeFactors",
            "details": "External spike-in size factors must be finite and > 0 for all contrast samples.",
        }, warnings_out
    # Re-center subset size factors to geometric mean = 1.0 matching R DESeq2 subset semantics
    sf_sub = sf_raw / float(np.exp(np.mean(np.log(sf_raw.to_numpy(dtype=float)))))

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        means = counts_sub.mean(axis=1).astype(float)
        variances = counts_sub.var(axis=1, ddof=1).astype(float)
        disp = ((variances - means) / (means ** 2)).replace(
            [float("inf"), float("-inf")], float("nan")
        ).dropna()
        dispersion_estimate = float(disp[disp > 0].median()) if (disp > 0).any() else 0.1
        mean_expression = float(means[means > 0].median()) if (means > 0).any() else 1.0
        power_estimate = bulk_power_estimate(
            n_per_group=(n_num, n_den),
            mean_expression=mean_expression,
            dispersion=dispersion_estimate,
            target_log2fc=lfc_thr,
            alpha=padj_thr,
        )

        design_formula = _build_design_formula(design_factor, usable_covariates)
        ref_kwargs = {"ref_level": [design_factor, denominator]} if lfc_shrink else {}
        par_kwargs = {"n_cpus": int(n_cpus)} if n_cpus else {}

        dds = DeseqDataSet(
            counts=counts_sub.T,
            metadata=meta_sub,
            design=design_formula,
            refit_cooks=True,
            quiet=True,
            **ref_kwargs,
            **par_kwargs,
        )

        wald_kwargs: dict[str, Any] = {
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

        with w.catch_warnings(record=True) as fit_warns:
            w.simplefilter("always")
            # Inject external heterologous spike-in size factors directly (bypassing
            # fit_size_factors() median-of-ratios so global unidirectional chromatin
            # shifts are preserved).
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

            stat_res = DeseqStats(dds, **wald_kwargs)
            stat_res.summary()
        warnings_out.extend(_serialize_fit_warnings(fit_warns))

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
                except Exception as sx:
                    shrink_reason = f"lfc_shrink failed: {str(sx)[:120]}"
        else:
            shrink_reason = "disabled"

        results_df = stat_res.results_df.dropna(subset=["padj"]).copy()
        results_df["log2FoldChange_raw"] = raw_lfc.reindex(results_df.index)
        if expose_convergence and "_LFC_converged" in dds.var.columns:
            results_df["lfc_converged"] = dds.var["_LFC_converged"].reindex(results_df.index)

        sig = results_df[results_df["padj"] < padj_thr].sort_values("padj")
        up_genes = list(sig[sig["log2FoldChange"] > 0].index)
        down_genes = list(sig[sig["log2FoldChange"] < 0].index)

        low_power = n_num <= 2 or n_den <= 2
        low_power_reason = (
            f"n={n_num} ({numerator}) vs n={n_den} ({denominator}): "
            f"dispersion estimation is unreliable with fewer than three "
            f"replicates per group. Interpret with caution."
        ) if low_power else None
        if low_power:
            warnings_out.append(low_power_reason)

        return {
            "status": "success",
            "results": results_df,
            "n_sig": len(sig),
            "n_up": len(up_genes),
            "n_down": len(down_genes),
            "sig_genes": list(sig.index),
            "up_genes": up_genes,
            "down_genes": down_genes,
            "n_replicates": {numerator: n_num, denominator: n_den},
            "low_power_warning": low_power,
            "low_power_reason": low_power_reason,
            "lfc_shrinkage": {"applied": shrink_applied, "reason": shrink_reason},
            "lfc_threshold_test": lfc_threshold_test,
            "fitted_design_formula": design_formula,
            "covariates_adjusted": usable_covariates,
            "covariates_dropped": [
                {"covariate": c, "reason": r} for c, r in dropped_covariates
            ],
            "power_estimate_at_lfc_min": power_estimate,
            "size_factors_used": {
                str(k): round(float(v), 6) for k, v in dds.obs["size_factors"].items()
            },
            "normalization_method": "heterologous_spikein_size_factors",
        }, warnings_out
    except Exception as exc:
        return {
            "status": "error",
            "error_type": type(exc).__name__,
            "details": str(exc)[:500],
        }, warnings_out
