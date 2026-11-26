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
        # An all-NaN count column must not shadow a populated scale column.
        if not sample_exposures and scale_col_candidates:
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
        if source_kind in {"explicit_spikein_counts", "explicit_spikein_scale_factors"}:
            return None, {
                "used": False,
                "source": source_kind,
                "explicit_requested": True,
                "invalid_reason": (
                    f"{source_kind} was explicitly supplied but all values were "
                    "None/missing; refusing silent fallback to conventional normalization."
                ),
            }
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
    """Run PyDESeq2 with pre-computed external (spike-in) size factors.

    Delegates to the canonical ``aria.scripts.rna_bulk.deseq2._run_deseq2`` core
    so formulaic ``Treatment(reference=denominator)`` anchoring, continuous /
    categorical covariate coercion, Wald ``lfcThreshold``, and apeGLM sign-flip
    checks never diverge from bulk RNA/ATAC DESeq2 (D8-001 / D2-BUG-6).
    """
    from aria.scripts.rna_bulk.deseq2 import _run_deseq2

    return _run_deseq2(
        counts,
        metadata,
        design_factor,
        numerator,
        denominator,
        padj_thr,
        lfc_thr,
        min_replicates_per_condition=min_replicates_per_condition,
        covariates=covariates,
        lfc_shrink=lfc_shrink,
        n_cpus=n_cpus,
        expose_convergence=expose_convergence,
        external_size_factors=external_size_factors,
    )
