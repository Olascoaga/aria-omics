"""D8 convergence regression tests for DESeq2 spike-in and pseudobulk helpers (D8-001, D2-BUG-6, CL-030)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_d8_001_deseq2_spikein_vs_bulk_reference_anchoring_and_covariate():
    """D8-001 / D2-BUG-6: spike-in DESeq2 must anchor Treatment(reference=denominator)
    for alphabetical contrasts (KO vs WT) and coerce numeric covariates like _run_deseq2."""
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk.deseq2 import _run_deseq2
    from aria.utils.epigenomic_spikein_de import run_deseq2_with_external_size_factors

    rng = np.random.default_rng(42)
    n_genes = 60
    base = rng.poisson(200, size=(n_genes, 6))
    base[:15, :3] = rng.poisson(900, size=(15, 3))
    base[15:30, 3:] = rng.poisson(900, size=(15, 3))
    samples = ["KO_1", "KO_2", "KO_3", "WT_1", "WT_2", "WT_3"]
    counts = pd.DataFrame(base, index=[f"PEAK_{i:03d}" for i in range(n_genes)], columns=samples)
    meta = pd.DataFrame(
        {"condition": ["KO"] * 3 + ["WT"] * 3, "age": ["20", "40", "60", "25", "45", "65"]},
        index=samples,
    )
    unit_sf = pd.Series(1.0, index=samples)

    r_bulk, _ = _run_deseq2(
        counts, meta[["condition"]], "condition", "KO", "WT", 0.05, 0.5,
        min_replicates_per_condition=3, lfc_shrink=True,
    )
    r_spike, _ = run_deseq2_with_external_size_factors(
        counts, meta[["condition"]], unit_sf, "condition", "KO", "WT", 0.05, 0.5,
        min_replicates_per_condition=3, lfc_shrink=True,
    )
    r_bulk_cov, _ = _run_deseq2(
        counts, meta, "condition", "KO", "WT", 0.05, 0.5,
        covariates=["age"], min_replicates_per_condition=3, lfc_shrink=True,
    )
    r_spike_cov, _ = run_deseq2_with_external_size_factors(
        counts, meta, unit_sf, "condition", "KO", "WT", 0.05, 0.5,
        covariates=["age"], min_replicates_per_condition=3, lfc_shrink=True,
    )

    assert r_bulk["status"] == "success" and r_bulk["lfc_shrinkage"]["applied"] is True
    assert r_bulk_cov["status"] == "success" and r_bulk_cov["lfc_shrinkage"]["applied"] is True
    assert r_spike["status"] == "success" and r_spike["lfc_shrinkage"]["applied"] is True
    assert r_spike_cov["status"] == "success" and r_spike_cov["lfc_shrinkage"]["applied"] is True
    assert r_spike["normalization_method"] == "heterologous_spikein_size_factors"


def test_cl_030_pseudobulk_shares_canonical_deseq2_formula_and_shrink_helpers():
    """CL-030: rna_pseudobulk_de must import formula/shrink helpers directly from rna_bulk.deseq2."""
    from aria.scripts import rna_pseudobulk_de as pb
    from aria.scripts.rna_bulk import deseq2 as bulk

    assert pb._quote_formula_term is bulk._quote_formula_term
    assert pb._can_wrap_in_formulaic_c is bulk._can_wrap_in_formulaic_c
    assert pb._build_design_formula is bulk._build_design_formula
    assert pb._shrink_coeff is bulk._shrink_coeff
    assert pb._detect_shrink_sign_flip is bulk._detect_shrink_sign_flip
