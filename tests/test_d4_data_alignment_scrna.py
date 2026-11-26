"""Regression tests for D4 silent data misalignment fixes in scRNA-seq."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def test_d4_006_pseudobulk_drops_nan_donors_and_rejects_noninvariant_covariates(tmp_path: Path):
    """D4-006: NaN donor/condition/group cells must be dropped BEFORE .astype(str)
    (never turned into a literal 'nan' replicate), and covariates varying within
    a pseudosample must be excluded with a warning rather than copying .iloc[0]."""
    ad = pytest.importorskip("anndata")
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    # Part 1: 2 real treated donors + 25 NaN-donor cells must NOT become 3 treated donors.
    rng = np.random.default_rng(11)
    n_cells = 150
    X = sparse.csr_matrix(rng.poisson(12, size=(n_cells, 40)).astype(np.int32))
    conds = ["treat"] * 75 + ["ctrl"] * 75
    reps = ["D1"] * 25 + ["D2"] * 25 + [np.nan] * 25 + ["C1"] * 25 + ["C2"] * 25 + ["C3"] * 25
    batches = ["B1"] * 12 + ["B2"] * 13 + ["B1"] * 25 + ["B2"] * 25 + ["B1"] * 25 + ["B2"] * 25 + ["B1"] * 25
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"leiden": ["0"] * n_cells, "condition": conds, "donor": reps, "batch": batches},
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=[f"GENE{i}" for i in range(40)]),
    )
    h5 = tmp_path / "pb_nan.h5ad"
    adata.write_h5ad(h5)
    out = rna_pseudobulk_de({
        "data_path": str(h5),
        "groupby": "leiden",
        "condition_col": "condition",
        "replicate_col": "donor",
        "covariates": ["batch"],
        "comparisons": [["treat", "ctrl"]],
        "min_cells_per_pseudosample": 10,
        "min_replicates_per_condition": 3,
        "output_dir": str(tmp_path / "out_nan"),
        "run_pathways": False,
    })
    g0 = out["per_group"]["0"]
    assert g0["n_pseudosamples"] == 5, f"Expected 5 pseudosamples (NaN dropped), got {g0}"
    assert g0.get("status") == "skipped" or (
        g0.get("per_comparison", {}).get("treat_vs_ctrl", {}).get("status") != "success"
    )
    assert any("missing" in w.lower() or "dropped" in w.lower() for w in out.get("warnings", []))

    # Part 2: 3 real donors per condition, but donor D1 spans two batches (B1 + B2).
    # 'batch' is non-invariant within pseudosample D1 and must be excluded with a warning,
    # never silently assigned from .iloc[0].
    n_cells2 = 150
    X2 = sparse.csr_matrix(rng.poisson(20, size=(n_cells2, 40)).astype(np.int32))
    reps2 = ["D1"] * 25 + ["D2"] * 25 + ["D3"] * 25 + ["C1"] * 25 + ["C2"] * 25 + ["C3"] * 25
    batches2 = ["B1"] * 12 + ["B2"] * 13 + ["B1"] * 25 + ["B2"] * 25 + ["B1"] * 25 + ["B2"] * 25 + ["B1"] * 25
    adata2 = ad.AnnData(
        X=X2,
        obs=pd.DataFrame(
            {"leiden": ["0"] * n_cells2, "condition": conds, "donor": reps2, "batch": batches2},
            index=[f"c2_{i}" for i in range(n_cells2)],
        ),
        var=pd.DataFrame(index=[f"GENE{i}" for i in range(40)]),
    )
    h5_2 = tmp_path / "pb_cov.h5ad"
    adata2.write_h5ad(h5_2)
    out2 = rna_pseudobulk_de({
        "data_path": str(h5_2),
        "groupby": "leiden",
        "condition_col": "condition",
        "replicate_col": "donor",
        "covariates": ["batch"],
        "comparisons": [["treat", "ctrl"]],
        "min_cells_per_pseudosample": 10,
        "min_replicates_per_condition": 3,
        "output_dir": str(tmp_path / "out_cov"),
        "run_pathways": False,
    })
    comp2 = out2["per_group"]["0"]["per_comparison"]["treat_vs_ctrl"]
    assert "batch" not in (comp2.get("covariates_adjusted") or []), (
        f"Non-invariant within-donor covariate 'batch' was silently fit via .iloc[0]: {comp2}"
    )
    assert any("batch" in w and "within" in w.lower() for w in out2.get("warnings", [])), (
        f"Expected warning about within-pseudosample non-invariant 'batch': {out2.get('warnings')}"
    )
