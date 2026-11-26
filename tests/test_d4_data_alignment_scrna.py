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


def test_d4_b905_scrna_helpers_strict_zip(tmp_path: Path):
    """B905 (rna_pseudobulk_de.py:295, rna_celltypist.py:207, rna_de_per_cluster.py:142):
    mismatched parallel arrays must raise ValueError (or flag degraded=True) instead
    of silently truncating."""
    from aria.scripts.rna_celltypist import _summarize_annotation_confidence
    from aria.scripts.rna_pseudobulk_de import _expressed_background

    # 1. _expressed_background: gene_names length != counts.shape[1] -> degraded=True.
    counts = np.ones((5, 3), dtype=np.int32)
    _, degraded = _expressed_background(counts, ["G1", "G2", "G3", "G4"])
    assert degraded is True, "Mismatched gene_names length silently truncated with degraded=False"

    # 2. _summarize_annotation_confidence: mismatched sequence lengths -> ValueError.
    with pytest.raises(ValueError):
        _summarize_annotation_confidence(
            ["0", "0", "1"],
            ["T", "T"],
            ["T", "T", "B"],
            [0.9, 0.8, 0.7],
        )

    # 3. rna_de_per_cluster: mismatched rank_genes_groups arrays -> ValueError.
    ad = pytest.importorskip("anndata")
    from aria.scripts.rna_de_per_cluster import rna_de_per_cluster

    adata = ad.AnnData(
        X=np.ones((6, 3), dtype=np.float32),
        obs=pd.DataFrame({"leiden": ["0", "0", "0", "1", "1", "1"]}, index=[f"c{i}" for i in range(6)]),
        var=pd.DataFrame(index=["G0", "G1", "G2"]),
    )
    adata.uns["rank_genes_groups"] = {
        "params": {"groupby": "leiden", "use_raw": False, "method": "wilcoxon"},
        "names": np.array([("G0", "G1"), ("G1", "G0"), ("G2", "G2")], dtype=[("0", "U4"), ("1", "U4")]),
        "logfoldchanges": np.array([(2.0, -2.0), (1.5, -1.5)], dtype=[("0", "f4"), ("1", "f4")]),
        "pvals_adj": np.array([(0.001, 0.001), (0.01, 0.01), (0.5, 0.5)], dtype=[("0", "f4"), ("1", "f4")]),
    }
    h5 = tmp_path / "rgg_mismatch.h5ad"
    adata.write_h5ad(h5)
    with pytest.raises(ValueError):
        rna_de_per_cluster({"data_path": str(h5), "groupby": "leiden", "output_dir": str(tmp_path / "out")})
