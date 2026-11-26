"""Regression tests for D4 silent data misalignment fixes in multiome / integration."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def test_d4_002_select_raw_counts_reindexes_permuted_var_and_obs_and_rejects_mismatch():
    """D4-002: _select_raw_counts_matrix must reindex raw.X when obs/var are permutations
    (and record it), accept exact matches directly, and reject non-matching indices."""
    ad = pytest.importorskip("anndata")
    from aria.utils.mudata_io import _select_raw_counts_matrix

    raw_mat = np.array(
        [[500, 0, 10], [600, 1, 15], [0, 800, 20], [1, 750, 12]], dtype=np.int32
    )
    # 1. Reordered var_names (permutation): must reindex raw.X and record it.
    adata = ad.AnnData(
        X=sparse.csr_matrix(raw_mat),
        obs=pd.DataFrame(index=["cell_1", "cell_2", "cell_3", "cell_4"]),
        var=pd.DataFrame(index=["TP53", "ACTB", "GAPDH"]),
    )
    adata.raw = adata.copy()
    adata = adata[:, ["ACTB", "GAPDH", "TP53"]].copy()
    adata.X = sparse.csr_matrix(np.log1p(adata.X.toarray() * 1.37))

    sel = _select_raw_counts_matrix(adata, modality_label="rna", require_raw_counts=True)
    assert sel["ok"] and sel["source"] == "raw.X"
    selected_col0 = np.asarray(sel["matrix"].toarray())[:, 0]
    assert np.array_equal(selected_col0, np.array([0, 1, 800, 750])), (
        f"Column 0 paired with var_names[0]='{adata.var_names[0]}' contains "
        f"{selected_col0.tolist()} (from '{adata.raw.var_names[0]}')"
    )
    assert sel.get("raw_reindexed") is True
    assert sel["inspections"]["raw.X"].get("var_reordered") is True
    assert sel["inspections"]["raw.X"].get("obs_reordered") is False

    # 2. Both obs_names and var_names permuted: must reindex both axes accurately.
    adata_both = ad.AnnData(
        X=sparse.csr_matrix(raw_mat),
        obs=pd.DataFrame(index=["cell_1", "cell_2", "cell_3", "cell_4"]),
        var=pd.DataFrame(index=["TP53", "ACTB", "GAPDH"]),
    )
    adata_both.raw = adata_both.copy()
    adata_both = adata_both[["cell_3", "cell_1", "cell_4", "cell_2"], ["GAPDH", "TP53", "ACTB"]].copy()
    adata_both.X = sparse.csr_matrix(np.log1p(adata_both.X.toarray() * 1.37))
    sel_both = _select_raw_counts_matrix(adata_both, modality_label="rna", require_raw_counts=True)
    assert sel_both["ok"] and sel_both["source"] == "raw.X"
    assert sel_both.get("raw_reindexed") is True
    expected_both = np.array(
        [[20, 0, 800], [10, 500, 0], [12, 1, 750], [15, 600, 1]], dtype=np.int32
    )
    assert np.array_equal(sel_both["matrix"].toarray(), expected_both)

    # 3. Same shape, different var_names (not a permutation): must reject raw.X with recorded reason.
    adata_bad = ad.AnnData(
        X=sparse.csr_matrix(raw_mat),
        obs=pd.DataFrame(index=["cell_1", "cell_2", "cell_3", "cell_4"]),
        var=pd.DataFrame(index=["TP53", "ACTB", "GAPDH"]),
    )
    adata_bad.raw = adata_bad.copy()
    adata_bad.var_names = pd.Index(["TP53", "ACTB", "OTHER_GENE"])
    adata_bad.X = sparse.csr_matrix(np.log1p(adata_bad.X.toarray() * 1.37))
    sel_bad = _select_raw_counts_matrix(adata_bad, modality_label="rna", require_raw_counts=True)
    assert not sel_bad["ok"]
    assert sel_bad["inspections"]["raw.X"]["reason"] == "raw_var_names_mismatch"
    assert "raw_var_names_mismatch" in sel_bad["details"]
