"""Regression tests for D4 silent data misalignment fixes in multiome / integration."""
from __future__ import annotations

from pathlib import Path

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


def test_d4_013_wnn_diagnostics_applies_fin_both_mask():
    """D4-013 (F841): compute_wnn_anti_degradation_diagnostics must apply fin_both_mask
    so RNA and ATAC modality weights are always evaluated on the exact same cells."""
    from scipy.spatial.distance import cdist
    from aria.utils.wnn_diagnostics import compute_wnn_anti_degradation_diagnostics

    rng = np.random.default_rng(42)
    n = 30
    emb = np.vstack([rng.normal(2.0, 0.2, (15, 6)), rng.normal(-2.0, 0.2, (15, 6))])
    D = cdist(emb, emb)
    rows, cols, vals = [], [], []
    for i in range(n):
        for j in np.argsort(D[i])[1:7]:
            rows.append(i)
            cols.append(j)
            vals.append(float(D[i, j]) + 1e-4)
    wnn_graph = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    clusters = ["0"] * 15 + ["1"] * 15

    # Case 1: Disjoint finite masks (0 cells have both finite) -> evaluable=False.
    diag_disjoint = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=emb,
        atac_embedding=emb,
        wnn_distances=wnn_graph,
        wnn_clusters=clusters,
        rna_clusters=clusters,
        atac_clusters=clusters,
        rna_weights=np.array([0.9] * 15 + [np.nan] * 15),
        atac_weights=np.array([np.nan] * 15 + [0.9] * 15),
        k=6,
        seed=0,
    )
    wa_disjoint = diag_disjoint["per_cluster_weight_audit"]
    assert wa_disjoint["evaluable"] is False
    assert wa_disjoint["reason"] == "non_finite_modality_weights"

    # Case 2: Partially overlapping finite masks -> means computed strictly on shared cells (sum = 1.0).
    rna_w = np.array([0.8] * 10 + [0.95] * 5 + [0.3] * 10 + [np.nan] * 5)
    atac_w = np.array([0.2] * 10 + [np.nan] * 5 + [0.7] * 10 + [0.95] * 5)
    diag_partial = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=emb,
        atac_embedding=emb,
        wnn_distances=wnn_graph,
        wnn_clusters=clusters,
        rna_clusters=clusters,
        atac_clusters=clusters,
        rna_weights=rna_w,
        atac_weights=atac_w,
        k=6,
        seed=0,
    )
    wa_partial = diag_partial["per_cluster_weight_audit"]
    assert wa_partial["evaluable"] is True
    assert abs(wa_partial["global_mean_rna_weight"] + wa_partial["global_mean_atac_weight"] - 1.0) < 1e-6
    for cl_info in wa_partial["clusters"].values():
        assert abs(cl_info["mean_rna_weight"] + cl_info["mean_atac_weight"] - 1.0) < 1e-6


def test_d4_018_mudata_io_h5_obs_length_mismatch_and_strict_zips(tmp_path: Path):
    """D4-018 / B905: Mismatched HDF5 /obs dataset lengths (values/mask or col/_index)
    must fail closed rather than silently truncating labels."""
    h5py = pytest.importorskip("h5py")
    from aria.utils.mudata_io import _read_h5_obs_series, join_obs_labels_from_h5ad

    # 1. Truncated /obs/<col> vs /obs/_index -> join_obs_labels_from_h5ad returns status='error'.
    h5_path = tmp_path / "ref_trunc.h5ad"
    with h5py.File(h5_path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        dt = h5py.string_dtype(encoding="utf-8")
        obs.create_dataset("_index", data=np.array(["BC1", "BC2", "BC3", "BC4"], dtype=object), dtype=dt)
        obs.create_dataset("cell_type", data=np.array(["T_cell", "B_cell"], dtype=object), dtype=dt)
    res = join_obs_labels_from_h5ad(["BC1", "BC2", "BC3", "BC4"], str(h5_path), label_col="cell_type")
    assert res["status"] == "error", f"Mismatched HDF5 obs length silently truncated: {res}"

    # 2. Nullable group values/mask length mismatch -> _read_h5_obs_series raises ValueError.
    h5_null = tmp_path / "ref_null.h5ad"
    with h5py.File(h5_null, "w") as f:
        obs = f.create_group("obs")
        grp = obs.create_group("cell_type")
        dt = h5py.string_dtype(encoding="utf-8")
        grp.create_dataset("values", data=np.array(["T_cell", "B_cell", "Mono"], dtype=object), dtype=dt)
        grp.create_dataset("mask", data=np.array([False, False], dtype=bool))
    with h5py.File(h5_null, "r") as f:
        with pytest.raises(ValueError):
            _read_h5_obs_series(f["obs"], "cell_type")

    # 3. Sample-scoped join with truncated sample_col -> status='error'.
    h5_sc = tmp_path / "ref_scoped.h5ad"
    with h5py.File(h5_sc, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        dt = h5py.string_dtype(encoding="utf-8")
        obs.create_dataset("_index", data=np.array(["BC1", "BC2", "BC3"], dtype=object), dtype=dt)
        obs.create_dataset("cell_type", data=np.array(["T", "B", "M"], dtype=object), dtype=dt)
        obs.create_dataset("sample", data=np.array(["S1"], dtype=object), dtype=dt)
    res_sc = join_obs_labels_from_h5ad(
        ["BC1", "BC2"], str(h5_sc), label_col="cell_type",
        sample_id="S1", reference_sample_col="sample",
    )
    assert res_sc["status"] == "error"


def test_d4_b905_wnn_count_discordant_clusters_strict_zip():
    """B905 (integration_wnn.py:373): _count_discordant_clusters must raise ValueError
    on mismatched wnn_labels / rna_labels lengths."""
    from aria.scripts.integration_wnn import _count_discordant_clusters

    with pytest.raises(ValueError):
        _count_discordant_clusters(["0", "0", "1"], ["0", "0"])
