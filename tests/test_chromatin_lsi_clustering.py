"""v4.6 scATAC step 2 — chromatin LSI + clustering guards.

Covers the honesty/structure contract for `aria/scripts/chromatin_lsi_clustering.py`:

- error routing is dependency-light and structured (no fabrication);
- the script is registered with an IPC contract;
- the TF-IDF helper is real linear algebra;
- the full TF-IDF -> SVD/LSI -> depth-drop -> neighbors/UMAP/Leiden pipeline
  forms clusters from a synthetic peak matrix, reports only computed quantities,
  and resumes from a valid prior output.

Synthetic labels only (ADR-011): no biological gene/peak names.
"""

import os
from pathlib import Path

import pytest


# ── Dependency-light error routing ───────────────────────────────────────────

def test_load_atac_adata_file_not_found(tmp_path):
    from aria.scripts.chromatin_lsi_clustering import _load_atac_adata
    res = _load_atac_adata(str(tmp_path / "missing.h5ad"))
    assert isinstance(res, dict)
    assert res["status"] == "error" and res["error_type"] == "FileNotFound"


def test_load_atac_adata_unsupported_suffix(tmp_path):
    from aria.scripts.chromatin_lsi_clustering import _load_atac_adata
    p = tmp_path / "peaks.txt"
    p.write_text("not an anndata")
    res = _load_atac_adata(str(p))
    assert isinstance(res, dict)
    assert res["status"] == "error" and res["error_type"] == "UnsupportedInput"


def test_h5mu_routing_blocks_without_tooling(tmp_path, monkeypatch):
    from aria.scripts.chromatin_lsi_clustering import _load_atac_adata
    from aria.utils import mudata_io
    p = tmp_path / "paired.h5mu"
    p.write_bytes(b"\x00")
    monkeypatch.setattr(mudata_io, "_import_mudata", lambda: None)
    res = _load_atac_adata(str(p))
    assert isinstance(res, dict)
    assert res["status"] == "error" and res["error_type"] == "MissingDependency"


# ── IPC contract registration ────────────────────────────────────────────────

def test_lsi_clustering_has_ipc_contract():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS
    key = "aria/scripts/chromatin_lsi_clustering.py"
    assert key in SCRIPT_CONTRACTS
    contract = SCRIPT_CONTRACTS[key]
    in_names = {f.name for f in contract.inputs}
    assert "data_path" in in_names


# ── TF-IDF helper is real linear algebra ─────────────────────────────────────

def test_tfidf_is_real_and_handles_empty_rows():
    pytest.importorskip("scipy")
    import numpy as np
    import scipy.sparse as sp
    from aria.scripts.chromatin_lsi_clustering import _tfidf

    counts = sp.csr_matrix(np.array([
        [2.0, 0.0, 1.0],
        [0.0, 0.0, 0.0],   # empty cell must not produce NaN/inf
        [1.0, 3.0, 0.0],
    ]))
    tfidf, cell_totals = _tfidf(counts)
    arr = tfidf.toarray()
    assert np.all(np.isfinite(arr))
    assert list(cell_totals) == [3.0, 0.0, 4.0]
    # a peak present in fewer cells gets a higher IDF weight than a ubiquitous one
    assert arr[2, 1] > 0.0


# ── Full pipeline on a synthetic peak matrix ─────────────────────────────────

def _make_synthetic_atac(tmp_path, seed=0):
    """Two cell populations with distinct accessible peak blocks plus per-cell
    depth variation, so an LSI depth axis and separable clusters both exist."""
    ad = pytest.importorskip("anndata")
    import numpy as np

    rng = np.random.default_rng(seed)
    n_per = 40
    n_peaks = 200
    blocks = []
    for grp in range(2):
        base = np.zeros((n_per, n_peaks))
        lo, hi = (0, 100) if grp == 0 else (100, 200)
        base[:, lo:hi] = 1.0
        # per-cell depth multiplier (creates the depth-correlated LSI axis)
        depth = rng.uniform(0.5, 3.0, size=(n_per, 1))
        counts = rng.poisson(base * depth * 2.0)
        blocks.append(counts)
    X = np.vstack(blocks).astype("float32")
    adata = ad.AnnData(X)
    p = tmp_path / "atac_peaks.h5ad"
    adata.write_h5ad(str(p))
    return p


def test_lsi_clustering_pipeline_structure(tmp_path):
    pytest.importorskip("scanpy")
    pytest.importorskip("sklearn")
    pytest.importorskip("scipy")
    pytest.importorskip("igraph")
    from aria.scripts.chromatin_lsi_clustering import chromatin_lsi_clustering

    data_path = _make_synthetic_atac(tmp_path)
    res = chromatin_lsi_clustering({
        "data_path": str(data_path),
        "n_components": 20,
        "resolution": 1.0,
        "seed": 0,
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert res["input_kind"] == "h5ad"
    assert res["n_cells_total"] == 80
    assert res["n_cells_after_doublet_filter"] <= res["n_cells_total"]
    assert res["n_peaks"] == 200
    assert res["doublets"]["ran"] is True
    assert "batch_qc" in res
    assert res["consensus_peaks"]["status"] in {"verified", "partial", "unverified"}
    assert res["rep_used"] == "X_lsi"
    # two planted populations should separate into >= 2 clusters
    assert res["n_clusters"] >= 2
    assert sum(res["cluster_sizes"].values()) == res["n_cells_used"]
    # honest plumbing: components accounting is consistent, no fabricated metrics
    assert res["n_components_used"] == (
        res["n_components_computed"] - len(res["dropped_components"])
    )
    assert len(res["depth_correlations"]) == res["n_components_computed"]
    assert all(0.0 <= r <= 1.0 for r in res["depth_correlations"])
    assert Path(res["output_path"]).exists()


def test_lsi_clustering_resumes_from_valid_output(tmp_path):
    pytest.importorskip("scanpy")
    pytest.importorskip("sklearn")
    pytest.importorskip("scipy")
    pytest.importorskip("igraph")
    from aria.scripts.chromatin_lsi_clustering import chromatin_lsi_clustering

    data_path = _make_synthetic_atac(tmp_path)
    params = {
        "data_path": str(data_path),
        "n_components": 20,
        "resolution": 1.0,
        "seed": 0,
        "output_dir": str(tmp_path / "out"),
    }
    first = chromatin_lsi_clustering(dict(params))
    assert first["status"] == "success"
    assert not first.get("resumed")

    second = chromatin_lsi_clustering(dict(params))
    assert second["status"] == "success"
    assert second.get("resumed") is True
    assert second["n_clusters"] == first["n_clusters"]
