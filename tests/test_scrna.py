"""ARIA scRNA-seq — native pytest (P1-11 follow-up).

Converted from the legacy script-style harness (top-level ok()/fail()/sys.exit,
which crashed pytest collection and was never run in CI). Source-text contract
checks run in every lane; the scanpy marker/DE checks are gated on scanpy (heavy
lane), the agent-import check on litellm, and the live PBMC check on a local
dataset — each lane runs the subset it supports and skips the rest.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).parent.parent
AGENT_FILE = REPO / "aria" / "agents" / "scrna_agent.py"


def _agent_src() -> str:
    parts = [AGENT_FILE.read_text(encoding="utf-8")]
    scrna_pkg = REPO / "aria" / "agents" / "scrna"
    if scrna_pkg.is_dir():
        for py_file in sorted(scrna_pkg.glob("*.py")):
            parts.append(py_file.read_text(encoding="utf-8"))
    return "\n".join(parts)


def make_mock_adata(n_cells=600, n_genes=500):
    """Three-population mock with scanpy rank_genes_groups precomputed."""
    import scanpy as sc
    import anndata as ad
    import pandas as pd
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(42)
    data = rng.negative_binomial(5, 0.5, (n_cells, n_genes)).astype(float)
    per = n_cells // 3
    data[:per, 1:11] *= 15
    data[per:2 * per, 101:111] *= 15
    data[2 * per:, 201:211] *= 15

    adata = ad.AnnData(
        X=csr_matrix(data),
        obs=pd.DataFrame(index=[f"C{i:04d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(n_genes)]),
    )
    adata.obs["leiden"] = pd.Categorical(
        ["0"] * per + ["1"] * per + ["2"] * (n_cells - 2 * per)
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=200, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=5, n_pcs=20)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                            key_added="rank_genes_groups")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                            key_added="de_wilcoxon")
    return adata


# ── agent import (needs litellm via aria.llm.provider) ───────────────────────

def test_scrna_agent_importable():
    pytest.importorskip("litellm")
    from aria.agents.scrna_agent import scRNAAgent  # noqa: F401


# ── source-text contract checks (no heavy import — run in every lane) ────────

def test_agent_delegates_qc_to_env_manager():
    # QC is delegated to aria/scripts/rna_qc.py via the EnvironmentManager
    # subprocess (no silent in-process QC fallback, by design).
    src = _agent_src()
    assert "self.env.run_in_stack" in src or "env_manager" in src
    assert "rna_qc.py" in src


def test_agent_has_doublet_detection():
    src = _agent_src().lower()
    assert "scrublet" in src or "doublet" in src


# ── scanpy-gated: ARIA rna_qc -> rna_clustering -> rna_de_per_cluster oracles ─

def _write_planted_scrna_h5ad(tmp_path: Path, n_per: int = 80, n_genes: int = 360):
    import anndata as ad
    import pandas as pd
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(42)
    n_cells = 3 * n_per
    X = rng.poisson(2.5, size=(n_cells, n_genes)).astype(np.float32)
    planted = {
        0: [f"GENE_{j:03d}" for j in range(1, 13)],
        1: [f"GENE_{j:03d}" for j in range(101, 113)],
        2: [f"GENE_{j:03d}" for j in range(201, 213)],
    }
    X[:n_per, 1:13] += rng.poisson(35.0, size=(n_per, 12))
    X[n_per:2 * n_per, 101:113] += rng.poisson(35.0, size=(n_per, 12))
    X[2 * n_per:, 201:213] += rng.poisson(35.0, size=(n_per, 12))

    # Add 4 damaged high-MT cells to verify rna_qc filtering
    mt_genes = [f"MT-G{j}" for j in range(5)]
    bad_X = np.zeros((4, n_genes + 5), dtype=np.float32)
    bad_X[:, :80] = 1.0
    bad_X[:, n_genes:] = 500.0
    full_X = np.vstack([np.hstack([X, np.zeros((n_cells, 5), dtype=np.float32)]), bad_X])
    adata = ad.AnnData(
        X=csr_matrix(full_X),
        obs=pd.DataFrame(index=[f"C{i:04d}" for i in range(full_X.shape[0])]),
        var=pd.DataFrame(index=[f"GENE_{j:03d}" for j in range(n_genes)] + mt_genes),
    )
    raw_path = tmp_path / "planted_scrna.h5ad"
    adata.write_h5ad(raw_path)
    return raw_path, planted, n_cells


def test_markers_distinct_per_cluster(tmp_path):
    """D3-013: invoke ARIA's real rna_qc + rna_clustering scripts and assert
    exact 3-cluster recovery and distinct planted top markers per cluster."""
    pytest.importorskip("scanpy")
    from aria.scripts.rna_clustering import rna_clustering
    from aria.scripts.rna_qc import rna_qc

    raw_path, planted, n_healthy = _write_planted_scrna_h5ad(tmp_path)
    qc_res = rna_qc({
        "data_path": str(raw_path),
        "organism": "Homo sapiens",
        "min_genes": 50,
        "min_cells": 3,
        "mt_threshold": 10.0,
        "run_scrublet": False,
        "output_dir": str(tmp_path / "qc"),
    })
    assert qc_res["status"] == "success", qc_res
    assert qc_res["n_cells_before"] == n_healthy + 4
    assert qc_res["n_cells_after"] == n_healthy

    cl_res = rna_clustering({
        "data_path": qc_res["output_path"],
        "resolution": 0.4,
        "n_hvg": 200,
        "n_pcs": 15,
        "n_neighbors": 12,
        "seed": 0,
        "output_dir": str(tmp_path / "clust"),
    })
    assert cl_res["status"] == "success", cl_res
    assert cl_res["n_clusters"] == 3, cl_res
    tops = {k: v[0] for k, v in cl_res["top_markers"].items() if v}
    assert len(set(tops.values())) == 3, f"expected 3 distinct cluster top markers: {tops}"
    matched = set()
    for cl, genes in cl_res["top_markers"].items():
        top12 = set(genes[:12])
        for pop_idx, exp_genes in planted.items():
            if len(top12 & set(exp_genes)) >= 8:
                matched.add(pop_idx)
    assert matched == {0, 1, 2}, cl_res["top_markers"]


def test_de_per_cluster_extraction_runs(tmp_path):
    """D3-013: invoke ARIA's real rna_de_per_cluster on clustered data and verify
    planted marker gene recovery with positive log2fc and controlled padj."""
    pytest.importorskip("scanpy")
    from aria.scripts.rna_clustering import rna_clustering
    from aria.scripts.rna_de_per_cluster import rna_de_per_cluster
    from aria.scripts.rna_qc import rna_qc

    raw_path, planted, _ = _write_planted_scrna_h5ad(tmp_path)
    qc_res = rna_qc({
        "data_path": str(raw_path),
        "organism": "Homo sapiens",
        "min_genes": 50,
        "min_cells": 3,
        "mt_threshold": 10.0,
        "run_scrublet": False,
        "output_dir": str(tmp_path / "qc"),
    })
    cl_res = rna_clustering({
        "data_path": qc_res["output_path"],
        "resolution": 0.4,
        "n_hvg": 200,
        "n_pcs": 15,
        "n_neighbors": 12,
        "seed": 0,
        "output_dir": str(tmp_path / "clust"),
    })
    de_res = rna_de_per_cluster({
        "data_path": cl_res["output_path"],
        "groupby": "leiden",
        "padj_max": 0.05,
        "lfc_min": 0.5,
        "top_n": 20,
        "output_dir": str(tmp_path / "de"),
    })
    assert de_res["status"] == "success", de_res
    de_genes = de_res["de_genes_by_cluster"]
    assert len(de_genes) == 3
    assert all(isinstance(v, list) for v in de_genes.values())
    matched = set()
    for cl, recs in de_genes.items():
        pos_hits = {r["gene"] for r in recs if r["log2fc"] > 1.0 and r["padj"] < 0.05}
        for pop_idx, exp_genes in planted.items():
            if len(pos_hits & set(exp_genes)) >= 8:
                matched.add(pop_idx)
    assert matched == {0, 1, 2}, de_genes


# ── live PBMC 3k (scanpy + local dataset) ────────────────────────────────────

def _find_pbmc():
    for p in (Path.home() / "aria-data" / "pbmc3k_test",
              Path.home() / "aria-data" / "pbmc3k_test" / "hg19"):
        if p.exists() and list(p.rglob("*.mtx*")):
            return p
    return None


def test_live_pbmc3k_markers():
    sc = pytest.importorskip("scanpy")
    pbmc = _find_pbmc()
    if pbmc is None:
        pytest.skip("PBMC 3k dataset not present (run install.sh to download)")
    adata = sc.read_10x_mtx(str(pbmc), var_names="gene_symbols", cache=False)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 5].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=40)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2,
                 directed=False)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    rgg = adata.uns["rank_genes_groups"]
    tops = {cl: rgg["names"][cl][0] for cl in rgg["names"].dtype.names}
    assert len(set(tops.values())) > 1, f"marker bug on PBMC: {tops}"
