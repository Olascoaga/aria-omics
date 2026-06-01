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
    return AGENT_FILE.read_text()


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


# ── scanpy-gated: per-cluster marker / DE extraction is distinct ─────────────

def test_markers_distinct_per_cluster():
    pytest.importorskip("scanpy")
    adata = make_mock_adata()
    rgg = adata.uns["rank_genes_groups"]
    labels = list(rgg["names"].dtype.names)
    markers = {cl: [g for g in rgg["names"][cl][:20] if g and str(g) != "nan"]
               for cl in labels}
    tops = {k: v[0] for k, v in markers.items() if v}
    # The historical bug returned names[0] for every cluster; distinct tops
    # confirm per-cluster extraction.
    assert len(set(tops.values())) > 1, f"all clusters share top marker: {tops}"


def test_de_per_cluster_extraction_runs():
    pytest.importorskip("scanpy")
    adata = make_mock_adata()
    rgg = adata.uns["de_wilcoxon"]
    de_genes = {}
    for cl in rgg["names"].dtype.names:
        names = rgg["names"][cl][:50]
        pvals = rgg["pvals_adj"][cl][:50]
        lfc = rgg["logfoldchanges"][cl][:50]
        sig = [str(g) for g, p, l in zip(names, pvals, lfc)
               if float(p) < 0.05 and abs(float(l)) > 0.5 and str(g) not in ("nan", "")]
        de_genes[str(cl)] = sig[:20]
    assert len(de_genes) >= 2                      # per-cluster, not collapsed
    assert all(isinstance(v, list) for v in de_genes.values())


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
    adata = sc.read_10x_mtx(str(pbmc), var_names="gene_symbols", cache=True)
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
