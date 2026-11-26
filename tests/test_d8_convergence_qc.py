"""D8 convergence regression tests for QC helpers (D8-002, D8-020)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def test_d8_002_mito_qc_ensembl_var_names_rna_qc_vs_integration_paired_qc(tmp_path: Path):
    """D8-002: rna_qc must detect MT genes in var symbol columns identically to integration_paired_qc."""
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    pytest.importorskip("scanpy")
    from scipy import sparse

    from aria.scripts.integration_paired_qc import integration_paired_qc
    from aria.scripts.rna_qc import rna_qc

    rng = np.random.default_rng(7)
    n_healthy, n_dead = 80, 20
    n_cells, n_nuc, n_mt = n_healthy + n_dead, 250, 5
    X = np.zeros((n_cells, n_nuc + n_mt), dtype=np.int32)
    X[:n_healthy, :n_nuc] = rng.poisson(4, size=(n_healthy, n_nuc)) + 1
    X[:n_healthy, n_nuc:] = rng.poisson(4, size=(n_healthy, n_mt))
    X[n_healthy:, :n_nuc] = rng.poisson(3, size=(n_dead, n_nuc)) + 1
    X[n_healthy:, n_nuc:] = rng.poisson(300, size=(n_dead, n_mt))
    barcodes = [f"CELL_{i:04d}-1" for i in range(n_cells)]
    gene_ids = [f"ENSG{100000 + i:011d}" for i in range(n_nuc + n_mt)]
    gene_syms = [f"GENE{i}" for i in range(n_nuc)] + [f"MT-ND{i + 1}" for i in range(n_mt)]
    rna = ad.AnnData(
        X=sparse.csr_matrix(X),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame({"gene_symbols": gene_syms}, index=gene_ids),
    )
    h5ad_path = tmp_path / "scrna_ensembl.h5ad"
    rna.write_h5ad(h5ad_path)
    atac = ad.AnnData(
        X=sparse.csr_matrix(rng.poisson(2, size=(n_cells, 250)).astype(np.int32) + 1),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=[f"chr1:{1000 + i * 500}-{1400 + i * 500}" for i in range(250)]),
    )
    h5mu_path = tmp_path / "paired_ensembl.h5mu"
    md.MuData({"rna": rna.copy(), "atac": atac}).write_h5mu(h5mu_path)

    qc_scrna = rna_qc({
        "data_path": str(h5ad_path),
        "output_dir": str(tmp_path / "scrna_out"),
        "organism": "Homo sapiens",
        "run_scrublet": False,
        "mt_threshold": 20.0,
        "min_genes": 50,
    })
    qc_paired = integration_paired_qc({
        "h5mu_path": str(h5mu_path),
        "output_dir": str(tmp_path / "paired_out"),
        "max_pct_mt": 20.0,
        "min_genes": 50,
        "min_rna_counts": 10,
        "min_peaks": 50,
        "min_atac_counts": 5,
        "resume": False,
    })
    assert qc_scrna["status"] == "success" and qc_paired["status"] == "success"
    paired_mt = qc_paired["rna_qc_summary"]["median_pct_mt_before"]
    scrna_mt_max = qc_scrna["mt_stats"]["max"]
    assert paired_mt is not None and paired_mt > 1.0 and qc_paired["n_cells_after"] == 80
    assert qc_scrna["n_cells_after"] == qc_paired["n_cells_after"]
    assert scrna_mt_max is not None and scrna_mt_max > 0.0
    assert qc_scrna["mt_measured"] is True
    assert qc_scrna["mt_gene_source"] == "var['gene_symbols']"
    assert qc_scrna["metrics_not_computed"] == []


def test_d8_002_rna_qc_zero_mt_genes_omits_filtering_and_records_uncomputed(tmp_path: Path):
    """D8-002 (ADR-002): 0 MT genes must not fabricate 0% or filter on mito."""
    ad = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    from scipy import sparse

    from aria.scripts.rna_qc import rna_qc

    rng = np.random.default_rng(11)
    n_cells, n_genes = 60, 260
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.int32) + 1
    rna = ad.AnnData(
        X=sparse.csr_matrix(X),
        obs=pd.DataFrame(index=[f"C{i:04d}-1" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"ENSG{200000 + i:011d}" for i in range(n_genes)]),
    )
    h5ad_path = tmp_path / "no_mt.h5ad"
    rna.write_h5ad(h5ad_path)

    res = rna_qc({
        "data_path": str(h5ad_path),
        "output_dir": str(tmp_path / "out"),
        "organism": "Homo sapiens",
        "run_scrublet": False,
        "mt_threshold": 0.0,  # would drop everything if fabricated pct_counts_mt were compared
        "min_genes": 50,
    })
    assert res["status"] == "success", res
    assert res["n_cells_after"] == n_cells
    assert res["mt_measured"] is False
    assert res["mito_fraction"] is None
    assert res["mt_threshold_used"] is None
    assert res["mt_stats"]["mean"] is None
    assert res["mt_stats"]["max"] is None
    assert res["mito_qc_degraded"] is True
    assert res.get("mito_qc_reason")
    assert any(m.startswith("mito_fraction") for m in res["metrics_not_computed"])

    filtered = ad.read_h5ad(res["output_path"])
    assert "pct_counts_mt" not in filtered.obs.columns


def test_d8_002_uncomputed_mito_surfaces_as_structured_signal_and_truthful_methods():
    """D8-002: uncomputed mito_fraction surfaces via warning_router and avoids false Methods claims."""
    from aria.utils.ambient_qc import detect_mt_mask

    class _DummyAnnData:
        def __init__(self, var_names, var_df):
            self.var_names = var_names
            self.var = var_df

    dummy = _DummyAnnData(
        ["ENSG00001", "ENSG00002"],
        pd.DataFrame({"gene_short_name": ["ACTB", "MT-CO1"]}),
    )
    mask, source = detect_mt_mask(dummy, preferred_prefix="MT-")
    assert mask.tolist() == [False, True]
    assert source == "var['gene_short_name']"

    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.registry import registry_with
    from aria.agents.narrative.scrna.text import build_scrna_methods

    qc_payload = {
        "status": "success",
        "n_cells_before": 100,
        "n_cells_after": 95,
        "pct_removed": 5.0,
        "mt_measured": False,
        "mt_threshold_used": None,
        "mito_fraction": None,
        "mito_qc_degraded": True,
        "mito_qc_reason": "no MT-/mt- mitochondrial genes matched in var_names or gene symbol columns",
        "metrics_not_computed": [
            "mito_fraction (no MT-/mt- mitochondrial genes matched in var_names or gene symbol columns)"
        ],
        "warnings": ["No mitochondrial genes matched; mito_fraction not computed."],
    }
    reg = registry_with([ScrnaNarrator()])
    blocks = reg.collect_blocks({"scrna_agent": {"findings": {"qc": qc_payload}}})
    qc_block = next(b for b in blocks if b.id == "scrna.qc")
    assert qc_block.confidence == "medium"
    sig_fields = {s["field"] for s in qc_block.metadata.get("degradation_signals", [])}
    assert "mito_qc_degraded" in sig_fields
    assert "metrics_not_computed" in sig_fields
    assert any("mito_fraction" in c.text for c in qc_block.caveats)

    methods = build_scrna_methods({"qc": qc_payload})
    assert "percent.mt" not in methods
    assert "mitochondrial filtering was skipped" in methods
