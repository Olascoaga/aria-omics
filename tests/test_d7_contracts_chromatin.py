"""Regression tests for D7 layer-contract fixes in chromatin & scATAC."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def test_d7_006_cross_modal_triangulation_reads_scatac_da_peaks_by_cluster():
    """D7-006: chromatin_diffacc._per_cluster_accessibility writes per-cluster DA
    peaks under 'da_peaks_by_cluster'. cross_modal_triangulation must read that
    canonical key and emit atac_da_association edges."""
    from aria.agents.narrative.synthesis.cross_modal_triangulation import (
        build_cross_modal_evidence_graph,
    )
    from aria.scripts.chromatin_diffacc import _per_cluster_accessibility

    rng = np.random.default_rng(5)
    n_cells = 60
    X = rng.uniform(0.0, 0.2, size=(n_cells, 10)).astype(np.float32)
    X[:30, 0] += 5.0
    X[:30, 1] += 4.5
    X[30:, 2] += 5.0
    peaks = [f"chr1:{1000 + i * 1000}-{1500 + i * 1000}" for i in range(10)]
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"leiden": ["0"] * 30 + ["1"] * 30}, index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=peaks),
    )
    pc_res = _per_cluster_accessibility(adata, groupby="leiden", padj_max=0.05, lfc_min=0.5, top_n=10, warnings=[])
    assert pc_res["ran"] is True
    assert pc_res["n_da_total"] >= 2
    assert any(pc_res["da_peaks_by_cluster"].values())

    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "scATAC": {
                    "findings": {
                        "differential_accessibility": {
                            "status": "success",
                            "per_cluster": pc_res,
                        }
                    }
                }
            },
        },
        "integration_agent": {
            "status": "done",
            "findings": {
                "peak_to_gene": {
                    "status": "success",
                    "top_links": [
                        {"peak": peaks[0], "gene": "GENE_A", "correlation": 0.72, "padj": 1e-6}
                    ],
                }
            },
        },
    }
    graph = build_cross_modal_evidence_graph(agent_results, exp_ctx={})
    atac_edges = [e for e in graph.get("edges", []) if e.get("relation") == "atac_da_association"]
    assert len(atac_edges) > 0, (
        f"Expected atac_da_association edges for {peaks[0]}, got 0 (da_peaks_by_cluster ignored)"
    )
