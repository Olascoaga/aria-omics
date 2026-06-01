"""P1-1 (a): unify bulk DE rigor with pseudobulk — apeGLM LFC shrinkage.

Pseudobulk already reports apeGLM-shrunken log2 fold changes (ADR-023) and gates
the effect-size threshold on the shrunken estimate while keeping the raw MLE.
Bulk DE used the raw MLE log2FC directly. This brings bulk to the same rigor:
the reported log2FoldChange is the apeGLM-shrunken value, log2FoldChange_raw
preserves the MLE, and p-values/padj are unchanged (apeGLM does not touch them).
"""

import numpy as np
import pytest


def _bulk_inputs(seed=3):
    import pandas as pd
    rng = np.random.default_rng(seed)
    samples = [f"ctrl{i}" for i in range(4)] + [f"treat{i}" for i in range(4)]
    meta = pd.DataFrame(
        {"condition": ["ctrl"] * 4 + ["treat"] * 4}, index=samples)
    n_genes = 60
    counts = rng.poisson(120, size=(n_genes, 8)).astype(float)
    counts[:8, 4:] *= rng.uniform(3.0, 5.0, size=(8, 4))   # strong DE, first 8
    df = pd.DataFrame(counts.round(),
                      index=[f"g{i}" for i in range(n_genes)], columns=samples)
    return df, meta


def test_bulk_applies_apeglm_and_preserves_raw():
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk_de import _run_deseq2

    counts, meta = _bulk_inputs()
    res, _w = _run_deseq2(counts, meta, "condition", "treat", "ctrl",
                          padj_thr=0.1, lfc_thr=0.5)
    assert res["status"] == "success", res
    assert res["lfc_shrinkage"]["applied"] is True
    assert res["lfc_shrinkage"]["method"] == "apeGLM"

    rdf = res["results"]
    assert "log2FoldChange_raw" in rdf.columns
    # apeGLM shrinks toward zero: |shrunk| <= |raw| (within fp tolerance).
    both = rdf.dropna(subset=["log2FoldChange", "log2FoldChange_raw"])
    assert (both["log2FoldChange"].abs()
            <= both["log2FoldChange_raw"].abs() + 1e-6).all()


def test_bulk_shrinkage_can_be_disabled():
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk_de import _run_deseq2

    counts, meta = _bulk_inputs()
    res, _w = _run_deseq2(counts, meta, "condition", "treat", "ctrl",
                          padj_thr=0.1, lfc_thr=0.5, lfc_shrink=False)
    assert res["status"] == "success", res
    assert res["lfc_shrinkage"]["applied"] is False
    rdf = res["results"]
    # With shrinkage off, the reported LFC equals the raw MLE.
    both = rdf.dropna(subset=["log2FoldChange", "log2FoldChange_raw"])
    assert np.allclose(both["log2FoldChange"], both["log2FoldChange_raw"])


def test_bulk_lfc_threshold_is_inside_wald_test():
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk_de import _run_deseq2

    counts, meta = _bulk_inputs()
    res0, _w0 = _run_deseq2(counts, meta, "condition", "treat", "ctrl",
                            padj_thr=0.1, lfc_thr=0.0)
    res1, _w1 = _run_deseq2(counts, meta, "condition", "treat", "ctrl",
                            padj_thr=0.1, lfc_thr=1.0)
    assert res0["status"] == "success", res0
    assert res1["status"] == "success", res1
    assert res0["lfc_threshold_test"]["applied"] is False
    assert res1["lfc_threshold_test"]["applied"] is True
    assert res1["lfc_threshold_test"]["lfc_null"] == 1.0
    assert res1["lfc_threshold_test"]["alt_hypothesis"] == "greaterAbs"

    # P1-1b: changing lfc_thr changes the Wald null and therefore p-values/padj;
    # significance is then simply padj-thresholded, with no second |LFC| gate.
    common = res0["results"].index.intersection(res1["results"].index)
    padj0 = res0["results"].loc[common, "padj"].dropna()
    padj1 = res1["results"].loc[padj0.index, "padj"].dropna()
    assert not np.allclose(padj0.loc[padj1.index], padj1)
    expected = set(res1["results"].index[res1["results"]["padj"] < 0.1])
    assert set(res1["sig_genes"]) == expected
