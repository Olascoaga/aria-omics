"""Stage 4 C4: apeGLM LFC shrinkage in pseudobulk DE.

Raw MLE log2 fold changes overestimate effect sizes for low-count / noisy genes.
pydeseq2's apeGLM lfc_shrink shrinks them toward zero, leaving p-values
unchanged. ARIA reports the shrunken estimate (and gates effect size on it) while
keeping the unshrunken MLE as log2fc_raw.
"""

import os
import tempfile

import numpy as np
import pytest


def _make_h5ad(path, seed=7):
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(seed)
    rows, obs = [], []
    for cond in ("ctrl", "treat"):
        for r in range(4):
            for _ in range(60):
                base = rng.poisson(80, 40).astype(float)
                if cond == "treat":
                    base[:6] *= rng.uniform(2.5, 4.0)   # strong DE, first 6 genes
                rows.append(base)
                obs.append((cond, f"{cond}_{r}", "A"))
    X = np.vstack(rows).round()
    O = pd.DataFrame(obs, columns=["condition", "replicate", "groupby"])
    A = ad.AnnData(X=X, obs=O)
    A.var_names = [f"g{i}" for i in range(40)]
    A.write_h5ad(path)


def _run(tmp, **overrides):
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de
    p = os.path.join(tmp, "pb.h5ad")
    _make_h5ad(p)
    params = {
        "data_path": p, "groupby": "groupby", "condition_col": "condition",
        "replicate_col": "replicate", "comparisons": [["treat", "ctrl"]],
        "min_replicates_per_condition": 3, "min_cells_per_pseudosample": 10,
    }
    params.update(overrides)
    return rna_pseudobulk_de(params)


def test_shrinkage_applied_and_raw_preserved():
    pytest.importorskip("pydeseq2")
    with tempfile.TemporaryDirectory() as tmp:
        res = _run(tmp)
    assert res["status"] == "success"
    assert res["lfc_shrinkage"]["requested"] is True
    blk = res["per_group"]["A"]["per_comparison"]["treat_vs_ctrl"]
    assert blk["lfc_shrinkage"]["applied"] is True
    assert blk["lfc_shrinkage"]["method"] == "apeGLM"
    sig = blk["all_sig"]
    assert sig, "expected significant DE genes"
    for g in sig:
        assert "log2fc_raw" in g and g["log2fc_raw"] is not None
        # apeGLM never increases magnitude: |shrunk| <= |raw| (+ rounding slack)
        assert abs(g["log2fc"]) <= abs(g["log2fc_raw"]) + 1e-3


def test_shrinkage_can_be_disabled():
    pytest.importorskip("pydeseq2")
    with tempfile.TemporaryDirectory() as tmp:
        res = _run(tmp, lfc_shrink=False)
    blk = res["per_group"]["A"]["per_comparison"]["treat_vs_ctrl"]
    assert blk["lfc_shrinkage"]["applied"] is False
    assert blk["lfc_shrinkage"]["reason"] == "disabled"
    # with shrinkage off, reported log2fc equals the raw MLE
    for g in blk["all_sig"]:
        assert abs(g["log2fc"] - g["log2fc_raw"]) < 1e-3


def test_shrinkage_clause_in_narrative():
    from aria.agents._narrative_scrna import _lfc_shrinkage_clause
    assert "apeGLM" in _lfc_shrinkage_clause({"lfc_shrinkage": {"requested": True}})
    assert _lfc_shrinkage_clause({"lfc_shrinkage": {"requested": False}}) == ""
    assert _lfc_shrinkage_clause({}) == ""
