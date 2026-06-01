"""P1-2: the multiple-testing family is pre-registered, not chosen post-hoc.

Five audits flagged FDR handling; the sharpest, most actionable part is the
anti-cherry-picking guarantee: the per-cluster vs global BH family must be fixed
from the analysis plan BEFORE p-values are seen, and the primary significance
call must depend ONLY on that declared strategy — never on which family yields
more significant genes.

(IHW and s-values are the heavier scientific additions; they are deferred with a
documented rationale — pydeseq2 0.5.4 exposes no s-values, and a faithful IHW
needs the validated Ignatiadis-Huber estimator rather than a hand-rolled
weighted-BH that could silently violate FDR control.)
"""

import os
import tempfile

import numpy as np
import pytest

from aria.utils.stats import (
    assert_fdr_family_not_post_hoc,
    preregister_fdr_family,
    primary_fdr_column,
)


# ── Pure helpers (no pydeseq2) ───────────────────────────────────────────────

def test_preregister_normalizes_and_declares():
    d = preregister_fdr_family("GLOBAL")
    assert d["fdr_strategy"] == "global"
    assert d["preregistered"] is True
    assert d["selected_before_results"] is True
    # Unknown / empty -> conservative per-cluster default.
    assert preregister_fdr_family("nonsense")["fdr_strategy"] == "per_cluster"
    assert preregister_fdr_family(None)["fdr_strategy"] == "per_cluster"


def test_primary_column_depends_only_on_strategy():
    assert primary_fdr_column("global") == "padj_global"
    assert primary_fdr_column("per_cluster") == "padj_local"
    assert primary_fdr_column("whatever") == "padj_local"


def test_post_hoc_switch_is_rejected():
    # Applied column consistent with the declared strategy -> ok.
    assert_fdr_family_not_post_hoc("global", "padj_global")
    assert_fdr_family_not_post_hoc("per_cluster", "padj_local")
    # Applied column that does not match the pre-registration -> integrity error.
    with pytest.raises(ValueError):
        assert_fdr_family_not_post_hoc("per_cluster", "padj_global")


# ── End-to-end: the pseudobulk result records the pre-registration ───────────

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
                    base[:6] *= rng.uniform(2.5, 4.0)
                rows.append(base)
                obs.append((cond, f"{cond}_{r}", "A"))
    X = np.vstack(rows).round()
    O = pd.DataFrame(obs, columns=["condition", "replicate", "groupby"])
    A = ad.AnnData(X=X, obs=O)
    A.var_names = [f"g{i}" for i in range(40)]
    A.write_h5ad(path)


def _run(tmp, strategy):
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de
    p = os.path.join(tmp, "pb.h5ad")
    _make_h5ad(p)
    return rna_pseudobulk_de({
        "data_path": p, "groupby": "groupby", "condition_col": "condition",
        "replicate_col": "replicate", "comparisons": [["treat", "ctrl"]],
        "min_replicates_per_condition": 3, "min_cells_per_pseudosample": 10,
        "fdr_strategy": strategy,
    })


def test_pseudobulk_records_preregistered_fdr_family():
    pytest.importorskip("pydeseq2")
    for strategy, family in (("per_cluster", "padj_local"),
                             ("global", "padj_global")):
        with tempfile.TemporaryDirectory() as tmp:
            res = _run(tmp, strategy)
        assert res["status"] == "success"
        mt = res["multiple_testing"]
        prereg = mt["fdr_preregistration"]
        assert prereg["preregistered"] is True
        assert prereg["fdr_strategy"] == strategy
        assert mt["fdr_strategy"] == strategy
        assert mt["primary_padj_column"] == family
