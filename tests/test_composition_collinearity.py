"""Stage 4 C3: the pseudobulk composition covariate is the cell type's own
log-proportion. When abundance shifts WITH condition it is collinear with the
contrast (variance inflation), so ARIA must drop it for that block and report
why — the shift is still covered by the differential-abundance layer.
"""

import numpy as np
import pytest


def test_abs_corr_detects_collinearity():
    from aria.scripts.rna_pseudobulk_de import _abs_corr
    cond = np.array([0, 0, 0, 1, 1, 1], dtype=float)
    collinear = np.array([0.1, 0.12, 0.09, 0.8, 0.82, 0.79])   # tracks cond
    orthogonal = np.array([0.5, 0.1, 0.9, 0.4, 0.6, 0.2])      # unrelated
    assert _abs_corr(collinear, cond) > 0.9
    assert _abs_corr(orthogonal, cond) < 0.8
    assert _abs_corr(np.ones(6), cond) is None                  # no variance


def _make_collinear_h5ad(path):
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(0)
    rows, obs = [], []
    # Type A abundance tracks condition (collinear); totals kept ~balanced by B.
    plan = {
        ("ctrl", "A"): 15, ("ctrl", "B"): 60,
        ("treat", "A"): 60, ("treat", "B"): 15,
    }
    for cond in ("ctrl", "treat"):
        for r in range(3):
            donor = f"{cond}_{r}"
            for ct in ("A", "B"):
                for _ in range(plan[(cond, ct)]):
                    rows.append(rng.poisson(80, 40))
                    obs.append((cond, donor, ct))
    X = np.vstack(rows).astype(float)
    O = pd.DataFrame(obs, columns=["condition", "replicate", "groupby"])
    A = ad.AnnData(X=X, obs=O)
    A.var_names = [f"g{i}" for i in range(40)]
    A.write_h5ad(path)


def test_composition_covariate_dropped_when_collinear(tmp_path):
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    p = tmp_path / "collinear.h5ad"
    _make_collinear_h5ad(str(p))
    res = rna_pseudobulk_de({
        "data_path": str(p),
        "groupby": "groupby",
        "condition_col": "condition",
        "replicate_col": "replicate",
        "comparisons": [["treat", "ctrl"]],
        "composition_covariate": True,           # request it
        "min_replicates_per_condition": 3,
        "min_cells_per_pseudosample": 10,
    })
    assert res["status"] == "success"
    block_a = res["per_group"]["A"]["per_comparison"]["treat_vs_ctrl"]
    assert block_a["status"] == "success"
    # C3: requested but dropped for collinearity, with an honest reason.
    assert block_a["corrected_for_composition"] is False
    assert "collinear" in str(block_a["composition_skipped_reason"])
