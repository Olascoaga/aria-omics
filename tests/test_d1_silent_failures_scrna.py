"""Permanent regression tests for D1 silent-failure fixes in scRNA pseudobulk DE (D1-005 part 2)."""
from __future__ import annotations

import pytest


def test_d1_005_part2_lfc_shrinkage_clause_none_partial_and_full():
    """D1-005 (part 2 narrator): _lfc_shrinkage_clause and ScrnaNarrator.methods() must be
    truthful for none, partial, and full apeGLM shrinkage application."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.scrna._common import _lfc_shrinkage_clause

    # Not requested -> empty clause
    assert _lfc_shrinkage_clause({"lfc_shrinkage": {"requested": False}}) == ""

    # None applied (all comparisons failed)
    pb_none = {
        "status": "success",
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "lfc_shrinkage": {
            "requested": True,
            "applied": False,
            "n_applied": 0,
            "n_failed": 2,
            "failed_reasons": ["treat_vs_ctrl: lfc_shrink failed: solver error"],
        },
        "per_group": {},
    }
    clause_none = _lfc_shrinkage_clause(pb_none)
    assert "are apeGLM-shrunken estimates" not in clause_none
    assert "unshrunken MLE" in clause_none
    assert "all 2 comparison(s)" in clause_none

    # Partial application (1 of 2 succeeded)
    pb_partial = {
        "status": "success",
        "lfc_shrinkage": {
            "requested": True,
            "applied": False,
            "n_applied": 1,
            "n_failed": 1,
            "failed_reasons": ["treat_vs_ctrl: lfc_shrink failed: solver error"],
        },
        "per_group": {},
    }
    clause_partial = _lfc_shrinkage_clause(pb_partial)
    assert "are apeGLM-shrunken estimates" not in clause_partial
    assert "1 of 2 comparison(s)" in clause_partial
    assert "unshrunken MLE" in clause_partial

    # Full application
    pb_full = {
        "status": "success",
        "lfc_shrinkage": {
            "requested": True,
            "applied": True,
            "n_applied": 2,
            "n_failed": 0,
            "failed_reasons": [],
        },
        "per_group": {},
    }
    clause_full = _lfc_shrinkage_clause(pb_full)
    assert "are apeGLM-shrunken estimates" in clause_full

    # Legacy dict without summary counts falls back to per_group inspection
    pb_legacy_failed = {
        "status": "success",
        "lfc_shrinkage": {"requested": True},
        "per_group": {
            "Tcell": {
                "per_comparison": {
                    "treat_vs_ctrl": {
                        "status": "success",
                        "lfc_shrinkage": {"applied": False, "method": None, "reason": "failed"},
                    }
                }
            }
        },
    }
    clause_legacy = _lfc_shrinkage_clause(pb_legacy_failed)
    assert "are apeGLM-shrunken estimates" not in clause_legacy
    assert "unshrunken MLE" in clause_legacy

    # End-to-end via ScrnaNarrator.methods()
    methods_none = " ".join(ScrnaNarrator().methods("scrna_agent", {"findings": {"pseudobulk_de": pb_none}}))
    assert "are apeGLM-shrunken estimates" not in methods_none
    assert "unshrunken MLE" in methods_none


def test_d1_005_part2_pseudobulk_script_records_none_partial_and_full_shrinkage(tmp_path, monkeypatch):
    """D1-005 (part 2 producer): rna_pseudobulk_de with real pydeseq2 must record realized
    lfc_shrinkage counts/effect truthfully for none (incl. broad IndexError), partial, and full."""
    pytest.importorskip("pydeseq2")
    ad = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd
    from pydeseq2.ds import DeseqStats
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    def _write_pb_h5ad(path, groups):
        rng = np.random.default_rng(42)
        donors = ["d1", "d2", "d3", "d4", "d5", "d6"]
        conds = ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]
        X_rows, obs_rows = [], []
        idx = 0
        for grp in groups:
            for d, c in zip(donors, conds):
                base = rng.poisson(30, size=(10, 20))
                if c == "treat":
                    base[:, :5] = rng.poisson(120, size=(10, 5))
                X_rows.append(base)
                for _ in range(10):
                    obs_rows.append({"cell_type": grp, "condition": c, "donor": d})
                    idx += 1
        ad.AnnData(
            X=np.vstack(X_rows).astype(np.int32),
            obs=pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(idx)]),
            var=pd.DataFrame(index=[f"G{i}" for i in range(20)]),
        ).write_h5ad(path)

    h5_two = tmp_path / "pb_two.h5ad"
    _write_pb_h5ad(h5_two, ["Tcell", "Bcell"])
    base_params = {
        "data_path": str(h5_two),
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["treat", "ctrl"]],
        "min_cells_per_pseudosample": 5,
        "min_replicates_per_condition": 2,
        "lfc_shrink": True,
    }

    # 1. NONE applied (inject broad IndexError to confirm handler stays broad & non-aborting)
    orig_shrink = DeseqStats.lfc_shrink
    monkeypatch.setattr(
        DeseqStats, "lfc_shrink",
        lambda self, coeff=None, **kw: (_ for _ in ()).throw(IndexError("apeGLM index failure")),
    )
    res_none = rna_pseudobulk_de({**base_params, "output_dir": str(tmp_path / "out_none")})
    assert res_none["status"] == "success"
    sk_none = res_none["lfc_shrinkage"]
    assert sk_none["requested"] is True
    assert sk_none["applied"] is False
    assert sk_none["n_applied"] == 0
    assert sk_none["n_failed"] == 2
    assert sk_none["method"] is None
    assert "unshrunken MLE" in sk_none["effect"]
    assert "reported log2fc is the shrunken estimate" not in sk_none["effect"]
    assert any("apeGLM index failure" in r for r in sk_none["failed_reasons"])

    # 2. PARTIAL applied (first group succeeds, second raises)
    calls = {"n": 0}

    def _partial_shrink(self, coeff=None, **kw):
        calls["n"] += 1
        if calls["n"] == 2:
            raise RuntimeError("second group solver failed")
        return orig_shrink(self, coeff=coeff, **kw)

    monkeypatch.setattr(DeseqStats, "lfc_shrink", _partial_shrink)
    res_partial = rna_pseudobulk_de({**base_params, "output_dir": str(tmp_path / "out_partial")})
    assert res_partial["status"] == "success"
    sk_partial = res_partial["lfc_shrinkage"]
    assert sk_partial["requested"] is True
    assert sk_partial["applied"] is False
    assert sk_partial["n_applied"] == 1
    assert sk_partial["n_failed"] == 1
    assert "1/2" in sk_partial["effect"]
    assert "unshrunken MLE" in sk_partial["effect"]
    assert "reported log2fc is the shrunken estimate" not in sk_partial["effect"]
    assert any("second group solver failed" in r for r in sk_partial["failed_reasons"])

    # 3. FULL applied (all succeed)
    monkeypatch.setattr(DeseqStats, "lfc_shrink", orig_shrink)
    h5_one = tmp_path / "pb_one.h5ad"
    _write_pb_h5ad(h5_one, ["Tcell"])
    res_full = rna_pseudobulk_de({
        **base_params, "data_path": str(h5_one), "output_dir": str(tmp_path / "out_full")
    })
    assert res_full["status"] == "success"
    sk_full = res_full["lfc_shrinkage"]
    assert sk_full["requested"] is True
    assert sk_full["applied"] is True
    assert sk_full["n_applied"] == 1
    assert sk_full["n_failed"] == 0
    assert sk_full["failed_reasons"] == []
    assert "reported log2fc is the shrunken estimate" in sk_full["effect"]
