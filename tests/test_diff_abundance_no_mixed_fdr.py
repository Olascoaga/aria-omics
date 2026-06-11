"""B1: differential abundance must not mix CLR and Fisher p-values in one FDR."""

import pytest


def test_fisher_fallback_is_excluded_from_donor_level_bh(tmp_path, monkeypatch):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    sm = pytest.importorskip("statsmodels.api")

    rows = []
    counts = {
        "ctrl": {"A": 50, "B": 45, "C": 1},
        "treat": {"A": 50, "B": 45, "C": 80},
    }
    for cond in ("ctrl", "treat"):
        for donor in range(1, 4):
            for cell_type, n_cells in counts[cond].items():
                for _ in range(n_cells):
                    rows.append({
                        "condition": cond,
                        "donor": f"{cond}_{donor}",
                        "cell_type": cell_type,
                    })
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    adata = ad.AnnData(
        X=np.zeros((len(obs), 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    h5ad = tmp_path / "mixed_fdr.h5ad"
    adata.write_h5ad(h5ad)

    calls = []

    class FakeOLS:
        def __init__(self, *args, **kwargs):
            self.call_index = len(calls)
            calls.append(self.call_index)

        def fit(self, *args, **kwargs):
            if self.call_index == 0:
                pval = 0.06
                coef = 0.25
            elif self.call_index == 1:
                pval = 0.90
                coef = -0.05
            else:
                raise np.linalg.LinAlgError("forced zero residual df")
            return type(
                "FakeFit",
                (),
                {
                    "params": pd.Series({"is_test": coef}),
                    "pvalues": pd.Series({"is_test": pval}),
                },
            )()

    monkeypatch.setattr(sm, "OLS", FakeOLS)

    from aria.scripts.rna_diff_abundance import rna_diff_abundance

    result = rna_diff_abundance({
        "data_path": str(h5ad),
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["treat", "ctrl"]],
        "significance_alpha": 0.10,
    })

    assert result["status"] == "success", result
    comp = result["per_comparison"]["treat_vs_ctrl"]
    by_name = {row["name"]: row for row in comp["per_cell_type"]}

    # A is a donor-level CLR test with p=0.06. BH over donor-level tests only
    # (A and B) gives 0.12, so it must not become significant just because C's
    # cell-level Fisher diagnostic p-value is tiny.
    assert by_name["A"]["model"] == "donor_clr_ols_hc3"
    assert by_name["A"]["fdr_included"] is True
    assert by_name["A"]["padj"] == pytest.approx(0.12)
    assert by_name["A"]["significant"] is False

    assert by_name["C"]["model"] == "fisher_exact_fallback"
    assert by_name["C"]["fdr_included"] is False
    assert by_name["C"]["pval_role"] == "cell_level_diagnostic"
    assert by_name["C"]["padj"] is None
    assert by_name["C"]["significant"] is False

    assert comp["fdr_family"] == "donor_level_clr_only"
    assert comp["degraded"] is True
    assert "Fisher" in comp["degradation_reason"]
    assert comp["n_significant"] == 0
    assert result["any_significant"] is False
