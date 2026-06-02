"""P1-3: donor-level compositional differential abundance tests."""

import pytest


def test_clr_transform_rows_sum_to_zero():
    pd = pytest.importorskip("pandas")
    from aria.scripts.rna_diff_abundance import _clr_transform_counts

    counts = pd.DataFrame(
        [[10, 20, 40], [5, 5, 90]],
        index=["donor1", "donor2"],
        columns=["A", "B", "C"],
    )
    clr = _clr_transform_counts(counts)

    assert list(clr.columns) == ["A", "B", "C"]
    assert abs(float(clr.loc["donor1"].sum())) < 1e-10
    assert abs(float(clr.loc["donor2"].sum())) < 1e-10
    assert clr.loc["donor2", "C"] > clr.loc["donor2", "A"]


def test_clr_design_adds_donor_fixed_effect_for_paired_design():
    pd = pytest.importorskip("pandas")
    from aria.scripts.rna_diff_abundance import _build_clr_design

    reps = ["d1__ctrl", "d1__treat", "d2__ctrl", "d2__treat"]
    cond = pd.Series(["ctrl", "treat", "ctrl", "treat"], index=reps)
    design = _build_clr_design(
        reps,
        cond,
        "treat",
        covariates=[],
        rep_to_covs={},
        rep_to_base={
            "d1__ctrl": "d1",
            "d1__treat": "d1",
            "d2__ctrl": "d2",
            "d2__treat": "d2",
        },
        paired_design=True,
    )

    assert "is_test" in design.columns
    assert any(col.startswith("_aria_donor_") for col in design.columns)


def test_diff_abundance_uses_paired_clr_model(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    rng = np.random.default_rng(7)
    rows = []
    for donor in range(1, 7):
        donor_scale = int(rng.normal(0, 3))
        for cond in ("ctrl", "treat"):
            plan = {
                "A": 70 + donor_scale,
                "B": 20 + donor_scale,
                "C": 70 - donor_scale,
            }
            if cond == "treat":
                plan["B"] = 90 + donor_scale
                plan["A"] = 45 + donor_scale
                plan["C"] = 45 - donor_scale
            for cell_type, n_cells in plan.items():
                for _ in range(max(1, n_cells)):
                    rows.append({
                        "condition": cond,
                        "donor": f"donor_{donor}",
                        "cell_type": cell_type,
                    })
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    adata = ad.AnnData(
        X=np.zeros((len(obs), 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    h5ad = tmp_path / "paired_composition.h5ad"
    adata.write_h5ad(h5ad)

    from aria.scripts.rna_diff_abundance import rna_diff_abundance

    result = rna_diff_abundance({
        "data_path": str(h5ad),
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["treat", "ctrl"]],
    })

    assert result["status"] == "success", result
    assert result["method"] == "donor_clr_ols_hc3"
    assert result["paired_design"] is True
    comp = result["per_comparison"]["treat_vs_ctrl"]
    b_row = next(row for row in comp["per_cell_type"] if row["name"] == "B")
    assert b_row["model"] == "donor_clr_ols_hc3"
    assert b_row["clr_log_ratio"] > 0
    assert b_row["significant"] is True, b_row
    assert result["any_significant"] is True
