"""Analytical oracle tests for every Benjamini-Hochberg / FDR implementation in ARIA.

Covers:
  - aria.utils.stats: bh_correct, pooled_bh_across_groups, contrast_family_significance,
    preregister_fdr_family, preregister_contrast_family, primary_fdr_column,
    assert_fdr_family_not_post_hoc, fdr_advanced_methods_disclosure
  - aria.utils.lobo_jackknife.bh_adjust
  - aria.utils.cross_modal_concordance.bh_adjust
  - aria.scripts.hic_topology._bh_fdr
  - aria.scripts.chromatin_footprint_tobias._bh_adjust
  - aria.scripts.rna_pseudobulk_de._global_bh
  - aria.benchmarks.synthetic_de._bh_adjust
  - aria.scripts.rna_diff_abundance (BH step inside rna_diff_abundance)

Compares numerically against statsmodels.stats.multitest.multipletests(..., method="fdr_bh")
and scipy.stats.false_discovery_control(..., method="bh"), plus closed-form hand-derived values.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import false_discovery_control

from aria.benchmarks.synthetic_de import _bh_adjust as synthetic_de_bh_adjust
from aria.scripts.chromatin_footprint_tobias import _bh_adjust as tobias_bh_adjust
from aria.scripts.hic_topology import _bh_fdr as hic_bh_fdr
from aria.scripts.rna_pseudobulk_de import _global_bh as pseudobulk_global_bh
from aria.utils.cross_modal_concordance import bh_adjust as cross_modal_bh_adjust
from aria.utils.lobo_jackknife import bh_adjust as lobo_bh_adjust
from aria.utils.stats import (
    assert_fdr_family_not_post_hoc,
    bh_correct,
    contrast_family_significance,
    fdr_advanced_methods_disclosure,
    pooled_bh_across_groups,
    preregister_contrast_family,
    preregister_fdr_family,
    primary_fdr_column,
)


def _reference_bh_nan_aware(pvals) -> np.ndarray:
    """Independent reference BH using SciPy's false_discovery_control on finite values."""
    arr = np.asarray(pvals, dtype=float)
    out = np.full_like(arr, np.nan, dtype=float)
    valid = np.isfinite(arr)
    if valid.any():
        clipped = np.clip(arr[valid], 0.0, 1.0)
        out[valid] = false_discovery_control(clipped, method="bh")
    return out


ALL_BH_IMPLS = [
    pytest.param(bh_correct, id="stats.bh_correct"),
    pytest.param(lobo_bh_adjust, id="lobo_jackknife.bh_adjust"),
    pytest.param(cross_modal_bh_adjust, id="cross_modal_concordance.bh_adjust"),
    pytest.param(hic_bh_fdr, id="hic_topology._bh_fdr"),
    pytest.param(tobias_bh_adjust, id="chromatin_footprint_tobias._bh_adjust"),
    pytest.param(pseudobulk_global_bh, id="rna_pseudobulk_de._global_bh"),
    pytest.param(synthetic_de_bh_adjust, id="synthetic_de._bh_adjust"),
]

NAN_SAFE_BH_IMPLS = [
    pytest.param(bh_correct, id="stats.bh_correct"),
    pytest.param(lobo_bh_adjust, id="lobo_jackknife.bh_adjust"),
    pytest.param(cross_modal_bh_adjust, id="cross_modal_concordance.bh_adjust"),
    pytest.param(hic_bh_fdr, id="hic_topology._bh_fdr"),
    pytest.param(tobias_bh_adjust, id="chromatin_footprint_tobias._bh_adjust"),
]


class TestAllBHImplementationsFiniteOracle:
    """Every local BH implementation must match closed-form BH, SciPy, and statsmodels."""

    @pytest.mark.parametrize("bh_fn", ALL_BH_IMPLS)
    def test_hand_computed_step_up_with_monotonicity_enforcement(self, bh_fn):
        # Unsorted p-values where raw p_i * m / i is non-monotone:
        # Sorted: [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.612, 0.900] (m=10)
        # Raw p*m/k:
        #   k=1: 0.0100, k=2: 0.0400, k=3: 0.1300, k=4: 0.1025, k=5: 0.0840,
        #   k=6: 0.1000, k=7: 0.105714..., k=8: 0.25625, k=9: 0.6800, k=10: 0.9000
        # Backward cummin fixes k=3,4 to 0.0840 (from k=5)!
        p_sorted = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.612, 0.900]
        expected_sorted = np.array([
            0.010, 0.040, 0.084, 0.084, 0.084, 0.100, 0.74 / 7.0, 0.25625, 0.680, 0.900
        ])
        perm = [4, 0, 9, 2, 7, 1, 5, 8, 3, 6]
        p_perm = [p_sorted[i] for i in perm]
        expected_perm = expected_sorted[perm]

        got = np.asarray(bh_fn(p_perm), dtype=float)
        np.testing.assert_allclose(got, expected_perm, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(
            got, false_discovery_control(p_perm, method="bh"), rtol=1e-12, atol=1e-12
        )

    @pytest.mark.parametrize("bh_fn", ALL_BH_IMPLS)
    def test_matches_statsmodels_multipletests_with_ties_zeros_and_ones(self, bh_fn):
        sm_multitest = pytest.importorskip("statsmodels.stats.multitest")
        pvals = [0.0, 0.0, 0.01, 0.01, 0.01, 0.04, 0.2, 0.2, 0.85, 1.0, 1.0]
        rng = np.random.default_rng(17)
        perm = rng.permutation(len(pvals))
        p_in = [pvals[i] for i in perm]

        _, sm_padj, _, _ = sm_multitest.multipletests(p_in, method="fdr_bh")
        got = np.asarray(bh_fn(p_in), dtype=float)

        np.testing.assert_allclose(got, sm_padj, rtol=1e-12, atol=1e-12)
        assert got[np.asarray(p_in) == 0.0].tolist() == [0.0, 0.0]
        assert np.all(got[np.asarray(p_in) == 1.0] == 1.0)

    @pytest.mark.parametrize("bh_fn", ALL_BH_IMPLS)
    def test_permutation_invariance_on_random_and_tied_vectors(self, bh_fn):
        rng = np.random.default_rng(123)
        base = np.concatenate([
            rng.uniform(0.0, 0.05, size=15),
            np.array([0.0, 0.02, 0.02, 0.02, 0.5, 1.0]),
            rng.uniform(0.05, 1.0, size=20),
        ])
        got_base = np.asarray(bh_fn(base.tolist()), dtype=float)
        ref_base = false_discovery_control(base, method="bh")
        np.testing.assert_allclose(got_base, ref_base, rtol=1e-12, atol=1e-12)

        for seed in (1, 2, 3):
            perm = np.random.default_rng(seed).permutation(base.size)
            got_perm = np.asarray(bh_fn(base[perm].tolist()), dtype=float)
            np.testing.assert_allclose(got_perm, got_base[perm], rtol=1e-12, atol=1e-12)

    @pytest.mark.parametrize("bh_fn", ALL_BH_IMPLS)
    def test_single_element_and_empty_edge_cases(self, bh_fn):
        assert len(bh_fn([])) == 0
        for val in (0.0, 0.037, 1.0):
            got = np.asarray(bh_fn([val]), dtype=float)
            assert got.shape == (1,)
            assert got[0] == pytest.approx(val, abs=1e-12)

    @pytest.mark.parametrize("bh_fn", ALL_BH_IMPLS)
    def test_clipping_at_one_when_raw_adjustment_exceeds_one(self, bh_fn):
        pvals = [0.6, 0.7, 0.8]
        # m=3: rank1=1.8->1.0, rank2=1.05->1.0, rank3=0.8
        # cummin from right -> all 0.8!
        got = np.asarray(bh_fn(pvals), dtype=float)
        np.testing.assert_allclose(got, [0.8, 0.8, 0.8], rtol=1e-12, atol=1e-12)

        pvals2 = [0.55, 0.95]
        # m=2: rank1=1.10, rank2=0.95 -> cummin gives [0.95, 0.95]
        got2 = np.asarray(bh_fn(pvals2), dtype=float)
        np.testing.assert_allclose(got2, [0.95, 0.95], rtol=1e-12, atol=1e-12)

        pvals3 = [0.6, 1.0]
        # m=2: rank1=1.20 -> clipped to 1.0, rank2=1.0
        got3 = np.asarray(bh_fn(pvals3), dtype=float)
        np.testing.assert_allclose(got3, [1.0, 1.0], rtol=1e-12, atol=1e-12)


class TestBHNanAndNonFiniteHandling:
    """NaN and +/-inf positions must be preserved as NaN without poisoning finite p-values."""

    @pytest.mark.parametrize("bh_fn", NAN_SAFE_BH_IMPLS)
    def test_nan_and_inf_preserved_without_inflating_m(self, bh_fn):
        pvals = [0.01, float("nan"), 0.04, float("inf"), 0.20, float("-inf"), 0.0]
        # Finite subset (m=4): [0.0, 0.01, 0.04, 0.20]
        # Exact BH (m=4): 0.0 -> 0.0, 0.01*4/2=0.02, 0.04*4/3=0.0533333..., 0.20*4/4=0.20
        expected = _reference_bh_nan_aware(pvals)
        got = np.asarray(bh_fn(pvals), dtype=float)

        assert np.isnan(got[1]) and np.isnan(got[3]) and np.isnan(got[5])
        finite_idx = [0, 2, 4, 6]
        np.testing.assert_allclose(got[finite_idx], expected[finite_idx], rtol=1e-12, atol=1e-12)

    @pytest.mark.parametrize("bh_fn", NAN_SAFE_BH_IMPLS)
    def test_all_nan_returns_all_nan(self, bh_fn):
        pvals = [float("nan"), float("nan"), float("inf")]
        got = np.asarray(bh_fn(pvals), dtype=float)
        assert got.shape == (3,)
        assert np.all(np.isnan(got))

    @pytest.mark.parametrize("bh_fn", NAN_SAFE_BH_IMPLS)
    def test_out_of_range_finite_pvalues_are_clipped_to_unit_interval(self, bh_fn):
        pvals = [-0.1, 0.02, 1.5]
        # Clipped to [0.0, 0.02, 1.0], m=3 -> [0.0, 0.03, 1.0]
        got = np.asarray(bh_fn(pvals), dtype=float)
        np.testing.assert_allclose(got, [0.0, 0.03, 1.0], rtol=1e-12, atol=1e-12)

    # D2-BUG-1 fixed by D8-016 on fix/d1-silent-failures; xfail marker removed.
    def test_pseudobulk_global_bh_nan_does_not_poison_finite_values(self):
        pytest.importorskip("statsmodels")
        pvals = [0.001, float("nan"), 0.02, 0.5]
        got = np.asarray(pseudobulk_global_bh(pvals), dtype=float)
        expected = _reference_bh_nan_aware(pvals)
        assert np.isnan(got[1])
        np.testing.assert_allclose(got[[0, 2, 3]], expected[[0, 2, 3]], rtol=1e-12, atol=1e-12)

    # D2-BUG-2 fixed by D8-016 on fix/d1-silent-failures; xfail marker removed.
    def test_synthetic_de_bh_adjust_nan_does_not_poison_finite_values(self):
        pytest.importorskip("statsmodels")
        pvals = [0.001, float("nan"), 0.02, 0.5]
        got = np.asarray(synthetic_de_bh_adjust(pvals), dtype=float)
        expected = _reference_bh_nan_aware(pvals)
        assert np.isnan(got[1])
        np.testing.assert_allclose(got[[0, 2, 3]], expected[[0, 2, 3]], rtol=1e-12, atol=1e-12)


class TestStatsHelpersAnalyticalOracle:
    """Exact oracle tests for pooled_bh_across_groups, contrast_family_significance, and FDR guards."""

    def test_pooled_bh_matches_flat_scipy_and_statsmodels_exactly(self):
        sm_multitest = pytest.importorskip("statsmodels.stats.multitest")
        groups = {
            "KO_vs_WT": {
                "G1": 0.001,
                "G2": 0.020,
                "G3": 0.400,
                "G_nan": float("nan"),
                "G_none": None,
                "G_bool": True,
            },
            "Drug_vs_WT": {
                "G1": 0.004,
                "G2": 0.030,
                "G4": 0.900,
            },
        }
        pooled = pooled_bh_across_groups(groups)

        # G_none and G_bool are omitted; G_nan is preserved as NaN.
        assert "G_none" not in pooled["KO_vs_WT"]
        assert "G_bool" not in pooled["KO_vs_WT"]
        assert math.isnan(pooled["KO_vs_WT"]["G_nan"])

        # Finite family has m = 6 tests: [0.001, 0.004, 0.020, 0.030, 0.400, 0.900]
        finite_keys = [
            ("KO_vs_WT", "G1"),
            ("Drug_vs_WT", "G1"),
            ("KO_vs_WT", "G2"),
            ("Drug_vs_WT", "G2"),
            ("KO_vs_WT", "G3"),
            ("Drug_vs_WT", "G4"),
        ]
        raw_p = [0.001, 0.004, 0.020, 0.030, 0.400, 0.900]
        # Closed-form BH (m=6):
        #   0.001*6/1 = 0.006
        #   0.004*6/2 = 0.012
        #   0.020*6/3 = 0.040
        #   0.030*6/4 = 0.045
        #   0.400*6/5 = 0.480
        #   0.900*6/6 = 0.900
        exact = [0.006, 0.012, 0.040, 0.045, 0.480, 0.900]
        _, sm_padj, _, _ = sm_multitest.multipletests(raw_p, method="fdr_bh")
        sp_padj = false_discovery_control(raw_p, method="bh")

        got = [pooled[c][g] for c, g in finite_keys]
        np.testing.assert_allclose(got, exact, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got, sm_padj, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got, sp_padj, rtol=1e-12, atol=1e-12)

    def test_pooled_bh_empty_and_group_order_invariance(self):
        assert pooled_bh_across_groups({}) == {}
        assert pooled_bh_across_groups({"A": {"g1": None, "g2": False}}) == {"A": {}}

        g_ab = {
            "A": {"g1": 0.002, "g2": 0.05},
            "B": {"g1": 0.01, "g3": 0.25},
        }
        g_ba = {
            "B": {"g3": 0.25, "g1": 0.01},
            "A": {"g2": 0.05, "g1": 0.002},
        }
        r_ab = pooled_bh_across_groups(g_ab)
        r_ba = pooled_bh_across_groups(g_ba)
        for grp in ("A", "B"):
            for gene in g_ab[grp]:
                assert r_ab[grp][gene] == pytest.approx(r_ba[grp][gene], abs=1e-12)

    def test_contrast_family_significance_exact_boundaries_and_invalid_lfc_filtering(self):
        # m = 5 finite p-values:
        #   c1/g_up:   p=0.005 -> adj = 0.025 (< 0.05), lfc=+2.0 (> 1.0) -> sig up
        #   c1/g_dn:   p=0.008 -> adj = 0.020 -> wait sorted:
        #     rank1: 0.005 * 5 / 1 = 0.025 -> cummin from rank2 (0.020) => 0.020
        #     rank2: 0.008 * 5 / 2 = 0.020 => 0.020
        #     rank3: 0.029 * 5 / 3 = 0.048333... (< 0.05) -> lfc=0.8 (<= 1.0) -> filtered when lfc_min=1.0, sig when lfc_min=None
        #     rank4: 0.040 * 5 / 4 = 0.050 (== padj_max, strict `<` excludes!)
        #     rank5: 0.500 * 5 / 5 = 0.500
        #   Plus invalid lfc cases (NaN, inf, bool, non-numeric) that MUST NOT be called sig.
        stats = {
            "c1": {
                "g_up": {"pvalue": 0.005, "log2fc": 2.0},
                "g_dn": {"pvalue": 0.008, "log2fc": -1.5},
                "g_small": {"pvalue": 0.029, "log2fc": 0.8},
                "g_boundary_p": {"pvalue": 0.040, "log2fc": 3.0},
                "g_nan_lfc": {"pvalue": 1e-6, "log2fc": float("nan")},
                "g_bool_lfc": {"pvalue": 1e-6, "log2fc": True},
                "g_str_lfc": {"pvalue": 1e-6, "log2fc": "bad"},
            },
            "c2": {
                "g_null": {"pvalue": 0.500, "log2fc": -2.5},
            },
        }
        # Note: m = 8 finite p-values (including the three 1e-6 rows!).
        # Sorted p-values (m=8):
        #   1e-6, 1e-6, 1e-6, 0.005, 0.008, 0.029, 0.040, 0.500
        # Raw p*8/rank:
        #   r1..3: <= 8e-6
        #   r4 (0.005): 0.005*8/4 = 0.010
        #   r5 (0.008): 0.008*8/5 = 0.0128
        #   r6 (0.029): 0.029*8/6 = 0.0386666...
        #   r7 (0.040): 0.040*8/7 = 0.0457142... (< 0.05)
        #   r8 (0.500): 0.500
        res_lfc = contrast_family_significance(stats, padj_max=0.05, lfc_min=1.0)
        c1 = res_lfc["c1"]
        assert set(c1["sig_genes"]) == {"g_up", "g_dn", "g_boundary_p"}
        assert c1["n_sig"] == 3
        assert c1["n_up"] == 2
        assert c1["n_down"] == 1
        assert c1["padj_family"]["g_up"] == pytest.approx(0.010, abs=1e-12)
        assert c1["padj_family"]["g_dn"] == pytest.approx(0.0128, abs=1e-12)
        assert c1["padj_family"]["g_small"] == pytest.approx(0.029 * 8 / 6, abs=1e-12)
        assert res_lfc["c2"]["n_sig"] == 0

        # When lfc_min=None, g_small passes while invalid-lfc entries remain excluded.
        res_nolim = contrast_family_significance(stats, padj_max=0.05, lfc_min=None)
        assert set(res_nolim["c1"]["sig_genes"]) == {"g_up", "g_dn", "g_small", "g_boundary_p"}
        assert res_nolim["c1"]["n_up"] == 3
        assert res_nolim["c1"]["n_down"] == 1

        # Exact strict `< padj_max` and `abs(lfc) > lfc_min` boundary check.
        strict_stats = {
            "c": {
                "at_padj": {"pvalue": 0.05, "log2fc": 2.0},
                "at_lfc": {"pvalue": 0.001, "log2fc": 1.0},
            }
        }
        res_strict = contrast_family_significance(strict_stats, padj_max=0.05, lfc_min=1.0)
        assert res_strict["c"]["sig_genes"] == []
        assert res_strict["c"]["n_sig"] == 0

    def test_preregistration_and_post_hoc_guards_exact_behavior(self):
        assert preregister_fdr_family("global")["fdr_strategy"] == "global"
        assert preregister_fdr_family(" GLOBAL ")["fdr_strategy"] == "global"
        assert preregister_fdr_family("per_cluster")["fdr_strategy"] == "per_cluster"
        assert preregister_fdr_family("unknown_value")["fdr_strategy"] == "per_cluster"
        assert preregister_fdr_family(None)["fdr_strategy"] == "per_cluster"
        assert preregister_fdr_family("global")["preregistered"] is True
        assert preregister_fdr_family("global")["selected_before_results"] is True

        assert preregister_contrast_family("global")["fdr_family"] == "global"
        assert preregister_contrast_family("per_contrast")["fdr_family"] == "per_contrast"
        assert preregister_contrast_family("bogus")["fdr_family"] == "per_contrast"
        assert preregister_contrast_family("")["fdr_family"] == "per_contrast"

        assert primary_fdr_column("global") == "padj_global"
        assert primary_fdr_column("per_cluster") == "padj_local"
        assert primary_fdr_column("other") == "padj_local"

        assert_fdr_family_not_post_hoc("global", "padj_global")
        assert_fdr_family_not_post_hoc("per_cluster", "padj_local")
        with pytest.raises(ValueError, match="post-hoc switch"):
            assert_fdr_family_not_post_hoc("global", "padj_local")
        with pytest.raises(ValueError, match="post-hoc switch"):
            assert_fdr_family_not_post_hoc("per_cluster", "padj_global")

        disc = fdr_advanced_methods_disclosure()
        assert disc["ihw"]["available"] is False
        assert disc["ihw"]["status"] == "not_implemented"
        assert disc["s_values"]["available"] is False
        assert disc["s_values"]["status"] == "not_available"
        assert "Benjamini-Hochberg" in disc["primary_method"]


class TestRnaDiffAbundanceBHIntegration:
    """Verify the BH adjustment inside aria.scripts.rna_diff_abundance matches SciPy/statsmodels."""

    def test_rna_diff_abundance_padj_matches_reference_bh(self, tmp_path):
        ad = pytest.importorskip("anndata")
        sm_multitest = pytest.importorskip("statsmodels.stats.multitest")
        import pandas as pd
        from aria.scripts.rna_diff_abundance import rna_diff_abundance

        # 4 cell types across 8 donors (4 KO, 4 WT). Cell type CT0 expands ~4x in KO.
        rows_obs = []
        counts_per_donor = {
            "WT": [
                {"CT0": 20, "CT1": 60, "CT2": 60, "CT3": 60},
                {"CT0": 22, "CT1": 58, "CT2": 62, "CT3": 58},
                {"CT0": 18, "CT1": 62, "CT2": 58, "CT3": 62},
                {"CT0": 20, "CT1": 59, "CT2": 61, "CT3": 60},
            ],
            "KO": [
                {"CT0": 80, "CT1": 40, "CT2": 40, "CT3": 40},
                {"CT0": 84, "CT1": 38, "CT2": 40, "CT3": 38},
                {"CT0": 78, "CT1": 42, "CT2": 39, "CT3": 41},
                {"CT0": 82, "CT1": 39, "CT2": 41, "CT3": 38},
            ],
        }
        for cond, donor_list in counts_per_donor.items():
            for d_idx, ct_map in enumerate(donor_list):
                donor = f"{cond}_d{d_idx}"
                for ct, n_cells in ct_map.items():
                    for _ in range(n_cells):
                        rows_obs.append({"condition": cond, "donor": donor, "cell_type": ct})
        obs = pd.DataFrame(rows_obs)
        obs.index = [f"cell_{i}" for i in range(len(obs))]
        adata = ad.AnnData(
            X=np.ones((len(obs), 3), dtype=np.float32),
            obs=obs,
            var=pd.DataFrame(index=["G0", "G1", "G2"]),
        )
        h5ad = tmp_path / "da.h5ad"
        adata.write_h5ad(h5ad)

        res = rna_diff_abundance({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "alpha": 0.05,
        })
        assert res["status"] == "success", res
        comp = res["per_comparison"]["KO_vs_WT"]
        primary = [r for r in comp["per_cell_type"] if r["fdr_included"]]
        assert len(primary) == 4

        raw_p = [r["pval"] for r in primary]
        got_padj = [r["padj"] for r in primary]
        _, sm_padj, _, _ = sm_multitest.multipletests(raw_p, method="fdr_bh")
        sp_padj = false_discovery_control(raw_p, method="bh")
        np.testing.assert_allclose(got_padj, sm_padj, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got_padj, sp_padj, rtol=1e-12, atol=1e-12)

        # CT0 MUST be significantly UP in KO vs WT, and reverse WT vs KO flips sign.
        ct0 = next(r for r in primary if r["name"] == "CT0")
        assert ct0["significant"] is True
        assert ct0["direction"] == "up"
        assert ct0["clr_log_ratio"] > 0

        res_rev = rna_diff_abundance({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["WT", "KO"]],
            "alpha": 0.05,
        })
        ct0_rev = next(
            r for r in res_rev["per_comparison"]["WT_vs_KO"]["per_cell_type"]
            if r["name"] == "CT0"
        )
        assert ct0_rev["direction"] == "down"
        assert ct0_rev["clr_log_ratio"] == pytest.approx(-ct0["clr_log_ratio"], rel=1e-9, abs=1e-9)
        assert ct0_rev["padj"] == pytest.approx(ct0["padj"], rel=1e-9, abs=1e-9)

    def test_stats_module_mutant_killers_preregister_disclosure_and_edges(self):
        from aria.utils.stats import (
            assert_fdr_family_not_post_hoc,
            fdr_advanced_methods_disclosure,
            preregister_contrast_family,
            preregister_fdr_family,
            primary_fdr_column,
        )

        # preregister_fdr_family & preregister_contrast_family normalization & exact fields
        pf_g = preregister_fdr_family("  GLOBAL ")
        assert pf_g["fdr_strategy"] == "global"
        assert pf_g["preregistered"] is True
        assert pf_g["selected_before_results"] is True
        assert "not chosen post-hoc" in pf_g["note"]
        assert preregister_fdr_family(None)["fdr_strategy"] == "per_cluster"
        assert preregister_fdr_family("invalid")["fdr_strategy"] == "per_cluster"

        pc_g = preregister_contrast_family("  GLOBAL ")
        assert pc_g["fdr_family"] == "global"
        assert pc_g["preregistered"] is True
        assert pc_g["selected_before_results"] is True
        assert "not chosen post-hoc" in pc_g["note"]
        assert preregister_contrast_family(None)["fdr_family"] == "per_contrast"
        assert preregister_contrast_family("bad")["fdr_family"] == "per_contrast"

        # primary_fdr_column & assert_fdr_family_not_post_hoc
        assert primary_fdr_column(" GLOBAL ") == "padj_global"
        assert primary_fdr_column("per_cluster") == "padj_local"
        assert_fdr_family_not_post_hoc("global", "padj_global")
        with pytest.raises(ValueError, match="post-hoc switch"):
            assert_fdr_family_not_post_hoc("global", "padj_local")

        # fdr_advanced_methods_disclosure exact semantic contract
        disc = fdr_advanced_methods_disclosure()
        assert disc["ihw"]["status"] == "not_implemented"
        assert disc["ihw"]["available"] is False
        assert "Ignatiadis-Huber" in disc["ihw"]["reason"] and "ADR-002" in disc["ihw"]["reason"]
        assert disc["s_values"]["status"] == "not_available"
        assert disc["s_values"]["available"] is False
        assert "pydeseq2 0.5.4" in disc["s_values"]["reason"]
        assert disc["primary_method"] == "pre-registered Benjamini-Hochberg (local/global dual)"

        # pooled_bh_across_groups: None/inf/empty-dict edges
        assert pooled_bh_across_groups(None) == {}
        inf_out = pooled_bh_across_groups({"c1": {"g_inf": float("inf"), "g_ok": 0.02}, "c2": None})
        assert math.isnan(inf_out["c1"]["g_inf"])
        assert inf_out["c1"]["g_ok"] == pytest.approx(0.02, abs=1e-12)
        assert inf_out["c2"] == {}

        # contrast_family_significance: lfc == 0.0 (neither up nor down), None lfc with lfc_min=None,
        # inf lfc excluded, missing pvalue defaults padj=1.0, empty/None contrast dict.
        cf = contrast_family_significance(
            {
                "c1": {
                    "g_zero_lfc": {"pvalue": 1e-4, "log2fc": 0.0},
                    "g_none_lfc": {"pvalue": 1e-4, "log2fc": None},
                    "g_inf_lfc": {"pvalue": 1e-4, "log2fc": float("inf")},
                    "g_no_p": {"pvalue": None, "log2fc": 3.0},
                    "g_bool_p": {"pvalue": True, "log2fc": 3.0},
                },
                "c_empty": None,
            },
            padj_max=0.05,
            lfc_min=None,
        )
        assert set(cf["c1"]["sig_genes"]) == {"g_zero_lfc", "g_none_lfc"}
        assert cf["c1"]["n_sig"] == 2
        assert cf["c1"]["n_up"] == 0 and cf["c1"]["n_down"] == 0
        assert cf["c_empty"]["n_sig"] == 0
        assert contrast_family_significance(None, padj_max=0.05, lfc_min=0.5) == {}
