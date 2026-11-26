"""Analytical oracle tests for aria/utils/lobo_jackknife.py.

Verifies every mathematical component against independent closed-form SciPy formulas:
  - Unmoderated (M < 20) exact Welch t-test and paired t-test fold stats
  - Moderated (M >= 20) Smyth (2004) empirical-Bayes variance moderation (d0=4, s0^2=median)
  - Sign flip on contrast reversal (KO vs WT <-> WT vs KO)
  - Alphabetical level-renaming invariance (KO vs KO_ctrl, A_trt vs Z_ctrl, Z_trt vs A_ctrl)
  - Permutation invariance over biological units and features
  - Exact single-unit outlier collapse detection (identifying the exact outlier unit)
  - Exact spike-in size-factor normalization (dividing raw / sf, NOT multiplying!)
  - Reordered pandas Series size_factors index alignment
  - Paired, partially paired, and unpaired design selection + biological unit enforcement
  - Edge cases: <3 replicates per arm, duplicate unit IDs, shape mismatch, zero/NaN values
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.stats import false_discovery_control

from aria.utils.lobo_jackknife import (
    _donor_id,
    _paired_fold_stats,
    _unpaired_fold_stats,
    compute_counts_lobo_jackknife,
    compute_lobo_jackknife,
)


def _closed_form_unpaired_moderated(mat: np.ndarray, n_idx: list[int], d_idx: list[int]):
    """Independent analytical Smyth-2004 moderated two-sample t-test with d0=4."""
    x1 = mat[:, n_idx]
    x2 = mat[:, d_idx]
    n1, n2 = x1.shape[1], x2.shape[1]
    lfc = x1.mean(axis=1) - x2.mean(axis=1)
    v1 = x1.var(axis=1, ddof=1)
    v2 = x2.var(axis=1, ddof=1)
    df = float(n1 + n2 - 2)
    s2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    s0_sq = float(np.median(s2[s2 > 1e-12]))
    d0 = 4.0
    s2_post = (d0 * s0_sq + df * s2) / (d0 + df)
    se = np.sqrt(np.maximum(s2_post * (1.0 / n1 + 1.0 / n2), 1e-15))
    t_stat = lfc / se
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df + d0)
    padj = false_discovery_control(pvals, method="bh")
    return lfc, pvals, padj


def _closed_form_paired_moderated(diff_mat: np.ndarray, pair_idx: list[int]):
    """Independent analytical Smyth-2004 moderated paired t-test with d0=4."""
    sub = diff_mat[:, pair_idx]
    k = sub.shape[1]
    lfc = sub.mean(axis=1)
    v = sub.var(axis=1, ddof=1)
    df = float(k - 1)
    s0_sq = float(np.median(v[v > 1e-12]))
    d0 = 4.0
    s2_post = (d0 * s0_sq + df * v) / (d0 + df)
    se = np.sqrt(np.maximum(s2_post / k, 1e-15))
    t_stat = lfc / se
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df + d0)
    padj = false_discovery_control(pvals, method="bh")
    return lfc, pvals, padj


class TestLoboAnalyticalFoldStats:
    def test_unpaired_unmoderated_matches_scipy_welch_exactly(self):
        # M = 5 (< 20 -> unmoderated exact Welch)
        mat = np.array([
            [10.0, 10.5, 9.8, 10.2, 4.0, 4.3, 3.9, 4.1],   # strong up (+6.05)
            [2.0, 2.2, 1.9, 2.1, 8.0, 8.3, 7.8, 7.9],       # strong down (-5.95)
            [5.0, 5.2, 4.8, 5.1, 5.1, 4.9, 5.0, 5.2],       # null
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],       # zero variance both arms
            [1.0, 1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 1.0],       # single outlier
        ])
        lfc, pvals, padj = _unpaired_fold_stats(mat, [0, 1, 2, 3], [4, 5, 6, 7], moderate=False)

        exp_lfc = mat[:, :4].mean(axis=1) - mat[:, 4:].mean(axis=1)
        _, raw_p = stats.ttest_ind(mat[:, :4], mat[:, 4:], axis=1, equal_var=False)
        exp_p = np.where(np.isfinite(raw_p), raw_p, 1.0)
        # Note: bh_adjust masks NaN before computing m, then NaNs become 1.0!
        valid = np.isfinite(raw_p)
        exp_padj = np.ones_like(raw_p)
        exp_padj[valid] = false_discovery_control(raw_p[valid], method="bh")

        np.testing.assert_allclose(lfc, exp_lfc, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(pvals, exp_p, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(padj, exp_padj, rtol=1e-12, atol=1e-12)

    def test_unpaired_moderated_matches_closed_form_smyth2004(self):
        rng = np.random.default_rng(42)
        mat = rng.normal(loc=5.0, scale=0.4, size=(25, 8))
        mat[0, :4] += 3.5   # strong up
        mat[1, :4] -= 3.5   # strong down
        lfc, pvals, padj = _unpaired_fold_stats(mat, [0, 1, 2, 3], [4, 5, 6, 7], moderate=True)
        exp_lfc, exp_p, exp_padj = _closed_form_unpaired_moderated(mat, [0, 1, 2, 3], [4, 5, 6, 7])
        np.testing.assert_allclose(lfc, exp_lfc, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(pvals, exp_p, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(padj, exp_padj, rtol=1e-12, atol=1e-12)

    def test_paired_unmoderated_and_moderated_match_closed_form(self):
        rng = np.random.default_rng(99)
        diff_mat = rng.normal(loc=0.0, scale=0.3, size=(24, 5))
        diff_mat[0, :] = [2.0, 2.1, 1.9, 2.05, 1.95]   # consistent up
        diff_mat[1, :] = [-2.5, -2.5, -2.5, -2.5, -2.5]  # constant nonzero diff (v=0)
        diff_mat[2, :] = [0.0, 0.0, 0.0, 0.0, 0.0]     # constant zero diff

        # Unmoderated on first 5 rows
        lfc_u, p_u, padj_u = _paired_fold_stats(diff_mat[:5], [0, 1, 2, 3, 4], moderate=False)
        exp_lfc_u = diff_mat[:5].mean(axis=1)
        np.testing.assert_allclose(lfc_u, exp_lfc_u, rtol=1e-12, atol=1e-12)
        # Row 0: exact scipy ttest_1samp
        _, p0 = stats.ttest_1samp(diff_mat[0], popmean=0.0)
        assert p_u[0] == pytest.approx(p0, rel=1e-12, abs=1e-12)
        # Row 1: v=0, lfc!=0 -> t=inf -> p=0.0
        assert p_u[1] == 0.0
        # Row 2: v=0, lfc==0 -> t=nan -> p=1.0, padj=1.0
        assert p_u[2] == 1.0
        assert padj_u[2] == 1.0

        # Moderated on all 24 rows
        lfc_m, p_m, padj_m = _paired_fold_stats(diff_mat, [0, 1, 2, 3, 4], moderate=True)
        exp_lfc_m, exp_p_m, exp_padj_m = _closed_form_paired_moderated(diff_mat, [0, 1, 2, 3, 4])
        np.testing.assert_allclose(lfc_m, exp_lfc_m, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(p_m, exp_p_m, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(padj_m, exp_padj_m, rtol=1e-12, atol=1e-12)


class TestLoboEndToEndOracles:
    def test_exact_stable_vs_single_outlier_driven_and_worst_unit_identification(self):
        # 25 features (>=20 -> moderated), 4 KO vs 4 WT.
        rng = np.random.default_rng(11)
        mat = rng.normal(loc=6.0, scale=0.25, size=(25, 8))
        # F_stable_up: every KO replicate is +3.0 above WT
        mat[0, :4] = [9.0, 9.1, 8.9, 9.05]
        mat[0, 4:] = [6.0, 5.9, 6.1, 6.0]
        # F_stable_dn: every KO replicate is -3.0 below WT
        mat[1, :4] = [3.0, 2.9, 3.1, 3.05]
        mat[1, 4:] = [6.0, 6.1, 5.9, 6.0]
        # F_outlier: ONLY KO_3 is inflated (+10); KO_0..KO_2 match WT (6.0)!
        # Full LFC = (16+6+6+6)/4 - 6 = +2.5; dropping KO_3 gives LFC = 0.0 (complete collapse!)
        mat[2, :4] = [6.0, 6.0, 6.0, 16.0]
        mat[2, 4:] = [6.0, 6.05, 5.95, 6.0]

        features = [f"F{i}" for i in range(25)]
        units = ["KO_0", "KO_1", "KO_2", "KO_3", "WT_0", "WT_1", "WT_2", "WT_3"]
        u2g = {u: ("KO" if u.startswith("KO") else "WT") for u in units}

        res = compute_lobo_jackknife(
            features,
            mat,
            units,
            u2g,
            numerator="KO",
            denominator="WT",
            alpha=0.05,
            target_features=["F0", "F1", "F2"],
            primary_significant_features=["F0", "F1", "F2"],
        )
        assert res["status"] == "measured"
        assert res["variance_moderation"] == "empirical_bayes_smyth2004"
        assert res["design_mode"] == "unpaired_biological_replicate"
        assert res["n_jackknife_folds"] == 8
        assert res["n_features_lobo_stable"] == 2
        assert res["n_single_unit_driven"] == 1
        assert res["lobo_stable_fraction"] == pytest.approx(2.0 / 3.0, abs=1e-4)

        f0 = res["per_feature"]["F0"]
        f1 = res["per_feature"]["F1"]
        f2 = res["per_feature"]["F2"]
        assert f0["full_lfc"] == pytest.approx(3.0125, abs=1e-6)
        assert f0["sign_consistent_fraction"] == 1.0
        assert f0["single_unit_driven"] is False

        assert f1["full_lfc"] == pytest.approx(-2.9875, abs=1e-6)
        assert f1["sign_consistent_fraction"] == 1.0
        assert f1["single_unit_driven"] is False

        assert f2["full_lfc"] == pytest.approx(2.5, abs=1e-6)
        assert f2["lfc_min"] == pytest.approx(0.0, abs=1e-6)
        assert f2["single_unit_driven"] is True
        assert f2["most_influential_unit"] == "KO_3"

    def test_contrast_reversal_and_alphabetical_level_renaming_invariance(self):
        rng = np.random.default_rng(21)
        mat = rng.normal(loc=5.0, scale=0.3, size=(22, 6))
        mat[0, :3] += 3.0
        mat[1, :3] -= 2.5
        features = [f"G{i}" for i in range(22)]
        units = ["u0", "u1", "u2", "u3", "u4", "u5"]

        u2g_kw = {u: ("KO" if i < 3 else "WT") for i, u in enumerate(units)}
        res_ko_vs_wt = compute_lobo_jackknife(
            features, mat, units, u2g_kw, numerator="KO", denominator="WT"
        )
        res_wt_vs_ko = compute_lobo_jackknife(
            features, mat, units, u2g_kw, numerator="WT", denominator="KO"
        )

        # Exact sign flip and padj equality across EVERY feature
        for f in features:
            r1 = res_ko_vs_wt["per_feature"][f]
            r2 = res_wt_vs_ko["per_feature"][f]
            assert r1["full_lfc"] == pytest.approx(-r2["full_lfc"], abs=1e-6)
            assert r1["lfc_min"] == pytest.approx(-r2["lfc_max"], abs=1e-6)
            assert r1["lfc_max"] == pytest.approx(-r2["lfc_min"], abs=1e-6)
            assert r1["full_padj"] == pytest.approx(r2["full_padj"], rel=1e-12, abs=1e-12)
            assert r1["max_jackknife_padj"] == pytest.approx(
                r2["max_jackknife_padj"], rel=1e-12, abs=1e-12
            )
            assert r1["single_unit_driven"] == r2["single_unit_driven"]

        # Alphabetical & prefix traps: (KO vs KO_ctrl), (Z_trt vs A_ctrl), (A_trt vs Z_ctrl)
        for num_lvl, den_lvl in [("KO", "KO_ctrl"), ("Z_trt", "A_ctrl"), ("A_trt", "Z_ctrl")]:
            u2g_alt = {u: (num_lvl if i < 3 else den_lvl) for i, u in enumerate(units)}
            res_alt = compute_lobo_jackknife(
                features, mat, units, u2g_alt, numerator=num_lvl, denominator=den_lvl
            )
            for f in features[:5]:
                assert res_alt["per_feature"][f]["full_lfc"] == pytest.approx(
                    res_ko_vs_wt["per_feature"][f]["full_lfc"], abs=1e-6
                )
                assert res_alt["per_feature"][f]["full_padj"] == pytest.approx(
                    res_ko_vs_wt["per_feature"][f]["full_padj"], rel=1e-12, abs=1e-12
                )

    def test_sample_and_feature_permutation_invariance(self):
        rng = np.random.default_rng(33)
        mat = rng.normal(loc=4.0, scale=0.4, size=(24, 8))
        mat[0, :4] += 2.8
        mat[1, :4] -= 2.2
        features = [f"P{i}" for i in range(24)]
        units = [f"S{i}" for i in range(8)]
        u2g = {u: ("TRT" if i < 4 else "CTRL") for i, u in enumerate(units)}

        base = compute_lobo_jackknife(features, mat, units, u2g, numerator="TRT", denominator="CTRL")

        f_perm = rng.permutation(len(features))
        u_perm = rng.permutation(len(units))
        mat_p = mat[np.ix_(f_perm, u_perm)]
        feat_p = [features[i] for i in f_perm]
        unit_p = [units[i] for i in u_perm]

        perm_res = compute_lobo_jackknife(
            feat_p, mat_p, unit_p, u2g, numerator="TRT", denominator="CTRL"
        )
        assert perm_res["n_features_lobo_stable"] == base["n_features_lobo_stable"]
        assert perm_res["n_single_unit_driven"] == base["n_single_unit_driven"]
        for f in features:
            b = base["per_feature"][f]
            p = perm_res["per_feature"][f]
            assert p["full_lfc"] == pytest.approx(b["full_lfc"], abs=1e-6)
            assert p["full_padj"] == pytest.approx(b["full_padj"], rel=1e-11, abs=1e-11)
            assert p["lfc_min"] == pytest.approx(b["lfc_min"], abs=1e-6)
            assert p["lfc_max"] == pytest.approx(b["lfc_max"], abs=1e-6)
            assert p["lfc_iqr"] == pytest.approx(b["lfc_iqr"], abs=1e-6)
            assert p["single_unit_driven"] == b["single_unit_driven"]

    def test_paired_donor_removes_donor_baseline_and_recovers_effect(self):
        # Huge donor baseline shift (+10, +20, +30, +40) hides a clean +1.5 treatment effect
        # in an unpaired fit, whereas paired_donor mode recovers exact LFC=+1.5 with tiny padj.
        donors = ["D1", "D2", "D3", "D4"]
        baselines = [10.0, 20.0, 30.0, 40.0]
        units, u2g, u2d = [], {}, {}
        rng = np.random.default_rng(77)
        mat = np.zeros((22, 8), dtype=float)
        for j, (d, base) in enumerate(zip(donors, baselines)):
            u_c, u_t = f"{d}_WT", f"{d}_KO"
            units.extend([u_c, u_t])
            u2g[u_c], u2g[u_t] = "WT", "KO"
            u2d[u_c], u2d[u_t] = d, d
            noise = rng.normal(0.0, 0.15, size=22)
            mat[:, 2 * j] = base + noise
            mat[:, 2 * j + 1] = base + noise + rng.normal(0.0, 0.05, size=22)
        # Plant exact +1.5 shift on feature 0 across all 4 donors
        for j in range(4):
            mat[0, 2 * j + 1] = mat[0, 2 * j] + 1.5

        features = [f"G{i}" for i in range(22)]
        res_paired = compute_lobo_jackknife(
            features, mat, units, u2g, numerator="KO", denominator="WT", unit_to_donor=u2d
        )
        res_unpaired = compute_lobo_jackknife(
            features, mat, units, u2g, numerator="KO", denominator="WT"
        )
        assert res_paired["design_mode"] == "paired_donor"
        assert res_paired["n_jackknife_folds"] == 4
        assert res_paired["n_unpaired_units_excluded"] == 0
        assert res_paired["per_feature"]["G0"]["full_lfc"] == pytest.approx(1.5, abs=1e-6)
        assert res_paired["per_feature"]["G0"]["full_padj"] < 1e-4
        assert res_paired["per_feature"]["G0"]["single_unit_driven"] is False
        # Unpaired is swamped by donor variance (padj >> paired padj)
        assert res_unpaired["per_feature"]["G0"]["full_padj"] > 0.5

    def test_edge_cases_replicates_duplicates_shapes_and_donor_id(self):
        # <3 replicates per arm -> not_testable_insufficient_replicates
        mat = np.ones((5, 4))
        res_2v2 = compute_lobo_jackknife(
            [f"F{i}" for i in range(5)],
            mat,
            ["a1", "a2", "b1", "b2"],
            {"a1": "KO", "a2": "KO", "b1": "WT", "b2": "WT"},
            numerator="KO",
            denominator="WT",
        )
        assert res_2v2["status"] == "not_testable_insufficient_replicates"
        assert res_2v2["minimum_required_per_arm"] == 3
        assert res_2v2["biological_unit_enforced"] is True
        assert res_2v2["numerator"] == "KO" and res_2v2["denominator"] == "WT"
        assert res_2v2["n_units_numerator"] == 2
        assert res_2v2["n_units_denominator"] == 2
        assert res_2v2["per_feature"] == {}

        # Single group -> 0 in denominator -> not_testable
        res_1grp = compute_lobo_jackknife(
            [f"F{i}" for i in range(5)],
            mat,
            ["a1", "a2", "a3", "a4"],
            {"a1": "KO", "a2": "KO", "a3": "KO", "a4": "KO"},
            numerator="KO",
            denominator="WT",
        )
        assert res_1grp["status"] == "not_testable_insufficient_replicates"
        assert res_1grp["n_units_denominator"] == 0

        # Paired mode with only 2 shared donors -> not_testable_insufficient_replicates
        res_2pair = compute_lobo_jackknife(
            [f"F{i}" for i in range(5)],
            mat,
            ["d1_k", "d1_w", "d2_k", "d2_w"],
            {"d1_k": "KO", "d1_w": "WT", "d2_k": "KO", "d2_w": "WT"},
            numerator="KO",
            denominator="WT",
            unit_to_donor={"d1_k": "D1", "d1_w": "D1", "d2_k": "D2", "d2_w": "D2"},
        )
        assert res_2pair["status"] == "not_testable_insufficient_replicates"
        assert res_2pair["design_mode"] == "paired_donor"

        # Duplicate biological unit IDs MUST raise ValueError
        with pytest.raises(ValueError, match="unique biological units"):
            compute_lobo_jackknife(
                ["F0", "F1"], np.ones((2, 2)), ["dup", "dup"], {"dup": "KO"},
                numerator="KO", denominator="WT",
            )

        # Shape mismatch MUST raise ValueError
        with pytest.raises(ValueError, match="shape"):
            compute_lobo_jackknife(
                ["F0", "F1"], np.ones((3, 2)), ["u1", "u2"], {"u1": "KO", "u2": "WT"},
                numerator="KO", denominator="WT",
            )

        # _donor_id normalization
        for bad in (None, True, False, float("nan"), float("inf"), "", "  ", "NA", "nan", "None", "null", "NULL", "<NA>"):
            assert _donor_id(bad) is None
        assert _donor_id(" Donor_01 ") == "Donor_01"
        assert _donor_id(42) == "42"

    def test_lobo_comprehensive_mutant_killers(self):
        from scipy import stats as sp_stats
        from aria.utils.lobo_jackknife import (
            _paired_fold_stats,
            _unpaired_fold_stats,
            bh_adjust,
        )

        # 1. bh_adjust edge cases: single finite element (m==1), p<0 / p>1 clipping, adj>1 clipping
        np.testing.assert_allclose(bh_adjust(np.array([0.04])), [0.04], atol=1e-12)
        np.testing.assert_allclose(bh_adjust(np.array([-0.5, 0.2])), [0.0, 0.2], atol=1e-12)
        np.testing.assert_allclose(bh_adjust(np.array([0.01, 1.5])), [0.02, 1.0], atol=1e-12)
        np.testing.assert_allclose(bh_adjust(np.array([0.6, 0.9])), [0.9, 0.9], atol=1e-12)
        np.testing.assert_allclose(bh_adjust(np.array([0.8, 0.9])), [0.9, 0.9], atol=1e-12)
        np.testing.assert_allclose(bh_adjust(np.array([0.8, 1.2])), [1.0, 1.0], atol=1e-12)

        # 2. _unpaired_fold_stats & _paired_fold_stats moderated s0_sq fallback (pos.size==0 & ==1)
        #    When all rows have zero variance, s0_sq MUST be 1e-4 -> exact analytical moderated p-value.
        mat_const = np.array([
            [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            [0.05, 0.05, 0.05, 0.0, 0.0, 0.0],
        ])
        lfc_c, p_c, q_c = _unpaired_fold_stats(mat_const, [0, 1, 2], [3, 4, 5], moderate=True)
        # df_obs = 4, d0 = 4, s0_sq = 1e-4 -> s2_post = 4e-4 / 8 = 5e-5, se = sqrt(5e-5 * 2/3)
        exp_se = float(np.sqrt(5e-5 * (1.0 / 3.0 + 1.0 / 3.0)))
        exp_p1 = float(2.0 * sp_stats.t.sf(0.05 / exp_se, df=8.0))
        assert p_c[1] == pytest.approx(exp_p1, rel=1e-9, abs=1e-12)

        # Single non-zero-variance row (pos.size == 1) uses that row's variance, NOT 1e-4
        mat_one = np.array([
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0],  # v1=1.0, v2=0 -> s2_pooled = 2*1/4 = 0.5
            [0.5, 0.5, 0.5, 0.0, 0.0, 0.0],  # v1=0, v2=0 -> uses s0_sq = 0.5
        ])
        _, p_one, _ = _unpaired_fold_stats(mat_one, [0, 1, 2], [3, 4, 5], moderate=True)
        se_row1 = float(np.sqrt((4.0 * 0.5 / 8.0) * (2.0 / 3.0)))
        assert p_one[1] == pytest.approx(float(2.0 * sp_stats.t.sf(0.5 / se_row1, df=8.0)), rel=1e-9)

        # Paired moderated zero-variance fallback (s0_sq = 1e-4) & unmoderated 0 < |lfc| < 1
        diff_small = np.array([[0.05, 0.05, 0.05, 0.05], [0.0, 0.0, 0.0, 0.0]])
        _, pp_mod, _ = _paired_fold_stats(diff_small, [0, 1, 2, 3], moderate=True)
        # df_obs = 3, d0 = 4, s2_post = 4e-4 / 7, se = sqrt(s2_post / 4)
        se_p = float(np.sqrt((4e-4 / 7.0) / 4.0))
        assert pp_mod[0] == pytest.approx(float(2.0 * sp_stats.t.sf(0.05 / se_p, df=7.0)), rel=1e-9)
        _, pp_unmod, _ = _paired_fold_stats(diff_small, [0, 1, 2, 3], moderate=False)
        assert pp_unmod[0] == 0.0 and pp_unmod[1] == 1.0

        # 3. Exact moderation boundary (n_features_universe == 20 vs 19) + full schema keys +
        #    padj_collapse, severe_effect_collapse, non-sig sign flip, exact lfc_iqr & rounding
        rng = np.random.default_rng(314)
        units = [f"U{i}" for i in range(8)]
        u2g = {u: ("KO" if i < 4 else "WT") for i, u in enumerate(units)}
        mat19 = rng.normal(5.0, 0.2, size=(19, 8))
        res19 = compute_lobo_jackknife(
            [f"G{i}" for i in range(19)], mat19, units, u2g, numerator="KO", denominator="WT"
        )
        assert res19["variance_moderation"] == "exact_welch_or_paired_t"
        assert res19["biological_unit_enforced"] is True
        assert res19["numerator"] == "KO" and res19["denominator"] == "WT"
        assert res19["n_features_universe"] == 19
        assert res19["design_mode"] == "unpaired_biological_replicate"
        assert res19["n_units_numerator"] == 4 and res19["n_units_denominator"] == 4

        mat20 = rng.normal(5.0, 0.2, size=(20, 8))
        # G0: rock-solid significant hit
        mat20[0, :4] = [8.0, 8.1, 7.9, 8.0]
        mat20[0, 4:] = [5.0, 5.1, 4.9, 5.0]
        # G1: non-significant with sign flip across folds (mean KO slightly > WT only via one unit)
        mat20[1, :4] = [4.9, 4.9, 4.9, 5.6]
        mat20[1, 4:] = [5.0, 5.0, 5.0, 5.0]
        # G2: severe_effect_collapse with non-significant full padj (abs_full=0.6 > 0.25,
        #     dropping unit 3 drops |lfc| to 0.1 < 0.35*0.6 = 0.21, sign_frac still 1.0!)
        mat20[2, :4] = [5.1, 5.1, 5.1, 7.1]
        mat20[2, 4:] = [5.0, 5.0, 5.0, 5.0]
        # G3: stable non-significant (sign_frac == 1.0, min_abs_lfc >= 0.35*abs_full) -> stable!
        mat20[3, :4] = [5.15, 5.18, 5.12, 5.15]
        mat20[3, 4:] = [5.0, 5.02, 4.98, 5.0]
        # G4: padj_collapse (borderline full padj < 0.05, dropping best unit pushes padj >= 0.10,
        #     while sign_frac == 1.0 and min_abs_lfc >= 0.35 * abs_full!)
        mat20[4, :4] = [5.38, 5.40, 5.42, 6.20]
        mat20[4, 4:] = [5.0, 5.05, 4.95, 5.0]

        res20 = compute_lobo_jackknife(
            [f"G{i}" for i in range(20)], mat20, units, u2g, numerator="KO", denominator="WT",
            target_features=["G0", "G1", "G2", "G3", "G4"],
        )
        assert res20["variance_moderation"] == "empirical_bayes_smyth2004"
        assert res20["n_features_universe"] == 20
        pf = res20["per_feature"]
        assert pf["G0"]["single_unit_driven"] is False
        assert pf["G0"]["sign_consistent_fraction"] == 1.0
        assert pf["G1"]["single_unit_driven"] is True
        assert pf["G1"]["sign_consistent_fraction"] == pytest.approx(7.0 / 8.0, abs=1e-4)
        assert pf["G2"]["single_unit_driven"] is True
        assert pf["G2"]["sign_consistent_fraction"] == 1.0
        assert pf["G3"]["single_unit_driven"] is False
        # Verify exact lfc_iqr against independent np.percentile([75, 25]) on leave-one-out fold LFCs
        fold_lfcs_g2 = []
        for drop in range(8):
            k_vals = [mat20[2, j] for j in range(4) if j != drop]
            w_vals = [mat20[2, j] for j in range(4, 8) if j != drop]
            fold_lfcs_g2.append(float(np.mean(k_vals) - np.mean(w_vals)))
        q75, q25 = np.percentile(fold_lfcs_g2, [75, 25])
        assert pf["G2"]["lfc_iqr"] == pytest.approx(round(float(q75 - q25), 6), abs=1e-9)
        assert res20["lobo_stable_fraction"] == pytest.approx(
            round(res20["n_features_lobo_stable"] / res20["n_features_evaluated"], 4), abs=1e-9
        )

        # Empty target_features overlap -> n_eval == 0 -> lobo_stable_fraction is None
        res_empty_eval = compute_lobo_jackknife(
            [f"G{i}" for i in range(20)], mat20, units, u2g, numerator="KO", denominator="WT",
            target_features=["NO_SUCH_GENE"],
        )
        assert res_empty_eval["n_features_evaluated"] == 0
        assert res_empty_eval["lobo_stable_fraction"] is None

        # Partially paired with >= 3 paired donors + 1 extra unpaired -> paired_donor with n_unpaired_units_excluded == 1
        u_p3 = ["d1_k", "d1_w", "d2_k", "d2_w", "d3_k", "d3_w", "extra_k"]
        g_p3 = {u: ("KO" if u.endswith("_k") else "WT") for u in u_p3}
        d_p3 = {"d1_k": "D1", "d1_w": "D1", "d2_k": "D2", "d2_w": "D2", "d3_k": "D3", "d3_w": "D3", "extra_k": "D4"}
        res_p3 = compute_lobo_jackknife(
            [f"F{i}" for i in range(4)], rng.normal(5.0, 0.2, size=(4, 7)),
            u_p3, g_p3, numerator="KO", denominator="WT", unit_to_donor=d_p3,
        )
        assert res_p3["design_mode"] == "paired_donor"
        assert res_p3["n_units_numerator"] == 3 and res_p3["n_units_denominator"] == 3
        assert res_p3["n_unpaired_units_excluded"] == 1

        # 4. compute_counts_lobo_jackknife: list/array size_factors, zero-count column, NaN/neg counts,
        #    target_features / max_reported_features / alpha forwarding, and zero size_factor fallback
        s6 = ["k1", "k2", "k3", "w1", "w2", "w3"]
        c6 = pd.DataFrame(
            {
                "k1": [100.0, 50.0, float("nan")],
                "k2": [105.0, 48.0, -5.0],
                "k3": [0.0, 0.0, 0.0],  # all-zero library -> lib_sizes floor 1.0
                "w1": [20.0, 50.0, 0.0],
                "w2": [22.0, 49.0, 0.0],
                "w3": [18.0, 51.0, 0.0],
            },
            index=["P0", "P1", "P_bad"],
        )
        m6 = pd.DataFrame({"condition": ["KO"] * 3 + ["WT"] * 3}, index=s6)
        res_cpm = compute_counts_lobo_jackknife(
            c6, m6, "condition", "KO", "WT",
            target_features=["P0", "P_bad"],
            max_reported_features=1,
            alpha=0.01,
        )
        assert res_cpm["normalization"] == "log2_cpm_plus_1"
        assert res_cpm["n_features_evaluated"] == 2
        assert res_cpm["per_feature_truncated"] is True
        assert list(res_cpm["per_feature"].keys()) == ["P0"]
        # P_bad had NaN and -5 -> clamped to 0.0 everywhere -> full_lfc == 0.0
        res_cpm_all = compute_counts_lobo_jackknife(c6, m6, "condition", "KO", "WT")
        assert res_cpm_all["per_feature"]["P_bad"]["full_lfc"] == 0.0

        # Sequence (list) size_factors path + exact 0.0 rejected to CPM fallback
        res_seq = compute_counts_lobo_jackknife(
            c6, m6, "condition", "KO", "WT", size_factors=[0.5, 0.5, 0.5, 2.0, 2.0, 2.0]
        )
        assert res_seq["normalization"] == "spikein_or_external_size_factor_log2"
        assert res_seq["per_feature"]["P0"]["full_lfc"] > res_cpm_all["per_feature"]["P0"]["full_lfc"]
        res_zero_sf = compute_counts_lobo_jackknife(
            c6, m6, "condition", "KO", "WT", size_factors=[0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        assert res_zero_sf["normalization"] == "log2_cpm_plus_1"


class TestCountsLoboSpikeInOracle:
    """Independent closed-form verification of compute_counts_lobo_jackknife & spike-in scaling.

    Directly kills D3-008 mutation (`raw / sf_vec` -> `raw * sf_vec`).
    """

    def test_spikein_normalization_matches_closed_form_log2_ratio_and_kills_multiply_mutant(self):
        # 6 biological units (3 KO, 3 WT), 4 features.
        # Global 4x chromatin loss in KO: KO has raw counts ~ 100, WT has raw counts ~ 400.
        # Spike-in size factors reflect equal input DNA per cell (sf_KO = 0.5, sf_WT = 2.0):
        # Wait: make raw_KO = 100 with sf_KO = 0.25 -> scaled = 100 / 0.25 = 400!
        # And raw_WT = 400 with sf_WT = 4.0 -> scaled = 400 / 4.0 = 100!
        # So under TRUE division (raw / sf):
        #   KO scaled = 400 -> log2(401) = 8.647458...
        #   WT scaled = 100 -> log2(101) = 6.658211...
        #   full_lfc = log2(401) - log2(101) = +1.989247 (POSITIVE!)
        # Under INVERTED multiplication (raw * sf):
        #   KO scaled = 25 -> log2(26) = 4.700440
        #   WT scaled = 1600 -> log2(1601) = 10.644757
        #   full_lfc = -5.944317 (NEGATIVE! Opposite sign & 8 log2 units off!)
        samples = ["KO_1", "KO_2", "KO_3", "WT_1", "WT_2", "WT_3"]
        counts = pd.DataFrame(
            {
                "KO_1": [100, 200, 50, 0],
                "KO_2": [104, 196, 52, 0],
                "KO_3": [96, 204, 48, 0],
                "WT_1": [400, 200, 200, 0],
                "WT_2": [416, 192, 208, 0],
                "WT_3": [384, 208, 192, 0],
            },
            index=["P_up_after_sf", "P_strong_up", "P_moderate", "P_allzero"],
        )
        meta = pd.DataFrame(
            {"condition": ["KO", "KO", "KO", "WT", "WT", "WT"]},
            index=samples,
        )
        sf_dict = {
            "KO_1": 0.25, "KO_2": 0.25, "KO_3": 0.25,
            "WT_1": 4.0, "WT_2": 4.0, "WT_3": 4.0,
        }

        res = compute_counts_lobo_jackknife(
            counts, meta, "condition", "KO", "WT", size_factors=sf_dict
        )
        assert res["status"] == "measured"
        assert res["normalization"] == "spikein_or_external_size_factor_log2"

        # Independent closed-form log2(raw / sf + 1)
        sf_arr = np.array([sf_dict[s] for s in samples])
        expected_log = np.log2(counts.to_numpy(dtype=float) / sf_arr[None, :] + 1.0)
        for idx, peak in enumerate(counts.index):
            exp_full_lfc = float(expected_log[idx, :3].mean() - expected_log[idx, 3:].mean())
            assert res["per_feature"][peak]["full_lfc"] == pytest.approx(exp_full_lfc, abs=1e-6)

        assert res["per_feature"]["P_up_after_sf"]["full_lfc"] > 1.9
        assert res["per_feature"]["P_up_after_sf"]["single_unit_driven"] is False
        assert res["per_feature"]["P_allzero"]["full_lfc"] == 0.0

    def test_cpm_normalization_matches_closed_form_without_size_factors(self):
        samples = ["KO_1", "KO_2", "KO_3", "WT_1", "WT_2", "WT_3"]
        counts = pd.DataFrame(
            {
                "KO_1": [800, 100, 100],
                "KO_2": [820, 90, 90],
                "KO_3": [780, 110, 110],
                "WT_1": [200, 400, 400],
                "WT_2": [210, 395, 395],
                "WT_3": [190, 405, 405],
            },
            index=["G_up", "G_dn1", "G_dn2"],
        )
        meta = pd.DataFrame(
            {"condition": ["KO", "KO", "KO", "WT", "WT", "WT"]},
            index=samples,
        )
        res = compute_counts_lobo_jackknife(counts, meta, "condition", "KO", "WT")
        assert res["normalization"] == "log2_cpm_plus_1"
        raw = counts.to_numpy(dtype=float)
        cpm = raw / raw.sum(axis=0, keepdims=True) * 1e6
        log_cpm = np.log2(cpm + 1.0)
        for i, g in enumerate(counts.index):
            exp_lfc = float(log_cpm[i, :3].mean() - log_cpm[i, 3:].mean())
            assert res["per_feature"][g]["full_lfc"] == pytest.approx(exp_lfc, abs=1e-6)

    def test_lobo_truncation_unmoderated_paired_partial_and_collapse_branches(self):
        from scipy import stats as sp_stats
        from aria.utils.lobo_jackknife import _paired_fold_stats, _unpaired_fold_stats

        # 1. Exact _unpaired_fold_stats & _paired_fold_stats (moderate=False & True)
        rng = np.random.default_rng(11)
        mat_small = rng.normal(5.0, 0.3, size=(6, 6))
        mat_small[0, :3] += 3.0  # strong signal
        mat_small[1, :] = 2.0    # zero variance -> p=1.0, padj=1.0
        lfc_u, p_u, q_u = _unpaired_fold_stats(mat_small, [0, 1, 2], [3, 4, 5], moderate=False)
        _, exp_p = sp_stats.ttest_ind(mat_small[:, :3], mat_small[:, 3:], axis=1, equal_var=False)
        exp_p = np.where(np.isfinite(exp_p), exp_p, 1.0)
        np.testing.assert_allclose(p_u, exp_p, rtol=1e-12, atol=1e-12)
        assert p_u[1] == 1.0 and q_u[1] == 1.0

        # Paired unmoderated zero-variance nonzero LFC -> p=0.0; zero-variance zero LFC -> p=1.0
        diff_mat = np.array([
            [2.0, 2.0, 2.0, 2.0],   # zero var, lfc=2 -> inf t -> p=0.0
            [0.0, 0.0, 0.0, 0.0],   # zero var, lfc=0 -> nan t -> p=1.0
            [1.0, 1.5, 2.0, 2.5],   # positive var
        ])
        lfc_p, p_p, q_p = _paired_fold_stats(diff_mat, [0, 1, 2, 3], moderate=False)
        assert p_p[0] == 0.0 and p_p[1] == 1.0
        _, exp_t_p = sp_stats.ttest_1samp(diff_mat[2], 0.0)
        assert p_p[2] == pytest.approx(exp_t_p, rel=1e-12, abs=1e-12)

        # Moderated fallback when all variances == 0 (s0_sq defaults to 1e-4)
        _, p_mod_zero, _ = _unpaired_fold_stats(np.zeros((4, 6)), [0, 1, 2], [3, 4, 5], moderate=True)
        np.testing.assert_allclose(p_mod_zero, np.ones(4))
        _, p_pmod_zero, _ = _paired_fold_stats(np.zeros((4, 4)), [0, 1, 2, 3], moderate=True)
        np.testing.assert_allclose(p_pmod_zero, np.ones(4))

        # 2. target_features deduplication, max_reported_features truncation, and primary_significant_features gating
        mat25 = rng.normal(4.0, 0.25, size=(25, 8))
        mat25[0, :4] += 3.0                      # rock-solid stable hit
        mat25[1, :4] = [4.0, 4.0, 4.0, 10.0]     # single-unit outlier (severe effect collapse!)
        mat25[1, 4:] = [4.0, 4.0, 4.0, 4.0]
        mat25[2, 0] = np.nan                     # invalid_stats -> single_unit_driven=True, None JSON values
        feats25 = [f"F{i}" for i in range(25)]
        units8 = [f"U{i}" for i in range(8)]
        u2g8 = {u: ("KO" if i < 4 else "WT") for i, u in enumerate(units8)}
        res_trunc = compute_lobo_jackknife(
            feats25, mat25, units8, u2g8,
            numerator="KO", denominator="WT",
            target_features=["F0", "F0", "F1", "F2", "F3", "MISSING"],
            primary_significant_features={"F0", "F1"},
            max_reported_features=3,
        )
        assert res_trunc["variance_moderation"] == "empirical_bayes_smyth2004"
        assert res_trunc["n_features_evaluated"] == 4  # deduped & MISSING skipped
        assert res_trunc["per_feature_truncated"] is True
        assert list(res_trunc["per_feature"].keys()) == ["F0", "F1", "F2"]
        assert res_trunc["per_feature"]["F0"]["single_unit_driven"] is False
        assert res_trunc["per_feature"]["F1"]["single_unit_driven"] is True
        assert res_trunc["per_feature"]["F2"]["single_unit_driven"] is True
        assert res_trunc["per_feature"]["F2"]["full_lfc"] is None

        # 3. Partially paired design: 2 paired donors + 2 unpaired per arm -> falls back to unpaired!
        u_part = ["d1_k", "d1_w", "d2_k", "d2_w", "u3_k", "u3_w"]
        g_part = {u: ("KO" if u.endswith("_k") else "WT") for u in u_part}
        d_part = {"d1_k": "D1", "d1_w": "D1", "d2_k": "D2", "d2_w": "D2", "u3_k": "D3", "u3_w": "D4"}
        res_part = compute_lobo_jackknife(
            [f"F{i}" for i in range(5)], rng.normal(5.0, 0.2, size=(5, 6)),
            u_part, g_part, numerator="KO", denominator="WT", unit_to_donor=d_part,
        )
        assert res_part["status"] == "measured"
        assert res_part["design_mode"] == "unpaired_biological_replicate"
        assert res_part["variance_moderation"] == "exact_welch_or_paired_t"

        # 4. compute_counts_lobo_jackknife with invalid size_factors falls back cleanly to CPM,
        #    and paired_unit_col wires paired_donor mode!
        samples6 = ["d1_k", "d2_k", "d3_k", "d1_w", "d2_w", "d3_w"]
        counts6 = pd.DataFrame(
            rng.poisson(100, size=(6, 6)),
            index=[f"G{i}" for i in range(6)],
            columns=samples6,
        )
        meta6 = pd.DataFrame(
            {"condition": ["KO"] * 3 + ["WT"] * 3, "donor": ["D1", "D2", "D3"] * 2},
            index=samples6,
        )
        res_bad_sf = compute_counts_lobo_jackknife(
            counts6, meta6, "condition", "KO", "WT",
            size_factors={"d1_k": -1.0, "d2_k": 1.0, "d3_k": 1.0, "d1_w": 1.0, "d2_w": 1.0, "d3_w": 1.0},
            paired_unit_col="donor",
        )
        assert res_bad_sf["normalization"] == "log2_cpm_plus_1"
        assert res_bad_sf["design_mode"] == "paired_donor"
        assert res_bad_sf["n_jackknife_folds"] == 3

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "D2-BUG-5: compute_counts_lobo_jackknife calls np.asarray(size_factors) "
            "on a pandas Series, stripping its sample index and misaligning size factors "
            "when the Series index order differs from counts_df.columns"
        ),
    )
    def test_reordered_pandas_series_size_factors_aligned_by_sample_name(self):
        samples = ["KO_1", "KO_2", "KO_3", "WT_1", "WT_2", "WT_3"]
        counts = pd.DataFrame(
            {
                "KO_1": [100, 50],
                "KO_2": [104, 52],
                "KO_3": [96, 48],
                "WT_1": [400, 200],
                "WT_2": [416, 208],
                "WT_3": [384, 192],
            },
            index=["P1", "P2"],
        )
        meta = pd.DataFrame(
            {"condition": ["KO", "KO", "KO", "WT", "WT", "WT"]},
            index=samples,
        )
        sf_dict = {
            "KO_1": 0.25, "KO_2": 0.25, "KO_3": 0.25,
            "WT_1": 4.0, "WT_2": 4.0, "WT_3": 4.0,
        }
        # Same sample->sf mapping, but Series index is reversed!
        sf_series_reordered = pd.Series(sf_dict).loc[samples[::-1]]

        res_dict = compute_counts_lobo_jackknife(
            counts, meta, "condition", "KO", "WT", size_factors=sf_dict
        )
        res_series = compute_counts_lobo_jackknife(
            counts, meta, "condition", "KO", "WT", size_factors=sf_series_reordered
        )
        assert res_series["per_feature"]["P1"]["full_lfc"] == pytest.approx(
            res_dict["per_feature"]["P1"]["full_lfc"], abs=1e-6
        )
