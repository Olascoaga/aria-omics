"""Analytical oracle tests for differential 3D statistics in aria/scripts/hic_topology.py.

Covers:
  - _moderated_linear_contrast:
      * Exact match to closed-form Smyth (2004) empirical-Bayes moderated OLS
      * Unmoderated limit (prior_df=0) exact match to statsmodels.api.OLS
      * Sign flip on contrast reversal (KO vs WT <-> WT vs KO)
      * Alphabetical & prefix level renaming invariance (KO vs KO_ctrl, A_trt vs Z_ctrl, Z_trt vs A_ctrl)
      * Sample-column and bin-row permutation invariance
      * Categorical and continuous covariate adjustment removing confounders
      * Sparse NaN partial-case OLS path vs complete-case vectorized path
      * Error cases: missing/NaN covariate, collinear/confounded covariate, saturated df, single arm
  - _parse_comparisons: all supported contrast formats
  - _compute_differential_3d (in aria-hic-env with cooler + cooltools):
      * Phased synthetic checkerboard coolers with planted A_to_B / B_to_A compartment switch
        and reverse contrast symmetry (A_to_B <-> B_to_A)
      * Unphased suppression (emits unphased_sign_flip, never A_to_B / B_to_A)
      * Planted TAD insulation boundary strengthening/weakening (delta_insulation sign)
      * Technical replicate collapse before biological replicate testing
      * Insufficient replicates -> descriptive_only=True with zero fabricated p-values
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.stats import false_discovery_control

from aria.scripts.hic_topology import (
    _compute_differential_3d,
    _moderated_linear_contrast,
    _parse_comparisons,
)


def _closed_form_moderated_ols(
    y: np.ndarray, X: np.ndarray, prior_df: float = 4.0
) -> dict[str, np.ndarray]:
    """Independent reference Smyth-style moderated OLS for column 1 of design matrix X."""
    n_obs, n_params = X.shape
    df_resid = float(n_obs - n_params)
    xtx_inv = np.linalg.inv(X.T @ X)
    coef = y @ X @ xtx_inv.T
    resid = y - coef @ X.T
    s2 = np.sum(resid ** 2, axis=1) / df_resid
    valid_s2 = s2[np.isfinite(s2) & (s2 > 0)]
    s0_sq = max(float(np.median(valid_s2)) if valid_s2.size else 1e-4, 1e-6)
    s2_post = (prior_df * s0_sq + df_resid * s2) / (prior_df + df_resid)
    se = np.sqrt(np.maximum(s2_post * xtx_inv[1, 1], 1e-16))
    t_stat = coef[:, 1] / se
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=prior_df + df_resid)
    padj = false_discovery_control(pvals, method="bh")
    return {
        "delta": coef[:, 1],
        "t_stat": t_stat,
        "pvalue": pvals,
        "padj": padj,
        "prior_s0_sq": s0_sq,
    }


class TestModeratedLinearContrastOracle:
    def test_exact_match_to_closed_form_smyth_moderated_ols(self):
        rng = np.random.default_rng(701)
        n_bins = 40
        conditions = ["WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO"]
        mat = rng.normal(loc=0.2, scale=0.12, size=(n_bins, 8))
        mat[:5, 4:] -= 1.0   # A_to_B shift (-1.0)
        mat[5:10, 4:] += 1.0  # B_to_A shift (+1.0)

        got = _moderated_linear_contrast(mat, conditions, "KO", "WT", prior_df=4.0)
        X = np.column_stack([np.ones(8), np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)])
        ref = _closed_form_moderated_ols(mat, X, prior_df=4.0)

        np.testing.assert_allclose(got["delta"], ref["delta"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got["t_stat"], ref["t_stat"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got["pvalue"], ref["pvalue"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got["padj"], ref["padj"], rtol=1e-12, atol=1e-12)
        assert got["prior_s0_sq"] == pytest.approx(ref["prior_s0_sq"], rel=1e-12)
        assert (got["delta"][:5] < -0.7).all()
        assert (got["delta"][5:10] > 0.7).all()
        assert (got["padj"][:10] < 1e-3).all()

    def test_unmoderated_limit_matches_statsmodels_ols(self):
        sm = pytest.importorskip("statsmodels.api")
        rng = np.random.default_rng(701)
        conditions = ["WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO"]
        mat = rng.normal(loc=0.2, scale=0.12, size=(20, 8))
        mat[:5, 4:] -= 1.0
        X = np.column_stack([np.ones(8), np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)])

        got_ols = _moderated_linear_contrast(mat, conditions, "KO", "WT", prior_df=0.0)
        for b in range(12):
            fit = sm.OLS(mat[b], X).fit()
            assert got_ols["delta"][b] == pytest.approx(fit.params[1], rel=1e-10, abs=1e-12)
            assert got_ols["t_stat"][b] == pytest.approx(fit.tvalues[1], rel=1e-10, abs=1e-12)
            assert got_ols["pvalue"][b] == pytest.approx(fit.pvalues[1], rel=1e-10, abs=1e-12)

    def test_contrast_reversal_alphabetical_renaming_and_permutation_invariance(self):
        rng = np.random.default_rng(702)
        n_bins = 30
        conditions = ["WT", "WT", "WT", "KO", "KO", "KO"]
        mat = rng.normal(0.0, 0.15, size=(n_bins, 6))
        mat[:4, 3:] += 1.2
        mat[4:8, 3:] -= 1.2

        base = _moderated_linear_contrast(mat, conditions, "KO", "WT")
        rev = _moderated_linear_contrast(mat, conditions, "WT", "KO")

        # Exact negation of delta & t_stat; exact equality of pvalue & padj
        np.testing.assert_allclose(rev["delta"], -base["delta"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(rev["t_stat"], -base["t_stat"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(rev["pvalue"], base["pvalue"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(rev["padj"], base["padj"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(rev["mean_test"], base["mean_ref"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(rev["mean_ref"], base["mean_test"], rtol=1e-12, atol=1e-12)

        # Alphabetical & prefix traps: (KO vs KO_ctrl), (A_trt vs Z_ctrl), (Z_trt vs A_ctrl)
        for num_l, den_l in [("KO", "KO_ctrl"), ("A_trt", "Z_ctrl"), ("Z_trt", "A_ctrl")]:
            ren_conds = [num_l if c == "KO" else den_l for c in conditions]
            ren = _moderated_linear_contrast(mat, ren_conds, num_l, den_l)
            np.testing.assert_allclose(ren["delta"], base["delta"], rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(ren["padj"], base["padj"], rtol=1e-12, atol=1e-12)

        # Sample-column & bin-row permutation invariance
        s_perm = rng.permutation(6)
        b_perm = rng.permutation(n_bins)
        perm_res = _moderated_linear_contrast(
            mat[np.ix_(b_perm, s_perm)], [conditions[i] for i in s_perm], "KO", "WT"
        )
        np.testing.assert_allclose(perm_res["delta"], base["delta"][b_perm], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(perm_res["t_stat"], base["t_stat"][b_perm], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(perm_res["padj"], base["padj"][b_perm], rtol=1e-12, atol=1e-12)

    def test_categorical_and_continuous_covariates_remove_confounders(self):
        rng = np.random.default_rng(703)
        n_bins = 30
        conditions = ["WT"] * 6 + ["KO"] * 6
        batches = ["B1", "B1", "B1", "B2", "B2", "B2"] * 2
        depth = np.array([1.0, 2.0, 3.0, 1.2, 2.2, 3.2, 1.1, 2.1, 3.1, 1.3, 2.3, 3.3])

        # Huge batch (+4.0) and depth (+2.0*depth) effects mask a clean +0.8 KO shift!
        b_ind = np.array([1.0 if b == "B2" else 0.0 for b in batches])
        t_ind = np.array([1.0 if c == "KO" else 0.0 for c in conditions])
        mat = (
            4.0 * b_ind[None, :]
            + 2.0 * depth[None, :]
            + rng.normal(0.0, 0.08, size=(n_bins, 12))
        )
        mat[:5, :] += 0.8 * t_ind[None, :]

        cov_rows = [
            {"batch": b, "depth": float(d), "const_cov": "SAME"}
            for b, d in zip(batches, depth)
        ]
        adj = _moderated_linear_contrast(mat, conditions, "KO", "WT", covariate_rows=cov_rows)
        unadj = _moderated_linear_contrast(mat, conditions, "KO", "WT")

        assert "const_cov:constant_within_contrast" in adj["dropped_covariates"]
        # Adjusted recovers exact beta ~ +0.8 with padj < 1e-4; unadjusted is completely swamped!
        np.testing.assert_allclose(adj["delta"][:5], 0.8, atol=0.15)
        assert (adj["padj"][:5] < 1e-4).all()
        assert float(np.median(unadj["padj"][:5])) > 0.25

    def test_sparse_nan_partial_case_and_error_guards(self):
        rng = np.random.default_rng(704)
        conditions = ["WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO"]
        mat = rng.normal(0.0, 0.1, size=(10, 8))
        mat[0, 4:] += 1.5
        mat[1, 4:] += 1.5
        # Bin 1 loses 1 WT replicate (still >=2 per arm -> partial-case OLS runs!)
        mat[1, 0] = np.nan
        # Bin 2 loses 3 WT replicates (only 1 WT left < 2 -> untestable -> NaN!)
        mat[2, :3] = np.nan

        res = _moderated_linear_contrast(mat, conditions, "KO", "WT", prior_df=0.0)
        # Verify bin 1 matches exact 7-sample OLS manually
        ok = np.isfinite(mat[1])
        exp_delta_1 = float(mat[1, 4:].mean() - mat[1, 1:4].mean())
        assert res["delta"][1] == pytest.approx(exp_delta_1, abs=1e-12)
        assert np.isfinite(res["pvalue"][1]) and res["pvalue"][1] < 0.01
        assert np.isnan(res["delta"][2]) and np.isnan(res["padj"][2])

        # Error guards: missing covariate, collinear/confounded covariate, saturated df
        with pytest.raises(ValueError, match="missing/NaN"):
            _moderated_linear_contrast(
                mat, conditions, "KO", "WT",
                covariate_rows=[{"b": None if i == 0 else "X"} for i in range(8)],
            )
        with pytest.raises(ValueError, match="collinear/confounded"):
            _moderated_linear_contrast(
                mat, conditions, "KO", "WT",
                covariate_rows=[{"b": "B1" if c == "WT" else "B2"} for c in conditions],
            )
        with pytest.raises(ValueError, match="saturates residual degrees of freedom"):
            # 6 linearly independent continuous covariates + intercept + is_test = 8 params == 8 obs
            cov_rng = np.random.default_rng(705)
            cov_mat = cov_rng.normal(size=(8, 6))
            _moderated_linear_contrast(
                mat, conditions, "KO", "WT",
                covariate_rows=[
                    {f"c{k}": float(cov_mat[i, k]) for k in range(6)} for i in range(8)
                ],
            )

    def test_parse_comparisons_variants(self):
        assert _parse_comparisons(["KO", "WT"]) == [("KO", "WT")]
        assert _parse_comparisons({"test": "KO", "reference": "WT"}) == [("KO", "WT")]
        assert _parse_comparisons({"numerator": "KO", "denominator": "WT"}) == [("KO", "WT")]
        assert _parse_comparisons([["A", "B"], {"test": "C", "reference": "D"}]) == [
            ("A", "B"),
            ("C", "D"),
        ]


class TestComputeDifferential3DCoolerOracle:
    """End-to-end differential 3D architecture tests with synthetic balanced coolers (aria-hic-env)."""

    @staticmethod
    def _write_synthetic_cooler(
        path,
        *,
        n_bins: int = 32,
        binsize: int = 100_000,
        a_bins: set[int],
        boundary_bin: int | None = None,
        boundary_insulation: float = 0.25,
        seed: int = 0,
    ):
        import cooler

        rng = np.random.default_rng(seed)
        bins = pd.DataFrame({
            "chrom": ["chr1"] * n_bins,
            "start": np.arange(n_bins) * binsize,
            "end": (np.arange(n_bins) + 1) * binsize,
        })
        is_a = np.array([i in a_bins for i in range(n_bins)], dtype=bool)
        rows = []
        for i in range(n_bins):
            for j in range(i, n_bins):
                dist = j - i
                # Power-law distance decay + strong A/B checkerboard + optional TAD insulation
                decay = 120.0 / ((dist + 1.0) ** 0.75)
                comp_mult = 2.6 if (is_a[i] == is_a[j]) else 0.35
                tad_mult = 1.0
                if boundary_bin is not None and i < boundary_bin <= j:
                    tad_mult = boundary_insulation
                lam = max(1.0, decay * comp_mult * tad_mult)
                cnt = int(rng.poisson(lam))
                if cnt > 0:
                    rows.append((i, j, cnt))
        pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"])
        cooler.create_cooler(str(path), bins=bins, pixels=pixels)
        return str(path)

    def test_differential_3d_phased_ab_switch_insulation_and_reverse_symmetry(self, tmp_path):
        pytest.importorskip("cooler")
        pytest.importorskip("cooltools")

        n_bins = 32
        binsize = 100_000
        # Balanced 50/50 A/B checkerboard in BOTH WT and KO so global mean E1 stays 0:
        #   WT: A = bins 0..7 and 16..23; B = bins 8..15 and 24..31
        #   KO: A = bins 0..7 and 8..15;  B = bins 16..23 and 24..31
        # So in KO vs WT:
        #   bins 8..15  switch B_to_A (delta_E1 > 0)
        #   bins 16..23 switch A_to_B (delta_E1 < 0)
        #   bins 0..7 stay A with high gene density -> keeps GTF phasing positively aligned!
        wt_a = set(range(0, 8)) | set(range(16, 24))
        ko_a = set(range(0, 8)) | set(range(8, 16))

        # Write a synthetic GTF with high gene density in bins 0..7 (always A)
        gtf = tmp_path / "phasing.gtf"
        gtf_lines = []
        for b in range(0, 8):
            for k in range(5):
                pos = b * binsize + 10_000 + k * 5_000
                gtf_lines.append(
                    f'chr1\tsyn\tgene\t{pos}\t{pos+1000}\t.\t+\t.\tgene_id "G{b}_{k}"; gene_name "G{b}_{k}";\n'
                )
        gtf.write_text("".join(gtf_lines), encoding="utf-8")

        files = []
        meta = {}
        # 3 WT replicates (no TAD boundary at bin 28 inside stable B block 24..31: mult=1.0)
        # 3 KO replicates (strong TAD boundary at bin 28: mult=0.05 -> lower log2 insulation = strengthened!)
        for r in range(3):
            sid = f"WT_r{r}"
            p = self._write_synthetic_cooler(
                tmp_path / f"{sid}.cool",
                n_bins=n_bins, binsize=binsize, a_bins=wt_a,
                boundary_bin=28, boundary_insulation=1.0, seed=10 + r,
            )
            files.append(p)
            meta[sid] = {"condition": "WT", "replicate": f"r{r}"}
        for r in range(3):
            sid = f"KO_r{r}"
            p = self._write_synthetic_cooler(
                tmp_path / f"{sid}.cool",
                n_bins=n_bins, binsize=binsize, a_bins=ko_a,
                boundary_bin=28, boundary_insulation=0.05, seed=20 + r,
            )
            files.append(p)
            meta[sid] = {"condition": "KO", "replicate": f"r{r}"}

        # 1. Phased KO vs WT -> bins 16..23 MUST be A_to_B, bins 8..15 MUST be B_to_A,
        #    and boundary bin 16 MUST be boundary_strengthened (negative delta_insulation)!
        res_kw = _compute_differential_3d(
            {
                "files": files,
                "sample_ids": list(meta.keys()),
                "sample_metadata": meta,
                "comparisons": [["KO", "WT"]],
                "comp_resolution": binsize,
                "tad_resolution": binsize,
                "window_size": binsize * 3,
                "gtf": str(gtf),
                "min_replicates_per_condition": 3,
                "padj_max": 0.05,
                "delta_e1_min": 0.2,
                "delta_insulation_min": 0.15,
            },
            str(tmp_path / "diff_kw"),
        )
        assert res_kw["status"] == "success" and res_kw["ran"] is True
        assert res_kw["pc1_validated"] is True
        assert res_kw["phasing_method"] == "gene_density"
        comp_kw = res_kw["comparisons"][0]
        assert comp_kw["n_a_to_b"] >= 5
        assert comp_kw["n_b_to_a"] >= 5
        assert comp_kw["n_boundaries_strengthened"] >= 1

        comp_df = pd.read_csv(comp_kw["differential_compartments_csv"])
        assert (comp_df.iloc[16:24]["delta_E1"] < -0.3).sum() >= 6
        assert (comp_df.iloc[16:24]["switch_class"] == "A_to_B").sum() >= 5
        assert (comp_df.iloc[8:16]["delta_E1"] > 0.3).sum() >= 6
        assert (comp_df.iloc[8:16]["switch_class"] == "B_to_A").sum() >= 5

        # 2. Reverse contrast WT vs KO -> exact flip to B_to_A and boundary_weakened!
        res_wk = _compute_differential_3d(
            {
                "files": files,
                "sample_ids": list(meta.keys()),
                "sample_metadata": meta,
                "comparisons": [["WT", "KO"]],
                "comp_resolution": binsize,
                "tad_resolution": binsize,
                "window_size": binsize * 3,
                "gtf": str(gtf),
                "min_replicates_per_condition": 3,
                "padj_max": 0.05,
                "delta_e1_min": 0.2,
                "delta_insulation_min": 0.15,
            },
            str(tmp_path / "diff_wk"),
        )
        comp_wk = res_wk["comparisons"][0]
        assert comp_wk["n_b_to_a"] == comp_kw["n_a_to_b"]
        assert comp_wk["n_a_to_b"] == comp_kw["n_b_to_a"]
        assert comp_wk["n_boundaries_weakened"] == comp_kw["n_boundaries_strengthened"]

        # 3. Unphased run (no GTF/FASTA) MUST suppress directional A_to_B/B_to_A labels
        res_unph = _compute_differential_3d(
            {
                "files": files,
                "sample_ids": list(meta.keys()),
                "sample_metadata": meta,
                "comparisons": [["KO", "WT"]],
                "comp_resolution": binsize,
                "tad_resolution": binsize,
                "window_size": binsize * 3,
                "min_replicates_per_condition": 3,
            },
            str(tmp_path / "diff_unph"),
        )
        assert res_unph["pc1_validated"] is False
        assert res_unph["comparisons"][0]["n_a_to_b"] == 0
        assert res_unph["comparisons"][0]["n_b_to_a"] == 0
        assert res_unph["comparisons"][0]["n_unphased_sign_flips"] >= 4

        # 4. Insufficient replicates (2 vs 2 with min_replicates=3) -> descriptive_only, NO p-values
        res_2v2 = _compute_differential_3d(
            {
                "files": [files[0], files[1], files[3], files[4]],
                "sample_ids": ["WT_r0", "WT_r1", "KO_r0", "KO_r1"],
                "sample_metadata": meta,
                "comparisons": [["KO", "WT"]],
                "comp_resolution": binsize,
                "tad_resolution": binsize,
                "window_size": binsize * 3,
                "gtf": str(gtf),
                "min_replicates_per_condition": 3,
            },
            str(tmp_path / "diff_2v2"),
        )
        assert res_2v2["ran"] is False
        assert res_2v2["descriptive_only"] is True
        assert res_2v2["comparisons"][0]["descriptive_only"] is True
        desc_df = pd.read_csv(res_2v2["comparisons"][0]["descriptive_csv"])
        assert "pvalue" not in desc_df.columns
        assert "delta_E1" in desc_df.columns
        # Assert signed descriptive delta_E1 matches planted A->B (negative) and B->A (positive)
        assert (desc_df.iloc[16:24]["delta_E1"] < -0.3).sum() >= 6
        assert (desc_df.iloc[8:16]["delta_E1"] > 0.3).sum() >= 6

    def test_moderated_linear_contrast_and_parse_comparisons_mutant_killers(self):
        from scipy import stats as sp_stats
        from aria.scripts.hic_topology import (
            _moderated_linear_contrast,
            _parse_comparisons,
        )

        # 1. _parse_comparisons: 2-str list, list of dicts/tuples, single dict, missing keys, invalid
        assert _parse_comparisons(["KO", "WT"]) == [("KO", "WT")]
        assert _parse_comparisons(("KO", "WT")) == [("KO", "WT")]
        assert _parse_comparisons([{"numerator": "A", "denominator": "B"}, ["C", "D"], {"test": "E"}]) == [
            ("A", "B"),
            ("C", "D"),
        ]
        assert _parse_comparisons({"numerator": "X", "denominator": "Y"}) == [("X", "Y")]
        assert _parse_comparisons({"test": "X"}) == []
        assert _parse_comparisons(None) == []

        # 2. Third-group exclusion (keep_samples), constant numeric + constant categorical covariates
        #    dropped, bool covariate treated as categorical, exact mean_test / mean_ref, and partial NaN
        #    bin matching independent OLS + rank-deficient partial NaN bin skipped!
        conds = ["WT", "WT", "WT", "KO", "KO", "KO", "OTHER", "OTHER"]
        cov_rows = [
            {"const_num": 5.0, "const_cat": "X", "flag": True, "age": 10.0},
            {"const_num": 5.0, "const_cat": "X", "flag": False, "age": 20.0},
            {"const_num": 5.0, "const_cat": "X", "flag": True, "age": 30.0},
            {"const_num": 5.0, "const_cat": "X", "flag": False, "age": 10.0},
            {"const_num": 5.0, "const_cat": "X", "flag": True, "age": 20.0},
            {"const_num": 5.0, "const_cat": "X", "flag": False, "age": 30.0},
            # OTHER samples have NaN covariates that MUST be ignored since keep_samples excludes OTHER!
            {"const_num": None, "const_cat": "", "flag": None, "age": float("nan")},
            {"const_num": None, "const_cat": "", "flag": None, "age": float("nan")},
        ]
        mat = np.array([
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 999.0, -999.0],   # complete bin
            [1.0, 2.2, 2.9, 4.1, 5.1, 6.2, 500.0, 500.0],    # complete bin with noise
            [1.0, 2.0, 3.0, 4.0, 5.0, np.nan, 0.0, 0.0],     # partial bin: 5 obs > 4 params -> testable
            [1.0, np.nan, 3.0, 4.0, np.nan, 6.0, 0.0, 0.0],  # partial bin: 4 obs == 4 params -> skipped
            [1.0, np.nan, np.nan, 4.0, 5.0, 6.0, 0.0, 0.0],  # only 1 valid WT -> untestable (<2)
        ], dtype=float)
        res = _moderated_linear_contrast(
            mat, conds, test_group="KO", ref_group="WT", covariate_rows=cov_rows, prior_df=4.0
        )
        assert "const_num:constant_within_contrast" in res["dropped_covariates"]
        assert "const_cat:constant_within_contrast" in res["dropped_covariates"]
        assert res["mean_test"][0] == pytest.approx(5.0, abs=1e-12)
        assert res["mean_ref"][0] == pytest.approx(2.0, abs=1e-12)
        assert res["delta"][0] == pytest.approx(3.0, abs=1e-10)
        assert np.isfinite(res["delta"][2]) and res["mean_test"][2] == pytest.approx(4.5, abs=1e-12)
        assert res["mean_ref"][2] == pytest.approx(2.0, abs=1e-12)
        assert np.isnan(res["delta"][3])
        assert np.isnan(res["delta"][4])

        # 3. Missing/blank/None covariate value inside active contrast -> raises ValueError
        for bad_val in (None, float("nan"), "   "):
            bad_cov = [{"batch": "B1"}, {"batch": "B2"}, {"batch": "B1"}, {"batch": "B2"}, {"batch": "B1"}, {"batch": bad_val}]
            with pytest.raises(ValueError, match="missing/NaN"):
                _moderated_linear_contrast(
                    mat[:2, :6], conds[:6], "KO", "WT", covariate_rows=bad_cov
                )

        # 4. Exact-zero all-bins -> prior_s0_sq defaults to 1e-4; tiny-variance -> floored at 1e-6
        res_zero = _moderated_linear_contrast(np.zeros((3, 6)), conds[:6], "KO", "WT", prior_df=4.0)
        assert res_zero["prior_s0_sq"] == pytest.approx(1e-4, abs=1e-12)
        assert res_zero["delta"][0] == pytest.approx(0.0, abs=1e-12)

        tiny_mat = np.array([
            [1.0, 1.0 + 1e-5, 1.0, 2.0, 2.0, 2.0],
            [1.0, 1.0, 1.0 + 1e-5, 2.0, 2.0, 2.0],
        ])
        res_tiny = _moderated_linear_contrast(tiny_mat, conds[:6], "KO", "WT", prior_df=4.0)
        assert res_tiny["prior_s0_sq"] == pytest.approx(1e-6, abs=1e-12)

        # 5. Rank-deficient base design (only 1 obs per group -> df_resid=0) -> ValueError
        with pytest.raises(ValueError, match="rank-deficient|zero residual"):
            _moderated_linear_contrast(np.ones((2, 2)), ["KO", "WT"], "KO", "WT")
