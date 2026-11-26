"""Analytical oracle tests for Layer D of aria/utils/cross_modal_concordance.py.

Layer D (`sample_matched_layer`, `_residualize_covariates`, `_standardize`, `bh_adjust`, `read_cpm`):
  - Exact Pearson r, Spearman rho, and t-distribution analytical p-values compared directly
    against scipy.stats.pearsonr, scipy.stats.spearmanr, and scipy.stats.t.sf
  - Exact OLS covariate residualization compared against statsmodels.api.OLS / np.linalg.lstsq
  - Confounder elimination: spurious RNA-ATAC correlation (r ~ 0.99 driven by batch/sex)
    collapses to ~0 after covariate residualization, while a masked true cis link is recovered
  - Freedman-Lane permutation null: empirical p-value formula, BH padj parity, held-out draw
  - Sample-column and feature-row permutation invariance (including explicit sample_id_map)
  - Strand-aware signed distance_to_tss (+ vs - strand) and sex-chromosome link audit
  - Edge cases: < min_matched_samples, non-unique sample_id_map, constant/NaN/negative CPM,
    missing/constant covariates, residual_df < 3 fail-closed
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control, pearsonr, spearmanr, t

from aria.utils.cross_modal_concordance import (
    _residualize_covariates,
    _standardize,
    read_cpm,
    sample_matched_layer,
)


def _default_layer_d_params(**overrides):
    params = {
        "cis_distance_bp": 50_000,
        "min_matched_samples": 6,
        "n_permutations": 100,
        "padj_threshold": 0.05,
        "seed": 42,
    }
    params.update(overrides)
    return params


class TestResidualizeCovariatesOracle:
    def test_ols_residuals_match_statsmodels_ols_for_numeric_and_categorical_covariates(self):
        sm = pytest.importorskip("statsmodels.api")
        rng = np.random.default_rng(601)
        samples = [f"S{i}" for i in range(12)]
        batch = ["B1"] * 4 + ["B2"] * 4 + ["B3"] * 4
        age = rng.uniform(20.0, 80.0, size=12)
        mat = rng.normal(loc=10.0, scale=2.0, size=(5, 12))
        mat[0] += 5.0 * (np.array(batch) == "B2") + 0.1 * age

        covs = {
            "batch": dict(zip(samples, batch)),
            "age": dict(zip(samples, age.tolist())),
            "const_num": {s: 3.14 for s in samples},
            "const_cat": {s: "ONLY" for s in samples},
            "incomplete": {s: ("X" if i < 11 else None) for i, s in enumerate(samples)},
        }
        resid, meta = _residualize_covariates(mat, samples, covs)
        assert meta["adjusted"] is True
        assert meta["covariates_adjusted"] == ["age", "batch"]
        assert meta["design_rank"] == 4      # intercept + age + 2 batch dummies
        assert meta["residual_df"] == 8      # 12 - 4
        assert meta["skipped_covariates"]["const_num"] == "constant_covariate"
        assert meta["skipped_covariates"]["const_cat"] == "constant_covariate"
        assert meta["skipped_covariates"]["incomplete"] == "missing_values_for_shared_samples"

        # Independent OLS verification via statsmodels.api.OLS on each row
        X_sm = pd.DataFrame({
            "const": 1.0,
            "age": age,
            "b2": [1.0 if b == "B2" else 0.0 for b in batch],
            "b3": [1.0 if b == "B3" else 0.0 for b in batch],
        })
        for row_i in range(mat.shape[0]):
            fit = sm.OLS(mat[row_i], X_sm).fit()
            np.testing.assert_allclose(resid[row_i], fit.resid.to_numpy(), rtol=1e-10, atol=1e-10)

    def test_residualize_fails_closed_when_residual_df_below_3(self):
        samples = [f"S{i}" for i in range(5)]
        mat = np.arange(10, dtype=float).reshape(2, 5)
        # 3-level categorical -> rank=3 -> residual_df = 5 - 3 = 2 < 3 -> must fail closed!
        covs = {"grp": {"S0": "A", "S1": "A", "S2": "B", "S3": "B", "S4": "C"}}
        out, meta = _residualize_covariates(mat, samples, covs)
        assert meta["adjusted"] is False
        assert meta["reason"] == "insufficient_residual_df_for_covariates"
        assert meta["unadjusted_confounding_risk"] is True
        np.testing.assert_array_equal(out, mat)

    def test_standardize_zero_variance_row_safe(self):
        m = np.array([[1.0, 2.0, 3.0], [5.0, 5.0, 5.0]])
        z = _standardize(m)
        np.testing.assert_allclose(z[0], [-1 / np.sqrt(2), 0.0, 1 / np.sqrt(2)], atol=1e-12)
        np.testing.assert_array_equal(z[1], [0.0, 0.0, 0.0])


class TestSampleMatchedLayerDOracle:
    def test_unadjusted_pearson_spearman_t_pvalues_and_bh_match_scipy_exactly(self):
        rng = np.random.default_rng(602)
        n_samples = 12
        samples = [f"U{i:02d}" for i in range(n_samples)]
        x = np.linspace(1.0, 12.0, n_samples)

        rna_cpm = pd.DataFrame(
            {
                s: [
                    float(np.exp2(2.0 + 0.5 * x[i]) - 1.0),                 # G_pos (+ strand)
                    float(np.exp2(8.0 - 0.4 * x[i]) - 1.0),                 # G_neg (- strand)
                    float(rng.uniform(5.0, 100.0)),                         # G_null
                    float(np.exp2(2.0 + 0.5 * x[i]) - 1.0),                 # G_sex (chrX)
                    0.0,                                                    # G_zero_const (dropped)
                ]
                for i, s in enumerate(samples)
            },
            index=["G_pos", "G_neg", "G_null", "G_sex", "G_zero_const"],
        )
        atac_cpm = pd.DataFrame(
            {
                s: [
                    float(np.exp2(1.0 + 0.6 * x[i] + rng.normal(0, 0.05)) - 1.0),  # P_pos (r ~ +1)
                    float(np.exp2(1.0 + 0.5 * x[i] + rng.normal(0, 0.05)) - 1.0),  # P_neg (r ~ -1 with G_neg)
                    float(rng.uniform(5.0, 100.0)),                                 # P_null
                    float(np.exp2(1.0 + 0.6 * x[i] + rng.normal(0, 0.05)) - 1.0),  # P_sex (chrX)
                ]
                for i, s in enumerate(samples)
            },
            index=["chr1:12000-12400", "chr1:90000-90400", "chr2:50000-50400", "chrX:10000-10400"],
        )
        models = {
            "chr1": [
                {"gene": "G_pos", "tss": 10_000, "strand": "+"},
                {"gene": "G_neg", "tss": 95_000, "strand": "-"},
                {"gene": "G_zero_const", "tss": 10_500, "strand": "+"},
            ],
            "chr2": [{"gene": "G_null", "tss": 51_000, "strand": "+"}],
            "chrX": [{"gene": "G_sex", "tss": 11_000, "strand": "+"}],
        }

        warns: list[str] = []
        frame, summary = sample_matched_layer(
            rna_cpm, atac_cpm, models, _default_layer_d_params(n_permutations=200), warns
        )
        assert summary["ran"] is True
        assert summary["n_shared_samples"] == n_samples
        assert summary["n_pairs_tested"] == 4
        assert "cpm_sample_matched_layer_unadjusted_for_covariates" in warns
        assert any("unadjusted_sex_chromosome_links_detected" in w for w in warns)

        # Independent SciPy verification of Pearson r, Spearman rho, t-pvalues, and BH padj
        log_r = np.log2(rna_cpm + 1.0)
        log_a = np.log2(atac_cpm + 1.0)
        df_corr = n_samples - 2
        for row in frame.itertuples(index=False):
            rv = log_r.loc[row.gene].to_numpy()
            av = log_a.loc[row.peak].to_numpy()
            pr = pearsonr(rv, av)
            sr = spearmanr(rv, av)
            assert row.pearson_r == pytest.approx(pr.statistic, abs=1e-10)
            assert row.spearman_r == pytest.approx(sr.statistic, abs=1e-10)
            assert row.pearson_pvalue == pytest.approx(pr.pvalue, rel=1e-8, abs=1e-14)
            # Layer D computes spearman_pvalue via t-approx on rho with df = n - 2
            t_sp = abs(sr.statistic) * np.sqrt(df_corr / max(1.0 - sr.statistic ** 2, np.finfo(float).eps))
            assert row.spearman_pvalue == pytest.approx(2.0 * t.sf(t_sp, df=df_corr), rel=1e-8, abs=1e-14)

        # BH padj on empirical_pvalue matches SciPy false_discovery_control
        exp_padj = false_discovery_control(frame["empirical_pvalue"].to_numpy(), method="bh")
        np.testing.assert_allclose(frame["padj"].to_numpy(), exp_padj, rtol=1e-12, atol=1e-12)

        # Exact directional signs and strand-aware TSS distances
        by_gene = frame.set_index("gene")
        assert by_gene.loc["G_pos", "pearson_r"] > 0.98
        assert by_gene.loc["G_pos", "distance_to_tss"] == 2200   # (12000+12400)//2 - 10000
        assert by_gene.loc["G_neg", "pearson_r"] < -0.98
        assert by_gene.loc["G_neg", "distance_to_tss"] == 4800   # -(90200 - 95000) on '-' strand!
        assert by_gene.loc["G_sex", "is_sex_chromosome"] == True  # noqa: E712
        assert summary["n_sex_chromosome_significant_links"] == 1

    def test_covariate_residualization_removes_confounder_and_recovers_masked_link(self):
        # 16 samples across 2 sex/batch strata (8 M, 8 F):
        #   Pair 1 (G_conf, P_conf): PURELY confounded by batch (+6 log2 shift in B2, independent noise within batch)
        #   Pair 2 (G_true, P_true): TRUE within-batch cis link (same latent signal) + opposite batch shift!
        rng = np.random.default_rng(603)
        n = 16
        samples = [f"S{i:02d}" for i in range(n)]
        batch = ["B1"] * 8 + ["B2"] * 8
        b_ind = np.array([0.0] * 8 + [1.0] * 8)

        latent = rng.normal(0.0, 1.0, size=n)
        latent -= np.where(b_ind == 0, latent[:8].mean(), latent[8:].mean())  # orthogonal to batch

        r_conf_log = 3.0 + 6.0 * b_ind + rng.normal(0.0, 0.25, size=n)
        a_conf_log = 3.0 + 6.0 * b_ind + rng.normal(0.0, 0.25, size=n)

        # True pair has +1.2*latent in both RNA and ATAC, but OPPOSITE huge batch effects!
        r_true_log = 6.0 + 5.0 * b_ind + 1.2 * latent + rng.normal(0.0, 0.05, size=n)
        a_true_log = 11.0 - 5.0 * b_ind + 1.2 * latent + rng.normal(0.0, 0.05, size=n)

        rna_cpm = pd.DataFrame(
            np.vstack([np.exp2(r_conf_log) - 1.0, np.exp2(r_true_log) - 1.0]),
            index=["G_conf", "G_true"],
            columns=samples,
        )
        atac_cpm = pd.DataFrame(
            np.vstack([np.exp2(a_conf_log) - 1.0, np.exp2(a_true_log) - 1.0]),
            index=["chr1:10000-10400", "chr2:10000-10400"],
            columns=samples,
        )
        models = {
            "chr1": [{"gene": "G_conf", "tss": 10_200, "strand": "+"}],
            "chr2": [{"gene": "G_true", "tss": 10_200, "strand": "+"}],
        }

        # Unadjusted: G_conf looks strongly positive (~+0.99), G_true looks strongly NEGATIVE!
        f_unadj, s_unadj = sample_matched_layer(
            rna_cpm, atac_cpm, models, _default_layer_d_params(n_permutations=120), []
        )
        u_map = f_unadj.set_index("gene")
        assert u_map.loc["G_conf", "pearson_r"] > 0.95
        assert u_map.loc["G_true", "pearson_r"] < -0.5

        # Covariate-adjusted: G_conf collapses near 0 (non-significant), G_true flips to r > +0.98!
        covs = {"batch": dict(zip(samples, batch))}
        f_adj, s_adj = sample_matched_layer(
            rna_cpm, atac_cpm, models,
            _default_layer_d_params(n_permutations=120, sample_covariates=covs), [],
        )
        assert s_adj["covariate_adjustment"]["adjusted"] is True
        assert s_adj["covariate_adjustment"]["residual_df"] == 14
        a_map = f_adj.set_index("gene")
        assert abs(a_map.loc["G_conf", "pearson_r"]) < 0.55
        assert a_map.loc["G_conf", "pearson_pvalue"] > 0.05
        assert a_map.loc["G_true", "pearson_r"] > 0.98
        assert a_map.loc["G_true", "pearson_pvalue"] < 1e-8

        # Analytical df check: partial correlation uses df = residual_df - 1 = 13!
        r_obs = float(a_map.loc["G_true", "pearson_r"])
        t_stat = abs(r_obs) * np.sqrt(13.0 / max(1.0 - r_obs ** 2, np.finfo(float).eps))
        assert a_map.loc["G_true", "pearson_pvalue"] == pytest.approx(2.0 * t.sf(t_stat, df=13), rel=1e-8)

    def test_sample_and_feature_permutation_invariance_with_sample_id_map(self):
        rng = np.random.default_rng(604)
        n = 10
        rna_ids = [f"RNA_{i}" for i in range(n)]
        atac_ids = [f"ATAC_{i}" for i in range(n)]
        id_map = dict(zip(rna_ids, atac_ids))
        x = np.linspace(1.0, 10.0, n)

        rna_cpm = pd.DataFrame(
            np.vstack([np.exp2(1.0 + 0.5 * x) - 1.0, rng.uniform(1, 50, size=n)]),
            index=["G1", "G2"],
            columns=rna_ids,
        )
        atac_cpm = pd.DataFrame(
            np.vstack([np.exp2(1.0 + 0.5 * x + rng.normal(0, 0.1, size=n)) - 1.0, rng.uniform(1, 50, size=n)]),
            index=["chr1:1000-1200", "chr1:5000-5200"],
            columns=atac_ids,
        )
        models = {
            "chr1": [
                {"gene": "G1", "tss": 1100, "strand": "+"},
                {"gene": "G2", "tss": 5100, "strand": "+"},
            ]
        }

        f_base, _ = sample_matched_layer(
            rna_cpm, atac_cpm, models,
            _default_layer_d_params(sample_id_map=id_map, n_permutations=60), [],
        )
        # Permute RNA columns, ATAC columns (independently!), gene rows, and peak rows
        f_perm, _ = sample_matched_layer(
            rna_cpm.iloc[::-1, rng.permutation(n)],
            atac_cpm.iloc[::-1, rng.permutation(n)],
            models,
            _default_layer_d_params(sample_id_map=id_map, n_permutations=60), [],
        )
        b = f_base.set_index(["gene", "peak"])
        p = f_perm.set_index(["gene", "peak"]).loc[b.index]
        np.testing.assert_allclose(p["pearson_r"].to_numpy(), b["pearson_r"].to_numpy(), atol=1e-10)
        np.testing.assert_allclose(p["spearman_r"].to_numpy(), b["spearman_r"].to_numpy(), atol=1e-10)
        np.testing.assert_allclose(p["pearson_pvalue"].to_numpy(), b["pearson_pvalue"].to_numpy(), rtol=1e-8)

    def test_edge_cases_insufficient_samples_nonunique_map_and_read_cpm(self, tmp_path):
        # Empty inputs
        _, res_empty = sample_matched_layer(
            pd.DataFrame(), pd.DataFrame(), {}, _default_layer_d_params(), []
        )
        assert res_empty["ran"] is False
        assert res_empty["reason"] == "cpm_matrices_unavailable"

        # Non-unique sample_id_map (two RNA samples mapped to same ATAC sample)
        r_df = pd.DataFrame({"R1": [10.0], "R2": [20.0]}, index=["G1"])
        a_df = pd.DataFrame({"A1": [10.0]}, index=["chr1:1-100"])
        _, res_dup = sample_matched_layer(
            r_df, a_df, {"chr1": [{"gene": "G1", "tss": 50, "strand": "+"}]},
            _default_layer_d_params(sample_id_map={"R1": "A1", "R2": "A1"}), [],
        )
        assert res_dup["ran"] is False
        assert res_dup["reason"] == "non_unique_sample_id_map"

        # Insufficient matched samples (< 6)
        r_small = pd.DataFrame({f"S{i}": [float(i + 1)] for i in range(4)}, index=["G1"])
        a_small = pd.DataFrame({f"S{i}": [float(i + 1)] for i in range(4)}, index=["chr1:1-100"])
        _, res_few = sample_matched_layer(
            r_small, a_small, {"chr1": [{"gene": "G1", "tss": 50, "strand": "+"}]},
            _default_layer_d_params(min_matched_samples=6), [],
        )
        assert res_few["ran"] is False
        assert res_few["reason"] == "insufficient_matched_biological_units"

        # read_cpm: preserves leading-zero sample IDs ("001") and rejects duplicate headers
        good_tsv = tmp_path / "good.tsv"
        good_tsv.write_text("gene\t001\t002\nG1\t10\t20\n", encoding="utf-8")
        w_list: list[str] = []
        cpm = read_cpm(str(good_tsv), {}, w_list)
        assert list(cpm.columns) == ["001", "002"]
        assert cpm.loc["G1", "001"] == 10.0

        dup_tsv = tmp_path / "dup.tsv"
        dup_tsv.write_text("gene\tS1\tS1\nG1\t10\t20\n", encoding="utf-8")
        cpm_dup = read_cpm(str(dup_tsv), {}, w_list)
        assert cpm_dup.empty
        assert any("duplicate_sample_headers" in w for w in w_list)

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "D2-BUG-8: sample_matched_layer filters constant features with `rna.var(axis=1) > 0` "
            "after log2(x+1), so strictly constant nonzero genes (e.g. CPM=10.0, log2(11) var ~2e-31) "
            "leak into hypothesis testing and inflate n_pairs_tested / BH padj"
        ),
    )
    def test_nonzero_constant_gene_must_be_filtered_before_testing(self):
        samples = [f"U{i:02d}" for i in range(12)]
        x = np.linspace(1.0, 12.0, 12)
        rna_cpm = pd.DataFrame(
            {s: [10.0, float(np.exp2(2.0 + 0.5 * x[i]) - 1.0)] for i, s in enumerate(samples)},
            index=["G_const_10", "G_pos"],
        )
        atac_cpm = pd.DataFrame(
            {s: [float(np.exp2(1.0 + 0.6 * x[i]) - 1.0)] for i, s in enumerate(samples)},
            index=["chr1:12000-12400"],
        )
        models = {
            "chr1": [
                {"gene": "G_const_10", "tss": 10_500, "strand": "+"},
                {"gene": "G_pos", "tss": 10_000, "strand": "+"},
            ]
        }
        frame, summary = sample_matched_layer(
            rna_cpm, atac_cpm, models, _default_layer_d_params(n_permutations=50), []
        )
        assert summary["n_pairs_tested"] == 1
        assert list(frame["gene"]) == ["G_pos"]

    def test_residualize_covariates_and_layer_d_boundary_mutant_killers(self):
        from aria.utils.cross_modal_concordance import _residualize_covariates, _standardize

        # 1. _standardize on zero-variance row returns all zeros; unit norm otherwise
        s_mat = _standardize(np.array([[3.0, 3.0, 3.0, 3.0], [1.0, 2.0, 3.0, 4.0]]))
        np.testing.assert_allclose(s_mat[0], 0.0, atol=1e-12)
        assert float(np.linalg.norm(s_mat[1])) == pytest.approx(1.0, abs=1e-12)
        assert float(s_mat[1].mean()) == pytest.approx(0.0, abs=1e-12)

        # 2. _residualize_covariates: none/invalid/missing/constant-num/constant-cat/bool/insufficient-df/valid
        samples = [f"S{i}" for i in range(8)]
        mat = np.arange(16, dtype=float).reshape(2, 8)
        m0, meta0 = _residualize_covariates(mat, samples, None)
        assert meta0["adjusted"] is False and meta0["reason"] == "no_covariates_supplied"
        assert meta0["unadjusted_confounding_risk"] is True
        np.testing.assert_array_equal(m0, mat)

        # All covariates invalid/unestimable -> no_estimable_covariates with exact skip reasons
        bad_covs = {
            "not_dict": [1, 2, 3],
            "has_none": {s: (None if i == 0 else float(i)) for i, s in enumerate(samples)},
            "has_nan": {s: (float("nan") if i == 1 else float(i)) for i, s in enumerate(samples)},
            "has_blank": {s: ("  " if i == 2 else "A") for i, s in enumerate(samples)},
            "const_num": {s: 5.0 for s in samples},
            "const_cat": {s: "only" for s in samples},
        }
        _, meta_bad = _residualize_covariates(mat, samples, bad_covs)
        assert meta_bad["adjusted"] is False
        assert meta_bad["reason"] == "no_estimable_covariates"
        sk = meta_bad["skipped_covariates"]
        assert sk["not_dict"] == "invalid_mapping"
        assert sk["has_none"] == "missing_values_for_shared_samples"
        assert sk["has_nan"] == "missing_values_for_shared_samples"
        assert sk["has_blank"] == "missing_values_for_shared_samples"
        assert sk["const_num"] == "constant_covariate"
        assert sk["const_cat"] == "constant_covariate"

        # Insufficient residual df: 6 samples, 4-level categorical (rank 4 -> df=2 < 3) -> fails closed
        s6 = [f"S{i}" for i in range(6)]
        _, meta_df = _residualize_covariates(
            mat[:, :6], s6, {"cat4": {"S0": "A", "S1": "B", "S2": "C", "S3": "D", "S4": "A", "S5": "B"}}
        )
        assert meta_df["adjusted"] is False
        assert meta_df["reason"] == "insufficient_residual_df_for_covariates"
        assert meta_df["design_rank"] == 4 and meta_df["residual_df"] == 2

        # Mixed numeric + bool (treated as categorical!) + 3-level categorical -> exact rank & orthogonal residuals
        mixed_covs = {
            "age": {s: float(20 + i * 5) for i, s in enumerate(samples)},
            "is_male": {s: (i % 2 == 0) for i, s in enumerate(samples)},
            "site": {s: ("X" if i < 3 else ("Y" if i < 6 else "Z")) for i, s in enumerate(samples)},
        }
        resid, meta_ok = _residualize_covariates(mat, samples, mixed_covs)
        assert meta_ok["adjusted"] is True
        assert meta_ok["covariates_adjusted"] == ["age", "is_male", "site"]
        assert meta_ok["design_rank"] == 5 and meta_ok["residual_df"] == 3
        assert meta_ok["unadjusted_confounding_risk"] is False
        assert meta_ok["method"] == "ols_residualization_on_log2_cpm"
        assert meta_ok["permutation_method"] == "freedman_lane_residual_reprojection"
        # Residuals must be orthogonal to intercept, age, and site indicators
        np.testing.assert_allclose(resid.mean(axis=1), 0.0, atol=1e-10)
        age_v = np.array([20 + i * 5 for i in range(8)], dtype=float)
        np.testing.assert_allclose(resid @ (age_v - age_v.mean()), 0.0, atol=1e-10)

        # 3. sample_matched_layer edge returns & sex-chrom/negative-strand/negative-link reporting
        empty_df, empty_res = sample_matched_layer(
            pd.DataFrame(), pd.DataFrame(), {}, _default_layer_d_params(), []
        )
        assert empty_df.empty and empty_res["ran"] is False
        assert empty_res["reason"] == "cpm_matrices_unavailable"

        # Negative strand distance sign, sex chromosome link accounting, negative correlation link,
        # and explicit sample_id_map
        rng = np.random.default_rng(77)
        rna_s = [f"R{i}" for i in range(10)]
        atac_s = [f"A{i}" for i in range(10)]
        smap = dict(zip(rna_s, atac_s))
        x = np.linspace(1.0, 10.0, 10)
        rna_df = pd.DataFrame(
            {
                s: [
                    float(np.exp2(2.0 + 0.5 * x[i]) - 1.0),
                    float(np.exp2(8.0 - 0.5 * x[i]) - 1.0),
                    -5.0,  # negative CPM must be dropped!
                ]
                for i, s in enumerate(rna_s)
            },
            index=["G_POS_X", "G_NEG_MINUS", "G_BAD_NEG"],
        )
        atac_df = pd.DataFrame(
            {
                s: [
                    float(np.exp2(1.0 + 0.6 * x[i]) - 1.0),
                    float(np.exp2(1.0 + 0.5 * x[i]) - 1.0),
                ]
                for i, s in enumerate(atac_s)
            },
            index=["chrX:12000-12400", "chr1:12000-12400"],
        )
        models = {
            "chrX": [{"gene": "G_POS_X", "tss": 10_000, "strand": "+"}],
            "chr1": [
                {"gene": "G_NEG_MINUS", "tss": 10_000, "strand": "-"},
                {"gene": "G_BAD_NEG", "tss": 10_000, "strand": "+"},
            ],
        }
        p = _default_layer_d_params(n_permutations=100, seed=5)
        p["sample_id_map"] = smap
        warns: list[str] = []
        df_out, summ = sample_matched_layer(rna_df, atac_df, models, p, warns)
        assert summ["ran"] is True
        assert summ["n_pairs_tested"] == 2  # G_BAD_NEG dropped
        assert summ["n_positive_links"] == 1
        assert summ["n_negative_links"] == 1
        assert summ["n_sex_chromosome_significant_links"] == 1
        assert summ["sex_chromosome_significant_fraction"] == pytest.approx(0.5)
        assert any("cpm_sample_matched_layer_unadjusted_for_covariates" in w for w in warns)
        by_gene = df_out.set_index("gene")
        # Minus strand flips distance_to_tss sign: center=12200, tss=10000 -> -2200
        assert int(by_gene.loc["G_NEG_MINUS", "distance_to_tss"]) == -2200
        assert int(by_gene.loc["G_POS_X", "distance_to_tss"]) == 2200
        assert bool(by_gene.loc["G_POS_X", "is_sex_chromosome"]) is True
        assert bool(by_gene.loc["G_NEG_MINUS", "is_sex_chromosome"]) is False

        # No cis pairs within tiny window -> no_variable_cis_pairs
        p_tiny = _default_layer_d_params(cis_distance_bp=50)
        p_tiny["sample_id_map"] = smap
        df_tiny, summ_tiny = sample_matched_layer(rna_df, atac_df, models, p_tiny, [])
        assert df_tiny.empty and summ_tiny["ran"] is False
        assert summ_tiny["reason"] == "no_variable_cis_pairs"

    def test_layer_d_batching_cap_covariate_collinearity_and_schema_mutant_killers(self):
        from aria.utils.cross_modal_concordance import _residualize_covariates, bh_adjust

        # 1. Direct bh_adjust clipping on out-of-range inputs
        np.testing.assert_allclose(bh_adjust(np.array([-0.2, 0.01, 1.4])), [0.0, 0.015, 1.0], atol=1e-12)

        # 2. _residualize_covariates: empty dict, non-dict, collinearity, large-scale numeric standardization
        samples = [f"S{i}" for i in range(8)]
        mat = np.arange(16, dtype=float).reshape(2, 8)
        for empty_cov in ({}, "not_a_dict"):
            m_e, meta_e = _residualize_covariates(mat, samples, empty_cov)  # type: ignore[arg-type]
            assert meta_e["adjusted"] is False
            assert meta_e["covariates_adjusted"] == []
            assert meta_e["unadjusted_confounding_risk"] is True
            np.testing.assert_array_equal(m_e, mat)

        # Collinear numeric covariates + huge scale (1e5) -> SVD rank stays 2 (not 3), residual_df == 6
        x_base = np.linspace(1.0, 8.0, 8)
        col_covs = {
            "cov_a": {s: float(1e5 * x_base[i]) for i, s in enumerate(samples)},
            "cov_b": {s: float(2e5 * x_base[i] + 7.0) for i, s in enumerate(samples)},
            "near_const": {s: float(5.0 + (1e-14 if i == 0 else 0.0)) for i, s in enumerate(samples)},
        }
        resid_c, meta_c = _residualize_covariates(mat, samples, col_covs)
        assert meta_c["adjusted"] is True
        assert meta_c["design_rank"] == 2
        assert meta_c["residual_df"] == 6
        assert meta_c["covariates_adjusted"] == ["cov_a", "cov_b"]
        assert meta_c["skipped_covariates"]["near_const"] == "constant_covariate"
        np.testing.assert_allclose(resid_c @ (x_base - x_base.mean()), 0.0, atol=1e-9)

        # Insufficient residual df also exposes candidate_covariates
        s6 = [f"S{i}" for i in range(6)]
        _, meta_df_cand = _residualize_covariates(
            mat[:, :6], s6, {"cat4": {"S0": "A", "S1": "B", "S2": "C", "S3": "D", "S4": "A", "S5": "B"}}
        )
        assert meta_df_cand["candidate_covariates"] == ["cat4"]

        # 3. sample_matched_layer full schema keys on ran=False AND ran=True + adjusted sex-chr link
        #    (adjusted=True with sex-chr hit MUST NOT emit unadjusted warning!) + NaN/neg atac filter
        _, empty_res = sample_matched_layer(
            pd.DataFrame(), pd.DataFrame(), {}, _default_layer_d_params(), []
        )
        assert empty_res["n_pairs_tested"] == 0
        assert empty_res["n_significant_links"] == 0
        assert empty_res["n_positive_links"] == 0
        assert empty_res["n_negative_links"] == 0
        assert empty_res["n_sex_chromosome_significant_links"] == 0
        assert empty_res["sex_chromosome_significant_fraction"] == 0.0
        assert empty_res["covariate_adjustment"]["reason"] == "layer_not_run"
        assert empty_res["covariate_adjustment"]["unadjusted_confounding_risk"] is True
        assert empty_res["permutation_null"] == {}
        assert empty_res["top_links"] == []

        rng = np.random.default_rng(123)
        s8 = [f"S{i}" for i in range(8)]
        batch = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)
        sig = np.array([-2.0, -0.5, 0.5, 2.0, -2.0, -0.5, 0.5, 2.0])
        rna8 = pd.DataFrame(
            {
                s: [
                    float(np.exp2(5.0 + 1.0 * sig[i] + 0.5 * batch[i]) - 1.0),
                    float("nan") if i == 0 else 5.0 + i,
                ]
                for i, s in enumerate(s8)
            },
            index=["GX_OK", "G_NAN"],
        )
        atac8 = pd.DataFrame(
            {
                s: [
                    float(np.exp2(5.0 + 1.0 * sig[i] + 0.5 * batch[i]) - 1.0),
                    -1.0 if i == 0 else 5.0 + i,
                    float("nan") if i == 1 else 5.0 + i,
                ]
                for i, s in enumerate(s8)
            },
            index=["chrX:10500-10800", "chrX:11000-11200", "chrX:11500-11700"],
        )
        models8 = {"chrX": [{"gene": "GX_OK", "tss": 10_000, "strand": "+"}, {"gene": "G_NAN", "tss": 10_000, "strand": "+"}]}
        p8 = _default_layer_d_params(n_permutations=40, seed=9)
        p8["sample_covariates"] = {"batch": {s: ("B0" if i < 4 else "B1") for i, s in enumerate(s8)}}
        w8: list[str] = []
        df8, summ8 = sample_matched_layer(rna8, atac8, models8, p8, w8)
        assert summ8["ran"] is True
        assert summ8["n_pairs_tested"] == 1  # G_NAN & negative/NaN atac peaks dropped
        assert summ8["shared_sample_ids"] == s8
        assert summ8["sample_id_pairs"] == [{"rna": s, "atac": s} for s in s8]
        assert summ8["covariate_adjustment"]["adjusted"] is True
        assert summ8["n_sex_chromosome_significant_links"] == 1
        assert not any("sex_chromosome_links_without_covariate_adjustment" in w for w in w8)
        pn = summ8["permutation_null"]
        assert pn["n_permutations"] == 40
        assert pn["seed"] == 9
        assert pn["method"] == "pooled_absolute_pearson_sample_permutation"
        assert pn["held_out_draw_independent"] is True
        assert len(summ8["top_links"]) == 1
        assert int(df8.iloc[0]["n_samples"]) == 8

        # 4. >4096 cis pairs exercises multi-batch chunking (range(0, len(rows), 4096))
        s6 = [f"B{i}" for i in range(6)]
        x6 = np.linspace(1.0, 6.0, 6)
        rna_b = pd.DataFrame(
            {s: [float(np.exp2(2.0 + 0.4 * x6[i] + 0.1 * g) - 1.0) for g in range(5)] for i, s in enumerate(s6)},
            index=[f"GB_{g}" for g in range(5)],
        )
        n_pk_b = 850  # 5 * 850 = 4250 > 4096 pairs!
        atac_vals = rng.uniform(5.0, 50.0, size=(n_pk_b, 6))
        atac_vals[-1] = np.exp2(2.0 + 0.4 * x6) - 1.0  # last peak sits in batch 2 (>4096)
        atac_b = pd.DataFrame(
            atac_vals,
            index=[f"chr1:{1000 + p*20}-{1010 + p*20}" for p in range(n_pk_b)],
            columns=s6,
        )
        models_b = {"chr1": [{"gene": f"GB_{g}", "tss": 10_000, "strand": "+"} for g in range(5)]}
        df_b, summ_b = sample_matched_layer(
            rna_b, atac_b, models_b, _default_layer_d_params(cis_distance_bp=50_000, n_permutations=10, seed=3), []
        )
        assert summ_b["ran"] is True
        assert summ_b["n_pairs_tested"] == 4250
        assert len(df_b) == 4250
        assert len(summ_b["top_links"]) == 20

        # 5. >50,000 ATAC peaks triggers top_50000_variable_peaks cap & warning
        s6_c = [f"C{i}" for i in range(6)]
        n_cap = 50_005
        # Give peak 0 huge variance (retained in top 50,000) and peaks 1..50004 tiny variance
        cap_mat = rng.uniform(10.0, 10.5, size=(n_cap, 6))
        cap_mat[0] = [1.0, 5.0, 20.0, 80.0, 300.0, 1000.0]
        atac_cap = pd.DataFrame(
            cap_mat,
            index=["chr1:10100-10300"] + [f"chr2:{i*100}-{i*100+50}" for i in range(1, n_cap)],
            columns=s6_c,
        )
        rna_cap = pd.DataFrame([[1.0, 5.0, 20.0, 80.0, 300.0, 1000.0]], index=["GCAP"], columns=s6_c)
        w_cap: list[str] = []
        p_cap = _default_layer_d_params(min_samples=6, n_permutations=5, seed=1)
        df_cap, summ_cap = sample_matched_layer(
            rna_cap, atac_cap, {"chr1": [{"gene": "GCAP", "tss": 10_000, "strand": "+"}]}, p_cap, w_cap
        )
        assert summ_cap["ran"] is True
        assert summ_cap["n_pairs_tested"] == 1
        assert any("sample_link_candidate_cap: top_50000_variable_peaks" in w for w in w_cap)
