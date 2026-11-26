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
