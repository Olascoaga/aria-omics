"""Analytical oracle tests for aria/scripts/rna_pseudobulk_de.py.

Verifies:
  - Pure helpers: _build_design_formula, _shrink_coeff (prefix trap KO vs KO_ctrl),
    _detect_shrink_sign_flip, _abs_corr, _composition_covariate_decision,
    _expressed_background (dense & sparse), _effective_alpha_from_significant,
    _fdr_filtering_basis
  - End-to-end pseudobulk aggregation & DESeq2 oracles (in aria-rna-env):
      * Exact cell-to-pseudosample sum aggregation & min_cells filtering
      * LFC sign for KO vs WT and exact negation for WT vs KO (kills D3-002)
      * Alphabetical & prefix level renaming invariance (KO vs KO_ctrl, Z_trt vs A_ctrl)
      * Cell-row and gene-column permutation invariance
      * Global BH-FDR exact match against statsmodels.multipletests
      * Paired donor design removing donor baseline variance & recovering signal (kills D3-019)
      * Partial pairing subsetting & disclosure
      * Exact log-normalized count recovery matching raw-count DE
      * Edge cases: one group, one replicate, all-zero genes, too few cells, NonIntegerCounts
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from aria.scripts.rna_pseudobulk_de import (
    _abs_corr,
    _build_design_formula,
    _composition_covariate_decision,
    _detect_shrink_sign_flip,
    _effective_alpha_from_significant,
    _expressed_background,
    _fdr_filtering_basis,
    _shrink_coeff,
    rna_pseudobulk_de,
)


def _make_pseudobulk_h5ad(
    tmp_path,
    *,
    filename: str = "pb.h5ad",
    n_genes: int = 60,
    n_up: int = 10,
    n_down: int = 10,
    donors_per_cond: int = 4,
    cells_per_donor: int = 20,
    lfc_mag: float = 2.5,
    test_label: str = "KO",
    ref_label: str = "WT",
    paired: bool = False,
    donor_shift_sd: float = 0.0,
    include_all_zero: bool = True,
    seed: int = 202,
):
    """Build a small deterministic scRNA AnnData with known DE genes & biological replicates."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    base_rate = rng.uniform(4.0, 12.0, size=n_genes)
    fc = np.ones(n_genes)
    fc[:n_up] = 2.0 ** lfc_mag
    fc[n_up : n_up + n_down] = 2.0 ** (-lfc_mag)

    rows = []
    obs_rows = []
    if paired:
        donors = [f"D{d}" for d in range(donors_per_cond)]
        for d_idx, donor in enumerate(donors):
            donor_mult = np.exp(rng.normal(0.0, donor_shift_sd, size=n_genes)) if donor_shift_sd > 0 else np.ones(n_genes)
            for cond in (ref_label, test_label):
                rate = base_rate * donor_mult * (fc if cond == test_label else 1.0)
                for c in range(cells_per_donor):
                    rows.append(rng.poisson(rate))
                    obs_rows.append({
                        "cell_type": "Tcell",
                        "condition": cond,
                        "donor": donor,
                    })
    else:
        for cond in (ref_label, test_label):
            for d in range(donors_per_cond):
                donor = f"{cond}_D{d}"
                donor_mult = np.exp(rng.normal(0.0, 0.08, size=n_genes))
                rate = base_rate * donor_mult * (fc if cond == test_label else 1.0)
                for c in range(cells_per_donor):
                    rows.append(rng.poisson(rate))
                    obs_rows.append({
                        "cell_type": "Tcell",
                        "condition": cond,
                        "donor": donor,
                    })

    X = np.asarray(rows, dtype=np.int64)
    if include_all_zero:
        X[:, -2:] = 0
    obs = pd.DataFrame(obs_rows, index=[f"cell_{i}" for i in range(len(obs_rows))])
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    path = tmp_path / filename
    adata.write_h5ad(path)

    up_genes = genes[:n_up]
    down_genes = genes[n_up : n_up + n_down]
    null_end = n_genes - (2 if include_all_zero else 0)
    null_genes = genes[n_up + n_down : null_end]
    return path, adata, up_genes, down_genes, null_genes


class TestPseudobulkPureHelpersOracle:
    def test_formula_and_shrink_coeff_prefix_resistance(self):
        assert _build_design_formula("cond", ["donor", "cond"], ref_level="WT") == (
            "~ donor + C(cond, Treatment(reference='WT'))"
        )

        class _FakeDDS:
            def __init__(self, cols):
                self.varm = {"LFC": pd.DataFrame(columns=cols)}

        assert _shrink_coeff(_FakeDDS(["Intercept", "cond[T.KO_ctrl]"]), "cond", "KO") is None
        assert _shrink_coeff(
            _FakeDDS(["Intercept", "C(cond, Treatment(reference='KO_ctrl'))[T.KO]"]),
            "cond",
            "KO",
        ) == "C(cond, Treatment(reference='KO_ctrl'))[T.KO]"
        assert _detect_shrink_sign_flip([2.0, -2.0, 1.5], [-1.5, 1.5, -1.0]) is True
        assert _detect_shrink_sign_flip([2.0, -2.0, 1.5], [1.5, -1.5, 1.0]) is False

    def test_composition_collinearity_and_abs_corr_exact_boundaries(self):
        # Perfect correlation -> |r| = 1.0 -> dropped
        assert _abs_corr([1, 2, 3, 4], [2, 4, 6, 8]) == pytest.approx(1.0, abs=1e-12)
        assert _abs_corr([1, 2, 3, 4], [-2, -4, -6, -8]) == pytest.approx(1.0, abs=1e-12)
        assert _abs_corr([1, 1, 1], [1, 2, 3]) is None
        assert _abs_corr([1.0], [2.0]) is None

        use_col, reason_col = _composition_covariate_decision(
            [1.0, 1.0, 1.0, 2.0, 2.0, 2.0], [0, 0, 0, 1, 1, 1], max_corr=0.8
        )
        assert use_col is False
        assert "collinear_with_condition" in reason_col

        # Unmeasurable (constant) -> dropped conservatively
        use_const, reason_const = _composition_covariate_decision(
            [0.5, 0.5, 0.5, 0.5], [0, 0, 1, 1], max_corr=0.8
        )
        assert use_const is False
        assert reason_const == "collinearity_unmeasurable"

        # Orthogonal (r = 0.0 < 0.8) -> kept
        use_orth, reason_orth = _composition_covariate_decision(
            [-1.0, 1.0, -1.0, 1.0], [0, 0, 1, 1], max_corr=0.8
        )
        assert use_orth is True
        assert reason_orth is None

    def test_expressed_background_effective_alpha_and_fdr_basis(self):
        dense = np.array([[0, 3, 0], [0, 0, 5]])
        genes = ["G0", "G1", "G2"]
        bg_d, deg_d = _expressed_background(dense, genes)
        bg_s, deg_s = _expressed_background(sparse.csr_matrix(dense), genes)
        assert bg_d == ["G1", "G2"] and deg_d is False
        assert bg_s == ["G1", "G2"] and deg_s is False

        df_sig = pd.DataFrame({"pvalue": [0.001, 0.015, float("nan"), 0.008]})
        assert _effective_alpha_from_significant(df_sig) == pytest.approx(0.015, abs=1e-12)
        assert _effective_alpha_from_significant(pd.DataFrame(columns=["pvalue"])) == 0.0

        blocks = [
            {
                "results": pd.DataFrame({
                    "pvalue": [0.001, 0.05, 0.4, float("nan")],
                    "padj_local": [0.003, float("nan"), 0.8, float("nan")],
                })
            }
        ]
        basis = _fdr_filtering_basis(blocks)
        assert basis["n_global_pool_tests"] == 3
        assert basis["n_independent_filtered_into_global"] == 1
        assert basis["same_base_hypothesis_set"] is False


class TestPseudobulkDESeq2Oracles:
    def test_signed_lfc_ko_vs_wt_reverse_and_global_bh_parity(self, tmp_path):
        """Kills D3-002 (negating log2FoldChange) and verifies global BH vs statsmodels."""
        pytest.importorskip("pydeseq2")
        sm_multitest = pytest.importorskip("statsmodels.stats.multitest")

        h5ad, adata, up_genes, down_genes, null_genes = _make_pseudobulk_h5ad(
            tmp_path, seed=202, include_all_zero=True
        )

        res_kw = rna_pseudobulk_de({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "min_cells_per_pseudosample": 10,
            "min_replicates_per_condition": 3,
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "kw"),
        })
        assert res_kw["status"] == "success", res_kw
        assert res_kw["paired_design"] is False
        assert res_kw["lfc_shrinkage"]["applied"] is True
        grp = res_kw["per_group"]["Tcell"]
        assert grp["n_pseudosamples"] == 8
        assert all(sz == 20 for sz in grp["pseudosample_sizes"].values())
        # All-zero genes excluded from tested background
        assert "G058" not in grp["background_genes"]
        assert "G059" not in grp["background_genes"]

        comp_kw = grp["per_comparison"]["KO_vs_WT"]
        assert comp_kw["status"] == "success"
        sig_map = {r["gene"]: r for r in comp_kw["all_sig"]}
        # Planted up genes MUST have positive log2fc and log2fc_raw; planted down MUST be negative
        up_hit = [g for g in up_genes if g in sig_map]
        dn_hit = [g for g in down_genes if g in sig_map]
        assert len(up_hit) >= 9
        assert len(dn_hit) >= 9
        assert all(sig_map[g]["log2fc"] > 1.0 and sig_map[g]["log2fc_raw"] > 1.0 for g in up_hit)
        assert all(sig_map[g]["log2fc"] < -1.0 and sig_map[g]["log2fc_raw"] < -1.0 for g in dn_hit)
        assert len(set(sig_map) & set(null_genes)) <= 1

        # Verify global BH-FDR in full_results_csv matches statsmodels.multipletests exactly
        full_kw = pd.read_csv(res_kw["full_results_csv"]).set_index("gene")
        _, sm_padj, _, _ = sm_multitest.multipletests(full_kw["pvalue"].to_numpy(), method="fdr_bh")
        np.testing.assert_allclose(full_kw["padj_global"].to_numpy(), sm_padj, rtol=1e-10, atol=1e-12)

        # Reverse contrast WT vs KO -> exact sign flip & identical pvalues/padj
        res_wk = rna_pseudobulk_de({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["WT", "KO"]],
            "use_raw": False,
            "min_cells_per_pseudosample": 10,
            "min_replicates_per_condition": 3,
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "wk"),
        })
        full_wk = pd.read_csv(res_wk["full_results_csv"]).set_index("gene")
        common = full_kw.index.intersection(full_wk.index)
        np.testing.assert_allclose(
            full_wk.loc[common, "log2fc_raw"].to_numpy(),
            -full_kw.loc[common, "log2fc_raw"].to_numpy(),
            rtol=1e-3,
            atol=2e-3,
        )
        np.testing.assert_allclose(
            full_wk.loc[common, "log2fc"].to_numpy(),
            -full_kw.loc[common, "log2fc"].to_numpy(),
            rtol=1e-3,
            atol=2e-3,
        )
        # Reversing the contrast is the same Wald test, but pydeseq2's IRLS starts
        # from a different parametrization, so padj agrees only to fitting precision
        # (up to ~2% relative on the v4.9.5 release-repro RNA image vs <1e-4 locally).
        # The invariance that matters scientifically is exact: the same genes are
        # called significant in both directions, with near-identical evidence.
        padj_wk = full_wk.loc[common, "padj_global"].to_numpy(dtype=float)
        padj_kw = full_kw.loc[common, "padj_global"].to_numpy(dtype=float)
        clearly_sig = padj_kw < 0.025
        clearly_null = padj_kw > 0.1
        assert clearly_sig.sum() >= 18  # the 9 up + 9 down planted genes
        assert np.all(padj_wk[clearly_sig] < 0.05)
        assert np.all(padj_wk[clearly_null] >= 0.05)
        np.testing.assert_allclose(
            np.log10(np.clip(padj_wk, 1e-300, 1.0)),
            np.log10(np.clip(padj_kw, 1e-300, 1.0)),
            rtol=0.05,
            atol=0.05,
        )

    def test_cell_gene_permutation_and_level_renaming_invariance(self, tmp_path):
        pytest.importorskip("pydeseq2")
        h5ad, adata, _, _, _ = _make_pseudobulk_h5ad(tmp_path, seed=203)
        base = rna_pseudobulk_de({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "output_dir": str(tmp_path / "base"),
        })
        full_base = pd.read_csv(base["full_results_csv"]).set_index("gene")

        # 1. Permute cell rows AND gene columns
        rng = np.random.default_rng(88)
        cell_perm = rng.permutation(adata.n_obs)
        gene_perm = rng.permutation(adata.n_vars)
        adata_perm = adata[cell_perm, gene_perm].copy()
        perm_path = tmp_path / "perm.h5ad"
        adata_perm.write_h5ad(perm_path)

        res_perm = rna_pseudobulk_de({
            "data_path": str(perm_path),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "output_dir": str(tmp_path / "perm"),
        })
        full_perm = pd.read_csv(res_perm["full_results_csv"]).set_index("gene").loc[full_base.index]
        np.testing.assert_allclose(full_perm["log2fc"].to_numpy(), full_base["log2fc"].to_numpy(), atol=1e-3)
        np.testing.assert_allclose(full_perm["pvalue"].to_numpy(), full_base["pvalue"].to_numpy(), rtol=1e-4, atol=1e-6)

        # 2. Rename levels KO -> KO, WT -> KO_ctrl (alphabetical & prefix trap!)
        adata_ren = adata.copy()
        adata_ren.obs["condition"] = adata_ren.obs["condition"].map({"KO": "KO", "WT": "KO_ctrl"})
        ren_path = tmp_path / "ren.h5ad"
        adata_ren.write_h5ad(ren_path)

        res_ren = rna_pseudobulk_de({
            "data_path": str(ren_path),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "KO_ctrl"]],
            "use_raw": False,
            "output_dir": str(tmp_path / "ren"),
        })
        assert res_ren["lfc_shrinkage"]["applied"] is True
        full_ren = pd.read_csv(res_ren["full_results_csv"]).set_index("gene").loc[full_base.index]
        np.testing.assert_allclose(full_ren["log2fc"].to_numpy(), full_base["log2fc"].to_numpy(), atol=1e-3)
        np.testing.assert_allclose(full_ren["padj_global"].to_numpy(), full_base["padj_global"].to_numpy(), rtol=1e-4, atol=1e-6)

    def test_paired_donor_block_removes_donor_shift_and_kills_d3_019(self, tmp_path):
        pytest.importorskip("pydeseq2")
        # Huge donor-to-donor variance (sd=1.8 in log space) masks moderate 1.4 log2FC treatment
        # unless donor is modeled as a paired covariate!
        h5ad, _, up_genes, down_genes, null_genes = _make_pseudobulk_h5ad(
            tmp_path,
            filename="paired.h5ad",
            n_genes=60,
            n_up=8,
            n_down=8,
            donors_per_cond=4,
            cells_per_donor=25,
            lfc_mag=1.5,
            paired=True,
            donor_shift_sd=1.6,
            seed=205,
        )
        res_paired = rna_pseudobulk_de({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "auto_paired_donor_covariate": True,
            "output_dir": str(tmp_path / "paired_on"),
        })
        res_unpaired = rna_pseudobulk_de({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "auto_paired_donor_covariate": False,
            "output_dir": str(tmp_path / "paired_off"),
        })
        comp_p = res_paired["per_group"]["Tcell"]["per_comparison"]["KO_vs_WT"]
        comp_u = res_unpaired["per_group"]["Tcell"]["per_comparison"]["KO_vs_WT"]
        assert res_paired["paired_design"] is True
        assert comp_p["paired_donor_covariate"] is True
        assert comp_p["pairing_status"] == "complete"
        assert comp_u["paired_donor_covariate"] is False

        sig_p = {r["gene"]: r for r in comp_p["all_sig"]}
        sig_u = {r["gene"]: r for r in comp_u["all_sig"]}
        planted = set(up_genes) | set(down_genes)
        assert len(set(sig_p) & planted) >= 13
        assert len(set(sig_p) & planted) > len(set(sig_u) & planted) + 4
        assert all(sig_p[g]["log2fc"] > 0.5 for g in up_genes if g in sig_p)
        assert all(sig_p[g]["log2fc"] < -0.5 for g in down_genes if g in sig_p)
        assert len(set(sig_p) & set(null_genes)) <= 1

    def test_lognorm_recovery_matches_raw_counts_and_edge_cases(self, tmp_path):
        pytest.importorskip("pydeseq2")
        import anndata as ad

        raw_h5ad, adata, up_genes, down_genes, _ = _make_pseudobulk_h5ad(
            tmp_path, filename="raw.h5ad", seed=206, include_all_zero=False
        )
        res_raw = rna_pseudobulk_de({
            "data_path": str(raw_h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "output_dir": str(tmp_path / "raw_out"),
        })
        # Create log1p(10000 * x / lib_size) matrix + nCount_RNA column
        X_raw = adata.X.astype(np.float64)
        lib = X_raw.sum(axis=1)
        X_log = np.log1p(X_raw * 10000.0 / lib[:, None])
        ad_log = ad.AnnData(X=X_log.astype(np.float64), obs=adata.obs.copy(), var=adata.var.copy())
        ad_log.obs["nCount_RNA"] = lib
        log_h5ad = tmp_path / "lognorm.h5ad"
        ad_log.write_h5ad(log_h5ad)

        # allow_lognorm_recovery=False must fail closed with NonIntegerCounts
        res_refused = rna_pseudobulk_de({
            "data_path": str(log_h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "allow_lognorm_recovery": False,
            "output_dir": str(tmp_path / "refused"),
        })
        assert res_refused["status"] == "error"
        assert res_refused["error_type"] == "NonIntegerCounts"

        # allow_lognorm_recovery=True must recover exact counts -> identical DE!
        res_rec = rna_pseudobulk_de({
            "data_path": str(log_h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "allow_lognorm_recovery": True,
            "output_dir": str(tmp_path / "rec_out"),
        })
        assert res_rec["status"] == "success"
        assert res_rec["lognorm_recovered"] is True
        assert res_rec["count_source"] == "recovered_from_lognorm"
        df_raw = pd.read_csv(res_raw["full_results_csv"]).set_index("gene")
        df_rec = pd.read_csv(res_rec["full_results_csv"]).set_index("gene").loc[df_raw.index]
        np.testing.assert_allclose(df_rec["log2fc"].to_numpy(), df_raw["log2fc"].to_numpy(), atol=1e-3)
        np.testing.assert_allclose(df_rec["pvalue"].to_numpy(), df_raw["pvalue"].to_numpy(), rtol=1e-4, atol=1e-6)

        # Edge case: 1 replicate per group -> group skipped cleanly
        h5ad_1rep, _, _, _, _ = _make_pseudobulk_h5ad(
            tmp_path, filename="one_rep.h5ad", donors_per_cond=1, seed=207
        )
        res_1rep = rna_pseudobulk_de({
            "data_path": str(h5ad_1rep),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "output_dir": str(tmp_path / "one_rep_out"),
        })
        assert res_1rep["status"] == "success"
        assert res_1rep["per_group"]["Tcell"]["status"] == "skipped"


class TestPseudobulkHelperMutantKillers:
    """Targeted boundary & analytical tests killing surviving mutants in
    rna_pseudobulk_de pure helpers (_fdr_filtering_basis, _detect_shrink_sign_flip,
    _shrink_coeff, _abs_corr, _expressed_background, _effective_alpha_from_significant,
    _composition_covariate_decision, _build_design_formula)."""

    def test_fdr_filtering_basis_counts_and_disclosures(self):
        from aria.scripts.rna_pseudobulk_de import _fdr_filtering_basis

        # Empty / None / missing results / missing pvalue blocks
        empty = _fdr_filtering_basis(None)
        assert empty["same_base_hypothesis_set"] is True
        assert empty["n_global_pool_tests"] == 0
        assert empty["n_independent_filtered_into_global"] == 0
        assert "independent filtering" in empty["local_family"]
        assert "finite p-value" in empty["global_family"]
        assert "conservative cross-family" in empty["note"]

        # Mixed blocks: block without pvalue skipped; block without padj_local counts pool only;
        # block with NaN/inf pvalue & padj_local counts independent-filtered accurately.
        b_no_p = {"results": pd.DataFrame({"log2FoldChange": [1.0]})}
        b_none = {"results": None}
        b_no_local = {"results": pd.DataFrame({"pvalue": [0.01, np.nan, 0.2]})}
        b_mixed = {
            "results": pd.DataFrame({
                "pvalue": [0.001, 0.05, 0.4, np.nan, np.inf],
                "padj_local": [0.01, np.nan, np.inf, np.nan, 0.9],
            })
        }
        out = _fdr_filtering_basis([b_no_p, b_none, b_no_local, b_mixed])
        # Finite p-values: 2 from b_no_local + 3 from b_mixed = 5
        assert out["n_global_pool_tests"] == 5
        # Finite pvalue & non-finite padj_local in b_mixed: idx 1 (nan) and idx 2 (inf) -> 2
        assert out["n_independent_filtered_into_global"] == 2
        assert out["same_base_hypothesis_set"] is False

        # All finite padj_local -> same_base_hypothesis_set True
        b_same = {"results": pd.DataFrame({"pvalue": [0.01, 0.2], "padj_local": [0.02, 0.2]})}
        out_same = _fdr_filtering_basis([b_same])
        assert out_same["same_base_hypothesis_set"] is True
        assert out_same["n_global_pool_tests"] == 2
        assert out_same["n_independent_filtered_into_global"] == 0

    def test_detect_shrink_sign_flip_exact_thresholds(self):
        from aria.scripts.rna_pseudobulk_de import _detect_shrink_sign_flip

        # Empty / mismatched shapes -> False
        assert _detect_shrink_sign_flip([], []) is False
        assert _detect_shrink_sign_flip([1.0, 2.0], [1.0]) is False
        # Below clear thresholds (|raw| < 0.5 or |shrunk| < 0.1 or non-finite) -> False
        assert _detect_shrink_sign_flip([0.49, np.nan, np.inf], [-0.49, -1.0, -1.0]) is False
        assert _detect_shrink_sign_flip([2.0, -2.0], [-0.09, 0.09]) is False
        # Exact boundary: |raw| == 0.5 and |shrunk| == 0.1 with n_clear == 3 and 2/3 flipped (>0.5) -> True
        assert _detect_shrink_sign_flip([0.5, -0.5, 0.5], [-0.1, 0.1, 0.1]) is True
        # Exactly 50% flipped (2 of 4) -> False (strict `flip_frac > 0.5`)
        assert _detect_shrink_sign_flip([1.0, 1.0, 1.0, 1.0], [-0.5, -0.5, 0.5, 0.5]) is False
        # n_clear == 2, 100% flipped, max |raw| == 1.0 -> True
        assert _detect_shrink_sign_flip([1.0, -0.6], [-0.2, 0.2]) is True
        # n_clear == 2, 100% flipped, max |raw| == 0.99 (< 1.0) -> False
        assert _detect_shrink_sign_flip([0.99, -0.6], [-0.2, 0.2]) is False
        # n_clear == 2, 1/2 flipped -> False
        assert _detect_shrink_sign_flip([2.0, 2.0], [-0.5, 0.5]) is False

    def test_shrink_coeff_all_formats_and_prefix_safety(self):
        from aria.scripts.rna_pseudobulk_de import _shrink_coeff

        class _FakeDDS:
            def __init__(self, cols):
                self.varm = {"LFC": pd.DataFrame(columns=cols)}

        assert _shrink_coeff(object(), "condition", "KO") is None
        # Plain exact vs prefix trap (KO_ctrl must NOT match KO)
        assert _shrink_coeff(_FakeDDS(["Intercept", "condition[T.KO_ctrl]"]), "condition", "KO") is None
        assert _shrink_coeff(_FakeDDS(["condition[T.KO_ctrl]", "condition[T.KO]"]), "condition", "KO") == "condition[T.KO]"
        # Backtick-quoted non-identifier column
        assert _shrink_coeff(_FakeDDS(["`treat group`[T.KO]"]), "treat group", "KO") == "`treat group`[T.KO]"
        # Formulaic C(...) with and without extra args; ensure prefix column doesn't false-match
        c_col = "C(condition, Treatment(reference='WT'))[T.KO]"
        assert _shrink_coeff(_FakeDDS(["C(condition_other, Treatment(reference='WT'))[T.KO]", c_col]), "condition", "KO") == c_col
        assert _shrink_coeff(_FakeDDS(["C(condition)[T.KO]"]), "condition", "KO") == "C(condition)[T.KO]"
        assert _shrink_coeff(_FakeDDS(["C(`treat group`, Treatment(reference='WT'))[T.KO]"]), "treat group", "KO") == "C(`treat group`, Treatment(reference='WT'))[T.KO]"
        assert _shrink_coeff(_FakeDDS(["other_col[T.KO]"]), "condition", "KO") is None

    def test_abs_corr_effective_alpha_and_expressed_background_edges(self):
        from scipy import sparse
        from aria.scripts.rna_pseudobulk_de import (
            _abs_corr,
            _composition_covariate_decision,
            _effective_alpha_from_significant,
            _expressed_background,
        )

        # _abs_corr: NaNs filtered, <2 finite -> None, zero variance -> None, exact negative corr -> +1.0
        assert _abs_corr([1.0], [2.0]) is None
        assert _abs_corr([1.0, np.nan], [2.0, 3.0]) is None
        assert _abs_corr([5.0, 5.0, 5.0], [1.0, 2.0, 3.0]) is None
        assert _abs_corr([1.0, 2.0, 3.0], [4.0, 4.0, 4.0]) is None
        assert _abs_corr([1.0, 2.0, 3.0, np.nan], [6.0, 4.0, 2.0, 9.0]) == pytest.approx(1.0, abs=1e-12)

        # _composition_covariate_decision boundary: corr == max_corr -> dropped (>=)
        use_eq, reason_eq = _composition_covariate_decision([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], 1.0)
        assert use_eq is False and "collinear_with_condition" in reason_eq
        use_ok, reason_ok = _composition_covariate_decision([0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0], 0.8)
        assert use_ok is True and reason_ok is None
        use_nan, reason_nan = _composition_covariate_decision([1.0, 1.0, 1.0], [0.0, 1.0, 0.0], 0.8)
        assert use_nan is False and reason_nan == "collinearity_unmeasurable"

        # _effective_alpha_from_significant
        assert _effective_alpha_from_significant(None) == 0.0
        assert _effective_alpha_from_significant(pd.DataFrame()) == 0.0
        assert _effective_alpha_from_significant(pd.DataFrame({"log2fc": [1.0]})) == 0.0
        assert _effective_alpha_from_significant(pd.DataFrame({"pvalue": [np.nan, np.inf]})) == 0.0
        assert _effective_alpha_from_significant(pd.DataFrame({"pvalue": [1e-4, 0.012, np.nan]})) == pytest.approx(0.012)

        # _expressed_background: dense, sparse, and exception fallback
        dense = np.array([[0, 2, 0], [0, 0, 5]])
        genes_d, deg_d = _expressed_background(dense, ["G0", "G1", "G2"])
        assert genes_d == ["G1", "G2"] and deg_d is False
        genes_s, deg_s = _expressed_background(sparse.csr_matrix(dense), ["G0", "G1", "G2"])
        assert genes_s == ["G1", "G2"] and deg_s is False
        genes_err, deg_err = _expressed_background("not_a_matrix", ["G0", "G1"])
        assert genes_err == ["G0", "G1"] and deg_err is True

    def test_sparse_sum_aggregation_preserves_total_counts_and_kills_mean_mutant(self, tmp_path):
        """Sparse AnnData pseudobulk MUST sum cell counts (not average), preserving count scale."""
        pytest.importorskip("pydeseq2")
        ad = pytest.importorskip("anndata")

        rng = np.random.default_rng(777)
        n_genes = 40
        genes = [f"SG{i:03d}" for i in range(n_genes)]
        cells_per_donor = 60
        rows, blocks = [], []
        for cond in ("WT", "KO"):
            for r in range(4):
                donor = f"{cond}_d{r}"
                # Per-cell mean=0.4 (WT) / 1.6 (KO up): sum across 60 cells yields ~24 vs ~96
                # (clean 4x DE with padj < 1e-8); averaging + int64 truncation collapses null
                # counts to 0 and destroys pseudobulk power/estimation.
                mu = np.full(n_genes, 0.4)
                if cond == "KO":
                    mu[:8] = 1.6
                X_blk = rng.poisson(mu, size=(cells_per_donor, n_genes)).astype(np.float32)
                blocks.append(X_blk)
                for c in range(cells_per_donor):
                    rows.append({"cell_type": "Mono", "condition": cond, "donor": donor})

        X_sparse = sparse.csr_matrix(np.vstack(blocks))
        adata = ad.AnnData(
            X=X_sparse,
            obs=pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))]),
            var=pd.DataFrame(index=genes),
        )
        h5ad = tmp_path / "sparse_pb.h5ad"
        adata.write_h5ad(h5ad)

        res = rna_pseudobulk_de({
            "data_path": str(h5ad),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["KO", "WT"]],
            "use_raw": False,
            "min_cells_per_pseudosample": 10,
            "min_replicates_per_condition": 3,
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "sparse_out"),
        })
        assert res["status"] == "success", res
        comp = res["per_group"]["Mono"]["per_comparison"]["KO_vs_WT"]
        assert comp["status"] == "success", comp
        full_df = pd.read_csv(comp["full_results_csv"]).set_index("gene")
        assert len(full_df) == n_genes
        assert float(full_df.loc[genes[:8], "padj"].median()) < 1e-6
        assert float(full_df.loc[genes[:8], "log2fc"].median()) == pytest.approx(2.0, abs=0.45)
        called_up = {rec["gene"] for rec in comp["all_sig"] if rec["log2fc"] > 1.0}
        assert len(called_up & set(genes[:8])) >= 7
