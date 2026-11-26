"""Analytical oracle tests for chromatin_diffacc.py, chromatin_bulk_diffacc.py,
and epigenomic_spikein_de.py.

Verifies:
  - Pure helpers: _gate_da_by_convergence, _binary_pct, _normalise_comparisons,
    _load_counts_matrix validation, _prepare_replicate_matrix technical-replicate summing,
    resolve_replicate_spikein_size_factors geometric-mean-1 math
  - Real pyDESeq2 / scanpy oracles (in aria-rna-env):
      * scATAC _pseudobulk_da: signed LFC & exact peak-ID recovery (kills D3-003a & D3-014),
        WT vs KO exact negation, KO vs KO_ctrl alphabetical trap, cell/peak permutation
        invariance, batch covariate & paired donor variance removal, edge cases
      * scATAC _per_cluster_accessibility: exact marker peak recovery & pct_in > pct_out
      * Bulk ATAC chromatin_bulk_diffacc: unmocked pyDESeq2 signed recovery (kills D3-003b,
        D3-003c, D3-007), technical replicate collapsing, sample/peak permutation invariance,
        batch covariate adjustment, edge cases
      * Heterologous spike-in size factors: global unidirectional chromatin loss recovery
        plus alphabetical KO vs WT reference-level check
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from aria.scripts.chromatin_bulk_diffacc import (
    _load_counts_matrix,
    _normalise_comparisons,
    _prepare_replicate_matrix,
    chromatin_bulk_diffacc,
)
from aria.scripts.chromatin_diffacc import (
    _binary_pct,
    _gate_da_by_convergence,
    _per_cluster_accessibility,
    _pseudobulk_da,
)
from aria.utils.epigenomic_spikein_de import resolve_replicate_spikein_size_factors


class TestChromatinDiffAccPureHelpersOracle:
    def test_gate_da_by_convergence_exact_filtering_and_sorting(self):
        df = pd.DataFrame(
            {
                "log2FoldChange": [2.0, -1.8, 0.01, 1.5, -2.1],
                "padj": [0.02, 0.001, 0.0005, 0.20, 0.01],
                "lfc_converged": [True, True, False, True, float("nan")],
            },
            index=["P1", "P2", "P_nonconv", "P_ns", "P_nanconv"],
        )
        sig, info = _gate_da_by_convergence(df, padj_max=0.05)
        assert info == {"convergence_available": True, "n_nonconverged_excluded": 2}
        assert list(sig.index) == ["P2", "P1"]

        # Without lfc_converged column -> no fabrication, keeps all padj < 0.05 sorted
        df_no = df.drop(columns=["lfc_converged"])
        sig_no, info_no = _gate_da_by_convergence(df_no, padj_max=0.05)
        assert info_no == {"convergence_available": False, "n_nonconverged_excluded": 0}
        assert list(sig_no.index) == ["P_nonconv", "P2", "P_nanconv", "P1"]

    def test_binary_pct_dense_and_sparse_parity(self):
        X = np.array([[0, 2], [3, 0], [1, 0], [0, 5]], dtype=float)
        mask = np.array([True, True, False, True])
        # Col 0 in mask: [0, 3, 0] -> 1/3; Col 1 in mask: [2, 0, 5] -> 2/3
        assert _binary_pct(X, mask, 0) == pytest.approx(1.0 / 3.0, abs=1e-12)
        assert _binary_pct(sparse.csr_matrix(X), mask, 0) == pytest.approx(1.0 / 3.0, abs=1e-12)
        assert _binary_pct(sparse.csr_matrix(X), mask, 1) == pytest.approx(2.0 / 3.0, abs=1e-12)

    def test_normalise_comparisons_all_input_forms(self):
        assert _normalise_comparisons(None) == []
        assert _normalise_comparisons("KO_vs_WT") == [{"test": "KO", "reference": "WT"}]
        assert _normalise_comparisons("KO vs WT") == [{"test": "KO", "reference": "WT"}]
        assert _normalise_comparisons({"numerator": "KO", "denominator": "WT"}) == [
            {"test": "KO", "reference": "WT"}
        ]
        assert _normalise_comparisons([["A", "B"], {"test": "C", "ref": "D"}]) == [
            {"test": "A", "reference": "B"},
            {"test": "C", "reference": "D"},
        ]

    def test_load_counts_matrix_rejects_invalid_inputs(self, tmp_path):
        good = tmp_path / "good.tsv"
        good.write_text("peak_id\ts1\ts2\nchr1:1-10\t5\t10\n", encoding="utf-8")
        df = _load_counts_matrix(good)
        assert df.loc["chr1:1-10", "s1"] == 5

        for name, content in [
            ("neg.tsv", "peak_id\ts1\ts2\nchr1:1-10\t-1\t5\n"),
            ("float.tsv", "peak_id\ts1\ts2\nchr1:1-10\t1.5\t5\n"),
            ("nan.tsv", "peak_id\ts1\ts2\nchr1:1-10\tnan\t5\n"),
            ("dup.tsv", "peak_id\ts1\ts2\nchr1:1-10\t1\t2\nchr1:1-10\t3\t4\n"),
            ("badcol.tsv", "gene\ts1\ts2\nG1\t1\t2\n"),
        ]:
            p = tmp_path / name
            p.write_text(content, encoding="utf-8")
            with pytest.raises(ValueError):
                _load_counts_matrix(p)

    def test_prepare_replicate_matrix_sums_technical_reps_and_validates_covariates(self):
        counts = pd.DataFrame(
            {"s1a": [10, 20], "s1b": [5, 15], "s2": [30, 40]},
            index=["P1", "P2"],
        )
        meta = pd.DataFrame(
            {
                "condition": ["KO", "KO", "WT"],
                "replicate": ["R1", "R1", "R2"],
                "batch": ["B1", "B1", "B2"],
                "inconsistent": ["X", "Y", "X"],
            },
            index=["s1a", "s1b", "s2"],
        )
        rep_c, rep_m, usable, dropped = _prepare_replicate_matrix(
            counts, meta, "condition", "replicate", ["batch", "inconsistent", "absent"]
        )
        assert rep_c.loc["P1", "KO::R1"] == 15
        assert rep_c.loc["P2", "KO::R1"] == 35
        assert rep_c.loc["P1", "WT::R2"] == 30
        assert rep_m.loc["KO::R1", "n_technical_samples"] == 2
        assert usable == ["batch"]
        dropped_names = {d["covariate"] for d in dropped}
        assert dropped_names == {"inconsistent", "absent"}

    def test_resolve_spikein_size_factors_exact_geometric_mean_and_tech_rep_sum(self):
        # Technical reps s1a (100 spike reads) + s1b (300 spike reads) -> rep R1 exposure = 400
        # Rep R2 exposure = 100
        # Geometric mean = sqrt(400 * 100) = 200 -> size factors = [2.0, 0.5]
        raw_meta = pd.DataFrame(
            {"condition": ["KO", "KO", "WT"], "replicate": ["R1", "R1", "R2"]},
            index=["s1a", "s1b", "s2"],
        )
        rep_meta = pd.DataFrame(
            {"condition": ["KO", "WT"], "replicate": ["R1", "R2"]},
            index=["KO::R1", "WT::R2"],
        )
        sf, info = resolve_replicate_spikein_size_factors(
            raw_meta, rep_meta, "condition", "replicate",
            spikein_counts={"s1a": 100.0, "s1b": 300.0, "s2": 100.0},
        )
        assert info["used"] is True
        assert sf.loc["KO::R1"] == pytest.approx(2.0, abs=1e-12)
        assert sf.loc["WT::R2"] == pytest.approx(0.5, abs=1e-12)
        assert float(np.exp(np.mean(np.log(sf.to_numpy())))) == pytest.approx(1.0, abs=1e-12)

        # Inverse coverage scale factors a_i = 1 / N_i MUST give identical size factors!
        sf_inv, _ = resolve_replicate_spikein_size_factors(
            raw_meta, rep_meta, "condition", "replicate",
            spikein_scale_factors={"s1a": 1.0 / 100.0, "s1b": 1.0 / 300.0, "s2": 1.0 / 100.0},
        )
        np.testing.assert_allclose(sf_inv.to_numpy(), sf.to_numpy(), rtol=1e-12, atol=1e-12)


class TestScAtacDiffAccOracles:
    def _make_scatac_adata(
        self,
        *,
        n_peaks: int = 90,
        n_up: int = 12,
        n_down: int = 12,
        n_reps: int = 4,
        cells_per_rep: int = 20,
        test_label: str = "KO",
        ref_label: str = "WT",
        batch_shift: bool = False,
        seed: int = 301,
    ):
        import anndata as ad

        rng = np.random.default_rng(seed)
        peaks = [f"chr1:{1000 + i * 500}-{1400 + i * 500}" for i in range(n_peaks)]
        base = rng.uniform(1.5, 4.0, size=n_peaks)
        fc = np.ones(n_peaks)
        fc[:n_up] = 2.0 ** (1.4 if batch_shift else 2.2)
        fc[n_up : n_up + n_down] = 2.0 ** (-1.4 if batch_shift else -2.2)
        batch_fc = np.where(rng.random(n_peaks) > 0.5, 6.0, 1.0 / 6.0) if batch_shift else np.ones(n_peaks)

        rows, obs_rows = [], []
        for cond in (ref_label, test_label):
            for r in range(n_reps):
                batch = "B1" if r < (n_reps // 2) else "B2"
                rep = f"{cond}_R{r}"
                rate = base * (fc if cond == test_label else 1.0) * (batch_fc if batch == "B2" else 1.0)
                for _ in range(cells_per_rep):
                    rows.append(rng.poisson(rate))
                    obs_rows.append({"condition": cond, "replicate": rep, "batch": batch, "leiden": "0"})
        X = sparse.csr_matrix(np.asarray(rows, dtype=np.int64))
        obs = pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(len(obs_rows))])
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=peaks))
        return adata, peaks[:n_up], peaks[n_up : n_up + n_down], peaks[n_up + n_down :]

    def test_pseudobulk_da_signed_recovery_reverse_renaming_and_permutation(self):
        """Kills D3-003a (negated log2FoldChange) and D3-014 (scrambled padj)."""
        pytest.importorskip("pydeseq2")
        adata, up_p, dn_p, null_p = self._make_scatac_adata(seed=301)

        res_kw = _pseudobulk_da(
            adata, "condition", "replicate", [{"test": "KO", "reference": "WT"}],
            padj_max=0.05, lfc_min=0.5, top_n=100, min_cells=10, min_reps=3,
            allow_mock=False, warnings=[],
        )
        assert res_kw["ran"] is True
        comp_kw = res_kw["comparisons"][0]
        assert comp_kw["status"] == "success"
        hit_map = {r["peak"]: r for r in comp_kw["top_peaks"]}

        up_found = [p for p in up_p if p in hit_map]
        dn_found = [p for p in dn_p if p in hit_map]
        assert len(up_found) >= 10
        assert len(dn_found) >= 10
        assert all(hit_map[p]["log2fc"] > 1.0 for p in up_found)
        assert all(hit_map[p]["log2fc"] < -1.0 for p in dn_found)
        assert len(set(hit_map) & set(null_p)) <= 2

        # Reverse contrast WT vs KO -> exact sign flip & swapped n_up/n_down
        res_wk = _pseudobulk_da(
            adata, "condition", "replicate", [{"test": "WT", "reference": "KO"}],
            padj_max=0.05, lfc_min=0.5, top_n=100, min_cells=10, min_reps=3,
            allow_mock=False, warnings=[],
        )
        comp_wk = res_wk["comparisons"][0]
        assert comp_wk["n_up"] == comp_kw["n_down"]
        assert comp_wk["n_down"] == comp_kw["n_up"]
        wk_map = {r["peak"]: r for r in comp_wk["top_peaks"]}
        for p in up_found:
            assert wk_map[p]["log2fc"] == pytest.approx(-hit_map[p]["log2fc"], abs=2e-3)

        # Level renaming prefix trap KO vs KO_ctrl + cell & peak permutation invariance
        rng = np.random.default_rng(77)
        ad_perm = adata[rng.permutation(adata.n_obs), rng.permutation(adata.n_vars)].copy()
        ad_perm.obs["condition"] = ad_perm.obs["condition"].map({"KO": "KO", "WT": "KO_ctrl"})
        res_ren = _pseudobulk_da(
            ad_perm, "condition", "replicate", [{"test": "KO", "reference": "KO_ctrl"}],
            padj_max=0.05, lfc_min=0.5, top_n=100, min_cells=10, min_reps=3,
            allow_mock=False, warnings=[],
        )
        comp_ren = res_ren["comparisons"][0]
        ren_map = {r["peak"]: r for r in comp_ren["top_peaks"]}
        assert set(ren_map) == set(hit_map)
        for p in hit_map:
            assert ren_map[p]["log2fc"] == pytest.approx(hit_map[p]["log2fc"], abs=2e-3)

    def test_pseudobulk_da_batch_covariate_recovers_masked_peaks(self):
        pytest.importorskip("pydeseq2")
        adata, up_p, dn_p, _ = self._make_scatac_adata(n_reps=6, batch_shift=True, seed=302)
        planted = set(up_p) | set(dn_p)

        res_adj = _pseudobulk_da(
            adata, "condition", "replicate", [{"test": "KO", "reference": "WT"}],
            padj_max=0.05, lfc_min=0.5, top_n=100, min_cells=10, min_reps=3,
            allow_mock=False, warnings=[], covariates=["batch"],
        )
        res_unadj = _pseudobulk_da(
            adata, "condition", "replicate", [{"test": "KO", "reference": "WT"}],
            padj_max=0.05, lfc_min=0.5, top_n=100, min_cells=10, min_reps=3,
            allow_mock=False, warnings=[], covariates=None,
        )
        adj_hits = {r["peak"] for r in res_adj["comparisons"][0]["top_peaks"]} & planted
        unadj_hits = {r["peak"] for r in res_unadj["comparisons"][0]["top_peaks"]} & planted
        assert len(adj_hits) >= 18
        assert len(adj_hits) > len(unadj_hits) + 8

    def test_per_cluster_wilcoxon_recovers_planted_cluster_markers_and_single_cluster_skips(self):
        pytest.importorskip("scanpy")
        import anndata as ad

        rng = np.random.default_rng(303)
        n_peaks = 40
        peaks = [f"chr1:{i*1000}-{i*1000+300}" for i in range(n_peaks)]
        # Cluster A (30 cells): peaks 0..4 strongly open; Cluster B (30 cells): peaks 5..9 strongly open
        X_a = rng.poisson(0.2, size=(30, n_peaks))
        X_b = rng.poisson(0.2, size=(30, n_peaks))
        X_a[:, :5] = rng.poisson(15.0, size=(30, 5))
        X_b[:, 5:10] = rng.poisson(15.0, size=(30, 5))
        X = np.vstack([X_a, X_b]).astype(np.float32)
        obs = pd.DataFrame({"leiden": ["0"] * 30 + ["1"] * 30}, index=[f"c{i}" for i in range(60)])
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=peaks))

        res = _per_cluster_accessibility(adata, "leiden", padj_max=0.05, lfc_min=0.5, top_n=20, warnings=[])
        assert res["ran"] is True
        c0_up = {r["peak"]: r for r in res["da_peaks_by_cluster"]["0"] if r["log2fc"] > 0}
        c1_up = {r["peak"]: r for r in res["da_peaks_by_cluster"]["1"] if r["log2fc"] > 0}
        assert set(peaks[:5]).issubset(c0_up)
        assert set(peaks[5:10]).issubset(c1_up)
        assert all(c0_up[p]["pct_in"] > c0_up[p]["pct_out"] + 0.5 for p in peaks[:5])

        # Single cluster skips cleanly
        adata_1c = adata.copy()
        adata_1c.obs["leiden"] = "0"
        res_1c = _per_cluster_accessibility(adata_1c, "leiden", 0.05, 0.5, 20, [])
        assert res_1c["ran"] is False


class TestBulkAtacDiffAccOracles:
    def _write_bulk_atac_fixture(
        self,
        tmp_path,
        *,
        tag: str = "base",
        n_peaks: int = 80,
        n_up: int = 12,
        n_down: int = 12,
        n_reps: int = 4,
        test_label: str = "KO",
        ref_label: str = "WT",
        seed: int = 401,
        permute: bool = False,
    ):
        rng = np.random.default_rng(seed)
        peaks = [f"chr1:{1000 + i * 500}-{1400 + i * 500}" for i in range(n_peaks)]
        base_mu = rng.uniform(80.0, 300.0, size=n_peaks)
        fc = np.ones(n_peaks)
        fc[:n_up] = 2.0 ** 2.3
        fc[n_up : n_up + n_down] = 2.0 ** (-2.3)

        # Include 2 technical replicates on KO_R0 to exercise technical-replicate summing
        sample_specs = []
        for cond in (ref_label, test_label):
            for r in range(n_reps):
                rep = f"R{r}"
                if cond == test_label and r == 0:
                    sample_specs.append((f"{cond}_{rep}_t1", cond, rep, 0.5))
                    sample_specs.append((f"{cond}_{rep}_t2", cond, rep, 0.5))
                else:
                    sample_specs.append((f"{cond}_{rep}", cond, rep, 1.0))

        alpha = 0.03
        r_param = 1.0 / alpha
        mat = {}
        meta_rows = []
        for sid, cond, rep, frac in sample_specs:
            mu = base_mu * (fc if cond == test_label else 1.0) * frac
            mat[sid] = rng.negative_binomial(r_param, r_param / (r_param + mu))
            meta_rows.append({"sample_id": sid, "condition": cond, "replicate": rep})

        counts = pd.DataFrame(mat, index=peaks)
        counts.iloc[-1, :] = 0  # all-zero peak
        meta = pd.DataFrame(meta_rows)

        if permute:
            p_rng = np.random.default_rng(seed + 999)
            counts = counts.iloc[p_rng.permutation(len(counts)), p_rng.permutation(counts.shape[1])]
            meta = meta.iloc[p_rng.permutation(len(meta))]

        c_path = tmp_path / f"{tag}_counts.tsv"
        m_path = tmp_path / f"{tag}_meta.tsv"
        counts_out = counts.copy()
        counts_out.insert(0, "peak_id", counts_out.index)
        counts_out.to_csv(c_path, sep="\t", index=False)
        meta.to_csv(m_path, sep="\t", index=False)
        return c_path, m_path, peaks[:n_up], peaks[n_up : n_up + n_down], peaks[n_up + n_down : -1]

    def test_bulk_atac_unmocked_signed_recovery_reverse_and_permutation(self, tmp_path):
        """Kills D3-003b, D3-003c (swap test/ref), and D3-007 (broken _run_deseq2)."""
        pytest.importorskip("pydeseq2")
        c_path, m_path, up_p, dn_p, null_p = self._write_bulk_atac_fixture(tmp_path, tag="kw", seed=401)

        res_kw = chromatin_bulk_diffacc({
            "counts_matrix_path": str(c_path),
            "sample_metadata_path": str(m_path),
            "comparisons": [{"test": "KO", "reference": "WT"}],
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "out_kw"),
        })
        assert res_kw["status"] == "success" and res_kw["ran"] is True
        assert res_kw["n_replicate_samples"] == 8  # 9 tech samples collapsed to 8 bio reps
        assert res_kw["n_peaks_tested"] == 79      # 1 all-zero peak dropped
        comp_kw = res_kw["comparisons"][0]
        full_kw = pd.read_csv(comp_kw["full_results_csv"]).set_index("peak")
        sig_kw = full_kw[full_kw["significant"]]

        assert len(set(sig_kw.index) & set(up_p)) >= 11
        assert len(set(sig_kw.index) & set(dn_p)) >= 11
        assert len(set(sig_kw.index) & set(null_p)) <= 2
        assert (full_kw.loc[up_p, "log2FoldChange"] > 1.0).all()
        assert (full_kw.loc[dn_p, "log2FoldChange"] < -1.0).all()

        # Reverse contrast WT vs KO -> exact sign flip
        res_wk = chromatin_bulk_diffacc({
            "counts_matrix_path": str(c_path),
            "sample_metadata_path": str(m_path),
            "comparisons": [{"test": "WT", "reference": "KO"}],
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "out_wk"),
        })
        full_wk = pd.read_csv(res_wk["comparisons"][0]["full_results_csv"]).set_index("peak")
        common = full_kw.index.intersection(full_wk.index)
        np.testing.assert_allclose(
            full_wk.loc[common, "log2FoldChange_raw"].to_numpy(),
            -full_kw.loc[common, "log2FoldChange_raw"].to_numpy(),
            rtol=1e-4,
            atol=1e-4,
        )
        np.testing.assert_allclose(
            full_wk.loc[common, "log2FoldChange"].to_numpy(),
            -full_kw.loc[common, "log2FoldChange"].to_numpy(),
            rtol=1e-2,
            atol=5e-3,
        )

        # Permuted sample & peak order -> identical LFCs & padj
        c_perm, m_perm, _, _, _ = self._write_bulk_atac_fixture(
            tmp_path, tag="perm", seed=401, permute=True
        )
        res_p = chromatin_bulk_diffacc({
            "counts_matrix_path": str(c_perm),
            "sample_metadata_path": str(m_perm),
            "comparisons": [{"test": "KO", "reference": "WT"}],
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "out_perm"),
        })
        full_p = pd.read_csv(res_p["comparisons"][0]["full_results_csv"]).set_index("peak").loc[full_kw.index]
        np.testing.assert_allclose(
            full_p["log2FoldChange"].to_numpy(),
            full_kw["log2FoldChange"].to_numpy(),
            rtol=1e-4,
            atol=1e-4,
        )

    def test_spikein_global_unidirectional_loss_recovered_for_treated_vs_ctrl(self, tmp_path):
        """With 80% of peaks globally reduced 4x in treated vs ctrl, spike-in size factors
        MUST recover negative LFCs (whereas median-of-ratios normalizes the global shift away)."""
        pytest.importorskip("pydeseq2")
        rng = np.random.default_rng(405)
        n_peaks = 60
        peaks = [f"chr1:{i*500}-{i*500+300}" for i in range(n_peaks)]
        base_mu = rng.uniform(200.0, 400.0, size=n_peaks)
        # 48 peaks lose 75% signal (0.25x = -2 log2FC); 12 peaks invariant (1.0x)
        fc = np.full(n_peaks, 0.25)
        fc[48:] = 1.0

        samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treated_1", "treated_2", "treated_3"]
        mat = {}
        meta_rows = []
        spikein = {}
        for s in samples:
            cond = "treated" if s.startswith("treated") else "ctrl"
            rep = s.split("_")[1]
            mu = base_mu * (fc if cond == "treated" else 1.0)
            mat[s] = rng.negative_binomial(50.0, 50.0 / (50.0 + mu))
            meta_rows.append({"sample_id": s, "condition": cond, "replicate": rep})
            spikein[s] = 1000.0  # Equal exogenous spike-in reads across all libraries

        c_path = tmp_path / "sp_counts.tsv"
        m_path = tmp_path / "sp_meta.tsv"
        df_c = pd.DataFrame(mat, index=peaks)
        df_c.insert(0, "peak_id", df_c.index)
        df_c.to_csv(c_path, sep="\t", index=False)
        pd.DataFrame(meta_rows).to_csv(m_path, sep="\t", index=False)

        res_sp = chromatin_bulk_diffacc({
            "data_type": "CUT_AND_RUN",
            "counts_matrix_path": str(c_path),
            "sample_metadata_path": str(m_path),
            "comparisons": [{"test": "treated", "reference": "ctrl"}],
            "spikein_counts": spikein,
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "sp_out"),
        })
        assert res_sp["status"] == "success" and res_sp["ran"] is True
        assert res_sp["normalization_mode"] == "heterologous_spikein"
        comp = res_sp["comparisons"][0]
        full = pd.read_csv(comp["full_results_csv"]).set_index("peak")
        # Under spike-in, the 48 lost peaks have LFC ~ -2.0 (negative) and the 12 invariant ~ 0.0!
        assert float(full.loc[peaks[:48], "log2FoldChange"].median()) == pytest.approx(-2.0, abs=0.35)
        assert float(full.loc[peaks[48:], "log2FoldChange"].median()) == pytest.approx(0.0, abs=0.35)
        assert comp["n_down"] >= 40
        assert comp["n_up"] == 0

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "D2-BUG-6: epigenomic_spikein_de.run_deseq2_with_external_size_factors omits "
            "Treatment(reference=denominator) in design formula; on pydeseq2>=0.5.4 when "
            "numerator < denominator alphabetically (KO vs WT), apeGLM shrinkage fails to "
            "find coeff [T.KO] (or uses inverted reference level)"
        ),
    )
    def test_spikein_alphabetical_trap_ko_vs_wt_applies_shrinkage_with_correct_sign(self, tmp_path):
        pytest.importorskip("pydeseq2")
        rng = np.random.default_rng(406)
        n_peaks = 60
        peaks = [f"chr1:{i*500}-{i*500+300}" for i in range(n_peaks)]
        base_mu = rng.uniform(200.0, 400.0, size=n_peaks)
        fc = np.full(n_peaks, 0.25)
        fc[48:] = 1.0

        samples = ["WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"]
        mat = {}
        meta_rows = []
        spikein = {}
        for s in samples:
            cond = "KO" if s.startswith("KO") else "WT"
            rep = s.split("_")[1]
            mu = base_mu * (fc if cond == "KO" else 1.0)
            mat[s] = rng.negative_binomial(50.0, 50.0 / (50.0 + mu))
            meta_rows.append({"sample_id": s, "condition": cond, "replicate": rep})
            spikein[s] = 1000.0

        c_path = tmp_path / "sp_kw_counts.tsv"
        m_path = tmp_path / "sp_kw_meta.tsv"
        df_c = pd.DataFrame(mat, index=peaks)
        df_c.insert(0, "peak_id", df_c.index)
        df_c.to_csv(c_path, sep="\t", index=False)
        pd.DataFrame(meta_rows).to_csv(m_path, sep="\t", index=False)

        res_kw = chromatin_bulk_diffacc({
            "data_type": "CUT_AND_RUN",
            "counts_matrix_path": str(c_path),
            "sample_metadata_path": str(m_path),
            "comparisons": [{"test": "KO", "reference": "WT"}],
            "spikein_counts": spikein,
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "output_dir": str(tmp_path / "sp_kw_out"),
        })
        assert res_kw["status"] == "success" and res_kw["ran"] is True
        comp = res_kw["comparisons"][0]
        # Shrinkage MUST apply for KO vs WT just as it did for treated vs ctrl!
        assert comp["lfc_shrinkage"]["applied"] is True
        full = pd.read_csv(comp["full_results_csv"]).set_index("peak")
        assert float(full.loc[peaks[:48], "log2FoldChange"].median()) == pytest.approx(-2.0, abs=0.35)
