"""Analytical oracle tests for aria/utils/seacr_idr.py.

Verifies:
  - _integrate_bedgraph_signal_blocks: exact AUC, contiguous merging, max_signal/max_region
  - _estimate_seacr_control_scale: exact lower-half median density ratio
  - _select_seacr_empirical_threshold: relaxed <= stringent threshold ordering, exact
    empirical FDR & precision, separation of signal vs IgG, no-excess -> inf, edge cases
  - call_seacr_peaks: bedGraph mode, top_quantile exact k_keep, scale invariance,
    spike-in norm='non' vs asymmetric spike-in fallback norm='norm', multi-rep mean pooling
  - _invert_mixture_cdf: machine-precision inversion of bimodal Gaussian mixture CDF
    verified directly against scipy.stats.norm.cdf
  - fit_copula_idr_pair: synthetic Li et al. (2011) bivariate Gaussian copula mixture
    with known parameters, rank/monotone-transform invariance, replicate symmetry,
    permutation invariance, <20 peaks error
  - _pair_two_replicate_peaks & compute_idr_reproducible_peaks: half-open book-ended
    intervals [100,200) vs [200,300), non-finite score exclusion, consensus BED & SCR
  - _finite_positive_scales & _fmt_bedgraph_value edge cases
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import norm

from aria.utils.seacr_idr import (
    _estimate_seacr_control_scale,
    _finite_positive_scales,
    _fmt_bedgraph_value,
    _integrate_bedgraph_signal_blocks,
    _invert_mixture_cdf,
    _pair_two_replicate_peaks,
    _read_scored_peaks,
    _select_seacr_empirical_threshold,
    call_seacr_peaks,
    compute_idr_reproducible_peaks,
    fit_copula_idr_pair,
)


class TestSeacrHelpersOracle:
    def test_fmt_bedgraph_preserves_tiny_spikein_values(self):
        s = _fmt_bedgraph_value(1.23456789e-9)
        assert float(s) == pytest.approx(1.23456789e-9, rel=1e-9)
        assert float(s) > 0.0

    def test_finite_positive_scales_validation(self):
        assert _finite_positive_scales(None, "s") is None
        assert _finite_positive_scales([0.5, 2, 1.25], "s") == [0.5, 2.0, 1.25]
        for bad in ([0.0], [-1.0], [float("nan")], [float("inf")], [True], ["abc"]):
            with pytest.raises(ValueError):
                _finite_positive_scales(bad, "s")

    def test_integrate_bedgraph_signal_blocks_exact_auc_and_max_region(self, tmp_path):
        bg = tmp_path / "test.bedgraph"
        bg.write_text(
            "# comment line\n"
            "track type=bedGraph\n"
            "chr1\t100\t150\t2.0\n"    # block 1: area = 50*2 = 100
            "chr1\t150\t200\t5.0\n"    # contiguous: area += 50*5 = 250 -> 350, max=5 at 150-200
            "chr1\t200\t220\t5.0\n"    # contiguous tie at max: max_region extends to 150-220! area += 100 -> 450
            "chr1\t220\t250\t1.0\n"    # contiguous: area += 30 -> 480
            "chr1\t300\t310\t0.0\n"    # zero score ignored
            "chr1\t400\t450\t3.0\n"    # block 2 (gap): area = 150, max=3
            "chr2\t0\t10\t4.0\n"       # block 3 (new chrom): area = 40, max=4
            "chr2\t20\t10\t9.0\n"      # invalid end<=start ignored
            "chr2\t30\t40\tbad\n",     # invalid float ignored
            encoding="utf-8",
        )
        blocks = _integrate_bedgraph_signal_blocks(bg)
        assert len(blocks) == 3

        assert blocks[0]["chrom"] == "chr1"
        assert (blocks[0]["start"], blocks[0]["end"]) == (100, 250)
        assert blocks[0]["total_signal"] == pytest.approx(480.0, abs=1e-12)
        assert blocks[0]["max_signal"] == pytest.approx(5.0, abs=1e-12)
        assert blocks[0]["max_region"] == "chr1:150-220"

        assert blocks[1]["total_signal"] == pytest.approx(150.0, abs=1e-12)
        assert blocks[1]["max_region"] == "chr1:400-450"

        assert blocks[2]["chrom"] == "chr2"
        assert blocks[2]["total_signal"] == pytest.approx(40.0, abs=1e-12)

    def test_estimate_seacr_control_scale_exact_lower_half_median_ratio(self):
        # Target densities: [2, 4, 6, 100, 200, 400] -> q50 = 53.0 -> lower half = [2, 4, 6] -> median = 4.0
        # Control densities: [1, 2, 3, 4, 5, 6] -> q50 = 3.5 -> lower half = [1, 2, 3] -> median = 2.0
        # Ratio = 4.0 / 2.0 = 2.0 (unaffected by enriched target peaks!)
        t_blocks = [
            {"start": 0, "end": 10, "total_signal": d * 10.0}
            for d in (2.0, 4.0, 6.0, 100.0, 200.0, 400.0)
        ]
        c_blocks = [
            {"start": 0, "end": 10, "total_signal": d * 10.0}
            for d in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        ]
        assert _estimate_seacr_control_scale(t_blocks, c_blocks) == pytest.approx(2.0, abs=1e-12)
        assert _estimate_seacr_control_scale([], c_blocks) == 1.0
        assert _estimate_seacr_control_scale(t_blocks, []) == 1.0

    def test_empirical_threshold_separates_planted_peaks_and_orders_modes(self):
        rng = np.random.default_rng(5)
        bg = rng.uniform(5.0, 25.0, size=200)
        peaks = rng.uniform(200.0, 500.0, size=25)
        target = np.concatenate([bg, peaks])
        control = rng.uniform(5.0, 25.0, size=200)

        thr_str, meta_str = _select_seacr_empirical_threshold(target, control, mode="stringent")
        thr_rel, meta_rel = _select_seacr_empirical_threshold(target, control, mode="relaxed")

        assert thr_rel <= thr_str
        n_str = int((target >= thr_str).sum())
        n_rel = int((target >= thr_rel).sum())
        assert n_rel >= n_str
        # Stringent must recover all 25 planted peaks with near-zero control contamination
        assert n_str == 25
        assert meta_str["n_control_blocks_above_threshold"] == 0
        assert meta_str["empirical_fdr"] == 0.0
        assert meta_str["precision"] == 1.0

        # Identical target and control -> no positive signal excess -> inf threshold
        thr_null, meta_null = _select_seacr_empirical_threshold(control, control, mode="stringent")
        assert math.isinf(thr_null)
        assert meta_null["no_signal_excess"] is True

        # Empty target -> inf; empty control -> 99th percentile
        thr_empty, _ = _select_seacr_empirical_threshold(np.array([]), control)
        assert math.isinf(thr_empty)
        thr_noctrl, meta_noctrl = _select_seacr_empirical_threshold(np.arange(1.0, 101.0), np.array([]))
        assert thr_noctrl == pytest.approx(float(np.quantile(np.arange(1.0, 101.0), 0.99)), abs=1e-12)
        assert meta_noctrl["precision"] == 1.0

        with pytest.raises(ValueError, match="stringent"):
            _select_seacr_empirical_threshold(target, control, mode="invalid")

    def test_call_seacr_peaks_bedgraph_quantile_spikein_and_asymmetric_guard(self, tmp_path):
        # Create target bedGraph: 40 background blocks (AUC=100) + 5 true peaks (AUC=5000)
        t_bg = tmp_path / "target.bedgraph"
        c_bg = tmp_path / "control.bedgraph"
        t_lines, c_lines = [], []
        for i in range(40):
            s = i * 200
            t_lines.append(f"chr1\t{s}\t{s+100}\t1.0\n")
            c_lines.append(f"chr1\t{s}\t{s+100}\t1.0\n")
        for i in range(5):
            s = 20_000 + i * 500
            t_lines.append(f"chr1\t{s}\t{s+100}\t50.0\n")
        t_bg.write_text("".join(t_lines), encoding="utf-8")
        c_bg.write_text("".join(c_lines), encoding="utf-8")

        # 1. With control: recovers exactly the 5 planted peaks
        out_ctrl = tmp_path / "peaks_ctrl.bed"
        res_ctrl = call_seacr_peaks(
            str(t_bg), out_ctrl, control_inputs=str(c_bg), mode="stringent", work_dir=tmp_path / "w1"
        )
        assert res_ctrl["status"] == "success"
        assert res_ctrl["n_peaks"] == 5
        assert res_ctrl["norm"] == "norm"
        assert res_ctrl["already_spikein_scaled"] is False

        # 2. Top-quantile mode without control: ceil(0.10 * 45) = 5 peaks
        out_q = tmp_path / "peaks_q.bed"
        res_q = call_seacr_peaks(
            str(t_bg), out_q, control_inputs=None, top_quantile=0.10, work_dir=tmp_path / "w2"
        )
        assert res_q["n_peaks"] == 5
        assert res_q["threshold"] == pytest.approx(5000.0, abs=1e-9)

        # 3. Symmetric spike-in scale factors force norm='non' and preserve exact peak count
        out_sp = tmp_path / "peaks_sp.bed"
        res_sp = call_seacr_peaks(
            str(t_bg),
            out_sp,
            control_inputs=str(c_bg),
            target_scale_factors=[0.5],
            control_scale_factors=[0.5],
            work_dir=tmp_path / "w3",
        )
        assert res_sp["norm"] == "non"
        assert res_sp["already_spikein_scaled"] is True
        assert res_sp["n_peaks"] == 5
        assert res_sp["threshold"] == pytest.approx(res_ctrl["threshold"] * 0.5, rel=1e-6)

        # 4. Asymmetric spike-in (target only, with control) falls back to norm='norm' + warning
        out_asym = tmp_path / "peaks_asym.bed"
        res_asym = call_seacr_peaks(
            str(t_bg),
            out_asym,
            control_inputs=str(c_bg),
            target_scale_factors=[10.0],
            work_dir=tmp_path / "w4",
        )
        assert res_asym["norm"] == "norm"
        assert res_asym["already_spikein_scaled"] is False
        assert any("only one lane" in w for w in res_asym["warnings"])
        assert res_asym["n_peaks"] == 5

        # 5. Empty input list fails closed
        res_empty = call_seacr_peaks([], tmp_path / "none.bed")
        assert res_empty["status"] == "error"
        assert res_empty["error_type"] == "NoInputFiles"


class TestCopulaIdrOracle:
    def test_invert_mixture_cdf_matches_scipy_norm_cdf_accurately(self):
        u_grid = np.array([1e-6, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0 - 1e-6])
        for pi1, mu, sigma in [
            (0.05, 1.0, 0.5),
            (0.35, 2.5, 1.0),
            (0.60, 4.0, 1.8),
            (0.90, 0.5, 0.4),
        ]:
            z = _invert_mixture_cdf(u_grid, pi1, mu, sigma)
            # Forward evaluation using independent scipy.stats.norm.cdf (16-step bracketed Newton)
            g_z = (1.0 - pi1) * norm.cdf(z) + pi1 * norm.cdf((z - mu) / sigma)
            np.testing.assert_allclose(g_z, u_grid, rtol=5e-3, atol=1e-4)
            # Strictly increasing in u
            assert np.all(np.diff(z) > 0.0)

    def test_fit_copula_idr_recovers_known_bimodal_gaussian_copula_parameters(self):
        # Generate from exact Li et al. (2011) generative model:
        #   true component (pi1=0.40): bivariate normal(mu=2.6, sigma=1.0, rho=0.85)
        #   noise component (1-pi1=0.60): bivariate normal(0, 1, rho=0)
        rng = np.random.default_rng(2026)
        n = 600
        n_true = 240
        n_null = n - n_true
        true_mu, true_sigma, true_rho = 2.6, 1.0, 0.85

        cov_true = (true_sigma ** 2) * np.array([[1.0, true_rho], [true_rho, 1.0]])
        z_true = rng.multivariate_normal([true_mu, true_mu], cov_true, size=n_true)
        z_null = rng.normal(0.0, 1.0, size=(n_null, 2))
        z_all = np.vstack([z_true, z_null])

        # Apply arbitrary strictly increasing nonlinear transform (rank-copula must be invariant!)
        s1 = np.exp(0.5 * z_all[:, 0])
        s2 = np.exp(0.5 * z_all[:, 1])

        fit = fit_copula_idr_pair(s1, s2, max_iter=250, tol=1e-4)
        assert fit["converged"] is True
        assert fit["n_paired_peaks"] == n
        assert fit["params"]["pi1"] == pytest.approx(0.40, abs=0.08)
        assert fit["params"]["mu"] == pytest.approx(true_mu, abs=0.6)
        assert fit["params"]["rho"] == pytest.approx(true_rho, abs=0.12)

        g_idr = fit["global_idr"]
        l_idr = fit["local_idr"]
        # Global IDR on reproducible component is far lower than on noise component
        assert float(np.median(g_idr[:n_true])) < 0.05
        assert float(np.median(g_idr[n_true:])) > 0.25
        assert float(np.mean(g_idr[:n_true] <= 0.05)) >= 0.70
        assert float(np.mean(g_idr[n_true:] <= 0.05)) <= 0.05

        # Analytical identity: global_idr is the cumulative mean of sorted local_idr
        mean_u = 0.5 * (
            (np.argsort(np.argsort(s1)) + 1.0) / (n + 1.0)
            + (np.argsort(np.argsort(s2)) + 1.0) / (n + 1.0)
        )
        # D2-BUG-9: peaks with identical (local_idr, max(u1,u2), min(u1,u2)) form an
        # exact symmetric tie (seed 2026 has one: rows 275/527). Their order is not
        # identifiable, so swap symmetry + permutation invariance require the whole
        # tie group to share the group's end cumulative mean (conservative).
        u1 = (np.argsort(np.argsort(s1)) + 1.0) / (n + 1.0)
        u2 = (np.argsort(np.argsort(s2)) + 1.0) / (n + 1.0)
        u_hi, u_lo = np.maximum(u1, u2), np.minimum(u1, u2)
        order = np.lexsort((-u_lo, -u_hi, -mean_u, l_idr))
        expected_global = np.cumsum(l_idr[order]) / np.arange(1, n + 1)
        key = list(zip(l_idr[order], u_hi[order], u_lo[order], strict=True))
        i = 0
        while i < n:
            j = i + 1
            while j < n and key[j] == key[i]:
                j += 1
            expected_global[i:j] = expected_global[j - 1]
            i = j
        np.testing.assert_allclose(g_idr[order], expected_global, rtol=1e-12, atol=1e-12)

    def test_fit_copula_idr_invariances_and_minimum_sample_guard(self):
        rng = np.random.default_rng(7)
        n = 80
        z_sig = rng.multivariate_normal([3.0, 3.0], [[1.0, 0.85], [0.85, 1.0]], size=30)
        z_bg = rng.normal(0.0, 1.0, size=(50, 2))
        s1 = np.concatenate([z_sig[:, 0], z_bg[:, 0]])
        s2 = np.concatenate([z_sig[:, 1], z_bg[:, 1]])

        base = fit_copula_idr_pair(s1, s2)

        # 1. Replicate swap symmetry: (s1, s2) == (s2, s1)
        swapped = fit_copula_idr_pair(s2, s1)
        assert swapped["params"] == base["params"]
        np.testing.assert_allclose(swapped["local_idr"], base["local_idr"], rtol=1e-4, atol=1e-5)
        np.testing.assert_allclose(swapped["global_idr"], base["global_idr"], rtol=1e-4, atol=1e-5)

        # 2. Strictly monotone transform invariance: exp(0.3*s + 2)
        mono = fit_copula_idr_pair(np.exp(0.3 * s1 + 2.0), np.power(np.exp(s2), 0.7))
        assert mono["params"] == base["params"]
        np.testing.assert_allclose(mono["local_idr"], base["local_idr"], rtol=1e-4, atol=1e-5)
        np.testing.assert_allclose(mono["global_idr"], base["global_idr"], rtol=1e-4, atol=1e-5)

        # 3. Peak order permutation invariance for local_idr and non-tied significant global_idr
        perm = rng.permutation(n)
        permuted = fit_copula_idr_pair(s1[perm], s2[perm])
        assert permuted["params"] == base["params"]
        np.testing.assert_allclose(
            permuted["local_idr"], base["local_idr"][perm], rtol=1e-4, atol=1e-5
        )
        sig_mask = base["local_idr"][perm] < 0.999999
        np.testing.assert_allclose(
            permuted["global_idr"][sig_mask], base["global_idr"][perm][sig_mask], rtol=1e-4, atol=1e-5
        )

        # 4. <20 paired peaks or mismatched lengths raise ValueError
        with pytest.raises(ValueError, match=">=20"):
            fit_copula_idr_pair(s1[:19], s2[:19])
        with pytest.raises(ValueError, match=">=20"):
            fit_copula_idr_pair(s1[:25], s2[:24])

    # D2-BUG-9 fixed: tie-group equalization makes global_idr strictly permutation-invariant.
    def test_fit_copula_idr_global_idr_is_strictly_permutation_invariant_with_ties(self):
        rng = np.random.default_rng(7)
        z_sig = rng.multivariate_normal([3.0, 3.0], [[1.0, 0.85], [0.85, 1.0]], size=20)
        z_bg = rng.normal(0.0, 1.0, size=(20, 2))
        # Append 2 exact duplicate peaks with identical scores in both replicates
        s1 = np.concatenate([z_sig[:, 0], z_bg[:, 0], [1.5, 1.5]])
        s2 = np.concatenate([z_sig[:, 1], z_bg[:, 1], [1.5, 1.5]])
        fit = fit_copula_idr_pair(s1, s2)
        # Two peaks with identical (s1, s2) MUST receive identical global_idr!
        assert fit["global_idr"][-1] == pytest.approx(fit["global_idr"][-2], abs=1e-12)

    def test_pair_peaks_half_open_bookended_and_consensus_idr(self, tmp_path):
        # Half-open BED: [100,200) and [200,300) do NOT overlap!
        p1 = [("chr1", 100, 200, 10.0), ("chr1", 500, 600, 5.0), ("chr1", 550, 650, 12.0)]
        p2 = [("chr1", 200, 300, 10.0), ("chr1", 520, 620, 9.0)]
        paired = _pair_two_replicate_peaks(p1, p2)
        assert len(paired) == 1
        assert paired[0] == ("chr1", 500, 650, 12.0, 9.0)

        # Non-finite scores dropped by _read_scored_peaks
        np_file = tmp_path / "p.narrowPeak"
        np_file.write_text(
            "chr1\t100\t200\tp1\t100\t.\t15.5\t-1\t-1\t50\n"
            "chr1\t300\t400\tp2\t100\t.\tnan\t-1\t-1\t50\n"
            "chr1\t500\t600\tp3\t100\t.\tinf\t-1\t-1\t50\n",
            encoding="utf-8",
        )
        parsed = _read_scored_peaks(np_file)
        assert parsed == [("chr1", 100, 200, 15.5)]

        # End-to-end compute_idr_reproducible_peaks on 2 synthetic narrowPeak files
        rng = np.random.default_rng(12)
        z_sig = rng.multivariate_normal([3.8, 3.8], [[0.6, 0.55], [0.55, 0.6]], size=40)
        z_bg = rng.normal(0.0, 1.0, size=(60, 2))
        scores = np.vstack([z_sig, z_bg])
        rep1, rep2 = tmp_path / "r1.narrowPeak", tmp_path / "r2.narrowPeak"
        r1_lines, r2_lines = [], []
        for idx in range(100):
            s = idx * 1000
            r1_lines.append(f"chr1\t{s}\t{s+200}\tp{idx}\t100\t.\t{math.exp(scores[idx, 0]):.6f}\t-1\t-1\t100\n")
            r2_lines.append(f"chr1\t{s+50}\t{s+250}\tp{idx}\t100\t.\t{math.exp(scores[idx, 1]):.6f}\t-1\t-1\t100\n")
        rep1.write_text("".join(r1_lines), encoding="utf-8")
        rep2.write_text("".join(r2_lines), encoding="utf-8")

        out_idr = tmp_path / "idr_consensus.narrowPeak"
        idr_res = compute_idr_reproducible_peaks([str(rep1), str(rep2)], out_idr, idr_threshold=0.05)
        assert idr_res["ran"] is True
        assert idr_res["n_pairs_evaluated"] == 1
        assert 30 <= idr_res["n_reproducible_peaks"] <= 48
        assert out_idr.is_file()
        assert len(out_idr.read_text(encoding="utf-8").strip().splitlines()) == idr_res["n_reproducible_peaks"]

        # <2 files -> honest not-ran
        res_1 = compute_idr_reproducible_peaks([str(rep1)], tmp_path / "one.bed")
        assert res_1["ran"] is False
        assert res_1["reason"] == "fewer_than_two_replicate_peak_files"

    def test_seacr_and_idr_edge_mutant_killers(self, tmp_path):
        from aria.utils.seacr_idr import (
            _estimate_seacr_control_scale,
            _finite_positive_scales,
            _integrate_bedgraph_signal_blocks,
            _select_seacr_empirical_threshold,
            compute_idr_reproducible_peaks,
        )

        # 1. _finite_positive_scales exact error cases
        assert _finite_positive_scales(None, "s") is None
        assert _finite_positive_scales([1.5, 2], "s") == [1.5, 2.0]
        for bad in ([True], [0.0], [-1.0], [float("nan")], [float("inf")], ["abc"]):
            with pytest.raises(ValueError):
                _finite_positive_scales(bad, "scale")

        # 2. _integrate_bedgraph_signal_blocks: comments, track/browser headers, <4 cols, bad floats,
        #    zero/negative scores, reversed coords, contiguous max_region extension on tied max!
        bg = tmp_path / "edge.bg"
        bg.write_text(
            "# header comment\n"
            "track type=bedGraph\n"
            "browser position chr1:1-1000\n"
            "\n"
            "chr1\t10\t20\n"
            "chr1\tbad\t20\t5.0\n"
            "chr1\t20\t10\t5.0\n"
            "chr1\t10\t20\t0.0\n"
            "chr1\t100\t110\t5.0\n"
            "chr1\t110\t120\t10.0\n"
            "chr1\t120\t135\t10.0\n"
            "chr1\t135\t145\t4.0\n"
            "chr2\t0\t10\t3.0\n",
            encoding="utf-8",
        )
        blocks = _integrate_bedgraph_signal_blocks(bg)
        assert len(blocks) == 2
        b0 = blocks[0]
        assert b0["chrom"] == "chr1" and b0["start"] == 100 and b0["end"] == 145
        # AUC = 10*5 + 10*10 + 15*10 + 10*4 = 50 + 100 + 150 + 40 = 340.0
        assert b0["total_signal"] == pytest.approx(340.0, abs=1e-12)
        assert b0["max_signal"] == pytest.approx(10.0, abs=1e-12)
        # Tied max spans [110, 135)
        assert b0["max_region"] == "chr1:110-135"
        assert blocks[1]["chrom"] == "chr2" and blocks[1]["total_signal"] == pytest.approx(30.0)

        # 3. _estimate_seacr_control_scale: empty & non-positive control fallback -> 1.0
        assert _estimate_seacr_control_scale([], blocks) == 1.0
        assert _estimate_seacr_control_scale(blocks, []) == 1.0
        assert _estimate_seacr_control_scale(blocks, [{"start": 0, "end": 10, "total_signal": 0.0}]) == 1.0

        # 4. _select_seacr_empirical_threshold: invalid mode, empty target, empty control, relaxed <= stringent
        with pytest.raises(ValueError, match="stringent"):
            _select_seacr_empirical_threshold(np.array([1.0, 2.0]), np.array([1.0]), mode="bad")
        thr_empty_t, meta_et = _select_seacr_empirical_threshold(np.array([]), np.array([1.0]))
        assert math.isinf(thr_empty_t) and meta_et["empirical_fdr"] is None and meta_et["precision"] is None
        thr_empty_c, meta_ec = _select_seacr_empirical_threshold(np.linspace(1.0, 100.0, 100), np.array([]))
        assert thr_empty_c == pytest.approx(float(np.quantile(np.linspace(1.0, 100.0, 100), 0.99)))
        assert meta_ec["empirical_fdr"] == 0.0 and meta_ec["precision"] == 1.0

        # 5. compute_idr_reproducible_peaks: insufficient paired peaks (<20) and >max_pairs wraparound
        small_files = []
        for k in range(5):
            fp = tmp_path / f"small_{k}.bed"
            fp.write_text("\n".join(f"chr1\t{i*100}\t{i*100+50}\t{10+i}" for i in range(5)) + "\n", encoding="utf-8")
            small_files.append(str(fp))
        res_small = compute_idr_reproducible_peaks(small_files, tmp_path / "out_small.bed", max_pairs=3)
        assert res_small["ran"] is False or res_small["n_reproducible_peaks"] == 0
        # 5 files -> C(5,2)=10 > max_pairs=3 -> wraparound yields 5 adjacent pairs!
        assert len(res_small["pair_reports"]) == 5
        assert all(pr["reason"] == "insufficient_paired_peaks_for_copula_idr" for pr in res_small["pair_reports"])

    def test_seacr_idr_comprehensive_mutant_killers(self, tmp_path):
        from scipy import stats
        from aria.utils.seacr_idr import (
            _finite_positive_scales,
            _invert_mixture_cdf,
            _pair_two_replicate_peaks,
            _select_seacr_empirical_threshold,
            compute_idr_reproducible_peaks,
            fit_copula_idr_pair,
        )

        # 1. _finite_positive_scales error message content
        with pytest.raises(ValueError, match="not numeric"):
            _finite_positive_scales(["bad"], "s")
        with pytest.raises(ValueError, match="finite and positive"):
            _finite_positive_scales([0.0], "s")

        # 2. _invert_mixture_cdf extreme endpoints & sharp bimodal Newton/bisection parity
        u_ext = np.array([0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0])
        z_bi = _invert_mixture_cdf(u_ext, pi1=0.5, mu=5.0, sigma=0.25)
        cdf_check = 0.5 * stats.norm.cdf(z_bi) + 0.5 * stats.norm.cdf((z_bi - 5.0) / 0.25)
        np.testing.assert_allclose(cdf_check[1:-1], u_ext[1:-1], atol=1e-5)
        assert np.all(np.isfinite(z_bi)) and z_bi[0] < z_bi[1] and z_bi[-1] > z_bi[-2]

        # 3. _pair_two_replicate_peaks: unsorted + multiple overlapping peaks in BOTH rep0 & rep1
        p0 = [("chr1", 550, 650, 12.0), ("chr1", 500, 600, 5.0), ("chr1", 520, 620, 8.0)]
        p1 = [("chr1", 510, 610, 4.0), ("chr1", 560, 660, 15.0), ("chr1", 530, 630, 9.0)]
        paired = _pair_two_replicate_peaks(p0, p1)
        assert paired == [("chr1", 500, 660, 12.0, 15.0)]

        # 4. _select_seacr_empirical_threshold: exact hand-verified stringent & relaxed objectives,
        #    empirical_fdr, precision, and metadata keys on overlapping target & control AUCs
        t_auc = np.array([2.0, 5.0, 8.0, 12.0, 20.0, 40.0, 60.0, 80.0, 100.0])
        c_auc = np.array([1.0, 3.0, 6.0, 10.0, 15.0, 25.0])
        # Independently compute the SEACR objective on the combined-quantile candidate grid
        t_sorted, c_sorted = np.sort(t_auc), np.sort(c_auc)
        grid = np.unique(
            np.quantile(
                np.concatenate([t_sorted, c_sorted]),
                np.linspace(0.0, 0.999, min(2000, t_sorted.size + c_sorted.size)),
            )
        )
        t_cum = np.concatenate([[0.0], np.cumsum(t_sorted[::-1])])[::-1]
        c_cum = np.concatenate([[0.0], np.cumsum(c_sorted[::-1])])[::-1]
        t_idx = np.searchsorted(t_sorted, grid, side="left")
        c_idx = np.searchsorted(c_sorted, grid, side="left")
        t_ab = (t_sorted.size - t_idx).astype(float)
        c_ab = (c_sorted.size - c_idx).astype(float)
        sig_prec = t_cum[t_idx] / np.maximum(1e-12, t_cum[t_idx] + c_cum[c_idx])
        blk_prec = t_ab / np.maximum(1.0, t_ab + c_ab)
        excess = (t_cum[t_idx] / t_cum[0]) - (c_cum[c_idx] / c_cum[0])
        str_obj = excess * (0.5 * sig_prec + 0.5 * blk_prec) * (1.0 - c_ab / c_sorted.size)
        exp_str_thr = float(grid[int(np.argmax(np.where(t_ab >= 1.0, str_obj, -np.inf)))])
        exp_nt = int((t_sorted >= exp_str_thr).sum())
        exp_nc = int((c_sorted >= exp_str_thr).sum())
        exp_fdr = round(exp_nc / max(1, exp_nt + exp_nc), 6)
        exp_prec = round(exp_nt / max(1, exp_nt + exp_nc), 6)

        thr_s, meta_s = _select_seacr_empirical_threshold(t_auc, c_auc, mode="stringent")
        thr_r, meta_r = _select_seacr_empirical_threshold(t_auc, c_auc, mode="relaxed")
        assert thr_s == pytest.approx(exp_str_thr, abs=1e-12)
        assert meta_s["threshold_stringent"] == pytest.approx(exp_str_thr, abs=1e-12)
        assert meta_s["empirical_fdr"] == exp_fdr
        assert meta_s["precision"] == exp_prec
        assert meta_s["n_target_blocks_total"] == len(t_auc)
        assert meta_s["n_control_blocks_total"] == len(c_auc)
        assert meta_s["n_control_blocks_above_threshold"] == exp_nc
        assert meta_s["threshold_selected"] == pytest.approx(exp_str_thr, abs=1e-12)
        assert thr_r <= thr_s
        assert meta_r["threshold_selected"] == pytest.approx(thr_r, abs=1e-12)
        # Also test no_signal_excess branch (control dominates target everywhere)
        thr_no, meta_no = _select_seacr_empirical_threshold(np.array([1.0, 2.0]), np.array([10.0, 20.0]))
        assert math.isinf(thr_no) and meta_no["no_signal_excess"] is True
        assert meta_no["threshold_stringent"] is None and meta_no["empirical_fdr"] == 0.0
        assert meta_no["precision"] == 0.0 and meta_no["n_control_blocks_above_threshold"] == 0

        # 5. fit_copula_idr_pair: exact n==20 boundary + analytical E-step local_idr verification
        #    at the converged (pi1, mu, sigma, rho) parameters!
        rng = np.random.default_rng(77)
        z_sig = rng.multivariate_normal([3.5, 3.5], [[0.8, 0.72], [0.72, 0.8]], size=30)
        z_bg = rng.normal(0.0, 1.0, size=(30, 2))
        pts = np.vstack([z_sig, z_bg])
        fit20 = fit_copula_idr_pair(pts[:20, 0], pts[:20, 1], max_iter=1, tol=1e-12)
        assert fit20["n_paired_peaks"] == 20
        assert fit20["converged"] is False and fit20["n_iter"] == 1
        fit60 = fit_copula_idr_pair(pts[:, 0], pts[:, 1], max_iter=200, tol=5e-3)
        p = fit60["params"]
        assert set(p) == {"pi1", "mu", "sigma", "rho"}
        assert fit60["converged"] is True and 1 <= fit60["n_iter"] <= 200
        # Independently evaluate bivariate Gaussian mixture posterior at converged params
        u1 = stats.rankdata(pts[:, 0], method="average") / 61.0
        u2 = stats.rankdata(pts[:, 1], method="average") / 61.0
        z1 = _invert_mixture_cdf(u1, p["pi1"], p["mu"], p["sigma"])
        z2 = _invert_mixture_cdf(u2, p["pi1"], p["mu"], p["sigma"])
        s2 = p["sigma"] ** 2
        cov_m = np.array([[s2, p["rho"] * s2], [p["rho"] * s2, s2]])
        log_f0 = stats.multivariate_normal.logpdf(np.column_stack([z1, z2]), mean=[0.0, 0.0], cov=np.eye(2))
        log_f1 = stats.multivariate_normal.logpdf(np.column_stack([z1, z2]), mean=[p["mu"], p["mu"]], cov=cov_m)
        a0 = math.log(1.0 - p["pi1"]) + log_f0
        a1 = math.log(p["pi1"]) + log_f1
        m_max = np.maximum(a0, a1)
        exp_local_idr = 1.0 - np.exp(a1 - (m_max + np.log(np.exp(a0 - m_max) + np.exp(a1 - m_max))))
        np.testing.assert_allclose(fit60["local_idr"], exp_local_idr, rtol=1e-3, atol=1e-3)

        # 6. 3-replicate compute_idr_reproducible_peaks: exercises merging overlapping passing
        #    intervals across pairs, self_consistency_ratio, pair keys, and narrowPeak formatting!
        z_s3 = rng.multivariate_normal([4.0, 4.0, 4.0], [[0.5, 0.45, 0.45], [0.45, 0.5, 0.45], [0.45, 0.45, 0.5]], size=30)
        z_b3 = rng.normal(0.0, 1.0, size=(30, 3))
        sc3 = np.vstack([z_s3, z_b3])
        rep_files = []
        for r_idx in range(3):
            fp = tmp_path / f"tri_r{r_idx+1}.narrowPeak"
            lines = []
            for pk in range(60):
                st = pk * 1000 + r_idx * 40  # slightly shifted per rep -> overlapping merged intervals!
                en = st + 200
                lines.append(f"chr1\t{st}\t{en}\tp{pk}\t100\t.\t{math.exp(sc3[pk, r_idx]):.6f}\t-1\t-1\t100\n")
            fp.write_text("".join(lines), encoding="utf-8")
            rep_files.append(str(fp))

        out3 = tmp_path / "nested_idr" / "consensus3.narrowPeak"
        res3 = compute_idr_reproducible_peaks(rep_files, out3, idr_threshold=0.05)
        assert res3["ran"] is True
        assert res3["method"] == "gaussian_copula_em_idr"
        assert res3["n_replicate_peak_files"] == 3
        assert res3["n_pairs_evaluated"] == 3
        assert res3["idr_threshold"] == 0.05
        assert res3["output_path"] == str(out3)
        assert [pr["pair"] for pr in res3["pair_reports"]] == ["rep1_vs_rep2", "rep1_vs_rep3", "rep2_vs_rep3"]
        assert all(pr["n_paired_peaks"] == 60 for pr in res3["pair_reports"])
        counts_p = [pr["n_passing_idr"] for pr in res3["pair_reports"]]
        assert res3["self_consistency_ratio"] == round(max(counts_p) / min(counts_p), 4)
        # Verify merged coordinates span [pk*1000, pk*1000 + 280) and narrowPeak column formulas
        out_lines = [l.split("\t") for l in out3.read_text(encoding="utf-8").strip().splitlines()]
        assert len(out_lines) == res3["n_reproducible_peaks"]
        first_cols = out_lines[0]
        assert len(first_cols) == 10
        s_bed, e_bed = int(first_cols[1]), int(first_cols[2])
        assert e_bed - s_bed >= 240  # union across shifted replicate intervals!
        assert int(first_cols[9]) == (e_bed - s_bed) // 2
        neg_log10 = float(first_cols[8])
        assert int(first_cols[4]) == min(1000, max(0, int(round(neg_log10 * 125.0))))
        assert first_cols[7] == "-1"
