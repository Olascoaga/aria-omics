"""Analytical oracle tests for aria/benchmarks/calibration.py and baseline scoring.

Tests:
  - Exact boundary conditions (recall == min_recall, empirical_fdr == max_empirical_fdr,
    false_positive_fraction == max_false_positive_fraction, n_called == 0)
  - Strict fail-closed on missing, None, non-numeric, NaN, inf, and bool metrics
  - Unknown benchmark rejection
  - Calibration baseline JSON schema & positive finite bounds
  - Analytical scoring of synthetic_de / synthetic_atac_da with planted truth
"""
from __future__ import annotations

import math

import pytest

from aria.benchmarks.calibration import (
    baseline_bounds,
    check_calibration,
    load_baseline,
)
from aria.benchmarks.synthetic_atac_da import (
    run_atac_da_benchmark,
    simulate_atac_da_dataset,
)


class TestCalibrationBaselineSchema:
    def test_baseline_loads_known_benchmarks_with_valid_bounds(self):
        raw = load_baseline()
        assert raw["baseline_version"] == "1"
        benchmarks = raw["benchmarks"]
        for name in ("synthetic_de_pseudobulk", "scatac_pseudobulk_da", "bulk_atac_da"):
            assert name in benchmarks
            b = baseline_bounds(name)
            assert 0.0 < b["min_recall"] <= 1.0
            assert 0.0 <= b["max_empirical_fdr"] < 1.0
            assert 0.0 <= b["max_false_positive_fraction"] < 1.0

    def test_unknown_benchmark_raises_key_error(self):
        with pytest.raises(KeyError, match="no calibration baseline"):
            baseline_bounds("nonexistent_benchmark_xyz")
        with pytest.raises(KeyError):
            check_calibration("nonexistent_benchmark_xyz", {"recall": 1.0, "empirical_fdr": 0.0})


class TestCheckCalibrationExactBoundaries:
    """Exact threshold checks on recall, empirical_fdr, and false_positive_fraction."""

    def test_exact_boundary_values_pass(self):
        b = baseline_bounds("synthetic_de_pseudobulk")
        res = check_calibration(
            "synthetic_de_pseudobulk",
            {
                "recall": b["min_recall"],
                "empirical_fdr": b["max_empirical_fdr"],
                "n_false_positive": 5,
                "n_called": 100,  # 5/100 == 0.05 == max_false_positive_fraction
            },
        )
        assert res["status"] == "pass"
        assert res["violations"] == []
        assert res["measured"]["recall"] == pytest.approx(b["min_recall"])
        assert res["measured"]["empirical_fdr"] == pytest.approx(b["max_empirical_fdr"])
        assert res["measured"]["false_positive_fraction"] == pytest.approx(0.05)
        assert "provenance" not in res["baseline"]
        assert res["baseline"]["min_recall"] == b["min_recall"]

    def test_epsilon_below_min_recall_fails(self):
        b = baseline_bounds("synthetic_de_pseudobulk")
        res = check_calibration(
            "synthetic_de_pseudobulk",
            {"recall": b["min_recall"] - 1e-6, "empirical_fdr": 0.0},
        )
        assert res["status"] == "fail"
        assert len(res["violations"]) == 1
        assert "recall" in res["violations"][0]

    def test_epsilon_above_max_empirical_fdr_fails(self):
        b = baseline_bounds("synthetic_de_pseudobulk")
        res = check_calibration(
            "synthetic_de_pseudobulk",
            {"recall": 1.0, "empirical_fdr": b["max_empirical_fdr"] + 1e-6},
        )
        assert res["status"] == "fail"
        assert len(res["violations"]) == 1
        assert "empirical_fdr" in res["violations"][0]

    def test_false_positive_fraction_exceeding_bound_fails_and_zero_called_safe(self):
        res_fail = check_calibration(
            "bulk_atac_da",
            {
                "recall": 0.9,
                "empirical_fdr": 0.06,
                "n_false_positive": 6,
                "n_called": 100,
            },
        )
        assert res_fail["status"] == "fail"
        assert any("false-positive fraction" in v for v in res_fail["violations"])
        assert res_fail["measured"]["false_positive_fraction"] == pytest.approx(0.06)

        # n_called == 0 with n_fp == 0 -> denom max(0, 1) = 1 -> frac = 0.0
        res_zero = check_calibration(
            "bulk_atac_da",
            {
                "recall": 0.8,
                "empirical_fdr": 0.0,
                "n_false_positive": 0,
                "n_called": 0,
            },
        )
        assert res_zero["status"] == "pass"
        assert res_zero["measured"]["false_positive_fraction"] == 0.0

        # Missing n_false_positive / n_called skips optional fp-fraction check
        res_opt = check_calibration(
            "bulk_atac_da",
            {"recall": 0.8, "empirical_fdr": 0.0},
        )
        assert res_opt["status"] == "pass"
        assert "false_positive_fraction" not in res_opt["measured"]

    def test_missing_none_and_non_numeric_metrics_fail_closed(self):
        res_missing = check_calibration("scatac_pseudobulk_da", {})
        assert res_missing["status"] == "fail"
        assert len(res_missing["violations"]) == 2
        assert any("missing metric 'recall'" in v for v in res_missing["violations"])
        assert any("missing metric 'empirical_fdr'" in v for v in res_missing["violations"])

        res_none = check_calibration(
            "scatac_pseudobulk_da", {"recall": None, "empirical_fdr": None}
        )
        assert res_none["status"] == "fail"
        assert len(res_none["violations"]) == 2

        res_str = check_calibration(
            "scatac_pseudobulk_da", {"recall": "not_a_number", "empirical_fdr": "bad"}
        )
        assert res_str["status"] == "fail"
        assert len(res_str["violations"]) == 2
        assert all("non-numeric metric" in v for v in res_str["violations"])

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "D2-BUG-3: calibration.check_calibration accepts NaN metrics as 'pass' "
            "because NaN < min_recall and NaN > max_empirical_fdr both evaluate False"
        ),
    )
    def test_nan_metrics_must_fail_closed(self):
        res = check_calibration(
            "synthetic_de_pseudobulk",
            {
                "recall": float("nan"),
                "empirical_fdr": float("nan"),
                "n_false_positive": float("nan"),
                "n_called": 100,
            },
        )
        assert res["status"] == "fail"
        assert len(res["violations"]) >= 2

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "D2-BUG-4: calibration.check_calibration accepts boolean True/False as "
            "numeric metrics (float(True)==1.0, float(False)==0.0) instead of failing closed"
        ),
    )
    def test_boolean_metrics_must_not_pass_as_numeric(self):
        res = check_calibration(
            "synthetic_de_pseudobulk",
            {"recall": True, "empirical_fdr": False},
        )
        assert res["status"] == "fail"


class TestSyntheticBenchmarkScoringOracle:
    """Analytical oracle checks on synthetic_atac_da scoring against exact TP/FP counts."""

    def test_atac_da_benchmark_exact_recall_and_fdr_math(self):
        ds = simulate_atac_da_dataset(
            n_peaks=100,
            n_da=20,
            n_cells_per_condition=30,
            n_replicates_per_condition=3,
            seed=7,
        )
        assert ds.adata.n_obs == 60
        assert ds.adata.n_vars == 100
        assert ds.n_da == 20
        assert len(ds.null_peaks) == 80
        assert set(ds.da_peaks).isdisjoint(set(ds.null_peaks))

        # Exact known caller: 15 of 20 true DA + 5 of 80 null -> recall=0.75, fdr=5/20=0.25
        da_list = sorted(ds.da_peaks)
        called_set = set(da_list[:15]) | set(ds.null_peaks[:5])
        res = run_atac_da_benchmark(
            ds,
            da_fn=lambda _d: called_set,
            min_recall=0.70,
            max_empirical_fdr=0.20,
        )
        assert res.n_true_de == 20
        assert res.n_called == 20
        assert res.n_true_positive == 15
        assert res.n_false_positive == 5
        assert res.recall == pytest.approx(0.75, abs=1e-12)
        assert res.empirical_fdr == pytest.approx(0.25, abs=1e-12)
        assert res.status == "fail"  # 0.25 > 0.20

        # Perfect caller: recall=1.0, fdr=0.0 -> pass
        res_perfect = run_atac_da_benchmark(
            ds,
            da_fn=lambda _d: set(ds.da_peaks),
            min_recall=0.5,
            max_empirical_fdr=0.2,
        )
        assert res_perfect.recall == 1.0
        assert res_perfect.empirical_fdr == 0.0
        assert res_perfect.status == "pass"

        # Empty caller: recall=0.0, fdr=0.0 -> fail
        res_empty = run_atac_da_benchmark(
            ds, da_fn=lambda _d: set(), min_recall=0.5, max_empirical_fdr=0.2
        )
        assert res_empty.recall == 0.0
        assert res_empty.empirical_fdr == 0.0
        assert res_empty.status == "fail"

        # Missing da_fn fails closed with NotImplementedError
        with pytest.raises(NotImplementedError):
            run_atac_da_benchmark(ds, da_fn=None)

    def test_atac_da_simulator_remainder_cell_distribution_exact(self):
        # 80 cells across 3 replicates -> 27, 27, 26 = 80 per condition, 160 total.
        ds = simulate_atac_da_dataset(
            n_peaks=30,
            n_da=5,
            n_cells_per_condition=80,
            n_replicates_per_condition=3,
            seed=19,
        )
        assert ds.adata.n_obs == 160
        cond_counts = ds.adata.obs["condition"].value_counts().to_dict()
        assert cond_counts == {"COND_A": 80, "COND_B": 80}
        assert math.isfinite(float(ds.adata.X.sum()))

    def test_check_calibration_optional_fp_fraction_and_baseline_sanitization(self):
        # Kills remaining mutants in check_calibration: n_called=0 clamping max(0,1)=1,
        # rounding to 4 decimals, optional skip when either n_fp or n_called is None,
        # and stripping 'provenance' from the returned baseline dict.
        bounds = baseline_bounds("bulk_atac_da")
        max_fp_frac = bounds["max_false_positive_fraction"]

        # Both counts omitted -> optional check skipped, measured lacks false_positive_fraction.
        r_opt = check_calibration(
            "bulk_atac_da",
            {"recall": bounds["min_recall"], "empirical_fdr": bounds["max_empirical_fdr"]},
        )
        assert r_opt["status"] == "pass"
        assert "false_positive_fraction" not in r_opt["measured"]
        assert "provenance" not in r_opt["baseline"]
        assert r_opt["benchmark"] == "bulk_atac_da"

        # Only ONE of n_false_positive / n_called supplied -> still skipped (requires BOTH).
        for partial in ({"n_false_positive": 5}, {"n_called": 10}):
            r_half = check_calibration(
                "bulk_atac_da",
                {
                    "recall": bounds["min_recall"],
                    "empirical_fdr": bounds["max_empirical_fdr"],
                    **partial,
                },
            )
            assert r_half["status"] == "pass"
            assert "false_positive_fraction" not in r_half["measured"]

        # n_called == 0 with n_false_positive == 0 -> denom clamped to 1.0, frac == 0.0
        r_zero = check_calibration(
            "bulk_atac_da",
            {
                "recall": bounds["min_recall"],
                "empirical_fdr": 0.0,
                "n_false_positive": 0,
                "n_called": 0,
            },
        )
        assert r_zero["status"] == "pass"
        assert r_zero["measured"]["false_positive_fraction"] == 0.0

        # Exact rounding to 4 decimals: 1/30 = 0.033333... -> 0.0333
        r_round = check_calibration(
            "bulk_atac_da",
            {
                "recall": bounds["min_recall"],
                "empirical_fdr": 0.0,
                "n_false_positive": 1,
                "n_called": 30,
            },
        )
        assert r_round["measured"]["false_positive_fraction"] == 0.0333
        assert (r_round["status"] == "pass") == (1 / 30 <= max_fp_frac)

        # Violation message contains exact formatted numbers
        r_fail = check_calibration(
            "bulk_atac_da",
            {
                "recall": bounds["min_recall"],
                "empirical_fdr": 0.0,
                "n_false_positive": 1,
                "n_called": 0,  # 1 / max(0, 1) = 1.0 > max_fp_frac
            },
        )
        assert r_fail["status"] == "fail"
        assert r_fail["measured"]["false_positive_fraction"] == 1.0
        assert any("false-positive fraction 1.0000 >" in v for v in r_fail["violations"])
