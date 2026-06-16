"""Guards for the multi-sample scATAC DA concordance building blocks (ADR-048
evidence): consensus-by-overlap pseudobulk assembly + DA-LFC concordance scoring.
Neutral labels only (ADR-011)."""
import numpy as np

from aria.benchmarks.consensus_pseudobulk import (
    merge_consensus_peaks,
    build_consensus_pseudobulk,
)
from aria.benchmarks.concordance_atac import score_da_lfc_concordance, _spearman


def test_merge_consensus_peaks_unions_overlap_keeps_disjoint():
    d1 = ["chr1:100-200", "chr1:500-600"]
    d2 = ["chr1:150-250", "chr2:10-20"]
    consensus = merge_consensus_peaks([d1, d2])
    # chr1:100-200 and chr1:150-250 overlap -> merged; others disjoint
    assert consensus == ["chr1:100-250", "chr1:500-600", "chr2:10-20"]


def test_merge_consensus_skips_unparseable():
    assert merge_consensus_peaks([["chr1:100-200", "not_a_peak", ""]]) == ["chr1:100-200"]


def test_build_consensus_pseudobulk_remaps_counts_by_midpoint():
    donors = {
        "d1": (["chr1:100-200", "chr1:500-600"], [10.0, 5.0]),
        "d2": (["chr1:150-250", "chr2:10-20"], [7.0, 3.0]),
    }
    mat, consensus, order = build_consensus_pseudobulk(donors)
    assert consensus == ["chr1:100-250", "chr1:500-600", "chr2:10-20"]
    assert order == ["d1", "d2"]
    # d1: chr1:100-200 (mid 150) -> consensus[0]; chr1:500-600 -> consensus[1]
    # d2: chr1:150-250 (mid 200) -> consensus[0]; chr2:10-20 -> consensus[2]
    np.testing.assert_array_equal(mat, np.array([[10.0, 7.0], [5.0, 0.0], [0.0, 3.0]]))


def test_spearman_monotone_and_ties():
    assert _spearman(np.array([1.0, 2.0, 3.0]), np.array([10.0, 20.0, 30.0])) == 1.0
    assert _spearman(np.array([1.0, 2.0, 3.0]), np.array([30.0, 20.0, 10.0])) == -1.0
    # constant vector -> no structure either way -> 1.0 (by convention)
    assert _spearman(np.array([1.0, 1.0, 1.0]), np.array([1.0, 2.0, 3.0])) == 1.0


def test_score_da_lfc_concordance_perfect_and_sets():
    aria = {"p1": 1.0, "p2": 2.0, "p3": -1.0, "p4": 0.5}
    ref = {"p1": 1.5, "p2": 2.5, "p3": -0.8, "p4": 0.3}  # same rank order + signs
    res = score_da_lfc_concordance(aria, ref, aria_sig=["p1", "p2"], ref_sig=["p2", "p3"])
    assert res["n_shared_tested"] == 4
    assert res["lfc_spearman"] == 1.0
    assert res["lfc_sign_agreement"] == 1.0
    assert res["sig_intersection"] == 1            # p2
    assert abs(res["sig_jaccard"] - 1 / 3) < 1e-9
    assert res["recall_aria_in_ref"] == 0.5
    assert res["recall_ref_in_aria"] == 0.5


def test_score_da_lfc_concordance_no_shared_is_honest_none():
    res = score_da_lfc_concordance({"p1": 1.0}, {"q1": 1.0})
    assert res["n_shared_tested"] == 0
    assert res["lfc_spearman"] is None
    assert res["lfc_sign_agreement"] is None
