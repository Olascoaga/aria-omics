"""D8 convergence regression tests for contrast normalization and comparison filtering (D8-012, D8-021)."""

from __future__ import annotations


def test_d8_012_contrast_normalizers_converge():
    """D8-012: all 4 contrast normalizers must accept ARIA's full key vocabulary and A_vs_B syntax."""
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.agents.scrna.pseudobulk import PseudobulkMixin
    from aria.scripts.chromatin_bulk_diffacc import _normalise_comparisons as bulk_atac_norm
    from aria.utils.design_power import _normalize_comparisons as power_norm

    groups = {"KO": ["k1", "k2", "k3"], "WT": ["w1", "w2", "w3"]}
    variants = [
        [{"test": "KO", "ref": "WT"}],
        [{"case": "KO", "control": "WT"}],
        [{"numerator": "KO", "denominator": "WT"}],
        [{"condition": "KO", "baseline": "WT"}],
        ["KO_vs_WT"],
        [["KO", "WT"]],
    ]
    for raw in variants:
        assert len(BulkRNAAgent._normalise_explicit_contrasts(raw, groups)) == 1, raw
        assert PseudobulkMixin._normalise_pseudobulk_comparisons(raw) == [["KO", "WT"]], raw
        assert bulk_atac_norm(raw) == [{"test": "KO", "reference": "WT"}], raw
        assert power_norm(raw, groups) == [("KO", "WT")], raw

    from aria.utils.design_matrix import normalize_comparisons

    assert bulk_atac_norm is normalize_comparisons
