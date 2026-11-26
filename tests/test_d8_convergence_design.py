"""D8 convergence regression tests for contrast normalization and comparison filtering (D8-012, D8-021)."""

from __future__ import annotations


def test_d8_012_contrast_normalizers_converge():
    """D8-012: all 4 contrast normalizers must accept ARIA's full key vocabulary and A_vs_B syntax."""
    import pytest

    pytest.importorskip("litellm")
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


def test_d8_021_comparison_row_filtering_converges():
    """D8-021: filter_comparison and _select_comparison_rows must agree on space-vs-underscore and single-contrast tables."""
    import pandas as pd

    from aria.scripts.chromatin_peak_annotation import _select_comparison_rows
    from aria.utils.cross_modal_concordance import filter_comparison

    df = pd.DataFrame({
        "peak": ["chr1:100-200", "chr1:500-600"],
        "comparison": ["treat vs ctrl", "treat vs ctrl"],
        "log2FoldChange": [1.5, -1.2],
        "padj": [0.001, 0.01],
    })
    comp = {"test": "treat", "reference": "ctrl"}
    r_ann = _select_comparison_rows(df, comp)
    r_xm = filter_comparison(df, comp)
    assert len(r_ann) == 2 and len(r_xm) == 2

    # Multi-contrast table must still strictly isolate the requested contrast.
    df_multi = pd.DataFrame({
        "peak": ["p1", "p2", "p3"],
        "comparison": ["treat vs ctrl", "other vs ctrl", "treat_vs_ctrl"],
    })
    assert len(filter_comparison(df_multi, comp)) == 2
    assert len(filter_comparison(df_multi, {"test": "missing", "reference": "ctrl"})) == 0
