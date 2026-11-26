"""D8 convergence regression tests for bool flags, categorical palette, and MEME counter (D8-023, D8-025, CL-028)."""

from __future__ import annotations

from pathlib import Path
import pytest


def test_d8_023_boolean_flag_parsers_converge():
    """D8-023: all 4 boolean-flag parsers must agree on 0.0, NaN, 'n', and 'on'."""
    pytest.importorskip("litellm")
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.scripts.chromatin_lsi_clustering import _as_bool as lsi_bool
    from aria.scripts.chromatin_peak_counts import _parse_bool_flag as cnt_bool
    from aria.scripts.rna_celltypist import _as_bool as ct_bool

    expected = [
        (0.0, False),
        (float("nan"), False),
        (float("inf"), False),
        ("n", False),
        ("on", True),
        ("yes", True),
        ("1", True),
        (True, True),
        (False, False),
        (None, False),
    ]
    for val, want in expected:
        outs = {
            "chromatin_lsi_clustering._as_bool": lsi_bool(val),
            "rna_celltypist._as_bool": ct_bool(val),
            "chromatin_peak_counts._parse_bool_flag": cnt_bool(val),
            "chromatin_agent._parse_bool_flag": ChromatinAgent._parse_bool_flag(val),
        }
        assert all(v is want for v in outs.values()), f"Bool parsers disagree on {val!r}: {outs}"

    from aria.utils.thresholds import parse_bool_flag

    assert lsi_bool is parse_bool_flag
    assert ct_bool is parse_bool_flag
    assert cnt_bool is parse_bool_flag


def test_d8_025_scrna_categorical_colors_converge_above_40():
    """D8-025: UMAP and PAGA categorical color palettes must match across >40 clusters."""
    from aria.scripts.rna_figure_paga import _categorical_colors as paga_colors
    from aria.scripts.rna_figure_umap import _categorical_colors as umap_colors

    c_umap = umap_colors(45)
    c_paga = paga_colors(45)
    assert len(set(c_umap)) == 45
    assert c_umap == c_paga
    assert paga_colors is umap_colors


def test_cl_028_meme_motif_counter_converges(tmp_path: Path):
    """CL-028: motifs._count_meme_motifs and chromatin_motifs._count_meme_motifs must share strict MOTIF matching."""
    from aria.scripts.chromatin_motifs import _count_meme_motifs as script_count
    from aria.utils.motifs import _count_meme_motifs as util_count

    meme = tmp_path / "test.meme"
    meme.write_text(
        "MEME version 4\n"
        "MOTIF_HEADER_COMMENT should not be counted\n"
        "MOTIF MA0001.1\tAGL3\n"
        "MOTIF\tMA0002.2 RUNX1\n",
        encoding="utf-8",
    )
    assert util_count(meme) == 2
    assert script_count(str(meme)) == 2
    assert util_count(tmp_path / "nonexistent.meme") is None
    assert script_count is util_count
