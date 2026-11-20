"""Wave 2C regression tests: Hi-C figures and Cross-Modal Starburst plots."""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest

from aria.scripts.hic_figure import (
    _format_resolution,
    _render_compartment_insulation_track,
    _render_differential_compartments,
    _render_ps_decay,
    hic_figure,
)
from aria.utils.cross_modal_concordance import (
    _adaptive_neglog10,
    render_starburst,
    summarize_tfs,
)


def _svg_text(svg_path: str | Path) -> str:
    return Path(svg_path).read_text(encoding="utf-8", errors="replace")


def test_format_resolution_subkb_and_noninteger_kb():
    assert _format_resolution(500) == "500 bp"
    assert _format_resolution(1000) == "1 kb"
    assert _format_resolution(2500) == "2.5 kb"
    assert _format_resolution(40000) == "40 kb"
    assert _format_resolution(1000000) == "1 Mb"


def test_hic_differential_volcano_no_300_spike_and_plots_unphased_sign_flip(tmp_path: Path):
    diff_csv = tmp_path / "diff_comp.csv"
    pd.DataFrame(
        {
            "delta_E1": [-0.8, 0.9, 0.4, 0.05, float("inf"), -0.7],
            "padj": [0.0, 1e-6, 0.01, 0.6, 0.001, float("nan")],
            "switch_class": [
                "unphased_sign_flip",
                "B_to_A",
                "A_to_B",
                "unchanged",
                "B_to_A",
                "A_to_B",
            ],
        }
    ).to_csv(diff_csv, index=False)

    bundle = _render_differential_compartments(str(diff_csv), tmp_path)
    assert bundle is not None
    assert Path(bundle["png"]).stat().st_size > 1000
    svg = _svg_text(bundle["svg"])
    assert "Sign flip (unphased)" in svg
    assert "padj = 0 capped" in svg
    # Must NOT have a y=300 axis tick.
    assert ">300<" not in svg


def test_hic_ps_decay_deduplicates_dist_and_handles_inf(tmp_path: Path):
    ps_tsv = tmp_path / "ps.tsv"
    dists = [10000, 10000, 20000, 40000, 80000, 160000, 320000, 640000, float("inf")]
    probs = [0.2, 0.18, 0.1, 0.05, 0.025, 0.012, 0.006, 0.003, 0.001]
    pd.DataFrame({"dist_bp": dists, "contact_prob": probs}).to_csv(
        ps_tsv, sep="\t", index=False
    )
    bundle = _render_ps_decay([("s1", str(ps_tsv))], tmp_path)
    assert bundle is not None
    assert Path(bundle["png"]).stat().st_size > 1000
    assert Path(bundle["svg"]).stat().st_size > 1000


def test_hic_compartment_track_single_panel_unphased_and_dense_boundaries(tmp_path: Path):
    comp_tsv = tmp_path / "comp.tsv"
    n = 120
    pd.DataFrame(
        {
            "chrom": ["chr1"] * n,
            "start": np.arange(n) * 5000,
            "end": (np.arange(n) + 1) * 5000,
            "E1": np.where(np.arange(n) % 2 == 0, 0.5, -0.4),
        }
    ).to_csv(comp_tsv, sep="\t", index=False)

    # Single-track unphased mode: should label +E1 (unphased), NOT A (+E1).
    out_single = tmp_path / "single"
    out_single.mkdir(exist_ok=True)
    bundle = _render_compartment_insulation_track(
        str(comp_tsv), None, out_single, pc1_validated=False
    )
    assert bundle is not None
    svg = _svg_text(bundle["svg"])
    assert "+E1 (unphased)" in svg
    assert "A (+E1)" not in svg

    # Dense boundaries (>80) + inf insulation values must not crash or draw wall.
    ins_tsv = tmp_path / "ins.tsv"
    ins_vals = np.sin(np.linspace(0, 6, n))
    ins_vals[0] = float("inf")
    pd.DataFrame(
        {
            "chrom": ["chr1"] * n,
            "start": np.arange(n) * 5000,
            "end": (np.arange(n) + 1) * 5000,
            "log2_insulation_score_50000": ins_vals,
            "is_boundary_50000": [True] * 100 + [False] * (n - 100),
        }
    ).to_csv(ins_tsv, sep="\t", index=False)
    out2 = tmp_path / "both"
    out2.mkdir()
    bundle2 = _render_compartment_insulation_track(
        str(comp_tsv), str(ins_tsv), out2, pc1_validated=True
    )
    assert bundle2 is not None
    svg2 = _svg_text(bundle2["svg"])
    assert "TAD boundaries (n=100)" in svg2


def test_hic_figure_isolates_corrupt_tsv_and_emits_figures_svg(tmp_path: Path):
    bad_ps = tmp_path / "empty_ps.tsv"
    bad_ps.write_text("", encoding="utf-8")
    good_diff = tmp_path / "diff.csv"
    pd.DataFrame(
        {"delta_E1": [0.5, -0.6], "padj": [0.001, 0.002], "switch_class": ["B_to_A", "A_to_B"]}
    ).to_csv(good_diff, index=False)

    res = hic_figure(
        {
            "output_dir": str(tmp_path / "figs"),
            "ps_curve_tsv": str(bad_ps),
            "differential_compartments_csv": str(good_diff),
        }
    )
    assert res["status"] == "success"
    assert res["ran"] is True
    assert "hic_differential_compartments" in res["figures"]
    assert "hic_differential_compartments" in res["figures_svg"]
    assert Path(res["figures_svg"]["hic_differential_compartments"]).is_file()
    assert "hic_contact_decay_ps" in res["skipped"]


def test_hic_contact_map_mouse_telocentric_gap_and_unbalanced_label(tmp_path: Path):
    cooler = pytest.importorskip("cooler")
    from aria.scripts.hic_figure import _render_contact_map

    # Simulate mouse chr1 with 0-3 Mb telocentric N-gap (all NaN weights and 0 contacts)
    # and mappable contacts starting at 3.0 Mb (bin 300 at 10 kb resolution).
    binsize = 10_000
    n_bins = 400
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * n_bins,
            "start": np.arange(n_bins) * binsize,
            "end": (np.arange(n_bins) + 1) * binsize,
            "weight": [np.nan] * 300 + [0.05] * 100,
        }
    )
    pix_rows = []
    for i in range(300, 395):
        for d in range(0, 6):
            pix_rows.append((i, i + d, 20 - d * 2))
    pixels = pd.DataFrame(pix_rows, columns=["bin1_id", "bin2_id", "count"])
    cool_path = tmp_path / "mouse_gap.cool"
    cooler.create_cooler(str(cool_path), bins, pixels)

    out_bal = tmp_path / "bal"
    out_bal.mkdir(exist_ok=True)
    bundle = _render_contact_map(str(cool_path), binsize, None, out_bal, max_bins=100)
    assert bundle is not None
    svg = _svg_text(bundle["svg"])
    assert "ICE-Balanced Contact Map" in svg
    assert "10 kb" in svg

    # Now test all-NaN weight fallback -> must NOT claim ICE-Balanced.
    bins_unbal = bins.copy()
    bins_unbal["weight"] = np.nan
    cool_unbal = tmp_path / "mouse_unbal.cool"
    cooler.create_cooler(str(cool_unbal), bins_unbal, pixels)
    out_raw = tmp_path / "raw"
    out_raw.mkdir(exist_ok=True)
    bundle_raw = _render_contact_map(str(cool_unbal), binsize, None, out_raw, max_bins=100)
    assert bundle_raw is not None
    svg_raw = _svg_text(bundle_raw["svg"])
    assert "Unbalanced Fallback" in svg_raw
    assert "ICE-Balanced" not in svg_raw


def test_cross_modal_starburst_adaptive_caps_casefold_and_empty_fallbacks(tmp_path: Path):
    vals, cap, capped = _adaptive_neglog10([0.0, 1e-4, 1e-2])
    assert capped is True
    assert 4.0 < cap < 15.0
    assert vals[0] == pytest.approx(cap)

    # Case-folded TF deduplication (FOXA1 and Foxa1 collapse to 1 evaluated TF).
    tfs = pd.DataFrame(
        [
            {
                "tf": "FOXA1",
                "chromatin_score": None,
                "chromatin_direction": "up",
                "chromatin_padj": 0.0,
                "rna_log2fc": 2.1,
                "classification": "concordant",
                "paralog_motif_ambiguity": False,
            },
            {
                "tf": "Foxa1",
                "chromatin_score": 3.0,
                "chromatin_direction": "up",
                "chromatin_padj": 1e-4,
                "rna_log2fc": 2.1,
                "classification": "concordant",
                "paralog_motif_ambiguity": False,
            },
        ]
    )
    tf_sum = summarize_tfs(tfs)
    assert tf_sum["n_tfs_evaluated"] == 1

    pathways = pd.DataFrame(
        [
            {
                "safe_label": "Circadian Rhythm",
                "rna_padj": 0.0,
                "atac_padj": 1e-5,
                "geometric_mean_padj": 0.0,
            }
        ]
    )
    # Non-empty de_da without finite LFCs + fisher_pvalue == 0.0.
    de_da_no_lfc = pd.DataFrame(
        [{"gene": "G1", "gene_key": "G1", "rna_log2fc": np.nan, "atac_log2fc": np.nan,
          "rna_padj": 0.01, "atac_padj": 0.01, "classification": "concordant_up"}]
    )
    png, svg = render_starburst(
        de_da_no_lfc,
        tfs,
        pathways,
        {"spearman_r": None, "fisher_pvalue": 0.0, "concordance_fraction": 1.0},
        tmp_path,
    )
    assert Path(png).is_file() and Path(svg).is_file()
    svg_txt = _svg_text(svg)
    assert "Fisher p &lt; 1e-300" in svg_txt or "Fisher p < 1e-300" in svg_txt
    assert "padj = 0 capped" in svg_txt
    assert "* Motif-family or dimer ambiguity" not in svg_txt
    assert ">300<" not in svg_txt

    # Empty pathways AND zero quadrant counts -> clean unavailable notice.
    out_empty = tmp_path / "empty"
    out_empty.mkdir()
    _, svg_empty = render_starburst(
        pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {}, out_empty
    )
    assert "No shared pathways or quadrant counts available" in _svg_text(svg_empty)
