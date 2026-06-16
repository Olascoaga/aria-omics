"""scATAC P0 W0.4 (T15) + W0.3: pseudobulk-DA convergence gate and DA/motif figures.

W0.4 / T15: a peak is reported significant only if its padj clears the threshold AND
its apeGLM LFC fit converged; non-converged peaks (LFC≈0 with a significant Wald p)
are excluded and counted. When the convergence flag is absent (older pydeseq2 / mock)
no gating is applied (honest).

W0.3: the pseudobulk-DA volcano/MA and the motif dotplot render inline from the full
DA table / structured motif findings; non-significant/non-converged peaks are never
drawn as real effects.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")


# ── W0.4 / T15 — convergence gate ────────────────────────────────────────────

def test_convergence_gate_excludes_nonconverged_significant_peaks():
    from aria.scripts.chromatin_diffacc import _gate_da_by_convergence

    df = pd.DataFrame(
        {
            "log2FoldChange": [2.0, -1.5, 0.01, 3.0],
            "padj": [0.001, 0.01, 0.001, 0.2],
            "lfc_converged": [True, True, False, True],
        },
        index=["peakA", "peakB", "peakC", "peakD"],
    )
    sig, info = _gate_da_by_convergence(df, padj_max=0.05)

    # peakC is padj-significant but NON-converged (apeGLM LFC≈0) -> excluded.
    assert list(sig.index) == ["peakA", "peakB"]
    assert info["convergence_available"] is True
    assert info["n_nonconverged_excluded"] == 1
    # peakD is converged but not padj-significant -> simply not in the set.
    assert "peakD" not in sig.index


def test_convergence_gate_no_gating_when_flag_absent():
    from aria.scripts.chromatin_diffacc import _gate_da_by_convergence

    df = pd.DataFrame(
        {"log2FoldChange": [2.0, 0.01], "padj": [0.001, 0.001]},
        index=["peakA", "peakC"],
    )
    sig, info = _gate_da_by_convergence(df, padj_max=0.05)

    # No convergence column -> do not fabricate a verdict; keep both sig peaks.
    assert set(sig.index) == {"peakA", "peakC"}
    assert info["convergence_available"] is False
    assert info["n_nonconverged_excluded"] == 0


def test_run_deseq2_exposes_convergence_only_when_requested():
    """The shared DESeq2 core must stay convergence-agnostic by default (bulk RNA
    unchanged); the scATAC lane opts in via expose_convergence=True."""
    import inspect
    from aria.scripts.rna_bulk_de import _run_deseq2

    sig = inspect.signature(_run_deseq2)
    assert sig.parameters["expose_convergence"].default is False


# ── W0.3 — DA volcano/MA + motif dotplot (inline) ────────────────────────────

def test_render_da_figures_writes_volcano_and_ma(tmp_path):
    from aria.agents._narrative_chromatin import render_da_figures

    csv = tmp_path / "da_full.csv"
    pd.DataFrame({
        "comparison": ["treat_vs_ctrl"] * 4,
        "peak": ["p1", "p2", "p3", "p4"],
        "log2fc": [2.0, -1.8, 0.02, 0.5],
        "padj": [0.001, 0.002, 0.001, 0.4],
        "base_mean": [100.0, 80.0, 50.0, 120.0],
        "lfc_converged": [True, True, False, True],
        "significant": [True, True, False, False],
    }).to_csv(csv, index=False)

    figs = render_da_figures(str(csv), tmp_path / "figures", padj_max=0.05)
    assert "da_volcano_treat_vs_ctrl" in figs
    assert "da_ma_treat_vs_ctrl" in figs
    assert Path(figs["da_volcano_treat_vs_ctrl"]).exists()
    # P4.1: line-art figures are dual PNG + publication-grade vector SVG. The dict
    # surfaces only the PNG; the SVG sits alongside as `<stem>.svg`.
    for key in ("da_volcano_treat_vs_ctrl", "da_ma_treat_vs_ctrl"):
        png = Path(figs[key])
        assert png.suffix == ".png"
        assert png.with_suffix(".svg").exists(), f"missing vector SVG for {key}"


def test_render_da_figures_honest_empty_without_csv(tmp_path):
    from aria.agents._narrative_chromatin import render_da_figures

    assert render_da_figures(None, tmp_path) == {}
    assert render_da_figures(str(tmp_path / "nope.csv"), tmp_path) == {}


def test_render_motif_dotplot(tmp_path):
    from aria.agents._narrative_chromatin import render_motif_dotplot

    motifs = {"per_group": {
        "0": {"top_motifs": [
            {"name": "GATA1", "log2_enrichment": 2.1, "padj": 1e-5},
            {"name": "SPI1", "log2_enrichment": 1.4, "padj": 1e-3},
        ]},
        "1": {"top_motifs": [
            {"name": "GATA1", "log2_enrichment": 0.3, "padj": 0.04},
        ]},
    }}
    path = render_motif_dotplot(motifs, tmp_path / "motif_dotplot.png")
    assert path and Path(path).exists()
    # P4.1: dual PNG + vector SVG.
    assert Path(path).with_suffix(".svg").exists()
    # No enriched motifs -> honest None.
    assert render_motif_dotplot({"per_group": {}}, tmp_path / "x.png") is None


def test_render_fragment_size_is_dual_png_and_svg(tmp_path):
    """P4.1: the fragment-size / nucleosome-banding figure is dual PNG + vector SVG,
    and stays an honest None when no size histogram was computed (fragment-less run)."""
    from aria.agents._narrative_chromatin import render_fragment_size_figure

    qc = {"fragment_sizes": {"size_histogram": {
        "bin_edges": [float(i) for i in range(0, 401, 5)],
        "counts": [i % 7 for i in range(80)],
    }}}
    path = render_fragment_size_figure(qc, tmp_path / "fragment_size.png")
    assert path and Path(path).exists()
    assert Path(path).with_suffix(".svg").exists()
    # No histogram (e.g. a peak-matrix run with no fragments) -> honest None.
    assert render_fragment_size_figure({}, tmp_path / "none.png") is None


def test_generate_figures_attaches_da_and_motif_inline(tmp_path):
    from aria.agents import _narrative_chromatin

    csv = tmp_path / "da_full.csv"
    pd.DataFrame({
        "comparison": ["treat_vs_ctrl"] * 2,
        "peak": ["p1", "p2"],
        "log2fc": [2.0, -1.8],
        "padj": [0.001, 0.002],
        "base_mean": [100.0, 80.0],
        "lfc_converged": [True, True],
        "significant": [True, True],
    }).to_csv(csv, index=False)

    findings = {
        "differential_accessibility": {
            "padj_max": 0.05,
            "pseudobulk": {"full_results_csv": str(csv)},
        },
        "motifs": {"per_group": {"0": {"top_motifs": [
            {"name": "GATA1", "log2_enrichment": 2.1, "padj": 1e-5}]}}},
    }
    # No h5ad / env_manager -> UMAP skipped, but DA + motif render inline.
    out = _narrative_chromatin.generate_figures(
        findings, None, tmp_path / "figures", env_manager=None)
    figs = out["figures"]
    assert "da_volcano_treat_vs_ctrl" in figs
    assert "motif_dotplot" in figs
    assert not any(k.startswith("umap_") for k in figs)
