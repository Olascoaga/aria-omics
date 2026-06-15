"""scATAC P0 W0.2: fragment-size/nucleosome figure (inline), marker-peak heatmap
and per-cluster QC depth violin (chromatin stack). No fabrication: a panel appears
only when its inputs exist.
"""

from __future__ import annotations

from pathlib import Path

import pytest


# ── Fragment-size / nucleosome banding (inline) ──────────────────────────────

def test_fragment_size_figure_renders_from_histogram(tmp_path):
    from aria.agents._narrative_chromatin import render_fragment_size_figure

    edges = [float(i) for i in range(0, 805, 5)]   # 0..800, 5 bp bins -> 160 bins
    counts = [0] * 160
    counts[30] = 500   # ~150 bp sub-nucleosome peak
    counts[40] = 300   # ~200 bp mono-nucleosome
    qc = {"fragment_sizes": {"size_histogram": {"bin_edges": edges,
                                                "counts": counts}}}
    path = render_fragment_size_figure(qc, tmp_path / "fragment_size.png")
    assert path and Path(path).exists()


def test_fragment_size_figure_honest_none_without_histogram(tmp_path):
    from aria.agents._narrative_chromatin import render_fragment_size_figure

    # .h5mu peak-matrix run: no fragments -> no histogram -> no figure.
    assert render_fragment_size_figure({}, tmp_path / "x.png") is None
    assert render_fragment_size_figure(
        {"fragment_sizes": {"median_size": 180}}, tmp_path / "y.png") is None
    assert render_fragment_size_figure(
        {"fragment_sizes": {"size_histogram": {"bin_edges": [0, 5],
                                               "counts": [0]}}},
        tmp_path / "z.png") is None


def test_qc_emits_size_histogram_in_fragment_sizes():
    """The QC fragment-size computation must emit a binned histogram (W0.2), not
    only summary stats, so the report can plot the banding pattern."""
    import inspect
    from aria.scripts import chromatin_qc

    src = inspect.getsource(chromatin_qc._compute_fragment_sizes)
    assert "size_histogram" in src


# ── Marker-peak heatmap extraction + generate_figures wiring ─────────────────

def test_marker_peaks_by_cluster_extraction():
    from aria.agents._narrative_chromatin import _marker_peaks_by_cluster

    findings = {"differential_accessibility": {"per_cluster": {
        "da_peaks_by_cluster": {
            "0": [{"peak": "chr1:1-2"}, {"peak": "chr1:3-4"}],
            "1": [{"peak": "chr2:5-6"}],
        }}}}
    out = _marker_peaks_by_cluster(findings)
    assert out == {"0": ["chr1:1-2", "chr1:3-4"], "1": ["chr2:5-6"]}


class _FakeEM:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, stack, script_path, params):
        self.calls.append({"stack": stack, "script_path": script_path,
                           "params": params})
        out = {}
        for name in ("marker_heatmap", "qc_depth_violin"):
            p = Path(params["output_dir"]) / f"{name}.png"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"\x89PNG")
            out[name] = str(p)
        return {"status": "success", "figures": out}


def test_generate_figures_renders_cluster_figures_in_chromatin_stack(tmp_path):
    from aria.agents import _narrative_chromatin

    em = _FakeEM()
    findings = {
        "lsi": {"output_path": "/clustered.h5ad"},
        "differential_accessibility": {"per_cluster": {"da_peaks_by_cluster": {
            "0": [{"peak": "chr1:1-2"}]}}},
    }
    out = _narrative_chromatin.generate_figures(
        findings, "/clustered.h5ad", tmp_path / "figures", env_manager=em)
    figs = out["figures"]
    assert "marker_heatmap" in figs
    assert "qc_depth_violin" in figs
    # The cluster-figure subprocess ran in the chromatin stack with the script
    # and the extracted marker peaks.
    cluster_call = next(c for c in em.calls
                        if c["script_path"].endswith("chromatin_figure_clusters.py"))
    assert cluster_call["stack"] == "chromatin"
    assert cluster_call["params"]["marker_peaks"] == {"0": ["chr1:1-2"]}


def test_generate_figures_skips_cluster_figures_without_env(tmp_path):
    from aria.agents import _narrative_chromatin

    findings = {"lsi": {"output_path": "/clustered.h5ad"}}
    out = _narrative_chromatin.generate_figures(
        findings, "/clustered.h5ad", tmp_path / "figures", env_manager=None)
    assert "marker_heatmap" not in out["figures"]
    assert "qc_depth_violin" not in out["figures"]
