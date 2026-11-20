"""Wave 2A regression tests: chromatin + multiome visualization bugs.

Covers (audit_chromatin #2, #5, #8, #9, #13, #19, #23; audit_multiome #5, #20):
motif dotplot size/colour/dedup, -log10(p) underflow caps (DA volcano, TOBIAS
volcano, p2g/eRegulon panels), correlation-heatmap text contrast, fragment-size
bands + subsample disclosure, NumPy-array depth truthiness, modality-aware FRiP /
bulk QC gates, scATAC peak annotation/ORA surfacing, singleton/constant cluster
violins, MOFA+ view-key resolution and confounder labels, WNN label column /
legend / IQR, falsy padj == 0.0, and empty-panel guards.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PathCollection, PolyCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


# ── helpers ────────────────────────────────────────────────────────────────

@pytest.fixture
def capture_chromatin(monkeypatch):
    """Capture figures passed to `_narrative_chromatin._save_dual` (still saved)."""
    import aria.agents._narrative_chromatin as nc

    captured: list = []
    real = nc._save_dual

    def _cap(fig, png_path):
        captured.append(fig)
        return real(fig, png_path)

    monkeypatch.setattr(nc, "_save_dual", _cap)
    return captured


@pytest.fixture
def capture_integration(monkeypatch):
    import aria.agents._narrative_integration as ni

    captured: list = []
    real = ni._save_dual

    def _cap(fig, png_path, dpi=220):
        captured.append(fig)
        return real(fig, png_path, dpi=60)

    monkeypatch.setattr(ni, "_save_dual", _cap)
    return captured


@pytest.fixture
def capture_closed(monkeypatch):
    """Record figures right before `plt.close(fig)` (for renderers that close
    internally without a save hook)."""
    captured: list = []
    real_close = plt.close

    def _close(fig=None):
        if fig is not None and hasattr(fig, "axes"):
            captured.append(fig)
        return real_close(fig)

    monkeypatch.setattr(plt, "close", _close)
    return captured


def _texts(fig) -> str:
    parts = [t.get_text() for t in fig.texts]
    for ax in fig.axes:
        parts += [ax.get_title(loc=loc) for loc in ("left", "center", "right")]
        parts += [ax.get_xlabel(), ax.get_ylabel()]
        parts += [t.get_text() for t in ax.texts]
        leg = ax.get_legend()
        if leg is not None:
            parts.append(leg.get_title().get_text())
            parts += [t.get_text() for t in leg.get_texts()]
        parts += [c.get_text() for c in ax.get_children()
                  if hasattr(c, "get_text") and c.__class__.__name__ == "Annotation"]
    return "\n".join(parts)


def _scatter_offsets(ax) -> np.ndarray:
    pts = [c.get_offsets() for c in ax.collections
           if isinstance(c, PathCollection) and len(c.get_offsets())]
    return np.vstack(pts) if pts else np.empty((0, 2))


# ── _narrative_chromatin: motif dotplot ─────────────────────────────────────

def test_motif_dotplot_bounded_sizes_filters_and_dedups(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_motif_dotplot

    motifs = {"per_group": {
        "0": {"top_motifs": [
            {"name": "SPI1", "log2_enrichment": 2.0, "padj": 0.0},       # underflow
            {"name": "SPI1", "log2_enrichment": 1.0, "padj": 1e-3},      # dup, weaker
            {"name": "CEBPA", "log2_enrichment": 1.5, "padj": 1e-20},
            {"name": "DEPL", "log2_enrichment": float("-inf"), "padj": 1e-50},
            {"name": "NEG", "log2_enrichment": -1.2, "padj": 1e-9},
            {"name": "INF", "log2_enrichment": float("inf"), "padj": 1e-4},
            {"name": "NANP", "log2_enrichment": 1.0, "padj": float("nan")},
            {"name": "BADP", "log2_enrichment": 1.0, "padj": 1.5},
        ]},
        "1": {"top_motifs": [
            {"name": "GATA1", "log2_enrichment": 0.8, "padj": 1e-5},
        ]},
    }}
    path = render_motif_dotplot(motifs, tmp_path / "motif.png")
    assert path and Path(path).exists() and Path(path).with_suffix(".svg").exists()
    fig = capture_chromatin[-1]
    ax = fig.axes[0]
    data_coll = [c for c in ax.collections if isinstance(c, PathCollection)
                 and len(c.get_offsets())][0]
    sizes = np.asarray(data_coll.get_sizes())
    assert sizes.min() >= 25.0 - 1e-9 and sizes.max() <= 280.0 + 1e-9
    # 3 valid (group, motif) rows: SPI1 dedup'd, CEBPA, GATA1; depleted/inf/nan dropped.
    assert len(data_coll.get_offsets()) == 3
    ylabels = [t.get_text() for t in ax.get_yticklabels()]
    assert set(ylabels) == {"SPI1", "CEBPA", "GATA1"}
    colors = np.asarray(data_coll.get_array(), dtype=float)
    assert np.all(np.isfinite(colors))
    text = _texts(fig)
    assert "-log10 padj" in text and "capped" in text   # size legend + cap disclosure
    # cap = max(20 * 1.15, 21) = 23 — never 300.
    assert "23.0" in text


def test_motif_dotplot_constant_color_and_all_invalid(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_motif_dotplot

    one = {"per_group": {"0": {"top_motifs": [
        {"name": "A", "log2_enrichment": 1.0, "padj": 0.01},
        {"name": "B", "log2_enrichment": 1.0, "padj": 0.02},
    ]}}}
    assert render_motif_dotplot(one, tmp_path / "m.png")
    coll = [c for c in capture_chromatin[-1].axes[0].collections
            if isinstance(c, PathCollection) and len(c.get_offsets())][0]
    norm = coll.norm
    assert norm.vmax > norm.vmin   # degenerate Normalize(v, v) padded

    bad = {"per_group": {"0": {"top_motifs": [
        {"name": "X", "log2_enrichment": float("-inf"), "padj": 1e-30},
        {"name": "Y", "log2_enrichment": -2.0, "padj": 0.0},
        {"name": "Z", "log2_enrichment": float("nan"), "padj": 0.01},
    ]}}}
    assert render_motif_dotplot(bad, tmp_path / "none.png") is None


def test_motif_dotplot_all_zero_padj_uses_sensible_floor(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_motif_dotplot

    motifs = {"per_group": {"0": {"top_motifs": [
        {"name": "A", "log2_enrichment": 1.0, "padj": 0.0},
        {"name": "B", "log2_enrichment": 2.0, "padj": 0.0},
    ]}}}
    assert render_motif_dotplot(motifs, tmp_path / "m.png")
    coll = [c for c in capture_chromatin[-1].axes[0].collections
            if isinstance(c, PathCollection) and len(c.get_offsets())][0]
    assert np.all(np.asarray(coll.get_sizes()) <= 280.0)


# ── _narrative_chromatin: DA volcano / MA ───────────────────────────────────

def _write_da_csv(path: Path) -> Path:
    rows = [
        # comparison, log2fc, padj, significant, base_mean
        ("A_vs_B", 2.0, 0.0, True, 50.0),
        ("A_vs_B", -1.5, 1e-10, True, 40.0),
        ("A_vs_B", 0.2, 0.5, False, 30.0),
        ("A_vs_B", float("inf"), 1e-3, True, 20.0),
        ("A_vs_B", 1.0, float("nan"), False, 10.0),
        ("A_vs_B", 1.0, -0.1, False, 10.0),
        ("A_vs_B", float("nan"), 0.01, False, 5.0),
    ]
    pd.DataFrame(rows, columns=["comparison", "log2fc", "padj", "significant",
                                "base_mean"]).to_csv(path, index=False)
    return path


def test_da_volcano_filters_nonfinite_and_caps_zero_padj(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_da_figures

    csv_path = _write_da_csv(tmp_path / "da.csv")
    figs = render_da_figures(str(csv_path), tmp_path / "out", padj_max=0.05)
    assert "da_volcano_A_vs_B" in figs and "da_ma_A_vs_B" in figs
    volcano = capture_chromatin[0]
    pts = _scatter_offsets(volcano.axes[0])
    assert np.all(np.isfinite(pts))
    assert len(pts) == 3                      # only the three valid rows survive
    ymax = float(pts[:, 1].max())
    # padj=0 placed at max(10*1.12, 10+1) = 11.2 (never 300).
    assert ymax == pytest.approx(11.2, rel=1e-6)
    assert ymax < 50
    text = _texts(volcano)
    assert "padj = 0" in text and "capped" in text
    ma = capture_chromatin[1]
    assert np.all(np.isfinite(_scatter_offsets(ma.axes[0])))


def test_da_volcano_all_invalid_rows_skip_cleanly(tmp_path):
    from aria.agents._narrative_chromatin import render_da_figures

    p = tmp_path / "bad.csv"
    pd.DataFrame({"comparison": ["c", "c"], "log2fc": [float("nan"), 1.0],
                  "padj": [0.1, float("inf")], "significant": [False, False]}
                 ).to_csv(p, index=False)
    figs = render_da_figures(str(p), tmp_path / "o")
    assert "da_volcano_c" not in figs


# ── _narrative_chromatin: correlation heatmap contrast ─────────────────────

def test_sample_correlation_text_contrast_uses_colormap_midpoint(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_sample_correlation

    rng = np.random.default_rng(0)
    base = rng.poisson(50, size=400).astype(float)
    df = pd.DataFrame({"peak_id": [f"p{i}" for i in range(400)]})
    for j in range(3):
        df[f"s{j}"] = base + rng.poisson(3, size=400)
    tsv = tmp_path / "counts.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    assert render_sample_correlation(tsv, tmp_path / "corr.png")
    ax = capture_chromatin[-1].axes[0]
    vals = [(float(t.get_text()), t.get_color()) for t in ax.texts]
    assert all(v > 0.5 for v, _ in vals)
    # Every cell is bright (>= midpoint of [0, 1]) -> black text, never white.
    assert all(c == "black" for _, c in vals)


# ── _narrative_chromatin: fragment size / TSS depth / FRiP ─────────────────

def _frag_qc(data_type=None, **fs_extra):
    edges = list(np.arange(0, 805, 5, dtype=float))
    counts = [int(100 + 50 * math.sin(i / 7.0) ** 2) for i in range(160)]
    qc = {"fragment_sizes": {"size_histogram": {"bin_edges": edges,
                                                "counts": counts}, **fs_extra}}
    if data_type:
        qc["data_type"] = data_type
    return qc


def _spans(ax):
    return [p for p in ax.patches if isinstance(p, (Rectangle,))
            or p.__class__.__name__ == "Polygon"]


def test_fragment_size_atac_bands_and_subsample_disclosure(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_fragment_size_figure

    qc = _frag_qc("scATAC", n_fragments=2_000_000, scan_truncated=True)
    assert render_fragment_size_figure(qc, tmp_path / "fs.png")
    fig = capture_chromatin[-1]
    text = _texts(fig)
    assert "subsample" in text and "2,000,000+" in text
    assert "mono-nucleosome (180-247 bp)" in text
    assert "di-nucleosome (315-473 bp)" in text
    # A true step line is drawn (not a piecewise-linear plot over a step fill).
    ax = fig.axes[0]
    assert any(line.get_drawstyle().startswith("steps") for line in ax.get_lines())


def test_fragment_size_chip_has_no_nucleosome_bands(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_fragment_size_figure

    assert render_fragment_size_figure(_frag_qc("ChIP", n_fragments_sampled=12345),
                                       tmp_path / "chip.png")
    fig = capture_chromatin[-1]
    text = _texts(fig)
    assert "nucleosome" not in text.lower().replace("fragment-size", "")
    assert "sampled" in text
    assert fig.axes[0].get_legend() is None


def test_fragment_size_cut_and_tag_uses_cutrun_bands(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_fragment_size_figure

    assert render_fragment_size_figure(_frag_qc(), tmp_path / "x.png",
                                       data_type="CUT_AND_TAG")
    text = _texts(capture_chromatin[-1])
    assert "sub-nucleosomal (<120 bp)" in text and "140-250" in text


def test_tss_depth_scatter_accepts_numpy_depth_and_masks_nonfinite(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_tss_depth_scatter

    qc = {"tss_depth_scatter": {
        "tsse": np.array([5.0, 8.0, float("nan"), 12.0]),
        "log10_depth": np.array([3.0, 3.5, 4.0, float("inf")]),
        "min_tss": 4.0, "warn_tss": 6.0}}
    path = render_tss_depth_scatter(qc, tmp_path / "tss.png")
    assert path and Path(path).exists()
    ax = capture_chromatin[-1].axes[0]
    pts = _scatter_offsets(ax)
    assert len(pts) == 2 and np.all(np.isfinite(pts))
    assert ax.get_xlabel() == "log10 fragment depth"


@pytest.mark.parametrize("data_type,gate,unit", [
    ("ChIP", 0.01, "samples"),
    ("CUT_AND_RUN", 0.05, "samples"),
    ("CUT_AND_TAG", 0.10, "samples"),
    ("bulk_ATAC", 0.20, "samples"),
    ("scATAC", 0.20, "barcodes"),
])
def test_frip_distribution_modality_gate_and_unit(tmp_path, capture_chromatin,
                                                  data_type, gate, unit):
    from aria.agents._narrative_chromatin import render_frip_distribution

    qc = {"frip_distribution": [0.03, 0.08, 0.12], "data_type": data_type}
    assert render_frip_distribution(qc, tmp_path / f"{data_type}.png")
    ax = capture_chromatin[-1].axes[0]
    assert ax.get_ylabel() == unit
    vlines = [ln.get_xdata()[0] for ln in ax.get_lines()]
    assert vlines and vlines[0] == pytest.approx(gate)


def test_frip_distribution_honors_qc_thresholds(tmp_path, capture_chromatin):
    from aria.agents._narrative_chromatin import render_frip_distribution

    qc = {"frip_distribution": [0.03, 0.08], "data_type": "ChIP",
          "thresholds": {"min_frip": 0.02, "require_tss": False}}
    assert render_frip_distribution(qc, tmp_path / "f.png")
    ax = capture_chromatin[-1].axes[0]
    assert ax.get_lines()[0].get_xdata()[0] == pytest.approx(0.02)
    assert "Per-sample" in ax.get_title()


# ── _narrative_chromatin: generate_figures wiring ──────────────────────────

def test_scatac_peak_annotation_and_ora_figures_are_surfaced(tmp_path):
    from aria.agents._narrative_chromatin import generate_figures

    ora_png = tmp_path / "scatac_peak_ora_GO.png"
    ora_png.write_bytes(b"\x89PNG")
    findings = {
        "qc": {"data_type": "scATAC", "status": "success"},
        "differential_accessibility": {"data_type": "scATAC"},
        "peak_annotation": {"ran": True, "data_type": "scATAC",
                            "feature_distribution_overall": {
                                "Promoter": 10, "Intronic": 20,
                                "Distal Intergenic": 15}},
        "peak_ora": {"ran": True, "figures": {"scatac_peak_ora_GO": str(ora_png)}},
    }
    out = generate_figures(findings, None, tmp_path / "figs", env_manager=None)
    figs = out["figures"]
    assert "scatac_peak_annotation" in figs and Path(figs["scatac_peak_annotation"]).exists()
    assert figs["scatac_peak_ora_GO"] == str(ora_png)
    assert "bulk_atac_peak_annotation" not in figs


def test_bulk_figure_params_pass_data_type_and_thresholds(tmp_path):
    from aria.agents._narrative_chromatin import _bulk_figure_params

    findings = {
        "qc": {"data_type": "CUT_AND_TAG", "frip": 0.12,
               "thresholds": {"min_frip": 0.10, "require_tss": False}},
        "peaks": {"data_type": "CUT_AND_TAG"},
    }
    params = _bulk_figure_params(findings, tmp_path)
    assert params["data_type"] == "CUT_AND_TAG"
    assert params["thresholds"]["min_frip"] == 0.10


# ── chromatin_figure_bulk ──────────────────────────────────────────────────

def test_bulk_tsse_threshold_aligned_with_chromatin_qc():
    from aria.scripts.chromatin_figure_bulk import ENCODE_TSSE_MIN
    from aria.scripts.chromatin_qc import _default_bulk_qc_thresholds

    assert ENCODE_TSSE_MIN == 6.0
    assert _default_bulk_qc_thresholds("bulk_ATAC")["min_tss"] == ENCODE_TSSE_MIN


def _bulk_qc_fig(monkeypatch, tmp_path, **kwargs):
    import aria.scripts.chromatin_figure_bulk as cfb

    captured = []
    real = cfb._save_dual

    def _cap(fig, png):
        captured.append(fig)
        return real(fig, png)

    monkeypatch.setattr(cfb, "_save_dual", _cap)
    qc = {"frip_per_sample": {"s1": 0.04, "s2": 0.06},
          "tss_enrichment_per_sample": {"s1": 3.0, "s2": 4.0},
          "dup_rate_per_sample": {"s1": 0.1, "s2": 0.2}}
    path, reason = cfb.render_qc_summary(qc, {}, tmp_path / "qc.png", **kwargs)
    assert path, reason
    return captured[-1]


def test_bulk_qc_summary_chip_uses_chip_gate_and_no_tss_gate(monkeypatch, tmp_path):
    fig = _bulk_qc_fig(monkeypatch, tmp_path, data_type="ChIP",
                       thresholds={"min_frip": 0.01, "require_tss": False})
    text = _texts(fig)
    assert "ChIP-seq" in text and "Bulk ATAC" not in text and "ENCODE" not in text
    frip_ax, tss_ax = fig.axes[0], fig.axes[1]
    frip_refs = [ln.get_ydata()[0] for ln in frip_ax.get_lines()]
    assert frip_refs == [pytest.approx(0.01)]
    assert tss_ax.get_lines() == []          # TSSe not a gate for ChIP


def test_bulk_qc_summary_atac_tss_gate_six_and_suptitle_room(monkeypatch, tmp_path):
    fig = _bulk_qc_fig(monkeypatch, tmp_path, data_type="bulk_ATAC")
    tss_ax = fig.axes[1]
    assert [ln.get_ydata()[0] for ln in tss_ax.get_lines()] == [pytest.approx(6.0)]
    # tight_layout(rect=[0,0,1,0.92]) leaves the top band to the suptitle.
    for ax in fig.axes:
        assert ax.get_position().y1 <= 0.92 + 1e-6


def test_bulk_figures_report_resolved_data_type(tmp_path):
    from aria.scripts.chromatin_figure_bulk import render_bulk_figures

    out = render_bulk_figures({"output_dir": str(tmp_path), "data_type": "CUT_AND_RUN",
                               "qc": {"frip": 0.07}})
    assert out["status"] == "success"
    assert out["data_type"] == "CUT_AND_RUN"
    assert out["thresholds"]["min_frip"] == 0.05
    assert "bulk_atac_qc_summary" in out["figures"]


def test_bulk_pca_title_uses_modality(monkeypatch, tmp_path):
    import aria.scripts.chromatin_figure_bulk as cfb

    captured = []
    monkeypatch.setattr(cfb, "_save_dual",
                        lambda fig, p: captured.append(fig) or str(p))
    rng = np.random.default_rng(1)
    mat = pd.DataFrame(rng.poisson(30, size=(50, 4)).astype(float),
                       columns=["a", "b", "c", "d"])
    path, _ = cfb.render_pca(mat, False, {}, tmp_path / "pca.png", data_type="ChIP")
    assert path
    assert captured[-1].axes[0].get_title(loc="left").startswith("ChIP-seq sample PCA")


# ── chromatin_figure_clusters ──────────────────────────────────────────────

def test_qc_depth_violin_singleton_and_constant_clusters(monkeypatch, tmp_path):
    pytest.importorskip("scipy")
    import aria.scripts.chromatin_figure_clusters as cfc

    captured = []
    real = cfc._save_dual
    monkeypatch.setattr(cfc, "_save_dual",
                        lambda fig, p: captured.append(fig) or real(fig, p))
    X = np.array([[10, 5], [20, 3], [7, 9],     # cluster "a": normal (3 cells)
                  [4, 4],                        # cluster "b": singleton
                  [5, 5], [5, 5]], dtype=float)  # cluster "c": constant depth
    obs = pd.DataFrame({"leiden": ["a", "a", "a", "b", "c", "c"]})
    adata = SimpleNamespace(X=X, obs=obs)
    path = cfc._qc_depth_violin(adata, "leiden", tmp_path / "violin.png")
    assert path and Path(path).exists()
    ax = captured[-1].axes[0]
    assert [t.get_text() for t in ax.get_xticklabels()] == ["a", "b", "c"]
    assert list(ax.get_xticks()) == [1, 2, 3]
    bodies = [c for c in ax.collections if isinstance(c, PolyCollection)]
    assert len(bodies) == 1                   # KDE only for the valid cluster


def test_qc_depth_violin_all_degenerate(monkeypatch, tmp_path):
    pytest.importorskip("scipy")
    import aria.scripts.chromatin_figure_clusters as cfc

    adata = SimpleNamespace(X=np.array([[3.0], [8.0]]),
                            obs=pd.DataFrame({"k": ["x", "y"]}))
    assert cfc._qc_depth_violin(adata, "k", tmp_path / "v.png")


# ── chromatin_footprint_tobias volcano ─────────────────────────────────────

def test_footprint_volcano_caps_zero_p_and_labels_largest_effects(tmp_path, capture_closed):
    from aria.scripts.chromatin_footprint_tobias import render_footprint_volcano

    res = tmp_path / "bindetect_results.txt"
    rows = [
        ("TF_SMALL", 0.05, 0.0), ("TF_BIG", 0.9, 0.0), ("TF_MID", 0.4, 0.0),
        ("TF_A", 0.3, 1e-8), ("TF_B", 0.2, 1e-6), ("TF_C", 0.15, 1e-4),
        ("TF_NEG", -0.5, 1e-5), ("TF_NS", 0.01, 0.6),
        ("TF_NAN", float("nan"), 0.01), ("TF_BADP", 0.3, 1.5),
    ]
    with open(res, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["name", "A_B_change", "A_B_pvalue"])
        for name, eff, p in rows:
            w.writerow([name, eff, p])
    out = render_footprint_volcano({"alpha": 0.05}, "A", "B", tmp_path / "fp",
                                   bindetect_results=str(res))
    assert out and Path(out["png"]).exists() and Path(out["svg"]).exists()
    fig = capture_closed[-1]
    ax = fig.axes[0]
    pts = _scatter_offsets(ax)
    assert np.all(np.isfinite(pts))
    assert len(pts) == 8                       # NaN effect / p>1 rows filtered
    ymax = float(pts[:, 1].max())
    assert ymax == pytest.approx(max(8 * 1.12, 9.0))
    assert ymax < 50
    annotated = [c.get_text() for c in ax.get_children()
                 if c.__class__.__name__ == "Annotation"]
    pos_labels = annotated[:5]
    # Tied p == 0 rows label by the largest |effect| first, deterministically.
    assert pos_labels[:3] == ["TF_BIG", "TF_MID", "TF_SMALL"]
    assert "capped" in _texts(fig)


# ── _narrative_integration: MOFA+ ───────────────────────────────────────────

def _mofa(view_rna="scRNA", view_atac="scATAC", r2=(0.12, 0.05), reasons=None):
    row = {"factor_id": 1, "category": "shared",
           "r2_by_view": {view_rna: r2[0], view_atac: r2[1]}}
    ex = {"factor_id": 2, "r2_by_view": {view_rna: 0.02, view_atac: 0.03},
          "excluded_from_biological_interpretation": True,
          "exclusion_reasons": reasons if reasons is not None else
          ["batch_or_donor_confounded:donor_id(eta2=0.412)"]}
    return {
        "status": "success",
        "views": {view_rna: {"kind": "rna"}, view_atac: {"kind": "atac"}},
        "top_factors": [row, ex], "excluded_factors": [ex],
        "variance_explained": {view_rna: 0.2, view_atac: 0.1},
        "variance_explained_additive": {view_rna: 0.22, view_atac: 0.11},
    }


@pytest.mark.parametrize("rna,atac", [("scRNA", "scATAC"), ("GEX", "Peaks"),
                                      ("RNA", "ATAC"), ("rna", "atac")])
def test_mofa_view_keys_resolved_case_insensitively(tmp_path, capture_integration, rna, atac):
    from aria.agents._narrative_integration import render_mofa_variance_explained

    assert render_mofa_variance_explained(_mofa(rna, atac), tmp_path)
    fig = capture_integration[-1]
    ax_hm = fig.axes[0]
    data = np.asarray(ax_hm.get_images()[0].get_array())
    assert data[0, 0] == pytest.approx(12.0) and data[0, 1] == pytest.approx(5.0)
    text = _texts(fig)
    assert "RNA=20.00%" in text and "ATAC=10.00%" in text
    assert "0.00%, ATAC=0.00%" not in text


def test_mofa_view_resolution_by_kind_only():
    from aria.agents._narrative_integration import _classify_view_kind, _view_value

    meta = {"modA": {"kind": "rna"}, "modB": {"kind": "atac"}}
    assert _classify_view_kind("modA", meta) == "rna"
    assert _view_value({"modA": 1, "modB": 2}, "atac", meta) == 2
    assert _classify_view_kind("Peaks") == "atac"
    assert _classify_view_kind("GEX") == "rna"


def test_mofa_excluded_reason_labels_are_accurate(tmp_path, capture_integration):
    from aria.agents._narrative_integration import (
        _exclusion_label, render_mofa_variance_explained,
    )

    assert _exclusion_label(["batch_or_donor_confounded:donor_id(eta2=0.412)"]) \
        == "excluded: batch/donor (donor_id)"
    assert _exclusion_label(["sequencing_depth_correlated:atac(|r|=0.812)"]) \
        == "excluded: depth |r|=0.81"
    assert _exclusion_label(["cell_cycle_correlated(|r|=0.613)",
                             "sequencing_depth_correlated:rna(|r|=0.75)"]) \
        == "excluded: cell-cycle |r|=0.61 + depth |r|=0.75"
    assert render_mofa_variance_explained(_mofa(), tmp_path)
    text = _texts(capture_integration[-1])
    assert "batch/donor (donor_id)" in text
    assert "depth/cell-cycle confounder" not in text
    assert "0.70" not in text


def test_mofa_all_zero_heatmap_has_nondegenerate_norm(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_mofa_variance_explained

    m = _mofa(r2=(0.0, 0.0))
    for row in m["top_factors"]:
        row["r2_by_view"] = {"scRNA": 0.0, "scATAC": float("nan")}
    m["excluded_factors"] = m["top_factors"][1:]
    assert render_mofa_variance_explained(m, tmp_path)
    im = capture_integration[-1].axes[0].get_images()[0]
    assert im.norm.vmax > im.norm.vmin
    assert np.all(np.isfinite(np.asarray(im.get_array())))


# ── _narrative_integration: regulatory triangulation ────────────────────────

def _p2g(**extra):
    base = {"status": "success", "n_links": 3, "top_links": [
        {"gene": "ZERO", "correlation": 0.8, "distance_bp": 5000, "padj": 0.0,
         "region_class": "distal_cis"},
        {"gene": "SMALL", "correlation": 0.4, "distance_bp": 1000, "padj": 1e-4,
         "region_class": "promoter_proximal"},
        {"gene": "NANP", "correlation": 0.5, "distance_bp": 2000, "padj": float("nan")},
        {"gene": "NANR", "correlation": float("nan"), "distance_bp": 2000, "padj": 0.01},
    ]}
    base.update(extra)
    return base


def test_parse_padj_keeps_zero_and_rejects_nonfinite():
    from aria.agents._narrative_integration import _parse_padj

    assert _parse_padj(0.0) == 0.0
    assert _parse_padj("0.0") == 0.0
    assert _parse_padj(None) == 0.05 and _parse_padj("") == 0.05
    assert _parse_padj(float("nan")) is None
    assert _parse_padj(float("inf")) is None
    assert _parse_padj(1.5) is None


def test_p2g_zero_padj_is_largest_bubble_not_shrunk(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    assert render_regulatory_triangulation(_p2g(), None, tmp_path)
    fig = capture_integration[-1]
    ax = fig.axes[0]
    colls = [c for c in ax.collections if isinstance(c, PathCollection)
             and len(c.get_offsets())]
    by_label = {c.get_label(): c for c in colls}
    distal = [c for k, c in by_label.items() if k.startswith("Distal")][0]
    prom = [c for k, c in by_label.items() if k.startswith("Promoter")][0]
    # padj == 0 (ZERO, distal) is drawn at the cap -> largest bubble.
    assert distal.get_sizes().max() > prom.get_sizes().max()
    assert len(distal.get_offsets()) + len(prom.get_offsets()) == 2  # NaN rows dropped
    text = _texts(fig)
    assert "ZERO" in text
    assert "-log10(padj)" in text and "capped" in text   # bubble-size legend


def test_p2g_truncation_is_disclosed(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    csv_path = tmp_path / "peak_gene_links.csv"
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "correlation", "distance_bp", "padj", "region_class"])
        for i in range(700):
            w.writerow([f"G{i}", 0.3 + (i % 50) / 100, i * 100, 1e-3, "distal_cis"])
    p2g = {"status": "success", "n_links": 700, "output_path": str(csv_path)}
    assert render_regulatory_triangulation(p2g, None, tmp_path)
    assert "top 600 of 700 shown" in _texts(capture_integration[-1])


def _ereg(mods, n=None):
    return {"status": "success", "n_eregulons": n if n is not None else len(mods),
            "top_eregulons": mods}


def test_eregulon_all_four_sign_classes_mapped(tmp_path, capture_integration):
    from aria.agents._narrative_integration import (
        _SIGN_CLASS_STYLE, render_regulatory_triangulation,
    )

    mods = [
        {"tf": "T1", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 1e-6},
        {"tf": "T2", "sign_class": "neg_tf_pos_p2g", "regulon_padj": 1e-5},
        {"tf": "T3", "sign_class": "pos_tf_neg_p2g", "regulon_padj": 1e-4},
        {"tf": "T4", "sign_class": "neg_tf_neg_p2g", "regulon_padj": 1e-3},
    ]
    assert render_regulatory_triangulation(None, _ereg(mods), tmp_path)
    fig = capture_integration[-1]
    ax = fig.axes[0]
    text = _texts(fig)
    for cls in _SIGN_CLASS_STYLE:
        assert cls in text
    bar_colors = [matplotlib.colors.to_hex(p.get_facecolor()) for p in ax.patches
                  if isinstance(p, Rectangle) and p.get_width() > 0]
    assert len(set(bar_colors)) == 4
    ticks = [t.get_text() for t in ax.get_yticklabels()]
    assert "T3 (TF+/p2g-)" in ticks and "T4 (TF-/p2g-)" in ticks


def test_eregulon_nan_padj_and_aucell_and_xlim_guards(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    mods = [
        {"tf": "NANP", "sign_class": "pos_tf_pos_p2g", "regulon_padj": float("nan")},
        {"tf": "ONE", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 1.0,
         "rna_atac_activity_corr": float("nan")},
        {"tf": "ZERO", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 0.0,
         "rna_atac_activity_corr": 0.42},
    ]
    assert render_regulatory_triangulation(None, _ereg(mods), tmp_path)
    ax = capture_integration[-1].axes[0]
    lo, hi = ax.get_xlim()
    assert math.isfinite(hi) and hi > lo
    ticks = [t.get_text() for t in ax.get_yticklabels()]
    assert not any(t.startswith("NANP") for t in ticks)
    text = _texts(capture_integration[-1])
    assert "AUCell r=nan" not in text and "AUCell r=NA" in text
    assert "padj=0 (capped)" in text


def test_eregulon_padj_one_only_has_positive_xlim(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    mods = [{"tf": "A", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 1.0}]
    assert render_regulatory_triangulation(None, _ereg(mods), tmp_path)
    lo, hi = capture_integration[-1].axes[0].get_xlim()
    assert hi >= 1.0


def test_eregulon_truncation_disclosed(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    mods = [{"tf": f"T{i}", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 1e-3}
            for i in range(20)]
    assert render_regulatory_triangulation(None, _ereg(mods, n=40), tmp_path)
    assert "top 12 of 40 shown" in _texts(capture_integration[-1])


def test_triangulation_no_plottable_rows_returns_none(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    p2g = {"status": "success", "n_links": 5, "top_links": [
        {"gene": "X", "correlation": float("nan"), "distance_bp": 1, "padj": 0.01}]}
    ereg = _ereg([{"tf": "Y", "regulon_padj": float("nan")}], n=3)
    n_before = len(plt.get_fignums())
    assert render_regulatory_triangulation(p2g, ereg, tmp_path) is None
    assert capture_integration == []
    assert len(plt.get_fignums()) == n_before   # no blank figure allocated


def test_triangulation_single_panel_when_other_side_empty(tmp_path, capture_integration):
    from aria.agents._narrative_integration import render_regulatory_triangulation

    ereg = _ereg([{"tf": "Y", "regulon_padj": float("nan")}], n=3)
    assert render_regulatory_triangulation(_p2g(), ereg, tmp_path)
    fig = capture_integration[-1]
    # Only the p2g panel exists — no blank eRegulon axes.
    assert len([a for a in fig.axes if a.get_title()]) == 1


# ── _narrative_integration: WNN ──────────────────────────────────────────────

def test_choose_label_column_prefers_leiden_wnn_over_donor_and_bool():
    from aria.agents._narrative_integration import _choose_label_column

    n = 12
    obs = {
        "batch": ["b1", "b2"] * 6,
        "donor_id": [f"d{i % 3}" for i in range(n)],
        "sample_id": [f"s{i % 4}" for i in range(n)],
        "pass_joint_qc": np.array([True, False] * 6),
        "is_good": ["True", "False"] * 6,
        "leiden_wnn": np.array([i % 3 for i in range(n)], dtype=np.int64),
    }
    key, labels = _choose_label_column(obs)
    assert key == "leiden_wnn"
    assert set(labels) == {"0", "1", "2"}


def test_choose_label_column_skips_non_state_generic_columns():
    from aria.agents._narrative_integration import _choose_label_column

    obs = {
        "donor_id": ["d1", "d2"] * 3,
        "mofa_training_subset": ["True", "False"] * 3,
        "atac_n_fragments": np.array([1, 2, 3, 4, 5, 6]),
        "state": ["x", "y", "z"] * 2,
    }
    assert _choose_label_column(obs)[0] == "state"
    # cell_type annotation still wins over leiden_wnn.
    obs2 = {"cell_type": ["T", "B"] * 3, "leiden_wnn": ["0", "1", "2"] * 2}
    assert _choose_label_column(obs2)[0] == "cell_type"


def _fake_bundle(n=60, n_groups=15, int_clusters=True):
    rng = np.random.default_rng(0)
    cl = np.array([i % n_groups for i in range(n)], dtype=np.int64)
    return {
        "barcodes": [f"c{i}" for i in range(n)],
        "wnn_umap": rng.normal(size=(n, 2)),
        "rna_umap": rng.normal(size=(n, 2)),
        "atac_umap": rng.normal(size=(n, 2)),
        "obs": {"leiden_wnn": cl if int_clusters else [str(x) for x in cl],
                "donor_id": [f"d{i % 2}" for i in range(n)],
                "rna_weight": rng.uniform(0.2, 0.8, size=n)},
    }


def test_wnn_umap_legend_uses_ncol(monkeypatch, tmp_path, capture_integration):
    import aria.agents._narrative_integration as ni

    monkeypatch.setattr(ni, "_read_wnn_h5mu", lambda p: _fake_bundle())
    wnn = {"status": "success", "output_path": str(tmp_path / "x.h5mu")}
    assert ni.render_wnn_umap_comparison(wnn, tmp_path)
    fig = capture_integration[-1]
    leg = fig.axes[2].get_legend()
    ncols = getattr(leg, "_ncols", getattr(leg, "_ncol", None))
    assert ncols == 2
    assert "leiden_wnn" in leg.get_title().get_text()


def test_wnn_weights_int_clusters_draw_iqr_with_legend(monkeypatch, tmp_path, capture_integration):
    import aria.agents._narrative_integration as ni

    monkeypatch.setattr(ni, "_read_wnn_h5mu", lambda p: _fake_bundle(n_groups=3))
    per_cl = {str(i): {"n_cells": 20, "mean_rna_weight": 0.5, "mean_atac_weight": 0.5}
              for i in range(3)}
    wnn = {"status": "success", "output_path": str(tmp_path / "x.h5mu"),
           "per_cluster_modality_weights": per_cl}
    assert ni.render_wnn_modality_weights(wnn, tmp_path)
    ax = capture_integration[-1].axes[0]
    # One IQR whisker line per cluster (in addition to the guide axvlines).
    whiskers = [ln for ln in ax.get_lines() if len(ln.get_xdata()) == 2
                and ln.get_ydata()[0] == ln.get_ydata()[1]
                and ln.get_xdata()[0] != ln.get_xdata()[1]]
    assert len(whiskers) == 3
    legend_text = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("Median & IQR" in t for t in legend_text)
