"""Wave 2B regressions: inspect plotted data as well as exported artifacts."""
from __future__ import annotations

from pathlib import Path
import importlib
import importlib.util
import sys
import types

# The science-only conda environment runs pytest without an editable ARIA install.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.text import Annotation
import numpy as np
import pandas as pd
import pytest

from aria.scripts.rna_bulk import plots as bulk
from aria.scripts import rna_pathway_viz as pathways
from aria.scripts import rna_figure_paga as paga
from aria.scripts import rna_figure_umap as umap


def _load_narrative_figures():
    """Load real leaf renderers without the optional agent/LLM startup stack."""
    if importlib.util.find_spec("litellm") is not None:
        return importlib.import_module("aria.agents.narrative.scrna.figures")
    # Namespace shells avoid executing heavyweight package __init__ files.
    # Every renderer and helper is still imported from the repository source.
    root = Path(__file__).resolve().parents[1]
    prefix = "aria.agents"
    previous = {name: module for name, module in sys.modules.items()
                if name == prefix or name.startswith(prefix + ".")}
    try:
        for name in ("aria.agents", "aria.agents.narrative", "aria.agents.narrative.scrna"):
            package = types.ModuleType(name)
            package.__path__ = [str(root.joinpath(*name.split(".")))]
            sys.modules[name] = package
        return importlib.import_module("aria.agents.narrative.scrna.figures")
    finally:
        for name in list(sys.modules):
            if name == prefix or name.startswith(prefix + "."):
                del sys.modules[name]
        sys.modules.update(previous)


narrative = _load_narrative_figures()


@pytest.fixture
def saved_figures(monkeypatch):
    """Keep artist references after close, while exercising real PNG/SVG export."""
    saved = {}
    original = Figure.savefig

    def capture(fig, filename, *args, **kwargs):
        saved[Path(filename).name] = fig
        return original(fig, filename, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", capture)
    yield saved
    plt.close("all")


def _volcano(tmp_path, frame, symbol_map=None):
    result = bulk._generate_plots(
        {"results": frame}, {"_symbol_map": symbol_map}, None,
        pd.DataFrame(), "condition", str(tmp_path), 0.05, 1.0,
    )
    assert "error" not in result
    assert Path(result["volcano"]).is_file()


def _offsets(ax):
    return np.concatenate([np.asarray(c.get_offsets()) for c in ax.collections])


def test_volcano_filters_invalid_values_and_aligns_symbol_annotations(tmp_path, saved_figures):
    frame = pd.DataFrame({
        "log2FoldChange": [20, -12, 2, np.inf, -np.inf, np.nan, 1, 1, 1, 1],
        "padj": [0, 1e-8, 0.5, 0.01, 0.01, 0.01, -0.1, 1.1, np.inf, np.nan],
    }, index=["ENSG00000123456.1", "ENSG00000123457", *[f"G{i}" for i in range(8)]])
    _volcano(tmp_path, frame, {"ENSG00000123456": "SYMBOL_A", "ENSG00000123457": "SYMBOL_B"})
    ax = saved_figures["volcano.svg"].axes[0]
    offsets = _offsets(ax)
    assert len(offsets) == 3
    assert np.isfinite(offsets).all()
    assert offsets[:, 1].max() == pytest.approx(9.2)
    annotations = [t for t in ax.texts if isinstance(t, Annotation)]
    assert {t.get_text() for t in annotations} == {"SYMBOL_A", "SYMBOL_B"}
    for annotation in annotations:
        assert any(np.allclose(annotation.xy, point) for point in offsets)
    assert any("underflow" in t.get_text() for t in ax.texts)


def test_volcano_keeps_distinct_full_ensembl_labels(tmp_path, saved_figures):
    ids = ["ENSG00000123456", "ENSG00000123457", "LONG_SYMBOL_" * 5]
    frame = pd.DataFrame({"log2FoldChange": [1, 2, 3], "padj": [0.001, 0.002, 0.003]}, index=ids)
    _volcano(tmp_path, frame)
    labels = [t.get_text() for t in saved_figures["volcano.svg"].axes[0].texts if isinstance(t, Annotation)]
    assert ids[0] in labels and ids[1] in labels
    assert any(label.endswith("…") and len(label) <= 40 for label in labels)


@pytest.mark.parametrize("padj,expected_max", [([0, 0], 10), ([1e-100, 0], 115)])
def test_volcano_zero_cap_is_adaptive_without_clipping_positive_pvalues(tmp_path, saved_figures, padj, expected_max):
    _volcano(tmp_path, pd.DataFrame({"log2FoldChange": [2, -2], "padj": padj}, index=["G1", "G2"]))
    points = _offsets(saved_figures["volcano.svg"].axes[0])
    assert points[:, 1].max() == pytest.approx(expected_max)
    if padj[0] > 0:
        assert 100 in points[:, 1]


def test_heatmap_failure_logs_and_returns_none(tmp_path, monkeypatch, caplog):
    def fail(*args, **kwargs):
        raise ValueError("synthetic rendering failure")

    monkeypatch.setattr(Axes, "imshow", fail)
    frame = pd.DataFrame([[1, 2]], index=["G1"])
    assert bulk._plot_heatmap(["G1"], frame, None, {}, "Top genes", "", str(tmp_path / "heatmap.svg")) is None
    assert "Heatmap failed" in caplog.text
    plt.close("all")


@pytest.mark.parametrize("n_samples,n_genes,design_factor,expected", [
    (2, 4, "batch", {"B0", "B1"}),
    (4, 2, None, {"0", "1", "2", "3"}),
    (3, 3, "absent", {"0", "1", "2"}),
])
def test_pca_mds_small_matrices_and_explicit_design(tmp_path, saved_figures, n_samples, n_genes, design_factor, expected):
    samples = [f"S{i}" for i in range(n_samples)]
    frame = pd.DataFrame(np.random.default_rng(1).normal(size=(n_genes, n_samples)), columns=samples)
    metadata = pd.DataFrame({"age": range(n_samples), "batch": [f"B{i % 2}" for i in range(n_samples)]}, index=samples)
    warnings = []
    paths = bulk._plot_pca_mds(frame, metadata.iloc[::-1], str(tmp_path), warnings, design_factor=design_factor)
    assert not warnings
    assert paths[2].shape == (n_samples, 2)
    assert len(paths[3]) == 2 and sum(paths[3]) <= 100.001
    for name in ("pca.svg", "mds.svg", "pca_mds.svg"):
        assert (tmp_path / name).is_file()
    assert {t.get_text() for t in saved_figures["pca.svg"].axes[0].get_legend().get_texts()} == expected


@pytest.mark.parametrize("n_samples,n_genes", [(1, 4), (2, 1), (5, 1), (0, 3), (3, 0)])
def test_pca_mds_skips_when_two_pcs_unavailable(tmp_path, n_samples, n_genes):
    samples = [f"S{i}" for i in range(n_samples)]
    frame = pd.DataFrame(np.ones((n_genes, n_samples)), columns=samples)
    warnings = []
    assert bulk._plot_pca_mds(frame, pd.DataFrame(index=samples), str(tmp_path), warnings) == (None, None, None, [])
    assert "skipped" in warnings[0] and "two" in warnings[0]


@pytest.mark.parametrize("design_factor,expected", [("batch", {"B0", "B1"}), (None, {"A0", "A1"}), ("absent", {"A0", "A1"})])
def test_sample_pca_grouping_and_percentages(tmp_path, saved_figures, design_factor, expected):
    metadata = pd.DataFrame({"age": ["A0", "A1"], "batch": ["B0", "B1"]}, index=["S0", "S1"])
    assert bulk._plot_sample_pca(np.array([[1, 2], [3, 4]]), metadata.index, metadata, [72.5, 20.0], str(tmp_path), design_factor)
    ax = saved_figures["sample_pca.svg"].axes[0]
    assert ax.get_xlabel() == "PC1 (72.5%)"
    assert ax.get_ylabel() == "PC2 (20.0%)"
    assert {t.get_text() for t in ax.get_legend().get_texts()} == expected


def test_sample_pca_skips_one_component(tmp_path):
    assert bulk._plot_sample_pca(np.ones((2, 1)), ["S0", "S1"], pd.DataFrame(index=["S0", "S1"]), [100], str(tmp_path)) is None


def test_gsea_signature_rejects_infinite_lfc_before_symbol_deduplication():
    frame = pd.DataFrame({"log2FoldChange": [np.inf, 2.5, -np.inf, np.nan, -1.5]}, index=["G1", "G1", "G2", "G3", "G4"])
    assert pathways._ranked_signature_frame(frame).to_dict("records") == [
        {"gene": "G1", "score": 2.5}, {"gene": "G4", "score": -1.5}]
    signature = pd.DataFrame({"gene": ["A", "B", "C", "D"], "score": [np.inf, -np.inf, 2, -1]})
    cleaned, raw, centered, ranks = pathways._prepare_gsea_signature(signature)
    assert list(cleaned["gene"]) == ["C", "D"]
    assert np.isfinite(raw).all() and np.isfinite(centered).all()
    assert ranks == {"C": 0, "D": 1}


def test_ora_adaptive_zero_cap_valid_odds_distinct_rows_and_svg(tmp_path, saved_figures, monkeypatch):
    import seaborn as sns

    seen = {}
    original = sns.scatterplot

    def capture(**kwargs):
        seen.update(kwargs)
        return original(**kwargs)

    monkeypatch.setattr(sns, "scatterplot", capture)
    prefix = "Shared pathway prefix " + "X" * 70
    rows = [
        {"term": prefix + " A", "padj": 0, "odds_ratio": 4, "overlap": "12/100"},
        {"term": prefix + " B", "padj": 1.234e-8, "odds_ratio": 2, "genes": ["G1", "G2"]},
        *[{"term": f"INVALID_{i}", "padj": 0.01, "odds_ratio": value} for i, value in enumerate([np.inf, -np.inf, np.nan, None, 0, -1])],
        *[{"term": f"INVALID_P_{i}", "padj": value, "odds_ratio": 2} for i, value in enumerate([np.inf, np.nan, -0.1, 1.1])],
    ]
    out = tmp_path / "ora.png"
    assert pathways.make_ora_dotplot(rows, "DB", "A vs B", str(out)) == str(out)
    assert out.stat().st_size > 1000 and out.with_suffix(".svg").stat().st_size > 1000
    data = seen["data"]
    assert len(data) == 2 and data["Term"].nunique() == 2
    true_neglog = -np.log10(1.234e-8)
    assert sorted(data["-log10(FDR)"]) == pytest.approx([true_neglog, true_neglog * 1.15])
    assert sorted(data["Genes"]) == [2, 12]
    assert np.isfinite(data["log2_OR"]).all()
    ax, cax = saved_figures["ora.png"].axes
    assert len(ax.get_yticklabels()) == 2
    assert "underflow" in cax.get_ylabel()
    assert seen["hue_norm"].vmax == pytest.approx(true_neglog * 1.15)


def test_ora_all_zero_fdr_has_sensible_ceiling(tmp_path, saved_figures):
    rows = [{"term": f"Term {i}", "padj": 0, "odds_ratio": i + 1} for i in range(2)]
    out = tmp_path / "zeros.png"
    assert pathways.make_ora_dotplot(rows, "DB", "A vs B", str(out))
    assert saved_figures["zeros.png"].axes[1].get_ylim()[1] < 20


def test_ora_all_invalid_odds_returns_none(tmp_path):
    assert pathways.make_ora_dotplot([{"term": "Term", "padj": 0.01, "odds_ratio": np.inf}], "DB", "A vs B", str(tmp_path / "empty.png")) is None


@pytest.mark.parametrize("size", [np.nan, np.inf, -np.inf, None, "invalid"])
def test_gsea_renderers_handle_unavailable_geneset_size(tmp_path, size):
    signature = pd.DataFrame({"gene": [f"G{i}" for i in range(12)], "score": np.linspace(2, -2, 12)})
    _, raw, centered, ranks = pathways._prepare_gsea_signature(signature)
    library = {"TERM_A": ["G0", "G1", "G2"]}
    result = pd.DataFrame({"nes": [2.0], "fdr": [0.01], "es": [0.5], "geneset_size": [size]}, index=["TERM_A"])
    assert pathways._geneset_size(result.iloc[0], 3) == 3
    pathways._render_gsea_running_sum(raw, centered, ranks, "TERM_A", library["TERM_A"], result.iloc[0], "A vs B", tmp_path / "running.png")
    pathways._render_gsea_top_table(raw, ranks, library, result, "A vs B", "DB", 1, tmp_path / "top.png", n=1)
    for name in ("running", "top"):
        for suffix in (".png", ".svg"):
            assert (tmp_path / (name + suffix)).stat().st_size > 1000


def test_de_export_accepts_none_symbol_map(tmp_path):
    frame = pd.DataFrame({"log2FoldChange": [2], "padj": [0.001]}, index=["GENE_A"])
    out = tmp_path / "de.tsv"
    assert pathways.export_de_table(frame, None, "A vs B", str(out)) == str(out)
    saved = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert saved.loc[0, "gene_symbol"] == ""
    assert saved.loc[0, "gene_id"] == "GENE_A"


@pytest.mark.parametrize("use_log", [False, True])
def test_paga_ignores_self_loops_and_legend_inverts_scaling(tmp_path, saved_figures, monkeypatch, use_log):
    conn = np.array([[100, 0.02], [0.02, 200]])
    seen = {}

    def layout(matrix):
        seen["conn"] = matrix.copy()
        return np.array([[0, 0], [1, 0]])

    monkeypatch.setattr(paga, "_spring_layout", layout)
    paga._draw_paga(conn, ["A", "B"], [10, 20], tmp_path / "paga.png", use_log=use_log)
    assert np.diag(seen["conn"]).tolist() == [0, 0]
    assert np.diag(conn).tolist() == [100, 200]  # input remains unchanged
    ax = saved_figures["paga.png"].axes[0]
    assert len(ax.lines) == 1 and ax.lines[0].get_linewidth() == pytest.approx(5)
    labels = [float(t.get_text().removeprefix("~")) for t in ax.get_legend().get_texts()]
    fractions = np.array([0.25, 0.5, 1.0])
    expected = (10 ** (fractions * np.log10(1 + 1000 * 0.02)) - 1) / 1000 if use_log else fractions * 0.02
    assert labels == pytest.approx(expected, abs=5.1e-5)


@pytest.mark.parametrize("use_log", [False, True])
def test_paga_without_offdiagonal_edges_has_no_width_legend(tmp_path, saved_figures, use_log):
    paga._draw_paga(np.eye(2), ["A", "B"], [10, 10], tmp_path / "empty_paga.png", use_log=use_log)
    ax = saved_figures["empty_paga.png"].axes[0]
    assert ax.get_legend() is None and not ax.lines


def test_paga_entrypoint_uses_offdiagonal_max_for_both_renderings(tmp_path, saved_figures, monkeypatch):
    import anndata as ad
    from scipy.sparse import csr_matrix
    from aria.utils import safe_h5ad

    adata = ad.AnnData(np.ones((4, 2)))
    adata.obs["leiden"] = pd.Categorical(["A", "A", "B", "B"])
    adata.uns["paga"] = {"groups": "leiden", "connectivities": csr_matrix([[100, 0.02], [0.02, 200]])}
    monkeypatch.setattr(safe_h5ad, "read_h5ad", lambda path: adata)
    result = paga.make_paga_figures({"h5ad_path": "unused.h5ad", "output_dir": str(tmp_path)})
    assert result["status"] == "success" and result["max_connectivity"] == 0.02
    assert set(result["figures"]) == {"paga_graph", "paga_log10_graph"}
    for name in ("paga_graph.png", "paga_graph_log10.png"):
        assert saved_figures[name].axes[0].lines[0].get_linewidth() == pytest.approx(5)
    np.testing.assert_array_equal(adata.uns["paga"]["connectivities"].diagonal(), [100, 200])


def test_continuous_umap_masks_nonfinite_and_keeps_finite_color_range(tmp_path, saved_figures):
    coords = np.arange(12).reshape(6, 2)
    values = np.array([0.1, 0.4, 0.8, np.inf, -np.inf, np.nan])
    umap._plot_one(coords, values, "dpt_pseudotime", tmp_path / "dpt.png", 2, "DPT", "UMAP")
    ax = saved_figures["dpt.png"].axes[0]
    missing, finite = ax.collections
    assert missing.get_zorder() < finite.get_zorder()
    np.testing.assert_array_equal(missing.get_offsets(), coords[3:])
    np.testing.assert_allclose(missing.get_facecolors()[0, :3], [0.82745098] * 3)
    np.testing.assert_array_equal(finite.get_offsets(), coords[:3])
    assert finite.get_clim() == pytest.approx((0.1, 0.8))
    assert "Non-finite" in ax.get_legend().get_texts()[0].get_text()


def test_continuous_umap_all_nonfinite_is_gray_without_colorbar(tmp_path, saved_figures):
    umap._plot_one(np.array([[0, 0], [1, 1]]), [np.nan, np.inf], "dpt_pseudotime", tmp_path / "missing.png", 2, "DPT", "UMAP")
    fig = saved_figures["missing.png"]
    assert len(fig.axes) == 1 and len(fig.axes[0].collections) == 1


@pytest.mark.parametrize("key,values", [
    ("codes", np.tile([0, 1], 20)),
    ("leiden", np.arange(40)),
    ("louvain", np.arange(40)),
    ("seurat_clusters", np.arange(40)),
    ("batch", np.arange(40, dtype=float)),
    ("sample_id", np.arange(40)),
    ("group", np.arange(40)),
    ("annotation", pd.Categorical(["A", "B"] * 20)),
])
def test_categorical_umap_uses_discrete_legend_and_deterministic_mixed_draw_order(tmp_path, saved_figures, key, values):
    coords = np.column_stack([np.arange(40), np.zeros(40)])
    for name in ("first.png", "second.png"):
        umap._plot_one(coords, values, key, tmp_path / name, 2, "Groups", "UMAP")
    first = saved_figures["first.png"]
    assert len(first.axes) == 1
    ax = first.axes[0]
    assert len(ax.collections) == 1
    expected = coords[np.random.default_rng(0).permutation(40)]
    np.testing.assert_array_equal(ax.collections[0].get_offsets(), expected)
    np.testing.assert_array_equal(saved_figures["second.png"].axes[0].collections[0].get_offsets(), expected)
    assert len(ax.get_legend().get_texts()) == pd.Series(values).nunique()
    palette = dict(zip(sorted(pd.Series(values).astype(str).unique()), umap._categorical_colors(pd.Series(values).nunique())))
    expected_colors = [palette[str(values[i])] for i in np.random.default_rng(0).permutation(40)]
    np.testing.assert_allclose(ax.collections[0].get_facecolors()[:, :3], expected_colors)


def test_high_cardinality_integer_measurement_remains_continuous(tmp_path, saved_figures):
    umap._plot_one(np.column_stack([np.arange(40)] * 2), np.arange(40), "n_counts", tmp_path / "counts.png", 2, "Counts", "UMAP")
    assert len(saved_figures["counts.png"].axes) == 2


@pytest.mark.parametrize("shape,expected", [((4, 100), 3), ((100, 4), 3), ((100, 100), 50), ((2, 5), 1)])
def test_fallback_umap_pca_is_clamped_to_both_dimensions(monkeypatch, shape, expected):
    import anndata as ad

    seen = {}
    def pca_call(adata, **kwargs):
        seen.update(kwargs)

    fake_scanpy = types.SimpleNamespace(
        pp=types.SimpleNamespace(normalize_total=lambda *a, **k: None,
                                 log1p=lambda *a, **k: None, pca=pca_call,
                                 neighbors=lambda *a, **k: None),
        tl=types.SimpleNamespace(umap=lambda *a, **k: None))
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    assert umap._compute_umap(ad.AnnData(np.ones(shape))) == "X_umap"
    assert seen["n_comps"] == expected


def test_fallback_umap_runs_real_pca_on_small_sample_matrix(tmp_path, monkeypatch):
    # The managed sandbox cannot write Numba's default cache under site-packages.
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))
    import anndata as ad
    import numba
    monkeypatch.setattr(numba.config, "CACHE_DIR", str(tmp_path / "numba-cache"))
    import scanpy as sc

    adata = ad.AnnData(np.random.default_rng(5).poisson(5, size=(4, 8)).astype(float))
    monkeypatch.setattr(sc.pp, "neighbors", lambda *a, **kw: None)
    monkeypatch.setattr(sc.tl, "umap", lambda *a, **kw: None)
    assert umap._compute_umap(adata) == "X_umap"
    assert adata.obsm["X_pca"].shape == (4, 3)
    assert np.isfinite(adata.obsm["X_pca"]).all()


def _interaction(score, rank, metric="specificity_rank"):
    return {"source": "GroupA", "target": "GroupB", "ligand": f"L{rank}", "receptor": "R",
            "rank": rank, "score": score, "rank_metric": metric}


def test_cellcomm_uses_real_pvalue_scores_and_orders_before_top_n(tmp_path, saved_figures):
    rows = [_interaction(1e-2, 1), _interaction(0, 2), _interaction(1e-6, 3),
            _interaction(np.inf, 4), _interaction(np.nan, 5), _interaction(-1, 6)]
    findings = {"cell_communication": {"method": "liana_rank_aggregate (specificity_rank)", "top_interactions": rows}}
    assert narrative.render_cellcomm_top_pairs_bar(findings, tmp_path / "pairs.png", top_n=3)
    ax = saved_figures["pairs.png"].axes[0]
    assert [p.get_width() for p in ax.patches] == pytest.approx([6.9, 6, 2])
    assert "L2-R" in ax.get_yticklabels()[0].get_text()
    assert "-log10(specificity_rank)" in ax.get_xlabel() and "underflow" in ax.get_xlabel()


@pytest.mark.parametrize("metric", ["lrscore", "magnitude", "specificity"])
def test_cellcomm_direct_scores_are_not_inverted(tmp_path, saved_figures, metric):
    findings = {"cell_communication": {"method": f"liana_rank_aggregate ({metric})", "top_interactions": [
        _interaction(0.2, 1, metric), _interaction(0.9, 2, metric), _interaction(0.5, 3, metric), _interaction(np.inf, 4, metric)]}}
    narrative.render_cellcomm_top_pairs_bar(findings, tmp_path / "direct.png")
    ax = saved_figures["direct.png"].axes[0]
    assert [p.get_width() for p in ax.patches] == pytest.approx([0.9, 0.5, 0.2])
    assert ax.get_xlabel() == metric and "L2-R" in ax.get_yticklabels()[0].get_text()


def test_cellcomm_all_zero_pvalues_get_annotated_finite_cap(tmp_path, saved_figures):
    findings = {"cell_communication": {"method": "liana", "top_interactions": [_interaction(0, 1), _interaction(0, 2)]}}
    narrative.render_cellcomm_top_pairs_bar(findings, tmp_path / "zero_pairs.png")
    ax = saved_figures["zero_pairs.png"].axes[0]
    assert [p.get_width() for p in ax.patches] == [10, 10]
    assert "underflow" in ax.get_xlabel()


def test_cellcomm_no_usable_scores_does_not_invent_rank_bars(tmp_path):
    findings = {"cell_communication": {"method": "liana", "top_interactions": [
        {"rank": 1}, _interaction(np.nan, 2), _interaction(np.inf, 3)]}}
    assert narrative.render_cellcomm_top_pairs_bar(findings, tmp_path / "invalid.png") is None


def test_cellcomm_can_plot_named_pvalues_without_ordinal_rank(tmp_path, saved_figures):
    findings = {"cell_communication": {"method": "liana (pval)", "top_interactions": [
        {"pval": 0, "source": "A", "target": "B", "ligand": "L0"},
        {"pval": 1e-4, "source": "A", "target": "B", "ligand": "L1"}]}}
    narrative.render_cellcomm_top_pairs_bar(findings, tmp_path / "pvalues.png")
    ax = saved_figures["pvalues.png"].axes[0]
    assert [patch.get_width() for patch in ax.patches] == pytest.approx([4.6, 4])
    assert "-log10(pval)" in ax.get_xlabel()


@pytest.mark.parametrize("maximum", [None, np.nan, np.inf, -np.inf])
def test_paga_html_handles_unavailable_max_connectivity_without_mutating_input(maximum):
    findings = {"trajectory": {"status": "success", "paga": {
        "max_connectivity": maximum, "top_connections": {"A-B": 0.01}, "n_strong": 0}}}
    section = narrative.build_scrna_html_section(findings)
    assert "Max connectivity = unavailable" in section
    assert "A-B" in section
    value = findings["trajectory"]["paga"]["max_connectivity"]
    assert value is maximum
