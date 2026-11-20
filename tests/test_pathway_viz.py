from __future__ import annotations


def test_ranked_signature_keeps_symbol_ids_without_symbol_map():
    pd = __import__("pandas")
    from aria.scripts.rna_pathway_viz import _ranked_signature_frame

    de = pd.DataFrame(
        {"log2FoldChange": [1.2, -2.5, 0.4]},
        index=["GENE_A", "GENE_B", "GENE_A"],
    )

    signature = _ranked_signature_frame(de, symbol_map={})

    assert signature.to_dict("records") == [
        {"gene": "GENE_B", "score": -2.5},
        {"gene": "GENE_A", "score": 1.2},
    ]


def test_ranked_signature_drops_unmapped_ensembl_ids():
    pd = __import__("pandas")
    from aria.scripts.rna_pathway_viz import _ranked_signature_frame

    de = pd.DataFrame(
        {"log2FoldChange": [1.2, -2.5]},
        index=["ENSG000001.1", "GENE_B"],
    )

    signature = _ranked_signature_frame(de, symbol_map={})

    assert signature.to_dict("records") == [{"gene": "GENE_B", "score": -2.5}]


def test_partition_gsea_results_excludes_nonfinite_and_sorts_deterministically():
    import numpy as np
    import pandas as pd
    from aria.scripts.rna_pathway_viz import _partition_gsea_results

    raw = pd.DataFrame(
        {
            "es": [0.17, 0.51, 0.52, -0.45, 0.30],
            "nes": [np.inf, 4.82, 4.94, -np.inf, np.nan],
            "pval": [0.0, 1.4e-6, 8.0e-7, 0.0, 0.01],
            "fdr": [0.0, 0.002, 0.002, 0.0, 0.05],
        },
        index=[
            "INF_POS_TERM",
            "RIBOSOME_B",
            "RIBOSOME_A",
            "INF_NEG_TERM",
            "NAN_TERM",
        ],
    )
    stable, unstable = _partition_gsea_results(raw)

    assert list(stable.index) == ["RIBOSOME_A", "RIBOSOME_B"]
    assert set(unstable.index) == {"INF_POS_TERM", "INF_NEG_TERM", "NAN_TERM"}


def test_make_gsea_running_sums_excludes_nonfinite_nes_from_plots_and_counts(
    tmp_path, monkeypatch
):
    import sys
    import types
    import numpy as np
    import pandas as pd
    import aria.scripts.rna_pathway_viz as viz

    calls = {"gsea_kwargs": None, "sig_genes": None, "running_terms": [], "top_table": None}

    fake_blitz = types.ModuleType("blitzgsea")
    fake_blitz.enrichr = types.SimpleNamespace(
        get_library=lambda lib: {
            "INF_TERM": ["G1"],
            "STABLE_POS": ["G1", "G2"],
            "STABLE_NEG": ["G3", "G4"],
        }
    )

    def _fake_gsea(signature, library, **kwargs):
        calls["gsea_kwargs"] = kwargs
        calls["sig_genes"] = list(signature["gene"])
        return pd.DataFrame(
            {
                "es": [0.17, 0.52, -0.48],
                "nes": [np.inf, 4.94, -4.50],
                "pval": [0.0, 8.0e-7, 5.0e-6],
                "sidak": [0.0, 0.004, 0.02],
                "fdr": [0.0, 0.002, 0.01],
                "geneset_size": [1542, 50, 48],
                "leading_edge": ["G1", "G1,G2", "G3,G4"],
            },
            index=pd.Index(["INF_TERM", "STABLE_POS", "STABLE_NEG"], name="Term"),
        )

    fake_blitz.gsea = _fake_gsea
    monkeypatch.setitem(sys.modules, "blitzgsea", fake_blitz)

    # Intercept renderers in light env (where matplotlib may be absent) while
    # validating exact arguments passed to them.
    def _spy_rs(**kwargs):
        calls["running_terms"].append(kwargs["term"])
        assert np.isfinite(float(kwargs["row"]["nes"]))
        kwargs["output_png"].write_bytes(b"PNG")
        return str(kwargs["output_png"])

    def _spy_tt(**kwargs):
        res = kwargs["result_stable"]
        calls["top_table"] = {
            "n": kwargs["n"],
            "terms": list(res.index),
            "all_finite": bool(np.isfinite(res["nes"]).all()),
            "gene_keys": sorted(kwargs["gene_to_rank"].keys()),
        }
        kwargs["output_png"].write_bytes(b"PNG")
        return str(kwargs["output_png"])

    monkeypatch.setattr(viz, "_render_gsea_running_sum", _spy_rs)
    monkeypatch.setattr(viz, "_render_gsea_top_table", _spy_tt)

    # Mixed-case mouse-like symbols must be upper-cased for library matching
    de = pd.DataFrame(
        {"log2FoldChange": [2.5, 1.4, -1.2, -2.8], "padj": [0.001, 0.01, 0.02, 0.001]},
        index=["g1", "G2", "g3", "G4"],
    )
    out = viz.make_gsea_running_sums(
        de_results_df=de,
        symbol_map={},
        contrast_name="REV_ERBa_KO vs WT",
        output_dir=str(tmp_path),
        top_n_plots=3,
    )

    assert calls["gsea_kwargs"] == {"min_size": 15, "max_size": 500, "seed": 0}
    assert calls["sig_genes"] == ["G1", "G2", "G3", "G4"]
    assert out["n_numeric_unstable"] == 1
    assert out["n_pathways"] == 2
    assert calls["running_terms"] == ["STABLE_POS", "STABLE_NEG"]
    assert "INF_TERM" not in calls["running_terms"]
    assert calls["top_table"]["all_finite"] is True
    assert calls["top_table"]["terms"] == ["STABLE_POS", "STABLE_NEG"]
    assert calls["top_table"]["n"] == 2
    assert calls["top_table"]["gene_keys"] == ["G1", "G2", "G3", "G4"]

    # CSV puts finite pathways first and keeps unstable row at tail for narrator disclosure
    saved = pd.read_csv(out["gsea_table"])
    assert list(saved["Term"]) == ["STABLE_POS", "STABLE_NEG", "INF_TERM"]


def test_load_gsea_library_respects_air_gap_and_uses_local_fallback(
    tmp_path, monkeypatch
):
    import types
    from aria.utils import ora
    from aria.scripts.rna_pathway_viz import _load_gsea_library

    # Stage only GO_Biological_Process_2021 locally; enable air-gap.
    lib_dir = tmp_path / "GO_Biological_Process_2021"
    lib_dir.mkdir(parents=True)
    (lib_dir / "GO_Biological_Process_2021.gmt").write_text(
        "My Pathway\t\tTrp53\tCdkn1a\n", encoding="utf-8"
    )
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")

    called = {"net": False}

    def _forbidden_get_library(name):
        called["net"] = True
        raise AssertionError("Network library fetch must not run under ARIA_AIR_GAPPED=1")

    fake_blitz = types.SimpleNamespace(
        enrichr=types.SimpleNamespace(get_library=_forbidden_get_library)
    )
    lib, label = _load_gsea_library("GO_Biological_Process_2023", fake_blitz)

    assert called["net"] is False
    assert label == "GO_Biological_Process_2021"
    assert lib == {"My Pathway": ["TRP53", "CDKN1A"]}


def test_native_gsea_renderers_produce_png_and_svg_for_pos_and_neg_nes(tmp_path):
    pytest = __import__("pytest")
    pytest.importorskip("matplotlib")
    import pandas as pd
    from aria.scripts.rna_pathway_viz import (
        _prepare_gsea_signature,
        _render_gsea_running_sum,
        _render_gsea_top_table,
    )

    genes = [f"G{i}" for i in range(60)]
    scores = [3.0 - 0.1 * i for i in range(60)]
    sig_df = pd.DataFrame({"gene": genes, "score": scores})
    _, scores_raw, scores_centered, gene_to_rank = _prepare_gsea_signature(sig_df)
    library = {
        "Up Pathway (GO:0000001)": genes[:12],
        "Down Pathway (GO:0000002)": genes[-12:],
    }
    res_df = pd.DataFrame(
        {
            "es": [0.72, -0.68],
            "nes": [4.25, -3.95],
            "pval": [1e-6, 3e-6],
            "fdr": [0.0012, 0.0024],
            "geneset_size": [12, 12],
        },
        index=pd.Index(
            ["Up Pathway (GO:0000001)", "Down Pathway (GO:0000002)"], name="Term"
        ),
    )

    rs_pos = tmp_path / "rs_pos.png"
    rs_neg = tmp_path / "rs_neg.png"
    tt_png = tmp_path / "top_table.png"

    _render_gsea_running_sum(
        scores_raw, scores_centered, gene_to_rank,
        "Up Pathway (GO:0000001)", library["Up Pathway (GO:0000001)"],
        res_df.iloc[0], "BMAL1_KO vs REV_ERBa_KO", rs_pos,
    )
    _render_gsea_running_sum(
        scores_raw, scores_centered, gene_to_rank,
        "Down Pathway (GO:0000002)", library["Down Pathway (GO:0000002)"],
        res_df.iloc[1], "BMAL1_KO vs REV_ERBa_KO", rs_neg,
    )
    _render_gsea_top_table(
        scores_raw, gene_to_rank, library, res_df,
        "BMAL1_KO vs REV_ERBa_KO", "GO_Biological_Process_2023", 2, tt_png, n=2,
    )

    for p in (rs_pos, rs_neg, tt_png):
        assert p.is_file() and p.stat().st_size > 2000
        assert p.with_suffix(".svg").is_file() and p.with_suffix(".svg").stat().st_size > 1000
