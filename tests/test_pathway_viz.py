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
    from aria.scripts.rna_pathway_viz import make_gsea_running_sums

    calls = {"gsea_kwargs": None, "running_terms": [], "top_table": None}

    class _FakeFig:
        def savefig(self, path, **kwargs):
            from pathlib import Path
            Path(path).write_bytes(b"PNG")

    fake_blitz = types.ModuleType("blitzgsea")
    fake_blitz.enrichr = types.SimpleNamespace(
        get_library=lambda lib: {"INF_TERM": ["G1"], "STABLE_A": ["G1"], "STABLE_B": ["G2"]}
    )

    def _fake_gsea(signature, library, **kwargs):
        calls["gsea_kwargs"] = kwargs
        return pd.DataFrame(
            {
                "es": [0.17, 0.52, 0.48],
                "nes": [np.inf, 4.94, 4.50],
                "pval": [0.0, 8.0e-7, 5.0e-6],
                "sidak": [0.0, 0.004, 0.02],
                "fdr": [0.0, 0.002, 0.01],
                "geneset_size": [1542, 50, 48],
                "leading_edge": ["G1", "G1", "G2"],
            },
            index=pd.Index(["INF_TERM", "STABLE_A", "STABLE_B"], name="Term"),
        )

    def _fake_running_sum(signature, term, library, result=None, compact=False):
        calls["running_terms"].append(term)
        assert np.isfinite(result["nes"]).all()
        return _FakeFig()

    def _fake_top_table(signature, library, result, n=15):
        calls["top_table"] = {
            "n": n,
            "terms": list(result.index),
            "sig_index": list(signature.index),
            "all_finite": bool(np.isfinite(result["nes"]).all()),
        }
        return _FakeFig()

    fake_blitz.gsea = _fake_gsea
    fake_blitz.plot = types.SimpleNamespace(
        running_sum=_fake_running_sum,
        top_table=_fake_top_table,
    )
    monkeypatch.setitem(sys.modules, "blitzgsea", fake_blitz)

    # Also stub matplotlib.pyplot.close if matplotlib is absent in light env
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        fake_mpl = types.ModuleType("matplotlib")
        fake_mpl.use = lambda *a, **k: None
        fake_plt = types.ModuleType("matplotlib.pyplot")
        fake_plt.close = lambda *a, **k: None
        monkeypatch.setitem(sys.modules, "matplotlib", fake_mpl)
        monkeypatch.setitem(sys.modules, "matplotlib.pyplot", fake_plt)

    de = pd.DataFrame(
        {"log2FoldChange": [2.5, -1.8], "padj": [0.001, 0.01]},
        index=["G1", "G2"],
    )
    out = make_gsea_running_sums(
        de_results_df=de,
        symbol_map={},
        contrast_name="REV_ERBa_KO vs WT",
        output_dir=str(tmp_path),
        top_n_plots=3,
    )

    assert calls["gsea_kwargs"] == {"min_size": 15, "max_size": 500, "seed": 0}
    assert out["n_numeric_unstable"] == 1
    assert out["n_pathways"] == 2
    assert calls["running_terms"] == ["STABLE_A", "STABLE_B"]
    assert "INF_TERM" not in calls["running_terms"]
    assert calls["top_table"]["all_finite"] is True
    assert calls["top_table"]["terms"] == ["STABLE_A", "STABLE_B"]
    assert calls["top_table"]["n"] == 2
    assert calls["top_table"]["sig_index"] == ["G1", "G2"]

    # CSV puts finite pathways first and keeps unstable row at tail for narrator disclosure
    saved = pd.read_csv(out["gsea_table"])
    assert list(saved["Term"]) == ["STABLE_A", "STABLE_B", "INF_TERM"]
