"""Guards for the A2 Kang/muscat external reference lane (reference_kang.py).

The scorer logic runs without pydeseq2 by patching ARIA's DE core; the real lane
is gated on a staged muscat export (ARIA_KANG_MUSCAT_EXPORT).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest


def _write_export(export: Path, effects: dict[str, float], clusters: list[str]):
    """Stage a synthetic muscat export: per-cluster pseudobulk + muscat tables
    whose logFC follows a planted per-gene effect, plus the sample table."""
    import numpy as np
    import pandas as pd

    export.mkdir(parents=True, exist_ok=True)
    a_cols = [f"stim{i}" for i in range(3)]
    b_cols = [f"ctrl{i}" for i in range(3)]
    pd.DataFrame({"sample": a_cols + b_cols,
                  "group": ["stim"] * 3 + ["ctrl"] * 3}).to_csv(
        export / "samples.tsv", sep="\t", index=False)
    rng = np.random.default_rng(1)
    genes = list(effects)
    for cl in clusters:
        counts = pd.DataFrame(
            rng.poisson(200, size=(len(genes), 6)), index=genes, columns=a_cols + b_cols)
        counts.insert(0, "gene", counts.index)
        counts.to_csv(export / f"pb_{cl}.tsv", sep="\t", index=False)
        mus = pd.DataFrame({
            "gene": genes,
            "logFC": [effects[g] + float(rng.normal(0, 0.1)) for g in genes],
            "p_val": [0.0001 if effects[g] != 0 else 0.6 for g in genes],
            "p_adj": [0.001 if effects[g] != 0 else 0.8 for g in genes],
        })
        mus.to_csv(export / f"muscat_{cl}.tsv", sep="\t", index=False)
    (export / "clusters.json").write_text(json.dumps({
        "clusters": clusters,
        "contrast": {"numerator": "stim", "denominator": "ctrl"},
        "dataset": "synthetic-test",
    }), encoding="utf-8")


def test_score_aria_vs_muscat_concordant(monkeypatch, tmp_path):
    import numpy as np
    import pandas as pd
    from aria.benchmarks import reference_kang
    import aria.scripts.rna_bulk_de as rbd

    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(60)]
    effects = {g: (float(rng.choice([-3, -2, 2, 3])) if i < 20 else 0.0)
               for i, g in enumerate(genes)}
    _write_export(tmp_path, effects, ["CD4_T", "B_cells"])

    def fake_run_deseq2(counts, meta, factor, num, den, **kw):
        idx = [str(g) for g in counts.index]
        lfc = [effects.get(g, 0.0) + float(rng.normal(0, 0.1)) for g in idx]
        padj = [0.001 if effects.get(g, 0.0) != 0 else 0.7 for g in idx]
        df = pd.DataFrame({"log2FoldChange": lfc, "pvalue": padj, "padj": padj},
                          index=idx)
        sig = [g for g in idx if effects.get(g, 0.0) != 0]
        return {"status": "success", "results": df, "sig_genes": sig}, []

    monkeypatch.setattr(rbd, "_run_deseq2", fake_run_deseq2)

    m = reference_kang.score_aria_vs_muscat(tmp_path)
    assert m["status"] == "pass", m
    s = m["summary"]
    assert s["n_clusters_scored"] == 2
    assert s["mean_lfc_pearson"] >= 0.9           # ARIA vs muscat logFC agree
    assert s["mean_lfc_spearman_sig"] >= 0.8      # rank agreement on DE genes
    assert s["mean_sig_jaccard"] >= 0.8           # significant sets overlap
    assert s["mean_shared_sig_direction_agreement"] >= 0.99
    # Per-cluster shape.
    cl = m["per_cluster"]["CD4_T"]
    assert cl["status"] == "success" and cl["n_muscat_sig"] == 20


def test_score_aria_vs_muscat_skips_without_export(tmp_path):
    from aria.benchmarks.reference_kang import score_aria_vs_muscat

    out = score_aria_vs_muscat(tmp_path / "absent")
    assert out["status"] == "skipped"
    assert out["benchmark"] == "A2_kang_muscat"


@pytest.mark.skipif(
    not os.environ.get("ARIA_KANG_MUSCAT_EXPORT"),
    reason="set ARIA_KANG_MUSCAT_EXPORT to a staged muscat/Kang export dir",
)
def test_a2_kang_muscat_real_export():
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.reference_kang import run_kang_muscat_benchmark

    m = run_kang_muscat_benchmark(os.environ["ARIA_KANG_MUSCAT_EXPORT"])
    assert m["status"] == "pass", m
    assert m["summary"]["n_clusters_scored"] >= 1


def test_score_aria_vs_muscat_real_deseq2_unmocked(tmp_path):
    """D3-006: run score_aria_vs_muscat against real unmocked pyDESeq2 on
    planted NB pseudobulk counts matching the muscat reference effects."""
    import numpy as np
    import pandas as pd

    pytest.importorskip("pydeseq2")
    from aria.benchmarks import reference_kang

    rng = np.random.default_rng(17)
    n_genes = 120
    genes = [f"G{i:03d}" for i in range(n_genes)]
    effects = {}
    for i, g in enumerate(genes):
        if i < 15:
            effects[g] = 2.5
        elif i < 30:
            effects[g] = -2.5
        else:
            effects[g] = 0.0

    a_cols = [f"stim{i}" for i in range(4)]
    b_cols = [f"ctrl{i}" for i in range(4)]
    pd.DataFrame({
        "sample": a_cols + b_cols,
        "group": ["stim"] * 4 + ["ctrl"] * 4,
    }).to_csv(tmp_path / "samples.tsv", sep="\t", index=False)

    clusters = ["CD4_T", "B_cells"]
    for cl in clusters:
        mat = np.zeros((n_genes, 8), dtype=int)
        for j, col in enumerate(a_cols + b_cols):
            mu = np.full(n_genes, 160.0, dtype=float)
            if col.startswith("stim"):
                for i, g in enumerate(genes):
                    mu[i] *= 2.0 ** effects[g]
            mat[:, j] = rng.negative_binomial(n=35, p=35.0 / (35.0 + mu))
        counts = pd.DataFrame(mat, index=genes, columns=a_cols + b_cols)
        counts.insert(0, "gene", counts.index)
        counts.to_csv(tmp_path / f"pb_{cl}.tsv", sep="\t", index=False)
        mus = pd.DataFrame({
            "gene": genes,
            "logFC": [effects[g] for g in genes],
            "p_val": [1e-6 if effects[g] != 0 else 0.7 for g in genes],
            "p_adj": [1e-4 if effects[g] != 0 else 0.8 for g in genes],
        })
        mus.to_csv(tmp_path / f"muscat_{cl}.tsv", sep="\t", index=False)

    (tmp_path / "clusters.json").write_text(json.dumps({
        "clusters": clusters,
        "contrast": {"numerator": "stim", "denominator": "ctrl"},
        "dataset": "synthetic-unmocked-deseq2",
    }), encoding="utf-8")

    m = reference_kang.score_aria_vs_muscat(tmp_path)
    assert m["status"] == "pass", m
    s = m["summary"]
    assert s["n_clusters_scored"] == 2
    assert s["mean_lfc_pearson"] >= 0.90
    assert s["mean_lfc_spearman_sig"] >= 0.85
    assert s["mean_sig_jaccard"] >= 0.80
    assert s["mean_shared_sig_direction_agreement"] >= 0.99
