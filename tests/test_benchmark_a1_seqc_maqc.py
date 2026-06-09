"""Guards for the A1 SEQC/MAQC reference-data lane (aria/benchmarks/reference_seqc.py).

Light unit tests (numpy/pandas) run everywhere; the scorer end-to-end test is
pydeseq2-gated (heavy lane); the real-bundle test runs only when
``ARIA_SEQC_MAQC_BUNDLE`` points at a staged reference bundle.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


def test_auc_separates_classes_and_handles_ties():
    from aria.benchmarks.reference_seqc import _auc

    assert _auc([3.0, 4.0, 5.0], [0.0, 1.0, 2.0]) == 1.0   # perfect
    assert _auc([0.0, 1.0, 2.0], [3.0, 4.0, 5.0]) == 0.0   # reversed
    assert _auc([1.0, 1.0], [1.0, 1.0]) == 0.5             # all ties
    assert _auc([1.0], []) is None                          # empty class


def test_load_bundle_returns_none_when_absent(tmp_path):
    from aria.benchmarks.reference_seqc import load_seqc_maqc_bundle

    assert load_seqc_maqc_bundle(tmp_path / "missing") is None


def _write_bundle(bundle_dir: Path, counts, samples, taqman):
    import pandas as pd

    bundle_dir.mkdir(parents=True, exist_ok=True)
    c = counts.copy()
    c.insert(0, "gene", c.index)
    c.to_csv(bundle_dir / "counts.tsv", sep="\t", index=False)
    pd.DataFrame(samples).to_csv(bundle_dir / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(taqman).to_csv(bundle_dir / "taqman.tsv", sep="\t", index=False)


def test_load_bundle_parses_staged_files(tmp_path):
    import pandas as pd
    from aria.benchmarks.reference_seqc import load_seqc_maqc_bundle

    counts = pd.DataFrame(
        {"A_1": [10, 20], "B_1": [5, 40]}, index=["GENE_1", "GENE_2"]
    )
    _write_bundle(
        tmp_path,
        counts,
        {"sample": ["A_1", "B_1"], "group": ["A", "B"]},
        {"gene": ["GENE_1", "GENE_2"], "log2_ab": [1.0, -3.0]},
    )
    bundle = load_seqc_maqc_bundle(tmp_path)
    assert bundle is not None
    assert list(bundle["counts"].index) == ["GENE_1", "GENE_2"]
    assert bundle["taqman_log2_ab"]["GENE_2"] == -3.0
    assert set(bundle["samples"]["group"]) == {"A", "B"}


def test_runner_skips_honestly_without_bundle(tmp_path):
    from aria.benchmarks.reference_seqc import run_seqc_maqc_a1_benchmark

    out = run_seqc_maqc_a1_benchmark(tmp_path / "absent")
    assert out["status"] == "skipped"
    assert "fabricat" not in str(out).lower() or "Nothing is" in str(out)
    assert out["expected_bundle"]["files"] == ["counts.tsv", "samples.tsv", "taqman.tsv"]


def test_score_seqc_maqc_recovers_reference_truth(tmp_path):
    """End-to-end: build a bundle from the simulated bulk dataset (treating its
    known true log2FC as the TaqMan truth) and confirm ARIA's real DE path scores
    high LFC concordance and AUC against it. Maps COND_B->A, COND_A->B so ARIA's
    log2(A/B) aligns with true_log2fc = log2(COND_B/COND_A)."""
    pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.reference_seqc import (
        load_seqc_maqc_bundle,
        score_seqc_maqc_a1,
    )
    from aria.benchmarks.synthetic_de import simulate_bulk_dataset

    ds = simulate_bulk_dataset(n_genes=400, n_de=60, replicates_per_condition=5, seed=7)
    cond = ds.metadata["condition"]
    group_of = {"COND_B": "A", "COND_A": "B"}
    samples = {
        "sample": list(ds.counts.columns),
        "group": [group_of[cond[s]] for s in ds.counts.columns],
    }
    taqman = {
        "gene": list(ds.true_log2fc),
        "log2_ab": [ds.true_log2fc[g] for g in ds.true_log2fc],
    }
    _write_bundle(tmp_path, ds.counts, samples, taqman)

    bundle = load_seqc_maqc_bundle(tmp_path)
    manifest = score_seqc_maqc_a1(bundle, numerator="A", denominator="B")

    conc = manifest["axes"]["lfc_concordance"]
    det = manifest["axes"]["taqman_de_detection"]
    assert manifest["status"] == "pass", manifest
    # SEQC-standard log-ratio concordance (Pearson) is the gated metric.
    assert conc["pearson"] >= 0.7
    # Rank concordance on real-signal genes is robust to the null-tie mass.
    assert conc["spearman_signal"] >= 0.5
    # TaqMan-DE detection AUC.
    assert det["auc"] >= 0.8
    # No titration mixtures in this bundle -> honestly not computed.
    assert manifest["axes"]["titration_monotonicity"]["status"] == "not_computed"


def test_multisite_cross_concordance_matrix(monkeypatch, tmp_path):
    """The cross-site lane must build a symmetric pairwise log2FC concordance
    matrix and summarize the off-diagonal. Patch the per-site DE so the test
    needs no pydeseq2: two near-identical sites + one decorrelated site."""
    import numpy as np
    import pandas as pd
    from aria.benchmarks import reference_seqc

    rng = np.random.default_rng(0)
    base = pd.Series(rng.normal(0, 2, size=300), index=[f"G{i}" for i in range(300)])
    site_lfc = {
        "S1": base,
        "S2": base + rng.normal(0, 0.1, size=300),          # ~identical to S1
        "S3": pd.Series(rng.normal(0, 2, size=300), index=base.index),  # unrelated
    }

    def fake_site_de_lfc(bundle, num, den, padj, lfc):
        site = bundle["__site__"]
        de = {"results": pd.DataFrame({"log2FoldChange": site_lfc[site],
                                       "pvalue": 0.01}, index=base.index)}
        return de, site_lfc[site], 5, 5

    def fake_summary(*a, **k):
        return {"pearson": 0.95, "spearman_signal": 0.95, "auc": 0.9, "n_overlap": 300}

    def fake_load(d):
        name = Path(d).name
        if name not in site_lfc:        # simulate an unstaged site
            return None
        return {"__site__": name, "taqman_log2_ab": {}}

    monkeypatch.setattr(reference_seqc, "_site_de_lfc", fake_site_de_lfc)
    monkeypatch.setattr(reference_seqc, "_taqman_summary", fake_summary)
    monkeypatch.setattr(reference_seqc, "load_seqc_maqc_bundle", fake_load)

    m = reference_seqc.run_seqc_maqc_multisite(
        {"S1": "S1", "S2": "S2", "S3": "S3", "S4_missing": "S4_missing"},
    )
    # S4 has no bundle -> honest skip; the other three score.
    assert m["per_site"]["S4_missing"]["status"] == "skipped"
    cs = m["cross_site"]
    assert cs["n_sites"] == 3 and cs["n_pairs"] == 3
    # Matrix symmetric, diagonal == 1.
    pm = cs["pearson_matrix"]
    assert pm["S1"]["S1"] == 1.0
    assert pm["S1"]["S2"] == pm["S2"]["S1"]
    # S1~S2 highly concordant; S3 decorrelated drags the min down.
    assert pm["S1"]["S2"] >= 0.99
    assert cs["min_offdiagonal_pearson"] < 0.5
    assert cs["mean_offdiagonal_pearson"] is not None


def test_ercc_dose_response_recovers_known_design(tmp_path):
    """ERCC dose-response uses only CPM (no pydeseq2): build an ERCC bundle whose
    A/B counts follow the known Mix1/Mix2 ratios across a concentration ladder,
    and confirm the scorer recovers the per-subgroup fold-changes and dynamic
    range."""
    import numpy as np
    import pandas as pd
    from aria.benchmarks.reference_seqc import (
        load_seqc_maqc_bundle,
        score_ercc_dose_response,
    )

    # 4 ERCC subgroups x known log2(Mix1/Mix2) and a concentration ladder.
    sub_log2 = {"A": 2.0, "B": 0.0, "C": -0.58, "D": -1.0}
    rng = np.random.default_rng(3)
    a_cols = [f"A_{i}" for i in range(1, 4)]
    b_cols = [f"B_{i}" for i in range(1, 4)]
    erows, truth_rows = {}, []
    for k in range(40):
        g = list(sub_log2)[k % 4]
        conc1 = float(10 ** rng.uniform(1, 5))            # 4-order ladder
        conc2 = conc1 / (2 ** sub_log2[g])
        scale = 30.0
        a = rng.poisson(conc1 * scale, size=3) + 5
        b = rng.poisson(conc2 * scale, size=3) + 5
        eid = f"ERCC-{k:05d}"
        erows[eid] = list(a) + list(b)
        truth_rows.append((eid, g, conc1, conc2, 2 ** sub_log2[g], sub_log2[g]))

    ercc = pd.DataFrame.from_dict(erows, orient="index", columns=a_cols + b_cols)
    ercc.insert(0, "ercc_id", ercc.index)
    # A large flat gene matrix so the library sizes are balanced across A/B.
    genes = pd.DataFrame(
        rng.poisson(100, size=(500, 6)), columns=a_cols + b_cols,
        index=[f"G{i}" for i in range(500)],
    )
    genes.insert(0, "gene", genes.index)

    bundle_dir = tmp_path
    genes.to_csv(bundle_dir / "counts.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": a_cols + b_cols,
                  "group": ["A"] * 3 + ["B"] * 3}).to_csv(
        bundle_dir / "samples.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": ["G0"], "log2_ab": [0.0]}).to_csv(
        bundle_dir / "taqman.tsv", sep="\t", index=False)
    ercc.to_csv(bundle_dir / "ercc_counts.tsv", sep="\t", index=False)
    pd.DataFrame(truth_rows, columns=[
        "ercc_id", "subgroup", "conc_mix1", "conc_mix2",
        "expected_fc", "log2_mix1_mix2",
    ]).to_csv(bundle_dir / "ercc_truth.tsv", sep="\t", index=False)

    bundle = load_seqc_maqc_bundle(bundle_dir)
    m = score_ercc_dose_response(bundle)

    assert m["status"] == "pass", m
    fc = m["axes"]["fold_change_recovery"]
    dr = m["axes"]["dynamic_range"]
    assert fc["pearson"] >= 0.6
    assert fc["slope_measured_vs_expected"] > 0.5          # A=Mix1 direction
    assert set(fc["by_subgroup"]) == {"A", "B", "C", "D"}
    # Per-subgroup measured tracks expected ordering A > B > C > D.
    means = {g: v["measured_log2_mean"] for g, v in fc["by_subgroup"].items()}
    assert means["A"] > means["B"] > means["D"]
    assert dr["pearson_log_cpm_vs_log_conc"] >= 0.9
    assert dr["dynamic_range_orders_of_magnitude"] >= 2.0


@pytest.mark.skipif(
    not os.environ.get("ARIA_SEQC_MAQC_BUNDLE"),
    reason="set ARIA_SEQC_MAQC_BUNDLE to a staged real reference bundle",
)
def test_a1_seqc_maqc_real_bundle():
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.reference_seqc import run_seqc_maqc_a1_benchmark

    manifest = run_seqc_maqc_a1_benchmark(os.environ["ARIA_SEQC_MAQC_BUNDLE"])
    assert manifest["status"] in ("pass", "fail"), manifest
    assert manifest["axes"]["lfc_concordance"]["n_genes_scored"] > 100
    assert manifest["axes"]["taqman_de_detection"]["auc"] is not None
