"""P1-1 (c): FDR across the bulk contrast family.

Bulk DE controlled FDR within each contrast (per-contrast BH). When several
contrasts are tested together, the family-wise error across the contrast family
is not controlled. This adds a pooled BH across ALL contrasts (the "contrast
family"), pre-registered like P1-2 (per_contrast default | global), recorded
additively. The core logic lives in pure stats helpers so it is testable without
pydeseq2; the W-CALIB bulk benchmark guards the DE math.
"""

import pytest

from aria.utils.stats import (
    contrast_family_significance,
    pooled_bh_across_groups,
    preregister_contrast_family,
)


def test_pooled_bh_is_more_conservative_than_per_group():
    import numpy as np
    from aria.utils.stats import bh_correct
    groups = {
        "c1": {"g1": 0.001, "g2": 0.02, "g3": 0.5},
        "c2": {"g1": 0.001, "g2": 0.04, "g3": 0.9},
    }
    pooled = pooled_bh_across_groups(groups)
    # Per-group BH on c1 alone vs pooled across both: pooling tests more
    # hypotheses, so each adjusted p-value is >= the per-group one.
    per_c1 = bh_correct(list(groups["c1"].values()))
    per_c1_map = dict(zip(groups["c1"], per_c1))
    for g in groups["c1"]:
        assert pooled["c1"][g] >= per_c1_map[g] - 1e-9
    assert pooled_bh_across_groups({}) == {}


def test_preregister_contrast_family_normalizes_and_declares():
    d = preregister_contrast_family("GLOBAL")
    assert d["fdr_family"] == "global"
    assert d["preregistered"] is True
    assert preregister_contrast_family("junk")["fdr_family"] == "per_contrast"
    assert preregister_contrast_family(None)["fdr_family"] == "per_contrast"


def test_contrast_family_significance_pools_and_gates_on_lfc():
    group_stats = {
        "B_vs_A": {
            "g1": {"pvalue": 1e-6, "log2fc": 3.0},   # strong, significant
            "g2": {"pvalue": 1e-6, "log2fc": 0.1},   # significant p but |lfc| too small
            "g3": {"pvalue": 0.5,  "log2fc": 2.0},   # not significant
        },
        "C_vs_A": {
            "g1": {"pvalue": 1e-6, "log2fc": -2.5},  # strong down
            "g4": {"pvalue": 0.9,  "log2fc": 0.0},
        },
    }
    fam = contrast_family_significance(group_stats, padj_max=0.05, lfc_min=1.0)
    assert fam["B_vs_A"]["sig_genes"] == ["g1"]       # g2 fails |lfc|, g3 fails padj
    assert fam["B_vs_A"]["n_up"] == 1 and fam["B_vs_A"]["n_down"] == 0
    assert fam["C_vs_A"]["sig_genes"] == ["g1"]
    assert fam["C_vs_A"]["n_down"] == 1
    # padj_family is the pooled adjustment (>= raw p-value).
    assert fam["B_vs_A"]["padj_family"]["g1"] >= 1e-6


def test_bulk_de_global_contrast_family_is_recorded(tmp_path):
    pytest.importorskip("pydeseq2")
    import numpy as np
    import pandas as pd
    from aria.scripts.rna_bulk_de import bulk_rna_de

    rng = np.random.default_rng(5)
    groups = {"A": 6, "B": 6, "C": 6}
    samples, conds = [], []
    for g, n in groups.items():
        for r in range(n):
            samples.append(f"{g}{r}")
            conds.append(g)
    n_genes = 200
    # Per-gene base (lognormal) so samples share structure and sample-QC does not
    # over-flag synthetic replicates; strong DE in distinct genes per condition.
    base = np.clip(np.exp(rng.normal(4.0, 0.8, size=n_genes)), 20.0, None)
    fc = np.ones((n_genes, len(samples)))
    fc[:10, 6:12] = 4.0    # DE up in B (samples 6..11)
    fc[10:20, 12:] = 4.0   # DE up in C (samples 12..17)
    counts = rng.poisson(base[:, None] * fc)
    counts_df = pd.DataFrame(np.asarray(counts).astype(int),
                             index=[f"g{i}" for i in range(n_genes)], columns=samples)
    counts_path = tmp_path / "counts.tsv"
    counts_df.to_csv(counts_path, sep="\t")
    meta_path = tmp_path / "meta.tsv"
    pd.DataFrame({"condition": conds}, index=samples).to_csv(meta_path, sep="\t")

    res = bulk_rna_de({
        "files": [str(counts_path)],
        "metadata_file": str(meta_path),
        "design_factor": "condition",
        "contrasts": [{"numerator": "B", "denominator": "A"},
                      {"numerator": "C", "denominator": "A"}],
        "fdr_family": "global",
        "run_pathways": False,
        "padj_threshold": 0.05,
        "lfc_threshold": 0.5,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success", res
    assert res["fdr_family"]["fdr_family"] == "global"
    assert res["fdr_family"]["preregistered"] is True
    succ = [c for c in res["contrasts"] if c.get("status") == "success"]
    assert succ, res
    for c in succ:
        assert "n_significant_contrast_family" in c
