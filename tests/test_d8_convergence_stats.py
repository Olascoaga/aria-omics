"""D8 convergence regression tests for statistical helpers (D8-016/D2-BUG-1/2, D8-014, D8-015, D8-024)."""

from __future__ import annotations

import numpy as np


def test_d8_016_bh_fdr_nan_poisoning_divergence():
    """D8-016 / D2-BUG-1 / D2-BUG-2: every BH helper must converge on stats.bh_correct
    and preserve NaNs without poisoning finite values."""
    from aria.benchmarks.synthetic_de import _bh_adjust as synth_bh
    from aria.scripts.chromatin_footprint_tobias import _bh_adjust as tobias_bh
    from aria.scripts.hic_topology import _bh_fdr as hic_bh
    from aria.scripts.rna_pseudobulk_de import _global_bh as pb_bh
    from aria.utils.cross_modal_concordance import bh_adjust as xmodal_bh
    from aria.utils.lobo_jackknife import bh_adjust as lobo_bh
    from aria.utils.stats import bh_correct as canonical_bh

    pvals = [0.001, 0.02, float("nan"), 0.4, 0.8]
    ref = canonical_bh(pvals)
    assert np.allclose(ref, [0.004, 0.04, np.nan, 0.5333333333333333, 0.8], equal_nan=True)
    copies = {
        "lobo_jackknife.bh_adjust": np.asarray(lobo_bh(np.asarray(pvals)), dtype=float),
        "cross_modal_concordance.bh_adjust": np.asarray(xmodal_bh(pvals), dtype=float),
        "hic_topology._bh_fdr": np.asarray(hic_bh(pvals), dtype=float),
        "chromatin_footprint_tobias._bh_adjust": np.asarray(tobias_bh(pvals), dtype=float),
        "rna_pseudobulk_de._global_bh": np.asarray(pb_bh(pvals), dtype=float),
        "synthetic_de._bh_adjust": np.asarray(synth_bh(pvals), dtype=float),
    }
    mismatches = {
        k: v.tolist()
        for k, v in copies.items()
        if not np.allclose(v, ref, equal_nan=True, atol=1e-12)
    }
    assert not mismatches, f"BH copies disagree with stats.bh_correct ({ref.tolist()}): {mismatches}"
    assert pb_bh is canonical_bh
    assert synth_bh is canonical_bh
    assert lobo_bh is canonical_bh
    assert xmodal_bh is canonical_bh
    assert hic_bh is canonical_bh


def test_d8_014_plot_neglog10_underflow_cap_converges(tmp_path):
    """D8-014: all -log10(p) underflow cappers must floor cap >= 5.0 (strictly > -log10(0.05)) and emit float NaN."""
    import math
    from pathlib import Path
    import pandas as pd
    import pytest

    pytest.importorskip("litellm")
    from aria.agents._narrative_chromatin import _neglog10_capped
    from aria.agents._narrative_integration import _neglog10_display
    from aria.scripts.rna_bulk.plots import _generate_plots
    from aria.utils.cross_modal_concordance import _adaptive_neglog10

    pvals = [0.0, 0.2, 0.8]
    _, cap_xm, _ = _adaptive_neglog10(pvals)
    _, _, cap_chr = _neglog10_capped(pvals)
    _, cap_int = _neglog10_display(pvals)
    assert cap_xm == 5.0 and cap_chr == 5.0 and cap_int == 5.0

    res_df = pd.DataFrame(
        {"log2FoldChange": [3.5, 0.2, -0.1], "padj": pvals},
        index=["GENE_ZERO", "GENE_B", "GENE_C"],
    )
    out = _generate_plots(
        {"results": res_df, "n_up": 1, "n_down": 0, "sig_genes": ["GENE_ZERO"]},
        {},
        None,
        None,
        "condition",
        str(tmp_path),
        padj_thr=0.05,
        lfc_thr=1.0,
    )
    svg = Path(out["volcano"]).read_text(encoding="utf-8")
    y_int_nan, _ = _neglog10_display([math.nan, 0.01])
    assert "shown at 5.00" in svg and "shown at 1.00" not in svg
    assert all(v is not None for v in y_int_nan)
    assert math.isnan(y_int_nan[0]) and y_int_nan[1] == pytest.approx(2.0)

    from aria.utils.stats import adaptive_neglog10

    assert _adaptive_neglog10 is adaptive_neglog10


def test_d8_015_cpm_and_logcpm_zero_depth_converges(tmp_path):
    """D8-015: CPM and log2(CPM+1) helpers must share full-library zero-depth handling without NaN crashes."""
    import pandas as pd
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.audit_agent import AuditAgent
    from aria.scripts.chromatin_figure_bulk import _pca as chrom_pca
    from aria.scripts.rna_bulk.transforms import _compute_cpm

    genes = [f"g{i}" for i in range(15)]
    counts = pd.DataFrame(
        {
            "ctrl_1": [100 + i * 5 for i in range(15)],
            "ctrl_2": [110 + i * 5 for i in range(15)],
            "treat_1": [400 - i * 5 for i in range(15)],
            "treat_2": [0] * 15,
        },
        index=genes,
    )
    tsv = tmp_path / "counts.tsv"
    counts.to_csv(tsv, sep="\t")

    cpm_bulk = _compute_cpm(counts, [])
    pca_res, pca_reason = chrom_pca(counts, is_cpm=False)
    assert cpm_bulk is not None and np.isfinite(cpm_bulk.to_numpy()).all()
    assert pca_res is None and pca_reason == "a sample has zero total counts"

    design = {
        "batch_factor": "batch",
        "groups": {"ctrl": ["ctrl_1", "ctrl_2"], "treat": ["treat_1", "treat_2"]},
        "batch_map": {"ctrl_1": "B1", "ctrl_2": "B2", "treat_1": "B1", "treat_2": "B2"},
    }
    agent = AuditAgent.__new__(AuditAgent)
    rep_findings = agent._check_replicate_correlation([str(tsv)], design)
    pca_findings = agent._check_pca_batch_dominance([str(tsv)], design)

    assert any("treat_2" in f.get("message", "") for f in rep_findings)
    assert isinstance(pca_findings, list)
    assert any("treat_2" in f.get("message", "") for f in pca_findings)
