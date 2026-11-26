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
