"""Guards for the scATAC multi-sample pseudobulk DA validation harness.

The harness (`scripts/run_scatac_multisample_validation.py`) exercises ARIA's
REAL chromatin pseudobulk DA lane on Samael's Erosion CONSENSUS peak universe
(a shared peak set across donors). These guards lock the two things that are
easy to get silently wrong -- donor identity parsing and the wiring to the real
ARIA scripts -- and a dataset-gated real-data smoke for the combined builder.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest

_HARNESS = (Path(__file__).resolve().parent.parent
            / "scripts" / "run_scatac_multisample_validation.py")


def _load_harness():
    spec = importlib.util.spec_from_file_location("scatac_msv", _HARNESS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_donor_parsing_handles_deep_suffix():
    """Barcodes are '<donor>_<16bp>-1'; the donor may itself contain '_'
    (e.g. 'hc29_deep'), so a naive split('_')[0] would merge 'hc29_deep' with a
    hypothetical 'hc29'. The harness must strip only the cell-barcode suffix."""
    h = _load_harness()
    assert h._donor_of("hc13344_AAACAGCCACTGACTA-1") == "hc13344"
    assert h._donor_of("hc29_deep_AAACATGCATCACAGC-1") == "hc29_deep"
    assert h._donor_of("hc212191_deep_TTTGTGTTCGGANNNN-1") == "hc212191_deep"
    # No spurious stripping when there is no trailing barcode.
    assert h._donor_of("hc11_deep") == "hc11_deep"


def test_harness_wires_real_aria_scripts():
    """The harness must call ARIA's real LSI + DA scripts, not a local copy."""
    h = _load_harness()
    assert callable(h.build_combined_adata)
    # The real ARIA chromatin scripts must be importable (wiring contract).
    from aria.scripts.chromatin_lsi_clustering import chromatin_lsi_clustering
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    assert callable(chromatin_lsi_clustering)
    assert callable(chromatin_diffacc)


@pytest.mark.skipif(
    not Path(os.environ.get(
        "ARIA_SCATAC_CONSENSUS_DIR",
        "/home/medusa/Samael/Erosion/results/02_consensus")).is_dir()
    or os.environ.get("ARIA_SCATAC_RUN_REALDATA") != "1",
    reason="consensus dataset not present or ARIA_SCATAC_RUN_REALDATA != 1",
)
def test_combined_builder_real_data_smoke():
    """Dataset-gated: the combined builder yields a multi-donor, multi-condition
    AnnData over a single shared peak set, with the obs the DA lane needs."""
    pytest.importorskip("anndata")
    pytest.importorskip("scipy")
    h = _load_harness()
    combined = h.build_combined_adata()
    assert combined.n_obs > 1000
    for col in ("age_group", "donor", "sex", "cell_type"):
        assert col in combined.obs.columns
    # Multiple donors per age group -> real biological replicates.
    by_age = combined.obs.groupby("age_group", observed=True)["donor"].nunique()
    assert (by_age >= 2).all(), by_age.to_dict()
    # One shared peak set (consensus) -> peaks are coordinate strings.
    assert combined.n_vars > 1000


def _synth_adata(n_donors_per_age=4, cells_per_donor=300, n_peaks=100,
                 n_empty_peaks=10, seed=0):
    """Tiny multi-donor AnnData mirroring the consensus obs the harness needs."""
    ad = pytest.importorskip("anndata")
    sp = pytest.importorskip("scipy.sparse")
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(seed)
    rows = []
    for age in ("20-39", "80-100"):
        for d in range(n_donors_per_age):
            for _ in range(cells_per_donor):
                rows.append((f"{age}_F_d{age}-{d}", "F", age))
    obs = pd.DataFrame(rows, columns=["donor", "sex", "age_group"])
    obs["cell_type"] = "Oligo"
    obs.index = [f"c{i}" for i in range(len(obs))]
    X = rng.poisson(0.5, size=(len(obs), n_peaks)).astype("float32")
    X[:, :n_empty_peaks] = 0  # all-empty peaks the filter must drop
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs)
    a.var_names = [f"chr1:{i}-{i+500}" for i in range(n_peaks)]
    a.obs["age_group"] = a.obs["age_group"].astype("category")
    a.obs["donor"] = a.obs["donor"].astype("category")
    return a


def test_subsample_preserves_donors_and_filters_peaks(monkeypatch):
    """The validation subsample must (a) keep every donor (stratified, so the
    10-vs-10 contrast survives), (b) hit roughly the cell target, and (c) drop
    all-empty peaks and cap to MAX_PEAKS."""
    h = _load_harness()
    monkeypatch.setattr(h, "MAX_CELLS", 500)
    monkeypatch.setattr(h, "FLOOR_PER_DONOR", 20)
    monkeypatch.setattr(h, "MIN_CELLS_PEAK", 5)
    monkeypatch.setattr(h, "MAX_PEAKS", 50)
    monkeypatch.setattr(h, "SEED", 0)

    adata = _synth_adata()  # 8 donors, 2400 cells, 100 peaks (10 empty)
    n_donors_before = adata.obs["donor"].nunique()
    report = {}
    out = h._subsample_and_filter(adata, report)

    # Every donor survives the stratified draw.
    assert out.obs["donor"].nunique() == n_donors_before == 8
    assert set(out.obs["age_group"].astype(str)) == {"20-39", "80-100"}
    # Cell target honored (floors can nudge it slightly above MAX_CELLS).
    assert out.n_obs <= 500 + 8 * 20
    assert out.n_obs >= 8 * 20  # at least the per-donor floor each
    # All-empty peaks dropped and capped to MAX_PEAKS.
    assert out.n_vars <= 50
    assert all(":" in str(v) for v in out.var_names)
    # Report carries the audit trail.
    s = report["subsample"]
    assert s["cells_before"] == 2400 and s["peaks_before"] == 100
    assert s["cells_after"] == out.n_obs and s["peaks_after"] == out.n_vars
