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


def test_presence_filter_drops_present_absent_and_rare_peaks(monkeypatch):
    """Lever B keeps peaks present in >=50% of donors in BOTH contrast groups,
    dropping (a) present/absent peaks (open in one age group, zero in the other)
    whose extreme |LFC| is a pseudocount artifact, and (b) rare peaks present in
    too few donors. A peak shared across both arms must survive."""
    ad = pytest.importorskip("anndata")
    sp = pytest.importorskip("scipy.sparse")
    import numpy as np
    import pandas as pd

    h = _load_harness()
    monkeypatch.setattr(h, "TEST_AGE", "80-100")
    monkeypatch.setattr(h, "REF_AGE", "20-39")
    monkeypatch.setattr(h, "MIN_DONOR_FRAC", 0.5)
    monkeypatch.setattr(h, "MIN_CELLS_DONOR_PRESENT", 1)

    cells_per_donor = 3
    rows, donor_of_cell, age_of_cell = [], [], []
    for age in ("80-100", "20-39"):
        for d in range(4):
            for _ in range(cells_per_donor):
                rows.append(f"{age}-d{d}")
                age_of_cell.append(age)
    n = len(rows)
    X = np.zeros((n, 3), dtype="float32")
    age_arr = np.array(age_of_cell)
    donor_arr = np.array(rows)
    X[:, 0] = 1                                   # peak0 shared -> kept
    X[age_arr == "80-100", 1] = 1                 # peak1 present/absent -> dropped
    X[(donor_arr == "80-100-d0") | (donor_arr == "20-39-d0"), 2] = 1  # rare -> dropped
    obs = pd.DataFrame({"donor": donor_arr, "age_group": age_arr})
    obs["donor"] = obs["donor"].astype("category")
    obs["age_group"] = obs["age_group"].astype("category")
    obs.index = [f"c{i}" for i in range(n)]
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs)
    a.var_names = ["shared", "present_absent", "rare"]

    out = h._presence_filter(a, {})
    assert list(out.var_names) == ["shared"]

    # Disabled (frac<=0) is a pass-through.
    monkeypatch.setattr(h, "MIN_DONOR_FRAC", 0.0)
    assert h._presence_filter(a, {}).n_vars == 3


def test_unify_peaks_merges_overlapping_and_sums_counts(monkeypatch):
    """Lever A collapses boundary-shifted per-donor duplicates of one region into
    a single column and sums their counts. Cell Ranger ARC calls peaks per sample,
    so the same region appears as chr1_100_200 in donor A and chr1_103_205 in
    donor B (exact-string disjoint, genomically overlapping); exact-string
    consensus matching left them as separate donor-specific columns (the block
    structure). Overlap-merge must (a) unite genomically overlapping peaks, (b)
    keep disjoint peaks separate, (c) drop non-canonical contigs, and (d) sum the
    counts of merged columns so the region's signal is donor-comparable."""
    ad = pytest.importorskip("anndata")
    sp = pytest.importorskip("scipy.sparse")
    import numpy as np

    h = _load_harness()
    monkeypatch.setattr(h, "UNIFY_PEAKS", True)

    # 2 cells x 4 peaks: two overlapping (same region, shifted) + one disjoint on
    # the same chrom + one on a non-canonical contig (must be dropped).
    peaks = ["chr1_100_200", "chr1_150_260", "chr1_5000_5100", "GL000009.2_10_90"]
    X = np.array([[1, 2, 5, 7],
                  [0, 3, 0, 9]], dtype="float32")
    a = ad.AnnData(X=sp.csr_matrix(X))
    a.var_names = peaks

    out = h._unify_peaks_by_overlap(a, {})
    # chr1_100_200 & chr1_150_260 overlap -> one merged interval chr1_100_260;
    # chr1_5000_5100 stays separate; GL000009.2 dropped (non-canonical).
    assert list(out.var_names) == ["chr1_100_260", "chr1_5000_5100"]
    # Merged column sums the two overlapping peaks per cell; disjoint untouched.
    merged = out.X.toarray() if sp.issparse(out.X) else np.asarray(out.X)
    assert merged[0, 0] == 3 and merged[1, 0] == 3   # (1+2), (0+3)
    assert merged[0, 1] == 5 and merged[1, 1] == 0   # chr1_5000_5100 unchanged

    # Disabled (UNIFY_PEAKS=0) is a pass-through.
    monkeypatch.setattr(h, "UNIFY_PEAKS", False)
    assert h._unify_peaks_by_overlap(a, {}).n_vars == 4


def test_diffacc_skip_per_cluster_runs_pseudobulk_only(monkeypatch, tmp_path):
    """skip_per_cluster=True must bypass the single-threaded descriptive
    per-cluster Wilcoxon lane (reporting it honestly as not-run) while the
    pseudobulk lane still dispatches. Default (absent) keeps per-cluster on.
    Both heavy lanes are stubbed so this stays a fast, deps-light gating guard."""
    ad = pytest.importorskip("anndata")
    sp = pytest.importorskip("scipy.sparse")
    import numpy as np
    import pandas as pd
    from aria.scripts import chromatin_diffacc as cd

    # Tiny clustered AnnData with the columns chromatin_diffacc requires.
    n = 8
    obs = pd.DataFrame({
        "leiden": ["0", "1"] * (n // 2),
        "age_group": ["80-100"] * (n // 2) + ["20-39"] * (n // 2),
        "donor": [f"d{i}" for i in range(n)],
    }, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=sp.csr_matrix(np.ones((n, 5), dtype="float32")), obs=obs)
    a.var_names = [f"chr1_{i}_100" for i in range(5)]
    data_path = tmp_path / "clustered.h5ad"
    a.write_h5ad(str(data_path))

    calls = {"per_cluster": 0, "pseudobulk": 0}

    def _spy_pc(*args, **kwargs):
        calls["per_cluster"] += 1
        return {"ran": True, "reason": None, "n_da_total": 0,
                "n_da_by_cluster": {}, "da_peaks_by_cluster": {},
                "_rows": [], "n_clusters": 2}

    def _stub_pb(*args, **kwargs):
        calls["pseudobulk"] += 1
        return {"ran": True, "reason": None, "comparisons": [], "_rows": []}

    monkeypatch.setattr(cd, "_per_cluster_accessibility", _spy_pc)
    monkeypatch.setattr(cd, "_pseudobulk_da", _stub_pb)

    base = {"data_path": str(data_path), "groupby": "leiden",
            "condition_col": "age_group", "replicate_col": "donor",
            "comparisons": [{"test": "80-100", "reference": "20-39"}],
            "output_dir": str(tmp_path)}

    # Skipped: per-cluster spy NOT called, pseudobulk still runs.
    out = cd.chromatin_diffacc({**base, "skip_per_cluster": True})
    assert out["status"] == "success"
    assert calls["per_cluster"] == 0 and calls["pseudobulk"] == 1
    assert out["per_cluster"]["ran"] is False
    assert "skip_per_cluster" in out["per_cluster"]["reason"]

    # Default (absent): per-cluster lane runs as before.
    out2 = cd.chromatin_diffacc(base)
    assert calls["per_cluster"] == 1 and calls["pseudobulk"] == 2
    assert out2["per_cluster"]["ran"] is True
