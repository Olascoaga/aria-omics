"""Preprint-readiness audit B3: pairing is classified, never silently dropped.

Partially-paired designs must not degrade to an independent model in silence, and
scATAC differential accessibility must be able to model the donor block.

Policy (option B, 2026-07-12): for a partially-paired contrast the within-subject
block is modeled on the COMPLETE-PAIRED SUBSET (donors present in both conditions),
which is always estimable and controls the type-I error contributed by donor
baseline shifts; the excluded unpaired samples are disclosed, never hidden.

  * complete    — every replicate spans both contrast conditions -> block modeled on all.
  * partial     — some replicates span both, some are nested -> block modeled on the
                  paired subset; unpaired samples disclosed.
  * independent — no replicate spans both conditions -> unpaired model, no block.

The deterministic classifier runs everywhere; the DESeq2 wiring / type-I checks are
gated on the scientific environment via importorskip.
"""
from __future__ import annotations

import pytest


# ── Deterministic pairing classifier ────────────────────────────────────────

def _rows(mapping):
    """mapping: {sample: (condition, replicate)} -> list-of-dict rows."""
    return [
        {"sample": s, "condition": c, "replicate": r}
        for s, (c, r) in mapping.items()
    ]


def test_classify_pairing_complete():
    from aria.utils.design_matrix import classify_pairing

    rows = _rows({
        "s1": ("ctrl", "d1"), "s2": ("stim", "d1"),
        "s3": ("ctrl", "d2"), "s4": ("stim", "d2"),
    })
    rec = classify_pairing(rows, "condition", "replicate", test="stim", ref="ctrl")
    assert rec["status"] == "complete"
    assert set(rec["paired_replicates"]) == {"d1", "d2"}
    assert rec["unpaired_replicates"] == []
    assert set(rec["paired_samples"]) == {"s1", "s2", "s3", "s4"}
    assert rec["unpaired_samples"] == []


def test_classify_pairing_partial():
    from aria.utils.design_matrix import classify_pairing

    rows = _rows({
        "s1": ("ctrl", "d1"), "s2": ("stim", "d1"),   # paired
        "s3": ("ctrl", "d2"), "s4": ("stim", "d2"),   # paired
        "s5": ("ctrl", "d3"),                          # nested
        "s6": ("stim", "d4"),                          # nested
    })
    rec = classify_pairing(rows, "condition", "replicate", test="stim", ref="ctrl")
    assert rec["status"] == "partial"
    assert set(rec["paired_replicates"]) == {"d1", "d2"}
    assert set(rec["unpaired_replicates"]) == {"d3", "d4"}
    assert set(rec["paired_samples"]) == {"s1", "s2", "s3", "s4"}
    assert set(rec["unpaired_samples"]) == {"s5", "s6"}


def test_classify_pairing_independent():
    from aria.utils.design_matrix import classify_pairing

    rows = _rows({
        "s1": ("ctrl", "d1"), "s2": ("ctrl", "d2"),
        "s3": ("stim", "d3"), "s4": ("stim", "d4"),
    })
    rec = classify_pairing(rows, "condition", "replicate", test="stim", ref="ctrl")
    assert rec["status"] == "independent"
    assert rec["paired_replicates"] == []
    assert set(rec["unpaired_samples"]) == {"s1", "s2", "s3", "s4"}


def test_classify_pairing_restricts_to_contrast_levels():
    """A replicate paired only across OTHER conditions is not paired for this
    contrast; the classifier must consider only the two contrast levels."""
    from aria.utils.design_matrix import classify_pairing

    rows = _rows({
        "s1": ("ctrl", "d1"), "s2": ("other", "d1"),  # d1 not in stim
        "s3": ("stim", "d2"), "s4": ("ctrl", "d2"),    # d2 paired ctrl/stim
    })
    rec = classify_pairing(rows, "condition", "replicate", test="stim", ref="ctrl")
    assert rec["paired_replicates"] == ["d2"]
    assert "d1" in rec["unpaired_replicates"]
    # the 'other'-condition sample is outside the contrast entirely
    assert "s2" not in rec["paired_samples"]
    assert "s2" not in rec["unpaired_samples"]


# ── RNA pseudobulk: partial pairing models the block on the paired subset ────

# Two partially-paired layouts. Condition NEVER enters the expression rate, so
# there is no true condition effect anywhere; only donor baselines drive counts.
_DISCLOSE_LAYOUT = [
    ("d1", "ctrl"), ("d1", "stim"),   # paired
    ("d2", "ctrl"), ("d2", "stim"),   # paired
    ("d3", "ctrl"),                    # nested
    ("d4", "stim"),                    # nested
]
_DISCLOSE_BASE = {"d1": 0.0, "d2": 2.5, "d3": 5.0, "d4": 5.0}

# Type-I layout: two balanced paired donors (d1,d2 — clean null) plus nested
# donors CONFOUNDED with condition (c* ctrl-only, t* stim-only). The stim-only
# donors carry a gene-subset signature (see _make_typei_scrna), so an unpaired
# fit reads the donor confound as a condition effect -> false positives; the
# paired subset (d1,d2) is null.
_TYPEI_LAYOUT = [
    ("d1", "ctrl"), ("d1", "stim"),
    ("d2", "ctrl"), ("d2", "stim"),
    ("c1", "ctrl"), ("c2", "ctrl"), ("c3", "ctrl"),
    ("t1", "stim"), ("t2", "stim"), ("t3", "stim"),
]


def _make_scrna(np, ad, pd, layout, donor_base, *, seed=0, cells_per_sample=60):
    rng = np.random.default_rng(seed)
    n_genes = 200
    Xrows, obs = [], []
    for donor, cond in layout:
        base = donor_base[donor]
        for _ in range(cells_per_sample):
            lam = np.exp(base / 5.0 + rng.normal(0, 0.1, n_genes))
            Xrows.append(rng.poisson(lam))
            obs.append({"cell_type": "T", "condition": cond, "donor": donor})
    X = np.asarray(Xrows, dtype=float)
    obs_df = pd.DataFrame(obs)
    obs_df.index = [f"cell{i}" for i in range(len(obs_df))]
    adata = ad.AnnData(X=X, obs=obs_df)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def _make_typei_scrna(np, ad, pd, *, seed=0, cells_per_sample=60):
    """Type-I fixture with a GENE-SPECIFIC donor confound (a uniform shift would
    be normalized away by DESeq2 size factors). Genes 0:100 carry a signature
    expressed only by the stim-only nested donors, so an unpaired fit reads it as
    a stim effect; the paired donors d1,d2 are identical across conditions (null)."""
    rng = np.random.default_rng(seed)
    n_genes = 200
    sig_genes = 100
    Xrows, obs = [], []
    for donor, cond in _TYPEI_LAYOUT:
        signature = donor.startswith("t")     # stim-only nested donors
        for _ in range(cells_per_sample):
            mu = np.full(n_genes, 3.0)
            if signature:
                mu[:sig_genes] += 5.0          # gene-subset confound
            lam = np.exp(mu / 5.0 + rng.normal(0, 0.1, n_genes))
            Xrows.append(rng.poisson(lam))
            obs.append({"cell_type": "T", "condition": cond, "donor": donor})
    X = np.asarray(Xrows, dtype=float)
    obs_df = pd.DataFrame(obs)
    obs_df.index = [f"cell{i}" for i in range(len(obs_df))]
    adata = ad.AnnData(X=X, obs=obs_df)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def _run_partial_pseudobulk(np, ad, pd, tmp_path, *, layout=None, donor_base=None,
                            adata=None, seed=0, auto_paired=True):
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    if adata is None:
        adata = _make_scrna(np, ad, pd, layout, donor_base, seed=seed)
    data_path = str(tmp_path / f"data_{seed}_{auto_paired}.h5ad")
    adata.write_h5ad(data_path)
    res = rna_pseudobulk_de({
        "data_path": data_path,
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["stim", "ctrl"]],
        "min_replicates_per_condition": 2,
        "min_cells_per_pseudosample": 10,
        "use_raw": False,
        "auto_paired_donor_covariate": auto_paired,
        "output_dir": str(tmp_path / f"out_{seed}_{auto_paired}"),
    })
    return res["per_group"]["T"]["per_comparison"]["stim_vs_ctrl"]


def test_rna_partial_pairing_discloses_and_models_block(tmp_path):
    np = pytest.importorskip("numpy")
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")

    entry = _run_partial_pseudobulk(
        np, ad, pd, tmp_path,
        layout=_DISCLOSE_LAYOUT, donor_base=_DISCLOSE_BASE,
    )
    assert entry.get("pairing_status") == "partial"
    assert entry.get("paired_block_modeled") is True
    # d3/d4 samples are disclosed as excluded, not silently used.
    excluded = entry.get("excluded_unpaired_samples") or []
    assert any("d3" in s for s in excluded)
    assert any("d4" in s for s in excluded)


def test_rna_partial_pairing_controls_type_i_vs_unpaired(tmp_path):
    """The paired-subset block must not produce MORE null false positives than the
    silent unpaired model it replaces (type-I control)."""
    np = pytest.importorskip("numpy")
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")

    blocked = _run_partial_pseudobulk(
        np, ad, pd, tmp_path, seed=1,
        adata=_make_typei_scrna(np, ad, pd, seed=1), auto_paired=True,
    )
    unpaired = _run_partial_pseudobulk(
        np, ad, pd, tmp_path, seed=1,
        adata=_make_typei_scrna(np, ad, pd, seed=1), auto_paired=False,
    )
    blocked_hits = int(blocked.get("n_significant", 0))
    unpaired_hits = int(unpaired.get("n_significant", 0))
    assert blocked.get("pairing_status") == "partial"
    # The silent unpaired fit reads the donor confound as a condition effect and
    # produces null false positives; the paired-subset block controls type-I.
    assert unpaired_hits > 0
    assert blocked_hits < unpaired_hits


# ── Chromatin DA: donor-aware pseudobulk ────────────────────────────────────

def test_chromatin_pseudobulk_da_is_donor_aware():
    np = pytest.importorskip("numpy")
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from aria.scripts.chromatin_diffacc import _pseudobulk_da

    rng = np.random.default_rng(0)
    layout = [
        ("d1", "ctrl"), ("d1", "stim"),
        ("d2", "ctrl"), ("d2", "stim"),
        ("d3", "ctrl"), ("d3", "stim"),
    ]
    donor_base = {"d1": 0.0, "d2": 3.0, "d3": 6.0}
    n_peaks = 120
    rows, obs = [], []
    for donor, cond in layout:
        for _ in range(40):
            lam = np.exp(donor_base[donor] / 6.0 + rng.normal(0, 0.1, n_peaks))
            rows.append(rng.poisson(lam))
            obs.append({"condition": cond, "donor": donor})
    X = np.asarray(rows, dtype=float)
    obs_df = pd.DataFrame(obs)
    obs_df.index = [f"c{i}" for i in range(len(obs_df))]
    adata = ad.AnnData(X=X, obs=obs_df)
    adata.var_names = [f"peak{i}" for i in range(n_peaks)]

    warnings: list = []
    res = _pseudobulk_da(
        adata, "condition", "donor",
        comparisons=[{"test": "stim", "reference": "ctrl"}],
        padj_max=0.1, lfc_min=0.0, top_n=50,
        min_cells=10, min_reps=2, allow_mock=False, warnings=warnings,
    )
    assert res["ran"] is True
    comp = res["comparisons"][0]
    # A complete-paired scATAC design must model the donor block, not ~condition.
    assert comp.get("pairing_status") == "complete"
    assert comp.get("paired_block_modeled") is True
