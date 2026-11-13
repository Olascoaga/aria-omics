"""V48 Slice I4 — Paired MOFA+ interpretable factors & external concordance tests.

Runs in two modes:
  1. `aria-env` (`pytest tests/test_v48_i4_mofa.py`): exercises all pure NumPy/SciPy
     concordance math (Hungarian 1-to-1 active-factor matching, reference-R^2
     weighting, rotation-invariant canonical subspace correlations, loading sign
     consistency, deterministic stratified subsampling with rare-population floor,
     hybrid absolute + view-relative ATAC/RNA activity criterion, IPC contract
     registration, IntegrationAgent confounder conflict surfacing, and AST
     no-fabrication / no-embedded-gene-panel guards).
  2. `aria-integration-env` (`python tests/test_v48_i4_mofa.py`): additionally
     executes full end-to-end `integration_mofa` with real `mofapy2` + `mudata`
     over planted paired `.h5mu` and `.h5ad` cohorts (shared cross-modal factor,
     modality-specific factor, sequencing-depth technical factor, `.obs` and
     caller-supplied S/G2M cell-cycle factor, stratified subsampling + exact
     closed-form posterior projection without `NaN`s, and `seed=0` determinism).
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.concordance_mofa import (
    is_factor_active_in_view,
    match_factors_hungarian,
    pearson_correlation,
    score_loading_concordance,
    score_mofa_concordance,
    score_mofa_seed_stability,
    spearman_correlation,
    subspace_canonical_concordance,
)
from aria.scripts.integration_mofa import (
    _deterministic_subsample,
    _mock_mofa,
    integration_mofa,
)


def _has_mofa_stack() -> bool:
    try:
        import anndata  # noqa: F401
        import mofapy2  # noqa: F401
        import mudata  # noqa: F401
        import scanpy  # noqa: F401

        return True
    except ImportError:
        return False


# ── Pure mathematical & governance tests (run in aria-env AND integration-env) ─


def test_hybrid_view_activity_protects_sparse_atac():
    # ATAC total R^2 = 0.08; a factor explaining 0.006 (7.5% of ATAC R^2) is active!
    assert is_factor_active_in_view(0.006, 0.08, active_r2_min=0.01, rel_r2_min=0.05)
    # Realistic sparse scATAC non-depth R^2 sum = 0.0037; a lineage factor with
    # r2_atac = 0.0010 (27% of biological ATAC variance, >= 0.0008 floor) is active!
    assert is_factor_active_in_view(0.0010, 0.0037, active_r2_min=0.01, rel_r2_min=0.05)
    # Noise below min_floor_r2=0.0008 stays inactive even if view total R^2 is small.
    assert not is_factor_active_in_view(0.0004, 0.0037, active_r2_min=0.01, rel_r2_min=0.05)
    # Standard absolute threshold >= 0.01 is always active.
    assert is_factor_active_in_view(0.025, 0.45, active_r2_min=0.01, rel_r2_min=0.05)


def test_weighted_pearson_and_spearman_properties():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    assert abs(pearson_correlation(x, y) - 1.0) < 1e-12
    assert abs(pearson_correlation(x, -y, weights=[1, 2, 3, 2, 1]) + 1.0) < 1e-12
    assert abs(spearman_correlation(x, np.exp(x)) - 1.0) < 1e-12
    assert pearson_correlation([3.0, 3.0, 3.0], [1.0, 2.0, 3.0]) is None


def test_hungarian_matching_permutation_sign_and_inactive_filtering():
    rng = np.random.default_rng(19)
    n = 200
    z1 = rng.normal(size=n)
    z2 = rng.normal(size=n)
    z3 = rng.normal(size=n)
    noise_a = rng.normal(scale=1e-4, size=n)
    noise_r = rng.normal(scale=1e-4, size=n)

    ref_scores = np.column_stack([z1, z2, z3, noise_r])
    # Permute [z3, -z1, z2, noise_a] with sign flip on z1
    aria_scores = np.column_stack([z3, -z1, z2, noise_a])

    ref_r2 = {"rna": [0.20, 0.10, 0.001, 0.0], "atac": [0.08, 0.0, 0.04, 0.0]}
    aria_r2 = {"rna": [0.001, 0.20, 0.10, 0.0], "atac": [0.04, 0.08, 0.0, 0.0]}

    res = match_factors_hungarian(
        aria_scores,
        ref_scores,
        aria_r2_per_view=aria_r2,
        ref_r2_per_view=ref_r2,
    )
    assert res["n_active_ref"] == 3
    assert res["n_active_aria"] == 3
    assert res["n_matched"] == 3
    assert res["n_unmatched_ref"] == 0
    assert res["mean_matched_abs_r"] > 0.999
    assert res["ref_r2_weighted_abs_r"] > 0.999

    # Verify ref Factor 1 (idx 0) matched aria Factor 2 (idx 1) with sign -1
    m_f1 = next(m for m in res["matches"] if m["ref_factor_id"] == 1)
    assert m_f1["aria_factor_id"] == 2
    assert m_f1["sign"] == -1


def test_subspace_canonical_concordance_is_rotation_invariant():
    rng = np.random.default_rng(42)
    z_ref = rng.normal(size=(250, 3))
    # Apply an arbitrary 3x3 orthogonal rotation Q (splits factor 1 & 2 at 45 deg!)
    q_orth, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    z_rot = z_ref @ q_orth

    # 1-to-1 Hungarian degrades under 45-deg rotation, whereas subspace stays 1.0!
    sub = subspace_canonical_concordance(z_rot, z_ref)
    assert sub["n_dims_aria"] == 3
    assert sub["n_dims_ref"] == 3
    assert sub["mean_canonical_corr"] == 1.0
    assert sub["min_top3_canonical_corr"] == 1.0


def test_loading_concordance_sign_alignment_and_overlap_guard():
    rng = np.random.default_rng(101)
    feats = [f"gene_{i}" for i in range(120)]
    w_ref = rng.normal(size=(120, 2))
    # Flip sign on factor 0 to match a sign=-1 score flip
    w_aria = np.column_stack([-w_ref[:, 0], w_ref[:, 1]])
    matches = [
        {"aria_factor_idx": 0, "aria_factor_id": 1, "ref_factor_idx": 0, "ref_factor_id": 1, "sign": -1},
        {"aria_factor_idx": 1, "aria_factor_id": 2, "ref_factor_idx": 1, "ref_factor_id": 2, "sign": 1},
    ]

    lc = score_loading_concordance(
        {"scRNA": (feats, w_aria)},
        {"rna": (feats, w_ref)},
        matches,
        min_shared_features=50,
        top_k=25,
    )
    assert lc["views"]["rna"]["status"] == "success"
    assert lc["views"]["rna"]["mean_pearson_r"] == 1.0
    assert lc["views"]["rna"]["mean_top_k_jaccard"] == 1.0
    assert lc["views"]["rna"]["sign_consistency_rate"] == 1.0

    # Insufficient feature overlap -> honest not_computed
    lc_low = score_loading_concordance(
        {"atac": ([f"chr1:{i}-{i+10}" for i in range(20)], w_aria[:20, :])},
        {"atac": ([f"chr1:{i}-{i+10}" for i in range(10, 30)], w_ref[:20, :])},
        matches,
        min_shared_features=50,
    )
    assert lc_low["views"]["atac"]["status"] == "insufficient_feature_overlap"
    assert "loading_concordance:atac" in lc_low["metrics_not_computed"]


def test_concordance_not_run_and_barcode_alignment():
    nr = score_mofa_concordance(
        aria_barcodes=["c1", "c2"],
        aria_scores=np.zeros((2, 2)),
        aria_r2_total={"rna": 0.3},
        ref_barcodes=None,
        ref_scores=None,
    )
    assert nr["status"] == "not_run"

    rng = np.random.default_rng(5)
    bcs = [f"BC_{i:03d}" for i in range(80)]
    z = rng.normal(size=(80, 2))
    perm = rng.permutation(80)
    bcs_perm = [bcs[i] for i in perm]
    z_perm = z[perm, :]

    full = score_mofa_concordance(
        aria_barcodes=bcs,
        aria_scores=z,
        aria_r2_total={"rna": 0.35, "atac": 0.14},
        ref_barcodes=bcs_perm,
        ref_scores=z_perm,
        ref_r2_total={"rna": 0.35, "atac": 0.14},
    )
    assert full["status"] == "success"
    assert full["n_shared_cells"] == 80
    assert full["factor_matching"]["mean_matched_abs_r"] == 1.0
    assert full["variance_explained_concordance"]["max_abs_diff"] == 0.0


def test_deterministic_stratified_subsampling_retains_rare_strata():
    import pandas as pd

    # 500 cells total: pop_major (492 cells) + pop_rare (8 cells = 1.6%)
    barcodes = [f"CELL_{i:04d}" for i in range(500)]
    labels = ["pop_major"] * 492 + ["pop_rare"] * 8
    obs_df = pd.DataFrame({"leiden_rna": labels}, index=barcodes)

    s1 = _deterministic_subsample(
        barcodes,
        [obs_df],
        max_cells=100,
        min_stratum_cells=5,
        explicit_strata_col=None,
        seed=0,
    )
    # Permute input row order and verify identical barcode set selected!
    perm = np.random.default_rng(99).permutation(500)
    bcs_perm = [barcodes[i] for i in perm]
    obs_perm = obs_df.iloc[perm].copy()
    s2 = _deterministic_subsample(
        bcs_perm,
        [obs_perm],
        max_cells=100,
        min_stratum_cells=5,
        explicit_strata_col=None,
        seed=0,
    )

    assert s1["subsampled"] is True
    assert s1["method"] == "stratified"
    assert s1["strata_col"] == "leiden_rna"
    assert len(s1["train_indices"]) == 100

    picked_bcs_1 = {barcodes[i] for i in s1["train_indices"]}
    picked_bcs_2 = {bcs_perm[i] for i in s2["train_indices"]}
    assert picked_bcs_1 == picked_bcs_2

    n_rare_picked = sum(1 for i in s1["train_indices"] if labels[i] == "pop_rare")
    assert n_rare_picked >= 5  # rare stratum floor enforced!
    # Rare cells have smaller inclusion weight than downsampled major cells
    w = s1["inclusion_weights"]
    assert w[0] > w[-1]


def test_no_hardcoded_biology_contract_and_agent_conflict_wiring():
    import aria.scripts.integration_mofa as mofa_mod

    assert not hasattr(mofa_mod, "CELL_CYCLE_GENES")
    mock_res = _mock_mofa(n_factors=8, reason="unit_test")
    assert mock_res["technical_factor_flag"] is False
    assert mock_res["technical_factor_id"] is None

    try:
        from aria.agents.integration_agent import IntegrationAgent
        from aria.utils.script_contracts import contract_for_script
        from tests.test_no_fabrication_guard import scan_source
    except ImportError:
        return  # Scientific-only env without pytest/litellm/pydantic; tested in aria-env.

    # No-fabrication static AST guard on new/updated V48 I4 modules
    for rel in (
        "aria/scripts/integration_mofa.py",
        "aria/benchmarks/concordance_mofa.py",
        "aria/scripts/benchmark_mofa_muon_ref.py",
    ):
        src = (ROOT / rel).read_text(encoding="utf-8")
        hits = scan_source(rel, src)
        assert hits == [], f"Fabrication guard hits in {rel}: {hits}"

    c = contract_for_script("aria/scripts/integration_mofa.py")
    assert c is not None
    assert c.validation_level == "production"
    assert c.validate_params({"h5mu_path": "/tmp/x.h5mu", "n_factors": 10, "seed": 0}) == []

    # Verify IntegrationAgent._resolve_conflicts surfaces excluded MOFA confounders
    agent = IntegrationAgent.__new__(IntegrationAgent)
    conflicts = agent._resolve_conflicts(
        "exp_i4",
        {"genome": "hg38"},
        {"user_question": "multiome factors"},
        {
            "mofa": {
                "status": "success",
                "technical_factor_flag": True,
                "technical_factor_id": 1,
                "excluded_factors": [
                    {"factor_id": 1, "exclusion_reasons": ["sequencing_depth_correlated"]}
                ],
            }
        },
    )
    assert any(
        x["type"] == "mofa_confounder_factor_excluded" and 1 in x["excluded_factor_ids"]
        for x in conflicts["conflicts"]
    )


# ── Full scientific runtime tests (executed when mofapy2 + mudata installed) ──


def test_end_to_end_mofa_planted_shared_specific_depth_and_cell_cycle(tmp_path: Path):
    if not _has_mofa_stack():
        return

    import anndata as ad
    import mudata as md
    import scipy.sparse as sp

    rng = np.random.default_rng(2026)
    n_cells = 180
    n_genes = 160
    n_peaks = 220
    barcodes = [f"BC_{i:04d}-1" for i in range(n_cells)]

    # Planted orthogonal latent drivers:
    # 1. z_shared: drives genes 0..39 AND peaks 0..49
    # 2. z_atac_only: drives peaks 50..99 ONLY
    # 3. z_cc: cell-cycle S/G2M driver on genes 40..69 (and recorded in obs S_score/G2M_score)
    # 4. depth_mult: 6x variation in global sequencing depth across cells!
    t = np.linspace(0, 2.0 * np.pi, n_cells, endpoint=False)
    z_shared = np.sin(t)
    z_atac_only = np.cos(2.0 * t)
    z_cc = np.sin(3.0 * t)
    depth_mult = np.exp(np.linspace(-0.9, 0.9, n_cells))

    base_rna = np.full((n_cells, n_genes), 1.5, dtype=np.float64)
    base_rna[:, 0:20] *= np.exp(1.4 * z_shared[:, None])
    base_rna[:, 20:40] *= np.exp(-1.3 * z_shared[:, None])
    base_rna[:, 40:55] *= np.exp(1.3 * z_cc[:, None])
    base_rna[:, 55:70] *= np.exp(-1.2 * z_cc[:, None])
    rna_counts = rng.poisson(base_rna * depth_mult[:, None]).astype(np.int32)

    base_atac = np.full((n_cells, n_peaks), 0.8, dtype=np.float64)
    base_atac[:, 0:25] *= np.exp(1.5 * z_shared[:, None])
    base_atac[:, 25:50] *= np.exp(-1.4 * z_shared[:, None])
    base_atac[:, 50:75] *= np.exp(1.5 * z_atac_only[:, None])
    base_atac[:, 75:100] *= np.exp(-1.4 * z_atac_only[:, None])
    atac_counts = rng.poisson(base_atac * depth_mult[:, None]).astype(np.int32)

    strata = ["rare_pop" if i < 8 else f"major_{i % 3}" for i in range(n_cells)]
    rna_ad = ad.AnnData(
        X=sp.csr_matrix(rna_counts),
        obs={
            "leiden_rna": strata,
            "S_score": z_cc + rng.normal(0.0, 0.05, size=n_cells),
            "G2M_score": 0.8 * z_cc + rng.normal(0.0, 0.05, size=n_cells),
        },
    )
    rna_ad.obs_names = barcodes
    rna_ad.var_names = [f"GENE_{j:03d}" for j in range(n_genes)]

    atac_ad = ad.AnnData(X=sp.csr_matrix(atac_counts))
    atac_ad.obs_names = barcodes
    atac_ad.var_names = [f"chr1:{j * 1000 + 100}-{j * 1000 + 600}" for j in range(n_peaks)]

    in_h5mu = tmp_path / "planted_paired.h5mu"
    md.MuData({"rna": rna_ad, "atac": atac_ad}).write_h5mu(str(in_h5mu))

    # Run 1 & Run 2 with seed=0 -> verify exact determinism + planted factor recovery
    res1 = integration_mofa(
        {
            "h5mu_path": str(in_h5mu),
            "output_dir": str(tmp_path / "out_r1"),
            "n_factors": 6,
            "n_top_genes": 160,
            "n_top_peaks": 220,
            "seed": 0,
        }
    )
    res2 = integration_mofa(
        {
            "h5mu_path": str(in_h5mu),
            "output_dir": str(tmp_path / "out_r2"),
            "n_factors": 6,
            "n_top_genes": 160,
            "n_top_peaks": 220,
            "seed": 0,
        }
    )
    assert res1["status"] == "success", res1
    assert res2["status"] == "success", res2
    assert res1["variance_explained"] == res2["variance_explained"]
    assert res1["n_cells"] == n_cells
    assert len(res1["shared_factors"]) >= 1
    assert res1["cell_cycle_evaluated"] is True
    assert res1["cell_cycle_factor_flag"] is True
    assert len(res1["cell_cycle_factor_ids"]) >= 1
    assert len(res1["excluded_factors"]) >= 1

    # Check exact score matrix equality across seed=0 runs
    m1 = md.read_h5mu(res1["h5mu_output_path"])
    m2 = md.read_h5mu(res2["h5mu_output_path"])
    z1 = np.asarray(m1.obsm["X_mofa"])
    z2 = np.asarray(m2.obsm["X_mofa"])
    assert np.max(np.abs(z1 - z2)) == 0.0
    assert np.all(np.isfinite(z1))
    assert "aria_mofa" in m1.uns

    # Verify planted shared factor and ATAC-specific factor are strongly recovered
    shared_corrs = [
        abs(r)
        for k in range(z1.shape[1])
        if (r := pearson_correlation(z1[:, k], z_shared)) is not None
    ]
    atac_corrs = [
        abs(r)
        for k in range(z1.shape[1])
        if (r := pearson_correlation(z1[:, k], z_atac_only)) is not None
    ]
    max_corr_shared = max(shared_corrs)
    max_corr_atac = max(atac_corrs)
    assert max_corr_shared > 0.85, f"Expected shared factor >0.85, got {max_corr_shared}"
    assert max_corr_atac > 0.80, f"Expected ATAC-only factor >0.80, got {max_corr_atac}"

    # Run 3: Subsampled (`max_cells=110 < 180`) + exact closed-form posterior projection
    res_sub = integration_mofa(
        {
            "h5mu_path": str(in_h5mu),
            "output_dir": str(tmp_path / "out_sub"),
            "n_factors": 6,
            "max_cells": 110,
            "n_top_genes": 160,
            "n_top_peaks": 220,
            "seed": 0,
        }
    )
    assert res_sub["status"] == "success", res_sub
    assert res_sub["subsampled"] is True
    assert res_sub["n_cells"] == n_cells
    assert res_sub["n_cells_trained"] == 110
    assert res_sub["subsampling"]["projection_min_train_abs_r"] >= 0.99

    m_sub = md.read_h5mu(res_sub["h5mu_output_path"])
    z_sub = np.asarray(m_sub.obsm["X_mofa"])
    assert z_sub.shape == (n_cells, res_sub["n_factors"])
    assert np.all(np.isfinite(z_sub))
    # Verify rare_pop retained >= 5 training cells!
    train_mask = np.asarray(m_sub.obs["mofa_training_subset"].values, dtype=bool)
    assert int(np.sum(train_mask[:8])) >= 5


def test_caller_supplied_cell_cycle_genes_and_missing_honesty(tmp_path: Path):
    if not _has_mofa_stack():
        return

    import anndata as ad
    import mudata as md
    import scipy.sparse as sp

    rng = np.random.default_rng(77)
    n_cells = 130
    barcodes = [f"C_{i:03d}" for i in range(n_cells)]
    cyc = np.sin(np.linspace(0, 4 * np.pi, n_cells))
    shared = np.cos(np.linspace(0, 2 * np.pi, n_cells))

    lam_r = np.full((n_cells, 100), 2.0)
    lam_r[:, 0:15] *= np.exp(1.5 * cyc[:, None])
    lam_r[:, 15:30] *= np.exp(1.4 * shared[:, None])
    rna_ad = ad.AnnData(X=sp.csr_matrix(rng.poisson(lam_r).astype(np.int32)))
    rna_ad.obs_names = barcodes
    rna_ad.var_names = [f"G_{j:03d}" for j in range(100)]

    lam_a = np.full((n_cells, 120), 1.0)
    lam_a[:, 0:25] *= np.exp(1.5 * shared[:, None])
    atac_ad = ad.AnnData(X=sp.csr_matrix(rng.poisson(lam_a).astype(np.int32)))
    atac_ad.obs_names = barcodes
    atac_ad.var_names = [f"chr1:{j*500+100}-{j*500+400}" for j in range(120)]

    p_h5mu = tmp_path / "no_obs_cc.h5mu"
    md.MuData({"rna": rna_ad, "atac": atac_ad}).write_h5mu(str(p_h5mu))

    # Without obs S_score/G2M_score or caller genes -> honestly recorded in metrics_not_computed
    r_none = integration_mofa(
        {
            "h5mu_path": str(p_h5mu),
            "output_dir": str(tmp_path / "cc_none"),
            "n_factors": 4,
            "seed": 0,
        }
    )
    assert r_none["status"] == "success"
    assert r_none["cell_cycle_evaluated"] is False
    assert "cell_cycle_factor_correlation" in r_none["metrics_not_computed"]

    # With caller-supplied phase gene sets -> dynamically scored and flagged!
    r_scored = integration_mofa(
        {
            "h5mu_path": str(p_h5mu),
            "output_dir": str(tmp_path / "cc_scored"),
            "n_factors": 4,
            "seed": 0,
            "cell_cycle_genes": {
                "s_genes": [f"G_{j:03d}" for j in range(0, 8)],
                "g2m_genes": [f"G_{j:03d}" for j in range(8, 15)],
            },
        }
    )
    assert r_scored["status"] == "success"
    assert r_scored["cell_cycle_evaluated"] is True
    assert r_scored["cell_cycle_source"] == "caller_cell_cycle_genes"
    assert r_scored["cell_cycle_factor_flag"] is True


if __name__ == "__main__":
    test_fns = [
        obj
        for name, obj in list(globals().items())
        if name.startswith("test_") and callable(obj)
    ]
    n_pass = 0
    for fn in test_fns:
        sig = inspect.signature(fn)
        if "tmp_path" in sig.parameters:
            with tempfile.TemporaryDirectory(prefix="test_v48_i4_") as td:
                fn(Path(td))
        else:
            fn()
        n_pass += 1
        print(f"  PASS {fn.__name__}")
    print(f"\nAll {n_pass}/{len(test_fns)} tests passed (mofa_stack={_has_mofa_stack()}).")
