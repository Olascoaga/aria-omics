"""
V48 Slice I3 — WNN Anti-Degradation Guardrail ("Lens, Not Replacement") Tests
-----------------------------------------------------------------------------
Deterministic unit and end-to-end tests for:
  - ``aria/utils/wnn_diagnostics.py`` (pure NumPy/SciPy 5-axis engine)
  - ``aria/scripts/integration_wnn_diagnostics.py`` (IPC script on chromatin stack)
  - ``aria/scripts/integration_wnn.py`` inline guardrail + HDF5 ``.uns`` persistence
  - ``aria/agents/integration_agent.py`` confidence downgrade & conflict wiring

All synthetic fixtures use strictly neutral identifiers (``pop_0``, ``rare_0``,
``gene_0``, ``chr1:100-600``) per ADR-011 and zero fabricated values per ADR-002.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.spatial import cKDTree

from aria.agents.integration_agent import IntegrationAgent
from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION
from aria.bus.message_bus import Confidence
from aria.scripts.integration_wnn_diagnostics import integration_wnn_diagnostics
from aria.utils.registry_integrity import (
    AGENT_SCRIPT_CONTRACTS,
    check_registry_integrity,
)
from aria.utils.script_contracts import SCRIPT_CONTRACTS
from aria.utils.wnn_diagnostics import (
    compute_wnn_anti_degradation_diagnostics,
    sanitize_for_hdf5_uns,
)
def _load_scan_source():
    guard_path = REPO_ROOT / "tests" / "test_no_fabrication_guard.py"
    spec = importlib.util.spec_from_file_location("_no_fab_guard", guard_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.scan_source


def _build_distance_csr_from_embedding(emb: np.ndarray, k: int = 12) -> sp.csr_matrix:
    """Build an exact off-diagonal Euclidean distance CSR graph."""
    n = emb.shape[0]
    normed = emb / np.maximum(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12)
    tree = cKDTree(normed)
    dists, idxs = tree.query(normed, k=min(n, k + 1))
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n):
        kept = 0
        for d, j in zip(dists[i], idxs[i]):
            if int(j) == i:
                continue
            rows.append(i)
            cols.append(int(j))
            vals.append(float(max(d, 1e-6)))
            kept += 1
            if kept >= k:
                break
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def _make_healthy_bimodal_fixture(
    seed: int = 42,
    n_per_major: int = 56,
    n_rare: int = 8,
    d_rna: int = 12,
    d_atac: int = 12,
    k: int = 12,
) -> dict[str, Any]:
    """Construct a well-separated 3-major + 1-rare bimodal manifold fixture.

    Each population occupies a multi-component direction across all dimensions
    so orthogonal dimension halves (D_build and D_eval) both carry clear signal.
    """
    rng = np.random.default_rng(seed)
    counts = [n_per_major, n_per_major, n_per_major, n_rare]
    names = ["pop_0", "pop_1", "pop_2", "rare_0"]

    # Orthogonal-ish multi-dimensional population signatures across all dims
    rna_centers = [
        np.array([2.5 if j % 4 == 0 else (-1.5 if j % 4 == 1 else 0.2) for j in range(d_rna)]),
        np.array([2.5 if j % 4 == 1 else (-1.5 if j % 4 == 2 else -0.2) for j in range(d_rna)]),
        np.array([2.5 if j % 4 == 2 else (-1.5 if j % 4 == 3 else 0.3) for j in range(d_rna)]),
        np.array([3.2 if j % 4 == 3 else (-2.0 if j % 4 == 0 else -0.4) for j in range(d_rna)]),
    ]
    atac_centers = [
        np.array([2.4 if j % 4 == 0 else (-1.4 if j % 4 == 2 else 0.2) for j in range(d_atac)]),
        np.array([2.4 if j % 4 == 1 else (-1.4 if j % 4 == 3 else -0.2) for j in range(d_atac)]),
        np.array([2.4 if j % 4 == 2 else (-1.4 if j % 4 == 0 else 0.3) for j in range(d_atac)]),
        np.array([3.0 if j % 4 == 3 else (-1.8 if j % 4 == 1 else -0.3) for j in range(d_atac)]),
    ]

    rna_blocks: list[np.ndarray] = []
    atac_blocks: list[np.ndarray] = []
    labels: list[str] = []
    for cnt, lbl, rc, ac in zip(counts, names, rna_centers, atac_centers):
        latent_cell = rng.normal(0.0, 0.20, size=(cnt, 1))
        rna_blocks.append(rc[None, :] + latent_cell + rng.normal(0.0, 0.22, size=(cnt, d_rna)))
        atac_blocks.append(ac[None, :] + latent_cell + rng.normal(0.0, 0.22, size=(cnt, d_atac)))
        labels.extend([lbl] * cnt)

    rna_emb = np.vstack(rna_blocks)
    atac_emb = np.vstack(atac_blocks)
    n_cells = rna_emb.shape[0]

    rna_l2 = rna_emb / np.linalg.norm(rna_emb, axis=1, keepdims=True)
    atac_l2 = atac_emb / np.linalg.norm(atac_emb, axis=1, keepdims=True)
    joint_concat = np.hstack([0.55 * rna_l2, 0.45 * atac_l2])
    wnn_csr = _build_distance_csr_from_embedding(joint_concat, k=k)

    rna_w = np.clip(0.55 + rng.normal(0.0, 0.04, size=n_cells), 0.35, 0.75)
    atac_w = 1.0 - rna_w

    return {
        "rna_embedding": rna_emb,
        "atac_embedding": atac_emb,
        "wnn_distances": wnn_csr,
        "wnn_clusters": list(labels),
        "rna_clusters": list(labels),
        "atac_clusters": list(labels),
        "rna_weights": rna_w,
        "atac_weights": atac_w,
        "k": k,
    }


# ── 1. Pure-Math 5-Axis Diagnostic Engine Tests ──────────────────────────────


def test_healthy_paired_fixture_passes_guardrail():
    fx = _make_healthy_bimodal_fixture(seed=11)
    diag = compute_wnn_anti_degradation_diagnostics(**fx, seed=0)

    assert diag["guardrail_status"] == "pass"
    assert diag["issue_codes"] == []
    assert diag["modality_preservation"]["mean_rna_recall_adjusted"] > 0.45
    assert diag["modality_preservation"]["mean_atac_recall_adjusted"] > 0.45

    ho = diag["held_out_check"]
    assert ho["evaluable"] is True
    assert ho["rna_manifold_informative"] is True
    assert ho["atac_manifold_informative"] is True
    assert ho["joint_bimodal_degraded"] is False
    assert ho["cross_modal_oof"]["rna_from_atac_excess_cosine"] > 0.25
    assert ho["cross_modal_oof"]["atac_from_rna_excess_cosine"] > 0.25

    rare_r = diag["rare_population_retention"]["rna_rare_clusters"]
    assert rare_r["n_rare_groups_detected"] == 1
    assert rare_r["n_retained"] == 1
    assert rare_r["n_degraded"] == 0
    assert diag["spurious_bridge_screen"]["n_flagged_clusters"] == 0


def test_pure_noise_atac_downweighted_warns_not_degraded():
    fx = _make_healthy_bimodal_fixture(seed=21)
    rng = np.random.default_rng(210)
    n_cells = fx["rna_embedding"].shape[0]
    noisy_atac = rng.normal(0.0, 1.0, size=fx["atac_embedding"].shape)
    wnn_from_rna = _build_distance_csr_from_embedding(fx["rna_embedding"], k=fx["k"])
    rna_w = np.full(n_cells, 0.92, dtype=np.float64)
    atac_w = np.full(n_cells, 0.08, dtype=np.float64)
    random_atac_clusters = [f"atac_noise_{i % 3}" for i in range(n_cells)]

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=noisy_atac,
        wnn_distances=wnn_from_rna,
        wnn_clusters=fx["wnn_clusters"],
        rna_clusters=fx["rna_clusters"],
        atac_clusters=random_atac_clusters,
        rna_weights=rna_w,
        atac_weights=atac_w,
        k=fx["k"],
        seed=0,
    )

    assert diag["held_out_check"]["atac_manifold_informative"] is False
    assert "atac_uninformative_downweighted" in diag["issue_codes"]
    assert "global_modality_imbalance" in diag["issue_codes"]
    assert "atac_structure_degraded" not in diag["issue_codes"]
    assert diag["guardrail_status"] == "warn"


def test_informative_atac_collapsed_flags_degraded():
    fx = _make_healthy_bimodal_fixture(seed=31)
    rng = np.random.default_rng(310)
    n_cells = fx["rna_embedding"].shape[0]

    # Build an independent ATAC manifold (shuffled cell order vs RNA) so RNA-only
    # graph completely ignores genuine ATAC local structure.
    perm = rng.permutation(n_cells)
    atac_indep = fx["atac_embedding"][perm]
    atac_clusters_indep = [fx["atac_clusters"][p] for p in perm]
    wnn_rna_only = _build_distance_csr_from_embedding(fx["rna_embedding"], k=fx["k"])

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=atac_indep,
        wnn_distances=wnn_rna_only,
        wnn_clusters=fx["wnn_clusters"],
        rna_clusters=fx["rna_clusters"],
        atac_clusters=atac_clusters_indep,
        rna_weights=np.full(n_cells, 0.55),
        atac_weights=np.full(n_cells, 0.45),
        k=fx["k"],
        seed=0,
    )

    assert diag["held_out_check"]["atac_manifold_informative"] is True
    assert "atac_structure_degraded" in diag["issue_codes"]
    assert diag["guardrail_status"] == "degraded"


def test_rare_population_scattered_flags_degraded():
    fx = _make_healthy_bimodal_fixture(seed=41, n_rare=8)
    n_cells = fx["rna_embedding"].shape[0]
    rare_indices = [
        i for i, lbl in enumerate(fx["rna_clusters"]) if lbl == "rare_0"
    ]
    major_indices = [
        i for i, lbl in enumerate(fx["rna_clusters"]) if lbl == "pop_0"
    ]

    # Corrupt WNN edges of rare_0 so all their joint neighbors point into pop_0
    lil = fx["wnn_distances"].tolil()
    for pos, r_idx in enumerate(rare_indices):
        lil.rows[r_idx] = [major_indices[(pos * 3 + m) % len(major_indices)] for m in range(fx["k"])]
        lil.data[r_idx] = [0.25 + 0.01 * m for m in range(fx["k"])]
    wnn_corrupted = lil.tocsr()

    wnn_labels_absorbed = [
        "pop_0" if lbl == "rare_0" else lbl for lbl in fx["wnn_clusters"]
    ]

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=fx["atac_embedding"],
        wnn_distances=wnn_corrupted,
        wnn_clusters=wnn_labels_absorbed,
        rna_clusters=fx["rna_clusters"],
        atac_clusters=fx["atac_clusters"],
        rna_weights=fx["rna_weights"],
        atac_weights=fx["atac_weights"],
        k=fx["k"],
        seed=0,
    )

    assert "rare_population_degraded" in diag["issue_codes"]
    assert diag["guardrail_status"] == "degraded"
    rare_grp = diag["rare_population_retention"]["rna_rare_clusters"]["groups"][0]
    assert rare_grp["status"] == "degraded"
    assert rare_grp["wnn_cohesion_raw"] == 0.0


def test_micro_rare_cluster_below_min_evaluable_is_not_false_alarm():
    fx = _make_healthy_bimodal_fixture(seed=51, n_rare=3)
    rare_indices = [i for i, lbl in enumerate(fx["rna_clusters"]) if lbl == "rare_0"]
    major_indices = [i for i, lbl in enumerate(fx["rna_clusters"]) if lbl == "pop_0"]

    lil = fx["wnn_distances"].tolil()
    for r_idx in rare_indices:
        lil.rows[r_idx] = major_indices[: fx["k"]]
        lil.data[r_idx] = [0.3] * fx["k"]

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=fx["atac_embedding"],
        wnn_distances=lil.tocsr(),
        wnn_clusters=fx["wnn_clusters"],
        rna_clusters=fx["rna_clusters"],
        atac_clusters=fx["atac_clusters"],
        rna_weights=fx["rna_weights"],
        atac_weights=fx["atac_weights"],
        k=fx["k"],
        min_evaluable_rare_cells=5,
        seed=0,
    )

    grp = diag["rare_population_retention"]["rna_rare_clusters"]["groups"][0]
    assert grp["n_cells"] == 3
    assert grp["evaluable"] is False
    assert grp["status"] == "not_evaluable"
    assert grp["reason"] == "below_min_evaluable_rare_cells"
    assert "rare_population_degraded" not in diag["issue_codes"]


def test_reference_audit_labels_rare_group_degraded():
    fx = _make_healthy_bimodal_fixture(seed=61, n_rare=8)
    audit_labels = list(fx["rna_clusters"])
    # Merge rare_0 into pop_0 in single-modality Leiden, leaving rare_0 visible
    # only via external audit_labels, then scatter its WNN neighborhood.
    merged_single = ["pop_0" if x == "rare_0" else x for x in fx["rna_clusters"]]
    rare_indices = [i for i, lbl in enumerate(audit_labels) if lbl == "rare_0"]
    major_indices = [i for i, lbl in enumerate(audit_labels) if lbl == "pop_0"]

    lil = fx["wnn_distances"].tolil()
    for r_idx in rare_indices:
        lil.rows[r_idx] = major_indices[: fx["k"]]
        lil.data[r_idx] = [0.3] * fx["k"]

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=fx["atac_embedding"],
        wnn_distances=lil.tocsr(),
        wnn_clusters=merged_single,
        rna_clusters=merged_single,
        atac_clusters=merged_single,
        rna_weights=fx["rna_weights"],
        atac_weights=fx["atac_weights"],
        audit_labels=audit_labels,
        k=fx["k"],
        seed=0,
    )

    assert "reference_rare_population_degraded" in diag["issue_codes"]
    assert diag["guardrail_status"] == "degraded"


def test_spurious_bridge_cluster_detected_vs_continuous_transition():
    fx = _make_healthy_bimodal_fixture(seed=71)
    # Create an artificial joint cluster 'bridge_art' mixing 10 cells of pop_0
    # and 10 cells of pop_1 that have zero native kNN edges between them.
    wnn_bridged = list(fx["wnn_clusters"])
    pop0_idx = [i for i, x in enumerate(wnn_bridged) if x == "pop_0"][:10]
    pop1_idx = [i for i, x in enumerate(wnn_bridged) if x == "pop_1"][:10]
    for idx in [*pop0_idx, *pop1_idx]:
        wnn_bridged[idx] = "bridge_art"

    diag_bad = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=fx["atac_embedding"],
        wnn_distances=fx["wnn_distances"],
        wnn_clusters=wnn_bridged,
        rna_clusters=fx["rna_clusters"],
        atac_clusters=fx["atac_clusters"],
        rna_weights=fx["rna_weights"],
        atac_weights=fx["atac_weights"],
        k=fx["k"],
        seed=0,
    )
    assert "bridge_art" in diag_bad["spurious_bridge_screen"]["flagged_clusters"]
    assert "spurious_bridge_candidate" in diag_bad["issue_codes"]
    assert diag_bad["guardrail_status"] == "degraded"

    # Now construct a genuine continuous transition along a 1D manifold where
    # boundary cells naturally share native kNN edges across pop_a and pop_b.
    rng = np.random.default_rng(72)
    t = np.linspace(-1.0, 1.0, 120)
    d = 12
    u_vec = np.array([1.0 if j % 2 == 0 else -0.5 for j in range(d)])
    v_vec = np.array([0.4 if j % 3 == 0 else 0.9 for j in range(d)])
    cont_rna = np.outer(np.cos(t), u_vec) + np.outer(np.sin(t), v_vec) + rng.normal(0.0, 0.08, (120, d))
    cont_atac = np.outer(np.cos(t), u_vec) + np.outer(np.sin(t), v_vec) + rng.normal(0.0, 0.08, (120, d))
    cont_csr = _build_distance_csr_from_embedding(np.hstack([cont_rna, cont_atac]), k=12)

    single_lbls = ["seg_left" if x < 0.0 else "seg_right" for x in t]
    joint_lbls = [
        "seg_mid" if abs(x) < 0.25 else ("seg_left" if x < 0.0 else "seg_right")
        for x in t
    ]
    diag_cont = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=cont_rna,
        atac_embedding=cont_atac,
        wnn_distances=cont_csr,
        wnn_clusters=joint_lbls,
        rna_clusters=single_lbls,
        atac_clusters=single_lbls,
        rna_weights=np.full(120, 0.5),
        atac_weights=np.full(120, 0.5),
        k=12,
        seed=0,
    )
    assert "seg_mid" not in diag_cont["spurious_bridge_screen"]["flagged_clusters"]
    assert "spurious_bridge_candidate" not in diag_cont["issue_codes"]


def test_orphan_wnn_edges_and_held_out_corruption_detected():
    fx = _make_healthy_bimodal_fixture(seed=81)
    rng = np.random.default_rng(810)
    n_cells = fx["rna_embedding"].shape[0]

    # Randomize WNN graph completely -> orphan edges + degraded held-out graph cosine
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n_cells):
        choices = rng.choice([j for j in range(n_cells) if j != i], size=fx["k"], replace=False)
        for c in choices:
            rows.append(i)
            cols.append(int(c))
            vals.append(0.5)
    random_csr = sp.csr_matrix((vals, (rows, cols)), shape=(n_cells, n_cells))

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=fx["atac_embedding"],
        wnn_distances=random_csr,
        wnn_clusters=fx["wnn_clusters"],
        rna_clusters=fx["rna_clusters"],
        atac_clusters=fx["atac_clusters"],
        rna_weights=fx["rna_weights"],
        atac_weights=fx["atac_weights"],
        k=fx["k"],
        seed=0,
    )
    assert diag["guardrail_status"] == "degraded"
    assert "spurious_bridge_candidate" in diag["issue_codes"]
    assert "held_out_joint_degraded" in diag["issue_codes"]


def test_seed_determinism_zero_norm_and_hdf5_uns_roundtrip(tmp_path: Path):
    import anndata as ad

    fx = _make_healthy_bimodal_fixture(seed=91)
    fx["rna_embedding"][0, :] = 0.0  # exercise zero-norm cell guard

    d1 = compute_wnn_anti_degradation_diagnostics(**fx, seed=13)
    d2 = compute_wnn_anti_degradation_diagnostics(**fx, seed=13)
    assert json.dumps(d1, sort_keys=True) == json.dumps(d2, sort_keys=True)
    assert d1["modality_preservation"]["n_zero_norm_rna_cells"] == 1

    # Also induce issues (list of dicts) so sanitize_for_hdf5_uns exercises list[dict]
    d_no_weights = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=fx["rna_embedding"],
        atac_embedding=fx["atac_embedding"],
        wnn_distances=fx["wnn_distances"],
        wnn_clusters=fx["wnn_clusters"],
        rna_clusters=fx["rna_clusters"],
        atac_clusters=fx["atac_clusters"],
        rna_weights=None,
        atac_weights=None,
        k=fx["k"],
        seed=13,
    )
    assert "missing_modality_weights" in d_no_weights["issue_codes"]

    clean_uns = sanitize_for_hdf5_uns(d_no_weights)
    adata = ad.AnnData(X=np.zeros((6, 4), dtype=np.float32))
    adata.uns["aria_wnn_diagnostics"] = clean_uns
    h5ad_file = tmp_path / "uns_roundtrip.h5ad"
    adata.write_h5ad(h5ad_file)
    reloaded = ad.read_h5ad(h5ad_file)
    assert reloaded.uns["aria_wnn_diagnostics"]["guardrail_status"] == d_no_weights["guardrail_status"]

    with pytest.raises(ValueError, match="At least 6 shared cells"):
        compute_wnn_anti_degradation_diagnostics(
            rna_embedding=np.ones((4, 6)),
            atac_embedding=np.ones((4, 6)),
            wnn_distances=sp.eye(4, format="csr"),
            wnn_clusters=["0"] * 4,
            rna_clusters=["0"] * 4,
            atac_clusters=["0"] * 4,
        )


# ── 2. AST No-Fabrication, Contract & Agent Wiring Tests ─────────────────────


def test_no_fabrication_ast_scan_on_i3_modules():
    scan_source = _load_scan_source()
    for rel in (
        "aria/utils/wnn_diagnostics.py",
        "aria/scripts/integration_wnn_diagnostics.py",
    ):
        src = (REPO_ROOT / rel).read_text(encoding="utf-8")
        violations = scan_source(src, rel)
        assert violations == [], f"Unexpected no-fabrication violations in {rel}: {violations}"


def test_script_contract_and_registry_integrity_for_i3():
    key = "aria/scripts/integration_wnn_diagnostics.py"
    assert key in SCRIPT_CONTRACTS
    assert SCRIPT_CONTRACTS[key].validation_level == "beta"
    assert key in AGENT_SCRIPT_CONTRACTS["integration_agent"]["required"]
    assert check_registry_integrity() == []
    assert INTEGRATION_VALIDATION["level"] == "beta"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is True
    assert INTEGRATION_VALIDATION["production_promoted"] is False


def test_integration_wnn_diagnostics_missing_file_fails_closed(tmp_path: Path):
    res = integration_wnn_diagnostics({"h5mu_path": str(tmp_path / "missing.h5mu")})
    assert res["status"] == "error"
    assert res["error_type"] == "MissingH5MU"


def test_integration_agent_guardrail_confidence_and_conflict_wiring():
    agent = IntegrationAgent.__new__(IntegrationAgent)
    captured_findings: list[dict[str, Any]] = []
    captured_dispatch: list[tuple[str, str, dict]] = []

    agent.publish_finding = lambda exp_id, payload, conf: captured_findings.append(
        {"payload": payload, "conf": conf}
    )

    class _StubEnv:
        def run_in_stack(self, stack, script_path, params, timeout=None):
            captured_dispatch.append((stack, script_path, params))
            return {"status": "success", "guardrail_status": "pass"}

    agent.env = _StubEnv()

    res = agent._run_wnn_diagnostics("/tmp/sample.h5mu", seed=7, audit_group_col="grp")
    assert res["status"] == "success"
    assert captured_dispatch[0][0] == "chromatin"
    assert captured_dispatch[0][1] == "aria/scripts/integration_wnn_diagnostics.py"
    assert captured_dispatch[0][2]["audit_group_col"] == "grp"

    # Guardrail 'warn' downgrades HIGH -> MEDIUM
    agent._evaluate_wnn_weights(
        "exp_w",
        {
            "status": "success",
            "mean_rna_weight": 0.58,
            "mean_atac_weight": 0.42,
            "anti_degradation": {
                "guardrail_status": "warn",
                "issue_codes": ["cluster_modality_collapse"],
            },
        },
        {},
    )
    assert captured_findings[-1]["conf"] == Confidence.MEDIUM
    assert captured_findings[-1]["payload"]["guardrail_status"] == "warn"

    # Guardrail 'degraded' downgrades HIGH -> LOW and surfaces conflicts
    degraded_wnn = {
        "status": "success",
        "mean_rna_weight": 0.58,
        "mean_atac_weight": 0.42,
        "n_discordant_clusters": 0,
        "anti_degradation": {
            "guardrail_status": "degraded",
            "issue_codes": [
                "rare_population_degraded",
                "spurious_bridge_candidate",
                "atac_structure_degraded",
            ],
            "spurious_bridge_screen": {"flagged_clusters": ["3"]},
        },
    }
    agent._evaluate_wnn_weights("exp_d", degraded_wnn, {})
    assert captured_findings[-1]["conf"] == Confidence.LOW
    assert "degraded" in captured_findings[-1]["payload"]["summary"]

    conflicts = agent._resolve_conflicts("exp_d", {}, {}, {"wnn": degraded_wnn})
    c_types = {c["type"] for c in conflicts["conflicts"]}
    assert "wnn_rare_population_degradation" in c_types
    assert "wnn_spurious_bridge_candidate" in c_types
    assert "wnn_modality_structure_degradation" in c_types


# ── 3. End-to-End .h5mu Integration Tests (aria-chromatin-env) ───────────────


def test_end_to_end_wnn_persists_guardrail_and_matches_standalone(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    pytest.importorskip("muon")

    from aria.scripts.integration_wnn import integration_wnn

    rng = np.random.default_rng(123)
    n_cells = 96
    n_genes = 120
    n_peaks = 160

    barcodes = [f"CELL_{i:03d}-1" for i in range(n_cells)]
    rna_counts = np.zeros((n_cells, n_genes), dtype=np.int32)
    atac_counts = np.zeros((n_cells, n_peaks), dtype=np.int32)

    # 3 distinct populations (32 cells each) with coherent RNA and ATAC blocks
    for g_idx in range(3):
        sl = slice(g_idx * 32, (g_idx + 1) * 32)
        g_genes = slice(g_idx * 40, (g_idx + 1) * 40)
        g_peaks = slice(g_idx * 50, (g_idx + 1) * 50)
        rna_counts[sl, :] = rng.poisson(1, size=(32, n_genes))
        rna_counts[sl, g_genes] += rng.poisson(12, size=(32, 40))
        atac_counts[sl, :] = rng.poisson(1, size=(32, n_peaks))
        atac_counts[sl, g_peaks] += rng.poisson(8, size=(32, 50))

    rna_ad = ad.AnnData(
        X=sp.csr_matrix(rna_counts),
        obs={"grp": [f"pop_{i // 32}" for i in range(n_cells)]},
    )
    rna_ad.obs_names = barcodes
    rna_ad.var_names = [f"gene_{j}" for j in range(n_genes)]

    atac_ad = ad.AnnData(X=sp.csr_matrix(atac_counts))
    atac_ad.obs_names = barcodes
    atac_ad.var_names = [f"chr1:{j * 1000 + 100}-{j * 1000 + 600}" for j in range(n_peaks)]

    mdata = md.MuData({"rna": rna_ad, "atac": atac_ad})
    in_h5mu = tmp_path / "paired_in.h5mu"
    mdata.write_h5mu(str(in_h5mu))

    wnn_out_dir = tmp_path / "wnn_out"
    wnn_res = integration_wnn(
        {
            "h5mu_path": str(in_h5mu),
            "output_dir": str(wnn_out_dir),
            "k": 12,
            "n_pcs": 12,
            "n_lsi": 12,
            "seed": 0,
        }
    )
    assert wnn_res["status"] == "success", wnn_res
    assert wnn_res["guardrail_status"] in ("pass", "warn")
    assert isinstance(wnn_res["anti_degradation"], dict)

    out_h5mu = Path(wnn_res["output_path"])
    saved = md.read_h5mu(str(out_h5mu))
    assert "aria_wnn_diagnostics" in saved.uns

    diag_json = tmp_path / "diag_report.json"
    standalone = integration_wnn_diagnostics(
        {
            "h5mu_path": str(out_h5mu),
            "k": 12,
            "seed": 0,
            "audit_group_col": "grp",
            "update_h5mu": True,
            "output_json": str(diag_json),
        }
    )
    assert standalone["status"] == "success", standalone
    assert diag_json.is_file()
    assert standalone["diagnostics"]["modality_preservation"] == wnn_res["anti_degradation"]["modality_preservation"]
    assert standalone["diagnostics"]["held_out_check"] == wnn_res["anti_degradation"]["held_out_check"]
