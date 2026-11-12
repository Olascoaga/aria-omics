"""V48 Slice I2 — validated deterministic WNN (paired RNA+ATAC joint state).

Group A (any env, no muon/mudata required):
  * pure NumPy concordance scorer (`aria/benchmarks/concordance_wnn.py`);
  * IPC contract level, orchestrator dispatch guardrail, IntegrationAgent WNN
    dispatch wiring and honest null-weight handling;
  * structured input errors raised before any heavy import.

Group B (`pytest.importorskip("muon")`, runs in `aria-chromatin-env`):
  * real WNN on a small synthetic paired `.h5mu` (3 separated populations);
  * strict seed determinism, `.h5mu` round trip + re-entrance, raw-count
    rejection, measured LSI depth-correlation filter;
  * cross-stack IPC E2E from `aria-env` into `aria-chromatin-env`.

Neutral labels only (ADR-011); honest `not_run` / null semantics (ADR-002).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aria.benchmarks.concordance_wnn import (
    knn_jaccard_overlap,
    modality_weight_concordance,
    pearson_correlation,
    score_seed_stability,
    score_wnn_concordance,
    spearman_correlation,
)
from aria.scripts import integration_wnn as wnn

WNN_SCRIPT = "aria/scripts/integration_wnn.py"
CHROMATIN_PY = Path("/home/medusa/anaconda3/envs/aria-chromatin-env/bin/python")


# ═══════════════════════════════════════════════════════════════════════════
# Group A.1 — concordance_wnn (pure NumPy)
# ═══════════════════════════════════════════════════════════════════════════

def test_pearson_perfect_and_negative():
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert pearson_correlation(a, [2.0, 4.0, 6.0, 8.0, 10.0]) == pytest.approx(1.0)
    assert pearson_correlation(a, [5.0, 4.0, 3.0, 2.0, 1.0]) == pytest.approx(-1.0)


def test_pearson_undefined_cases_return_none():
    assert pearson_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None
    assert pearson_correlation([1.0], [2.0]) is None
    assert pearson_correlation([1.0, 2.0], [1.0, 2.0, 3.0]) is None
    assert pearson_correlation([np.nan, 1.0, np.nan], [1.0, 2.0, 3.0]) is None


def test_pearson_filters_nans():
    a = [1.0, 2.0, np.nan, 4.0, 5.0]
    b = [2.0, 4.0, 100.0, 8.0, np.inf]
    # Only (1,2), (2,4), (4,8) remain -> perfect linear relation.
    assert pearson_correlation(a, b) == pytest.approx(1.0)


def test_spearman_monotonic_nonlinear_and_ties():
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert spearman_correlation(a, [1.0, 8.0, 27.0, 64.0, 125.0]) == pytest.approx(1.0)
    assert spearman_correlation(a, [10.0, 5.0, 1.0, 0.5, 0.1]) == pytest.approx(-1.0)

    # Average-tie ranks: x ranks [1,2.5,2.5,4], y ranks [1,2,3,4].
    x = [1.0, 2.0, 2.0, 3.0]
    y = [1.0, 2.0, 3.0, 4.0]
    expected = float(np.corrcoef([1.0, 2.5, 2.5, 4.0], [1.0, 2.0, 3.0, 4.0])[0, 1])
    assert spearman_correlation(x, y) == pytest.approx(expected)


def test_spearman_constant_and_nan_handling():
    assert spearman_correlation([3.0, 3.0, 3.0], [1.0, 2.0, 3.0]) is None
    assert spearman_correlation([1.0], [1.0]) is None
    assert spearman_correlation(
        [1.0, np.nan, 2.0, 3.0], [1.0, 5.0, 4.0, 9.0]
    ) == pytest.approx(1.0)


def test_modality_weight_concordance_perfect_agreement():
    w = [0.1, 0.4, 0.6, 0.9]
    out = modality_weight_concordance(w, w)
    assert out["n_cells_compared"] == 4
    assert out["pearson_r"] == pytest.approx(1.0)
    assert out["spearman_rho"] == pytest.approx(1.0)
    assert out["mae"] == 0.0
    assert out["median_abs_diff"] == 0.0
    assert out["dominant_modality_agreement"] == 1.0
    assert out["aria_mean_rna_weight"] == out["ref_mean_rna_weight"] == 0.5


def test_modality_weight_concordance_negative_and_mae():
    aria = [0.2, 0.4, 0.6, 0.8]
    ref = [0.8, 0.6, 0.4, 0.2]
    out = modality_weight_concordance(aria, ref)
    assert out["pearson_r"] == pytest.approx(-1.0)
    assert out["spearman_rho"] == pytest.approx(-1.0)
    # |diffs| = [0.6, 0.2, 0.2, 0.6]
    assert out["mae"] == pytest.approx(0.4)
    assert out["median_abs_diff"] == pytest.approx(0.4)
    # Dominant (w>=0.5): aria [F,F,T,T] vs ref [T,T,F,F] -> 0 agreement.
    assert out["dominant_modality_agreement"] == 0.0


def test_modality_weight_concordance_partial_dominance_and_errors():
    out = modality_weight_concordance([0.6, 0.6, 0.3, 0.3], [0.7, 0.4, 0.2, 0.55])
    # aria dom [T,T,F,F]; ref dom [T,F,F,T] -> 2/4 agree.
    assert out["dominant_modality_agreement"] == 0.5
    with pytest.raises(ValueError):
        modality_weight_concordance([0.1, 0.2], [0.1])
    empty = modality_weight_concordance([np.nan], [np.nan])
    assert empty["n_cells_compared"] == 0
    assert empty["pearson_r"] is None
    assert empty["mae"] is None


def test_knn_jaccard_excludes_self():
    aria = [[0, 1, 2], [1, 0, 2]]
    ref = [[1, 2], [0, 2]]
    out = knn_jaccard_overlap(aria, ref)
    assert out["n_cells_compared"] == 2
    assert out["mean_jaccard"] == 1.0
    assert out["mean_recall"] == 1.0


def test_knn_jaccard_match_effective_k():
    cell_ids = ["c0", "c1"]
    # ARIA returns k=3 non-self neighbors; reference returns k-1=2 (self in slot 1).
    aria = [["c0", "n1", "n2", "n3"], ["c1", "m1", "m2", "m3"]]
    ref = [["c0", "n1", "n2"], ["c1", "m1", "m2"]]

    matched = knn_jaccard_overlap(aria, ref, cell_ids=cell_ids, match_effective_k=True)
    assert matched["mean_jaccard"] == 1.0
    assert matched["mean_recall"] == 1.0

    unmatched = knn_jaccard_overlap(
        aria, ref, cell_ids=cell_ids, match_effective_k=False
    )
    assert unmatched["mean_jaccard"] == pytest.approx(round(2.0 / 3.0, 4))
    assert unmatched["mean_recall"] == 1.0


def test_knn_jaccard_errors_and_empty():
    with pytest.raises(ValueError):
        knn_jaccard_overlap([[1]], [[1], [2]])
    assert knn_jaccard_overlap([], [])["mean_jaccard"] is None
    assert knn_jaccard_overlap([[0]], [[0]])["n_cells_compared"] == 0


def test_seed_stability_identical_runs():
    labels = ["0", "0", "1", "1", "2", "2"]
    weights = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    out = score_seed_stability([labels, list(labels)], [weights, list(weights)])
    assert out["status"] == "success"
    assert out["n_runs"] == 2
    assert out["mean_pairwise_ari"] == 1.0
    assert out["min_pairwise_ari"] == 1.0
    assert out["max_weight_abs_diff"] == 0.0


def test_seed_stability_relabelled_partition_and_weight_drift():
    a = ["0", "0", "1", "1"]
    b = ["x", "x", "y", "y"]
    out = score_seed_stability([a, b], [[0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.75]])
    assert out["mean_pairwise_ari"] == 1.0
    assert out["max_weight_abs_diff"] == pytest.approx(0.25)


def test_seed_stability_not_run_with_single_run():
    out = score_seed_stability([["0", "1"]])
    assert out["status"] == "not_run"
    assert out["n_runs"] == 1
    assert out["mean_pairwise_ari"] is None
    assert out["max_weight_abs_diff"] is None


def test_score_wnn_concordance_not_run_without_reference():
    out = score_wnn_concordance(
        aria_barcodes=["a", "b"], aria_clusters=["0", "1"],
    )
    assert out["status"] == "not_run"
    assert "reason" in out
    assert "cluster_ari" not in out
    assert out["seed_stability"]["status"] == "not_run"

    empty_ref = score_wnn_concordance(
        aria_barcodes=["a", "b"], aria_clusters=["0", "1"],
        ref_barcodes=["a", "b"], ref_clusters=[],
    )
    assert empty_ref["status"] == "not_run"


def test_score_wnn_concordance_not_run_with_zero_shared_barcodes():
    out = score_wnn_concordance(
        aria_barcodes=["a", "b"], aria_clusters=["0", "1"],
        ref_barcodes=["x", "y"], ref_clusters=["0", "1"],
    )
    assert out["status"] == "not_run"
    assert out["n_shared_cells"] == 0
    assert out["n_aria_cells"] == 2
    assert out["n_ref_cells"] == 2
    assert "cluster_ari" not in out


def test_score_wnn_concordance_length_mismatch_raises():
    with pytest.raises(ValueError):
        score_wnn_concordance(
            aria_barcodes=["a", "b"], aria_clusters=["0"],
            ref_barcodes=["a"], ref_clusters=["0"],
        )


def _concordance_fixture():
    aria_bc = [f"cell{i}" for i in range(9)]
    aria_cl = ["0", "0", "0", "1", "1", "1", "2", "2", "2"]
    aria_w = [0.2, 0.3, 0.25, 0.6, 0.7, 0.65, 0.5, 0.45, 0.55]
    labels = ["group_a"] * 3 + ["group_b"] * 3 + ["group_c"] * 3
    perm = [4, 8, 0, 6, 2, 7, 1, 5, 3]
    ref_bc = [aria_bc[i] for i in perm]
    # Reference uses different label names but the same partition.
    ref_cl = [{"0": "A", "1": "B", "2": "C"}[aria_cl[i]] for i in perm]
    ref_w = [aria_w[i] for i in perm]
    return aria_bc, aria_cl, aria_w, labels, ref_bc, ref_cl, ref_w


def test_score_wnn_concordance_permuted_reference_alignment():
    aria_bc, aria_cl, aria_w, labels, ref_bc, ref_cl, ref_w = _concordance_fixture()
    aria_knn = {bc: [bc] + [x for x in aria_bc if x != bc][:2] for bc in aria_bc}
    ref_knn = {bc: aria_knn[bc][1:] for bc in aria_bc}

    out = score_wnn_concordance(
        aria_barcodes=aria_bc, aria_clusters=aria_cl, aria_rna_weights=aria_w,
        ref_barcodes=ref_bc, ref_clusters=ref_cl, ref_rna_weights=ref_w,
        aria_knn=aria_knn, ref_knn=ref_knn,
        seed_cluster_runs=[aria_cl, list(aria_cl)],
        seed_weight_runs=[aria_w, list(aria_w)],
    )
    assert out["status"] == "success"
    assert out["n_shared_cells"] == 9
    assert out["n_aria_clusters"] == out["n_ref_clusters"] == 3
    assert out["cluster_ari"] == 1.0
    assert out["cluster_nmi"] == pytest.approx(1.0)
    mw = out["modality_weights"]
    assert mw["status"] == "success"
    assert mw["pearson_r"] == pytest.approx(1.0)
    assert mw["mae"] == 0.0
    assert mw["dominant_modality_agreement"] == 1.0
    assert out["knn_overlap"]["status"] == "success"
    assert out["knn_overlap"]["mean_jaccard"] == 1.0
    assert out["seed_stability"]["mean_pairwise_ari"] == 1.0
    assert out["label_recovery"]["status"] == "not_run"

    # Deterministic: identical inputs give identical scores.
    again = score_wnn_concordance(
        aria_barcodes=aria_bc, aria_clusters=aria_cl, aria_rna_weights=aria_w,
        ref_barcodes=ref_bc, ref_clusters=ref_cl, ref_rna_weights=ref_w,
        aria_knn=aria_knn, ref_knn=ref_knn,
        seed_cluster_runs=[aria_cl, list(aria_cl)],
        seed_weight_runs=[aria_w, list(aria_w)],
    )
    assert again == out


def test_score_wnn_concordance_label_recovery_neutral_labels():
    aria_bc, aria_cl, aria_w, labels, ref_bc, ref_cl, ref_w = _concordance_fixture()
    rna_only = ["0", "0", "0", "1", "1", "1", "1", "1", "1"]  # merges b+c
    atac_only = list(aria_cl)
    labels_with_missing = list(labels)
    labels_with_missing[0] = None
    labels_with_missing[1] = "nan"

    out = score_wnn_concordance(
        aria_barcodes=aria_bc, aria_clusters=aria_cl,
        ref_barcodes=ref_bc, ref_clusters=ref_cl,
        aria_rna_only_clusters=rna_only,
        aria_atac_only_clusters=atac_only,
        reference_labels=labels_with_missing,
    )
    lr = out["label_recovery"]
    assert lr["status"] == "success"
    assert lr["n_labeled_cells"] == 7
    assert lr["n_unique_labels"] == 3
    assert lr["aria_wnn_ari"] == 1.0
    assert lr["ref_wnn_ari"] == 1.0
    assert lr["aria_atac_only_ari"] == 1.0
    assert lr["aria_rna_only_ari"] < 1.0


# ═══════════════════════════════════════════════════════════════════════════
# Group A.2 — contracts, governance, IntegrationAgent wiring, input errors
# ═══════════════════════════════════════════════════════════════════════════

def _require_orchestration_deps():
    pytest.importorskip("pydantic")
    pytest.importorskip("litellm")


def _wnn_success_payload(out_path: str) -> dict:
    return {
        "status": "success",
        "output_path": out_path,
        "n_cells": 75,
        "n_joint_clusters": 3,
        "n_rna_only_clusters": 3,
        "n_atac_only_clusters": 3,
        "n_discordant_clusters": 0,
        "per_cluster_modality_weights": {"0": {"n_cells": 25}},
        "warnings": [],
    }


def test_wnn_contract_is_beta_and_validates_output_path(tmp_path):
    pytest.importorskip("pydantic")
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    contract = SCRIPT_CONTRACTS[WNN_SCRIPT]
    assert contract.validation_level == "beta"

    h5mu = tmp_path / "paired.h5mu"
    h5mu.write_bytes(b"")
    assert contract.validate_params({"h5mu_path": str(h5mu)}) == []
    assert contract.validate_params({"data_path": str(h5mu)}) == []

    out = tmp_path / "wnn_integrated.h5mu"
    out.write_bytes(b"")
    assert contract.validate_result(_wnn_success_payload(str(out))) == []

    ghost = _wnn_success_payload(str(tmp_path / "missing.h5mu"))
    assert contract.validate_result(ghost) != []

    incomplete = {
        k: v for k, v in _wnn_success_payload(str(out)).items()
        if k != "n_discordant_clusters"
    }
    assert contract.validate_result(incomplete) != []


def test_integration_dispatch_guardrail_stays_scaffold():
    _require_orchestration_deps()
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION

    assert INTEGRATION_VALIDATION["level"] == "beta"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is True
    assert INTEGRATION_VALIDATION["production_promoted"] is False


class _RecordingEnv:
    def __init__(self, result: dict):
        self.calls: list[dict] = []
        self.result = result

    def run_in_stack(self, stack, script_path, params, timeout=None):
        self.calls.append({"stack": stack, "script_path": script_path,
                           "params": dict(params)})
        return dict(self.result)


class _Advisor:
    def advise_wnn_k(self, **kwargs):
        return SimpleNamespace(
            decision_id="wnn-k-test", chosen_value=20,
            chosen_by="advisor", approved_by_user=False,
        )


def _make_agent(env_result: dict):
    from aria.agents.integration_agent import IntegrationAgent

    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.env = _RecordingEnv(env_result)
    agent.advisor = _Advisor()
    agent._findings = []
    agent.publish_status = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None
    agent.publish_finding = (
        lambda exp_id, payload, conf: agent._findings.append((payload, conf))
    )
    agent.publish_blocking_escalation = lambda **kw: (None, {"choice": "accept"})
    agent.checkpoint_choice_text = lambda decision: "Use recommended k=20"
    agent.numeric_checkpoint_override = lambda decision, cast: None
    agent._format_wnn_checkpoint = lambda decision: "WNN k checkpoint"
    agent._get_output_dir = lambda files: "/tmp/aria_wnn_test_out"
    return agent


def test_run_wnn_dispatches_paired_h5mu_to_chromatin_stack():
    _require_orchestration_deps()
    agent = _make_agent({"status": "error", "error_type": "X", "details": "d"})
    h5mu = "/data/sample1_paired.h5mu"
    result = agent._run_wnn(
        "exp-1", {"seed": 7}, {}, {"scRNA": [h5mu], "scATAC": [h5mu]},
    )
    assert result["status"] == "error"
    assert len(agent.env.calls) == 1
    call = agent.env.calls[0]
    assert call["stack"] == "chromatin"
    assert call["script_path"] == WNN_SCRIPT
    assert call["params"]["h5mu_path"] == h5mu
    assert call["params"]["seed"] == 7
    assert call["params"]["k"] == 20


def test_run_wnn_default_seed_is_zero():
    _require_orchestration_deps()
    agent = _make_agent({"status": "error", "error_type": "X", "details": "d"})
    agent._run_wnn("exp-1", {}, {}, {"scRNA": ["/d/p.h5mu"], "scATAC": ["/d/p.h5mu"]})
    assert agent.env.calls[0]["params"]["seed"] == 0


def test_evaluate_wnn_weights_null_is_low_confidence_not_fabricated():
    _require_orchestration_deps()
    from aria.bus.message_bus import Confidence

    agent = _make_agent({})
    agent._evaluate_wnn_weights(
        "exp-1", {"status": "success", "mean_rna_weight": None,
                  "mean_atac_weight": None}, {},
    )
    assert len(agent._findings) == 1
    payload, conf = agent._findings[0]
    assert conf == Confidence.LOW
    assert payload["rna_weight"] is None
    assert payload["atac_weight"] is None
    text = str(payload)
    assert "0.6" not in text and "0.4" not in text


def test_error_missing_paired_h5mu(tmp_path):
    out = wnn.integration_wnn({
        "h5mu_path": str(tmp_path / "absent.h5mu"),
        "output_dir": str(tmp_path / "out"),
    })
    assert out["status"] == "error"
    assert out["error_type"] == "MissingPairedH5mu"


def test_error_missing_rna(tmp_path):
    atac = tmp_path / "atac.h5ad"
    atac.write_bytes(b"")
    out = wnn.integration_wnn({
        "rna_files": [str(tmp_path / "absent_rna.h5ad")],
        "atac_files": [str(atac)],
        "output_dir": str(tmp_path / "out"),
    })
    assert out["status"] == "error"
    assert out["error_type"] == "MissingRNA"


def test_error_missing_atac(tmp_path):
    rna = tmp_path / "rna.h5ad"
    rna.write_bytes(b"")
    out = wnn.integration_wnn({
        "rna_files": [str(rna)],
        "atac_files": [],
        "output_dir": str(tmp_path / "out"),
    })
    assert out["status"] == "error"
    assert out["error_type"] == "MissingATAC"


def test_error_mismatched_paired_h5mu(tmp_path):
    a = tmp_path / "a.h5mu"
    b = tmp_path / "b.h5mu"
    a.write_bytes(b"")
    b.write_bytes(b"")
    out = wnn.integration_wnn({
        "rna_files": [str(a)],
        "atac_files": [str(b)],
        "output_dir": str(tmp_path / "out"),
    })
    assert out["status"] == "error"
    assert out["error_type"] == "MismatchedPairedH5mu"


def test_count_discordant_clusters_uses_contingency():
    wnn_lbl = ["0"] * 10 + ["1"] * 10
    rna_lbl = ["a"] * 10 + ["b"] * 6 + ["c"] * 4
    n_disc, details = wnn._count_discordant_clusters(wnn_lbl, rna_lbl)
    assert n_disc == 1
    assert details["0"]["discordant"] is False
    assert details["1"]["modal_rna_purity"] == pytest.approx(0.6)
    assert details["1"]["discordant"] is True


def test_l2_normalize_rows_unit_norm_and_zero_safe():
    out = wnn._l2_normalize_rows(np.array([[3.0, 4.0], [0.0, 0.0]]))
    assert np.allclose(out[0], [0.6, 0.8])
    assert np.allclose(out[1], [0.0, 0.0])


# ═══════════════════════════════════════════════════════════════════════════
# Group B — real WNN on a synthetic paired .h5mu (aria-chromatin-env)
# ═══════════════════════════════════════════════════════════════════════════

N_PER_GROUP = 25
N_GENES = 60
N_PEAKS = 80
GROUPS = ("group_a", "group_b", "group_c")
FAST_PARAMS = {"k": 10, "n_pcs": 10, "n_lsi_comps": 12, "seed": 0}


def _synthetic_counts(seed: int = 0, *, dominant_depth_peak: bool = False):
    rng = np.random.default_rng(seed)
    n_cells = N_PER_GROUP * len(GROUPS)
    labels = np.repeat(np.arange(len(GROUPS)), N_PER_GROUP)

    rna_lam = np.full((n_cells, N_GENES), 0.5)
    atac_lam = np.full((n_cells, N_PEAKS), 0.3)
    gene_block = N_GENES // len(GROUPS)
    peak_block = N_PEAKS // len(GROUPS)
    for g in range(len(GROUPS)):
        rows = labels == g
        rna_lam[np.ix_(rows, np.arange(g * gene_block, (g + 1) * gene_block))] = 15.0
        atac_lam[np.ix_(rows, np.arange(g * peak_block, (g + 1) * peak_block))] = 6.0

    rna = rng.poisson(rna_lam).astype(np.float32)
    atac = rng.poisson(atac_lam).astype(np.float32)
    if dominant_depth_peak:
        # Last peak dominates total counts and scales with a per-cell depth
        # factor unrelated to group identity -> a depth-driven LSI axis.
        depth = rng.uniform(1.0, 40.0, size=n_cells)
        atac[:, -1] = np.round(depth * 25.0).astype(np.float32)
    barcodes = [f"cell{i:03d}-1" for i in range(n_cells)]
    group_labels = [GROUPS[g] for g in labels]
    return barcodes, rna, atac, group_labels


def _write_synthetic_h5mu(path: Path, *, seed: int = 0, rna_x=None,
                          dominant_depth_peak: bool = False) -> Path:
    import anndata as ad
    import mudata as md
    import pandas as pd
    import scipy.sparse as sp

    barcodes, rna, atac, group_labels = _synthetic_counts(
        seed, dominant_depth_peak=dominant_depth_peak
    )
    rna_mat = rna if rna_x is None else rna_x
    rna_ad = ad.AnnData(
        X=sp.csr_matrix(rna_mat),
        obs=pd.DataFrame(index=pd.Index(barcodes)),
        var=pd.DataFrame(index=pd.Index([f"gene{i:02d}" for i in range(N_GENES)])),
    )
    atac_ad = ad.AnnData(
        X=sp.csr_matrix(atac),
        obs=pd.DataFrame(index=pd.Index(barcodes)),
        var=pd.DataFrame(
            index=pd.Index([f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(N_PEAKS)])
        ),
    )
    mdata = md.MuData({"rna": rna_ad, "atac": atac_ad})
    mdata.obs["group_label"] = pd.Categorical(group_labels)
    mdata.write_h5mu(str(path))
    return path


@pytest.fixture(scope="module")
def muon_stack():
    pytest.importorskip("muon")
    return pytest.importorskip("mudata")


@pytest.fixture(scope="module")
def wnn_runs(muon_stack, tmp_path_factory):
    """Run WNN twice (seed=0) on the same synthetic .h5mu."""
    base = tmp_path_factory.mktemp("wnn_i2")
    h5mu = _write_synthetic_h5mu(base / "synthetic_paired.h5mu")
    res1 = wnn.integration_wnn({
        "h5mu_path": str(h5mu), "output_dir": str(base / "run1"), **FAST_PARAMS,
    })
    res2 = wnn.integration_wnn({
        "h5mu_path": str(h5mu), "output_dir": str(base / "run2"), **FAST_PARAMS,
    })
    return {"input": h5mu, "base": base, "r1": res1, "r2": res2}


def test_wnn_runs_successfully_on_small_matrix(wnn_runs):
    res = wnn_runs["r1"]
    assert res["status"] == "success", res
    assert res["input_kind"] == "paired_h5mu"
    assert res["n_cells"] == N_PER_GROUP * len(GROUPS)
    assert res["n_genes"] == N_GENES
    assert res["n_peaks"] == N_PEAKS
    assert res["seed"] == 0
    # <200 cells is warned about, not crashed on.
    assert any("200" in w for w in res["warnings"])
    assert res["n_joint_clusters"] >= 2
    assert res["discordance_method"] == "contingency_modal_purity_lt_0.70"
    assert isinstance(res["n_discordant_clusters"], int)
    assert 0 <= res["n_discordant_clusters"] <= res["n_joint_clusters"]
    assert set(res["cluster_concordance"]) == set(res["joint_cluster_sizes"])
    assert res["mean_rna_weight"] is not None
    assert res["mean_atac_weight"] is not None
    assert res["mean_rna_weight"] + res["mean_atac_weight"] == pytest.approx(1.0, abs=1e-3)
    assert set(res["per_cluster_modality_weights"]) == set(res["joint_cluster_sizes"])
    for entry in res["per_cluster_modality_weights"].values():
        assert entry["dominant_modality"] in {"rna", "atac", "balanced"}
    for key in ("wnn_vs_rna_ari", "wnn_vs_atac_ari", "rna_vs_atac_ari",
                "wnn_vs_rna_nmi", "wnn_vs_atac_nmi", "rna_vs_atac_nmi"):
        assert -1.0 <= res[key] <= 1.0
    # Three well-separated populations: joint partition recovers them.
    from aria.benchmarks.concordance_atac import adjusted_rand_index
    barcodes, _, _, group_labels = _synthetic_counts(0)
    joint = [res["joint_clusters"][bc] for bc in barcodes]
    assert adjusted_rand_index(joint, group_labels) > 0.8
    assert "atac_frip" in res["metrics_not_computed"]
    assert res["atac_frip"] is None


def test_wnn_ipc_result_satisfies_beta_contract(wnn_runs):
    pytest.importorskip("pydantic")
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    assert SCRIPT_CONTRACTS[WNN_SCRIPT].validate_result(wnn_runs["r1"]) == []


def test_wnn_is_strictly_deterministic(wnn_runs):
    import mudata as md
    from aria.benchmarks.concordance_atac import adjusted_rand_index

    r1, r2 = wnn_runs["r1"], wnn_runs["r2"]
    assert r1["status"] == r2["status"] == "success"
    bcs = sorted(r1["joint_clusters"])
    assert bcs == sorted(r2["joint_clusters"])
    l1 = [r1["joint_clusters"][b] for b in bcs]
    l2 = [r2["joint_clusters"][b] for b in bcs]
    assert adjusted_rand_index(l1, l2) == 1.0
    assert r1["per_cluster_modality_weights"] == r2["per_cluster_modality_weights"]

    m1 = md.read_h5mu(r1["output_path"])
    m2 = md.read_h5mu(r2["output_path"])
    assert list(m1.obs_names) == list(m2.obs_names)
    for col in ("rna_weight", "atac_weight"):
        assert np.allclose(
            np.asarray(m1.obs[col], dtype=float),
            np.asarray(m2.obs[col], dtype=float),
            rtol=0.0, atol=1e-12,
        )
    assert np.allclose(
        np.asarray(m1.obsm["X_wnn_umap"]), np.asarray(m2.obsm["X_wnn_umap"]),
        rtol=0.0, atol=1e-6,
    )
    stab = score_seed_stability(
        [l1, l2],
        [np.asarray(m1.obs["rna_weight"], dtype=float),
         np.asarray(m2.obs["rna_weight"], dtype=float)],
    )
    assert stab["mean_pairwise_ari"] == 1.0
    assert stab["max_weight_abs_diff"] == 0.0


def test_wnn_h5mu_round_trip_preserves_counts_and_annotations(wnn_runs):
    import mudata as md
    import scipy.sparse as sp

    res = wnn_runs["r1"]
    out_path = Path(res["output_path"])
    assert out_path.name == "wnn_integrated.h5mu"
    assert out_path.exists()
    assert not list(out_path.parent.glob(".wnn_integrated.tmp.*"))

    mdata = md.read_h5mu(str(out_path))
    barcodes, rna, atac, group_labels = _synthetic_counts(0)
    assert list(mdata.obs_names) == barcodes

    rna_mod = mdata.mod["rna"]
    atac_mod = mdata.mod["atac"]
    assert rna_mod.n_vars == N_GENES
    assert atac_mod.n_vars == N_PEAKS
    rna_counts = rna_mod.layers["counts"]
    rna_counts = rna_counts.toarray() if sp.issparse(rna_counts) else np.asarray(rna_counts)
    atac_counts = atac_mod.layers["counts"]
    atac_counts = atac_counts.toarray() if sp.issparse(atac_counts) else np.asarray(atac_counts)
    assert np.array_equal(rna_counts, rna)
    assert np.array_equal(atac_counts, atac)
    assert np.all(rna_counts == np.round(rna_counts))
    assert "highly_variable" in rna_mod.var.columns

    assert [str(x) for x in mdata.obs["group_label"]] == group_labels
    for col in ("leiden_wnn", "leiden_rna", "leiden_atac", "rna_weight", "atac_weight"):
        assert col in mdata.obs.columns, col
    assert [str(x) for x in mdata.obs["leiden_wnn"]] == [
        res["joint_clusters"][bc] for bc in barcodes
    ]
    assert mdata.obsm["X_wnn_umap"].shape == (len(barcodes), 2)
    assert "wnn_connectivities" in mdata.obsp
    assert mdata.obsp["wnn_connectivities"].shape == (len(barcodes), len(barcodes))
    assert "X_pca_l2" in rna_mod.obsm
    assert "X_lsi_l2" in atac_mod.obsm
    assert np.allclose(np.linalg.norm(rna_mod.obsm["X_pca_l2"], axis=1), 1.0, atol=1e-5)
    assert np.allclose(np.linalg.norm(atac_mod.obsm["X_lsi_l2"], axis=1), 1.0, atol=1e-5)

    meta = mdata.uns["aria_wnn"]
    assert meta["method"] == "muon_wnn_l2"
    assert int(meta["seed"]) == 0
    assert int(meta["k"]) == res["k"]
    assert int(meta["n_multineighbors"]) <= len(barcodes) - 2
    assert int(meta["n_bandwidth_neighbors"]) <= len(barcodes) - 2


def test_wnn_output_is_reentrant(wnn_runs):
    res = wnn_runs["r1"]
    again = wnn.integration_wnn({
        "h5mu_path": res["output_path"],
        "output_dir": str(wnn_runs["base"] / "reentrant"),
        **FAST_PARAMS,
    })
    assert again["status"] == "success", again
    assert again["rna_counts_source"] == "layers['counts']"
    assert again["atac_counts_source"] in {"X", "layers['counts']"}
    assert again["n_genes"] == N_GENES
    assert again["joint_clusters"] == res["joint_clusters"]


def test_wnn_rejects_fractional_counts_without_counts_layer(muon_stack, tmp_path):
    _, rna, _, _ = _synthetic_counts(0)
    h5mu = _write_synthetic_h5mu(tmp_path / "fractional.h5mu", rna_x=rna + 0.25)
    out = wnn.integration_wnn({
        "h5mu_path": str(h5mu), "output_dir": str(tmp_path / "out"), **FAST_PARAMS,
    })
    assert out["status"] == "error"
    assert out["error_type"] == "NonIntegerCounts"
    assert not (tmp_path / "out" / "wnn_integrated.h5mu").exists()


def _assert_depth_filter_consistent(res: dict, cutoff: float):
    corrs = res["lsi_depth_correlations"]
    dropped = res["dropped_lsi_components"]
    assert len(corrs) == res["n_lsi_components_raw"]
    assert res["n_lsi_components_used"] == res["n_lsi_components_raw"] - len(dropped)
    for idx, r_val in enumerate(corrs):
        # Correlations are rounded to 4 decimals in the payload.
        if idx in dropped:
            assert r_val >= cutoff - 1e-4
        else:
            assert r_val < cutoff + 1e-4


def test_lsi_depth_filter_drops_measured_component(muon_stack, tmp_path):
    h5mu = _write_synthetic_h5mu(
        tmp_path / "depth.h5mu", seed=1, dominant_depth_peak=True
    )
    cutoff = 0.5
    res = wnn.integration_wnn({
        "h5mu_path": str(h5mu), "output_dir": str(tmp_path / "out"),
        "depth_corr_cutoff": cutoff, **FAST_PARAMS,
    })
    assert res["status"] == "success", res
    assert res["dropped_lsi_components"], res["lsi_depth_correlations"]
    _assert_depth_filter_consistent(res, cutoff)
    assert any("depth-associated" in w for w in res["warnings"])


def test_lsi_depth_filter_does_not_blindly_drop_first_component(muon_stack, tmp_path):
    h5mu = _write_synthetic_h5mu(tmp_path / "nodrop.h5mu", seed=2)
    cutoff = 0.9999
    res = wnn.integration_wnn({
        "h5mu_path": str(h5mu), "output_dir": str(tmp_path / "out"),
        "depth_corr_cutoff": cutoff, **FAST_PARAMS,
    })
    assert res["status"] == "success", res
    _assert_depth_filter_consistent(res, cutoff)
    if res["lsi_depth_correlations"][0] < cutoff:
        assert 0 not in res["dropped_lsi_components"]


# ═══════════════════════════════════════════════════════════════════════════
# Group B.4 — cross-stack IPC E2E (aria-env -> aria-chromatin-env)
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.skipif(
    not CHROMATIN_PY.exists(),
    reason="aria-chromatin-env interpreter not available on this machine",
)
def test_cross_stack_ipc_wnn_on_paired_h5ad(tmp_path):
    pytest.importorskip("pydantic")
    ad = pytest.importorskip("anndata")
    import pandas as pd
    import scipy.sparse as sp

    from aria.utils.environment_manager import EnvironmentManager
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    barcodes, rna, atac, _ = _synthetic_counts(0)
    rna_path = tmp_path / "rna.h5ad"
    atac_path = tmp_path / "atac.h5ad"
    ad.AnnData(
        X=sp.csr_matrix(rna),
        obs=pd.DataFrame(index=pd.Index(barcodes)),
        var=pd.DataFrame(index=pd.Index([f"gene{i:02d}" for i in range(N_GENES)])),
    ).write_h5ad(str(rna_path))
    ad.AnnData(
        X=sp.csr_matrix(atac),
        obs=pd.DataFrame(index=pd.Index(barcodes)),
        var=pd.DataFrame(
            index=pd.Index([f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(N_PEAKS)])
        ),
    ).write_h5ad(str(atac_path))

    env = EnvironmentManager(workspace_dir=str(tmp_path / "workspace"))
    result = env.run_in_stack(
        stack="chromatin",
        script_path=WNN_SCRIPT,
        params={
            "rna_files": [str(rna_path)],
            "atac_files": [str(atac_path)],
            "output_dir": str(tmp_path / "wnn_out"),
            **FAST_PARAMS,
        },
    )
    assert result["status"] == "success", result
    assert result["input_kind"] == "paired_h5ad"
    assert Path(result["output_path"]).exists()
    assert SCRIPT_CONTRACTS[WNN_SCRIPT].validate_result(result) == []


# ═══════════════════════════════════════════════════════════════════════════
# Group C — committed external concordance artifact invariants
# ═══════════════════════════════════════════════════════════════════════════

def test_wnn_external_concordance_manifest_invariants():
    import json

    manifest_path = (
        Path(__file__).resolve().parents[1]
        / "docs"
        / "benchmark_results"
        / "multiome_concordance"
        / "v48_i2_wnn_concordance.json"
    )
    if not manifest_path.is_file():
        pytest.skip("WNN external concordance manifest not present")

    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert data["benchmark"] == "V48_I2_wnn_external_concordance"
    assert data["seurat_status"]["status"] == "success"

    conc = data["concordance"]
    assert conc["status"] == "success"
    assert conc["n_shared_cells"] >= 2000
    assert conc["seed_stability"]["mean_pairwise_ari"] == 1.0
    assert conc["seed_stability"]["max_weight_abs_diff"] == 0.0
    assert conc["cluster_ari"] >= 0.85
    assert conc["cluster_nmi"] >= 0.85
    assert conc["modality_weights"]["spearman_rho"] >= 0.85
    assert conc["modality_weights"]["mae"] <= 0.08
    assert conc["knn_overlap"]["mean_jaccard"] >= 0.65
    assert conc["label_recovery"]["status"] == "success"
    assert conc["label_recovery"]["aria_wnn_nmi"] >= 0.80

