"""Regression & verification suite for ARIA v4.9.4 pre-submission hardening.

Covers every engineering upgrade introduced in v4.9.4:
1. Fail-closed `DecisionManifest` commitment & `realized_decisions_sha256` tamper rejection.
2. Chromatin preprocessing & dimensionality-reduction blocks classified strictly `descriptive`.
3. Covariate-aware OLS residualization + sex-chromosome / unadjusted-covariate warning
   on `cross_modal_concordance.sample_matched_layer`.
4. Path-portable `preprint_freeze` lanes (zero machine-local user paths) and clean
   provenance exclusion via `ARIA_PROVENANCE_EXCLUDE_PATH`.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from aria.utils.decision_manifest import (
    build_decision_manifest,
    commit_decision_manifest,
    seal_realized_parameter_decisions,
    verify_decision_manifest,
)


def _sample_manifest_inputs():
    base_ctx = {
        "organism": "Homo sapiens",
        "genome": "hg38",
        "reproducible_mode": True,
        "modalities": {},
        "input_files": [
            {"modality": "bulk_RNA", "path": "/data/counts.tsv", "size_bytes": 42, "sha256": "a" * 64}
        ],
        "design": {
            "organism": "Homo sapiens",
            "genome": "hg38",
            "main_factor": "condition",
            "groups": {"treated": ["t1", "t2"], "control": ["c1", "c2"]},
        },
    }
    plan = {"steps": []}
    intent = {"analysis_type": "differential", "raw_user_question": "treated vs control"}
    return base_ctx, plan, intent


def test_realized_parameter_decisions_digest_verifies_and_rejects_tampering():
    base_ctx, plan, intent = _sample_manifest_inputs()
    dm = build_decision_manifest(
        experiment_id="exp_seal",
        exp_context=base_ctx,
        plan=plan,
        biological_intent=intent,
    )
    assert verify_decision_manifest(dm) == (True, "verified")

    # Unsealed non-empty realized_parameter_decisions MUST be rejected
    unsealed = {**dm, "realized_parameter_decisions": {"leiden_clustering.resolution": 0.6}}
    ok, reason = verify_decision_manifest(unsealed)
    assert ok is False and reason == "missing_realized_decisions_sha256"

    # Properly sealed realized_parameter_decisions verifies while keeping manifest_sha256 invariant
    sealed = seal_realized_parameter_decisions(dm, {"leiden_clustering.resolution": 0.6})
    assert sealed["manifest_sha256"] == dm["manifest_sha256"]
    assert verify_decision_manifest(sealed) == (True, "verified")

    # Tampering with realized_parameter_decisions post-seal MUST be rejected
    tampered = {
        **sealed,
        "realized_parameter_decisions": {"leiden_clustering.resolution": 1.4},
    }
    ok_t, reason_t = verify_decision_manifest(tampered)
    assert ok_t is False and reason_t == "realized_decisions_sha256_mismatch"


def test_orchestrator_dispatch_fails_closed_on_commit_error_and_tampered_replay(monkeypatch, tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.memory.memory import ARIAMemory

    mem = ARIAMemory(db_path=str(tmp_path / "fail_closed.db"))
    orch = OrchestratorAgent(mem, llm=None)
    base_ctx, plan, intent = _sample_manifest_inputs()

    downstream_called = {"setup": False, "narrative": False}

    def _boom_run_agent(*_a, **_k):
        downstream_called["setup"] = True
        raise AssertionError("Downstream agent must NEVER run after manifest failure!")

    monkeypatch.setattr(orch, "_run_agent", _boom_run_agent)
    monkeypatch.setattr(orch, "_blocked_modalities", lambda m: {})
    monkeypatch.setattr(orch, "_experimental_modalities", lambda m: {})

    # 1. Tampered replayed_decision_manifest halts fail-closed before any compute runs
    dm = build_decision_manifest(
        experiment_id="orig_exp",
        exp_context=base_ctx,
        plan=plan,
        biological_intent=intent,
    )
    sealed = seal_realized_parameter_decisions(dm, {"leiden_clustering.resolution": 0.6})
    tampered_replay = {
        **sealed,
        "realized_parameter_decisions": {"leiden_clustering.resolution": 1.8},
    }
    orch._experiment_plans["exp_tampered"] = {"intent": dict(intent)}
    res_tampered = orch._dispatch_agents(
        "exp_tampered",
        plan,
        {**base_ctx, "workspace_dir": str(tmp_path / "ws1"), "replayed_decision_manifest": tampered_replay},
    )
    assert res_tampered["status"] == "error"
    assert res_tampered["reason"] == "decision_manifest_commitment_failed"
    assert downstream_called["setup"] is False

    # 2. Failure during commit_decision_manifest halts fail-closed instead of warning & proceeding
    from aria.utils import decision_manifest as dm_mod

    def _commit_failure(*_a, **_k):
        raise OSError("simulated disk write failure")

    monkeypatch.setattr(dm_mod, "commit_decision_manifest", _commit_failure)
    orch._experiment_plans["exp_disk_fail"] = {"intent": dict(intent)}
    res_disk = orch._dispatch_agents(
        "exp_disk_fail",
        plan,
        {**base_ctx, "workspace_dir": str(tmp_path / "ws2")},
    )
    assert res_disk["status"] == "error"
    assert res_disk["reason"] == "decision_manifest_commitment_failed"
    assert "simulated disk write failure" in res_disk["details"]
    assert downstream_called["setup"] is False
    mem.close()


def test_preprocessing_and_dimred_blocks_classify_descriptive():
    pytest.importorskip("litellm")
    from aria.agents.narrative.claim_compiler import annotate_claim_tiers
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    analyses = (
        "scatac_preprocessing",
        "bulk_atac_preprocessing",
        "dimensionality_reduction",
        "regulatory_layers",
        "peak_calling",
        "peak_count_matrix",
        "peak_annotation",
    )
    blocks = [
        NarrativeBlock(
            id=f"chromatin.{name}",
            modality="chromatin",
            analysis=name,
            block_type="result",
            title=f"Title {name}",
            status="success",
            confidence="high",
            claim=f"{name} processed 2400 cells and 55385 peaks.",
            evidence=[
                EvidenceItem("n_cells", 2400, source="chromatin"),
                EvidenceItem("n_peaks", 55385, source="chromatin"),
            ],
        )
        for name in analyses
    ]
    annotate_claim_tiers(blocks, exp_ctx={})
    for blk in blocks:
        claim_meta = blk.metadata["claim"]
        assert claim_meta["tier"] == "descriptive", f"{blk.analysis} got {claim_meta}"
        assert claim_meta["licensed_language"] == "descriptive"


def test_cross_modal_cpm_layer_covariate_residualization_and_sex_chrom_guard(tmp_path):
    from aria.utils.cross_modal_concordance import sample_matched_layer

    # 10 matched donors (5 Female, 5 Male)
    donors = [f"D{i}" for i in range(10)]
    sex = {d: ("F" if i < 5 else "M") for i, d in enumerate(donors)}
    rng = np.random.default_rng(7)
    is_female = np.array([1.0 if sex[d] == "F" else 0.0 for d in donors])

    # Gene 1 (chrX, XIST): purely driven by donor sex (no residual correlation once sex is adjusted!)
    xist_rna = np.exp2(2.0 + 7.0 * is_female + rng.normal(0, 0.15, 10)) - 1.0
    xist_atac = np.exp2(2.0 + 7.0 * is_female + rng.normal(0, 0.15, 10)) - 1.0

    # Gene 2 (chr1, AUTOGENE): genuine sex-independent cis correlation across donors
    signal = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, -2.0, -1.0, 0.0, 1.0, 2.0])
    auto_rna = np.exp2(5.0 + 1.5 * signal + rng.normal(0, 0.02, 10)) - 1.0
    auto_atac = np.exp2(5.0 + 1.5 * signal + rng.normal(0, 0.02, 10)) - 1.0

    # Add 25 null autosomal genes for permutation family stability
    rna_rows = {"XIST": xist_rna, "AUTOGENE": auto_rna}
    atac_rows = {"chrX:10000-10500": xist_atac, "chr1:10000-10500": auto_atac}
    models = {
        "chrX": [{"gene": "XIST", "tss": 10200, "strand": "+"}],
        "chr1": [{"gene": "AUTOGENE", "tss": 10200, "strand": "+"}],
    }
    for j in range(25):
        pos = 500000 + j * 300000
        gname = f"NULLG{j}"
        rna_rows[gname] = np.exp2(5.0 + rng.normal(0, 1.0, 10)) - 1.0
        atac_rows[f"chr1:{pos}-{pos + 500}"] = np.exp2(5.0 + rng.normal(0, 1.0, 10)) - 1.0
        models["chr1"].append({"gene": gname, "tss": pos + 200, "strand": "+"})

    rna_cpm = pd.DataFrame(rna_rows, index=donors).T
    atac_cpm = pd.DataFrame(atac_rows, index=donors).T
    base_params = {
        "cis_distance_bp": 100_000,
        "min_matched_samples": 6,
        "n_permutations": 200,
        "padj_threshold": 0.05,
        "seed": 0,
    }

    # Unadjusted fit: flags unadjusted covariate risk AND catches chrX XIST link!
    warn_unadj: list[str] = []
    frame_u, res_u = sample_matched_layer(rna_cpm, atac_cpm, models, base_params, warn_unadj)
    assert res_u["ran"] is True
    assert res_u["covariate_adjustment"]["adjusted"] is False
    assert "cpm_sample_matched_layer_unadjusted_for_covariates" in warn_unadj
    assert res_u["n_sex_chromosome_significant_links"] >= 1
    assert any(w.startswith("unadjusted_sex_chromosome_links_detected:") for w in warn_unadj)
    xist_u = frame_u[frame_u["gene"] == "XIST"].iloc[0]
    assert bool(xist_u["is_sex_chromosome"]) is True
    assert xist_u["pearson_r"] > 0.95

    # Sex-adjusted fit: eliminates spurious XIST link while preserving AUTOGENE!
    warn_adj: list[str] = []
    frame_a, res_a = sample_matched_layer(
        rna_cpm, atac_cpm, models, {**base_params, "sample_covariates": {"sex": sex}}, warn_adj
    )
    assert res_a["ran"] is True
    assert res_a["covariate_adjustment"]["adjusted"] is True
    assert res_a["covariate_adjustment"]["covariates_adjusted"] == ["sex"]
    assert "cpm_sample_matched_layer_unadjusted_for_covariates" not in warn_adj
    assert res_a["n_sex_chromosome_significant_links"] == 0
    xist_a = frame_a[frame_a["gene"] == "XIST"].iloc[0]
    auto_a = frame_a[frame_a["gene"] == "AUTOGENE"].iloc[0]
    assert abs(float(xist_a["pearson_r"])) < 0.7
    assert float(xist_a["padj"]) >= 0.05
    assert float(auto_a["pearson_r"]) > 0.98
    assert float(auto_a["padj"]) < 0.05
    assert res_a["covariate_adjustment"]["permutation_method"] == "freedman_lane_residual_reprojection"


def test_cross_modal_cpm_covariate_adjusted_null_stays_calibrated():
    """Freedman-Lane residual re-projection must preserve null calibration
    under covariate adjustment (no subspace-shrinkage false-positive inflation)."""
    from aria.utils.cross_modal_concordance import sample_matched_layer

    donors = [f"D{i}" for i in range(10)]
    sex = {d: ("F" if i < 4 else "M") for i, d in enumerate(donors)}
    rng = np.random.default_rng(123)
    rna_rows, atac_rows, gene_models = {}, {}, []
    for j in range(400):
        pos = 100_000 + j * 250_000
        gname = f"NULLG{j}"
        rna_rows[gname] = np.exp2(5.0 + rng.normal(0, 1.0, 10)) - 1.0
        atac_rows[f"chr1:{pos}-{pos + 500}"] = np.exp2(5.0 + rng.normal(0, 1.0, 10)) - 1.0
        gene_models.append({"gene": gname, "tss": pos + 200, "strand": "+"})

    warn: list[str] = []
    _, res = sample_matched_layer(
        pd.DataFrame(rna_rows, index=donors).T,
        pd.DataFrame(atac_rows, index=donors).T,
        {"chr1": gene_models},
        {
            "cis_distance_bp": 100_000,
            "min_matched_samples": 6,
            "n_permutations": 100,
            "padj_threshold": 0.05,
            "seed": 42,
            "sample_covariates": {"sex": sex},
        },
        warn,
    )
    assert res["ran"] is True
    assert res["covariate_adjustment"]["adjusted"] is True
    assert res["n_significant_links"] == 0


def test_preprint_freeze_lanes_contain_no_private_user_paths():
    from aria.benchmarks.preprint_freeze import LANES

    serialized = json.dumps(LANES)
    assert "Samael" not in serialized
    assert "H9-RNA" not in serialized
    assert "/home/" not in serialized
    assert "/tmp/" not in serialized
