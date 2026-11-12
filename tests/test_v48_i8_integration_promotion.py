"""V48 Slice I8 (ADR-071) — IntegrationAgent scaffold -> beta promotion,
pre-dispatch readiness/registry/lock gate, CP3.5 acknowledgement enforcement,
and post-QC MuData propagation tests.
"""

from __future__ import annotations

import ast
import copy
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from aria.agents.audit_agent import AuditAgent
from aria.agents.data_audit_agent import DataAuditAgent
from aria.agents.integration_agent import IntegrationAgent
from aria.agents.modality_audit import build_capability_matrix
from aria.agents.orchestrator_agent import (
    INTEGRATION_VALIDATION,
    MODALITY_VALIDATION,
    MODALITY_VALIDATION_LEVELS,
    OrchestratorAgent,
)
from aria.agents.setup_agent import SetupAgent
from aria.memory.memory import ARIAMemory
from aria.runtime.experiment_session import ExperimentSession
from aria.utils.environment_audit import audit_environment
from aria.utils.environment_specs import ENVIRONMENT_SPECS, lock_environment_names
from aria.utils.multiome_contracts import (
    enforce_same_cell_integration_gate,
    infer_multiome_contract,
    validate_multiome_contract,
)
from aria.utils.registry_integrity import check_registry_integrity

ROOT = Path(__file__).resolve().parents[1]

PAIRED_H5MU_CTX: dict[str, Any] = {
    "genome": "hg38",
    "organism": "Homo sapiens",
    "modalities": {"scATAC": ["/virtual/hc11_paired.h5mu"]},
    "assay_detections": [
        {
            "path": "/virtual/hc11_paired.h5mu",
            "modality": "scATAC",
            "confidence": "high",
            "evidence": {
                "format": "h5mu",
                "modalities": ["rna", "atac"],
                "paired_rna": True,
            },
        }
    ],
}


# ── 1. Registry & ADR-056 / ADR-071 Tier Governance ─────────────────────────


def test_integration_validation_is_beta_ack_gated_and_never_production():
    assert INTEGRATION_VALIDATION["level"] == "beta"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is True
    assert INTEGRATION_VALIDATION["scope"] == "paired_same_cell_rna_atac"
    assert INTEGRATION_VALIDATION["production_promoted"] is False
    assert INTEGRATION_VALIDATION["promotion_policy_adr"] == "ADR-056"
    assert INTEGRATION_VALIDATION["promotion_adr"] == "ADR-071"

    # Per-assay modality table remains strictly 9 per-assay keys (1-to-1 with docs)
    assert set(MODALITY_VALIDATION) == {
        "scRNA",
        "bulk_RNA",
        "bulk_RNA_raw",
        "scATAC",
        "bulk_ATAC",
        "ChIP",
        "CUT_AND_RUN",
        "CUT_AND_TAG",
        "HiC",
    }
    assert MODALITY_VALIDATION_LEVELS["scATAC"] == "production"
    assert MODALITY_VALIDATION_LEVELS["bulk_ATAC"] == "beta"
    assert "integration" not in MODALITY_VALIDATION


def test_registry_integrity_enforces_adr071_invariants_as_hard_errors(monkeypatch):
    assert [i for i in check_registry_integrity() if i.severity == "error"] == []

    bad_cases = [
        ({"level": "scaffold", "dispatch_enabled": True}, "scaffold_dispatch_enabled"),
        (
            {
                "level": "beta",
                "dispatch_enabled": True,
                "requires_ack": False,
                "production_promoted": False,
                "scope": "paired_same_cell_rna_atac",
            },
            "integration_beta_missing_ack_gate",
        ),
        (
            {
                "level": "beta",
                "dispatch_enabled": True,
                "requires_ack": True,
                "production_promoted": True,
                "scope": "paired_same_cell_rna_atac",
            },
            "integration_beta_production_flag_invalid",
        ),
        (
            {
                "level": "beta",
                "dispatch_enabled": True,
                "requires_ack": True,
                "production_promoted": False,
                "scope": "unpaired_any",
            },
            "integration_beta_scope_invalid",
        ),
        (
            {
                "level": "production",
                "dispatch_enabled": True,
            },
            "integration_premature_production_tier",
        ),
    ]
    for mutated, expected_code in bad_cases:
        monkeypatch.setattr(
            "aria.agents.orchestrator_agent.INTEGRATION_VALIDATION",
            mutated,
        )
        issues = check_registry_integrity()
        errors = [i for i in issues if i.severity == "error"]
        assert any(i.code == expected_code for i in errors), (
            f"Expected hard error {expected_code!r}, got {issues}"
        )


# ── 2. Paired .h5mu Detection, Multi-Container Ambiguity & Plan Routing ──────


def test_atac_only_h5mu_is_not_promoted_to_paired_mudata_by_filename_fallback():
    atac_only_ctx = {
        "genome": "hg38",
        "modalities": {"scATAC": ["/virtual/atac_only.h5mu"]},
        "assay_detections": [
            {
                "path": "/virtual/atac_only.h5mu",
                "modality": "scATAC",
                "confidence": "high",
                "evidence": {
                    "format": "h5mu",
                    "modalities": ["atac"],
                    "paired_rna": False,
                },
            }
        ],
    }
    contract = infer_multiome_contract(atac_only_ctx)
    assert contract["object_type"] != "paired_mudata"
    assert contract["same_cell"] is False
    gate = enforce_same_cell_integration_gate(atac_only_ctx)
    assert gate["allowed"] is False


def test_multiple_unmerged_paired_mudata_containers_fail_closed():
    multi_ctx = {
        "genome": "hg38",
        "modalities": {
            "scATAC": ["/virtual/donor1.h5mu", "/virtual/donor2.h5mu"],
        },
        "assay_detections": [
            {
                "path": "/virtual/donor1.h5mu",
                "modality": "scATAC",
                "evidence": {"format": "h5mu", "paired_rna": True},
            },
            {
                "path": "/virtual/donor2.h5mu",
                "modality": "scATAC",
                "evidence": {"format": "h5mu", "paired_rna": True},
            },
        ],
    }
    gate = enforce_same_cell_integration_gate(multi_ctx)
    assert gate["allowed"] is False
    assert any(
        b["check"] == "multiome_multiple_mudata_unsupported"
        for b in gate["blockers"]
    )


def test_orchestrator_fallback_and_llm_plans_set_integration_needed_for_paired_h5mu():
    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch._sessions = {}
    orch._experiment_plans = {"exp_paired_plan": {"intent": {}}}

    fallback = orch._fallback_plan(PAIRED_H5MU_CTX)
    assert fallback["integration_needed"] is True
    assert fallback["integration_type"] == "WNN"

    # Even if an LLM omits/disables integration_needed on single-key modalities={"scATAC": [.h5mu]},
    # _design_analysis_plan deterministically enforces integration_needed=True.
    orch.think_structured = lambda *a, **kw: {
        "steps": [
            {
                "order": 1,
                "agent": "chromatin_agent",
                "analysis": "scATAC QC",
                "depends_on": [],
                "can_parallel": False,
            }
        ],
        "contrasts": [],
        "integration_needed": False,
        "integration_type": "none",
        "estimated_complexity": "medium",
        "rationale": "LLM omitted integration",
    }
    plan = orch._design_analysis_plan("exp_paired_plan", dict(PAIRED_H5MU_CTX))
    assert plan["integration_needed"] is True
    assert plan["integration_type"] == "WNN"


# ── 3. Capability Matrix, AuditAgent & CP3.5 Acknowledgement Lifecycle ──────


def test_capability_matrix_and_audit_agent_require_cp35_ack_for_paired_multiome(tmp_path):
    # Explicit same-cell multiome with scRNA-only modality card (green) still triggers
    # integration_readiness requires_ack_for_dispatch -> AuditAgent status == "blocking"
    explicit_same_cell_ctx = {
        "genome": "hg38",
        "modalities": {
            "scRNA": ["/virtual/rna.h5ad"],
        },
        "design": {
            "groups": {"ctrl": ["c1", "c2", "c3"], "stim": ["s1", "s2", "s3"]},
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col": "cell_type",
            },
        },
        "multiome_contract": {
            "object_type": "explicit_same_cell",
            "same_cell": True,
            "cell_namespace": "shared_10x_barcodes",
            "atac_feature_space": "peak_matrix",
        },
    }
    matrix = build_capability_matrix(
        explicit_same_cell_ctx,
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert matrix["dispatch"]["requires_ack"] == []
    assert matrix["dispatch"]["blocked"] == []
    int_ready = matrix["integration_readiness"]
    assert int_ready["tier"] == "beta"
    assert int_ready["requires_ack"] is True
    assert int_ready["requires_ack_for_dispatch"] is True
    assert int_ready["production_promoted"] is False
    assert int_ready["governance_blockers"]
    assert any(
        f["check"] == "integration_beta_ack_required" for f in matrix["findings"]
    )

    mem = ARIAMemory(db_path=str(tmp_path / "mem.sqlite"))
    try:
        auditor = AuditAgent(memory=mem, llm=object())
        auditor.publish_status = lambda *a, **kw: None
        audit_res = auditor.run_audit(
            explicit_same_cell_ctx,
            "exp_int_ack",
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )
        assert audit_res["status"] == "blocking"
    finally:
        mem.close()


def test_cp35_ack_lifecycle_strict_affirmative_memory_logging_and_reaudit_reset(
    tmp_path, monkeypatch
):
    mem = ARIAMemory(db_path=str(tmp_path / "cp35.sqlite"))
    try:
        mem.create_wing("exp_cp35", name="hc11", organism="Homo sapiens", genome="hg38")
        orch = OrchestratorAgent.__new__(OrchestratorAgent)
        orch.memory = mem
        orch.llm = object()
        orch._sessions = {}
        orch._pending_checkpoints = {}
        orch._pending_dispatch = {}
        orch._experiment_plans = {"exp_cp35": {"intent": {}}}
        orch._agent_results = {}
        orch._log_handlers = {}
        escalations: list[dict[str, Any]] = []
        orch.publish_status = lambda *a, **kw: None
        orch.publish_escalation = lambda **kw: escalations.append(kw)

        dispatched: list[tuple[str, dict, dict]] = []

        class _NoopThread:
            def __init__(self, target, args, daemon):
                self.args = args

            def start(self):
                dispatched.append(self.args)

        monkeypatch.setattr("aria.agents.orchestrator_agent.threading.Thread", _NoopThread)

        plan = {"steps": [], "integration_needed": True}
        gate_res = orch._gate_and_dispatch(
            "exp_cp35",
            plan,
            {**copy.deepcopy(PAIRED_H5MU_CTX), "decision_actor": "human_reviewer"},
        )
        assert gate_res["status"] == "audit_blocking"
        assert len(escalations) == 1
        assert escalations[0]["checkpoint"] == 3.5
        assert "integration_agent is beta" in escalations[0]["question"]
        assert orch._get_session("exp_cp35").cp35_acknowledged is False

        # Ambiguous / non-affirmative response must fail closed (cancelled, never dispatched)
        ambiguous = orch._after_audit_checkpoint("exp_cp35", "maybe later", None)
        assert ambiguous["status"] == "cancelled"
        assert dispatched == []
        assert orch._get_session("exp_cp35").cp35_acknowledged is False

        # Re-gate and acknowledge affirmatively with "Proceed anyway"
        gate_res_2 = orch._gate_and_dispatch(
            "exp_cp35",
            plan,
            {**copy.deepcopy(PAIRED_H5MU_CTX), "decision_actor": "human_reviewer"},
        )
        assert gate_res_2["status"] == "audit_blocking"
        ack_res = orch._after_audit_checkpoint("exp_cp35", "Proceed anyway", None)
        assert ack_res["status"] == "analysis_running"
        assert len(dispatched) == 1
        session = orch._get_session("exp_cp35")
        assert session.cp35_acknowledged is True
        assert session.cp35_decision == "Proceed anyway"
        assert session.cp35_decision_actor == "human_reviewer"

        decisions = mem.get_decisions("exp_cp35")
        cp35_rows = [d for d in decisions if str(d["checkpoint"]) == "3.5"]
        assert len(cp35_rows) == 1
        assert cp35_rows[0]["decision"] == "Proceed anyway"
        assert cp35_rows[0]["made_by"] == "human_reviewer"

        # Re-running _gate_and_dispatch resets stale cp35_acknowledged back to False
        orch._gate_and_dispatch("exp_cp35", plan, copy.deepcopy(PAIRED_H5MU_CTX))
        assert orch._get_session("exp_cp35").cp35_acknowledged is False
    finally:
        mem.close()


# ── 4. Orchestrator _dispatch_agents Gate Order & Single-Lane Preservation ──


def _make_stub_orchestrator():
    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch._sessions = {}
    orch._experiment_plans = {}
    orch._agent_results = {}
    orch._log_handlers = {}
    orch.memory = object()
    orch.llm = object()
    orch.publish_status = lambda *a, **kw: None
    orch.publish_escalation = lambda *a, **kw: None
    findings: list[dict[str, Any]] = []
    orch.publish_finding = lambda exp, payload, conf: findings.append(payload)
    calls: list[tuple[str, dict[str, Any]]] = []

    def _fake_run_agent(agent_name: str, experiment_id: str, context: dict):
        calls.append((agent_name, context))
        return {"status": "done", "agent": agent_name}

    orch._run_agent = _fake_run_agent
    return orch, calls, findings


@pytest.mark.parametrize(
    "exp_id,modalities,expected_lanes",
    [
        (
            "exp_split_unpaired",
            {"scRNA": ["/virtual/rna.h5ad"], "scATAC": ["/virtual/atac.h5ad"]},
            ["scrna_agent", "chromatin_agent"],
        ),
        (
            "exp_scrna_bulkrna",
            {"scRNA": ["/virtual/rna.h5ad"], "bulk_RNA": ["/virtual/counts.tsv"]},
            ["scrna_agent", "bulk_rna_agent"],
        ),
        (
            "exp_bulk_rna_atac",
            {"bulk_RNA": ["/virtual/counts.tsv"], "bulk_ATAC": ["/virtual/peaks.bam"]},
            ["bulk_rna_agent", "chromatin_agent"],
        ),
    ],
)
def test_dispatch_blocks_non_same_cell_integration_while_running_single_lanes(
    exp_id, modalities, expected_lanes
):
    orch, calls, _ = _make_stub_orchestrator()
    orch._experiment_plans[exp_id] = {"intent": {}}
    ctx = {"genome": "hg38", "modalities": modalities, "integration_ack_confirmed": True}
    plan = orch._fallback_plan(ctx)

    orch._dispatch_agents(exp_id, plan, ctx)

    called_names = [name for name, _ in calls]
    for lane in expected_lanes:
        assert lane in called_names
    assert "integration_agent" not in called_names
    assert "narrative_agent" in called_names

    int_res = orch._agent_results[exp_id]["integration_agent"]
    assert int_res["status"] == "skipped"
    assert int_res["reason"] == "multiome_same_cell_gate_blocked"
    assert int_res["validation_level"] == "beta"


def test_dispatch_requires_boolean_cp35_ack_before_running_paired_integration():
    orch, calls, _ = _make_stub_orchestrator()

    # 1. Unacknowledged paired .h5mu (including truthy non-bool string "true") fails closed
    for idx, bad_ctx in enumerate(
        [
            copy.deepcopy(PAIRED_H5MU_CTX),
            {**copy.deepcopy(PAIRED_H5MU_CTX), "integration_ack_confirmed": "true"},
            {**copy.deepcopy(PAIRED_H5MU_CTX), "cp35_acknowledged": "yes"},
        ]
    ):
        exp_id = f"exp_no_ack_{idx}"
        orch._experiment_plans[exp_id] = {"intent": {}}
        plan = orch._fallback_plan(bad_ctx)
        orch._dispatch_agents(exp_id, plan, bad_ctx)
        int_res = orch._agent_results[exp_id]["integration_agent"]
        assert int_res["status"] == "skipped"
        assert int_res["reason"] == "integration_ack_required"

    assert "integration_agent" not in [name for name, _ in calls]

    # 2. Session-level CP3.5 acknowledgement dispatches IntegrationAgent cleanly
    orch._experiment_plans["exp_ack_ok"] = {"intent": {"user_question": "joint states"}}
    orch._get_session("exp_ack_ok").cp35_acknowledged = True
    plan_ok = orch._fallback_plan(PAIRED_H5MU_CTX)
    orch._dispatch_agents("exp_ack_ok", plan_ok, copy.deepcopy(PAIRED_H5MU_CTX))

    called_ok = [name for name, _ in calls]
    assert "chromatin_agent" in called_ok
    assert "integration_agent" in called_ok
    assert orch._agent_results["exp_ack_ok"]["integration_agent"]["status"] == "done"


# ── 5. IntegrationAgent Internal Propagation, Strategy & Setup Fixes ────────


def test_integration_agent_blocks_non_wnn_strategies_without_same_cell_contract():
    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.publish_finding = lambda *a, **kw: None
    agent.publish_status = lambda *a, **kw: None
    agent.publish_escalation = lambda *a, **kw: None
    calls: list[str] = []
    agent.env = SimpleNamespace(
        run_in_stack=lambda stack, script_path, params, timeout=None: calls.append(script_path)
    )

    # Bulk RNA + Bulk ATAC selects 'peak_to_gene_only' and MUST still fail closed!
    out_p2g = agent.run(
        "exp_bulk_bypass",
        {
            "exp_context": {
                "genome": "hg38",
                "modalities": {
                    "bulk_RNA": ["/virtual/counts.tsv"],
                    "bulk_ATAC": ["/virtual/atac.bam"],
                },
            },
            "biological_intent": {"user_question": "link peaks to genes"},
        },
    )
    assert out_p2g["strategy"] == "peak_to_gene_only"
    assert out_p2g["status"] == "blocked"
    assert out_p2g["reason"] == "multiome_same_cell_gate_blocked"
    assert calls == []

    # scRNA + ChIP with MOFA keyword selects 'mofa' and MUST also fail closed!
    out_mofa = agent.run(
        "exp_mofa_bypass",
        {
            "exp_context": {
                "genome": "hg38",
                "modalities": {
                    "scRNA": ["/virtual/rna.h5ad"],
                    "ChIP": ["/virtual/chip.bam"],
                },
            },
            "biological_intent": {"user_question": "run MOFA factor decomposition"},
        },
    )
    assert out_mofa["strategy"] == "mofa"
    assert out_mofa["status"] == "blocked"
    assert out_mofa["reason"] == "multiome_same_cell_gate_blocked"
    assert calls == []


def test_post_qc_validated_h5mu_overrides_conflicting_raw_paths_across_all_stages(
    tmp_path,
):
    raw_h5mu = tmp_path / "raw_unfiltered.h5mu"
    raw_h5mu.write_bytes(b"RAW")
    post_qc_h5mu = tmp_path / "paired_qc_filtered.h5mu"
    post_qc_h5mu.write_bytes(b"POST_QC")
    wnn_out_h5mu = tmp_path / "wnn_integrated.h5mu"
    wnn_out_h5mu.write_bytes(b"WNN")

    recorded: list[tuple[str, str, dict[str, Any]]] = []

    class _Env:
        def run_in_stack(self, stack: str, script_path: str, params: dict, timeout=None):
            name = Path(script_path).name
            recorded.append((stack, name, copy.deepcopy(params)))
            if name == "integration_paired_qc.py":
                return {
                    "status": "success",
                    "output_path": str(post_qc_h5mu),
                    "n_cells_before": 3143,
                    "n_cells_after": 2517,
                    "n_genes": 3000,
                    "n_peaks": 10000,
                }
            if name == "integration_wnn.py":
                return {
                    "status": "success",
                    "output_path": str(wnn_out_h5mu),
                    "n_cells": 2517,
                    "mean_rna_weight": 0.55,
                    "mean_atac_weight": 0.45,
                    "n_joint_clusters": 11,
                    "n_rna_only_clusters": 11,
                    "n_discordant_clusters": 2,
                    "warnings": [],
                }
            if name == "integration_peak2gene.py":
                return {
                    "status": "success",
                    "output_path": str(tmp_path / "p2g.csv"),
                    "n_links": 10,
                    "n_positive_correlations": 8,
                    "n_negative_correlations": 2,
                    "top_links": [],
                    "warnings": [],
                }
            if name == "integration_eregulons.py":
                return {
                    "status": "success",
                    "output_path": str(tmp_path / "ereg.csv"),
                    "triplets_path": str(tmp_path / "triplets.csv"),
                    "n_eregulons": 3,
                    "n_triplets": 25,
                    "n_unique_tfs": 3,
                    "n_pos_tf_pos_p2g": 3,
                    "n_neg_tf_pos_p2g": 0,
                    "n_paralog_ambiguous": 0,
                    "n_depth_confounded": 0,
                    "top_eregulons": [],
                    "associative_only": True,
                    "perturbation_validated": False,
                    "chromatin_conformation_validated": False,
                    "warnings": [],
                }
            if name == "integration_mofa.py":
                return {
                    "status": "success",
                    "output_path": str(tmp_path / "mofa.hdf5"),
                    "n_factors": 8,
                    "n_cells": 2517,
                    "variance_explained": {"scRNA": 0.25, "scATAC": 0.14},
                    "top_factors": [],
                    "factor1_top_features": [],
                    "technical_factor_flag": False,
                    "technical_factor_id": None,
                    "warnings": [],
                }
            raise AssertionError(f"Unexpected script {name}")

    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.env = _Env()
    agent.memory = SimpleNamespace(
        store_decision=lambda *a, **kw: None,
        create_tunnel=lambda *a, **kw: None,
    )
    agent.llm = SimpleNamespace(
        complete=lambda *a, **kw: "ok",
        complete_heavy=lambda *a, **kw: "ok",
        complete_medium=lambda *a, **kw: "ok",
        complete_light=lambda *a, **kw: "ok",
    )
    agent.advisor = SimpleNamespace(
        advise_wnn_k=lambda **kw: SimpleNamespace(
            decision_id="wnn_k",
            chosen_value=20,
            justification="ok",
            candidates=[],
            chosen_by="advisor",
            approved_by_user=False,
        )
    )
    agent.publish_status = lambda *a, **kw: None
    agent.publish_finding = lambda *a, **kw: None
    agent.publish_blocking_escalation = lambda **kw: (
        "msg",
        {"choice": kw["options"][0]},
    )

    exp_ctx = {
        "genome": "hg38",
        "organism": "Homo sapiens",
        "modalities": {"scATAC": [str(raw_h5mu)]},
        "assay_detections": [
            {
                "path": str(raw_h5mu),
                "modality": "scATAC",
                "evidence": {"format": "h5mu", "paired_rna": True},
            }
        ],
        # Conflicting raw h5mu_path and per-stage overrides MUST NOT bypass post-QC container!
        "h5mu_path": str(raw_h5mu),
        "enable_mofa": True,
        "genome_config": {
            "gtf": "/ref/genes.gtf",
            "fasta": "/ref/genome.fa",
        },
        "wnn_params": {"h5mu_path": str(raw_h5mu)},
        "peak2gene_params": {"h5mu_path": str(raw_h5mu)},
        "eregulon_params": {"h5mu_path": str(raw_h5mu)},
        "mofa_params": {"h5mu_path": str(raw_h5mu), "max_cells": 2517},
    }

    res = agent.run(
        "exp_prop",
        {
            "exp_context": exp_ctx,
            "biological_intent": {"user_question": "multiome cell states"},
        },
    )
    assert res["status"] == "done"
    assert res["strategy"] == "wnn_then_mofa"
    # WNN output_path preserved for IntegrationNarrator Figure 1
    assert res["findings"]["wnn"]["output_path"] == str(wnn_out_h5mu)

    by_script = {name: params for _, name, params in recorded}
    assert by_script["integration_paired_qc.py"]["h5mu_path"] == str(raw_h5mu)
    for downstream in (
        "integration_wnn.py",
        "integration_peak2gene.py",
        "integration_eregulons.py",
        "integration_mofa.py",
    ):
        assert by_script[downstream]["h5mu_path"] == str(post_qc_h5mu), downstream
    assert by_script["integration_mofa.py"]["modalities"]["scRNA"] == [str(post_qc_h5mu)]
    assert by_script["integration_mofa.py"]["modalities"]["scATAC"] == [str(post_qc_h5mu)]
    assert by_script["integration_mofa.py"]["max_cells"] == 2517
    assert by_script["integration_peak2gene.py"]["gtf_path"] == "/ref/genes.gtf"
    assert by_script["integration_peak2gene.py"]["genome_fasta"] == "/ref/genome.fa"
    assert by_script["integration_eregulons.py"]["gtf_path"] == "/ref/genes.gtf"
    assert by_script["integration_eregulons.py"]["genome_fasta"] == "/ref/genome.fa"


def test_paired_qc_missing_output_file_fails_closed(tmp_path):
    raw_h5mu = tmp_path / "raw.h5mu"
    raw_h5mu.write_bytes(b"RAW")
    missing_post_qc = tmp_path / "does_not_exist.h5mu"

    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.env = SimpleNamespace(
        run_in_stack=lambda stack, script_path, params, timeout=None: {
            "status": "success",
            "output_path": str(missing_post_qc),
        }
    )
    agent.publish_status = lambda *a, **kw: None
    agent.publish_finding = lambda *a, **kw: None

    out = agent.run(
        "exp_missing_qc_out",
        {
            "exp_context": {
                "genome": "hg38",
                "modalities": {"scATAC": [str(raw_h5mu)]},
                "assay_detections": [
                    {
                        "path": str(raw_h5mu),
                        "modality": "scATAC",
                        "evidence": {"format": "h5mu", "paired_rna": True},
                    }
                ],
            },
            "biological_intent": {"user_question": "cell types"},
        },
    )
    assert out["status"] == "failed"
    assert out["reason"] == "paired_qc_missing_output"


def test_select_strategy_enable_mofa_precedence_and_setup_agent_envs():
    agent = IntegrationAgent.__new__(IntegrationAgent)

    # enable_mofa=False strictly overrides keywords, n_modalities >= 3, and mofa_params
    assert (
        agent._select_strategy(
            has_scrna=True,
            has_scatac=True,
            has_hic=True,
            n_modalities=3,
            intent={"user_question": "Run MOFA latent factor decomposition"},
            exp_ctx={"enable_mofa": False, "mofa_params": {"max_cells": 1000}},
        )
        == "wnn"
    )
    # enable_mofa=True promotes 2-modality paired run without keywords
    assert (
        agent._select_strategy(
            has_scrna=True,
            has_scatac=True,
            has_hic=False,
            n_modalities=2,
            intent={"user_question": "Identify cell types"},
            exp_ctx={"enable_mofa": True},
        )
        == "wnn_then_mofa"
    )

    setup = SetupAgent.__new__(SetupAgent)
    envs_h5mu = setup._needed_envs({"scATAC": ["hc11_paired.h5mu"]}, has_fastq=False)
    assert "aria-chromatin-env" in envs_h5mu
    assert "aria-integration-env" in envs_h5mu

    envs_split = setup._needed_envs(
        {"scRNA": ["rna.h5ad"], "scATAC": ["atac.h5ad"]}, has_fastq=False
    )
    assert "aria-integration-env" in envs_split


def test_data_audit_forwards_multiome_passthrough_keys(tmp_path):
    h5mu = tmp_path / "paired.h5mu"
    h5mu.write_bytes(b"")
    mem = ARIAMemory(db_path=str(tmp_path / "da.sqlite"))
    try:
        da = DataAuditAgent(mem, llm=object())
        da._last_scan_report = {"truncated": False}
        da._scan_directory = lambda *a, **kw: [h5mu]
        da._classify_files = lambda *a, **kw: {"scATAC": [str(h5mu)]}
        da._last_assay_detections = [
            {
                "path": str(h5mu),
                "modality": "scATAC",
                "confidence": "high",
                "evidence": {"format": "h5mu", "paired_rna": True},
            }
        ]
        da._detect_organism = lambda *a, **kw: ("hg38", "Homo sapiens")
        da._validate_design = lambda *a, **kw: []
        da._inspect_h5ad_metadata = lambda *a, **kw: {}
        da.publish_status = lambda *a, **kw: None
        da.publish_escalation = lambda *a, **kw: None

        res = da.run(
            "exp_da",
            {
                "data_dir": str(tmp_path),
                "user_question": "multiome",
                "enable_mofa": True,
                "gtf_path": "/ref/gencode.gtf",
                "genome_fasta": "/ref/hg38.fa",
                "decision_actor": "headless_policy",
            },
        )
        ctx = res["exp_context"]
        assert ctx["enable_mofa"] is True
        assert ctx["gtf_path"] == "/ref/gencode.gtf"
        assert ctx["genome_fasta"] == "/ref/hg38.fa"
        assert ctx["decision_actor"] == "headless_policy"
        assert ctx["multiome_contract"]["same_cell"] is True
    finally:
        mem.close()


# ── 6. Hermetic Lock Contract & Anti-Fabrication Guard on Slice I8 Files ────


def test_integration_env_is_locked_active_scientific_runtime():
    spec = ENVIRONMENT_SPECS["integration"]
    assert spec.env_name == "aria-integration-env"
    assert spec.lock_required is True
    assert spec.lock_policy == "active_scientific_runtime"
    assert "aria-integration-env" in lock_environment_names()

    conda_lock = ROOT / "envs" / "aria-integration-env.linux-64.lock"
    pip_lock = ROOT / "envs" / "aria-integration-env.pip.lock"
    assert conda_lock.is_file()
    assert pip_lock.is_file()

    report = audit_environment(
        "aria-integration-env",
        env_name="aria-integration-env",
    )
    if report.get("status") != "missing_conda":
        assert report["status"] == "pass", report
        assert report["conda_exact"] is True
        assert report["pip_exact"] is True


def test_slice_i8_modified_files_have_zero_fabrication_violations():
    import importlib.util

    fab_path = ROOT / "tests" / "test_no_fabrication_guard.py"
    spec = importlib.util.spec_from_file_location("_fab_guard_i8", fab_path)
    assert spec and spec.loader
    fab_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fab_mod)

    targets = [
        "aria/agents/orchestrator_agent.py",
        "aria/agents/integration_agent.py",
        "aria/agents/modality_audit.py",
        "aria/agents/audit_agent.py",
        "aria/agents/setup_agent.py",
        "aria/runtime/experiment_session.py",
        "aria/utils/multiome_contracts.py",
        "aria/utils/registry_integrity.py",
        "aria/utils/environment_specs.py",
    ]
    for rel in targets:
        src = (ROOT / rel).read_text(encoding="utf-8")
        ast.parse(src)
        violations = fab_mod.scan_source(rel, src)
        assert violations == [], f"Fabrication guard flagged {rel}: {violations}"


def test_integration_agent_output_dir_does_not_nest_recursively(tmp_path: Path):
    """Post-QC `.h5mu` lives inside `<input>/aria_integration/`; downstream stages
    must reuse that directory rather than nesting `aria_integration/aria_integration`."""
    agent = IntegrationAgent.__new__(IntegrationAgent)
    raw = tmp_path / "hc11_paired.h5mu"
    out1 = agent._get_output_dir([str(raw)])
    assert out1 == str(tmp_path / "aria_integration")
    post_qc = Path(out1) / "paired_qc_hc11.h5mu"
    out2 = agent._get_output_dir([str(post_qc)])
    assert out2 == out1
