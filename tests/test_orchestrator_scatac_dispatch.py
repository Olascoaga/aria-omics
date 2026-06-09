"""Guards for the chromatin-only (scATAC) orchestrator dispatch path.

Failing-first reproduction of the live-validation blocker (2026-06-08): a
chromatin-only run (no RNA modality, no inferred DE groups) made the
DE-oriented DesignAgent return ``failed/no_samples``, which the orchestrator
turned into ``cancelled`` BEFORE the ChromatinAgent could be dispatched. The
fake-env unit tests called the agent directly and never exercised this gate;
only the live orchestrator/bus path did.

Fix (Option A): DesignAgent skips the DE design phase for chromatin-only runs
and synthesizes a minimal ``no_de_design`` design; the orchestrator proceeds
straight to plan/CP2 -> audit gate -> CP3.5 ack -> dispatch.

Per ADR-011 these tests use neutral synthetic labels only.
"""

import pytest

from aria.agents.design_agent import DesignAgent
from aria.memory.memory import ARIAMemory


def _agent(tmp_path):
    return DesignAgent(memory=ARIAMemory(db_path=str(tmp_path / "mem.db")), llm=None)


def test_design_skips_de_phase_for_chromatin_only_run(tmp_path):
    """scATAC-only: start_design must SKIP (not fail) with a minimal design."""
    agent = _agent(tmp_path)
    exp_context = {
        "modalities": {"scATAC": ["/synthetic/atac.h5mu"]},
        "organism": "Homo sapiens",
        "genome": "hg38",
    }
    result = agent.start_design(
        experiment_id="exp_scatac",
        exp_context=exp_context,
        biological_intent={"question": "characterize chromatin accessibility"},
    )

    assert result["status"] == "skipped", result
    assert result["reason"] == "no_de_design_needed"
    design = result["design"]
    assert design["design_type"] == "no_de_design"
    assert design["organism"] == "Homo sapiens"
    assert design["genome"] == "hg38"
    assert design["groups"] == {}
    assert design["main_factor"] is None
    assert design["modalities_without_de"] == ["scATAC"]


@pytest.mark.parametrize("modality", ["bulk_ATAC", "ChIP", "CUT_AND_RUN"])
def test_design_skips_de_phase_for_other_chromatin_modalities(tmp_path, modality):
    agent = _agent(tmp_path)
    result = agent.start_design(
        experiment_id="exp_chrom",
        exp_context={"modalities": {modality: ["/synthetic/x.bam"]}},
        biological_intent={"question": "accessibility"},
    )
    assert result["status"] == "skipped", result
    assert result["design"]["design_type"] == "no_de_design"


def test_design_still_fails_when_no_modalities_at_all(tmp_path):
    """Regression guard: a truly empty run must still fail with no_samples."""
    agent = _agent(tmp_path)
    result = agent.start_design(
        experiment_id="exp_empty",
        exp_context={"modalities": {}},
        biological_intent={"question": "nothing"},
    )
    assert result["status"] == "failed"
    assert result["reason"] == "no_samples"


def test_rna_run_still_uses_de_design_phase(tmp_path):
    """An scRNA run with inferred groups must NOT take the skip path; it must
    publish the groups checkpoint (awaiting_user), keeping the DE path intact."""
    agent = _agent(tmp_path)
    inferred_design = {
        "source": "h5ad_obs",
        "confidence": "high",
        "organism": "Homo sapiens",
        "genome": "hg38",
        "groups": {"COND_A": ["alpha", "beta", "gamma"],
                   "COND_B": ["alpha", "beta", "gamma"]},
        "main_factor": "grp",
        "condition_col": "grp",
        "replicate_col": "rep",
    }
    result = agent.start_design(
        experiment_id="exp_rna",
        exp_context={"modalities": {"scRNA": ["/synthetic/rna.h5ad"]},
                     "inferred_design": inferred_design,
                     "organism": "Homo sapiens", "genome": "hg38"},
        biological_intent={"question": "compare A vs B"},
    )
    assert result["status"] == "awaiting_user", result
    assert result["step"] == "groups"


def test_chromatin_qc_finding_handles_uncomputed_frip_tss(tmp_path):
    """Live-validation guard: a pre-called peak-matrix (.h5mu) QC honestly
    reports FRiP/TSS as None. Publishing that finding must NOT crash on the
    None contract (it did: 'NoneType' < float), and must report dimensions and
    'not computed' instead of a fabricated value.
    """
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.bus.message_bus import bus, MessageType

    memory = ARIAMemory(db_path=str(tmp_path / "chrom.db"))
    experiment_id = "exp_qc_none"
    memory.create_wing(experiment_id, name="scatac", organism="Homo sapiens",
                       genome="hg38")
    agent = ChromatinAgent(memory=memory, llm=None)

    qc_result = {  # shape returned by chromatin_qc on the real HC11 .h5mu
        "status": "success",
        "data_type": "scATAC",
        "n_cells": 3143,
        "n_peaks": 60990,
        "frip": None,
        "tss_enrichment": None,
        "warnings": ["frip (requires called peaks — run chromatin_peaks first)"],
    }

    # Must not raise.
    agent._publish_qc_finding(experiment_id, qc_result, "scATAC")

    findings = [m for m in bus.get_log()
                if m.experiment_id == experiment_id
                and m.type == MessageType.FINDING]
    assert findings, "no QC finding published"
    payload = findings[-1].payload
    assert payload["frip"] is None and payload["tss_score"] is None
    assert "not computed" in payload["summary"]
    assert "3,143 cells" in payload["summary"]
    assert "60,990 peaks" in payload["summary"]
    memory.close()


def _nested_chromatin_agent_result():
    """The shape ChromatinAgent.run() actually returns (per-modality wrapper
    with analysis findings nested one level deeper)."""
    return {
        "status": "done",
        "findings": {
            "scATAC": {
                "status": "done",
                "findings": {
                    "qc": {"status": "success", "data_type": "scATAC",
                           "n_cells": 3143, "n_peaks": 60990,
                           "frip": None, "tss_enrichment": None,
                           "metrics_not_computed": ["frip", "tss_enrichment"]},
                    "lsi": {"status": "success", "n_clusters": 8,
                            "n_cells_used": 3143, "n_peaks": 60990,
                            "dropped_components": [0]},
                    "differential_accessibility": {
                        "status": "success",
                        "per_cluster": {"n_da_peaks": 13294, "n_clusters": 8}},
                    "motifs": {"status": "success", "n_enriched": 5,
                               "collection": "JASPAR2024_CORE_vertebrates"},
                },
            }
        },
    }


def test_unwrap_chromatin_findings_flattens_modality_wrapper():
    from aria.agents.narrative.narrators.chromatin import (
        unwrap_chromatin_findings,
    )
    flat = unwrap_chromatin_findings(_nested_chromatin_agent_result())
    assert set(flat) >= {"qc", "lsi", "differential_accessibility", "motifs"}
    # A flat dict (older test / direct feed) is returned unchanged.
    already_flat = {"findings": {"qc": {"status": "success"}}}
    assert "qc" in unwrap_chromatin_findings(already_flat)


def test_chromatin_narrator_emits_blocks_from_nested_agent_result():
    """Live-validation guard: the narrator must produce LSI/DA/motif blocks from
    the REAL nested agent return, not only from a flat findings feed."""
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    narr = ChromatinNarrator()
    result = _nested_chromatin_agent_result()
    assert narr.accepts("chromatin_agent", result)
    blocks = narr.collect("chromatin_agent", result)
    analyses = {b.analysis for b in blocks}
    assert "qc" in analyses
    assert "clustering" in analyses or "dimensionality_reduction" in analyses, analyses
    assert any("motif" in a for a in analyses), analyses


def test_run_ledger_reconciles_nested_chromatin_findings_as_ran():
    """Live-validation guard: the run ledger must mark LSI/DA/motifs as 'ran'
    (not 'planned but not run') from the real nested agent return. Peak calling
    stays not-run honestly for a pre-called peak matrix (.h5mu)."""
    from aria.agents.narrative.run_ledger import build_run_ledger

    agent_results = {"chromatin_agent": _nested_chromatin_agent_result()}
    exp_ctx = {"user_question": "characterize chromatin accessibility via LSI "
                                "clustering and motif enrichment"}
    ledger = build_run_ledger(exp_ctx, agent_results)
    by_analysis = {e["analysis"]: e for e in ledger["entries"]
                   if e["modality"] == "chromatin"}
    assert by_analysis["qc"]["status"] == "ran", by_analysis["qc"]
    assert by_analysis["dimensionality_reduction"]["status"] == "ran"
    assert by_analysis["differential_accessibility"]["status"] == "ran"
    assert by_analysis["motif_enrichment"]["status"] == "ran"


def test_executive_summary_concrete_block_includes_chromatin_outputs(tmp_path):
    """Live-validation guard: the anti-hallucination CONCRETE block fed to the
    executive-summary LLM must report the scATAC LSI/DA/motif outputs, so the
    summary cannot claim 'no dimensionality reduction outputs were recorded'."""
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent(memory=ARIAMemory(db_path=str(tmp_path / "n.db")),
                           llm=None)
    result = _nested_chromatin_agent_result()
    # enrich DA/motif findings to the shapes the summarizer reads
    da = result["findings"]["scATAC"]["findings"]["differential_accessibility"]
    da["per_cluster"] = {"ran": True, "n_da_total": 9957,
                         "n_da_by_cluster": {str(i): 1 for i in range(8)}}
    mot = result["findings"]["scATAC"]["findings"]["motifs"]
    mot["ran"] = True
    mot["motif_source"] = {"n_motifs": 879,
                           "collection": "JASPAR2024_CORE_vertebrates"}
    mot["per_group"] = {str(i): {"n_enriched": (1 if i < 5 else 0)}
                        for i in range(8)}

    concrete = agent._summarize_agent_results_for_llm(
        {"chromatin_agent": result})
    assert "CHROMATIN" in concrete
    assert "LSI/clustering RAN" in concrete
    assert "differential accessibility RAN" in concrete
    assert "9,957" in concrete
    assert "motif enrichment RAN" in concrete
    memory = agent.memory
    memory.close()


def test_orchestrator_dispatches_chromatin_only_run_after_cp1(tmp_path):
    """Wiring guard: resolving CP1 on a chromatin-only run must reach plan_ready
    (CP2 published), not cancel. Drives _after_checkpoint_1 directly with a
    constructed CP1 message; llm=None falls back to the deterministic plan.
    """
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.bus.message_bus import Message, MessageType, bus

    memory = ARIAMemory(db_path=str(tmp_path / "orch.db"))
    experiment_id = "exp_orch_scatac"
    memory.create_wing(experiment_id, name="scatac", organism="Homo sapiens",
                       genome="hg38")
    orch = OrchestratorAgent(memory)

    exp_context = {
        "modalities": {"scATAC": ["/synthetic/atac.h5mu"]},
        "organism": "Homo sapiens",
        "genome": "hg38",
        "user_question": "characterize chromatin accessibility",
    }
    orch._experiment_plans[experiment_id] = {
        "intent": {"analysis_type": "cell_type"},
        "context": {"user_question": exp_context["user_question"]},
    }
    session = orch._get_session(experiment_id)
    session.intent = {"analysis_type": "cell_type"}

    cp1 = Message(
        sender="data_audit_agent",
        type=MessageType.ESCALATION,
        checkpoint=1,
        experiment_id=experiment_id,
        payload={
            "checkpoint": 1,
            "question": "Is this correct?",
            "options": ["Confirm and continue", "Cancel"],
            "context": {"exp_context": exp_context},
        },
    )

    result = orch._after_checkpoint_1(experiment_id, "Confirm and continue", cp1)

    assert result["status"] == "plan_ready", result
    # The CP2 plan must include the chromatin agent for scATAC.
    agents_in_plan = {s.get("agent") for s in result["plan"].get("steps", [])}
    assert "chromatin_agent" in agents_in_plan, result["plan"]
    # A CP2 escalation must now be pending for this experiment.
    pend = [m for m in bus.get_pending_checkpoints(experiment_id=experiment_id)
            if m.payload.get("checkpoint") == 2]
    assert pend, "CP2 was not published for the chromatin-only run"
    memory.close()
