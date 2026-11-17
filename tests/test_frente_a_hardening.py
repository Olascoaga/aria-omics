"""Regression guards for Frente A readiness in ARIA:
1. Mouse canonical default is ``mm39`` across ``DataAuditAgent`` and ``SetupAgent``,
   while explicit ``mm10`` / ``GRCm38`` / ``hg19`` declarations are never silently
   overwritten by generic species heuristics.
2. Per-project ``memory.db`` isolation via ``ARIA_MEMORY_DB`` and
   ``run_headless(..., memory_db_path=...)``, plus ``ParameterAdvisor``
   historical-bonus suppression in ``reproducible_mode`` and full 32-hex UUIDs.
3. ``aria.doctor._check_llm`` executes ``ARIA_DOCTOR_LLM_PROBE=1`` cleanly with
   ``privacy.egress_allowed()`` (zero-argument contract).
"""

from __future__ import annotations

from pathlib import Path
import pytest

from aria.agents.data_audit_agent import (
    DataAuditAgent,
    apply_metadata_corrections,
    default_genome_for_organism,
)
from aria.memory.memory import ARIAMemory, resolve_memory_db_path


def _bare_audit_agent() -> DataAuditAgent:
    return DataAuditAgent.__new__(DataAuditAgent)


def test_mouse_defaults_to_mm39_and_preserves_explicit_mm10_or_mm39():
    da = _bare_audit_agent()
    assert default_genome_for_organism("Mus musculus") == "mm39"

    # Generic mouse cohort -> mm39
    g_gen, org_gen = da._infer_genome_organism(
        [Path("sample1_R1.fastq.gz")],
        Path("/data/mouse_cohort"),
        "Compare WT vs KO mouse cortex",
    )
    assert org_gen == "Mus musculus"
    assert g_gen == "mm39"

    # Explicit mm39 alongside word "mouse" -> mm39 (never shadowed by mm10)
    g_39, org_39 = da._infer_genome_organism(
        [Path("grad_rep1.bam")],
        Path("/data/CUTANDRUN/mm39_gradientes"),
        "Analyze mouse CUT&RUN gradients",
    )
    assert org_39 == "Mus musculus"
    assert g_39 == "mm39"

    # Explicit mm10 / GRCm38 alongside word "mouse" -> mm10 preserved
    g_10, org_10 = da._infer_genome_organism(
        [Path("legacy_GRCm38_rep1.bam")],
        Path("/data/mouse_legacy"),
        "Analyze mouse chromatin",
    )
    assert org_10 == "Mus musculus"
    assert g_10 == "mm10"

    # Explicit hg19 alongside word "human" -> hg19 preserved
    g_hg19, org_hg19 = da._infer_genome_organism(
        [Path("patient1_hg19.bam")],
        Path("/data/human_cohort"),
        "Analyze human ChIP-seq",
    )
    assert org_hg19 == "Homo sapiens"
    assert g_hg19 == "hg19"


def test_apply_metadata_corrections_fills_default_mm39_when_genome_unknown():
    ctx = {"modalities": {"CUT_AND_RUN": ["/a/s1.bam"]}, "organism": "unknown", "genome": "unknown"}
    out = apply_metadata_corrections(ctx, {"organism": "Mus musculus"})
    assert out["organism"] == "Mus musculus"
    assert out["genome"] == "mm39"


def test_setup_agent_honors_audited_genome_and_refuses_cross_assembly_substitution(tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.setup_agent import GENOME_REGISTRY, SetupAgent

    assert "mm39" in GENOME_REGISTRY
    assert "mm10" in GENOME_REGISTRY

    sa = SetupAgent.__new__(SetupAgent)
    sa.genomes_dir = tmp_path / "genomes"
    sa.genomes_dir.mkdir(parents=True, exist_ok=True)

    # Generic mouse without audited genome -> mm39
    assert sa._detect_genome("Mus musculus", {"summary": "WT vs KO"}, audited_genome="unknown") == "mm39"
    # Audited mm39 -> mm39
    assert sa._detect_genome("Mus musculus", {"summary": "WT vs KO"}, audited_genome="mm39") == "mm39"
    # Audited mm10 / GRCm38 -> mm10 (never overwritten with mm39)
    assert sa._detect_genome("Mus musculus", {"summary": "WT vs KO"}, audited_genome="mm10") == "mm10"
    assert sa._detect_genome("Mus musculus", {"summary": "WT vs KO"}, audited_genome="GRCm38") == "mm10"
    # Audited hg19 (not in auto-download registry) -> None rather than silently substituting hg38!
    assert sa._detect_genome("Homo sapiens", {"summary": "human cohort"}, audited_genome="hg19") is None

    warnings: list[str] = []
    cfg = sa._local_staged_genome_cfg("hg19", warnings)
    assert cfg["genome_key"] == "hg19"
    assert cfg["fasta_ready"] is False
    assert any("refusing to substitute" in w for w in warnings)


def test_aria_memory_db_isolation_env_and_headless(monkeypatch, tmp_path):
    proj_a_db = tmp_path / "proj_a" / "memory.db"
    proj_b_db = tmp_path / "proj_b" / "memory.db"

    monkeypatch.setenv("ARIA_MEMORY_DB", str(proj_a_db))
    assert resolve_memory_db_path() == str(proj_a_db)

    mem_a = ARIAMemory()
    mem_a.create_wing("exp_a", "Project A", organism="Mus musculus", genome="mm39")
    mem_a.store_decision("dec_a", "exp_a", 3, "leiden_clustering resolution", "0.8")
    mem_a.close()

    # Explicit db_path overrides ARIA_MEMORY_DB cleanly
    mem_b = ARIAMemory(db_path=str(proj_b_db))
    assert mem_b.list_wings() == []
    mem_b.close()

    # Verify run_headless forwards memory_db_path to ARIAMemory
    pytest.importorskip("litellm")
    import aria.headless as headless
    import aria.agents.orchestrator_agent as orch_mod

    captured: dict[str, str] = {}

    class _StubOrch:
        def __init__(self, memory, llm=None):
            captured["db_path"] = memory.db_path

        def run(self, experiment_id, ctx):
            return {"status": "not_started"}

    monkeypatch.setattr(orch_mod, "OrchestratorAgent", _StubOrch)
    res = headless.run_headless(
        str(tmp_path),
        "Analyze Oligos cohort",
        memory_db_path=str(proj_b_db),
    )
    assert res.status == "failed"
    assert captured["db_path"] == str(proj_b_db)


def test_parameter_advisor_disables_historical_bonus_in_reproducible_mode_and_uses_full_uuid():
    pytest.importorskip("litellm")
    from aria.llm.parameter_advisor import ParameterAdvisor, ParameterCandidate

    adv = ParameterAdvisor.__new__(ParameterAdvisor)
    a = ParameterCandidate(value=1.0, score=0.50)
    b = ParameterCandidate(value=0.4, score=0.49)
    hist = [{"decision": "0.4"}, {"decision": "0.4"}]

    # Interactive mode preserves bounded tie-breaking
    assert adv._choose_best([a, b], hist, disable_historical_bonus=False).value == 0.4
    # Reproducible / batch mode suppresses cross-run historical drift
    assert ParameterAdvisor._historical_bonus_disabled({"reproducible_mode": True}) is True
    assert adv._choose_best([a, b], hist, disable_historical_bonus=True).value == 1.0

    # Full 32-hex UUID check via advise_wnn_k with stubbed memory + LLM
    class _Mem:
        def list_wings(self):
            return []

        def store_decision(self, **kw):
            self.last = kw

    class _LLM:
        def complete(self, **kw):
            return "Objective balance favored k."

    mem = _Mem()
    full_adv = ParameterAdvisor(mem, _LLM())
    full_adv.MetricEvaluator = None
    monkeypatch_eval = lambda adata, k: {"rna_weight": 0.5, "atac_weight": 0.5, "knn_overlap": 0.7}
    from aria.llm import parameter_advisor as pa_mod
    orig_wnn = pa_mod.MetricEvaluator.wnn_k
    pa_mod.MetricEvaluator.wnn_k = staticmethod(monkeypatch_eval)
    try:
        dec = full_adv.advise_wnn_k(object(), "exp_uuid", {"reproducible_mode": True})
    finally:
        pa_mod.MetricEvaluator.wnn_k = orig_wnn
    assert len(dec.decision_id) == 32
    int(dec.decision_id, 16)  # valid 128-bit hex UUID


def test_doctor_llm_probe_respects_zero_arg_egress_allowed(monkeypatch):
    from aria.doctor import _check_llm

    monkeypatch.setenv("ARIA_DOCTOR_LLM_PROBE", "1")
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    issues_air = _check_llm()
    assert any(i.code == "llm_latency" and "skipped" in i.message for i in issues_air)

    # Also test non-air-gapped path with no commercial keys configured (zero network traffic, no TypeError)
    monkeypatch.setenv("ARIA_AIR_GAPPED", "0")
    for k in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GEMINI_API_KEY", "GOOGLE_API_KEY"):
        monkeypatch.delenv(k, raising=False)
    issues_open = _check_llm()
    assert any(i.code == "llm_mode" for i in issues_open)
