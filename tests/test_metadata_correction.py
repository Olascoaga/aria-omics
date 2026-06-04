"""Real-run bug (2026-06-04): CHECKPOINT 1 option [3] "Correct metadata" did
nothing — the orchestrator had no branch for it and proceeded to design with the
unchanged (possibly wrong) modality/organism. The user must be able to correct
the modality (RNA/ATAC/any ARIA modality) and the species from there.
"""

from aria.agents.data_audit_agent import (
    apply_metadata_corrections,
    SUPPORTED_MODALITIES,
    default_genome_for_organism,
)


def test_correction_rekeys_modality_and_sets_species():
    ctx = {
        "modalities": {"bulk_ATAC": ["/a/B1_1.fq.gz", "/a/B1_2.fq.gz"]},
        "organism": "Homo sapiens", "genome": "hg38",
    }
    out = apply_metadata_corrections(
        ctx, {"modality": "bulk_RNA_raw", "organism": "Mus musculus", "genome": "mm10"})
    assert set(out["modalities"]) == {"bulk_RNA_raw"}
    assert out["modalities"]["bulk_RNA_raw"] == ["/a/B1_1.fq.gz", "/a/B1_2.fq.gz"]
    assert out["organism"] == "Mus musculus"
    assert out["genome"] == "mm10"


def test_correction_collapses_all_files_to_chosen_modality():
    ctx = {"modalities": {"bulk_ATAC": ["/a/x.fq.gz"], "unknown": ["/a/y.fq.gz"]}}
    out = apply_metadata_corrections(ctx, {"modality": "bulk_RNA_raw"})
    assert set(out["modalities"]) == {"bulk_RNA_raw"}
    assert sorted(out["modalities"]["bulk_RNA_raw"]) == ["/a/x.fq.gz", "/a/y.fq.gz"]


def test_correction_partial_only_organism_leaves_modality():
    ctx = {"modalities": {"bulk_RNA_raw": ["/a/x"]}, "organism": "unknown",
           "genome": "unknown"}
    out = apply_metadata_corrections(ctx, {"organism": "Mus musculus"})
    assert out["organism"] == "Mus musculus"
    assert set(out["modalities"]) == {"bulk_RNA_raw"}


def test_correction_none_is_noop():
    ctx = {"modalities": {"bulk_RNA_raw": ["/a/x"]}, "organism": "Homo sapiens"}
    assert apply_metadata_corrections(ctx, None) is ctx
    assert apply_metadata_corrections(ctx, {}) is ctx


def test_supported_modalities_offers_rna_and_atac_and_chromatin():
    for m in ("bulk_RNA_raw", "bulk_RNA", "scRNA", "scATAC", "bulk_ATAC"):
        assert m in SUPPORTED_MODALITIES


def test_default_genome_for_organism():
    assert default_genome_for_organism("Homo sapiens") == "hg38"
    assert default_genome_for_organism("Mus musculus") == "mm10"
    assert default_genome_for_organism("Totally unknown organism") is None


def test_orchestrator_after_checkpoint_1_applies_corrections(monkeypatch):
    # The orchestrator must actually apply CP1 corrections to exp_context before
    # design (the bug was that "Correct metadata" silently proceeded unchanged).
    import pytest
    pytest.importorskip("litellm")
    from aria.agents.orchestrator_agent import OrchestratorAgent
    import aria.agents.design_agent as da

    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch.llm = object()
    orch._experiment_plans = {"exp1": {"intent": {}}}
    orch.publish_finding = lambda *a, **k: None
    orch.memory = type("M", (), {"store_decision": lambda self, **k: None})()

    # Stop right after corrections are applied (design stubbed to "failed").
    monkeypatch.setattr(da.DesignAgent, "__init__", lambda self, memory, llm: None)
    monkeypatch.setattr(da.DesignAgent, "start_design",
                        lambda self, **k: {"status": "failed", "reason": "stub"})

    exp_context = {"modalities": {"bulk_ATAC": ["/a/x.fq.gz", "/a/y.fq.gz"]},
                   "organism": "Homo sapiens", "genome": "hg38"}
    msg = type("Msg", (), {"payload": {
        "context": {"exp_context": exp_context}, "question": "q"}})()

    orch._after_checkpoint_1(
        "exp1", "Correct metadata", msg,
        corrections={"modality": "bulk_RNA_raw", "organism": "Mus musculus",
                     "genome": "mm10"})

    assert set(exp_context["modalities"]) == {"bulk_RNA_raw"}
    assert exp_context["modalities"]["bulk_RNA_raw"] == ["/a/x.fq.gz", "/a/y.fq.gz"]
    assert exp_context["organism"] == "Mus musculus"
    assert exp_context["genome"] == "mm10"
