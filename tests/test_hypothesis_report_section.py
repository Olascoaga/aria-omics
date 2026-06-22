"""S9 guards: ranking + SPECULATIVE report section + active non-promotion wall.

The visible face of the tier. Ranking orders rival hypotheses by independent
converging evidence (never narrative elegance); the section gathers audited
evidence from the applicable adapters, runs the agent, and renders an explicitly
walled section; the enforcer is now active in the real claim-compilation path.
"""

from __future__ import annotations

import json

import pytest

from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    LLMProposer,
    SpeculativePromotionError,
    assert_no_speculative_promotion,
    build_speculative_section,
    gather_evidence,
    rank_hypotheses,
    render_speculative_section_html,
)


def _exp() -> DiscriminatingExperiment:
    return DiscriminatingExperiment(
        perturbation="KD", readout="qPCR",
        predicted_direction="decrease", refuting_outcome="no change",
    )


# ── ranking ─────────────────────────────────────────────────────────────────

def _signals() -> list[EvidenceSignal]:
    return [
        EvidenceSignal(entity="GATA1", entity_kind="gene", modality="scRNA",
                       measure="log2fc",
                       audited_node_ref="ledger://scRNA/pseudobulk_de",
                       value=2.0),
        EvidenceSignal(entity="KLF1", entity_kind="tf_motif", modality="scATAC",
                       measure="motif_enrich",
                       audited_node_ref="ledger://chromatin/motif_enrichment",
                       value=3.0),
    ]


def test_ranking_prefers_more_independent_lines():
    two_lines = Hypothesis(
        id="two", mechanism="GATA1 and KLF1 may interact",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"], experiment=_exp(),
    )
    one_line = Hypothesis(
        id="one", mechanism="GATA1 may act alone",
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"], experiment=_exp(),
    )
    ranked = rank_hypotheses([one_line, two_lines], _signals())
    assert [h.id for h in ranked] == ["two", "one"]
    assert ranked[0].rank_evidence["n_independent_lines"] == 2


def test_ranking_populates_rank_evidence():
    h = Hypothesis(id="h", mechanism="m", entities=["GATA1"],
                   observation_refs=[], experiment=_exp())
    rank_hypotheses([h], _signals())
    assert h.rank_evidence["mean_abs_effect"] == 2.0


# ── evidence gathering by modality ──────────────────────────────────────────

def test_gather_evidence_detects_present_modalities():
    agent_results = {
        "bulk_rna_agent": {"findings": {"contrasts": [
            {"top_genes": [{"symbol": "MYC", "log2fc": 1.5}]}]}},
    }
    signals = gather_evidence(agent_results)
    assert any(s.entity == "MYC" for s in signals)


def test_gather_evidence_empty_without_agents():
    assert gather_evidence({}) == []


# ── section build + render ──────────────────────────────────────────────────

def _bulk_results() -> dict:
    return {"bulk_rna_agent": {"findings": {"low_power_warning": True, "contrasts": [
        {"top_genes": [{"symbol": "GATA1", "log2fc": 2.3},
                       {"symbol": "KLF1", "log2fc": 1.9}]}]}}}


def _ledger() -> dict:
    return {"entries": [
        {"node_id": "ledger://bulk/differential_expression", "status": "ran"},
        {"node_id": "ledger://bulk/pathway_enrichment", "status": "ran"}]}


def _good_json() -> str:
    return json.dumps([{
        "id": "g1",
        "mechanism": "GATA1 and KLF1 may share an erythroid program",
        "entities": ["GATA1", "KLF1"],
        "observation_refs": ["ledger://bulk/differential_expression"],
        "experiment": {"perturbation": "GATA1 KD", "readout": "KLF1 qPCR",
                       "predicted_direction": "decrease",
                       "refuting_outcome": "no change"},
        "devils_advocate": {"simpler_explanation": "shared stimulus",
                            "confounds": ["low replication"]},
    }])


def test_build_section_returns_none_without_evidence():
    assert build_speculative_section({}, None, {}) is None


def test_build_and_render_section():
    proposer = LLMProposer(lambda p, s: _good_json())
    section = build_speculative_section(
        _bulk_results(), _ledger(), {}, proposer=proposer
    )
    assert section["ran"] is True
    assert section["header"].startswith("Machine-generated hypotheses")
    html = render_speculative_section_html(section)
    assert "SPECULATIVE" in html
    assert "not part of the audited claim manifest" in html
    assert "erythroid program" in html
    assert "Discriminating experiment" in html
    assert "hypothesis://g1" in html


def test_render_empty_for_none_or_not_ran():
    assert render_speculative_section_html(None) == ""
    assert render_speculative_section_html({"ran": False}) == ""


def test_section_honest_null_renders_note():
    # Default (null) proposer -> no hypotheses -> honest-null note, no crash.
    section = build_speculative_section(_bulk_results(), _ledger(), {})
    html = render_speculative_section_html(section)
    assert "honest-null" in html


def test_section_renders_generation_failure_note_not_honest_null():
    # C: a truncated/unparseable proposer response must render as a generation
    # issue, never as "no defensible hypothesis (honest-null)".
    truncated = '```json\n[\n  {\n    "id": "v0",\n    "mechanism": "co-regulation of the markers may'
    section = build_speculative_section(
        _bulk_results(), _ledger(), {},
        proposer=LLMProposer(lambda p, s: truncated),
    )
    assert section["null_reason"] == "parse_error"
    html = render_speculative_section_html(section)
    assert "could not be parsed" in html
    assert "honest-null" not in html


# ── active non-promotion wall ───────────────────────────────────────────────

def test_enforcer_is_a_noop_on_clean_claims():
    assert_no_speculative_promotion([
        {"claim_id": "c1", "tier": "associative",
         "node_id": "ledger://scRNA/pseudobulk_de"},
    ])


def test_enforcer_raises_on_speculative_leak():
    with pytest.raises(SpeculativePromotionError):
        assert_no_speculative_promotion([
            {"claim_id": "x", "tier": "SPECULATIVE"},
        ])
