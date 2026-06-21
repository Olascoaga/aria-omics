"""S1 guards for the HypothesisAgent grounding wall (ADR-057).

The core invariant: a hypothesis may be free over the *connection* it proposes,
but every *fact* it names must resolve to real audited evidence. A hypothesis
citing an entity absent from the audited evidence — or arising from an analysis
the run did not execute — is REJECTED, never caveated into the output.
"""

from __future__ import annotations

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    build_evidence_index,
    verify_hypothesis_grounding,
)


def _signals() -> list[EvidenceSignal]:
    return [
        EvidenceSignal(
            entity="GATA1",
            entity_kind="gene",
            modality="scRNA",
            measure="LFC",
            audited_node_ref="ledger://scRNA/pseudobulk_de",
            value=2.1,
            direction="up",
        ),
        EvidenceSignal(
            entity="KLF1",
            entity_kind="tf_motif",
            modality="scATAC",
            measure="motif_enrich",
            audited_node_ref="ledger://chromatin/motif_enrichment",
            value=3.4,
            direction="up",
        ),
    ]


def _ran_ledger() -> dict:
    return {
        "entries": [
            {
                "node_id": "ledger://scRNA/pseudobulk_de",
                "status": "ran",
                "modality": "scRNA",
                "analysis": "pseudobulk_de",
            },
            {
                "node_id": "ledger://chromatin/motif_enrichment",
                "status": "ran",
                "modality": "chromatin",
                "analysis": "motif_enrichment",
            },
            {
                "node_id": "ledger://chromatin/regulatory_layers",
                "status": "skipped",
                "reason": "no_peak2gene_inputs",
                "modality": "chromatin",
                "analysis": "regulatory_layers",
            },
        ]
    }


def _experiment() -> DiscriminatingExperiment:
    return DiscriminatingExperiment(
        perturbation="KLF1 knockdown",
        readout="GATA1 expression by qPCR",
        predicted_direction="down",
        refuting_outcome="GATA1 unchanged",
    )


def test_grounded_hypothesis_passes():
    hyp = Hypothesis(
        id="h1",
        mechanism="KLF1 motif accessibility may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=[
            "ledger://scRNA/pseudobulk_de",
            "ledger://chromatin/motif_enrichment",
        ],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is True
    assert result.missing_entities == []
    assert result.not_run_refs == []


def test_hypothesis_citing_absent_entity_is_rejected():
    # FOXP3 is never measured in the audited evidence -> invented fact.
    hyp = Hypothesis(
        id="h2",
        mechanism="FOXP3 drives the observed accessibility shift",
        entities=["GATA1", "FOXP3"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert "FOXP3" in result.missing_entities
    assert "GATA1" not in result.missing_entities
    assert "FOXP3" in (result.reason or "")


def test_hypothesis_from_not_run_analysis_is_rejected():
    # Entities are real, but the cited observation analysis was skipped.
    hyp = Hypothesis(
        id="h3",
        mechanism="peak-to-gene linkage explains the co-occurrence",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://chromatin/regulatory_layers"],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert result.missing_entities == []
    assert any(
        n["node_id"] == "ledger://chromatin/regulatory_layers"
        and n["status"] == "skipped"
        for n in result.not_run_refs
    )


def test_observation_ref_with_no_ledger_node_is_rejected():
    hyp = Hypothesis(
        id="h4",
        mechanism="some connection",
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/velocity"],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert any(
        n["status"] == "no_ledger_node" for n in result.not_run_refs
    )


def test_entity_grounding_is_case_insensitive():
    hyp = Hypothesis(
        id="h5",
        mechanism="connection",
        entities=["gata1", "Klf1"],
        observation_refs=[],
        experiment=_experiment(),
    )
    # No run_ledger -> only entity grounding is enforced.
    result = verify_hypothesis_grounding(hyp, _signals(), None)
    assert result.grounded is True


def test_build_evidence_index_skips_blank_and_nonsignals():
    signals = _signals() + [
        EvidenceSignal(
            entity="", entity_kind="gene", modality="scRNA",
            measure="LFC", audited_node_ref="ledger://scRNA/pseudobulk_de",
        ),
        "not-a-signal",  # type: ignore[list-item]
    ]
    index = build_evidence_index(signals)  # type: ignore[arg-type]
    assert set(index) == {"gata1", "klf1"}


def test_agent_accepts_grounded_rejects_ungrounded():
    grounded = Hypothesis(
        id="ok",
        mechanism="KLF1 accessibility sustains GATA1",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=_experiment(),
    )
    invented = Hypothesis(
        id="bad",
        mechanism="MYB invented",
        entities=["MYB"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=_experiment(),
    )

    def proposer(signals, exp_ctx):
        return [grounded, invented]

    agent = HypothesisAgent(proposer=proposer)
    out = agent.generate(_signals(), _ran_ledger())
    assert out["ran"] is True
    assert out["requires_ack"] is True
    assert [h["id"] for h in out["hypotheses"]] == ["ok"]
    assert out["rejected"][0]["hypothesis_id"] == "bad"
    assert "MYB" in out["rejected"][0]["grounding"]["missing_entities"]
    assert out["honest_null"] is False


def test_agent_honest_null_with_default_proposer():
    agent = HypothesisAgent()
    out = agent.generate(_signals(), _ran_ledger())
    assert out["ran"] is True
    assert out["hypotheses"] == []
    assert out["honest_null"] is True


def test_causal_gate_blocks_when_verification_failed():
    def proposer(signals, exp_ctx):
        raise AssertionError("proposer must not run when the gate is closed")

    agent = HypothesisAgent(proposer=proposer)
    out = agent.generate(_signals(), _ran_ledger(), w_ledger_passed=False)
    assert out["ran"] is False
    assert out["reason"] == "verification_gate_not_passed"
    assert out["hypotheses"] == []


def test_agent_does_not_mutate_inputs():
    signals = _signals()
    ledger = _ran_ledger()
    n_signals = len(signals)
    n_entries = len(ledger["entries"])
    HypothesisAgent().generate(signals, ledger)
    assert len(signals) == n_signals
    assert len(ledger["entries"]) == n_entries
