"""S2 guards for the SPECULATIVE publication gates (ADR-057 rails #8 + #11).

A hypothesis must EARN publication: a complete, concrete discriminating
experiment (falsifiability) and a hedged, non-assertive register (language lint).
A candidate failing either is REJECTED, never caveated; when all candidates are
rejected the agent reports an honest-null with aggregated per-gate reasons.
"""

from __future__ import annotations

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    check_falsifiability,
    check_language,
)


def _signals() -> list[EvidenceSignal]:
    return [
        EvidenceSignal(
            entity="GATA1", entity_kind="gene", modality="scRNA",
            measure="LFC", audited_node_ref="ledger://scRNA/pseudobulk_de",
        ),
        EvidenceSignal(
            entity="KLF1", entity_kind="tf_motif", modality="scATAC",
            measure="motif_enrich",
            audited_node_ref="ledger://chromatin/motif_enrichment",
        ),
    ]


def _ledger() -> dict:
    return {
        "entries": [
            {"node_id": "ledger://scRNA/pseudobulk_de", "status": "ran"},
            {"node_id": "ledger://chromatin/motif_enrichment", "status": "ran"},
        ]
    }


def _complete_experiment() -> DiscriminatingExperiment:
    return DiscriminatingExperiment(
        perturbation="CRISPRi knockdown of KLF1",
        readout="GATA1 expression by RT-qPCR",
        predicted_direction="decrease",
        refuting_outcome="GATA1 expression unchanged after KLF1 knockdown",
    )


def _observed(*entities) -> list[dict]:
    """H15: faithful observed_claims citing real signal_ids."""
    by_ent = {s.entity.lower(): s for s in _signals()}
    return [
        {
            "signal_id": by_ent[e.lower()].signal_id,
            "stated_direction": by_ent[e.lower()].direction or "na",
        }
        for e in entities
    ]


def _valid_hypothesis(**overrides) -> Hypothesis:
    base = dict(
        id="h",
        mechanism="KLF1 accessibility may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=_observed("GATA1"),
        experiment=_complete_experiment(),
        devils_advocate={
            "simpler_explanation": "co-regulation by a shared upstream factor",
            "confounds": [],
        },
    )
    base.update(overrides)
    return Hypothesis(**base)


# ── H7: word-boundary lint hardening ────────────────────────────────────────

def test_hedge_marker_matches_before_punctuation():
    # H7: the old "may " (literal trailing space) missed a hedge before
    # punctuation; word boundaries fix it.
    assert check_language(
        _valid_hypothesis(mechanism="KLF1 may sustain GATA1.")
    ).passed is True
    assert check_language(
        _valid_hypothesis(mechanism="this may, in turn, sustain GATA1")
    ).passed is True


def test_hedge_word_not_matched_inside_a_larger_word():
    # H7: "may" inside "mayonnaise" is not a hedge — needs a real marker.
    assert check_language(
        _valid_hypothesis(mechanism="the mayonnaise mapping shifts GATA1")
    ).passed is False


def test_direction_token_not_matched_inside_a_word():
    # H7: "up" must NOT satisfy the direction check inside "upstream".
    hyp = _valid_hypothesis(experiment=DiscriminatingExperiment(
        perturbation="KLF1 knockdown", readout="GATA1 by qPCR",
        predicted_direction="altered upstream regulation",
        refuting_outcome="no effect"))
    assert check_falsifiability(hyp).passed is False


def test_stem_marker_still_matches_inflections():
    # H7: stem markers keep catching inflections (suggests/proposed).
    assert check_language(
        _valid_hypothesis(mechanism="the data suggests a shared GATA1 program")
    ).passed is True


# ── falsifiability ──────────────────────────────────────────────────────────

def test_falsifiability_passes_for_complete_experiment():
    assert check_falsifiability(_valid_hypothesis()).passed is True


def test_falsifiability_rejects_incomplete_experiment():
    hyp = _valid_hypothesis(
        experiment=DiscriminatingExperiment(
            perturbation="KLF1 knockdown", readout="", predicted_direction="",
            refuting_outcome="",
        )
    )
    result = check_falsifiability(hyp)
    assert result.passed is False
    assert "incomplete" in (result.reason or "")


def test_falsifiability_rejects_vacuous_experiment():
    hyp = _valid_hypothesis(
        experiment=DiscriminatingExperiment(
            perturbation="further studies are needed",
            readout="more experiments",
            predicted_direction="decrease",
            refuting_outcome="unclear",
        )
    )
    result = check_falsifiability(hyp)
    assert result.passed is False
    assert "non-specific" in (result.reason or "")


def test_falsifiability_rejects_directionless_prediction():
    hyp = _valid_hypothesis(
        experiment=DiscriminatingExperiment(
            perturbation="KLF1 knockdown",
            readout="GATA1 qPCR",
            predicted_direction="something happens",
            refuting_outcome="nothing happens",
        )
    )
    result = check_falsifiability(hyp)
    assert result.passed is False
    assert "direction" in (result.reason or "")


# ── language lint ───────────────────────────────────────────────────────────

def test_language_passes_for_hedged_mechanism():
    assert check_language(_valid_hypothesis()).passed is True


def test_language_rejects_causal_verb():
    hyp = _valid_hypothesis(mechanism="KLF1 accessibility drives GATA1 expression")
    result = check_language(hyp)
    assert result.passed is False
    assert "assertive/causal" in (result.reason or "")


def test_language_rejects_assertive_finding():
    hyp = _valid_hypothesis(
        mechanism="we demonstrate KLF1 regulates GATA1"
    )
    assert check_language(hyp).passed is False


def test_language_rejects_unhedged_neutral_mechanism():
    hyp = _valid_hypothesis(mechanism="KLF1 accessibility and GATA1 expression")
    result = check_language(hyp)
    assert result.passed is False
    assert "speculative marker" in (result.reason or "")


# ── agent wiring + honest-null aggregation ──────────────────────────────────

def test_agent_accepts_fully_valid_hypothesis():
    agent = HypothesisAgent(proposer=lambda s, c: [_valid_hypothesis(id="ok")])
    out = agent.generate(_signals(), _ledger(), w_claim_passed=True, w_ledger_passed=True)
    assert [h["id"] for h in out["hypotheses"]] == ["ok"]
    assert out["honest_null"] is False
    assert "null_summary" not in out


def test_agent_rejects_unfalsifiable_and_reports_failure():
    bad = _valid_hypothesis(
        id="nofalsify",
        experiment=DiscriminatingExperiment(
            perturbation="KLF1 knockdown", readout="", predicted_direction="",
            refuting_outcome="",
        ),
    )
    agent = HypothesisAgent(proposer=lambda s, c: [bad])
    out = agent.generate(_signals(), _ledger(), w_claim_passed=True, w_ledger_passed=True)
    assert out["hypotheses"] == []
    assert out["honest_null"] is True
    gates = {f["gate"] for f in out["rejected"][0]["failures"]}
    assert "falsifiability" in gates
    assert out["null_summary"].get("falsifiability") == 1


def test_agent_honest_null_aggregates_reasons_across_gates():
    causal = _valid_hypothesis(
        id="causal", mechanism="KLF1 causes GATA1 expression"
    )
    invented = _valid_hypothesis(
        id="invented", mechanism="MYB may explain it", entities=["MYB"]
    )
    agent = HypothesisAgent(proposer=lambda s, c: [causal, invented])
    out = agent.generate(_signals(), _ledger(), w_claim_passed=True, w_ledger_passed=True)
    assert out["hypotheses"] == []
    assert out["honest_null"] is True
    summary = out["null_summary"]
    assert summary.get("language") == 1
    assert summary.get("grounding") == 1


def test_a_single_candidate_can_fail_multiple_gates():
    # Ungrounded entity AND assertive language in one shot.
    hyp = _valid_hypothesis(
        id="double",
        mechanism="FOXP3 drives the program",
        entities=["FOXP3"],
    )
    agent = HypothesisAgent(proposer=lambda s, c: [hyp])
    out = agent.generate(_signals(), _ledger(), w_claim_passed=True, w_ledger_passed=True)
    gates = {f["gate"] for f in out["rejected"][0]["failures"]}
    assert {"grounding", "language"} <= gates
