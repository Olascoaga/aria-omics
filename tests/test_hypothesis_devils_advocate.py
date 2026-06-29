"""S3 guards for the adversarial devils_advocate gate (ADR-057 rail #5).

Even a speculation must declare its alternatives: a hypothesis must offer a
simpler/competing explanation and must own every confound the audited evidence
it cites already flags (carried in EvidenceSignal.caveats_inherited). Hiding a
known confound, or offering no alternative, is REJECTED — never caveated.
"""

from __future__ import annotations

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    build_evidence_index,
    build_signals_by_entity,
    check_devils_advocate,
    visible_confounds,
)


def test_visible_confounds_unions_across_contexts():
    # H4: GATA1 is measured in two contexts; only ONE flags batch. The gate must
    # still see batch — confounds are unioned across all of the entity's signals,
    # so a confound flagged in any context cannot be hidden by another.
    sigs = [
        EvidenceSignal(entity="GATA1", entity_kind="gene", modality="scRNA",
                       measure="log2fc",
                       audited_node_ref="ledger://scRNA/pseudobulk_de",
                       context="A", caveats_inherited=[]),
        EvidenceSignal(entity="GATA1", entity_kind="gene", modality="scRNA",
                       measure="log2fc",
                       audited_node_ref="ledger://scRNA/pseudobulk_de",
                       context="B",
                       caveats_inherited=["residual batch effect (integration)"]),
    ]
    hyp = Hypothesis(id="h", mechanism="m", entities=["GATA1"])
    assert "batch" in visible_confounds(hyp, build_signals_by_entity(sigs))


def _experiment() -> DiscriminatingExperiment:
    return DiscriminatingExperiment(
        perturbation="CRISPRi knockdown of KLF1",
        readout="GATA1 expression by RT-qPCR",
        predicted_direction="decrease",
        refuting_outcome="GATA1 expression unchanged",
    )


def _clean_signals() -> list[EvidenceSignal]:
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


def _batch_confounded_signals() -> list[EvidenceSignal]:
    # The audited DE result for GATA1 already flagged a batch confound.
    return [
        EvidenceSignal(
            entity="GATA1", entity_kind="gene", modality="scRNA",
            measure="LFC", audited_node_ref="ledger://scRNA/pseudobulk_de",
            caveats_inherited=[
                "condition is partially confounded with batch",
            ],
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


def _hyp(**overrides) -> Hypothesis:
    base = dict(
        id="h",
        mechanism="KLF1 accessibility may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=_experiment(),
        devils_advocate={
            "simpler_explanation": "shared upstream regulator",
            "confounds": [],
        },
    )
    base.update(overrides)
    return Hypothesis(**base)


# ── unit: the gate ──────────────────────────────────────────────────────────

def test_passes_with_simpler_explanation_and_no_visible_confounds():
    index = build_evidence_index(_clean_signals())
    assert check_devils_advocate(_hyp(), index).passed is True


def test_rejects_missing_simpler_explanation():
    index = build_evidence_index(_clean_signals())
    hyp = _hyp(devils_advocate={"confounds": []})
    result = check_devils_advocate(hyp, index)
    assert result.passed is False
    assert "simpler" in (result.reason or "")


def test_rejects_unacknowledged_visible_confound():
    index = build_evidence_index(_batch_confounded_signals())
    # The hypothesis cites GATA1 (batch-confounded) but does not own "batch".
    result = check_devils_advocate(_hyp(), index)
    assert result.passed is False
    assert "batch" in (result.reason or "")


def test_passes_when_visible_confound_is_acknowledged():
    index = build_evidence_index(_batch_confounded_signals())
    hyp = _hyp(
        devils_advocate={
            "simpler_explanation": "shared upstream regulator",
            "confounds": ["batch effect not fully removed"],
        }
    )
    assert check_devils_advocate(hyp, index).passed is True


def test_visible_confounds_reads_inherited_caveats():
    index = build_evidence_index(_batch_confounded_signals())
    assert "batch" in visible_confounds(_hyp(), index)
    assert visible_confounds(_hyp(), build_evidence_index(_clean_signals())) == set()


# ── agent wiring ────────────────────────────────────────────────────────────

def test_agent_rejects_hypothesis_hiding_a_known_confound():
    agent = HypothesisAgent(proposer=lambda s, c: [_hyp(id="hide")])
    out = agent.generate(_batch_confounded_signals(), _ledger())
    assert out["hypotheses"] == []
    assert out["honest_null"] is True
    gates = {f["gate"] for f in out["rejected"][0]["failures"]}
    assert "devils_advocate" in gates
    assert out["null_summary"].get("devils_advocate") == 1


def test_agent_accepts_hypothesis_that_owns_the_confound():
    good = _hyp(
        id="owns",
        devils_advocate={
            "simpler_explanation": "shared upstream regulator",
            "confounds": ["residual batch effect"],
        },
    )
    agent = HypothesisAgent(proposer=lambda s, c: [good])
    out = agent.generate(_batch_confounded_signals(), _ledger())
    assert [h["id"] for h in out["hypotheses"]] == ["owns"]
    assert out["honest_null"] is False
