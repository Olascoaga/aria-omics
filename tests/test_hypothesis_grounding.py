"""S1 guards for the HypothesisAgent grounding wall (ADR-057).

The core invariant: a hypothesis may be free over the *connection* it proposes,
but every *fact* it names must resolve to real audited evidence. A hypothesis
citing an entity absent from the audited evidence — or arising from an analysis
the run did not execute — is REJECTED, never caveated into the output.
"""

from __future__ import annotations

import pytest

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    build_evidence_index,
    build_signals_by_entity,
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


def _observed(*entities) -> list[dict]:
    """H15: faithful observed_claims citing the real signal_id + audited direction."""
    by_ent = {s.entity.lower(): s for s in _signals()}
    return [
        {
            "signal_id": by_ent[e.lower()].signal_id,
            "stated_direction": by_ent[e.lower()].direction,
        }
        for e in entities
    ]


def test_grounded_hypothesis_passes():
    hyp = Hypothesis(
        id="h1",
        mechanism="KLF1 motif accessibility may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=[
            "ledger://scRNA/pseudobulk_de",
            "ledger://chromatin/motif_enrichment",
        ],
        observed_claims=_observed("GATA1", "KLF1"),
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


def test_mechanism_prose_naming_undeclared_entity_is_rejected():
    # The structured entities are all grounded and the cited observation ran, but
    # the mechanism PROSE smuggles in SPI1 — never measured, never declared. The
    # wall must guard the rendered prose, not only the structured entities field;
    # otherwise an LLM evades grounding by naming an invented entity in the text
    # the reader actually sees.
    hyp = Hypothesis(
        id="h_prose",
        mechanism="GATA1 accessibility may be co-opted by SPI1 at shared loci",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert "SPI1" in result.ungrounded_prose_entities
    # The declared entities are all grounded -> the structured-only check is clean.
    assert result.missing_entities == []
    assert "SPI1" in (result.reason or "")


def test_experiment_smuggling_undeclared_entity_is_rejected():
    # H1 bug 2: the mechanism is clean and the entities are grounded, but the
    # discriminating experiment perturbs TP53 — never measured, never declared.
    # The wall must scan EVERY generated field, not only the mechanism.
    hyp = Hypothesis(
        id="h_exp",
        mechanism="GATA1 accessibility may sustain the erythroid program",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=DiscriminatingExperiment(
            perturbation="TP53 knockout",
            readout="GATA1 expression by qPCR",
            predicted_direction="down",
            refuting_outcome="GATA1 unchanged",
        ),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert "TP53" in result.ungrounded_prose_entities
    assert result.missing_entities == []


def test_vacuous_hypothesis_naming_nothing_is_rejected():
    # H1 bug 1: no declared entities and no cited observation -> anchored to
    # nothing. "Grounded by naming nothing" must not be accepted.
    hyp = Hypothesis(
        id="h_void",
        mechanism="the observed shift may reflect a regulatory rewiring",
        entities=[],
        observation_refs=[],
        experiment=DiscriminatingExperiment(
            "perturb the system", "measure a readout", "up", "no change"
        ),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert result.vacuous is True
    assert "vacuous" in (result.reason or "")


def test_hypothesis_with_entities_but_no_observation_is_vacuous():
    # Grounded entities but zero cited observations is still anchored to no
    # audited result.
    hyp = Hypothesis(
        id="h_noref",
        mechanism="KLF1 may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=[],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert result.vacuous is True


def test_mechanism_prose_with_only_grounded_entities_passes():
    # A mechanism that names only measured entities (and hedge/connective words)
    # must NOT be flagged: the prose check targets invented entities, not prose.
    hyp = Hypothesis(
        id="h_ok_prose",
        mechanism="KLF1 motif accessibility may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=_observed("GATA1"),
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is True
    assert result.ungrounded_prose_entities == []


def test_citing_any_audited_node_is_allowed_but_fabricated_node_is_rejected():
    # Round-4 (softens H10): citing an audited node that produced OTHER entities
    # (the pathway/motif node alongside the DE node) is legitimate provenance, not
    # misattribution — a hypothesis reasonably cites the DE node for its genes AND
    # the enriched-pathway node for context. What STAYS rejected is a ref that
    # points at NO audited node in the run (a fabricated citation).
    ok = Hypothesis(
        id="ok", mechanism="GATA1 may act", entities=["GATA1"],
        observation_refs=[
            "ledger://scRNA/pseudobulk_de",       # GATA1's node
            "ledger://chromatin/motif_enrichment",  # an audited node (KLF1's)
        ],
        observed_claims=_observed("GATA1"), experiment=_experiment(),
    )
    r_ok = verify_hypothesis_grounding(ok, _signals(), _ran_ledger())
    assert r_ok.misattributed_refs == []  # both are audited nodes in the run

    fabricated = Hypothesis(
        id="fab", mechanism="GATA1 may act", entities=["GATA1"],
        observation_refs=["ledger://bulk/made_up_node"],
        observed_claims=_observed("GATA1"), experiment=_experiment(),
    )
    r_fab = verify_hypothesis_grounding(fabricated, _signals(), _ran_ledger())
    assert r_fab.grounded is False
    assert "ledger://bulk/made_up_node" in r_fab.misattributed_refs


def test_readout_smuggling_an_entity_is_rejected():
    # H11 (F2): the experiment readout is scanned too — TP53 in the readout, with
    # TP53 absent from the evidence, is a smuggled fact (reverses the H1 carve-out).
    hyp = Hypothesis(
        id="ro", mechanism="GATA1 may act", entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=DiscriminatingExperiment(
            "GATA1 knockdown", "TP53 protein abundance by Western blot",
            "down", "no change"),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert "TP53" in result.ungrounded_prose_entities


def test_readout_assay_vocabulary_is_not_flagged():
    # H11: honest assay descriptions in the readout (RT-qPCR / FACS / GFP / DAPI)
    # must NOT be mistaken for smuggled entities.
    hyp = Hypothesis(
        id="ra", mechanism="GATA1 may act", entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=_observed("GATA1"),
        experiment=DiscriminatingExperiment(
            "GATA1 knockdown", "GATA1 by RT-qPCR, FACS and GFP reporter with DAPI",
            "down", "no change"),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is True
    assert result.ungrounded_prose_entities == []


def test_citation_to_the_producing_node_is_accepted():
    # The same entity cited to the node that actually produced it is fine.
    hyp = Hypothesis(
        id="ok",
        mechanism="GATA1 may act",
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=_observed("GATA1"),
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is True
    assert result.misattributed_refs == []


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
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=_observed("GATA1"),
        experiment=_experiment(),
    )
    # No run_ledger -> the ledger-node check is skipped; entity grounding,
    # prose grounding and non-vacuity are still enforced.
    result = verify_hypothesis_grounding(hyp, _signals(), None)
    assert result.grounded is True


def _ctx_signal(entity, context, value, direction):
    return EvidenceSignal(
        entity=entity, entity_kind="gene", modality="bulk_RNA",
        measure="log2fc", audited_node_ref="ledger://bulk/differential_expression",
        value=value, direction=direction, context=context,
    )


def test_signal_id_distinguishes_context():
    # H4: the SAME entity in two contexts is two distinct signals.
    a = _ctx_signal("GATA1", "old_vs_young", 2.0, "up")
    b = _ctx_signal("GATA1", "treated_vs_ctrl", -1.0, "down")
    assert a.signal_id and b.signal_id
    assert a.signal_id != b.signal_id


def test_build_signals_by_entity_keeps_every_context():
    a = _ctx_signal("GATA1", "old_vs_young", 2.0, "up")
    b = _ctx_signal("GATA1", "treated_vs_ctrl", -1.0, "down")
    by_entity = build_signals_by_entity([a, b])
    assert len(by_entity["gata1"]) == 2
    # build_evidence_index still yields one deterministic representative.
    assert build_evidence_index([a, b])["gata1"] is a


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
        mechanism="KLF1 motif accessibility may sustain GATA1 expression",
        entities=["GATA1", "KLF1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=_observed("GATA1"),
        experiment=_experiment(),
        devils_advocate={
            "simpler_explanation": "co-regulation by a shared upstream factor",
            "confounds": [],
        },
    )
    invented = Hypothesis(
        id="bad",
        mechanism="MYB may explain the shift",
        entities=["MYB"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        experiment=_experiment(),
    )

    def proposer(signals, exp_ctx):
        return [grounded, invented]

    agent = HypothesisAgent(proposer=proposer)
    out = agent.generate(_signals(), _ran_ledger(), w_claim_passed=True, w_ledger_passed=True)
    assert out["ran"] is True
    assert out["requires_ack"] is True
    assert [h["id"] for h in out["hypotheses"]] == ["ok"]
    bad = out["rejected"][0]
    assert bad["hypothesis_id"] == "bad"
    grounding_fail = next(
        f for f in bad["failures"] if f["gate"] == "grounding"
    )
    assert "MYB" in grounding_fail["missing_entities"]
    assert out["honest_null"] is False


def test_agent_honest_null_with_default_proposer():
    agent = HypothesisAgent()
    out = agent.generate(_signals(), _ran_ledger(), w_claim_passed=True, w_ledger_passed=True)
    assert out["ran"] is True
    assert out["hypotheses"] == []
    assert out["honest_null"] is True


def test_generate_requires_explicit_verification_flags():
    # H13 (F5): the W-CLAIM/W-LEDGER flags are fail-closed — a caller MUST pass
    # them explicitly. There is no permissive default for a future caller to
    # forget, so the gate cannot be silently bypassed.
    agent = HypothesisAgent()
    with pytest.raises(TypeError):
        agent.generate(_signals(), _ran_ledger())


def test_causal_gate_blocks_when_verification_failed():
    def proposer(signals, exp_ctx):
        raise AssertionError("proposer must not run when the gate is closed")

    agent = HypothesisAgent(proposer=proposer)
    out = agent.generate(_signals(), _ran_ledger(), w_claim_passed=True, w_ledger_passed=False)
    assert out["ran"] is False
    assert out["reason"] == "verification_gate_not_passed"
    assert out["hypotheses"] == []


def test_agent_does_not_mutate_inputs():
    signals = _signals()
    ledger = _ran_ledger()
    n_signals = len(signals)
    n_entries = len(ledger["entries"])
    HypothesisAgent().generate(signals, ledger, w_claim_passed=True, w_ledger_passed=True)
    assert len(signals) == n_signals
    assert len(ledger["entries"]) == n_entries


# ── H15 (round-3, Codex blocker 2): signal-level direction/context grounding ──

def test_observed_claim_contradicting_audited_direction_is_rejected():
    # GATA1 is audited UP; a hypothesis that CITES that exact signal but reads it
    # as "increased"... is faithful. The contradiction is reading it "decreased".
    sigs = _signals()
    gata1 = next(s for s in sigs if s.entity == "GATA1")  # direction up
    hyp = Hypothesis(
        id="contra",
        mechanism="GATA1 accessibility may sustain the program",
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=[
            {"signal_id": gata1.signal_id, "stated_direction": "decreased"}
        ],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, sigs, _ran_ledger())
    assert result.grounded is False
    assert result.contradicting_claims
    assert result.contradicting_claims[0]["entity"] == "GATA1"
    assert "contradicts" in (result.reason or "")


def test_observed_claim_unknown_signal_id_is_rejected():
    hyp = Hypothesis(
        id="unk",
        mechanism="GATA1 may act",
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=[
            {"signal_id": "sig_does_not_exist", "stated_direction": "up"}
        ],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert "sig_does_not_exist" in result.unknown_signals


def test_observed_claim_for_a_real_signal_is_allowed_even_if_entity_unlisted():
    # Round-4 (softens H15's misattributed_signals): the evidence uses official
    # symbols while a model cites in common names, so requiring the cited signal's
    # gene to be in `entities` rejected FAITHFUL citations of real measurements
    # (NR1D1 cited while the hypothesis names REV-ERBα). A cited REAL signal is
    # evidence — invention stays blocked because the signal must EXIST and its
    # direction must not be contradicted.
    sigs = _signals()
    klf1 = next(s for s in sigs if s.entity == "KLF1")
    hyp = Hypothesis(
        id="cite", mechanism="GATA1 may act with a partner", entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=[
            {"signal_id": klf1.signal_id, "stated_direction": klf1.direction or "na"}
        ],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, sigs, _ran_ledger())
    assert result.misattributed_signals == []  # no longer rejected
    assert result.grounded is True

    # But citing a NON-existent signal is still rejected as invention.
    bad = Hypothesis(
        id="badcite", mechanism="GATA1 may act", entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=[{"signal_id": "sig_does_not_exist", "stated_direction": "up"}],
        experiment=_experiment(),
    )
    rbad = verify_hypothesis_grounding(bad, sigs, _ran_ledger())
    assert rbad.grounded is False
    assert "sig_does_not_exist" in rbad.unknown_signals


def _design_signals() -> list[EvidenceSignal]:
    """Evidence from a real KO study: official symbols + contrast-label contexts."""
    return [
        EvidenceSignal(
            entity="NR1D1", entity_kind="gene", modality="bulk_RNA",
            measure="log2fc", audited_node_ref="ledger://bulk/differential_expression",
            value=-2.1, direction="down", context="rev_erba_ko_vs_wt",
        ),
        EvidenceSignal(
            entity="TP53", entity_kind="gene", modality="bulk_RNA",
            measure="log2fc", audited_node_ref="ledger://bulk/differential_expression",
            value=1.5, direction="up", context="rev_erba_ko_vs_wt",
        ),
        EvidenceSignal(
            entity="apoptosis", entity_kind="pathway", modality="bulk_RNA",
            measure="ora", audited_node_ref="ledger://bulk/pathway_enrichment",
            value=None, direction="na", context="rev_erba_ko_vs_wt",
        ),
    ]


def _design_ledger() -> dict:
    return {"entries": [
        {"node_id": "ledger://bulk/differential_expression", "status": "ran"},
        {"node_id": "ledger://bulk/pathway_enrichment", "status": "ran"},
    ]}


def test_common_name_of_a_perturbation_target_grounds_via_design():
    # Round-4 (calibration): the model names the KO'd gene by its common name
    # (REV-ERBα = NR1D1, declared by the rev_erba_ko contrast). It grounds via the
    # design target, cites NR1D1's real signal, mentions TP53 (measured), a residue
    # (Ser15, dropped) and the pathway node. A run that used to emit NOTHING emits.
    sigs = _design_signals()
    nr1d1 = next(s for s in sigs if s.entity == "NR1D1")
    hyp = Hypothesis(
        id="h4",
        mechanism="REV-ERBα loss may de-repress TP53 and promote apoptosis via Ser15",
        entities=["REV-ERBα"],
        observation_refs=[
            "ledger://bulk/pathway_enrichment",
            "ledger://bulk/differential_expression",
        ],
        observed_claims=[{"signal_id": nr1d1.signal_id, "stated_direction": "down"}],
        experiment=DiscriminatingExperiment(
            "REV-ERBα knockdown", "TP53 target expression by qPCR",
            "increase", "no change"),
        devils_advocate={"simpler_explanation": "shared stress response",
                         "confounds": []},
    )
    result = verify_hypothesis_grounding(hyp, sigs, _design_ledger())
    assert result.grounded is True, result.reason


def test_assay_and_concept_prose_tokens_are_not_flagged_as_entities():
    # Round-4: ICP-MS / ATAC / ECM / SASP are assays/concepts, not genes; a residue
    # (Ser15) and a labelled condition (BMAL1_KO -> the design target) must not read
    # as ungrounded entities.
    sigs = _design_signals()
    nr1d1 = next(s for s in sigs if s.entity == "NR1D1")
    hyp = Hypothesis(
        id="prose",
        mechanism="NR1D1 loss may reshape the ECM and SASP, probed by ATAC and ICP-MS",
        entities=["NR1D1"],
        observation_refs=["ledger://bulk/differential_expression"],
        observed_claims=[{"signal_id": nr1d1.signal_id, "stated_direction": "down"}],
        experiment=DiscriminatingExperiment(
            "NR1D1 knockdown", "chromatin accessibility by ATAC-seq",
            "increase", "no change"),
        devils_advocate={"simpler_explanation": "batch", "confounds": []},
    )
    result = verify_hypothesis_grounding(hyp, sigs, _design_ledger())
    assert result.ungrounded_prose_entities == []
    assert result.grounded is True, result.reason


def test_invented_gene_still_rejected_after_calibration():
    # The calibration must NOT reopen the invention hole: a gene that is neither
    # measured nor a design target is still rejected.
    sigs = _design_signals()
    nr1d1 = next(s for s in sigs if s.entity == "NR1D1")
    hyp = Hypothesis(
        id="inv", mechanism="SOX2 may drive the program", entities=["SOX2"],
        observation_refs=["ledger://bulk/differential_expression"],
        observed_claims=[{"signal_id": nr1d1.signal_id, "stated_direction": "down"}],
        experiment=_experiment(),
        devils_advocate={"simpler_explanation": "x", "confounds": []},
    )
    result = verify_hypothesis_grounding(hyp, sigs, _design_ledger())
    assert result.grounded is False
    assert "SOX2" in result.missing_entities


def test_missing_observed_claims_is_rejected():
    # Entities + refs are fine, but the hypothesis declares no signal it reads.
    hyp = Hypothesis(
        id="no_obs",
        mechanism="GATA1 may act",
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=[],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, _signals(), _ran_ledger())
    assert result.grounded is False
    assert result.missing_observed_claims is True


def test_downstream_speculation_opposite_direction_is_allowed():
    # The OBSERVED claim is faithful (GATA1 up); the mechanism freely speculates a
    # DOWNSTREAM repression in the opposite direction. That speculation must NOT be
    # rejected — only a false restatement of the audited signal is (no over-reject).
    sigs = _signals()
    gata1 = next(s for s in sigs if s.entity == "GATA1")  # up
    hyp = Hypothesis(
        id="downstream",
        mechanism=(
            "elevated GATA1 may in turn repress its downstream targets, "
            "lowering their accessibility"
        ),
        entities=["GATA1"],
        observation_refs=["ledger://scRNA/pseudobulk_de"],
        observed_claims=[
            {"signal_id": gata1.signal_id, "stated_direction": "up"}
        ],
        experiment=_experiment(),
    )
    result = verify_hypothesis_grounding(hyp, sigs, _ran_ledger())
    assert result.grounded is True
    assert result.contradicting_claims == []
