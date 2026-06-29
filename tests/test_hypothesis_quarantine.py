"""S4 guards for the ledger quarantine (ADR-057 rail #6).

A machine hypothesis lives in its own ``hypothesis://`` namespace, linked to the
audited evidence it cites but mechanically marked non-promotable: NOTHING
downstream (the audited claim manifest, RO-Crate export, aria diff) may promote
it to a claim. Zero feedback loop — a speculation never becomes a premise.
"""

from __future__ import annotations

import pytest

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.claim_compiler import TIER_ORDER
from aria.agents.narrative.hypothesis import (
    SPECULATIVE_TIER,
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    SpeculativePromotionError,
    assert_no_speculative_promotion,
    is_hypothesis_node,
    quarantine_hypotheses,
    quarantine_node_id,
)
from aria.agents.narrative.run_ledger import (
    _NON_DESCRIPTIVE_TIERS,
    node_id_for,
)


def _experiment() -> DiscriminatingExperiment:
    return DiscriminatingExperiment(
        perturbation="CRISPRi knockdown of KLF1",
        readout="GATA1 expression by RT-qPCR",
        predicted_direction="decrease",
        refuting_outcome="GATA1 expression unchanged",
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


def _hyp(**overrides) -> Hypothesis:
    base = dict(
        id="h1",
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


# ── namespace + manifest ────────────────────────────────────────────────────

def test_quarantine_node_id_uses_hypothesis_scheme():
    assert quarantine_node_id(_hyp(id="abc")) == "hypothesis://abc"
    assert is_hypothesis_node("hypothesis://abc") is True
    assert is_hypothesis_node("ledger://scRNA/pseudobulk_de") is False


def test_quarantine_stamps_node_and_builds_non_promotable_manifest():
    hyp = _hyp(id="h7")
    manifest = quarantine_hypotheses([hyp])
    assert hyp.ledger_node == "hypothesis://h7"
    entry = manifest[0]
    assert entry["hypothesis_node"] == "hypothesis://h7"
    assert entry["tier"] == SPECULATIVE_TIER
    assert entry["promotable"] is False
    assert "ledger://scRNA/pseudobulk_de" in entry["cites"]


def test_agent_emits_quarantine_for_accepted_hypotheses():
    agent = HypothesisAgent(proposer=lambda s, c: [_hyp(id="ok")])
    out = agent.generate(_signals(), _ledger(), w_claim_passed=True, w_ledger_passed=True)
    assert out["hypotheses"][0]["ledger_node"] == "hypothesis://ok"
    assert out["quarantine"][0]["hypothesis_node"] == "hypothesis://ok"
    assert out["quarantine"][0]["promotable"] is False


# ── the wall is by construction ─────────────────────────────────────────────

def test_speculative_is_not_a_claim_tier():
    assert SPECULATIVE_TIER not in TIER_ORDER
    assert SPECULATIVE_TIER not in _NON_DESCRIPTIVE_TIERS


def test_claim_node_ids_never_collide_with_hypothesis_namespace():
    for modality in ("scRNA", "bulk", "chromatin", "report"):
        nid = node_id_for(modality, "pseudobulk_de")
        assert not is_hypothesis_node(nid)


# ── the active enforcer ─────────────────────────────────────────────────────

def test_enforcer_passes_clean_claim_manifest():
    claims = [
        {"claim_id": "c1", "tier": "associative",
         "node_id": "ledger://scRNA/pseudobulk_de"},
        {"claim_id": "c2", "tier": "descriptive"},
    ]
    assert_no_speculative_promotion(claims)  # no raise


def test_enforcer_rejects_speculative_tier_in_claims():
    claims = [{"claim_id": "leak", "tier": SPECULATIVE_TIER}]
    with pytest.raises(SpeculativePromotionError):
        assert_no_speculative_promotion(claims)


def test_enforcer_rejects_hypothesis_node_in_claims():
    claims = [
        {"claim_id": "leak", "tier": "associative",
         "ledger_node_id": "hypothesis://h1"},
    ]
    with pytest.raises(SpeculativePromotionError):
        assert_no_speculative_promotion(claims)


def test_duplicate_ids_get_distinct_quarantine_nodes():
    # H1 bug 3: two hypotheses with the same id must NOT collapse onto one
    # quarantine node, or one would silently inherit the other's cites/provenance.
    a = _hyp(id="dup", entities=["GATA1"])
    b = _hyp(id="dup", entities=["KLF1"])
    manifest = quarantine_hypotheses([a, b])
    nodes = [m["hypothesis_node"] for m in manifest]
    assert a.ledger_node != b.ledger_node
    assert len(set(nodes)) == len(nodes)
    assert nodes[0] == "hypothesis://dup"


def test_unsafe_id_is_slugged():
    # An LLM-authored id cannot inject scheme/path characters into the namespace.
    hyp = _hyp(id="../evil id://x")
    node = quarantine_node_id(hyp)
    assert node.startswith("hypothesis://")
    assert "://x" not in node[len("hypothesis://"):]
    assert " " not in node


def test_quarantined_hypothesis_cannot_pass_as_a_claim():
    # Simulate a leak: build a "claim" from the quarantine manifest and prove the
    # enforcer catches it before it could be emitted as audited.
    hyp = _hyp(id="sneaky")
    quarantine_hypotheses([hyp])
    forged_claim = {
        "claim_id": "sneaky",
        "tier": SPECULATIVE_TIER,
        "node_id": hyp.ledger_node,
    }
    with pytest.raises(SpeculativePromotionError):
        assert_no_speculative_promotion([forged_claim])
