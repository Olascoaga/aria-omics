"""Guards for Benchmark B2 (claim/narrative governance, governance_b2.py)."""

from __future__ import annotations


def test_b2_governance_reduces_false_narrative_and_causal_overreach():
    from aria.benchmarks.governance_b2 import score_b2

    m = score_b2()
    s = m["summary"]
    # ARIA governance reduces false narrative vs the ungoverned arm.
    assert s["governed_false_narrative_rate"] < s["naive_false_narrative_rate"]
    # ARIA never licenses causal/mechanistic conclusions from observational omics.
    assert s["governed_causal_overreach_rate"] == 0.0
    assert s["naive_causal_overreach_rate"] == 1.0
    assert m["status"] == "pass"


def test_b2_corpus_is_labelled_with_failure_taxonomy():
    from aria.benchmarks.governance_b2 import build_claim_corpus

    labels = {b.metadata["b2_label"] for b in build_claim_corpus()}
    assert {"clean", "unsupported", "fabricated", "overclaim",
            "causal_inflation"} <= labels
