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
            "causal_inflation", "missing_caveat"} <= labels


def test_b2_corpus_spans_multiple_modalities_and_is_sized():
    from aria.benchmarks.governance_b2 import build_claim_corpus

    corpus = build_claim_corpus()
    assert len(corpus) >= 48
    assert len({b.modality for b in corpus}) >= 3
    assert len({b.analysis for b in corpus}) >= 5


def test_b2_four_arm_ablation_is_coherent():
    from aria.benchmarks.governance_b2 import score_b2, ARMS

    m = score_b2()
    a = m["arms"]
    assert set(a) == set(ARMS)
    # false-narrative is monotone: full system never worse than its ablation,
    # and structure already beats the ungoverned LLM.
    assert (a["governed"]["false_narrative_rate"]
            <= a["guards_off"]["false_narrative_rate"]
            <= a["naive_llm"]["false_narrative_rate"])
    # the guards are what kill causal overreach (structure alone does not).
    assert a["naive_llm"]["causal_overreach_rate"] == 1.0
    assert a["guards_off"]["causal_overreach_rate"] == 1.0
    assert a["governed"]["causal_overreach_rate"] == 0.0
    assert a["template_only"]["causal_overreach_rate"] == 0.0
    # template-only is safe but rigid; governance preserves real interpretation.
    assert a["template_only"]["false_narrative_rate"] == 0.0
    assert (a["governed"]["informative_coverage"]
            > a["template_only"]["informative_coverage"])
    assert m["axis_pass"]["preserves_informative_coverage"] is True


def test_b2_labeling_rows_do_not_leak_gold():
    from aria.benchmarks.governance_b2 import corpus_rows_for_labeling

    rows = corpus_rows_for_labeling()
    assert len(rows) >= 48
    for r in rows:
        assert set(r) == {"claim_id", "modality", "analysis",
                          "claim_text", "evidence_summary"}
        assert "label" not in r  # gold is never exported to the blind sheet
