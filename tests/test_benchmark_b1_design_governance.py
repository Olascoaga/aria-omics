"""Guards for Benchmark B1 (adversarial design governance, governance_b1.py).

Pure-Python (pandas/numpy only): drives ARIA's real readiness audit +
design-matrix validator. Runs in every CI lane.
"""

from __future__ import annotations

import pytest


def test_corpus_is_labelled_and_covers_categories():
    pytest.importorskip("pandas")
    from aria.benchmarks.governance_b1 import build_corpus

    corpus = build_corpus()
    assert len(corpus) >= 20
    cats = {c.category for c in corpus}
    assert {"defensible", "under_replicated", "batch_confounded",
            "continuous_as_categorical", "degenerate"} <= cats
    assert all(c.gold in ("infer", "escalate", "block") for c in corpus)


def test_per_category_decisions_are_correct():
    pytest.importorskip("pandas")
    from aria.benchmarks.governance_b1 import build_corpus, aria_decision

    for case in build_corpus():
        res = aria_decision(case)
        if case.category in ("under_replicated", "batch_confounded",
                             "continuous_as_categorical", "degenerate"):
            assert res["decision"] == "block", (case.case_id, res)
        elif case.category == "defensible":
            assert res["decision"] == "infer", (case.case_id, res)
        elif case.category == "borderline_two_reps":
            assert res["decision"] == "escalate", (case.case_id, res)


def test_b1_headline_unsafe_execution_is_zero():
    pytest.importorskip("pandas")
    from aria.benchmarks.governance_b1 import score_b1

    m = score_b1()
    s = m["summary"]
    # The headline: ARIA never runs inferential analysis on an indefensible design.
    assert s["unsafe_execution_rate"] == 0.0, m["confusion_matrix"]
    assert s["correct_refusal_rate"] >= 0.9
    assert s["correct_inference_rate"] >= 0.9
    assert m["status"] == "pass"
    # No defensible design is wrongly blocked (governance is not trivially "block all").
    assert m["confusion_matrix"]["infer"]["block"] == 0
