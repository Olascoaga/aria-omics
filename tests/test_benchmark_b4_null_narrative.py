"""Guards for Benchmark B4 (null-narrative governance, governance_b4.py)."""

from __future__ import annotations

import pytest


def test_null_corpus_has_honest_and_fabricated():
    from aria.benchmarks.governance_b4 import build_null_corpus

    corpus = build_null_corpus()
    kinds = {b.metadata.get("b4_kind") for b in corpus}
    assert kinds == {"honest", "fabricated"}
    assert sum(b.metadata["b4_kind"] == "fabricated" for b in corpus) >= 5


def test_no_fabricated_biology_slips_through_on_null():
    from aria.benchmarks.governance_b4 import score_b4

    m = score_b4()
    s = m["summary"]
    # Headline: no fabricated biology becomes a licensed conclusion on null data.
    assert s["fabricated_narrative_rate"] == 0.0, m["per_block"]
    assert s["fabricated_caught_rate"] == 1.0
    # Governance is not trivially "reject everything": honest "no signal" passes.
    assert s["honest_false_flag_rate"] <= 0.25
    assert m["status"] == "pass"
    # Every fabricated block is either hard-rejected or tier-capped to descriptive.
    assert s["fabricated_hard_rejected"] + s["fabricated_tier_capped_descriptive"] == s["n_fabricated"]
