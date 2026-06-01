"""P1-14 — ParameterAdvisor honesty.

Three fixes, each a failing-first guard:
  1. Historical recall is conditioned on `bio_context` (organism), not just an
     `analysis_type` substring — a mouse-brain decision is not recalled for a
     human run.
  2. The historical-approval bonus is BOUNDED, so it can only break a near-tie,
     never override a real objective-score difference.
  3. When clustering metrics cannot be measured (subprocess failure), candidates
     are honest `not measured` (no fabricated silhouette/modularity) — never a
     comparable substitute.
"""

from __future__ import annotations

import pytest

# parameter_advisor -> aria.llm.provider imports litellm eagerly; skip cleanly
# where litellm is absent (the light CI lane) instead of erroring at collection.
pytest.importorskip("litellm")

from aria.llm.parameter_advisor import (   # noqa: E402
    ParameterAdvisor,
    ParameterCandidate,
    ParameterDecision,
)


def _advisor():
    """ParameterAdvisor without __init__ (only memory + class consts needed)."""
    return ParameterAdvisor.__new__(ParameterAdvisor)


class _FakeMemory:
    def __init__(self, wings, decisions_by_wing):
        self._wings = wings
        self._decisions = decisions_by_wing

    def list_wings(self):
        return self._wings

    def get_decisions(self, wing_id):
        return self._decisions.get(wing_id, [])


# ── 1. recall conditioned on organism ────────────────────────────────────────

def test_recall_skips_other_organism():
    mem = _FakeMemory(
        wings=[
            {"id": "w_human", "organism": "Homo sapiens"},
            {"id": "w_mouse", "organism": "Mus musculus"},
        ],
        decisions_by_wing={
            "w_human": [{"question": "leiden_clustering resolution",
                         "decision": "1.0", "made_at": "2026-01-02"}],
            "w_mouse": [{"question": "leiden_clustering resolution",
                         "decision": "0.4", "made_at": "2026-01-03"}],
        },
    )
    adv = _advisor()
    adv.memory = mem
    hist = adv._recall_similar_decisions(
        "exp", "leiden_clustering", {"organism": "Homo sapiens"}
    )
    decisions = {h["decision"] for h in hist}
    assert decisions == {"1.0"}            # mouse decision excluded


def test_recall_unknown_organism_falls_back_to_type_match():
    mem = _FakeMemory(
        wings=[{"id": "w1", "organism": "Mus musculus"}],
        decisions_by_wing={
            "w1": [{"question": "leiden_clustering resolution",
                    "decision": "0.4", "made_at": "2026-01-03"}],
        },
    )
    adv = _advisor()
    adv.memory = mem
    hist = adv._recall_similar_decisions("exp", "leiden_clustering", {})
    assert [h["decision"] for h in hist] == ["0.4"]   # no organism gate applied


# ── 2. bounded historical bonus ──────────────────────────────────────────────

def test_historical_bonus_cannot_override_real_score_gap():
    adv = _advisor()
    strong = ParameterCandidate(value=1.0, score=0.90)   # clearly better, no history
    weak = ParameterCandidate(value=0.4, score=0.50)     # worse, approved 10x
    hist = [{"decision": "0.4"}] * 10
    best = adv._choose_best([strong, weak], hist)
    assert best.value == 1.0          # cap (0.05) cannot bridge a 0.40 gap


def test_historical_bonus_breaks_near_tie():
    adv = _advisor()
    a = ParameterCandidate(value=1.0, score=0.50)        # no history
    b = ParameterCandidate(value=0.4, score=0.49)        # approved, near-tie
    best = adv._choose_best([a, b], [{"decision": "0.4"}, {"decision": "0.4"}])
    assert best.value == 0.4          # bounded bonus legitimately breaks the tie


# ── 3. not-measured fallback (no fabricated comparable metrics) ───────────────

def test_unmeasured_candidates_have_no_fabricated_metrics():
    cands = ParameterAdvisor._unmeasured_leiden_candidates([0.2, 0.5, 0.8])
    assert len(cands) == 3
    for c in cands:
        assert c.metrics.get("measured") is False
        assert "silhouette" not in c.metrics
        assert "modularity" not in c.metrics
        assert any("NOT measured" in f for f in c.flags)
    # neutral mid-range prior ranks the middle value highest
    assert max(cands, key=lambda c: c.score).value == 0.5


def test_score_leiden_uses_prior_when_unmeasured():
    adv = _advisor()
    unmeasured = {"measured": False, "prior_score": 0.42}
    assert adv._score_leiden(unmeasured, {}) == 0.42
    # a measured candidate is still scored on its metrics
    measured = {"silhouette": 0.5, "modularity": None,
                "n_singleton_clusters": 0, "min_cluster_size": 100}
    assert adv._score_leiden(measured, {}) > 0


def test_subprocess_failure_yields_not_measured(monkeypatch):
    adv = _advisor()

    class _Env:
        def run_in_stack(self, **_kw):
            return {"status": "error", "error_type": "FileUnreadable",
                    "details": "boom"}

    adv._env = _Env()
    cands = adv._evaluate_via_subprocess("/x.h5ad", [0.2, 0.5, 0.8], {})
    assert cands and all(c.metrics.get("measured") is False for c in cands)
    assert all("silhouette" not in c.metrics for c in cands)


def test_checkpoint_format_shows_not_measured():
    adv = _advisor()
    adv.memory = _FakeMemory(wings=[], decisions_by_wing={})
    cand = ParameterAdvisor._unmeasured_leiden_candidates([0.5])[0]
    cand.recommended = True
    decision = ParameterDecision(
        decision_id="d1", experiment_id="exp",
        analysis_type="leiden_clustering", parameter_name="resolution",
        candidates=[cand], chosen_value=0.5, chosen_by="advisor",
        biological_context={}, justification="prior only", warnings=cand.flags,
    )
    out = adv.format_for_checkpoint(decision)
    assert "metrics not measured" in out
    assert "NOT measured" in out          # warning surfaced
