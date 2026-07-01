"""H20 guards: the blind-evaluation harness + promotion go/no-go gate.

Exercises the deterministic mechanics with FAKE proposers (no LLM, no network):
factuality isolates invention, the grading sheet is arm-blind but de-anonymisable,
stability is run-to-run overlap, and the promotion gate stays non-promotable until
mechanical AND human criteria both pass. The real-LLM run is a separate process
step that fills the human panel scores.
"""

from __future__ import annotations

import json

from aria.agents.narrative.hypothesis import parse_hypotheses
from aria.benchmarks.hypothesis_blind_eval import (
    build_grading_sheet,
    evaluate_promotion_gate,
    factuality,
    run_blind_eval,
    scenarios,
    stability,
)


def _governed_proposer(signals, exp_ctx):
    """A faithful fake: one grounded, gate-passing hypothesis per scenario."""
    by = {s.entity.lower(): s for s in signals}
    ents = [s.entity for s in signals[:2]]
    oc = [
        {"signal_id": by[e.lower()].signal_id,
         "stated_direction": by[e.lower()].direction}
        for e in ents
    ]
    caveats = sorted({c for s in signals for c in s.caveats_inherited})
    item = {
        "id": "g1",
        "mechanism": f"{ents[0]} and {ents[1]} may act in a shared program",
        "entities": ents,
        "observation_refs": [signals[0].audited_node_ref],
        "observed_claims": oc,
        "experiment": {
            "perturbation": f"{ents[0]} knockdown",
            "readout": f"{ents[1]} by qPCR",
            "predicted_direction": "decrease",
            "refuting_outcome": "no change",
        },
        "devils_advocate": {
            "simpler_explanation": "a shared upstream driver",
            "confounds": caveats,
        },
    }
    return parse_hypotheses(json.dumps([item]))


def _inventing_proposer(signals, exp_ctx):
    """An ungoverned fake that fabricates an entity absent from the evidence."""
    item = {
        "id": "b1",
        "mechanism": "FOXP3 secretly drives the whole program",
        "entities": ["FOXP3"],
        "observation_refs": [],
        "observed_claims": [],
        "experiment": {"perturbation": "x", "readout": "y",
                       "predicted_direction": "up", "refuting_outcome": "z"},
        "devils_advocate": {"simpler_explanation": "", "confounds": []},
    }
    return parse_hypotheses(json.dumps([item]))


# ── factuality ───────────────────────────────────────────────────────────────

def test_factuality_scores_governed_and_baseline_apart():
    scen = scenarios()["senescence_bulk_rna"]
    gov = _governed_proposer(scen["signals"], scen["exp_ctx"])
    base = _inventing_proposer(scen["signals"], scen["exp_ctx"])
    fg = factuality(gov, scen["signals"], scen["run_ledger"])
    fb = factuality(base, scen["signals"], scen["run_ledger"])
    assert fg["rate"] == 1.0 and fg["invented"] == []
    assert fb["rate"] == 0.0
    assert "FOXP3" in fb["invented"][0]["fabricated_entities"]


def test_stability_is_one_for_a_deterministic_arm():
    scen = scenarios()["senescence_bulk_rna"]
    stab = stability(lambda: _governed_proposer(scen["signals"], scen["exp_ctx"]))
    assert stab["mean_pairwise_jaccard"] == 1.0


# ── blind grading sheet ──────────────────────────────────────────────────────

def test_grading_sheet_is_arm_blind_and_deanonymisable():
    scen = scenarios()["senescence_bulk_rna"]
    outputs = {
        "senescence_bulk_rna": {
            "governed": _governed_proposer(scen["signals"], scen["exp_ctx"]),
            "ungoverned": _inventing_proposer(scen["signals"], scen["exp_ctx"]),
        }
    }
    sheet = build_grading_sheet(outputs, seed=3)
    # rows never leak the arm or scenario.
    for row in sheet["rows"]:
        assert "arm" not in row and "scenario" not in row
        assert row["blind_id"].startswith("hyp_")
    # the SEPARATE key de-anonymises every row.
    for row in sheet["rows"]:
        assert row["blind_id"] in sheet["key"]
    arms = {v["arm"] for v in sheet["key"].values()}
    assert arms == {"governed", "ungoverned"}
    # deterministic given the seed.
    again = build_grading_sheet(outputs, seed=3)
    assert [r["blind_id"] for r in sheet["rows"]] == [r["blind_id"] for r in again["rows"]]


# ── promotion gate ───────────────────────────────────────────────────────────

def _report(gov=1.0, base=0.0, evasions=0, stab=1.0, key=None):
    return {
        "governed": {"factuality_overall": {"rate": gov},
                     "stability_overall": {"mean_pairwise_jaccard": stab}},
        "ungoverned": {"factuality_overall": {"rate": base}},
        "redteam": {"evasions": evasions},
        "grading_key": key or {},
    }


def test_gate_awaits_human_review_when_mechanical_passes():
    gate = evaluate_promotion_gate(_report())
    assert gate["mechanical"]["passed"] is True
    assert gate["promotable"] is False
    assert gate["reason"] == "awaiting_human_review"


def test_gate_fails_on_redteam_evasion():
    gate = evaluate_promotion_gate(_report(evasions=1))
    assert gate["mechanical"]["criteria"]["redteam_zero_evasions"] is False
    assert gate["mechanical"]["passed"] is False
    assert gate["promotable"] is False


def test_gate_fails_when_governance_does_not_beat_baseline():
    gate = evaluate_promotion_gate(_report(gov=1.0, base=1.0))
    assert gate["mechanical"]["criteria"]["governance_beats_baseline"] is False
    assert gate["promotable"] is False


def test_gate_fails_on_governed_invention():
    gate = evaluate_promotion_gate(_report(gov=0.9))
    assert gate["mechanical"]["criteria"]["governed_zero_invention"] is False


def test_gate_promotes_only_when_mechanical_and_human_pass():
    key = {"hyp_a": {"arm": "governed"}, "hyp_b": {"arm": "governed"}}
    human = {
        "hyp_a": {"plausibility": 4, "novelty": 4, "experimental_utility": 4},
        "hyp_b": {"plausibility": 5, "novelty": 3, "experimental_utility": 4},
    }
    gate = evaluate_promotion_gate(_report(key=key), human_scores=human)
    assert gate["human"]["passed"] is True
    assert gate["promotable"] is True
    assert gate["reason"] == "all_criteria_passed"


def test_gate_blocks_promotion_on_weak_human_panel():
    key = {"hyp_a": {"arm": "governed"}}
    human = {"hyp_a": {"plausibility": 2, "novelty": 2, "experimental_utility": 2}}
    gate = evaluate_promotion_gate(_report(key=key), human_scores=human)
    assert gate["human"]["passed"] is False
    assert gate["promotable"] is False
    assert gate["reason"] == "human_panel_below_threshold"


# ── full run ─────────────────────────────────────────────────────────────────

def test_run_blind_eval_end_to_end_with_fakes():
    report = run_blind_eval(_governed_proposer, _inventing_proposer, seed=7)
    assert report["governed"]["factuality_overall"]["rate"] == 1.0
    assert report["ungoverned"]["factuality_overall"]["rate"] == 0.0
    assert report["redteam"]["passed"] is True
    assert report["promotion_gate"]["mechanical"]["passed"] is True
    assert report["promotion_gate"]["promotable"] is False  # awaits humans
    # every scenario produced both arms; the sheet has no key embedded.
    assert set(report["scenarios"]) == set(scenarios())
    assert "key" not in report["grading_sheet"]
