"""H20 guards: the red-team battery must report ZERO evasions.

The versioned, deterministic proof that every wall built across rounds 1-3
(H1-H19) holds under adversarial input. Zero evasions is a hard promotion-gate
criterion — if any case here starts evading, the SPECULATIVE tier cannot be
promoted.
"""

from __future__ import annotations

from aria.benchmarks.hypothesis_redteam import (
    build_cases,
    evaluate_redteam,
    run_case,
    run_redteam,
)

# Which agent gate should catch each "gate"-kind attack.
_EXPECTED_GATE = {
    "invented_entity_structured": "grounding",
    "invented_entity_in_mechanism": "grounding",
    "invented_entity_in_readout": "grounding",
    "directional_contradiction": "grounding",
    "unknown_signal_id": "grounding",
    "no_observed_claims": "grounding",
    "dropped_confound": "devils_advocate",
    "non_falsifiable_experiment": "falsifiability",
    "assertive_causal_language": "language",
}


def test_redteam_reports_zero_evasions():
    summary = evaluate_redteam()
    assert summary["passed"] is True
    assert summary["evasions"] == 0
    assert summary["evaded_cases"] == []
    assert summary["n_cases"] >= 14


def test_every_case_holds():
    for result in run_redteam():
        assert result["evaded"] is False, f"{result['name']} evaded: {result['detail']}"


def test_gate_attacks_are_caught_by_the_expected_wall():
    by_name = {r["name"]: r for r in run_redteam()}
    for name, gate in _EXPECTED_GATE.items():
        caught = by_name[name]["detail"].get("caught_by", [])
        assert gate in caught, f"{name} not caught by {gate} (got {caught})"


def test_battery_covers_all_walls():
    walls = {c.wall for c in build_cases()}
    for expected in (
        "grounding", "grounding_h15", "governance", "quarantine_h17",
        "verification_h14", "devils_advocate_h16", "falsifiability",
        "language", "ranking_h18",
    ):
        assert expected in walls


def test_verification_absence_is_fail_closed():
    result = next(
        r for r in run_redteam() if r["name"] == "verification_evidence_absent"
    )
    assert result["evaded"] is False
    assert result["detail"]["ran"] is False
    assert result["detail"]["reason"] == "verification_evidence_absent"


def test_nested_quarantine_node_raises():
    case = next(c for c in build_cases() if c.kind == "quarantine")
    result = run_case(case)
    assert result["evaded"] is False
    assert "hypothesis://" in result["detail"]["raised"]
