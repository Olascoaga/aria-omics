"""Plumbing guards for the c3_blind_multifactorial_corpus preprint-freeze lane.

These validate the lane registration, the held-out multifactorial corpus, the
blindness invariant and the blind scoring boundary. The lane is human-gold: its
receipt cannot be produced until an INDEPENDENT human authors design_gold.csv
(the protocol forbids synthesizing it), so the lane stays blocked on that
resource. The synthetic gold here is a test fixture that exercises the scorer,
never freeze evidence.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from aria.benchmarks.governance_b1 import build_corpus
from aria.benchmarks.design_blind_gold import (
    DECISIONS, SHEET_COLUMNS, build_multifactorial_corpus, export_design_sheet,
    load_design_gold, score_blind_design_gold,
)
from aria.benchmarks.preprint_freeze import LANES


REPO_ROOT = Path(__file__).resolve().parents[1]
FROZEN_TEMPLATE = (
    REPO_ROOT / "docs/benchmark_results/preprint_v1/human/design_gold_TEMPLATE.csv"
)


def _lane():
    return next(x for x in LANES if x["lane_id"] == "c3_blind_multifactorial_corpus")


def _synthetic_gold():
    """A fixture human gold spanning every scenario (NOT freeze evidence)."""
    rule = {
        "factorial_confounded": "block",
        "factorial_continuous": "block",
        "factorial_under_replicated": "block",
        "factorial_incomplete": "escalate",
        "factorial_low_power": "escalate",
        "factorial_degenerate": "escalate",
        "factorial_complete": "infer",
        "nested_design": "infer",
    }
    return {c.case_id: rule[c.category] for c in build_multifactorial_corpus()}


def test_lane_is_registered_and_executable():
    lane = _lane()
    assert lane["claims"] == ["claim_3"]
    assert lane["implementation"] == "available"
    assert lane["evidence_kind"] == "human_gold"
    assert "run_c3_design_blind_gold.py score" in lane["command"]
    assert "human/design_gold.csv" in lane["command"]


def test_lane_declares_manifest_and_human_resource():
    lane = _lane()
    assert tuple(lane["expected_artifacts"]) == (
        "claim_3/blind_design_gold.json",
    )
    assert "human:design_gold" in lane["resources"]


def test_corpus_is_multifactorial_and_held_out_from_b1():
    mf = build_multifactorial_corpus()
    assert len(mf) >= 8
    b1_ids = {c.case_id for c in build_corpus()}
    mf_ids = {c.case_id for c in mf}
    # Held-out: ARIA's primitives were not tuned on these ids.
    assert not (mf_ids & b1_ids)
    # Every scenario carries a real second factor / covariate structure.
    assert all(c.covariates for c in mf)


def test_corpus_ships_no_gold_labels_blindness_invariant():
    # The correct decision is authored by a human, never in code.
    assert all(c.gold == "" for c in build_multifactorial_corpus())


def test_export_sheet_is_blind_and_complete():
    sheet = export_design_sheet()
    rows = list(csv.DictReader(io.StringIO(sheet)))
    assert [*rows[0].keys()] == list(SHEET_COLUMNS)
    assert len(rows) == len(build_multifactorial_corpus())
    # No leaked labels: every gold_decision cell is empty.
    assert all(r["gold_decision"] == "" for r in rows)


def test_load_design_gold_rejects_out_of_vocab():
    good = "case_id,gold_decision\nmf2x2_balanced_n3,infer\n"
    assert load_design_gold(good) == {"mf2x2_balanced_n3": "infer"}
    with pytest.raises(ValueError, match="gold_decision"):
        load_design_gold("case_id,gold_decision\nx,maybe\n")


def test_score_blind_gold_runs_aria_and_reports_agreement():
    manifest = score_blind_design_gold(_synthetic_gold())
    assert manifest["status"] == "pass"
    assert manifest["corpus_disjoint_from_b1"] is True
    # Real 3x3 governance confusion matrix over the decision vocabulary.
    conf = manifest["confusion_matrix"]
    assert set(conf) == set(DECISIONS)
    assert manifest["summary"]["n_scored"] == len(build_multifactorial_corpus())
    assert "kappa" in manifest["cohen_kappa"]
    # ARIA's decision is computed for every scored case (blind).
    assert all(pc["aria_decision"] in DECISIONS for pc in manifest["per_case"])


def test_partial_gold_is_incomplete_not_fabricated():
    partial = {"mf2x2_balanced_n3": "infer"}
    manifest = score_blind_design_gold(partial)
    assert manifest["status"] == "incomplete"
    assert len(manifest["unscored_scenarios"]) > 0


def test_frozen_template_is_committed_and_unlabeled():
    assert FROZEN_TEMPLATE.is_file()
    rows = list(csv.DictReader(io.StringIO(
        FROZEN_TEMPLATE.read_text(encoding="utf-8"))))
    assert rows and all(r["gold_decision"] == "" for r in rows)
