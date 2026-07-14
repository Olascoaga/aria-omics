"""Plumbing guards for the c4_report_e2e_false_narrative preprint-freeze lane.

These validate the lane registration and the E2E false-narrative logic (real
compiler boundary, prose-level leakage, blind faithfulness sheet, human scoring)
WITHOUT the aria-rna-env DESeq2 subprocess — a fixture DE result stands in for
the real analysis. The receipt is gated on the independent human faithfulness
gold (human:b2_annotators), which is never synthesized.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from aria.benchmarks import report_false_narrative as rfn
from aria.benchmarks.preprint_freeze import LANES


REPO_ROOT = Path(__file__).resolve().parents[1]
FROZEN_TEMPLATE = (
    REPO_ROOT / "docs/benchmark_results/preprint_v1/human/"
    "report_faithfulness_TEMPLATE.csv"
)

# A fixture DE result in the aria_pseudobulk_da_from_tsv output shape — stands in
# for the real aria-rna-env DESeq2 run so the guard stays fast and hermetic.
FIXTURE_DE = {
    "status": "success", "n_tested": 60, "n_sig": 12, "n_up": 6, "n_down": 6,
    "sig_peaks": [f"GENE{i:03d}" for i in range(12)],
    "lfc_by_peak": {f"GENE{i:03d}": (2.5 if i % 2 else -2.5) for i in range(12)},
    "padj_by_peak": {f"GENE{i:03d}": 0.01 for i in range(12)},
}


def _lane():
    return next(x for x in LANES if x["lane_id"] == "c4_report_e2e_false_narrative")


def _e2e():
    agent_results = rfn.de_to_agent_results(FIXTURE_DE)
    legit, ledger, exp = rfn.build_legit_blocks_and_ledger(agent_results)
    return rfn.compile_e2e(legit, rfn.false_narrative_blocks(), exp, ledger)


def test_lane_is_registered_and_executable():
    lane = _lane()
    assert lane["claims"] == ["claim_4"]
    assert lane["implementation"] == "available"
    assert lane["evidence_kind"] == "e2e_human_gold"
    assert "run_c4_report_e2e_false_narrative.py score" in lane["command"]
    assert "human/report_faithfulness.csv" in lane["command"]


def test_lane_binds_both_envs_and_the_human_resource():
    lane = _lane()
    assert lane["resources"] == [
        "env:aria-env", "env:aria-rna-env", "human:b2_annotators"
    ]
    assert tuple(lane["expected_artifacts"]) == (
        "claim_4/report_e2e/report_e2e_human_gold.json",
    )


def test_injected_false_blocks_cover_every_mechanism():
    labels = {b.metadata.get("b2_label") for b in rfn.false_narrative_blocks()}
    assert labels == set(rfn.FALSE_LABELS)


def test_real_compiler_lets_no_false_narrative_reach_the_report():
    e2e = _e2e()
    # Governance is safe without collapsing informativeness.
    assert e2e["n_false_leaked"] == 0
    assert e2e["n_legit_emitted_with_injection"] >= 1
    assert e2e["safe"] is True
    # Every injected false narrative is either withheld or emitted only as safe,
    # evidence-derived prose (neutralized); none leaks its adversarial claim.
    assert (e2e["n_false_withheld"] + e2e["n_false_neutralized"]
            == e2e["n_false_injected"])
    outcomes = {o["outcome"] for o in e2e["false_outcomes"]}
    assert outcomes <= {"withheld", "neutralized"}


def test_no_emitted_narrative_contains_an_adversarial_claim():
    e2e = _e2e()
    rendered = rfn._normalize(
        " ".join(n["rendered_claim"] for n in e2e["emitted_narratives"]))
    for block in rfn.false_narrative_blocks():
        assert rfn._normalize(block.claim) not in rendered


def test_faithfulness_sheet_is_blind_and_covers_emitted_narratives():
    e2e = _e2e()
    sheet = rfn.export_faithfulness_sheet(e2e)
    rows = list(csv.DictReader(io.StringIO(sheet)))
    assert [*rows[0].keys()] == list(rfn.FAITHFULNESS_COLUMNS)
    assert len(rows) == len(e2e["emitted_narratives"])
    assert all(r["human_verdict"] == "" for r in rows)


def test_load_faithfulness_gold_rejects_out_of_vocab():
    good = "block_id,human_verdict\nb1,faithful\n"
    assert rfn.load_faithfulness_gold(good) == {"b1": "faithful"}
    with pytest.raises(ValueError, match="human_verdict"):
        rfn.load_faithfulness_gold("block_id,human_verdict\nb1,maybe\n")


def test_score_needs_complete_human_gold_and_no_flagged_false():
    e2e = _e2e()
    ids = [n["block_id"] for n in e2e["emitted_narratives"]]
    # A complete, all-faithful human gold plus safe automation -> pass.
    full = {bid: "faithful" for bid in ids}
    m = rfn.score_against_human_gold(full, e2e)
    assert m["status"] == "pass"
    assert m["automated"]["safe"] is True
    # A single human "false" verdict fails the joint claim.
    flagged = dict(full); flagged[ids[0]] = "false"
    assert rfn.score_against_human_gold(flagged, e2e)["status"] == "fail"
    # A partial gold is incomplete, never silently passed.
    assert rfn.score_against_human_gold(
        {ids[0]: "faithful"}, e2e)["status"] == "incomplete"


def test_frozen_faithfulness_template_is_committed_and_unlabeled():
    assert FROZEN_TEMPLATE.is_file()
    rows = list(csv.DictReader(io.StringIO(
        FROZEN_TEMPLATE.read_text(encoding="utf-8"))))
    assert rows and all(r["human_verdict"] == "" for r in rows)
