"""FASE 8: fail-closed inventory and clean-publication guards."""

from __future__ import annotations

import json
from pathlib import Path

from aria.benchmarks.preprint_freeze import (
    CLAIM_IDS,
    LANES,
    _sanitize_public_json_artifact,
    build_inventory,
    execute_lane,
    write_inventory,
)


def _all_available():
    resources = {r for lane in LANES for r in lane["resources"]}
    return {resource: True for resource in resources}


def test_inventory_covers_every_claim_with_unique_required_lanes(tmp_path):
    payload = build_inventory(Path.cwd(), tmp_path, resource_overrides=_all_available())

    assert tuple(c["claim_id"] for c in payload["claims"]) == CLAIM_IDS
    ids = [lane["lane_id"] for lane in payload["lanes"]]
    assert len(ids) == len(set(ids))
    assert all(lane["required_for_freeze"] is True for lane in payload["lanes"])
    assert all(lane["claims"] for lane in payload["lanes"])


def test_historical_artifacts_never_make_the_freeze_ready(tmp_path):
    # A historical-looking file outside the dedicated root has no effect. Even
    # with every resource available, no current-commit receipts means blocked.
    historical = tmp_path.parent / "a1_bulk_de_v4.5.5.json"
    historical.write_text('{"status":"pass"}', encoding="utf-8")
    payload = build_inventory(Path.cwd(), tmp_path, resource_overrides=_all_available())

    assert payload["freeze_gate"]["ready"] is False
    assert payload["freeze_gate"]["n_verified_lanes"] == 0
    assert payload["freeze_gate"]["tag_action_authorized"] is False
    assert all(claim["status"] == "blocked" for claim in payload["claims"])


def test_missing_implementation_and_human_inputs_are_explicit(tmp_path):
    payload = build_inventory(Path.cwd(), tmp_path, resource_overrides=_all_available())
    by_id = {lane["lane_id"]: lane for lane in payload["lanes"]}

    assert by_id["c1_h9_fastq_e2e"]["status"] == "blocked_missing_implementation"
    assert by_id["c3_blind_multifactorial_corpus"]["status"] == "blocked_missing_implementation"
    assert by_id["c4_b2_multi_annotator_gold"]["status"] == "ready_to_run"
    assert by_id["c4_report_e2e_false_narrative"]["status"] == "blocked_missing_implementation"


def test_inventory_is_path_portable_and_atomic(tmp_path):
    payload = build_inventory(Path.cwd(), tmp_path, resource_overrides=_all_available())
    out = write_inventory(payload, tmp_path / "inventory.json")
    text = out.read_text(encoding="utf-8")
    loaded = json.loads(text)

    assert loaded["output_root"] == "docs/benchmark_results/preprint_v1"
    assert '"/home/' not in text
    assert '"/tmp/' not in text
    assert not (tmp_path / "inventory.json.tmp").exists()


def test_receipt_must_match_current_commit_to_verify(tmp_path):
    payload = build_inventory(Path.cwd(), tmp_path, resource_overrides=_all_available())
    write_inventory(payload, tmp_path / "inventory.json")
    receipts = tmp_path / "receipts"
    receipts.mkdir()
    lane_id = "c5_semantic_negation"
    write_inventory({
        "status": "pass",
        "git_commit": "stale-commit",
    }, receipts / f"{lane_id}.json")

    refreshed = build_inventory(
        Path.cwd(), tmp_path, resource_overrides=_all_available()
    )
    by_id = {lane["lane_id"]: lane for lane in refreshed["lanes"]}
    assert by_id[lane_id]["status"] == "stale_receipt"
    assert refreshed["freeze_gate"]["ready"] is False


def test_executor_refuses_unimplemented_lane(tmp_path):
    try:
        execute_lane(Path.cwd(), tmp_path, "c1_h9_fastq_e2e")
    except RuntimeError as exc:
        assert "no executable implementation" in str(exc)
    else:
        raise AssertionError("unimplemented lane was executed")


def test_public_json_sanitizer_removes_only_conda_prefix(tmp_path):
    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps({
        "provenance": {
            "environment": {
                "conda_prefix": "/machine/local/envs/aria-env",
                "env_name": "aria-env",
            },
            "unexpected_path": "/machine/local/input.tsv",
        },
    }), encoding="utf-8")

    _sanitize_public_json_artifact(artifact)

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    assert payload["provenance"]["environment"]["conda_prefix"] is None
    assert payload["provenance"]["environment"]["env_name"] == "aria-env"
    assert payload["provenance"]["unexpected_path"] == "/machine/local/input.tsv"
