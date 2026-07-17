"""FASE 8: fail-closed inventory and clean-publication guards."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

from aria.benchmarks.preprint_freeze import (
    CLAIM_IDS,
    LANES,
    _receipt_status,
    _sanitize_public_json_artifact,
    _source_snapshot_hash,
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

    assert by_id["c1_h9_fastq_e2e"]["status"] == "ready_to_run"
    assert by_id["c2_scatac_donor_aware"]["status"] == "ready_to_run"
    assert by_id["c3_blind_multifactorial_corpus"]["status"] == "ready_to_run"
    assert by_id["c4_b2_multi_annotator_gold"]["status"] == "ready_to_run"
    assert by_id["c4_report_e2e_false_narrative"]["status"] == "ready_to_run"
    assert by_id["c5_multimodal_null_permutations"]["status"] == "ready_to_run"
    # FASE 8 completion invariant: every required lane now has an executable
    # implementation; nothing remains blocked for want of code.
    assert all(lane["implementation"] == "available" for lane in payload["lanes"])
    assert not any(lane["status"] == "blocked_missing_implementation"
                   for lane in payload["lanes"])


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


def test_executor_rejects_unknown_lane(tmp_path):
    # With every registered lane now implemented, the durable executor guard is
    # the unknown-lane rejection.
    import pytest
    with pytest.raises(KeyError):
        execute_lane(Path.cwd(), tmp_path, "not_a_registered_lane")


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


def test_receipt_revalidates_artifact_hash(tmp_path):
    lane = next(item for item in LANES if item["lane_id"] == "c1_a1_synthetic_bulk_de")
    artifact_paths = lane["expected_artifacts"]
    artifacts = []
    for relative in artifact_paths:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative, encoding="utf-8")
        artifacts.append({
            "path": relative,
            "size_bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    provenance = {
        "git_commit": "commit",
        "git_tree_sha": "tree",
        "source_snapshot_hash": "source",
        "workflow_hash": "workflow",
        "aria_version": "4.7.0",
    }
    receipt = {
        "schema_version": "aria.preprint_freeze.receipt.v1",
        "lane_id": lane["lane_id"],
        "claims": lane["claims"],
        "status": "pass",
        **provenance,
        "environment": lane["environment"],
        "command": lane["command"],
        "returncode": 0,
        "artifacts": artifacts,
    }
    receipts = tmp_path / "receipts"
    receipts.mkdir()
    write_inventory(receipt, receipts / f"{lane['lane_id']}.json")

    assert _receipt_status(lane, tmp_path, provenance) == "verified"
    (tmp_path / artifact_paths[0]).write_text("tampered", encoding="utf-8")
    assert _receipt_status(lane, tmp_path, provenance) == "artifact_mismatch"


def test_external_comparator_receipt_covers_supporting_files():
    lane = next(item for item in LANES if item["lane_id"] == "c1_a1_external_comparators")
    expected = set(lane["expected_artifacts"])

    assert lane["command"].startswith("PYTHONPATH=. ")
    assert "claim_1/a1_external/a1_external_comparators.json" in expected
    assert "claim_1/a1_external/inputs/a1_counts.tsv" in expected
    assert "claim_1/a1_external/r_outputs/deseq2.tsv" in expected
    assert "claim_1/a1_external/r_outputs/edgeR_QLF.tsv" in expected
    assert "claim_1/a1_external/r_outputs/limma_voom.tsv" in expected


def test_source_snapshot_hash_excludes_only_output_root(tmp_path):
    repo = tmp_path / "repo"
    output = repo / "docs/benchmark_results/preprint_v1"
    output.mkdir(parents=True)
    source = repo / "aria/source.py"
    source.parent.mkdir()
    source.write_text("VALUE = 1\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    baseline = _source_snapshot_hash(repo, output)

    artifact = output / "inventory.json"
    artifact.write_text("{}\n", encoding="utf-8")
    subprocess.run(["git", "add", str(artifact)], cwd=repo, check=True)
    assert _source_snapshot_hash(repo, output) == baseline

    source.write_text("VALUE = 2\n", encoding="utf-8")
    subprocess.run(["git", "add", str(source)], cwd=repo, check=True)
    assert _source_snapshot_hash(repo, output) != baseline
