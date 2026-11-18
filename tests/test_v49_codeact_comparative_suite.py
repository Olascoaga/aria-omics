"""v4.9 Track B: 40-capsule governed vs. uncontracted CodeAct comparative benchmark (ADR-083 Tier 1).

Verifies that the suite executes deterministically over 40 capsules (8 x 5
archetypes), that the governed ARIA arms meet the frozen-seed acceptance gate
while the uncontracted baseline exhibits a measured silent-error rate, that the
locked-DecisionManifest replay reproduces every governed decision, and that the
committed JSON receipt / vignette match live execution and stay path-portable
and free of external benchmark identifiers.
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from aria.benchmarks.adapter import scan_codebase_for_benchmark_leakage
from aria.benchmarks.codeact_comparative_suite import (
    ARCHETYPES,
    ARMS,
    CAPSULES_PER_ARCHETYPE,
    ROBUSTNESS_SWEEP_SEEDS,
    SCHEMA_VERSION,
    build_capsules,
    capsule_manifest,
    evaluate_claim_requests_in_process,
    median_of_ratios_size_factors,
    payload_fingerprint,
    run_codeact_comparative_suite,
    run_seed_robustness_sweep,
    GOVERNED_POLICY,
    _a3_claim_request,
    _a3_proposal,
    _a3_stats,
)
from aria.utils.decision_manifest import verify_decision_manifest
from aria.utils.provenance import hash_params

REPO_ROOT = Path(__file__).resolve().parents[1]
RECEIPT = REPO_ROOT / "docs" / "benchmark_results" / "v49_comparative_benchmark" / "comparative_benchmark_v49.json"
VIGNETTE = REPO_ROOT / "docs" / "vignettes" / "13_governed_vs_uncontracted_benchmark.md"

FORBIDDEN_BENCHMARK_TOKENS = (
    "DE01_pseudobulk_de",
    "NRM01_sparse_normalization",
    "dr_05_pca_preprocessing_sentinels",
    "gse297023_q2_senescence_not_mechanism",
    "tapestri_ccus",
    "compbiobench-data-v1",
)
ABSOLUTE_PREFIXES = ("/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/")


@pytest.fixture(scope="module")
def live(tmp_path_factory):
    return run_codeact_comparative_suite(scratch_dir=tmp_path_factory.mktemp("codeact_cmp"))


@pytest.fixture(scope="module")
def receipt():
    assert RECEIPT.is_file(), "run scripts/run_v49_codeact_comparative_benchmark.py to emit the receipt"
    return json.loads(RECEIPT.read_text(encoding="utf-8"))


def _iter_strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield str(k)
            yield from _iter_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _iter_strings(v)


# ── capsule construction ─────────────────────────────────────────────────────


def test_forty_deterministic_capsules_across_five_archetypes():
    caps = build_capsules()
    assert len(caps) == 40 == len(ARCHETYPES) * CAPSULES_PER_ARCHETYPE
    assert len({c.capsule_id for c in caps}) == 40
    for arch in ARCHETYPES:
        assert sum(c.archetype == arch for c in caps) == 8

    # Fixed composition of every archetype.
    variants = [c.variant for c in caps]
    assert variants.count("perfect_batch_confounding") == 2
    assert variants.count("technical_replicate_inflation") == 2
    assert sum(v.endswith("batch_balanced_4v4") for v in variants) == 4
    assert sum(v.startswith("null_donor_variance_") for v in variants) == 4
    assert variants.count("planted_donor_consistent_4v4") == 4
    assert variants.count("observational_causal_overreach") == 4
    assert variants.count("speculative_hypothesis_promotion") == 2
    assert variants.count("randomized_interventional_licensed") == 2
    assert sum(v.startswith("raw_sequencer_") for v in variants) == 2
    assert variants.count("unpaired_split_diagonal_imputation") == 2
    assert variants.count("n1_vs_n1_unreplicable_contrast") == 2
    assert variants.count("paired_same_cell_barcode_bijection") == 2
    assert sum(v.startswith("invalid_") for v in variants) == 2
    assert sum(c.should_abstain for c in caps) == 18

    # Null single-cell capsules span 480-640 cells; scrubbed IDs only.
    for c in caps:
        if c.variant.startswith("null_donor_variance_"):
            assert c.payload["counts"].shape[0] in (480, 640)
            assert c.expected["up"] == [] and c.expected["down"] == []
        if c.archetype == ARCHETYPES[4]:
            assert all(f.startswith("FEAT_") for f in c.payload["tables"]["de"]["feature_id"])

    # Byte-identical regeneration from seeds.
    again = build_capsules()
    assert [payload_fingerprint(a.payload) for a in caps] == [payload_fingerprint(b.payload) for b in again]


# ── arm metrics ──────────────────────────────────────────────────────────────


def test_governed_arm_meets_gate_and_baseline_has_measured_silent_errors(live):
    assert live["schema_version"] == SCHEMA_VERSION
    assert live["n_capsules"] == 40
    assert set(live["metrics"]) == set(ARMS)
    for arm in ("aria_governed", "aria_locked_manifest", "aria_zero_llm_deterministic"):
        m = live["metrics"][arm]
        assert m["total_suite_size"] == 40
        assert m["silent_error_rate"] == 0.0
        assert m["conditional_accuracy_on_supported"] == 1.0
        assert m["abstention_precision"] == 1.0
        assert m["abstention_recall"] == 1.0
        assert m["unconditional_end_to_end_accuracy"] == 1.0
    base = live["metrics"]["uncontracted_codeact_baseline"]
    assert base["silent_error_rate"] > 0.35
    assert base["n_confident_wrong"] >= 15
    assert live["silent_error_gap"]["absolute_gap"] > 0.35


def test_trap_outcomes_per_capsule(live):
    by_id = {c["capsule_id"]: c for c in live["capsules"]}
    for rec in live["capsules"]:
        gov = rec["arms"]["aria_governed"]
        base = rec["arms"]["uncontracted_codeact_baseline"]
        if rec["should_abstain"]:
            assert gov["abstained"] and gov["outcome"] == "correct_abstention", rec["capsule_id"]
            assert not base["abstained"] and base["outcome"] == "confident_wrong", rec["capsule_id"]
        else:
            assert gov["outcome"] == "correct_answer", rec["capsule_id"]

    def reason(cid):
        return by_id[cid]["arms"]["aria_governed"]["abstain_reason"]

    assert reason("A1-01") == "design_blocking:condition_covariate_confounding"
    assert reason("A1-03") == "design_blocking:insufficient_condition_replicates"
    assert by_id["A1-03"]["arms"]["aria_governed"]["evidence"]["raw_unit_check"].startswith("rejected:")
    assert reason("A3-01").startswith("causal_claim_not_licensed:")
    assert reason("A3-05") == "speculative_promotion_blocked"
    assert reason("A4-01") == "unsupported_sequencer_native_input:BCL"
    assert reason("A4-02") == "unsupported_sequencer_native_input:POD5"
    assert reason("A4-03").startswith("same_cell_gate_blocked:")
    assert reason("A5-07") == reason("A5-08") == "derived_query_violation"

    # Null donor capsules: zero governed calls, massive cell-level false positives.
    for cid in ("A2-01", "A2-02", "A2-03", "A2-04"):
        gov_ans = by_id[cid]["arms"]["aria_governed"]["answer"]
        assert gov_ans == {"up": [], "down": []}
        assert by_id[cid]["arms"]["uncontracted_codeact_baseline"]["evidence"]["n_significant"] >= 10

    # Randomized capsules are causally licensed and published.
    for cid in ("A3-07", "A3-08"):
        ev = by_id[cid]["arms"]["aria_governed"]["evidence"]
        assert ev["causal_license"]["licensed"] is True
        assert "causal_experimental" in ev["published_tiers"]


def test_locked_manifest_and_zero_llm_reproduce_every_governed_decision(live):
    lc = live["concordance"]["aria_locked_manifest_vs_aria_governed"]
    assert lc["n_tasks"] == 40 and lc["n_identical_decisions"] == 40
    assert lc["decision_concordance"] == 1.0 and lc["prediction_concordance"] == 1.0
    zc = live["concordance"]["aria_zero_llm_deterministic_vs_aria_governed"]
    assert zc["decision_concordance"] == 1.0 and zc["mismatched_capsules"] == []
    dmv = live["decision_manifest_verification"]
    assert dmv["n_manifests_verified"] == 40
    assert dmv["n_tampered_manifests_detected"] == 40
    assert dmv["n_input_fingerprints_matched"] == 40
    iso = live["zero_llm_isolation"]
    assert iso["llm_handle"] is None and iso["n_capsules_regenerated_identically"] == 40
    assert iso["isolated_env_flags"]["ARIA_LLM_CACHE"] == "0"


def test_capsule_manifest_detects_policy_tampering():
    cap = build_capsules()[0]
    manifest = capsule_manifest(cap, GOVERNED_POLICY)
    assert verify_decision_manifest(manifest) == (True, "verified")
    tampered = json.loads(json.dumps(manifest))
    tampered["locked_parameter_decisions"]["policy"]["alpha"] = 0.2
    ok, reason = verify_decision_manifest(tampered)
    assert not ok and reason == "manifest_sha256_mismatch"


def test_stack_isolated_claim_gate_ipc_matches_in_process():
    caps = [c for c in build_capsules() if c.archetype == ARCHETYPES[2]]
    reqs = []
    for c in caps:
        s = _a3_stats(c, GOVERNED_POLICY)
        reqs.append(_a3_claim_request(c, s, _a3_proposal(c, s, None)))
    in_proc = evaluate_claim_requests_in_process(reqs)
    proc = subprocess.run(
        [sys.executable, "-m", "aria.benchmarks.codeact_comparative_suite", "--claim-gate-ipc"],
        input=json.dumps(reqs), capture_output=True, text=True, cwd=str(REPO_ROOT), check=True,
    )
    ipc = json.loads(proc.stdout.split("\n@@ARIA_CLAIM_GATE_IPC@@\n", 1)[1])
    assert ipc == json.loads(json.dumps(in_proc))


def test_median_of_ratios_resists_compositional_shift():
    rng = np.random.default_rng(7)
    base = rng.integers(50, 500, size=(100, 1)).astype(float)
    mat = np.repeat(base, 6, axis=1)
    mat[:5, :3] *= 20.0  # a few features massively induced in one arm
    sf = median_of_ratios_size_factors(mat)
    assert np.allclose(sf[:3], sf[3:].mean(), rtol=0.05)
    lib = mat.sum(axis=0)
    assert lib[:3].mean() / lib[3:].mean() > 1.5  # CPM scaling would be biased


# ── committed receipt & vignette ─────────────────────────────────────────────


def test_receipt_matches_live_execution(live, receipt):
    assert receipt["schema_version"] == SCHEMA_VERSION
    assert receipt["receipt_sha256"] == live["receipt_sha256"]
    assert receipt["decisions_sha256"] == live["decisions_sha256"]
    for arm in ARMS:
        assert receipt["metrics"][arm] == live["metrics"][arm]
    live_digests = {
        c["capsule_id"]: {a: c["arms"][a]["decision_digest"] for a in ARMS} for c in live["capsules"]
    }
    rec_digests = {
        c["capsule_id"]: {a: c["arms"][a]["decision_digest"] for a in ARMS} for c in receipt["capsules"]
    }
    assert rec_digests == live_digests
    body = {k: v for k, v in receipt.items() if k not in {"runtime", "seed_robustness_sweep", "receipt_sha256"}}
    assert hash_params(body) == receipt["receipt_sha256"]


def test_receipt_seed_sweep_is_reproducible(receipt):
    sweep = receipt["seed_robustness_sweep"]
    assert sweep["n_seeds"] == len(ROBUSTNESS_SWEEP_SEEDS) == 20
    body = {k: v for k, v in sweep.items() if k != "sweep_sha256"}
    assert hash_params(body) == sweep["sweep_sha256"]
    # Re-run the first two seeds live and compare per-seed rows exactly.
    partial = run_seed_robustness_sweep(ROBUSTNESS_SWEEP_SEEDS[:2])
    assert partial["per_seed"] == sweep["per_seed"][:2]
    # Governance gates never fail across re-draws: no trap answered, no sign error.
    for d in sweep["aria_governed_failure_details"]:
        assert d["should_abstain"] is False
        assert d.get("n_sign_errors", 0) == 0
    assert sweep["summary"]["n_seeds_governed_abstention_precision_and_recall_1"] == 20
    assert sweep["summary"]["baseline_silent_error_rate"]["min"] > 0.35


def test_vignette_documents_all_metrics_and_limitations(receipt):
    text = VIGNETTE.read_text(encoding="utf-8")
    assert receipt["receipt_sha256"] in text
    for key in (
        "native_contract_coverage",
        "conditional_accuracy_on_supported",
        "unconditional_end_to_end_accuracy",
        "abstention_precision",
        "abstention_recall",
        "silent_error_rate",
    ):
        assert key in text
    for arch in ARCHETYPES:
        assert arch in text
    for i in range(1, 9):
        for a in range(1, 6):
            assert f"`A{a}-{i:02d}`" in text
    assert "Honest scope & limitations" in text
    assert "Seed-robustness sweep" in text


def test_public_artifacts_are_path_portable(receipt):
    for s in _iter_strings(receipt):
        assert not s.startswith(ABSOLUTE_PREFIXES), s
    for path in (RECEIPT, VIGNETTE):
        text = path.read_text(encoding="utf-8")
        for marker in ("/home/", "/tmp/", "/Users/"):
            assert marker not in text, f"{path.name} leaks {marker}"


def test_zero_benchmark_leakage_in_package_and_artifacts():
    assert scan_codebase_for_benchmark_leakage(REPO_ROOT / "aria", FORBIDDEN_BENCHMARK_TOKENS) == []
    for path in (RECEIPT, VIGNETTE, REPO_ROOT / "scripts" / "run_v49_codeact_comparative_benchmark.py"):
        text = path.read_text(encoding="utf-8")
        for tok in FORBIDDEN_BENCHMARK_TOKENS:
            assert tok not in text, f"{tok} leaked into {path.name}"
