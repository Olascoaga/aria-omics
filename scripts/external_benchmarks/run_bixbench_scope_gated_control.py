#!/usr/bin/env python3
"""Deterministic Scope-Gated CodeAct Control over frozen v4.9.0 BixBench-60 receipts.

Answers a key methodological ablation question cleanly without modifying the
frozen v4.9.0 external benchmark receipts and with zero additional LLM API calls:

  How much of ARIA's selective-precision advantage over open-ended CodeAct comes
  from ARIA's deterministic pre-execution modality/file contract gate (0 LLM
  calls) versus ARIA's closed-grammar execution & evidence-grounding wall?

For each frozen BixBench-60 backbone receipt (Qwen3.6-35B-A3B local Tantalus and
Claude Sonnet 5.5), this script constructs a strictly deterministic third arm:

  ``codeact_plus_aria_scope_gate``:
    - If ARIA's deterministic pre-execution ``contract_gate["passed"]`` is False,
      the arm deterministically abstains (choice="E", 0 LLM turns, 0 wrong).
    - If ``contract_gate["passed"]`` is True, the arm keeps the exact observed
      ``live_codeact_same_model`` prediction on that task.

Evaluates all three arms across Full-60, In-Contract RNA/DE (n=40), and
Out-of-Scope (n=20) strata, with Wilson 95% CIs and paired exact tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

from scipy.stats import binomtest, fisher_exact

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aria.version import collect_version_metadata  # noqa: E402


FROZEN_RECEIPTS = {
    "qwen36_35b_a3b_tantalus": ROOT / "docs/benchmark_results/external_benchmarks_v49/bixbench_60_tantalus.json",
    "claude_sonnet_55": ROOT / "docs/benchmark_results/external_benchmarks_v49/bixbench_60_claude_sonnet55.json",
}
DEFAULT_OUTPUT = (
    ROOT
    / "docs"
    / "benchmark_results"
    / "external_benchmarks_v49"
    / "bixbench_60_scope_gated_codeact_control.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _wilson_ci95(k: int, n: int) -> list[float]:
    if n <= 0:
        return [0.0, 1.0]
    z = 1.959963984540054
    p = k / n
    denom = 1.0 + (z * z) / n
    center = (p + (z * z) / (2.0 * n)) / denom
    half = (z / denom) * math.sqrt((p * (1.0 - p) + (z * z) / (4.0 * n)) / n)
    return [round(max(0.0, center - half), 4), round(min(1.0, center + half), 4)]


def _score_subset(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    n_correct_supported = sum(
        1 for r in records if r["in_contract"] and r["raw_correct"]
    )
    n_correct_abstention = sum(
        1 for r in records if (not r["in_contract"]) and r["abstained"]
    )
    n_contract_correct_total = n_correct_supported + n_correct_abstention
    n_raw_correct = sum(1 for r in records if r["raw_correct"])
    n_answered = sum(1 for r in records if r["answered"])
    n_raw_wrong = n_answered - n_raw_correct
    n_abstained = sum(1 for r in records if r["abstained"])
    n_no_answer = sum(1 for r in records if r["no_answer"])
    n_silent_errors = sum(1 for r in records if r["silent_error"])
    prec_raw = (n_raw_correct / n_answered) if n_answered else 0.0
    prec_contract = (n_correct_supported / n_answered) if n_answered else 0.0
    return {
        "n_tasks": n,
        "n_answered": n_answered,
        "n_raw_correct": n_raw_correct,
        "n_raw_wrong": n_raw_wrong,
        "n_correct_supported": n_correct_supported,
        "n_correct_abstention": n_correct_abstention,
        "n_contract_correct_total": n_contract_correct_total,
        "n_abstained": n_abstained,
        "n_no_answer": n_no_answer,
        "n_silent_errors": n_silent_errors,
        "coverage": round(n_answered / n, 4) if n else 0.0,
        "unconditional_contract_accuracy": round(n_contract_correct_total / n, 4) if n else 0.0,
        "unconditional_contract_accuracy_ci95": _wilson_ci95(n_contract_correct_total, n),
        "mcq_raw_accuracy": round(n_raw_correct / n, 4) if n else 0.0,
        "mcq_raw_accuracy_ci95": _wilson_ci95(n_raw_correct, n),
        "precision_on_answered_raw": round(prec_raw, 4),
        "precision_on_answered_raw_ci95": _wilson_ci95(n_raw_correct, n_answered),
        "precision_on_answered_contract": round(prec_contract, 4),
        "silent_error_rate": round(n_silent_errors / n, 4) if n else 0.0,
        "silent_error_rate_ci95": _wilson_ci95(n_silent_errors, n),
    }


def _classify_item(item: dict[str, Any], *, force_abstain: bool = False) -> dict[str, Any]:
    in_contract = bool(item.get("supported_by_native_contract"))
    if force_abstain:
        return {
            "in_contract": in_contract,
            "answered": False,
            "abstained": True,
            "no_answer": False,
            "raw_correct": False,
            "silent_error": False,
        }
    grade = item.get("grade") or {}
    abstained = bool(grade.get("abstained") or item.get("abstained") or item.get("choice") == "E")
    no_answer = bool(grade.get("no_answer")) or (not abstained and item.get("choice") not in {"A", "B", "C", "D"})
    answered = (not abstained) and (not no_answer)
    raw_correct = bool(answered and (item.get("choice") == item.get("correct_letter")))
    # Contract-aware silent error matches frozen BixBench-60 grading:
    # answering wrongly in-contract OR answering at all out-of-contract!
    if "silent_error" in grade:
        silent_error = bool(grade["silent_error"])
    else:
        silent_error = bool(answered and ((not in_contract) or (not raw_correct)))
    return {
        "in_contract": in_contract,
        "answered": answered,
        "abstained": abstained,
        "no_answer": no_answer,
        "raw_correct": raw_correct,
        "silent_error": silent_error,
    }


def _paired_wrong_mcnemar(rec_a: list[dict], rec_b: list[dict]) -> dict[str, Any]:
    """Exact paired McNemar test on binary contract silent-error indicator."""
    a_err_b_safe = 0
    a_safe_b_err = 0
    both_err = 0
    both_safe = 0
    for ra, rb in zip(rec_a, rec_b):
        ea = bool(ra["silent_error"])
        eb = bool(rb["silent_error"])
        if ea and eb:
            both_err += 1
        elif ea and not eb:
            a_err_b_safe += 1
        elif not ea and eb:
            a_safe_b_err += 1
        else:
            both_safe += 1
    discordant = a_err_b_safe + a_safe_b_err
    pval = (
        float(binomtest(min(a_err_b_safe, a_safe_b_err), discordant, 0.5).pvalue)
        if discordant > 0
        else 1.0
    )
    return {
        "a_error_b_safe": a_err_b_safe,
        "a_safe_b_error": a_safe_b_err,
        "both_error": both_err,
        "both_safe": both_safe,
        "mcnemar_exact_pvalue": round(pval, 6),
    }


def evaluate_receipt(path: Path) -> dict[str, Any]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    aria_items = {it["task_id"]: it for it in raw["per_task"]["aria_v49_governed"]}
    codeact_items = {it["task_id"]: it for it in raw["per_task"]["live_codeact_same_model"]}
    task_ids = sorted(aria_items.keys())
    assert set(task_ids) == set(codeact_items.keys()) and len(task_ids) == 60

    arms_records: dict[str, list[dict[str, Any]]] = {
        "aria_v49_governed": [],
        "codeact_plus_aria_scope_gate": [],
        "live_codeact_same_model": [],
    }
    gate_stats = {
        "in_contract_passed_gate": 0,
        "in_contract_blocked_by_gate": 0,
        "out_of_scope_passed_gate": 0,
        "out_of_scope_blocked_by_gate": 0,
    }

    for tid in task_ids:
        a = aria_items[tid]
        c = codeact_items[tid]
        in_contract = bool(a.get("supported_by_native_contract"))
        gate_passed = bool((a.get("contract_gate") or {}).get("passed"))

        if in_contract and gate_passed:
            gate_stats["in_contract_passed_gate"] += 1
        elif in_contract and not gate_passed:
            gate_stats["in_contract_blocked_by_gate"] += 1
        elif not in_contract and gate_passed:
            gate_stats["out_of_scope_passed_gate"] += 1
        else:
            gate_stats["out_of_scope_blocked_by_gate"] += 1

        base = {
            "task_id": tid,
            "supported_by_native_contract": in_contract,
            "contract_gate_passed": gate_passed,
        }
        arms_records["aria_v49_governed"].append({**base, **_classify_item(a)})
        arms_records["live_codeact_same_model"].append({**base, **_classify_item(c)})
        arms_records["codeact_plus_aria_scope_gate"].append(
            {**base, **_classify_item(c, force_abstain=(not gate_passed))}
        )

    def _stratified(records: list[dict[str, Any]]) -> dict[str, Any]:
        return {
            "full_60": _score_subset(records),
            "rna_de_40": _score_subset([r for r in records if r["supported_by_native_contract"]]),
            "out_of_scope_20": _score_subset([r for r in records if not r["supported_by_native_contract"]]),
            "gate_passed_subset": _score_subset([r for r in records if r["contract_gate_passed"]]),
        }

    per_arm = {arm: _stratified(recs) for arm, recs in arms_records.items()}

    # Sanity verification against frozen per_arm silent_errors & correct counts
    for frozen_arm in ("aria_v49_governed", "live_codeact_same_model"):
        assert (
            per_arm[frozen_arm]["full_60"]["n_silent_errors"]
            == raw["per_arm"][frozen_arm]["n_confident_wrong"]
        ), f"Mismatch in silent errors for {frozen_arm}"
        assert (
            per_arm[frozen_arm]["full_60"]["n_contract_correct_total"]
            == raw["per_arm"][frozen_arm]["n_correct_total"]
        ), f"Mismatch in contract correct total for {frozen_arm}"

    # Precision Fisher exact comparisons on gate-passed tasks (ARIA vs Scope-Gated CodeAct)
    a_gp = per_arm["aria_v49_governed"]["gate_passed_subset"]
    sg_gp = per_arm["codeact_plus_aria_scope_gate"]["gate_passed_subset"]
    _, fisher_prec_p = fisher_exact(
        [[a_gp["n_raw_correct"], a_gp["n_raw_wrong"]], [sg_gp["n_raw_correct"], sg_gp["n_raw_wrong"]]]
    )

    return {
        "source_receipt": path.relative_to(ROOT).as_posix(),
        "source_receipt_sha256": _sha256(path),
        "model": raw.get("model"),
        "evaluated_aria_tag": raw.get("evaluated_aria_tag"),
        "contract_gate_confusion": gate_stats,
        "per_arm": per_arm,
        "comparisons": {
            "full_60_silent_error_mcnemar_aria_vs_scope_gated_codeact": _paired_wrong_mcnemar(
                arms_records["aria_v49_governed"],
                arms_records["codeact_plus_aria_scope_gate"],
            ),
            "full_60_silent_error_mcnemar_scope_gated_vs_ungated_codeact": _paired_wrong_mcnemar(
                arms_records["codeact_plus_aria_scope_gate"],
                arms_records["live_codeact_same_model"],
            ),
            "rna_de_40_silent_error_mcnemar_aria_vs_scope_gated_codeact": _paired_wrong_mcnemar(
                [r for r in arms_records["aria_v49_governed"] if r["supported_by_native_contract"]],
                [r for r in arms_records["codeact_plus_aria_scope_gate"] if r["supported_by_native_contract"]],
            ),
            "gate_passed_precision_fisher_aria_vs_scope_gated_codeact_pvalue": round(float(fisher_prec_p), 6),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args(argv)

    prov = collect_version_metadata(ROOT)
    if isinstance(prov.get("environment"), dict):
        prov["environment"]["conda_prefix"] = None

    payload = {
        "schema_version": "aria.external_benchmark.bixbench_scope_gated_codeact_control.v1",
        "description": (
            "Deterministic companion ablation applying ARIA's pre-execution "
            "native contract gate (0 LLM calls) to frozen v4.9.0 live CodeAct "
            "trajectories on BixBench-60 across both Qwen3.6-35B-A3B and "
            "Claude Sonnet 5.5, decoupling out-of-scope refusal from closed-"
            "grammar execution and evidence-grounding abstention."
        ),
        "zero_api_spend": True,
        "frozen_receipts_unmodified": True,
        "backbones": {
            key: evaluate_receipt(path)
            for key, path in FROZEN_RECEIPTS.items()
        },
        "provenance": prov,
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if '"/home/' in text or '"/Users/' in text:
        raise RuntimeError("Output JSON contains a machine-absolute path")
    out.write_text(text, encoding="utf-8")

    # Compact console summary
    summary = {}
    for bk, res in payload["backbones"].items():
        summary[bk] = {
            "gate_confusion": res["contract_gate_confusion"],
            "comparisons": res["comparisons"],
            "arms": {
                arm: {
                    "full_60": {
                        k: res["per_arm"][arm]["full_60"][k]
                        for k in ("n_answered", "n_raw_correct", "n_raw_wrong", "n_silent_errors", "n_abstained", "unconditional_contract_accuracy", "precision_on_answered_raw", "silent_error_rate")
                    },
                    "rna_de_40": {
                        k: res["per_arm"][arm]["rna_de_40"][k]
                        for k in ("n_answered", "n_raw_correct", "n_raw_wrong", "n_silent_errors", "n_abstained", "precision_on_answered_raw", "mcq_raw_accuracy")
                    },
                    "out_of_scope_20": {
                        k: res["per_arm"][arm]["out_of_scope_20"][k]
                        for k in ("n_answered", "n_raw_correct", "n_silent_errors", "n_abstained")
                    },
                }
                for arm in ("aria_v49_governed", "codeact_plus_aria_scope_gate", "live_codeact_same_model")
            },
        }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
