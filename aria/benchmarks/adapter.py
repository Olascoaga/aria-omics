"""Zero-Leakage External Benchmark Adapter & Full-Denominator Scoring (Section 5.3).

Enforces strict scientific honesty across external agent benchmarks (`GenoTEX`,
`BixBench`, `scBench`, `scBench-Long`, `CompBioBench`):
1. Per-task isolation (`isolated_benchmark_environment`: fresh `ARIA_MEMORY_DB`,
   `ARIA_LLM_CACHE=0`, `ARIA_REPRODUCIBLE=1`, separate scratch workspace).
2. Full-denominator accounting reporting ALL 6 mandatory metrics:
   - `total_suite_size`
   - `native_contract_coverage` (`n_supported / total_suite_size`)
   - `conditional_accuracy_on_supported` (`n_correct_supported / n_supported`)
   - `unconditional_end_to_end_accuracy` (`n_correct_total / total_suite_size`)
   - `abstention_precision` & `abstention_recall`
   - `silent_error_rate` (`n_confident_wrong / total_suite_size` — ARIA's core
     safety metric against open-ended `CodeAct` baselines).
3. Static anti-leakage codebase scanner (`scan_codebase_for_benchmark_leakage`)
   ensuring zero benchmark task IDs, accessions, or answer tokens ever enter `aria/`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

BENCHMARK_ADAPTER_SCHEMA = "aria.benchmark_adapter.v1"


@dataclass(frozen=True)
class BenchmarkTaskSpec:
    benchmark_name: str
    task_id: str
    supported_by_native_contract: bool
    should_abstain: bool = False


@dataclass(frozen=True)
class BenchmarkPrediction:
    task_id: str
    abstained: bool
    is_correct: bool
    abstain_reason: str | None = None
    ledger_verified: bool = True


def isolated_benchmark_environment(
    scratch_dir: str | Path,
    *,
    base_env: dict[str, str] | None = None,
) -> dict[str, str]:
    """Construct a zero-cache, isolated-memory, reproducible env for one benchmark task."""
    root = Path(scratch_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    env = dict(base_env or {})
    env.update({
        "ARIA_MEMORY_DB": str(root / "isolated_memory.db"),
        "ARIA_FORK_LEDGER_PATH": str(root / "isolated_fork_ledger.jsonl"),
        "ARIA_LLM_CACHE": "0",
        "ARIA_REPRODUCIBLE": "1",
        "ARIA_DISABLE_HISTORICAL_BONUS": "1",
        "PYTHONHASHSEED": "0",
    })
    return env


def evaluate_benchmark_predictions(
    suite_name: str,
    tasks: Iterable[BenchmarkTaskSpec],
    predictions: Iterable[BenchmarkPrediction],
) -> dict[str, Any]:
    """Compute full-denominator coverage, conditional/unconditional accuracy, abstention, and silent-error rate."""
    task_list = list(tasks)
    if not task_list:
        raise ValueError("Benchmark task list cannot be empty")
    pred_by_id = {p.task_id: p for p in predictions}

    total = len(task_list)
    n_supported = 0
    n_correct_supported = 0
    n_correct_total = 0
    n_abstained = 0
    n_true_abstain = 0
    n_should_abstain = 0
    n_confident_wrong = 0

    for t in task_list:
        must_abstain = bool(t.should_abstain or not t.supported_by_native_contract)
        if t.supported_by_native_contract and not t.should_abstain:
            n_supported += 1
        if must_abstain:
            n_should_abstain += 1

        pred = pred_by_id.get(t.task_id)
        if pred is None:
            # Missing prediction counts as neither answered nor correctly answered
            continue

        if pred.abstained:
            n_abstained += 1
            if must_abstain:
                n_true_abstain += 1
                if pred.is_correct:
                    n_correct_total += 1
        else:
            if pred.is_correct and pred.ledger_verified:
                n_correct_total += 1
                if t.supported_by_native_contract and not t.should_abstain:
                    n_correct_supported += 1
            else:
                n_confident_wrong += 1

    coverage = n_supported / total
    cond_acc = (n_correct_supported / n_supported) if n_supported > 0 else None
    uncond_acc = n_correct_total / total
    abst_prec = (n_true_abstain / n_abstained) if n_abstained > 0 else None
    abst_rec = (n_true_abstain / n_should_abstain) if n_should_abstain > 0 else None
    silent_err = n_confident_wrong / total

    return {
        "schema_version": BENCHMARK_ADAPTER_SCHEMA,
        "suite_name": suite_name,
        "total_suite_size": total,
        "n_supported_tasks": n_supported,
        "n_should_abstain_tasks": n_should_abstain,
        "n_abstained": n_abstained,
        "n_correct_supported": n_correct_supported,
        "n_correct_total": n_correct_total,
        "n_confident_wrong": n_confident_wrong,
        "native_contract_coverage": round(coverage, 6),
        "conditional_accuracy_on_supported": (
            round(cond_acc, 6) if cond_acc is not None else None
        ),
        "unconditional_end_to_end_accuracy": round(uncond_acc, 6),
        "abstention_precision": (
            round(abst_prec, 6) if abst_prec is not None else None
        ),
        "abstention_recall": (
            round(abst_rec, 6) if abst_rec is not None else None
        ),
        "silent_error_rate": round(silent_err, 6),
    }


def scan_codebase_for_benchmark_leakage(
    aria_package_dir: str | Path,
    forbidden_tokens: Iterable[str],
) -> list[dict[str, Any]]:
    """Scan `aria/` Python source files to prove zero hardcoded benchmark task IDs or answers."""
    pkg = Path(aria_package_dir).resolve()
    tokens = [str(t).strip() for t in forbidden_tokens if str(t).strip() and len(str(t).strip()) >= 5]
    hits: list[dict[str, Any]] = []
    if not tokens:
        return hits
    for py_file in sorted(pkg.rglob("*.py")):
        try:
            text = py_file.read_text(encoding="utf-8")
        except OSError:
            continue
        for tok in tokens:
            if tok in text:
                hits.append({
                    "file": py_file.relative_to(pkg.parent).as_posix(),
                    "token": tok,
                })
    return hits
