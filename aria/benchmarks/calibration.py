"""Scientific calibration gate (S8, pre-integration audit).

ARIA's recovery benchmarks already compute recall + empirical FDR on synthetic
ground truth. Their pass/fail thresholds used to live hardcoded and scattered in
each test, so there was no single, versioned statement of the minimum scientific
performance ARIA must keep — a refactor could quietly erode recall and every test
might still be green against its own loose number.

This module makes that an invariant: `calibration_baseline.json` is the single
source of the floors/ceilings, and `check_calibration` compares a benchmark's
measured metrics against its entry and returns a structured pass/fail with the
exact violations. The benchmark run-functions accept `min_recall` /
`max_empirical_fdr`, so callers pull the bounds from here (`baseline_bounds`)
instead of inventing numbers.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_BASELINE_PATH = Path(__file__).resolve().parent / "calibration_baseline.json"


def load_baseline() -> dict[str, Any]:
    """Load the versioned calibration baseline (the single source of bounds)."""
    return json.loads(_BASELINE_PATH.read_text(encoding="utf-8"))


def baseline_bounds(benchmark: str) -> dict[str, float]:
    """Return the {min_recall, max_empirical_fdr, ...} bounds for a benchmark."""
    benchmarks = load_baseline()["benchmarks"]
    if benchmark not in benchmarks:
        raise KeyError(
            f"no calibration baseline for '{benchmark}'; "
            f"known: {sorted(benchmarks)}"
        )
    return benchmarks[benchmark]


def check_calibration(benchmark: str, metrics: dict[str, Any]) -> dict[str, Any]:
    """Compare measured metrics against the versioned baseline for `benchmark`.

    Checks (only those present in both the baseline and the metrics):
      - recall            >= min_recall
      - empirical_fdr     <= max_empirical_fdr
      - false-positive fraction (n_false_positive / max(n_called, 1))
                          <= max_false_positive_fraction

    Returns {benchmark, status: "pass"|"fail", violations: [...], measured,
    baseline}. A measured metric that is missing is itself a violation — a gate
    must not pass on absent evidence.
    """
    bounds = baseline_bounds(benchmark)
    violations: list[str] = []
    measured: dict[str, Any] = {}

    def _need(key: str) -> float | None:
        if key not in metrics or metrics[key] is None:
            violations.append(f"missing metric '{key}'")
            return None
        try:
            return float(metrics[key])
        except (TypeError, ValueError):
            violations.append(f"non-numeric metric '{key}'={metrics[key]!r}")
            return None

    if "min_recall" in bounds:
        recall = _need("recall")
        measured["recall"] = recall
        if recall is not None and recall < bounds["min_recall"]:
            violations.append(
                f"recall {recall:.4f} < min_recall {bounds['min_recall']}")

    if "max_empirical_fdr" in bounds:
        fdr = _need("empirical_fdr")
        measured["empirical_fdr"] = fdr
        if fdr is not None and fdr > bounds["max_empirical_fdr"]:
            violations.append(
                f"empirical_fdr {fdr:.4f} > max_empirical_fdr "
                f"{bounds['max_empirical_fdr']}")

    if "max_false_positive_fraction" in bounds:
        n_fp = metrics.get("n_false_positive")
        n_called = metrics.get("n_called")
        if n_fp is None or n_called is None:
            # Optional check: only enforced when both counts are reported.
            pass
        else:
            frac = float(n_fp) / max(float(n_called), 1.0)
            measured["false_positive_fraction"] = round(frac, 4)
            if frac > bounds["max_false_positive_fraction"]:
                violations.append(
                    f"false-positive fraction {frac:.4f} > "
                    f"{bounds['max_false_positive_fraction']}")

    return {
        "benchmark": benchmark,
        "status": "pass" if not violations else "fail",
        "violations": violations,
        "measured": measured,
        "baseline": {k: v for k, v in bounds.items() if k != "provenance"},
    }
