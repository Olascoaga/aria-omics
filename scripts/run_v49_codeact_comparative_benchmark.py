#!/usr/bin/env python3
"""Run the v4.9 40-capsule governed vs. uncontracted CodeAct comparative benchmark.

Emits the JSON receipt and the vignette:

  docs/benchmark_results/v49_comparative_benchmark/comparative_benchmark_v49.json
  docs/vignettes/13_governed_vs_uncontracted_benchmark.md

Usage (scientific stack; the narrative C4/C7 claim gate is dispatched to the
``aria-env`` stack automatically when not importable here):

  conda run -n aria-rna-env python scripts/run_v49_codeact_comparative_benchmark.py
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aria.benchmarks.codeact_comparative_suite import (  # noqa: E402
    ROBUSTNESS_SWEEP_SEEDS,
    render_vignette_markdown,
    run_codeact_comparative_suite,
    run_seed_robustness_sweep,
)
from aria.utils.provenance import hash_params  # noqa: E402

DEFAULT_RECEIPT = "docs/benchmark_results/v49_comparative_benchmark/comparative_benchmark_v49.json"
DEFAULT_VIGNETTE = "docs/vignettes/13_governed_vs_uncontracted_benchmark.md"
_ABSOLUTE_PREFIXES = ("/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/")


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _assert_portable(text: str, label: str) -> None:
    for prefix in _ABSOLUTE_PREFIXES:
        if f'"{prefix}' in text or f"`{prefix}" in text or f" {prefix}" in text:
            raise SystemExit(f"FAIL: machine-absolute path prefix {prefix!r} leaked into {label}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--receipt", default=DEFAULT_RECEIPT, help="Repo-relative JSON receipt path.")
    parser.add_argument("--vignette", default=DEFAULT_VIGNETTE, help="Repo-relative vignette markdown path.")
    parser.add_argument(
        "--sweep-seeds", type=int, default=len(ROBUSTNESS_SWEEP_SEEDS),
        help="Number of alternative base seeds for the robustness sweep (0 disables).",
    )
    args = parser.parse_args()

    receipt = run_codeact_comparative_suite()
    if args.sweep_seeds > 0:
        sweep = run_seed_robustness_sweep(ROBUSTNESS_SWEEP_SEEDS[: args.sweep_seeds])
        sweep_runtime = sweep.pop("runtime", {})
        sweep["sweep_sha256"] = hash_params(sweep)
        receipt["seed_robustness_sweep"] = sweep
        receipt["runtime"]["sweep_claim_gate_execution"] = sweep_runtime.get("claim_gate_execution")

    receipt_text = json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n"
    _assert_portable(receipt_text, args.receipt)
    receipt_link = Path(os.path.relpath(REPO_ROOT / args.receipt, (REPO_ROOT / args.vignette).parent)).as_posix()
    vignette_text = render_vignette_markdown(
        receipt, receipt_relpath=args.receipt, receipt_link=receipt_link,
    )
    _assert_portable(vignette_text, args.vignette)

    _atomic_write(REPO_ROOT / args.receipt, receipt_text)
    _atomic_write(REPO_ROOT / args.vignette, vignette_text)

    m = receipt["metrics"]
    conc = receipt["concordance"]["aria_locked_manifest_vs_aria_governed"]
    print(f"Receipt : {args.receipt}")
    print(f"Vignette: {args.vignette}")
    print(f"claim gate: {receipt['runtime']['claim_gate_execution']}")
    for arm in receipt["arms"]:
        row = m[arm]
        print(
            f"  {arm:32s} silent_error={row['silent_error_rate']:.3f} "
            f"cond_acc={row['conditional_accuracy_on_supported']} "
            f"abst_P={row['abstention_precision']} abst_R={row['abstention_recall']}"
        )
    print(f"  locked-manifest decision concordance: {conc['n_identical_decisions']}/{conc['n_tasks']}")
    if "seed_robustness_sweep" in receipt:
        sm = receipt["seed_robustness_sweep"]["summary"]
        print(
            f"  sweep ({receipt['seed_robustness_sweep']['n_seeds']} seeds): governed SER mean="
            f"{sm['aria_governed_silent_error_rate']['mean']} max={sm['aria_governed_silent_error_rate']['max']}; "
            f"baseline SER mean={sm['baseline_silent_error_rate']['mean']}"
        )

    gov = m["aria_governed"]
    ok = (
        gov["silent_error_rate"] == 0.0
        and gov["conditional_accuracy_on_supported"] == 1.0
        and gov["abstention_precision"] == 1.0
        and gov["abstention_recall"] == 1.0
        and conc["decision_concordance"] == 1.0
    )
    if not ok:
        print("FAIL: governed arm did not meet the frozen-seed acceptance gate", file=sys.stderr)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
