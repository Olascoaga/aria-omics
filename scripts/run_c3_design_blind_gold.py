#!/usr/bin/env python3
"""C3 blind multifactorial design-gold kit CLI — the bridge to the human gold.

Subcommands:
  export   Write the BLIND labeling sheet (one row per held-out multifactorial
           design scenario, empty ``gold_decision`` column). An independent human
           fills each ``gold_decision`` in {infer, escalate, block} from design
           principles, WITHOUT seeing ARIA's output, and saves it as
           ``design_gold.csv``.
  score    Load the filled human ``design_gold.csv``, run ARIA's REAL governance
           decision blind on each scenario, and write a confusion-matrix /
           agreement (Cohen's kappa) manifest with provenance.

No fabrication: the gold comes only from the supplied human sheet; this CLI never
invents a decision and never runs ARIA against a peeked label. The scenarios are
held out from the B1 corpus so ARIA's primitives were not tuned on them.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aria.benchmarks.design_blind_gold import (  # noqa: E402
    export_design_sheet, load_design_gold, score_blind_design_gold,
)


def _export(args) -> int:
    sheet = export_design_sheet()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(sheet, encoding="utf-8")
    n = sheet.count("\n") - 1
    print(f"wrote blind design-gold sheet ({n} scenarios) -> {out}")
    print("An INDEPENDENT human fills `gold_decision` in {infer, escalate, block} "
          "per scenario WITHOUT seeing ARIA output, saves as design_gold.csv, "
          "then run `score`.")
    return 0


def _score(args) -> int:
    gold_path = Path(args.gold)
    if not gold_path.is_file():
        print(f"design gold not found: {gold_path}. Provide an independent human "
              f"design_gold.csv (see `export`). Refusing to fabricate a gold.",
              file=sys.stderr)
        return 1
    human_gold = load_design_gold(gold_path.read_text(encoding="utf-8"))
    if not human_gold:
        print(f"design gold {gold_path} has no filled gold_decision rows.",
              file=sys.stderr)
        return 1

    manifest = score_blind_design_gold(human_gold)

    from aria.version import __version__, collect_version_metadata
    manifest["aria_version"] = __version__
    manifest["provenance"] = collect_version_metadata()
    manifest["gold_source"] = gold_path.name

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    s = manifest["summary"]
    print(f"scored {s['n_scored']}/{s['n_scenarios']} scenarios: "
          f"agreement={s['agreement_rate']} unsafe_execution_rate="
          f"{s['unsafe_execution_rate']} kappa={manifest['cohen_kappa'].get('kappa')}")
    print(f"wrote {out}")
    if manifest["status"] != "pass":
        print(f"status={manifest['status']} (unscored={manifest['unscored_scenarios']} "
              f"unknown={manifest['unknown_gold_case_ids']})", file=sys.stderr)
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    pe = sub.add_parser("export", help="write the blind labeling sheet")
    pe.add_argument(
        "--out",
        default="docs/benchmark_results/preprint_v1/human/design_gold_TEMPLATE.csv",
    )
    pe.set_defaults(func=_export)

    ps = sub.add_parser("score", help="score ARIA blind vs the human design gold")
    ps.add_argument("gold", help="path to the filled human design_gold.csv")
    ps.add_argument(
        "--out",
        default="docs/benchmark_results/preprint_v1/claim_3/blind_design_gold.json",
    )
    ps.set_defaults(func=_score)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
