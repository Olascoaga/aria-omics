#!/usr/bin/env python3
"""B2 multi-annotator kit CLI — the bridge from the scaffold to the human gold.

Subcommands:
  export   Write the BLIND labeling sheet (one row per claim, empty label column).
           Give the same sheet to >=2 annotators; they fill `label` independently.
  score    Load the filled annotator sheets, compute inter-annotator agreement
           (Cohen's kappa for 2, Krippendorff's alpha for >=2), adjudicate into a
           final gold (with an optional adjudicator sheet), and write a summary.

No fabrication: agreement/gold come only from the supplied human sheets; this CLI
never invents labels and never runs ARIA governance.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.governance_b2 import corpus_rows_for_labeling  # noqa: E402
from aria.benchmarks.b2_annotation import (  # noqa: E402
    export_labeling_sheet, load_annotations, cohen_kappa, krippendorff_alpha,
    adjudicate,
)


def _export(args) -> int:
    sheet = export_labeling_sheet(corpus_rows_for_labeling())
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(sheet, encoding="utf-8")
    n = sheet.count("\n") - 1
    print(f"wrote blind labeling sheet ({n} claims) -> {out}")
    print("Give this SAME sheet to >=2 annotators; they fill the `label` column "
          "independently against B2_RUBRIC.md, then run `score`.")
    return 0


def _score(args) -> int:
    anns = [load_annotations(Path(p).read_text(encoding="utf-8")) for p in args.sheets]
    overrides = (load_annotations(Path(args.adjudicator).read_text(encoding="utf-8"))
                 if args.adjudicator else None)

    pairwise = {f"{i + 1}v{j + 1}": cohen_kappa(anns[i], anns[j])
                for i, j in combinations(range(len(anns)), 2)}
    alpha = krippendorff_alpha(anns)
    gold = adjudicate(anns, overrides=overrides)
    n_gold = sum(1 for g in gold.values() if g["gold"] is not None)
    n_open = sum(1 for g in gold.values() if g["gold"] is None)

    summary = {
        "n_annotators": len(anns),
        "n_claims_labelled": {f"annotator_{i + 1}": len(a) for i, a in enumerate(anns)},
        "cohen_kappa_pairwise": pairwise,
        "krippendorff_alpha": alpha,
        "n_gold_resolved": n_gold,
        "n_needs_adjudication": n_open,
        "gold": gold,
    }
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(summary, indent=2, sort_keys=True),
                                  encoding="utf-8")
    print(json.dumps({k: summary[k] for k in
                      ("n_annotators", "cohen_kappa_pairwise", "krippendorff_alpha",
                       "n_gold_resolved", "n_needs_adjudication")},
                     indent=2, sort_keys=True))
    if n_open:
        print(f"\n{n_open} claim(s) need adjudication; supply --adjudicator to resolve.")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    pe = sub.add_parser("export", help="write the blind labeling sheet")
    pe.add_argument("--out", default="docs/benchmark_results/b2_claim/b2_labeling_sheet.csv")
    pe.set_defaults(func=_export)
    ps = sub.add_parser("score", help="agreement + adjudication from filled sheets")
    ps.add_argument("sheets", nargs="+", help="filled annotator CSVs (>=2)")
    ps.add_argument("--adjudicator", default=None, help="adjudicator CSV for ties")
    ps.add_argument("--out", default=None, help="write the full summary JSON here")
    ps.set_defaults(func=_score)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
