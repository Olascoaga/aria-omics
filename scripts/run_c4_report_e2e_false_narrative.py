#!/usr/bin/env python3
"""C4 report-level false-narrative E2E kit CLI.

Subcommands:
  run     Execute the E2E now: real pyDESeq2 (dispatched to aria-rna-env) ->
          legitimate NarrativeBlocks + real run ledger -> inject the B2 false
          narratives -> ARIA's real compile_public_claims -> render report.html.
          Emit the automated withhold/emit manifest and the BLIND human
          faithfulness sheet (one row per emitted narrative). No human needed.
  score   Combine the automated result with an INDEPENDENT human faithfulness
          gold (report_faithfulness.csv) into the final receipt manifest. The
          human gold is never synthesized.

The `run` step exercises the bus/compiler/report boundary and proves that
injected false narratives are withheld while legitimate claims are emitted; the
`score` step (the freeze lane command) is gated on the independent human gold.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aria.benchmarks import report_false_narrative as rfn  # noqa: E402


def _run(args) -> int:
    de = rfn.run_real_de(args.work_dir, rna_env=args.rna_env)
    agent_results = rfn.de_to_agent_results(de)
    legit, ledger, exp_ctx = rfn.build_legit_blocks_and_ledger(agent_results)
    e2e = rfn.compile_e2e(legit, rfn.false_narrative_blocks(), exp_ctx, ledger)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rfn.render_report_html(e2e, out_dir / "report.html")

    from aria.version import __version__, collect_version_metadata
    auto = {
        "benchmark": "C4_report_e2e_false_narrative_auto",
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "de_summary": {"n_significant": de.get("n_sig"),
                       "n_tested": de.get("n_tested")},
        "result": {k: v for k, v in e2e.items() if k != "emitted_narratives"},
        "emitted_narratives": e2e["emitted_narratives"],
    }
    (out_dir / "false_narrative_auto.json").write_text(
        json.dumps(auto, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    sheet_path = Path(args.sheet_out)
    sheet_path.parent.mkdir(parents=True, exist_ok=True)
    sheet_path.write_text(rfn.export_faithfulness_sheet(e2e), encoding="utf-8")

    print(f"E2E: {e2e['n_false_withheld']}/{e2e['n_false_injected']} false "
          f"withheld, {e2e['n_legit_emitted_with_injection']} legit emitted, "
          f"safe={e2e['safe']}")
    print(f"wrote report.html + false_narrative_auto.json -> {out_dir}")
    print(f"wrote blind faithfulness sheet -> {sheet_path}")
    if not e2e["safe"]:
        print("UNSAFE: a false narrative leaked or legit emission collapsed",
              file=sys.stderr)
        return 1
    return 0


def _score(args) -> int:
    gold_path = Path(args.gold)
    if not gold_path.is_file():
        print(f"human faithfulness gold not found: {gold_path}. Provide an "
              f"independent report_faithfulness.csv (see `run`). Refusing to "
              f"fabricate a gold.", file=sys.stderr)
        return 1
    # Re-run the deterministic E2E so the receipt binds to the current source.
    de = rfn.run_real_de(args.work_dir, rna_env=args.rna_env)
    agent_results = rfn.de_to_agent_results(de)
    legit, ledger, exp_ctx = rfn.build_legit_blocks_and_ledger(agent_results)
    e2e = rfn.compile_e2e(legit, rfn.false_narrative_blocks(), exp_ctx, ledger)

    human_gold = rfn.load_faithfulness_gold(gold_path.read_text(encoding="utf-8"))
    if not human_gold:
        print(f"gold {gold_path} has no filled human_verdict rows.",
              file=sys.stderr)
        return 1
    manifest = rfn.score_against_human_gold(human_gold, e2e)

    from aria.version import __version__, collect_version_metadata
    manifest["aria_version"] = __version__
    manifest["provenance"] = collect_version_metadata()
    manifest["gold_source"] = gold_path.name

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print(f"status={manifest['status']} automated.safe={manifest['automated']['safe']} "
          f"human.flagged_false={manifest['human']['n_flagged_false_by_human']}")
    print(f"wrote {out}")
    return 0 if manifest["status"] == "pass" else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("run", help="execute the E2E and export the blind sheet")
    pr.add_argument("--output-dir",
                    default="docs/benchmark_results/preprint_v1/claim_4/report_e2e")
    pr.add_argument(
        "--sheet-out",
        default="docs/benchmark_results/preprint_v1/human/"
                "report_faithfulness_TEMPLATE.csv")
    pr.add_argument("--work-dir", default="/tmp/c4_report_e2e")
    pr.add_argument("--rna-env", default="aria-rna-env")
    pr.set_defaults(func=_run)

    ps = sub.add_parser("score", help="score against the human faithfulness gold")
    ps.add_argument("gold")
    ps.add_argument(
        "--out",
        default="docs/benchmark_results/preprint_v1/claim_4/"
                "report_e2e/report_e2e_human_gold.json")
    ps.add_argument("--work-dir", default="/tmp/c4_report_e2e")
    ps.add_argument("--rna-env", default="aria-rna-env")
    ps.set_defaults(func=_score)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
