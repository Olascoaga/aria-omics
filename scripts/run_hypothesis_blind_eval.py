#!/usr/bin/env python
"""H20: run the SPECULATIVE HypothesisAgent blind-evaluation benchmark.

Runs two arms over the four real scenarios with a REAL LLM:
  - governed  : the full HypothesisAgent (all walls H14-H19)
  - ungoverned: the same model, a naive prompt, NO gates (the baseline)

then the red-team battery and the promotion go/no-go gate. Writes the report +
the anonymised blind grading sheet to ``docs/benchmark_results/hypothesis/``. The
human panel fills the sheet; re-run ``evaluate_promotion_gate`` with the scores to
get the final promotion verdict. Deterministic parts (red-team, factuality,
gate logic) are guarded by tests and need no network.

Usage:
    python scripts/run_hypothesis_blind_eval.py [--seed N] [--out DIR]
    python scripts/run_hypothesis_blind_eval.py --redteam-only   # no LLM needed
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from aria.benchmarks.hypothesis_blind_eval import run_blind_eval
from aria.benchmarks.hypothesis_redteam import evaluate_redteam
from aria.version import collect_version_metadata

_OUT_DEFAULT = Path("docs/benchmark_results/hypothesis")


def _public_provenance() -> dict:
    """Version metadata with machine-absolute paths stripped (public artifact).

    ``collect_version_metadata`` records ``environment.conda_prefix`` (an absolute
    home path); a committed, preprint-citable benchmark JSON must stay relocatable
    and leak no author path (guarded by
    ``tests/test_public_artifacts_no_absolute_paths.py``). Keep the env name + lock
    sha256 provenance, drop the absolute prefix.
    """
    prov = collect_version_metadata()
    env = prov.get("environment")
    if isinstance(env, dict):
        env.pop("conda_prefix", None)
    return prov


# A naive baseline system prompt: hypothesis generation with NO governance rules.
_NAIVE_SYSTEM = (
    "You are a creative molecular-biology hypothesis generator. Given a list of "
    "measurements, propose interesting mechanistic hypotheses that connect them. "
    "Return a JSON array of objects with keys id, mechanism, entities, "
    "observation_refs, observed_claims, experiment, devils_advocate."
)


def _build_proposers():
    """Wire a governed and an ungoverned proposer over the real LLMProvider."""
    from aria.llm.provider import LLMProvider
    from aria.agents.narrative.hypothesis import LLMProposer
    from aria.agents.narrative.hypothesis import build_proposer_prompt, parse_hypotheses

    llm = LLMProvider()
    governed = LLMProposer.from_provider(llm)

    def ungoverned(signals, exp_ctx):
        # Same evidence, naive prompt, NO gates downstream.
        from aria.llm.provider import TaskTier
        prompt = build_proposer_prompt(signals, exp_ctx, 4)
        raw = llm.complete(prompt=prompt, system=_NAIVE_SYSTEM,
                           tier=TaskTier.HEAVY, max_tokens=8192)
        return parse_hypotheses(raw)

    return governed, ungoverned


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=_OUT_DEFAULT)
    ap.add_argument("--redteam-only", action="store_true",
                    help="Run only the deterministic red-team battery (no LLM).")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    if args.redteam_only:
        redteam = evaluate_redteam()
        redteam["provenance"] = _public_provenance()
        path = args.out / "hypothesis_redteam.json"
        path.write_text(json.dumps(redteam, indent=2, default=str))
        print(f"red-team: passed={redteam['passed']} "
              f"evasions={redteam['evasions']}/{redteam['n_cases']}")
        print(f"wrote {path}")
        return 0 if redteam["passed"] else 1

    governed, ungoverned = _build_proposers()
    report = run_blind_eval(governed, ungoverned, seed=args.seed)
    report["provenance"] = _public_provenance()

    report_path = args.out / "hypothesis_blind_eval.json"
    sheet_path = args.out / "hypothesis_grading_sheet.json"
    report_path.write_text(json.dumps(report, indent=2, default=str))
    # The blind sheet is written WITHOUT the de-anonymisation key.
    sheet_path.write_text(json.dumps(report["grading_sheet"], indent=2, default=str))

    gate = report["promotion_gate"]
    print(f"governed factuality: {report['governed']['factuality_overall']['rate']}")
    print(f"baseline factuality: {report['ungoverned']['factuality_overall']['rate']}")
    print(f"red-team: passed={report['redteam']['passed']} "
          f"evasions={report['redteam']['evasions']}")
    print(f"gate: mechanical={gate['mechanical']['passed']} "
          f"promotable={gate['promotable']} reason={gate['reason']}")
    print(f"wrote {report_path}")
    print(f"wrote {sheet_path}  (blind — hand to the panel)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
