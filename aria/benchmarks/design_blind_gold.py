"""C3 blind multifactorial design gold — independent complement to B1.

``governance_b1`` scores ARIA's design governance against an ARIA-authored
adversarial corpus (the labels ship in code). That measures internal
consistency but cannot, by construction, rule out the "graded its own homework"
objection for Claim 3.

This module supplies the INDEPENDENT, BLIND half:

  * a frozen HELD-OUT corpus of MULTIFACTORIAL design scenarios whose case ids are
    disjoint from the B1 corpus (ARIA's primitives were not tuned on them);
  * the scenarios ship WITHOUT any decision label — the correct decision is
    authored by an independent human into ``design_gold.csv`` (one row per case),
    never in code (blindness invariant);
  * a scoring boundary that runs ARIA's REAL governance decision
    (``governance_b1.aria_decision``) on each scenario WITHOUT ever reading the
    human gold, then reports the confusion matrix, the same three protocol rates
    as B1 (correct inference / refusal / escalation and the headline unsafe
    execution rate) and ARIA-vs-human agreement (Cohen's kappa).

No fabrication: the gold comes only from the supplied human sheet; this module
never invents a decision and never runs ARIA against a peeked label.
"""

from __future__ import annotations

import csv
import io
from typing import Any, Mapping, Sequence

from aria.benchmarks.b2_annotation import cohen_kappa
from aria.benchmarks.governance_b1 import DesignCase, aria_decision

__all__ = [
    "DECISIONS", "SHEET_COLUMNS", "build_multifactorial_corpus",
    "export_design_sheet", "load_design_gold", "score_blind_design_gold",
]

DECISIONS = ("infer", "escalate", "block")
SHEET_COLUMNS = (
    "case_id", "category", "n_conditions", "n_samples", "factors",
    "batch_structure", "description", "gold_decision",
)


def _factorial(cells: Mapping[tuple[str, str], int], factor2: str,
               *, batch: str | None = None) -> list[dict[str, Any]]:
    """Build samples for a two-factor design.

    ``cells`` maps ``(condition, level2) -> n_replicates``. ``factor2`` is the
    metadata column name of the second factor. ``batch`` selects the batch
    structure: ``"confounded_f2"`` aliases batch to the second factor,
    ``"confounded_cond"`` aliases batch to condition, ``"balanced"`` alternates,
    ``None`` omits batch.
    """
    samples: list[dict[str, Any]] = []
    levels2 = sorted({lvl for _, lvl in cells})
    for (cond, lvl2), n in cells.items():
        for r in range(n):
            s = {
                "sample": f"{cond}_{lvl2}_{r + 1}",
                "condition": cond,
                factor2: lvl2,
                "donor": f"{cond}_{lvl2}_{r + 1}",
                "cell_type": "ctype0",
            }
            if batch == "confounded_f2":
                s["batch"] = f"b_{lvl2}"
            elif batch == "confounded_cond":
                s["batch"] = f"b_{cond}"
            elif batch == "balanced":
                s["batch"] = f"b{(levels2.index(lvl2) + r) % 2}"
            samples.append(s)
    return samples


def _case(case_id: str, category: str, description: str, factors: str,
          batch_structure: str, samples: list[dict[str, Any]],
          covariates: tuple[str, ...]) -> DesignCase:
    # gold is deliberately empty: the decision is authored by an independent
    # human, never in code. aria_decision() never reads it.
    return DesignCase(
        case_id=case_id, category=category, description=description,
        samples=samples, gold="", condition_col="condition",
        covariates=covariates,
    )


def build_multifactorial_corpus() -> list[DesignCase]:
    """Frozen held-out multifactorial design scenarios (blind; no gold in code)."""
    corpus: list[DesignCase] = []

    # Balanced complete 2x2 factorials with adequate replication.
    corpus.append(_case(
        "mf2x2_balanced_n3", "factorial_complete",
        "2x2 condition x genotype, 3 reps/cell, batch balanced across cells",
        "condition x genotype", "balanced",
        _factorial({("treated", "WT"): 3, ("treated", "KO"): 3,
                    ("control", "WT"): 3, ("control", "KO"): 3},
                   "genotype", batch="balanced"),
        covariates=("genotype", "batch")))
    corpus.append(_case(
        "mf3x2_balanced_n3", "factorial_complete",
        "3x2 condition(A,B,C) x timepoint(t0,t1), 3 reps/cell, no batch",
        "condition x timepoint", "none",
        _factorial({("A", "t0"): 3, ("A", "t1"): 3, ("B", "t0"): 3,
                    ("B", "t1"): 3, ("C", "t0"): 3, ("C", "t1"): 3},
                   "timepoint"),
        covariates=("timepoint",)))

    # Low-power complete factorial.
    corpus.append(_case(
        "mf2x2_two_reps", "factorial_low_power",
        "2x2 condition x genotype with only 2 reps/cell (low power)",
        "condition x genotype", "balanced",
        _factorial({("treated", "WT"): 2, ("treated", "KO"): 2,
                    ("control", "WT"): 2, ("control", "KO"): 2},
                   "genotype", batch="balanced"),
        covariates=("genotype", "batch")))

    # Incomplete factorial: a missing cell breaks the interaction contrast.
    corpus.append(_case(
        "mf2x2_missing_cell", "factorial_incomplete",
        "2x2 with the (treated,KO) cell absent (unestimable interaction)",
        "condition x genotype", "none",
        _factorial({("treated", "WT"): 3, ("control", "WT"): 3,
                    ("control", "KO"): 3}, "genotype"),
        covariates=("genotype",)))

    # Singleton cell: one factorial cell has a single replicate.
    corpus.append(_case(
        "mf2x2_singleton_cell", "factorial_under_replicated",
        "2x2 where the (treated,KO) cell has a single replicate",
        "condition x genotype", "none",
        _factorial({("treated", "WT"): 3, ("treated", "KO"): 1,
                    ("control", "WT"): 3, ("control", "KO"): 3}, "genotype"),
        covariates=("genotype",)))

    # Second factor perfectly confounded with batch.
    corpus.append(_case(
        "mf2x2_genotype_batch_confounded", "factorial_confounded",
        "2x2 where genotype is perfectly confounded with processing batch",
        "condition x genotype", "confounded_f2",
        _factorial({("treated", "WT"): 3, ("treated", "KO"): 3,
                    ("control", "WT"): 3, ("control", "KO"): 3},
                   "genotype", batch="confounded_f2"),
        covariates=("genotype", "batch")))

    # Two factors aliased with each other (no independent variation).
    corpus.append(_case(
        "mf_condition_genotype_aliased", "factorial_confounded",
        "genotype aliased with condition (all treated=KO, all control=WT)",
        "condition x genotype", "none",
        _factorial({("treated", "KO"): 4, ("control", "WT"): 4}, "genotype"),
        covariates=("genotype",)))

    # Condition perfectly confounded with batch in a factorial layout.
    corpus.append(_case(
        "mf2x2_condition_batch_confounded", "factorial_confounded",
        "2x2 where condition is perfectly confounded with batch",
        "condition x genotype", "confounded_cond",
        _factorial({("treated", "WT"): 3, ("treated", "KO"): 3,
                    ("control", "WT"): 3, ("control", "KO"): 3},
                   "genotype", batch="confounded_cond"),
        covariates=("genotype", "batch")))

    # Degenerate second factor: only one level (collapses to single factor).
    corpus.append(_case(
        "mf2x1_single_second_level", "factorial_degenerate",
        "declared 2-factor design where the second factor has one level only",
        "condition x genotype", "none",
        _factorial({("treated", "WT"): 3, ("control", "WT"): 3}, "genotype"),
        covariates=("genotype",)))

    # Continuous dose disguised as a categorical second factor.
    dose_samples: list[dict[str, Any]] = []
    for ci, cond in enumerate(("treated", "control")):
        for i in range(4):
            dose_samples.append({
                "sample": f"{cond}_d{i}", "condition": cond,
                "dose": f"{0.5 + i * 0.7:.1f}", "donor": f"{cond}_{i}",
                "cell_type": "ctype0",
            })
    corpus.append(_case(
        "mf_continuous_dose_second_factor", "factorial_continuous",
        "second factor is a continuous dose (each value unique per cell)",
        "condition x dose(continuous)", "none",
        dose_samples, covariates=("dose",)))

    # Nested design: donors nested within condition, adequate replication.
    nested: list[dict[str, Any]] = []
    for ci, cond in enumerate(("treated", "control")):
        for d in range(3):
            for c in range(2):
                nested.append({
                    "sample": f"{cond}_d{d}_c{c}", "condition": cond,
                    "donor": f"{cond}_donor{d}", "cell_type": "ctype0",
                    "batch": f"b{d % 2}",
                })
    corpus.append(_case(
        "mf_nested_donor_replicates", "nested_design",
        "donors nested within condition, 2 pseudo-units per donor, 3 donors/arm",
        "condition / donor(nested)", "balanced",
        nested, covariates=("batch",)))

    return corpus


def export_design_sheet(corpus: Sequence[DesignCase] | None = None) -> str:
    """Write the BLIND labeling sheet (one row per scenario, empty gold column)."""
    corpus = list(corpus) if corpus is not None else build_multifactorial_corpus()
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(SHEET_COLUMNS))
    writer.writeheader()
    for case in corpus:
        n_conditions = len({str(s["condition"]) for s in case.samples})
        writer.writerow({
            "case_id": case.case_id,
            "category": case.category,
            "n_conditions": n_conditions,
            "n_samples": len(case.samples),
            "factors": ";".join(case.covariates),
            "batch_structure": "batch" if any("batch" in s for s in case.samples)
                               else "none",
            "description": case.description,
            "gold_decision": "",
        })
    return buf.getvalue()


def load_design_gold(csv_text: str) -> dict[str, str]:
    """Parse a filled human sheet into ``{case_id: decision}``.

    Only non-empty, in-vocabulary decisions are returned; malformed rows raise so
    a partial or mislabeled gold can never be scored silently.
    """
    reader = csv.DictReader(io.StringIO(csv_text))
    gold: dict[str, str] = {}
    for row in reader:
        case_id = (row.get("case_id") or "").strip()
        decision = (row.get("gold_decision") or "").strip().lower()
        if not case_id or not decision:
            continue
        if decision not in DECISIONS:
            raise ValueError(
                f"case {case_id!r} has gold_decision {decision!r}; "
                f"expected one of {DECISIONS}"
            )
        gold[case_id] = decision
    return gold


def score_blind_design_gold(
    human_gold: Mapping[str, str],
    corpus: Sequence[DesignCase] | None = None,
) -> dict[str, Any]:
    """Score ARIA's blind governance decisions against an independent human gold."""
    corpus = list(corpus) if corpus is not None else build_multifactorial_corpus()
    by_id = {case.case_id: case for case in corpus}

    missing = sorted(set(by_id) - set(human_gold))
    extra = sorted(set(human_gold) - set(by_id))

    confusion = {g: {d: 0 for d in DECISIONS} for g in DECISIONS}
    aria_labels: dict[str, str] = {}
    per_case: list[dict[str, Any]] = []
    for case_id in sorted(set(by_id) & set(human_gold)):
        case = by_id[case_id]
        res = aria_decision(case)          # blind: gold is never passed in
        decision = res["decision"]
        gold = human_gold[case_id]
        aria_labels[case_id] = decision
        confusion[gold][decision] += 1
        per_case.append({
            "case_id": case_id, "category": case.category,
            "human_gold": gold, "aria_decision": decision,
            "agree": gold == decision,
            "card_status": res["card_status"], "vdm_status": res["vdm_status"],
        })

    n_scored = len(per_case)
    n_block_gold = sum(confusion["block"].values())
    n_infer_gold = sum(confusion["infer"].values())
    n_esc_gold = sum(confusion["escalate"].values())
    n_agree = sum(confusion[g][g] for g in DECISIONS)
    kappa = cohen_kappa(dict(human_gold), aria_labels) if aria_labels else {}

    summary = {
        "n_scenarios": len(corpus),
        "n_scored": n_scored,
        "n_agree": n_agree,
        "agreement_rate": round(n_agree / max(n_scored, 1), 4),
        "correct_inference_rate": round(
            confusion["infer"]["infer"] / max(n_infer_gold, 1), 4),
        "correct_refusal_rate": round(
            confusion["block"]["block"] / max(n_block_gold, 1), 4),
        "correct_escalation_rate": round(
            confusion["escalate"]["escalate"] / max(n_esc_gold, 1), 4),
        "unsafe_execution_rate": round(
            confusion["block"]["infer"] / max(n_block_gold, 1), 4),
    }
    status = "pass" if (not missing and not extra and n_scored > 0) else "incomplete"
    return {
        "status": status,
        "benchmark": "C3_blind_multifactorial_design_gold",
        "benchmark_version": "v1",
        "scope": "independent_blind_gold_multifactorial_held_out_designs",
        "method_under_test": (
            "ARIA ScRNAAuditAgent readiness + validate_design_matrix "
            "(governance_b1.aria_decision)"
        ),
        "corpus_disjoint_from_b1": True,
        "gold_authorship": "independent_human_design_gold.csv",
        "summary": summary,
        "cohen_kappa": kappa,
        "confusion_matrix": confusion,
        "per_case": per_case,
        "unscored_scenarios": missing,
        "unknown_gold_case_ids": extra,
        "caveats": [
            "The gold is a single independent human design annotation; ARIA's "
            "decision is computed blind (the gold is never read before scoring).",
            "Gold authorship independence is recorded at run time in provenance; a "
            "self-authored gold must disclose that caveat.",
            "Corpus case ids are disjoint from governance_b1.build_corpus() so "
            "ARIA's primitives were not tuned on these scenarios.",
        ],
    }
