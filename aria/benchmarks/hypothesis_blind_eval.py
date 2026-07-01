"""H20 blind-evaluation harness for the SPECULATIVE HypothesisAgent (ADR-057).

Reproducible comparison of two arms over real-shaped scenarios:

  - **governed** — the full ``HypothesisAgent`` with every wall (H14-H19): the
    accepted, grounded, gated set.
  - **ungoverned** — the SAME LLM proposer output with NO gates: every candidate
    it emits (the "LLM without governance" baseline).

The harness computes the MECHANICAL, deterministic metrics (factuality =
invention rate, stability = run-to-run top-k overlap, plus the red-team battery),
anonymises both arms into a blind grading sheet for the human panel, and evaluates
the promotion go/no-go gate. Human scoring (plausibility / novelty / experimental
utility, 1-5) is PROCESS, not code: the harness produces the blind sheet and folds
the returned scores into the gate; it never fabricates them.

The tier stays experimental / opt-in / non-promotable until this gate passes
(zero invention in the governed set, governance beats the baseline, zero red-team
evasions, acceptable stability, AND a passing human panel).
"""

from __future__ import annotations

import hashlib
import random
import statistics
from typing import Any, Callable

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import EvidenceSignal
from aria.agents.narrative.hypothesis.grounding import verify_hypothesis_grounding
from aria.agents.narrative.hypothesis.types import Hypothesis
from aria.benchmarks.hypothesis_redteam import evaluate_redteam

# A proposer: audited evidence -> candidate hypotheses (the LLM lives here).
Proposer = Callable[[list, dict], list]

GATE_THRESHOLDS = {
    "governed_factuality_min": 1.0,  # zero invention in the accepted set
    "stability_min": 0.5,            # top-k entity-signature Jaccard across reruns
    "human_min": 3.5,                # median panel score (1-5) on the human axes
    "top_k": 3,
}


# ── scenarios (real-shaped, one per single-modality adapter) ─────────────────

def _sig(entity, kind, modality, measure, node, value, direction, ctx, caveats):
    return EvidenceSignal(
        entity=entity, entity_kind=kind, modality=modality, measure=measure,
        audited_node_ref=node, value=value, direction=direction, context=ctx,
        caveats_inherited=list(caveats),
    )


def scenarios() -> dict[str, dict]:
    """The four real-validated evidence scenarios, one per modality."""
    de = "ledger://bulk/differential_expression"
    pb = "ledger://scRNA/pseudobulk_de"
    mo = "ledger://chromatin/motif_enrichment"
    return {
        "senescence_bulk_rna": {
            "run_ledger": {"entries": [{"node_id": de, "status": "ran"}]},
            "exp_ctx": {"biological_question": "drivers of replicative senescence?"},
            "signals": [
                _sig("CDKN1A", "gene", "bulk_RNA", "log2fc", de, 2.4, "up",
                     "senescent_vs_control", ["low_replication"]),
                _sig("LMNB1", "gene", "bulk_RNA", "log2fc", de, -2.1, "down",
                     "senescent_vs_control", ["low_replication"]),
                _sig("IL6", "gene", "bulk_RNA", "log2fc", de, 1.8, "up",
                     "senescent_vs_control", ["low_replication"]),
            ],
        },
        "cd8_exhaustion_scrna": {
            "run_ledger": {"entries": [{"node_id": pb, "status": "ran"}]},
            "exp_ctx": {"biological_question": "CD8 T-cell exhaustion program?"},
            "signals": [
                _sig("TOX", "gene", "scRNA", "log2fc", pb, 1.9, "up",
                     "exhausted_vs_effector", ["batch"]),
                _sig("PDCD1", "gene", "scRNA", "log2fc", pb, 2.2, "up",
                     "exhausted_vs_effector", ["batch"]),
                _sig("TCF7", "gene", "scRNA", "log2fc", pb, -1.7, "down",
                     "exhausted_vs_effector", ["batch"]),
            ],
        },
        "k562_gm12878_bulk_atac": {
            "run_ledger": {"entries": [{"node_id": mo, "status": "ran"}]},
            "exp_ctx": {"biological_question": "lineage TF programs, K562 vs GM12878?"},
            "signals": [
                _sig("KLF1", "tf_motif", "bulk_ATAC", "motif_enrich", mo, 2.6, "up",
                     "K562_vs_GM12878", ["motif_not_binding", "low_replication"]),
                _sig("GATA1", "tf_motif", "bulk_ATAC", "motif_enrich", mo, 2.1, "up",
                     "K562_vs_GM12878", ["motif_not_binding", "low_replication"]),
                _sig("SPI1", "tf_motif", "bulk_ATAC", "motif_enrich", mo, 2.3, "up",
                     "GM12878_vs_K562", ["motif_not_binding", "low_replication"]),
            ],
        },
        "exhaustion_scatac": {
            "run_ledger": {"entries": [{"node_id": mo, "status": "ran"}]},
            "exp_ctx": {"biological_question": "chromatin drivers of CD8 exhaustion?"},
            "signals": [
                _sig("TOX", "tf_motif", "scATAC", "motif_enrich", mo, 2.0, "up",
                     "exhausted_cluster", ["motif_not_binding"]),
                _sig("NR4A1", "tf_motif", "scATAC", "motif_enrich", mo, 1.8, "up",
                     "exhausted_cluster", ["motif_not_binding"]),
            ],
        },
    }


# ── arms ─────────────────────────────────────────────────────────────────────

def run_governed(proposer: Proposer, scenario: dict) -> list[dict]:
    """The governed arm: the accepted, gated, grounded hypothesis set."""
    agent = HypothesisAgent(proposer=proposer)
    out = agent.generate(
        scenario["signals"], scenario.get("run_ledger"),
        scenario.get("exp_ctx"),
        w_claim_passed=True, w_ledger_passed=True,
    )
    return out["hypotheses"]


def run_ungoverned(proposer: Proposer, scenario: dict) -> list[dict]:
    """The baseline arm: EVERY candidate the proposer emits, no gate applied."""
    candidates = proposer(scenario["signals"], scenario.get("exp_ctx") or {}) or []
    return [
        c.to_dict() if isinstance(c, Hypothesis) else dict(c)
        for c in candidates
    ]


# ── mechanical metrics ───────────────────────────────────────────────────────

def _as_hypothesis(hyp: Any) -> Hypothesis:
    return hyp if isinstance(hyp, Hypothesis) else Hypothesis.from_dict(dict(hyp))


def factuality(hyps: list, signals: list, run_ledger: dict | None) -> dict:
    """Invention rate: the fraction of hypotheses that fabricate NO fact.

    A hypothesis is factual iff it names no entity absent from the audited
    evidence (structured OR prose), cites no unknown/misattributed signal, and
    states no direction that contradicts an audited signal. This isolates
    INVENTION from falsifiability/hedging so it scores both arms on the same axis,
    reusing the grounding verifier's invention checks only.
    """
    total = 0
    factual = 0
    invented: list[dict] = []
    for raw in hyps or []:
        total += 1
        hyp = _as_hypothesis(raw)
        g = verify_hypothesis_grounding(hyp, signals, run_ledger)
        fabricated = (
            list(g.missing_entities)
            + list(g.ungrounded_prose_entities)
            + list(g.unknown_signals)
            + list(g.misattributed_signals)
        )
        contradicts = bool(g.contradicting_claims)
        if not fabricated and not contradicts:
            factual += 1
        else:
            invented.append({
                "id": getattr(hyp, "id", None),
                "fabricated_entities": sorted(set(fabricated)),
                "contradicts": contradicts,
            })
    return {
        "n": total,
        "factual": factual,
        "rate": (factual / total) if total else 1.0,
        "invented": invented,
    }


def _top_k_signature(hyps: list, k: int) -> frozenset:
    """A stable set signature of the top-k hypotheses (entity sets)."""
    sig = []
    for raw in list(hyps)[:k]:
        hyp = _as_hypothesis(raw)
        sig.append(frozenset(str(e).strip().lower() for e in (hyp.entities or [])))
    return frozenset(sig)


def _jaccard(a: frozenset, b: frozenset) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    return len(a & b) / len(union) if union else 1.0


def stability(run_fn: Callable[[], list], k: int = 3, reruns: int = 3) -> dict:
    """Mean pairwise top-k overlap across reruns (run-to-run consistency).

    Deterministic proposers score 1.0; a real stochastic LLM reveals how stable
    the gated set is across generations. Advisory: it measures consistency, never
    correctness.
    """
    runs = [_top_k_signature(run_fn(), k) for _ in range(max(2, reruns))]
    pairs = [
        _jaccard(runs[i], runs[j])
        for i in range(len(runs))
        for j in range(i + 1, len(runs))
    ]
    return {
        "reruns": len(runs),
        "top_k": k,
        "mean_pairwise_jaccard": round(statistics.mean(pairs), 4) if pairs else 1.0,
    }


# ── blind grading sheet ──────────────────────────────────────────────────────

def _blind_id(scenario: str, arm: str, idx: int, seed: int) -> str:
    basis = f"{seed}|{scenario}|{arm}|{idx}"
    return "hyp_" + hashlib.sha1(basis.encode("utf-8")).hexdigest()[:10]


def build_grading_sheet(
    arm_outputs: dict[str, dict[str, list]], seed: int = 0
) -> dict:
    """Anonymise both arms into a shuffled, arm-blind grading sheet + a hidden key.

    ``arm_outputs`` is ``{scenario: {"governed": [...], "ungoverned": [...]}}``.
    The sheet carries only the readable hypothesis content under an opaque
    ``blind_id``; the ``key`` (kept separate, NEVER shown to graders) maps each
    ``blind_id`` back to its ``(scenario, arm)`` so scores can be de-anonymised
    after grading. Deterministic given ``seed``.
    """
    rows: list[dict] = []
    key: dict[str, dict] = {}
    for scenario, arms in sorted(arm_outputs.items()):
        for arm, hyps in sorted(arms.items()):
            for idx, raw in enumerate(hyps or []):
                hyp = _as_hypothesis(raw)
                bid = _blind_id(scenario, arm, idx, seed)
                exp = hyp.experiment.to_dict()
                rows.append({
                    "blind_id": bid,
                    "mechanism": hyp.mechanism,
                    "entities": list(hyp.entities or []),
                    "experiment": {
                        "perturbation": exp.get("perturbation"),
                        "readout": exp.get("readout"),
                        "predicted_direction": exp.get("predicted_direction"),
                        "refuting_outcome": exp.get("refuting_outcome"),
                    },
                    "devils_advocate": hyp.devils_advocate,
                })
                key[bid] = {"scenario": scenario, "arm": arm}
    rng = random.Random(seed)
    rng.shuffle(rows)
    return {
        "schema": "aria.hypothesis_blind_sheet.v1",
        "seed": seed,
        "rubric": {
            "scale": "1-5",
            "axes": ["plausibility", "novelty", "experimental_utility"],
            "note": "Grade blind. Do not attempt to infer which system produced a row.",
        },
        "rows": rows,
        "key": key,  # SEPARATE from the sheet shown to graders.
    }


# ── promotion go/no-go gate ──────────────────────────────────────────────────

def _median_by_axis(human_scores: dict | None, key: dict, arm: str) -> dict:
    """Median human score per axis for a given arm, using the de-anonymisation key."""
    if not human_scores:
        return {}
    axes = ("plausibility", "novelty", "experimental_utility")
    collected: dict[str, list[float]] = {a: [] for a in axes}
    for bid, scores in human_scores.items():
        if (key.get(bid) or {}).get("arm") != arm:
            continue
        for a in axes:
            if isinstance(scores.get(a), (int, float)):
                collected[a].append(float(scores[a]))
    return {a: round(statistics.median(v), 2) for a, v in collected.items() if v}


def evaluate_promotion_gate(
    report: dict,
    human_scores: dict | None = None,
    thresholds: dict | None = None,
) -> dict:
    """Decide go/no-go for promoting the SPECULATIVE tier from a benchmark report.

    Mechanical criteria (hard, code-checked): the governed set invents nothing,
    governance beats the ungoverned baseline on factuality, the red-team battery
    has zero evasions, and the governed set is stable enough. The human criterion
    (median plausibility + experimental utility over the panel) folds in only when
    scores are supplied. The tier is promotable ONLY when mechanical AND human
    both pass; absent human scores, the gate reports the mechanical verdict and
    stays non-promotable pending review.
    """
    th = {**GATE_THRESHOLDS, **(thresholds or {})}
    gov = report["governed"]["factuality_overall"]["rate"]
    base = report["ungoverned"]["factuality_overall"]["rate"]
    evasions = report["redteam"]["evasions"]
    stab = report["governed"]["stability_overall"]["mean_pairwise_jaccard"]

    criteria = {
        "governed_zero_invention": gov >= th["governed_factuality_min"],
        "governance_beats_baseline": gov > base,
        "redteam_zero_evasions": evasions == 0,
        "stability_ok": stab >= th["stability_min"],
    }
    mechanical_passed = all(criteria.values())

    human_passed: bool | None = None
    human_medians: dict = {}
    if human_scores:
        human_medians = _median_by_axis(
            human_scores, report.get("grading_key", {}), "governed"
        )
        human_passed = bool(human_medians) and all(
            human_medians.get(a, 0.0) >= th["human_min"]
            for a in ("plausibility", "experimental_utility")
        )

    if human_passed is None:
        reason = "awaiting_human_review"
    elif not mechanical_passed:
        reason = "mechanical_criteria_failed"
    elif not human_passed:
        reason = "human_panel_below_threshold"
    else:
        reason = "all_criteria_passed"

    return {
        "schema": "aria.hypothesis_promotion_gate.v1",
        "thresholds": th,
        "mechanical": {"passed": mechanical_passed, "criteria": criteria,
                       "governed_factuality": gov, "baseline_factuality": base,
                       "redteam_evasions": evasions, "governed_stability": stab},
        "human": {"passed": human_passed, "medians": human_medians},
        "promotable": bool(mechanical_passed and human_passed),
        "reason": reason,
    }


# ── full run ─────────────────────────────────────────────────────────────────

def run_blind_eval(
    governed_proposer: Proposer,
    ungoverned_proposer: Proposer,
    *,
    seed: int = 0,
    stability_reruns: int = 3,
    thresholds: dict | None = None,
) -> dict:
    """Run both arms over all scenarios, the red-team battery, and the gate.

    ``governed_proposer`` and ``ungoverned_proposer`` are injected so the harness
    is testable with deterministic fakes and driven by a real ``LLMProposer`` from
    the CLI. Returns the full report + the anonymised grading sheet.
    """
    th = {**GATE_THRESHOLDS, **(thresholds or {})}
    scens = scenarios()
    arm_outputs: dict[str, dict[str, list]] = {}
    per_scenario: dict[str, dict] = {}

    gov_fact_tot = {"n": 0, "factual": 0}
    base_fact_tot = {"n": 0, "factual": 0}
    gov_stab: list[float] = []

    for name, scen in scens.items():
        governed = run_governed(governed_proposer, scen)
        ungoverned = run_ungoverned(ungoverned_proposer, scen)
        arm_outputs[name] = {"governed": governed, "ungoverned": ungoverned}

        gov_f = factuality(governed, scen["signals"], scen.get("run_ledger"))
        base_f = factuality(ungoverned, scen["signals"], scen.get("run_ledger"))
        stab = stability(
            lambda s=scen: run_governed(governed_proposer, s),
            k=th["top_k"], reruns=stability_reruns,
        )
        gov_fact_tot["n"] += gov_f["n"]
        gov_fact_tot["factual"] += gov_f["factual"]
        base_fact_tot["n"] += base_f["n"]
        base_fact_tot["factual"] += base_f["factual"]
        gov_stab.append(stab["mean_pairwise_jaccard"])
        per_scenario[name] = {
            "governed_factuality": gov_f,
            "ungoverned_factuality": base_f,
            "governed_stability": stab,
            "n_governed": len(governed),
            "n_ungoverned": len(ungoverned),
        }

    def _rate(t):
        return (t["factual"] / t["n"]) if t["n"] else 1.0

    sheet = build_grading_sheet(arm_outputs, seed=seed)
    redteam = evaluate_redteam()

    report = {
        "schema": "aria.hypothesis_blind_eval.v1",
        "seed": seed,
        "scenarios": sorted(scens),
        "per_scenario": per_scenario,
        "governed": {
            "factuality_overall": {**gov_fact_tot, "rate": _rate(gov_fact_tot)},
            "stability_overall": {
                "mean_pairwise_jaccard": round(
                    statistics.mean(gov_stab), 4) if gov_stab else 1.0
            },
        },
        "ungoverned": {
            "factuality_overall": {**base_fact_tot, "rate": _rate(base_fact_tot)},
        },
        "redteam": redteam,
        "grading_key": sheet["key"],
    }
    gate = evaluate_promotion_gate(report, human_scores=None, thresholds=th)
    report["promotion_gate"] = gate
    report["grading_sheet"] = {k: v for k, v in sheet.items() if k != "key"}
    return report
