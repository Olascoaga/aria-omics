"""Benchmark B1: DesignAgent governance on an adversarial design corpus.

Benchmark B is ARIA's contribution: governance, not statistical accuracy. B1
measures whether ARIA infers when a design is defensible and blocks/escalates
when it is not. It drives ARIA's REAL governance primitives — the per-modality
readiness audit (`ScRNAAuditAgent`, ADR P2-6) and the design-matrix validator
(`validate_design_matrix`, X7) — over a labelled adversarial corpus, and reports
the three components the protocol requires separately:

- correct inference rate (defensible designs ARIA proceeds on);
- correct refusal rate (indefensible designs ARIA blocks);
- unsafe execution rate (HEADLINE, must be ~0): an indefensible design ARIA
  would run an inferential analysis on.

Decision classes map the real primitive outputs:
- green/ok      -> ``infer``;
- yellow/warning-> ``escalate`` (proceed only with explicit acknowledgement);
- red/blocking  -> ``block`` (removed from inferential dispatch).

The decision is the worst (most conservative) of the two primitives, mirroring
how ARIA gates dispatch. Coverage boundary: the explicit-reference-direction
case (alphabetical denominator, P0-5) is governed by a separate gate
(`tests/test_explicit_contrast_gate.py`) and is represented here through the
incomplete-design path (escalate), not re-implemented.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

__all__ = ["DesignCase", "build_corpus", "aria_decision", "score_b1",
           "run_b1_benchmark", "write_b1_figure"]

_DECISIONS = ("infer", "escalate", "block")
_RANK = {"infer": 0, "escalate": 1, "block": 2}
_CARD_TO_RANK = {"green": 0, "yellow": 1, "red": 2}
_VDM_TO_RANK = {"ok": 0, "warning": 1, "blocking": 2}


@dataclass
class DesignCase:
    case_id: str
    category: str
    description: str
    samples: list[dict[str, Any]]          # each: sample, condition, batch, donor, ...
    gold: str                              # infer | escalate | block
    condition_col: str = "condition"
    covariates: tuple[str, ...] = ()
    pseudobulk_cols: dict[str, str] = field(default_factory=lambda: {
        "condition_col": "condition", "replicate_col": "donor",
        "groupby_col": "cell_type"})
    drop_pseudobulk_cols: tuple[str, ...] = ()


def _exp_context(case: DesignCase) -> dict[str, Any]:
    groups: dict[str, list[str]] = {}
    batch_map: dict[str, str] = {}
    for s in case.samples:
        cond = str(s[case.condition_col])
        groups.setdefault(cond, []).append(str(s["sample"]))
        if "batch" in s and s["batch"] is not None:
            batch_map[str(s["sample"])] = str(s["batch"])
    pb = {k: v for k, v in case.pseudobulk_cols.items()
          if k not in case.drop_pseudobulk_cols}
    design: dict[str, Any] = {"groups": groups, "pseudobulk": pb}
    if batch_map:
        design["batch_map"] = batch_map
    return {"design": design}


def _metadata_df(case: DesignCase):
    import pandas as pd

    df = pd.DataFrame(case.samples)
    if "sample" in df.columns:
        df = df.set_index("sample")
    return df


def aria_decision(case: DesignCase) -> dict[str, Any]:
    """Run ARIA's real readiness audit + design-matrix validator on a case."""
    from aria.agents.modality_audit import ScRNAAuditAgent
    from aria.utils.design_matrix import validate_design_matrix

    card = ScRNAAuditAgent().audit(_exp_context(case))
    vdm = validate_design_matrix(
        _metadata_df(case), case.condition_col, list(case.covariates))

    card_rank = _CARD_TO_RANK.get(card.get("status", "green"), 0)
    vdm_rank = _VDM_TO_RANK.get(vdm.get("status", "ok"), 0)
    worst = max(card_rank, vdm_rank)
    decision = _DECISIONS[worst]
    return {
        "decision": decision,
        "card_status": card.get("status"),
        "vdm_status": vdm.get("status"),
        "dispatch_policy": card.get("dispatch_policy"),
    }


def _case(case_id, category, description, samples, gold, **kw) -> DesignCase:
    return DesignCase(case_id=case_id, category=category, description=description,
                      samples=samples, gold=gold, **kw)


def _grp(conditions: dict[str, int], *, batch=None, donors=True):
    """Build a sample list: {condition: n_reps}. batch="confounded" maps each
    condition to its own batch; batch="balanced" alternates batches across
    conditions; batch=None omits batch."""
    samples = []
    conds = list(conditions)
    for ci, (cond, n) in enumerate(conditions.items()):
        for r in range(n):
            s = {"sample": f"{cond}_{r+1}", "condition": cond,
                 "donor": f"d{ci}_{r+1}", "cell_type": "ctype0"}
            if batch == "confounded":
                s["batch"] = f"b{ci}"
            elif batch == "balanced":
                s["batch"] = f"b{r % 2}"
            samples.append(s)
    return samples


def build_corpus() -> list[DesignCase]:
    """~30 labelled adversarial + defensible design cases."""
    corpus: list[DesignCase] = []

    # 1. Defensible: >=3 replicates, complete design, balanced batch -> infer.
    for i, n in enumerate((3, 4, 5, 6)):
        corpus.append(_case(
            f"defensible_n{n}", "defensible",
            f"{n} balanced replicates/condition, complete design, balanced batch",
            _grp({"treated": n, "control": n}, batch="balanced"), "infer"))
    corpus.append(_case(
        "defensible_3cond", "defensible",
        "three conditions, 3 reps each, balanced batch",
        _grp({"A": 3, "B": 3, "C": 3}, batch="balanced"), "infer"))
    corpus.append(_case(
        "defensible_nobatch", "defensible",
        "3 reps/condition, no batch column (nothing to confound)",
        _grp({"treated": 3, "control": 3}), "infer"))

    # 2. Under-replicated (n<2 in a condition) -> block.
    corpus.append(_case(
        "under_rep_singleton", "under_replicated",
        "one condition has a single replicate",
        _grp({"treated": 1, "control": 3}), "block"))
    corpus.append(_case(
        "under_rep_both_one", "under_replicated",
        "both conditions have a single replicate",
        _grp({"treated": 1, "control": 1}), "block"))
    corpus.append(_case(
        "under_rep_three_groups", "under_replicated",
        "three conditions, one with a single replicate",
        _grp({"A": 3, "B": 3, "C": 1}), "block"))

    # 3. Borderline two replicates -> escalate (low power, ack required).
    for i in range(3):
        corpus.append(_case(
            f"two_reps_{i}", "borderline_two_reps",
            "exactly two replicates per condition (low power)",
            _grp({"treated": 2, "control": 2}, batch="balanced"), "escalate"))

    # 4. Batch perfectly confounded with condition -> block.
    for i in range(4):
        corpus.append(_case(
            f"batch_confounded_{i}", "batch_confounded",
            "each condition maps to a distinct batch",
            _grp({"treated": 3, "control": 3}, batch="confounded"),
            "block", covariates=("batch",)))

    # 5. Incomplete design (missing replicate/groupby column) -> escalate.
    corpus.append(_case(
        "incomplete_no_replicate", "incomplete_design",
        "no biological replicate/donor column confirmed",
        _grp({"treated": 3, "control": 3}), "escalate",
        drop_pseudobulk_cols=("replicate_col",)))
    corpus.append(_case(
        "incomplete_no_groupby", "incomplete_design",
        "no cell-group column confirmed",
        _grp({"treated": 3, "control": 3}), "escalate",
        drop_pseudobulk_cols=("groupby_col",)))
    corpus.append(_case(
        "incomplete_ambiguous_ref", "incomplete_design",
        "labels A/B with no confirmed contrast (ambiguous reference, P0-5)",
        _grp({"A": 3, "B": 3}), "escalate",
        drop_pseudobulk_cols=("replicate_col", "groupby_col")))

    # 6. Continuous covariate disguised as a categorical group -> block.
    cont = [{"sample": f"s{i}", "condition": f"{2.0 + i * 0.5}",
             "donor": f"d{i}", "cell_type": "ctype0"} for i in range(8)]
    corpus.append(_case(
        "continuous_as_categorical", "continuous_as_categorical",
        "condition is a continuous dose (each value unique, 1 sample)",
        cont, "block"))
    cont2 = [{"sample": f"t{i}", "condition": f"{3.1 + i}",
              "donor": f"d{i}", "cell_type": "ctype0"} for i in range(6)]
    corpus.append(_case(
        "continuous_as_categorical_2", "continuous_as_categorical",
        "continuous time disguised as a group", cont2, "block"))

    # 7. Degenerate: single condition level -> block.
    corpus.append(_case(
        "single_level", "degenerate",
        "only one condition level (no contrast possible)",
        _grp({"treated": 5}), "block"))
    corpus.append(_case(
        "single_level_2", "degenerate",
        "all samples share one condition", _grp({"only": 4}), "block"))

    return corpus


def score_b1(corpus: list[DesignCase] | None = None) -> dict[str, Any]:
    """Score ARIA's governance decisions against the adversarial corpus."""
    corpus = corpus if corpus is not None else build_corpus()
    confusion = {g: {d: 0 for d in _DECISIONS} for g in _DECISIONS}
    per_case = []
    for case in corpus:
        res = aria_decision(case)
        confusion[case.gold][res["decision"]] += 1
        per_case.append({
            "case_id": case.case_id, "category": case.category,
            "gold": case.gold, "decision": res["decision"],
            "correct": case.gold == res["decision"],
            "card_status": res["card_status"], "vdm_status": res["vdm_status"],
        })

    n_infer_gold = sum(confusion["infer"].values())
    n_block_gold = sum(confusion["block"].values())
    n_esc_gold = sum(confusion["escalate"].values())
    correct_inference_rate = confusion["infer"]["infer"] / max(n_infer_gold, 1)
    correct_refusal_rate = confusion["block"]["block"] / max(n_block_gold, 1)
    # Headline: indefensible design (gold=block) on which ARIA would still run
    # an inferential analysis (decision=infer).
    unsafe_execution_rate = confusion["block"]["infer"] / max(n_block_gold, 1)
    correct_escalation_rate = confusion["escalate"]["escalate"] / max(n_esc_gold, 1)
    n_correct = sum(confusion[g][g] for g in _DECISIONS)

    summary = {
        "n_cases": len(corpus),
        "n_correct": n_correct,
        "overall_accuracy": round(n_correct / max(len(corpus), 1), 4),
        "correct_inference_rate": round(correct_inference_rate, 4),
        "correct_refusal_rate": round(correct_refusal_rate, 4),
        "correct_escalation_rate": round(correct_escalation_rate, 4),
        "unsafe_execution_rate": round(unsafe_execution_rate, 4),
        "n_block_gold": n_block_gold,
    }
    axis_pass = {
        "unsafe_execution_near_zero": unsafe_execution_rate <= 0.0,
        "refusal": correct_refusal_rate >= 0.9,
        "inference": correct_inference_rate >= 0.9,
    }
    status = "pass" if all(axis_pass.values()) else "fail"
    return {
        "status": status,
        "benchmark": "B1_design_governance",
        "benchmark_version": "v1",
        "scope": "governance_adversarial_design_corpus",
        "method_under_test": "ARIA ScRNAAuditAgent readiness + validate_design_matrix",
        "axis_pass": axis_pass,
        "summary": summary,
        "confusion_matrix": confusion,
        "per_case": per_case,
        "messages": [
            "B1 drives ARIA's real readiness audit + design-matrix validator over "
            "an adversarial design corpus; unsafe execution rate is the headline "
            "and must be ~0."
        ],
    }


def write_b1_figure(manifest: dict[str, Any], path: str) -> str:
    from pathlib import Path

    s = manifest.get("summary", {})
    vals = [
        ("Correct inference", s.get("correct_inference_rate", 0.0)),
        ("Correct refusal", s.get("correct_refusal_rate", 0.0)),
        ("Correct escalation", s.get("correct_escalation_rate", 0.0)),
        ("Unsafe execution (↓)", s.get("unsafe_execution_rate", 0.0)),
    ]
    width, left, bar_w, bar_h, gap, top = 760, 250, 380, 40, 24, 70
    rows = []
    for i, (label, v) in enumerate(vals):
        y = top + i * (bar_h + gap)
        vv = max(0.0, min(1.0, float(v)))
        good = vv <= 0.05 if "Unsafe" in label else vv >= 0.9
        fill = "#2F6F73" if good else "#A23B3B"
        rows.append(
            f'<text x="20" y="{y + 26}" font-size="16" fill="#1f2933">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * vv:.1f}" height="{bar_h}" fill="{fill}"/>'
            f'<text x="{left + bar_w + 16}" y="{y + 26}" font-size="15" fill="#1f2933">{vv:.2f}</text>'
        )
    status = manifest.get("status", "unknown").upper()
    height = top + len(vals) * (bar_h + gap) + 30
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="20" y="34" font-size="20" font-weight="700" fill="#111827">'
        'Fig. 3: B1 DesignAgent governance (adversarial corpus)</text>'
        f'<text x="620" y="34" font-size="16" font-weight="700" fill="#2F6F73">{status}</text>'
        f'<text x="20" y="{height - 12}" font-size="12" fill="#4b5563">'
        f'{s.get("n_cases", 0)} cases; unsafe execution rate = inferential dispatch on '
        f'an indefensible design (target 0).</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_b1_benchmark(
    *,
    output_dir: str | None = None,
    manifest_name: str = "b1_design_governance_v4.5.5.json",
    figure_name: str = "fig3_b1_design_governance_v4.5.5.svg",
) -> dict[str, Any]:
    import json
    from pathlib import Path

    manifest = score_b1()
    if output_dir:
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        figure_path = outdir / figure_name
        manifest_path = outdir / manifest_name
        write_b1_figure(manifest, str(figure_path))
        manifest["artifacts"] = {"manifest_json": str(manifest_path),
                                 "figure_svg": str(figure_path)}
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
