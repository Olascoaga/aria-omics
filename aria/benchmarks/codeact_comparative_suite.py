"""40-Capsule Governed vs. Uncontracted CodeAct Comparative Benchmark (ADR-083 Tier 1).

Builds ``40`` deterministic, seeded multi-omics benchmark capsules (``8`` capsules x
``5`` archetypes) with a known ground truth and scores four analysis arms through
the full-denominator adapter (``aria.benchmarks.adapter``):

1. ``aria_governed`` — live ARIA governance and donor-aware primitives:
   ``resolve_design_contract`` + ``validate_design_matrix`` (design / confounding /
   replicate gates), ``compute_lobo_jackknife`` / ``compute_counts_lobo_jackknife``
   (biological-unit-enforced, Smyth-moderated inference + ``bh_adjust``),
   ``resolve_causal_estimand_license`` + ``compile_public_claims`` (C4 causal
   licensing), ``assert_no_speculative_promotion`` (C7 quarantine wall),
   ``detect_unsupported_inputs`` (raw sequencer boundary),
   ``enforce_same_cell_integration_gate`` (V48 same-cell gate) and
   ``execute_derived_query`` (typed derived-query AST).
2. ``aria_locked_manifest`` — every capsule is first committed to a
   ``DecisionManifest`` (``build_decision_manifest``), serialized, re-loaded through
   ``load_decision_manifest`` (``verify_decision_manifest``) and replayed with the
   analysis policy read ONLY from the verified manifest. A tampered copy of every
   manifest must fail verification.
3. ``aria_zero_llm_deterministic`` — capsules are regenerated from their seeds
   inside an ``isolated_benchmark_environment`` and solved with ``llm=None``
   (deterministic template proposals instead of a proposer).
4. ``uncontracted_codeact_baseline`` — a deterministic, scripted emulation of
   common unguarded analysis code: cell-level ``scipy.stats.ranksums``, Welch tests
   that ignore batch confounding and technical-replicate structure, unvalidated
   ``pd.merge``, positional cell pairing and causal prose on correlation.

The decision path of every ARIA arm is LLM-free by construction (the LLM may only
draft claim prose; code licenses it). No live LLM is called in any arm.

The narrative claim compiler lives in the ``aria-env`` (agents) stack. When this
module runs in a scientific stack without that dependency set (e.g.
``aria-rna-env``), the epistemic gate is executed in a stack-isolated subprocess
on the narrative interpreter via JSON IPC; the results are identical.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Callable

import numpy as np
import pandas as pd

from aria.benchmarks.adapter import (
    BenchmarkPrediction,
    BenchmarkTaskSpec,
    evaluate_benchmark_predictions,
    isolated_benchmark_environment,
)
from aria.benchmarks.derived_query_ast import (
    DerivedQueryPlan,
    DerivedQueryViolation,
    FilterOp,
    GroupSummaryOp,
    ProjectOp,
    RankOp,
    TableContract,
    ValidatedJoinOp,
    execute_derived_query,
)
from aria.utils.decision_manifest import (
    build_decision_manifest,
    load_decision_manifest,
    verify_decision_manifest,
)
from aria.utils.design_matrix import (
    factors_confounded_with_condition,
    resolve_design_contract,
    validate_design_matrix,
)
from aria.utils.lobo_jackknife import (
    bh_adjust,
    compute_counts_lobo_jackknife,
    compute_lobo_jackknife,
)
from aria.utils.multiome_contracts import enforce_same_cell_integration_gate
from aria.utils.provenance import hash_params
from aria.utils.sequencer_formats import detect_unsupported_inputs

SCHEMA_VERSION = "aria.codeact_comparative_benchmark.v1"
SUITE_NAME = "aria_v49_codeact_comparative_40"
SUITE_BASE_SEED = 490_083
CAPSULES_PER_ARCHETYPE = 8

ARCHETYPES: tuple[str, ...] = (
    "design_and_confounding_governance",
    "donor_pseudobulk_vs_cell_pseudoreplication",
    "epistemic_claim_and_causal_licensing",
    "sequencer_and_modality_boundary_abstention",
    "scrubbed_id_derived_query_reasoning",
)

ARMS: tuple[str, ...] = (
    "aria_governed",
    "aria_locked_manifest",
    "aria_zero_llm_deterministic",
    "uncontracted_codeact_baseline",
)

GOVERNED_POLICY: dict[str, Any] = {
    "alpha": 0.05,
    "bulk_lfc_min_log2": 1.0,
    "sc_lfc_min_log2": 0.5,
    "min_replicates_per_condition": 2,
    "biological_unit_aggregation": "mean_log2_over_technical_replicates",
    "single_cell_unit": "donor_pseudobulk_sum_counts",
    "single_cell_normalization": "deseq2_median_of_ratios_size_factors",
    "single_cell_test": "smyth_moderated_t_on_log2cpm_donor_pseudobulk",
    "bulk_test": "ols_condition_plus_confirmed_covariates_smyth_moderated",
    "recovery_min_precision": 0.9,
    "recovery_min_recall_floor": 0.5,
}

NARRATIVE_STACK_ENV = "aria-env"
NARRATIVE_PYTHON_ENV_VAR = "ARIA_NARRATIVE_PYTHON"
_IPC_SENTINEL = "\n@@ARIA_CLAIM_GATE_IPC@@\n"

_REPO_ROOT = Path(__file__).resolve().parents[2]


# ─────────────────────────────────────────────────────────────────────────────
# Capsule model
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class Capsule:
    capsule_id: str
    archetype: str
    variant: str
    seed: int
    question: str
    should_abstain: bool
    supported_by_native_contract: bool
    expected: dict[str, Any]
    design: dict[str, Any] = field(default_factory=dict)
    modalities: dict[str, list[str]] = field(default_factory=dict)
    payload: dict[str, Any] = field(default_factory=dict)

    def public_record(self) -> dict[str, Any]:
        return {
            "capsule_id": self.capsule_id,
            "archetype": self.archetype,
            "variant": self.variant,
            "seed": self.seed,
            "question": self.question,
            "should_abstain": self.should_abstain,
            "supported_by_native_contract": self.supported_by_native_contract,
            "expected": self.expected,
            "input_fingerprint": payload_fingerprint(self.payload),
        }


def _feature_ids(n: int, prefix: str = "FEAT") -> list[str]:
    return [f"{prefix}_{i:04d}" for i in range(1, n + 1)]


def _capsule_seed(archetype_idx: int, k: int, base_seed: int) -> int:
    return int(base_seed) + 1000 * (archetype_idx + 1) + k


def _digest_obj(obj: Any, h: "hashlib._Hash") -> None:
    if isinstance(obj, pd.DataFrame):
        h.update(b"DF")
        h.update(obj.to_csv(float_format="%.10g").encode("utf-8"))
    elif isinstance(obj, pd.Series):
        h.update(b"SR")
        h.update(obj.to_csv(float_format="%.10g").encode("utf-8"))
    elif isinstance(obj, np.ndarray):
        h.update(b"ND" + str(obj.shape).encode() + str(obj.dtype).encode())
        h.update(np.ascontiguousarray(obj).tobytes())
    elif isinstance(obj, dict):
        h.update(b"{")
        for key in sorted(obj, key=str):
            h.update(str(key).encode("utf-8"))
            _digest_obj(obj[key], h)
        h.update(b"}")
    elif isinstance(obj, (list, tuple)):
        h.update(b"[")
        for item in obj:
            _digest_obj(item, h)
        h.update(b"]")
    else:
        h.update(json.dumps(obj, sort_keys=True, default=str).encode("utf-8"))


def payload_fingerprint(payload: dict[str, Any]) -> str:
    """Deterministic SHA-256 over a capsule's synthetic input payload."""
    h = hashlib.sha256()
    _digest_obj(payload, h)
    return h.hexdigest()


def _round(value: Any, nd: int = 6) -> Any:
    if isinstance(value, (float, np.floating)):
        v = float(value)
        return round(v, nd) if np.isfinite(v) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, dict):
        return {str(k): _round(v, nd) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round(v, nd) for v in value]
    return value


# ─────────────────────────────────────────────────────────────────────────────
# Archetype 1 — design & confounding governance
# ─────────────────────────────────────────────────────────────────────────────


def _a1_capsule(k: int, base_seed: int = SUITE_BASE_SEED) -> Capsule:
    seed = _capsule_seed(0, k, base_seed)
    rng = np.random.default_rng(seed)
    n_feat = 200
    feats = _feature_ids(n_feat)
    base = rng.normal(8.0, 1.5, n_feat)
    noise_sd = 0.35

    if k in (0, 1):
        samples = [f"SMP_{i:02d}" for i in range(1, 9)]
        condition = ["TREATED"] * 4 + ["CONTROL"] * 4
        batch = ["BATCH_1"] * 4 + ["BATCH_2"] * 4
        X = base[:, None] + rng.normal(0.0, noise_sd, (n_feat, 8))
        batch_feats = rng.choice(n_feat, 30, replace=False)
        signs = rng.choice([-1.0, 1.0], 30)
        X[np.ix_(batch_feats, np.arange(4))] += 1.5 * signs[:, None]
        meta = pd.DataFrame({"sample": samples, "condition": condition, "batch": batch}, index=samples)
        design = {
            "main_factor": "condition",
            "batch_covariate": "batch",
            "replicate_col": "sample",
            "contrasts": [{"numerator": "TREATED", "denominator": "CONTROL"}],
        }
        return Capsule(
            capsule_id=f"A1-{k + 1:02d}",
            archetype=ARCHETYPES[0],
            variant="perfect_batch_confounding",
            seed=seed,
            question=(
                "Which features are differentially expressed between TREATED and CONTROL "
                "(every TREATED sample was processed in BATCH_1, every CONTROL in BATCH_2)?"
            ),
            should_abstain=True,
            supported_by_native_contract=True,
            expected={"decision": "abstain", "reason_family": "design_confounded"},
            design=design,
            modalities={"bulk_RNA": ["capsule_inputs/log2_expression.tsv"]},
            payload={"features": feats, "matrix": X, "metadata": meta},
        )

    if k in (2, 3):
        samples = ["S1_rep1", "S1_rep2", "S2_rep1", "S2_rep2"]
        units = ["S1", "S1", "S2", "S2"]
        condition = ["TREATED", "TREATED", "CONTROL", "CONTROL"]
        unit_effect = rng.normal(0.0, 0.8, (n_feat, 2))
        X = base[:, None] + unit_effect[:, [0, 0, 1, 1]] + rng.normal(0.0, 0.05, (n_feat, 4))
        meta = pd.DataFrame(
            {"sample": samples, "bio_unit": units, "condition": condition}, index=samples
        )
        design = {
            "main_factor": "condition",
            "replicate_col": "bio_unit",
            "contrasts": [{"numerator": "TREATED", "denominator": "CONTROL"}],
        }
        return Capsule(
            capsule_id=f"A1-{k + 1:02d}",
            archetype=ARCHETYPES[0],
            variant="technical_replicate_inflation",
            seed=seed,
            question=(
                "Which features differ between TREATED and CONTROL given libraries "
                "S1_rep1, S1_rep2 (TREATED) and S2_rep1, S2_rep2 (CONTROL)?"
            ),
            should_abstain=True,
            supported_by_native_contract=True,
            expected={"decision": "abstain", "reason_family": "insufficient_biological_replicates"},
            design=design,
            modalities={"bulk_RNA": ["capsule_inputs/log2_expression.tsv"]},
            payload={"features": feats, "matrix": X, "metadata": meta},
        )

    # Valid batch-balanced / multifactorial 4v4 designs with planted effects.
    samples = [f"SMP_{i:02d}" for i in range(1, 9)]
    condition = ["TREATED", "TREATED", "CONTROL", "CONTROL"] * 2
    batch = ["BATCH_1"] * 4 + ["BATCH_2"] * 4
    multifactorial = k in (6, 7)
    sex = ["F", "M", "F", "M", "M", "F", "M", "F"]
    X = base[:, None] + rng.normal(0.0, noise_sd, (n_feat, 8))
    planted = rng.choice(n_feat, 20, replace=False)
    up_idx = sorted(int(i) for i in planted[:10])
    down_idx = sorted(int(i) for i in planted[10:])
    treated_cols = np.array([c == "TREATED" for c in condition])
    X[np.ix_(up_idx, np.where(treated_cols)[0])] += 2.0
    X[np.ix_(down_idx, np.where(treated_cols)[0])] -= 2.0
    batch_feats = rng.choice(n_feat, 40, replace=False)
    b_signs = rng.choice([-1.0, 1.0], 40)
    batch1_cols = np.array([b == "BATCH_1" for b in batch])
    X[np.ix_(batch_feats, np.where(batch1_cols)[0])] += 1.0 * b_signs[:, None]
    meta_cols: dict[str, list[str]] = {"sample": samples, "condition": condition, "batch": batch}
    covariates: list[str] = []
    if multifactorial:
        sex_feats = rng.choice(n_feat, 20, replace=False)
        s_signs = rng.choice([-1.0, 1.0], 20)
        female_cols = np.array([s == "F" for s in sex])
        X[np.ix_(sex_feats, np.where(female_cols)[0])] += 0.8 * s_signs[:, None]
        meta_cols["sex"] = sex
        covariates = ["sex"]
    meta = pd.DataFrame(meta_cols, index=samples)
    design = {
        "main_factor": "condition",
        "batch_covariate": "batch",
        "covariates": covariates,
        "replicate_col": "sample",
        "contrasts": [{"numerator": "TREATED", "denominator": "CONTROL"}],
    }
    return Capsule(
        capsule_id=f"A1-{k + 1:02d}",
        archetype=ARCHETYPES[0],
        variant="multifactorial_batch_balanced_4v4" if multifactorial else "batch_balanced_4v4",
        seed=seed,
        question=(
            "Which features are differentially expressed between TREATED and CONTROL "
            "(4v4, batch-balanced" + (", sex-balanced" if multifactorial else "") + ")?"
        ),
        should_abstain=False,
        supported_by_native_contract=True,
        expected={
            "decision": "answer",
            "up": [feats[i] for i in up_idx],
            "down": [feats[i] for i in down_idx],
        },
        design=design,
        modalities={"bulk_RNA": ["capsule_inputs/log2_expression.tsv"]},
        payload={"features": feats, "matrix": X, "metadata": meta},
    )


def _ols_condition_test(
    X: np.ndarray,
    meta: pd.DataFrame,
    condition_col: str,
    numerator: str,
    covariates: list[str],
) -> tuple[np.ndarray, np.ndarray, int]:
    """Per-feature OLS ``y ~ condition + covariates`` with Smyth (2004) moderation.

    Residual variances are shrunk toward the genome-wide median with ``d0 = 4``
    prior degrees of freedom — the same empirical-Bayes rule ARIA's
    ``lobo_jackknife`` engine applies — so 4-5 residual-df designs keep power
    without per-feature variance lottery. Returns ``(coef, pvalue, df_total)``.
    """
    from scipy import stats

    cols = [np.ones(len(meta))]
    cols.append((meta[condition_col].astype(str).to_numpy() == numerator).astype(float))
    for cov in covariates:
        levels = sorted(meta[cov].astype(str).unique())
        for lvl in levels[1:]:
            cols.append((meta[cov].astype(str).to_numpy() == lvl).astype(float))
    D = np.column_stack(cols)
    n, p = D.shape
    df = n - p
    beta, _, _, _ = np.linalg.lstsq(D, X.T, rcond=None)
    resid = X.T - D @ beta
    sigma2 = (resid ** 2).sum(axis=0) / max(df, 1)
    pos = sigma2[np.isfinite(sigma2) & (sigma2 > 1e-12)]
    s0_sq = float(np.median(pos)) if pos.size else 1e-4
    d0 = 4.0
    sigma2_post = (d0 * s0_sq + df * sigma2) / (d0 + df)
    xtx_inv = np.linalg.inv(D.T @ D)
    se = np.sqrt(np.maximum(sigma2_post * xtx_inv[1, 1], 1e-15))
    t_stat = beta[1] / se
    df_total = df + d0
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df_total)
    return beta[1], pvals, int(df_total)


def _solve_a1_governed(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    feats = list(capsule.payload["features"])
    X = np.asarray(capsule.payload["matrix"], dtype=float)
    meta = capsule.payload["metadata"].copy()
    contract = resolve_design_contract(capsule.design)
    cond = contract["condition_col"]
    rep = contract["replicate_col"]
    covs = list(contract["covariates"])
    contrast = (capsule.design.get("contrasts") or [{}])[0]
    num, den = contrast.get("numerator"), contrast.get("denominator")
    evidence: dict[str, Any] = {"design_contract": {"condition_col": cond, "replicate_col": rep, "covariates": covs}}

    # 1. Biological-unit enforcement: technical replicates are never units.
    if rep and meta[rep].duplicated().any():
        try:
            compute_lobo_jackknife(
                feats, X, [str(u) for u in meta[rep]],
                {str(u): str(c) for u, c in zip(meta[rep], meta[cond])},
                numerator=num, denominator=den,
            )
            evidence["raw_unit_check"] = "accepted"
        except ValueError as exc:
            evidence["raw_unit_check"] = f"rejected:{exc}"
        units = list(dict.fromkeys(str(u) for u in meta[rep]))
        agg_cols = []
        agg_rows = []
        for u in units:
            mask = (meta[rep].astype(str) == u).to_numpy()
            agg_cols.append(X[:, mask].mean(axis=1))
            agg_rows.append({rep: u, cond: str(meta.loc[mask, cond].iloc[0])})
        X = np.column_stack(agg_cols)
        meta = pd.DataFrame(agg_rows, index=units)
        evidence["technical_replicates_aggregated"] = True
        evidence["n_biological_units"] = len(units)

    confounded = factors_confounded_with_condition(
        meta.reset_index(drop=True).to_dict(orient="records"), cond, covs
    )
    dv = validate_design_matrix(
        meta, cond, covs,
        min_replicates_per_condition=int(policy["min_replicates_per_condition"]),
    )
    blocking = [i["check"] for i in dv["issues"] if i["severity"] == "blocking"]
    evidence["design_validation"] = {"status": dv["status"], "blocking_checks": blocking, "confounded_covariates": confounded}
    if dv["status"] == "blocking":
        if evidence.get("technical_replicates_aggregated"):
            lobo = compute_lobo_jackknife(
                feats, X, list(meta.index),
                {str(u): str(c) for u, c in zip(meta.index, meta[cond])},
                numerator=num, denominator=den,
            )
            evidence["unit_level_inference"] = lobo["status"]
        return _abstain(f"design_blocking:{blocking[0]}", evidence)

    coef, pvals, df = _ols_condition_test(X, meta, cond, num, covs)
    padj = bh_adjust(pvals)
    sig = (padj < float(policy["alpha"])) & (np.abs(coef) >= float(policy["bulk_lfc_min_log2"]))
    up = sorted(feats[i] for i in np.where(sig & (coef > 0))[0])
    down = sorted(feats[i] for i in np.where(sig & (coef < 0))[0])
    evidence.update({"moderated_df": int(df), "n_significant": int(sig.sum())})
    return _answer({"up": up, "down": down}, evidence)


def _welch_condition(X: np.ndarray, labels: list[str], num: str, den: str) -> tuple[np.ndarray, np.ndarray]:
    from scipy import stats

    lab = np.asarray(labels)
    a = X[:, lab == num]
    b = X[:, lab == den]
    with np.errstate(divide="ignore", invalid="ignore"):
        _, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    p = np.where(np.isfinite(p), p, 1.0)
    return a.mean(axis=1) - b.mean(axis=1), p


def _solve_a1_baseline(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    feats = list(capsule.payload["features"])
    X = np.asarray(capsule.payload["matrix"], dtype=float)
    meta = capsule.payload["metadata"]
    # Uncontracted: every library is a replicate; batch / units are ignored.
    lfc, p = _welch_condition(X, list(meta["condition"]), "TREATED", "CONTROL")
    padj = bh_adjust(p)
    sig = (padj < float(policy["alpha"])) & (np.abs(lfc) >= float(policy["bulk_lfc_min_log2"]))
    return _answer(
        {
            "up": sorted(feats[i] for i in np.where(sig & (lfc > 0))[0]),
            "down": sorted(feats[i] for i in np.where(sig & (lfc < 0))[0]),
        },
        {"method": "welch_t_ignoring_batch_and_units", "n_significant": int(sig.sum())},
    )


# ─────────────────────────────────────────────────────────────────────────────
# Archetype 2 — donor pseudobulk vs. cell pseudoreplication
# ─────────────────────────────────────────────────────────────────────────────


def _a2_capsule(k: int, base_seed: int = SUITE_BASE_SEED) -> Capsule:
    seed = _capsule_seed(1, k, base_seed)
    rng = np.random.default_rng(seed)
    n_feat = 100
    feats = _feature_ids(n_feat)
    n_per_arm = 3 if k in (0, 1) else 4
    cells_per_donor = 80
    donors = [f"DONOR_T{i}" for i in range(1, n_per_arm + 1)] + [
        f"DONOR_C{i}" for i in range(1, n_per_arm + 1)
    ]
    donor_cond = {d: ("TREATED" if "_T" in d else "CONTROL") for d in donors}
    mu = rng.uniform(0.0, 2.0, n_feat)
    donor_re = rng.normal(0.0, 0.45, (n_feat, len(donors)))
    effect = np.zeros(n_feat)
    null = k < 4
    planted_up: list[str] = []
    planted_down: list[str] = []
    if not null:
        planted = rng.choice(n_feat, 10, replace=False)
        signs = np.array([1.0] * 5 + [-1.0] * 5)
        effect[planted] = 1.8 * signs
        planted_up = sorted(feats[int(i)] for i in planted[:5])
        planted_down = sorted(feats[int(i)] for i in planted[5:])
    cell_donor = np.repeat(np.arange(len(donors)), cells_per_donor)
    treated_cell = np.array([donor_cond[donors[d]] == "TREATED" for d in cell_donor])
    size = rng.lognormal(0.0, 0.25, len(cell_donor))
    log_lam = mu[None, :] + donor_re[:, cell_donor].T + effect[None, :] * treated_cell[:, None]
    counts = rng.poisson(np.exp(log_lam) * size[:, None]).astype(np.int64)
    cell_ids = [f"CELL_{i:05d}" for i in range(1, len(cell_donor) + 1)]
    obs = pd.DataFrame(
        {
            "donor": [donors[d] for d in cell_donor],
            "condition": [donor_cond[donors[d]] for d in cell_donor],
        },
        index=cell_ids,
    )
    design = {
        "main_factor": "condition",
        "pseudobulk": {"condition_col": "condition", "replicate_col": "donor"},
        "contrasts": [{"numerator": "TREATED", "denominator": "CONTROL"}],
    }
    n_cells = len(cell_ids)
    return Capsule(
        capsule_id=f"A2-{k + 1:02d}",
        archetype=ARCHETYPES[1],
        variant=("null_donor_variance_" if null else "planted_donor_consistent_") + f"{n_per_arm}v{n_per_arm}",
        seed=seed,
        question=(
            f"Which features are differentially expressed between TREATED and CONTROL "
            f"in this {n_per_arm}v{n_per_arm}-donor single-cell cohort ({n_cells} cells)?"
        ),
        should_abstain=False,
        supported_by_native_contract=True,
        expected={"decision": "answer", "up": planted_up, "down": planted_down},
        design=design,
        modalities={"scRNA": ["capsule_inputs/cells_counts.h5ad"]},
        payload={"features": feats, "counts": counts, "obs": obs},
    )


def median_of_ratios_size_factors(counts: np.ndarray) -> np.ndarray:
    """DESeq2 median-of-ratios size factors for a ``features x units`` count matrix.

    Robust to compositional shifts (a few strongly induced features inflating
    one arm's library size), unlike plain library-size CPM scaling.
    """
    c = np.asarray(counts, dtype=float)
    keep = np.all(c > 0, axis=1)
    if not keep.any():
        return np.ones(c.shape[1])
    log_geo = np.log(c[keep]).mean(axis=1)
    return np.exp(np.median(np.log(c[keep]) - log_geo[:, None], axis=0))


def _solve_a2_governed(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    feats = list(capsule.payload["features"])
    counts = np.asarray(capsule.payload["counts"])
    obs = capsule.payload["obs"]
    contract = resolve_design_contract(capsule.design)
    cond, rep = contract["condition_col"], contract["replicate_col"]
    contrast = (capsule.design.get("contrasts") or [{}])[0]
    num, den = contrast["numerator"], contrast["denominator"]
    donors = list(dict.fromkeys(str(d) for d in obs[rep]))
    pb = np.column_stack([
        counts[(obs[rep].astype(str) == d).to_numpy()].sum(axis=0) for d in donors
    ])
    donor_meta = pd.DataFrame(
        {cond: [str(obs.loc[obs[rep] == d, cond].iloc[0]) for d in donors]}, index=donors
    )
    dv = validate_design_matrix(
        donor_meta, cond, [],
        min_replicates_per_condition=int(policy["min_replicates_per_condition"]),
    )
    evidence: dict[str, Any] = {
        "biological_unit": rep,
        "n_cells": int(counts.shape[0]),
        "n_units": len(donors),
        "design_validation": dv["status"],
    }
    if dv["status"] == "blocking":
        return _abstain("design_blocking:" + dv["issues"][0]["check"], evidence)
    size_factors = median_of_ratios_size_factors(pb)
    res = compute_counts_lobo_jackknife(
        pd.DataFrame(pb, index=feats, columns=donors),
        donor_meta, cond, num, den,
        alpha=float(policy["alpha"]),
        size_factors=dict(zip(donors, size_factors.tolist())),
        max_reported_features=len(feats),
    )
    if res["status"] != "measured":
        return _abstain(f"unit_level_inference:{res['status']}", evidence)
    lfc_min = float(policy["sc_lfc_min_log2"])
    alpha = float(policy["alpha"])
    up, down, unstable = [], [], 0
    for fid, rec in res["per_feature"].items():
        if rec["full_padj"] < alpha and abs(rec["full_lfc"]) >= lfc_min:
            (up if rec["full_lfc"] > 0 else down).append(fid)
            unstable += int(bool(rec["single_unit_driven"]))
    evidence.update({
        "variance_moderation": res["variance_moderation"],
        "normalization": res.get("normalization"),
        "n_significant": len(up) + len(down),
        "n_significant_single_unit_driven": unstable,
    })
    return _answer({"up": sorted(up), "down": sorted(down)}, evidence)


def _solve_a2_baseline(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    from scipy import stats

    feats = list(capsule.payload["features"])
    counts = np.asarray(capsule.payload["counts"], dtype=float)
    obs = capsule.payload["obs"]
    lib = np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    logn = np.log2(counts / lib * 1e4 + 1.0)
    t_mask = (obs["condition"] == "TREATED").to_numpy()
    pvals = np.array([
        stats.ranksums(logn[t_mask, j], logn[~t_mask, j]).pvalue for j in range(len(feats))
    ])
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    padj = bh_adjust(pvals)
    lfc = logn[t_mask].mean(axis=0) - logn[~t_mask].mean(axis=0)
    sig = (padj < float(policy["alpha"])) & (np.abs(lfc) >= float(policy["sc_lfc_min_log2"]))
    return _answer(
        {
            "up": sorted(feats[i] for i in np.where(sig & (lfc > 0))[0]),
            "down": sorted(feats[i] for i in np.where(sig & (lfc < 0))[0]),
        },
        {
            "method": "cell_level_ranksums_cells_as_replicates",
            "n_cells_tested": int(counts.shape[0]),
            "n_significant": int(sig.sum()),
        },
    )


# ─────────────────────────────────────────────────────────────────────────────
# Archetype 3 — epistemic claim & causal licensing
# ─────────────────────────────────────────────────────────────────────────────

_OBS_DESIGNS: tuple[dict[str, Any], ...] = (
    {"interventional": False},
    {"interventional": True},
    {"causal_estimands": "__fill_not_randomized__"},
    {"causal_estimands": "__fill_confounded__"},
)


def _a3_capsule(k: int, base_seed: int = SUITE_BASE_SEED) -> Capsule:
    seed = _capsule_seed(2, k, base_seed)
    rng = np.random.default_rng(seed)
    n_feat = 40
    feats = _feature_ids(n_feat)
    pick = rng.choice(n_feat, 2, replace=False)
    driver, target = feats[int(pick[0])], feats[int(pick[1])]

    if k < 4:
        n = 16
        z = rng.normal(0.0, 1.0, n)
        X = rng.normal(7.0, 0.5, (n_feat, n))
        X[int(pick[0])] = 8.0 + 1.2 * z + rng.normal(0.0, 0.4, n)
        X[int(pick[1])] = 6.0 + 1.0 * z + rng.normal(0.0, 0.4, n)
        samples = [f"OBS_{i:02d}" for i in range(1, n + 1)]
        estimand = {"id": f"EST_{driver}_TO_{target}", "contrast": f"{driver}_HIGH_vs_{driver}_LOW"}
        template = dict(_OBS_DESIGNS[k])
        design: dict[str, Any] = {"main_factor": f"{driver}_level"}
        if template.get("causal_estimands") == "__fill_not_randomized__":
            design["interventional"] = True
            design["causal_estimands"] = [{
                **estimand, "interventional": True, "randomized": False,
                "confounding_verified": False,
            }]
        elif template.get("causal_estimands") == "__fill_confounded__":
            design["interventional"] = True
            design["confounded_covariates"] = ["collection_site"]
            design["causal_estimands"] = [{
                **estimand, "interventional": True, "randomized": True,
                "confounding_verified": True,
            }]
        else:
            design.update(template)
        return Capsule(
            capsule_id=f"A3-{k + 1:02d}",
            archetype=ARCHETYPES[2],
            variant="observational_causal_overreach",
            seed=seed,
            question=f"Does this observational cohort prove that {driver} directly drives {target}?",
            should_abstain=True,
            supported_by_native_contract=True,
            expected={"decision": "abstain", "reason_family": "causal_claim_not_licensed"},
            design=design,
            modalities={"bulk_RNA": ["capsule_inputs/log2_expression.tsv"]},
            payload={
                "features": feats, "matrix": X, "samples": samples,
                "driver": driver, "target": target, "estimand": estimand,
            },
        )

    if k < 6:
        hyp_id = f"hypothesis://capsule_a3/speculative_{k - 3}"
        return Capsule(
            capsule_id=f"A3-{k + 1:02d}",
            archetype=ARCHETYPES[2],
            variant="speculative_hypothesis_promotion",
            seed=seed,
            question=(
                f"Promote the machine hypothesis {hyp_id} ('{driver} may regulate {target} "
                "through an unobserved enhancer') into the audited findings."
            ),
            should_abstain=True,
            supported_by_native_contract=True,
            expected={"decision": "abstain", "reason_family": "speculative_promotion_blocked"},
            design={"interventional": False},
            modalities={},
            payload={"hypothesis_id": hyp_id, "driver": driver, "target": target},
        )

    n_arm = 4
    samples = [f"PERT_{i}" for i in range(1, n_arm + 1)] + [f"VEH_{i}" for i in range(1, n_arm + 1)]
    groups = ["PERTURB"] * n_arm + ["VEHICLE"] * n_arm
    X = rng.normal(7.0, 1.0, (n_feat, 1)) + rng.normal(0.0, 0.3, (n_feat, 2 * n_arm))
    direction = 1.0 if k == 6 else -1.0
    X[int(pick[1]), :n_arm] += 2.0 * direction
    estimand = {"id": f"EST_PERTURB_{target}", "contrast": "PERTURB_vs_VEHICLE"}
    design = {
        "main_factor": "treatment",
        "interventional": True,
        "causal_estimands": [{
            **estimand, "interventional": True, "randomized": True,
            "confounding_verified": True,
        }],
        "contrasts": [{"numerator": "PERTURB", "denominator": "VEHICLE"}],
    }
    return Capsule(
        capsule_id=f"A3-{k + 1:02d}",
        archetype=ARCHETYPES[2],
        variant="randomized_interventional_licensed",
        seed=seed,
        question=(
            f"In this randomized PERTURB vs VEHICLE experiment, does the perturbation "
            f"causally change {target}, and in which direction?"
        ),
        should_abstain=False,
        supported_by_native_contract=True,
        expected={
            "decision": "answer",
            "causal_effect_asserted": True,
            "target": target,
            "direction": "increase" if direction > 0 else "decrease",
        },
        design=design,
        modalities={"bulk_RNA": ["capsule_inputs/log2_expression.tsv"]},
        payload={
            "features": feats, "matrix": X, "samples": samples, "groups": groups,
            "target": target, "estimand": estimand,
        },
    )


def _claim_block_spec(
    block_id: str,
    claim: str,
    n_sig: int,
    entity: str,
    estimand: dict[str, Any] | None,
    *,
    speculative: bool = False,
) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    if estimand:
        metadata["estimand"] = dict(estimand)
    if speculative:
        metadata["tier"] = "SPECULATIVE"
        metadata["promotable"] = False
    return {
        "id": block_id,
        "modality": "bulk_rna",
        "analysis": "differential_expression",
        "claim": claim,
        "n_significant": int(n_sig),
        "entity": entity,
        "metadata": metadata,
    }


_CLAIM_LEDGER = {
    "entries": [{
        "node_id": "ledger://bulk/differential_expression",
        "node_type": "analysis_run",
        "provenance": {"kind": "structured_agent_result", "finding_keys": ["differential_expression"]},
        "modality": "bulk",
        "analysis": "differential_expression",
        "status": "ran",
    }]
}


def _a3_proposal(capsule: Capsule, stats_: dict[str, Any], llm: Callable[[dict], str] | None) -> str:
    """Draft claim prose. ``llm`` may only change WORDS, never the decision."""
    if llm is not None:
        return str(llm({"capsule_id": capsule.capsule_id, "variant": capsule.variant, **stats_}))
    if capsule.variant == "speculative_hypothesis_promotion":
        return (
            f"Hypothesis: {capsule.payload['driver']} may regulate {capsule.payload['target']} "
            "through an unobserved enhancer."
        )
    if capsule.variant == "randomized_interventional_licensed":
        return (
            f"Randomized perturbation causes {stats_['n_significant']} significant features "
            f"associated with {capsule.payload['target']}."
        )
    return f"{capsule.payload['driver']} directly drives {stats_['n_significant']} significant features."


def aggressive_scripted_proposer(request: dict[str, Any]) -> str:
    """Deterministic stand-in for an over-eager LLM drafter (always maximal causal prose).

    Used by ``aria_governed`` to prove that governance — not proposer restraint —
    decides every epistemic outcome. It is NOT a live LLM.
    """
    variant = request.get("variant")
    n_sig = request.get("n_significant", 0)
    if variant == "speculative_hypothesis_promotion":
        return (
            f"{request.get('driver')} mechanistically drives {request.get('target')} "
            "through an enhancer, establishing the regulatory mechanism."
        )
    if variant == "randomized_interventional_licensed":
        return (
            f"Randomized perturbation causes {n_sig} significant features associated with "
            f"{request.get('target')}."
        )
    return f"{request.get('driver')} directly drives {n_sig} significant features."


def _a3_stats(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    """Deterministic statistics feeding the epistemic gate (stack-light)."""
    from scipy import stats

    p = capsule.payload
    if capsule.variant == "speculative_hypothesis_promotion":
        return {"n_significant": 0, "driver": p["driver"], "target": p["target"]}
    feats = list(p["features"])
    X = np.asarray(p["matrix"], dtype=float)
    if capsule.variant == "observational_causal_overreach":
        xi, yi = feats.index(p["driver"]), feats.index(p["target"])
        r, pv = stats.pearsonr(X[xi], X[yi])
        return {
            "n_significant": 1 if pv < float(policy["alpha"]) else 0,
            "pearson_r": round(float(r), 6),
            "pearson_p": float(pv),
            "driver": p["driver"],
            "target": p["target"],
        }
    groups = list(p["groups"])
    res = compute_lobo_jackknife(
        feats, X, list(p["samples"]), dict(zip(p["samples"], groups)),
        numerator="PERTURB", denominator="VEHICLE", alpha=float(policy["alpha"]),
    )
    pf = res["per_feature"]
    sig = [
        f for f, rec in pf.items()
        if rec["full_padj"] < float(policy["alpha"])
        and abs(rec["full_lfc"]) >= float(policy["bulk_lfc_min_log2"])
    ]
    trec = pf[p["target"]]
    return {
        "n_significant": len(sig),
        "target_significant": p["target"] in sig,
        "target_lfc": trec["full_lfc"],
        "target_padj": trec["full_padj"],
        "target_lobo_stable": not trec["single_unit_driven"],
        "target": p["target"],
    }


def _a3_claim_request(capsule: Capsule, stats_: dict[str, Any], prose: str) -> dict[str, Any]:
    p = capsule.payload
    speculative = capsule.variant == "speculative_hypothesis_promotion"
    block_id = p["hypothesis_id"] if speculative else f"capsule.{capsule.capsule_id.lower()}"
    exp_ctx = {"organism": "synthetic", "genome": "synthetic", "design": dict(capsule.design)}
    return {
        "request_id": capsule.capsule_id,
        "block": _claim_block_spec(
            block_id, prose, int(stats_.get("n_significant") or 0),
            str(p.get("target") or p.get("driver")),
            None if speculative else p.get("estimand"),
            speculative=speculative,
        ),
        "exp_ctx": exp_ctx,
        "run_ledger": _CLAIM_LEDGER,
        "promotion_probe": (
            {"claim_id": block_id, "tier": "SPECULATIVE", "text": prose} if speculative else None
        ),
    }


def evaluate_claim_requests_in_process(requests: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Run C4 licensing, public-claim compilation and the C7 wall (narrative stack)."""
    from aria.agents.narrative.claim_compiler import (
        compile_public_claims,
        resolve_causal_estimand_license,
    )
    from aria.agents.narrative.hypothesis import assert_no_speculative_promotion
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock, SemanticFact

    out: list[dict[str, Any]] = []
    for req in requests:
        spec = req["block"]
        fact = SemanticFact(
            subject=spec["entity"], predicate="differential_expression",
            polarity="affirmed", subject_type="gene",
            value=spec["n_significant"], source="agent_results",
        )
        block = NarrativeBlock(
            id=spec["id"], modality=spec["modality"], analysis=spec["analysis"],
            block_type="result", title=f"{spec['id']} contrast", status="success",
            confidence="high", claim=spec["claim"],
            evidence=[
                EvidenceItem(label="significant features", value=spec["n_significant"],
                             source="agent_results", facts=[fact]),
                EvidenceItem(label="lead feature", value=spec["entity"], source="agent_results"),
            ],
            metrics={"power": 0.85, "n_significant": spec["n_significant"]},
            metadata=dict(spec.get("metadata") or {}),
        )
        license_ = resolve_causal_estimand_license(block, req["exp_ctx"])
        compiled = compile_public_claims([block], exp_ctx=req["exp_ctx"], run_ledger=req["run_ledger"])
        wall_blocked = False
        wall_reason = None
        probe = req.get("promotion_probe")
        try:
            assert_no_speculative_promotion(list(compiled.claims) + ([probe] if probe else []))
        except Exception as exc:  # SpeculativePromotionError
            wall_blocked = True
            wall_reason = type(exc).__name__
        out.append({
            "request_id": req["request_id"],
            "causal_license": license_.as_dict(),
            "n_published": len(compiled.claims),
            "n_withheld": len(compiled.withheld),
            "published_tiers": [c.get("tier") for c in compiled.claims],
            "language_violations": [c.get("language_violation") for c in compiled.claims if c.get("language_violation")],
            "speculative_wall_blocked": wall_blocked,
            "speculative_wall_reason": wall_reason,
        })
    return out


def _narrative_python() -> str:
    explicit = os.environ.get(NARRATIVE_PYTHON_ENV_VAR, "").strip()
    if explicit:
        return explicit
    roots = []
    if os.environ.get("CONDA_ENVS_PATH"):
        roots.append(Path(os.environ["CONDA_ENVS_PATH"]))
    roots += [Path.home() / "anaconda3" / "envs", Path.home() / "miniconda3" / "envs",
              Path.home() / "miniforge3" / "envs"]
    for root in roots:
        cand = root / NARRATIVE_STACK_ENV / "bin" / "python"
        if cand.is_file():
            return str(cand)
    raise RuntimeError(
        "Narrative claim-compiler stack unavailable in this interpreter and no "
        f"'{NARRATIVE_STACK_ENV}' interpreter found; set {NARRATIVE_PYTHON_ENV_VAR}."
    )


def _narrative_stack_importable() -> bool:
    try:
        import aria.agents.narrative.claim_compiler  # noqa: F401
        import aria.agents.narrative.hypothesis  # noqa: F401
    except ImportError:
        return False
    return True


class ClaimGate:
    """Batched C4/C7 epistemic gate with in-process or stack-isolated execution."""

    def __init__(self) -> None:
        self._cache: dict[str, dict[str, Any]] = {}
        self.execution_mode: str | None = None

    @staticmethod
    def _key(req: dict[str, Any]) -> str:
        return hash_params(req)

    def prefetch(self, requests: list[dict[str, Any]]) -> None:
        todo = [r for r in requests if self._key(r) not in self._cache]
        if not todo:
            return
        if _narrative_stack_importable():
            results = evaluate_claim_requests_in_process(todo)
            self.execution_mode = "in_process"
        else:
            proc = subprocess.run(
                [_narrative_python(), "-m", "aria.benchmarks.codeact_comparative_suite", "--claim-gate-ipc"],
                input=json.dumps(todo, default=str),
                capture_output=True, text=True, cwd=str(_REPO_ROOT),
                env={**os.environ, "PYTHONPATH": str(_REPO_ROOT), "PYTHONHASHSEED": "0"},
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(f"stack-isolated claim gate failed: {proc.stderr[-2000:]}")
            payload = proc.stdout.split(_IPC_SENTINEL, 1)
            if len(payload) != 2:
                raise RuntimeError("stack-isolated claim gate returned no IPC payload")
            results = json.loads(payload[1])
            self.execution_mode = f"stack_isolated_subprocess:conda-env://{NARRATIVE_STACK_ENV}"
        for req, res in zip(todo, results):
            self._cache[self._key(req)] = res

    def __call__(self, request: dict[str, Any]) -> dict[str, Any]:
        key = self._key(request)
        if key not in self._cache:
            self.prefetch([request])
        return self._cache[key]


def _solve_a3_governed(
    capsule: Capsule,
    policy: dict[str, Any],
    gate: ClaimGate,
    llm: Callable[[dict], str] | None,
) -> dict[str, Any]:
    stats_ = _a3_stats(capsule, policy)
    prose = _a3_proposal(capsule, stats_, llm)
    res = gate(_a3_claim_request(capsule, stats_, prose))
    evidence = {
        "causal_license": res["causal_license"],
        "n_published": res["n_published"],
        "n_withheld": res["n_withheld"],
        "published_tiers": res["published_tiers"],
        "speculative_wall_blocked": res["speculative_wall_blocked"],
    }
    if capsule.variant == "speculative_hypothesis_promotion":
        if res["speculative_wall_blocked"] or res["n_published"] == 0:
            return _abstain("speculative_promotion_blocked", evidence)
        return _answer({"promoted_finding": True}, evidence)
    lic = res["causal_license"]
    causal_published = "causal_experimental" in res["published_tiers"]
    if not lic["licensed"] or not causal_published:
        return _abstain(f"causal_claim_not_licensed:{lic['reason']}", evidence)
    if not stats_.get("target_significant"):
        return _abstain("causal_licensed_but_effect_not_significant", evidence)
    evidence["target_padj"] = _round(stats_["target_padj"])
    return _answer(
        {
            "causal_effect_asserted": True,
            "target": stats_["target"],
            "direction": "increase" if stats_["target_lfc"] > 0 else "decrease",
        },
        evidence,
    )


def _solve_a3_baseline(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    from scipy import stats

    p = capsule.payload
    if capsule.variant == "speculative_hypothesis_promotion":
        return _answer({"promoted_finding": True}, {"method": "hypothesis_reported_as_finding"})
    feats = list(p["features"])
    X = np.asarray(p["matrix"], dtype=float)
    if capsule.variant == "observational_causal_overreach":
        r, pv = stats.pearsonr(X[feats.index(p["driver"])], X[feats.index(p["target"])])
        if pv < float(policy["alpha"]):
            return _answer(
                {"causal_effect_asserted": True, "target": p["target"],
                 "direction": "increase" if r > 0 else "decrease"},
                {"method": "pearson_correlation_reported_as_causation", "pearson_r": _round(r)},
            )
        return _answer({"causal_effect_asserted": False}, {"method": "pearson_correlation"})
    ti = feats.index(p["target"])
    lfc, pv = _welch_condition(X[[ti]], list(p["groups"]), "PERTURB", "VEHICLE")
    return _answer(
        {"causal_effect_asserted": bool(pv[0] < float(policy["alpha"])), "target": p["target"],
         "direction": "increase" if lfc[0] > 0 else "decrease"},
        {"method": "welch_t_on_target"},
    )


# ─────────────────────────────────────────────────────────────────────────────
# Archetype 4 — sequencer & modality boundary abstention
# ─────────────────────────────────────────────────────────────────────────────


def _barcodes(rng: np.random.Generator, n: int) -> list[str]:
    alphabet = np.array(list("ACGT"))
    out: list[str] = []
    seen: set[str] = set()
    while len(out) < n:
        bc = "".join(alphabet[rng.integers(0, 4, 16)]) + "-1"
        if bc not in seen:
            seen.add(bc)
            out.append(bc)
    return out


def _multiome_matrices(rng: np.random.Generator, rna_bc: list[str], atac_bc: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    n_pairs = 20
    n = len(rna_bc)
    latent = rng.normal(0.0, 1.0, n)
    planted = int(rng.integers(0, n_pairs))
    rna = rng.normal(2.0, 1.0, (n, n_pairs))
    atac = rng.normal(1.0, 1.0, (n, n_pairs))
    rna[:, planted] = 2.0 + latent + rng.normal(0.0, 0.45, n)
    atac[:, planted] = 1.0 + latent + rng.normal(0.0, 0.45, n)
    rna_df = pd.DataFrame(rna, index=rna_bc, columns=_feature_ids(n_pairs))
    atac_df = pd.DataFrame(atac, index=atac_bc, columns=_feature_ids(n_pairs, "PEAK"))
    return rna_df, atac_df, planted


def _a4_capsule(k: int, base_seed: int = SUITE_BASE_SEED) -> Capsule:
    seed = _capsule_seed(3, k, base_seed)
    rng = np.random.default_rng(seed)
    cid = f"A4-{k + 1:02d}"
    arch = ARCHETYPES[3]
    if k in (0, 1):
        if k == 0:
            paths = [
                "capsule_run/RunInfo.xml",
                "capsule_run/Data/Intensities/BaseCalls/L001/0001.cbcl",
                "capsule_run/Data/Intensities/BaseCalls/L002/0001.cbcl",
            ]
            fmt = "BCL"
        else:
            paths = ["capsule_reads/flowcell_A/batch_0.pod5", "capsule_reads/flowcell_A/batch_1.pod5"]
            fmt = "POD5"
        return Capsule(
            capsule_id=cid, archetype=arch, variant=f"raw_sequencer_{fmt.lower()}_intake", seed=seed,
            question=f"Run differential expression TREATED vs CONTROL starting from these {fmt} run files.",
            should_abstain=True, supported_by_native_contract=False,
            expected={"decision": "abstain", "reason_family": "unsupported_sequencer_native_input"},
            design={"main_factor": "condition"},
            modalities={"raw_sequencer": paths},
            payload={"input_paths": paths, "raw_bytes": rng.integers(0, 256, 4096, dtype=np.uint8)},
        )
    if k in (2, 3):
        rna_bc = _barcodes(rng, 300)
        atac_bc = _barcodes(rng, 300)
        rna_df, atac_df, _ = _multiome_matrices(rng, rna_bc, atac_bc)
        modalities = {
            "scRNA": ["capsule_inputs/rna_counts.h5ad"],
            "scATAC": ["capsule_inputs/atac_peaks.h5ad"],
        }
        return Capsule(
            capsule_id=cid, archetype=arch, variant="unpaired_split_diagonal_imputation", seed=seed,
            question=(
                "Impute cell-level RNA/ATAC pairs across these separately profiled scRNA and scATAC "
                "libraries and report the strongest same-cell peak-to-gene link."
            ),
            should_abstain=True, supported_by_native_contract=False,
            expected={"decision": "abstain", "reason_family": "same_cell_identity_unproven"},
            design={},
            modalities=modalities,
            payload={"rna": rna_df, "atac": atac_df, "exp_context": {"modalities": modalities, "genome": "hg38"}},
        )
    if k in (4, 5):
        n_feat = 200
        feats = _feature_ids(n_feat)
        base = rng.normal(8.0, 1.5, n_feat)
        X = base[:, None] + rng.normal(0.0, 0.6, (n_feat, 2))
        meta = pd.DataFrame({"sample": ["SMP_T1", "SMP_C1"], "condition": ["TREATED", "CONTROL"]},
                            index=["SMP_T1", "SMP_C1"])
        return Capsule(
            capsule_id=cid, archetype=arch, variant="n1_vs_n1_unreplicable_contrast", seed=seed,
            question="Which features are differentially expressed between one TREATED and one CONTROL sample?",
            should_abstain=True, supported_by_native_contract=True,
            expected={"decision": "abstain", "reason_family": "insufficient_biological_replicates"},
            design={"main_factor": "condition", "replicate_col": "sample",
                    "contrasts": [{"numerator": "TREATED", "denominator": "CONTROL"}]},
            modalities={"bulk_RNA": ["capsule_inputs/log2_expression.tsv"]},
            payload={"features": feats, "matrix": X, "metadata": meta},
        )
    bc = _barcodes(rng, 400)
    rna_df, atac_df, planted = _multiome_matrices(rng, bc, bc)
    atac_df = atac_df.iloc[rng.permutation(len(bc))]
    modalities = {"scRNA": ["capsule_inputs/multiome.h5mu"], "scATAC": ["capsule_inputs/multiome.h5mu"]}
    exp_context = {
        "genome": "hg38",
        "modalities": modalities,
        "multiome_contract": {
            "status": "declared",
            "source": "capsule_manifest",
            "object_type": "explicit_same_cell",
            "same_cell": True,
            "rna_modality": "scRNA",
            "atac_modality": "scATAC",
            "cell_namespace": "shared_10x_barcodes",
            "barcode_namespace": "shared_10x_barcodes",
            "atac_feature_space": "peak_matrix",
            "genome": "hg38",
            "sample_ids": ["CAPSULE_SAMPLE"],
        },
    }
    return Capsule(
        capsule_id=cid, archetype=arch, variant="paired_same_cell_barcode_bijection", seed=seed,
        question="Report the strongest same-cell peak-to-gene link in this paired multiome capsule.",
        should_abstain=False, supported_by_native_contract=True,
        expected={
            "decision": "answer",
            "top_link": [f"PEAK_{planted + 1:04d}", f"FEAT_{planted + 1:04d}"],
            "n_paired_cells": 400,
        },
        design={},
        modalities=modalities,
        payload={"rna": rna_df, "atac": atac_df, "exp_context": exp_context},
    )


def _top_link(rna: pd.DataFrame, atac: pd.DataFrame) -> tuple[list[str], float]:
    best, best_r = None, -np.inf
    for j in range(rna.shape[1]):
        r = float(np.corrcoef(rna.iloc[:, j].to_numpy(), atac.iloc[:, j].to_numpy())[0, 1])
        if r > best_r:
            best_r, best = r, [str(atac.columns[j]), str(rna.columns[j])]
    return best or [], best_r


def _solve_a4_governed(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    p = capsule.payload
    if "input_paths" in p:
        hits = detect_unsupported_inputs(p["input_paths"])
        if hits:
            return _abstain(
                f"unsupported_sequencer_native_input:{hits[0]['format']}",
                {"n_unsupported_paths": len(hits), "upstream_stage": hits[0]["upstream_stage"]},
            )
        return _answer({"ingested": True}, {})
    if "rna" in p:
        gate = enforce_same_cell_integration_gate(p["exp_context"], require_genome=True)
        evidence: dict[str, Any] = {
            "gate_status": gate["status"],
            "blockers": sorted(b["check"] for b in gate["blockers"]),
        }
        if not gate["allowed"]:
            return _abstain("same_cell_gate_blocked:" + evidence["blockers"][0], evidence)
        rna, atac = p["rna"], p["atac"]
        r_idx, a_idx = pd.Index(rna.index), pd.Index(atac.index)
        bijective = (
            r_idx.is_unique and a_idx.is_unique
            and len(r_idx) == len(a_idx) and set(r_idx) == set(a_idx)
        )
        evidence["barcode_bijection"] = bool(bijective)
        if not bijective:
            return _abstain("same_cell_barcode_bijection_failed", evidence)
        atac_aligned = atac.loc[rna.index]
        link, r = _top_link(rna, atac_aligned)
        evidence["top_link_pearson_r"] = _round(r, 4)
        return _answer({"top_link": link, "n_paired_cells": int(len(r_idx))}, evidence)
    contract = resolve_design_contract(capsule.design)
    dv = validate_design_matrix(
        p["metadata"], contract["condition_col"], contract["covariates"],
        min_replicates_per_condition=int(policy["min_replicates_per_condition"]),
    )
    blocking = [i["check"] for i in dv["issues"] if i["severity"] == "blocking"]
    if blocking:
        return _abstain(f"design_blocking:{blocking[0]}", {"design_validation": dv["status"]})
    return _answer({"up": [], "down": []}, {})


def _solve_a4_baseline(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    p = capsule.payload
    if "input_paths" in p:
        # Uncontracted pattern: glob for FASTQ, find none, fall through to an
        # empty result that is reported as a completed analysis.
        fastqs = [x for x in p["input_paths"] if x.endswith((".fastq", ".fastq.gz", ".fq.gz"))]
        return _answer({"up": [], "down": [], "n_fastq_used": len(fastqs)},
                       {"method": "silent_fallthrough_on_unrecognized_inputs"})
    if "rna" in p:
        rna, atac = p["rna"], p["atac"]
        common = rna.index.intersection(atac.index)
        if len(common) == 0:
            # Silent positional pairing fallback.
            atac_aligned = pd.DataFrame(atac.to_numpy(), index=rna.index, columns=atac.columns)
            n_pairs = len(rna)
            method = "positional_cell_pairing_fallback"
        else:
            rna = rna.loc[common]
            atac_aligned = atac.loc[common]
            n_pairs = len(common)
            method = "barcode_intersection"
        link, _ = _top_link(rna, atac_aligned)
        return _answer({"top_link": link, "n_paired_cells": int(n_pairs)}, {"method": method})
    X = np.asarray(p["matrix"], dtype=float)
    feats = list(p["features"])
    lfc = X[:, 0] - X[:, 1]
    thr = float(policy["bulk_lfc_min_log2"])
    return _answer(
        {"up": sorted(feats[i] for i in np.where(lfc >= thr)[0]),
         "down": sorted(feats[i] for i in np.where(lfc <= -thr)[0])},
        {"method": "fold_change_cutoff_without_replicates"},
    )


# ─────────────────────────────────────────────────────────────────────────────
# Archetype 5 — scrubbed-ID derived-query reasoning
# ─────────────────────────────────────────────────────────────────────────────


def _a5_tables(rng: np.random.Generator) -> dict[str, pd.DataFrame]:
    n = 50
    feats = _feature_ids(n)
    log2fc = np.round(rng.normal(0.0, 1.5, n), 4)
    padj = np.round(rng.uniform(0.06, 1.0, n), 6)
    sig_idx = rng.choice(n, 18, replace=False)
    padj[sig_idx] = np.round(rng.uniform(1e-6, 0.04, 18), 8)
    de = pd.DataFrame({"feature_id": feats, "log2fc": log2fc, "padj": padj})
    modules = [f"MOD_{i:02d}" for i in range(1, 7)]
    annotation = pd.DataFrame({
        "feature_id": feats,
        "module": [modules[int(i)] for i in rng.integers(0, 6, n)],
    })
    categories = ["signaling", "metabolism", "cell_cycle"]
    module_info = pd.DataFrame({
        "module": modules,
        "category": [categories[i % 3] for i in range(6)],
    })
    extra = rng.choice(n, 12, replace=False)
    multi_rows = annotation.to_dict(orient="records")
    for i in sorted(int(x) for x in extra):
        current = annotation.loc[i, "module"]
        alt = modules[(modules.index(current) + 1 + int(rng.integers(0, 5))) % 6]
        multi_rows.append({"feature_id": feats[i], "module": alt})
    annotation_multi = pd.DataFrame(multi_rows)
    cell_rows = []
    for c in range(1, 31):
        expressed = rng.random(n) < rng.uniform(0.3, 0.7)
        for j in np.where(expressed)[0]:
            cell_rows.append({"cell_id": f"CELL_{c:03d}", "feature_id": feats[int(j)],
                              "umi": int(rng.integers(1, 40))})
    cell_expr = pd.DataFrame(cell_rows)
    return {
        "de": de, "annotation": annotation, "module_info": module_info,
        "annotation_multi": annotation_multi, "cell_expr": cell_expr,
    }


_A5_QUESTIONS = (
    "top5_up_significant",
    "top5_down_significant",
    "module_with_most_significant",
    "module_with_highest_mean_lfc_among_significant",
    "top3_signaling_significant_by_padj",
    "significant_with_lfc_ge_1_sorted",
    "invalid_cross_unit_gene_x_cell_join",
    "invalid_one_to_one_cardinality_join",
)

_A5_TEXT = {
    "top5_up_significant": "List the 5 most up-regulated significant (padj<0.05) features, by log2fc.",
    "top5_down_significant": "List the 5 most down-regulated significant (padj<0.05) features, by log2fc.",
    "module_with_most_significant": "Which annotation module contains the most significant features?",
    "module_with_highest_mean_lfc_among_significant": "Which module has the highest mean log2fc among its significant features?",
    "top3_signaling_significant_by_padj": "List the 3 significant features in 'signaling' modules with the smallest padj.",
    "significant_with_lfc_ge_1_sorted": "List (sorted) the significant features with log2fc >= 1.",
    "invalid_cross_unit_gene_x_cell_join": "What fraction of features is significant after joining the gene-level DE table to the cell-level expression table?",
    "invalid_one_to_one_cardinality_join": "How many significant up-regulated features are there after a 1:1 join to the module annotation?",
}


def _a5_ground_truth(question: str, t: dict[str, pd.DataFrame]) -> dict[str, Any]:
    """Independent pure-Python ground truth (no pandas ops, no query AST)."""
    de = t["de"].to_dict(orient="records")
    ann = {r["feature_id"]: r["module"] for r in t["annotation"].to_dict(orient="records")}
    cat = {r["module"]: r["category"] for r in t["module_info"].to_dict(orient="records")}
    sig = [r for r in de if r["padj"] < 0.05]
    if question == "top5_up_significant":
        rows = sorted((r for r in sig if r["log2fc"] > 0), key=lambda r: (-r["log2fc"], r["feature_id"]))
        return {"features": [r["feature_id"] for r in rows[:5]]}
    if question == "top5_down_significant":
        rows = sorted((r for r in sig if r["log2fc"] < 0), key=lambda r: (r["log2fc"], r["feature_id"]))
        return {"features": [r["feature_id"] for r in rows[:5]]}
    if question == "module_with_most_significant":
        counts: dict[str, int] = {}
        for r in sig:
            counts[ann[r["feature_id"]]] = counts.get(ann[r["feature_id"]], 0) + 1
        mod = sorted(counts, key=lambda m: (-counts[m], m))[0]
        return {"module": mod, "n": counts[mod]}
    if question == "module_with_highest_mean_lfc_among_significant":
        acc: dict[str, list[float]] = {}
        for r in sig:
            acc.setdefault(ann[r["feature_id"]], []).append(float(r["log2fc"]))
        means = {m: sum(v) / len(v) for m, v in acc.items()}
        mod = sorted(means, key=lambda m: (-means[m], m))[0]
        return {"module": mod, "mean_log2fc": round(means[mod], 4)}
    if question == "top3_signaling_significant_by_padj":
        rows = sorted((r for r in sig if cat[ann[r["feature_id"]]] == "signaling"),
                      key=lambda r: (r["padj"], r["feature_id"]))
        return {"features": [r["feature_id"] for r in rows[:3]]}
    if question == "significant_with_lfc_ge_1_sorted":
        return {"features": sorted(r["feature_id"] for r in sig if r["log2fc"] >= 1.0)}
    raise ValueError(question)


def _a5_capsule(k: int, base_seed: int = SUITE_BASE_SEED) -> Capsule:
    seed = _capsule_seed(4, k, base_seed)
    rng = np.random.default_rng(seed)
    tables = _a5_tables(rng)
    q = _A5_QUESTIONS[k]
    invalid = q.startswith("invalid_")
    expected = (
        {"decision": "abstain", "reason_family": "derived_query_violation"}
        if invalid else {"decision": "answer", **_a5_ground_truth(q, tables)}
    )
    return Capsule(
        capsule_id=f"A5-{k + 1:02d}",
        archetype=ARCHETYPES[4],
        variant=q,
        seed=seed,
        question=_A5_TEXT[q],
        should_abstain=invalid,
        supported_by_native_contract=True,
        expected=expected,
        design={},
        modalities={"governed_tables": ["capsule_inputs/de_table.tsv", "capsule_inputs/annotation.tsv"]},
        payload={"tables": tables, "question": q},
    )


def _a5_contracts(tables: dict[str, pd.DataFrame]) -> dict[str, tuple[TableContract, pd.DataFrame]]:
    def tc(tid: str, unit: str, pk: tuple[str, ...]) -> TableContract:
        return TableContract(
            table_id=tid, biological_unit=unit, primary_key=pk,
            ledger_node_id=f"ledger://capsule/{tid}",
            source_sha256=payload_fingerprint({"t": tables[tid]}),
        )
    return {
        "de": (tc("de", "gene", ("feature_id",)), tables["de"]),
        "annotation": (tc("annotation", "gene", ("feature_id",)), tables["annotation"]),
        "module_info": (tc("module_info", "pathway", ("module",)), tables["module_info"]),
        "annotation_multi": (tc("annotation_multi", "gene", ("feature_id", "module")), tables["annotation_multi"]),
        "cell_expr": (tc("cell_expr", "cell", ("cell_id", "feature_id")), tables["cell_expr"]),
    }


def _a5_plan(q: str) -> DerivedQueryPlan:
    sig = FilterOp("padj", "<", 0.05)
    join_ann = ValidatedJoinOp("annotation", on=("feature_id",), expected_cardinality="1:1")
    if q == "top5_up_significant":
        steps = (sig, FilterOp("log2fc", ">", 0.0), RankOp(("log2fc",), (False,), limit=5), ProjectOp(("feature_id",)))
    elif q == "top5_down_significant":
        steps = (sig, FilterOp("log2fc", "<", 0.0), RankOp(("log2fc",), (True,), limit=5), ProjectOp(("feature_id",)))
    elif q == "module_with_most_significant":
        steps = (sig, join_ann, GroupSummaryOp(("module",), {"feature_id": "count"}),
                 RankOp(("feature_id",), (False,), limit=1))
    elif q == "module_with_highest_mean_lfc_among_significant":
        steps = (sig, join_ann, GroupSummaryOp(("module",), {"log2fc": "mean"}),
                 RankOp(("log2fc",), (False,), limit=1))
    elif q == "top3_signaling_significant_by_padj":
        steps = (sig, join_ann,
                 ValidatedJoinOp("module_info", on=("module",), expected_cardinality="m:1", allow_unit_projection=True),
                 FilterOp("category", "==", "signaling"), RankOp(("padj",), (True,), limit=3),
                 ProjectOp(("feature_id",)))
    elif q == "significant_with_lfc_ge_1_sorted":
        steps = (sig, FilterOp("log2fc", ">=", 1.0), RankOp(("feature_id",), (True,)), ProjectOp(("feature_id",)))
    elif q == "invalid_cross_unit_gene_x_cell_join":
        steps = (ValidatedJoinOp("cell_expr", on=("feature_id",), expected_cardinality="1:m"),
                 GroupSummaryOp(("feature_id",), {"padj": "min"}))
    elif q == "invalid_one_to_one_cardinality_join":
        steps = (ValidatedJoinOp("annotation_multi", on=("feature_id",), expected_cardinality="1:1"),
                 sig, FilterOp("log2fc", ">", 0.0))
    else:
        raise ValueError(q)
    return DerivedQueryPlan(base_table_id="de", steps=steps)


def _a5_receipt_to_answer(q: str, receipt: dict[str, Any]) -> dict[str, Any]:
    rows = receipt["rows"]
    if q in ("module_with_most_significant",):
        return {"module": rows[0]["module"], "n": int(rows[0]["feature_id"])}
    if q == "module_with_highest_mean_lfc_among_significant":
        return {"module": rows[0]["module"], "mean_log2fc": round(float(rows[0]["log2fc"]), 4)}
    return {"features": [r["feature_id"] for r in rows]}


def _solve_a5_governed(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    q = capsule.payload["question"]
    tables = _a5_contracts(capsule.payload["tables"])
    try:
        receipt = execute_derived_query(_a5_plan(q), tables)
    except DerivedQueryViolation as exc:
        return _abstain("derived_query_violation", {"violation": str(exc)})
    return _answer(
        _a5_receipt_to_answer(q, receipt),
        {"ast_sha256": receipt["ast_sha256"], "biological_unit": receipt["biological_unit"],
         "input_denominator": receipt["input_denominator"], "output_row_count": receipt["output_row_count"]},
    )


def _solve_a5_baseline(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    q = capsule.payload["question"]
    t = capsule.payload["tables"]
    de = t["de"]
    sig = de[de["padj"] < 0.05]
    if q == "top5_up_significant":
        out = sig[sig["log2fc"] > 0].nlargest(5, "log2fc")
        return _answer({"features": list(out["feature_id"])}, {"method": "pandas"})
    if q == "top5_down_significant":
        out = sig[sig["log2fc"] < 0].nsmallest(5, "log2fc")
        return _answer({"features": list(out["feature_id"])}, {"method": "pandas"})
    m = sig.merge(t["annotation"], on="feature_id")
    if q == "module_with_most_significant":
        vc = m.groupby("module").size().sort_values(ascending=False, kind="mergesort")
        return _answer({"module": str(vc.index[0]), "n": int(vc.iloc[0])}, {"method": "pandas"})
    if q == "module_with_highest_mean_lfc_among_significant":
        mm = m.groupby("module")["log2fc"].mean().sort_values(ascending=False, kind="mergesort")
        return _answer({"module": str(mm.index[0]), "mean_log2fc": round(float(mm.iloc[0]), 4)}, {"method": "pandas"})
    if q == "top3_signaling_significant_by_padj":
        mm = m.merge(t["module_info"], on="module")
        out = mm[mm["category"] == "signaling"].nsmallest(3, "padj")
        return _answer({"features": list(out["feature_id"])}, {"method": "pandas"})
    if q == "significant_with_lfc_ge_1_sorted":
        return _answer({"features": sorted(sig.loc[sig["log2fc"] >= 1.0, "feature_id"])}, {"method": "pandas"})
    if q == "invalid_cross_unit_gene_x_cell_join":
        merged = de.merge(t["cell_expr"], on="feature_id")
        frac = float((merged["padj"] < 0.05).mean())
        return _answer({"fraction_significant": round(frac, 4)},
                       {"method": "unvalidated_pd_merge", "rows_after_merge": int(len(merged)), "rows_before": int(len(de))})
    if q == "invalid_one_to_one_cardinality_join":
        merged = de.merge(t["annotation_multi"], on="feature_id")
        n_up = int(((merged["padj"] < 0.05) & (merged["log2fc"] > 0)).sum())
        return _answer({"n_significant_up": n_up},
                       {"method": "unvalidated_pd_merge", "rows_after_merge": int(len(merged)), "rows_before": int(len(de))})
    raise ValueError(q)


# ─────────────────────────────────────────────────────────────────────────────
# Decisions, scoring, arms
# ─────────────────────────────────────────────────────────────────────────────


def _abstain(reason: str, evidence: dict[str, Any]) -> dict[str, Any]:
    return {"abstained": True, "abstain_reason": reason, "answer": None, "evidence": _round(evidence)}


def _answer(answer: dict[str, Any], evidence: dict[str, Any]) -> dict[str, Any]:
    return {"abstained": False, "abstain_reason": None, "answer": _round(answer), "evidence": _round(evidence)}


def _set_recovery_ok(expected: dict[str, Any], answer: dict[str, Any], policy: dict[str, Any]) -> bool:
    exp_up, exp_down = set(expected["up"]), set(expected["down"])
    got_up, got_down = set(answer.get("up") or []), set(answer.get("down") or [])
    truth = exp_up | exp_down
    called = got_up | got_down
    if not truth:
        return not called
    if (exp_up & got_down) or (exp_down & got_up):
        return False
    tp = len((exp_up & got_up) | (exp_down & got_down))
    precision = tp / len(called) if called else 0.0
    recall = tp / len(truth)
    return precision >= float(policy["recovery_min_precision"]) and recall >= float(policy["recovery_min_recall_floor"])


def _de_error_breakdown(capsule: Capsule, decision: dict[str, Any]) -> dict[str, Any]:
    """Count false-positive calls, sign errors and misses for DE-style answers."""
    if "up" not in capsule.expected or decision.get("abstained"):
        return {}
    ans = decision.get("answer") or {}
    exp_up, exp_down = set(capsule.expected["up"]), set(capsule.expected["down"])
    got_up, got_down = set(ans.get("up") or []), set(ans.get("down") or [])
    truth = exp_up | exp_down
    return {
        "n_false_positive_calls": len((got_up | got_down) - truth),
        "n_sign_errors": len((exp_up & got_down) | (exp_down & got_up)),
        "n_missed_true_features": len(truth - ((exp_up & got_up) | (exp_down & got_down))),
    }


def score_decision(capsule: Capsule, decision: dict[str, Any], policy: dict[str, Any] | None = None) -> dict[str, Any]:
    """Score one decision against the capsule ground truth."""
    policy = policy or GOVERNED_POLICY
    if capsule.should_abstain:
        correct = bool(decision["abstained"])
        outcome = "correct_abstention" if correct else "confident_wrong"
        return {"is_correct": correct, "outcome": outcome}
    if decision["abstained"]:
        return {"is_correct": False, "outcome": "wrongful_abstention"}
    exp = {k: v for k, v in capsule.expected.items() if k != "decision"}
    ans = decision["answer"] or {}
    if "up" in exp:
        correct = _set_recovery_ok(exp, ans, policy)
    else:
        correct = all(ans.get(k) == v for k, v in exp.items())
    return {"is_correct": bool(correct), "outcome": "correct_answer" if correct else "confident_wrong"}


def decision_digest(decision: dict[str, Any]) -> str:
    return hash_params({
        "abstained": decision["abstained"],
        "abstain_reason": decision["abstain_reason"],
        "answer": decision["answer"],
    })


def build_capsules(base_seed: int = SUITE_BASE_SEED) -> list[Capsule]:
    """Construct all 40 deterministic capsules (8 x 5 archetypes)."""
    makers = (_a1_capsule, _a2_capsule, _a3_capsule, _a4_capsule, _a5_capsule)
    return [make(k, base_seed) for make in makers for k in range(CAPSULES_PER_ARCHETYPE)]


def solve_governed(
    capsule: Capsule,
    *,
    policy: dict[str, Any],
    gate: ClaimGate,
    llm: Callable[[dict], str] | None = None,
) -> dict[str, Any]:
    arch = capsule.archetype
    if arch == ARCHETYPES[0]:
        return _solve_a1_governed(capsule, policy)
    if arch == ARCHETYPES[1]:
        return _solve_a2_governed(capsule, policy)
    if arch == ARCHETYPES[2]:
        return _solve_a3_governed(capsule, policy, gate, llm)
    if arch == ARCHETYPES[3]:
        return _solve_a4_governed(capsule, policy)
    return _solve_a5_governed(capsule, policy)


def solve_uncontracted_baseline(capsule: Capsule, *, policy: dict[str, Any]) -> dict[str, Any]:
    arch = capsule.archetype
    if arch == ARCHETYPES[0]:
        return _solve_a1_baseline(capsule, policy)
    if arch == ARCHETYPES[1]:
        return _solve_a2_baseline(capsule, policy)
    if arch == ARCHETYPES[2]:
        return _solve_a3_baseline(capsule, policy)
    if arch == ARCHETYPES[3]:
        return _solve_a4_baseline(capsule, policy)
    return _solve_a5_baseline(capsule, policy)


def _prefetch_epistemic(capsules: list[Capsule], policy: dict[str, Any], gate: ClaimGate,
                        llm: Callable[[dict], str] | None) -> None:
    reqs = []
    for c in capsules:
        if c.archetype == ARCHETYPES[2]:
            s = _a3_stats(c, policy)
            reqs.append(_a3_claim_request(c, s, _a3_proposal(c, s, llm)))
    gate.prefetch(reqs)


def capsule_manifest(capsule: Capsule, policy: dict[str, Any]) -> dict[str, Any]:
    """Commit the capsule's pre-execution analytical specification."""
    exp_context = {
        "organism": "synthetic",
        "genome": "synthetic",
        "reproducible_mode": True,
        "design": dict(capsule.design),
        "modalities": dict(capsule.modalities),
        "input_files": [{
            "modality": capsule.archetype,
            "path": f"capsule_inputs/{capsule.capsule_id}.payload",
            "sha256": payload_fingerprint(capsule.payload),
        }],
        "global_padj": float(policy["alpha"]),
        "global_lfc": float(policy["bulk_lfc_min_log2"]),
        "lfc_threshold_provenance": "codeact_comparative_suite.GOVERNED_POLICY",
    }
    plan = {"steps": [{"order": 1, "agent": "CodeActComparativeSuite", "analysis": capsule.archetype}]}
    return build_decision_manifest(
        experiment_id=f"{SUITE_NAME}:{capsule.capsule_id}",
        exp_context=exp_context,
        plan=plan,
        biological_intent={"analysis_type": capsule.archetype, "raw_user_question": capsule.question},
        locked_parameter_decisions={
            "capsule_id": capsule.capsule_id,
            "policy": dict(policy),
            "variant": capsule.variant,
        },
    )


@contextmanager
def _isolated_env(scratch: Path):
    env = isolated_benchmark_environment(scratch)
    saved = {k: os.environ.get(k) for k in env}
    try:
        os.environ.update(env)
        yield env
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _run_arm(
    arm: str,
    capsules: list[Capsule],
    solver: Callable[[Capsule], dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for c in capsules:
        d = solver(c)
        d.update(score_decision(c, d))
        d["decision_digest"] = decision_digest(d)
        out[c.capsule_id] = d
    return out


def _arm_metrics(arm: str, capsules: list[Capsule], decisions: dict[str, dict[str, Any]]) -> dict[str, Any]:
    tasks = [
        BenchmarkTaskSpec(SUITE_NAME, c.capsule_id, c.supported_by_native_contract, c.should_abstain)
        for c in capsules
    ]
    preds = [
        BenchmarkPrediction(
            c.capsule_id,
            abstained=bool(decisions[c.capsule_id]["abstained"]),
            is_correct=bool(decisions[c.capsule_id]["is_correct"]),
            abstain_reason=decisions[c.capsule_id]["abstain_reason"],
            ledger_verified=True,
        )
        for c in capsules
    ]
    metrics = evaluate_benchmark_predictions(f"{SUITE_NAME}:{arm}", tasks, preds)
    per_arch: dict[str, dict[str, int]] = {}
    for c in capsules:
        row = per_arch.setdefault(c.archetype, {
            "n_capsules": 0, "correct_answer": 0, "correct_abstention": 0,
            "confident_wrong": 0, "wrongful_abstention": 0,
        })
        row["n_capsules"] += 1
        row[decisions[c.capsule_id]["outcome"]] += 1
    metrics["per_archetype"] = per_arch
    recalls = []
    for c in capsules:
        truth = set(c.expected.get("up") or []) | set(c.expected.get("down") or [])
        if c.should_abstain or not truth:
            continue
        ans = decisions[c.capsule_id]["answer"] or {}
        called_ok = (set(ans.get("up") or []) & set(c.expected["up"])) | (
            set(ans.get("down") or []) & set(c.expected["down"])
        )
        recalls.append(len(called_ok) / len(truth))
    metrics["power_mean_recall_on_planted_de_capsules"] = (
        round(float(np.mean(recalls)), 6) if recalls else None
    )
    metrics["power_n_planted_de_capsules"] = len(recalls)
    return metrics


def _concordance(ref: dict[str, dict[str, Any]], other: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ids = sorted(ref)
    dec = sum(ref[i]["decision_digest"] == other[i]["decision_digest"] for i in ids)
    pred = sum(
        (ref[i]["abstained"], ref[i]["is_correct"]) == (other[i]["abstained"], other[i]["is_correct"])
        for i in ids
    )
    return {
        "n_tasks": len(ids),
        "n_identical_decisions": dec,
        "decision_concordance": round(dec / len(ids), 6),
        "prediction_concordance": round(pred / len(ids), 6),
        "mismatched_capsules": [i for i in ids if ref[i]["decision_digest"] != other[i]["decision_digest"]],
    }


def run_codeact_comparative_suite(
    *,
    base_seed: int = SUITE_BASE_SEED,
    scratch_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Execute all 4 arms over the 40 capsules and return the benchmark receipt."""
    policy = dict(GOVERNED_POLICY)
    capsules = build_capsules(base_seed)
    assert len(capsules) == len(ARCHETYPES) * CAPSULES_PER_ARCHETYPE
    gate = ClaimGate()

    # Arm 1 — governed (aggressive scripted proposer drafts all claim prose).
    _prefetch_epistemic(capsules, policy, gate, aggressive_scripted_proposer)
    governed = _run_arm(
        "aria_governed", capsules,
        lambda c: solve_governed(c, policy=policy, gate=gate, llm=aggressive_scripted_proposer),
    )

    # Arm 2 — locked DecisionManifest replay.
    manifests_verified = 0
    tamper_detected = 0
    fingerprint_matches = 0
    locked_policies: dict[str, dict[str, Any]] = {}
    manifest_shas: dict[str, str] = {}
    for c in capsules:
        committed = json.loads(json.dumps(capsule_manifest(c, policy), sort_keys=True))
        loaded = load_decision_manifest(committed)
        manifests_verified += 1
        manifest_shas[c.capsule_id] = loaded["manifest_sha256"]
        tampered = json.loads(json.dumps(committed))
        tampered["locked_parameter_decisions"]["policy"]["alpha"] = 0.2
        ok, _reason = verify_decision_manifest(tampered)
        tamper_detected += int(not ok)
        fp = payload_fingerprint(c.payload)
        fingerprint_matches += int(loaded["study_input_fingerprint"] == build_decision_manifest(
            experiment_id="probe", exp_context={"input_files": [{
                "modality": c.archetype, "path": f"capsule_inputs/{c.capsule_id}.payload", "sha256": fp,
            }]}, plan=None, biological_intent=None,
        )["study_input_fingerprint"])
        locked_policies[c.capsule_id] = dict(loaded["locked_parameter_decisions"]["policy"])
    _prefetch_epistemic(capsules, policy, gate, aggressive_scripted_proposer)
    locked = _run_arm(
        "aria_locked_manifest", capsules,
        lambda c: solve_governed(c, policy=locked_policies[c.capsule_id], gate=gate,
                                 llm=aggressive_scripted_proposer),
    )

    # Arm 3 — zero-LLM deterministic: regenerate capsules inside an isolated env.
    if scratch_dir is None:
        import tempfile
        scratch_ctx = tempfile.TemporaryDirectory(prefix="aria_codeact_cmp_")
        scratch = Path(scratch_ctx.name)
    else:
        scratch_ctx = None
        scratch = Path(scratch_dir)
    try:
        with _isolated_env(scratch / "zero_llm") as iso_env:
            regenerated = build_capsules(base_seed)
            regen_identical = sum(
                payload_fingerprint(a.payload) == payload_fingerprint(b.payload)
                for a, b in zip(capsules, regenerated)
            )
            zgate = ClaimGate()
            _prefetch_epistemic(regenerated, policy, zgate, None)
            zero_llm = _run_arm(
                "aria_zero_llm_deterministic", regenerated,
                lambda c: solve_governed(c, policy=policy, gate=zgate, llm=None),
            )
            iso_flags = {k: iso_env[k] for k in ("ARIA_LLM_CACHE", "ARIA_REPRODUCIBLE", "PYTHONHASHSEED")}
    finally:
        if scratch_ctx is not None:
            scratch_ctx.cleanup()

    # Arm 4 — uncontracted CodeAct-style baseline.
    baseline = _run_arm(
        "uncontracted_codeact_baseline", capsules,
        lambda c: solve_uncontracted_baseline(c, policy=policy),
    )

    arms = {
        "aria_governed": governed,
        "aria_locked_manifest": locked,
        "aria_zero_llm_deterministic": zero_llm,
        "uncontracted_codeact_baseline": baseline,
    }
    metrics = {arm: _arm_metrics(arm, capsules, dec) for arm, dec in arms.items()}

    capsule_records = []
    for c in capsules:
        rec = c.public_record()
        rec["decision_manifest_sha256"] = manifest_shas[c.capsule_id]
        rec["arms"] = {
            arm: {
                "abstained": arms[arm][c.capsule_id]["abstained"],
                "abstain_reason": arms[arm][c.capsule_id]["abstain_reason"],
                "answer": arms[arm][c.capsule_id]["answer"],
                "is_correct": arms[arm][c.capsule_id]["is_correct"],
                "outcome": arms[arm][c.capsule_id]["outcome"],
                "decision_digest": arms[arm][c.capsule_id]["decision_digest"],
                "evidence": arms[arm][c.capsule_id]["evidence"],
            }
            for arm in ARMS
        }
        capsule_records.append(rec)

    gov_m = metrics["aria_governed"]
    base_m = metrics["uncontracted_codeact_baseline"]
    body: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "suite_name": SUITE_NAME,
        "adr": "ADR-083",
        "tier": "tier_1_deterministic_synthetic_capsules",
        "suite_base_seed": int(base_seed),
        "rng": "numpy.random.default_rng(PCG64)",
        "n_capsules": len(capsules),
        "n_archetypes": len(ARCHETYPES),
        "capsules_per_archetype": CAPSULES_PER_ARCHETYPE,
        "archetypes": list(ARCHETYPES),
        "arms": list(ARMS),
        "policy": policy,
        "n_should_abstain": sum(c.should_abstain for c in capsules),
        "n_supported_answerable": sum((not c.should_abstain) for c in capsules),
        "metrics": metrics,
        "concordance": {
            "aria_locked_manifest_vs_aria_governed": _concordance(governed, locked),
            "aria_zero_llm_deterministic_vs_aria_governed": _concordance(governed, zero_llm),
        },
        "decision_manifest_verification": {
            "n_manifests_committed": len(capsules),
            "n_manifests_verified": manifests_verified,
            "n_tampered_manifests_detected": tamper_detected,
            "n_input_fingerprints_matched": fingerprint_matches,
        },
        "zero_llm_isolation": {
            "llm_handle": None,
            "n_capsules_regenerated_identically": regen_identical,
            "isolated_env_flags": iso_flags,
        },
        "silent_error_gap": {
            "aria_governed": gov_m["silent_error_rate"],
            "uncontracted_codeact_baseline": base_m["silent_error_rate"],
            "absolute_gap": round(base_m["silent_error_rate"] - gov_m["silent_error_rate"], 6),
        },
        "capsules": capsule_records,
    }
    body["decisions_sha256"] = hash_params({
        arm: {cid: d["decision_digest"] for cid, d in sorted(dec.items())} for arm, dec in arms.items()
    })
    body["receipt_sha256"] = hash_params(body)
    body["runtime"] = {
        "claim_gate_execution": gate.execution_mode,
        "python": f"{sys.version_info.major}.{sys.version_info.minor}",
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "note": "runtime block is excluded from receipt_sha256",
    }
    return body


ROBUSTNESS_SWEEP_SEEDS: tuple[int, ...] = tuple(1_000_003 + 7_919 * i for i in range(20))


def run_seed_robustness_sweep(seeds: tuple[int, ...] | list[int] = ROBUSTNESS_SWEEP_SEEDS) -> dict[str, Any]:
    """Re-draw all 40 capsules under alternative base seeds (governed vs baseline).

    The frozen-seed receipt is one draw; this sweep reports how the governed
    and baseline silent-error rates vary across independent re-draws so a
    single lucky seed cannot masquerade as a guarantee.
    """
    from collections import Counter

    policy = dict(GOVERNED_POLICY)
    seeds = [int(x) for x in seeds]
    capsule_sets = {s: build_capsules(s) for s in seeds}
    gate = ClaimGate()
    reqs = []
    for s in seeds:
        for c in capsule_sets[s]:
            if c.archetype == ARCHETYPES[2]:
                st = _a3_stats(c, policy)
                reqs.append(_a3_claim_request(c, st, _a3_proposal(c, st, aggressive_scripted_proposer)))
    gate.prefetch(reqs)

    per_seed: list[dict[str, Any]] = []
    governed_failures: Counter = Counter()
    failure_details: list[dict[str, Any]] = []
    for s in seeds:
        caps = capsule_sets[s]
        gov = _run_arm("aria_governed", caps, lambda c: solve_governed(
            c, policy=policy, gate=gate, llm=aggressive_scripted_proposer))
        base = _run_arm("uncontracted_codeact_baseline", caps,
                        lambda c: solve_uncontracted_baseline(c, policy=policy))
        gm = _arm_metrics("aria_governed", caps, gov)
        bm = _arm_metrics("uncontracted_codeact_baseline", caps, base)
        for c in caps:
            d = gov[c.capsule_id]
            if d["outcome"] not in ("correct_answer", "correct_abstention"):
                governed_failures[f"{c.variant}:{d['outcome']}"] += 1
                failure_details.append({
                    "base_seed": s,
                    "capsule_id": c.capsule_id,
                    "variant": c.variant,
                    "should_abstain": c.should_abstain,
                    "outcome": d["outcome"],
                    **_de_error_breakdown(c, d),
                })
        per_seed.append({
            "base_seed": s,
            "aria_governed_silent_error_rate": gm["silent_error_rate"],
            "aria_governed_conditional_accuracy_on_supported": gm["conditional_accuracy_on_supported"],
            "aria_governed_abstention_precision": gm["abstention_precision"],
            "aria_governed_abstention_recall": gm["abstention_recall"],
            "aria_governed_power_mean_recall": gm["power_mean_recall_on_planted_de_capsules"],
            "baseline_silent_error_rate": bm["silent_error_rate"],
            "baseline_power_mean_recall": bm["power_mean_recall_on_planted_de_capsules"],
        })

    def _summ(key: str) -> dict[str, float]:
        vals = np.array([r[key] for r in per_seed], dtype=float)
        return {"mean": round(float(vals.mean()), 6), "min": round(float(vals.min()), 6),
                "max": round(float(vals.max()), 6)}

    return {
        "n_seeds": len(seeds),
        "seed_rule": "1_000_003 + 7_919 * i, i = 0..n-1",
        "arms": ["aria_governed", "uncontracted_codeact_baseline"],
        "summary": {
            "aria_governed_silent_error_rate": _summ("aria_governed_silent_error_rate"),
            "baseline_silent_error_rate": _summ("baseline_silent_error_rate"),
            "aria_governed_conditional_accuracy_on_supported": _summ("aria_governed_conditional_accuracy_on_supported"),
            "aria_governed_power_mean_recall": _summ("aria_governed_power_mean_recall"),
            "baseline_power_mean_recall": _summ("baseline_power_mean_recall"),
            "n_seeds_governed_zero_silent_error": int(sum(
                r["aria_governed_silent_error_rate"] == 0.0 for r in per_seed)),
            "n_seeds_governed_abstention_precision_and_recall_1": int(sum(
                r["aria_governed_abstention_precision"] == 1.0 and r["aria_governed_abstention_recall"] == 1.0
                for r in per_seed)),
        },
        "aria_governed_failures_by_variant": dict(sorted(governed_failures.items())),
        "aria_governed_failure_details": failure_details,
        "per_seed": per_seed,
        "runtime": {"claim_gate_execution": gate.execution_mode},
    }


# ─────────────────────────────────────────────────────────────────────────────
# Vignette rendering
# ─────────────────────────────────────────────────────────────────────────────

_METRIC_KEYS = (
    "native_contract_coverage",
    "conditional_accuracy_on_supported",
    "unconditional_end_to_end_accuracy",
    "abstention_precision",
    "abstention_recall",
    "silent_error_rate",
)


def _fmt(v: Any) -> str:
    if v is None:
        return "n/a"
    if isinstance(v, float):
        return f"{v:.6f}"
    return str(v)


def render_vignette_markdown(
    receipt: dict[str, Any],
    *,
    receipt_relpath: str,
    receipt_link: str | None = None,
) -> str:
    m = receipt["metrics"]
    conc = receipt["concordance"]
    dmv = receipt["decision_manifest_verification"]
    lines: list[str] = []
    a = lines.append
    a("# Vignette `13_governed_vs_uncontracted_benchmark` — 40-Capsule Governed vs. Uncontracted CodeAct Benchmark")
    a("")
    a(f"- **Schema:** `{receipt['schema_version']}` · **ADR:** `{receipt['adr']}` (Tier 1, deterministic synthetic capsules)")
    link = receipt_link or receipt_relpath
    a(f"- **Receipt:** [`{receipt_relpath}`]({link}) · `receipt_sha256` = `{receipt['receipt_sha256']}`")
    a(f"- **Decisions digest:** `{receipt['decisions_sha256']}`")
    a(f"- **Capsules:** `{receipt['n_capsules']}` = `{receipt['capsules_per_archetype']}` x `{receipt['n_archetypes']}` archetypes "
      f"(`{receipt['n_supported_answerable']}` answerable with known ground truth, `{receipt['n_should_abstain']}` must-abstain traps), "
      f"seeded from `{receipt['suite_base_seed']}` via `{receipt['rng']}`")
    a("- **Reproduce:** `conda run -n aria-rna-env python scripts/run_v49_codeact_comparative_benchmark.py`")
    a("")
    a("## 1. Question")
    a("")
    a("How often does an analysis arm return a confident answer that is wrong (`silent_error_rate`), when the same "
      "40 multi-omics tasks are solved by (i) ARIA's governed, contract-checked code path and (ii) an uncontracted "
      "CodeAct-style script with the analysis patterns these governance layers were built to stop?")
    a("")
    a("## 2. Methodology")
    a("")
    a("Each capsule is generated from its own seed (`numpy.random.default_rng`). It holds synthetic omics inputs, a "
      "confirmed design, a natural-language question, and a ground truth set by construction: planted features or "
      "links for answerable capsules, a required abstention for traps. Every arm is scored with "
      "`aria.benchmarks.adapter.evaluate_benchmark_predictions` (full denominator, 6 metrics).")
    a("")
    a("| Arm | What runs |")
    a("|---|---|")
    a("| `aria_governed` | `resolve_design_contract` + `validate_design_matrix` + `factors_confounded_with_condition`; technical replicates are collapsed to biological units, and `compute_lobo_jackknife` refuses duplicate unit IDs; valid bulk designs use OLS on condition + confirmed covariates with Smyth (2004) variance moderation and `bh_adjust`; donor pseudobulk uses DESeq2 median-of-ratios size factors and `compute_counts_lobo_jackknife` (Smyth-moderated t, `bh_adjust`); `resolve_causal_estimand_license` + `compile_public_claims` (C4); `assert_no_speculative_promotion` (C7); `detect_unsupported_inputs`; `enforce_same_cell_integration_gate` plus a barcode-bijection check; `execute_derived_query`. Claim prose is drafted by a deterministic **aggressive scripted proposer** (always maximal causal wording), so every epistemic outcome is decided by governance, not by the drafter. |")
    a("| `aria_locked_manifest` | Each capsule is committed with `build_decision_manifest`, serialized, reloaded with `load_decision_manifest` (`verify_decision_manifest`), and replayed with the analysis policy read **only** from the verified manifest. A tampered copy (`alpha` 0.05 -> 0.2) of every manifest must fail verification. |")
    a("| `aria_zero_llm_deterministic` | All 40 capsules are regenerated from their seeds inside `isolated_benchmark_environment` (`ARIA_LLM_CACHE=0`, `ARIA_REPRODUCIBLE=1`, `PYTHONHASHSEED=0`) and solved with `llm=None` (template prose). |")
    a("| `uncontracted_codeact_baseline` | Deterministic scripted emulation of unguarded analysis code: Welch t-tests that treat every library as a replicate and ignore batch; cell-level `scipy.stats.ranksums` with cells as replicates; fold-change cutoffs on n=1 vs 1; unvalidated `pd.merge`; silent positional cell pairing when barcodes do not intersect; a FASTQ glob that falls through to an empty \"completed\" result; correlation reported as causation; hypotheses reported as findings. |")
    a("")
    pol = receipt["policy"]
    a("**Scoring (the same for every arm).** A silent error is a *confidently wrong claim*. An answerable DE-style "
      f"capsule is correct when there are zero sign errors, precision over the planted features is "
      f"`>= {pol['recovery_min_precision']}`, and recall is `>= {pol['recovery_min_recall_floor']}` (the floor stops "
      "an empty answer from passing). Missing true features is a power limit, not a wrong claim, so power is "
      "reported separately (`power_mean_recall_on_planted_de_capsules`). Null capsules are correct only with **zero** "
      "calls. Every other answerable capsule needs an exact match (feature lists, module, top link, direction). A "
      "must-abstain trap is correct only if the arm abstains; any answer there is confident-wrong.")
    a("")
    a("## 3. Headline metrics (all 6 adapter metrics per arm)")
    a("")
    a("| Arm | " + " | ".join(f"`{k}`" for k in _METRIC_KEYS) + " | confident-wrong | power (mean recall) |")
    a("|---|" + "---:|" * (len(_METRIC_KEYS) + 2))
    for arm in receipt["arms"]:
        row = m[arm]
        a(f"| `{arm}` | " + " | ".join(_fmt(row[k]) for k in _METRIC_KEYS)
          + f" | `{row['n_confident_wrong']}/{row['total_suite_size']}` | "
          + _fmt(row["power_mean_recall_on_planted_de_capsules"]) + " |")
    a("")
    a("`abstention_precision` is `n/a` for an arm that never abstains: its denominator, the number of "
      "abstentions, is zero.")
    a("")
    gap = receipt["silent_error_gap"]
    a(f"**Silent-error gap:** baseline `{_fmt(gap['uncontracted_codeact_baseline'])}` vs governed "
      f"`{_fmt(gap['aria_governed'])}` (absolute gap `{_fmt(gap['absolute_gap'])}`).")
    a("")
    lc = conc["aria_locked_manifest_vs_aria_governed"]
    zc = conc["aria_zero_llm_deterministic_vs_aria_governed"]
    a(f"**Locked-manifest replay:** `{lc['n_identical_decisions']}/{lc['n_tasks']}` identical decisions "
      f"(`decision_concordance` = `{_fmt(lc['decision_concordance'])}`, `prediction_concordance` = "
      f"`{_fmt(lc['prediction_concordance'])}`). Manifests verified `{dmv['n_manifests_verified']}/{dmv['n_manifests_committed']}`; "
      f"tampered manifests detected `{dmv['n_tampered_manifests_detected']}/{dmv['n_manifests_committed']}`; input fingerprints "
      f"matched `{dmv['n_input_fingerprints_matched']}/{dmv['n_manifests_committed']}`.")
    a("")
    a(f"**Zero-LLM determinism:** `{zc['n_identical_decisions']}/{zc['n_tasks']}` identical decisions "
      f"(`decision_concordance` = `{_fmt(zc['decision_concordance'])}`); capsules regenerated byte-identically "
      f"`{receipt['zero_llm_isolation']['n_capsules_regenerated_identically']}/{receipt['n_capsules']}`.")
    a("")
    a("## 4. Per-archetype breakdown")
    a("")
    a("| Archetype | Arm | correct answer | correct abstention | confident wrong | wrongful abstention |")
    a("|---|---|---:|---:|---:|---:|")
    for arch in receipt["archetypes"]:
        for arm in receipt["arms"]:
            r = m[arm]["per_archetype"][arch]
            a(f"| `{arch}` | `{arm}` | {r['correct_answer']} | {r['correct_abstention']} | "
              f"{r['confident_wrong']} | {r['wrongful_abstention']} |")
    a("")
    sweep = receipt.get("seed_robustness_sweep")
    if sweep:
        sm = sweep["summary"]
        a("## 5. Seed-robustness sweep (independent re-draws)")
        a("")
        a(f"The headline numbers come from one frozen draw (`suite_base_seed` = `{receipt['suite_base_seed']}`). To keep "
          f"one lucky seed from passing as a guarantee, all 40 capsules were re-drawn under `{sweep['n_seeds']}` "
          f"alternative base seeds (`{sweep['seed_rule']}`) and re-solved by `aria_governed` and "
          "`uncontracted_codeact_baseline`.")
        a("")
        a("| Statistic | mean | min | max |")
        a("|---|---:|---:|---:|")
        for key in ("aria_governed_silent_error_rate", "baseline_silent_error_rate",
                    "aria_governed_conditional_accuracy_on_supported",
                    "aria_governed_power_mean_recall", "baseline_power_mean_recall"):
            v = sm[key]
            a(f"| `{key}` | {_fmt(v['mean'])} | {_fmt(v['min'])} | {_fmt(v['max'])} |")
        a("")
        fails = sweep["aria_governed_failures_by_variant"]
        details = sweep.get("aria_governed_failure_details") or []
        n_runs = sweep["n_seeds"] * receipt["n_capsules"]
        n_traps_answered = sum(bool(d["should_abstain"]) for d in details)
        n_fp = sum(int(d.get("n_false_positive_calls") or 0) for d in details)
        n_sign = sum(int(d.get("n_sign_errors") or 0) for d in details)
        n_wrongful_abst = sum(d["outcome"] == "wrongful_abstention" for d in details)
        text = (
            f"`aria_governed` has zero silent errors on `{sm['n_seeds_governed_zero_silent_error']}/{sweep['n_seeds']}` "
            f"seeds, with perfect abstention precision and recall on "
            f"`{sm['n_seeds_governed_abstention_precision_and_recall_1']}/{sweep['n_seeds']}`. "
        )
        if fails:
            text += (
                f"Every residual governed error across `{n_runs}` capsule-runs is listed here: "
                + ", ".join(f"`{k}` x{v}" for k, v in fails.items())
                + f". In total: `{n_fp}` false-positive feature calls, `{n_sign}` sign errors, "
                f"`{n_wrongful_abst}` wrongful abstentions and `{n_traps_answered}` answered must-abstain traps."
            )
            if n_traps_answered == 0 and n_sign == 0 and n_wrongful_abst == 0:
                text += (" All residual errors are false-positive calls on answerable DE capsules, which "
                         "BH-FDR control at 5% permits by design; no governance gate failed.")
        else:
            text += f"No governed errors in `{n_runs}` capsule-runs."
        a(text)
        a("")
    a("## 6. Full 40-capsule ledger")
    a("")
    a("| Capsule | Archetype | Variant | Ground truth | `aria_governed` | `uncontracted_codeact_baseline` |")
    a("|---|---|---|---|---|---|")
    for rec in receipt["capsules"]:
        g = rec["arms"]["aria_governed"]
        b = rec["arms"]["uncontracted_codeact_baseline"]

        def cell(d: dict[str, Any]) -> str:
            if d["abstained"]:
                return f"{d['outcome']} (`{d['abstain_reason']}`)"
            return d["outcome"]

        gt = "abstain" if rec["should_abstain"] else "answer"
        a(f"| `{rec['capsule_id']}` | `{rec['archetype']}` | `{rec['variant']}` | {gt} | {cell(g)} | {cell(b)} |")
    a("")
    a("## 7. Honest scope & limitations")
    a("")
    a("- **Synthetic Tier 1 only.** Capsules are seeded simulations with ground truth set by construction. They test "
      "governance boundaries (design, units, licensing, intake, same-cell identity, join semantics) and show nothing "
      "about accuracy on real cohorts. Public-data validation lives in vignettes `01`-`12`.")
    a("- **The baseline is a scripted emulation, not a live LLM agent.** `uncontracted_codeact_baseline` runs real "
      "code (`scipy.stats.ranksums`, Welch t-tests, `pd.merge`, Pearson correlation), but which unguarded pattern it "
      "applies to each capsule was fixed by the authors. Two behaviors are modeled outright: the FASTQ glob that "
      "silently falls through on BCL/POD5 intake, and the positional pairing fallback for unpaired cells. Its "
      "`silent_error_rate` is the error rate *of those patterns on these capsules*. It is not an estimate for any "
      "specific CodeAct product or model; a careful human or agent could avoid every one of these traps.")
    a("- **No live LLM in any arm.** The `aria_governed` claim drafter is a deterministic aggressive stand-in. The "
      "zero-LLM arm uses template prose. Decisions are identical across both by construction of the governance path, "
      "which is the claim being tested; this says nothing about end-to-end LLM planning quality.")
    a(f"- **The capsule mix sets the rates.** `{receipt['n_should_abstain']}/{receipt['n_capsules']}` capsules are "
      "must-abstain traps. Both the baseline's silent-error rate and ARIA's coverage depend on this designed ratio, "
      "not on real-world prevalence.")
    a("- **A frozen-seed 0.0 is not a guarantee.** BH-FDR control permits occasional false-positive calls, and the "
      "seed sweep above shows them. The frozen receipt is a regression anchor, not a bound.")
    a("- **Archetype 5 tests unit and cardinality guards.** It does not test free-form querying: ARIA executes a "
      "fixed typed AST per question, and the baseline runs equivalent pandas code. Where both are valid they agree.")
    a("- **The narrative claim compiler runs on the `aria-env` stack.** When the runner executes in "
      "`aria-rna-env`, the C4/C7 gate runs in a stack-isolated subprocess (JSON IPC); the receipt records "
      "the execution mode in its `runtime` block, which is excluded from `receipt_sha256`.")
    a("")
    return "\n".join(lines)


def _main(argv: list[str]) -> int:
    if "--claim-gate-ipc" in argv:
        requests = json.loads(sys.stdin.read())
        results = evaluate_claim_requests_in_process(requests)
        sys.stdout.write(_IPC_SENTINEL + json.dumps(results, default=str))
        return 0
    print(json.dumps(run_codeact_comparative_suite()["metrics"], indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
