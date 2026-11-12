"""ADR-077 bulk ATAC multi-replicate external R concordance + null benchmark.

Deterministic, fully specified benchmark for the bulk ATAC DA lane
(``aria/scripts/chromatin_bulk_diffacc.py``):

* **Design.** ``n = 4`` biological replicates per condition (``treated`` vs
  ``control``), each sequenced as 2 technical libraries (16 libraries that ARIA
  must collapse to 8 biological replicates, ADR-077 B1), with a balanced
  2-level batch (``batch1``/``batch2``, 2 biological replicates per condition per
  batch) fitted as ``~ batch + condition``.
* **Truth.** 600 consensus peaks: 500 null + 100 planted differential peaks
  (50 gained, 50 lost accessibility; |log2FC| ~ U(1, 2)). A per-peak batch shift
  is applied to every peak so an unadjusted model would be misspecified.
* **Generative model.** Biological replicate abundance is negative-binomial
  (Gamma-Poisson with a mean-dependent dispersion trend); each technical library
  is a Poisson draw around its biological replicate with its own size factor,
  so summing technical libraries is the statistically correct aggregation.

Comparators (Bioconductor, run in ``aria-bench-env``) consume ARIA's own
replicate-level matrix so every method sees identical input: DESeq2 (same
``lfcThreshold`` / ``greaterAbs`` hypothesis as ARIA), edgeR ``glmQLFTest`` and
limma-voom (point-null tests, called at FDR + |log2FC|), plus the
threshold-matched edgeR ``glmTreat`` / limma ``treat`` variants.

The global null is probed with condition-label permutations over the 500 null
peaks (library-level, so technical aggregation is exercised on every run).
Permutations keep the batch balance (2 "treated" biological replicates per
batch), exclude the true labelling and its complement, and move technical
libraries together with their biological replicate.

Pure helpers live here; orchestration across conda envs is in
``scripts/run_bulk_atac_concordance.py``. No fabrication: a failed ARIA or R
step is recorded as such and never imputed.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# ── fixed benchmark specification ────────────────────────────────────────────

SIM_SPEC: dict[str, Any] = {
    "seed": 20261112,
    "conditions": {"reference": "control", "test": "treated"},
    "n_bio_per_condition": 4,
    "n_tech_per_bio": 2,
    "batches": ["batch1", "batch2"],
    "n_null_peaks": 500,
    "n_up_peaks": 50,
    "n_down_peaks": 50,
    "abs_log2fc_range": [1.0, 2.0],
    "base_mean_lognormal": {"meanlog": math.log(150.0), "sdlog": 0.8,
                            "clip": [20.0, 3000.0]},
    "dispersion_trend": {"asymptotic": 0.03, "extra_poisson": 2.0},
    "batch_log2_shift_sd": 0.4,
    "library_size_factor_range": [0.6, 1.4],
}

THRESHOLDS = {"padj": 0.05, "lfc": 0.5}
N_PERMUTATIONS = 20

# Pre-registered acceptance gates (set before the benchmark was executed).
GATES = {
    "lfc_spearman_vs_r_deseq2_min": 0.95,
    "sign_agreement_vs_r_deseq2_min": 0.95,
    "significant_jaccard_vs_r_deseq2_min": 0.80,
    "aria_truth_fdp_max": 0.10,
    "aria_truth_sensitivity_min": 0.70,
    "permutation_fraction_with_any_discovery_max": 0.10,
}


# ── simulation ────────────────────────────────────────────────────────────────

def simulate_multireplicate_bulk_atac(spec: dict | None = None):
    """Return ``(library_counts, library_metadata, truth)``.

    ``library_counts`` is peaks x 16 libraries (``peak_id`` index),
    ``library_metadata`` has ``sample_id/condition/replicate/batch/technical``
    and ``truth`` has ``peak_id/true_log2fc/is_differential/direction/
    base_mean/batch_log2_shift``.
    """
    spec = {**SIM_SPEC, **(spec or {})}
    rng = np.random.default_rng(int(spec["seed"]))
    n_null, n_up, n_down = (int(spec["n_null_peaks"]), int(spec["n_up_peaks"]),
                            int(spec["n_down_peaks"]))
    n_peaks = n_null + n_up + n_down
    peak_ids = [f"chr{1 + i // 100}:{100_000 + (i % 100) * 2_000}-"
                f"{100_000 + (i % 100) * 2_000 + 500}" for i in range(n_peaks)]

    bm = spec["base_mean_lognormal"]
    base_mean = np.clip(rng.lognormal(bm["meanlog"], bm["sdlog"], n_peaks),
                        *bm["clip"])
    dt = spec["dispersion_trend"]
    dispersion = dt["asymptotic"] + dt["extra_poisson"] / base_mean

    lfc = np.zeros(n_peaks)
    order = rng.permutation(n_peaks)
    up_idx, down_idx = order[:n_up], order[n_up:n_up + n_down]
    lo, hi = spec["abs_log2fc_range"]
    lfc[up_idx] = rng.uniform(lo, hi, n_up)
    lfc[down_idx] = -rng.uniform(lo, hi, n_down)
    batch_shift = rng.normal(0.0, float(spec["batch_log2_shift_sd"]), n_peaks)

    ref, test = spec["conditions"]["reference"], spec["conditions"]["test"]
    n_bio, n_tech = int(spec["n_bio_per_condition"]), int(spec["n_tech_per_bio"])
    batches = list(spec["batches"])
    per_batch = n_bio // len(batches)
    sf_lo, sf_hi = spec["library_size_factor_range"]

    columns: dict[str, np.ndarray] = {}
    meta_rows: list[dict[str, Any]] = []
    bio_counter = 0
    for condition in (ref, test):
        for j in range(n_bio):
            bio_counter += 1
            bio_id = f"bio{bio_counter}"
            batch = batches[min(j // per_batch, len(batches) - 1)]
            mean = base_mean * np.power(2.0, lfc * (condition == test))
            mean = mean * np.power(2.0, batch_shift * (batch == batches[-1]))
            shape = 1.0 / dispersion
            abundance = rng.gamma(shape, mean / shape)
            for t in range(n_tech):
                lib = f"{condition}_{bio_id}_tech{t + 1}"
                size_factor = rng.uniform(sf_lo, sf_hi)
                columns[lib] = rng.poisson(abundance * size_factor)
                meta_rows.append({"sample_id": lib, "condition": condition,
                                  "replicate": bio_id, "batch": batch,
                                  "technical": f"tech{t + 1}"})

    counts = pd.DataFrame(columns, index=pd.Index(peak_ids, name="peak_id"))
    meta = pd.DataFrame(meta_rows)
    direction = np.where(lfc > 0, "up", np.where(lfc < 0, "down", "null"))
    truth = pd.DataFrame({
        "peak_id": peak_ids, "true_log2fc": lfc,
        "is_differential": lfc != 0, "direction": direction,
        "base_mean": base_mean, "batch_log2_shift": batch_shift,
    })
    return counts.astype(int), meta, truth


def collapse_technical(library_counts: pd.DataFrame,
                       library_metadata: pd.DataFrame,
                       condition_col: str = "condition",
                       replicate_col: str = "replicate") -> pd.DataFrame:
    """Independent reference implementation of the technical-replicate sum
    (columns ``<condition>::<replicate>``), used to verify ARIA's B1 step."""
    meta = library_metadata.set_index("sample_id")
    keys = (meta[condition_col].astype(str) + "::"
            + meta[replicate_col].astype(str))
    out = library_counts.T.groupby(keys.reindex(library_counts.columns)).sum().T
    return out


def write_counts_tsv(counts: pd.DataFrame, path: Path,
                     id_col: str = "peak_id") -> None:
    out = counts.copy()
    out.insert(0, id_col, out.index.astype(str))
    out.to_csv(path, sep="\t", index=False)


def matrix_sha256(counts: pd.DataFrame) -> str:
    payload = counts.to_csv(sep="\t").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


# ── permutations ──────────────────────────────────────────────────────────────

def batch_balanced_permutations(library_metadata: pd.DataFrame,
                                n_permutations: int = N_PERMUTATIONS,
                                seed: int = SIM_SPEC["seed"],
                                ) -> list[dict[str, str]]:
    """Deterministic batch-balanced condition relabellings of the biological
    replicates. Each keeps the per-batch condition balance, and the true
    labelling and its complement are excluded. Returns ``[{bio_id: label}]``."""
    bio = (library_metadata.drop_duplicates("replicate")
           [["replicate", "condition", "batch"]])
    ref = SIM_SPEC["conditions"]["reference"]
    test = SIM_SPEC["conditions"]["test"]
    true_test = frozenset(bio.loc[bio["condition"] == test, "replicate"])
    true_ref = frozenset(bio.loc[bio["condition"] == ref, "replicate"])
    per_batch_choices = []
    for _batch, rows in bio.groupby("batch", sort=True):
        reps = sorted(rows["replicate"])
        k = int((rows["condition"] == test).sum())
        per_batch_choices.append(list(itertools.combinations(reps, k)))
    candidates = []
    for combo in itertools.product(*per_batch_choices):
        test_set = frozenset(itertools.chain.from_iterable(combo))
        if test_set in (true_test, true_ref):
            continue
        candidates.append(test_set)
    candidates.sort(key=lambda s: sorted(s))
    rng = np.random.default_rng(int(seed) + 1)
    if n_permutations > len(candidates):
        raise ValueError(
            f"requested {n_permutations} permutations but only "
            f"{len(candidates)} batch-balanced relabellings exist")
    picks = rng.choice(len(candidates), size=n_permutations, replace=False)
    out = []
    for idx in sorted(int(i) for i in picks):
        test_set = candidates[idx]
        out.append({rep: (test if rep in test_set else ref)
                    for rep in bio["replicate"]})
    return out


def overlap_with_truth(test_set: frozenset, bio_metadata: pd.DataFrame) -> int:
    test = SIM_SPEC["conditions"]["test"]
    true_test = set(bio_metadata.loc[bio_metadata["condition"] == test,
                                     "replicate"])
    return len(set(test_set) & true_test)


# ── scoring ───────────────────────────────────────────────────────────────────

def _corr(a: pd.Series, b: pd.Series, method: str) -> float | None:
    mask = np.isfinite(a.to_numpy(float)) & np.isfinite(b.to_numpy(float))
    if mask.sum() < 3:
        return None
    value = a[mask].corr(b[mask], method=method)
    return None if value is None or not np.isfinite(value) else round(float(value), 6)


def _neglog10(p: pd.Series) -> pd.Series:
    return -np.log10(p.astype(float).clip(lower=1e-300))


def _set_metrics(called: set, reference: set) -> dict[str, Any]:
    inter = len(called & reference)
    union = len(called | reference)
    return {
        "n_called": len(called), "n_reference": len(reference),
        "n_intersection": inter,
        "jaccard": round(inter / union, 6) if union else None,
        "sensitivity": round(inter / len(reference), 6) if reference else None,
        "precision": round(inter / len(called), 6) if called else None,
    }


def aria_significant(aria: pd.DataFrame) -> set:
    if "significant" in aria.columns:
        sig = aria["significant"].astype(str).str.lower().isin({"true", "1"})
        return set(aria.index[sig])
    return set(aria.index[aria["padj"] < THRESHOLDS["padj"]])


def reference_significant(ref: pd.DataFrame, threshold_matched: bool) -> set:
    padj_ok = ref["padj"].astype(float) < THRESHOLDS["padj"]
    if threshold_matched:
        return set(ref.index[padj_ok])
    lfc_ok = ref["log2FoldChange"].astype(float).abs() >= THRESHOLDS["lfc"]
    return set(ref.index[padj_ok & lfc_ok])


def score_against_reference(aria: pd.DataFrame, ref: pd.DataFrame,
                            threshold_matched: bool) -> dict[str, Any]:
    """Concordance of ARIA (peak-indexed full DA table) vs one R method."""
    common = aria.index.intersection(ref.index)
    a, r = aria.loc[common], ref.loc[common]
    a_lfc = a["log2FoldChange_raw"] if "log2FoldChange_raw" in a else a[
        "log2FoldChange"]
    r_lfc = r["log2FoldChange"].astype(float)
    finite = np.isfinite(a_lfc.to_numpy(float)) & np.isfinite(r_lfc.to_numpy(float))
    sign_all = float((np.sign(a_lfc[finite]) == np.sign(r_lfc[finite])).mean()) \
        if finite.any() else None
    a_sig = aria_significant(a)
    r_sig = reference_significant(r, threshold_matched)
    either = sorted(a_sig | r_sig)
    sign_sig = (float((np.sign(a_lfc.loc[either]) == np.sign(
        r_lfc.loc[either])).mean()) if either else None)
    return {
        "n_common_peaks": int(len(common)),
        "lfc_pearson": _corr(a_lfc.astype(float), r_lfc, "pearson"),
        "lfc_spearman": _corr(a_lfc.astype(float), r_lfc, "spearman"),
        "lfc_spearman_shrunk_aria": _corr(
            a["log2FoldChange"].astype(float), r_lfc, "spearman"),
        "neglog10_padj_spearman": _corr(_neglog10(a["padj"]),
                                        _neglog10(r["padj"]), "spearman"),
        "sign_agreement_all_peaks": (round(sign_all, 6)
                                     if sign_all is not None else None),
        "sign_agreement_significant_union": (round(sign_sig, 6)
                                             if sign_sig is not None else None),
        "significant_overlap": _set_metrics(a_sig, r_sig),
        "reference_significance_rule": (
            f"padj < {THRESHOLDS['padj']} (lfcThreshold-matched test)"
            if threshold_matched else
            f"padj < {THRESHOLDS['padj']} and |log2FC| >= {THRESHOLDS['lfc']}"),
    }


def score_against_truth(called: set, est_lfc: pd.Series,
                        truth: pd.DataFrame) -> dict[str, Any]:
    t = truth.set_index("peak_id")
    true_de = set(t.index[t["is_differential"]])
    tp = called & true_de
    fp = called - true_de
    direction_ok = sum(
        1 for p in tp
        if np.sign(float(est_lfc.get(p, 0.0))) == np.sign(t.loc[p, "true_log2fc"]))
    return {
        "n_called": len(called), "n_true_differential": len(true_de),
        "true_positives": len(tp), "false_positives": len(fp),
        "false_negatives": len(true_de - called),
        "sensitivity": round(len(tp) / len(true_de), 6) if true_de else None,
        "precision": round(len(tp) / len(called), 6) if called else None,
        "false_discovery_proportion": (round(len(fp) / len(called), 6)
                                       if called else 0.0),
        "direction_correct_among_tp": (round(direction_ok / len(tp), 6)
                                       if tp else None),
    }


def summarise_permutations(per_perm: list[dict[str, Any]],
                           n_null_peaks: int) -> dict[str, Any]:
    ran = [p for p in per_perm if p.get("status") == "success"]
    n_sig = [int(p.get("n_sig", 0)) for p in ran]
    pvals = np.concatenate([np.asarray(p.get("pvalues") or [], dtype=float)
                            for p in ran]) if ran else np.array([])
    pvals = pvals[np.isfinite(pvals)]
    total_tests = sum(int(p.get("n_tested", 0)) for p in ran)
    any_disc = sum(1 for n in n_sig if n > 0)
    return {
        "n_permutations_requested": len(per_perm),
        "n_permutations_completed": len(ran),
        "n_null_peaks": n_null_peaks,
        "total_tests": total_tests,
        "total_false_discoveries": int(sum(n_sig)),
        "mean_false_discoveries_per_permutation": (
            round(float(np.mean(n_sig)), 6) if n_sig else None),
        "max_false_discoveries_in_one_permutation": max(n_sig) if n_sig else None,
        "fraction_permutations_with_any_discovery": (
            round(any_disc / len(ran), 6) if ran else None),
        "per_test_false_positive_rate": (
            round(sum(n_sig) / total_tests, 8) if total_tests else None),
        "pooled_raw_pvalue_fraction_below_0.05": (
            round(float((pvals < 0.05).mean()), 6) if pvals.size else None),
        "pooled_raw_pvalue_fraction_below_0.01": (
            round(float((pvals < 0.01).mean()), 6) if pvals.size else None),
        "n_pooled_pvalues": int(pvals.size),
        "note": ("Raw p-values come from the |log2FC| > lfc_min (greaterAbs) "
                 "Wald test, which is conservative under the point null, so "
                 "pooled fractions are expected to fall below the nominal "
                 "level."),
    }


def evaluate_gates(aria_vs_r: dict[str, Any], aria_truth: dict[str, Any],
                   null_summary: dict[str, Any],
                   aggregation_exact: bool) -> dict[str, Any]:
    deseq = aria_vs_r.get("deseq2") or {}
    overlap = deseq.get("significant_overlap") or {}

    def _ge(value, threshold):
        return value is not None and value >= threshold

    def _le(value, threshold):
        return value is not None and value <= threshold

    checks = {
        "technical_aggregation_exact": bool(aggregation_exact),
        "lfc_spearman_vs_r_deseq2": _ge(
            deseq.get("lfc_spearman"), GATES["lfc_spearman_vs_r_deseq2_min"]),
        "sign_agreement_vs_r_deseq2": _ge(
            deseq.get("sign_agreement_all_peaks"),
            GATES["sign_agreement_vs_r_deseq2_min"]),
        "significant_jaccard_vs_r_deseq2": _ge(
            overlap.get("jaccard"), GATES["significant_jaccard_vs_r_deseq2_min"]),
        "aria_truth_fdp": _le(aria_truth.get("false_discovery_proportion"),
                              GATES["aria_truth_fdp_max"]),
        "aria_truth_sensitivity": _ge(aria_truth.get("sensitivity"),
                                      GATES["aria_truth_sensitivity_min"]),
        "permutation_null_control": _le(
            null_summary.get("fraction_permutations_with_any_discovery"),
            GATES["permutation_fraction_with_any_discovery_max"]),
    }
    return {"thresholds": dict(GATES), "checks": checks,
            "all_passed": all(checks.values())}


# ── public-artifact hygiene ───────────────────────────────────────────────────

def scrub_paths(obj: Any, repo_root: Path) -> Any:
    """Make every path-like string repo-relative (or a basename when outside
    the repo) so the public receipt carries no machine-absolute path."""
    root = str(Path(repo_root).resolve())
    if isinstance(obj, dict):
        return {k: scrub_paths(v, repo_root) for k, v in obj.items()}
    if isinstance(obj, list):
        return [scrub_paths(v, repo_root) for v in obj]
    if isinstance(obj, str) and obj.startswith("/") and "/" in obj[1:]:
        if obj.startswith(root + "/"):
            return obj[len(root) + 1:]
        return Path(obj).name
    return obj


# ── embedded Bioconductor comparator script (~ batch + condition) ────────────

R_COMPARATORS = r"""
suppressPackageStartupMessages({
  library(jsonlite); library(DESeq2); library(edgeR); library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]; meta_path <- args[2]; out_dir <- args[3]
padj_thr <- as.numeric(args[4]); lfc_thr <- as.numeric(args[5])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cd <- read.delim(counts_path, check.names = FALSE)
mat <- as.matrix(cd[, -1, drop = FALSE]); rownames(mat) <- cd[[1]]
mode(mat) <- "integer"
meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
meta <- meta[colnames(mat), , drop = FALSE]
meta$condition <- relevel(factor(meta$condition), ref = "control")
meta$batch <- factor(meta$batch)

write_res <- function(peak, lfc, p, padj, name) {
  df <- data.frame(peak = peak, log2FoldChange = lfc, pvalue = p, padj = padj)
  path <- file.path(out_dir, paste0(name, ".tsv"))
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(status = "success", results_tsv = basename(path), n_tested = nrow(df))
}
methods <- list()
run <- function(name, expr) {
  methods[[name]] <<- tryCatch(expr, error = function(e)
    list(status = "error", error_type = class(e)[1], details = conditionMessage(e)))
}

run("deseq2", {
  dds <- DESeqDataSetFromMatrix(mat, meta, design = ~ batch + condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treated", "control"),
                 lfcThreshold = lfc_thr, altHypothesis = "greaterAbs",
                 alpha = padj_thr)
  write_res(rownames(res), res$log2FoldChange, res$pvalue, res$padj, "deseq2")
})

y <- DGEList(counts = mat)
y <- calcNormFactors(y)
design <- model.matrix(~ batch + condition, data = meta)
coef <- "conditiontreated"
run("edgeR_QLF", {
  yd <- estimateDisp(y, design)
  fit <- glmQLFit(yd, design)
  q <- glmQLFTest(fit, coef = coef)
  tab <- topTags(q, n = Inf, sort.by = "none")$table
  write_res(rownames(tab), tab$logFC, tab$PValue, tab$FDR, "edgeR_QLF")
})
run("edgeR_glmTreat", {
  yd <- estimateDisp(y, design)
  fit <- glmQLFit(yd, design)
  q <- glmTreat(fit, coef = coef, lfc = lfc_thr)
  tab <- topTags(q, n = Inf, sort.by = "none")$table
  write_res(rownames(tab), tab$logFC, tab$PValue, tab$FDR, "edgeR_glmTreat")
})
v <- voom(y, design, plot = FALSE)
vfit <- lmFit(v, design)
run("limma_voom", {
  fit <- eBayes(vfit)
  tab <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  write_res(rownames(tab), tab$logFC, tab$P.Value, tab$adj.P.Val, "limma_voom")
})
run("limma_treat", {
  fit <- treat(vfit, lfc = lfc_thr)
  tab <- topTreat(fit, coef = coef, number = Inf, sort.by = "none")
  write_res(rownames(tab), tab$logFC, tab$P.Value, tab$adj.P.Val, "limma_treat")
})

write_json(list(
  status = "success", methods = methods,
  design = "~ batch + condition",
  versions = list(
    R = paste(R.version$major, R.version$minor, sep = "."),
    DESeq2 = as.character(packageVersion("DESeq2")),
    edgeR = as.character(packageVersion("edgeR")),
    limma = as.character(packageVersion("limma")))
), file.path(out_dir, "r_comparators.json"), auto_unbox = TRUE, pretty = TRUE,
   null = "null", digits = NA)
"""

# Which R methods test the same |log2FC| > lfc_min hypothesis as ARIA.
R_METHOD_THRESHOLD_MATCHED = {
    "deseq2": True,
    "edgeR_QLF": False,
    "edgeR_glmTreat": True,
    "limma_voom": False,
    "limma_treat": True,
}
