"""ADR-078 differential TF footprinting: external R + live TOBIAS concordance.

Deterministic, fully specified benchmark for the replicate-level footprint
inference that promotes TOBIAS differential TF footprinting to production
(``aria/scripts/chromatin_footprint_tobias.py``):

* **Suite 1 — independent Welch.** ``n = 4`` vs ``4`` biological replicates over
  40 motifs (10 planted true effects, 30 nulls). ARIA's
  ``infer_replicate_footprints`` vs R ``stats::t.test(var.equal = FALSE)`` +
  ``stats::p.adjust(method = "BH")``.
* **Suite 2 — paired donors.** 4 donors measured in both conditions, with large
  donor offsets that only a paired analysis removes. ARIA vs R
  ``stats::t.test(paired = TRUE)`` + BH.
* **Suite 3 — batch-adjusted OLS.** 4 vs 4 replicates balanced over 2 batches
  with a large +/-2.5 batch offset. ARIA vs R ``stats::lm(score ~ batch + group)``
  (``groupA`` coefficient) + BH, plus a fully confounded negative control
  (batch1 = all A, batch2 = all B) that must stay ``descriptive_only``.
* **Suite 4 — null calibration.** 20 pure-noise cohorts (4 vs 4, 50 motifs):
  the number of BH discoveries at FDR <= 0.05 is reported per cohort.
* **Suite 5 — live TOBIAS 0.17.3.** A deterministic 40 kb mini-genome with 20
  peaks carrying planted motif instances; Tn5 cut sites are simulated with a
  protected, depth-neutral CEBPB footprint in 3 independent ``Treated``
  replicates vs 3 independent ``Control`` replicates. The production bulk
  driver runs ATACorrect -> ScoreBigwig -> BINDetect in ``aria-tobias-env`` and
  R re-tests the per-replicate ``*_mean_score`` columns of the live
  ``bindetect_results.txt``.

The acceptance gates in ``GATES`` were fixed before the benchmark was executed;
all seeds are fixed in ``SPEC``. No fabrication: a failed ARIA, R or TOBIAS step
is recorded as such and the gate fails. The returned receipt contains no
absolute filesystem paths (only names relative to the run directory).
"""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aria.scripts.chromatin_footprint_tobias import (  # noqa: E402
    infer_replicate_footprints,
)

BENCHMARK_ID = "adr078_tobias_replicate_concordance"
ALPHA = 0.05

SPEC: dict[str, Any] = {
    "seed": 20261112,
    "welch": {"n_per_group": 4, "n_motifs": 40, "n_planted": 10,
              "planted_effect": 1.5, "noise_sd": 0.25},
    "paired": {"n_donors": 4, "n_motifs": 40, "n_planted": 10,
               "planted_effect": 1.0, "noise_sd": 0.15, "donor_sd": 2.0},
    "batch": {"n_per_group": 4, "n_motifs": 40, "n_planted": 10,
              "planted_effect": 1.5, "noise_sd": 0.25, "batch_offset": 2.5,
              "batches": ["batch1", "batch2"]},
    "null": {"n_cohorts": 20, "n_per_group": 4, "n_motifs": 50, "noise_sd": 1.0},
    "live_tobias": {"genome_length": 40_000, "n_peaks": 20, "peak_width": 500,
                    "n_replicates_per_condition": 3, "reads_per_replicate": 60_000,
                    "read_length": 36, "groups": ["Treated", "Control"],
                    "protected_tf": "CEBPB"},
}

# Pre-registered acceptance gates (fixed before execution).
GATES: dict[str, Any] = {
    "effect_pearson_r_min": 0.999999,
    "max_abs_padj_diff_max": 1e-10,
    "significant_jaccard_min": 1.0,
    "null_false_positive_cohorts_max": 0,
    "live_max_abs_padj_diff_max": 1e-10,
}

_ABS_PREFIXES = ("/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/",
                 "/data/", "/var/", "/opt/")

# ── environment helpers ──────────────────────────────────────────────────────


def _conda_env_prefix(env: str) -> Path | None:
    """Resolve a conda env prefix by name (no activation needed)."""
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if not conda:
        return None
    try:
        out = subprocess.run([conda, "env", "list", "--json"], check=True,
                             capture_output=True, text=True, timeout=120)
        for prefix in json.loads(out.stdout).get("envs", []):
            if Path(prefix).name == env:
                return Path(prefix)
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    return None


def _env_tool(env: str, tool: str) -> list[str] | None:
    """Command prefix running ``tool`` from ``env`` (direct bin, else conda run)."""
    prefix = _conda_env_prefix(env)
    if prefix is not None and (prefix / "bin" / tool).is_file():
        return [str(prefix / "bin" / tool)]
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if conda:
        return [conda, "run", "--no-capture-output", "-n", env, tool]
    return None


def _run_tool(cmd: list[str], *, cwd: Path | None = None,
              timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, capture_output=True, text=True,
                          cwd=str(cwd) if cwd else None, timeout=timeout)


# ── shared statistics helpers ────────────────────────────────────────────────

_R_SCRIPT = r"""
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; infile <- args[2]; designfile <- args[3]; outfile <- args[4]
m <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)
d <- read.delim(designfile, stringsAsFactors = FALSE, colClasses = "character")
a <- d[d$group == "A", ]; b <- d[d$group == "B", ]
if (mode == "paired") { a <- a[order(a$pair), ]; b <- b[order(b$pair), ] }
est <- numeric(nrow(m)); p <- numeric(nrow(m))
for (i in seq_len(nrow(m))) {
  xa <- as.numeric(m[i, a$sample]); xb <- as.numeric(m[i, b$sample])
  if (mode == "welch") {
    tt <- stats::t.test(xa, xb, var.equal = FALSE)
    est[i] <- mean(xa) - mean(xb); p[i] <- tt$p.value
  } else if (mode == "paired") {
    tt <- stats::t.test(xa, xb, paired = TRUE)
    est[i] <- unname(tt$estimate); p[i] <- tt$p.value
  } else if (mode == "lm") {
    df <- data.frame(score = c(xa, xb),
                     group = factor(c(rep("A", length(xa)), rep("B", length(xb))),
                                    levels = c("B", "A")),
                     batch = factor(c(a$batch, b$batch)))
    co <- summary(stats::lm(score ~ batch + group, data = df))$coefficients
    est[i] <- co["groupA", "Estimate"]; p[i] <- co["groupA", "Pr(>|t|)"]
  } else stop(paste("unknown mode", mode))
}
padj <- stats::p.adjust(p, method = "BH")
out <- data.frame(motif = m$motif, estimate = sprintf("%.17g", est),
                  pvalue = sprintf("%.17g", p), padj = sprintf("%.17g", padj))
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
"""


def _r_version(rscript: list[str]) -> str | None:
    try:
        out = _run_tool([*rscript, "-e", "cat(R.version.string)"], timeout=300)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def _run_r(rscript: list[str], mode: str, matrix: dict[str, dict[str, float]],
           design: list[dict[str, str]], workdir: Path, tag: str
           ) -> dict[str, dict[str, float]]:
    """Run the R comparator; returns ``{motif: {estimate, pvalue, padj}}``."""
    workdir.mkdir(parents=True, exist_ok=True)
    script = workdir / "footprint_r_comparator.R"
    script.write_text(_R_SCRIPT, encoding="utf-8")
    samples = [row["sample"] for row in design]
    mat_path = workdir / f"{tag}_matrix.tsv"
    with open(mat_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["motif", *samples])
        for motif, scores in matrix.items():
            writer.writerow([motif, *[repr(float(scores[s])) for s in samples]])
    design_path = workdir / f"{tag}_design.tsv"
    with open(design_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "group", "pair", "batch"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in design:
            writer.writerow({k: row.get(k, "") for k in writer.fieldnames})
    out_path = workdir / f"{tag}_r_results.tsv"
    _run_tool([*rscript, str(script), mode, str(mat_path), str(design_path),
               str(out_path)], timeout=1800)
    results = {}
    with open(out_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            results[row["motif"]] = {k: float(row[k])
                                     for k in ("estimate", "pvalue", "padj")}
    return results


def _pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 2:
        return None
    xa, ya = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if np.std(xa) == 0 or np.std(ya) == 0:
        return None
    return float(np.corrcoef(xa, ya)[0, 1])


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def _compare(aria_rows: dict[str, dict[str, float]],
             r_rows: dict[str, dict[str, float]],
             planted: set[str] | None = None) -> dict[str, Any]:
    """Concordance metrics between ARIA and R per-motif results."""
    shared = sorted(set(aria_rows) & set(r_rows))
    eff_a = [aria_rows[m]["change"] for m in shared]
    eff_r = [r_rows[m]["estimate"] for m in shared]
    padj_diff = [abs(aria_rows[m]["padj"] - r_rows[m]["padj"]) for m in shared]
    sig_a = {m for m in shared if aria_rows[m]["padj"] <= ALPHA}
    sig_r = {m for m in shared if r_rows[m]["padj"] <= ALPHA}
    out: dict[str, Any] = {
        "n_motifs_aria": len(aria_rows),
        "n_motifs_r": len(r_rows),
        "n_motifs_shared": len(shared),
        "effect_pearson_r": _pearson(eff_a, eff_r),
        "max_abs_effect_diff": max((abs(a - b) for a, b in zip(eff_a, eff_r)),
                                   default=None),
        "max_abs_padj_diff": max(padj_diff, default=None),
        "n_significant_aria": len(sig_a),
        "n_significant_r": len(sig_r),
        "significant_jaccard": _jaccard(sig_a, sig_r),
    }
    if planted is not None:
        out["n_planted"] = len(planted)
        out["planted_recovered_aria"] = len(sig_a & planted)
        out["false_discoveries_aria"] = len(sig_a - planted)
    return out


def _concordance_gates(metrics: dict[str, Any]) -> dict[str, bool]:
    r = metrics.get("effect_pearson_r")
    d = metrics.get("max_abs_padj_diff")
    return {
        "effect_pearson_r": r is not None and r >= GATES["effect_pearson_r_min"],
        "max_abs_padj_diff": d is not None and d < GATES["max_abs_padj_diff_max"],
        "significant_jaccard": (metrics.get("significant_jaccard", 0.0)
                                >= GATES["significant_jaccard_min"]),
        "shared_all_motifs": (metrics.get("n_motifs_shared")
                              == metrics.get("n_motifs_aria")
                              == metrics.get("n_motifs_r")),
    }


def _aria_rows(summary: dict[str, Any]) -> dict[str, dict[str, float]]:
    return {str(row["motif_id"]): {"change": float(row["change"]),
                                    "padj": float(row["padj"]),
                                    "pvalue": float(row["pvalue"])}
            for row in summary.get("all_results") or []}


# ── suites 1-4: synthetic replicate score tables ─────────────────────────────


def _planted_table(rng: np.random.Generator, reps_a: list[str], reps_b: list[str],
                   n_motifs: int, n_planted: int, effect: float, noise_sd: float,
                   offsets: dict[str, float] | None = None,
                   ) -> tuple[dict[str, dict[str, Any]], set[str]]:
    table: dict[str, dict[str, Any]] = {}
    planted: set[str] = set()
    for idx in range(n_motifs):
        motif = f"MOTIF_{idx:03d}"
        sign = 1.0 if idx % 2 == 0 else -1.0
        delta = sign * effect if idx < n_planted else 0.0
        if idx < n_planted:
            planted.add(motif)
        base = float(rng.normal(0.0, 0.5))
        scores = {}
        for rep in reps_a:
            scores[rep] = (base + delta + (offsets or {}).get(rep, 0.0)
                           + float(rng.normal(0.0, noise_sd)))
        for rep in reps_b:
            scores[rep] = (base + (offsets or {}).get(rep, 0.0)
                           + float(rng.normal(0.0, noise_sd)))
        table[motif] = {"tf": f"TF{idx:03d}", "motif_id": motif, "scores": scores,
                        "n_sites": 100 + idx}
    return table, planted


def _matrix(table: dict[str, dict[str, Any]]) -> dict[str, dict[str, float]]:
    return {motif: dict(rec["scores"]) for motif, rec in table.items()}


def _suite_welch(rscript: list[str], workdir: Path) -> dict[str, Any]:
    spec = SPEC["welch"]
    rng = np.random.default_rng(SPEC["seed"] + 1)
    n = spec["n_per_group"]
    reps_a = [f"A_r{i + 1}" for i in range(n)]
    reps_b = [f"B_r{i + 1}" for i in range(n)]
    table, planted = _planted_table(rng, reps_a, reps_b, spec["n_motifs"],
                                    spec["n_planted"], spec["planted_effect"],
                                    spec["noise_sd"])
    groups = {**{r: "A" for r in reps_a}, **{r: "B" for r in reps_b}}
    summary = infer_replicate_footprints(table, groups, "A", "B", alpha=ALPHA,
                                         top_n=1000)
    r_rows = _run_r(rscript, "welch", _matrix(table),
                    [{"sample": r, "group": groups[r]} for r in reps_a + reps_b],
                    workdir, "welch")
    metrics = _compare(_aria_rows(summary), r_rows, planted)
    metrics["aria_test"] = (summary.get("inference") or {}).get("test")
    metrics["r_call"] = "stats::t.test(x_a, x_b, var.equal = FALSE); p.adjust(BH)"
    gates = _concordance_gates(metrics)
    gates["aria_inference_success"] = (
        (summary.get("inference") or {}).get("status") == "success")
    return {"metrics": metrics, "gates": gates, "pass": all(gates.values())}


def _suite_paired(rscript: list[str], workdir: Path) -> dict[str, Any]:
    spec = SPEC["paired"]
    rng = np.random.default_rng(SPEC["seed"] + 2)
    donors = [f"donor{i + 1}" for i in range(spec["n_donors"])]
    reps_a = [f"A_{d}" for d in donors]
    reps_b = [f"B_{d}" for d in donors]
    donor_offset = {d: float(rng.normal(0.0, spec["donor_sd"])) for d in donors}
    offsets = {**{f"A_{d}": v for d, v in donor_offset.items()},
               **{f"B_{d}": v for d, v in donor_offset.items()}}
    table, planted = _planted_table(rng, reps_a, reps_b, spec["n_motifs"],
                                    spec["n_planted"], spec["planted_effect"],
                                    spec["noise_sd"], offsets=offsets)
    groups = {**{r: "A" for r in reps_a}, **{r: "B" for r in reps_b}}
    replicate_ids = {f"{g}_{d}": d for g in ("A", "B") for d in donors}
    summary = infer_replicate_footprints(table, groups, "A", "B", alpha=ALPHA,
                                         top_n=1000, replicate_ids=replicate_ids)
    design = ([{"sample": f"A_{d}", "group": "A", "pair": d} for d in donors]
              + [{"sample": f"B_{d}", "group": "B", "pair": d} for d in donors])
    r_rows = _run_r(rscript, "paired", _matrix(table), design, workdir, "paired")
    metrics = _compare(_aria_rows(summary), r_rows, planted)
    inference = summary.get("inference") or {}
    metrics["aria_test"] = inference.get("test")
    metrics["aria_pairing_status"] = inference.get("pairing_status")
    metrics["r_call"] = "stats::t.test(x_a, x_b, paired = TRUE); p.adjust(BH)"
    gates = _concordance_gates(metrics)
    gates["aria_inference_success"] = inference.get("status") == "success"
    gates["aria_used_paired_test"] = str(inference.get("test", "")).startswith("paired")
    return {"metrics": metrics, "gates": gates, "pass": all(gates.values())}


def _suite_batch(rscript: list[str], workdir: Path) -> dict[str, Any]:
    spec = SPEC["batch"]
    rng = np.random.default_rng(SPEC["seed"] + 3)
    n = spec["n_per_group"]
    b1, b2 = spec["batches"]
    reps_a = [f"A_r{i + 1}" for i in range(n)]
    reps_b = [f"B_r{i + 1}" for i in range(n)]
    batches = {rep: (b1 if i < n // 2 else b2)
               for reps in (reps_a, reps_b) for i, rep in enumerate(reps)}
    offsets = {rep: (spec["batch_offset"] if batches[rep] == b1
                     else -spec["batch_offset"]) for rep in batches}
    table, planted = _planted_table(rng, reps_a, reps_b, spec["n_motifs"],
                                    spec["n_planted"], spec["planted_effect"],
                                    spec["noise_sd"], offsets=offsets)
    groups = {**{r: "A" for r in reps_a}, **{r: "B" for r in reps_b}}
    summary = infer_replicate_footprints(table, groups, "A", "B", alpha=ALPHA,
                                         top_n=1000, replicate_batches=batches)
    design = [{"sample": r, "group": groups[r], "batch": batches[r]}
              for r in reps_a + reps_b]
    r_rows = _run_r(rscript, "lm", _matrix(table), design, workdir, "batch")
    metrics = _compare(_aria_rows(summary), r_rows, planted)
    inference = summary.get("inference") or {}
    metrics["aria_test"] = inference.get("test")
    metrics["r_call"] = ("summary(stats::lm(score ~ batch + group))"
                         "$coefficients['groupA', ]; p.adjust(BH)")

    confounded_batches = {**{r: b1 for r in reps_a}, **{r: b2 for r in reps_b}}
    confounded = infer_replicate_footprints(
        table, groups, "A", "B", alpha=ALPHA, replicate_batches=confounded_batches)
    conf_inf = confounded.get("inference") or {}
    metrics["confounded_control"] = {"status": conf_inf.get("status"),
                                      "reason": conf_inf.get("reason")}
    gates = _concordance_gates(metrics)
    gates["aria_inference_success"] = inference.get("status") == "success"
    gates["aria_batch_adjusted"] = bool(inference.get("batch_adjusted"))
    gates["confounded_is_descriptive_only"] = (
        conf_inf.get("status") == "descriptive_only"
        and conf_inf.get("reason") == "batch_condition_confounded")
    return {"metrics": metrics, "gates": gates, "pass": all(gates.values())}


def _suite_null() -> dict[str, Any]:
    spec = SPEC["null"]
    n = spec["n_per_group"]
    reps_a = [f"A_r{i + 1}" for i in range(n)]
    reps_b = [f"B_r{i + 1}" for i in range(n)]
    groups = {**{r: "A" for r in reps_a}, **{r: "B" for r in reps_b}}
    per_cohort = []
    for cohort in range(spec["n_cohorts"]):
        rng = np.random.default_rng(SPEC["seed"] + 1000 + cohort)
        table, _ = _planted_table(rng, reps_a, reps_b, spec["n_motifs"], 0, 0.0,
                                  spec["noise_sd"])
        summary = infer_replicate_footprints(table, groups, "A", "B", alpha=ALPHA)
        per_cohort.append(int(summary.get("n_significant") or 0))
    fp_cohorts = sum(1 for k in per_cohort if k > 0)
    metrics = {
        "n_cohorts": spec["n_cohorts"],
        "n_motifs_per_cohort": spec["n_motifs"],
        "false_discoveries_per_cohort": per_cohort,
        "total_false_discoveries": int(sum(per_cohort)),
        "false_positive_cohorts": fp_cohorts,
        "nominal_cohort_fwer_under_global_null": ALPHA,
    }
    gates = {"false_positive_cohorts": (
        fp_cohorts <= GATES["null_false_positive_cohorts_max"])}
    return {"metrics": metrics, "gates": gates, "pass": all(gates.values())}


# ── suite 5: live TOBIAS on a synthetic mini-genome ──────────────────────────

_MOTIFS = {
    # name: consensus; the PWM puts 0.97 on the consensus base.
    "CEBPB": "ATTGCGCAAT",
    "NFKB_BG": "GGGACTTTCC",
    "AP1_BG": "TGACTCAT",
}


def _write_meme(path: Path) -> None:
    lines = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies",
             "A 0.25 C 0.25 G 0.25 T 0.25", ""]
    for idx, (name, consensus) in enumerate(_MOTIFS.items(), start=1):
        lines += [f"MOTIF MA{9000 + idx}.1 {name}", "",
                  f"letter-probability matrix: alength= 4 w= {len(consensus)} "
                  f"nsites= 20 E= 0", ]
        for base in consensus:
            lines.append("  " + "  ".join(
                "0.970" if base == b else "0.010" for b in "ACGT"))
        lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _synthesize_genome(workdir: Path, rng: np.random.Generator
                       ) -> tuple[str, list[tuple[int, int]], dict[str, list[int]]]:
    spec = SPEC["live_tobias"]
    length, width = spec["genome_length"], spec["peak_width"]
    seq = list(rng.choice(list("ACGT"), size=length))
    peaks, sites = [], {name: [] for name in _MOTIFS}
    spacing = length // (spec["n_peaks"] + 1)
    for k in range(spec["n_peaks"]):
        start = spacing * (k + 1) - width // 2
        peaks.append((start, start + width))
        center = start + width // 2
        for name, offset in (("CEBPB", 0), ("NFKB_BG", 120), ("AP1_BG", -150)):
            motif = _MOTIFS[name]
            pos = center + offset - len(motif) // 2
            seq[pos:pos + len(motif)] = list(motif)
            sites[name].append(pos + len(motif) // 2)
    fasta = workdir / "mini_genome.fa"
    body = "".join(seq)
    with open(fasta, "w", encoding="utf-8") as fh:
        fh.write(">chrS\n")
        for i in range(0, length, 60):
            fh.write(body[i:i + 60] + "\n")
    (workdir / "mini_genome.fa.fai").write_text(
        f"chrS\t{length}\t6\t60\t61\n", encoding="utf-8")
    with open(workdir / "peaks.bed", "w", encoding="utf-8") as fh:
        for i, (s, e) in enumerate(peaks):
            fh.write(f"chrS\t{s}\t{e}\tpeak{i + 1}\n")
    return body, peaks, sites


def _cut_weights(length: int, peaks: list[tuple[int, int]],
                 cebpb_sites: list[int], protection: float, flank_gain: float
                 ) -> np.ndarray:
    weights = np.full(length, 0.02)
    for s, e in peaks:
        weights[s:e] = 1.0
    for c in cebpb_sites:
        weights[c - 45:c - 15] *= flank_gain
        weights[c + 15:c + 45] *= flank_gain
        weights[c - 15:c + 15] *= protection
    return weights / weights.sum()


def _write_replicate_bam(samtools: list[str], body: str, weights: np.ndarray,
                         rng: np.random.Generator, n_reads: int, dest: Path) -> None:
    read_len = SPEC["live_tobias"]["read_length"]
    length = len(body)
    cuts = rng.choice(length, size=n_reads, p=weights)
    reverse = rng.random(n_reads) < 0.5
    rows = []
    for i, (cut, rev) in enumerate(zip(cuts.tolist(), reverse.tolist())):
        start = cut - read_len + 1 if rev else cut
        if start < 0 or start + read_len > length:
            continue
        rows.append((start, rev, i))
    rows.sort()
    sam = dest.with_suffix(".sam")
    with open(sam, "w", encoding="utf-8") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n")
        fh.write(f"@SQ\tSN:chrS\tLN:{length}\n")
        for start, rev, i in rows:
            fh.write(f"r{i}\t{16 if rev else 0}\tchrS\t{start + 1}\t60\t{read_len}M\t*"
                     f"\t0\t0\t{body[start:start + read_len]}\t*\n")
    _run_tool([*samtools, "view", "-b", "-o", str(dest), str(sam)])
    _run_tool([*samtools, "index", str(dest)])
    sam.unlink()


def _rel(path: Any, base: Path) -> str | None:
    if not path:
        return None
    try:
        return str(Path(str(path)).resolve().relative_to(base.resolve()))
    except ValueError:
        return Path(str(path)).name


def _suite_live_tobias(rscript: list[str], workdir: Path,
                       tobias_env: str) -> dict[str, Any]:
    spec = SPEC["live_tobias"]
    samtools = _env_tool(tobias_env, "samtools")
    tobias = _env_tool(tobias_env, "TOBIAS")
    if samtools is None or tobias is None:
        return {"status": "not_run", "reason": f"{tobias_env} not resolvable",
                "gates": {"ran": False}, "pass": False}
    try:
        version = _run_tool([*tobias, "--version"], timeout=300).stdout.strip()
    except (OSError, subprocess.SubprocessError) as exc:
        return {"status": "not_run", "reason": f"TOBIAS unavailable: {exc}",
                "gates": {"ran": False}, "pass": False}

    live = workdir / "live_tobias"
    live.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SPEC["seed"] + 5)
    body, peaks, sites = _synthesize_genome(live, rng)
    _write_meme(live / "motifs.meme")
    treated, control = spec["groups"]
    replicate_bams: dict[str, dict[str, list[str]]] = {treated: {}, control: {}}
    for group in (treated, control):
        for r in range(spec["n_replicates_per_condition"]):
            if group == treated:
                protection = float(rng.uniform(0.15, 0.35))
            else:
                protection = float(rng.uniform(0.9, 1.1))
            # Depth-neutral footprint: the cuts removed from the 30 bp core are
            # redistributed over the 60 bp of flanks, so coverage at the
            # unprotected background motifs is unchanged between conditions.
            flank_gain = 1.0 + (1.0 - protection) / 2.0
            weights = _cut_weights(len(body), peaks, sites["CEBPB"],
                                   protection, flank_gain)
            n_reads = int(spec["reads_per_replicate"] * rng.uniform(0.85, 1.15))
            # Independent biological replicates carry condition-unique IDs; a
            # shared ID would declare a donor-matched (paired) design.
            rep_id = f"{group[0]}{r + 1}"
            bam = live / f"{group}_{rep_id}.bam"
            _write_replicate_bam(samtools, body, weights, rng, n_reads, bam)
            replicate_bams[group][rep_id] = [str(bam)]

    from aria.utils.environment_manager import EnvironmentManager

    out_dir = live / "footprinting"
    params = {
        "mode": "bulk", "data_type": "bulk_ATAC",
        "replicate_bams": replicate_bams,
        "genome_fasta": str(live / "mini_genome.fa"),
        "genome_fai": str(live / "mini_genome.fa.fai"),
        "peaks_bed": str(live / "peaks.bed"),
        "motif_meme": str(live / "motifs.meme"),
        "group_a": treated, "group_b": control,
        "output_dir": str(out_dir), "cores": 2,
        "min_replicates_per_condition": 3, "footprint_fdr": ALPHA,
        "make_aggregate_plots": False, "make_volcano_plot": True,
        "resume": False,
    }
    manager = EnvironmentManager(workspace_dir=str(workdir / "aria_workspace"))
    result = manager.run_in_stack(
        stack="tobias", script_path="aria/scripts/chromatin_footprint_tobias.py",
        params=params)
    summary = (result or {}).get("differential_summary") or {}
    inference = summary.get("inference") or {}
    volcano = (result or {}).get("volcano_plot") or {}
    manifest = (result or {}).get("stage_manifest_path")
    metrics: dict[str, Any] = {
        "tobias_version": version,
        "ran": bool((result or {}).get("ran")),
        "status": (result or {}).get("status"),
        "reason": (result or {}).get("reason"),
        "method": (result or {}).get("method"),
        "inference_status": inference.get("status"),
        "inference_test": inference.get("test"),
        "pairing_status": inference.get("pairing_status"),
        "protected_tf": spec["protected_tf"],
        "replicates_per_condition": inference.get("replicates_per_condition"),
        "n_motifs_tested": summary.get("n_motifs_tested"),
        "n_significant": summary.get("n_significant"),
        "volcano_png": _rel(volcano.get("png"), workdir),
        "volcano_svg": _rel(volcano.get("svg"), workdir),
        "stage_manifest": _rel(manifest, workdir),
        "results_table": _rel(summary.get("results_table"), workdir),
        "bindetect_results": _rel((result or {}).get("bindetect_results"), workdir),
    }
    gates = {
        "ran": metrics["ran"],
        "inference_success": inference.get("status") == "success",
        "independent_replicates_use_welch": str(inference.get("test", "")).startswith(
            "welch"),
        "volcano_png_svg": bool(volcano.get("png") and Path(volcano["png"]).is_file()
                                and volcano.get("svg")
                                and Path(volcano["svg"]).is_file()),
        "stage_manifest": bool(manifest and Path(manifest).is_file()),
    }
    results_table = summary.get("results_table")
    bindetect = (result or {}).get("bindetect_results")
    replicate_groups = (result or {}).get("replicate_groups") or {}
    if (gates["ran"] and results_table and Path(results_table).is_file()
            and bindetect and Path(bindetect).is_file() and replicate_groups):
        matrix: dict[str, dict[str, float]] = {}
        with open(bindetect, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                motif_id = str(row.get("motif_id") or row.get("name"))
                try:
                    scores = {sig: float(row[f"{sig}_mean_score"])
                              for sig in replicate_groups}
                except (KeyError, TypeError, ValueError):
                    continue
                if all(math.isfinite(v) for v in scores.values()):
                    matrix[motif_id] = scores
        design = [{"sample": sig, "group": "A" if grp == treated else "B"}
                  for sig, grp in replicate_groups.items()]
        r_rows = _run_r(rscript, "welch", matrix, design, workdir, "live_tobias")
        aria_rows, tf_names = {}, {}
        with open(results_table, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                aria_rows[str(row["motif_id"])] = {
                    "change": float(row["change"]), "padj": float(row["padj"]),
                    "pvalue": float(row["pvalue"])}
                tf_names[str(row["motif_id"])] = row.get("tf")
        concordance = _compare(aria_rows, r_rows)
        metrics["r_concordance"] = concordance
        metrics["aria_rows"] = {m: {"tf": tf_names.get(m), **v}
                                for m, v in sorted(aria_rows.items())}
        metrics["recovery_note"] = (
            "Not gated: this suite certifies that ARIA's replicate statistics on "
            "live TOBIAS output match external R. On this small synthetic "
            "genome, unprotected background motifs can also differ because "
            "TOBIAS footprint scores are not depth-invariant; recovery is "
            "reported, not claimed.")
        metrics["significant_tfs"] = sorted(
            str(tf_names.get(m)) for m, v in aria_rows.items() if v["padj"] <= ALPHA)
        diff = concordance.get("max_abs_padj_diff")
        gates["r_max_abs_padj_diff"] = (diff is not None and diff
                                        < GATES["live_max_abs_padj_diff_max"])
        gates["r_significant_jaccard"] = concordance["significant_jaccard"] >= 1.0
        gates["r_shared_all_motifs"] = (concordance["n_motifs_shared"]
                                        == concordance["n_motifs_aria"] > 0)
    else:
        gates["r_max_abs_padj_diff"] = False
    return {"status": "ran" if gates["ran"] else "failed", "metrics": metrics,
            "gates": gates, "pass": all(gates.values())}


# ── receipt hygiene + entry point ────────────────────────────────────────────


def _absolute_strings(obj: Any) -> list[str]:
    found: list[str] = []
    if isinstance(obj, str):
        if obj.startswith(_ABS_PREFIXES) or (len(obj) > 2 and obj[1:3] == ":\\"):
            found.append(obj)
    elif isinstance(obj, dict):
        for value in obj.values():
            found.extend(_absolute_strings(value))
    elif isinstance(obj, (list, tuple)):
        for value in obj:
            found.extend(_absolute_strings(value))
    return found


def _scrub(obj: Any) -> Any:
    """Replace any absolute path that slipped into free text (e.g. a reason)."""
    if isinstance(obj, str):
        if obj.startswith(_ABS_PREFIXES) or (len(obj) > 2 and obj[1:3] == ":\\"):
            return Path(obj.replace("\\", "/")).name
        return obj
    if isinstance(obj, dict):
        return {k: _scrub(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_scrub(v) for v in obj]
    return obj


def _suite_error(exc: Exception) -> dict[str, Any]:
    return {"status": "error", "error": f"{type(exc).__name__}",
            "gates": {"ran": False}, "pass": False}


def run_benchmark(workdir: str | Path | None = None, *,
                  bench_env: str = "aria-bench-env",
                  tobias_env: str = "aria-tobias-env",
                  run_live_tobias: bool = True) -> dict[str, Any]:
    """Run all suites and return a path-clean JSON-serializable receipt."""
    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory(prefix="aria_adr078_")
        base = Path(tmp.name)
    else:
        base = Path(workdir)
        base.mkdir(parents=True, exist_ok=True)
    try:
        rscript = _env_tool(bench_env, "Rscript")
        r_version = _r_version(rscript) if rscript else None
        suites: dict[str, Any] = {}
        for name, fn in (("welch_independent", _suite_welch),
                         ("paired_donor", _suite_paired),
                         ("batch_adjusted_ols", _suite_batch)):
            if rscript is None or r_version is None:
                suites[name] = {"status": "not_run",
                                "reason": f"Rscript unavailable in {bench_env}",
                                "gates": {"ran": False}, "pass": False}
                continue
            try:
                suites[name] = {"status": "ran", **fn(rscript, base / "r")}
            except Exception as exc:  # recorded, never imputed
                suites[name] = _suite_error(exc)
        try:
            suites["null_calibration"] = {"status": "ran", **_suite_null()}
        except Exception as exc:
            suites["null_calibration"] = _suite_error(exc)
        if not run_live_tobias:
            suites["live_tobias"] = {"status": "not_run",
                                     "reason": "live TOBIAS suite disabled",
                                     "gates": {"ran": False}, "pass": False}
        elif rscript is None or r_version is None:
            suites["live_tobias"] = {"status": "not_run",
                                     "reason": f"Rscript unavailable in {bench_env}",
                                     "gates": {"ran": False}, "pass": False}
        else:
            try:
                suites["live_tobias"] = _suite_live_tobias(rscript, base, tobias_env)
            except Exception as exc:
                suites["live_tobias"] = _suite_error(exc)

        try:
            import aria

            aria_version = getattr(aria, "__version__", None)
        except Exception:
            aria_version = None
        receipt: dict[str, Any] = {
            "benchmark": BENCHMARK_ID,
            "adr": "ADR-078",
            "schema_version": 1,
            "description": ("Replicate-level differential TF footprint inference "
                            "(Welch / paired / batch-adjusted OLS + BH across "
                            "motifs) vs external R stats, null calibration, and a "
                            "live TOBIAS ATACorrect->ScoreBigwig->BINDetect run."),
            "spec": SPEC,
            "gates": GATES,
            "alpha": ALPHA,
            "environment": {"aria_version": aria_version, "r_version": r_version,
                            "r_env": bench_env, "tobias_env": tobias_env},
            "suites": suites,
        }
        receipt = _scrub(json.loads(json.dumps(receipt, default=str)))
        receipt["paths_clean"] = not _absolute_strings(receipt)
        receipt["all_gates_pass"] = bool(
            receipt["paths_clean"]
            and all(s.get("pass") for s in receipt["suites"].values()))
        return receipt
    finally:
        if tmp is not None:
            tmp.cleanup()
