#!/usr/bin/env python3
"""ADR-077: bulk ATAC multi-replicate R concordance + permutation-null benchmark.

Orchestrates (via ``conda run``) the deterministic benchmark specified in
``aria/benchmarks/bulk_atac_concordance.py``:

  1. (this env) simulate 16 libraries = 8 biological replicates x 2 technical,
     balanced 2-level batch, 600 peaks (500 null + 100 planted DA).
  2. (--rna-env) ARIA ``aria/scripts/chromatin_bulk_diffacc.py`` through its real
     IPC entry point (``python <script> in.json out.json``), ``~ batch +
     condition``, technical libraries collapsed by ARIA.
  3. (this env) verify ARIA's replicate matrix equals an independent
     technical-replicate sum (B1).
  4. (--bench-env) Bioconductor DESeq2 / edgeR glmQLFTest / limma-voom (+ the
     threshold-matched glmTreat / treat) on ARIA's replicate matrix.
  5. (--rna-env, worker mode) 20 batch-balanced label permutations over the 500
     null peaks through ``chromatin_bulk_diffacc``.
  6. score + write a path-scrubbed receipt JSON.

Usage:
  python scripts/run_bulk_atac_concordance.py \
      [--output docs/benchmark_results/bulk_atac/adr077_multireplicate_r_concordance.json]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

from aria.benchmarks import bulk_atac_concordance as bac  # noqa: E402

DEFAULT_OUTPUT = (ROOT / "docs" / "benchmark_results" / "bulk_atac"
                  / "adr077_multireplicate_r_concordance.json")
DIFFACC_SCRIPT = ROOT / "aria" / "scripts" / "chromatin_bulk_diffacc.py"


def _conda_run(env: str, cmd: list[str]) -> None:
    child_env = dict(os.environ)
    child_env["PYTHONPATH"] = (str(ROOT) + os.pathsep
                               + child_env.get("PYTHONPATH", "")).rstrip(os.pathsep)
    subprocess.run(["conda", "run", "--no-capture-output", "-n", env, *cmd],
                   check=True, cwd=str(ROOT), env=child_env)


def _diffacc_params(counts_tsv: Path, meta_tsv: Path, out_dir: Path) -> dict:
    return {
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts_tsv),
        "sample_metadata_path": str(meta_tsv),
        "condition_col": "condition",
        "replicate_col": "replicate",
        "covariates": ["batch"],
        "required_covariates": ["batch"],
        "comparisons": [{"test": "treated", "reference": "control"}],
        "padj_max": bac.THRESHOLDS["padj"],
        "lfc_min": bac.THRESHOLDS["lfc"],
        "min_replicates_per_condition": 3,
        "output_dir": str(out_dir),
        "allow_mock": False,
    }


# ── worker mode (runs INSIDE the rna env) ─────────────────────────────────────

def _worker_permutations(spec_path: Path) -> int:
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc
    import pydeseq2

    spec = json.loads(spec_path.read_text())
    work = Path(spec["work_dir"])
    lib_meta = pd.read_csv(spec["library_metadata"], sep="\t")
    results = []
    for i, labels in enumerate(spec["permutations"], start=1):
        perm_dir = work / f"perm_{i:02d}"
        perm_dir.mkdir(parents=True, exist_ok=True)
        meta = lib_meta.copy()
        meta["condition"] = meta["replicate"].map(labels)
        meta_tsv = perm_dir / "library_metadata.tsv"
        meta.to_csv(meta_tsv, sep="\t", index=False)
        t0 = time.time()
        res = chromatin_bulk_diffacc(_diffacc_params(
            Path(spec["null_counts"]), meta_tsv, perm_dir / "da"))
        row = {"permutation": i, "status": "error", "labels": labels,
               "seconds": round(time.time() - t0, 2)}
        comp = (res.get("comparisons") or [{}])[0]
        if res.get("status") == "success" and res.get("ran") \
                and comp.get("status") == "success":
            full = pd.read_csv(comp["full_results_csv"])
            row.update({
                "status": "success",
                "n_tested": int(len(full)),
                "n_sig": int(comp.get("n_sig", 0)),
                "fitted_design_formula": comp.get("fitted_design_formula"),
                "pvalues": [float(p) for p in full["pvalue"].tolist()],
            })
        else:
            row["reason"] = res.get("reason") or comp.get("reason")
        results.append(row)
        print(f"perm {i:02d}: {row['status']} n_sig={row.get('n_sig')}",
              flush=True)
    Path(spec["output_json"]).write_text(json.dumps({
        "pydeseq2": pydeseq2.__version__, "results": results}))
    return 0


# ── orchestrator ──────────────────────────────────────────────────────────────

def _git_head() -> dict:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, check=True)
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                               capture_output=True, text=True, check=True)
        return {"commit": head.stdout.strip(),
                "working_tree_dirty": bool(dirty.stdout.strip())}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "working_tree_dirty": None}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--work-dir", default=None,
                    help="scratch dir (default: a fresh temp dir)")
    ap.add_argument("--rna-env", default="aria-rna-env")
    ap.add_argument("--bench-env", default="aria-bench-env")
    ap.add_argument("--n-permutations", type=int, default=bac.N_PERMUTATIONS)
    ap.add_argument("--worker-permutations", default=None,
                    help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if args.worker_permutations:
        return _worker_permutations(Path(args.worker_permutations))

    t_start = time.time()
    work = Path(args.work_dir or tempfile.mkdtemp(prefix="aria_adr077_bulk_atac_"))
    work.mkdir(parents=True, exist_ok=True)
    print(f"work dir: {work}", flush=True)

    # 1) simulate
    lib_counts, lib_meta, truth = bac.simulate_multireplicate_bulk_atac()
    lib_counts_tsv = work / "library_peak_counts.tsv"
    lib_meta_tsv = work / "library_metadata.tsv"
    bac.write_counts_tsv(lib_counts, lib_counts_tsv)
    lib_meta.to_csv(lib_meta_tsv, sep="\t", index=False)
    truth.to_csv(work / "truth.tsv", sep="\t", index=False)

    # 2) ARIA via the real script IPC entry point (rna env)
    aria_in, aria_out = work / "aria_in.json", work / "aria_out.json"
    aria_in.write_text(json.dumps(_diffacc_params(
        lib_counts_tsv, lib_meta_tsv, work / "aria_da")))
    _conda_run(args.rna_env, ["python", str(DIFFACC_SCRIPT),
                              str(aria_in), str(aria_out)])
    aria = json.loads(aria_out.read_text())
    comp = (aria.get("comparisons") or [{}])[0]
    if aria.get("status") != "success" or not aria.get("ran") \
            or comp.get("status") != "success":
        print("ARIA bulk ATAC DA did not run:", aria.get("reason"),
              comp.get("reason"))
        return 1
    aria_full = pd.read_csv(comp["full_results_csv"]).set_index("peak")

    # 3) B1 technical aggregation check
    aria_rep = pd.read_csv(aria["replicate_counts_path"], sep="\t") \
        .set_index("peak_id")
    expected = bac.collapse_technical(lib_counts, lib_meta)
    expected = expected.loc[aria_rep.index, sorted(expected.columns)]
    aggregation_exact = bool(
        list(sorted(aria_rep.columns)) == list(expected.columns)
        and (aria_rep[expected.columns].to_numpy()
             == expected.to_numpy()).all())

    # 4) R comparators on ARIA's replicate-level matrix
    r_counts_tsv = work / "replicate_counts_for_r.tsv"
    aria_rep.reset_index().to_csv(r_counts_tsv, sep="\t", index=False)
    bio_meta = lib_meta.drop_duplicates("replicate").copy()
    bio_meta["sample"] = bio_meta["condition"] + "::" + bio_meta["replicate"]
    r_meta_tsv = work / "replicate_metadata_for_r.tsv"
    bio_meta[["sample", "condition", "replicate", "batch"]].to_csv(
        r_meta_tsv, sep="\t", index=False)
    r_script = work / "adr077_bulk_atac_r_comparators.R"
    r_script.write_text(bac.R_COMPARATORS)
    r_out = work / "r_out"
    _conda_run(args.bench_env, [
        "Rscript", str(r_script), str(r_counts_tsv), str(r_meta_tsv),
        str(r_out), str(bac.THRESHOLDS["padj"]), str(bac.THRESHOLDS["lfc"])])
    rj = json.loads((r_out / "r_comparators.json").read_text())

    aria_vs_r: dict[str, dict] = {}
    truth_scores: dict[str, dict] = {}
    aria_sig = bac.aria_significant(aria_full)
    aria_lfc = aria_full.get("log2FoldChange_raw", aria_full["log2FoldChange"])
    truth_scores["aria_chromatin_bulk_diffacc"] = bac.score_against_truth(
        aria_sig, aria_lfc, truth)
    for name, meta in (rj.get("methods") or {}).items():
        if meta.get("status") != "success":
            aria_vs_r[name] = {"status": meta.get("status"),
                               "error_type": meta.get("error_type"),
                               "details": meta.get("details")}
            continue
        ref = pd.read_csv(r_out / meta["results_tsv"], sep="\t") \
            .set_index("peak").dropna(subset=["padj"])
        matched = bac.R_METHOD_THRESHOLD_MATCHED.get(name, False)
        aria_vs_r[name] = {"status": "success",
                           "threshold_matched_hypothesis": matched,
                           **bac.score_against_reference(aria_full, ref, matched)}
        truth_scores[f"r_{name}"] = bac.score_against_truth(
            bac.reference_significant(ref, matched),
            ref["log2FoldChange"], truth)
        s = aria_vs_r[name]
        print(f"{name}: lfc_spearman={s['lfc_spearman']} "
              f"sign={s['sign_agreement_all_peaks']} "
              f"jaccard={s['significant_overlap']['jaccard']}", flush=True)

    # 5) permutation null over the 500 null peaks (worker in rna env)
    null_ids = truth.loc[~truth["is_differential"], "peak_id"]
    null_counts_tsv = work / "null_library_peak_counts.tsv"
    bac.write_counts_tsv(lib_counts.loc[null_ids], null_counts_tsv)
    perms = bac.batch_balanced_permutations(lib_meta, args.n_permutations)
    spec_path = work / "perm_spec.json"
    perm_out = work / "perm_results.json"
    spec_path.write_text(json.dumps({
        "work_dir": str(work / "permutations"),
        "null_counts": str(null_counts_tsv),
        "library_metadata": str(lib_meta_tsv),
        "permutations": perms,
        "output_json": str(perm_out),
    }))
    _conda_run(args.rna_env, ["python", str(Path(__file__).resolve()),
                              "--worker-permutations", str(spec_path)])
    perm = json.loads(perm_out.read_text())
    per_perm = perm["results"]
    null_summary = bac.summarise_permutations(per_perm, int(len(null_ids)))
    bio_lookup = lib_meta.drop_duplicates("replicate")
    permutation_rows = []
    for row in per_perm:
        test_set = frozenset(r for r, lab in row["labels"].items()
                             if lab == "treated")
        permutation_rows.append({
            "permutation": row["permutation"],
            "status": row["status"],
            "permuted_treated_replicates": sorted(test_set),
            "n_true_treated_in_permuted_treated":
                bac.overlap_with_truth(test_set, bio_lookup),
            "n_tested": row.get("n_tested"),
            "n_false_discoveries": row.get("n_sig"),
            "fitted_design_formula": row.get("fitted_design_formula"),
            "reason": row.get("reason"),
        })

    gates = bac.evaluate_gates(aria_vs_r,
                               truth_scores["aria_chromatin_bulk_diffacc"],
                               null_summary, aggregation_exact)

    from aria.version import __version__
    receipt = {
        "benchmark": "ADR-077_bulk_ATAC_multireplicate_external_R_concordance",
        "adr": "ADR-077",
        "scope": ("bulk ATAC replicate-level differential accessibility "
                  "(aria/scripts/chromatin_bulk_diffacc.py): technical-replicate "
                  "aggregation, ~ batch + condition fit, concordance with "
                  "Bioconductor DESeq2 / edgeR / limma-voom, recovery of planted "
                  "differential peaks, and label-permutation null control"),
        "aria_version": __version__,
        "git": _git_head(),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - t_start, 1),
        "reproduce": ("python scripts/run_bulk_atac_concordance.py "
                      "(requires conda envs aria-rna-env + aria-bench-env)"),
        "data": {
            "kind": "deterministic simulation (no external data)",
            "generator": "aria.benchmarks.bulk_atac_concordance."
                         "simulate_multireplicate_bulk_atac",
            "spec": bac.SIM_SPEC,
            "n_libraries": int(lib_counts.shape[1]),
            "n_biological_replicates": int(lib_meta["replicate"].nunique()),
            "n_peaks": int(lib_counts.shape[0]),
            "library_counts_sha256": bac.matrix_sha256(lib_counts),
            "design_table": bio_meta[["sample", "condition", "replicate",
                                      "batch"]].to_dict(orient="records"),
        },
        "aria": {
            "script": "aria/scripts/chromatin_bulk_diffacc.py",
            "environment_stack": "rna",
            "invocation": "script IPC (python <script> in.json out.json)",
            "validation_level": aria.get("validation_level"),
            "method": aria.get("method"),
            "fitted_design_formula": comp.get("fitted_design_formula"),
            "covariates_adjusted": aria.get("covariates_adjusted"),
            "lfc_threshold_test": comp.get("lfc_threshold_test"),
            "lfc_shrinkage": comp.get("lfc_shrinkage"),
            "n_replicate_samples": aria.get("n_replicate_samples"),
            "n_peaks_tested": aria.get("n_peaks_tested"),
            "n_sig": comp.get("n_sig"), "n_up": comp.get("n_up"),
            "n_down": comp.get("n_down"),
            "n_nonconverged_excluded": comp.get("n_nonconverged_excluded"),
            "pydeseq2": perm.get("pydeseq2"),
            "thresholds": bac.THRESHOLDS,
        },
        "technical_replicate_aggregation": {
            "n_libraries_in": int(lib_counts.shape[1]),
            "n_biological_replicates_out": int(aria_rep.shape[1]),
            "matches_independent_sum_exactly": aggregation_exact,
        },
        "external_comparators": {
            "environment": "aria-bench-env (Bioconductor)",
            "design": rj.get("design"),
            "versions": rj.get("versions"),
            "input": "ARIA's replicate-level count matrix (identical input)",
            "aria_lfc_compared": ("raw MLE log2FC (log2FoldChange_raw); "
                                  "apeGLM-shrunk Spearman reported separately"),
            "aria_vs_method": aria_vs_r,
        },
        "planted_truth_recovery": truth_scores,
        "permutation_null": {
            "design": ("batch-balanced relabelling of the 8 biological "
                       "replicates (2 'treated' per batch), true labelling "
                       "and its complement excluded, technical libraries move "
                       "with their biological replicate; 500 null peaks only"),
            "summary": null_summary,
            "per_permutation": permutation_rows,
        },
        "gates": gates,
        "caveats": [
            "Simulated data: the generative model (Gamma-Poisson biological "
            "variation, Poisson technical sampling, per-peak batch shift) is "
            "stated in full; real-data concordance is covered by the K562 vs "
            "GM12878 artifacts in this directory.",
            "edgeR glmQLFTest and limma-voom test the point null (log2FC = 0) "
            "and are called at FDR plus |log2FC|; ARIA and R DESeq2 test "
            "|log2FC| > lfc_min. glmTreat / treat are the threshold-matched "
            "edgeR / limma variants.",
            "Under the global null every discovery is false, so the fraction "
            "of permutations with any discovery is the empirical FDR.",
        ],
    }
    receipt = bac.scrub_paths(receipt, ROOT)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2, sort_keys=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}")
    print("gates:", json.dumps(gates["checks"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
