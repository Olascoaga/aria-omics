"""
ARIA scRNA End-to-End Validation Harness (Fase 3)
--------------------------------------------------
Runs the scRNA pipeline scripts (S1–S8) in sequence on a real dataset using
EnvironmentManager directly. No LLM, no MessageBus, no Orchestrator — just
subprocess calls to validate that every step works on real data.

Stages:
  1. rna_qc.py              — QC + Scrublet (S4)
  2. rna_advise_resolution  — evaluate Leiden candidates, pick best silhouette
  3. rna_clustering.py      — final clustering (idempotent, S1+S2+S5)
  4. rna_celltypist.py      — database-backed annotation (S3)
  5. rna_de_per_cluster.py  — per-cluster DE
  6. rna_pathway_per_cluster.py — per-cluster ORA (S8)

Usage:
  conda activate aria-env
  python tests/test_scrna_e2e.py [--data DIR_OR_H5AD] [--organism "Homo sapiens"]
                                 [--tissue pbmc] [--workspace /tmp/aria_e2e]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from aria.utils.environment_manager import env_manager


GRN = "\033[92m"; RED = "\033[91m"; YLW = "\033[93m"
CYN = "\033[96m"; DIM = "\033[2m"; RST = "\033[0m"; BLD = "\033[1m"


def banner(t: str) -> None:
    print(f"\n{BLD}{CYN}━━ {t} {'━' * max(0, 56 - len(t))}{RST}")


def ok(msg: str, detail: str = "") -> None:
    d = f"  {DIM}{detail}{RST}" if detail else ""
    print(f"  {GRN}✓{RST} {msg}{d}")


def warn(msg: str, detail: str = "") -> None:
    d = f"  {DIM}{detail}{RST}" if detail else ""
    print(f"  {YLW}!{RST} {msg}{d}")


def fail(msg: str, detail: str = "") -> None:
    d = f"\n    {DIM}{detail}{RST}" if detail else ""
    print(f"  {RED}✗{RST} {msg}{d}")


def step(name: str, fn) -> dict:
    banner(name)
    t0 = time.perf_counter()
    result = fn()
    dt = time.perf_counter() - t0
    status = result.get("status", "?")
    if status == "success":
        ok(f"{name} OK", f"{dt:.1f}s")
    elif status == "skipped":
        warn(f"{name} skipped", result.get("reason", ""))
    else:
        fail(f"{name} FAILED",
             f"{result.get('error_type', '?')}: "
             f"{str(result.get('details', ''))[:300]}")
    return result


def _sample_id_from_path(path: Path) -> str:
    stem = path.stem
    for suffix in ("_raw_feature_bc_matrix",
                   "_filtered_feature_bc_matrix",
                   "_feature_bc_matrix",
                   "_matrix"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem or "sample"


def _expand_inputs(data_arg: str) -> list[Path]:
    """
    Resolve --data into a list of input paths.

    Accepts:
      - a single .h5ad / .h5 file path
      - a MEX directory (matrix.mtx + barcodes + genes)
      - a directory containing multiple .h5/.h5ad files
      - a comma-separated list of any of the above
    """
    if "," in data_arg:
        parts = [Path(p.strip()).resolve() for p in data_arg.split(",")
                 if p.strip()]
        return parts

    p = Path(data_arg).resolve()
    if p.is_file():
        return [p]

    if p.is_dir():
        # MEX dir? — Scanpy needs matrix.mtx(.gz) + barcodes + genes/features
        if any((p / n).exists() for n in ("matrix.mtx", "matrix.mtx.gz")):
            return [p]
        # Otherwise: collect .h5/.h5ad children (sorted for determinism).
        children = sorted(
            list(p.glob("*.h5")) + list(p.glob("*.h5ad")),
            key=lambda x: x.name,
        )
        if not children:
            raise SystemExit(f"No .h5 / .h5ad / MEX inputs found under {p}")
        return children

    raise SystemExit(f"--data path does not exist: {p}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/tmp/aria_e2e_pbmc3k/input",
                    help=("Single file (MEX dir, .h5ad, .h5), "
                          "comma-separated list, or a directory of .h5/.h5ad."))
    ap.add_argument("--organism", default="Homo sapiens")
    ap.add_argument("--tissue", default="pbmc",
                    help="CellTypist tissue hint (pbmc, brain, immune, ...)")
    ap.add_argument("--resolutions", default="0.2,0.5,0.8,1.2",
                    help="comma-separated list of leiden resolutions to try")
    ap.add_argument("--workspace", default="/tmp/aria_e2e_pbmc3k/out",
                    help="output dir for intermediate artifacts")
    # Pseudobulk mode — when --pseudobulk-condition is set, skip QC/cluster/DE
    # and run between-condition DE directly on a preprocessed h5ad input.
    ap.add_argument("--pseudobulk-condition", default=None,
                    help=("obs column with experimental condition (e.g. "
                          "'age_group'); enables pseudobulk DE mode"))
    ap.add_argument("--pseudobulk-replicate", default="orig.ident",
                    help="obs column with biological replicate ID")
    ap.add_argument("--pseudobulk-groupby",   default="subclass",
                    help="obs column for cell-type stratification")
    ap.add_argument("--pseudobulk-compare",   default=None,
                    help=("comma-separated pairs as test:ref (semicolons to "
                          "chain). e.g. '80-100:20-39;60-79:20-39'"))
    ap.add_argument("--pseudobulk-covariates", default="",
                    help="comma-separated obs columns for the design formula")
    # Trajectory mode — preprocessed h5ad with clusters/celltypes; run
    # PAGA + DPT and skip QC / integration / clustering / pathways.
    ap.add_argument("--trajectory-h5ad", default=None,
                    help=("path to preprocessed h5ad (already QC'd + "
                          "clustered) — when set, enters trajectory-only "
                          "mode"))
    ap.add_argument("--trajectory-groupby", default="cell_type",
                    help="obs column for PAGA grouping (default: cell_type)")
    ap.add_argument("--trajectory-root", default=None,
                    help="cell type label to use as DPT root (optional)")
    ap.add_argument("--emit-html", action="store_true",
                    help=("After stages succeed, render the NarrativeAgent "
                          "HTML report into <workspace>/report/"))
    args = ap.parse_args()

    is_trajectory = args.trajectory_h5ad is not None
    if is_trajectory:
        # In trajectory mode, --data is ignored; the input is the
        # preprocessed h5ad. Set inputs to that single file for
        # consistency in `report["data"]`.
        inputs = [Path(args.trajectory_h5ad).resolve()]
    else:
        inputs = _expand_inputs(args.data)
    workspace = Path(args.workspace).resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    is_multi  = len(inputs) >= 2 and not is_trajectory
    is_pseudobulk = args.pseudobulk_condition is not None

    report: dict = {"data":        [str(p) for p in inputs],
                    "organism":    args.organism,
                    "tissue_hint": args.tissue,
                    "workspace":   str(workspace),
                    "n_samples":   len(inputs),
                    "stages":      {}}

    print(f"{BLD}ARIA scRNA E2E — "
          f"{len(inputs)} sample(s){RST}")
    if is_multi:
        for p in inputs:
            print(f"    {DIM}- {_sample_id_from_path(p)}  ({p}){RST}")
    print(f"{DIM}organism={args.organism} | tissue={args.tissue}"
          f" | workspace={workspace}{RST}")

    # ── Trajectory mode — early dispatch for preprocessed h5ad ──────────
    if is_trajectory:
        return _run_trajectory_flow(args, inputs[0], workspace, report)

    # ── Pseudobulk mode — early dispatch for preprocessed h5ad ──────────
    if is_pseudobulk:
        if len(inputs) != 1 or inputs[0].suffix != ".h5ad":
            return _save_and_exit(
                {**report, "error": "pseudobulk needs a single .h5ad input"},
                workspace, 1)
        return _run_pseudobulk_flow(args, inputs[0], workspace, report)

    # ── 1. QC + Scrublet ─────────────────────────────────────────────────
    if not is_multi:
        # Single-sample fast-path — identical to v4.3.2 behaviour.
        qc = step("1. rna_qc.py (QC + Scrublet)",
                  lambda: env_manager.run_in_stack(
                      stack="rna",
                      script_path="aria/scripts/rna_qc.py",
                      params={
                          "data_path": str(inputs[0]),
                          "organism":  args.organism,
                          "biological_context": {"summary": "e2e validation"},
                      },
                  ))
        report["stages"]["qc"] = qc
        if qc.get("status") != "success":
            return _save_and_exit(report, workspace, 1)
        print(f"    {DIM}cells: {qc.get('n_cells_before')} → "
              f"{qc.get('n_cells_after')} "
              f"({qc.get('pct_removed', 0):.1f}% removed) | "
              f"MT≤{qc.get('mt_threshold_used', '?')}%{RST}")
        scrub = qc.get("scrublet", {})
        if scrub.get("ran"):
            print(f"    {DIM}scrublet: {scrub.get('n_doublets')} doublets "
                  f"@ thr={scrub.get('threshold_used'):.3f} "
                  f"(rate={scrub.get('doublet_rate'):.2%}){RST}")
        qc_h5ad   = qc["output_path"]
        batch_col = None
    else:
        # Multi-sample: per-sample QC, then rna_concat.
        manifest      = []
        per_sample    = []
        qc_workspace  = workspace / "per_sample_qc"
        qc_workspace.mkdir(parents=True, exist_ok=True)
        banner("1. rna_qc.py per sample")
        t0 = time.perf_counter()
        for p in inputs:
            sid = _sample_id_from_path(p)
            r = env_manager.run_in_stack(
                stack="rna",
                script_path="aria/scripts/rna_qc.py",
                params={
                    "data_path":          str(p),
                    "organism":           args.organism,
                    "biological_context": {"summary": "e2e validation"},
                    "sample_id":          sid,
                    "output_dir":         str(qc_workspace),
                },
            )
            if r.get("status") != "success":
                fail(f"QC failed on {sid}",
                     f"{r.get('error_type', '?')}: "
                     f"{str(r.get('details', ''))[:200]}")
                report["stages"]["qc"] = {"status": "error",
                                          "failed_sample": sid,
                                          "details": r.get("details", "")}
                return _save_and_exit(report, workspace, 1)
            per_sample.append({
                "sample_id":     sid,
                "n_cells_after": r.get("n_cells_after", 0),
                "pct_removed":   r.get("pct_removed", 0),
                "scrublet":      r.get("scrublet", {}),
            })
            manifest.append({"path": r["output_path"], "sample_id": sid})
            print(f"    {DIM}{sid:<20}  "
                  f"{r.get('n_cells_after', 0):>6} cells  "
                  f"({r.get('pct_removed', 0):.1f}% removed)  "
                  f"doublets={r.get('scrublet', {}).get('n_doublets', 0)}{RST}")
        ok(f"per-sample QC: {len(manifest)} samples",
           f"{time.perf_counter()-t0:.1f}s")
        report["stages"]["qc"] = {"status": "success",
                                  "per_sample": per_sample,
                                  "n_samples":  len(manifest)}

        concat = step("1b. rna_concat.py",
                      lambda: env_manager.run_in_stack(
                          stack="rna",
                          script_path="aria/scripts/rna_concat.py",
                          params={"samples":    manifest,
                                  "output_dir": str(workspace),
                                  "join":       "inner"},
                      ))
        report["stages"]["concat"] = concat
        if concat.get("status") != "success":
            return _save_and_exit(report, workspace, 1)
        print(f"    {DIM}{concat.get('n_cells_total')} cells × "
              f"{concat.get('n_genes_shared')} shared genes "
              f"across {concat.get('n_samples')} samples{RST}")
        qc_h5ad   = concat["output_path"]
        batch_col = concat.get("batch_col", "batch")

    # ── 1c. Integration (Harmony) — only when multi-batch ────────────────
    if batch_col:
        integ = step("1c. rna_integration.py (Harmony)",
                     lambda: env_manager.run_in_stack(
                         stack="rna",
                         script_path="aria/scripts/rna_integration.py",
                         params={"data_path": qc_h5ad,
                                 "batch_col": batch_col,
                                 "output_dir": str(workspace)},
                     ))
        report["stages"]["integration"] = integ
        if integ.get("status") == "success":
            print(f"    {DIM}n_batches={integ.get('n_batches')}, "
                  f"silhouette {integ.get('silhouette_before'):.3f} → "
                  f"{integ.get('silhouette_after'):.3f} "
                  f"(Δ={integ.get('batch_correction_delta'):+.3f}, "
                  f"lower=better){RST}")
            qc_h5ad = integ["output_path"]
        elif integ.get("status") == "skipped":
            print(f"    {DIM}skipped: {integ.get('reason', '?')}{RST}")
        else:
            return _save_and_exit(report, workspace, 1)

    # ── 2. Advise Leiden resolution ──────────────────────────────────────
    resolutions = [float(x.strip()) for x in args.resolutions.split(",") if x.strip()]
    advice = step("2. rna_advise_resolution.py",
                  lambda: env_manager.run_in_stack(
                      stack="rna",
                      script_path="aria/scripts/rna_advise_resolution.py",
                      params={"data_path": qc_h5ad,
                              "resolutions": resolutions},
                  ))
    report["stages"]["advise"] = advice
    if advice.get("status") != "success":
        return _save_and_exit(report, workspace, 1)
    candidates = advice.get("candidates", [])
    for c in candidates:
        print(f"    {DIM}res={c['resolution']:.2f}: "
              f"{c['n_clusters']:2d} clusters, "
              f"silhouette={c['silhouette']:.3f}, "
              f"min={c['min_cluster_size']}, "
              f"singletons={c['n_singleton_clusters']}{RST}")
    # Pick: highest silhouette but require ≥3 clusters AND ≤1 singleton.
    eligible = [c for c in candidates
                if c["n_clusters"] >= 3 and c["n_singleton_clusters"] <= 1]
    pool = eligible or candidates
    best = max(pool, key=lambda c: c["silhouette"])
    chosen_res = best["resolution"]
    ok(f"Chose resolution={chosen_res} (silhouette={best['silhouette']:.3f}, "
       f"{best['n_clusters']} clusters)")

    # ── 3. Clustering ────────────────────────────────────────────────────
    clust = step("3. rna_clustering.py", lambda: env_manager.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_clustering.py",
        params={"data_path": qc_h5ad, "resolution": chosen_res},
    ))
    report["stages"]["clustering"] = clust
    if clust.get("status") != "success":
        return _save_and_exit(report, workspace, 1)
    print(f"    {DIM}{clust.get('n_clusters')} clusters, "
          f"rep={clust.get('rep_used', '?')}{RST}")
    clust_h5ad = clust["output_path"]
    top_markers = clust.get("top_markers", {})

    # ── 4. CellTypist annotation ─────────────────────────────────────────
    cty = step("4. rna_celltypist.py", lambda: env_manager.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_celltypist.py",
        params={
            "data_path":   clust_h5ad,
            "organism":    args.organism,
            "tissue_hint": args.tissue,
            "cluster_col": "leiden",
            "majority_voting": True,
        },
    ))
    report["stages"]["celltypist"] = cty
    annotated_h5ad = cty.get("output_path", clust_h5ad)
    if cty.get("status") == "success":
        per_cluster = cty.get("per_cluster", {})
        print(f"    {DIM}model: {cty.get('model_used', '?')}{RST}")
        for cl_id, info in list(per_cluster.items())[:8]:
            print(f"    {DIM}  cluster {cl_id}: {info.get('label')} "
                  f"({info.get('frequency', 0)*100:.0f}%){RST}")
        if len(per_cluster) > 8:
            print(f"    {DIM}  … +{len(per_cluster)-8} more{RST}")

    # ── 5. DE per cluster ────────────────────────────────────────────────
    de = step("5. rna_de_per_cluster.py", lambda: env_manager.run_in_stack(
        stack="rna",
        script_path="aria/scripts/rna_de_per_cluster.py",
        params={
            "data_path": annotated_h5ad,
            "groupby":   "leiden",
            "padj_max":  0.05,
            "lfc_min":   0.5,
            "top_n":     50,
            "output_dir": str(workspace),
        },
    ))
    report["stages"]["de"] = de
    if de.get("status") != "success":
        return _save_and_exit(report, workspace, 1)
    de_by_cluster = de.get("de_genes_by_cluster") or {}
    n_sig_total = sum(len(v) for v in de_by_cluster.values())
    print(f"    {DIM}n_significant_genes total: "
          f"{de.get('n_significant_genes', n_sig_total)} "
          f"across {len(de_by_cluster)} clusters{RST}")

    # ── 6. Pathway ORA per cluster ───────────────────────────────────────
    if not de_by_cluster:
        warn("Skipping pathway ORA — no DE genes")
        report["stages"]["pathways"] = {"status": "skipped",
                                         "reason": "no_de_genes"}
    else:
        pw = step("6. rna_pathway_per_cluster.py",
                  lambda: env_manager.run_in_stack(
                      stack="rna",
                      script_path="aria/scripts/rna_pathway_per_cluster.py",
                      params={
                          "de_genes_by_cluster":   de_by_cluster,
                          "organism":              args.organism,
                          "top_genes_per_cluster": 200,
                          "padj_db_max":           0.05,
                          "output_dir":            str(workspace),
                      },
                  ))
        report["stages"]["pathways"] = pw
        if pw.get("status") == "success":
            per = pw.get("per_cluster", {})
            ranked = sorted(per.items(),
                            key=lambda kv: kv[1].get("n_significant", 0),
                            reverse=True)[:5]
            for cl, info in ranked:
                results = info.get("results") or {}
                top_terms = []
                for db_label, entries in results.items():
                    if entries:
                        top_terms.append(f"{db_label}: {entries[0]['term'][:45]}")
                top_str = " | ".join(top_terms[:3])
                print(f"    {DIM}  cluster {cl}: "
                      f"{info.get('n_significant', 0)} sig pathways "
                      f"({top_str}){RST}")

    return _save_and_exit(report, workspace, 0, args=args)


def _run_trajectory_flow(args, h5ad: Path, workspace: Path,
                          report: dict) -> int:
    """
    Preprocessed-input fast path: skip QC/integration/clustering and run
    PAGA + DPT (+ scVelo if spliced/unspliced layers exist) on the
    grouping column provided.
    """
    report["mode"] = "trajectory"
    report["trajectory_inputs"] = {
        "h5ad":    str(h5ad),
        "groupby": args.trajectory_groupby,
        "root":    args.trajectory_root,
    }
    print(f"{DIM}mode=trajectory | groupby={args.trajectory_groupby} | "
          f"root={args.trajectory_root or 'auto'}{RST}")

    traj = step("1. rna_trajectory.py (PAGA + DPT)",
                lambda: env_manager.run_in_stack(
                    stack="rna",
                    script_path="aria/scripts/rna_trajectory.py",
                    params={
                        "data_path":      str(h5ad),
                        "root_cell_type": args.trajectory_root,
                        "cell_type_col":  args.trajectory_groupby,
                        "output_dir":     str(workspace),
                    },
                ))
    report["stages"]["trajectory"] = traj
    if traj.get("status") != "success":
        return _save_and_exit(report, workspace, 1, args=args)

    paga = traj.get("paga", {}) or {}
    pt = traj.get("pseudotime", {}) or {}
    print(f"    {DIM}PAGA: {paga.get('n_connections', 0)} edges, "
          f"max={paga.get('max_connectivity', 0):.4f}, "
          f"{paga.get('n_strong', 0)} strong (>{paga.get('strong_threshold', 0.05)}){RST}")
    if pt.get("computed"):
        pt_by = pt.get("pseudotime_by_group", {}) or {}
        ordered = sorted(pt_by.items(), key=lambda kv: kv[1])
        print(f"    {DIM}DPT root={pt.get('root_used')}, order: "
              f"{' → '.join(g for g, _ in ordered)}{RST}")
    vel = traj.get("velocity", {}) or {}
    if vel.get("computed"):
        print(f"    {DIM}scVelo: {vel.get('method', '?')}{RST}")
    else:
        print(f"    {DIM}scVelo: skipped ({vel.get('reason', '?')[:60]}){RST}")

    return _save_and_exit(report, workspace, 0, args=args)


def _run_pseudobulk_flow(args, h5ad: Path, workspace: Path,
                        report: dict) -> int:
    """
    Preprocessed-input fast path: skip QC/integration/clustering and run
    between-condition DE on the cell types already in obs, followed by
    pathway ORA per (group, comparison).
    """
    # Parse comparisons (e.g. "80-100:20-39;60-79:20-39")
    comparisons = []
    for pair in args.pseudobulk_compare.split(";") if args.pseudobulk_compare else []:
        pair = pair.strip()
        if not pair:
            continue
        if ":" not in pair:
            return _save_and_exit(
                {**report,
                 "error": f"--pseudobulk-compare expects test:ref, got '{pair}'"},
                workspace, 1)
        test, ref = pair.split(":", 1)
        comparisons.append([test.strip(), ref.strip()])
    if not comparisons:
        return _save_and_exit(
            {**report, "error": "--pseudobulk-compare is required in pseudobulk mode"},
            workspace, 1)

    covariates = [c.strip() for c in args.pseudobulk_covariates.split(",") if c.strip()]
    report["mode"]              = "pseudobulk"
    report["pseudobulk_inputs"] = {
        "h5ad":          str(h5ad),
        "condition":     args.pseudobulk_condition,
        "replicate":     args.pseudobulk_replicate,
        "groupby":       args.pseudobulk_groupby,
        "comparisons":   comparisons,
        "covariates":    covariates,
    }
    print(f"{DIM}mode=pseudobulk | groupby={args.pseudobulk_groupby} | "
          f"condition={args.pseudobulk_condition} | "
          f"replicate={args.pseudobulk_replicate}{RST}")
    print(f"{DIM}comparisons: "
          f"{', '.join(f'{a}_vs_{b}' for a,b in comparisons)}{RST}")
    if covariates:
        print(f"{DIM}covariates: {', '.join(covariates)}{RST}")

    pb = step("1. rna_pseudobulk_de.py",
              lambda: env_manager.run_in_stack(
                  stack="rna",
                  script_path="aria/scripts/rna_pseudobulk_de.py",
                  params={
                      "data_path":     str(h5ad),
                      "groupby":       args.pseudobulk_groupby,
                      "condition_col": args.pseudobulk_condition,
                      "replicate_col": args.pseudobulk_replicate,
                      "comparisons":   comparisons,
                      "covariates":    covariates,
                      "min_cells_per_pseudosample":   10,
                      "min_replicates_per_condition": 2,
                      "padj_max":      0.05,
                      "lfc_min":       0.5,
                      "top_n":         50,
                      "output_dir":    str(workspace),
                  },
              ))
    report["stages"]["pseudobulk"] = pb
    if pb.get("status") != "success":
        return _save_and_exit(report, workspace, 1)

    per_group = pb.get("per_group", {})
    print()
    print(f"    {DIM}{'group':<22}  {'pseudo':>6}  {'comp':>16}  "
          f"{'sig':>4}  top up gene{RST}")
    de_for_ora: dict = {}   # key: f"{group}_{comp}"  →  list of gene records
    for group in sorted(per_group):
        info = per_group[group]
        if info.get("status") == "skipped":
            print(f"    {DIM}{group:<22}  SKIP  {info.get('reason','')[:40]}{RST}")
            continue
        n_ps = info.get("n_pseudosamples", 0)
        for comp_key, comp in info.get("per_comparison", {}).items():
            if comp.get("status") != "success":
                print(f"    {DIM}{group:<22}  {n_ps:>6}  {comp_key:>16}  "
                      f"SKIP  {comp.get('reason', comp.get('error_type',''))[:30]}{RST}")
                continue
            n_sig = comp["n_significant"]
            top = (comp.get("top_genes") or [{}])[0]
            top_gene = top.get("gene", "-")
            print(f"    {DIM}{group:<22}  {n_ps:>6}  {comp_key:>16}  "
                  f"{n_sig:>4}  {top_gene} (LFC={top.get('log2fc',0):+.2f}){RST}")
            if comp.get("all_sig"):
                # Convert pseudobulk DE records to the format that
                # rna_pathway_per_cluster expects: {gene, log2fc, padj}.
                de_for_ora[f"{group}::{comp_key}"] = [
                    {"gene":   r["gene"],
                     "log2fc": r["log2fc"],
                     "padj":   r["padj"]}
                    for r in comp["all_sig"]
                ]

    if not de_for_ora:
        warn("No significant pseudobulk hits to pathway-enrich.")
        return _save_and_exit(report, workspace, 0, args=args)

    pw = step("2. rna_pathway_per_cluster.py (per group×comparison)",
              lambda: env_manager.run_in_stack(
                  stack="rna",
                  script_path="aria/scripts/rna_pathway_per_cluster.py",
                  params={
                      "de_genes_by_cluster":   de_for_ora,
                      "organism":              args.organism,
                      "top_genes_per_cluster": 200,
                      "padj_db_max":           0.05,
                      "output_dir":            str(workspace),
                  },
              ))
    report["stages"]["pathways"] = pw
    if pw.get("status") == "success":
        per = pw.get("per_cluster", {})
        ranked = sorted(per.items(),
                        key=lambda kv: kv[1].get("n_significant", 0),
                        reverse=True)[:8]
        for cl, info in ranked:
            results = info.get("results") or {}
            top_terms = []
            for db_label, entries in results.items():
                if entries:
                    top_terms.append(f"{db_label}: {entries[0]['term'][:40]}")
            top_str = " | ".join(top_terms[:3])
            print(f"    {DIM}  {cl}: "
                  f"{info.get('n_significant', 0)} sig pathways "
                  f"({top_str}){RST}")

    return _save_and_exit(report, workspace, 0, args=args)


def _emit_narrative_report(report: dict, workspace: Path, args) -> None:
    """
    Render the NarrativeAgent HTML report from the freshly-completed
    e2e_report. Stub the LLM + bus so the harness stays subprocess-only.
    """
    from aria.scripts.rna_narrative_adapter import adapt
    from aria.agents import _narrative_scrna
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.llm.provider import TaskTier
    from aria.utils.environment_manager import env_manager as _envm

    banner("7. NarrativeAgent — HTML report")
    report_dir = workspace / "report"
    figures_dir = report_dir / "figures"
    report_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    bundle = adapt(report, workspace=workspace)
    scrna_envelope = bundle["agent_results"]["scrna_agent"]
    findings = scrna_envelope["findings"]["scRNA"]["findings"]
    h5ad_path = scrna_envelope.get("output_h5ad")
    print(f"    {DIM}h5ad for figures: {h5ad_path}{RST}")

    t0 = time.perf_counter()
    _narrative_scrna.generate_figures(
        findings,
        h5ad_path=h5ad_path,
        output_dir=figures_dir,
        env_manager=_envm,
    )
    figs = findings.get("figures") or {}
    print(f"    {DIM}figures: "
          f"{sum(1 for k in figs if k.startswith('umap_'))} UMAPs, "
          f"{1 if 'per_celltype_de_bar' in figs else 0} DE bar, "
          f"{len(figs.get('pathway_dotplots') or {})} pathway blocks "
          f"({time.perf_counter()-t0:.1f}s){RST}")

    class _StubLLM:
        def complete(self, prompt, system="", tier=TaskTier.MEDIUM,
                     max_tokens=400, **kw):
            raise RuntimeError("LLM stubbed in harness")

    class _StubMemory:
        def get_decisions(self, experiment_id):
            return []
        def __getattr__(self, name):
            return lambda *a, **kw: None

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.name        = "narrative_agent"
    agent.memory      = _StubMemory()
    agent.llm         = _StubLLM()
    agent.reports_dir = report_dir

    from aria.bus.message_bus import bus as _bus
    _bus.register(agent.name, agent)

    exp_id = f"e2e_{int(time.time())}"
    result = agent.run(exp_id, {
        "exp_context":       bundle["exp_context"],
        "biological_intent": bundle["intent"],
        "agent_results":     bundle["agent_results"],
        "findings":          bundle["findings_list"],
    })

    ok(f"HTML report: {result.get('report_path')}")
    report["narrative"] = {
        "status":      result.get("status"),
        "report_path": result.get("report_path"),
        "figures":     figs,
    }


def _save_and_exit(report: dict, workspace: Path, code: int,
                   args=None) -> int:
    out = workspace / "e2e_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2, default=str)
    banner("Summary")
    stages = report.get("stages", {})
    for name, r in stages.items():
        st = r.get("status", "?")
        color = GRN if st == "success" else (YLW if st == "skipped" else RED)
        print(f"  {color}{st:>9}{RST}  {name}")
    print(f"\n{DIM}Report: {out}{RST}")
    if code == 0:
        print(f"{GRN}{BLD}✓ E2E validation PASSED{RST}\n")
        if args is not None and getattr(args, "emit_html", False):
            try:
                _emit_narrative_report(report, workspace, args)
                # Re-save now that report["narrative"] is populated
                with open(out, "w") as f:
                    json.dump(report, f, indent=2, default=str)
            except Exception as e:
                import traceback
                fail("NarrativeAgent emission failed",
                     f"{type(e).__name__}: {e}\n{traceback.format_exc()}")
                return 2
    else:
        print(f"{RED}{BLD}✗ E2E validation FAILED — first failure above{RST}\n")
    return code


if __name__ == "__main__":
    sys.exit(main())
