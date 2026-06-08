"""Benchmark A1 external R comparator runner.

This IPC script runs inside ``aria-bench-env`` via
``EnvironmentManager.run_in_stack(stack="benchmark", ...)``. It exports the same
neutral synthetic bulk matrix used by the preliminary A1 lane, calls the R
reference comparators (DESeq2, edgeR-QLF, limma-voom), and scores their result
tables against the frozen synthetic truth.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from aria.benchmarks.synthetic_de import (
    SyntheticBulkDEDataset,
    _rank_biased_overlap,
    simulate_bulk_dataset,
)
from aria.scripts._base import run_script


def _sim_kwargs(quick: bool) -> dict[str, int]:
    return (
        {"n_genes": 600, "n_de": 60, "replicates_per_condition": 6}
        if quick else
        {"n_genes": 1000, "n_de": 120, "replicates_per_condition": 6}
    )


def _write_inputs(dataset: SyntheticBulkDEDataset, out_dir: Path) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    counts_path = out_dir / "a1_counts.tsv"
    metadata_path = out_dir / "a1_metadata.tsv"
    truth_path = out_dir / "a1_truth.tsv"

    counts = dataset.counts.copy()
    counts.insert(0, "gene", counts.index)
    counts.to_csv(counts_path, sep="\t", index=False)

    metadata = dataset.metadata.copy()
    metadata.insert(0, "sample", metadata.index)
    metadata.to_csv(metadata_path, sep="\t", index=False)

    import pandas as pd
    truth = pd.DataFrame({
        "gene": list(dataset.true_log2fc),
        "true_log2fc": [dataset.true_log2fc[g] for g in dataset.true_log2fc],
        "is_true_de": [g in dataset.de_genes for g in dataset.true_log2fc],
        "direction": [
            dataset.de_genes.get(g, "null") for g in dataset.true_log2fc
        ],
    })
    truth.to_csv(truth_path, sep="\t", index=False)
    return {
        "counts_tsv": str(counts_path),
        "metadata_tsv": str(metadata_path),
        "truth_tsv": str(truth_path),
    }


def _finite(value: Any, default: float = 0.0) -> float:
    try:
        val = float(value)
        return val if math.isfinite(val) else default
    except Exception:
        return default


def _score_table(
    dataset: SyntheticBulkDEDataset,
    table_path: str,
    *,
    padj_max: float,
    lfc_min: float,
    top_k: int,
) -> dict[str, Any]:
    import pandas as pd

    df = pd.read_csv(table_path, sep="\t")
    if "gene" not in df.columns:
        raise ValueError(f"Comparator table lacks 'gene': {table_path}")
    df["gene"] = df["gene"].astype(str)
    df = df.set_index("gene", drop=False)
    for col in ("log2FoldChange", "pvalue", "padj"):
        if col not in df.columns:
            df[col] = float("nan")
        df[col] = pd.to_numeric(df[col], errors="coerce")

    true_de = set(dataset.de_genes)
    null = set(dataset.null_genes)
    called = set(
        df.loc[
            (df["padj"] < padj_max) & (df["log2FoldChange"].abs() >= lfc_min),
            "gene",
        ].astype(str)
    )
    tp = called & true_de
    fp = called & null

    truth = pd.Series(dataset.true_log2fc, name="true_log2fc", dtype=float)
    lfc_frame = pd.DataFrame({
        "estimated": df["log2FoldChange"],
        "true": truth,
    }).dropna()
    lfc_de = lfc_frame[lfc_frame["true"].abs() > 0]

    ranked_est = (
        df.assign(
            _padj=df["padj"].fillna(1.0),
            _abs_lfc=df["log2FoldChange"].abs().fillna(0.0),
        )
        .sort_values(["_padj", "_abs_lfc"], ascending=[True, False])
        .index.astype(str)
        .tolist()
    )
    ranked_truth = (
        pd.Series(dataset.true_log2fc, dtype=float)
        .abs()
        .sort_values(ascending=False)
        .index.astype(str)
        .tolist()
    )
    est_top = set(ranked_est[:top_k])
    truth_top = set(ranked_truth[:top_k])
    top_intersection = est_top & truth_top

    recall = len(tp) / max(len(true_de), 1)
    empirical_fdr = len(fp) / max(len(called), 1)
    precision = len(tp) / max(len(called), 1)
    spearman_de = (
        _finite(lfc_de["estimated"].corr(lfc_de["true"], method="spearman"))
        if not lfc_de.empty else 0.0
    )
    pearson_de = (
        _finite(lfc_de["estimated"].corr(lfc_de["true"], method="pearson"))
        if not lfc_de.empty else 0.0
    )
    top_k_jaccard = len(top_intersection) / max(len(est_top | truth_top), 1)
    return {
        "lfc_concordance": {
            "pearson_true_de": round(pearson_de, 4),
            "spearman_true_de": round(spearman_de, 4),
            "n_true_de_scored": int(len(lfc_de)),
        },
        "ranking_concordance": {
            "top_k": int(top_k),
            "top_k_jaccard": round(top_k_jaccard, 4),
            "top_k_truth_recall": round(len(top_intersection) / max(len(truth_top), 1), 4),
            "rank_biased_overlap_p0.9": round(
                _rank_biased_overlap(ranked_est[:top_k], ranked_truth[:top_k], p=0.9),
                4,
            ),
        },
        "significant_call_concordance": {
            "recall": round(recall, 4),
            "precision": round(precision, 4),
            "empirical_fdr": round(empirical_fdr, 4),
            "n_true_de": len(true_de),
            "n_called": len(called),
            "n_true_positive": len(tp),
            "n_false_positive": len(fp),
        },
    }


def _score_methods(
    dataset: SyntheticBulkDEDataset,
    r_payload: dict[str, Any],
    *,
    padj_max: float,
    lfc_min: float,
    top_k: int,
) -> dict[str, Any]:
    methods = {}
    for key, rec in (r_payload.get("methods") or {}).items():
        method = dict(rec or {})
        if method.get("status") != "success":
            methods[key] = method
            continue
        try:
            method["scores"] = _score_table(
                dataset,
                str(method["results_tsv"]),
                padj_max=padj_max,
                lfc_min=lfc_min,
                top_k=top_k,
            )
        except Exception as exc:
            method = {
                "status": "error",
                "method": method.get("method", key),
                "error_type": type(exc).__name__,
                "details": str(exc),
                "results_tsv": method.get("results_tsv"),
            }
        methods[key] = method
    return methods


def benchmark_a1_external_comparators(params: dict[str, Any]) -> dict[str, Any]:
    seed = int(params.get("seed", 11))
    quick = bool(params.get("quick", False))
    padj_max = float(params.get("padj_max", 0.05))
    lfc_min = float(params.get("lfc_min", 0.5))
    top_k = int(params.get("top_k", 0) or 0)
    rscript_bin = str(params.get("rscript_bin", "Rscript"))

    output_dir = Path(params.get("output_dir") or tempfile.mkdtemp(prefix="aria_a1_external_"))
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = output_dir / "inputs"
    r_output_dir = output_dir / "r_outputs"
    r_output_dir.mkdir(parents=True, exist_ok=True)
    r_json = r_output_dir / "r_comparators.json"

    dataset = simulate_bulk_dataset(seed=seed, **_sim_kwargs(quick))
    if top_k <= 0:
        top_k = min(100, max(1, dataset.n_de))
    artifacts = _write_inputs(dataset, input_dir)

    resolved_rscript = shutil.which(rscript_bin)
    if not resolved_rscript:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"Rscript not found on PATH: {rscript_bin}",
            "benchmark": "A1_bulk_de_external_comparators",
            "artifacts": artifacts,
        }

    r_script = Path(__file__).with_suffix(".R")
    cmd = [
        resolved_rscript,
        str(r_script),
        "--counts", artifacts["counts_tsv"],
        "--metadata", artifacts["metadata_tsv"],
        "--output-dir", str(r_output_dir),
        "--output-json", str(r_json),
        "--padj", str(padj_max),
        "--lfc", str(lfc_min),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=int(params.get("timeout", 14400)))
    if proc.returncode != 0:
        return {
            "status": "error",
            "error_type": "RComparatorFailed",
            "details": proc.stderr[-2000:] or proc.stdout[-2000:],
            "command": cmd,
            "benchmark": "A1_bulk_de_external_comparators",
            "artifacts": artifacts,
        }
    if not r_json.exists():
        return {
            "status": "error",
            "error_type": "MissingOutput",
            "details": f"R comparator runner did not write {r_json}",
            "command": cmd,
            "benchmark": "A1_bulk_de_external_comparators",
            "artifacts": artifacts,
        }
    r_payload = json.loads(r_json.read_text(encoding="utf-8"))
    if r_payload.get("status") != "success":
        return {
            "status": "error",
            "error_type": r_payload.get("error_type", "RComparatorError"),
            "details": r_payload.get("details", ""),
            "command": cmd,
            "benchmark": "A1_bulk_de_external_comparators",
            "artifacts": artifacts | {"r_json": str(r_json)},
        }

    methods = _score_methods(
        dataset, r_payload, padj_max=padj_max, lfc_min=lfc_min, top_k=top_k
    )
    method_status = {
        key: rec.get("status") == "success" and "scores" in rec
        for key, rec in methods.items()
    }
    manifest = {
        "status": "success",
        "external_comparator_status": (
            "complete" if method_status and all(method_status.values()) else "error"
        ),
        "benchmark": "A1_bulk_de_external_comparators",
        "benchmark_version": "v1",
        "artifact_version": "v4.5.5",
        "scope": "external_r_comparator_synthetic_truth",
        "comparison": {"numerator": "COND_B", "denominator": "COND_A"},
        "dataset": dataset.params,
        "thresholds": {"padj": padj_max, "lfc_threshold": lfc_min},
        "methods": methods,
        "method_status": method_status,
        "r_versions": r_payload.get("versions", {}),
        "artifacts": artifacts | {"r_json": str(r_json)},
        "messages": [
            "A1 external comparator runner executed DESeq2, edgeR-QLF, and "
            "limma-voom in aria-bench-env through JSON IPC.",
            "This validates the comparator execution path on synthetic truth; "
            "SEQC/MAQC/ERCC external references remain a separate lane.",
        ],
    }
    manifest_path = output_dir / str(
        params.get("manifest_name", "a1_external_comparators_v4.5.5.json")
    )
    manifest["artifacts"]["manifest_json"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


if __name__ == "__main__":
    run_script(benchmark_a1_external_comparators)
