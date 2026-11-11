#!/usr/bin/env python3
"""V48 Slice I2 — Multi-env WNN external concordance runner (ARIA vs Seurat v5).

Orchestrates:
  1. (`--chromatin-env`, default `aria-chromatin-env`):
     `aria/scripts/benchmark_wnn_export.py` -> runs V48 I1 paired QC (optional)
     + two seeded executions of `aria/scripts/integration_wnn.py` + exports
     shared `X_pca` and depth-filtered `X_lsi` embeddings.
  2. (`--seurat-env`, default `aria-bench-signac-env`):
     `aria/scripts/benchmark_wnn_seurat.R` -> runs `Seurat::FindMultiModalNeighbors`
     (`l2.norm = TRUE`) + `FindClusters(graph.name = "wsnn")` on the IDENTICAL
     embeddings and cells.
  3. (caller env):
     `aria.benchmarks.concordance_wnn.score_wnn_concordance` -> cluster ARI/NMI,
     per-cell modality weight Spearman/Pearson/MAE, weighted-kNN Jaccard/recall,
     seed determinism, and optional neutral reference-label alignment, emitted
     with repo-relative paths only (guarded by
     `tests/test_public_artifacts_no_absolute_paths.py`).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.concordance_wnn import score_wnn_concordance  # noqa: E402


def _resolve_conda_bin() -> str:
    found = shutil.which("conda")
    if found:
        return found
    candidate = Path(sys.executable).resolve().parents[2] / "bin" / "conda"
    if candidate.is_file():
        return str(candidate)
    return "conda"


def _conda_run(env: str, cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [_resolve_conda_bin(), "run", "--no-capture-output", "-n", env, *cmd],
        check=False,
    )


def _sanitize_public_value(val: Any) -> Any:
    """Strip machine-absolute filesystem prefixes from public manifest values."""
    if isinstance(val, dict):
        return {str(k): _sanitize_public_value(v) for k, v in val.items()}
    if isinstance(val, list):
        return [_sanitize_public_value(v) for v in val]
    if isinstance(val, str):
        p = Path(val)
        if val.startswith(("/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/")):
            try:
                return str(p.resolve().relative_to(ROOT.resolve()))
            except Exception:
                return p.name
        return val
    return val


def _read_aria_cells(path: Path) -> dict[str, list[Any]]:
    cols: dict[str, list[Any]] = {
        "barcodes": [],
        "wnn1": [],
        "wnn2": [],
        "rna_only": [],
        "atac_only": [],
        "rna_w1": [],
        "rna_w2": [],
        "ref_label": [],
    }
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cols["barcodes"].append(row["barcode"])
            cols["wnn1"].append(row["leiden_wnn"])
            cols["wnn2"].append(row["leiden_wnn_run2"])
            cols["rna_only"].append(row["leiden_rna"])
            cols["atac_only"].append(row["leiden_atac"])
            cols["rna_w1"].append(float(row["rna_weight"]))
            cols["rna_w2"].append(float(row["rna_weight_run2"]))
            cols["ref_label"].append(row.get("reference_label", ""))
    return cols


def _read_seurat_cells(path: Path) -> dict[str, list[Any]]:
    cols: dict[str, list[Any]] = {
        "barcodes": [],
        "clusters": [],
        "rna_w": [],
    }
    if not path.is_file():
        return cols
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cols["barcodes"].append(row["barcode"])
            cols["clusters"].append(row["seurat_cluster"])
            cols["rna_w"].append(float(row["rna_weight"]))
    return cols


def _read_knn_map(path: Path) -> dict[str, list[str]]:
    if not path.is_file():
        return {}
    out: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            raw = row.get("neighbors", "") or ""
            out[row["barcode"]] = [x for x in raw.split(",") if x]
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5mu", required=True, help="Paired RNA+ATAC .h5mu")
    p.add_argument("--work-dir", default="/tmp/v48_i2_wnn_bench")
    p.add_argument("--output-dir", default="docs/benchmark_results/multiome_concordance")
    p.add_argument("--manifest-name", default="v48_i2_wnn_concordance.json")
    p.add_argument("--chromatin-env", default="aria-chromatin-env")
    p.add_argument("--seurat-env", default="aria-bench-signac-env")
    p.add_argument("--run-paired-qc", action="store_true")
    p.add_argument("--reference-h5ad", default=None)
    p.add_argument("--reference-label-col", default=None)
    p.add_argument("--reference-sample-col", default=None)
    p.add_argument("--sample-id", default=None)
    p.add_argument("--label-col", default=None)
    p.add_argument("--k", type=int, default=20)
    p.add_argument("--resolution", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--dataset", default=None)
    p.add_argument("--skip-export", action="store_true")
    p.add_argument("--skip-seurat", action="store_true")
    args = p.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    aria_summary_json = work / "aria_wnn_summary.json"
    aria_cells_tsv = work / "aria_wnn_cells.tsv"
    aria_knn_tsv = work / "aria_wnn_knn.tsv"
    seurat_cells_tsv = work / "seurat_wnn_cells.tsv"
    seurat_knn_tsv = work / "seurat_wnn_knn.tsv"
    seurat_json = work / "seurat_status.json"

    if not args.skip_export:
        export_cmd = [
            "python",
            str(ROOT / "aria" / "scripts" / "benchmark_wnn_export.py"),
            "--h5mu", str(args.h5mu),
            "--out-dir", str(work),
            "--k", str(args.k),
            "--resolution", str(args.resolution),
            "--seed", str(args.seed),
        ]
        if args.run_paired_qc:
            export_cmd.append("--run-paired-qc")
        if args.reference_h5ad:
            export_cmd.extend(["--reference-h5ad", str(args.reference_h5ad)])
        if args.reference_label_col:
            export_cmd.extend(["--reference-label-col", str(args.reference_label_col)])
        if args.reference_sample_col:
            export_cmd.extend(["--reference-sample-col", str(args.reference_sample_col)])
        if args.sample_id:
            export_cmd.extend(["--sample-id", str(args.sample_id)])
        if args.label_col:
            export_cmd.extend(["--label-col", str(args.label_col)])

        proc_exp = _conda_run(args.chromatin_env, export_cmd)
        if proc_exp.returncode != 0:
            raise RuntimeError(
                f"benchmark_wnn_export.py failed with exit code {proc_exp.returncode}"
            )

    if not args.skip_seurat:
        seurat_cmd = [
            "Rscript",
            str(ROOT / "aria" / "scripts" / "benchmark_wnn_seurat.R"),
            "--in-dir", str(work),
            "--out-cells", str(seurat_cells_tsv),
            "--out-knn", str(seurat_knn_tsv),
            "--out-json", str(seurat_json),
            "--k", str(args.k),
            "--resolution", str(args.resolution),
            "--seed", str(args.seed),
        ]
        _conda_run(args.seurat_env, seurat_cmd)

    aria_summary = (
        json.loads(aria_summary_json.read_text(encoding="utf-8"))
        if aria_summary_json.is_file()
        else {"status": "not_run", "reason": "missing aria_wnn_summary.json"}
    )
    seurat_status = (
        json.loads(seurat_json.read_text(encoding="utf-8"))
        if seurat_json.is_file()
        else {"status": "not_run", "reason": "missing seurat_status.json"}
    )

    if aria_summary.get("status") != "success" or not aria_cells_tsv.is_file():
        concordance = {
            "status": "not_run",
            "reason": "ARIA WNN export did not complete successfully.",
        }
    else:
        aria_cols = _read_aria_cells(aria_cells_tsv)
        aria_knn = _read_knn_map(aria_knn_tsv)
        has_ref_labels = any(bool(str(x).strip()) for x in aria_cols["ref_label"])

        if seurat_status.get("status") == "success" and seurat_cells_tsv.is_file():
            seurat_cols = _read_seurat_cells(seurat_cells_tsv)
            seurat_knn = _read_knn_map(seurat_knn_tsv)
            concordance = score_wnn_concordance(
                aria_barcodes=aria_cols["barcodes"],
                aria_clusters=aria_cols["wnn1"],
                aria_rna_weights=aria_cols["rna_w1"],
                ref_barcodes=seurat_cols["barcodes"],
                ref_clusters=seurat_cols["clusters"],
                ref_rna_weights=seurat_cols["rna_w"],
                aria_knn=aria_knn,
                ref_knn=seurat_knn,
                aria_rna_only_clusters=aria_cols["rna_only"],
                aria_atac_only_clusters=aria_cols["atac_only"],
                reference_labels=aria_cols["ref_label"] if has_ref_labels else None,
                seed_cluster_runs=[aria_cols["wnn1"], aria_cols["wnn2"]],
                seed_weight_runs=[aria_cols["rna_w1"], aria_cols["rna_w2"]],
            )
        else:
            concordance = score_wnn_concordance(
                aria_barcodes=aria_cols["barcodes"],
                aria_clusters=aria_cols["wnn1"],
                aria_rna_weights=aria_cols["rna_w1"],
                ref_barcodes=None,
                ref_clusters=None,
                seed_cluster_runs=[aria_cols["wnn1"], aria_cols["wnn2"]],
                seed_weight_runs=[aria_cols["rna_w1"], aria_cols["rna_w2"]],
            )
            concordance["reason"] = (
                f"Seurat reference unavailable: "
                f"{seurat_status.get('reason') or seurat_status.get('status')}"
            )

    from aria.version import __version__, collect_version_metadata

    prov = collect_version_metadata()
    # Strip any host-specific conda prefix before writing public artifact.
    if isinstance(prov, dict):
        prov.pop("conda_prefix", None)

    manifest = _sanitize_public_value({
        "benchmark": "V48_I2_wnn_external_concordance",
        "scope": "paired_multiome_wnn_vs_seurat_findmultimodalneighbors",
        "aria_version": __version__,
        "provenance": prov,
        "dataset": args.dataset or Path(args.h5mu).stem,
        "design": {
            "aria_method": (
                "aria/scripts/integration_wnn.py (muon.pp.neighbors over unit-L2 "
                "X_pca + depth-filtered Signac TF-IDF/SVD X_lsi, igraph Leiden on "
                "joint wnn_connectivities)"
            ),
            "reference_method": (
                "Seurat::FindMultiModalNeighbors (l2.norm=TRUE, k.nn=k) + "
                "Seurat::FindClusters (wsnn)"
            ),
            "feature_space": "identical shared barcodes, X_pca, and post-DepthCor X_lsi",
            "k": int(args.k),
            "resolution": float(args.resolution),
            "seed": int(args.seed),
        },
        "aria_summary": aria_summary,
        "seurat_status": seurat_status,
        "concordance": concordance,
        "caveats": [
            "muon.pp.neighbors and Seurat::FindMultiModalNeighbors implement the "
            "Hao et al. 2021 WNN kernel independently (differences in local "
            "bandwidth smoothing, UMAP vs SNN edge pruning, and igraph Leiden vs "
            "Seurat Louvain/SLM on wsnn produce expected non-identical partitions).",
            "Joint WNN embeddings are an exploratory lens over paired cells and "
            "never replace per-modality count-based DE/DA inference.",
        ],
    })

    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / args.manifest_name
    out_file.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {out_file.relative_to(ROOT) if out_file.is_relative_to(ROOT) else out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
