#!/usr/bin/env python3
"""V48 Slice I2 — ARIA WNN execution + export stage for external concordance.

Runs inside `aria-chromatin-env` (muon/mudata/scanpy). Executes:
  1. Optional joint RNA+ATAC QC + `.h5ad` `/obs` label join (`integration_paired_qc`)
     so raw single-donor `.h5mu` inputs undergo the exact V48 I1->I2 pipeline.
  2. Two independent executions of `integration_wnn` with identical `seed` to
     measure empirical seed determinism (`ARI == 1.0`, weight delta `0.0`).
  3. Export of shared `barcodes.tsv`, `X_pca.tsv`, `X_lsi.tsv`, ARIA cell-level
     clusters/weights (`aria_wnn_cells.tsv`), and distance-ordered non-self WNN
     neighbor barcodes (`aria_wnn_knn.tsv`) so `benchmark_wnn_seurat.R` runs
     `Seurat::FindMultiModalNeighbors` on the IDENTICAL cell and feature space.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def _extract_ordered_neighbors(dist_csr, barcodes: list[str], k: int) -> list[list[str]]:
    """Return distance-sorted non-self neighbor barcodes per cell from CSR graph."""
    import numpy as np

    indptr = dist_csr.indptr
    indices = dist_csr.indices
    data = dist_csr.data
    out: list[list[str]] = []
    for row_i in range(dist_csr.shape[0]):
        start, end = int(indptr[row_i]), int(indptr[row_i + 1])
        cols = indices[start:end]
        vals = data[start:end]
        order = np.argsort(vals, kind="mergesort")
        nbs = [
            barcodes[int(cols[pos])]
            for pos in order
            if int(cols[pos]) != row_i
        ][:k]
        out.append(nbs)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5mu", required=True, help="Paired RNA+ATAC .h5mu input")
    p.add_argument("--out-dir", required=True, help="Working directory for exports")
    p.add_argument("--run-paired-qc", action="store_true",
                   help="Run integration_paired_qc (I1) before integration_wnn (I2)")
    p.add_argument("--reference-h5ad", default=None,
                   help="Optional reference .h5ad for obs_barcode_join label transfer")
    p.add_argument("--reference-label-col", default=None)
    p.add_argument("--reference-sample-col", default=None)
    p.add_argument("--sample-id", default=None)
    p.add_argument("--label-col", default=None,
                   help="Column in mdata.obs carrying external cell labels")
    p.add_argument("--k", type=int, default=20)
    p.add_argument("--resolution", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)

    import mudata as md
    import numpy as np

    from aria.scripts.integration_paired_qc import integration_paired_qc
    from aria.scripts.integration_wnn import integration_wnn

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    wnn_input_h5mu = str(args.h5mu)
    qc_summary: dict | None = None
    inferred_label_col = args.label_col

    if args.run_paired_qc:
        qc_dir = out / "paired_qc"
        qc_params: dict = {
            "h5mu_path": str(args.h5mu),
            "output_dir": str(qc_dir),
        }
        if args.sample_id:
            qc_params["sample_id"] = args.sample_id
        if args.reference_h5ad and args.reference_label_col:
            qc_params["reference_h5ad"] = str(args.reference_h5ad)
            qc_params["reference_label_col"] = str(args.reference_label_col)
            if args.reference_sample_col:
                qc_params["reference_sample_col"] = str(args.reference_sample_col)
            if not inferred_label_col:
                inferred_label_col = str(args.reference_label_col)

        qc_res = integration_paired_qc(qc_params)
        if qc_res.get("status") != "success":
            (out / "aria_wnn_summary.json").write_text(
                json.dumps({"status": "error", "stage": "paired_qc", "result": qc_res}, indent=2),
                encoding="utf-8",
            )
            return 1
        wnn_input_h5mu = str(qc_res["output_path"])
        qc_summary = {
            "n_cells_before": qc_res.get("n_cells_before"),
            "n_cells_after": qc_res.get("n_cells_after"),
            "n_genes": qc_res.get("n_genes"),
            "n_peaks": qc_res.get("n_peaks"),
            "label_join": qc_res.get("label_join"),
        }

    # Run 1 (primary) and Run 2 (determinism check) with identical seed.
    run1_dir = out / "wnn_run1"
    run2_dir = out / "wnn_run2"
    base_params = {
        "h5mu_path": wnn_input_h5mu,
        "k": int(args.k),
        "resolution": float(args.resolution),
        "seed": int(args.seed),
    }
    res1 = integration_wnn({**base_params, "output_dir": str(run1_dir)})
    if res1.get("status") != "success":
        (out / "aria_wnn_summary.json").write_text(
            json.dumps({"status": "error", "stage": "wnn_run1", "result": res1}, indent=2),
            encoding="utf-8",
        )
        return 1

    res2 = integration_wnn({**base_params, "output_dir": str(run2_dir)})
    if res2.get("status") != "success":
        (out / "aria_wnn_summary.json").write_text(
            json.dumps({"status": "error", "stage": "wnn_run2", "result": res2}, indent=2),
            encoding="utf-8",
        )
        return 1

    mdata1 = md.read_h5mu(str(res1["output_path"]))
    mdata2 = md.read_h5mu(str(res2["output_path"]))

    barcodes = [str(bc) for bc in mdata1.obs_names]
    pca_mat = np.asarray(mdata1.mod["rna"].obsm["X_pca"], dtype=np.float64)
    lsi_mat = np.asarray(mdata1.mod["atac"].obsm["X_lsi"], dtype=np.float64)

    (out / "barcodes.tsv").write_text("\n".join(barcodes) + "\n", encoding="utf-8")
    np.savetxt(out / "X_pca.tsv", pca_mat, fmt="%.8g", delimiter="\t")
    np.savetxt(out / "X_lsi.tsv", lsi_mat, fmt="%.8g", delimiter="\t")

    labels_col = None
    if inferred_label_col and inferred_label_col in mdata1.obs.columns:
        labels_col = inferred_label_col

    with (out / "aria_wnn_cells.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "barcode",
            "leiden_wnn",
            "leiden_wnn_run2",
            "leiden_rna",
            "leiden_atac",
            "rna_weight",
            "rna_weight_run2",
            "atac_weight",
            "reference_label",
        ])
        w_wnn1 = [str(x) for x in mdata1.obs["leiden_wnn"].values]
        w_wnn2 = [str(x) for x in mdata2.obs["leiden_wnn"].values]
        w_rna = [str(x) for x in mdata1.obs["leiden_rna"].values]
        w_atac = [str(x) for x in mdata1.obs["leiden_atac"].values]
        rw1 = np.asarray(mdata1.obs["rna_weight"].values, dtype=np.float64)
        rw2 = np.asarray(mdata2.obs["rna_weight"].values, dtype=np.float64)
        aw1 = np.asarray(mdata1.obs["atac_weight"].values, dtype=np.float64)
        ref_lbls = (
            ["" if x is None else str(x) for x in mdata1.obs[labels_col].values]
            if labels_col
            else [""] * len(barcodes)
        )
        for idx, bc in enumerate(barcodes):
            writer.writerow([
                bc,
                w_wnn1[idx],
                w_wnn2[idx],
                w_rna[idx],
                w_atac[idx],
                f"{rw1[idx]:.8f}",
                f"{rw2[idx]:.8f}",
                f"{aw1[idx]:.8f}",
                ref_lbls[idx],
            ])

    knn_lists = _extract_ordered_neighbors(
        mdata1.obsp["wnn_distances"], barcodes, int(res1["k"])
    )
    with (out / "aria_wnn_knn.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["barcode", "neighbors"])
        for bc, nbs in zip(barcodes, knn_lists):
            writer.writerow([bc, ",".join(nbs)])

    umap_max_diff = float(
        np.max(
            np.abs(
                np.asarray(mdata1.obsm["X_wnn_umap"], dtype=np.float64)
                - np.asarray(mdata2.obsm["X_wnn_umap"], dtype=np.float64)
            )
        )
    )

    summary = {
        "status": "success",
        "paired_qc": qc_summary,
        "aria_wnn": {
            k_name: v_val
            for k_name, v_val in res1.items()
            if k_name not in ("joint_clusters", "output_path")
        },
        "determinism_umap_max_abs_diff": round(umap_max_diff, 8),
        "label_column_used": labels_col,
    }
    (out / "aria_wnn_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
