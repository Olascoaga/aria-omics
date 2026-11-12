#!/usr/bin/env python3
"""External SCENIC+ 1.0.0 reference runner for V48 Slice I6 concordance.

Executed inside `scenicplus_env` (`scenicplus==1.0.0`, `arboreto`, `pyranges`).
Runs the unmodified upstream SCENIC+ enhancer-GRN workflow offline on the exact
post-QC paired feature/cell universe and fixed-500bp JASPAR2024 motif cistromes:
  1. `SCENICPLUS` object construction (`X_ACC` regions x cells, `X_EXP` cells x genes,
     `uns['Cistromes']['Unfiltered']` PyRanges per expressed TF).
  2. `get_search_space` (`pr_annot`, `pr_chromsizes`, `upstream=[1000, distance_bp]`,
     `downstream=[1000, distance_bp]`, zero network/BioMart calls).
  3. `calculate_regions_to_genes_relationships` (`GBM` + Spearman `SR`).
  4. `calculate_TFs_to_genes_relationships` (`arboreto` `GBM` + Spearman `rho`).
  5. `gsea_approach.build_grn` (`cistromes_key='Unfiltered'`) + `format_egrns`.

Propagates `--seed` into both the driver process and every Ray worker via a
temporary `PYTHONPATH` `sitecustomize.py` hook so `seed=0` vs `seed=1` measures
SCENIC+'s genuine end-to-end GBM + GSEA stochasticity ceiling without modifying
any installed environment files.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import tempfile
import time
from typing import Any


def _install_ray_seed_sitecustomize(seed: int, temp_dir: Path) -> None:
    sc_dir = temp_dir / "sitecustomize_seed"
    sc_dir.mkdir(parents=True, exist_ok=True)
    (sc_dir / "sitecustomize.py").write_text(
        "import os\n"
        "import random\n"
        "import numpy as np\n"
        "_s = int(os.environ.get('ARIA_SCENICPLUS_SEED', '0'))\n"
        "random.seed(_s)\n"
        "np.random.seed(_s)\n"
        "try:\n"
        "    import arboreto.core as _ac\n"
        "    _ac.DEMON_SEED = _s\n"
        "    _ac.fit_model.__defaults__ = (_ac.EARLY_STOP_WINDOW_LENGTH, _s)\n"
        "except Exception:\n"
        "    pass\n",
        encoding="utf-8",
    )
    os.environ["ARIA_SCENICPLUS_SEED"] = str(int(seed))
    prev_pp = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = f"{sc_dir}:{prev_pp}" if prev_pp else str(sc_dir)


def _compute_aucell(
    counts_csc,
    feature_names: list[str],
    gene_sets: dict[str, list[str]],
    top_frac: float = 0.05,
) -> tuple[list[str], Any]:
    import numpy as np

    n_cells, n_feat = counts_csc.shape
    keys = sorted(gene_sets.keys())
    if n_cells == 0 or n_feat == 0 or not keys:
        return keys, np.zeros((n_cells, len(keys)), dtype=np.float32)

    k_top = max(1, int(round(float(top_frac) * float(n_feat))))
    feat_to_idx = {f: i for i, f in enumerate(feature_names)}
    dense = counts_csc.toarray().astype(np.float64)
    # Deterministic tie-breaking by feature column index (matching ARIA AUCell)
    order = np.argsort(-dense, axis=1, kind="mergesort")
    top_idx = order[:, :k_top]  # (n_cells, k_top)
    weights = (float(k_top) - np.arange(k_top, dtype=np.float64))  # k, k-1, ..., 1

    scores = np.zeros((n_cells, len(keys)), dtype=np.float32)
    for m_i, key in enumerate(keys):
        g_idxs = sorted({feat_to_idx[g] for g in gene_sets[key] if g in feat_to_idx})
        if not g_idxs:
            continue
        m_sz = len(g_idxs)
        n_max = min(m_sz, k_top)
        max_auc = float(n_max) * float(k_top) - float(n_max * (n_max - 1)) / 2.0
        if max_auc <= 0.0:
            continue
        in_set = np.isin(top_idx, np.asarray(g_idxs, dtype=top_idx.dtype))
        auc = (in_set * weights[None, :]).sum(axis=1) / max_auc
        scores[:, m_i] = auc.astype(np.float32)
    return keys, scores


def _sign_from_eregulon_context(ctx: frozenset[str] | set[str]) -> str:
    pos_tf = "positive tf2g" in ctx
    pos_r2g = "positive r2g" in ctx
    if pos_tf and pos_r2g:
        return "pos_tf_pos_p2g"
    if (not pos_tf) and pos_r2g:
        return "neg_tf_pos_p2g"
    if pos_tf and (not pos_r2g):
        return "pos_tf_neg_p2g"
    return "neg_tf_neg_p2g"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ray-cpus", type=int, default=16)
    ap.add_argument("--distance-bp", type=int, default=500_000)
    ap.add_argument("--min-target-genes", type=int, default=5)
    ap.add_argument("--gsea-n-perm", type=int, default=1000)
    ap.add_argument("--out-triplets", required=True)
    ap.add_argument("--out-aucell-npz", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args(argv)

    t_start = time.time()
    in_dir = Path(args.in_dir)

    with tempfile.TemporaryDirectory(prefix=f"scplus_seed{args.seed}_") as tmp_raw:
        tmp_dir = Path(tmp_raw)
        ray_tmp = tmp_dir / "ray_tmp"
        ray_tmp.mkdir(parents=True, exist_ok=True)
        _install_ray_seed_sitecustomize(int(args.seed), tmp_dir)

        import numpy as np
        import pandas as pd
        import pyranges as pr
        import scipy.io as sio
        import scipy.sparse as sp

        import arboreto.core as _ac
        import scenicplus
        from scenicplus.TF_to_gene import calculate_TFs_to_genes_relationships
        from scenicplus.enhancer_to_gene import (
            SGBM_KWARGS as E2G_SGBM_KWARGS,
            calculate_regions_to_genes_relationships,
            get_search_space,
        )
        from scenicplus.grn_builder.gsea_approach import build_grn
        from scenicplus.scenicplus_class import SCENICPLUS
        from scenicplus.utils import format_egrns

        np.random.seed(int(args.seed))
        _ac.DEMON_SEED = int(args.seed)
        _ac.fit_model.__defaults__ = (_ac.EARLY_STOP_WINDOW_LENGTH, int(args.seed))

        barcodes = [
            x.strip()
            for x in (in_dir / "barcodes.tsv").read_text(encoding="utf-8").splitlines()
            if x.strip()
        ]
        genes = [
            x.strip()
            for x in (in_dir / "genes.tsv").read_text(encoding="utf-8").splitlines()
            if x.strip()
        ]
        regions = [
            x.strip()
            for x in (in_dir / "regions.tsv").read_text(encoding="utf-8").splitlines()
            if x.strip()
        ]

        rna_csc = sp.csc_matrix(sio.mmread(str(in_dir / "rna_counts.mtx")))
        atac_csc = sp.csc_matrix(sio.mmread(str(in_dir / "atac_counts.mtx")))
        if rna_csc.shape != (len(barcodes), len(genes)):
            raise ValueError(
                f"RNA shape mismatch: {rna_csc.shape} vs ({len(barcodes)}, {len(genes)})"
            )
        if atac_csc.shape != (len(barcodes), len(regions)):
            raise ValueError(
                f"ATAC shape mismatch: {atac_csc.shape} vs ({len(barcodes)}, {len(regions)})"
            )

        # Standard CP10k log1p normalization for SCENIC+ X_EXP (cells x genes)
        # and X_ACC (regions x cells)
        rna_dense = rna_csc.toarray().astype(np.float32)
        rna_lib = np.maximum(rna_dense.sum(axis=1, keepdims=True), 1.0)
        x_exp = np.log1p((rna_dense / rna_lib) * 1e4).astype(np.float32)

        atac_csr = atac_csc.tocsr().astype(np.float32)
        atac_lib = np.maximum(np.asarray(atac_csr.sum(axis=1)).ravel(), 1.0)
        inv_atac = sp.diags((1e4 / atac_lib).astype(np.float32))
        x_acc_cells_x_reg = inv_atac @ atac_csr
        x_acc_cells_x_reg.data = np.log1p(x_acc_cells_x_reg.data).astype(np.float32)
        x_acc = x_acc_cells_x_reg.T.tocsr()  # (n_regions, n_cells)

        meta_regions = pd.DataFrame(index=pd.Index(regions, name="Region"))
        meta_genes = pd.DataFrame(index=pd.Index(genes, name="Gene"))
        meta_cells = pd.DataFrame(index=pd.Index(barcodes, name="Cell"))

        scplus_obj = SCENICPLUS(
            X_ACC=x_acc,
            X_EXP=x_exp,
            metadata_regions=meta_regions,
            metadata_genes=meta_genes,
            metadata_cell=meta_cells,
            menr={},
        )

        # Load shared JASPAR2024 cistromes (`{tf}_direct` -> PyRanges)
        cist_df = pd.read_csv(in_dir / "cistromes.tsv", sep="\t")
        cistromes_unfiltered: dict[str, pr.PyRanges] = {}
        for tf, grp in cist_df.groupby("tf", sort=True):
            tf_str = str(tf).strip()
            if not tf_str or "_" in tf_str:
                continue
            cistromes_unfiltered[f"{tf_str}_direct"] = pr.PyRanges(
                chromosomes=grp["Chromosome"].astype(str).to_numpy(),
                starts=grp["Start"].astype(np.int32).to_numpy(),
                ends=grp["End"].astype(np.int32).to_numpy(),
            )
        scplus_obj.uns["Cistromes"] = {"Unfiltered": cistromes_unfiltered}

        annot_df = pd.read_csv(in_dir / "annot.tsv", sep="\t")
        chrom_df = pd.read_csv(in_dir / "chromsizes.tsv", sep="\t")
        pr_annot = pr.PyRanges(annot_df)
        pr_chromsizes = pr.PyRanges(chrom_df)

        t0 = time.time()
        get_search_space(
            scplus_obj,
            species=None,
            assembly=None,
            pr_annot=pr_annot,
            pr_chromsizes=pr_chromsizes,
            use_gene_boundaries=False,
            upstream=[1000, int(args.distance_bp)],
            downstream=[1000, int(args.distance_bp)],
            remove_promoters=False,
            inplace=True,
        )
        t_search = round(time.time() - t0, 2)

        ray_kwargs = {
            "_temp_dir": str(ray_tmp),
            "include_dashboard": False,
            "ignore_reinit_error": True,
        }

        t1 = time.time()
        calculate_regions_to_genes_relationships(
            scplus_obj,
            search_space_key="search_space",
            mask_expr_dropout=False,
            importance_scoring_method="GBM",
            importance_scoring_kwargs=dict(E2G_SGBM_KWARGS),
            correlation_scoring_method="SR",
            ray_n_cpu=int(args.ray_cpus),
            add_distance=True,
            inplace=True,
            **ray_kwargs,
        )
        t_r2g = round(time.time() - t1, 2)

        t2 = time.time()
        calculate_TFs_to_genes_relationships(
            scplus_obj,
            tf_file=str(in_dir / "tfs.txt"),
            method="GBM",
            ray_n_cpu=int(args.ray_cpus),
            key="TF2G_adj",
            **ray_kwargs,
        )
        t_tf2g = round(time.time() - t2, 2)

        t3 = time.time()
        build_grn(
            scplus_obj,
            adj_key="TF2G_adj",
            cistromes_key="Unfiltered",
            region_to_gene_key="region_to_gene",
            order_regions_to_genes_by="importance",
            order_TFs_to_genes_by="importance",
            gsea_n_perm=int(args.gsea_n_perm),
            quantiles=(0.85, 0.90),
            top_n_regionTogenes_per_gene=(5, 10, 15),
            top_n_regionTogenes_per_region=(),
            binarize_using_basc=False,
            min_regions_per_gene=0,
            rho_dichotomize_tf2g=True,
            rho_dichotomize_r2g=True,
            rho_dichotomize_eregulon=True,
            keep_only_activating=False,
            rho_threshold=0.05,
            NES_thr=0.0,
            adj_pval_thr=0.05,
            min_target_genes=int(args.min_target_genes),
            inplace=True,
            key_added="eRegulons",
            ray_n_cpu=int(args.ray_cpus),
            merge_eRegulons=True,
            keep_extended_motif_annot=False,
            disable_tqdm=True,
            **ray_kwargs,
        )
        t_grn = round(time.time() - t3, 2)

        ereg_list = list(scplus_obj.uns.get("eRegulons") or [])
        if ereg_list:
            format_egrns(
                scplus_obj,
                eregulons_key="eRegulons",
                TF2G_key="TF2G_adj",
                key_added="eRegulon_metadata",
            )
            meta_df = scplus_obj.uns["eRegulon_metadata"]
        else:
            meta_df = pd.DataFrame()

        # Extract triplet rows + per-module gene sets directly from `eRegulons`
        triplet_rows: list[dict[str, Any]] = []
        module_gene_sets: dict[str, set[str]] = {}
        n_self_loops_dropped = 0
        seen_trip: set[tuple[str, str, str, str]] = set()

        for ereg in ereg_list:
            tf = str(ereg.transcription_factor).strip()
            sc = _sign_from_eregulon_context(set(ereg.context))
            mod_key = f"{tf}__{sc}"
            for r2g in ereg.regions2genes:
                region = str(r2g[0]).strip()
                gene = str(r2g[1]).strip()
                r2g_imp = float(r2g[2])
                r2g_rho = float(r2g[3])
                if tf.upper() == gene.upper():
                    n_self_loops_dropped += 1
                    continue
                module_gene_sets.setdefault(mod_key, set()).add(gene)
                key4 = (tf, region, gene, sc)
                if key4 in seen_trip:
                    continue
                seen_trip.add(key4)
                triplet_rows.append(
                    {
                        "tf": tf,
                        "peak": region,
                        "gene": gene,
                        "sign_class": sc,
                        "r2g_importance": round(r2g_imp, 6),
                        "r2g_rho": round(r2g_rho, 6),
                        "gsea_nes": (
                            round(float(ereg.gsea_enrichment_score), 6)
                            if ereg.gsea_enrichment_score is not None
                            else None
                        ),
                        "gsea_adj_pval": (
                            float(ereg.gsea_adj_pval)
                            if ereg.gsea_adj_pval is not None
                            else None
                        ),
                    }
                )

        # Enforce `min_target_genes` AFTER dropping SCENIC+'s injected TF self-loop
        valid_mod_keys = {
            k for k, gs in module_gene_sets.items() if len(gs) >= int(args.min_target_genes)
        }
        triplet_rows = [
            r
            for r in triplet_rows
            if f"{r['tf']}__{r['sign_class']}" in valid_mod_keys
        ]
        filtered_gene_sets = {
            k: sorted(module_gene_sets[k]) for k in sorted(valid_mod_keys)
        }

        out_trip_p = Path(args.out_triplets)
        out_trip_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_trip_p, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "tf",
                    "peak",
                    "gene",
                    "sign_class",
                    "r2g_importance",
                    "r2g_rho",
                    "gsea_nes",
                    "gsea_adj_pval",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for r in triplet_rows:
                writer.writerow(r)

        mod_keys, aucell_mat = _compute_aucell(
            rna_csc, genes, filtered_gene_sets, top_frac=0.05
        )
        np.savez_compressed(
            str(args.out_aucell_npz),
            module_keys=np.asarray(mod_keys, dtype=object),
            aucell=aucell_mat,
        )

        summary = {
            "status": "success",
            "scenicplus_version": str(getattr(scenicplus, "__version__", "1.0.0")),
            "seed": int(args.seed),
            "n_cells": len(barcodes),
            "n_genes": len(genes),
            "n_regions": len(regions),
            "n_cistromes_unfiltered": len(cistromes_unfiltered),
            "n_search_space_pairs": int(len(scplus_obj.uns["search_space"])),
            "n_eregulons_raw": len(ereg_list),
            "n_eregulons_post_selfloop_filter": len(valid_mod_keys),
            "n_unique_tfs": len({r["tf"] for r in triplet_rows}),
            "n_triplets": len(triplet_rows),
            "n_self_loops_dropped": int(n_self_loops_dropped),
            "n_metadata_rows_raw": int(len(meta_df)),
            "timings_sec": {
                "search_space": t_search,
                "regions_to_genes_gbm": t_r2g,
                "tfs_to_genes_gbm": t_tf2g,
                "build_grn_gsea": t_grn,
                "total": round(time.time() - t_start, 2),
            },
        }
        Path(args.out_json).write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
