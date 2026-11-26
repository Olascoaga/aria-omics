#!/usr/bin/env python3
"""Real-data donor-label permutation proof of Contract 1 on Kang et al. (GSE96583).

Evaluates pseudoreplication inflation on REAL unstimulated human PBMC single-cell
UMI counts from 8 biological donors (GEO GSE96583, Kang et al. Nat Biotechnol
2018; staged at ``~/.aria/public_atlas_datasets/03_scrna_multidonor``).

Because all 8 donor libraries come from the identical unstimulated control arm
(zero true perturbation effect; real inter-donor biological heterogeneity), any
balanced 4-vs-4 donor split defines an exact real-data global null. Across all
35 unique balanced 4-vs-4 donor partitions (C(8,4)/2 = 35) and major PBMC cell
types, this script compares on the EXACT SAME cells and expressed genes:

1. Naive cell-level Wilcoxon rank-sum test (Scanpy ``rank_genes_groups``, BH FDR < 0.05)
2. Naive cell-level Welch t-test (BH FDR < 0.05)
3. ARIA donor-aware pseudobulk pyDESeq2 (BH FDR < 0.05)

Writes a path-portable, provenance-stamped JSON receipt (zero machine-absolute
paths) suitable for direct publication verification.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
import math
import os
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ttest_ind

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aria.version import collect_version_metadata  # noqa: E402


DEFAULT_H5AD = (
    Path.home()
    / ".aria"
    / "public_atlas_datasets"
    / "03_scrna_multidonor"
    / "inputs"
    / "kang2018_pbmc_8donors_ifnb.h5ad"
)
DEFAULT_OUTPUT = (
    ROOT
    / "docs"
    / "benchmark_results"
    / "v494_hardening"
    / "kang_real_pseudoreplication_null.json"
)


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    out = np.ones(n, dtype=float)
    valid = np.isfinite(p)
    if not np.any(valid):
        return out
    pv = np.clip(p[valid], 0.0, 1.0)
    m = pv.size
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    q = ranked * (m / np.arange(1, m + 1, dtype=float))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    inv = np.empty(m, dtype=float)
    inv[order] = q
    out[valid] = inv
    return out


def _unique_balanced_splits(donors: list[str]) -> list[tuple[list[str], list[str]]]:
    """Return all C(n, n/2)/2 unique unordered balanced splits (pinned first donor)."""
    donors = sorted(donors)
    half = len(donors) // 2
    anchor = donors[0]
    rest = donors[1:]
    splits: list[tuple[list[str], list[str]]] = []
    for combo in combinations(rest, half - 1):
        grp_a = sorted([anchor, *combo])
        grp_b = sorted([d for d in donors if d not in grp_a])
        splits.append((grp_a, grp_b))
    return splits


def _summarize_counts(vals: list[int]) -> dict[str, float | int]:
    arr = np.asarray(vals, dtype=float)
    return {
        "n_tests": int(arr.size),
        "total_false_positives": int(np.sum(arr)),
        "mean_false_positives": round(float(np.mean(arr)), 2),
        "median_false_positives": round(float(np.median(arr)), 2),
        "iqr_false_positives": [
            round(float(np.percentile(arr, 25)), 2),
            round(float(np.percentile(arr, 75)), 2),
        ],
        "max_false_positives": int(np.max(arr)),
        "min_false_positives": int(np.min(arr)),
        "fwer_any_fp_rate": round(float(np.mean(arr > 0)), 4),
        "n_splits_with_ge1_fp": int(np.sum(arr > 0)),
    }


def run_benchmark(
    h5ad_path: Path,
    *,
    alpha: float = 0.05,
    min_cells_per_donor: int = 10,
    min_donors_per_arm: int = 4,
    min_total_counts: int = 20,
    min_donors_detected: int = 3,
    max_splits: int | None = None,
) -> dict:
    import anndata as ad
    import scanpy as sc
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    warnings.filterwarnings("ignore")

    adata_full = ad.read_h5ad(str(h5ad_path))
    ctrl = adata_full[adata_full.obs["condition"].astype(str) == "control"].copy()
    donors_all = sorted(ctrl.obs["donor_id"].astype(str).unique().tolist())
    assert len(donors_all) == 8, f"Expected 8 control donors, found {len(donors_all)}"

    all_splits = _unique_balanced_splits(donors_all)
    if max_splits is not None and max_splits < len(all_splits):
        all_splits = all_splits[:max_splits]

    # Eligible cell types: all 8 donors have >= min_cells_per_donor control cells
    ct_counts = pd.crosstab(
        ctrl.obs["author_cell_type"].astype(str),
        ctrl.obs["donor_id"].astype(str),
    )
    eligible_cts = sorted(
        ct for ct in ct_counts.index
        if (ct_counts.loc[ct] >= min_cells_per_donor).all()
    )

    per_cell_type: dict[str, dict] = {}
    all_wilcoxon_fps: list[int] = []
    all_welch_fps: list[int] = []
    all_pseudobulk_fps: list[int] = []
    all_genes_tested: list[int] = []

    for ct in eligible_cts:
        sub = ctrl[ctrl.obs["author_cell_type"].astype(str) == ct].copy()
        X_raw = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
        X_raw = np.rint(X_raw).astype(np.int64)
        obs_donors = sub.obs["donor_id"].astype(str).to_numpy()
        gene_names = np.asarray(sub.var_names.astype(str))

        # Aggregate exact per-donor pseudobulk raw counts
        pb_rows = {}
        for d in donors_all:
            pb_rows[d] = X_raw[obs_donors == d].sum(axis=0)
        pb_df = pd.DataFrame(pb_rows, index=gene_names).T  # donors x genes

        # Filter to genes expressed in >= min_donors_detected donors and total >= min_total_counts
        keep_mask = (
            ((pb_df > 0).sum(axis=0) >= min_donors_detected)
            & (pb_df.sum(axis=0) >= min_total_counts)
        ).to_numpy()
        pb_filt = pb_df.loc[:, keep_mask].copy()
        X_filt = X_raw[:, keep_mask]
        genes_filt = gene_names[keep_mask]
        n_genes = int(genes_filt.size)
        n_cells = int(X_filt.shape[0])

        # Precompute library-size log1p-normalized expression for cell-level tests
        lib = X_raw.sum(axis=1, keepdims=True)
        lib = np.where(lib > 0, lib, 1.0)
        X_log = np.log1p(X_filt * (1e4 / lib))

        # Build reusable AnnData for Scanpy Wilcoxon
        ad_ct = ad.AnnData(
            X=X_log.astype(np.float32),
            obs=pd.DataFrame({"donor_id": obs_donors}, index=sub.obs_names),
            var=pd.DataFrame(index=genes_filt),
        )

        ct_wilcoxon: list[int] = []
        ct_welch: list[int] = []
        ct_pb: list[int] = []

        for split_idx, (grp_a, grp_b) in enumerate(all_splits):
            grp_a_set = set(grp_a)
            cell_grp = np.where(np.isin(obs_donors, grp_a), "pseudo_A", "pseudo_B")

            # 1) Naive cell-level Scanpy Wilcoxon
            ad_ct.obs["pseudo_cond"] = pd.Categorical(cell_grp, categories=["pseudo_A", "pseudo_B"])
            sc.tl.rank_genes_groups(
                ad_ct,
                groupby="pseudo_cond",
                groups=["pseudo_A"],
                reference="pseudo_B",
                method="wilcoxon",
                pts=False,
                use_raw=False,
                tie_correct=True,
            )
            w_df = sc.get.rank_genes_groups_df(ad_ct, group="pseudo_A")
            n_fp_wilcox = int((w_df["pvals_adj"].fillna(1.0) < alpha).sum())

            # 2) Naive cell-level Welch t-test (BH-adjusted)
            mask_a = cell_grp == "pseudo_A"
            mask_b = ~mask_a
            _, p_welch = ttest_ind(
                X_log[mask_a], X_log[mask_b], axis=0, equal_var=False, nan_policy="omit"
            )
            q_welch = _bh_adjust(np.where(np.isfinite(p_welch), p_welch, 1.0))
            n_fp_welch = int((q_welch < alpha).sum())

            # 3) ARIA donor-aware pseudobulk pyDESeq2
            meta = pd.DataFrame(
                {"condition": ["pseudo_A" if d in grp_a_set else "pseudo_B" for d in pb_filt.index]},
                index=pb_filt.index,
            )
            dds = DeseqDataSet(
                counts=pb_filt,
                metadata=meta,
                design="~condition",
                refit_cooks=True,
                n_cpus=1,
                quiet=True,
            )
            dds.deseq2()
            ds = DeseqStats(
                dds,
                contrast=["condition", "pseudo_A", "pseudo_B"],
                n_cpus=1,
                quiet=True,
            )
            ds.summary()
            res = ds.results_df
            padj = res["padj"].fillna(1.0).to_numpy(dtype=float)
            n_fp_pb = int((padj < alpha).sum())

            ct_wilcoxon.append(n_fp_wilcox)
            ct_welch.append(n_fp_welch)
            ct_pb.append(n_fp_pb)
            all_wilcoxon_fps.append(n_fp_wilcox)
            all_welch_fps.append(n_fp_welch)
            all_pseudobulk_fps.append(n_fp_pb)
            all_genes_tested.append(n_genes)

        per_cell_type[ct] = {
            "n_cells": n_cells,
            "n_genes_tested": n_genes,
            "n_donors": len(donors_all),
            "donor_cell_counts": {d: int(ct_counts.loc[ct, d]) for d in donors_all},
            "naive_cell_wilcoxon": _summarize_counts(ct_wilcoxon),
            "naive_cell_welch": _summarize_counts(ct_welch),
            "aria_pseudobulk_deseq2": _summarize_counts(ct_pb),
        }

    total_gene_tests = int(sum(all_genes_tested))
    tot_w = int(sum(all_wilcoxon_fps))
    tot_t = int(sum(all_welch_fps))
    tot_pb = int(sum(all_pseudobulk_fps))

    prov = collect_version_metadata(ROOT)
    if isinstance(prov.get("environment"), dict):
        prov["environment"]["conda_prefix"] = None

    return {
        "schema_version": "aria.benchmark.kang_real_pseudoreplication_null.v1",
        "dataset": {
            "name": "Kang et al. 2018 Human Lupus PBMC (unstimulated control arm)",
            "accession": "GEO GSE96583 / Bioconductor EH2259",
            "doi": "10.1038/nbt.4042",
            "input_ref": "public-atlas://03_scrna_multidonor/inputs/kang2018_pbmc_8donors_ifnb.h5ad",
            "n_donors": len(donors_all),
            "n_control_cells_total": int(ctrl.n_obs),
            "n_eligible_cell_types": len(eligible_cts),
            "eligible_cell_types": eligible_cts,
        },
        "design": {
            "null_construction": "all_unique_balanced_4v4_donor_splits_on_unstimulated_control_pbmcs",
            "n_splits": len(all_splits),
            "total_celltype_split_runs": len(all_pseudobulk_fps),
            "alpha_bh_fdr": alpha,
            "min_cells_per_donor": min_cells_per_donor,
            "min_total_counts": min_total_counts,
            "min_donors_detected": min_donors_detected,
        },
        "overall": {
            "total_gene_tests": total_gene_tests,
            "naive_cell_wilcoxon": {
                **_summarize_counts(all_wilcoxon_fps),
                "per_gene_false_discovery_fraction": round(tot_w / max(total_gene_tests, 1), 6),
            },
            "naive_cell_welch": {
                **_summarize_counts(all_welch_fps),
                "per_gene_false_discovery_fraction": round(tot_t / max(total_gene_tests, 1), 6),
            },
            "aria_pseudobulk_deseq2": {
                **_summarize_counts(all_pseudobulk_fps),
                "per_gene_false_discovery_fraction": round(tot_pb / max(total_gene_tests, 1), 6),
            },
            "reduction_factor_vs_wilcoxon": (
                round(tot_w / max(tot_pb, 1), 1) if tot_pb > 0 else float("inf")
            ),
            "reduction_factor_vs_welch": (
                round(tot_t / max(tot_pb, 1), 1) if tot_pb > 0 else float("inf")
            ),
        },
        "per_cell_type": per_cell_type,
        "provenance": prov,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--max-splits", type=int, default=None)
    args = parser.parse_args(argv)

    result = run_benchmark(Path(args.h5ad).expanduser(), max_splits=args.max_splits)
    # Replace non-finite floats before JSON serialization
    for k in ("reduction_factor_vs_wilcoxon", "reduction_factor_vs_welch"):
        v = result["overall"].get(k)
        if isinstance(v, float) and not math.isfinite(v):
            result["overall"][k] = None
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if '"/home/' in text or '"/Users/' in text:
        raise RuntimeError("Output JSON contains a machine-absolute path")
    out_path.write_text(text, encoding="utf-8")
    print(json.dumps({
        "output": str(out_path.relative_to(ROOT)) if out_path.is_relative_to(ROOT) else str(out_path),
        "eligible_cell_types": result["dataset"]["eligible_cell_types"],
        "n_splits": result["design"]["n_splits"],
        "overall": result["overall"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
