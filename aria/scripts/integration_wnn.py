"""
ARIA WNN Integration Script (V48 Slice I2)
------------------------------------------
Builds a deterministic joint RNA+ATAC embedding using Weighted Nearest Neighbors
(Hao et al. 2021, Cell; implemented via muon/scanpy in aria-chromatin-env).

Pipeline over same-cell paired multiome data (.h5mu preferred, or paired .h5ad):
  1. Verify non-negative integer raw counts for RNA and ATAC via mudata_io
     (never double-log-normalizing pre-normalized matrices).
  2. Deterministic RNA preprocessing (stamped seed):
     normalize_total(1e4) -> log1p -> HVG -> scale -> PCA -> unit L2 norm -> kNN.
  3. Deterministic ATAC preprocessing (stamped seed):
     Signac log-TF-log-IDF -> TruncatedSVD -> drop LSI components correlated with
     log10 sequencing depth (|Pearson r| >= depth_corr_cutoff) -> unit L2 norm -> kNN.
  4. Joint WNN graph (muon.pp.neighbors with clamped neighbor bounds), Leiden
     clustering directly on the joint WNN connectivity graph (plus RNA-only and
     ATAC-only Leiden on the identical cell universe), and joint WNN UMAP.
  5. Honest extraction of per-cell and per-cluster modality weights, partition
     agreement (ARI/NMI), contingency-measured discordant cluster count, and
     atomic .h5mu output write.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.utils.wnn_diagnostics import (
    compute_wnn_anti_degradation_diagnostics,
    sanitize_for_hdf5_uns,
)

_FRIP_COL_HINTS = ("frip", "fraction_reads_in_peaks", "pct_reads_in_peaks")


def _l2_normalize_rows(mat):
    """Normalize each cell embedding row to unit L2 norm (Hao et al. 2021 / muon)."""
    import numpy as np

    arr = np.asarray(mat, dtype=np.float64)
    norms = np.linalg.norm(arr, axis=1, keepdims=True)
    safe_norms = np.where(norms > 1e-12, norms, 1.0)
    return (arr / safe_norms).astype(np.float32)


def _find_measured_frip(obs_dfs: list[Any]) -> tuple[float | None, str | None]:
    """Return median FRiP only if an explicit per-cell FRiP column exists."""
    import numpy as np

    for df in obs_dfs:
        if df is None:
            continue
        for col in df.columns:
            low = str(col).lower()
            if any(h in low for h in _FRIP_COL_HINTS):
                try:
                    vals = np.asarray(df[col].values, dtype=float)
                    vals = vals[np.isfinite(vals)]
                    if vals.size:
                        med = float(np.median(vals))
                        # Convert 0..100 percentage columns to 0..1 fraction if needed.
                        if med > 1.0 and med <= 100.0 and "pct" in low:
                            med = med / 100.0
                        return round(med, 4), str(col)
                except Exception:
                    continue
    return None, None


def _count_discordant_clusters(
    wnn_labels: list[str],
    rna_labels: list[str],
    *,
    purity_threshold: float = 0.70,
) -> tuple[int, dict[str, dict[str, Any]]]:
    """Count WNN clusters whose modal RNA-only cluster share is below threshold.

    Uses the exact contingency distribution rather than a synthetic formula.
    """
    from collections import Counter, defaultdict

    buckets: dict[str, list[str]] = defaultdict(list)
    for w_lbl, r_lbl in zip(wnn_labels, rna_labels):
        buckets[str(w_lbl)].append(str(r_lbl))

    n_discordant = 0
    details: dict[str, dict[str, Any]] = {}
    for w_lbl in sorted(buckets.keys()):
        members = buckets[w_lbl]
        counts = Counter(members)
        modal_rna, modal_n = counts.most_common(1)[0]
        purity = float(modal_n) / float(len(members)) if members else 1.0
        is_disc = bool(purity < purity_threshold)
        if is_disc:
            n_discordant += 1
        details[w_lbl] = {
            "n_cells": int(len(members)),
            "modal_rna_cluster": str(modal_rna),
            "modal_rna_purity": round(purity, 4),
            "discordant": is_disc,
        }
    return n_discordant, details


def _cluster_graph_leiden(
    *,
    obs_df,
    connectivities,
    distances,
    neighbors_meta: dict,
    resolution: float,
    seed: int,
) -> list[str]:
    """Run deterministic igraph Leiden directly on a precomputed neighbor graph."""
    import anndata as ad
    import scanpy as sc

    graph_adata = ad.AnnData(
        obs=obs_df[[]].copy(),
        obsp={
            "connectivities": connectivities,
            "distances": distances,
        },
        uns={
            "neighbors": {
                "connectivities_key": "connectivities",
                "distances_key": "distances",
                "params": dict(neighbors_meta.get("params", {})),
            }
        },
    )
    sc.tl.leiden(
        graph_adata,
        resolution=float(resolution),
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=int(seed),
        key_added="leiden",
    )
    return [str(x) for x in graph_adata.obs["leiden"].values]


def integration_wnn(params: dict) -> dict:
    import numpy as np

    h5mu_param = params.get("h5mu_path") or params.get("data_path")
    rna_files = list(params.get("rna_files") or [])
    atac_files = list(params.get("atac_files") or [])
    genome = str(params.get("genome") or "hg38")
    organism = str(params.get("organism") or "Homo sapiens")
    k = int(params.get("k", 20))
    resolution = float(params.get("resolution", 0.5))
    n_top_genes = int(params.get("n_top_genes", 3000))
    n_pcs = int(params.get("n_pcs", 50))
    n_lsi_comps = int(params.get("n_lsi_comps", 50))
    depth_corr_cutoff = float(params.get("depth_corr_cutoff", 0.75))
    min_cells = max(10, int(params.get("min_cells", 50)))
    seed = int(params.get("seed", 0))
    output_dir = str(params.get("output_dir") or "/tmp/aria_wnn")
    warnings: list[str] = []

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Resolve whether input is a paired .h5mu or split modality files.
    resolved_h5mu: str | None = None
    if h5mu_param is not None and str(h5mu_param).strip():
        cand = Path(str(h5mu_param))
        if not cand.exists():
            return {
                "status": "error",
                "error_type": "MissingPairedH5mu",
                "details": f"Paired .h5mu file does not exist: {h5mu_param}",
            }
        resolved_h5mu = str(cand)
    else:
        valid_rna = [str(f) for f in rna_files if Path(str(f)).exists()]
        valid_atac = [str(f) for f in atac_files if Path(str(f)).exists()]
        if not valid_rna:
            return {
                "status": "error",
                "error_type": "MissingRNA",
                "details": "No valid RNA files found for WNN integration.",
            }
        if not valid_atac:
            return {
                "status": "error",
                "error_type": "MissingATAC",
                "details": "No valid ATAC files found for WNN integration.",
            }
        if Path(valid_rna[0]).suffix.lower() == ".h5mu":
            if Path(valid_atac[0]).suffix.lower() == ".h5mu":
                if Path(valid_rna[0]).resolve() != Path(valid_atac[0]).resolve():
                    return {
                        "status": "error",
                        "error_type": "MismatchedPairedH5mu",
                        "details": (
                            "rna_files and atac_files point to different .h5mu "
                            f"containers: {valid_rna[0]} vs {valid_atac[0]}."
                        ),
                    }
            resolved_h5mu = valid_rna[0]

    try:
        import anndata as ad
        import mudata as md
        import muon as mu
        import scanpy as sc
        import scipy.sparse as sp
        from sklearn.decomposition import TruncatedSVD

        from aria.benchmarks.concordance_atac import (
            adjusted_rand_index,
            normalized_mutual_info,
        )
        from aria.scripts.chromatin_lsi_clustering import _tfidf
        from aria.utils.mudata_io import (
            _select_raw_counts_matrix,
            read_h5mu_paired,
        )

        global_obs_preserved = None
        input_kind = "paired_h5ad"

        if resolved_h5mu is not None:
            input_kind = "paired_h5mu"
            paired = read_h5mu_paired(resolved_h5mu, require_raw_counts=True)
            if paired.get("status") != "success":
                return {
                    "status": "error",
                    "error_type": paired.get("error_type", "H5muReadFailed"),
                    "details": paired.get("details", "Could not read paired .h5mu."),
                    "warnings": warnings,
                }
            common_cells = list(paired["shared_barcodes"])
            src_mdata = paired["mdata"]
            src_rna = paired["rna_adata"]
            src_atac = paired["atac_adata"]
            rna_raw_full = paired["rna_counts_matrix"]
            atac_raw_full = paired["atac_counts_matrix"]
            rna_counts_source = str(paired["rna_counts_source"])
            atac_counts_source = str(paired["atac_counts_source"])

            # Align raw count matrices to common_cells in deterministic RNA order.
            rna_idx_map = {str(bc): i for i, bc in enumerate(src_rna.obs_names)}
            atac_idx_map = {str(bc): i for i, bc in enumerate(src_atac.obs_names)}
            r_take = [rna_idx_map[bc] for bc in common_cells]
            a_take = [atac_idx_map[bc] for bc in common_cells]

            rna_raw = sp.csr_matrix(rna_raw_full[r_take, :])
            atac_raw = sp.csr_matrix(atac_raw_full[a_take, :])

            rna_adata = ad.AnnData(
                X=rna_raw.copy(),
                obs=src_rna.obs.iloc[r_take].copy(),
                var=src_rna.var.copy(),
            )
            atac_adata = ad.AnnData(
                X=atac_raw.copy(),
                obs=src_atac.obs.iloc[a_take].copy(),
                var=src_atac.var.copy(),
            )
            if src_mdata.obs is not None and len(src_mdata.obs.columns) > 0:
                m_idx_map = {str(bc): i for i, bc in enumerate(src_mdata.obs_names)}
                if all(bc in m_idx_map for bc in common_cells):
                    m_take = [m_idx_map[bc] for bc in common_cells]
                    global_obs_preserved = src_mdata.obs.iloc[m_take].copy()
        else:
            # Split modality files: check ATAC first so raw fragment blockers
            # surface explicitly as NotImplemented (P0-8 / ADR-002).
            atac_loaded = _load_atac(valid_atac, genome, organism)
            if atac_loaded is None:
                return {
                    "status": "error",
                    "error_type": "ATACLoadFailed",
                    "details": "Could not load ATAC peak count matrix.",
                }
            rna_loaded = _load_rna(valid_rna[0])
            if rna_loaded is None:
                return {
                    "status": "error",
                    "error_type": "RNALoadFailed",
                    "details": f"Could not load RNA from {valid_rna[0]}",
                }

            rna_check = _select_raw_counts_matrix(
                rna_loaded, modality_label="rna", require_raw_counts=True
            )
            if not rna_check["ok"]:
                return {
                    "status": "error",
                    "error_type": rna_check["error_type"],
                    "details": rna_check["details"],
                }
            atac_check = _select_raw_counts_matrix(
                atac_loaded, modality_label="atac", require_raw_counts=True
            )
            if not atac_check["ok"]:
                return {
                    "status": "error",
                    "error_type": atac_check["error_type"],
                    "details": atac_check["details"],
                }

            rna_counts_source = str(rna_check["source"])
            atac_counts_source = str(atac_check["source"])

            atac_obs_set = {str(bc) for bc in atac_loaded.obs_names}
            common_cells = [
                str(bc) for bc in rna_loaded.obs_names if str(bc) in atac_obs_set
            ]
            rna_idx_map = {str(bc): i for i, bc in enumerate(rna_loaded.obs_names)}
            atac_idx_map = {str(bc): i for i, bc in enumerate(atac_loaded.obs_names)}
            r_take = [rna_idx_map[bc] for bc in common_cells]
            a_take = [atac_idx_map[bc] for bc in common_cells]

            rna_raw = sp.csr_matrix(rna_check["matrix"][r_take, :])
            atac_raw = sp.csr_matrix(atac_check["matrix"][a_take, :])

            rna_adata = ad.AnnData(
                X=rna_raw.copy(),
                obs=rna_loaded.obs.iloc[r_take].copy(),
                var=rna_loaded.var.copy(),
            )
            atac_adata = ad.AnnData(
                X=atac_raw.copy(),
                obs=atac_loaded.obs.iloc[a_take].copy(),
                var=atac_loaded.var.copy(),
            )

        n_cells = len(common_cells)
        if n_cells < 200:
            warnings.append(
                f"Only {n_cells} cells shared between RNA and ATAC. "
                f"WNN recommends at least 200 shared cells for stable weights."
            )
        if n_cells < min_cells:
            return {
                "status": "error",
                "error_type": "InsufficientSharedCells",
                "details": (
                    f"{n_cells} shared cells is below minimum ({min_cells}) for WNN. "
                    f"Verify that RNA and ATAC data come from the same cells."
                ),
                "warnings": warnings,
            }

        if rna_adata.n_vars < 3 or atac_adata.n_vars < 3:
            return {
                "status": "error",
                "error_type": "InsufficientFeatures",
                "details": (
                    f"Feature space too small for WNN: {rna_adata.n_vars} genes, "
                    f"{atac_adata.n_vars} peaks."
                ),
            }

        # Clamp neighbor search hyperparameters safely below n_cells - 1 so
        # muon.pp.neighbors (which queries n_multineighbors + 1) never overflows.
        max_valid_k = max(2, n_cells - 2)
        k_eff = min(max(2, k), max_valid_k)
        if k_eff < k:
            warnings.append(
                f"Clamped WNN k from {k} to {k_eff} for {n_cells} cells."
            )
        n_bw = min(max(2, int(params.get("n_bandwidth_neighbors", 20))), max_valid_k)
        n_multi = min(
            max(k_eff, int(params.get("n_multineighbors", 200))),
            max_valid_k,
        )

        # Preserve exact raw integer counts in layers["counts"] for downstream
        # steps (I3 diagnostics, I4 MOFA+, I5 peak2gene, I6 eRegulons).
        rna_adata.layers["counts"] = rna_raw.copy()
        atac_adata.layers["counts"] = atac_raw.copy()

        # ── 1. Deterministic RNA preprocessing ────────────────────────────
        sc.pp.normalize_total(rna_adata, target_sum=1e4)
        sc.pp.log1p(rna_adata)

        # Identify HVGs on non-constant genes without dropping non-HVG genes from
        # the output MuData container.
        gene_sums = np.asarray(rna_raw.sum(axis=0)).ravel()
        nonzero_gene_mask = gene_sums > 0
        if int(nonzero_gene_mask.sum()) < 3:
            return {
                "status": "error",
                "error_type": "EmptyRNAMatrix",
                "details": "Fewer than 3 expressed genes across shared cells.",
            }

        n_hvg_target = min(max(2, n_top_genes), int(nonzero_gene_mask.sum()))
        rna_work = rna_adata[:, nonzero_gene_mask].copy()
        if n_hvg_target < rna_work.n_vars:
            sc.pp.highly_variable_genes(
                rna_work,
                n_top_genes=n_hvg_target,
                flavor="seurat",
                subset=True,
            )
        hvg_names = set(map(str, rna_work.var_names))
        rna_adata.var["highly_variable"] = [
            str(g) in hvg_names for g in rna_adata.var_names
        ]
        n_hvg_used = int(rna_work.n_vars)

        sc.pp.scale(rna_work, max_value=10)
        n_pcs_eff = min(max(2, n_pcs), n_cells - 1, n_hvg_used - 1)
        if n_pcs_eff < n_pcs:
            warnings.append(
                f"Clamped RNA n_pcs from {n_pcs} to {n_pcs_eff} "
                f"(n_cells={n_cells}, n_hvg={n_hvg_used})."
            )
        sc.tl.pca(
            rna_work,
            n_comps=n_pcs_eff,
            svd_solver="arpack",
            random_state=seed,
        )
        rna_adata.obsm["X_pca"] = np.asarray(
            rna_work.obsm["X_pca"], dtype=np.float32
        )
        rna_adata.obsm["X_pca_l2"] = _l2_normalize_rows(rna_adata.obsm["X_pca"])

        sc.pp.neighbors(
            rna_adata,
            n_neighbors=k_eff,
            use_rep="X_pca_l2",
            metric="euclidean",
            random_state=seed,
        )
        sc.tl.umap(rna_adata, random_state=seed)
        rna_adata.obs["leiden_rna"] = _cluster_graph_leiden(
            obs_df=rna_adata.obs,
            connectivities=rna_adata.obsp["connectivities"],
            distances=rna_adata.obsp["distances"],
            neighbors_meta=rna_adata.uns["neighbors"],
            resolution=resolution,
            seed=seed,
        )

        # ── 2. Deterministic ATAC LSI + measured depth-correlation drop ───
        peak_sums = np.asarray(atac_raw.sum(axis=0)).ravel()
        nonzero_peak_mask = peak_sums > 0
        if int(nonzero_peak_mask.sum()) < 3:
            return {
                "status": "error",
                "error_type": "EmptyATACMatrix",
                "details": "Fewer than 3 accessible peaks across shared cells.",
            }

        tfidf_mat, atac_cell_totals = _tfidf(atac_raw[:, nonzero_peak_mask])
        n_lsi_eff = min(
            max(2, n_lsi_comps),
            n_cells - 1,
            int(nonzero_peak_mask.sum()) - 1,
        )
        if n_lsi_eff < n_lsi_comps:
            warnings.append(
                f"Clamped ATAC n_lsi_comps from {n_lsi_comps} to {n_lsi_eff}."
            )

        svd = TruncatedSVD(n_components=n_lsi_eff, random_state=seed)
        lsi_raw = svd.fit_transform(tfidf_mat)

        log_depth = np.log10(np.maximum(atac_cell_totals, 1.0))
        depth_corr: list[float] = []
        for idx in range(lsi_raw.shape[1]):
            comp = lsi_raw[:, idx]
            if np.std(comp) == 0.0 or np.std(log_depth) == 0.0:
                depth_corr.append(0.0)
            else:
                r_val = float(np.corrcoef(comp, log_depth)[0, 1])
                depth_corr.append(abs(r_val) if np.isfinite(r_val) else 0.0)

        dropped_lsi = [
            idx for idx, r_val in enumerate(depth_corr)
            if r_val >= depth_corr_cutoff
        ]
        keep_lsi = [idx for idx in range(lsi_raw.shape[1]) if idx not in dropped_lsi]
        if not keep_lsi:
            warnings.append(
                f"All {lsi_raw.shape[1]} LSI components exceeded depth "
                f"correlation cutoff {depth_corr_cutoff}; kept all components."
            )
            keep_lsi = list(range(lsi_raw.shape[1]))
            dropped_lsi = []
        elif dropped_lsi:
            warnings.append(
                f"Dropped depth-associated ATAC LSI component(s) {dropped_lsi} "
                f"(|r| vs log10 depth >= {depth_corr_cutoff})."
            )

        atac_adata.obsm["X_lsi"] = np.asarray(
            lsi_raw[:, keep_lsi], dtype=np.float32
        )
        atac_adata.obsm["X_lsi_l2"] = _l2_normalize_rows(atac_adata.obsm["X_lsi"])
        atac_adata.obs["atac_total_counts"] = atac_cell_totals.astype(np.float64)

        sc.pp.neighbors(
            atac_adata,
            n_neighbors=k_eff,
            use_rep="X_lsi_l2",
            metric="euclidean",
            random_state=seed,
        )
        sc.tl.umap(atac_adata, random_state=seed)
        atac_adata.obs["leiden_atac"] = _cluster_graph_leiden(
            obs_df=atac_adata.obs,
            connectivities=atac_adata.obsp["connectivities"],
            distances=atac_adata.obsp["distances"],
            neighbors_meta=atac_adata.uns["neighbors"],
            resolution=resolution,
            seed=seed,
        )

        # ── 3. Build paired MuData & run WNN ──────────────────────────────
        mdata = md.MuData({"rna": rna_adata, "atac": atac_adata})
        if global_obs_preserved is not None:
            for col in global_obs_preserved.columns:
                mdata.obs[str(col)] = global_obs_preserved[col].values

        mu.pp.neighbors(
            mdata,
            n_neighbors=k_eff,
            n_bandwidth_neighbors=n_bw,
            n_multineighbors=n_multi,
            key_added="wnn",
            weight_key="mod_weight",
            random_state=seed,
        )
        mu.tl.umap(mdata, neighbors_key="wnn", random_state=seed)
        mdata.obsm["X_wnn_umap"] = np.asarray(
            mdata.obsm["X_umap"], dtype=np.float32
        ).copy()

        # Cluster directly on the joint WNN graph (never multiplex per-modality).
        wnn_labels = _cluster_graph_leiden(
            obs_df=mdata.obs,
            connectivities=mdata.obsp["wnn_connectivities"],
            distances=mdata.obsp["wnn_distances"],
            neighbors_meta=mdata.uns["wnn"],
            resolution=resolution,
            seed=seed,
        )
        rna_labels = [str(x) for x in rna_adata.obs["leiden_rna"].values]
        atac_labels = [str(x) for x in atac_adata.obs["leiden_atac"].values]

        mdata.obs["leiden_wnn"] = wnn_labels
        mdata.obs["leiden_rna"] = rna_labels
        mdata.obs["leiden_atac"] = atac_labels

        # ── 4. Extract real modality weights (no fabrication) ─────────────
        rna_w_series = mdata.obs.get("rna:mod_weight")
        atac_w_series = mdata.obs.get("atac:mod_weight")

        per_cluster_weights: dict[str, dict[str, Any]] = {}
        if rna_w_series is not None and atac_w_series is not None:
            rna_w = np.asarray(rna_w_series.values, dtype=np.float64)
            atac_w = np.asarray(atac_w_series.values, dtype=np.float64)
            mdata.obs["rna_weight"] = rna_w
            mdata.obs["atac_weight"] = atac_w
            mean_rna_w = float(np.mean(rna_w))
            mean_atac_w = float(np.mean(atac_w))
            median_rna_w = float(np.median(rna_w))
            median_atac_w = float(np.median(atac_w))

            for cl in sorted(set(wnn_labels)):
                mask = np.asarray([lbl == cl for lbl in wnn_labels], dtype=bool)
                cl_rna = float(np.mean(rna_w[mask]))
                cl_atac = float(np.mean(atac_w[mask]))
                dom = (
                    "rna"
                    if cl_rna > cl_atac + 0.05
                    else ("atac" if cl_atac > cl_rna + 0.05 else "balanced")
                )
                per_cluster_weights[str(cl)] = {
                    "n_cells": int(mask.sum()),
                    "mean_rna_weight": round(cl_rna, 4),
                    "mean_atac_weight": round(cl_atac, 4),
                    "dominant_modality": dom,
                }
        else:
            mean_rna_w = None
            mean_atac_w = None
            median_rna_w = None
            median_atac_w = None
            warnings.append(
                "WNN modality weights were not exposed by muon and were not "
                "computed; reported as null (never estimated)."
            )

        if mean_rna_w is not None and mean_rna_w > 0.85:
            warnings.append(
                f"RNA dominates WNN (weight={mean_rna_w:.2f}). "
                f"Verify ATAC data quality (FRiP, TSS enrichment)."
            )
        if mean_atac_w is not None and mean_atac_w > 0.85:
            warnings.append(
                f"ATAC dominates WNN (weight={mean_atac_w:.2f}). "
                f"Verify RNA library complexity and UMI depth."
            )

        # ── 5. Cross-modal cluster agreement & contingency discordance ────
        wnn_vs_rna_ari = adjusted_rand_index(wnn_labels, rna_labels)
        wnn_vs_atac_ari = adjusted_rand_index(wnn_labels, atac_labels)
        rna_vs_atac_ari = adjusted_rand_index(rna_labels, atac_labels)
        wnn_vs_rna_nmi = normalized_mutual_info(wnn_labels, rna_labels)
        wnn_vs_atac_nmi = normalized_mutual_info(wnn_labels, atac_labels)
        rna_vs_atac_nmi = normalized_mutual_info(rna_labels, atac_labels)

        n_discordant, cluster_concordance = _count_discordant_clusters(
            wnn_labels, rna_labels, purity_threshold=0.70
        )

        joint_cluster_sizes = {
            str(cl): int(sum(1 for x in wnn_labels if x == cl))
            for cl in sorted(set(wnn_labels))
        }
        n_joint_clusters = len(joint_cluster_sizes)
        n_rna_clusters = len(set(rna_labels))
        n_atac_clusters = len(set(atac_labels))

        atac_frip, frip_col = _find_measured_frip(
            [atac_adata.obs, global_obs_preserved]
        )
        metrics_not_computed: list[str] = []
        if atac_frip is None:
            metrics_not_computed.append("atac_frip")

        # Stamp HDF5-safe provenance metadata (no None values in .uns).
        mdata.uns["aria_wnn"] = {
            "method": "muon_wnn_l2",
            "input_kind": str(input_kind),
            "rna_counts_source": str(rna_counts_source),
            "atac_counts_source": str(atac_counts_source),
            "seed": int(seed),
            "k": int(k_eff),
            "n_bandwidth_neighbors": int(n_bw),
            "n_multineighbors": int(n_multi),
            "resolution": float(resolution),
            "n_hvg": int(n_hvg_used),
            "n_pcs": int(n_pcs_eff),
            "n_lsi_components_raw": int(n_lsi_eff),
            "n_lsi_components_used": int(len(keep_lsi)),
            "dropped_lsi_components": [int(c) for c in dropped_lsi],
            "depth_corr_cutoff": float(depth_corr_cutoff),
            "n_joint_clusters": int(n_joint_clusters),
            "wnn_vs_rna_ari": round(float(wnn_vs_rna_ari), 4),
            "wnn_vs_atac_ari": round(float(wnn_vs_atac_ari), 4),
        }

        # ── 5b. V48 Slice I3 anti-degradation guardrail ("Lens, not replacement")
        run_guardrail = bool(params.get("run_guardrail", True))
        anti_degradation: dict[str, Any] | None = None
        guardrail_status = "skipped"
        if run_guardrail and n_cells >= 6:
            audit_group_col = params.get("audit_group_col")
            audit_labels: list[str] | None = None
            if audit_group_col:
                for df in (global_obs_preserved, rna_adata.obs, atac_adata.obs):
                    if df is not None and str(audit_group_col) in df.columns:
                        audit_labels = [
                            str(v) for v in df[str(audit_group_col)].tolist()
                        ]
                        break

            wnn_graph_for_diag = (
                mdata.obsp["wnn_distances"]
                if "wnn_distances" in mdata.obsp
                else mdata.obsp["wnn_connectivities"]
            )
            is_conn_diag = "wnn_distances" not in mdata.obsp
            anti_degradation = compute_wnn_anti_degradation_diagnostics(
                rna_embedding=np.asarray(rna_adata.obsm["X_pca_l2"], dtype=np.float64),
                atac_embedding=np.asarray(
                    atac_adata.obsm["X_lsi_l2"], dtype=np.float64
                ),
                wnn_distances=wnn_graph_for_diag,
                wnn_clusters=wnn_labels,
                rna_clusters=rna_labels,
                atac_clusters=atac_labels,
                rna_weights=(
                    rna_w
                    if (rna_w_series is not None and atac_w_series is not None)
                    else None
                ),
                atac_weights=(
                    atac_w
                    if (rna_w_series is not None and atac_w_series is not None)
                    else None
                ),
                audit_labels=audit_labels,
                graph_is_connectivity=is_conn_diag,
                k=int(k_eff),
                seed=int(seed),
            )
            guardrail_status = str(anti_degradation["guardrail_status"])
            mdata.uns["aria_wnn_diagnostics"] = sanitize_for_hdf5_uns(
                anti_degradation
            )
            for issue in anti_degradation.get("issues", []):
                if issue.get("severity") == "degraded":
                    warnings.append(
                        f"WNN guardrail [{issue['code']}]: {issue['detail']}"
                    )

        # ── 6. Atomic .h5mu output write ──────────────────────────────────
        out_path = Path(output_dir) / "wnn_integrated.h5mu"
        tmp_path = Path(output_dir) / f".wnn_integrated.tmp.{os.getpid()}.h5mu"
        try:
            mdata.write_h5mu(str(tmp_path))
            os.replace(str(tmp_path), str(out_path))
        finally:
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

        joint_clusters_map: dict[str, str] = {}
        if n_cells <= 50_000:
            joint_clusters_map = {
                str(bc): str(lbl) for bc, lbl in zip(common_cells, wnn_labels)
            }
        else:
            warnings.append(
                f"Omitted inline joint_clusters map in JSON for {n_cells} cells "
                f"(persisted in mdata.obs['leiden_wnn'])."
            )

        return {
            "status": "success",
            "input_kind": input_kind,
            "seed": int(seed),
            "k": int(k_eff),
            "resolution": float(resolution),
            "n_cells": int(n_cells),
            "n_genes": int(rna_adata.n_vars),
            "n_hvg": int(n_hvg_used),
            "n_pcs": int(n_pcs_eff),
            "n_peaks": int(atac_adata.n_vars),
            "n_lsi_components_raw": int(n_lsi_eff),
            "n_lsi_components_used": int(len(keep_lsi)),
            "dropped_lsi_components": [int(c) for c in dropped_lsi],
            "lsi_depth_correlations": [round(float(r), 4) for r in depth_corr],
            "rna_counts_source": rna_counts_source,
            "atac_counts_source": atac_counts_source,
            "mean_rna_weight": (
                round(float(mean_rna_w), 4) if mean_rna_w is not None else None
            ),
            "mean_atac_weight": (
                round(float(mean_atac_w), 4) if mean_atac_w is not None else None
            ),
            "median_rna_weight": (
                round(float(median_rna_w), 4) if median_rna_w is not None else None
            ),
            "median_atac_weight": (
                round(float(median_atac_w), 4) if median_atac_w is not None else None
            ),
            "n_joint_clusters": int(n_joint_clusters),
            "n_rna_only_clusters": int(n_rna_clusters),
            "n_atac_only_clusters": int(n_atac_clusters),
            "joint_cluster_sizes": joint_cluster_sizes,
            "per_cluster_modality_weights": per_cluster_weights,
            "wnn_vs_rna_ari": round(float(wnn_vs_rna_ari), 4),
            "wnn_vs_atac_ari": round(float(wnn_vs_atac_ari), 4),
            "rna_vs_atac_ari": round(float(rna_vs_atac_ari), 4),
            "wnn_vs_rna_nmi": round(float(wnn_vs_rna_nmi), 4),
            "wnn_vs_atac_nmi": round(float(wnn_vs_atac_nmi), 4),
            "rna_vs_atac_nmi": round(float(rna_vs_atac_nmi), 4),
            "n_discordant_clusters": int(n_discordant),
            "discordance_method": "contingency_modal_purity_lt_0.70",
            "cluster_concordance": cluster_concordance,
            "joint_clusters": joint_clusters_map,
            "atac_frip": atac_frip,
            "atac_frip_source": frip_col,
            "metrics_not_computed": metrics_not_computed,
            "guardrail_status": guardrail_status,
            "anti_degradation": anti_degradation,
            "output_path": str(out_path),
            "warnings": warnings,
        }

    except NotImplementedError as e:
        # P0-8: an unimplemented scaffold step (e.g. raw fragment peak-matrix build)
        # surfaces as an explicit structural blocker, never a fabricated result.
        return {
            "status": "error",
            "error_type": "NotImplemented",
            "details": str(e)[:500],
            "validation_level": "scaffold",
            "warnings": warnings,
        }

    except ImportError as e:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"WNN integration dependencies are required: {e}",
            "warnings": warnings,
        }

    except Exception as e:
        return {
            "status": "error",
            "error_type": "WNNFailed",
            "details": f"{type(e).__name__}: {str(e)[:500]}",
            "warnings": warnings,
        }


def _load_rna(path: str):
    """Load scRNA-seq count matrix from .h5ad, 10x .h5, or 10x MEX directory."""
    try:
        import scanpy as sc

        p = Path(path)
        if p.suffix.lower() == ".h5ad":
            return sc.read_h5ad(str(p))
        if p.suffix.lower() == ".h5":
            return sc.read_10x_h5(str(p))
        if p.is_dir():
            return sc.read_10x_mtx(str(p), var_names="gene_symbols", cache=True)
    except Exception:
        return None
    return None


def _load_atac(files: list, genome: str, organism: str):
    """Load scATAC peak-count matrix (.h5ad) or raise explicit blocker on fragments.

    P0-8 (ADR-002/ADR-025): raw fragment peak-matrix construction
    (`count_fragments_in_peaks`) is not implemented here; ARIA raises an explicit
    `NotImplementedError` rather than fabricating a placeholder peak matrix.
    Precomputed peak-by-cell `.h5ad` matrices are loaded directly and verified
    for raw counts upstream.
    """
    frag_files = [
        str(f) for f in files
        if "fragment" in str(f).lower() or str(f).lower().endswith(".tsv.gz")
    ]
    if frag_files:
        raise NotImplementedError(
            "scATAC peak-matrix construction from raw fragments for WNN is not "
            "implemented (needs snapatac2/episcanpy count_fragments_in_peaks). "
            "Supply a paired .h5mu or peak-count .h5ad; ARIA refuses to "
            "fabricate a peak matrix."
        )

    h5ad_files = [str(f) for f in files if Path(str(f)).suffix.lower() == ".h5ad"]
    if h5ad_files:
        try:
            import scanpy as sc

            return sc.read_h5ad(h5ad_files[0])
        except Exception:
            return None
    return None


if __name__ == "__main__":
    run_script(integration_wnn)
