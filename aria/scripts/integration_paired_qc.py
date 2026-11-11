"""
ARIA Paired Multiome (.h5mu) Ingestion, Contract & Joint QC Script (V48 I1)
---------------------------------------------------------------------------
Validates a same-cell paired RNA+ATAC `.h5mu` input, computes measured-only
per-cell RNA and ATAC quality-control metrics over deterministically aligned
cell barcodes, filters both modalities jointly over the exact same cells,
optionally joins cell-type annotations from a reference `.h5ad` `/obs` table
via `h5py` without loading the reference expression matrix into memory, and
atomically writes a validated `.h5mu` ready for downstream WNN / MOFA+ /
peak-to-gene integration (`I2`+).

Executed inside `aria-chromatin-env` (`stack="chromatin"`) where `mudata` and
`muon` are locked and verified, while `INTEGRATION_VALIDATION` remains
`scaffold` (`dispatch_enabled=False`) until explicit tier promotion in Slice I8.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from aria.scripts._base import run_script


_VAR_SYMBOL_COLS = ("gene_symbols", "gene_symbol", "gene_name", "feature_name", "symbol")


def _detect_mt_mask(rna_adata: Any) -> tuple[Any, str | None]:
    """Identify mitochondrial genes (`MT-` / `mt-`) without fabricating 0% on Ensembl IDs."""
    import numpy as np

    var_names = [str(g) for g in rna_adata.var_names]
    mask = np.array([g.upper().startswith("MT-") for g in var_names], dtype=bool)
    if int(np.sum(mask)) > 0:
        return mask, "var_names"

    var_df = getattr(rna_adata, "var", None)
    if var_df is not None:
        for col in _VAR_SYMBOL_COLS:
            if col in var_df.columns:
                vals = [str(v) for v in var_df[col].values]
                col_mask = np.array(
                    [v.upper().startswith("MT-") for v in vals],
                    dtype=bool,
                )
                if int(np.sum(col_mask)) > 0:
                    return col_mask, f"var['{col}']"

    return np.zeros(len(var_names), dtype=bool), None


def _row_totals_and_detected(matrix: Any) -> tuple[Any, Any]:
    """Compute per-cell total counts and detected features without densifying sparse CSR/CSC."""
    import numpy as np
    try:
        from scipy import sparse as sp
        is_sparse = sp.issparse(matrix)
    except ImportError:
        is_sparse = False

    if is_sparse:
        csr = matrix.tocsr()
        totals = np.asarray(csr.sum(axis=1)).ravel().astype(np.float64)
        detected = np.asarray((csr > 0).sum(axis=1)).ravel().astype(np.int64)
        return totals, detected

    arr = np.asarray(matrix, dtype=np.float64)
    totals = np.sum(arr, axis=1)
    detected = np.sum(arr > 0, axis=1).astype(np.int64)
    return totals, detected


def _subset_matrix_rows(matrix: Any, row_indices: list[int]) -> Any:
    import numpy as np
    idx = np.asarray(row_indices, dtype=np.int64)
    return matrix[idx, :]


def integration_paired_qc(params: dict[str, Any]) -> dict[str, Any]:
    """Run paired `.h5mu` contract verification, joint RNA+ATAC QC, and optional `/obs` label join."""
    import numpy as np
    from aria.utils.mudata_io import (
        _import_mudata,
        join_obs_labels_from_h5ad,
        read_h5mu_paired,
    )

    h5mu_path = str(params.get("h5mu_path") or params.get("data_path") or "").strip()
    if not h5mu_path:
        return {
            "status": "error",
            "error_type": "MissingInput",
            "details": "integration_paired_qc requires 'h5mu_path' (or 'data_path').",
        }

    require_raw_counts = bool(params.get("require_raw_counts", True))
    paired = read_h5mu_paired(h5mu_path, require_raw_counts=require_raw_counts)
    if paired.get("status") != "success":
        if paired.get("error_type") == "MuDataToolingMissing":
            return {
                "status": "error",
                "error_type": "MissingDependency",
                "details": paired["details"],
            }
        return paired

    min_shared_frac = float(params.get("min_shared_barcode_fraction", 0.50))
    allow_partial = bool(params.get("allow_partial_overlap", False))
    min_observed_frac = min(
        float(paired["shared_fraction_rna"]),
        float(paired["shared_fraction_atac"]),
    )
    if not allow_partial and min_observed_frac < min_shared_frac:
        return {
            "status": "error",
            "error_type": "LowBarcodeOverlap",
            "details": (
                f"Shared cell-barcode fraction ({min_observed_frac:.3f}) is "
                f"below min_shared_barcode_fraction={min_shared_frac:.3f} "
                f"(n_rna={paired['n_rna_cells']}, n_atac={paired['n_atac_cells']}, "
                f"n_shared={paired['n_shared_cells']}). Pass "
                f"allow_partial_overlap=True only if partial pairing is intentional."
            ),
        }

    mdata_in = paired["mdata"]
    rna_mod = str(paired["rna_modality"])
    atac_mod = str(paired["atac_modality"])
    rna_adata = paired["rna_adata"]
    atac_adata = paired["atac_adata"]
    shared_barcodes: list[str] = list(paired["shared_barcodes"])

    rna_pos = {str(bc): idx for idx, bc in enumerate(rna_adata.obs_names)}
    atac_pos = {str(bc): idx for idx, bc in enumerate(atac_adata.obs_names)}
    rna_rows = [rna_pos[bc] for bc in shared_barcodes]
    atac_rows = [atac_pos[bc] for bc in shared_barcodes]

    rna_mat = _subset_matrix_rows(paired["rna_counts_matrix"], rna_rows)
    atac_mat = _subset_matrix_rows(paired["atac_counts_matrix"], atac_rows)

    # ── Measured-only RNA QC ──────────────────────────────────────────────
    rna_totals, rna_n_genes = _row_totals_and_detected(rna_mat)
    mt_mask, mt_source = _detect_mt_mask(rna_adata)
    n_mt_genes = int(np.sum(mt_mask))
    mt_measured = n_mt_genes > 0

    warnings: list[str] = []
    metrics_not_computed: list[str] = [
        "tss_enrichment (fragment-level metric; requires fragments.tsv.gz + TSS reference)",
        "frip (fragment-level metric; requires fragments.tsv.gz + called peaks)",
    ]

    if mt_measured:
        try:
            from scipy import sparse as sp
            if sp.issparse(rna_mat):
                mt_totals = np.asarray(rna_mat.tocsr()[:, mt_mask].sum(axis=1)).ravel().astype(np.float64)
            else:
                mt_totals = np.sum(np.asarray(rna_mat, dtype=np.float64)[:, mt_mask], axis=1)
        except ImportError:
            mt_totals = np.sum(np.asarray(rna_mat, dtype=np.float64)[:, mt_mask], axis=1)
        pct_counts_mt = np.divide(
            100.0 * mt_totals,
            rna_totals,
            out=np.zeros_like(rna_totals, dtype=np.float64),
            where=rna_totals > 0,
        )
    else:
        pct_counts_mt = None
        metrics_not_computed.append(
            "pct_counts_mt (no MT-/mt- mitochondrial gene symbols detected in RNA var)"
        )
        warnings.append(
            "No MT-/mt- mitochondrial gene symbols detected in RNA modality; "
            "pct_counts_mt left uncomputed (None) rather than fabricated as 0%."
        )

    # ── Measured-only ATAC QC ─────────────────────────────────────────────
    atac_totals, atac_n_peaks = _row_totals_and_detected(atac_mat)
    if paired["atac_is_binary"]:
        warnings.append(
            "ATAC matrix is binary (0/1); total_peak_counts equals n_peaks_by_counts."
        )

    atac_frag_col = paired.get("atac_fragment_count_col")
    atac_mito_col = paired.get("atac_mito_col")
    if atac_frag_col is None:
        metrics_not_computed.append(
            "atac_fragments_per_cell (no per-cell fragment obs column in ATAC modality)"
        )
    if atac_mito_col is None:
        metrics_not_computed.append(
            "atac_mito_fraction (no per-cell mitochondrial obs column in ATAC modality)"
        )

    # ── Explicit Joint QC Thresholds ──────────────────────────────────────
    min_genes = int(params.get("min_genes", 200))
    min_rna_counts = int(params.get("min_rna_counts", 500))
    max_pct_mt = float(params.get("max_pct_mt", 20.0))
    min_peaks = int(params.get("min_peaks", 200))
    min_atac_counts = int(params.get("min_atac_counts", 500))
    max_atac_mito_frac = (
        None if params.get("max_atac_mito_frac") is None and atac_mito_col is None
        else float(params.get("max_atac_mito_frac", 0.10))
    )

    pass_rna = (rna_totals >= min_rna_counts) & (rna_n_genes >= min_genes)
    if mt_measured and pct_counts_mt is not None:
        pass_rna = pass_rna & (pct_counts_mt <= max_pct_mt)

    pass_atac = (atac_totals >= min_atac_counts) & (atac_n_peaks >= min_peaks)
    atac_mito_scale = None
    if atac_mito_col is not None and max_atac_mito_frac is not None:
        raw_atac_mito = np.asarray(
            atac_adata.obs.iloc[atac_rows][atac_mito_col].values,
            dtype=np.float64,
        )
        finite_mito = raw_atac_mito[np.isfinite(raw_atac_mito)]
        if finite_mito.size > 0:
            if float(np.max(finite_mito)) > 1.0:
                atac_mito_frac_01 = raw_atac_mito / 100.0
                atac_mito_scale = "percent_converted_to_fraction"
            else:
                atac_mito_frac_01 = raw_atac_mito
                atac_mito_scale = "fraction"
            pass_atac = pass_atac & (
                ~np.isfinite(atac_mito_frac_01) | (atac_mito_frac_01 <= max_atac_mito_frac)
            )

    pass_joint = pass_rna & pass_atac
    n_shared = len(shared_barcodes)
    n_pass_rna = int(np.sum(pass_rna))
    n_pass_atac = int(np.sum(pass_atac))
    n_after = int(np.sum(pass_joint))

    if n_after == 0:
        return {
            "status": "error",
            "error_type": "ZeroCellsAfterQC",
            "details": (
                f"Joint RNA+ATAC QC retained 0 / {n_shared} shared cells "
                f"(n_pass_rna={n_pass_rna}, n_pass_atac={n_pass_atac}; "
                f"thresholds: min_genes={min_genes}, min_rna_counts={min_rna_counts}, "
                f"max_pct_mt={max_pct_mt if mt_measured else None}, "
                f"min_peaks={min_peaks}, min_atac_counts={min_atac_counts})."
            ),
        }

    kept_indices = [i for i, keep in enumerate(pass_joint) if bool(keep)]
    kept_barcodes = [shared_barcodes[i] for i in kept_indices]

    # ── Optional Deterministic Reference `/obs` Barcode Join ──────────────
    reference_h5ad = params.get("reference_h5ad")
    sample_id_raw = params.get("sample_id")
    sample_id = str(sample_id_raw).strip() if sample_id_raw is not None and str(sample_id_raw).strip() else None
    label_join_summary: dict[str, Any] = {
        "requested": bool(reference_h5ad),
        "ran": False,
        "mechanism": "obs_barcode_join",
    }
    transferred_labels: list[str] | None = None
    transferred_mask: list[bool] | None = None
    target_label_col: str | None = None

    if reference_h5ad:
        ref_label_col = str(params.get("reference_label_col") or "subclass").strip()
        target_label_col = str(params.get("target_label_col") or ref_label_col).strip()
        ref_sample_col = params.get("reference_sample_col")
        ref_sample_col = str(ref_sample_col).strip() if ref_sample_col else None
        require_ref_match = bool(params.get("require_reference_match", True))
        min_warn_match_frac = float(params.get("min_warn_reference_match_fraction", 0.50))

        join_res = join_obs_labels_from_h5ad(
            kept_barcodes,
            str(reference_h5ad),
            label_col=ref_label_col,
            sample_id=sample_id,
            reference_sample_col=ref_sample_col,
        )
        if join_res.get("status") != "success":
            return join_res

        if require_ref_match and int(join_res["n_matched"]) == 0:
            return {
                "status": "error",
                "error_type": "NoReferenceBarcodesMatched",
                "details": (
                    f"Zero post-QC cells matched reference .h5ad '{reference_h5ad}' "
                    f"(label_col={ref_label_col!r}, sample_id={sample_id!r}, "
                    f"sample_query={join_res.get('sample_query_barcodes')}, "
                    f"sample_ref={join_res.get('sample_reference_barcodes')}). "
                    f"Pass sample_id / reference_sample_col if barcodes are donor-prefixed, "
                    f"or require_reference_match=False to proceed with unassigned cells."
                ),
            }

        if float(join_res["match_fraction"]) < min_warn_match_frac:
            warnings.append(
                f"Reference /obs barcode join matched {join_res['n_matched']}/{join_res['n_query']} "
                f"({join_res['match_fraction']:.1%}) post-QC cells (< {min_warn_match_frac:.0%}); "
                f"unmatched cells are marked 'unassigned'."
            )

        transferred_labels = list(join_res["labels"])
        transferred_mask = list(join_res["matched_mask"])
        label_join_summary = {
            "requested": True,
            "ran": True,
            "mechanism": "obs_barcode_join",
            "reference_h5ad": str(reference_h5ad),
            "reference_label_col": ref_label_col,
            "target_label_col": target_label_col,
            "sample_id": sample_id,
            "reference_sample_col": ref_sample_col,
            "match_strategy": join_res["match_strategy"],
            "n_query": int(join_res["n_query"]),
            "n_matched": int(join_res["n_matched"]),
            "n_unmatched": int(join_res["n_unmatched"]),
            "match_fraction": round(float(join_res["match_fraction"]), 4),
            "label_counts": dict(join_res["label_counts"]),
        }

    # ── Build Filtered Paired MuData & Persist Atomically ─────────────────
    md = _import_mudata()
    assert md is not None

    rna_kept_rows = [rna_pos[bc] for bc in kept_barcodes]
    atac_kept_rows = [atac_pos[bc] for bc in kept_barcodes]
    rna_filt = rna_adata[rna_kept_rows, :].copy()
    atac_filt = atac_adata[atac_kept_rows, :].copy()

    kept_idx_arr = np.asarray(kept_indices, dtype=np.int64)
    rna_filt.obs["total_counts"] = rna_totals[kept_idx_arr]
    rna_filt.obs["n_genes_by_counts"] = rna_n_genes[kept_idx_arr]
    if mt_measured and pct_counts_mt is not None:
        rna_filt.obs["pct_counts_mt"] = pct_counts_mt[kept_idx_arr]

    atac_filt.obs["total_peak_counts"] = atac_totals[kept_idx_arr]
    atac_filt.obs["n_peaks_by_counts"] = atac_n_peaks[kept_idx_arr]

    if sample_id is not None:
        rna_filt.obs["sample_id"] = sample_id
        atac_filt.obs["sample_id"] = sample_id

    if transferred_labels is not None and transferred_mask is not None and target_label_col:
        rna_filt.obs[target_label_col] = transferred_labels
        rna_filt.obs[f"{target_label_col}_transferred"] = transferred_mask
        atac_filt.obs[target_label_col] = transferred_labels
        atac_filt.obs[f"{target_label_col}_transferred"] = transferred_mask

    mdata_out = md.MuData({rna_mod: rna_filt, atac_mod: atac_filt})
    mdata_out.update()

    # Preserve shared top-level mdata.obs columns aligned to kept_barcodes.
    if getattr(mdata_in, "obs", None) is not None and len(mdata_in.obs.columns) > 0:
        shared_in_obs = mdata_in.obs.reindex(kept_barcodes)
        for col in shared_in_obs.columns:
            if col not in mdata_out.obs.columns:
                mdata_out.obs[col] = shared_in_obs[col].values

    mdata_out.obs["rna_total_counts"] = rna_totals[kept_idx_arr]
    mdata_out.obs["rna_n_genes_by_counts"] = rna_n_genes[kept_idx_arr]
    if mt_measured and pct_counts_mt is not None:
        mdata_out.obs["rna_pct_counts_mt"] = pct_counts_mt[kept_idx_arr]
    mdata_out.obs["atac_total_peak_counts"] = atac_totals[kept_idx_arr]
    mdata_out.obs["atac_n_peaks_by_counts"] = atac_n_peaks[kept_idx_arr]
    mdata_out.obs["pass_joint_qc"] = True
    if sample_id is not None:
        mdata_out.obs["sample_id"] = sample_id
    if transferred_labels is not None and transferred_mask is not None and target_label_col:
        mdata_out.obs[target_label_col] = transferred_labels
        mdata_out.obs[f"{target_label_col}_transferred"] = transferred_mask

    in_path = Path(h5mu_path)
    out_dir = Path(params.get("output_dir") or in_path.parent)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem_tag = sample_id or in_path.stem
    out_name = str(params.get("output_filename") or f"paired_qc_{stem_tag}.h5mu")
    if not out_name.lower().endswith(".h5mu"):
        out_name = f"{out_name}.h5mu"
    out_path = out_dir / out_name
    tmp_path = out_dir / f".{out_Path_safe(out_name)}.tmp.{os.getpid()}.h5mu"

    try:
        mdata_out.write_h5mu(str(tmp_path))
        os.replace(str(tmp_path), str(out_path))
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass

    pct_removed = round(100.0 * float(n_shared - n_after) / float(n_shared), 2)

    return {
        "status": "success",
        "input_path": str(in_path),
        "output_path": str(out_path),
        "rna_modality": rna_mod,
        "atac_modality": atac_mod,
        "rna_counts_source": paired["rna_counts_source"],
        "atac_counts_source": paired["atac_counts_source"],
        "rna_is_raw_counts": paired["rna_is_raw_counts"],
        "atac_is_raw_counts": paired["atac_is_raw_counts"],
        "atac_is_binary": paired["atac_is_binary"],
        "exact_same_input_order": paired["exact_same_order"],
        "n_rna_cells_input": int(paired["n_rna_cells"]),
        "n_atac_cells_input": int(paired["n_atac_cells"]),
        "n_cells_before": int(n_shared),
        "n_pass_rna": n_pass_rna,
        "n_pass_atac": n_pass_atac,
        "n_cells_after": n_after,
        "pct_removed": pct_removed,
        "n_genes": int(rna_filt.n_vars),
        "n_peaks": int(atac_filt.n_vars),
        "thresholds_used": {
            "min_genes": min_genes,
            "min_rna_counts": min_rna_counts,
            "max_pct_mt": max_pct_mt if mt_measured else None,
            "min_peaks": min_peaks,
            "min_atac_counts": min_atac_counts,
            "max_atac_mito_frac": max_atac_mito_frac if atac_mito_col else None,
            "atac_mito_scale": atac_mito_scale,
            "min_shared_barcode_fraction": min_shared_frac,
        },
        "rna_qc_summary": {
            "mt_gene_source": mt_source,
            "n_mt_genes": n_mt_genes,
            "median_total_counts_before": round(float(np.median(rna_totals)), 2),
            "median_total_counts_after": round(float(np.median(rna_totals[kept_idx_arr])), 2),
            "median_genes_before": round(float(np.median(rna_n_genes)), 2),
            "median_genes_after": round(float(np.median(rna_n_genes[kept_idx_arr])), 2),
            "median_pct_mt_before": (
                round(float(np.median(pct_counts_mt)), 3)
                if mt_measured and pct_counts_mt is not None else None
            ),
            "median_pct_mt_after": (
                round(float(np.median(pct_counts_mt[kept_idx_arr])), 3)
                if mt_measured and pct_counts_mt is not None else None
            ),
        },
        "atac_qc_summary": {
            "median_peak_counts_before": round(float(np.median(atac_totals)), 2),
            "median_peak_counts_after": round(float(np.median(atac_totals[kept_idx_arr])), 2),
            "median_peaks_before": round(float(np.median(atac_n_peaks)), 2),
            "median_peaks_after": round(float(np.median(atac_n_peaks[kept_idx_arr])), 2),
            "frip": None,
            "tss_enrichment": None,
            "fragment_count_col": atac_frag_col,
            "median_fragments_per_cell": paired.get("median_atac_fragments_per_cell"),
            "mito_obs_col": atac_mito_col,
            "median_mito_obs": paired.get("median_atac_mito_fraction"),
        },
        "label_join": label_join_summary,
        "metrics_not_computed": metrics_not_computed,
        "warnings": warnings,
    }


def out_Path_safe(filename: str) -> str:
    return filename[:-5] if filename.lower().endswith(".h5mu") else filename


if __name__ == "__main__":
    run_script(integration_paired_qc)
