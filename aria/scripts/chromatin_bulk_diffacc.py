"""Bulk ATAC differential accessibility over peak-count matrices.

Runs replicate-gated DESeq2 on a peak x sample count matrix. This script requires
explicit condition, replicate, and comparison metadata; it never infers biology
from filenames.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from aria.scripts._base import mocks_allowed, run_script


def _error(error_type: str, details: str, **extra) -> dict:
    out = {"status": "error", "error_type": error_type, "details": details}
    out.update(extra)
    return out


def _skipped(reason: str, **extra) -> dict:
    out = {
        "status": "success",
        "ran": False,
        "reason": reason,
        "data_type": "bulk_ATAC",
        "analysis": "differential_accessibility",
        "validation_level": "beta",
        "comparisons": [],
        "warnings": [],
    }
    out.update(extra)
    return out


def _normalise_comparisons(raw) -> list[dict[str, str]]:
    if raw in (None, "", []):
        return []
    if isinstance(raw, str):
        text = raw.strip()
        for sep in ("_vs_", " vs ", "_VS_", " VS "):
            if sep in text:
                left, right = text.split(sep, 1)
                if left.strip() and right.strip():
                    return [{"test": left.strip(), "reference": right.strip()}]
        return []
    if isinstance(raw, dict):
        test = raw.get("test") or raw.get("numerator")
        ref = raw.get("reference") or raw.get("ref") or raw.get("denominator")
        return ([{"test": str(test), "reference": str(ref)}]
                if test and ref else [])
    out = []
    if isinstance(raw, list):
        for item in raw:
            if isinstance(item, dict):
                out.extend(_normalise_comparisons(item))
            elif isinstance(item, (list, tuple)) and len(item) >= 2:
                out.append({"test": str(item[0]), "reference": str(item[1])})
            elif isinstance(item, str):
                out.extend(_normalise_comparisons(item))
    return out


def _load_counts_matrix(path: str | Path):
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    if df.empty or df.shape[1] < 2:
        raise ValueError("counts matrix must contain peak_id plus sample columns.")
    peak_col = df.columns[0]
    if str(peak_col) not in {"peak_id", "peak", "feature_id"}:
        raise ValueError("first counts matrix column must be peak_id/peak.")
    if df[peak_col].isna().any() or df[peak_col].astype(str).duplicated().any():
        raise ValueError("peak IDs must be non-empty and unique.")
    counts = df.set_index(peak_col)
    for col in counts.columns:
        counts[col] = pd.to_numeric(counts[col], errors="raise")
    if (counts < 0).any().any():
        raise ValueError("counts matrix contains negative values.")
    return counts.astype(int)


def _load_sample_metadata(path: str | Path):
    import pandas as pd

    meta = pd.read_csv(path, sep="\t")
    if "sample_id" not in meta.columns:
        raise ValueError("sample metadata must contain a sample_id column.")
    if meta["sample_id"].isna().any():
        raise ValueError("sample metadata contains empty sample_id values.")
    meta["sample_id"] = meta["sample_id"].astype(str)
    if meta["sample_id"].duplicated().any():
        raise ValueError("sample metadata sample_id values must be unique.")
    return meta.set_index("sample_id")


def _prepare_replicate_matrix(counts, metadata, condition_col: str,
                              replicate_col: str, covariates: list[str]):
    import pandas as pd

    samples = list(counts.columns)
    missing_meta = [s for s in samples if s not in metadata.index]
    if missing_meta:
        raise ValueError(
            "sample metadata is missing count-matrix sample(s): "
            f"{missing_meta}"
        )
    meta = metadata.loc[samples].copy()
    for col in (condition_col, replicate_col):
        if col not in meta.columns:
            raise ValueError(f"required metadata column missing: {col}")
        if meta[col].isna().any() or (meta[col].astype(str).str.strip() == "").any():
            raise ValueError(f"metadata column {col} contains empty values.")
        meta[col] = meta[col].astype(str)

    requested_covariates = list(dict.fromkeys(c for c in covariates if c))
    covariates_dropped: list[dict[str, Any]] = []
    for cov in requested_covariates:
        if cov not in meta.columns:
            covariates_dropped.append({
                "covariate": cov,
                "reason": "not present in the sample metadata",
            })
    present_covariates = [c for c in requested_covariates if c in meta.columns]
    meta["_replicate_sample"] = (
        meta[condition_col].astype(str) + "::" + meta[replicate_col].astype(str)
    )

    agg_counts = []
    agg_meta_rows = []
    bad_covariate_reps: dict[str, dict[str, list[str]]] = {
        cov: {"missing": [], "non_constant": []}
        for cov in present_covariates
    }
    for rep_sample, rows in meta.groupby("_replicate_sample", sort=False):
        sample_ids = list(rows.index)
        agg_counts.append(counts[sample_ids].sum(axis=1).rename(rep_sample))
        row = {
            "_sample": rep_sample,
            condition_col: str(rows[condition_col].iloc[0]),
            replicate_col: str(rows[replicate_col].iloc[0]),
            "n_technical_samples": int(len(sample_ids)),
        }
        for cov in present_covariates:
            vals = [
                str(v).strip()
                for v in rows[cov].dropna().astype(str).unique()
                if str(v).strip()
            ]
            if len(vals) == 1:
                row[cov] = vals[0]
            elif len(vals) == 0:
                bad_covariate_reps[cov]["missing"].append(str(rep_sample))
            else:
                bad_covariate_reps[cov]["non_constant"].append(str(rep_sample))
        agg_meta_rows.append(row)

    agg_counts_df = pd.concat(agg_counts, axis=1)
    agg_meta = pd.DataFrame(agg_meta_rows).set_index("_sample")

    dropped_covs = {d["covariate"] for d in covariates_dropped}
    for cov, problems in bad_covariate_reps.items():
        if problems["non_constant"]:
            covariates_dropped.append({
                "covariate": cov,
                "reason": (
                    "not constant within biological replicate after "
                    "technical-sample aggregation"
                ),
                "affected_replicates": problems["non_constant"],
                "n_affected_replicates": len(problems["non_constant"]),
            })
            dropped_covs.add(cov)
        elif problems["missing"]:
            covariates_dropped.append({
                "covariate": cov,
                "reason": (
                    "missing within biological replicate after "
                    "technical-sample aggregation"
                ),
                "affected_replicates": problems["missing"],
                "n_affected_replicates": len(problems["missing"]),
            })
            dropped_covs.add(cov)

    usable_covariates = [
        c for c in present_covariates
        if c not in dropped_covs
        and c in agg_meta.columns
        and not agg_meta[c].isna().any()
    ]
    return agg_counts_df, agg_meta, usable_covariates, covariates_dropped


def _write_replicate_counts(path: Path, counts) -> None:
    out = counts.copy()
    out.insert(0, "peak_id", out.index)
    out.to_csv(path, sep="\t", index=False)


def _write_comparison_table(path: Path, results_df, sig_df, comp_key: str) -> None:
    table = results_df.copy()
    table.insert(0, "peak", table.index)
    table["comparison"] = comp_key
    table["significant"] = table.index.isin(sig_df.index)
    table.to_csv(path, index=False)


def _write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = sorted({k for row in rows for k in row})
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def chromatin_bulk_diffacc(params: dict) -> dict:
    if params.get("data_type", "bulk_ATAC") != "bulk_ATAC":
        return _error(
            "UnsupportedDataType",
            "chromatin_bulk_diffacc currently supports bulk_ATAC only.",
            data_type=params.get("data_type"),
        )

    counts_path = params.get("counts_matrix_path")
    metadata_path = params.get("sample_metadata_path")
    if not counts_path or not Path(str(counts_path)).exists():
        return _error("FileNotFound", f"counts_matrix_path not found: {counts_path}")
    if not metadata_path or not Path(str(metadata_path)).exists():
        return _error("FileNotFound", f"sample_metadata_path not found: {metadata_path}")

    condition_col = str(params.get("condition_col") or "condition")
    replicate_col = str(params.get("replicate_col") or "replicate")
    comparisons = _normalise_comparisons(
        params.get("comparisons") or params.get("comparison"))
    if not comparisons:
        return _skipped(
            "no explicit comparison supplied; ARIA does not invent a "
            "bulk ATAC reference contrast",
        )

    from aria.utils.thresholds import AnalysisThresholds
    thr = AnalysisThresholds.from_exp_context(
        params.get("exp_context") or {},
        modality="bulk_ATAC",
        log2fc_default=0.5,
        min_replicates=int(params.get("min_replicates_per_condition", 3)),
    )
    padj_max = float(params.get("padj_max", thr.padj))
    lfc_min = float(params.get("lfc_min", thr.log2fc))
    min_reps = int(params.get("min_replicates_per_condition",
                              thr.min_replicates))
    top_n = int(params.get("top_n", 20))
    covariates = [str(c) for c in (params.get("covariates") or []) if str(c)]
    required_covariates = {
        str(c) for c in (
            params.get("required_covariates")
            or params.get("essential_covariates")
            or []
        )
        if str(c)
    }
    n_cpus = params.get("n_cpus")
    allow_mock = mocks_allowed(params)
    warnings: list[str] = []

    try:
        counts = _load_counts_matrix(counts_path)
        metadata = _load_sample_metadata(metadata_path)
        if replicate_col not in metadata.columns and replicate_col == "replicate":
            if "replicate_id" in metadata.columns:
                replicate_col = "replicate_id"
        rep_counts, rep_meta, usable_covariates, covariates_dropped = (
            _prepare_replicate_matrix(
                counts, metadata, condition_col, replicate_col, covariates)
        )
    except Exception as exc:
        return _skipped(str(exc))

    for dropped in covariates_dropped:
        warnings.append(
            f"Covariate '{dropped.get('covariate')}' was not adjusted for in "
            f"bulk ATAC DA: {dropped.get('reason')}."
        )

    dropped_required = [
        d for d in covariates_dropped
        if d.get("covariate") in required_covariates
    ]
    if dropped_required:
        return _skipped(
            "required covariate(s) cannot be represented after "
            "technical-sample aggregation",
            covariates_requested=covariates,
            covariates_adjusted=usable_covariates,
            covariates_dropped=covariates_dropped,
            warnings=warnings,
        )

    nonzero = rep_counts.sum(axis=1) > 0
    n_dropped_zero = int((~nonzero).sum())
    if n_dropped_zero:
        warnings.append(f"Dropped {n_dropped_zero} all-zero peak(s) before DESeq2.")
    rep_counts = rep_counts.loc[nonzero]
    if rep_counts.empty:
        return _skipped("all peaks had zero counts after replicate aggregation")

    output_dir = Path(str(params.get("output_dir") or
                          Path(str(counts_path)).parent / "bulk_atac_da"))
    output_dir.mkdir(parents=True, exist_ok=True)
    replicate_counts_path = output_dir / "bulk_atac_replicate_peak_counts.tsv"
    _write_replicate_counts(replicate_counts_path, rep_counts)

    from aria.scripts.chromatin_diffacc import _gate_da_by_convergence
    from aria.scripts.rna_bulk_de import _run_deseq2

    comparison_results = []
    summary_rows = []
    n_success = 0
    n_sig_total = 0
    levels = set(rep_meta[condition_col].astype(str))

    for comp in comparisons:
        test = str(comp["test"])
        ref = str(comp["reference"])
        comp_key = f"{test}_vs_{ref}"
        if test not in levels or ref not in levels:
            row = {
                "test": test,
                "reference": ref,
                "status": "skipped",
                "reason": f"contrast levels not present in metadata: {sorted(levels)}",
            }
            comparison_results.append(row)
            summary_rows.append(row)
            continue

        de, de_warnings = _run_deseq2(
            rep_counts,
            rep_meta,
            condition_col,
            test,
            ref,
            padj_max,
            lfc_min,
            allow_mock=allow_mock,
            min_replicates_per_condition=min_reps,
            covariates=usable_covariates,
            n_cpus=int(n_cpus) if n_cpus else None,
            expose_convergence=True,
        )
        warnings.extend(de_warnings)
        if de.get("status") != "success":
            row = {
                "test": test,
                "reference": ref,
                "status": "error",
                "error_type": de.get("error_type"),
                "reason": de.get("details", "")[:300],
            }
            comparison_results.append(row)
            summary_rows.append(row)
            continue

        results_df = de["results"]
        sig_df, conv_info = _gate_da_by_convergence(results_df, padj_max)
        if conv_info["n_nonconverged_excluded"]:
            warnings.append(
                f"[{comp_key}] excluded "
                f"{conv_info['n_nonconverged_excluded']} non-converged peak(s) "
                "from the DA significant set."
            )
        n_up = int((sig_df["log2FoldChange"] > 0).sum())
        n_down = int((sig_df["log2FoldChange"] < 0).sum())
        n_sig = int(len(sig_df))
        n_sig_total += n_sig
        n_success += 1

        full_csv = output_dir / f"bulk_atac_da_{comp_key}_full.csv"
        _write_comparison_table(full_csv, results_df, sig_df, comp_key)
        top_peaks = []
        for peak, row_df in sig_df.head(top_n).iterrows():
            top_peaks.append({
                "peak": str(peak),
                "log2fc": round(float(row_df["log2FoldChange"]), 3),
                "padj": round(float(row_df["padj"]), 6),
            })
        result_row = {
            "test": test,
            "reference": ref,
            "status": "success",
            "n_sig": n_sig,
            "n_up": n_up,
            "n_down": n_down,
            "n_replicates": de.get("n_replicates"),
            "top_peaks": top_peaks,
            "full_results_csv": str(full_csv),
            "low_power_warning": de.get("low_power_warning", False),
            "low_power_reason": de.get("low_power_reason"),
            "convergence_available": conv_info["convergence_available"],
            "n_nonconverged_excluded": conv_info["n_nonconverged_excluded"],
            "lfc_shrinkage": de.get("lfc_shrinkage"),
            "lfc_threshold_test": de.get("lfc_threshold_test"),
            "fitted_design_formula": de.get("fitted_design_formula"),
            "covariates_adjusted": de.get(
                "covariates_adjusted", usable_covariates),
            "covariates_dropped": covariates_dropped + (
                de.get("covariates_dropped") or []),
            "power_estimate_at_lfc_min": de.get("power_estimate_at_lfc_min"),
        }
        comparison_results.append(result_row)
        summary_rows.append({
            "test": test,
            "reference": ref,
            "status": "success",
            "n_sig": n_sig,
            "n_up": n_up,
            "n_down": n_down,
            "full_results_csv": str(full_csv),
        })

    summary_csv = output_dir / "bulk_atac_da_summary.csv"
    _write_summary(summary_csv, summary_rows)

    if n_success == 0:
        return _skipped(
            "no bulk ATAC DA comparison completed successfully",
            comparisons=comparison_results,
            warnings=warnings,
            output_csv=str(summary_csv),
            replicate_counts_path=str(replicate_counts_path),
            padj_max=padj_max,
            lfc_min=lfc_min,
            covariates_requested=covariates,
            covariates_adjusted=usable_covariates,
            covariates_dropped=covariates_dropped,
        )

    return {
        "status": "success",
        "ran": True,
        "reason": None,
        "data_type": "bulk_ATAC",
        "analysis": "differential_accessibility",
        "validation_level": "beta",
        "method": "replicate-level DESeq2 over bulk ATAC peak counts",
        "counts_matrix_path": str(counts_path),
        "sample_metadata_path": str(metadata_path),
        "replicate_counts_path": str(replicate_counts_path),
        "output_csv": str(summary_csv),
        "condition_col": condition_col,
        "replicate_col": replicate_col,
        "covariates_requested": covariates,
        "covariates_adjusted": usable_covariates,
        "covariates_dropped": covariates_dropped,
        "padj_max": padj_max,
        "lfc_min": lfc_min,
        "min_replicates_per_condition": min_reps,
        "n_peaks_tested": int(rep_counts.shape[0]),
        "n_replicate_samples": int(rep_counts.shape[1]),
        "n_comparisons": len(comparison_results),
        "n_comparisons_success": n_success,
        "n_sig_total": n_sig_total,
        "comparisons": comparison_results,
        "warnings": warnings,
    }


if __name__ == "__main__":
    run_script(chromatin_bulk_diffacc)
