"""Bulk RNA-seq sample QC (PCA-based outlier detection) and outlier sensitivity
(A7 split of rna_bulk_de.py; bodies verbatim).

Re-exported from aria.scripts.rna_bulk_de."""
from __future__ import annotations
import warnings

from aria.utils.stats import contrast_family_significance
from aria.scripts.rna_bulk.transforms import _run_vst, _select_variable_genes
from aria.scripts.rna_bulk.deseq2 import _run_deseq2
from aria.scripts.rna_bulk.plots import _plot_pca_mds


def _sample_qc(counts, metadata, output_dir: str,
               warnings: list,
               biotype_map: dict | None = None) -> dict:
    """
    Sample-level QC for bulk RNA-seq.

    Two complementary checks:
      (a) Pairwise replicate concordance — Spearman correlation between
          samples within the same biological group. A sample that
          correlates < 0.85 with its replicates is technically suspect.
      (b) PCA-based global outlier detection — samples > 2.5 SD from the
          centroid in PC1-PC2 are flagged. These are dataset-level
          outliers (could be wrong condition assignment, batch effect).

    A sample is reported as "outlier" if it fails BOTH (more conservative).
    Primary DE retains flagged outliers; removal is a sensitivity analysis only.
    """
    try:
        import pandas as pd
        import numpy as np
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        # ── Library size stats ────────────────────────────────────────────
        lib_sizes  = counts.sum(axis=0)
        size_ratio = float(lib_sizes.max() / max(lib_sizes.min(), 1))
        if size_ratio > 10:
            warnings.append(
                f"Library size range: {size_ratio:.1f}× "
                f"(max={int(lib_sizes.max()):,}, "
                f"min={int(lib_sizes.min()):,}). "
                f"DESeq2 normalization handles modest variation; "
                f"consider re-sequencing low-depth samples if >100×."
            )

        # ── (a) Replicate concordance ─────────────────────────────────────
        # Use log2(counts+1) so highly-expressed genes don't dominate.
        # Spearman is rank-based → robust to outlier genes.
        log_counts = np.log2(counts.astype(float) + 1)
        corr_full  = log_counts.corr(method="spearman")

        # Identify the condition column from metadata
        condition_col = None
        for c in metadata.columns:
            if c not in ("sample", "batch", "replicate"):
                condition_col = c
                break

        replicate_outliers = []
        replicate_report   = {}

        if condition_col:
            for group, samples in metadata.groupby(condition_col):
                sample_ids = list(samples.index)
                if len(sample_ids) < 3:
                    # Need ≥3 reps to assess concordance reliably
                    continue
                sub_corr = corr_full.loc[sample_ids, sample_ids]
                # Mean correlation of each sample with the others in its group
                # (exclude self-correlation = 1.0)
                for s in sample_ids:
                    others = [x for x in sample_ids if x != s]
                    mean_r = sub_corr.loc[s, others].mean()
                    replicate_report[s] = round(float(mean_r), 4)
                    if mean_r < 0.85:
                        replicate_outliers.append(s)
                        warnings.append(
                            f"Sample '{s}' has mean Spearman r={mean_r:.3f} "
                            f"with other {group} replicates (<0.85). "
                            f"Possible technical issue or sample swap."
                        )

        # ── (b) VST-based dimensionality reduction (PCA + MDS) ────────────
        # Methodology (v3.8):
        #   • VST instead of log2(raw+1) → homoscedastic, corrects lib size
        #   • Top variable protein_coding genes → remove noise, focus signal
        #   • NO StandardScaler (VST already stabilizes variance)
        # See methodology docstring at top of this file for full justification.
        vst_matrix = _run_vst(counts, metadata, warnings)
        vst_variable, n_pc, n_var = _select_variable_genes(
            vst_matrix, n_top=2000,
            biotype_map=biotype_map,
            warnings=warnings,
        )
        if n_pc > 0:
            warnings.append(
                f"DR input: {n_var} most-variable protein_coding genes "
                f"(from {n_pc} total protein_coding)."
            )
        else:
            warnings.append(
                f"DR input: {n_var} most-variable genes "
                f"(no biotype filtering — GTF annotation lacked gene_biotype)."
            )

        # Run PCA + MDS on same matrix (apples to apples)
        pca_plot, mds_plot, coords, var_exp = _plot_pca_mds(
            vst_variable, metadata, output_dir, warnings,
        )
        if coords is None:
            # DR failed entirely — empty fallback
            coords = np.zeros((counts.shape[1], 2))
            var_exp = [0, 0]

        # 2.5 SD outlier threshold (same logic as before, now on VST PCA)
        pc12        = coords[:, :2]
        mean_pc     = pc12.mean(axis=0)
        std_pc      = pc12.std(axis=0) + 1e-8
        z_scores    = np.abs((pc12 - mean_pc) / std_pc)
        pca_outlier_mask = z_scores.max(axis=1) > 2.5

        pca_outliers = [
            counts.columns[i]
            for i in range(len(pca_outlier_mask))
            if pca_outlier_mask[i]
        ]

        # ── Combine: a sample is an outlier if it fails BOTH checks ──────
        # OR if its replicate correlation is so low (<0.70) that it's
        # almost certainly broken regardless of PCA.
        critical_replicate = [
            s for s in replicate_outliers
            if replicate_report.get(s, 1.0) < 0.70
        ]
        confirmed_outliers = sorted(set(
            critical_replicate +
            [s for s in pca_outliers if s in replicate_outliers]
        ))

        if pca_outliers and not confirmed_outliers:
            warnings.append(
                f"PCA flagged {pca_outliers} as off-centroid (>2.5 SD), "
                f"but their replicate correlation is good. Likely real "
                f"biological variation — keeping in analysis."
            )

        return {
            "n_samples":           int(counts.shape[1]),
            "outliers":            confirmed_outliers,
            "pca_outliers":        pca_outliers,
            "replicate_outliers":  replicate_outliers,
            "replicate_correlations": replicate_report,
            "pca_variance":        var_exp,
            "lib_size_range":      [int(lib_sizes.min()), int(lib_sizes.max())],
            "size_ratio":          round(float(size_ratio), 1),
            "pca_plot":            pca_plot,
            "mds_plot":            mds_plot,
            "vst_matrix":          vst_matrix,     # for downstream heatmaps
            "n_genes_dr":          int(n_var),
            "n_protein_coding":    int(n_pc),
        }

    except Exception as e:
        warnings.append(f"Sample QC failed: {e}")
        try:
            lib_sizes = counts.sum(axis=0)
            lib_range = [int(lib_sizes.min()), int(lib_sizes.max())]
            size_ratio = round(
                float(lib_sizes.max() / max(lib_sizes.min(), 1)), 1
            )
        except Exception:
            lib_range = []
            size_ratio = None
        return {"n_samples": int(counts.shape[1]), "outliers": [],
                "pca_variance": [], "lib_size_range": lib_range,
                "size_ratio": size_ratio, "error": str(e)}


def _prune_outliers_for_design(outliers: list, metadata,
                               design_factor: str,
                               warnings: list) -> list:
    """
    Keep QC from destroying the statistical design.

    Sample-level QC can be noisy in very small or synthetic datasets. Removing
    every low-concordance sample is worse than reporting the warning and then
    running the requested contrast with explicit caveats. Only remove a sample
    when every affected group still has at least two samples afterward.
    """
    if not outliers:
        return []

    outliers = [s for s in outliers if s in metadata.index]
    if not outliers or design_factor not in metadata.columns:
        return outliers

    kept = []
    skipped = []
    for sample in outliers:
        trial = set(kept + [sample])
        remaining = metadata.drop(index=list(trial), errors="ignore")
        group_sizes = remaining[design_factor].value_counts()
        if (group_sizes < 2).any() or len(group_sizes) < 2:
            skipped.append(sample)
        else:
            kept.append(sample)

    if skipped:
        warnings.append(
            "QC flagged samples as potential outliers but kept them because "
            "removal would leave at least one group with <2 replicates: "
            f"{skipped}."
        )

    return kept


def _run_outlier_sensitivity(
    *,
    counts,
    metadata,
    design_factor: str,
    contrasts_in: list,
    flagged_outliers: list,
    removable_outliers: list,
    primary_contrasts: list,
    padj_thr: float,
    lfc_thr: float,
    allow_mock: bool,
    min_reps: int,
    covariates: list,
    lfc_shrink: bool,
    fdr_family: dict,
    warnings: list,
) -> dict:
    """Run P1-5 outlier sensitivity without changing primary DE results."""
    flagged = [str(s) for s in (flagged_outliers or []) if s in metadata.index]
    removable = [str(s) for s in (removable_outliers or []) if s in metadata.index]
    policy = (
        "primary_includes_all_samples; sensitivity_removes_flagged_outliers_when_design_allows"
    )
    summary = {
        "status": "not_applicable",
        "policy": policy,
        "flagged_samples": flagged,
        "removed_samples": removable,
        "contrasts": [],
    }
    if not flagged:
        summary["reason"] = "no_qc_outliers_flagged"
        return summary
    if not removable:
        summary["status"] = "skipped"
        summary["reason"] = "flagged_outliers_protected_by_design"
        return summary

    counts_s = counts.drop(columns=removable, errors="ignore")
    metadata_s = metadata.drop(index=removable, errors="ignore")
    min_samples = max(2, metadata_s.shape[0] // 4)
    keep = (counts_s > 10).sum(axis=1) >= min_samples
    counts_s_filt = counts_s[keep]
    if counts_s_filt.empty:
        summary["status"] = "error"
        summary["reason"] = "no_genes_after_sensitivity_filter"
        warnings.append(
            "Outlier sensitivity skipped: no genes remained after low-count "
            "filtering on the pruned sample set."
        )
        return summary

    primary_by_name = {
        c.get("name"): c for c in primary_contrasts
        if c.get("status") == "success"
    }
    sensitivity_results = []
    sensitivity_family_stats: dict = {}

    for contrast in contrasts_in:
        num = str(contrast.get("numerator", "")).strip()
        den = str(contrast.get("denominator", "")).strip()
        name = contrast.get("name", f"{num} vs {den}")
        primary = primary_by_name.get(name)
        if not primary:
            continue

        available_groups = list(metadata_s[design_factor].unique())
        if num not in available_groups or den not in available_groups:
            sensitivity_results.append({
                "name": name,
                "status": "skipped",
                "reason": "contrast_group_removed_by_outlier_sensitivity",
            })
            continue

        de_result, de_warn = _run_deseq2(
            counts_s_filt,
            metadata_s,
            design_factor,
            num,
            den,
            padj_thr,
            lfc_thr,
            allow_mock=allow_mock,
            min_replicates_per_condition=min_reps,
            covariates=covariates,
            lfc_shrink=lfc_shrink,
        )
        warnings.extend([f"[{name} sensitivity] {w}" for w in de_warn])
        if de_result.get("status") == "error":
            sensitivity_results.append({
                "name": name,
                "status": "error",
                "reason": de_result.get("details", ""),
                "removed_samples": removable,
            })
            continue

        sig_ids = [str(g) for g in (de_result.get("sig_genes", []) or [])]
        sensitivity_results.append({
            "name": name,
            "status": "success",
            "removed_samples": removable,
            "n_genes_tested": int(counts_s_filt.shape[0]),
            "n_significant_sensitivity": int(de_result.get("n_sig", 0)),
            "n_upregulated_sensitivity": int(de_result.get("n_up", 0)),
            "n_downregulated_sensitivity": int(de_result.get("n_down", 0)),
            "sig_gene_ids_sensitivity": sig_ids,
            "fitted_design_formula": de_result.get("fitted_design_formula"),
            "covariates_adjusted": de_result.get("covariates_adjusted", []),
            "covariates_dropped": de_result.get("covariates_dropped", []),
        })

        rdf = de_result.get("results")
        if rdf is not None and len(rdf) > 0:
            sensitivity_family_stats[name] = {
                str(g): {"pvalue": float(row["pvalue"]),
                         "log2fc": float(row["log2FoldChange"])}
                for g, row in rdf.dropna(subset=["pvalue"]).iterrows()
            }

    if fdr_family.get("fdr_family") == "global" and sensitivity_family_stats:
        fam = contrast_family_significance(
            sensitivity_family_stats, padj_max=padj_thr, lfc_min=None
        )
        for sens in sensitivity_results:
            if sens.get("status") != "success":
                continue
            f = fam.get(sens["name"])
            if not f:
                continue
            sig_ids = [str(g) for g in f["sig_genes"]]
            sens["n_significant_sensitivity"] = int(f["n_sig"])
            sens["n_upregulated_sensitivity"] = int(f["n_up"])
            sens["n_downregulated_sensitivity"] = int(f["n_down"])
            sens["sig_gene_ids_sensitivity"] = sig_ids

    robust_values = []
    for sens in sensitivity_results:
        if sens.get("status") != "success":
            continue
        primary = primary_by_name.get(sens["name"], {})
        primary_ids = set(
            str(g) for g in (
                primary.get("all_sig_gene_ids")
                or primary.get("all_sig_genes")
                or []
            )
        )
        sensitivity_ids = set(sens.get("sig_gene_ids_sensitivity") or [])
        primary_only = sorted(primary_ids - sensitivity_ids)
        sensitivity_only = sorted(sensitivity_ids - primary_ids)
        robust = not primary_only and not sensitivity_only
        robust_values.append(robust)
        sens.update({
            "n_significant_primary": int(primary.get("n_significant", 0)),
            "sig_gene_ids_primary": sorted(primary_ids),
            "overlap_n": len(primary_ids & sensitivity_ids),
            "primary_only_n": len(primary_only),
            "sensitivity_only_n": len(sensitivity_only),
            "primary_only_gene_ids": primary_only[:50],
            "sensitivity_only_gene_ids": sensitivity_only[:50],
            "conclusion_robust": robust,
            "interpretation": (
                "Outlier removal did not change the significant gene set."
                if robust else
                "Outlier removal changed the significant gene set; primary "
                "results should be interpreted with this sensitivity caveat."
            ),
        })
        primary["outlier_sensitivity"] = {
            k: sens[k] for k in (
                "status",
                "removed_samples",
                "n_significant_primary",
                "n_significant_sensitivity",
                "overlap_n",
                "primary_only_n",
                "sensitivity_only_n",
                "conclusion_robust",
                "interpretation",
            ) if k in sens
        }

    summary["status"] = "success" if sensitivity_results else "skipped"
    summary["reason"] = (
        "completed" if sensitivity_results else "no_primary_contrast_to_compare"
    )
    summary["contrasts"] = sensitivity_results
    summary["conclusion_robust"] = all(robust_values) if robust_values else None
    return summary


# ── DESeq2 ────────────────────────────────────────────────────────────────────

