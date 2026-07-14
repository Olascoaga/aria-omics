"""
ARIA Bulk RNA-seq Differential Expression Script
-------------------------------------------------
Full bulk RNA-seq pipeline executed inside aria-rna-env by EnvironmentManager.

Fixes vs old inline implementation:
  1. Design factor extracted from biological intent (not hardcoded "sample")
  2. Robust metadata parsing: detects groups from column names automatically
  3. Sample outlier detection (PCA-based) before running DESeq2; primary DE
     retains all samples and outlier removal is reported as sensitivity only
  4. Pathway enrichment via local ORA after DE (GO BP, KEGG, Reactome)
  5. Visualizations saved as SVG (volcano, sample PCA, heatmap)
  6. Runs in aria-rna-env (isolated) not in aria-env (base)

Input params:
    files:           list  — count matrix files (.tsv, .csv, .txt)
                            OR list of per-sample count files
    metadata_file:   str   — explicit metadata TSV with
                            sample, condition, batch columns
    design_factor:   str   — column in metadata to test (e.g. "condition")
    comparison:      dict  — {"numerator": "treated", "denominator": "control"}
    organism:        str   — for pathway databases (human/mouse)
    genome:          str
    output_dir:      str
    run_pathways:    bool  (default: True)
    padj_threshold:  float (default: 0.05)
    lfc_threshold:   float (default: 1.0)

Output:
    {
      "status":          "success",
      "n_genes_tested":  int,
      "n_significant":   int,
      "n_upregulated":   int,
      "n_downregulated": int,
      "top_genes":       [{"gene", "log2fc", "padj", "direction"}],
      "sample_qc":       {"n_samples", "outliers", "pca_variance"},
      "pathways": {
        "GO_BP":    [{"term", "genes", "padj", "NES"}],
        "KEGG":     [...],
        "Reactome": [...],
      },
      "plots":           {"volcano": str, "pca": str, "heatmap": str},
      "design_used":     str,
      "comparison_used": dict,
      "warnings":        [str]
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script
from aria.utils.count_classifier import classify_matrix, validate_raw_count_matrix
from aria.utils.stats import (
    contrast_family_significance,
    fdr_advanced_methods_disclosure,
    preregister_contrast_family,
)
from pathlib import Path

# P2-8: helpers extracted to the aria.scripts.rna_bulk subpackage and
# re-exported here so the public surface (and tests) are unchanged.
from aria.scripts.rna_bulk.gtf_io import (  # noqa: E402,F401
    _load_symbol_map, _locate_gtf, _gtf_to_symbol_map,
    _load_gene_annotation, _parse_gtf_biotype_and_length, _to_symbols,
)
from aria.scripts.rna_bulk.ora import (  # noqa: E402,F401
    _enrichr_enrichment, _run_pathway_enrichment, _get_gene_sets,
    _gseapy_organism, _mock_pathways,
)
from aria.scripts.rna_bulk.plots import (  # noqa: E402,F401
    _plot_pca_mds, _save_single_dr_plot, _generate_plots,
    _plot_heatmap, _plot_sample_pca,
)
# A7: the remaining counts/contrasts/transforms/QC/DESeq2 helper groups were
# extracted to the aria.scripts.rna_bulk subpackage and are re-exported here so
# the public surface (and every caller/test) is unchanged.
from aria.scripts.rna_bulk.contrasts import (  # noqa: E402,F401
    _slugify, _format_top_genes, _suggest_contrasts, _auto_contrasts,
    _contrast_overlap,
)
from aria.scripts.rna_bulk.deseq2 import (  # noqa: E402,F401
    _BENIGN_WARNING_CATEGORIES, _build_design_formula, _resolve_covariates,
    _shrink_coeff, _serialize_fit_warnings, _run_deseq2, _mock_de_result,
)
from aria.scripts.rna_bulk.transforms import (  # noqa: E402,F401
    _run_vst, _select_variable_genes, _compute_tpm,
)
from aria.scripts.rna_bulk.counts import (  # noqa: E402,F401
    _load_counts, _enforce_metadata_correspondence, _metadata_inference_allowed,
    _load_or_infer_metadata, _aggregate_technical_replicates, _infer_groups,
    _resolve_comparison,
)
from aria.scripts.rna_bulk.qc import (  # noqa: E402,F401
    _sample_qc, _prune_outliers_for_design, _run_outlier_sensitivity,
)

# The paper plot theme (`_P`) now lives in aria/scripts/rna_bulk/plots.py.


def bulk_rna_de(params: dict) -> dict:
    from pathlib import Path
    import numpy as np
    import warnings as warn_mod
    # F8 (preprint audit): do NOT globally silence warnings — that hid pydeseq2
    # convergence/numeric warnings outside the audit trail. Suppress only benign
    # third-party API-churn categories; numeric/convergence warnings during the fit
    # are captured and surfaced by _run_deseq2 (see _serialize_fit_warnings).
    for _benign in _BENIGN_WARNING_CATEGORIES:
        warn_mod.filterwarnings("ignore", category=_benign)

    files          = params.get("files", [])
    metadata_file  = params.get("metadata_file", "")
    design_factor  = params.get("design_factor", "condition")
    # P0-4: covariates confirmed at DesignAgent CHECKPOINT 2.4 (e.g. batch).
    covariates     = params.get("covariates", []) or []
    technical_replicate_col = str(
        params.get("technical_replicate_col", "") or ""
    )
    # P1-1/ADR-023: apeGLM LFC shrinkage on by default (bulk = pseudobulk rigor).
    lfc_shrink     = bool(params.get("lfc_shrink", True))
    # P1-1c: pre-register the contrast-FDR family before any p-values are seen.
    fdr_family     = preregister_contrast_family(params.get("fdr_family", "per_contrast"))

    metadata_inference_allowed = _metadata_inference_allowed(params)

    # v3: accept list of contrasts OR single comparison (backward compat)
    contrasts_in   = params.get("contrasts", [])
    if not contrasts_in:
        single_comp = params.get("comparison", {})
        if single_comp:
            contrasts_in = [{
                "numerator":   single_comp.get("numerator", ""),
                "denominator": single_comp.get("denominator", ""),
                "name":        f"{single_comp.get('numerator','?')} vs "
                               f"{single_comp.get('denominator','?')}",
            }]

    organism       = params.get("organism", "Homo sapiens")
    output_dir     = params.get("output_dir", "/tmp/aria_bulk_de")
    run_pathways   = bool(params.get("run_pathways", True))
    padj_thr       = float(params.get("padj_threshold", 0.05))
    lfc_thr        = float(params.get("lfc_threshold", 1.0))
    min_reps       = int(params.get("min_replicates_per_condition", 3))
    # B5: explicit, audited exclusion of count columns that lack metadata. Only
    # genuine orphan columns may be named here; everything else needs metadata.
    excluded_samples = params.get("excluded_samples", []) or []
    allow_mock     = mocks_allowed(params)
    warnings       = []

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if not metadata_file and not metadata_inference_allowed:
        return {
            "status":     "error",
            "error_type": "MetadataRequired",
            "details": (
                "Production bulk RNA DE requires an explicit metadata_file "
                "aligned to the count matrix. Column-name metadata inference is "
                "available only with an explicit legacy/dev opt-in."
            ),
        }
    if metadata_file and not Path(metadata_file).exists() and not metadata_inference_allowed:
        return {
            "status":     "error",
            "error_type": "MetadataRequired",
            "details": (
                f"metadata_file does not exist: {metadata_file}. Production "
                "bulk RNA DE will not fall back to column-name inference."
            ),
        }

    # ── 1. Load counts matrix ─────────────────────────────────────────────
    allow_nonraw = bool(params.get("allow_nonraw_counts", False))
    counts, load_warn, count_meta = _load_counts(files, allow_nonraw=allow_nonraw)
    warnings.extend(load_warn)
    if counts is None:
        if count_meta.get("refused"):
            # Raw-count guard (B10): hard-refuse a non-raw matrix rather than
            # coerce it into pseudo-counts for DESeq2.
            return {
                "status":       "error",
                "error_type":   count_meta.get("error_type", "NonRawCounts"),
                "details":      count_meta.get("details", ""),
                "count_source": count_meta.get("kind"),
            }
        return {
            "status":     "error",
            "error_type": "CountsLoadFailed",
            "details":    "Could not load count matrix from provided files.",
        }
    count_source = count_meta.get("count_source", "raw_counts")

    # ── 2. Load or infer metadata ─────────────────────────────────────────
    metadata, meta_warn = _load_or_infer_metadata(
        counts, metadata_file, design_factor,
        allow_inference=metadata_inference_allowed,
    )
    warnings.extend(meta_warn)
    if metadata is None:
        return {
            "status":     "error",
            "error_type": "MetadataFailed",
            "details": (
                "Could not construct sample metadata. "
                "Provide a metadata TSV with 'sample' and condition columns, "
                "aligned to the count matrix. Column-name inference is a "
                "legacy/dev opt-in, not production behavior."
            ),
        }
    if design_factor not in metadata.columns:
        return {
            "status":     "error",
            "error_type": "DesignFactorMissing",
            "details": (
                f"Column '{design_factor}' not found in metadata. "
                f"Available columns: {list(metadata.columns)}"
            ),
        }

    # ── 2a. Total metadata correspondence (B5) ────────────────────────────
    # Partial metadata must not silently reduce n. Every count column needs a
    # metadata row; genuine orphans may only be dropped via an explicit, audited
    # excluded_samples list. This runs before QC/aggregation/DE.
    try:
        counts, metadata_correspondence = _enforce_metadata_correspondence(
            counts, metadata, excluded_samples,
        )
    except ValueError as exc:
        return {
            "status":     "error",
            "error_type": "MetadataCorrespondenceError",
            "details":    str(exc),
        }
    if metadata_correspondence["n_excluded"]:
        warnings.append(
            "Explicit audited exclusion (B5): dropped count columns without "
            f"metadata before QC/DE: {metadata_correspondence['excluded_samples']}."
        )

    technical_replicate_aggregation = {
        "ran": False,
        "reason": "not_declared",
    }
    if technical_replicate_col:
        try:
            counts, metadata, technical_replicate_aggregation = (
                _aggregate_technical_replicates(
                    counts,
                    metadata,
                    design_factor=design_factor,
                    unit_col=technical_replicate_col,
                    covariates=covariates,
                )
            )
        except ValueError as exc:
            return {
                "status": "error",
                "error_type": "TechnicalReplicateContractError",
                "details": str(exc),
            }
        warnings.append(
            "Technical libraries were summed within biological units before "
            "QC and DE; inferential replicate counts use biological units."
        )

    # ── 2b. Load gene annotation maps (biotype, length) from GTF ────────
    # Used for DR filtering (protein_coding) and TPM computation.
    # Graceful fallback if GTF not locatable → warnings explain consequences.
    gene_annotation = _load_gene_annotation(files, warnings)
    biotype_map     = gene_annotation.get("biotype", {})
    gene_lengths    = gene_annotation.get("length", {})

    # ── 3. Sample QC (shared across contrasts) ───────────────────────────
    sample_qc = _sample_qc(counts, metadata, output_dir, warnings,
                           biotype_map=biotype_map)
    flagged_outliers = list(sample_qc.get("outliers", []) or [])
    outlier_samples = _prune_outliers_for_design(
        sample_qc.get("outliers", []), metadata, design_factor, warnings
    )
    sample_qc["outliers"] = flagged_outliers
    sample_qc["candidate_outliers"] = flagged_outliers
    sample_qc["outliers_removed_primary"] = []
    sample_qc["sensitivity_outliers_removed"] = outlier_samples
    sample_qc["outlier_policy"] = (
        "primary_includes_all_samples; sensitivity_removes_flagged_outliers_when_design_allows"
    )
    if outlier_samples:
        warnings.append(
            f"Sample outliers retained in primary DE and removed only in "
            f"sensitivity analysis when design support allows: {outlier_samples}."
        )

    # ── 3b. TPM supplementary table ──────────────────────────────────────
    # TPM is NOT used internally by ARIA (DESeq2 uses raw, DR uses VST).
    # Computed here as a supplementary output for downstream tools that
    # require TPM (ssGSEA, deconvolution methods, external sharing).
    try:
        tpm_matrix = _compute_tpm(counts, gene_lengths, warnings)
        if tpm_matrix is not None:
            tpm_path = Path(output_dir) / "counts_tpm.tsv"
            tpm_matrix.to_csv(tpm_path, sep="\t")
            warnings.append(
                f"TPM table written ({len(tpm_matrix)} genes) — "
                f"supplementary only; not used for DE or PCA."
            )
    except Exception as e:
        warnings.append(f"TPM export failed (non-fatal): {e}")

    # ── 4. Filter low-count genes (shared across contrasts) ──────────────
    min_samples = max(2, metadata.shape[0] // 4)
    keep        = (counts > 10).sum(axis=1) >= min_samples
    counts_filt = counts[keep]
    n_filtered  = int((~keep).sum())
    if n_filtered:
        warnings.append(
            f"{n_filtered} low-count genes removed "
            f"(< 10 counts in < {min_samples} samples)."
        )

    # ── 5. Require explicit contrasts ───────────────────────────────────
    if not contrasts_in:
        suggestions = _suggest_contrasts(metadata, design_factor)
        return {
            "status":     "error",
            "error_type": "ExplicitContrastRequired",
            "details": (
                "Differential expression requires an explicit contrast with "
                "both numerator/test level and denominator/reference level. "
                "ARIA will not choose the reference level from group order."
            ),
            "design_factor": design_factor,
            "available_groups": sorted(str(g) for g in metadata[design_factor].unique()),
            "suggested_contrasts": suggestions,
            "warnings": warnings,
        }

    # ── 6. Run each contrast ─────────────────────────────────────────────
    contrast_results = []
    # P1-1c: per-contrast (gene -> pvalue, log2fc) for the pooled contrast-family BH.
    family_stats: dict = {}

    for contrast in contrasts_in:
        num  = str(contrast.get("numerator",   "")).strip()
        den  = str(contrast.get("denominator", "")).strip()
        name = contrast.get("name", f"{num} vs {den}")

        if not num or not den:
            warnings.append(
                f"Skipping contrast '{name}': numerator and denominator "
                "must both be explicit."
            )
            continue

        # Validate groups exist in metadata
        available_groups = list(metadata[design_factor].unique())
        if num not in available_groups or den not in available_groups:
            warnings.append(
                f"Skipping contrast '{name}': "
                f"groups {num!r}, {den!r} not in {available_groups}"
            )
            continue

        # Run DE for this contrast
        de_result, de_warn = _run_deseq2(
            counts_filt, metadata, design_factor,
            num, den, padj_thr, lfc_thr,
            allow_mock=allow_mock,
            min_replicates_per_condition=min_reps,
            covariates=covariates,
            lfc_shrink=lfc_shrink,
        )

        if de_result.get("status") == "error":
            contrast_results.append({
                "name":         name,
                "numerator":    num,
                "denominator":  den,
                "status":       "error",
                "error":        de_result.get("details", ""),
                "warnings":     de_warn,
            })
            warnings.extend([f"[{name}] {w}" for w in de_warn])
            continue

        warnings.extend([f"[{name}] {w}" for w in de_warn])

        # Load gene symbol map (used both for pathway enrichment and
        # for annotating top_genes with HGNC symbols)
        symbol_map = _load_symbol_map(files, warnings)
        background_symbols = _to_symbols(list(counts_filt.index), symbol_map) \
            if symbol_map else [str(g) for g in counts_filt.index]
        background_symbols = sorted({
            str(g) for g in background_symbols
            if g and str(g).lower() != "nan"
        })

        # Pathway enrichment per contrast
        pathways = {}
        ora_meta = {"method": "none", "gene_set_versions": {}}
        if run_pathways and de_result.get("sig_genes"):
            pathways, pw_warn, ora_meta = _run_pathway_enrichment(
                sig_genes=de_result["sig_genes"],
                up_genes=de_result.get("up_genes", []),
                down_genes=de_result.get("down_genes", []),
                organism=organism,
                output_dir=output_dir,
                symbol_map=symbol_map,
                background_genes=background_symbols,
                allow_mock=allow_mock,
            )
            warnings.extend([f"[{name}] {w}" for w in pw_warn])

        # Plots per contrast — one volcano and one heatmap each
        contrast_dir = Path(output_dir) / _slugify(name)
        contrast_dir.mkdir(exist_ok=True)
        figures_dir = contrast_dir / "figures"
        tables_dir  = contrast_dir / "tables"
        figures_dir.mkdir(exist_ok=True)
        tables_dir.mkdir(exist_ok=True)

        # Inject symbol_map into sample_qc so heatmap/DR plots can label
        # rows with HGNC symbols instead of bare Ensembl IDs.
        sample_qc["_symbol_map"] = symbol_map

        plots = _generate_plots(
            de_result=de_result,
            sample_qc=sample_qc,
            counts_filt=counts_filt,
            metadata=metadata,
            design_factor=design_factor,
            output_dir=str(figures_dir),
            padj_thr=padj_thr,
            lfc_thr=lfc_thr,
            title_suffix=name,
        )

        # ── Pathway visualizations (ORA dotplots + GSEA running sums) ──
        # Lazy import: avoid pulling matplotlib/blitzgsea at module load
        try:
            from aria.scripts.rna_pathway_viz import (
                make_ora_dotplot,
                make_gsea_running_sums,
                export_de_table,
                export_pathways_table,
            )

            # ORA dotplots — one per database
            ora_dotplots = {}
            for db, terms in (pathways or {}).items():
                if not isinstance(terms, list) or not terms:
                    continue
                dot_path = figures_dir / f"pathway_dotplot_{_slugify(db)}.png"
                result_path = make_ora_dotplot(
                    pathways_list=terms,
                    db_name=db,
                    contrast_name=name,
                    output_path=str(dot_path),
                )
                if result_path:
                    ora_dotplots[db] = result_path
            plots["ora_dotplots"] = ora_dotplots

            # GSEA running sums (uses ranked log2FC, complementary to ORA)
            results_df = de_result.get("results")
            if results_df is not None and len(results_df) > 0:
                gsea_out = make_gsea_running_sums(
                    de_results_df=results_df,
                    symbol_map=symbol_map,
                    contrast_name=name,
                    output_dir=str(figures_dir),
                    organism=_gseapy_organism(organism),
                )
                plots["gsea_running_sums"] = gsea_out.get("running_sums", [])
                plots["gsea_top_table"]    = gsea_out.get("top_table_fig")
                plots["gsea_table"]        = gsea_out.get("gsea_table")
                if gsea_out.get("n_pathways"):
                    warnings.append(
                        f"[{name}] GSEA: {gsea_out['n_pathways']} pathways "
                        f"at FDR<0.25"
                    )

            # Export supplementary tables (DE genes + pathways)
            de_tsv = export_de_table(
                de_results_df=results_df,
                symbol_map=symbol_map,
                contrast_name=name,
                output_path=str(tables_dir / "de_genes.tsv"),
                padj_thr=padj_thr,
            )
            pw_tsv = export_pathways_table(
                pathways_dict=pathways,
                contrast_name=name,
                output_path=str(tables_dir / "pathways.tsv"),
            )
            plots["tables"] = {
                "de_genes": de_tsv,
                "pathways": pw_tsv,
            }

        except Exception as e:
            warnings.append(
                f"[{name}] Pathway visualization failed: {str(e)[:150]}"
            )

        # Format top genes (annotated with HGNC symbols when known)
        top_genes = _format_top_genes(de_result, symbol_map=symbol_map)

        # Convert all sig_genes to symbols (or keep as IDs) for overlap computation
        all_sig = de_result.get("sig_genes", []) or []
        all_sig_symbols = _to_symbols(all_sig, symbol_map) if symbol_map \
                            else list(all_sig)

        contrast_results.append({
            "name":              name,
            "numerator":         num,
            "denominator":       den,
            "status":            "success",
            "n_genes_tested":    int(counts_filt.shape[0]),
            "n_significant":     int(de_result.get("n_sig", 0)),
            "n_upregulated":     int(de_result.get("n_up", 0)),
            "n_downregulated":   int(de_result.get("n_down", 0)),
            "n_replicates":      de_result.get("n_replicates"),
            "low_power_warning": bool(de_result.get("low_power_warning", False)),
            "low_power_reason":  de_result.get("low_power_reason"),
            "dispersion_estimate": de_result.get("dispersion_estimate"),
            "mean_expression_estimate": de_result.get("mean_expression_estimate"),
            "power_estimate_at_lfc_min": de_result.get("power_estimate_at_lfc_min"),
            "design_check":       de_result.get("design_check"),
            # P0-4: the design actually fitted by DESeq2 (covariate-adjusted),
            # surfaced so the report Methods can state it verbatim.
            "fitted_design_formula": de_result.get("fitted_design_formula"),
            "covariates_adjusted":   de_result.get("covariates_adjusted", []),
            "covariates_dropped":    de_result.get("covariates_dropped", []),
            # P1-1/ADR-023: reported log2fc is the apeGLM-shrunken estimate; the
            # raw MLE is preserved per gene as log2fc_raw and in the DE table.
            "lfc_shrinkage":     de_result.get("lfc_shrinkage"),
            # P1-1b: padj comes from a Wald test against |LFC| > lfc_threshold,
            # so significance no longer applies a second post-hoc LFC filter.
            "lfc_threshold_test": de_result.get("lfc_threshold_test"),
            "top_genes":         top_genes,
            # Full DE list for cross-contrast overlap (in symbols when available,
            # else Ensembl IDs — both work for set intersection)
            "all_sig_genes":     all_sig_symbols,
            "all_sig_gene_ids":   [str(g) for g in all_sig],
            # Direction split (gene IDs) so the BiologicalSynthesisAgent can score
            # direction concordance of shared genes across contrasts.
            "up_gene_ids":       [str(g) for g in (de_result.get("up_genes") or [])],
            "down_gene_ids":     [str(g) for g in (de_result.get("down_genes") or [])],
            "pathways":          pathways,
            "pathway_background": {
                "background_size": len(background_symbols),
                "background_source": "dataset_expressed_genes",
            },
            # P1-7/W-PRIV: ORA engine + the exact versioned gene-set release per
            # database, so methodology.json records reproducible provenance and
            # discloses any database skipped for privacy/availability reasons.
            "pathway_ora": ora_meta,
            "plots":             plots,
            "contrast_dir":      str(contrast_dir),
            "figures_dir":       str(figures_dir),
            "tables_dir":        str(tables_dir),
        })

        # P1-1c: stash this contrast's per-gene p-value + shrunken log2fc for the
        # pooled contrast-family BH computed once all contrasts have run.
        rdf = de_result.get("results")
        if rdf is not None and len(rdf) > 0:
            family_stats[name] = {
                str(g): {"pvalue": float(row["pvalue"]),
                         "log2fc": float(row["log2FoldChange"])}
                for g, row in rdf.dropna(subset=["pvalue"]).iterrows()
            }

    if not contrast_results:
        return {
            "status":     "error",
            "error_type": "AllContrastsFailed",
            "details":    "No contrasts produced valid results.",
            "warnings":   warnings,
        }

    # ── 6b. Contrast-family FDR (P1-1c) ──────────────────────────────────
    # Always compute the pooled BH across the contrast family for audit; when the
    # pre-registered family is "global", the primary significance call follows it.
    n_tests_family = sum(len(v) for v in family_stats.values())
    fam = contrast_family_significance(
        family_stats, padj_max=padj_thr, lfc_min=None)
    for cr in contrast_results:
        if cr.get("status") != "success":
            continue
        f = fam.get(cr["name"])
        if not f:
            continue
        cr["n_significant_contrast_family"] = f["n_sig"]
        cr["sig_genes_contrast_family"] = [str(g) for g in f["sig_genes"]]
        if fdr_family["fdr_family"] == "global":
            # Primary follows the pre-registered global contrast-family BH.
            cr["n_significant"]   = f["n_sig"]
            cr["n_upregulated"]   = f["n_up"]
            cr["n_downregulated"] = f["n_down"]
            cr["all_sig_genes"]   = (
                _to_symbols(f["sig_genes"], symbol_map) if symbol_map
                else [str(g) for g in f["sig_genes"]]
            )
            cr["all_sig_gene_ids"] = [str(g) for g in f["sig_genes"]]

    # ── 6c. Primary + outlier sensitivity (P1-5) ────────────────────────
    # Primary DE above deliberately used all samples. If QC flagged outliers and
    # removal does not break the design, rerun DE on the pruned matrix and record
    # whether the conclusion changes. Sensitivity never replaces primary calls.
    outlier_sensitivity = _run_outlier_sensitivity(
        counts=counts,
        metadata=metadata,
        design_factor=design_factor,
        contrasts_in=contrasts_in,
        flagged_outliers=flagged_outliers,
        removable_outliers=outlier_samples,
        primary_contrasts=contrast_results,
        padj_thr=padj_thr,
        lfc_thr=lfc_thr,
        allow_mock=allow_mock,
        min_reps=min_reps,
        covariates=covariates,
        lfc_shrink=lfc_shrink,
        fdr_family=fdr_family,
        warnings=warnings,
    )

    # ── 7. Aggregate summary ─────────────────────────────────────────────
    total_sig  = sum(c.get("n_significant",   0) for c in contrast_results)
    total_up   = sum(c.get("n_upregulated",   0) for c in contrast_results)
    total_down = sum(c.get("n_downregulated", 0) for c in contrast_results)

    # Cross-contrast overlap (shared DE genes)
    overlap_info = _contrast_overlap(contrast_results)

    # Strip non-JSON-serializable items before returning
    # (vst_matrix is a DataFrame, symbol_map is big and already elsewhere)
    sample_qc_clean = {k: v for k, v in sample_qc.items()
                        if k not in ("vst_matrix", "_symbol_map")}

    # ── Methodology decisions record (NEW in v3.8) ──────────────────────
    # Explicit record of normalization choices, gene filters, and their
    # justifications. Rendered as a table in the HTML report so reviewers
    # and collaborators can audit methodology without reading the code.
    methodology = {
        "decisions": [
            {
                "step":           "Differential expression (DESeq2)",
                "input":          "Raw integer counts",
                "normalization":  "DESeq2 median-of-ratios (internal)",
                "gene_filter":    f"≥10 counts in ≥{max(2, metadata.shape[0] // 4)} samples",
                "justification":  "DESeq2 handles library-size normalization internally; external normalization would break the negative-binomial model assumptions.",
            },
            {
                "step":           "Sample outliers",
                "input":          "All count-matrix samples",
                "normalization":  "Primary DE retains all samples",
                "gene_filter":    "Design-safe removal is evaluated only as sensitivity",
                "justification":  "Automatic pre-DE pruning can turn a null primary result into a significant result; ARIA reports primary and outlier-removal sensitivity separately.",
            },
            {
                "step":           "PCA + MDS (sample-level structure)",
                "input":          "VST-transformed counts",
                "normalization":  "Variance-Stabilizing Transformation (pydeseq2)",
                "gene_filter":    f"Top 2000 most-variable protein_coding genes",
                "justification":  "VST produces homoscedastic values suitable for Euclidean-based methods (DESeq2 authors' recommendation). Protein-coding filter removes pseudogene/ncRNA noise. Variable-gene filter focuses on informative signal.",
            },
            {
                "step":           "Heatmap (padj top 50)",
                "input":          "VST (or log2(counts+1) fallback)",
                "normalization":  "Row z-score",
                "gene_filter":    "50 most significant DE genes (sorted by padj)",
                "justification":  "Statistically confident signal; symbol-annotated rows.",
            },
            {
                "step":           "Heatmap (|log2FC| top 50)",
                "input":          "VST (or log2(counts+1) fallback)",
                "normalization":  "Row z-score",
                "gene_filter":    "50 DE genes with largest effect sizes",
                "justification":  "Complementary view — surfaces largest effect sizes, which may not have smallest padj for low-count genes.",
            },
            {
                "step":           "Pathway enrichment (ORA)",
                "input":          "DE gene symbols",
                "normalization":  "Local hypergeometric test against versioned GMT libraries (Enrichr is opt-in)",
                "gene_filter":    "padj<0.05 from Wald lfcThreshold test; ORA universe = genes retained after expression filtering",
                "justification":  "Standard over-representation for discrete gene lists, computed locally against the dataset-expressed background so the gene list never leaves the machine (W-PRIV); the exact gene-set release is recorded in pathway_ora.gene_set_versions.",
            },
            {
                "step":           "GSEA (pre-ranked)",
                "input":          "All DE-tested genes ranked by log2FC",
                "normalization":  "blitzgsea running sum",
                "gene_filter":    "None (uses full ranking)",
                "justification":  "More sensitive than ORA to coordinated small effects; complementary to cutoff-based enrichment.",
            },
            {
                "step":           "TPM (supplementary export)",
                "input":          "Raw counts + exon-union gene lengths",
                "normalization":  "TPM (transcripts per million)",
                "gene_filter":    "Genes with length annotation in GTF",
                "justification":  "Supplementary table for downstream tools that require TPM (ssGSEA, deconvolution). NOT used by ARIA for DE or PCA — both use more appropriate methods.",
            },
        ],
        "n_genes_dr":         sample_qc.get("n_genes_dr", 0),
        "n_protein_coding":   sample_qc.get("n_protein_coding", 0),
        "technical_replicate_aggregation": technical_replicate_aggregation,
        "metadata_correspondence": metadata_correspondence,
    }

    return {
        "status":           "success",
        "n_contrasts":      len(contrast_results),
        "contrasts":        contrast_results,
        "n_significant":    total_sig,       # legacy total
        "n_upregulated":    total_up,
        "n_downregulated":  total_down,
        "sample_qc":        sample_qc_clean,
        "outlier_sensitivity": outlier_sensitivity,
        "design_used":      f"~{design_factor}",
        "padj_threshold":   padj_thr,
        "lfc_threshold":    lfc_thr,
        # P1-1c: pre-registered contrast-FDR family + the pooled-BH family size.
        "fdr_family":       {**fdr_family, "n_tests_family": n_tests_family},
        # P1-2 closure (ADR-027): IHW + s-values honestly disclosed as not
        # implemented, never faked. Primary FDR stays pre-registered BH.
        "fdr_advanced_methods": fdr_advanced_methods_disclosure(),
        "overlap":          overlap_info,
        "methodology":      methodology,
        "count_source":     count_source,
        "technical_replicate_aggregation": technical_replicate_aggregation,
        "metadata_correspondence": metadata_correspondence,
        "warnings":         warnings,
    }


if __name__ == "__main__":
    run_script(bulk_rna_de)
