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
    metadata_file:   str   (optional) — explicit metadata TSV with
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
from aria.utils.count_classifier import classify_matrix
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

# The paper plot theme (`_P`) now lives in aria/scripts/rna_bulk/plots.py.


def bulk_rna_de(params: dict) -> dict:
    from pathlib import Path
    import numpy as np
    import warnings as warn_mod
    warn_mod.filterwarnings("ignore")

    files          = params.get("files", [])
    metadata_file  = params.get("metadata_file", "")
    design_factor  = params.get("design_factor", "condition")
    # P0-4: covariates confirmed at DesignAgent CHECKPOINT 2.4 (e.g. batch).
    covariates     = params.get("covariates", []) or []
    # P1-1/ADR-023: apeGLM LFC shrinkage on by default (bulk = pseudobulk rigor).
    lfc_shrink     = bool(params.get("lfc_shrink", True))
    # P1-1c: pre-register the contrast-FDR family before any p-values are seen.
    fdr_family     = preregister_contrast_family(params.get("fdr_family", "per_contrast"))

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
    allow_mock     = mocks_allowed(params)
    warnings       = []

    Path(output_dir).mkdir(parents=True, exist_ok=True)

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
        counts, metadata_file, design_factor
    )
    warnings.extend(meta_warn)
    if metadata is None:
        return {
            "status":     "error",
            "error_type": "MetadataFailed",
            "details": (
                "Could not construct sample metadata. "
                "Provide a metadata TSV with 'sample' and condition columns, "
                "or name samples as: condition_replicate (e.g. ctrl_1, treat_1)."
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
        "warnings":         warnings,
    }


def _slugify(s: str) -> str:
    """Make a string safe for use as a directory name."""
    import re
    return re.sub(r"[^a-z0-9_]+", "_", str(s).lower()).strip("_") or "contrast"


def _format_top_genes(de_result: dict, symbol_map: dict = None) -> list:
    """Format top genes from a DE result. Adds symbol if mapping available."""
    results = de_result.get("results")
    if results is None or len(results) == 0:
        return []
    sm = symbol_map or {}
    top = []
    for g in de_result.get("sig_genes", [])[:30]:
        if g not in results.index:
            continue
        row = results.loc[g]
        try:
            clean_id = str(g).split(".")[0]
            symbol   = sm.get(clean_id, "")
            top.append({
                "gene":      g,                    # Ensembl ID (or whatever)
                "symbol":    symbol or g,          # HGNC symbol if known, else fallback to ID
                "log2fc":    round(float(row["log2FoldChange"]), 3),
                "padj":      float(row["padj"]),
                "direction": "up" if row["log2FoldChange"] > 0 else "down",
            })
        except Exception:
            continue
    return top


def _suggest_contrasts(metadata, design_factor: str) -> list[dict]:
    """Suggest candidate contrasts without authorizing execution.

    P0-5: suggestions are display-only. The caller must pass one or more of
    them back explicitly as ``contrasts`` before DE can run.
    """
    groups = sorted(metadata[design_factor].unique())

    if len(groups) < 2:
        return []

    # Identify control
    ctrl_keywords = ["wt", "wildtype", "control", "ctrl",
                     "vehicle", "dmso", "untreated", "mock",
                     "normal", "healthy", "baseline", "scramble"]
    control = None
    for kw in ctrl_keywords:
        for g in groups:
            if g.lower() == kw:
                control = g
                break
        if control:
            break
    if not control:
        for kw in ctrl_keywords:
            for g in groups:
                if kw in g.lower():
                    control = g
                    break
            if control:
                break

    contrasts = []
    if control:
        for g in groups:
            if g == control:
                continue
            contrasts.append({
                "numerator":   g,
                "denominator": control,
                "name":        f"{g} vs {control}",
            })
    else:
        # Pairwise suggestions only. Do not claim a reference was selected.
        ref = groups[0]
        for g in groups[1:]:
            contrasts.append({
                "numerator":   g,
                "denominator": ref,
                "name":        f"{g} vs {ref}",
            })

    return contrasts


def _auto_contrasts(metadata, design_factor: str) -> tuple:
    """Backward-compatible wrapper returning suggestions, not executable DE.

    Kept for legacy diagnostics that import the helper. Production execution
    must call ``bulk_rna_de`` with explicit ``contrasts``.
    """
    suggestions = _suggest_contrasts(metadata, design_factor)
    warning = (
        "Automatic contrast generation is disabled for production DE. "
        "Use these suggestions only after explicit user confirmation."
    )
    return suggestions, [warning] if suggestions else []


def _contrast_overlap(contrast_results: list) -> dict:
    """
    Compute DE gene overlap between contrasts.
    Uses the FULL list of significant DE genes per contrast (not just
    the top 30 used for display) — otherwise overlap counts are
    misleadingly small.
    """
    successful = [c for c in contrast_results if c.get("status") == "success"]
    if len(successful) < 2:
        return {}

    # Prefer all_sig_genes (full DE list) over top_genes (display top 30)
    gene_sets = {}
    for c in successful:
        if c.get("all_sig_genes"):
            gene_sets[c["name"]] = set(c["all_sig_genes"])
        else:
            # Fallback if all_sig_genes not present
            gene_sets[c["name"]] = set(g["gene"] for g in c.get("top_genes", []))

    names = list(gene_sets.keys())
    overlaps = {}
    for i, a in enumerate(names):
        for b in names[i+1:]:
            shared = gene_sets[a] & gene_sets[b]
            n_a, n_b = len(gene_sets[a]), len(gene_sets[b])
            # Hypergeometric expectation if independent: rough sanity check
            # (assumes ~30k expressed genes universe; can be refined with the
            # actual n_genes_tested if passed in)
            jaccard = len(shared) / max(len(gene_sets[a] | gene_sets[b]), 1)
            overlaps[f"{a} ∩ {b}"] = {
                "n_shared":      len(shared),
                "n_in_first":    n_a,
                "n_in_second":   n_b,
                "jaccard":       round(jaccard, 3),
                "shared_genes":  sorted(shared)[:50],   # cap for serialization
            }
    return overlaps


# ── Data loading ──────────────────────────────────────────────────────────────

def _load_counts(files: list, allow_nonraw: bool = False) -> tuple:
    """
    Load counts matrix from various formats.
    Returns (DataFrame, warnings_list, count_meta).
    genes × samples orientation enforced.

    Raw-count guard (audit 2026-05-29, B10 / P-RAWCLASS): DESeq2 requires raw
    integer counts. The loaded matrix is classified before rounding so that
    TPM/CPM/FPKM/log-normalized/scaled inputs are NOT silently coerced into
    pseudo-counts. Non-raw matrices are hard-refused unless ``allow_nonraw`` is
    set, in which case they are coerced at explicit low confidence. ``count_meta``
    always carries ``count_source`` / ``kind`` (or ``refused`` details).
    """
    import pandas as pd
    warnings = []

    valid = [f for f in files if Path(f).exists()]
    if not valid:
        return None, ["No valid count files found."], {}

    # Detect format
    count_files = [f for f in valid
                   if any(f.endswith(x) for x in
                          [".tsv", ".csv", ".txt", ".counts",
                           ".featureCounts", ".htseq"])]
    if not count_files:
        # Try all valid files
        count_files = valid

    # Try to load as single matrix or merge multiple per-sample files
    if len(count_files) == 1:
        sep = "\t" if count_files[0].endswith(".tsv") else ","
        try:
            counts = pd.read_csv(count_files[0], sep=sep, index_col=0,
                                  comment="#")
            # Drop non-count columns (featureCounts format)
            drop_cols = [c for c in counts.columns
                         if c in ("Chr","Start","End","Strand","Length")]
            counts = counts.drop(columns=drop_cols, errors="ignore")
            # Keep only numeric columns
            counts = counts.select_dtypes(include="number")
            if counts.empty:
                return None, ["Count matrix has no numeric columns."], {}
        except Exception as e:
            return None, [f"Failed to load {count_files[0]}: {e}"], {}
    else:
        # Multiple per-sample files — merge by gene ID
        frames = []
        for f in count_files:
            sep = "\t" if f.endswith(".tsv") else ","
            try:
                df = pd.read_csv(f, sep=sep, index_col=0, comment="#",
                                  header=None)
                df.columns = [Path(f).stem]
                frames.append(df.select_dtypes(include="number"))
            except Exception as e:
                warnings.append(f"Skipping {f}: {e}")
        if not frames:
            return None, ["Could not load any count files."], {}
        counts = pd.concat(frames, axis=1).fillna(0)

    # Ensure genes × samples (more genes than samples in typical experiments)
    if counts.shape[1] > counts.shape[0]:
        warnings.append(
            f"Transposing matrix: detected {counts.shape[1]} rows × "
            f"{counts.shape[0]} cols — expected genes as rows."
        )
        counts = counts.T

    # Raw-count guard (B10 / P-RAWCLASS): classify BEFORE rounding so a
    # normalized matrix is not silently turned into pseudo-counts.
    source_hint = ";".join(str(f) for f in count_files)
    info = classify_matrix(
        counts.values,
        gene_ids=list(counts.index),
        source_hint=source_hint,
    )
    if info["is_raw_counts"]:
        count_source = "raw_counts"
    elif allow_nonraw:
        count_source = "coerced_nonraw"
        warnings.append(
            f"Count matrix does not look like raw counts "
            f"(kind={info['kind']}, score={info.get('raw_count_score', 0):.2f}, "
            f"max={info['max']:.2f}); coercing to "
            f"integers because allow_nonraw_counts=True. DESeq2 results are "
            f"LOW CONFIDENCE — supply raw counts for a valid negative-binomial "
            f"fit."
        )
    else:
        return None, warnings, {
            "refused":    True,
            "error_type": "NonRawCounts",
            "kind":       info["kind"],
            "raw_count_score": info.get("raw_count_score"),
            "confidence": info.get("confidence"),
            "sub_scores": info.get("sub_scores", {}),
            "score_basis": info.get("score_basis", {}),
            "details": (
                f"Count matrix does not look like raw counts "
                f"(kind={info['kind']}, score={info.get('raw_count_score', 0):.2f}, "
                f"max={info['max']:.2f}, "
                f"min={info['min']:.2f}). DESeq2 requires raw integer counts; "
                f"TPM/CPM/FPKM/log-normalized/scaled inputs are invalid. Supply "
                f"a raw-count matrix, or set allow_nonraw_counts=True to coerce "
                f"at low confidence."
            ),
        }

    # Round to integers (required by DESeq2)
    counts = counts.round().astype(int)

    return counts, warnings, {
        "count_source": count_source,
        "kind": info["kind"],
        "raw_count_score": info.get("raw_count_score"),
        "confidence": info.get("confidence"),
        "sub_scores": info.get("sub_scores", {}),
        "score_basis": info.get("score_basis", {}),
    }


def _load_or_infer_metadata(counts, metadata_file: str,
                              design_factor: str) -> tuple:
    """
    Load explicit metadata or infer groups from column names.

    Supported naming patterns:
      ctrl_1, ctrl_2, treat_1, treat_2     → condition = {ctrl, treat}
      WT_rep1, KO_rep1                      → condition = {WT, KO}
      sample_A_1, sample_B_1               → condition = {A, B}
    """
    import pandas as pd
    import re
    warnings = []

    # Try loading explicit metadata file
    if metadata_file and Path(metadata_file).exists():
        try:
            meta = pd.read_csv(metadata_file, sep="\t", index_col=0)
            # Align to count matrix samples
            common = [s for s in counts.columns if s in meta.index]
            if len(common) < 2:
                warnings.append(
                    f"Metadata file has {len(common)} matching samples. "
                    f"Falling back to automatic detection."
                )
            else:
                return meta.loc[common], warnings
        except Exception as e:
            warnings.append(f"Metadata file load failed: {e}. "
                            f"Attempting automatic detection.")

    # Automatic group detection from column names
    samples = list(counts.columns)
    groups  = _infer_groups(samples)

    if groups is None or len(set(groups.values())) < 2:
        return None, warnings + [
            "Could not infer experimental groups from sample names. "
            "Please provide a metadata TSV file with columns: "
            "sample, condition (and optionally: batch, replicate)."
        ]

    meta = pd.DataFrame({
        "sample":      samples,
        design_factor: [groups[s] for s in samples],
    }, index=samples)

    n_groups = len(set(groups.values()))
    warnings.append(
        f"Experimental groups inferred from sample names: "
        f"{dict(set((v,sum(1 for x in groups.values() if x==v)) for v in set(groups.values())))}. "
        f"If this is incorrect, provide a metadata file."
    )

    return meta, warnings


def _infer_groups(samples: list) -> dict | None:
    """
    Detect condition groups from sample names using regex patterns.
    Returns {sample_name: group_label} or None if detection fails.
    """
    import re

    # Pattern 1: condition_replicate (ctrl_1, treat_1, ctrl_2, treat_2)
    p1 = re.compile(r'^([A-Za-z][A-Za-z0-9]+)[_\-](\d+)$')
    matches = {s: m.group(1) for s in samples
               if (m := p1.match(s))}
    if len(matches) == len(samples) and len(set(matches.values())) >= 2:
        return matches

    # Pattern 2: prefix_suffix with letters (WT_rep1, KO_rep1)
    p2 = re.compile(r'^([A-Za-z][A-Za-z0-9]+)[_\-]([Rr]ep\d+|[A-Za-z]\d*)$')
    matches = {s: m.group(1) for s in samples
               if (m := p2.match(s))}
    if len(matches) == len(samples) and len(set(matches.values())) >= 2:
        return matches

    # Pattern 3: split by last underscore, use prefix as group
    if all("_" in s for s in samples):
        groups = {s: "_".join(s.split("_")[:-1]) for s in samples}
        if len(set(groups.values())) >= 2:
            return groups

    # Pattern 4: alphabetical prefix before first digit
    p4 = re.compile(r'^([A-Za-z][A-Za-z0-9\-]*?)(\d.*)$')
    matches = {s: m.group(1).rstrip("_-") for s in samples
               if (m := p4.match(s))}
    if len(matches) == len(samples) and len(set(matches.values())) >= 2:
        return matches

    return None


def _resolve_comparison(metadata, design_factor: str,
                          comparison: dict) -> tuple:
    """
    Resolve which groups to compare.
    P0-5: never infer or substitute numerator/reference levels.
    """
    warnings = []
    groups   = sorted(metadata[design_factor].unique())

    if len(groups) < 2:
        return comparison, [
            f"Only one group found in '{design_factor}': {groups}. "
            f"Cannot run differential expression."
        ]

    num = comparison.get("numerator",   "")
    den = comparison.get("denominator", "")

    if not num or not den:
        warnings.append(
            "Explicit numerator and denominator are required; no comparison "
            "was inferred from group names."
        )
        return {"numerator": "", "denominator": ""}, warnings

    if num not in groups:
        warnings.append(f"Numerator '{num}' not in {groups}.")
    if den not in groups:
        warnings.append(f"Denominator '{den}' not in {groups}.")

    return {"numerator": num, "denominator": den}, warnings


# ── Sample QC ─────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# METHODOLOGY LAYER — explicit normalization & dimensionality reduction
# ══════════════════════════════════════════════════════════════════════════════
#
# Decisions baked into this module (all justified in methods section of report):
#
#   1. DESeq2 DE test        → raw counts (DESeq2 normalizes internally)
#   2. PCA / MDS             → VST + top N variable protein_coding genes
#   3. Heatmap (padj top)    → log2(counts+1) + row z-score
#   4. Heatmap (|log2FC| top)→ log2(counts+1) + row z-score (NEW in v3.8)
#   5. TPM (supplementary)   → gene-length × library-size normalized (NEW)
#
# Why VST over log2(raw+1) + StandardScaler for PCA:
#   - VST is the DESeq2 authors' recommendation (Love, Huber, Anders 2014).
#   - Produces homoscedastic values — equal variance across expression levels.
#   - StandardScaler z-scores per gene over-weight low-variance genes and
#     compress high-variance (biologically meaningful) genes.
#   - log2(raw+1) doesn't correct for library size differences across samples.
#
# Why top N variable protein_coding only for DR:
#   - Pseudogenes/rRNAs dominate raw variance without biological meaning.
#   - Low-variance genes are noise; top 2000 captures informative signal.
#   - Matches the standard DESeq2 vignette workflow.


def _run_vst(counts_raw, metadata, warnings: list):
    """
    Variance-Stabilizing Transformation via pydeseq2.

    Returns a DataFrame (genes × samples) of VST-transformed values.
    Falls back to log2(normed+1) if pydeseq2 VST is unavailable.

    Used for: PCA, MDS, heatmaps at the sample level (NOT for DE testing).
    DESeq2 DE still receives raw counts — it has its own internal normalization.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        import pandas as pd
        import numpy as np

        # pydeseq2 expects samples × genes; use intercept-only design since
        # this is a normalization step, not a test for DE.
        dds = DeseqDataSet(
            counts=counts_raw.T.astype(int),
            metadata=metadata,
            design="~1",
            refit_cooks=False,
            quiet=True,
        )
        dds.fit_size_factors()

        # Try VST first (preferred — fast, handles large datasets).
        # If not available in the installed version, fall back to rlog,
        # then to log2(normalized_counts+1).
        try:
            dds.vst_fit(use_design=False)
            vst = dds.vst_transform()
        except AttributeError:
            try:
                dds.deseq2()
                vst = np.log2(dds.layers["normed_counts"] + 1)
                warnings.append(
                    "pydeseq2 VST unavailable — using log2(size-factor normalized + 1)."
                )
            except Exception as e:
                warnings.append(
                    f"VST/rlog failed ({e}); falling back to raw log2. "
                    f"PCA/MDS may be affected by library size differences."
                )
                vst = np.log2(counts_raw.T.astype(float) + 1)

        # vst is samples × genes; transpose back to genes × samples
        vst_df = pd.DataFrame(
            np.asarray(vst).T,
            index=counts_raw.index,
            columns=counts_raw.columns,
        )
        return vst_df

    except ImportError:
        warnings.append(
            "pydeseq2 not available for VST — using log2(counts+1) + lib-size scaling. "
            "PCA may be dominated by library size differences."
        )
        import pandas as pd
        import numpy as np
        lib = counts_raw.sum(axis=0)
        scale_factor = lib.median() / lib
        normed = counts_raw * scale_factor
        return pd.DataFrame(
            np.log2(normed.astype(float) + 1),
            index=counts_raw.index,
            columns=counts_raw.columns,
        )


def _select_variable_genes(matrix, n_top: int = 2000,
                             biotype_map: dict | None = None,
                             warnings: list | None = None):
    """
    Select the top-N most variable genes from a VST-transformed matrix.

    If biotype_map is provided, restrict to protein_coding genes first.
    Returns (filtered_matrix, n_protein_coding, n_after_variance_filter).

    Args:
        matrix:      DataFrame (genes × samples) — should be VST-transformed.
        n_top:       number of most-variable genes to keep.
        biotype_map: optional {ensembl_id_no_version: biotype_string}.
        warnings:    list to append advisory messages to.
    """
    if warnings is None:
        warnings = []

    # Strip Ensembl version suffix from matrix index (for biotype lookup)
    if biotype_map:
        def _lookup_biotype(gid):
            return biotype_map.get(str(gid).split(".")[0], "unknown")
        biotypes = matrix.index.map(_lookup_biotype)
        n_pc     = int((biotypes == "protein_coding").sum())
        n_total  = len(matrix)
        pc_frac  = n_pc / max(n_total, 1)

        if n_pc < 500:
            warnings.append(
                f"Only {n_pc} protein_coding genes found in matrix "
                f"({pc_frac:.0%} of {n_total}). Falling back to all biotypes "
                f"for DR — results may be affected by pseudogenes/ncRNAs."
            )
        elif pc_frac < 0.70 and biotype_map:
            warnings.append(
                f"GTF annotation has unusually low protein_coding fraction "
                f"({pc_frac:.0%}). Expected ~70-85% for human/mouse."
            )
            matrix_pc = matrix[biotypes == "protein_coding"]
            matrix    = matrix_pc
        else:
            matrix = matrix[biotypes == "protein_coding"]
    else:
        n_pc = 0

    # Top-N most variable
    variance = matrix.var(axis=1)
    n_keep   = min(n_top, len(matrix))
    top_idx  = variance.nlargest(n_keep).index

    return matrix.loc[top_idx], n_pc, n_keep


def _compute_tpm(counts_raw, gene_lengths: dict, warnings: list):
    """
    Compute TPM (Transcripts Per Million) from raw counts.

    TPM = (reads_per_gene / gene_length_kb) / (sum_of_all / 1e6)

    TPM is NOT used by ARIA for DE or PCA (both use better methods).
    It's computed as a supplementary table for downstream tools that
    require TPM input (e.g., ssGSEA, single-sample deconvolution).

    Args:
        counts_raw:   DataFrame (genes × samples), raw integer counts.
        gene_lengths: {ensembl_id: length_bp} — from GTF exon sum.
                     If missing, returns None with a warning.
    """
    try:
        import pandas as pd
        import numpy as np

        if not gene_lengths:
            warnings.append(
                "Gene lengths not available from GTF — cannot compute TPM. "
                "DE analysis is unaffected (uses raw counts)."
            )
            return None

        # Map matrix rows (possibly versioned IDs) to lengths
        def _get_length(gid):
            clean = str(gid).split(".")[0]
            return gene_lengths.get(clean)

        lengths = counts_raw.index.map(_get_length)
        has_len = ~pd.isna(lengths)
        n_with_len = int(has_len.sum())

        if n_with_len < len(counts_raw) * 0.5:
            warnings.append(
                f"TPM: only {n_with_len}/{len(counts_raw)} genes have lengths "
                f"in GTF — TPM table will be incomplete."
            )

        # Drop genes without lengths (can't TPM-normalize them)
        mat      = counts_raw.loc[has_len]
        lens_kb  = pd.Series(
            lengths[has_len].astype(float) / 1000.0,
            index=mat.index,
        )

        # RPK = reads per kilobase per gene
        rpk = mat.div(lens_kb, axis=0)

        # Scaling factor = sum of RPK per sample / 1e6
        scale = rpk.sum(axis=0) / 1e6

        tpm = rpk.div(scale, axis=1)
        tpm = tpm.round(3)
        return tpm

    except Exception as e:
        warnings.append(f"TPM computation failed: {e}")
        return None






# ══════════════════════════════════════════════════════════════════════════════


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

def _build_design_formula(design_factor: str, covariates: list) -> str:
    """Build a DESeq2 design with covariates first and the factor of interest
    last (DESeq2/pydeseq2 convention). Covariates are deduped and never repeat
    the factor of interest. P0-4: honors the confirmed `~ batch + condition`
    design instead of a hardcoded `~ condition`."""
    terms = [c for c in dict.fromkeys(covariates or []) if c and c != design_factor]
    return "~ " + " + ".join(terms + [design_factor])


def _resolve_covariates(metadata, design_factor: str, covariates: list) -> tuple:
    """Keep only covariates usable in this contrast subset: present as a column,
    distinct from the factor of interest, and varying (>=2 non-null levels).
    Returns (usable, dropped) where dropped is a list of (name, reason). A
    confirmed-but-unusable covariate is DISCLOSED, never silently ignored."""
    usable: list = []
    dropped: list = []
    for cov in dict.fromkeys(covariates or []):
        if not cov or cov == design_factor:
            continue
        if cov not in metadata.columns:
            dropped.append((cov, "not present in the sample metadata"))
            continue
        if metadata[cov].dropna().nunique() < 2:
            dropped.append(
                (cov, "constant within the contrast (no levels to adjust for)")
            )
            continue
        usable.append(cov)
    return usable, dropped


def _shrink_coeff(dds, condition_col: str, test_lvl: str):
    """Find the apeGLM LFC coefficient column for test-vs-ref (P1-1 / ADR-023).

    pydeseq2 names coefficients patsy-style, e.g. ``condition[T.treat]``. With the
    dds reference fixed to the contrast's ref level, the only non-reference level
    is ``test_lvl``, so that column is exactly the test-vs-ref effect.
    """
    try:
        cols = [str(c) for c in dds.varm["LFC"].columns]
    except Exception:
        return None
    exact = f"{condition_col}[T.{test_lvl}]"
    if exact in cols:
        return exact
    for c in cols:
        if c.startswith(str(condition_col)) and str(test_lvl) in c:
            return c
    return None


def _run_deseq2(counts, metadata, design_factor: str,
                numerator: str, denominator: str,
                padj_thr: float, lfc_thr: float,
                allow_mock: bool = False,
                min_replicates_per_condition: int = 3,
                covariates: list = None,
                lfc_shrink: bool = True) -> tuple:
    """
    Run DESeq2 via pydeseq2 with correct design factor.
    Returns (result_dict, warnings_list).

    min_replicates_per_condition: production floor. Default 3. n=2 is
    permitted only when the caller explicitly lowers this (e.g. for
    pilot/legacy datasets); the resulting block is flagged with
    low_power_warning in the return dict.
    """
    import pandas as pd
    warnings = []

    # Filter to comparison groups only
    mask      = metadata[design_factor].isin([numerator, denominator])
    meta_sub  = metadata[mask].copy()
    counts_sub = counts[meta_sub.index].copy()

    group_sizes = meta_sub[design_factor].value_counts().to_dict()
    n_num = int(group_sizes.get(numerator, 0))
    n_den = int(group_sizes.get(denominator, 0))

    if n_num < min_replicates_per_condition or n_den < min_replicates_per_condition:
        return {
            "status":     "error",
            "error_type": "InsufficientReplicates",
            "details":    (
                f"DESeq2 requires at least {min_replicates_per_condition} "
                f"replicates per group. Found {group_sizes}. Lower "
                f"min_replicates_per_condition explicitly (e.g. to 2) to run "
                f"with a low_power_warning instead."
            ),
        }, warnings

    # P0-4: honor confirmed covariates (e.g. batch). Keep only those usable in
    # this subset; disclose any confirmed covariate we cannot adjust for.
    usable_covariates, dropped_covariates = _resolve_covariates(
        meta_sub, design_factor, covariates
    )
    for cov, reason in dropped_covariates:
        warnings.append(
            f"Confirmed covariate '{cov}' was NOT adjusted for: {reason}. "
            f"The fitted model does not control for it."
        )

    from aria.utils.design_matrix import validate_design_matrix
    design_check = validate_design_matrix(
        meta_sub,
        condition_col=design_factor,
        covariates=usable_covariates,
        min_replicates_per_condition=min_replicates_per_condition,
    )
    warnings.extend(
        f"Design matrix {issue['severity']}: {issue['message']}"
        for issue in design_check.get("issues", [])
    )
    if design_check.get("status") == "blocking":
        return {
            "status": "error",
            "error_type": "InvalidDesignMatrix",
            "details": "; ".join(
                issue["message"] for issue in design_check.get("issues", [])
                if issue.get("severity") == "blocking"
            )[:500],
            "design_check": design_check,
        }, warnings

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        from aria.utils.power_estimation import bulk_power_estimate
        import warnings as w
        w.filterwarnings("ignore")

        means = counts_sub.mean(axis=1).astype(float)
        variances = counts_sub.var(axis=1, ddof=1).astype(float)
        disp = ((variances - means) / (means ** 2)).replace(
            [float("inf"), float("-inf")], float("nan")
        )
        disp = disp.dropna()
        dispersion_estimate = float(disp[disp > 0].median()) \
            if (disp > 0).any() else 0.1
        mean_expression = float(means[means > 0].median()) \
            if (means > 0).any() else 1.0
        power_estimate = bulk_power_estimate(
            n_per_group=(n_num, n_den),
            mean_expression=mean_expression,
            dispersion=dispersion_estimate,
            target_log2fc=lfc_thr,
            alpha=padj_thr,
        )

        # P0-4: fit the covariate-adjusted design (covariates first, factor of
        # interest last). The contrast still names `design_factor` explicitly.
        design_formula = _build_design_formula(design_factor, usable_covariates)

        # P1-1/ADR-023: fix the reference to the contrast's denominator so the
        # apeGLM coefficient (design_factor[T.numerator]) is exactly the
        # numerator-vs-denominator effect we shrink. Only when shrinkage is on,
        # so the disabled path is byte-identical to the legacy behavior.
        ref_kwargs = {"ref_level": [design_factor, denominator]} if lfc_shrink else {}

        # pydeseq2 expects samples × genes. The public API changed from
        # design_factors=... to design="~ factor"; support both.
        try:
            dds = DeseqDataSet(
                counts=counts_sub.T,
                metadata=meta_sub,
                design=design_formula,
                refit_cooks=True,
                quiet=True,
                **ref_kwargs,
            )
        except TypeError:
            dds = DeseqDataSet(
                counts=counts_sub.T,
                metadata=meta_sub,
                design_factors=usable_covariates + [design_factor],
                refit_cooks=True,
                **ref_kwargs,
            )
        dds.deseq2()

        # P1-1b: put the biological effect-size threshold inside the Wald test
        # instead of filtering |log2FC| after BH. This makes pvalue/padj answer
        # the null |LFC| <= lfc_thr directly, matching DESeq2's lfcThreshold.
        wald_kwargs = {
            "contrast": [design_factor, numerator, denominator],
            "alpha": padj_thr,
            "quiet": True,
        }
        lfc_threshold_test = {
            "requested": float(lfc_thr),
            "applied": False,
            "alt_hypothesis": None,
            "lfc_null": 0.0,
            "significance_rule": "padj < threshold",
        }
        if float(lfc_thr) > 0:
            wald_kwargs.update({
                "lfc_null": float(lfc_thr),
                "alt_hypothesis": "greaterAbs",
            })
            lfc_threshold_test.update({
                "applied": True,
                "alt_hypothesis": "greaterAbs",
                "lfc_null": float(lfc_thr),
            })
        stat_res = DeseqStats(dds, **wald_kwargs)
        stat_res.summary()

        # P1-1/ADR-023: apeGLM LFC shrinkage. Keep the raw MLE LFC, shrink in
        # place, and fall back to raw if the coefficient is unavailable or
        # shrinkage raises (never break the DE). apeGLM leaves p-values unchanged.
        shrink_applied = False
        shrink_reason = None
        raw_lfc = stat_res.results_df["log2FoldChange"].copy()
        if lfc_shrink:
            coeff = _shrink_coeff(dds, design_factor, numerator)
            if coeff is None:
                shrink_reason = "apeGLM coefficient not found"
            else:
                try:
                    stat_res.lfc_shrink(coeff=coeff)
                    shrink_applied = True
                except Exception as _sx:
                    shrink_reason = f"lfc_shrink failed: {str(_sx)[:120]}"
        else:
            shrink_reason = "disabled"

        results_df = stat_res.results_df.dropna(subset=["padj"])
        results_df["log2FoldChange_raw"] = raw_lfc.reindex(results_df.index)

        sig = results_df[results_df["padj"] < padj_thr]
        # Sort by padj ascending (most significant first) so that
        # downstream "top N" slicing returns the most reliable hits,
        # not the most extreme log2FC (which can be noisy at low counts).
        sig = sig.sort_values("padj")

        up_genes   = list(sig[sig["log2FoldChange"] > 0].index)
        down_genes = list(sig[sig["log2FoldChange"] < 0].index)

        if len(sig) == 0:
            warnings.append(
                f"No significant genes found at padj < {padj_thr} "
                f"and |log2FC| > {lfc_thr}. "
                f"Consider relaxing thresholds or checking data quality."
            )

        low_power = (n_num <= 2 or n_den <= 2)
        low_power_reason = (
            f"n={n_num} ({numerator}) vs n={n_den} ({denominator}): "
            f"dispersion estimation is unreliable with fewer than three "
            f"replicates per group. DESeq2 produced results, but effect-size "
            f"estimates and FDR are noisy. Interpret with caution."
        ) if low_power else None
        if low_power:
            warnings.append(low_power_reason)

        return {
            "status":            "success",
            "results":           results_df,
            "n_sig":             len(sig),
            "n_up":              len(up_genes),
            "n_down":            len(down_genes),
            "sig_genes":         list(sig.index),
            "up_genes":          up_genes,
            "down_genes":        down_genes,
            "n_replicates":      {"test": n_num, "ref": n_den},
            "low_power_warning": low_power,
            "low_power_reason":  low_power_reason,
            "dispersion_estimate": dispersion_estimate,
            "mean_expression_estimate": mean_expression,
            "power_estimate_at_lfc_min": power_estimate,
            "design_check":      design_check,
            "fitted_design_formula": design_formula,
            "covariates_adjusted":   usable_covariates,
            "covariates_dropped":    [
                {"covariate": c, "reason": r} for c, r in dropped_covariates
            ],
            "lfc_shrinkage": {
                "requested": bool(lfc_shrink),
                "applied":   shrink_applied,
                "method":    "apeGLM" if shrink_applied else None,
                "reason":    shrink_reason,
            },
            "lfc_threshold_test": lfc_threshold_test,
        }, warnings

    except ImportError:
        if allow_mock:
            warnings.append(
                "pydeseq2 not installed — returning mock DE results "
                "(explicit mock mode)."
            )
            return _mock_de_result(padj_thr, lfc_thr), warnings
        return {
            "status":     "error",
            "error_type": "MissingDependency",
            "details":    (
                "pydeseq2 is required for bulk RNA differential expression. "
                "Install it in the RNA environment or rerun with explicit "
                "allow_mock=true for development only."
            ),
        }, warnings

    except Exception as e:
        return {
            "status":     "error",
            "error_type": "DESeq2Failed",
            "details":    str(e)[:500],
        }, warnings + [f"DESeq2 error: {e}"]


def _mock_de_result(padj_thr, lfc_thr) -> dict:
    """Mock DE result for environments without pydeseq2."""
    import pandas as pd
    import numpy as np
    rng   = np.random.default_rng(42)
    genes = [f"GENE_{i:04d}" for i in range(100)]
    mock_df = pd.DataFrame({
        "log2FoldChange": rng.normal(0, 2, 100),
        "padj":           rng.uniform(0, 0.2, 100),
        "pvalue":         rng.uniform(0, 0.1, 100),
        "baseMean":       rng.exponential(500, 100),
    }, index=genes)

    sig  = mock_df[(mock_df["padj"] < padj_thr) &
                   (mock_df["log2FoldChange"].abs() > lfc_thr)]
    up   = list(sig[sig["log2FoldChange"] > 0].index)
    down = list(sig[sig["log2FoldChange"] < 0].index)

    return {
        "status":    "success",
        "results":   mock_df,
        "n_sig":     len(sig),
        "n_up":      len(up),
        "n_down":    len(down),
        "sig_genes": list(sig.index),
        "up_genes":  up,
        "down_genes": down,
        "note":      "mock — pydeseq2 not installed",
    }


# ── Pathway enrichment ────────────────────────────────────────────────────────























# ── Visualizations ────────────────────────────────────────────────────────────







if __name__ == "__main__":
    run_script(bulk_rna_de)
