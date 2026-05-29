"""
ARIA Bulk RNA-seq Differential Expression Script
-------------------------------------------------
Full bulk RNA-seq pipeline executed inside aria-rna-env by EnvironmentManager.

Fixes vs old inline implementation:
  1. Design factor extracted from biological intent (not hardcoded "sample")
  2. Robust metadata parsing: detects groups from column names automatically
  3. Sample outlier detection (PCA-based) before running DESeq2
  4. Pathway enrichment via gseapy after DE (GO BP, KEGG, Reactome)
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
from pathlib import Path

# ── Paper theme — applied to every plot in this module ────────────────────────
_P = dict(
    fig     = "white",
    ax      = "white",
    text    = "#1e293b",
    muted   = "#475569",
    dim     = "#94a3b8",
    border  = "#e2e8f0",
    title   = "#0f172a",
    up      = "#dc2626",   # upregulated  (red  — convention)
    down    = "#2563eb",   # downregulated (blue — convention)
    ns      = "#d1d5db",   # non-significant (light gray)
    ref     = "#94a3b8",   # threshold reference lines
    annot   = "#374151",   # gene label annotations
    palette = ["#1d4ed8", "#dc2626", "#059669", "#d97706",
               "#7c3aed", "#db2777", "#0891b2"],
)
# ──────────────────────────────────────────────────────────────────────────────


def bulk_rna_de(params: dict) -> dict:
    from pathlib import Path
    import numpy as np
    import warnings as warn_mod
    warn_mod.filterwarnings("ignore")

    files          = params.get("files", [])
    metadata_file  = params.get("metadata_file", "")
    design_factor  = params.get("design_factor", "condition")

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
    counts, load_warn = _load_counts(files)
    warnings.extend(load_warn)
    if counts is None:
        return {
            "status":     "error",
            "error_type": "CountsLoadFailed",
            "details":    "Could not load count matrix from provided files.",
        }

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
    outlier_samples = _prune_outliers_for_design(
        sample_qc.get("outliers", []), metadata, design_factor, warnings
    )
    sample_qc["outliers"] = outlier_samples
    if outlier_samples:
        counts   = counts.drop(columns=outlier_samples, errors="ignore")
        metadata = metadata.drop(index=outlier_samples, errors="ignore")
        warnings.append(
            f"Sample outliers removed: {outlier_samples}. "
            f"Clustered away from group centroid in PCA."
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

    # ── 5. If no contrasts provided, auto-generate from groups ──────────
    if not contrasts_in:
        contrasts_in, auto_warn = _auto_contrasts(
            metadata, design_factor
        )
        warnings.extend(auto_warn)
        if not contrasts_in:
            return {
                "status":     "error",
                "error_type": "NoContrastsBuildable",
                "details":    f"Could not build contrasts from "
                              f"{design_factor}={list(metadata[design_factor].unique())}",
            }

    # ── 6. Run each contrast ─────────────────────────────────────────────
    contrast_results = []

    for contrast in contrasts_in:
        num  = contrast.get("numerator",   "")
        den  = contrast.get("denominator", "")
        name = contrast.get("name", f"{num} vs {den}")

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
        if run_pathways and de_result.get("sig_genes"):
            pathways, pw_warn = _run_pathway_enrichment(
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
            "top_genes":         top_genes,
            # Full DE list for cross-contrast overlap (in symbols when available,
            # else Ensembl IDs — both work for set intersection)
            "all_sig_genes":     all_sig_symbols,
            "pathways":          pathways,
            "pathway_background": {
                "background_size": len(background_symbols),
                "background_source": "dataset_expressed_genes",
            },
            "plots":             plots,
            "contrast_dir":      str(contrast_dir),
            "figures_dir":       str(figures_dir),
            "tables_dir":        str(tables_dir),
        })

    if not contrast_results:
        return {
            "status":     "error",
            "error_type": "AllContrastsFailed",
            "details":    "No contrasts produced valid results.",
            "warnings":   warnings,
        }

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
                "normalization":  "Enrichr Fisher's exact test",
                "gene_filter":    "padj<0.05 and |log2FC|>threshold; ORA universe = genes retained after expression filtering",
                "justification":  "Standard over-representation for discrete gene lists using the dataset-expressed background rather than Enrichr's default all-database universe.",
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
        "design_used":      f"~{design_factor}",
        "padj_threshold":   padj_thr,
        "lfc_threshold":    lfc_thr,
        "overlap":          overlap_info,
        "methodology":      methodology,
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


def _auto_contrasts(metadata, design_factor: str) -> tuple:
    """Generate contrasts when none provided — each group vs control."""
    warnings = []
    groups = sorted(metadata[design_factor].unique())

    if len(groups) < 2:
        return [], [
            f"Only one group in '{design_factor}': {groups}. "
            f"Cannot run differential expression."
        ]

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
        # Alphabetical pairwise (only first vs rest to limit explosion)
        ref = groups[0]
        for g in groups[1:]:
            contrasts.append({
                "numerator":   g,
                "denominator": ref,
                "name":        f"{g} vs {ref}",
            })
        warnings.append(
            f"No control group identified in {groups}. "
            f"Using '{ref}' as reference."
        )

    return contrasts, warnings


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

def _load_counts(files: list) -> tuple:
    """
    Load counts matrix from various formats.
    Returns (DataFrame, warnings_list).
    genes × samples orientation enforced.
    """
    import pandas as pd
    warnings = []

    valid = [f for f in files if Path(f).exists()]
    if not valid:
        return None, ["No valid count files found."]

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
                return None, ["Count matrix has no numeric columns."]
        except Exception as e:
            return None, [f"Failed to load {count_files[0]}: {e}"]
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
            return None, ["Could not load any count files."]
        counts = pd.concat(frames, axis=1).fillna(0)

    # Ensure genes × samples (more genes than samples in typical experiments)
    if counts.shape[1] > counts.shape[0]:
        warnings.append(
            f"Transposing matrix: detected {counts.shape[1]} rows × "
            f"{counts.shape[0]} cols — expected genes as rows."
        )
        counts = counts.T

    # Round to integers (required by DESeq2)
    counts = counts.round().astype(int)

    return counts, warnings


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
    If comparison not specified, infer from available groups.
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

    if num and den:
        if num not in groups:
            warnings.append(
                f"Numerator '{num}' not in {groups}. "
                f"Using '{groups[-1]}' instead."
            )
            num = groups[-1]
        if den not in groups:
            warnings.append(
                f"Denominator '{den}' not in {groups}. "
                f"Using '{groups[0]}' instead."
            )
            den = groups[0]
    else:
        # Infer: common patterns
        TREATED_KEYWORDS = {"treat", "treated", "mut", "mutant", "ko",
                             "knockdown", "kd", "overexpression", "oe",
                             "infected", "stimulated", "high", "disease"}
        CTRL_KEYWORDS    = {"ctrl", "control", "wt", "wildtype",
                             "scramble", "vehicle", "untreated",
                             "healthy", "low", "normal"}

        def score(g, keywords):
            return sum(1 for k in keywords if k in g.lower())

        ctrl_scores = {g: score(g, CTRL_KEYWORDS)   for g in groups}
        trt_scores  = {g: score(g, TREATED_KEYWORDS) for g in groups}

        best_ctrl = max(ctrl_scores, key=ctrl_scores.get)
        best_trt  = max(trt_scores,  key=trt_scores.get)

        if best_ctrl == best_trt or (ctrl_scores[best_ctrl] == 0
                                     and trt_scores[best_trt] == 0):
            # Fallback: use alphabetical order
            den, num = sorted(groups)[:2]
            warnings.append(
                f"Could not infer comparison direction from group names {groups}. "
                f"Using alphabetical order: {num} vs {den}. "
                f"Specify comparison explicitly if this is wrong."
            )
        else:
            num, den = best_trt, best_ctrl
            warnings.append(
                f"Comparison inferred: {num} (treated) vs {den} (control). "
                f"Verify this matches your experimental design."
            )

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


def _plot_pca_mds(vst_variable, metadata, output_dir: str,
                    warnings: list) -> tuple:
    """
    Generate PCA + MDS plots from a VST-transformed, variable-filtered matrix.

    Both plots operate on the SAME data matrix → apples-to-apples comparison.
    PCA: linear orthogonal components capturing variance.
    MDS: preserves Euclidean distances (non-linear relationships).

    Returns (pca_path, mds_path, pca_coords, pca_variance_pct).
    """
    try:
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
        from sklearn.manifold import MDS
        from scipy.spatial.distance import pdist, squareform

        # Samples × genes for PCA/MDS (standard orientation)
        X = vst_variable.T.values

        # Identify the condition column
        condition_col = None
        for c in metadata.columns:
            if c not in ("sample", "batch", "replicate"):
                condition_col = c
                break
        conditions = (metadata[condition_col].astype(str)
                       if condition_col else
                       pd.Series(["all"] * X.shape[0], index=metadata.index))

        # ── PCA ──────────────────────────────────────────────────────────
        # NO StandardScaler — VST already stabilizes variance.
        # PCA on genes × samples with centered (not z-scored) values.
        pca     = PCA(n_components=min(5, X.shape[0] - 1))
        coords  = pca.fit_transform(X)
        var_pct = pca.explained_variance_ratio_ * 100

        # ── MDS (classical / metric) ─────────────────────────────────────
        # Euclidean distance on VST — same convention as edgeR's plotMDS.
        dist_matrix = squareform(pdist(X, metric="euclidean"))
        # Use future defaults explicitly to avoid FutureWarnings on sklearn 1.9+
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore", FutureWarning)
            try:
                # sklearn 1.9+ signature
                mds = MDS(n_components=2, metric=True,
                           normalized_stress="auto", random_state=42,
                           n_init=4, dissimilarity="precomputed")
            except TypeError:
                mds = MDS(n_components=2, dissimilarity="precomputed",
                           random_state=42)
            mds_coords = mds.fit_transform(dist_matrix)

        # ── Plot both side by side in the dark theme ─────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.patch.set_facecolor(_P["fig"])

        groups = sorted(conditions.unique())
        color_map = {g: _P["palette"][i % len(_P["palette"])]
                     for i, g in enumerate(groups)}

        for ax, coords_plot, title, xlbl, ylbl in [
            (axes[0], coords[:, :2], "PCA (VST, top variable PC genes)",
             f"PC1 ({var_pct[0]:.1f}%)", f"PC2 ({var_pct[1]:.1f}%)"),
            (axes[1], mds_coords,    "MDS (VST, Euclidean distance)",
             "MDS1", "MDS2"),
        ]:
            ax.set_facecolor(_P["ax"])
            for g in groups:
                mask = (conditions.values == g)
                ax.scatter(coords_plot[mask, 0], coords_plot[mask, 1],
                           s=120, c=color_map[g], label=g,
                           edgecolor="white", linewidth=0.8, alpha=0.9)
            # Sample labels
            for i, s in enumerate(metadata.index):
                ax.annotate(s, (coords_plot[i, 0], coords_plot[i, 1]),
                            fontsize=8, color=_P["annot"],
                            xytext=(5, 5), textcoords="offset points")
            ax.set_title(title, color=_P["title"], fontsize=11, fontweight="bold")
            ax.set_xlabel(xlbl, color=_P["muted"])
            ax.set_ylabel(ylbl, color=_P["muted"])
            ax.tick_params(colors=_P["muted"])
            ax.legend(facecolor=_P["fig"], edgecolor=_P["border"],
                      labelcolor=_P["text"], fontsize=9, loc="best")
            for spine in ax.spines.values():
                spine.set_edgecolor(_P["border"])

        plt.tight_layout()
        combined_path = str(Path(output_dir) / "pca_mds.svg")
        plt.savefig(combined_path, format="svg",
                     facecolor=_P["fig"])
        plt.close(fig)

        # Also save individual SVGs for manuscript use
        pca_path = _save_single_dr_plot(
            coords[:, :2], conditions, metadata.index, color_map,
            f"PCA (VST, top variable PC genes)",
            f"PC1 ({var_pct[0]:.1f}%)", f"PC2 ({var_pct[1]:.1f}%)",
            str(Path(output_dir) / "pca.svg"),
        )
        mds_path = _save_single_dr_plot(
            mds_coords, conditions, metadata.index, color_map,
            "MDS (VST, Euclidean distance)",
            "MDS1", "MDS2",
            str(Path(output_dir) / "mds.svg"),
        )

        return pca_path, mds_path, coords, [round(float(v), 3)
                                               for v in var_pct[:2]]

    except Exception as e:
        warnings.append(f"PCA/MDS plotting failed: {e}")
        return None, None, None, []


def _save_single_dr_plot(coords, conditions, sample_names, color_map,
                          title, xlbl, ylbl, output_path):
    """Helper — single-axis DR plot for manuscript figures."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 5))
        fig.patch.set_facecolor(_P["fig"])
        ax.set_facecolor(_P["ax"])
        for g in sorted(color_map.keys()):
            mask = (conditions.values == g)
            ax.scatter(coords[mask, 0], coords[mask, 1],
                       s=140, c=color_map[g], label=g,
                       edgecolor="white", linewidth=0.8, alpha=0.9)
        for i, s in enumerate(sample_names):
            ax.annotate(s, (coords[i, 0], coords[i, 1]),
                        fontsize=8, color=_P["annot"],
                        xytext=(5, 5), textcoords="offset points")
        ax.set_title(title, color=_P["title"], fontsize=12, fontweight="bold")
        ax.set_xlabel(xlbl, color=_P["muted"])
        ax.set_ylabel(ylbl, color=_P["muted"])
        ax.tick_params(colors=_P["muted"])
        ax.legend(facecolor=_P["fig"], edgecolor=_P["border"],
                  labelcolor=_P["text"], fontsize=9, loc="best")
        for spine in ax.spines.values():
            spine.set_edgecolor(_P["border"])
        plt.tight_layout()
        plt.savefig(output_path, format="svg",
                     facecolor=_P["fig"])
        plt.close(fig)
        return output_path
    except Exception:
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
    Outliers are removed from downstream DE.
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
        return {"n_samples": int(counts.shape[1]), "outliers": [],
                "pca_variance": [], "error": str(e)}


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


# ── DESeq2 ────────────────────────────────────────────────────────────────────

def _run_deseq2(counts, metadata, design_factor: str,
                numerator: str, denominator: str,
                padj_thr: float, lfc_thr: float,
                allow_mock: bool = False,
                min_replicates_per_condition: int = 3) -> tuple:
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

    from aria.utils.design_matrix import validate_design_matrix
    design_check = validate_design_matrix(
        meta_sub,
        condition_col=design_factor,
        covariates=[],
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

        # pydeseq2 expects samples × genes. The public API changed from
        # design_factors=... to design="~ factor"; support both.
        try:
            dds = DeseqDataSet(
                counts=counts_sub.T,
                metadata=meta_sub,
                design=f"~ {design_factor}",
                refit_cooks=True,
                quiet=True,
            )
        except TypeError:
            dds = DeseqDataSet(
                counts=counts_sub.T,
                metadata=meta_sub,
                design_factors=design_factor,
                refit_cooks=True,
            )
        dds.deseq2()

        stat_res = DeseqStats(
            dds,
            contrast=[design_factor, numerator, denominator],
        )
        stat_res.summary()

        results_df = stat_res.results_df.dropna(subset=["padj"])

        sig = results_df[
            (results_df["padj"]             < padj_thr) &
            (results_df["log2FoldChange"].abs() > lfc_thr)
        ]
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

def _load_symbol_map(files: list, warnings: list,
                       gtf_hint: str = None) -> dict:
    """
    Load Ensembl-ID → HGNC-symbol mapping.

    Strategy:
      1. Look for counts_with_symbols.tsv next to the counts file.
      2. If missing, auto-regenerate from a GTF (search standard locations
         + hint provided by the agent).
      3. If still no GTF found, fall back to {} and warn loudly.

    The auto-regeneration matters because resume logic (v3.1+) skips
    featureCounts when counts_matrix.tsv exists — but if the matrix was
    written by an older version of ARIA (before v3.3), the symbols file
    won't exist. This function patches that gap without re-running
    featureCounts.
    """
    if not files:
        return {}
    try:
        import pandas as pd
        counts_path = Path(files[0])
        sym_path    = counts_path.parent / "counts_with_symbols.tsv"

        # ── Path A: symbols file already exists ──────────────────────
        if sym_path.exists():
            df = pd.read_csv(sym_path, sep="\t", index_col=0)
            if "gene_symbol" not in df.columns:
                warnings.append(
                    "counts_with_symbols.tsv missing 'gene_symbol' "
                    "column — regenerating."
                )
            else:
                clean = {}
                for k, v in df["gene_symbol"].items():
                    if isinstance(k, str) and isinstance(v, str):
                        clean[k.split(".")[0]] = v
                return clean

        # ── Path B: regenerate from GTF ──────────────────────────────
        gtf_path = _locate_gtf(counts_path, gtf_hint)
        if not gtf_path:
            warnings.append(
                "counts_with_symbols.tsv missing AND no GTF found in "
                "standard locations (~/.aria/genomes/*/annotation.gtf*). "
                "Pathway enrichment will use Ensembl IDs (likely 0 matches)."
            )
            return {}

        warnings.append(
            f"Auto-generating gene-symbol map from GTF: {gtf_path.name}"
        )
        mapping = _gtf_to_symbol_map(str(gtf_path))
        if not mapping:
            warnings.append(
                "GTF parse returned 0 mappings — check gene_name "
                "annotations are present."
            )
            return {}

        # Persist the regenerated symbols file for next time
        try:
            counts_df = pd.read_csv(counts_path, sep="\t", index_col=0)
            counts_with_sym = counts_df.copy()
            counts_with_sym.insert(
                0, "gene_symbol",
                [mapping.get(str(g).split(".")[0], str(g))
                 for g in counts_with_sym.index]
            )
            counts_with_sym.to_csv(str(sym_path), sep="\t")
            warnings.append(
                f"Wrote {sym_path.name} ({len(mapping):,} symbol mappings) "
                f"for future runs."
            )
        except Exception as e:
            warnings.append(f"Could not persist symbols file: {e}")

        return mapping

    except Exception as e:
        warnings.append(f"Symbol map load failed: {e}")
        return {}


def _locate_gtf(counts_path: Path, hint: str = None) -> Path | None:
    """
    Find a GTF file in standard locations.

    Search order:
      1. Explicit hint from the caller
      2. Sibling of the counts file (~/.../counts/annotation.gtf)
      3. ~/.aria/genomes/*/annotation.gtf{,.gz}
      4. ~/.aria/genomes/hg38/annotation.gtf etc.
    """
    candidates = []
    if hint:
        candidates.append(Path(hint))

    # Sibling of counts file
    candidates.extend([
        counts_path.parent / "annotation.gtf",
        counts_path.parent / "annotation.gtf.gz",
        counts_path.parent.parent / "annotation.gtf",
    ])

    # ~/.aria/genomes/*/annotation.gtf
    aria_genomes = Path.home() / ".aria" / "genomes"
    if aria_genomes.exists():
        for genome_dir in aria_genomes.iterdir():
            if genome_dir.is_dir():
                candidates.extend([
                    genome_dir / "annotation.gtf",
                    genome_dir / "annotation.gtf.gz",
                ])

    for c in candidates:
        if c.exists() and c.stat().st_size > 100_000:
            return c
    return None


def _gtf_to_symbol_map(gtf_file: str) -> dict:
    """
    Parse a GTF for {ensembl_id (no version): gene_symbol}.
    Streaming, low memory. Same logic as rna_quantify._build_ensembl_to_symbol_map.
    """
    import gzip as _gzip
    import re as _re

    mapping = {}
    opener = _gzip.open if str(gtf_file).endswith(".gz") else open
    try:
        with opener(gtf_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attrs = fields[8]
                gid_match = _re.search(r'gene_id "([^"]+)"', attrs)
                sym_match = _re.search(r'gene_name "([^"]+)"', attrs)
                if gid_match and sym_match:
                    gid = gid_match.group(1).split(".")[0]
                    mapping[gid] = sym_match.group(1)
    except Exception:
        return {}

    return mapping


def _load_gene_annotation(files: list, warnings: list) -> dict:
    """
    Load {biotype_map, length_map} from the GTF — used for DR filtering
    (protein_coding only) and TPM computation.

    biotype_map: {ensembl_id_no_version: biotype_string}
    length_map:  {ensembl_id_no_version: union_exon_length_bp}

    Exon-union length is the gene-model-appropriate length for TPM
    (sums non-overlapping exon regions). Falls back to gene coordinate
    length if exon parsing fails.

    Returns empty dicts gracefully if GTF not locatable.
    """
    if not files:
        return {"biotype": {}, "length": {}}
    try:
        counts_path = Path(files[0])
        gtf_path    = _locate_gtf(counts_path)
        if not gtf_path:
            warnings.append(
                "GTF not found for biotype/length annotation — "
                "DR will use all biotypes; TPM will not be computed."
            )
            return {"biotype": {}, "length": {}}

        biotype_map, length_map = _parse_gtf_biotype_and_length(str(gtf_path))
        warnings.append(
            f"Gene annotation loaded from {gtf_path.name}: "
            f"{len(biotype_map):,} biotypes, {len(length_map):,} lengths."
        )
        return {"biotype": biotype_map, "length": length_map}
    except Exception as e:
        warnings.append(f"Gene annotation load failed (non-fatal): {e}")
        return {"biotype": {}, "length": {}}


def _parse_gtf_biotype_and_length(gtf_file: str) -> tuple:
    """
    Single pass through GTF to extract both biotype (from gene lines)
    and exon-union length (sum of non-overlapping exon spans per gene).

    Returns (biotype_map, length_map).
    """
    import gzip as _gzip
    import re as _re
    from collections import defaultdict

    biotype_map = {}
    exons_by_gene = defaultdict(list)  # gid -> list of (start, end)

    opener = _gzip.open if str(gtf_file).endswith(".gz") else open
    try:
        with opener(gtf_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                feature = fields[2]
                if feature not in ("gene", "exon"):
                    continue

                attrs = fields[8]
                gid_match = _re.search(r'gene_id "([^"]+)"', attrs)
                if not gid_match:
                    continue
                gid = gid_match.group(1).split(".")[0]

                if feature == "gene":
                    bio_match = (_re.search(r'gene_biotype "([^"]+)"', attrs)
                                  or _re.search(r'gene_type "([^"]+)"', attrs))
                    if bio_match:
                        biotype_map[gid] = bio_match.group(1)

                elif feature == "exon":
                    try:
                        start = int(fields[3])
                        end   = int(fields[4])
                        exons_by_gene[gid].append((start, end))
                    except ValueError:
                        continue
    except Exception:
        return biotype_map, {}

    # Compute union-exon length per gene (merge overlapping intervals)
    length_map = {}
    for gid, intervals in exons_by_gene.items():
        if not intervals:
            continue
        # Merge overlapping intervals
        intervals.sort()
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        length_map[gid] = sum(e - s + 1 for s, e in merged)

    return biotype_map, length_map


def _to_symbols(gene_ids: list, symbol_map: dict) -> list:
    """
    Convert Ensembl IDs to HGNC symbols. Genes without a symbol are
    dropped (Enrichr won't match them anyway). Returns deduplicated list.
    """
    if not symbol_map:
        return list(gene_ids)
    seen = set()
    out  = []
    for gid in gene_ids:
        if not isinstance(gid, str):
            continue
        clean = gid.split(".")[0]   # strip version
        sym = symbol_map.get(clean)
        # If not in map, check if input might already be a symbol
        if not sym:
            if not clean.startswith(("ENSG", "ENSMUSG", "ENS")):
                sym = clean   # already a symbol
            else:
                continue       # unmapped Ensembl ID — skip
        if sym not in seen:
            seen.add(sym)
            out.append(sym)
    return out


def _run_pathway_enrichment(sig_genes: list, up_genes: list,
                              down_genes: list, organism: str,
                              output_dir: str,
                              symbol_map: dict = None,
                              background_genes: list | None = None,
                              allow_mock: bool = False) -> tuple:
    """
    Run pathway enrichment via gseapy.
    Tests GO Biological Process, KEGG, and Reactome.

    If symbol_map is provided, Ensembl IDs are converted to HGNC symbols
    before sending to Enrichr (which requires symbols).
    """
    warnings  = []
    pathways  = {}
    gene_sets = _get_gene_sets(organism)

    if not sig_genes:
        return {}, ["No significant genes for pathway enrichment."]

    # Convert IDs → symbols if we have a mapping
    sig_symbols  = _to_symbols(sig_genes,  symbol_map or {})
    up_symbols   = _to_symbols(up_genes,   symbol_map or {})
    down_symbols = _to_symbols(down_genes, symbol_map or {})
    background_symbols = sorted({
        str(g) for g in (background_genes or [])
        if g and str(g).lower() != "nan"
    })
    background_set = set(background_symbols)
    if background_set:
        sig_symbols = [g for g in sig_symbols if g in background_set]
        up_symbols = [g for g in up_symbols if g in background_set]
        down_symbols = [g for g in down_symbols if g in background_set]

    if symbol_map and len(sig_symbols) < len(sig_genes) * 0.3:
        warnings.append(
            f"Symbol mapping yielded only {len(sig_symbols)}/{len(sig_genes)} "
            f"genes (< 30%). Pathway enrichment may be underpowered."
        )

    if not sig_symbols:
        warnings.append(
            "No valid gene symbols to enrich — check that GTF gene_name "
            "annotations are present."
        )
        return {}, warnings

    try:
        import gseapy as gp
        import time as _time

        for i, (db, db_name) in enumerate(gene_sets.items()):
            # Rate-limit: Enrichr blocks aggressive callers. Pattern from
            # production scripts: ~5-15 sec between calls (we use 8 — enough
            # to avoid throttling but tolerable for 3 databases).
            if i > 0:
                _time.sleep(8)

            try:
                enrichr_kwargs = {
                    "gene_list": sig_symbols[:500],
                    "gene_sets": db_name,
                    "organism": _gseapy_organism(organism),
                    "outdir": None,
                    "verbose": False,
                }
                if background_symbols:
                    enrichr_kwargs["background"] = background_symbols
                try:
                    enr = gp.enrichr(**enrichr_kwargs)
                except TypeError as exc:
                    if "background" not in str(exc):
                        raise
                    enrichr_kwargs.pop("background", None)
                    enr = gp.enrichr(**enrichr_kwargs)
                results = enr.results

                if results is not None and not results.empty:
                    sig_pw = results[results["Adjusted P-value"] < 0.05]
                    sig_pw = sig_pw.sort_values("Adjusted P-value")

                    pathways[db] = [
                        {
                            "term":    row["Term"],
                            "padj":    round(float(row["Adjusted P-value"]), 5),
                            "overlap": row.get("Overlap", ""),
                            "odds_ratio": round(float(row.get("Odds Ratio", 1)), 2),
                            "combined_score": round(float(row.get("Combined Score", 0)), 1),
                            "genes":   row.get("Genes", "").split(";")[:10],
                        }
                        for _, row in sig_pw.head(20).iterrows()
                    ]

                    if not pathways[db]:
                        warnings.append(
                            f"No significant {db} pathways at padj < 0.05 "
                            f"({len(results)} terms tested)."
                        )
                else:
                    warnings.append(
                        f"{db} returned empty results from Enrichr "
                        f"(0 of {len(sig_symbols)} symbols matched)."
                    )
            except Exception as e:
                # Capture the REAL error so the LLM doesn't need to
                # invent a cause (and bad explanations won't poison
                # the report).
                err_msg = str(e)[:200]
                warnings.append(
                    f"{db} enrichment ERROR: {err_msg}. "
                    f"(Sent {len(sig_symbols)} symbols, organism="
                    f"{_gseapy_organism(organism)}.)"
                )

    except ImportError:
        if allow_mock:
            warnings.append(
                "gseapy not installed — returning mock pathway enrichment "
                "(explicit mock mode)."
            )
            pathways = _mock_pathways(sig_genes)
        else:
            warnings.append(
                "gseapy not installed — pathway enrichment skipped. "
                "No simulated pathway terms were generated."
            )

    # Directional enrichment: up vs down separately (uses symbols)
    if up_symbols and down_symbols:
        try:
            import gseapy as gp
            for direction, genes in [("up", up_symbols[:200]),
                                      ("down", down_symbols[:200])]:
                enrichr_kwargs = {
                    "gene_list": genes,
                    "gene_sets": "KEGG_2021_Human"
                                 if "sapiens" in organism.lower()
                                 else "KEGG_2019_Mouse",
                    "organism": _gseapy_organism(organism),
                    "outdir": None,
                    "verbose": False,
                }
                if background_symbols:
                    enrichr_kwargs["background"] = background_symbols
                try:
                    enr = gp.enrichr(**enrichr_kwargs)
                except TypeError as exc:
                    if "background" not in str(exc):
                        raise
                    enrichr_kwargs.pop("background", None)
                    enr = gp.enrichr(**enrichr_kwargs)
                if enr.results is not None and not enr.results.empty:
                    sig_pw = enr.results[
                        enr.results["Adjusted P-value"] < 0.1
                    ].head(10)
                    pathways[f"KEGG_{direction}"] = [
                        {"term": row["Term"],
                         "padj": round(float(row["Adjusted P-value"]), 5)}
                        for _, row in sig_pw.iterrows()
                    ]
        except Exception:
            pass

    return pathways, warnings


def _get_gene_sets(organism: str) -> dict:
    """Return gene set databases appropriate for organism."""
    if "sapiens" in organism.lower():
        return {
            "GO_BP":    "GO_Biological_Process_2021",
            "KEGG":     "KEGG_2021_Human",
            "Reactome": "Reactome_2022",
        }
    elif "musculus" in organism.lower():
        return {
            "GO_BP":    "GO_Biological_Process_2021",
            "KEGG":     "KEGG_2019_Mouse",
            "Reactome": "Reactome_2022",
        }
    else:
        return {"GO_BP": "GO_Biological_Process_2021"}


def _gseapy_organism(organism: str) -> str:
    """
    Return the organism string in the format gseapy/Enrichr expects.

    gseapy 1.1.13+ explicitly validates against a lowercase whitelist:
        human, hsapiens, h. sapiens, hs, mouse, mus musculus, ...
    Capitalized strings like "Human" or "Homo sapiens" are REJECTED with:
        "Invalid organism 'Human'. Valid options are: ..."
    """
    s = organism.lower()
    if "sapiens" in s or s in ("human", "hs", "hsapiens", "h. sapiens"):
        return "human"
    if "musculus" in s or s in ("mouse", "mm"):
        return "mouse"
    if "rerio" in s or s in ("zebrafish", "danio rerio"):
        return "fish"
    if "melanogaster" in s or s in ("fly", "drosophila"):
        return "fly"
    if "elegans" in s:
        return "worm"
    if "cerevisiae" in s or s == "yeast":
        return "yeast"
    return "human"   # safest default for the most common use case


def _mock_pathways(sig_genes: list) -> dict:
    """Mock pathways when gseapy unavailable."""
    return {
        "GO_BP": [
            {"term": "immune response", "padj": 0.001,
             "genes": sig_genes[:5], "overlap": "12/234"},
            {"term": "inflammatory response", "padj": 0.003,
             "genes": sig_genes[:3], "overlap": "8/156"},
        ],
        "KEGG": [
            {"term": "T cell receptor signaling pathway", "padj": 0.005,
             "genes": sig_genes[:4], "overlap": "6/87"},
        ],
        "note": "mock — install gseapy for real enrichment",
    }


# ── Visualizations ────────────────────────────────────────────────────────────

def _generate_plots(de_result: dict, sample_qc: dict,
                    counts_filt, metadata, design_factor: str,
                    output_dir: str, padj_thr: float,
                    lfc_thr: float,
                    title_suffix: str = "") -> dict:
    """Generate volcano, heatmap, and sample PCA plots."""
    plots = {}

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        results_df = de_result.get("results")
        if results_df is not None and not results_df.empty:
            # ── Volcano plot ───────────────────────────────────────────
            fig, ax = plt.subplots(figsize=(8, 6))
            fig.patch.set_facecolor(_P["fig"])
            ax.set_facecolor(_P["ax"])

            lfc = results_df["log2FoldChange"].clip(-8, 8)
            neg_log_p = -np.log10(results_df["padj"].clip(1e-300, 1) + 1e-300)

            # Color by significance
            sig_mask = (
                (results_df["padj"] < padj_thr) &
                (results_df["log2FoldChange"].abs() > lfc_thr)
            )
            up_mask   = sig_mask & (results_df["log2FoldChange"] > 0)
            down_mask = sig_mask & (results_df["log2FoldChange"] < 0)

            ax.scatter(lfc[~sig_mask], neg_log_p[~sig_mask],
                       color=_P["ns"], alpha=0.5, s=8, linewidths=0)
            ax.scatter(lfc[up_mask], neg_log_p[up_mask],
                       color=_P["up"], alpha=0.8, s=12, linewidths=0)
            ax.scatter(lfc[down_mask], neg_log_p[down_mask],
                       color=_P["down"], alpha=0.8, s=12, linewidths=0)

            # Reference lines
            ax.axhline(-np.log10(padj_thr), color=_P["ref"],
                       linestyle="--", alpha=0.6, linewidth=0.8)
            ax.axvline( lfc_thr, color=_P["ref"],
                       linestyle="--", alpha=0.6, linewidth=0.8)
            ax.axvline(-lfc_thr, color=_P["ref"],
                       linestyle="--", alpha=0.6, linewidth=0.8)

            # Labels for top genes
            top_sig = results_df[sig_mask].nsmallest(10, "padj")
            for gene, row in top_sig.iterrows():
                ax.annotate(
                    str(gene)[:12],
                    xy=(float(row["log2FoldChange"]),
                        float(-np.log10(row["padj"] + 1e-300))),
                    xytext=(3, 3), textcoords="offset points",
                    fontsize=6, color=_P["annot"], alpha=0.9,
                )

            ax.set_xlabel("log₂ Fold Change", color=_P["muted"])
            ax.set_ylabel("-log₁₀ adjusted p-value", color=_P["muted"])
            n_up   = de_result.get("n_up",   0)
            n_down = de_result.get("n_down", 0)
            title_text = (
                f"Differential Expression"
                + (f" — {title_suffix}" if title_suffix else "")
                + f"\n{n_up} up (red)  {n_down} down (blue)"
            )
            ax.set_title(title_text, color=_P["title"], fontsize=11,
                         fontweight="bold")
            ax.tick_params(colors=_P["muted"])
            for spine in ax.spines.values():
                spine.set_edgecolor(_P["border"])

            volcano_path = str(Path(output_dir) / "volcano.svg")
            plt.tight_layout()
            plt.savefig(volcano_path, format="svg",
                        facecolor=_P["fig"])
            plt.close()
            plots["volcano"] = volcano_path

        # ── Heatmaps of top DE genes (two complementary views) ──────────
        # (a) Top 50 by padj — most statistically confident genes
        # (b) Top 50 by |log2FC| — largest effect sizes (user-requested v3.8)
        #
        # Input data: VST if available (sample_qc provides it), else
        # log2(counts+1). Then row z-score for visualization.
        # Symbol annotation: use symbol_map to replace ENSG with HGNC.
        de_results_df = de_result.get("results")
        vst_matrix    = sample_qc.get("vst_matrix")   # may be None
        symbol_map    = sample_qc.get("_symbol_map", {})  # passed via QC dict

        if de_results_df is not None and not de_results_df.empty:
            # Intersect with gene matrix (only plot genes we have data for)
            available_all = [g for g in de_results_df.index
                              if (vst_matrix is not None and g in vst_matrix.index)
                              or (counts_filt is not None and g in counts_filt.index)]

            if available_all:
                # (a) Top 50 by padj — already sig_genes-ordered by padj asc
                top_padj = [g for g in de_result.get("sig_genes", [])[:50]
                            if g in available_all]
                hm_a = _plot_heatmap(
                    genes=top_padj,
                    vst_matrix=vst_matrix,
                    counts_filt=counts_filt,
                    symbol_map=symbol_map,
                    title_prefix="Top 50 DE genes by padj",
                    title_suffix=title_suffix,
                    output_path=str(Path(output_dir) / "heatmap_padj.svg"),
                )
                if hm_a:
                    plots["heatmap"]      = hm_a   # backward compat
                    plots["heatmap_padj"] = hm_a

                # (b) Top 50 by |log2FC| (NEW in v3.8)
                lfc_sorted = de_results_df.dropna(subset=["padj","log2FoldChange"])
                lfc_sorted = lfc_sorted[lfc_sorted["padj"] < padj_thr]
                lfc_sorted = lfc_sorted.assign(
                    abs_lfc=lfc_sorted["log2FoldChange"].abs()
                ).sort_values("abs_lfc", ascending=False)
                top_lfc = [g for g in lfc_sorted.index[:50]
                            if g in available_all]
                hm_b = _plot_heatmap(
                    genes=top_lfc,
                    vst_matrix=vst_matrix,
                    counts_filt=counts_filt,
                    symbol_map=symbol_map,
                    title_prefix="Top 50 DE genes by |log2FC|",
                    title_suffix=title_suffix,
                    output_path=str(Path(output_dir) / "heatmap_lfc.svg"),
                )
                if hm_b:
                    plots["heatmap_lfc"] = hm_b

    except Exception as e:
        plots["error"] = str(e)

    # PCA and MDS were already generated in _sample_qc (v3.8: VST-based)
    if sample_qc.get("pca_plot"):
        plots["pca"] = sample_qc["pca_plot"]
    if sample_qc.get("mds_plot"):
        plots["mds"] = sample_qc["mds_plot"]

    return plots


def _plot_heatmap(genes: list, vst_matrix, counts_filt,
                    symbol_map: dict,
                    title_prefix: str, title_suffix: str,
                    output_path: str) -> str | None:
    """
    Plot a row-z-score heatmap of the given genes.

    Data source priority:
      1. VST matrix if provided (homoscedastic, lib-size corrected)
      2. log2(counts_filt+1) as fallback

    Row labels are replaced with HGNC symbols when available.
    Dark theme matching the rest of the report.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        if not genes:
            return None

        if vst_matrix is not None:
            available = [g for g in genes if g in vst_matrix.index]
            if not available:
                return None
            hm_data = vst_matrix.loc[available]
        elif counts_filt is not None:
            available = [g for g in genes if g in counts_filt.index]
            if not available:
                return None
            hm_data = np.log2(counts_filt.loc[available].astype(float) + 1)
        else:
            return None

        # Row z-score across samples
        row_mean = hm_data.mean(axis=1)
        row_std  = hm_data.std(axis=1).replace(0, 1)
        hm_z     = ((hm_data.T - row_mean) / row_std).T

        # Gene labels: symbols when available
        row_labels = []
        for g in available:
            clean = str(g).split(".")[0]
            sym   = (symbol_map or {}).get(clean, "")
            row_labels.append(sym if sym else clean)

        n_genes = len(available)
        fig, ax = plt.subplots(
            figsize=(
                max(6, len(hm_data.columns) * 0.6),
                max(5, n_genes * 0.18),
            )
        )
        fig.patch.set_facecolor(_P["fig"])

        im = ax.imshow(hm_z, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
        ax.set_xticks(range(len(hm_data.columns)))
        ax.set_xticklabels(hm_data.columns, rotation=45,
                            ha="right", fontsize=9, color=_P["text"])
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(row_labels, fontsize=7, color=_P["text"])
        ax.set_facecolor(_P["ax"])

        cbar = plt.colorbar(im, ax=ax, label="Z-score (VST or log2)",
                             fraction=0.025, pad=0.04)
        cbar.ax.tick_params(colors=_P["muted"])
        cbar.set_label("Z-score (VST or log2)", color=_P["text"])

        title = title_prefix + (f" — {title_suffix}" if title_suffix else "")
        ax.set_title(title, color=_P["title"], fontsize=11,
                     fontweight="bold", pad=10)

        plt.tight_layout()
        plt.savefig(output_path, format="svg",
                     facecolor=_P["fig"])
        plt.close(fig)
        return output_path
    except Exception as e:
        log.warning(f"Heatmap failed ({title_prefix}): {e}")
        return None


def _plot_sample_pca(coords, samples, metadata, var_exp: list,
                      output_dir: str) -> str | None:
    """Save sample PCA plot colored by condition."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=(7, 5))
        fig.patch.set_facecolor(_P["fig"])
        ax.set_facecolor(_P["ax"])

        COLORS = _P["palette"]

        # Get condition for each sample
        condition_col = None
        for col in metadata.columns:
            if col not in ("sample", "batch", "replicate"):
                condition_col = col
                break

        groups = metadata[condition_col].values if condition_col else \
                 ["unknown"] * len(samples)
        unique = sorted(set(groups))

        for i, grp in enumerate(unique):
            mask = [g == grp for g in groups]
            ax.scatter(
                coords[mask, 0], coords[mask, 1],
                label=grp, color=COLORS[i % len(COLORS)],
                s=80, alpha=0.85, edgecolors=_P["fig"], linewidths=0.5,
            )

        for j, s in enumerate(samples):
            ax.annotate(s[:10], (coords[j, 0], coords[j, 1]),
                        xytext=(4, 4), textcoords="offset points",
                        fontsize=6, color=_P["muted"])

        pct1 = round(float(var_exp[0]) * 100, 1) if len(var_exp) > 0 else 0
        pct2 = round(float(var_exp[1]) * 100, 1) if len(var_exp) > 1 else 0
        ax.set_xlabel(f"PC1 ({pct1}%)", color=_P["muted"])
        ax.set_ylabel(f"PC2 ({pct2}%)", color=_P["muted"])
        ax.set_title("Sample PCA", color=_P["title"], fontweight="bold")
        ax.tick_params(colors=_P["muted"])
        for spine in ax.spines.values():
            spine.set_edgecolor(_P["border"])
        ax.legend(fontsize=8, facecolor=_P["fig"],
                  labelcolor=_P["text"], edgecolor=_P["border"])

        pca_path = str(Path(output_dir) / "sample_pca.svg")
        plt.tight_layout()
        plt.savefig(pca_path, format="svg",
                    facecolor=_P["fig"])
        plt.close()
        return pca_path

    except Exception:
        return None


if __name__ == "__main__":
    run_script(bulk_rna_de)
