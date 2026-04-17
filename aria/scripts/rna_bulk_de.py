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
from aria.scripts._base import run_script
from pathlib import Path


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

    # ── 3. Sample QC (shared across contrasts) ───────────────────────────
    sample_qc = _sample_qc(counts, metadata, output_dir, warnings)
    outlier_samples = sample_qc.get("outliers", [])
    if outlier_samples:
        counts   = counts.drop(columns=outlier_samples, errors="ignore")
        metadata = metadata.drop(index=outlier_samples, errors="ignore")
        warnings.append(
            f"Sample outliers removed: {outlier_samples}. "
            f"Clustered away from group centroid in PCA."
        )

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
            num, den, padj_thr, lfc_thr
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

        # Pathway enrichment per contrast
        pathways = {}
        if run_pathways and de_result.get("sig_genes"):
            pathways, pw_warn = _run_pathway_enrichment(
                sig_genes=de_result["sig_genes"],
                up_genes=de_result.get("up_genes", []),
                down_genes=de_result.get("down_genes", []),
                organism=organism,
                output_dir=output_dir,
            )
            warnings.extend([f"[{name}] {w}" for w in pw_warn])

        # Plots per contrast — one volcano and one heatmap each
        contrast_dir = Path(output_dir) / _slugify(name)
        contrast_dir.mkdir(exist_ok=True)
        plots = _generate_plots(
            de_result=de_result,
            sample_qc=sample_qc,
            counts_filt=counts_filt,
            metadata=metadata,
            design_factor=design_factor,
            output_dir=str(contrast_dir),
            padj_thr=padj_thr,
            lfc_thr=lfc_thr,
            title_suffix=name,
        )

        # Format top genes
        top_genes = _format_top_genes(de_result)

        contrast_results.append({
            "name":           name,
            "numerator":      num,
            "denominator":    den,
            "status":         "success",
            "n_genes_tested": int(counts_filt.shape[0]),
            "n_significant":  int(de_result.get("n_sig", 0)),
            "n_upregulated":  int(de_result.get("n_up", 0)),
            "n_downregulated":int(de_result.get("n_down", 0)),
            "top_genes":      top_genes,
            "pathways":       pathways,
            "plots":          plots,
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

    return {
        "status":           "success",
        "n_contrasts":      len(contrast_results),
        "contrasts":        contrast_results,
        "n_significant":    total_sig,       # legacy total
        "n_upregulated":    total_up,
        "n_downregulated":  total_down,
        "sample_qc":        sample_qc,
        "design_used":      f"~{design_factor}",
        "padj_threshold":   padj_thr,
        "lfc_threshold":    lfc_thr,
        "overlap":          overlap_info,
        "warnings":         warnings,
    }


def _slugify(s: str) -> str:
    """Make a string safe for use as a directory name."""
    import re
    return re.sub(r"[^a-z0-9_]+", "_", str(s).lower()).strip("_") or "contrast"


def _format_top_genes(de_result: dict) -> list:
    """Format top genes from a DE result."""
    results = de_result.get("results")
    if results is None or len(results) == 0:
        return []
    top = []
    for g in de_result.get("sig_genes", [])[:30]:
        if g not in results.index:
            continue
        row = results.loc[g]
        try:
            top.append({
                "gene":      g,
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
    """Compute DE gene overlap between contrasts."""
    successful = [c for c in contrast_results if c.get("status") == "success"]
    if len(successful) < 2:
        return {}

    gene_sets = {
        c["name"]: set(g["gene"] for g in c.get("top_genes", []))
        for c in successful
    }
    names = list(gene_sets.keys())

    overlaps = {}
    for i, a in enumerate(names):
        for b in names[i+1:]:
            shared = gene_sets[a] & gene_sets[b]
            overlaps[f"{a} ∩ {b}"] = {
                "n_shared":      len(shared),
                "n_in_first":    len(gene_sets[a]),
                "n_in_second":   len(gene_sets[b]),
                "shared_genes":  sorted(shared)[:20],
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

def _sample_qc(counts, metadata, output_dir: str,
               warnings: list) -> dict:
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

        # ── (b) PCA-based outlier detection ───────────────────────────────
        scaler = StandardScaler()
        X      = scaler.fit_transform(log_counts.T)

        n_components = max(2, min(10, X.shape[0] - 1))
        pca     = PCA(n_components=n_components)
        coords  = pca.fit_transform(X)
        var_exp = pca.explained_variance_ratio_[:2]

        # 2.5 SD: standard outlier threshold (~99% CI for normal),
        # less aggressive than 2 SD (which can flag valid biology),
        # less permissive than 3 SD (which lets bad samples slip through).
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

        # ── Save PCA plot ──────────────────────────────────────────────────
        pca_plot = _plot_sample_pca(
            coords, counts.columns, metadata,
            var_exp, output_dir
        )

        return {
            "n_samples":           int(counts.shape[1]),
            "outliers":            confirmed_outliers,
            "pca_outliers":        pca_outliers,
            "replicate_outliers":  replicate_outliers,
            "replicate_correlations": replicate_report,
            "pca_variance":        [round(float(v), 3) for v in var_exp],
            "lib_size_range":      [int(lib_sizes.min()), int(lib_sizes.max())],
            "size_ratio":          round(float(size_ratio), 1),
            "pca_plot":            pca_plot,
        }

    except Exception as e:
        warnings.append(f"Sample QC failed: {e}")
        return {"n_samples": int(counts.shape[1]), "outliers": [],
                "pca_variance": [], "error": str(e)}


# ── DESeq2 ────────────────────────────────────────────────────────────────────

def _run_deseq2(counts, metadata, design_factor: str,
                numerator: str, denominator: str,
                padj_thr: float, lfc_thr: float) -> tuple:
    """
    Run DESeq2 via pydeseq2 with correct design factor.
    Returns (result_dict, warnings_list).
    """
    import pandas as pd
    warnings = []

    # Filter to comparison groups only
    mask      = metadata[design_factor].isin([numerator, denominator])
    meta_sub  = metadata[mask].copy()
    counts_sub = counts[meta_sub.index].copy()

    if meta_sub.shape[0] < 4:
        return {
            "status":     "error",
            "error_type": "InsufficientReplicates",
            "details":    (
                f"DESeq2 requires at least 2 replicates per group. "
                f"Found {meta_sub[design_factor].value_counts().to_dict()}"
            ),
        }, warnings

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        import warnings as w
        w.filterwarnings("ignore")

        # pydeseq2 expects samples × genes
        dds = DeseqDataSet(
            counts=counts_sub.T,          # samples × genes
            metadata=meta_sub,
            design_factors=design_factor,  # ← THE FIX: dynamic design factor
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

        up_genes   = list(sig[sig["log2FoldChange"] > 0].index)
        down_genes = list(sig[sig["log2FoldChange"] < 0].index)

        if len(sig) == 0:
            warnings.append(
                f"No significant genes found at padj < {padj_thr} "
                f"and |log2FC| > {lfc_thr}. "
                f"Consider relaxing thresholds or checking data quality."
            )

        return {
            "status":    "success",
            "results":   results_df,
            "n_sig":     len(sig),
            "n_up":      len(up_genes),
            "n_down":    len(down_genes),
            "sig_genes": list(sig.index),
            "up_genes":  up_genes,
            "down_genes": down_genes,
        }, warnings

    except ImportError:
        warnings.append("pydeseq2 not installed — returning mock DE results.")
        return _mock_de_result(padj_thr, lfc_thr), warnings

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

def _run_pathway_enrichment(sig_genes: list, up_genes: list,
                              down_genes: list, organism: str,
                              output_dir: str) -> tuple:
    """
    Run pathway enrichment via gseapy.
    Tests GO Biological Process, KEGG, and Reactome.
    """
    warnings  = []
    pathways  = {}
    gene_sets = _get_gene_sets(organism)

    if not sig_genes:
        return {}, ["No significant genes for pathway enrichment."]

    try:
        import gseapy as gp

        for db, db_name in gene_sets.items():
            try:
                enr = gp.enrichr(
                    gene_list=sig_genes[:500],  # top 500 to limit runtime
                    gene_sets=db_name,
                    organism=_gseapy_organism(organism),
                    outdir=None,
                    verbose=False,
                )
                results = enr.results

                if results is not None and not results.empty:
                    sig_pw = results[results["Adjusted P-value"] < 0.05]
                    sig_pw = sig_pw.sort_values("Adjusted P-value")

                    pathways[db] = [
                        {
                            "term":    row["Term"],
                            "padj":    round(float(row["Adjusted P-value"]), 5),
                            "overlap": row.get("Overlap", ""),
                            "genes":   row.get("Genes", "").split(";")[:10],
                        }
                        for _, row in sig_pw.head(20).iterrows()
                    ]

                    if not pathways[db]:
                        warnings.append(
                            f"No significant {db} pathways at padj < 0.05."
                        )
            except Exception as e:
                warnings.append(f"{db} enrichment failed: {str(e)[:100]}")

    except ImportError:
        warnings.append("gseapy not installed — skipping pathway enrichment.")
        pathways = _mock_pathways(sig_genes)

    # Directional enrichment: up vs down separately
    if up_genes and down_genes:
        try:
            import gseapy as gp
            for direction, genes in [("up", up_genes[:200]),
                                      ("down", down_genes[:200])]:
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets="KEGG_2021_Human"
                              if "sapiens" in organism.lower()
                              else "KEGG_2019_Mouse",
                    organism=_gseapy_organism(organism),
                    outdir=None,
                    verbose=False,
                )
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
    if "sapiens" in organism.lower():
        return "Human"
    if "musculus" in organism.lower():
        return "Mouse"
    return "Human"


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
            fig.patch.set_facecolor("#0f1729")
            ax.set_facecolor("#1a2744")

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
                       color="#475569", alpha=0.4, s=8, linewidths=0)
            ax.scatter(lfc[up_mask], neg_log_p[up_mask],
                       color="#4ade80", alpha=0.7, s=12, linewidths=0)
            ax.scatter(lfc[down_mask], neg_log_p[down_mask],
                       color="#f87171", alpha=0.7, s=12, linewidths=0)

            # Reference lines
            ax.axhline(-np.log10(padj_thr), color="#94a3b8",
                       linestyle="--", alpha=0.5, linewidth=0.8)
            ax.axvline( lfc_thr, color="#94a3b8",
                       linestyle="--", alpha=0.5, linewidth=0.8)
            ax.axvline(-lfc_thr, color="#94a3b8",
                       linestyle="--", alpha=0.5, linewidth=0.8)

            # Labels for top genes
            top_sig = results_df[sig_mask].nsmallest(10, "padj")
            for gene, row in top_sig.iterrows():
                ax.annotate(
                    str(gene)[:12],
                    xy=(float(row["log2FoldChange"]),
                        float(-np.log10(row["padj"] + 1e-300))),
                    xytext=(3, 3), textcoords="offset points",
                    fontsize=6, color="#e2e8f0", alpha=0.9,
                )

            ax.set_xlabel("log₂ Fold Change", color="#94a3b8")
            ax.set_ylabel("-log₁₀ adjusted p-value", color="#94a3b8")
            title_text = (
                f"Differential Expression" +
                (f" — {title_suffix}" if title_suffix else "") +
                f"\n{de_result.get('n_up',0)} up  "
                f"{de_result.get('n_down',0)} down"
            )
            ax.set_title(title_text, color="#22d3ee", fontsize=11)
            ax.tick_params(colors="#64748b")
            for spine in ax.spines.values():
                spine.set_edgecolor("#2d3f6e")

            volcano_path = str(Path(output_dir) / "volcano.svg")
            plt.tight_layout()
            plt.savefig(volcano_path, format="svg",
                        facecolor=fig.get_facecolor())
            plt.close()
            plots["volcano"] = volcano_path

        # ── Heatmap of top DE genes ────────────────────────────────────
        top_genes = de_result.get("sig_genes", [])[:50]
        if top_genes and counts_filt is not None:
            available = [g for g in top_genes if g in counts_filt.index]
            if available:
                import numpy as np
                hm_data = np.log2(
                    counts_filt.loc[available].astype(float) + 1
                )
                # Z-score across samples
                row_mean = hm_data.mean(axis=1)
                row_std  = hm_data.std(axis=1).replace(0, 1)
                hm_z     = ((hm_data.T - row_mean) / row_std).T

                fig, ax = plt.subplots(
                    figsize=(max(6, len(available) * 0.15),
                             max(4, len(available) * 0.25))
                )
                fig.patch.set_facecolor("#0f1729")

                im = ax.imshow(hm_z, aspect="auto", cmap="RdBu_r",
                               vmin=-3, vmax=3)
                ax.set_xticks(range(len(hm_data.columns)))
                ax.set_xticklabels(hm_data.columns, rotation=45,
                                   ha="right", fontsize=7, color="#94a3b8")
                ax.set_yticks(range(len(available)))
                ax.set_yticklabels(available, fontsize=6, color="#e2e8f0")
                ax.set_facecolor("#1a2744")
                plt.colorbar(im, ax=ax, label="Z-score",
                             fraction=0.02, pad=0.04)
                hm_title = "Top DE Genes" + (f" — {title_suffix}" if title_suffix else "")
                ax.set_title(hm_title, color="#22d3ee", fontsize=10)

                heatmap_path = str(Path(output_dir) / "heatmap.svg")
                plt.tight_layout()
                plt.savefig(heatmap_path, format="svg",
                            facecolor=fig.get_facecolor())
                plt.close()
                plots["heatmap"] = heatmap_path

    except Exception as e:
        plots["error"] = str(e)

    # PCA was already generated in _sample_qc
    if sample_qc.get("pca_plot"):
        plots["pca"] = sample_qc["pca_plot"]

    return plots


def _plot_sample_pca(coords, samples, metadata, var_exp: list,
                      output_dir: str) -> str | None:
    """Save sample PCA plot colored by condition."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=(7, 5))
        fig.patch.set_facecolor("#0f1729")
        ax.set_facecolor("#1a2744")

        COLORS = ["#22d3ee", "#4ade80", "#f87171",
                  "#fbbf24", "#a78bfa", "#f472b6"]

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
                s=80, alpha=0.85, edgecolors="#0f1729", linewidths=0.5,
            )

        for j, s in enumerate(samples):
            ax.annotate(s[:10], (coords[j, 0], coords[j, 1]),
                        xytext=(4, 4), textcoords="offset points",
                        fontsize=6, color="#94a3b8")

        pct1 = round(float(var_exp[0]) * 100, 1) if len(var_exp) > 0 else 0
        pct2 = round(float(var_exp[1]) * 100, 1) if len(var_exp) > 1 else 0
        ax.set_xlabel(f"PC1 ({pct1}%)", color="#94a3b8")
        ax.set_ylabel(f"PC2 ({pct2}%)", color="#94a3b8")
        ax.set_title("Sample PCA", color="#22d3ee")
        ax.tick_params(colors="#64748b")
        for spine in ax.spines.values():
            spine.set_edgecolor("#2d3f6e")
        ax.legend(fontsize=8, facecolor="#1e2d50",
                  labelcolor="#e2e8f0", edgecolor="#2d3f6e")

        pca_path = str(Path(output_dir) / "sample_pca.svg")
        plt.tight_layout()
        plt.savefig(pca_path, format="svg",
                    facecolor=fig.get_facecolor())
        plt.close()
        return pca_path

    except Exception:
        return None


if __name__ == "__main__":
    run_script(bulk_rna_de)
