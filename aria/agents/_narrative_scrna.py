"""
ARIA NarrativeAgent — scRNA / pseudobulk extension
---------------------------------------------------
Helpers that produce the scRNA-specific text, figures, and HTML blocks for
the NarrativeAgent. Kept in a separate module so the main agent stays
focused on cross-modal orchestration and the scRNA logic can grow without
bloating it.

The functions here are pure (no LLM, no bus). They:

  - read the same `agent_results["scrna_agent"]["findings"]["scRNA"]["findings"]`
    dict that NarrativeAgent already consumes
  - produce text summaries, methods blocks, and HTML cards
  - optionally call out to env_manager to render UMAPs in the rna stack
  - render pathway dotplots in-process via aria.scripts.rna_pathway_viz

Shape consumed (see aria/scripts/rna_narrative_adapter.py for construction):

    findings: {
      qc, integration, clustering, clustering_decision, cell_types,
      differential_expression, pathways,
      pseudobulk_de:        { groupby, condition_col, replicate_col,
                              covariates, thresholds, n_groups, per_group },
      pseudobulk_pathways:  { organism, databases, per_cluster },
      figures:              { umap_<key>: png_path, pathway_dotplots: {...},
                              per_celltype_de_bar: png_path },
    }
"""

from __future__ import annotations

import base64
import csv
import html
import logging
from pathlib import Path
from typing import Optional

log = logging.getLogger("aria.narrative.scrna")


# ── Shape normalisation ───────────────────────────────────────────────────

def unwrap_scrna_findings(agent_result: dict) -> dict:
    """
    Return the scRNA findings dict from a scrna_agent envelope, robust to
    both shapes that exist in the codebase:

        - Adapter / multimodal-wrapped:  {findings: {scRNA: {findings: {...}}}}
        - scrna_agent.run() direct:      {findings: {qc, clustering, ...}}

    Without this helper the TUI / Orchestrator path silently returns empty
    findings (the inner scRNA wrapper does not exist on the direct emit).
    """
    f = agent_result.get("findings", {}) or {}
    wrapped = (f.get("scRNA", {}) or {}).get("findings", {}) or {}
    return wrapped or f


# ── Text summaries ────────────────────────────────────────────────────────

def _fmt_int(value) -> str:
    try:
        return f"{int(value):,}"
    except Exception:
        return str(value)


def _fdr_primary_clause(pb: dict) -> str:
    """Name the adjusted-p column that defined significance, per fdr_strategy.

    Legacy results without an explicit strategy used global pooled BH.
    """
    strategy = ((pb or {}).get("multiple_testing", {}) or {}).get(
        "fdr_strategy", "global"
    )
    return "padj_local (per-cluster BH)" if strategy == "per_cluster" else "padj_global"


def _lfc_shrinkage_clause(pb: dict) -> str:
    """Disclose apeGLM LFC shrinkage when it was requested (C4)."""
    shrink = (pb or {}).get("lfc_shrinkage") or {}
    if not shrink.get("requested"):
        return ""
    return (
        "Reported log2 fold changes are apeGLM-shrunken estimates (pydeseq2 "
        "lfc_shrink); the effect-size threshold is applied to the shrunken "
        "value while p-values are unchanged, and the unshrunken MLE is kept as "
        "log2fc_raw. "
    )


def _fmt_stat(value) -> str:
    """Compact numeric display that preserves very small nonzero values."""
    if not isinstance(value, (int, float)):
        return str(value)
    if value == 0:
        return "0"
    value = float(value)
    if abs(value) < 1e-3:
        return f"{value:.2e}"
    return f"{value:.4g}"


def _group_label(groupby: str | None, n: int | None = None) -> str:
    """Human-readable label for an obs grouping column."""
    if not groupby:
        return "groups"
    if groupby in {"cell_type", "celltype", "cell_type_celltypist"}:
        singular = "cell type"
    elif groupby == "leiden":
        singular = "Leiden cluster"
    else:
        singular = f"{groupby} group"
    return singular if n == 1 else f"{singular}s"


def _label_cell_type(value) -> str:
    if isinstance(value, dict):
        return (value.get("cell_type")
                or value.get("celltypist_label")
                or value.get("label")
                or "")
    return str(value) if value else ""


def _annotation_state(findings: dict) -> dict:
    ct_block = findings.get("cell_types") or {}
    labels = [
        _label_cell_type(v)
        for v in (ct_block.get("cell_types", {}) or {}).values()
    ]
    labels = [x for x in labels if x]
    invalid = {"annotation_failed", "failed", "unknown", "nan", "none"}
    valid = [
        x for x in labels
        if x.strip().lower() not in invalid
        and not x.strip().lower().startswith("unresolved cluster")
    ]
    source = "unknown"
    for v in (ct_block.get("cell_types", {}) or {}).values():
        if isinstance(v, dict) and v.get("annotation_source"):
            source = str(v.get("annotation_source"))
            break
    label_col = ct_block.get("label_col")
    return {
        "has_valid": bool(valid),
        "labels": valid,
        "n_unique": len(set(valid)),
        "source": source,
        "label_col": label_col,
        "is_marker_fallback": source in {
            "marker_fallback", "unresolved_marker_fallback",
        },
    }


def _top_de_blocks(pb: dict, limit: int = 5) -> list[tuple[str, str, dict]]:
    rows = []
    for group, info in (pb.get("per_group", {}) or {}).items():
        for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
            if comp.get("status") == "success":
                rows.append((str(group), str(comp_key), comp))
    rows.sort(key=lambda row: row[2].get("n_significant", 0), reverse=True)
    return rows[:limit]


def _top_pathway_blocks(pwp: dict, limit: int = 3) -> list[tuple[str, dict]]:
    blocks = list((pwp.get("per_cluster", {}) or {}).items())
    blocks.sort(key=lambda kv: kv[1].get("n_significant", 0), reverse=True)
    return blocks[:limit]


def _gene_name(rec: dict) -> str:
    return str(rec.get("symbol") or rec.get("gene") or rec.get("name") or "?")


def _gene_brief(rec: dict) -> str:
    gene = _gene_name(rec)
    lfc = rec.get("log2fc", rec.get("log2FoldChange"))
    padj = rec.get("padj_global", rec.get("padj", rec.get("padj_local")))
    details = []
    if isinstance(lfc, (int, float)):
        details.append(f"log2FC={lfc:+.2f}")
    if isinstance(padj, (int, float)):
        details.append(f"FDR={_fmt_stat(padj)}")
    return f"{gene} ({', '.join(details)})" if details else gene


def _top_directional_genes(comp: dict, direction: str,
                           limit: int = 3) -> list[str]:
    records = comp.get("top_genes") or comp.get("all_sig") or []
    if direction == "up":
        rows = [r for r in records if isinstance(r, dict)
                and r.get("log2fc", 0) > 0]
        rows.sort(key=lambda r: r.get("log2fc", 0), reverse=True)
    else:
        rows = [r for r in records if isinstance(r, dict)
                and r.get("log2fc", 0) < 0]
        rows.sort(key=lambda r: r.get("log2fc", 0))
    return [_gene_brief(r) for r in rows[:limit]]


def _find_pathway_block(pwp: dict, group: str, comp_key: str) -> tuple[str, dict]:
    per_cluster = pwp.get("per_cluster", {}) or {}
    if not per_cluster:
        return "", {}
    group_s = str(group)
    comp_s = str(comp_key)
    candidates = [
        f"{group_s}::{comp_s}",
        f"{group_s}__{comp_s}",
        f"{group_s} {comp_s}",
        f"{group_s}_{comp_s}",
    ]
    for key in candidates:
        if key in per_cluster:
            return key, per_cluster[key]
    for key, block in per_cluster.items():
        key_s = str(key)
        if group_s in key_s and comp_s in key_s:
            return key_s, block
    return "", {}


def _top_pathway_terms(block: dict, limit: int = 3) -> list[str]:
    terms = []
    for db_name, rows in (block.get("results", {}) or {}).items():
        for row in rows or []:
            if not isinstance(row, dict):
                continue
            term = row.get("term") or row.get("Term")
            if not term:
                continue
            padj = _term_value(
                row, "adjusted_p", "Adjusted P-value", "adj_p", "padj",
                default=None,
            )
            label = f"{db_name}: {term}"
            if isinstance(padj, (int, float)):
                label += f" (FDR={_fmt_stat(padj)})"
            terms.append(label)
            if len(terms) >= limit:
                return terms
    return terms


def _describe_top_de_blocks(findings: dict, limit: int = 3) -> list[str]:
    """Detailed, grounded interpretation for the highest-priority DE blocks."""
    pb = findings.get("pseudobulk_de") or {}
    pwp = findings.get("pseudobulk_pathways") or {}
    lines = []
    for rank, (group, comp_key, comp) in enumerate(
        _top_de_blocks(pb, limit=limit), 1
    ):
        n_global = comp.get("n_significant_global", comp.get("n_significant", 0))
        n_local = comp.get("n_significant_local", comp.get("n_significant", 0))
        n_up = comp.get("n_up_global", comp.get("n_up", 0))
        n_down = comp.get("n_down_global", comp.get("n_down", 0))
        power = comp.get("power_estimate_at_lfc_min")
        power_txt = (
            f"; approximate power={power:.0%}"
            if isinstance(power, (int, float)) else ""
        )
        correction_txt = (
            " composition-corrected"
            if comp.get("corrected_for_composition") else
            " not composition-corrected"
        )
        fdr_txt = ""
        if n_local != n_global:
            fdr_txt = (
                f" Local FDR found {_fmt_int(n_local)} genes, but global FDR "
                f"retained {_fmt_int(n_global)}; prioritize the global-FDR "
                f"set for cross-cell-type claims."
            )
        up = _top_directional_genes(comp, "up")
        down = _top_directional_genes(comp, "down")
        gene_txt = []
        if up:
            gene_txt.append("top up: " + ", ".join(up))
        if down:
            gene_txt.append("top down: " + ", ".join(down))
        gene_clause = " " + "; ".join(gene_txt) + "." if gene_txt else ""
        _, pw_block = _find_pathway_block(pwp, group, comp_key)
        terms = _top_pathway_terms(pw_block)
        if terms:
            pathway_clause = (
                " ORA support: " + "; ".join(terms) + "."
            )
        elif pw_block:
            pathway_clause = (
                " ORA was run for this block but did not yield a compact "
                "top-term signal."
            )
        else:
            pathway_clause = (
                " No matched ORA block was available, so interpretation rests "
                "on the DE statistics alone."
            )
        caveats = []
        if comp.get("low_power_warning"):
            caveats.append("low replicate support")
        if not comp.get("corrected_for_composition"):
            skip_reason = str(comp.get("composition_skipped_reason") or "")
            if "collinear" in skip_reason:
                # C3: the covariate was deliberately dropped because it was
                # collinear with the condition; the shift is in the DA layer.
                caveats.append(
                    "composition covariate dropped (collinear with condition; "
                    "abundance shift reported separately)"
                )
            else:
                caveats.append("no composition covariate")
        caveat_txt = (
            " Caveat: " + ", ".join(caveats) + "."
            if caveats else ""
        )
        lines.append(
            f"DE block {rank}: {group} {comp_key} had "
            f"{_fmt_int(n_global)} global-FDR DE genes "
            f"({_fmt_int(n_up)} up, {_fmt_int(n_down)} down; "
            f"{_fmt_int(n_local)} local;{correction_txt}{power_txt})."
            f"{gene_clause}{pathway_clause}{fdr_txt}{caveat_txt}"
        )
    return lines


def _describe_pathway_support(findings: dict, limit: int = 3) -> list[str]:
    pwp = findings.get("pseudobulk_pathways") or {}
    if not (pwp.get("per_cluster") or {}):
        return []
    lines = []
    for block_key, block in _top_pathway_blocks(pwp, limit=limit):
        terms = _top_pathway_terms(block, limit=3)
        if not terms:
            continue
        # C2 (audit 2026-05-29): prefer this cluster's own ORA universe size
        # (genes tested in its pseudobulk); fall back to the aggregate only for
        # legacy results that lack a per-cluster background.
        bg_size = block.get("background_size") or pwp.get("background_size")
        bg = (
            f" against {_fmt_int(bg_size)} expressed genes in this cell type"
            if bg_size else ""
        )
        lines.append(
            f"Pathway support for {str(block_key).replace('::', ' ')}: "
            f"{_fmt_int(block.get('n_significant', len(terms)))} enriched "
            f"term(s){bg}; strongest terms were {', '.join(terms)}."
        )
    return lines


def _describe_abundance_de_relationship(findings: dict) -> list[str]:
    da = findings.get("differential_abundance") or {}
    pb = findings.get("pseudobulk_de") or {}
    if not da or not pb:
        return []
    shifted = set()
    for comp_info in (da.get("per_comparison") or {}).values():
        if comp_info.get("status") != "success":
            continue
        for row in comp_info.get("per_cell_type", []) or []:
            if row.get("significant"):
                shifted.add(str(row.get("name")))
    if not shifted:
        return []
    corrected = []
    uncorrected = []
    for group, info in (pb.get("per_group", {}) or {}).items():
        for comp in (info.get("per_comparison", {}) or {}).values():
            if comp.get("status") != "success":
                continue
            target = corrected if comp.get("corrected_for_composition") else uncorrected
            if str(group) in shifted:
                target.append(str(group))
    lines = []
    if corrected:
        lines.append(
            "Composition-aware interpretation: abundance shifts overlapped "
            f"with DE blocks for {', '.join(sorted(set(corrected)))}; those "
            "models included a log-proportion covariate, so within-cell-type "
            "expression effects are less confounded by changing cell-type mix."
        )
    if uncorrected:
        lines.append(
            "Composition caveat: abundance shifts were detected for "
            f"{', '.join(sorted(set(uncorrected)))}, but the matched DE block "
            "was not marked composition-corrected. Treat within-cell-type "
            "effect sizes cautiously."
        )
    return lines


def _describe_cellcomm_context(findings: dict, limit: int = 5) -> list[str]:
    ccc = findings.get("cell_communication") or {}
    if ccc.get("status") not in ("done", "success"):
        return []
    top = ccc.get("top_interactions") or []
    if not top:
        return []
    rows = []
    for ia in top[:limit]:
        if not isinstance(ia, dict):
            continue
        metric = ia.get("rank_metric") or (
            (ccc.get("method", "").split("(")[-1].rstrip(")").strip())
            if ccc.get("method") else "score"
        )
        rank = ia.get("rank")
        rank_txt = f"rank #{int(rank)}" if isinstance(rank, (int, float)) else "top-ranked"
        rows.append(
            f"{ia.get('source', '?')} -> {ia.get('target', '?')} "
            f"{ia.get('ligand', '?')}-{ia.get('receptor', '?')} "
            f"({rank_txt}, metric={metric or 'score'})"
        )
    if not rows:
        return []
    return [
        "Communication interpretation: the leading non-autocrine pairs were "
        + "; ".join(rows)
        + ". These are transcript-supported interaction candidates and need "
        "manual ligand/receptor and cell-label review before functional claims."
    ]


def _describe_trajectory_context(findings: dict) -> list[str]:
    traj = findings.get("trajectory") or {}
    if traj.get("status") not in ("done", "success"):
        return []
    paga = traj.get("paga", {}) or {}
    pt = traj.get("pseudotime", {}) or {}
    lines = []
    max_conn = paga.get("max_connectivity")
    thr = paga.get("strong_threshold", 0.05)
    if isinstance(max_conn, (int, float)):
        if max_conn < thr:
            lines.append(
                f"Trajectory depth: maximum PAGA connectivity was "
                f"{max_conn:.4f}, below the configured strong-edge threshold "
                f"({thr}); use this as neighborhood context, not as evidence "
                "for active state transitions."
            )
        else:
            lines.append(
                f"Trajectory depth: PAGA found "
                f"{_fmt_int(paga.get('n_strong', 0))} strong edge(s) at "
                f"threshold {thr}; the ranking of connected groups is useful "
                "for hypothesis generation but remains non-causal without "
                "velocity or time-course support."
            )
    if pt.get("computed") and pt.get("pseudotime_by_group"):
        ordered = sorted(
            pt.get("pseudotime_by_group", {}).items(),
            key=lambda kv: kv[1],
        )
        first = ordered[0][0]
        last = ordered[-1][0]
        lines.append(
            f"DPT context: groups span pseudotime from {first} to {last} "
            f"(root={pt.get('root_used', 'auto')}). Interpret this as an "
            "ordering on the observed manifold, not elapsed biological time."
        )
    return lines


def summarize_scrna_text(findings: dict) -> str:
    """Multi-line text summary for findings_sections['scrna']."""
    lines = []

    qc = findings.get("qc") or {}
    if qc:
        n_b = qc.get("n_cells_before")
        n_a = qc.get("n_cells_after")
        if n_b and n_a:
            n_samples = qc.get("n_samples")
            sample_txt = (f" across {n_samples} samples"
                          if n_samples not in (None, "", "?") else "")
            lines.append(
                f"Data quality and representation: after QC, "
                f"{n_a:,} of {n_b:,} cells were retained "
                f"({qc.get('pct_removed', 0)}% removed{sample_txt})."
            )
        elif n_a:
            lines.append(
                f"Data quality and representation: after QC, "
                f"{n_a:,} cells were retained."
            )

    integ = findings.get("integration") or {}
    if integ.get("status") in ("done", "success"):
        s_b = integ.get("silhouette_before")
        s_a = integ.get("silhouette_after")
        method = integ.get("method", "harmony")
        if s_b is not None and s_a is not None:
            lines.append(
                f"Batch correction ({method}) across "
                f"{integ.get('n_batches', '?')} batches: "
                f"silhouette {s_b:+.3f} → {s_a:+.3f} "
                f"(Δ={integ.get('batch_correction_delta', 0):+.3f}; "
                f"lower silhouette indicates better mixing)."
            )
    elif integ.get("status") == "skipped":
        lines.append(
            f"Batch correction was not applied: {integ.get('reason', 'skipped')}. "
            "Interpret sample/batch-colored UMAPs as diagnostic context for "
            "residual integration structure."
        )

    clu = findings.get("clustering") or {}
    if clu.get("n_clusters"):
        if clu.get("predef_clusters"):
            groupby = clu.get("groupby", "input annotation")
            lines.append(
                f"Cell-state structure: reused {clu['n_clusters']} "
                f"{_group_label(groupby, clu.get('n_clusters'))} from "
                f"input obs['{groupby}']; Leiden clustering was skipped."
            )
        else:
            lines.append(
                f"Cell-state structure: Leiden clustering identified "
                f"{clu['n_clusters']} clusters "
                f"at resolution {clu.get('resolution', '?')}."
            )

    ct = (findings.get("cell_types") or {}).get("cell_types", {}) or {}
    if ct:
        ann = _annotation_state(findings)
        unique = sorted(set(ann["labels"]))
        if unique:
            qualifier = "unresolved fallback " if ann["is_marker_fallback"] else ""
            lines.append(
                f"Cell-type annotation produced {len(unique)} "
                f"{qualifier}labels "
                f"(top: {', '.join(unique[:5])}{'…' if len(unique) > 5 else ''})."
            )
        else:
            lines.append(
                "Cell-type annotation did not produce usable biological labels; "
                "downstream sections should be interpreted at Leiden-cluster "
                "resolution."
            )

    de = findings.get("differential_expression") or {}
    de_status = de.get("status")
    if de_status and de_status != "success":
        error_type = de.get("error_type") or "Error"
        details = (de.get("details") or "").strip()
        details_clause = f": {details[:160]}" if details else ""
        lines.append(
            f"Per-cluster marker discovery (Wilcoxon) did not complete "
            f"({error_type}{details_clause}). Cluster-marker TSV is "
            "unavailable for this run; cell-type identities and "
            "between-condition results below come from independent paths "
            "(input obs labels and/or pseudobulk DE) and remain valid."
        )

    # T1.1: Cell-type abundance (differential abundance) goes BEFORE
    # pseudobulk DE because composition shifts confound within-cell-type
    # contrasts. The pseudobulk DE header below also reports whether each
    # block was composition-corrected.
    da = findings.get("differential_abundance") or {}
    if da and da.get("per_comparison") is not None:
        method = da.get("method") or "unknown"
        method_desc = (
            "donor-level centered log-ratio OLS with HC3 robust standard errors"
            if method == "donor_clr_ols_hc3" else method
        )
        if method == "donor_clr_ols_hc3" and (da.get("model") or {}).get(
            "paired_donor_fixed_effects"
        ):
            method_desc += " and donor fixed effects for the paired design"
        alpha = da.get("significance_alpha", 0.10)
        for comp_key, comp_info in (da.get("per_comparison") or {}).items():
            if comp_info.get("status") != "success":
                lines.append(
                    f"Cell-type abundance ({comp_key}): not computed "
                    f"({comp_info.get('reason', 'unknown')})."
                )
                continue
            rows = comp_info.get("per_cell_type", []) or []
            n_sig = comp_info.get("n_significant", 0)
            n_reps = comp_info.get("n_replicates", {})
            lines.append(
                f"Cell-type abundance ({comp_key}): {method_desc} on "
                f"{len(rows)} cell types, n={n_reps.get('test', '?')} vs "
                f"n={n_reps.get('ref', '?')} replicates. "
                f"{n_sig} cell type(s) shift significantly at padj < {alpha}."
            )
            if n_sig:
                top_shifts = [r for r in rows if r.get("significant")]
                top_shifts = sorted(
                    top_shifts,
                    key=lambda r: abs(r.get("log2_fold_change", 0.0)),
                    reverse=True,
                )[:5]
                desc = [
                    f"{r['name']} ({r['direction']}, "
                    f"log2FC={r['log2_fold_change']:.2f}, "
                    f"padj={r['padj']:.2g})"
                    for r in top_shifts
                ]
                if desc:
                    lines.append("Largest abundance shifts: " + "; ".join(desc) + ".")

    pb = findings.get("pseudobulk_de") or {}
    if pb:
        n_groups = pb.get("n_groups", 0)
        per_group = pb.get("per_group", {}) or {}
        n_success = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success"
        )
        n_skipped = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "skipped"
        )
        n_with_de_local = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success"
            and c.get("n_significant_local", c.get("n_significant", 0)) > 0
        )
        n_with_de_global = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success"
            and c.get("n_significant_global", c.get("n_significant", 0)) > 0
        )
        n_corrected = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success" and c.get("corrected_for_composition")
        )
        thr = pb.get("thresholds", {}) or {}
        mt = pb.get("multiple_testing", {}) or {}
        n_tests_global = mt.get("n_tests_global")
        cond_col = pb.get("condition_col") or "condition"
        cond_label = cond_col.replace("_", " ").strip()
        cond_label = cond_label[:1].upper() + cond_label[1:] if cond_label else "Condition"
        composition_clause = (
            f" {n_corrected}/{n_success} blocks were composition-corrected "
            f"(a continuous log-proportion covariate was added to the "
            f"DESeq2 design because rna_diff_abundance flagged significant "
            f"shifts)."
            if n_corrected else
            (" No block used a composition covariate "
             "(rna_diff_abundance found no significant shifts).")
            if da else ""
        )
        lines.append(
            f"{cond_label}-associated expression programs: pseudobulk DE "
            f"(DESeq2 on pseudosamples) ran across {n_groups} "
            f"{_group_label(pb.get('groupby'), n_groups)} and "
            f"{n_success} analyzable "
            f"group x comparison blocks"
            f"{f' ({n_skipped} skipped for replicate support)' if n_skipped else ''}. "
            f"{n_with_de_local} blocks yielded locally significant DE and "
            f"{n_with_de_global} blocks were significant under pooled global "
            f"BH correction"
            f"{f' across {_fmt_int(n_tests_global)} gene-block tests' if n_tests_global else ''} "
            f"at FDR < {thr.get('padj_max', 0.05)} and "
            f"|log2FC| > {thr.get('lfc_min', 0.5)}."
            f"{composition_clause}"
        )
        n_low_power = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success" and c.get("low_power_warning")
        )
        if n_low_power:
            lines.append(
                f"Caveat: {n_low_power} of {n_success} analyzable blocks ran "
                f"with n<=2 replicates on at least one side. Dispersion is "
                f"poorly estimated, effect-size estimates are noisy, and FDR "
                f"is unreliable for those blocks. Interpret with caution and "
                f"prefer n>=3 designs where possible."
            )
        top = _top_de_blocks(pb, limit=5)
        if top:
            desc = []
            for group, comp_key, comp in top:
                tag = " [low power]" if comp.get("low_power_warning") else ""
                desc.append(
                    f"{group} {comp_key}{tag}: "
                    f"{_fmt_int(comp.get('n_significant_global', comp.get('n_significant', 0)))} "
                    f"global-FDR DE genes "
                    f"({_fmt_int(comp.get('n_up_global', comp.get('n_up', 0)))} up, "
                    f"{_fmt_int(comp.get('n_down_global', comp.get('n_down', 0)))} down; "
                    f"{_fmt_int(comp.get('n_significant_local', comp.get('n_significant', 0)))} local)"
                )
            lines.append("Largest DE blocks: " + "; ".join(desc) + ".")
        lines.extend(_describe_abundance_de_relationship(findings))
        lines.extend(_describe_top_de_blocks(findings, limit=3))

    pwp = findings.get("pseudobulk_pathways") or {}
    if pwp.get("per_cluster"):
        n_blocks = len(pwp["per_cluster"])
        n_sig_blocks = sum(
            1 for b in pwp["per_cluster"].values()
            if b.get("n_significant", 0) > 0
        )
        if pwp.get("background_source") == "per_cluster_expressed_genes":
            bg_summary = (" using each cell type's own pseudobulk-expressed "
                          "genes as the ORA universe")
        elif pwp.get("background_size"):
            bg_summary = (f" using {_fmt_int(pwp.get('background_size'))} "
                          f"dataset-expressed genes as background")
        else:
            bg_summary = ""
        lines.append(
            f"Pathway over-representation (Enrichr) on top-200 DE genes per "
            f"(group × comparison): {n_sig_blocks}/{n_blocks} blocks "
            f"with significant enrichment"
            f"{bg_summary}."
        )
        examples = []
        for block_key, block in _top_pathway_blocks(pwp, limit=3):
            results = block.get("results", {}) or {}
            first_term = None
            for terms in results.values():
                if terms:
                    first = terms[0]
                    first_term = first.get("term") or first.get("Term")
                    break
            if first_term:
                examples.append(f"{block_key.replace('::', ' ')}: {first_term}")
        if examples:
            lines.append("Top enriched examples: " + "; ".join(examples) + ".")
        lines.extend(_describe_pathway_support(findings, limit=3))

    ccc = findings.get("cell_communication") or {}
    if ccc.get("status") in ("done", "success"):
        method = ccc.get("method", "?")
        n_int = ccc.get("n_interactions", 0)
        n_ct = ccc.get("n_cell_types", 0)
        n_auto = ccc.get("n_autocrine_dropped", 0)
        top_pairs = ccc.get("top_pairs", [])[:3]
        pair_str = (
            f" Top sender→receiver pairs: {'; '.join(top_pairs)}."
            if top_pairs else ""
        )
        lines.append(
            f"Cell-cell communication landscape ({method}): "
            f"{n_int} significant L-R interactions across {n_ct} cell types "
            f"({n_auto} autocrine pairs excluded).{pair_str}"
        )
        lines.extend(_describe_cellcomm_context(findings, limit=5))

    traj = findings.get("trajectory") or {}
    if traj.get("status") in ("done", "success"):
        paga = traj.get("paga", {}) or {}
        pt = traj.get("pseudotime", {}) or {}
        max_conn = paga.get("max_connectivity", 0) or 0
        n_strong = paga.get("n_strong", 0)
        thr = paga.get("strong_threshold", 0.05)
        connectivity_note = (
            f"; max edge {max_conn:.4f} — weak, consistent with "
            f"mature / non-developmental populations"
            if max_conn < thr else
            f"; {n_strong} edges above {thr} threshold "
            f"(exploratory manifold connectivity; not proof of active "
            f"differentiation without velocity or time-course data)"
        )
        traj_line = (
            f"Trajectory context: PAGA on "
            f"{paga.get('n_connections', 0)} "
            f"cluster pairs{connectivity_note}."
        )
        if pt.get("computed"):
            root = pt.get("root_used", "auto")
            pt_by = pt.get("pseudotime_by_group", {}) or {}
            if pt_by:
                ordered = sorted(pt_by.items(), key=lambda kv: kv[1])
                order_str = " → ".join(g for g, _ in ordered)
                traj_line += (
                    f" DPT pseudotime (root: {root}) orders groups "
                    f"{order_str}."
                )
        vel = traj.get("velocity", {}) or {}
        if not vel.get("computed") and vel.get("reason"):
            traj_line += (
                f" RNA velocity skipped ({vel['reason'][:60]})."
            )
        lines.append(traj_line)
        lines.extend(_describe_trajectory_context(findings))

    return "\n".join(lines) if lines else (
        "scRNA analysis completed. See findings table for details."
    )


def build_scrna_integrated_interpretation(findings: dict,
                                          intent: Optional[dict] = None) -> str:
    """
    Deterministic final interpretation from structured scRNA outputs.
    This avoids letting a generic LLM decide which completed analyses exist.
    """
    intent = intent or {}
    parts = []

    qc = findings.get("qc") or {}
    pb = findings.get("pseudobulk_de") or {}
    pwp = findings.get("pseudobulk_pathways") or {}
    ccc = findings.get("cell_communication") or {}
    traj = findings.get("trajectory") or {}
    ann = _annotation_state(findings)
    resolution_word = "cell-type" if ann["has_valid"] else "cluster"

    question = _concise_question(intent)
    if qc.get("n_cells_after"):
        parts.append(
            f"Integrated interpretation: ARIA had enough retained cells "
            f"({_fmt_int(qc.get('n_cells_after'))}) to address {question} "
            f"at {resolution_word} resolution, with the main inferential weight "
            f"coming from donor-level pseudobulk contrasts rather than "
            f"cell-level tests."
        )
    else:
        parts.append(
            f"Integrated interpretation: ARIA addressed {question} using "
            f"the structured scRNA outputs available in this run."
        )

    top_de = _top_de_blocks(pb, limit=3)
    if top_de:
        de_txt = []
        for group, comp_key, comp in top_de:
            de_txt.append(
                f"{group} {comp_key} "
                f"({_fmt_int(comp.get('n_significant', 0))} DE genes)"
            )
        parts.append(
            "The strongest between-condition transcriptional shifts were "
            + "; ".join(de_txt)
            + ". These blocks should be treated as the primary candidates "
            "for biological follow-up because they combine cell-type "
            "specificity with replicate-aware differential expression."
        )
        parts.extend(_describe_abundance_de_relationship(findings))
        parts.extend(_describe_top_de_blocks(findings, limit=3))

    if pwp.get("per_cluster"):
        n_blocks = len(pwp.get("per_cluster", {}) or {})
        n_sig = sum(
            1 for b in (pwp.get("per_cluster", {}) or {}).values()
            if b.get("n_significant", 0) > 0
        )
        parts.append(
            f"Pathway enrichment supported the DE results in {n_sig}/{n_blocks} "
            f"group x comparison blocks, giving a functional layer for "
            f"prioritising the largest pseudobulk signals."
        )
        parts.extend(_describe_pathway_support(findings, limit=3))

    if ccc.get("status") in ("done", "success"):
        top_pairs = ccc.get("top_pairs", [])[:3]
        pair_txt = f" Top ranked non-autocrine pairs were {', '.join(top_pairs)}." \
                   if top_pairs else ""
        unit = "cell types" if ann["has_valid"] else "Leiden clusters"
        parts.append(
            f"LIANA added a communication layer with "
            f"{_fmt_int(ccc.get('n_interactions', 0))} ligand-receptor "
            f"interactions across {_fmt_int(ccc.get('n_cell_types', 0))} "
            f"{unit} after excluding "
            f"{_fmt_int(ccc.get('n_autocrine_dropped', 0))} autocrine pairs."
            f"{pair_txt}"
        )
        if ann["is_marker_fallback"]:
            parts.append(
                "Because CellTypist was unavailable, communication labels are "
                "unresolved fallback labels and should be manually curated before "
                "being treated as cell identities."
            )
        parts.extend(_describe_cellcomm_context(findings, limit=5))

    if traj.get("status") in ("done", "success"):
        paga = traj.get("paga", {}) or {}
        pt = traj.get("pseudotime", {}) or {}
        traj_txt = (
            f"Trajectory analysis placed these findings in a lineage context: "
            f"PAGA evaluated {_fmt_int(paga.get('n_connections', 0))} "
            f"cluster pairs and found {_fmt_int(paga.get('n_strong', 0))} "
            f"edge(s) above the configured connectivity threshold."
        )
        if pt.get("computed"):
            traj_txt += (
                " DPT pseudotime was computed, but it should be interpreted "
                "as an ordering on the observed manifold, not as proof of "
                "active differentiation without velocity or time-course data."
            )
        if ann["is_marker_fallback"]:
            traj_txt += (
                " Group names come from unresolved fallback labels, so "
                "the trajectory section should be read as a hypothesis for "
                "manual curation."
            )
        parts.append(traj_txt)
        parts.extend(_describe_trajectory_context(findings))

    if not parts:
        return ""
    return "\n".join(parts)


def _concise_question(intent: Optional[dict]) -> str:
    raw = str((intent or {}).get("summary") or "").strip()
    if not raw:
        return "the submitted single-cell RNA-seq question"
    cut_markers = [
        "\n", " Use ", " Reuse ", " Run all ", " Do not ", " Interpret ",
        " 1.", " 2.", " 3.",
    ]
    end = len(raw)
    for marker in cut_markers:
        idx = raw.find(marker)
        if idx > 0:
            end = min(end, idx)
    concise = raw[:end].strip(" .")
    if len(concise) > 180:
        concise = concise[:177].rsplit(" ", 1)[0].rstrip(" .,") + "..."
    return concise or "the submitted single-cell RNA-seq question"


# ── Methods block ─────────────────────────────────────────────────────────

def build_scrna_methods(findings: dict) -> str:
    """Methods section text for the report (scRNA + optional pseudobulk)."""
    lines = []

    qc = findings.get("qc") or {}
    if qc:
        mt = qc.get("mt_threshold")
        mt_str = f" with mt-fraction cap at {mt}%" if mt else ""
        lines.append(
            f"Raw count matrices were processed using scanpy. "
            f"Cells were filtered by adaptive MAD-based thresholds on "
            f"total_counts, n_genes, and percent.mt{mt_str}. Doublets were "
            f"flagged with Scrublet. Counts were normalised to 10,000 "
            f"per cell and log1p-transformed."
        )

    integ = findings.get("integration") or {}
    if integ.get("status") in ("done", "success"):
        lines.append(
            f"Batch correction was performed with "
            f"{integ.get('method', 'Harmony')} on the {integ.get('rep_used', 'X_pca')} "
            f"representation across the '{integ.get('batch_col', 'batch')}' "
            f"covariate. Mixing quality was assessed by silhouette score on "
            f"the corrected embedding."
        )

    clu = findings.get("clustering") or {}
    cdec = findings.get("clustering_decision") or {}
    if clu.get("n_clusters"):
        embedding_label = findings.get("embedding_label") or "UMAP"
        # scrna_agent emits {resolution, justification, n_clusters} —
        # earlier code looked for {recommended, n_candidates} which never
        # existed, so Methods printed "resolution=? across ? candidates".
        # Fall back across both shapes for safety.
        res = (cdec.get("resolution")
               or cdec.get("recommended")
               or clu.get("resolution"))
        n_cand = cdec.get("n_candidates")
        cand_str = (f" (selected by silhouette across {n_cand} candidates)"
                    if n_cand else "")
        if clu.get("predef_clusters") or cdec.get("predef_clusters"):
            groupby = clu.get("groupby") or cdec.get("groupby") or "input annotation"
            lines.append(
                f"Dimensionality reduction used PCA (50 components) followed "
                f"by k-NN graph construction (k=15) and {embedding_label} visualisation. "
                f"ARIA reused pre-existing obs['{groupby}'] labels as "
                f"{_group_label(groupby)} and skipped Leiden clustering, "
                f"yielding {clu['n_clusters']} groups."
            )
        else:
            lines.append(
                f"Dimensionality reduction used PCA (50 components) followed by "
                f"k-NN graph construction (k=15) and {embedding_label} visualisation. "
                f"Leiden clustering at resolution={res if res is not None else '?'}"
                f"{cand_str} yielded {clu['n_clusters']} clusters."
            )

    ct = findings.get("cell_types") or {}
    if ct.get("model_used"):
        lines.append(
            f"Cell-type annotation used CellTypist with model "
            f"'{ct['model_used']}', assigning a majority label per cluster."
        )

    pb = findings.get("pseudobulk_de") or {}
    if pb:
        thr = pb.get("thresholds", {}) or {}
        mt = pb.get("multiple_testing", {}) or {}
        n_tests_clause = (
            f" (n={mt.get('n_tests_global')})"
            if mt.get("n_tests_global") else ""
        )
        powers = [
            c.get("power_estimate_at_lfc_min")
            for g in (pb.get("per_group", {}) or {}).values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success"
            and isinstance(c.get("power_estimate_at_lfc_min"), (int, float))
        ]
        power_clause = (
            f" Approximate power to detect |log2FC|>{thr.get('lfc_min', 0.5)} "
            f"ranged from {min(powers):.0%} to {max(powers):.0%} across "
            f"analyzable blocks."
            if powers else ""
        )
        cov = ", ".join(pb.get("covariates", []) or []) or "none"
        paired = bool(pb.get("paired_design"))
        paired_cov = any(
            c.get("paired_donor_covariate")
            for g in (pb.get("per_group", {}) or {}).values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success"
        )
        sample_unit = (
            f"{pb.get('replicate_col', 'replicate')} × "
            f"{pb.get('condition_col', 'condition')}"
            if paired else pb.get("replicate_col", "replicate")
        )
        design_terms = [pb.get("condition_col", "condition")]
        if paired_cov:
            design_terms.append(pb.get("replicate_col", "replicate"))
        if cov != "none":
            design_terms.append(cov)
        lines.append(
            f"Between-condition differential expression was performed by "
            f"pseudobulk aggregation: raw counts were summed per "
            f"({pb.get('groupby', 'cell_type')} × "
            f"{sample_unit}) and fitted with "
            f"pyDESeq2 (design ~ {' + '.join(design_terms)}). "
            f"Pseudosamples with < "
            f"{thr.get('min_cells_per_pseudosample', 10)} cells were dropped; "
            f"groups requiring ≥ "
            f"{thr.get('min_replicates_per_condition', 2)} replicates per "
            f"condition. Local BH correction was computed within each "
            f"cell-type × comparison block; global BH correction was computed "
            f"across all gene × block tests"
            f"{n_tests_clause}. "
            f"Significance for narrative summaries and ORA input used "
            f"{_fdr_primary_clause(pb)} &lt; {thr.get('padj_max', 0.05)} and "
            f"|log2FC| &gt; {thr.get('lfc_min', 0.5)}. "
            f"{_lfc_shrinkage_clause(pb)}"
            f"For Seurat-derived h5ads with log-normalised raw.X, counts were "
            f"recovered as expm1(x) × nCount_RNA / 10000 prior to aggregation."
            f"{power_clause}"
        )

    pwp = findings.get("pseudobulk_pathways") or {}
    if pwp.get("per_cluster"):
        dbs = list((pwp.get("databases") or {}).keys()) or [
            "GO_BP", "KEGG", "Reactome"
        ]
        if pwp.get("background_source") == "per_cluster_expressed_genes":
            bg_clause = (
                " The ORA universe was cell-type-specific: for each cell type "
                "the genes detected in its own pseudobulk, which avoids the "
                "enrichment inflation a single global background causes "
                "(per-cluster sizes are reported per block)."
            )
        elif pwp.get("background_size"):
            bg_clause = (
                f" The ORA background was {pwp.get('background_size')} genes "
                f"detected in the analyzed dataset."
            )
        else:
            bg_clause = (
                " The ORA background was Enrichr's default universe because no "
                "dataset-expressed background was available."
            )
        lines.append(
            f"Over-representation analysis (gseapy / Enrichr endpoint) was "
            f"run on the top-200 DE genes per (group × comparison) against "
            f"{', '.join(dbs)}. Significance: adjusted p &lt; 0.05."
            f"{bg_clause}"
        )

    ccc = findings.get("cell_communication") or {}
    if ccc.get("status") in ("done", "success"):
        method = ccc.get("method", "LIANA rank_aggregate")
        if "liana" in method.lower():
            lines.append(
                f"Cell-cell communication was inferred with LIANA "
                f"(rank_aggregate). The rank metric used was "
                f"'{method.split('(')[-1].rstrip(')').strip() or 'specificity_rank'}' "
                f"(lower rank = more specific or stronger interaction; "
                f"when LIANA emits only NaN magnitude_rank values the "
                f"pipeline falls back to specificity_rank to keep results "
                f"deterministic and ranked). Autocrine pairs "
                f"(source == target) are excluded a priori because they "
                f"trivially overlap and dominate the score distribution."
            )
        else:
            lines.append(
                f"Cell-cell communication was scored by mean-expression "
                f"product over a curated set of high-confidence "
                f"ligand-receptor pairs (LIANA unavailable). Autocrine "
                f"pairs are excluded."
            )

    traj = findings.get("trajectory") or {}
    if traj.get("status") in ("done", "success"):
        pt = traj.get("pseudotime", {}) or {}
        groupby = traj.get("groupby", "cell type")
        method_parts = [
            f"Trajectory inference was performed in scanpy. A k-nearest "
            f"neighbour graph (k=15) was reused from the clustering step "
            f"or recomputed on the Harmony-corrected PCA representation. "
            f"PAGA (Partition-based Graph Abstraction; Wolf et al. 2019) "
            f"was applied on the '{groupby}' grouping to estimate "
            f"cluster-level connectivity."
        ]
        if pt.get("computed"):
            root = pt.get("root_used", "auto")
            method_parts.append(
                f"Diffusion pseudotime (DPT; Haghverdi et al. 2016) was "
                f"computed on the diffusion map embedding with root cell "
                f"selected via {root}."
            )
        vel = traj.get("velocity", {}) or {}
        if not vel.get("computed"):
            method_parts.append(
                f"RNA velocity was not computed: "
                f"{vel.get('reason', 'spliced / unspliced layers absent')}. "
                f"Re-quantification with velocyto / kb-python `nac` mode "
                f"would be required to enable scVelo."
            )
        else:
            method_parts.append(
                f"RNA velocity was estimated with scVelo "
                f"({vel.get('method', 'stochastic')} model) on spliced / "
                f"unspliced layers."
            )
        lines.append(" ".join(method_parts))

    return "\n\n".join(lines)


# ── Pseudobulk DE table (HTML rows) ───────────────────────────────────────

def extract_pseudobulk_de_table(findings: dict,
                                  top_genes_per_row: int = 4) -> str:
    """
    Build an HTML <tbody> with one row per (group × comparison) summarising
    pseudobulk DE. Returns the inner rows (caller wraps in <table>).
    """
    pb = findings.get("pseudobulk_de") or {}
    per_group = pb.get("per_group", {}) or {}
    if not per_group:
        return ""

    rows = []
    # Sort groups by max n_significant across their comparisons (desc)
    def _max_sig(group_info):
        return max(
            (c.get("n_significant", 0)
             for c in (group_info.get("per_comparison", {}) or {}).values()),
            default=-1,
        )
    ordered = sorted(per_group.items(),
                     key=lambda kv: _max_sig(kv[1]),
                     reverse=True)
    for group, info in ordered:
        n_ps = info.get("n_pseudosamples", "?")
        comps = info.get("per_comparison", {}) or {}
        if not comps:
            rows.append(
                f"<tr><td>{html.escape(str(group))}</td>"
                f"<td>{n_ps}</td><td colspan='4'>"
                f"<em>{html.escape(str(info.get('reason', 'no comparison')))}</em>"
                f"</td></tr>"
            )
            continue
        for comp_key, comp in comps.items():
            if comp.get("status") == "skipped":
                rows.append(
                    f"<tr><td>{html.escape(str(group))}</td>"
                    f"<td>{n_ps}</td><td>{html.escape(str(comp_key))}</td>"
                    f"<td colspan='3' style='color:var(--muted)'>"
                    f"<em>skipped: "
                    f"{html.escape(str(comp.get('reason', '')))}</em></td></tr>"
                )
                continue
            if comp.get("status") != "success":
                continue
            n_sig = comp.get("n_significant_global", comp.get("n_significant", 0))
            n_sig_local = comp.get("n_significant_local", n_sig)
            n_up   = comp.get("n_up_global", comp.get("n_up", 0))
            n_down = comp.get("n_down_global", comp.get("n_down", 0))
            top_genes = comp.get("top_genes", []) or []
            up_tops = [g["gene"] for g in top_genes
                       if g.get("log2fc", 0) > 0][:top_genes_per_row]
            dn_tops = [g["gene"] for g in top_genes
                       if g.get("log2fc", 0) < 0][:top_genes_per_row]
            rows.append(
                f"<tr>"
                f"<td><strong>{html.escape(str(group))}</strong></td>"
                f"<td>{n_ps}</td>"
                f"<td>{html.escape(str(comp_key))}</td>"
                f"<td><strong>{n_sig}</strong> "
                f"<span style='color:var(--muted);font-size:0.85em'>"
                f"global / {n_sig_local} local</span></td>"
                f"<td style='color:var(--red)'>{n_up} ↑ "
                f"<span style='color:var(--muted);font-size:0.85em'>"
                f"{html.escape(', '.join(up_tops))}</span></td>"
                f"<td style='color:var(--blue)'>{n_down} ↓ "
                f"<span style='color:var(--muted);font-size:0.85em'>"
                f"{html.escape(', '.join(dn_tops))}</span></td>"
                f"</tr>"
            )
    return "\n".join(rows)


# ── Figure rendering ──────────────────────────────────────────────────────

def render_pathway_dotplots(findings: dict,
                              output_dir: Path,
                              top_n_blocks: int = 8,
                              top_n_terms: int = 12) -> dict:
    """
    Render ORA dotplots for the top-K (group × comparison) blocks of the
    pseudobulk_pathways stage. Uses aria.scripts.rna_pathway_viz directly
    (in-process; needs matplotlib + seaborn + pandas).

    Returns dict mapping block_key → list[png_path] (one per database).
    """
    from aria.scripts import rna_pathway_viz as pviz

    # Prefer pseudobulk_pathways (between-condition); fall back to per-cluster
    # marker pathways from the standard scRNA pipeline. Both share schema.
    pwp = (findings.get("pseudobulk_pathways")
           or findings.get("pathways") or {})
    per_cluster = pwp.get("per_cluster", {}) or {}
    if not per_cluster:
        return {}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Rank blocks by total n_significant pathways across all databases
    ranked = sorted(
        per_cluster.items(),
        key=lambda kv: kv[1].get("n_significant", 0),
        reverse=True,
    )[:top_n_blocks]

    figures: dict = {}
    for block_key, block in ranked:
        results = block.get("results", {}) or {}
        per_db = []
        for db_name, terms in results.items():
            if not terms:
                continue
            safe_block = (
                str(block_key).replace("::", "__").replace(" ", "_")
                .replace("/", "_")
            )
            out_path = output_dir / f"pathway_{safe_block}__{db_name}.png"
            png = pviz.make_ora_dotplot(
                pathways_list=terms,
                db_name=db_name,
                contrast_name=block_key.replace("::", " — "),
                output_path=str(out_path),
                top_n=top_n_terms,
            )
            if png:
                per_db.append(png)
        if per_db:
            figures[block_key] = per_db
    return figures


def render_cellcomm_heatmap(findings: dict,
                              output_path: Path) -> Optional[str]:
    """
    Heatmap of n_interactions per (source × target). Top_interactions only.
    """
    ccc = findings.get("cell_communication") or {}
    top = ccc.get("top_interactions") or []
    if not top:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df = pd.DataFrame(top)
    if df.empty or "source" not in df.columns or "target" not in df.columns:
        return None

    counts = (df.groupby(["source", "target"]).size()
              .reset_index(name="n"))
    cell_types = sorted(set(counts["source"]).union(counts["target"]))
    mat = pd.DataFrame(0, index=cell_types, columns=cell_types, dtype=int)
    for _, r in counts.iterrows():
        mat.loc[r["source"], r["target"]] = int(r["n"])

    n = len(cell_types)
    fig, ax = plt.subplots(figsize=(max(5, 0.55 * n + 2),
                                     max(4, 0.55 * n + 1.5)),
                           dpi=160)
    im = ax.imshow(mat.values, cmap="magma_r", aspect="auto",
                   vmin=0, vmax=max(1, mat.values.max()))
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(cell_types, rotation=55, ha="right", fontsize=8)
    ax.set_yticklabels(cell_types, fontsize=8)
    ax.set_xlabel("Receiver", fontsize=9)
    ax.set_ylabel("Sender", fontsize=9)
    ax.set_title(
        f"Cell-cell communication — interactions among top {len(top)} "
        f"(autocrine excluded)",
        fontsize=10, fontweight="bold",
    )
    # Cell-level labels
    for i in range(n):
        for j in range(n):
            v = mat.values[i, j]
            if v > 0:
                ax.text(j, i, str(v),
                        ha="center", va="center", fontsize=7,
                        color="white" if v >= mat.values.max() * 0.5
                              else "#1e293b")
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02,
                 label="n interactions").ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


def render_cellcomm_top_pairs_bar(findings: dict,
                                    output_path: Path,
                                    top_n: int = 15) -> Optional[str]:
    """Horizontal barplot of top-N L-R interactions by rank/score."""
    ccc = findings.get("cell_communication") or {}
    top = ccc.get("top_interactions") or []
    if not top:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Ranks: lower = better for spec/mag. Prefer explicit rank order when
    # present so tied/underflowed LIANA scores do not render as all-zero bars.
    rows = []
    for ia in top[:top_n]:
        label = (f"{ia.get('source', '?')[:18]} → "
                 f"{ia.get('target', '?')[:18]}  "
                 f"({ia.get('ligand', '?')}-{ia.get('receptor', '?')})")
        rows.append((label, float(ia.get("rank", ia.get("score", 0)))))
    if not rows:
        return None
    labels = [r[0] for r in rows]
    scores = np.array([r[1] for r in rows], dtype=float)

    is_rank = (ccc.get("method") or "").startswith("liana")
    if is_rank:
        # Invert: higher bar = lower rank = better.
        max_score = scores.max() if scores.max() > 0 else 1.0
        plot_vals = max_score - scores + (max_score * 0.05)
        xlabel = f"strength (inverted rank; metric: {ccc.get('method')})"
    else:
        plot_vals = scores
        xlabel = "score"

    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(rows) + 1.2)),
                           dpi=160)
    y = np.arange(len(rows))
    ax.barh(y, plot_vals, color="#0d9488", edgecolor="#0f172a",
            linewidth=0.4, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()  # first = top
    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_title("Top ligand-receptor interactions",
                 fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


def extract_cellcomm_table(findings: dict, top_n: int = 20) -> str:
    """HTML <tbody> rows for top L-R interactions."""
    ccc = findings.get("cell_communication") or {}
    top = (ccc.get("top_interactions") or [])[:top_n]
    if not top:
        return ""
    rows = []
    for ia in top:
        src = html.escape(str(ia.get("source", "")))
        tgt = html.escape(str(ia.get("target", "")))
        lig = html.escape(str(ia.get("ligand", "")))
        rec = html.escape(str(ia.get("receptor", "")))
        score = ia.get("score", "?")
        rank = ia.get("rank")
        metric = ia.get("rank_metric") or (
            (ccc.get("method", "").split("(")[-1].rstrip(")").strip())
            if ccc.get("method") else ""
        )
        pval = ia.get("cellphone_pval")
        pval_str = (f"<code>{_fmt_stat(pval)}</code>"
                    if isinstance(pval, (int, float)) and pval > 0 else "—")
        rank_str = f"#{int(rank)}" if isinstance(rank, (int, float)) else "—"
        score_str = _fmt_stat(score)
        metric_str = html.escape(str(metric or "score"))
        rows.append(
            f"<tr><td>{src}</td><td>{tgt}</td>"
            f"<td><strong>{lig}</strong></td><td>{rec}</td>"
            f"<td>{rank_str}</td><td><code>{score_str}</code><br>"
            f"<span style='color:var(--muted);font-size:0.82em'>"
            f"{metric_str}</span></td><td>{pval_str}</td></tr>"
        )
    return "\n".join(rows)


def extract_trajectory_tables(findings: dict) -> dict:
    """
    Return {paga_rows, pseudotime_rows} HTML strings for the trajectory
    section. Each is a <tbody> inner snippet (caller wraps in <table>).
    """
    traj = findings.get("trajectory") or {}
    paga = traj.get("paga", {}) or {}
    pt = traj.get("pseudotime", {}) or {}

    paga_rows = ""
    top_conn = paga.get("top_connections", {}) or {}
    if top_conn:
        max_c = paga.get("max_connectivity") or max(
            (v for v in top_conn.values() if isinstance(v, (int, float))),
            default=0,
        )
        thr = paga.get("strong_threshold", 0.05)
        rows = []
        for edge, val in top_conn.items():
            strong = isinstance(val, (int, float)) and val >= thr
            badge = (
                '<span style="background:#dcfce7;color:var(--green);'
                'padding:2px 6px;border-radius:3px;font-size:0.75em;'
                'font-weight:600">strong</span>'
                if strong else
                '<span style="color:var(--muted);font-size:0.85em">weak</span>'
            )
            rel = (val / max_c) if max_c else 0
            bar_w = max(2, int(rel * 100))
            rows.append(
                f"<tr><td>{html.escape(str(edge))}</td>"
                f"<td><code>{val}</code></td>"
                f"<td>{badge}</td>"
                f"<td><div style='background:#e2e8f0;width:120px;"
                f"height:8px;border-radius:3px;overflow:hidden'>"
                f"<div style='background:var(--teal);width:{bar_w}px;"
                f"height:100%;'></div></div></td></tr>"
            )
        paga_rows = "\n".join(rows)

    pseudotime_rows = ""
    pt_by = pt.get("pseudotime_by_group", {}) or {}
    if pt_by:
        ordered = sorted(pt_by.items(), key=lambda kv: kv[1])
        max_pt = max(pt_by.values()) if pt_by else 1
        rows = []
        for rank, (group, val) in enumerate(ordered, 1):
            rel = (val / max_pt) if max_pt else 0
            bar_w = max(2, int(rel * 140))
            rows.append(
                f"<tr><td>{rank}</td>"
                f"<td><strong>{html.escape(str(group))}</strong></td>"
                f"<td><code>{val:.4f}</code></td>"
                f"<td><div style='background:#e2e8f0;width:160px;"
                f"height:8px;border-radius:3px;overflow:hidden'>"
                f"<div style='background:var(--blue);width:{bar_w}px;"
                f"height:100%'></div></div></td></tr>"
            )
        pseudotime_rows = "\n".join(rows)

    return {"paga_rows": paga_rows,
            "pseudotime_rows": pseudotime_rows}


# ── Supplementary table export ────────────────────────────────────────────

def _write_tsv(path: Path, rows: list[dict]) -> Optional[str]:
    """Write rows as TSV. Returns path on success, None if no rows."""
    if not rows:
        return None
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def _term_value(term: dict, *keys, default=""):
    for key in keys:
        if key in term:
            return term.get(key)
    return default


def export_supplementary_tables(findings: dict, output_dir: Path) -> dict:
    """
    Materialize scRNA result objects into report/tables/*.tsv.

    The analytical scripts often return rich in-memory structures for the
    NarrativeAgent but the report staging layer only copied bulk RNA tables.
    This exporter keeps the report directory self-contained for scRNA runs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables: dict[str, str] = {}

    qc = findings.get("qc") or {}
    qc_rows = []
    for sample in qc.get("per_sample", []) or []:
        if isinstance(sample, dict):
            qc_rows.append(sample)
    p = _write_tsv(output_dir / "scrna_qc_per_sample.tsv", qc_rows)
    if p:
        tables["qc_per_sample"] = p

    ct = (findings.get("cell_types") or {}).get("cell_types", {}) or {}
    ct_rows = []
    for cluster, value in ct.items():
        row = {"cluster": cluster, "label": _label_cell_type(value)}
        if isinstance(value, dict):
            for key, val in value.items():
                if isinstance(val, (str, int, float, bool)) or val is None:
                    row[key] = val
                elif isinstance(val, list):
                    row[key] = ", ".join(map(str, val))
        ct_rows.append(row)
    p = _write_tsv(output_dir / "scrna_cell_types.tsv", ct_rows)
    if p:
        tables["cell_types"] = p

    # Standard per-cluster marker DE.
    de = findings.get("differential_expression") or {}
    marker_rows = []
    for cluster, genes in (de.get("de_genes_by_cluster", {}) or {}).items():
        for gene in genes or []:
            if isinstance(gene, dict):
                row = {"cluster": cluster}
                row.update(gene)
            else:
                row = {"cluster": cluster, "gene": gene}
            marker_rows.append(row)
    p = _write_tsv(output_dir / "scrna_cluster_markers.tsv", marker_rows)
    if p:
        tables["cluster_markers"] = p

    pb = findings.get("pseudobulk_de") or {}
    pb_summary_rows = []
    pb_gene_rows = []
    for group, info in (pb.get("per_group", {}) or {}).items():
        n_ps = info.get("n_pseudosamples")
        for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
            pb_summary_rows.append({
                "group": group,
                "comparison": comp_key,
                "status": comp.get("status"),
                "n_pseudosamples": n_ps,
                "n_significant": comp.get("n_significant", 0),
                "n_significant_local": comp.get(
                    "n_significant_local", comp.get("n_significant", 0)
                ),
                "n_significant_global": comp.get(
                    "n_significant_global", comp.get("n_significant", 0)
                ),
                "n_up": comp.get("n_up", 0),
                "n_up_local": comp.get("n_up_local", comp.get("n_up", 0)),
                "n_up_global": comp.get("n_up_global", comp.get("n_up", 0)),
                "n_down": comp.get("n_down", 0),
                "n_down_local": comp.get(
                    "n_down_local", comp.get("n_down", 0)
                ),
                "n_down_global": comp.get(
                    "n_down_global", comp.get("n_down", 0)
                ),
                "reason": comp.get("reason", ""),
            })
            records = comp.get("all_sig") or comp.get("top_genes") or []
            for rec in records:
                if not isinstance(rec, dict):
                    continue
                row = {
                    "group": group,
                    "comparison": comp_key,
                    "gene": rec.get("gene"),
                    "log2fc": rec.get("log2fc"),
                    "padj": rec.get("padj"),
                    "padj_local": rec.get("padj_local"),
                    "padj_global": rec.get("padj_global"),
                    "pvalue": rec.get("pvalue"),
                }
                for key, val in rec.items():
                    if key not in row:
                        row[key] = val
                pb_gene_rows.append(row)
    p = _write_tsv(output_dir / "scrna_pseudobulk_de_summary.tsv",
                   pb_summary_rows)
    if p:
        tables["pseudobulk_de_summary"] = p
    p = _write_tsv(output_dir / "scrna_pseudobulk_de_genes.tsv",
                   pb_gene_rows)
    if p:
        tables["pseudobulk_de_genes"] = p

    def _pathway_rows(container: dict, mode: str) -> list[dict]:
        rows = []
        for block_key, block in (container.get("per_cluster", {}) or {}).items():
            results = block.get("results", {}) or {}
            for db_name, terms in results.items():
                for term in terms or []:
                    if not isinstance(term, dict):
                        continue
                    rows.append({
                        "mode": mode,
                        "block": block_key,
                        "database": db_name,
                        "term": _term_value(term, "term", "Term"),
                        "adjusted_p": _term_value(
                            term, "adjusted_p", "Adjusted P-value",
                            "adj_p", "padj",
                        ),
                        "p_value": _term_value(term, "p_value", "P-value"),
                        "overlap": _term_value(term, "overlap", "Overlap"),
                        "odds_ratio": _term_value(
                            term, "odds_ratio", "Odds Ratio"
                        ),
                        "combined_score": _term_value(
                            term, "combined_score", "Combined Score"
                        ),
                        "genes": _term_value(
                            term, "genes", "Genes", "lead_genes"
                        ),
                    })
        return rows

    pathway_rows = []
    pathway_rows.extend(_pathway_rows(findings.get("pathways") or {},
                                      "cluster_markers"))
    pathway_rows.extend(_pathway_rows(findings.get("pseudobulk_pathways") or {},
                                      "pseudobulk_de"))
    p = _write_tsv(output_dir / "scrna_pathway_enrichment.tsv",
                   pathway_rows)
    if p:
        tables["pathway_enrichment"] = p

    ccc = findings.get("cell_communication") or {}
    cc_rows = []
    for rec in ccc.get("top_interactions", []) or []:
        if isinstance(rec, dict):
            cc_rows.append(rec)
    p = _write_tsv(output_dir / "scrna_cellcomm_interactions.tsv", cc_rows)
    if p:
        tables["cellcomm_interactions"] = p

    traj = findings.get("trajectory") or {}
    paga = traj.get("paga", {}) or {}
    paga_rows = []
    for edge, val in (paga.get("top_connections", {}) or {}).items():
        if "->" in str(edge):
            source, target = str(edge).split("->", 1)
        elif "→" in str(edge):
            source, target = str(edge).split("→", 1)
        else:
            source, target = "", ""
        paga_rows.append({
            "edge": edge,
            "source": source.strip(),
            "target": target.strip(),
            "connectivity": val,
            "strong_threshold": paga.get("strong_threshold"),
            "is_strong": (
                isinstance(val, (int, float))
                and val >= (paga.get("strong_threshold", 0.05) or 0.05)
            ),
        })
    p = _write_tsv(output_dir / "scrna_paga_connections.tsv", paga_rows)
    if p:
        tables["paga_connections"] = p

    pt = traj.get("pseudotime", {}) or {}
    pt_rows = [
        {"group": group, "mean_dpt": val}
        for group, val in (pt.get("pseudotime_by_group", {}) or {}).items()
    ]
    pt_rows.sort(key=lambda r: r["mean_dpt"])
    for i, row in enumerate(pt_rows, 1):
        row["rank"] = i
    p = _write_tsv(output_dir / "scrna_pseudotime_by_group.tsv", pt_rows)
    if p:
        tables["pseudotime_by_group"] = p

    if tables:
        findings["tables"] = tables
    return tables


def render_per_celltype_de_bar(findings: dict, output_path: Path) -> Optional[str]:
    """
    Stacked bar chart of n_up / n_down DE genes per (group × comparison)
    from the pseudobulk stage. Returns PNG path on success, None if no data.
    """
    pb = findings.get("pseudobulk_de") or {}
    per_group = pb.get("per_group", {}) or {}
    rows = []
    for group, info in per_group.items():
        for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
            if comp.get("status") != "success":
                continue
            rows.append((
                f"{group}\n({comp_key})",
                int(comp.get("n_up", 0)),
                int(comp.get("n_down", 0)),
            ))
    if not rows:
        return None

    rows.sort(key=lambda r: r[1] + r[2], reverse=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [r[0] for r in rows]
    n_up   = np.array([r[1] for r in rows])
    n_down = np.array([r[2] for r in rows])

    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(rows) + 2), 4.2),
                           dpi=160)
    x = np.arange(len(rows))
    ax.bar(x, n_up,   color="#991b1b", label="up",   width=0.7)
    ax.bar(x, -n_down, color="#1d4ed8", label="down", width=0.7)
    ax.axhline(0, color="#1e293b", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("DE genes  (up ↑ / down ↓)", fontsize=9)
    ax.set_title("Pseudobulk DE — per cell type",
                 fontsize=10, fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


# ── HTML embedding ────────────────────────────────────────────────────────

def _embed_png(path: str) -> str:
    """Inline a PNG as a base64 data URI (returns '' on failure)."""
    try:
        p = Path(path)
        if not p.exists():
            return ""
        data = p.read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        return f"data:image/png;base64,{b64}"
    except Exception as e:
        log.warning(f"_embed_png failed for {path}: {e}")
        return ""


def build_scrna_html_section(findings: dict,
                              max_pathway_blocks: int = 8) -> str:
    """
    Build the inner-HTML for the scRNA findings card: UMAPs, per-cell-type
    DE bar, DE table, and pathway dotplot mosaic. Returns concatenated
    HTML snippet (no outer <div class="card"> — caller wraps).
    """
    parts: list[str] = []
    figs = findings.get("figures") or {}
    ann = _annotation_state(findings)
    label_unit = "cell type" if ann["has_valid"] else "cluster"
    if ann["is_marker_fallback"]:
        parts.append(
            '<div class="warning">'
            'Cell labels in this report are unresolved fallback labels because '
            'CellTypist did not complete. Treat UMAP, '
            'trajectory, and communication labels as curation targets, not '
            'final cell identities.'
            '</div>'
        )
    elif not ann["has_valid"]:
        parts.append(
            '<div class="warning">'
            'Cell-type annotation did not produce usable biological labels. '
            'Embedding, trajectory, and communication results are reported at '
            'Leiden-cluster resolution.'
            '</div>'
        )

    # 1. UMAP figures ─────────────────────────────────────────────────────
    # Exclude trajectory-specific UMAPs (e.g. dpt_pseudotime) — those are
    # rendered inside the Trajectory section to keep the narrative
    # contiguous.
    umaps = {
        k: v for k, v in figs.items()
        if k.startswith("umap_") and k != "umap_dpt_pseudotime"
    }
    if umaps:
        embedding_label = findings.get("embedding_label") or "UMAP"
        parts.append('<h4 style="margin-top:1rem">Embedding</h4>')
        parts.append('<div style="display:flex;flex-wrap:wrap;gap:1rem">')
        for key, path in sorted(umaps.items()):
            uri = _embed_png(path)
            if not uri:
                continue
            pretty = key.replace("umap_", "")
            pretty = {
                "cell_type_marker": "marker-based cell label",
                "cell_type_celltypist": "CellTypist cell label",
                "leiden": "Leiden cluster",
                "batch": "batch",
                "sample_id": "sample",
            }.get(pretty, pretty)
            caption = html.escape(f"{embedding_label} — {pretty}")
            parts.append(
                f'<figure style="flex:1 1 320px;min-width:300px;max-width:480px">'
                f'<img src="{uri}" alt="{caption}">'
                f'<figcaption>{caption}</figcaption>'
                f'</figure>'
            )
        parts.append('</div>')

    # 2. Per-cell-type DE summary bar ─────────────────────────────────────
    de_bar = figs.get("per_celltype_de_bar")
    if de_bar:
        uri = _embed_png(de_bar)
        if uri:
            parts.append(f'<h4>Pseudobulk DE — counts per {label_unit}</h4>')
            parts.append(
                f'<figure><img src="{uri}" '
                f'alt="Per cell-type DE bar"></figure>'
            )

    # 3. Pseudobulk DE table ──────────────────────────────────────────────
    table_rows = extract_pseudobulk_de_table(findings)
    if table_rows:
        parts.append(f'<h4>DE summary by ({label_unit} × comparison)</h4>')
        parts.append(
            '<table style="width:100%;font-size:0.85em">'
            '<thead><tr>'
                f'<th>{html.escape(label_unit.title())}</th>'
            '<th>n<sub>pseudo</sub></th>'
            '<th>Comparison</th>'
            '<th>Sig.</th>'
            '<th>Up (top genes)</th>'
            '<th>Down (top genes)</th>'
            '</tr></thead>'
            f'<tbody>{table_rows}</tbody>'
            '</table>'
        )

    # 4. Trajectory section (PAGA + DPT) ──────────────────────────────────
    traj = findings.get("trajectory") or {}
    if traj.get("status") in ("done", "success"):
        parts.append('<h4 style="margin-top:1.4rem">'
                     f'Trajectory — PAGA + DPT by {label_unit}</h4>')
        parts.append(
            '<p style="color:var(--muted);font-size:0.88em">'
            'This section reports graph connectivity and DPT ordering for '
            f'{html.escape(label_unit)} groups. It is an exploratory manifold '
            'summary, not causal evidence of differentiation by itself.'
            '</p>'
        )

        # PAGA + DPT-coloured UMAP figures, side by side
        traj_figs = []
        for fkey in ("paga_graph", "paga_log10_graph",
                     "umap_dpt_pseudotime"):
            p = figs.get(fkey)
            if p:
                uri = _embed_png(p)
                if uri:
                    cap = html.escape({
                        "paga_graph":         "PAGA — group connectivity",
                        "paga_log10_graph":   "PAGA — log-scaled edges",
                        "umap_dpt_pseudotime": "UMAP — DPT pseudotime",
                    }[fkey])
                    traj_figs.append(
                        f'<figure style="flex:1 1 280px;min-width:260px;'
                        f'max-width:430px">'
                        f'<img src="{uri}" alt="{cap}">'
                        f'<figcaption>{cap}</figcaption>'
                        f'</figure>'
                    )
        if traj_figs:
            parts.append('<div style="display:flex;flex-wrap:wrap;'
                         'gap:1rem">')
            parts.extend(traj_figs)
            parts.append('</div>')

        # Tables: PAGA top connections + DPT pseudotime by group
        tables = extract_trajectory_tables(findings)
        if tables["paga_rows"]:
            paga_meta = traj.get("paga", {}) or {}
            max_c = paga_meta.get("max_connectivity", 0)
            n_str = paga_meta.get("n_strong", 0)
            thr_  = paga_meta.get("strong_threshold", 0.05)
            note = (
                f'<p style="color:var(--muted);font-size:0.82em;'
                f'margin-top:0.4rem;font-style:italic">'
                f'Max connectivity = {max_c:.4f}. '
                f'{n_str} edge(s) above the {thr_} threshold. '
                f'In mature / non-developmental populations, absolute '
                f'connectivities are typically &lt; 0.01 — interpret '
                f'rankings rather than absolute magnitudes.</p>'
            )
            parts.append(
                '<h4 style="margin-top:1.2rem">PAGA — top connections</h4>'
                '<table style="width:100%;font-size:0.88em">'
                '<thead><tr><th>Edge</th><th>Connectivity</th>'
                '<th>Strength</th><th>Visual</th></tr></thead>'
                f'<tbody>{tables["paga_rows"]}</tbody>'
                '</table>'
                + note
            )
        if tables["pseudotime_rows"]:
            pt = traj.get("pseudotime", {}) or {}
            root_str = html.escape(str(pt.get("root_used", "auto")))
            parts.append(
                f'<h4 style="margin-top:1.2rem">DPT pseudotime by '
                f'{html.escape(label_unit)} group '
                f'(root: {root_str})</h4>'
                '<table style="width:100%;font-size:0.88em">'
                '<thead><tr><th>Rank</th><th>Group</th>'
                '<th>Mean DPT</th><th>Visual</th></tr></thead>'
                f'<tbody>{tables["pseudotime_rows"]}</tbody>'
                '</table>'
            )

    # 4b. Cell-cell communication section ─────────────────────────────────
    ccc = findings.get("cell_communication") or {}
    if ccc.get("status") in ("done", "success"):
        parts.append('<h4 style="margin-top:1.4rem">'
                     f'Cell-cell communication by {label_unit}</h4>')
        parts.append(
            '<p style="color:var(--muted);font-size:0.88em">'
            'Ligand-receptor scores are summarized between observed '
            f'{html.escape(label_unit)} groups. These results require '
            'manual review of sender and receiver labels before biological '
            'interpretation.'
            '</p>'
        )

        ccc_figs = []
        for fkey, caption in (
            ("cellcomm_heatmap",   "Sender → receiver interaction count"),
            ("cellcomm_top_pairs", "Top ligand-receptor interactions"),
        ):
            p = figs.get(fkey)
            if p:
                uri = _embed_png(p)
                if uri:
                    ccc_figs.append(
                        f'<figure style="flex:1 1 320px;min-width:300px;'
                        f'max-width:520px"><img src="{uri}" '
                        f'alt="{html.escape(caption)}">'
                        f'<figcaption>{html.escape(caption)}</figcaption>'
                        f'</figure>'
                    )
        if ccc_figs:
            parts.append('<div style="display:flex;flex-wrap:wrap;'
                         'gap:1rem">')
            parts.extend(ccc_figs)
            parts.append('</div>')

        cc_rows = extract_cellcomm_table(findings)
        if cc_rows:
            method = html.escape(str(ccc.get("method", "?")))
            n_ct = ccc.get("n_cell_types", "?")
            n_int = ccc.get("n_interactions", "?")
            n_auto = ccc.get("n_autocrine_dropped", 0)
            parts.append(
                f'<h4 style="margin-top:1rem">'
                f'Top L-R interactions  '
                f'<span style="color:var(--muted);font-weight:400;'
                f'font-size:0.85em">'
                f'({method} · {n_int} interactions across {n_ct} '
                f'{html.escape(label_unit)} groups · {n_auto} autocrine '
                f'pairs excluded)</span></h4>'
                '<table style="width:100%;font-size:0.85em">'
                '<thead><tr><th>Sender</th><th>Receiver</th>'
                '<th>Ligand</th><th>Receptor</th>'
                '<th>Rank</th><th>Metric value</th><th>CellPhone p</th>'
                '</tr></thead>'
                f'<tbody>{cc_rows}</tbody></table>'
            )

    # 5. Pathway dotplots ─────────────────────────────────────────────────
    pw_figs = figs.get("pathway_dotplots") or {}
    if pw_figs:
        parts.append(
            '<h4 style="margin-top:1.4rem">'
            f'Pathway enrichment — top {label_unit} groups</h4>'
        )
        # Render up to N blocks, two-column grid
        n = 0
        parts.append(
            '<div style="display:grid;'
            'grid-template-columns:repeat(auto-fit, minmax(340px, 1fr));'
            'gap:1rem">'
        )
        for block_key, png_list in pw_figs.items():
            if n >= max_pathway_blocks:
                break
            for png in png_list:
                uri = _embed_png(png)
                if not uri:
                    continue
                # Database label = filename suffix between __ and .png
                db_label = (
                    Path(png).stem.rsplit("__", 1)[-1]
                    if "__" in Path(png).stem else "ORA"
                )
                caption = html.escape(
                    f"{block_key.replace('::', ' — ')}  ·  {db_label}"
                )
                parts.append(
                    f'<figure><img src="{uri}" alt="{caption}">'
                    f'<figcaption>{caption}</figcaption></figure>'
                )
            n += 1
        parts.append('</div>')

    # 6. Supplementary table links ────────────────────────────────────────
    table_links = findings.get("tables") or {}
    if table_links:
        labels = {
            "qc_per_sample": "QC per sample",
            "cell_types": "Cell types",
            "cluster_markers": "Cluster markers",
            "pseudobulk_de_summary": "Pseudobulk DE summary",
            "pseudobulk_de_genes": "Pseudobulk DE genes",
            "pathway_enrichment": "Pathway enrichment",
            "cellcomm_interactions": "Cell-cell communication",
            "paga_connections": "PAGA connections",
            "pseudotime_by_group": "Pseudotime by group",
        }
        links = []
        for key, path in table_links.items():
            name = labels.get(key, key.replace("_", " ").title())
            rel = f"tables/{Path(path).name}"
            links.append(
                f'<a href="{html.escape(rel)}" style="color:var(--blue);'
                f'text-decoration:underline">{html.escape(name)}</a>'
            )
        if links:
            parts.append(
                '<h4 style="margin-top:1.4rem">Supplementary tables</h4>'
                '<p style="font-size:0.85rem;color:var(--muted)">'
                + " &middot; ".join(links)
                + '</p>'
            )

    return "\n".join(parts)


# ── Top-level orchestration helper ────────────────────────────────────────

def generate_figures(findings: dict,
                      h5ad_path: Optional[str],
                      output_dir: Path,
                      env_manager=None,
                      umap_color_keys: Optional[list] = None) -> dict:
    """
    Generate all scRNA figures and write their paths back into
    findings['figures']. Mutates `findings` and returns it.

    Args:
        findings:        scRNA findings dict (mutated in place).
        h5ad_path:       path to AnnData with UMAP. Required for UMAP figs.
        output_dir:      where to write PNGs.
        env_manager:     ARIA EnvironmentManager (so UMAP runs in rna stack).
                         If None, UMAP rendering is skipped.
        umap_color_keys: obs columns to color the UMAP by. If None, picks
                         sensible defaults from the design.

    Returns:
        findings dict with `figures` populated:
            {umap_<key>: png_path, per_celltype_de_bar: png_path,
             pathway_dotplots: {block_key: [png_paths]}}
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figs = findings.setdefault("figures", {})

    # 1. UMAP figures via the rna stack ─────────────────────────────────
    if h5ad_path and env_manager is not None:
        if umap_color_keys is None:
            keys: list = []
            ann = _annotation_state(findings)
            if ann.get("label_col"):
                keys.append(ann["label_col"])
            pb = findings.get("pseudobulk_de") or {}
            if pb.get("groupby") and pb.get("groupby") not in keys:
                keys.append(pb["groupby"])
            if pb.get("condition_col") and pb.get("condition_col") not in keys:
                keys.append(pb["condition_col"])
            if not keys:
                # Standard mode: pick a sensible set in priority order. The
                # rna_figure_umap script silently skips missing columns, so
                # listing redundant fallbacks is safe.
                ct = findings.get("cell_types") or {}
                for candidate in (
                    "cell_type_celltypist",
                    ct.get("label_col"),
                    "leiden",
                    findings.get("integration", {}).get("batch_col"),
                    "batch",
                    "sample_id",
                ):
                    if candidate and candidate not in keys:
                        keys.append(candidate)
            else:
                for candidate in (
                    "leiden",
                    findings.get("integration", {}).get("batch_col"),
                    "batch",
                    "sample_id",
                ):
                    if candidate and candidate not in keys:
                        keys.append(candidate)
            umap_color_keys = keys
        if umap_color_keys:
            try:
                res = env_manager.run_in_stack(
                    stack="rna",
                    script_path="aria/scripts/rna_figure_umap.py",
                    params={
                        "h5ad_path":  str(h5ad_path),
                        "color_by":   umap_color_keys,
                        "output_dir": str(output_dir),
                    },
                )
                if res.get("status") == "success":
                    if res.get("embedding_label"):
                        findings["embedding_label"] = res.get("embedding_label")
                    if res.get("embedding_key"):
                        findings["embedding_key"] = res.get("embedding_key")
                    if res.get("embedding_was_computed"):
                        findings["embedding_was_computed"] = True
                    for key, path in (res.get("figures") or {}).items():
                        figs[f"umap_{key}"] = path
                else:
                    log.warning(
                        f"UMAP figure generation failed: "
                        f"{res.get('error_type')} — {res.get('details', '')[:200]}"
                    )
            except Exception as e:
                log.warning(f"UMAP figure subprocess crashed: {e}")

    # 2. Per-cell-type DE summary bar ───────────────────────────────────
    bar_path = render_per_celltype_de_bar(
        findings, output_dir / "pseudobulk_de_per_celltype_bar.png"
    )
    if bar_path:
        figs["per_celltype_de_bar"] = bar_path

    # 3. Pathway dotplots ───────────────────────────────────────────────
    pw_figs = render_pathway_dotplots(findings, output_dir / "pathways")
    if pw_figs:
        figs["pathway_dotplots"] = pw_figs

    # 3b. Cell-cell communication figures ───────────────────────────────
    heat_path = render_cellcomm_heatmap(
        findings, output_dir / "cellcomm_heatmap.png"
    )
    if heat_path:
        figs["cellcomm_heatmap"] = heat_path
    bar_path_ccc = render_cellcomm_top_pairs_bar(
        findings, output_dir / "cellcomm_top_pairs.png"
    )
    if bar_path_ccc:
        figs["cellcomm_top_pairs"] = bar_path_ccc

    # 4. Trajectory figures (PAGA graph + DPT UMAP) ─────────────────────
    traj = findings.get("trajectory") or {}
    if traj.get("status") in ("done", "success") and env_manager is not None:
        traj_h5ad = traj.get("output_path") or h5ad_path
        if traj_h5ad:
            try:
                paga_res = env_manager.run_in_stack(
                    stack="rna",
                    script_path="aria/scripts/rna_figure_paga.py",
                    params={
                        "h5ad_path":  str(traj_h5ad),
                        "output_dir": str(output_dir / "trajectory"),
                        "groupby":    traj.get("groupby"),
                    },
                )
                if paga_res.get("status") == "success":
                    figs.update(paga_res.get("figures") or {})
                else:
                    log.warning(
                        f"PAGA figure failed: "
                        f"{paga_res.get('error_type')} — "
                        f"{paga_res.get('details', '')[:200]}"
                    )
            except Exception as e:
                log.warning(f"PAGA figure subprocess crashed: {e}")

            # DPT-coloured UMAP — only if dpt_pseudotime obs col exists.
            pt = traj.get("pseudotime", {}) or {}
            if pt.get("computed"):
                try:
                    dpt_res = env_manager.run_in_stack(
                        stack="rna",
                        script_path="aria/scripts/rna_figure_umap.py",
                        params={
                            "h5ad_path":  str(traj_h5ad),
                            "color_by":   ["dpt_pseudotime"],
                            "output_dir": str(output_dir / "trajectory"),
                        },
                    )
                    if dpt_res.get("status") == "success":
                        path = (dpt_res.get("figures") or {}).get(
                            "dpt_pseudotime"
                        )
                        if path:
                            figs["umap_dpt_pseudotime"] = path
                except Exception as e:
                    log.warning(f"DPT UMAP subprocess crashed: {e}")

    return findings
