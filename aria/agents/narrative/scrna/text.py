"""scRNA narrative text: descriptions, summary, integrated interpretation, methods.

Extracted verbatim from aria/agents/_narrative_scrna.py (A7); behavior pinned by
tests/test_narrative_scrna_contract.py.
"""
from __future__ import annotations

from typing import Optional

from aria.agents.narrative.scrna._common import *  # noqa: F401,F403


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
        # T9: disclose the compositional-dependence assumption of the per-type
        # CLR + BH approach (CLR values are sum-to-zero; per-type tests are not
        # independent). The runtime method is unchanged.
        if method == "donor_clr_ols_hc3":
            lines.append(
                "Compositional-dependence assumption: per-cell-type CLR tests "
                "are fit independently and BH-corrected together, but CLR values "
                "are sum-to-zero, so per-type abundance directions are not "
                "independent and should be read jointly as one compositional "
                "shift. A joint compositional model (e.g. scCODA or propeller) "
                "would model this dependence directly."
            )
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
        amb = qc.get("ambient_correction") or {}
        if amb.get("ran"):
            lines.append(
                f"Ambient-RNA decontamination was applied "
                f"({amb.get('method', 'SoupX/decontX')}) before downstream "
                f"analysis."
            )
        else:
            lines.append(
                "Ambient-RNA decontamination (SoupX/decontX) is available as an "
                "optional step and was not applied in this run; potential "
                "ambient contamination was instead screened by a cross-cluster "
                "marker-ubiquity check."
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
    ct_typist = ct.get("celltypist") or {}
    model_used = ct.get("model_used") or ct_typist.get("model_used")
    if model_used:
        line = (
            f"Cell-type annotation used CellTypist with model "
            f"'{model_used}', assigning a majority label per cluster."
        )
        # N-ANNO3: disclose a silent immune-default fallback in Methods.
        if (ct_typist.get("model_source") == "default_immune_fallback"
                or ct_typist.get("model_warning")):
            line += (
                " No tissue or model hint was supplied, so this immune-default "
                "model may not match the tissue; cell-type labels are reported "
                "as model-derived and potentially mismatched."
            )
        lines.append(line)

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

