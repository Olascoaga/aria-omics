"""Deterministic prose composition for narrative blocks."""

from __future__ import annotations

import math

from aria.agents.narrative.types import NarrativeBlock


def compose_block_prose(block: NarrativeBlock) -> str:
    """Return reviewer-readable prose for one structured narrative block."""
    if block.status != "success":
        return _compose_non_success(block)

    if block.analysis == "differential_expression":
        return _compose_de(block)
    if block.analysis == "pathway_enrichment":
        return _compose_pathway(block)
    if block.analysis == "gsea_preranked":
        return _compose_gsea(block)
    if block.analysis == "differential_abundance":
        return _compose_composition(block)
    if block.analysis == "pseudobulk_de":
        return _compose_pseudobulk(block)
    if block.analysis == "pseudobulk_pathways":
        return _compose_pathway(block)
    if block.analysis == "cell_communication":
        return _compose_cellcomm(block)
    if block.analysis == "trajectory":
        return _compose_trajectory(block)
    if block.analysis == "paired_qc":
        return _compose_paired_qc(block)
    if block.analysis == "wnn":
        return _compose_wnn(block)
    if block.analysis == "mofa":
        return _compose_mofa(block)
    if block.analysis == "peak_to_gene":
        return _compose_peak_to_gene(block)
    if block.analysis == "eregulons":
        return _compose_eregulons(block)
    if block.analysis == "cross_modal_triangulation":
        return _compose_triangulation(block)
    if block.analysis in {"qc", "sample_qc"} or block.block_type == "qc":
        return _compose_qc(block)
    if block.analysis == "power":
        return _compose_power(block)
    return _sentence(block.claim)


def _compose_non_success(block: NarrativeBlock) -> str:
    if block.error:
        return (
            f"{block.title} did not produce an interpretable result: "
            f"{block.error}."
        )
    if block.claim:
        return _sentence(block.claim)
    return f"{block.title} did not produce an interpretable result."


def _compose_qc(block: NarrativeBlock) -> str:
    samples = _ev(block, "samples")
    before = _ev(block, "cells before QC")
    after = _ev(block, "cells after QC")
    removed = _ev(block, "outliers removed")
    flagged = _ev(block, "outliers flagged")
    removed_primary = _ev(block, "outliers removed primary")
    removed_sensitivity = _ev(block, "outliers removed sensitivity")
    if before is not None and after is not None:
        return (
            f"{block.title} retained {after} of {before} cells for downstream "
            "analysis, so later sections should be read as post-QC results."
        )
    if samples is not None:
        if flagged is not None:
            outlier_text = (
                f", flagged {flagged} outlier sample(s), removed "
                f"{removed_primary if removed_primary is not None else 0} "
                "in the primary analysis"
            )
            if removed_sensitivity is not None:
                outlier_text += (
                    f", and removed {removed_sensitivity} in sensitivity"
                )
        else:
            outlier_text = (
                f" and removed {removed} outlier sample(s)"
                if removed is not None else ""
            )
        return (
            f"{block.title} evaluated {samples} sample(s){outlier_text} before "
            "differential analysis."
        )
    return _sentence(block.claim)


def _compose_de(block: NarrativeBlock) -> str:
    n_sig = _ev(block, "DE genes", block.metrics.get("n_significant"))
    up = _ev(block, "up genes", block.metrics.get("n_upregulated"))
    down = _ev(block, "down genes", block.metrics.get("n_downregulated"))
    genes = _top_gene_phrases(block)
    direction = ""
    if up is not None and down is not None:
        direction = f" ({up} higher and {down} lower in the numerator group)"
    gene_text = f" The strongest reported genes include {genes}." if genes else ""
    return (
        f"{block.title} identified {n_sig} differentially expressed genes"
        f"{direction}.{gene_text}"
    )


def _compose_pseudobulk(block: NarrativeBlock) -> str:
    n_global = _ev(
        block, "global-FDR DE genes", block.metrics.get("n_significant_global")
    )
    n_local = _ev(
        block, "local-FDR DE genes", block.metrics.get("n_significant_local")
    )
    up = _ev(block, "up genes")
    down = _ev(block, "down genes")
    genes = _top_gene_phrases(block)
    local_text = (
        f"; {n_local} passed the per-block local FDR threshold"
        if n_local is not None and n_local != n_global else ""
    )
    direction = ""
    if up is not None and down is not None:
        direction = f" ({up} higher, {down} lower)"
    gene_text = f" Top genes in this block include {genes}." if genes else ""
    correction = (
        " The DE model included a composition covariate."
        if block.metrics.get("corrected_for_composition") else ""
    )
    return (
        f"{block.title} contributed {n_global} global-FDR DE genes{direction}"
        f"{local_text}.{gene_text}{correction}"
    )


def _compose_composition(block: NarrativeBlock) -> str:
    tested = _ev(block, "cell types tested")
    shifted = _ev(block, "significant shifts", block.metrics.get("n_significant"))
    shifts = []
    for ev in block.evidence:
        if "abundance direction" in ev.label and ev.value:
            name = ev.label.replace(" abundance direction", "")
            shifts.append(f"{name} {ev.value}")
    examples = f" The largest listed shifts were {', '.join(shifts[:5])}." \
        if shifts else ""
    if tested is not None:
        return (
            f"{block.title} detected abundance changes in {shifted} of "
            f"{tested} tested cell type(s).{examples}"
        )
    return _sentence(block.claim)


def _compose_pathway(block: NarrativeBlock) -> str:
    n_terms = _ev(block, "enriched terms", block.metrics.get("n_terms"))
    terms = _term_values(block)
    term_text = f" Leading terms include {', '.join(terms[:5])}." if terms else ""
    return (
        f"{block.title} found {n_terms} enriched pathway or gene-set term(s)."
        f"{term_text} These terms summarize the DE gene list rather than "
        "establishing causal pathway activity."
    )


def _compose_gsea(block: NarrativeBlock) -> str:
    n_terms = _ev(block, "FDR<0.25 pathways", block.metrics.get("n_pathways"))
    top = block.metrics.get("top_pathways") or []
    pieces = []
    for row in top[:3]:
        if row.get("numeric_unstable"):
            continue
        term = row.get("term")
        if not term:
            continue
        nes = row.get("nes")
        fdr = row.get("fdr")
        stats = []
        if _is_finite_number(nes):
            # Use the same 6-significant-figure normalization the strict evidence
            # verifier applies, so the displayed NES matches the value on the
            # evidence card (`:.3g` rounded to 5.43 while the card held 5.4321,
            # which made the verifier reject the prose and abort the report).
            stats.append(f"NES={float(nes):.6g}")
        if _is_finite_number(fdr):
            if float(fdr) <= 0:
                stats.append("FDR below display precision")
            else:
                stats.append(f"FDR={float(fdr):.6g}")
        pieces.append(f"{term} ({', '.join(stats)})" if stats else str(term))
    top_text = f" The top ranked signals were {', '.join(pieces)}." if pieces else ""
    if block.metrics.get("n_numeric_unstable"):
        top_text += (
            " One or more GSEA rows had non-finite numerical values and were "
            "excluded from the ranked-signal prose."
        )
    return (
        f"{block.title} used the full log2FC-ranked gene list and found "
        f"{n_terms} pathway(s) at the standard exploratory GSEA threshold "
        f"of FDR < 0.25.{top_text}"
    )


def _is_finite_number(value) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _compose_cellcomm(block: NarrativeBlock) -> str:
    interactions = _ev(block, "interactions", block.metrics.get("n_interactions"))
    cell_types = _ev(block, "cell types")
    pairs = []
    for ev in block.evidence:
        if "->" in ev.label and ev.value:
            pairs.append(f"{ev.label} {ev.value}")
    pair_text = f" Top candidates include {', '.join(pairs[:5])}." if pairs else ""
    if cell_types is not None:
        return (
            f"{block.title} reported {interactions} non-autocrine "
            f"ligand-receptor candidate(s) across {cell_types} cell type(s)."
            f"{pair_text}"
        )
    return _sentence(block.claim)


def _compose_trajectory(block: NarrativeBlock) -> str:
    edges = _ev(block, "PAGA connections")
    strong = _ev(block, "strong PAGA edges")
    order = _ev(block, "DPT ordered groups")
    order_text = f" DPT ordered groups as {order}." if order else ""
    return (
        f"{block.title} summarized manifold connectivity with {edges} PAGA "
        f"connection(s), including {strong} strong edge(s).{order_text}"
    )


def _compose_power(block: NarrativeBlock) -> str:
    min_power = block.metrics.get("min_power")
    max_power = block.metrics.get("max_power")
    if isinstance(min_power, (int, float)) and isinstance(max_power, (int, float)):
        return (
            f"{block.title} estimated approximate power between "
            f"{min_power:.0%} and {max_power:.0%} across analyzable contrasts."
        )
    return _sentence(block.claim)


def _compose_paired_qc(block: NarrativeBlock) -> str:
    before = _ev(block, "cells before QC", block.metrics.get("n_cells_before"))
    after = _ev(block, "cells after QC", block.metrics.get("n_cells_after"))
    pct = _ev(block, "removed fraction pct", block.metrics.get("pct_removed"))
    pct_str = str(pct) if str(pct).endswith("%") else f"{pct}%"
    genes = _ev(block, "RNA genes", block.metrics.get("n_genes"))
    peaks = _ev(block, "ATAC peaks", block.metrics.get("n_peaks"))
    matched = _ev(block, "label join matched cells")
    join_txt = (
        f" Reference barcode label join matched {matched} shared cells."
        if matched is not None
        else ""
    )
    return (
        f"{block.title} retained {after} of {before} verified shared barcodes "
        f"({pct_str} removed) spanning {genes} RNA genes and {peaks} ATAC peaks."
        f"{join_txt}"
    )


def _compose_wnn(block: NarrativeBlock) -> str:
    n_cells = _ev(block, "shared cells", block.metrics.get("n_cells"))
    n_joint = _ev(block, "WNN joint clusters", block.metrics.get("n_joint_clusters"))
    n_rna = _ev(block, "RNA-only clusters", block.metrics.get("n_rna_only_clusters"))
    n_atac = _ev(block, "ATAC-only clusters", block.metrics.get("n_atac_only_clusters"))
    ari_r = _ev(block, "WNN vs RNA ARI", block.metrics.get("wnn_vs_rna_ari"))
    ari_a = _ev(block, "WNN vs ATAC ARI", block.metrics.get("wnn_vs_atac_ari"))
    nmi_r = _ev(block, "WNN vs RNA NMI", block.metrics.get("wnn_vs_rna_nmi"))
    nmi_a = _ev(block, "WNN vs ATAC NMI", block.metrics.get("wnn_vs_atac_nmi"))
    n_disc = _ev(block, "discordant clusters", block.metrics.get("n_discordant_clusters"))
    guard = _ev(block, "WNN guardrail status", block.metrics.get("guardrail_status"))
    rw = _ev(block, "mean RNA weight", block.metrics.get("mean_rna_weight"))
    aw = _ev(block, "mean ATAC weight", block.metrics.get("mean_atac_weight"))
    w_txt = (
        f" Mean modality weights were {rw} for RNA and {aw} for ATAC."
        if (rw is not None and aw is not None)
        else ""
    )
    return (
        f"{block.title} integrated {n_cells} shared cells into {n_joint} joint "
        f"clusters compared with {n_rna} RNA-only and {n_atac} ATAC-only clusters "
        f"(WNN vs RNA ARI {ari_r}, NMI {nmi_r}; WNN vs ATAC ARI {ari_a}, "
        f"NMI {nmi_a}; {n_disc} discordant cluster(s); guardrail status {guard})."
        f"{w_txt}"
    )


def _compose_mofa(block: NarrativeBlock) -> str:
    n_cells = _ev(block, "MOFA cells evaluated", block.metrics.get("n_cells"))
    n_fac = _ev(block, "fitted factors", block.metrics.get("n_factors"))
    n_act = _ev(block, "active factors", block.metrics.get("n_active_factors"))
    n_sh = _ev(block, "shared factors", block.metrics.get("n_shared_factors"))
    n_sp = _ev(
        block,
        "modality-specific factors",
        block.metrics.get("n_modality_specific_factors"),
    )
    n_ex = _ev(
        block,
        "excluded confounder factors",
        block.metrics.get("n_excluded_factors"),
    )
    rna_r2 = _ev(block, "RNA non-additive R2 pct", block.metrics.get("rna_r2_pct"))
    atac_r2 = _ev(block, "ATAC non-additive R2 pct", block.metrics.get("atac_r2_pct"))
    rna_str = str(rna_r2) if str(rna_r2).endswith("%") else f"{rna_r2}%"
    atac_str = str(atac_r2) if str(atac_r2).endswith("%") else f"{atac_r2}%"
    return (
        f"{block.title} evaluated {n_cells} shared cells across {n_fac} fitted "
        f"factors, retaining {n_act} active factors ({n_sh} shared and {n_sp} "
        f"modality-specific) while quarantining {n_ex} technical or cell-cycle "
        f"confounder factor(s); total non-additive R2 was {rna_str} in RNA and "
        f"{atac_str} in ATAC."
    )


def _compose_peak_to_gene(block: NarrativeBlock) -> str:
    n_cells = _ev(block, "shared cells", block.metrics.get("n_cells"))
    dist_kb = _ev(block, "cis window kb", block.metrics.get("distance_kb"))
    n_tested = _ev(block, "tested peak-gene pairs", block.metrics.get("n_tested"))
    n_links = _ev(block, "BH-FDR cis links", block.metrics.get("n_links"))
    n_pos = _ev(block, "positive cis links", block.metrics.get("n_positive_correlations"))
    n_neg = _ev(block, "negative cis links", block.metrics.get("n_negative_correlations"))
    n_prom = _ev(block, "promoter-proximal links", block.metrics.get("n_promoter_proximal"))
    n_dist = _ev(block, "distal cis links", block.metrics.get("n_distal_cis"))
    candidates = []
    for ev in block.evidence:
        if ev.label.startswith("candidate cis link ") and ev.value is not None:
            gene = ev.label[len("candidate cis link "):]
            candidates.append(f"{gene} (r={ev.value})")
    cand_txt = (
        f" Top associative cis candidates include {', '.join(candidates[:4])}."
        if candidates
        else ""
    )
    return (
        f"{block.title} tested {n_tested} cis pairs within {dist_kb} kb across "
        f"{n_cells} shared cells and identified {n_links} BH-FDR candidate cis "
        f"associations ({n_pos} positive, {n_neg} negative; {n_prom} "
        f"promoter-proximal, {n_dist} distal cis).{cand_txt}"
    )


def _compose_eregulons(block: NarrativeBlock) -> str:
    n_cells = _ev(block, "shared cells", block.metrics.get("n_shared_cells"))
    n_tested = _ev(
        block,
        "candidate modules tested",
        block.metrics.get("n_candidate_modules_tested"),
    )
    n_ereg = _ev(
        block,
        "BH-FDR candidate eRegulons",
        _ev(block, "candidate eRegulons", block.metrics.get("n_eregulons")),
    )
    n_trip = _ev(block, "TF-enhancer-gene triplets", block.metrics.get("n_triplets"))
    n_tfs = _ev(block, "unique TFs", block.metrics.get("n_unique_tfs"))
    n_pos = _ev(
        block,
        "concordant pos_tf_pos_p2g modules",
        block.metrics.get("n_pos_tf_pos_p2g"),
    )
    n_neg = _ev(
        block,
        "inverse neg_tf_pos_p2g modules",
        block.metrics.get("n_neg_tf_pos_p2g"),
    )
    n_par = _ev(
        block,
        "paralog-ambiguous modules",
        block.metrics.get("n_paralog_ambiguous"),
    )
    n_dep = _ev(
        block,
        "depth-confounded modules",
        block.metrics.get("n_depth_confounded"),
    )
    mods = []
    for ev in block.evidence:
        if ev.label.startswith("candidate eRegulon ") and ev.value is not None:
            tf = ev.label[len("candidate eRegulon "):]
            mods.append(f"{tf} (r={ev.value})")
    mod_txt = (
        f" Leading associative candidate modules include {', '.join(mods[:4])}."
        if mods
        else ""
    )
    return (
        f"{block.title} evaluated {n_tested} motif-linked modules across "
        f"{n_cells} shared cells and retained {n_ereg} BH-FDR candidate "
        f"eRegulons ({n_pos} pos_tf_pos_p2g, {n_neg} neg_tf_pos_p2g) spanning "
        f"{n_tfs} unique TFs and {n_trip} triplets, with {n_par} "
        f"paralog-ambiguous and {n_dep} depth-confounded module(s) flagged."
        f"{mod_txt}"
    )


def _compose_triangulation(block: NarrativeBlock) -> str:
    n_stages = _ev(block, "integration stages executed", block.metrics.get("n_stages"))
    n_nodes = _ev(block, "evidence graph nodes", block.metrics.get("n_nodes"))
    n_edges = _ev(block, "evidence graph edges", block.metrics.get("n_edges"))
    n_axes = _ev(block, "triangulated axes", block.metrics.get("n_axes_total"))
    n_conv = _ev(block, "convergent axes", block.metrics.get("n_convergent_total"))
    n_assoc = _ev(
        block,
        "associative cis candidates",
        block.metrics.get("associative_cis_candidate"),
    )
    n_disc = _ev(block, "discordant axes", block.metrics.get("discordant_signal"))
    n_cav = _ev(block, "caveated ambiguous axes", block.metrics.get("caveated_ambiguous"))
    guard = _ev(
        block,
        "WNN guardrail status",
        block.metrics.get("wnn_guardrail_status"),
    )
    return (
        f"{block.title} synthesized {n_stages} executed integration stage(s) "
        f"into {n_nodes} nodes and {n_edges} non-causal associative edges across "
        f"{n_axes} axes ({n_conv} convergent, {n_assoc} associative cis "
        f"candidates, {n_disc} discordant, and {n_cav} caveated ambiguous; "
        f"WNN guardrail {guard})."
    )


def _ev(block: NarrativeBlock, label: str, default=None):
    for ev in block.evidence:
        if ev.label == label:
            return ev.value
    return default


def _top_gene_phrases(block: NarrativeBlock) -> str:
    genes = []
    for ev in block.evidence:
        if not ev.label.startswith("top gene "):
            continue
        gene = ev.label.replace("top gene ", "")
        if ev.value is None:
            genes.append(gene)
        else:
            genes.append(f"{gene} (log2FC {ev.value})")
    return ", ".join(genes[:6])


def _term_values(block: NarrativeBlock) -> list[str]:
    terms = []
    for ev in block.evidence:
        if ev.label == "enriched terms":
            continue
        label = ev.label.lower()
        if " term" in label and ev.value:
            terms.append(str(ev.value))
    return terms


def _sentence(text: str) -> str:
    text = str(text or "").strip()
    if not text:
        return ""
    return text if text.endswith((".", "!", "?")) else f"{text}."
