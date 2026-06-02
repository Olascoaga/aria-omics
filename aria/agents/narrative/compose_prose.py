"""Deterministic prose composition for narrative blocks."""

from __future__ import annotations

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
        term = row.get("term")
        if not term:
            continue
        nes = row.get("nes")
        fdr = row.get("fdr")
        stats = []
        if nes is not None:
            stats.append(f"NES={nes}")
        if fdr is not None:
            stats.append(f"FDR={fdr}")
        pieces.append(f"{term} ({', '.join(stats)})" if stats else str(term))
    top_text = f" The top ranked signals were {', '.join(pieces)}." if pieces else ""
    return (
        f"{block.title} used the full log2FC-ranked gene list and found "
        f"{n_terms} pathway(s) at the standard exploratory GSEA threshold "
        f"of FDR < 0.25.{top_text}"
    )


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
