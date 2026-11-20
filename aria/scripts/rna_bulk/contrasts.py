"""Bulk RNA-seq contrast helpers: slug, top-gene formatting, contrast suggestion
and overlap (A7 split of rna_bulk_de.py; bodies verbatim, behavior-preserving).

Re-exported from aria.scripts.rna_bulk_de."""
from __future__ import annotations
import re
import math

from aria.scripts.rna_bulk.gtf_io import _strip_gene_version


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
            lfc, padj = float(row["log2FoldChange"]), float(row["padj"])
            if not math.isfinite(lfc) or not math.isfinite(padj):
                continue
            clean_id = _strip_gene_version(str(g))
            symbol   = sm.get(clean_id, "")
            top.append({
                "gene":      g,                    # Ensembl ID (or whatever)
                "symbol":    symbol or g,          # HGNC symbol if known, else fallback to ID
                "log2fc":    round(lfc, 3),
                "padj":      padj,
                "direction": "up" if lfc > 0 else "down",
            })
        except Exception:
            continue
    return top


def _suggest_contrasts(metadata, design_factor: str) -> list[dict]:
    """Suggest candidate contrasts without authorizing execution.

    P0-5: suggestions are display-only. The caller must pass one or more of
    them back explicitly as ``contrasts`` before DE can run.
    """
    groups = sorted({str(g) for g in metadata[design_factor].unique()})

    if len(groups) < 2:
        return []

    # Identify control
    ctrl_keywords = ["wt", "wildtype", "control", "ctrl",
                     "vehicle", "dmso", "untreated", "mock",
                     "normal", "healthy", "baseline", "scramble"]
    control = None
    for kw in ctrl_keywords:
        for g in groups:
            if str(g).lower() == kw:
                control = g
                break
        if control:
            break
    if not control:
        for kw in ctrl_keywords:
            for g in groups:
                if kw in str(g).lower():
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

    # Use full IDs when every contrast supplies them, keeping unmapped genes
    # and a consistent namespace. Legacy records fall back to full symbols.
    gene_sets = {}
    use_ids = all("all_sig_gene_ids" in c for c in successful)
    for c in successful:
        if use_ids:
            gene_sets[c["name"]] = set(c["all_sig_gene_ids"] or [])
        elif "all_sig_genes" in c:
            gene_sets[c["name"]] = set(c["all_sig_genes"] or [])
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

