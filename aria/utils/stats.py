"""Small statistical helpers shared by ARIA scripts."""

from __future__ import annotations


def bh_correct(pvals):
    """Benjamini-Hochberg correction without requiring statsmodels."""
    import numpy as np

    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = adj
    return out


# ── Pre-registered FDR family (P1-2) ─────────────────────────────────────────
# The per-cluster vs global BH family must be fixed from the analysis plan BEFORE
# any p-values are computed, so the choice cannot be presented as a post-hoc,
# discovery-maximizing decision. These helpers make that pre-registration
# explicit and assert that the applied family never diverges from the declared
# strategy.

_VALID_FDR_STRATEGIES = ("per_cluster", "global")


def preregister_fdr_family(strategy) -> dict:
    """Normalize and pre-register a multiple-testing family.

    Returns a declaration dict (recorded in `methodology.json`) stating that the
    family was chosen before results were seen. Unknown/empty strategies fall
    back to the conservative `per_cluster` default.
    """
    s = str(strategy or "per_cluster").strip().lower()
    if s not in _VALID_FDR_STRATEGIES:
        s = "per_cluster"
    return {
        "fdr_strategy": s,
        "preregistered": True,
        "selected_before_results": True,
        "note": (
            "The multiple-testing family (per-cluster vs global BH) is fixed from "
            "the analysis plan before any p-values are computed; it is not chosen "
            "post-hoc to maximize the number of significant genes."
        ),
    }


def primary_fdr_column(strategy) -> str:
    """Map a (pre-registered) strategy to its primary adjusted-p-value column.

    Deterministic in the strategy ONLY — never in the data or result counts.
    `global` -> `padj_global`; anything else -> `padj_local`.
    """
    return "padj_global" if str(strategy).strip().lower() == "global" \
        else "padj_local"


def assert_fdr_family_not_post_hoc(declared_strategy, applied_column) -> None:
    """Integrity guard: the applied primary padj column must be the one the
    pre-registered strategy maps to. Raises ValueError on a post-hoc switch."""
    expected = primary_fdr_column(declared_strategy)
    if applied_column != expected:
        raise ValueError(
            f"FDR family integrity violation: pre-registered "
            f"'{declared_strategy}' maps to '{expected}', but the applied "
            f"primary column is '{applied_column}' (post-hoc switch)."
        )


# ── FDR across the bulk contrast family (P1-1c) ──────────────────────────────
# Bulk DE controls FDR within each contrast. When several contrasts are tested
# together, the family-wise error across the contrast family is not controlled.
# A pooled BH across all contrasts does, with the family pre-registered.

_VALID_CONTRAST_FAMILIES = ("per_contrast", "global")


def preregister_contrast_family(strategy) -> dict:
    """Pre-register the bulk contrast-FDR family before results are seen (P1-1c):
    `per_contrast` (BH within each contrast) or `global` (pooled BH across the
    contrast family). Unknown/empty falls back to the conservative-by-default
    `per_contrast`."""
    s = str(strategy or "per_contrast").strip().lower()
    if s not in _VALID_CONTRAST_FAMILIES:
        s = "per_contrast"
    return {
        "fdr_family": s,
        "preregistered": True,
        "selected_before_results": True,
        "note": (
            "The bulk FDR family (per-contrast vs global BH across the contrast "
            "family) is fixed from the analysis plan before any p-values are "
            "computed; it is not chosen post-hoc to maximize significant genes."
        ),
    }


def pooled_bh_across_groups(group_pvalues: dict) -> dict:
    """Pool p-values across groups (e.g. DE contrasts) and apply ONE BH
    correction over the whole family, mapping the adjusted values back per
    group. Controls FDR across the family rather than within each group.

    `group_pvalues`: {group: {gene: pvalue}}. Returns the same shape with
    pooled-BH adjusted p-values.
    """
    keys: list = []
    pvals: list = []
    for grp, genes in (group_pvalues or {}).items():
        for gene, p in (genes or {}).items():
            if p is None:
                continue
            keys.append((grp, gene))
            pvals.append(float(p))
    out: dict = {grp: {} for grp in (group_pvalues or {})}
    if not pvals:
        return out
    adj = bh_correct(pvals)
    for (grp, gene), a in zip(keys, adj):
        out[grp][gene] = float(a)
    return out


def contrast_family_significance(group_stats: dict, *, padj_max: float,
                                 lfc_min: float | None) -> dict:
    """Derive each contrast's family-significant set under pooled BH (P1-1c).

    `group_stats`: {contrast: {gene: {"pvalue": p, "log2fc": l}}}. Applies one
    pooled BH across the family, then calls a gene significant when its
    family-adjusted p-value < `padj_max` and, when `lfc_min` is not None,
    |log2fc| > `lfc_min`. Use `lfc_min=None` when the upstream Wald p-value
    already tested against an effect-size null (P1-1b). Returns
    {contrast: {"padj_family", "sig_genes", "n_sig", "n_up", "n_down"}}.
    """
    pvals = {
        c: {g: s.get("pvalue") for g, s in (genes or {}).items()
            if s.get("pvalue") is not None}
        for c, genes in (group_stats or {}).items()
    }
    family_padj = pooled_bh_across_groups(pvals)

    out: dict = {}
    for c, genes in (group_stats or {}).items():
        fp = family_padj.get(c, {})
        sig, n_up, n_down = [], 0, 0
        for g, s in (genes or {}).items():
            lfc = s.get("log2fc")
            passes_lfc = True if lfc_min is None else (
                lfc is not None and abs(lfc) > lfc_min
            )
            if fp.get(g, 1.0) < padj_max and passes_lfc:
                sig.append(g)
                if lfc is not None and lfc > 0:
                    n_up += 1
                elif lfc is not None and lfc < 0:
                    n_down += 1
        out[c] = {
            "padj_family": fp,
            "sig_genes": sig,
            "n_sig": len(sig),
            "n_up": n_up,
            "n_down": n_down,
        }
    return out
