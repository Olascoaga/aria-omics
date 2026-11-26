"""Small statistical helpers shared by ARIA scripts."""

from __future__ import annotations


def bh_correct(pvals):
    """Benjamini-Hochberg correction without requiring statsmodels."""
    import numpy as np

    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan)
    valid = np.isfinite(pvals)
    n = int(valid.sum())
    if n == 0:
        return out
    values = np.clip(pvals[valid], 0.0, 1.0)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty(n, dtype=float)
    restored[order] = adj
    out[valid] = restored
    return out


def adaptive_neglog10(
    pvals,
    *,
    headroom: float = 1.15,
    min_cap: float = 5.0,
    all_zero_floor: float = 10.0,
):
    """Data-adaptive ``-log10(p)`` with underflow capping for ``p == 0``.

    Returns ``(values, cap, has_zero_cap)``. Strictly positive finite ``p`` in
    ``(0, 1]`` map to ``-log10(p)``; float64-underflowed ``p == 0.0`` map to
    ``max(max_finite * headroom, max_finite + 1.0, min_cap)`` (or
    ``all_zero_floor`` when no strictly positive finite ``p`` exists), ensuring
    capped ultra-significant values always lie strictly above ``-log10(0.05)``.
    Non-finite or out-of-range inputs map to ``np.nan`` (never Python ``None``).
    """
    import numpy as np

    arr = np.asarray(pvals, dtype=float)
    valid = np.isfinite(arr) & (arr >= 0.0) & (arr <= 1.0)
    pos = arr[valid & (arr > 0.0)]
    if pos.size > 0:
        max_fin = float(-np.log10(np.min(pos)))
        cap = max(max_fin * float(headroom), max_fin + 1.0, float(min_cap))
    else:
        cap = float(all_zero_floor)
    zero_mask = valid & (arr == 0.0)
    out = np.full(arr.shape, np.nan, dtype=float)
    out[valid & (arr > 0.0)] = -np.log10(arr[valid & (arr > 0.0)])
    out[zero_mask] = cap
    return out, cap, bool(np.any(zero_mask))


def compute_cpm(counts, *, zero_policy: str = "zero"):
    """Compute Counts Per Million (CPM) using full-library totals per sample.

    Accepts a 1D sequence/array (single sample) or a 2D features × samples
    ``numpy.ndarray`` / ``pandas.DataFrame``. ``zero_policy`` controls samples
    with non-positive library totals:
      - ``"zero"`` (default): emit finite ``0.0`` for zero-depth columns;
      - ``"nan"``: emit ``np.nan`` for zero-depth columns;
      - ``"raise"``: raise ``ValueError("a sample has zero total counts")``.
    """
    import numpy as np

    is_df = hasattr(counts, "to_numpy") and hasattr(counts, "columns") and hasattr(counts, "index")
    arr = np.asarray(counts.to_numpy(dtype=float) if is_df else counts, dtype=float)
    if arr.ndim == 1:
        total = float(np.sum(arr))
        if not np.isfinite(total) or total <= 0.0:
            if zero_policy == "raise":
                raise ValueError("a sample has zero total counts")
            fill = np.nan if zero_policy == "nan" else 0.0
            return np.full_like(arr, fill, dtype=float)
        return (arr / total) * 1e6

    lib_sizes = np.sum(arr, axis=0)
    zero_mask = (~np.isfinite(lib_sizes)) | (lib_sizes <= 0.0)
    if np.any(zero_mask) and zero_policy == "raise":
        raise ValueError("a sample has zero total counts")
    safe_lib = np.where(zero_mask, 1.0, lib_sizes)
    cpm = (arr / safe_lib[None, :]) * 1e6
    if np.any(zero_mask):
        fill = np.nan if zero_policy == "nan" else 0.0
        cpm[:, zero_mask] = fill
    if is_df:
        import pandas as pd

        return pd.DataFrame(cpm, index=counts.index, columns=counts.columns)
    return cpm


def compute_log2_cpm(counts, *, zero_policy: str = "zero"):
    """Compute ``log2(CPM + 1)`` using full-library totals per sample."""
    import numpy as np

    cpm = compute_cpm(counts, zero_policy=zero_policy)
    if hasattr(cpm, "to_numpy") and hasattr(cpm, "columns"):
        import pandas as pd

        return pd.DataFrame(
            np.log2(np.maximum(cpm.to_numpy(dtype=float), 0.0) + 1.0),
            index=cpm.index,
            columns=cpm.columns,
        )
    return np.log2(np.maximum(cpm, 0.0) + 1.0)


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


def fdr_advanced_methods_disclosure() -> dict:
    """Honest disclosure for IHW and apeGLM s-values (P1-2 closure, ADR-027).

    ARIA does NOT implement covariate-weighted IHW or s-values, and discloses
    that rather than faking them. There is no validated pure-Python IHW
    (Ignatiadis-Huber is R/Bioconductor only) and pydeseq2 0.5.4 exposes no
    s-values; a hand-rolled covariate-weighted BH could silently break FDR
    control, which the no-fabrication policy (ADR-002) forbids. The primary FDR
    therefore stays the pre-registered Benjamini-Hochberg family. Recorded under
    `methodology.json["multiple_testing"]["advanced_methods"]`.
    """
    return {
        "ihw": {
            "status": "not_implemented",
            "available": False,
            "reason": (
                "No validated pure-Python IHW estimator; Independent Hypothesis "
                "Weighting (Ignatiadis-Huber) is R/Bioconductor only. A "
                "hand-rolled covariate-weighted BH could silently break FDR "
                "control (ADR-002), so it is not shipped."
            ),
        },
        "s_values": {
            "status": "not_available",
            "available": False,
            "reason": "pydeseq2 0.5.4 (apeGLM) exposes no s-values.",
        },
        "primary_method": "pre-registered Benjamini-Hochberg (local/global dual)",
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
    import math

    keys: list = []
    pvals: list = []
    for grp, genes in (group_pvalues or {}).items():
        for gene, p in (genes or {}).items():
            if p is None or isinstance(p, bool):
                continue
            keys.append((grp, gene))
            value = float(p)
            pvals.append(value if math.isfinite(value) else float("nan"))
    out: dict = {grp: {} for grp in (group_pvalues or {})}
    if not pvals:
        return out
    adj = bh_correct(pvals)
    for (grp, gene), a in zip(keys, adj, strict=True):
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
    import math

    pvals = {
        c: {g: s.get("pvalue") for g, s in (genes or {}).items()
            if s.get("pvalue") is not None and not isinstance(s.get("pvalue"), bool)}
        for c, genes in (group_stats or {}).items()
    }
    family_padj = pooled_bh_across_groups(pvals)

    out: dict = {}
    for c, genes in (group_stats or {}).items():
        fp = family_padj.get(c, {})
        sig, n_up, n_down = [], 0, 0
        for g, s in (genes or {}).items():
            raw_lfc = s.get("log2fc")
            lfc: float | None = None
            lfc_invalid = False
            if raw_lfc is not None:
                if isinstance(raw_lfc, bool):
                    lfc_invalid = True
                else:
                    try:
                        val = float(raw_lfc)
                        if math.isfinite(val):
                            lfc = val
                        else:
                            lfc_invalid = True
                    except (TypeError, ValueError):
                        lfc_invalid = True
            if lfc_invalid:
                continue
            passes_lfc = True if lfc_min is None else (
                lfc is not None and abs(lfc) > lfc_min
            )
            padj_val = fp.get(g, 1.0)
            if math.isfinite(padj_val) and padj_val < padj_max and passes_lfc:
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
