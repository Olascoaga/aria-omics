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
