"""Shared analysis-threshold value object (P0-7).

The user confirms significance thresholds once (CP3); the orchestrator records
them in `exp_context["global_padj"]` / `["global_lfc"]`. Scripts must NOT carry
their own loose threshold literals — that silently overrides the confirmed
decision. `AnalysisThresholds` is the single, additive vehicle that resolves the
confirmed thresholds (with documented fallbacks) and hands them to a script's
parameter dict.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AnalysisThresholds:
    """Significance / inclusion thresholds for a DE-style analysis."""

    padj: float = 0.05
    log2fc: float = 0.5
    min_cells: int = 10
    min_replicates: int = 3
    modality: str = "generic"

    @classmethod
    def from_exp_context(cls, exp_ctx: dict | None, modality: str = "generic", *,
                         log2fc_default: float = 0.5,
                         min_cells: int = 10,
                         min_replicates: int = 3) -> "AnalysisThresholds":
        """Resolve thresholds from the confirmed experiment context.

        `global_padj` / `global_lfc` are the user-confirmed CP3 values. When they
        are absent (direct/legacy invocation), fall back to the documented
        defaults (`log2fc_default` lets each modality keep its prior default —
        e.g. 0.5 for scRNA pseudobulk, 1.0 for bulk).
        """
        ctx = exp_ctx or {}

        def _num(key: str, default: float) -> float:
            try:
                return float(ctx.get(key, default))
            except (TypeError, ValueError):
                return float(default)

        return cls(
            padj=_num("global_padj", 0.05),
            log2fc=_num("global_lfc", log2fc_default),
            min_cells=int(min_cells),
            min_replicates=int(min_replicates),
            modality=str(modality),
        )

    def as_pseudobulk_params(self) -> dict:
        """Threshold keys for `rna_pseudobulk_de.py` (additive — merge into the
        existing params dict; no loose literals in the agent)."""
        return {
            "padj_max": self.padj,
            "lfc_min": self.log2fc,
            "min_cells_per_pseudosample": self.min_cells,
            "min_replicates_per_condition": self.min_replicates,
        }
