"""Deterministic robustness summaries for methodology provenance."""

from __future__ import annotations

from typing import Any

from aria.agents import _narrative_scrna


def build_robustness_multiverse(agent_results: dict[str, Any] | None) -> dict:
    """Summarize available multiverse checks without hidden reruns.

    P-MULTIVERSE originally proposed FDR strategy x composition-covariate
    reruns. ARIA already computes both local and global BH families in each
    pseudobulk block; this manifest records the genes stable across those
    families and states the realized composition-covariate state explicitly.
    """
    sc = (agent_results or {}).get("scrna_agent", {})
    findings = _narrative_scrna.unwrap_scrna_findings(sc)
    pb = findings.get("pseudobulk_de") or {}
    entries = []
    for group, info in (pb.get("per_group", {}) or {}).items():
        for comparison, comp in (info.get("per_comparison", {}) or {}).items():
            if comp.get("status") != "success":
                continue
            mv = comp.get("robustness_multiverse") or {}
            if not mv:
                local = int(comp.get("n_significant_local", 0) or 0)
                global_ = int(comp.get("n_significant_global", 0) or 0)
                stable = min(local, global_)
            else:
                local = int(mv.get("n_local", comp.get("n_significant_local", 0)) or 0)
                global_ = int(mv.get("n_global", comp.get("n_significant_global", 0)) or 0)
                stable = int(mv.get("stable_significant_genes", min(local, global_)) or 0)
            entries.append({
                "group": group,
                "comparison": comparison,
                "stable_significant_genes": stable,
                "n_local_fdr": local,
                "n_global_fdr": global_,
                "composition_covariate": (
                    "included" if comp.get("corrected_for_composition")
                    else "not_included"
                ),
                "composition_axis_rerun": False,
            })
    return {
        "status": "available" if entries else "not_available",
        "method": "FDR-family stability over local/global BH pseudobulk calls",
        "entries": entries,
        "n_entries": len(entries),
        "note": (
            "Composition on/off is not rerun implicitly; the manifest records "
            "the realized composition-covariate state for each block."
        ),
    }
