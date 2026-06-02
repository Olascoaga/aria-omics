"""Deterministic robustness summaries for methodology provenance."""

from __future__ import annotations

from typing import Any

from aria.agents import _narrative_scrna


def build_robustness_multiverse(agent_results: dict[str, Any] | None) -> dict:
    """Summarize available multiverse checks without hidden reruns.

    P-MULTIVERSE originally proposed FDR strategy x composition-covariate
    reruns. ARIA computes both local and global BH families in each pseudobulk
    block from the same p-value table; this manifest records the true gene-ID
    intersection when the block-level submanifest is available and states the
    realized composition-covariate state explicitly. It never substitutes
    ``min(local, global)`` as if it were an intersection.
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
            local = int(mv.get("n_local", comp.get("n_significant_local", 0)) or 0)
            global_ = int(mv.get("n_global", comp.get("n_significant_global", 0)) or 0)
            has_stability = (
                bool(mv)
                and mv.get("stability_basis") == "gene_id_intersection"
                and mv.get("stable_significant_genes") is not None
            )
            stable = (
                int(mv.get("stable_significant_genes") or 0)
                if has_stability else None
            )
            entries.append({
                "group": group,
                "comparison": comparison,
                "stable_significant_genes": stable,
                "stability_status": (
                    "computed" if has_stability else "not_computed"
                ),
                "stability_basis": (
                    mv.get("stability_basis") if has_stability else None
                ),
                "stable_gene_ids": list(mv.get("stable_gene_ids") or []),
                "stable_gene_ids_truncated": bool(
                    mv.get("stable_gene_ids_truncated", False)
                ),
                "fdr_axis_evaluated": bool(
                    mv.get("fdr_axis_evaluated", has_stability)
                ),
                "fdr_family_variants": mv.get("fdr_family_variants", {}),
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
        "method": (
            "gene-ID intersection across local/global BH pseudobulk FDR families"
        ),
        "entries": entries,
        "n_entries": len(entries),
        "note": (
            "FDR-family stability is reported only when a block records the "
            "actual gene-ID intersection. Composition on/off is not rerun "
            "implicitly; the manifest records the realized composition-covariate "
            "state for each block."
        ),
    }
