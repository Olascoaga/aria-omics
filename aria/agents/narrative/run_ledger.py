"""Planned-vs-Run ledger (P-LEDGER, audit 2026-05-29).

A deterministic per-run manifest that reconciles what the analysis plan called
for against what execution actually produced. Its job is to make any
plan -> execution divergence visible — specifically the dispatch-integrity gap
that produced the PBMC thin report, where pseudobulk DE was planned and approved
but silently skipped by a free-text keyword gate.

This is pure bookkeeping over the structured plan (DesignIntelligence) and the
structured agent results; it asserts no biology and runs no LLM.
"""

from __future__ import annotations

from typing import Any


# Canonical analyses, the plan phrases that imply them, and the agent-result
# finding keys that prove they ran. Technical/process vocabulary only (no
# biological content — ADR-011).
_SCRNA_ANALYSES: list[dict[str, Any]] = [
    {"key": "qc", "label": "Quality control",
     "plan_kw": ["qc", "quality control"], "finding_keys": ["qc"]},
    {"key": "clustering", "label": "Clustering",
     "plan_kw": ["leiden", "cluster"], "finding_keys": ["clustering"]},
    {"key": "annotation", "label": "Cell-type annotation",
     "plan_kw": ["celltypist", "annotat", "cell type", "cell-type", "marker"],
     "finding_keys": ["cell_types"]},
    {"key": "pseudobulk_de", "label": "Pseudobulk differential expression",
     "plan_kw": ["pseudobulk", "deseq", "donor-level", "donor level"],
     "finding_keys": ["pseudobulk_de"]},
    {"key": "differential_abundance", "label": "Differential abundance",
     "plan_kw": ["abundance", "composition"],
     "finding_keys": ["differential_abundance"]},
    {"key": "pathway_enrichment", "label": "Pathway enrichment (ORA/GSEA)",
     "plan_kw": ["pathway", "ora", "enrichment", "gsea"],
     "finding_keys": ["pseudobulk_pathways", "pathways"]},
    {"key": "cell_communication", "label": "Cell-cell communication",
     "plan_kw": ["liana", "ligand", "communication", "cell-cell", "cell cell"],
     "finding_keys": ["cell_communication"]},
    {"key": "trajectory", "label": "Trajectory (PAGA/DPT)",
     "plan_kw": ["trajectory", "paga", "dpt", "pseudotime"],
     "finding_keys": ["trajectory"]},
    {"key": "rna_velocity", "label": "RNA velocity",
     "plan_kw": ["velocity", "scvelo"], "finding_keys": ["velocity"]},
]


def _plan_phrases(exp_ctx: dict) -> list[str]:
    di = (exp_ctx or {}).get("design_intelligence", {}) or {}
    phrases: list[str] = []
    for bucket in ("recommended", "optional"):
        for item in di.get(bucket, []) or []:
            phrases.append(str(item).lower())
    return phrases


def _is_planned(spec: dict, phrases: list[str]) -> bool:
    return any(kw in phrase for phrase in phrases for kw in spec["plan_kw"])


def _status_from_finding(val: Any) -> tuple[str, str | None]:
    """Map an agent-result finding value to (status, reason)."""
    if val is None:
        return "not_run", None
    if isinstance(val, dict):
        st = str(val.get("status") or "").lower()
        if st == "skipped":
            return "skipped", val.get("reason")
        if st == "error":
            return "error", (val.get("error_type") or val.get("details"))
        if st in ("success", "ok"):
            return "ran", None
        # Structured dict with content but no explicit status -> it ran.
        return ("ran", None) if val else ("not_run", None)
    if isinstance(val, (list, tuple)):
        return ("ran", None) if len(val) else ("not_run", None)
    return "ran", None


def _scrna_findings(agent_results: dict) -> dict:
    sc = (agent_results or {}).get("scrna_agent", {})
    try:
        from aria.agents import _narrative_scrna
        return _narrative_scrna.unwrap_scrna_findings(sc) or {}
    except Exception:
        if isinstance(sc, dict):
            return sc.get("findings", sc) or {}
        return {}


def build_run_ledger(exp_ctx: dict, agent_results: dict) -> dict:
    """Reconcile planned vs executed analyses into an auditable manifest.

    Returns ``{"entries": [...], "divergences": [...], "n_divergences": int,
    "modalities": [...]}``. A divergence is an analysis the plan called for that
    did not run or was skipped — the signal that would have caught the PBMC thin
    report before a rerun.
    """
    phrases = _plan_phrases(exp_ctx)
    findings = _scrna_findings(agent_results)
    has_scrna = bool(findings) or "scrna_agent" in (agent_results or {})

    entries: list[dict] = []
    if has_scrna:
        for spec in _SCRNA_ANALYSES:
            planned = _is_planned(spec, phrases)
            val = None
            for fk in spec["finding_keys"]:
                if fk in findings:
                    val = findings[fk]
                    break
            status, reason = _status_from_finding(val)
            ran = status == "ran"
            divergence = bool(planned and not ran)
            entries.append({
                "modality": "scRNA",
                "analysis": spec["key"],
                "label": spec["label"],
                "planned": planned,
                "status": status,
                "reason": reason,
                "divergence": divergence,
            })

    divergences = [e for e in entries if e["divergence"]]
    return {
        "entries": entries,
        "divergences": divergences,
        "n_divergences": len(divergences),
        "modalities": ["scRNA"] if has_scrna else [],
    }
