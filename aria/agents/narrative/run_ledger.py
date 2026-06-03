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


# P3-1 (pre-4.6 polish): the same planned-vs-run pattern for chromatin, wired in
# from day one so the moment the v4.6 scATAC stack lands a thin chromatin report
# (QC ran but LSI/peaks/motifs did not) is reconciled as a divergence rather than
# vanishing silently. Technical/process vocabulary only (no biology — ADR-011).
# finding_keys mirror ChromatinAgent's structured `findings` keys.
_CHROMATIN_ANALYSES: list[dict[str, Any]] = [
    {"key": "qc", "label": "Quality control",
     "plan_kw": ["qc", "quality control"], "finding_keys": ["qc"]},
    {"key": "dimensionality_reduction",
     "label": "Dimensionality reduction (LSI/TF-IDF)",
     "plan_kw": ["lsi", "latent semantic", "tf-idf", "tfidf",
                 "dimensionality", "svd"],
     "finding_keys": ["lsi_params", "lsi"]},
    {"key": "peak_calling", "label": "Peak calling",
     "plan_kw": ["peak", "macs"], "finding_keys": ["peaks"]},
    {"key": "differential_accessibility", "label": "Differential accessibility",
     "plan_kw": ["differential accessibility", "differentially accessible",
                 "differential peak"],
     "finding_keys": ["differential_accessibility"]},
    {"key": "motif_enrichment", "label": "TF motif enrichment",
     "plan_kw": ["motif", "transcription factor", "chromvar", "tf enrichment"],
     "finding_keys": ["motifs"]},
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


def _chromatin_findings(agent_results: dict) -> dict:
    ch = (agent_results or {}).get("chromatin_agent", {})
    if isinstance(ch, dict):
        return ch.get("findings", ch) or {}
    return {}


def _entries_for(modality: str, specs: list[dict], findings: dict,
                 phrases: list[str]) -> list[dict]:
    """Reconcile one modality's analyses against its structured findings."""
    entries: list[dict] = []
    for spec in specs:
        planned = _is_planned(spec, phrases)
        val = None
        for fk in spec["finding_keys"]:
            if fk in findings:
                val = findings[fk]
                break
        status, reason = _status_from_finding(val)
        ran = status == "ran"
        entries.append({
            "modality": modality,
            "analysis": spec["key"],
            "label": spec["label"],
            "planned": planned,
            "status": status,
            "reason": reason,
            "divergence": bool(planned and not ran),
        })
    return entries


def build_run_ledger(exp_ctx: dict, agent_results: dict) -> dict:
    """Reconcile planned vs executed analyses into an auditable manifest.

    Returns ``{"entries": [...], "divergences": [...], "n_divergences": int,
    "modalities": [...]}``. A divergence is an analysis the plan called for that
    did not run or was skipped — the signal that would have caught the PBMC thin
    report before a rerun. Covers scRNA and (P3-1) chromatin from day one.
    """
    phrases = _plan_phrases(exp_ctx)
    entries: list[dict] = []
    modalities: list[str] = []

    sc_findings = _scrna_findings(agent_results)
    if bool(sc_findings) or "scrna_agent" in (agent_results or {}):
        entries += _entries_for("scRNA", _SCRNA_ANALYSES, sc_findings, phrases)
        modalities.append("scRNA")

    ch_findings = _chromatin_findings(agent_results)
    if bool(ch_findings) or "chromatin_agent" in (agent_results or {}):
        entries += _entries_for("chromatin", _CHROMATIN_ANALYSES,
                                ch_findings, phrases)
        modalities.append("chromatin")

    divergences = [e for e in entries if e["divergence"]]
    return {
        "entries": entries,
        "divergences": divergences,
        "n_divergences": len(divergences),
        "modalities": modalities,
    }
