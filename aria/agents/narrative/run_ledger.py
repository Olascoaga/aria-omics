"""Planned-vs-Run ledger (P-LEDGER, audit 2026-05-29) + claim linkage (W-LEDGER).

A deterministic per-run manifest that reconciles what the analysis plan called
for against what execution actually produced. Its job is to make any
plan -> execution divergence visible — specifically the dispatch-integrity gap
that produced the PBMC thin report, where pseudobulk DE was planned and approved
but silently skipped by a free-text keyword gate.

W-LEDGER (pre-4.6 Tier-A) closes the loop: every ledger entry gets a stable
``node_id`` and every report claim is linked to its node, so a claim is traceable
to both an evidence card (W-CLAIM) and a run-ledger node. The linkage is actively
verified — a claim of associative-or-stronger tier may not cite a node the run
marked not-run/skipped/error (the thin-report contradiction).

This is pure bookkeeping over the structured plan (DesignIntelligence) and the
structured agent results; it asserts no biology and runs no LLM.
"""

from __future__ import annotations

from typing import Any


class LedgerLinkageError(ValueError):
    """A report claim cites a ledger analysis the run did not execute."""


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


# Bulk RNA-seq analyses. Bulk pathway/GSEA results live nested under each
# contrast (not a clean top-level finding), so the bulk findings are normalized
# (see _bulk_findings_normalized) before reconciliation. Process vocabulary only
# (no biology — ADR-011).
_BULK_ANALYSES: list[dict[str, Any]] = [
    {"key": "qc", "label": "Sample QC",
     "plan_kw": ["qc", "quality control"], "finding_keys": ["sample_qc"]},
    {"key": "differential_expression", "label": "Differential expression",
     "plan_kw": ["deseq", "differential expression", "differentially expressed"],
     "finding_keys": ["contrasts"]},
    {"key": "pathway_enrichment", "label": "Pathway enrichment (ORA/GSEA)",
     "plan_kw": ["pathway", "ora", "enrichment", "gsea"],
     "finding_keys": ["pathway_enrichment"]},
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


def _bulk_findings_normalized(agent_results: dict) -> dict:
    """Bulk findings reshaped so the generic reconciler can read them.

    Bulk pathway/GSEA results are nested per contrast, so synthesize a single
    ``pathway_enrichment`` value that is truthy iff at least one contrast carries
    ORA terms or a GSEA table — exactly how the bulk narrator decides a pathway
    block exists. Never fabricates: returns ``{}`` for absent inputs.
    """
    bk = (agent_results or {}).get("bulk_rna_agent", {})
    findings = bk.get("findings", bk) if isinstance(bk, dict) else {}
    if not isinstance(findings, dict):
        return {}
    out = dict(findings)
    contrasts = findings.get("contrasts") or []
    pathway_ran = False
    for c in contrasts:
        if not isinstance(c, dict):
            continue
        if c.get("pathways"):
            pathway_ran = True
            break
        plots = c.get("plots") or {}
        if isinstance(plots, dict) and plots.get("gsea_table"):
            pathway_ran = True
            break
    # Only assert the node when there is evidence it ran; otherwise leave it
    # absent so reconciliation reports honest not_run (no fabrication).
    if pathway_ran:
        out["pathway_enrichment"] = {"status": "success"}
    return out


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
            "node_id": f"ledger://{modality}/{spec['key']}",
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

    bk_findings = _bulk_findings_normalized(agent_results)
    if bool(bk_findings) or "bulk_rna_agent" in (agent_results or {}):
        entries += _entries_for("bulk", _BULK_ANALYSES, bk_findings, phrases)
        modalities.append("bulk")

    divergences = [e for e in entries if e["divergence"]]
    return {
        "entries": entries,
        "divergences": divergences,
        "n_divergences": len(divergences),
        "modalities": modalities,
    }


# ── W-LEDGER: claim ↔ ledger-node linkage + active verification ──────────────

# Claim tiers (Claim Compiler) at or above "associative" assert a real result;
# a descriptive claim only restates a measured quantity and is never a linkage
# violation even if its analysis did not run.
_NON_DESCRIPTIVE_TIERS = {
    "associative", "weak_mechanistic", "strong_mechanistic",
    "causal_experimental",
}
# Statuses that mean the analysis did not actually produce results.
_NOT_RUN_STATUSES = {"not_run", "skipped", "error"}

# A narrative block's ``analysis`` label can be a narrator-specific variant of a
# canonical ledger analysis key. Map the variants so the claim resolves to the
# same node the ledger reconciled (modality disambiguates: bulk vs scRNA).
_LEDGER_KEY_ALIASES = {
    "sample_qc": "qc",
    "data_quality": "qc",
    "pseudobulk_pathways": "pathway_enrichment",
    "gsea_preranked": "pathway_enrichment",
    "contrast": "differential_expression",
}


def _ledger_modality_group(modality: str) -> str:
    """Map a block/claim modality label to a ledger modality group."""
    m = str(modality or "").lower()
    if "bulk" in m:
        return "bulk"
    if "atac" in m or "chromatin" in m:
        return "chromatin"
    return "scRNA"


def node_id_for(modality: str, analysis: str) -> str:
    """Deterministic ledger node id for a claim's (modality, analysis)."""
    group = _ledger_modality_group(modality)
    key = _LEDGER_KEY_ALIASES.get(str(analysis or ""), str(analysis or ""))
    return f"ledger://{group}/{key}"


def _node_index(run_ledger: dict | None) -> dict[str, dict]:
    return {
        e["node_id"]: e
        for e in (run_ledger or {}).get("entries", []) or []
        if isinstance(e, dict) and e.get("node_id")
    }


def link_claims_to_ledger(claims: list[dict],
                          run_ledger: dict | None) -> dict[str, Any]:
    """Link each claim manifest to its ledger node, in place.

    Adds ``ledger_node_id``, ``ledger_status``, and ``ledger_linked`` to every
    claim dict. A claim whose analysis has no matching node records
    ``ledger_status="no_ledger_node"`` honestly (coverage gap, not fabrication).
    A non-descriptive claim linked to a not-run/skipped/error node is a
    violation. Returns a serializable linkage summary for methodology.json.
    """
    index = _node_index(run_ledger)
    linked = 0
    unlinked = 0
    violations: list[dict] = []
    for claim in claims or []:
        if not isinstance(claim, dict):
            continue
        nid = node_id_for(claim.get("modality"), claim.get("analysis"))
        node = index.get(nid)
        if node is None:
            claim["ledger_node_id"] = None
            claim["ledger_status"] = "no_ledger_node"
            claim["ledger_linked"] = False
            unlinked += 1
            continue
        status = node.get("status")
        claim["ledger_node_id"] = nid
        claim["ledger_status"] = status
        claim["ledger_linked"] = True
        linked += 1
        if claim.get("tier") in _NON_DESCRIPTIVE_TIERS and status in _NOT_RUN_STATUSES:
            violations.append({
                "claim_id": claim.get("claim_id"),
                "node_id": nid,
                "ledger_status": status,
                "reason": node.get("reason"),
                "tier": claim.get("tier"),
            })
    return {
        "linked": linked,
        "unlinked": unlinked,
        "violations": violations,
        "n_violations": len(violations),
    }


def verify_blocks_against_ledger(blocks: list,
                                 run_ledger: dict | None,
                                 *, strict: bool = True) -> dict[str, Any]:
    """Actively verify rendered blocks against the run ledger.

    A block whose claim tier is associative-or-stronger may not cite a ledger
    node the run marked not-run/skipped/error. When ``strict`` is true the first
    such contradiction raises :class:`LedgerLinkageError` so the report cannot be
    rendered with a claim about an analysis that did not run. Blocks without a
    matching node (coverage gap) and descriptive claims are never violations.
    """
    index = _node_index(run_ledger)
    violations: list[dict] = []
    for block in blocks or []:
        tier = (getattr(block, "metadata", {}) or {}).get("claim", {}) or {}
        tier = tier.get("tier")
        if tier not in _NON_DESCRIPTIVE_TIERS:
            continue
        nid = node_id_for(getattr(block, "modality", None),
                          getattr(block, "analysis", None))
        node = index.get(nid)
        if node is None:
            continue
        status = node.get("status")
        if status in _NOT_RUN_STATUSES:
            violations.append({
                "claim_id": getattr(block, "id", None),
                "node_id": nid,
                "ledger_status": status,
                "reason": node.get("reason"),
                "tier": tier,
            })
    manifest = {
        "status": "ok" if not violations else "violated",
        "n_violations": len(violations),
        "violations": violations,
    }
    if strict and violations:
        v = violations[0]
        detail = f" ({v['reason']})" if v.get("reason") else ""
        raise LedgerLinkageError(
            f"{v['claim_id']}: claim of tier {v['tier']!r} cites ledger node "
            f"{v['node_id']} whose run status is {v['ledger_status']!r}{detail}"
        )
    return manifest
