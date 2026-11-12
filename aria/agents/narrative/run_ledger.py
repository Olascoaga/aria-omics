"""Planned-vs-Run ledger (P-LEDGER, audit 2026-05-29) + claim linkage (W-LEDGER).

A deterministic per-run manifest that reconciles what the analysis plan called
for against what execution actually produced. Its job is to make any
plan -> execution divergence visible — specifically the dispatch-integrity gap
that produced the PBMC thin report, where pseudobulk DE was planned and approved
but silently skipped by a free-text keyword gate.

W-LEDGER (pre-4.6 Tier-A) closes the loop: every ledger entry gets a stable
``node_id`` and every report claim is linked to its node, so a claim is traceable
to both an evidence card (W-CLAIM) and a run-ledger node. C3 additionally requires
both nodes to be typed and provenance-bearing and fails public eligibility when a
node is absent or its analysis did not run.

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
    {"key": "marker_discovery", "label": "Per-cluster marker discovery",
     "plan_kw": ["marker", "wilcoxon", "cluster-vs-rest"],
     "finding_keys": ["differential_expression"]},
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
    {"key": "scatac_preprocessing",
     "label": "scATAC raw preprocessing & matrix construction",
     "plan_kw": ["chromap", "fastq to fragments", "fragments to matrix",
                 "barcode whitelist"],
     "finding_keys": ["fastq_to_fragments", "fragments_to_matrix",
                      "cohort_matrix_union"]},
    # ADR-077: raw bulk ATAC FASTQ -> filtered BAM front-end (fastp + bwa-mem2 +
    # ATAC filtering). Reconciled from the bulk lane's ``preprocessing`` finding
    # by ``_reconcile_chromatin`` (the scATAC node above is untouched).
    {"key": "bulk_atac_preprocessing",
     "label": "Raw bulk ATAC alignment and filtering",
     "plan_kw": ["bwa-mem2", "bwa mem", "bulk atac alignment",
                 "atac alignment", "atac_align"],
     "finding_keys": ["bulk_atac_preprocessing"]},
    {"key": "dimensionality_reduction",
     "label": "Dimensionality reduction (LSI/TF-IDF)",
     "plan_kw": ["lsi", "latent semantic", "tf-idf", "tfidf",
                 "dimensionality", "svd"],
     "finding_keys": ["lsi_params", "lsi", "lsi_clustering"]},
    {"key": "peak_calling", "label": "Peak calling",
     "plan_kw": ["peak", "macs"], "finding_keys": ["peaks"]},
    {"key": "peak_count_matrix", "label": "Peak-count matrix",
     "plan_kw": ["peak count", "count matrix", "peak-by-sample",
                 "peak x sample"],
     "finding_keys": ["peak_counts"]},
    {"key": "differential_accessibility", "label": "Differential accessibility",
     "plan_kw": ["differential accessibility", "differentially accessible",
                 "differential peak"],
     "finding_keys": ["differential_accessibility"]},
    {"key": "motif_enrichment", "label": "TF motif enrichment",
     "plan_kw": ["motif enrichment", "transcription factor enrichment",
                 "tf enrichment"],
     "finding_keys": ["motifs"]},
    {"key": "regulatory_layers", "label": "Regulatory layers",
     "plan_kw": ["chromvar", "footprint", "peak-to-gene", "peak to gene",
                 "gene activity", "gene score", "label transfer"],
     "finding_keys": ["regulatory"]},
    {"key": "matrix_pipeline", "label": "scATAC matrix pipeline",
     "plan_kw": ["matrix pipeline", "fragment matrix", "peak matrix"],
     "finding_keys": ["matrix_pipeline"]},
    {"key": "peak_annotation", "label": "Genomic peak annotation",
     "plan_kw": ["peak annotation", "nearest tss", "genomic annotation"],
     "finding_keys": ["peak_annotation"]},
    {"key": "peak_ora", "label": "DA-peak functional enrichment",
     "plan_kw": ["peak ora", "peak pathway", "functional enrichment"],
     "finding_keys": ["peak_ora"]},
    {"key": "differential_tf_footprinting",
     "label": "Differential TF footprinting",
     "plan_kw": ["footprint", "tobias", "bindetect"],
     "finding_keys": ["footprinting"]},
]


# Bulk RNA-seq analyses. Bulk pathway/GSEA results live nested under each
# contrast (not a clean top-level finding), so the bulk findings are normalized
# (see _bulk_findings_normalized) before reconciliation. Process vocabulary only
# (no biology — ADR-011).
_BULK_ANALYSES: list[dict[str, Any]] = [
    {"key": "preprocessing", "label": "Raw FASTQ preprocessing & alignment",
     "plan_kw": ["fastq", "fastp", "star align", "starsolo", "featurecounts",
                 "umi dedup"],
     "finding_keys": ["preprocessing"]},
    {"key": "qc", "label": "Sample QC",
     "plan_kw": ["qc", "quality control"], "finding_keys": ["sample_qc"]},
    {"key": "differential_expression", "label": "Differential expression",
     "plan_kw": ["deseq", "differential expression", "differentially expressed"],
     "finding_keys": ["contrasts"]},
    {"key": "pathway_enrichment", "label": "Pathway enrichment (ORA/GSEA)",
     "plan_kw": ["pathway", "ora", "enrichment", "gsea"],
     "finding_keys": ["pathway_enrichment"]},
]


# V48 Slice I7 (ADR-070): paired multiome integration analyses reconciled
# directly from IntegrationAgent structured outputs. Process vocabulary only.
_INTEGRATION_ANALYSES: list[dict[str, Any]] = [
    {"key": "paired_qc", "label": "Paired multiome QC",
     "plan_kw": ["paired qc", "multiome qc", "shared barcode"],
     "finding_keys": ["paired_qc"]},
    {"key": "wnn", "label": "Weighted Nearest Neighbors (WNN)",
     "plan_kw": ["wnn", "weighted nearest", "joint manifold"],
     "finding_keys": ["wnn"]},
    {"key": "mofa", "label": "MOFA+ factor decomposition",
     "plan_kw": ["mofa", "factor analysis", "variance decomposition"],
     "finding_keys": ["mofa"]},
    {"key": "peak_to_gene", "label": "Cis peak-to-gene linkage",
     "plan_kw": ["peak2gene", "peak-to-gene", "peak to gene", "cis link"],
     "finding_keys": ["peak_to_gene"]},
    {"key": "eregulons", "label": "Candidate eRegulon inference",
     "plan_kw": ["eregulon", "scenic", "cistrome", "regulatory triplet"],
     "finding_keys": ["eregulons"]},
    {"key": "cross_modal_triangulation",
     "label": "Cross-modal evidence triangulation",
     "plan_kw": ["triangulation", "cross-modal evidence"],
     "finding_keys": ["cross_modal_triangulation"]},
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
        st = str(val.get("status") or "").lower().strip()
        if st == "not_run":
            return "not_run", val.get("reason")
        if st in ("skipped", "blocked"):
            return "skipped", (val.get("reason") or val.get("details"))
        if st in ("error", "failed"):
            return "error", (
                val.get("error_type") or val.get("details") or val.get("reason")
            )
        if st in ("success", "ok", "done"):
            return "ran", None
        if val.get("ran") is False:
            return "skipped", val.get("reason")
        if val.get("ran") is True:
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
    if not isinstance(ch, dict):
        return {}
    # ChromatinAgent nests analysis findings under a per-modality wrapper
    # (findings -> scATAC -> findings -> {qc, lsi, ...}); flatten to the analysis
    # keys this reconciler expects (mirrors the narrator's unwrap).
    from aria.agents.narrative.narrators.chromatin import (
        unwrap_chromatin_findings,
    )
    return unwrap_chromatin_findings(ch) or {}


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
    if not out.get("sample_qc"):
        preprocessing = findings.get("preprocessing") or {}
        if isinstance(preprocessing, dict) and preprocessing.get("qc"):
            out["sample_qc"] = preprocessing["qc"]
    contrasts = findings.get("contrasts") or []
    pathway_ran = False
    for c in contrasts:
        if not isinstance(c, dict):
            continue
        if c.get("pathways"):
            pathway_ran = True
            break
        # Mirror the bulk narrator's GSEA block-creation condition exactly
        # (narrators/bulk_rna.py: gsea_table OR running_sums OR top_table_fig),
        # so a GSEA-only contrast does not read as a not-run pathway node.
        plots = c.get("plots") or {}
        if isinstance(plots, dict) and (
            plots.get("gsea_table")
            or plots.get("gsea_running_sums")
            or plots.get("gsea_top_table")
        ):
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
            "node_type": "analysis_run",
            "provenance": {
                "kind": "structured_agent_result",
                "modality": modality,
                "finding_keys": list(spec["finding_keys"]),
            },
            "planned": planned,
            "status": status,
            "reason": reason,
            "divergence": bool(planned and not ran),
        })
    return entries


def _is_bulk_atac_preprocessing(preprocessing: Any) -> bool:
    return isinstance(preprocessing, dict) and (
        preprocessing.get("data_type") == "bulk_ATAC"
        or "bam_files" in preprocessing
    )


def _reconcile_chromatin(findings: dict, phrases: list[str]) -> list[dict]:
    """Reconcile chromatin findings, including the bulk ATAC raw front-end.

    The bulk lane stores its FASTQ -> BAM step under ``preprocessing``; it maps
    to ``ledger://chromatin/bulk_atac_preprocessing`` with script/stack/input/
    artifact provenance. ``executed`` requires a successful status AND at least
    one produced BAM — never inferred from the request alone.
    """
    findings = dict(findings or {})
    preprocessing = findings.get("preprocessing")
    bulk_pre = _is_bulk_atac_preprocessing(preprocessing)
    if bulk_pre:
        findings["bulk_atac_preprocessing"] = preprocessing
    entries = _entries_for("chromatin", _CHROMATIN_ANALYSES, findings, phrases)
    if not bulk_pre:
        return entries
    entry = next(e for e in entries
                 if e["analysis"] == "bulk_atac_preprocessing")
    status = str(preprocessing.get("status") or "").lower().strip()
    artifacts: list[str] = []
    for rec in preprocessing.get("bam_files") or []:
        if isinstance(rec, dict):
            if rec.get("status", "success") == "success" and rec.get("bam"):
                artifacts.append(str(rec["bam"]))
        elif rec:
            artifacts.append(str(rec))
    inputs = [str(f) for f in (preprocessing.get("fastq_files") or []) if f]
    entry["provenance"].update({
        "script": "aria/scripts/atac_align.py",
        "environment_stack": "atacseq",
        "inputs": inputs,
        "artifacts": artifacts,
        "validation_level": preprocessing.get("validation_level",
                                              "production"),
    })
    if status in ("success", "done") and artifacts:
        entry["status"] = "executed"
        entry["reason"] = None
        entry["divergence"] = False
    elif entry["status"] == "ran":
        # A success stamp without a produced BAM is not an executed alignment.
        entry["status"] = "skipped"
        entry["reason"] = "no filtered BAM produced"
        entry["divergence"] = bool(entry["planned"])
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
        entries += _reconcile_chromatin(ch_findings, phrases)
        modalities.append("chromatin")

    bk_findings = _bulk_findings_normalized(agent_results)
    if bool(bk_findings) or "bulk_rna_agent" in (agent_results or {}):
        entries += _entries_for("bulk", _BULK_ANALYSES, bk_findings, phrases)
        modalities.append("bulk")

    int_findings = _integration_findings_normalized(agent_results, exp_ctx)
    if bool(int_findings) or "integration_agent" in (agent_results or {}):
        entries += _entries_for(
            "integration", _INTEGRATION_ANALYSES, int_findings, phrases
        )
        modalities.append("integration")

    divergences = [e for e in entries if e["divergence"]]
    return {
        "entries": entries,
        "divergences": divergences,
        "n_divergences": len(divergences),
        "modalities": modalities,
    }


def _integration_findings_normalized(
    agent_results: dict,
    exp_ctx: dict | None = None,
) -> dict:
    """Normalize IntegrationAgent outputs for ledger reconciliation.

    Synthesizes ``cross_modal_triangulation`` as executed iff the deterministic
    evidence graph genuinely ran over structured outputs (mirrors
    ``IntegrationNarrator.collect``; never fabricates on absent inputs).
    """
    integ = (agent_results or {}).get("integration_agent")
    if not isinstance(integ, dict):
        return {}
    from aria.agents.narrative.narrators.integration import (
        unwrap_integration_findings,
    )
    findings = unwrap_integration_findings(integ)
    if not isinstance(findings, dict):
        return {}
    out = dict(findings)
    try:
        from aria.agents.narrative.synthesis.cross_modal_triangulation import (
            build_cross_modal_evidence_graph,
        )
        graph = build_cross_modal_evidence_graph(agent_results, exp_ctx or {})
        if graph.get("ran"):
            out["cross_modal_triangulation"] = {"status": "success"}
    except Exception:
        pass
    return out


# ── W-LEDGER: claim ↔ ledger-node linkage + active verification ──────────────

# Audited claim tiers retained as the ADR-057 quarantine boundary. C3 linkage
# now applies to every non-empty result claim, not only to this tier subset.
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
    "power": "differential_expression",
    "peak2gene": "regulatory_layers",
    "integrated_signal": "biological_synthesis",
    "convergent_evidence": "biological_synthesis",
    "divergent_evidence": "biological_synthesis",
    "integration_limitations": "biological_synthesis",
    "scrna_integrated_signal": "biological_synthesis",
    "scrna_context_layers": "biological_synthesis",
    "scrna_integration_limitations": "biological_synthesis",
}

_REPORT_DERIVED_ANALYSES = {
    "executive_summary",
    "integrated_signal",
    "convergent_evidence",
    "divergent_evidence",
    "integration_limitations",
    "scrna_integrated_signal",
    "scrna_context_layers",
    "scrna_integration_limitations",
}


def _ledger_modality_group(modality: str) -> str:
    """Map a block/claim modality label to a ledger modality group."""
    m = str(modality or "").lower().strip()
    if "report" in m or "integrated" in m:
        return "report"
    if m == "integration" or "multiome" in m:
        return "integration"
    if "bulk" in m:
        return "bulk"
    if "atac" in m or "chromatin" in m:
        return "chromatin"
    return "scRNA"


def node_id_for(modality: str, analysis: str) -> str:
    """Deterministic ledger node id for a claim's (modality, analysis)."""
    group = _ledger_modality_group(modality)
    raw = str(analysis or "")
    if group == "integration":
        # Modality-aware mapping: chromatin's single-assay `peak2gene` maps to
        # `regulatory_layers`, whereas paired integration owns `peak_to_gene`.
        int_aliases = {
            "peak2gene": "peak_to_gene",
            "eregulon": "eregulons",
            "triangulation": "cross_modal_triangulation",
        }
        key = int_aliases.get(raw, raw)
    else:
        key = _LEDGER_KEY_ALIASES.get(raw, raw)
    return f"ledger://{group}/{key}"


def _node_index(run_ledger: dict | None) -> dict[str, dict]:
    return {
        e["node_id"]: e
        for e in (run_ledger or {}).get("entries", []) or []
        if isinstance(e, dict) and e.get("node_id")
    }


def ensure_report_ledger_nodes(run_ledger: dict, blocks: list) -> None:
    """Add canonical nodes for report-derived analyses that actually ran.

    Primary assay nodes must come from :func:`build_run_ledger`; this allowlist
    only covers deterministic report products whose execution is represented by
    a typed ``NarrativeBlock`` rather than an agent finding. Unknown analyses are
    deliberately not synthesized so public compilation can fail closed.
    """
    if not isinstance(run_ledger, dict):
        return
    entries = run_ledger.setdefault("entries", [])
    index = _node_index(run_ledger)
    by_node: dict[str, list] = {}
    for block in blocks or []:
        analysis = str(getattr(block, "analysis", "") or "")
        if analysis not in _REPORT_DERIVED_ANALYSES:
            continue
        nid = node_id_for(getattr(block, "modality", None), analysis)
        by_node.setdefault(nid, []).append(block)
    for nid, node_blocks in by_node.items():
        if nid in index:
            continue
        canonical_analysis = nid.rsplit("/", 1)[-1]
        entries.append({
            "modality": "report",
            "analysis": canonical_analysis,
            "label": (
                "Executive summary" if canonical_analysis == "executive_summary"
                else "Integrated biological synthesis"
            ),
            "node_id": nid,
            "node_type": "derived_report_analysis",
            "provenance": {
                "kind": "compiled_narrative_block",
                "block_ids": sorted({
                    str(getattr(block, "id", "")) for block in node_blocks
                }),
            },
            "planned": True,
            "status": "ran",
            "reason": "compiled from governed narrative blocks",
            "divergence": False,
        })
    if by_node:
        modalities = run_ledger.setdefault("modalities", [])
        if "report" not in modalities:
            modalities.append("report")
    run_ledger["divergences"] = [
        entry for entry in entries
        if isinstance(entry, dict) and entry.get("divergence")
    ]
    run_ledger["n_divergences"] = len(run_ledger["divergences"])


def link_claims_to_ledger(claims: list[dict],
                          run_ledger: dict | None) -> dict[str, Any]:
    """Link each claim manifest to its ledger node, in place.

    Adds ``ledger_node_id``, ``ledger_status``, and ``ledger_linked`` to every
    claim dict. C3 makes missing, untyped or provenance-free nodes violations,
    as well as claims linked to analyses that did not run. The caller can then
    withhold those manifests before they reach a public surface.
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
            claim["ledger_node_type"] = None
            claim["ledger_provenance"] = None
            claim["ledger_status"] = "no_ledger_node"
            claim["ledger_linked"] = False
            unlinked += 1
            violations.append({
                "claim_id": claim.get("claim_id"),
                "node_id": nid,
                "ledger_status": "no_ledger_node",
                "reason": "no ledger node exists for the public claim",
                "tier": claim.get("tier"),
            })
            continue
        status = node.get("status")
        claim["ledger_node_id"] = nid
        claim["ledger_node_type"] = node.get("node_type")
        claim["ledger_provenance"] = node.get("provenance")
        claim["ledger_status"] = status
        claim["ledger_linked"] = True
        linked += 1
        reason = None
        if not claim.get("evidence_card_id"):
            reason = "public claim has no evidence-card id"
        elif not node.get("node_type"):
            reason = "ledger node has no node_type"
        elif not node.get("provenance"):
            reason = "ledger node has no provenance"
        elif status in _NOT_RUN_STATUSES:
            reason = "public claim cites an analysis that did not run"
        if reason:
            violations.append({
                "claim_id": claim.get("claim_id"),
                "node_id": nid,
                "ledger_status": status,
                "reason": reason,
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

    Every non-empty public result claim must resolve to a typed, provenance-bearing
    ledger node and an evidence-card id. It may not cite an analysis marked
    not-run/skipped/error. Diagnostic blocks with no claim remain visible and are
    outside this result-claim invariant.
    """
    index = _node_index(run_ledger)
    violations: list[dict] = []
    for block in blocks or []:
        if not str(getattr(block, "claim", "") or "").strip():
            continue
        metadata = getattr(block, "metadata", {}) or {}
        tier = (metadata.get("claim", {}) or {}).get("tier")
        nid = node_id_for(getattr(block, "modality", None),
                          getattr(block, "analysis", None))
        node = index.get(nid)
        if node is None:
            violations.append({
                "claim_id": getattr(block, "id", None),
                "node_id": nid,
                "ledger_status": "no_ledger_node",
                "reason": "no ledger node exists for the public claim",
                "tier": tier,
            })
            continue
        status = node.get("status")
        verification = metadata.get("claim_verification") or {}
        card = verification.get("evidence_card") or {}
        reason = None
        if not card.get("evidence_card_id"):
            reason = "public claim has no evidence-card id"
        elif not card.get("node_type"):
            reason = "evidence node has no node_type"
        elif not card.get("provenance"):
            reason = "evidence node has no provenance"
        elif not node.get("node_type"):
            reason = "ledger node has no node_type"
        elif not node.get("provenance"):
            reason = "ledger node has no provenance"
        elif status in _NOT_RUN_STATUSES:
            reason = "public claim cites an analysis that did not run"
        if reason:
            violations.append({
                "claim_id": getattr(block, "id", None),
                "node_id": nid,
                "ledger_status": status,
                "reason": reason,
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
        if v["ledger_status"] == "no_ledger_node":
            raise LedgerLinkageError(
                f"{v['claim_id']}: no ledger node exists for public claim "
                f"{v['node_id']}{detail}"
            )
        raise LedgerLinkageError(
            f"{v['claim_id']}: claim of tier {v['tier']!r} cites ledger node "
            f"{v['node_id']} whose run status is {v['ledger_status']!r}{detail}"
        )
    return manifest
