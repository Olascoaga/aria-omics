"""
ARIA NarrativeAgent
--------------------
Synthesizes all agent findings into a scientific report.

Receives from the Orchestrator:
  - findings:         list of Finding dicts from all agents (via MessageBus)
  - agent_results:    raw result dicts from each agent
  - exp_context:      experiment metadata (organism, genome, question)
  - biological_intent: parsed intent from OrchestratorAgent

Produces:
  - HTML report  (~/.aria/reports/{experiment_id}.html)
  - Methods section (plain text, copy-paste ready for manuscript)
  - Executive summary (1 paragraph, PI-readable)

Report structure:
  1. Executive Summary        — what was found, in plain language
  2. Data & QC                — what data, how many cells/reads passed
  3. Findings by modality     — RNA, Chromatin, 3D genome
  4. Integration findings     — cross-modal insights (if multimodal)
  5. Conflicts & Limitations  — where data is ambiguous (honest uncertainty)
  6. Methods section          — reproducible, exportable for manuscript
  7. Appendix                 — parameter decisions log from ARIAMemory

Design principles:
  - Every claim carries its confidence level (HIGH / MEDIUM / LOW)
  - LOW confidence findings are reported but explicitly flagged
  - INSUFFICIENT findings are listed as "data insufficient to conclude"
  - DebateCouncil verdicts and limitations are always included
  - Methods section uses exact parameter values from ParameterAdvisor decisions
  - No hallucinated biology — only what the data actually showed
"""

from __future__ import annotations

import json
import logging
import uuid
import html as _html
import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from aria.utils.provenance import collect_llm_usage, collect_provenance

from aria import __version__ as ARIA_VERSION
from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, MessageType
from aria.llm.provider import LLMProvider, TaskTier
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.narrative")


NARRATIVE_SYSTEM = """
You are ARIA's NarrativeAgent — a scientific writer specializing in
multi-omics research reports.

Your job is to synthesize computational findings into clear, accurate,
scientifically defensible prose.

Rules you never break:
1. Only report what the data showed. Never extrapolate beyond findings.
2. Always include confidence levels. "HIGH confidence" means multiple
   independent lines of evidence. "MEDIUM confidence" means one line
   of evidence with caveats. "LOW confidence" means preliminary signal
   requiring validation.
3. Distinguish correlation from causation. Peak-to-gene links are
   correlational, not regulatory, unless validated by orthogonal data.
4. When DebateCouncil revised a claim, use the REVISED version, not
   the original overclaim.
5. The limitations section is not optional. Every finding has limits.
6. Methods must be reproducible: include exact parameter values,
   tool versions, and decision rationales.
7. If data was insufficient to conclude, say so directly. Do not
   soften "insufficient evidence" into "preliminary findings suggest."

Tone: precise, measured, scientifically appropriate.
Audience: biology PhD / MD-PhD level.
""".strip()


class NarrativeAgent(BaseAgent):

    name        = "narrative_agent"
    description = (
        "Synthesizes all findings into a scientific HTML report "
        "with reproducible methods section."
    )

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.reports_dir = Path.home() / ".aria" / "reports"
        self.reports_dir.mkdir(parents=True, exist_ok=True)

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        self.publish_status(experiment_id,
                            "NarrativeAgent: synthesizing findings...", 0.0)

        exp_ctx       = context.get("exp_context", {})
        intent        = context.get("biological_intent", {})
        agent_results = context.get("agent_results", {})
        raw_findings  = context.get("findings", [])

        # Retrieve decisions log from memory
        decisions = self._get_decisions_log(experiment_id)

        # Group findings by confidence
        grouped = self._group_findings(raw_findings)

        # Build the report directory up front so figure subprocesses can
        # write into report_dir/figures/ before the HTML sections are
        # composed (the TUI path used to skip this — only the harness
        # called generate_figures, so reports came out at ~30 KB with
        # zero embedded PNGs).
        report_dir = self._build_report_dir(experiment_id, intent, exp_ctx)
        self._generate_scrna_figures(agent_results, report_dir)

        self.publish_status(experiment_id,
                            "Writing executive summary...", 0.2)

        # Generate each report section
        executive_summary = self._write_executive_summary(
            exp_ctx, intent, grouped, agent_results
        )

        self.publish_status(experiment_id,
                            "Writing findings sections...", 0.4)

        findings_sections = self._write_findings_sections(
            grouped, agent_results, exp_ctx
        )

        self.publish_status(experiment_id,
                            "Writing methods section...", 0.7)

        methods = self._write_methods_section(
            exp_ctx, agent_results, decisions
        )

        self.publish_status(experiment_id,
                            "Rendering HTML report...", 0.85)

        # Render HTML
        report_path = self._render_html_report(
            experiment_id=experiment_id,
            exp_ctx=exp_ctx,
            intent=intent,
            executive_summary=executive_summary,
            findings_sections=findings_sections,
            grouped_findings=grouped,
            methods=methods,
            decisions=decisions,
            agent_results=agent_results,
            report_dir=report_dir,
        )

        self.publish_status(experiment_id,
                            f"Report saved: {report_path}", 1.0)

        return {
            "status":            "done",
            "report_path":       str(report_path),
            "executive_summary": executive_summary,
            "n_findings":        len(raw_findings),
            "n_high":            len(grouped["high"]),
            "n_medium":          len(grouped["medium"]),
            "n_low":             len(grouped["low"]),
        }

    # ── Section writers ───────────────────────────────────────────────────

    def _write_executive_summary(self, exp_ctx: dict, intent: dict,
                                  grouped: dict,
                                  agent_results: dict) -> str:
        """
        One-paragraph executive summary for the PI.
        Uses deterministic structured summaries when scRNA outputs contain
        enough detail; otherwise falls back to the LLM.

        Feeds the LLM BOTH the findings bus AND concrete agent results
        (DE counts, pathways, contrast details) so it can't hallucinate
        that analyses "didn't run" when their outputs are clearly present.
        """
        deterministic = self._deterministic_executive_summary(
            exp_ctx, intent, grouped, agent_results
        )
        if deterministic:
            return deterministic

        high_findings = grouped["high"][:5]
        med_findings  = grouped["medium"][:3]
        n_low         = len(grouped["low"])
        n_insuff      = len(grouped["insufficient"])

        finding_summaries = "\n".join([
            f"- [HIGH] {f.get('summary', str(f))[:200]}"
            for f in high_findings
        ] + [
            f"- [MEDIUM] {f.get('summary', str(f))[:200]}"
            for f in med_findings
        ])

        # Build a CONCRETE results block from agent_results (anti-hallucination)
        concrete = self._summarize_agent_results_for_llm(agent_results)

        prompt = f"""
You are writing the executive summary of a bioinformatics report for a
PI scientist. Be direct and specific. No marketing language.

Experiment: {exp_ctx.get('organism', 'unknown')} / {exp_ctx.get('genome', 'unknown')}
Biological question: {intent.get('summary', exp_ctx.get('user_question', 'unknown'))}
Analysis type: {intent.get('analysis_type', 'unknown')}

═══════════════════════════════════════════════════════════════════
CONCRETE RESULTS FROM PIPELINE (use these numbers, don't invent):
═══════════════════════════════════════════════════════════════════
{concrete}
═══════════════════════════════════════════════════════════════════

Additional findings on the confidence bus:
{finding_summaries}

Note: {n_low} LOW-confidence findings, {n_insuff} insufficient analyses.

Write 3-4 sentences. STRICT requirements:
1. Lead with the biological question and whether the data answers it
2. Quote specific numbers FROM THE CONCRETE RESULTS SECTION ABOVE
3. State confidence honestly (HIGH/MEDIUM/LOW)
4. End with the most important next step or actionable limitation

CRITICAL: If CONCRETE RESULTS shows DE genes, pathways, or contrasts,
the pipeline DID run those analyses — do NOT claim they're incomplete.
Only report incompleteness if CONCRETE RESULTS says "not run" or is empty.

FORBIDDEN phrases: "comprehensive analysis", "robust evidence",
"foundational", "provides insights", "reveals important",
"this study", "elucidate", "leverage". Use plain scientific English.

Example good tone: "Bulk RNA-seq comparing treatment vs control
(3v3 replicates) identified 2,232 DE genes with medium confidence.
Top genes and pathway enrichment suggest a coherent state shift, but
the result should be treated as differential-expression evidence rather
than proof of direct regulation. Sample size (n=3/group) limits power
for small-effect genes."
"""
        try:
            return self.llm.complete(
                prompt=prompt,
                system=NARRATIVE_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=400,
            )
        except Exception as e:
            log.warning(f"Executive summary LLM failed: {e}")
            return self._fallback_executive_summary(grouped, intent)

    def _deterministic_executive_summary(self, exp_ctx: dict,
                                         intent: dict,
                                         grouped: dict,
                                         agent_results: dict) -> str:
        """Build a grounded summary directly from structured outputs."""
        sc = (agent_results or {}).get("scrna_agent") or \
             (agent_results or {}).get("rna_agent") or {}
        if sc.get("status") != "done":
            return ""

        from aria.agents import _narrative_scrna
        sc_f = _narrative_scrna.unwrap_scrna_findings(sc)
        qc = sc_f.get("qc", {}) or {}
        clu = sc_f.get("clustering", {}) or {}
        pb = sc_f.get("pseudobulk_de", {}) or {}
        pwp = sc_f.get("pseudobulk_pathways", {}) or {}
        ccc = sc_f.get("cell_communication", {}) or {}
        traj = sc_f.get("trajectory", {}) or {}
        if not any((pb, pwp, ccc, traj)):
            return ""

        question = intent.get("summary", exp_ctx.get("user_question", ""))
        n_cells = qc.get("n_cells_after")
        n_before = qc.get("n_cells_before")
        n_clusters = clu.get("n_clusters")
        groupby = clu.get("groupby") or (sc_f.get("clustering_decision") or {}).get("groupby")
        predef_clusters = bool(
            clu.get("predef_clusters")
            or (sc_f.get("clustering_decision") or {}).get("predef_clusters")
        )
        qc_clause = ""
        if n_cells and n_before:
            qc_clause = (
                f"after QC {n_cells:,}/{n_before:,} cells were retained"
            )
        elif n_cells:
            qc_clause = f"after QC {n_cells:,} cells were retained"
        if n_clusters:
            cluster_label = (
                f"{n_clusters} {_narrative_scrna._group_label(groupby, n_clusters)} "
                f"from obs['{groupby}']"
                if predef_clusters and groupby else
                f"{n_clusters} Leiden clusters"
            )
            qc_clause = (
                f"{qc_clause}, with {cluster_label}"
                if qc_clause else cluster_label
            )

        parts = []
        if question and qc_clause:
            parts.append(
                f"ARIA addressed the question '{question}' with "
                f"MEDIUM-confidence scRNA evidence: {qc_clause}."
            )
        elif qc_clause:
            parts.append(
                f"ARIA produced MEDIUM-confidence scRNA results: {qc_clause}."
            )

        de = sc_f.get("differential_expression") or {}
        de_status = de.get("status")
        if de_status and de_status != "success":
            parts.append(
                f"Per-cluster marker discovery did not complete "
                f"({de.get('error_type', 'Error')}); pseudobulk DE "
                f"below is the primary differential-expression signal."
            )

        if pb:
            per_group = pb.get("per_group", {}) or {}
            n_success = sum(
                1 for g in per_group.values()
                for c in (g.get("per_comparison", {}) or {}).values()
                if c.get("status") == "success"
            )
            n_with_de = sum(
                1 for g in per_group.values()
                for c in (g.get("per_comparison", {}) or {}).values()
                if c.get("status") == "success"
                and c.get("n_significant", 0) > 0
            )
            top_blocks = _narrative_scrna._top_de_blocks(pb, limit=2)
            top_txt = ""
            if top_blocks:
                top_txt = "; strongest blocks: " + "; ".join(
                    f"{g} {cmp} ({c.get('n_significant', 0):,} DE genes)"
                    for g, cmp, c in top_blocks
                )
            parts.append(
                f"Pseudobulk DE did run across {pb.get('n_groups', 0)} "
                f"{_narrative_scrna._group_label(pb.get('groupby'), pb.get('n_groups', 0))}: "
                f"{n_success} "
                f"group x comparison blocks were analyzable and "
                f"{n_with_de} had significant DE{top_txt}."
            )

        if pwp.get("per_cluster"):
            n_blocks = len(pwp["per_cluster"])
            n_sig = sum(
                1 for b in pwp["per_cluster"].values()
                if b.get("n_significant", 0) > 0
            )
            parts.append(
                f"Pathway ORA supported {n_sig}/{n_blocks} DE blocks; "
                f"LIANA found {ccc.get('n_interactions', 0)} non-autocrine "
                f"L-R interactions across {ccc.get('n_cell_types', 0)} "
                f"cell types."
            )
        elif ccc.get("status") in ("done", "success"):
            parts.append(
                f"LIANA found {ccc.get('n_interactions', 0)} non-autocrine "
                f"L-R interactions across {ccc.get('n_cell_types', 0)} "
                f"cell types."
            )

        if traj.get("status") in ("done", "success"):
            paga = traj.get("paga", {}) or {}
            pt = traj.get("pseudotime", {}) or {}
            parts.append(
                f"Trajectory analysis evaluated "
                f"{paga.get('n_connections', 0)} PAGA cluster pairs with "
                f"{paga.get('n_strong', 0)} edges above threshold; "
                f"DPT pseudotime computed={bool(pt.get('computed'))}. "
                f"The actionable next step is biological curation of the "
                f"largest between-condition DE/pathway blocks and validation of "
                f"the top communication pairs."
            )

        return " ".join(parts[:4])

    def _summarize_agent_results_for_llm(self, agent_results: dict) -> str:
        """
        Build a concise, concrete summary of agent outputs for the LLM.

        This is the anti-hallucination substrate — the LLM sees raw numbers
        from the actual pipeline and cannot claim an analysis is missing
        when its outputs are present in this block.
        """
        if not agent_results:
            return "(no agent results available)"

        lines = []

        # Bulk RNA
        bulk = agent_results.get("bulk_rna_agent", {})
        if bulk.get("status") == "done":
            findings = bulk.get("findings", {}) or {}
            contrasts = findings.get("contrasts", []) or []
            successful = [c for c in contrasts if c.get("status") == "success"]
            if successful:
                lines.append(
                    f"BULK RNA-seq: {len(successful)} contrast(s) ran "
                    f"successfully (DESeq2 completed, not just QC/alignment)."
                )
                for c in successful:
                    name = c.get("name", "?")
                    n_sig = c.get("n_significant", 0)
                    n_up  = c.get("n_upregulated", 0)
                    n_dn  = c.get("n_downregulated", 0)
                    tops  = c.get("top_genes", [])[:6]
                    top_str = ", ".join(
                        f"{g.get('symbol') or g.get('gene')}"
                        f"({'↑' if g.get('log2fc',0)>0 else '↓'}"
                        f"{abs(g.get('log2fc',0)):.1f})"
                        for g in tops
                    )
                    pw_dict = c.get("pathways", {}) or {}
                    n_pw = sum(len(v) if isinstance(v, list) else 0
                               for v in pw_dict.values())
                    pw_top = []
                    for db in list(pw_dict.keys())[:3]:
                        terms = pw_dict.get(db) or []
                        if terms and isinstance(terms, list):
                            pw_top.append(f"{db}: {terms[0].get('term','?')[:60]}")
                    lines.append(
                        f"  • {name}: {n_sig} DE genes ({n_up} up, {n_dn} down). "
                        f"Top: {top_str}."
                    )
                    if pw_top:
                        lines.append(
                            f"    Enriched pathways ({n_pw} total): "
                            + "; ".join(pw_top)
                        )
                # Overlap
                overlap = findings.get("overlap", {}) or {}
                for pair, info in list(overlap.items())[:2]:
                    lines.append(
                        f"  • Shared DE genes [{pair}]: "
                        f"{info.get('n_shared',0)} "
                        f"(Jaccard={info.get('jaccard',0)})"
                    )
                # Sample QC
                sqc = findings.get("sample_qc", {}) or {}
                if sqc.get("n_samples"):
                    lines.append(
                        f"  • Sample QC: {sqc['n_samples']} samples, "
                        f"library-size range {sqc.get('size_ratio', 0)}×, "
                        f"{len(sqc.get('outliers', []))} outliers removed."
                    )
            else:
                lines.append("BULK RNA-seq: ran but no successful contrasts.")

        # scRNA
        sc = agent_results.get("scrna_agent", {})
        if sc.get("status") == "done":
            from aria.agents import _narrative_scrna
            sc_f = _narrative_scrna.unwrap_scrna_findings(sc)
            # Keep f around for legacy fallback keys (n_cells_after_qc etc.)
            f    = sc.get("findings", {}) or {}
            qc   = sc_f.get("qc", {}) or {}
            clus = sc_f.get("clustering", {}) or \
                sc_f.get("clustering_decision", {}) or {}
            n_cells    = qc.get("n_cells_after") or f.get("n_cells_after_qc", "?")
            n_clusters = clus.get("n_clusters") or f.get("n_clusters", "?")
            ct = sc_f.get("cell_types", {}) or {}
            ct_raw = ct.get("cell_types", {}) or {}
            ct_list = list({
                v.get("cell_type", "") if isinstance(v, dict) else str(v)
                for v in ct_raw.values()
            } - {""})[:5]
            ct_str = f" Cell types: {', '.join(ct_list)}." if ct_list else ""

            # Integration
            integ = sc_f.get("integration", {}) or {}
            integ_str = ""
            if integ.get("status") == "done":
                integ_str = (
                    f" Harmony batch correction: {integ.get('n_batches','?')} batches, "
                    f"silhouette {integ.get('silhouette_before','?')} → "
                    f"{integ.get('silhouette_after','?')}."
                )

            # Trajectory
            traj = sc_f.get("trajectory", {}) or {}
            traj_str = ""
            if traj.get("status") in ("done", "success"):
                pt = traj.get("pseudotime", {}) or {}
                paga = traj.get("paga", {}) or {}
                n_conn = paga.get("n_connections") or \
                    len(paga.get("top_connections", {}) or {})
                vel = traj.get("velocity", {}) or {}
                traj_str = (
                    f" PAGA: {n_conn} cluster-pair connections, "
                    f"{paga.get('n_strong', 0)} above threshold. "
                    f"DPT pseudotime computed: {pt.get('computed', False)}."
                    + (" RNA velocity computed." if vel.get("computed") else "")
                )

            # Cell-cell communication
            ccc = sc_f.get("cell_communication", {}) or {}
            ccc_str = ""
            if ccc.get("status") in ("done", "success"):
                top_p = ccc.get("top_pairs", [])[:3]
                ccc_str = (
                    f" Cell-cell comm ({ccc.get('method','?')}): "
                    f"{ccc.get('n_interactions','?')} L-R interactions. "
                    f"Top pairs: {', '.join(top_p)}."
                ) if top_p else (
                    f" Cell-cell comm: {ccc.get('n_interactions','?')} interactions."
                )

            # Pseudobulk DE and pathway enrichment
            pb = sc_f.get("pseudobulk_de", {}) or {}
            pb_str = ""
            if pb:
                per_group = pb.get("per_group", {}) or {}
                n_success = sum(
                    1 for g in per_group.values()
                    for c in (g.get("per_comparison", {}) or {}).values()
                    if c.get("status") == "success"
                )
                n_with_de = sum(
                    1 for g in per_group.values()
                    for c in (g.get("per_comparison", {}) or {}).values()
                    if c.get("status") == "success"
                    and c.get("n_significant", 0) > 0
                )
                top_blocks = _narrative_scrna._top_de_blocks(pb, limit=3)
                top_str = "; ".join(
                    f"{g} {cmp}: {c.get('n_significant', 0)} DE"
                    for g, cmp, c in top_blocks
                )
                pb_str = (
                    f" Pseudobulk DE: {pb.get('n_groups', '?')} "
                    f"{pb.get('groupby', 'group')} groups, {n_success} "
                    f"successful group x comparison blocks, {n_with_de} "
                    f"with significant DE."
                    + (f" Largest blocks: {top_str}." if top_str else "")
                )

            pwp = sc_f.get("pseudobulk_pathways", {}) or {}
            pw_str = ""
            if pwp.get("per_cluster"):
                n_blocks = len(pwp["per_cluster"])
                n_sig = sum(
                    1 for b in pwp["per_cluster"].values()
                    if b.get("n_significant", 0) > 0
                )
                pw_str = (
                    f" Pathway ORA: {n_sig}/{n_blocks} DE blocks "
                    f"with significant enrichment."
                )

            lines.append(
                f"scRNA-seq: {n_cells} cells after QC, "
                f"{n_clusters} clusters identified.{ct_str}"
                f"{integ_str}{pb_str}{pw_str}{traj_str}{ccc_str}"
            )

        # Chromatin
        chrom = agent_results.get("chromatin_agent", {})
        if chrom.get("status") == "done":
            chrom_f = chrom.get("findings", {}) or {}
            chrom_lines = []
            for assay in ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN"):
                af = chrom_f.get(assay, {}).get("findings", {})
                peaks = af.get("peaks", {})
                if peaks.get("n_peaks"):
                    frip = peaks.get("frip", "?")
                    frip_str = (f"{frip:.2f}" if isinstance(frip, float)
                                else str(frip))
                    chrom_lines.append(
                        f"{assay}: {peaks['n_peaks']:,} peaks "
                        f"(FRiP={frip_str})"
                    )
            if chrom_lines:
                lines.append("CHROMATIN: " + "; ".join(chrom_lines) + ".")

        # 3D genome / HiC
        hic = agent_results.get("genome_arch_agent", {})
        if hic.get("status") == "done":
            hic_f  = hic.get("findings", {}) or {}
            topo   = hic_f.get("topology", {}) or {}
            tads   = topo.get("tads", {})
            comp   = topo.get("compartments", {})
            loops  = topo.get("loops", {})
            hic_parts = []
            if tads.get("n_tads"):
                hic_parts.append(f"{tads['n_tads']:,} TADs")
            if comp.get("pct_A"):
                hic_parts.append(
                    f"compartments {comp['pct_A']:.1f}% A / "
                    f"{comp.get('pct_B', 0):.1f}% B"
                )
            if loops.get("n_loops"):
                hic_parts.append(f"{loops['n_loops']:,} loops")
            if hic_parts:
                lines.append("3D GENOME: " + ", ".join(hic_parts) + ".")

        return "\n".join(lines) if lines else "(no concrete outputs recorded)"

    def _write_findings_sections(self, grouped: dict,
                                  agent_results: dict,
                                  exp_ctx: dict) -> dict:
        """
        Write structured findings for each modality.
        Returns dict of section_name -> prose text.
        """
        sections = {}

        # Bulk RNA (v3: multi-contrast aware)
        bulk = agent_results.get("bulk_rna_agent", {})
        if bulk.get("status") == "done":
            sections["bulk_rna"] = self._summarize_bulk_rna(bulk, grouped)

        # scRNA
        scrna = agent_results.get("scrna_agent", {})
        if scrna.get("status") == "done":
            sections["scrna"] = self._summarize_rna(scrna, grouped)

        # Legacy "rna_agent" key (backward compatibility)
        rna = agent_results.get("rna_agent", {})
        if rna.get("status") == "done" and "scrna" not in sections:
            sections["scrna"] = self._summarize_rna(rna, grouped)

        # Chromatin findings
        chrom = agent_results.get("chromatin_agent", {})
        if chrom.get("status") == "done":
            sections["chromatin"] = self._summarize_chromatin(chrom, grouped)

        # 3D genome findings
        hic = agent_results.get("genome_arch_agent", {})
        if hic.get("status") == "done":
            sections["hic"] = self._summarize_hic(hic, grouped)

        # Integration findings
        integ = agent_results.get("integration_agent", {})
        if integ.get("status") == "done":
            sections["integration"] = self._summarize_integration(
                integ, grouped
            )

        # Final structured synthesis for single-cell reports
        sc = agent_results.get("scrna_agent") or agent_results.get("rna_agent", {})
        if sc.get("status") == "done":
            from aria.agents import _narrative_scrna
            sc_findings = _narrative_scrna.unwrap_scrna_findings(sc)
            synthesis = _narrative_scrna.build_scrna_integrated_interpretation(
                sc_findings,
                {"summary": exp_ctx.get("user_question", "")},
            )
            if synthesis:
                sections["synthesis"] = synthesis

        # Conflicts and limitations
        sections["conflicts"] = self._summarize_conflicts(
            agent_results, grouped
        )

        return sections

    def _write_methods_section(self, exp_ctx: dict,
                                agent_results: dict,
                                decisions: list) -> str:
        """
        Write a reproducible methods section.
        Uses exact parameter values from ParameterAdvisor decisions.
        This section should be copy-pasteable into a manuscript Methods.
        """
        lines = []

        organism = exp_ctx.get("organism", "unknown organism")
        genome   = exp_ctx.get("genome", "unknown assembly")

        lines.append(f"**Data processing**\n")
        lines.append(
            f"Genomic analyses were performed on {organism} data "
            f"aligned to the {genome} reference assembly."
        )

        # Bulk RNA methods
        bulk = agent_results.get("bulk_rna_agent", {})
        if bulk.get("status") == "done":
            findings_bulk = bulk.get("findings", {})
            preprocessing = findings_bulk.get("preprocessing", {})
            contrasts     = findings_bulk.get("contrasts", [])
            lfc_thr       = findings_bulk.get("lfc_threshold",  1.0)
            padj_thr_used = findings_bulk.get("padj_threshold", 0.05)

            lines.append(f"\n**Bulk RNA-seq**\n")

            # Preprocessing methods
            if preprocessing:
                qc_info = preprocessing.get("qc", {})
                lines.append(
                    f"Raw paired-end FASTQ reads were trimmed with fastp "
                    f"(default parameters, min read length 36 bp, "
                    f"Phred quality ≥ 20), "
                    f"retaining an average of "
                    f"{self._avg_pct_passed(qc_info):.1f}% of reads across "
                    f"{qc_info.get('n_samples', '?')} samples. "
                    f"Per-sample MultiQC reports were generated."
                )
                lines.append(
                    f"Reads were aligned to the {genome} reference genome "
                    f"(Ensembl release 112) using STAR (2-pass mode) with "
                    f"--outSAMtype BAM SortedByCoordinate and --quantMode "
                    f"GeneCounts enabled. Gene-level counts were computed "
                    f"with featureCounts (subread) using the Ensembl GTF "
                    f"annotation; paired-end reads were counted as "
                    f"fragments."
                )

            # DE methods
            lines.append(
                f"Differential expression was performed with DESeq2 "
                f"(via pydeseq2). Low-count genes "
                f"(<10 counts in fewer than 25% of samples) were "
                f"filtered. Cook's distance was used to refit "
                f"outlier counts. Sample QC used PCA-based outlier "
                f"detection (>2 SD from group centroid in PC1–PC2)."
            )

            if contrasts:
                contrast_names = [c.get("name", "?") for c in contrasts
                                   if c.get("status") == "success"]
                design_used = findings_bulk.get("design_used", "~condition")
                if len(contrast_names) > 1:
                    lines.append(
                        f"Contrasts tested: {'; '.join(contrast_names)}. "
                        f"Each contrast was analyzed independently using "
                        f"the Wald test with the design {design_used}. "
                        f"Significance thresholds: adjusted p-value < {padj_thr_used} "
                        f"(Benjamini–Hochberg) and |log2 fold-change| > "
                        f"{lfc_thr}."
                    )
                else:
                    lines.append(
                        f"Significance thresholds: adjusted p-value < {padj_thr_used} "
                        f"and |log2 fold-change| > {lfc_thr}."
                    )

                if lfc_thr < 0.8:
                    lines.append(
                        f"Note: a lower log2FC threshold of {lfc_thr} "
                        f"(≈1.5× fold change) was used because the "
                        f"biological question involves transcription "
                        f"factor perturbation, where direct target genes "
                        f"typically show modest effect sizes."
                    )
                powers = [
                    c.get("power_estimate_at_lfc_min")
                    for c in contrasts
                    if c.get("status") == "success"
                    and isinstance(c.get("power_estimate_at_lfc_min"), (int, float))
                ]
                if powers:
                    lines.append(
                        f"Approximate power to detect |log2FC| > {lfc_thr} "
                        f"at adjusted alpha {padj_thr_used} ranged from "
                        f"{min(powers):.0%} to {max(powers):.0%} across "
                        f"contrasts, using a negative-binomial Wald "
                        f"approximation from replicate count, mean expression, "
                        f"and dispersion."
                    )

            # Pathway enrichment methods
            has_pathways = any(
                c.get("pathways") for c in contrasts
                if c.get("status") == "success"
            )
            if has_pathways:
                bg_sizes = [
                    ((c.get("pathway_background") or {}).get("background_size"))
                    for c in contrasts
                    if (c.get("pathway_background") or {}).get("background_size")
                ]
                bg_clause = (
                    f" ORA used the dataset-expressed background "
                    f"({min(bg_sizes)}-{max(bg_sizes)} genes across contrasts) "
                    f"instead of Enrichr's default universe."
                    if bg_sizes else
                    " ORA used Enrichr's default universe because no "
                    "dataset-expressed background was available."
                )
                lines.append(
                    f"Over-representation analysis was performed with "
                    f"gseapy (Enrichr endpoint) against GO Biological "
                    f"Process, KEGG, and Reactome gene sets. "
                    f"Significance was defined as adjusted p-value < 0.05."
                    f"{bg_clause}"
                )

        # Single-cell RNA methods (delegated to _narrative_scrna for full
        # pseudobulk/standard/trajectory/cellcomm handling)
        scrna = agent_results.get("scrna_agent",
                                    agent_results.get("rna_agent", {}))
        if scrna.get("status") == "done":
            from aria.agents import _narrative_scrna
            sc_findings = _narrative_scrna.unwrap_scrna_findings(scrna)
            scrna_methods = _narrative_scrna.build_scrna_methods(sc_findings)
            if scrna_methods:
                lines.append(f"\n**Single-cell RNA-seq**\n")
                lines.append(scrna_methods)

        # Chromatin methods
        chrom = agent_results.get("chromatin_agent", {})
        if chrom.get("status") == "done":
            lines.append(f"\n**Chromatin accessibility**\n")
            lines.append(
                f"ATAC-seq peaks were called using MACS3 with parameters "
                f"appropriate for the assay type (--nomodel --extsize 200 "
                f"for scATAC-seq; --nolambda for CUT&RUN/CUT&TAG). "
                f"Quality control thresholds: FRiP > 0.20, "
                f"TSS enrichment score > 4.0 (ENCODE standards). "
                f"For single-cell ATAC-seq, dimensionality reduction used "
                f"LSI (TF-IDF normalization + SVD), discarding the first "
                f"component which correlates with sequencing depth."
            )

        # HiC methods
        hic = agent_results.get("genome_arch_agent", {})
        if hic.get("status") == "done":
            lines.append(f"\n**3D genome organization**\n")
            lines.append(
                f"Hi-C contact matrices were balanced using ICE normalization "
                f"(cooler balance). Compartment A/B analysis was performed "
                f"via eigenvector decomposition (PC1) at 100 kb resolution. "
                f"Topologically associating domains were identified using "
                f"the Insulation Score method (Crane et al., 2015). "
                f"Window size was calibrated on chromosome 1 as a proxy "
                f"for full-genome optimization."
            )

        # Integration methods
        integ = agent_results.get("integration_agent", {})
        if integ.get("status") == "done":
            strategy = integ.get("strategy", "")
            lines.append(f"\n**Multi-omics integration**\n")
            if "wnn" in strategy:
                lines.append(
                    f"RNA and ATAC data from the same cells were integrated "
                    f"using Weighted Nearest Neighbors (WNN), implemented "
                    f"in muon. Modality weights were computed per cell based "
                    f"on local prediction accuracy."
                )
            if "mofa" in strategy:
                lines.append(
                    f"Multi-omics factor analysis (MOFA+) was applied to "
                    f"decompose variance across modalities into shared and "
                    f"modality-specific factors."
                )

        # Parameter decisions appendix reference
        if decisions:
            lines.append(
                f"\n**Reproducibility**\n"
                f"All {len(decisions)} parameter decisions are logged in "
                f"ARIA's institutional memory database with full justification. "
                f"See Supplementary Appendix for the complete decisions log."
            )

        return "\n".join(lines)

    # ── Per-agent summarizers ─────────────────────────────────────────────

    @staticmethod
    def _avg_pct_passed(qc_info: dict) -> float:
        """Average percentage of reads passing QC across samples."""
        samples = qc_info.get("samples", [])
        if not samples:
            return 0.0
        pct_list = [s.get("pct_passed", 0) for s in samples]
        return sum(pct_list) / max(len(pct_list), 1)

    def _summarize_bulk_rna(self, bulk_result: dict,
                              grouped: dict) -> str:
        """
        Summarize bulk RNA-seq findings. Handles multi-contrast format (v3).

        Returns prose describing:
          - Sample QC (libraries, outliers)
          - Per-contrast results (DE counts, top genes)
          - Cross-contrast overlap (shared vs unique biology)
          - LLM interpretation (if available)
        """
        findings = bulk_result.get("findings", {})
        parts    = []

        # Preprocessing summary (if raw FASTQs were processed)
        preprocessing = findings.get("preprocessing", {})
        if preprocessing:
            qc = preprocessing.get("qc", {})
            align = preprocessing.get("alignment", {})
            if qc.get("n_samples"):
                parts.append(
                    f"Raw FASTQ processing: "
                    f"{qc.get('n_samples', '?')} samples trimmed with fastp."
                )
            if align.get("n_aligned"):
                parts.append(
                    f"Reads aligned with STAR "
                    f"({align.get('n_aligned', '?')} samples, 2-pass mode)."
                )

        # Sample QC
        sqc = findings.get("sample_qc", {})
        if sqc:
            n_samples = sqc.get("n_samples", "?")
            outliers  = sqc.get("outliers", [])
            lib_range = sqc.get("size_ratio", 1)
            if outliers:
                parts.append(
                    f"PCA-based QC on {n_samples} samples identified "
                    f"outliers: {outliers}, which were excluded from "
                    f"differential expression."
                )
            else:
                parts.append(
                    f"PCA-based QC on {n_samples} samples passed; "
                    f"library size range {lib_range:.1f}× "
                    f"(no outliers removed)."
                )

        # Multi-contrast results
        contrasts = findings.get("contrasts", [])
        lfc_thr   = findings.get("lfc_threshold", 1.0)

        if contrasts:
            n_ok = sum(1 for c in contrasts if c.get("status") == "success")
            padj_str = f"padj < 0.05, |log2FC| > {lfc_thr}"

            if len(contrasts) == 1:
                c = contrasts[0]
                lp_tag = " [LOW POWER]" if c.get("low_power_warning") else ""
                parts.append(
                    f"Differential expression ({c.get('name', '?')}{lp_tag}, "
                    f"DESeq2, {padj_str}): "
                    f"{c.get('n_significant', 0)} significant genes "
                    f"({c.get('n_upregulated', 0)} upregulated, "
                    f"{c.get('n_downregulated', 0)} downregulated)."
                )
                if c.get("low_power_warning") and c.get("low_power_reason"):
                    parts.append(f"Caveat: {c['low_power_reason']}")
            else:
                parts.append(
                    f"Differential expression was performed across "
                    f"{n_ok} contrasts (DESeq2, {padj_str}):"
                )
                for c in contrasts:
                    if c.get("status") != "success":
                        continue
                    lp_tag = " [LOW POWER]" if c.get("low_power_warning") else ""
                    parts.append(
                        f"  • {c.get('name', '?')}{lp_tag}: "
                        f"{c.get('n_significant', 0)} DE genes "
                        f"({c.get('n_upregulated', 0)} up, "
                        f"{c.get('n_downregulated', 0)} down)."
                    )
                n_low_power = sum(
                    1 for c in contrasts if c.get("low_power_warning")
                )
                if n_low_power:
                    parts.append(
                        f"Caveat: {n_low_power} of {n_ok} contrasts ran with "
                        f"n<=2 replicates on at least one side. Dispersion is "
                        f"poorly estimated and FDR is unreliable for those "
                        f"contrasts. Interpret with caution."
                    )

            # Cross-contrast overlap (uses FULL DE list now, not top 30)
            overlap = findings.get("overlap", {})
            if overlap and len(contrasts) > 1:
                for pair_name, info in list(overlap.items())[:3]:
                    shared  = info.get("n_shared", 0)
                    n_a     = info.get("n_in_first", 0)
                    n_b     = info.get("n_in_second", 0)
                    jaccard = info.get("jaccard", 0)
                    if shared > 0:
                        # Show first few shared genes (already in symbols
                        # if symbol_map was available)
                        examples = info.get("shared_genes", [])[:8]
                        ex_str   = (f" Examples: {', '.join(examples)}."
                                     if examples else "")
                        parts.append(
                            f"Shared DE genes [{pair_name}]: "
                            f"{shared} genes (Jaccard={jaccard}, "
                            f"{n_a} and {n_b} per contrast).{ex_str}"
                        )

            # Top genes per contrast — by HGNC symbol when available
            for c in contrasts:
                if c.get("status") != "success":
                    continue
                tops = c.get("top_genes", [])[:8]
                if not tops:
                    continue
                gene_strs = []
                for g in tops:
                    sym = g.get("symbol") or g.get("gene", "?")
                    lfc = g.get("log2fc", 0)
                    arrow = "↑" if lfc > 0 else "↓"
                    gene_strs.append(f"{sym} ({arrow}{abs(lfc):.1f})")
                parts.append(
                    f"Top DE genes [{c.get('name')}]: "
                    + ", ".join(gene_strs[:6]) + "."
                )

            # Top pathways across ALL contrasts (not just first)
            for c in contrasts:
                if c.get("status") != "success":
                    continue
                pw_dict = c.get("pathways", {})
                if not pw_dict:
                    continue
                top_dbs = list(pw_dict.keys())[:3]
                examples = []
                for db in top_dbs:
                    terms = pw_dict.get(db, [])
                    if isinstance(terms, list) and terms:
                        examples.append(
                            f"{db}: {terms[0].get('term', '?')[:60]}"
                        )
                if examples:
                    parts.append(
                        f"Top enriched pathways [{c.get('name')}]: "
                        + "; ".join(examples) + "."
                    )

        else:
            # Legacy single-contrast fallback
            n_sig = findings.get("n_significant", 0)
            if n_sig:
                comp = findings.get("comparison_used", {})
                parts.append(
                    f"Differential expression "
                    f"({comp.get('numerator','?')} vs "
                    f"{comp.get('denominator','?')}): "
                    f"{n_sig} significant genes "
                    f"({findings.get('n_upregulated', 0)} up, "
                    f"{findings.get('n_downregulated', 0)} down)."
                )

        # LLM interpretation (already generated by BulkRNAAgent)
        interpretation = findings.get("interpretation", "")
        if interpretation and isinstance(interpretation, str):
            parts.append("\n" + self._guard_bulk_interpretation(
                interpretation.strip()
            ))

        return "\n".join(parts) if parts else \
               "Bulk RNA-seq analysis completed. See findings table for details."

    # ── Legacy scRNA summarizer ───────────────────────────────────────────

    def _summarize_rna(self, rna_result: dict,
                        grouped: dict) -> str:
        from aria.agents import _narrative_scrna
        # Delegate to the scRNA module: it knows about pseudobulk + standard
        # + trajectory + cellcomm. unwrap_scrna_findings is robust to both
        # the adapter (wrapped) and direct scrna_agent.run() (unwrapped) shapes.
        sc = _narrative_scrna.unwrap_scrna_findings(rna_result)
        return _narrative_scrna.summarize_scrna_text(sc)

    def _summarize_chromatin(self, chrom_result: dict,
                              grouped: dict) -> str:
        findings = chrom_result.get("findings", {})
        parts    = []

        for assay in ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN"):
            assay_f = (findings.get(assay) or {}).get("findings", {}) or {}
            peaks   = assay_f.get("peaks", {})
            if peaks.get("n_peaks"):
                parts.append(
                    f"{assay}: {peaks['n_peaks']:,} peaks called "
                    f"(FRiP={peaks.get('frip', '?'):.2f})."
                )
            motifs = assay_f.get("motifs", {})
            if motifs.get("top_motifs"):
                top = motifs["top_motifs"][:3]
                conf = motifs.get("confidence", "medium")
                parts.append(
                    f"TF motif enrichment ({conf} confidence): "
                    f"{', '.join(top)}."
                )

        return " ".join(parts) if parts else \
               "Chromatin analysis completed. See findings table for details."

    def _summarize_hic(self, hic_result: dict,
                        grouped: dict) -> str:
        findings = hic_result.get("findings", {})
        topo     = findings.get("topology", {})
        parts    = []

        comp = topo.get("compartments", {})
        if comp.get("pct_A"):
            parts.append(
                f"A/B compartment analysis: "
                f"{comp['pct_A']:.1f}% A-compartment, "
                f"{comp['pct_B']:.1f}% B-compartment."
            )
            n_sw = comp.get("n_ab_switches", 0)
            if n_sw:
                parts.append(f"{n_sw} compartment switches detected.")

        tads = topo.get("tads", {})
        if tads.get("n_tads"):
            parts.append(
                f"{tads['n_tads']:,} TADs identified "
                f"(median size: {tads.get('median_size_kb', '?')} kb)."
            )

        loops = topo.get("loops", {})
        if loops.get("n_loops"):
            parts.append(
                f"{loops['n_loops']:,} chromatin loops called."
            )

        return " ".join(parts) if parts else \
               "3D genome analysis completed. See findings table for details."

    def _summarize_integration(self, integ_result: dict,
                                grouped: dict) -> str:
        findings = integ_result.get("findings", {})
        parts    = []

        wnn = findings.get("wnn", {})
        if wnn.get("status") == "success":
            rna_w = wnn.get("mean_rna_weight", 0)
            atac_w = wnn.get("mean_atac_weight", 0)
            disc  = wnn.get("n_discordant_clusters", 0)
            parts.append(
                f"WNN integration: RNA weight={rna_w:.2f}, "
                f"ATAC weight={atac_w:.2f}."
            )
            if disc:
                parts.append(
                    f"{disc} clusters differed between RNA-only "
                    f"and joint embedding, indicating chromatin-driven "
                    f"cell identity information."
                )

        p2g = findings.get("peak_to_gene", {})
        if p2g.get("n_links"):
            n_neg = p2g.get("n_negative_correlations", 0)
            parts.append(
                f"Peak-to-gene analysis: {p2g['n_links']} regulatory links. "
            )
            if n_neg:
                parts.append(
                    f"{n_neg} open-chromatin/silent-gene pairs identified "
                    f"(candidate poised regulatory elements)."
                )

        mofa = findings.get("mofa", {})
        if mofa.get("n_factors"):
            if mofa.get("cell_cycle_factor"):
                fid = mofa.get("cell_cycle_factor_id", "?")
                parts.append(
                    f"MOFA+ Factor {fid} captures cell cycle variation "
                    f"(excluded from biological interpretation)."
                )

        return " ".join(parts) if parts else \
               "Integration analysis completed. See findings table for details."

    def _summarize_conflicts(self, agent_results: dict,
                              grouped: dict) -> str:
        """Honest accounting of conflicts and limitations."""
        conflicts = agent_results.get("integration_agent", {}).get(
            "findings", {}
        ).get("conflicts", {}).get("conflicts", [])

        low_findings = grouped["low"]
        insuff       = grouped["insufficient"]

        parts = []

        if conflicts:
            parts.append(
                f"**Cross-modal conflicts ({len(conflicts)}):**"
            )
            for c in conflicts[:3]:
                parts.append(f"- {c.get('description', str(c))[:150]}")

        if low_findings:
            parts.append(
                f"\n**Findings requiring validation ({len(low_findings)}):**"
            )
            for f in low_findings[:3]:
                parts.append(f"- {f.get('summary', str(f))[:120]}")

        if insuff:
            parts.append(
                f"\n**Analyses with insufficient data "
                f"({len(insuff)} findings):**"
            )
            parts.append(
                "The following analyses could not be completed due to "
                "insufficient data or failed quality control. "
                "No conclusions should be drawn from these:"
            )
            for f in insuff[:3]:
                parts.append(f"- {f.get('summary', str(f))[:120]}")

        if not parts:
            return "No cross-modal conflicts identified."

        return "\n".join(parts)

    # ── HTML rendering ────────────────────────────────────────────────────

    def _build_report_dir(self, experiment_id: str,
                           intent: dict, exp_ctx: dict) -> Path:
        """
        Construct (and mkdir) ~/.aria/reports/aria_<ts>_<slug>_<suffix>/
        with figures/ + tables/ subdirs. Idempotent: call once early in
        run() so figure subprocesses have a stable target.
        """
        reproducible = bool((exp_ctx or {}).get("reproducible_mode"))
        ts        = (
            "reproducible"
            if reproducible else
            datetime.now().strftime("%Y%m%d_%H%M%S")
        )
        slug      = self._build_slug(intent, exp_ctx)
        if reproducible:
            first_hash = ""
            for rec in (exp_ctx or {}).get("input_files", []) or []:
                if rec.get("sha256") and rec.get("sha256") != "unavailable":
                    first_hash = str(rec["sha256"])[:12]
                    break
            suffix = first_hash or "nohash"
        else:
            suffix    = (experiment_id[-4:] if len(experiment_id) >= 4
                          else experiment_id)
        report_name = (f"aria_{ts}_{slug}_{suffix}" if slug
                        else f"aria_{ts}_{suffix}")
        report_dir = self.reports_dir / report_name
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "figures").mkdir(exist_ok=True)
        (report_dir / "tables").mkdir(exist_ok=True)
        return report_dir

    def _generate_scrna_figures(self, agent_results: dict,
                                  report_dir: Path) -> None:
        """
        Orchestrate scRNA figure rendering (UMAPs + DE bar + pathway dotplots
        + cellcomm + trajectory). Mutates agent_results so that the in-memory
        scrna_agent envelope carries the figure paths into the section
        builders that follow.
        """
        sc = (agent_results or {}).get("scrna_agent") or {}
        if sc.get("status") != "done":
            return
        try:
            from aria.agents import _narrative_scrna
            from aria.utils.environment_manager import env_manager
        except Exception as e:
            log.warning(f"Cannot import figure helpers: {e}")
            return
        findings = _narrative_scrna.unwrap_scrna_findings(sc)
        h5ad = sc.get("output_h5ad") or sc.get("output_path")
        try:
            _narrative_scrna.generate_figures(
                findings=findings,
                h5ad_path=h5ad,
                output_dir=report_dir / "figures",
                env_manager=env_manager,
            )
            # `findings` may be a reference to sc["findings"] (legacy shape)
            # or a wrapped child (multimodal shape). Either way the mutation
            # above already wrote into the dict that downstream readers see.
        except Exception as e:
            log.warning(f"scRNA figure generation failed: {e}", exc_info=True)

    def _render_html_report(self, experiment_id: str,
                             exp_ctx: dict,
                             intent: dict,
                             executive_summary: str,
                             findings_sections: dict,
                             grouped_findings: dict,
                             methods: str,
                             decisions: list,
                             agent_results: dict = None,
                             report_dir: Optional[Path] = None) -> Path:
        """
        Render the full HTML report.
        Self-contained: CSS embedded, no external dependencies.
        """
        organism = exp_ctx.get("organism", "Unknown organism")
        genome   = exp_ctx.get("genome", "Unknown assembly")
        question = intent.get("summary",
                               exp_ctx.get("user_question", ""))
        agent_results = agent_results or {}
        reproducible = bool(exp_ctx.get("reproducible_mode"))
        date_str = (
            "&lt;timestamp redacted for byte-identity&gt;"
            if reproducible else
            datetime.now().strftime("%B %d, %Y")
        )
        exp_short = experiment_id[:8]

        # Build/stage report artifacts before composing section HTML. scRNA
        # tables are generated from in-memory results here, and the findings
        # object is annotated with links used by build_scrna_html_section().
        if report_dir is None:
            report_dir = self._build_report_dir(experiment_id, intent, exp_ctx)
        report_dir.mkdir(parents=True, exist_ok=True)
        self._stage_artifacts(agent_results, report_dir)

        # Findings table rows
        findings_rows = self._build_findings_table(grouped_findings)

        # Decisions log rows
        decisions_rows = self._build_decisions_table(decisions)

        version = _html.escape(str(ARIA_VERSION))
        conflicts_html = self._plain_text_to_html(
            findings_sections.get("conflicts", "No conflicts detected.")
        )
        methods_html = self._plain_text_to_html(methods)
        provenance = exp_ctx.get("provenance") or collect_provenance()
        llm_usage_since = provenance.get("timestamp_utc")
        llm_usage = collect_llm_usage(llm_usage_since)
        if reproducible:
            provenance = dict(provenance)
            provenance["timestamp_utc"] = "<timestamp redacted for byte-identity>"
        provenance_html = self._build_provenance_section(
            provenance=provenance,
            input_files=exp_ctx.get("input_files", []),
            agent_results=agent_results,
            llm_usage=llm_usage,
        )

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ARIA Report — {exp_short}</title>
<style>
  :root {{
    --bg:     #ffffff;
    --bg-alt: #f8fafc;
    --panel:  #f1f5f9;
    --text:   #1e293b;
    --muted:  #475569;
    --dim:    #94a3b8;
    --border: #e2e8f0;
    --navy:   #0f172a;
    --blue:   #1d4ed8;
    --teal:   #0d9488;
    --green:  #15803d;
    --amber:  #92400e;
    --red:    #991b1b;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.65;
    padding: 3rem 2rem;
    max-width: 900px;
    margin: 0 auto;
    font-size: 15px;
  }}
  h1 {{
    color: var(--navy);
    font-size: 1.7rem;
    font-weight: 700;
    margin-bottom: 0.25rem;
  }}
  h2 {{
    color: var(--navy);
    font-size: 0.85rem;
    font-weight: 700;
    margin: 2.5rem 0 0.8rem;
    border-bottom: 2px solid var(--navy);
    padding-bottom: 0.35rem;
    text-transform: uppercase;
    letter-spacing: 0.06em;
  }}
  h3 {{
    color: var(--muted);
    font-size: 0.78rem;
    font-weight: 700;
    margin: 1.2rem 0 0.4rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
  }}
  h4 {{
    color: var(--navy);
    font-size: 0.95rem;
    font-weight: 600;
    margin: 1.4rem 0 0.4rem;
  }}
  p {{ margin: 0.6rem 0; }}
  .meta {{
    color: var(--muted);
    font-size: 0.85rem;
    margin-bottom: 2.5rem;
    padding-bottom: 1rem;
    border-bottom: 1px solid var(--border);
  }}
  .card {{
    background: var(--bg-alt);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 1.25rem 1.5rem;
    margin: 0.8rem 0;
  }}
  .badge {{
    display: inline-block;
    padding: 0.2rem 0.65rem;
    border-radius: 4px;
    font-size: 0.7rem;
    font-weight: 700;
    letter-spacing: 0.05em;
    margin-right: 0.4rem;
    text-transform: uppercase;
  }}
  .high   {{ background: #dcfce7; color: var(--green);  border: 1px solid #86efac; }}
  .medium {{ background: #fef3c7; color: var(--amber);  border: 1px solid #fcd34d; }}
  .low    {{ background: #fee2e2; color: var(--red);    border: 1px solid #fca5a5; }}
  .insuff {{ background: var(--panel); color: var(--dim); border: 1px solid var(--border); }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.875rem;
    margin: 0.8rem 0;
  }}
  th {{
    background: var(--navy);
    color: #ffffff;
    text-align: left;
    padding: 0.55rem 0.75rem;
    font-weight: 600;
    font-size: 0.78rem;
    text-transform: uppercase;
    letter-spacing: 0.04em;
  }}
  td {{
    border-bottom: 1px solid var(--border);
    padding: 0.5rem 0.75rem;
    vertical-align: top;
  }}
  tr:nth-child(even) td {{ background: var(--bg-alt); }}
  pre {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-left: 4px solid var(--teal);
    border-radius: 4px;
    padding: 1rem 1.25rem;
    font-size: 0.82rem;
    white-space: pre-wrap;
    color: var(--text);
    margin: 0.8rem 0;
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', monospace;
  }}
  .warning {{
    border-left: 4px solid var(--amber);
    padding: 0.75rem 1rem;
    background: #fffbeb;
    color: var(--amber);
    font-size: 0.9rem;
    margin: 0.8rem 0;
    border-radius: 0 4px 4px 0;
  }}
  footer {{
    margin-top: 4rem;
    padding-top: 1.5rem;
    border-top: 2px solid var(--border);
    color: var(--dim);
    font-size: 0.8rem;
    text-align: center;
  }}
  a {{ color: var(--blue); text-decoration: underline; }}
  figure {{ margin: 1.5rem 0; text-align: center; }}
  figcaption {{
    color: var(--muted);
    font-size: 0.82rem;
    margin-top: 0.4rem;
    font-style: italic;
  }}
  figure img {{
    max-width: 100%;
    height: auto;
    border-radius: 4px;
    border: 1px solid var(--border);
  }}
  @media print {{
    body {{ padding: 1rem; font-size: 11pt; }}
    h2 {{ page-break-after: avoid; }}
    .card {{ border: 1px solid #ccc; break-inside: avoid; }}
    footer {{ margin-top: 2rem; }}
  }}
</style>
</head>
<body>

<h1>ARIA Analysis Report</h1>
<p class="meta">
  Experiment: <strong>{exp_short}</strong> &nbsp;|&nbsp;
  {organism} / {genome} &nbsp;|&nbsp;
  {date_str} &nbsp;|&nbsp;
  Generated by ARIA v{version}
</p>

<div class="card">
  <h3>Biological Question</h3>
  <p><em>{_html.escape(str(question))}</em></p>
</div>

<h2>Provenance</h2>
{provenance_html}

<h2>Executive Summary</h2>
<div class="card">
  <p>{self._plain_text_to_html(executive_summary)}</p>
</div>

<h2>Quality Control Summary</h2>
{self._build_qc_section(grouped_findings, agent_results, exp_ctx)}

<h2>Findings</h2>
{self._build_findings_section(findings_sections, agent_results)}

<h2>All Findings ({sum(len(v) for v in grouped_findings.values())} total)</h2>
<table>
  <tr>
    <th>Confidence</th>
    <th>Finding</th>
    <th>Agent</th>
  </tr>
  {findings_rows}
</table>

<h2>Cross-modal Conflicts &amp; Limitations</h2>
<div class="card">
  <pre>{conflicts_html}</pre>
</div>

<h2>Methods</h2>
<div class="card">
  <pre>{methods_html}</pre>
</div>

<h2>Parameter Decisions Log</h2>
<p style="color: var(--muted); font-size: 0.875rem;">
  All parameter decisions made during this analysis, with justifications.
  Export this section as the basis for your manuscript Methods section.
</p>
<table>
  <tr>
    <th>Checkpoint</th>
    <th>Parameter</th>
    <th>Decision</th>
    <th>Rationale</th>
  </tr>
  {decisions_rows}
</table>

<footer>
  Generated by ARIA (Agentic Research Intelligence for -omics Analysis) &nbsp;|&nbsp;
  github.com/Olascoaga/aria-omics &nbsp;|&nbsp;
  All findings subject to expert review before publication.
</footer>

</body>
</html>"""

        # Build a per-experiment report directory:
        #   ~/.aria/reports/aria_YYYYMMDD_HHMMSS_slug_uuid4/
        #     ├── report.html              (this file)
        #     ├── figures/                 (copied from contrast_dir/figures/)
        #     │   ├── pca_all_samples.svg
        #     │   ├── bmal1_vs_wt/
        #     │   │   ├── volcano.svg
        #     │   │   ├── pathway_dotplot_GO_BP.png
        #     │   │   └── ...
        #     │   └── rev_erba_vs_wt/...
        #     └── tables/                  (copied from contrast_dir/tables/)
        #         ├── bmal1_vs_wt_de_genes.tsv
        #         └── ...
        report_path = report_dir / "report.html"
        report_path.write_text(html, encoding="utf-8")
        methodology_path = report_dir / "methodology.json"
        methodology_path.write_text(
            json.dumps(
                self._build_methodology_json(
                    provenance=provenance,
                    exp_ctx=exp_ctx,
                    agent_results=agent_results,
                    decisions=decisions,
                    llm_usage=llm_usage,
                ),
                indent=2,
                sort_keys=True,
                default=str,
            ),
            encoding="utf-8",
        )
        if reproducible:
            self._write_memory_snapshot(report_dir)
        log.info(f"Report written to {report_path}")
        return report_path

    def _write_memory_snapshot(self, report_dir: Path) -> None:
        try:
            import shutil
            db_path = getattr(self.memory, "db_path", "")
            if db_path and db_path != ":memory:" and Path(db_path).exists():
                shutil.copy2(db_path, report_dir / "memory_snapshot.sqlite")
            else:
                (report_dir / "memory_snapshot.sqlite").write_bytes(b"")
        except Exception as exc:
            log.warning(f"Could not write memory snapshot: {exc}")

    def _build_methodology_json(self, provenance: dict, exp_ctx: dict,
                                agent_results: dict, decisions: list,
                                llm_usage: dict | None = None) -> dict:
        thresholds = {}
        bulk = (agent_results or {}).get("bulk_rna_agent", {})
        bulk_findings = bulk.get("findings", bulk) if isinstance(bulk, dict) else {}
        if isinstance(bulk_findings, dict):
            if bulk_findings.get("padj_threshold") is not None:
                thresholds["bulk_padj"] = bulk_findings.get("padj_threshold")
            if bulk_findings.get("lfc_threshold") is not None:
                thresholds["bulk_lfc_min"] = bulk_findings.get("lfc_threshold")
        sc = (agent_results or {}).get("scrna_agent", {})
        try:
            from aria.agents import _narrative_scrna
            sc_f = _narrative_scrna.unwrap_scrna_findings(sc)
            pb = sc_f.get("pseudobulk_de") or {}
            thresholds["scrna_pseudobulk"] = pb.get("thresholds", {})
            thresholds["scrna_multiple_testing"] = pb.get("multiple_testing", {})
        except Exception:
            pass
        tools = {}
        try:
            from importlib.metadata import version, PackageNotFoundError
            for pkg in ("scanpy", "anndata", "pydeseq2", "gseapy", "numpy", "pandas"):
                try:
                    tools[pkg] = version(pkg)
                except PackageNotFoundError:
                    tools[pkg] = "not installed"
        except Exception:
            tools = {}
        return {
            "provenance": provenance,
            "inputs": exp_ctx.get("input_files", []),
            "design": exp_ctx.get("design", {}),
            "design_intelligence": exp_ctx.get("design_intelligence", {}),
            "thresholds": thresholds,
            "seeds": {
                "global": 0,
                "scanpy": 0,
                "harmony": 0,
            },
            "tools": tools,
            "llm_usage": llm_usage or collect_llm_usage(
                provenance.get("timestamp_utc")
            ),
            "decisions": decisions or [],
        }

    def _build_provenance_section(self, provenance: dict,
                                  input_files: list,
                                  agent_results: dict,
                                  llm_usage: dict | None = None) -> str:
        rows = []
        for key in [
            "aria_version", "git_sha", "git_dirty", "python_version",
            "platform", "conda_env", "timestamp_utc",
        ]:
            rows.append(
                "<tr>"
                f"<td>{_html.escape(key)}</td>"
                f"<td><code>{_html.escape(str(provenance.get(key, '')))}</code></td>"
                "</tr>"
            )
        input_rows = []
        for rec in input_files or []:
            input_rows.append(
                "<tr>"
                f"<td>{_html.escape(str(rec.get('modality', '')))}</td>"
                f"<td><code>{_html.escape(str(rec.get('path', '')))}</code></td>"
                f"<td>{_html.escape(str(rec.get('size_bytes', '')))}</td>"
                f"<td><code>{_html.escape(str(rec.get('sha256', '')))}</code></td>"
                "</tr>"
            )
        if not input_rows:
            input_rows.append(
                "<tr><td colspan='4'><em>No input hashes recorded.</em></td></tr>"
            )
        param_rows = []
        for label, digest in self._collect_param_hashes(agent_results):
            param_rows.append(
                "<tr>"
                f"<td>{_html.escape(str(label))}</td>"
                f"<td><code>{_html.escape(str(digest))}</code></td>"
                "</tr>"
            )
        if not param_rows:
            param_rows.append(
                "<tr><td colspan='2'><em>No per-stage parameter hashes recorded.</em></td></tr>"
            )
        llm_usage = llm_usage or {}
        llm_rows = []
        for key in (
            "calls", "cache_hits", "prompt_tokens", "completion_tokens",
            "total_tokens", "estimated_cost_usd",
        ):
            llm_rows.append(
                "<tr>"
                f"<td>{_html.escape(key)}</td>"
                f"<td><code>{_html.escape(str(llm_usage.get(key, 0)))}</code></td>"
                "</tr>"
            )
        return (
            "<div class='card'>"
            "<h3>Runtime</h3>"
            "<table><tr><th>Field</th><th>Value</th></tr>"
            + "".join(rows)
            + "</table>"
            + "<h3>Inputs</h3>"
            + "<table><tr><th>Modality</th><th>Path</th><th>Bytes</th><th>SHA-256</th></tr>"
            + "".join(input_rows)
            + "</table>"
            + "<h3>Stage Parameter Hashes</h3>"
            + "<table><tr><th>Stage</th><th>params_sha256</th></tr>"
            + "".join(param_rows)
            + "</table>"
            + "<h3>LLM Usage</h3>"
            + "<table><tr><th>Field</th><th>Value</th></tr>"
            + "".join(llm_rows)
            + "</table>"
            + "<h3>Conda Lockfiles</h3>"
            + self._build_lockfile_section()
            + "</div>"
        )

    @staticmethod
    def _collect_param_hashes(obj, prefix: str = "") -> list[tuple[str, str]]:
        rows: list[tuple[str, str]] = []
        if isinstance(obj, dict):
            digest = obj.get("params_sha256")
            if digest:
                rows.append((prefix or "root", str(digest)))
            for key, val in obj.items():
                if key == "params_sha256":
                    continue
                child_prefix = f"{prefix}.{key}" if prefix else str(key)
                rows.extend(NarrativeAgent._collect_param_hashes(val, child_prefix))
        elif isinstance(obj, list):
            for idx, val in enumerate(obj):
                rows.extend(
                    NarrativeAgent._collect_param_hashes(val, f"{prefix}[{idx}]")
                )
        seen = set()
        unique = []
        for label, digest in rows:
            pair = (label, digest)
            if pair not in seen:
                seen.add(pair)
                unique.append(pair)
        return unique

    @staticmethod
    def _build_lockfile_section() -> str:
        root = Path(__file__).resolve().parents[2]
        env_dir = root / "envs"
        # Pair every conda lockfile with its sibling pip lock (if present).
        # `conda create --name X --file <env>.linux-64.lock` reproduces the
        # conda side byte-by-byte; `pip install -r <env>.pip.lock` covers
        # the pip side. Both are needed for full peer-reviewable reproduction.
        lockfiles = sorted(env_dir.glob("*.linux-64.lock"))
        pip_locks = sorted(env_dir.glob("*.pip.lock"))
        if not lockfiles and not pip_locks:
            return (
                "<div class='warning'>No conda lockfiles found in "
                "<code>envs/*.linux-64.lock</code>. Run "
                "<code>scripts/generate_locks.sh</code> before tagging v4.4.</div>"
            )
        blocks = []
        for lock in lockfiles:
            try:
                text = lock.read_text(encoding="utf-8")[:20000]
            except Exception as exc:
                text = f"Could not read lockfile: {exc}"
            blocks.append(
                f"<details><summary>{_html.escape(lock.name)}</summary>"
                f"<pre>{_html.escape(text)}</pre></details>"
            )
        for pip_lock in pip_locks:
            try:
                text = pip_lock.read_text(encoding="utf-8")[:20000]
            except Exception as exc:
                text = f"Could not read pip lockfile: {exc}"
            blocks.append(
                f"<details><summary>{_html.escape(pip_lock.name)} "
                f"(pip side; use after the conda lock)</summary>"
                f"<pre>{_html.escape(text)}</pre></details>"
            )
        # If a conda lock has no sibling pip lock, note that pip was empty
        # (rather than missing) so a reviewer doesn't suspect omission.
        for lock in lockfiles:
            sibling_pip = env_dir / lock.name.replace(".linux-64.lock", ".pip.lock")
            if not sibling_pip.exists():
                blocks.append(
                    f"<div><em>{_html.escape(lock.stem)}: no pip packages "
                    f"in env (pip lock not emitted).</em></div>"
                )
        return "".join(blocks)

    @staticmethod
    def _build_slug(intent: dict, exp_ctx: dict) -> str:
        """
        Build a short URL-safe slug from biological entities or the question.
        Examples:
          {entities: ['GeneA', 'GeneB']} -> "genea_geneb"
          {summary: "Effect of treatment on cells"} -> "effect_treatment_cells"
        """
        import re as _re
        # Prefer biological entities (most specific)
        entities = intent.get("biological_entities", []) or []
        if entities:
            parts = [_re.sub(r"[^a-z0-9]+", "", str(e).lower())
                     for e in entities[:3]]
            parts = [p for p in parts if p and p not in
                     ("cells", "cell", "h9", "h1", "wildtype", "wt")]
            if parts:
                return "_".join(parts)[:40]

        # Fallback: first few words of the summary or user_question
        text = (str(intent.get("summary", "")) or
                str(exp_ctx.get("user_question", ""))).lower()
        if text:
            words = _re.findall(r"[a-z0-9]+", text)
            stop  = {"the", "a", "an", "of", "to", "in", "and", "or", "for",
                     "with", "is", "are", "what", "how", "this", "that",
                     "vs", "versus", "between", "from", "on"}
            keep  = [w for w in words if w not in stop and len(w) > 2][:3]
            if keep:
                return "_".join(keep)[:40]
        return ""

    # ── HTML helpers ──────────────────────────────────────────────────────

    def _build_qc_section(self, grouped: dict,
                           agent_results: dict = None,
                           exp_ctx: dict = None) -> str:
        high  = len(grouped["high"])
        med   = len(grouped["medium"])
        low   = len(grouped["low"])
        ins   = len(grouped["insufficient"])
        total = high + med + low + ins

        # Build per-modality QC rows from agent_results
        qc_rows = ""
        if agent_results:
            rows = []

            # Bulk RNA
            bulk = agent_results.get("bulk_rna_agent", {})
            if bulk.get("status") == "done":
                sqc = bulk.get("findings", {}).get("sample_qc", {}) or {}
                n   = sqc.get("n_samples", "?")
                out = sqc.get("outliers", [])
                ratio = sqc.get("size_ratio", "")
                ratio_str = (f" · library-size range {ratio:.1f}×"
                             if isinstance(ratio, float) else "")
                out_str = (f" · <span style='color:var(--amber)'>"
                           f"{len(out)} outlier(s) removed</span>"
                           if out else " · no outliers removed")
                rows.append(
                    f"<tr><td>Bulk RNA-seq</td>"
                    f"<td>{n} samples{ratio_str}{out_str}</td></tr>"
                )

            # scRNA
            sc = agent_results.get("scrna_agent", {})
            if sc.get("status") == "done":
                from aria.agents import _narrative_scrna as _ns
                sc_f = _ns.unwrap_scrna_findings(sc)
                qc   = sc_f.get("qc", {}) or {}
                n_b  = qc.get("n_cells_before", "?")
                n_a  = qc.get("n_cells_after", "?")
                pct  = qc.get("pct_removed", "")
                pct_str = (f" · {pct:.1f}% removed"
                           if isinstance(pct, (int, float)) else "")
                rows.append(
                    f"<tr><td>scRNA-seq</td>"
                    f"<td>{n_b} → {n_a} cells{pct_str}</td></tr>"
                )

            # Chromatin
            chrom = agent_results.get("chromatin_agent", {})
            if chrom.get("status") == "done":
                chrom_f = chrom.get("findings", {}) or {}
                for assay in ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN"):
                    af    = chrom_f.get(assay, {}).get("findings", {}) or {}
                    peaks = af.get("peaks", {}) or {}
                    if peaks.get("n_peaks"):
                        frip = peaks.get("frip", "?")
                        frip_str = (f"{frip:.2f}"
                                    if isinstance(frip, float) else str(frip))
                        tss  = peaks.get("tss_enrichment", "")
                        tss_str = (f" · TSS enrichment {tss:.1f}"
                                   if isinstance(tss, float) else "")
                        rows.append(
                            f"<tr><td>{assay}</td>"
                            f"<td>{peaks['n_peaks']:,} peaks · "
                            f"FRiP={frip_str}{tss_str}</td></tr>"
                        )

            # HiC
            hic = agent_results.get("genome_arch_agent", {})
            if hic.get("status") == "done":
                hic_f = hic.get("findings", {}) or {}
                bal   = hic_f.get("balancing", {}) or {}
                if bal.get("status"):
                    rows.append(
                        f"<tr><td>Hi-C</td>"
                        f"<td>ICE balancing: {bal.get('status', 'done')}</td></tr>"
                    )

            if rows:
                qc_rows = f"""
<table style="margin-top:1rem">
  <tr><th>Modality</th><th>QC outcome</th></tr>
  {''.join(rows)}
</table>"""

        # Audit findings panel (v4.1)
        audit_html = ""
        audit_findings = (exp_ctx or {}).get("audit_findings", [])
        if audit_findings:
            rows_html = []
            for f in audit_findings:
                sev   = f.get("severity", "warning")
                color = "var(--red)" if sev == "blocking" else "var(--amber)"
                label = "BLOCKING" if sev == "blocking" else "WARNING"
                msg   = (f.get("message", "")
                          .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
                rec   = (f.get("recommendation", "")
                          .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
                rows_html.append(
                    f"<tr>"
                    f"<td style='color:{color};font-weight:600'>{label}</td>"
                    f"<td>{f.get('check','').replace('_',' ').title()}</td>"
                    f"<td>{msg}<br><em style='color:var(--muted);font-size:0.8em'>"
                    f"Recommendation: {rec}</em></td>"
                    f"</tr>"
                )
            audit_html = f"""
<h4 style="margin-top:1.2rem;color:var(--red)">Pre-analysis Quality Audit</h4>
<table style="margin-top:0.5rem">
  <tr><th>Severity</th><th>Check</th><th>Detail</th></tr>
  {''.join(rows_html)}
</table>"""

        if total == 0 and not qc_rows and not audit_html:
            return '<div class="card"><p>No QC data available.</p></div>'

        return f"""
<div class="card">
  <span class="badge high">HIGH {high}</span>
  <span class="badge medium">MEDIUM {med}</span>
  <span class="badge low">LOW {low}</span>
  <span class="badge insuff">INSUFFICIENT {ins}</span>
  <p style="margin-top:0.8rem; color: var(--muted); font-size:0.875rem;">
    {total} findings across all modalities. LOW and INSUFFICIENT findings
    require validation before inclusion in publications.
  </p>
  {qc_rows}
  {audit_html}
</div>"""

    def _build_findings_section(self, sections: dict,
                                  agent_results: dict = None) -> str:
        """Build the findings cards. When agent_results provided, embed plots."""
        parts = []
        section_labels = {
            "bulk_rna":    ("Bulk RNA-seq", "var(--green)"),
            "scrna":       ("Single-cell RNA-seq", "var(--green)"),
            "rna":         ("RNA-seq", "var(--green)"),
            "chromatin":   ("Chromatin", "var(--teal)"),
            "hic":         ("3D Genome", "#a78bfa"),
            "integration": ("Integration", "#f472b6"),
            "synthesis":   ("Integrated Interpretation", "var(--navy)"),
        }
        for key, (label, color) in section_labels.items():
            text = sections.get(key, "")
            if not text:
                continue

            # Build plot embeds if we have bulk RNA results
            plot_html = ""
            if key == "bulk_rna" and agent_results:
                plot_html = self._build_bulk_rna_plots(
                    agent_results.get("bulk_rna_agent", {})
                )
            elif key == "scrna" and agent_results:
                from aria.agents import _narrative_scrna
                sc_envelope = agent_results.get("scrna_agent") or \
                              agent_results.get("rna_agent", {})
                sc_findings = _narrative_scrna.unwrap_scrna_findings(
                    sc_envelope
                )
                plot_html = _narrative_scrna.build_scrna_html_section(
                    sc_findings
                )

            # Escape HTML-breaking newlines into <br> for readability
            body_html = self._plain_text_to_html(text)
            parts.append(f"""
<div class="card">
  <h3 style="color:{color}">{label}</h3>
  <p>{body_html}</p>
  {plot_html}
</div>""")
        return "\n".join(parts) if parts else \
               '<div class="card"><p>No modality findings available.</p></div>'

    @staticmethod
    def _plain_text_to_html(text: str) -> str:
        """Escape plain report text and preserve line breaks."""
        escaped = _html.escape(str(text or ""))
        return escaped.replace("\n", "<br>")

    def _build_methodology_table(self, bulk_result: dict) -> str:
        """
        Render the methodology decisions table.

        Each row documents a pipeline step: what input data was used,
        what normalization was applied, what gene filter, and why.
        This lets reviewers audit the analysis without reading the code.
        """
        findings = bulk_result.get("findings", {})
        methodology = findings.get("methodology", {})
        decisions = methodology.get("decisions", [])
        if not decisions:
            return ""

        rows = []
        for d in decisions:
            rows.append(
                f'<tr>'
                f'<td style="color:var(--navy);font-weight:600">{d.get("step", "")}</td>'
                f'<td>{d.get("input", "")}</td>'
                f'<td>{d.get("normalization", "")}</td>'
                f'<td style="color:var(--muted);font-size:0.9em">{d.get("gene_filter", "")}</td>'
                f'<td style="color:var(--muted);font-size:0.88em;font-style:italic">'
                f'{d.get("justification", "")}</td>'
                f'</tr>'
            )

        return f"""
<details style="margin-top:1.5rem">
  <summary style="cursor:pointer;color:var(--navy);font-weight:600;
                   padding:0.6rem 0;border-bottom:1px solid var(--border)">
    Methods &amp; Decisions (click to expand)
  </summary>
  <div style="margin-top:0.8rem;overflow-x:auto">
    <table style="width:100%;border-collapse:collapse;font-size:0.9em">
      <thead>
        <tr style="border-bottom:2px solid var(--border);color:var(--navy)">
          <th style="text-align:left;padding:0.4rem;width:22%">Step</th>
          <th style="text-align:left;padding:0.4rem;width:15%">Input</th>
          <th style="text-align:left;padding:0.4rem;width:18%">Normalization</th>
          <th style="text-align:left;padding:0.4rem;width:15%">Gene filter</th>
          <th style="text-align:left;padding:0.4rem;width:30%">Justification</th>
        </tr>
      </thead>
      <tbody>
        {''.join(f'<tr style="border-bottom:1px solid var(--border)">{row.replace("<tr>", "").replace("</tr>", "")}</tr>' for row in rows)}
      </tbody>
    </table>
    <p style="font-size:0.82em;color:var(--muted);margin-top:0.6rem;font-style:italic">
      All methodology choices are explicit and auditable.
      Raw counts, VST matrix, and TPM are all exported as supplementary TSV
      tables for downstream analysis.
    </p>
  </div>
</details>
"""

    def _build_bulk_rna_plots(self, bulk_result: dict) -> str:
        """
        Embed all bulk RNA visualizations in the HTML report:
          - Sample PCA (shared across contrasts)
          - Per contrast: volcano, heatmap, ORA dotplots (per database),
            GSEA running sums, GSEA top table
          - Links to supplementary TSV tables in tables/
        """
        findings  = bulk_result.get("findings", {})
        contrasts = findings.get("contrasts", [])
        sqc       = findings.get("sample_qc", {})

        import base64
        from pathlib import Path

        def _embed(path: str) -> str:
            """Inline a file as a data URI (auto-detect SVG/PNG)."""
            try:
                p = Path(path)
                if not p.exists():
                    return ""
                data = p.read_bytes()
                b64  = base64.b64encode(data).decode("ascii")
                if p.suffix.lower() == ".svg":
                    return f"data:image/svg+xml;base64,{b64}"
                if p.suffix.lower() == ".png":
                    return f"data:image/png;base64,{b64}"
                if p.suffix.lower() in (".jpg", ".jpeg"):
                    return f"data:image/jpeg;base64,{b64}"
                return ""
            except Exception:
                return ""

        def _figure(src: str, caption: str, alt: str = "") -> str:
            if not src:
                return ""
            return (
                f'<figure style="margin:1.2rem 0;text-align:center">'
                f'<img src="{src}" alt="{alt or caption}" '
                f'style="max-width:100%;height:auto;border-radius:6px;'
                f'background:var(--bg-alt);border:1px solid var(--border)">'
                f'<figcaption style="color:var(--muted);font-size:0.82rem;'
                f'margin-top:0.4rem;font-style:italic">{caption}</figcaption>'
                f'</figure>'
            )

        html_parts = []

        # ── Shared PCA + MDS (top of bulk section) ─────────────────────
        # Both are VST-based with protein_coding-filtered variable genes (v3.8)
        n_genes_dr     = sqc.get("n_genes_dr", 0)
        n_pc           = sqc.get("n_protein_coding", 0)
        dr_basis       = (f"VST, top {n_genes_dr:,} variable protein_coding genes"
                          if n_pc > 0 else
                          f"VST, top {n_genes_dr:,} most-variable genes")

        pca_path = sqc.get("pca_plot")
        if pca_path:
            src = _embed(pca_path)
            if src:
                html_parts.append(_figure(
                    src,
                    f"Sample PCA — {dr_basis}. "
                    f"Variance-stabilized counts (homoscedastic, library-size corrected).",
                    "Sample PCA",
                ))

        mds_path = sqc.get("mds_plot")
        if mds_path:
            src = _embed(mds_path)
            if src:
                html_parts.append(_figure(
                    src,
                    f"Sample MDS — {dr_basis}, Euclidean distance. "
                    f"Non-linear sample-to-sample distances; complements linear PCA.",
                    "Sample MDS",
                ))

        # ── Methods & Decisions table (v3.8) ─────────────────────────
        # Explicit record of every methodology choice. Collapsible to keep
        # the report focused on results but auditable on demand.
        methods_html = self._build_methodology_table(bulk_result)
        if methods_html:
            html_parts.append(methods_html)

        # ── Per contrast section ─────────────────────────────────────
        for c in contrasts:
            if c.get("status") != "success":
                continue
            cname  = c.get("name", "contrast")
            plots  = c.get("plots", {})

            html_parts.append(
                f'<h4 style="color:var(--navy);margin-top:1.5rem;'
                f'border-bottom:2px solid var(--navy);padding-bottom:0.3rem">'
                f'{cname}</h4>'
            )

            # Volcano (existing)
            src = _embed(plots.get("volcano", ""))
            if src:
                html_parts.append(_figure(src,
                    f"Volcano plot — {cname}. "
                    f"Red: significant DE genes (padj<0.05 & |log2FC| above threshold).",
                    f"Volcano {cname}"))

            # Heatmap by padj (most significant, NEW v3.8)
            src = _embed(plots.get("heatmap_padj", plots.get("heatmap", "")))
            if src:
                html_parts.append(_figure(src,
                    f"Heatmap of top 50 DE genes by padj — {cname}. "
                    f"Input: VST-transformed counts, row z-scored. "
                    f"Rows: HGNC symbols when available.",
                    f"Heatmap padj {cname}"))

            # Heatmap by |log2FC| (largest effect sizes, NEW v3.8)
            src = _embed(plots.get("heatmap_lfc", ""))
            if src:
                html_parts.append(_figure(src,
                    f"Heatmap of top 50 DE genes by |log2FC| — {cname}. "
                    f"Complementary view: largest effect sizes (may differ from padj-ranked).",
                    f"Heatmap lfc {cname}"))

            # ORA dotplots (per database)
            ora = plots.get("ora_dotplots", {}) or {}
            for db, dot_path in ora.items():
                src = _embed(dot_path)
                if src:
                    html_parts.append(_figure(src,
                        f"ORA dotplot ({db}) — top 15 enriched terms. "
                        f"X: log₂(Odds Ratio); color: −log₁₀(FDR); "
                        f"size: gene count.",
                        f"ORA {db} {cname}"))

            # GSEA running sum plots (top 3)
            gsea_runs = plots.get("gsea_running_sums", []) or []
            for i, gpath in enumerate(gsea_runs):
                src = _embed(gpath)
                if src:
                    html_parts.append(_figure(src,
                        f"GSEA running sum #{i+1} (ranked by log₂FC) — {cname}.",
                        f"GSEA running sum {i+1} {cname}"))

            # GSEA top table summary
            tt = _embed(plots.get("gsea_top_table", ""))
            if tt:
                html_parts.append(_figure(tt,
                    f"GSEA top 15 enriched gene sets — {cname}. "
                    f"Sorted by FDR; shows NES distribution.",
                    f"GSEA top table {cname}"))

            # Supplementary table links (relative paths within report dir)
            tables = plots.get("tables", {}) or {}
            de_tsv = tables.get("de_genes")
            pw_tsv = tables.get("pathways")
            link_parts = []
            if de_tsv:
                # Use relative path: report_dir/tables/<contrast>_de_genes.tsv
                rel = f"tables/{Path(de_tsv).parent.parent.name}_de_genes.tsv"
                link_parts.append(
                    f'<a href="{rel}" style="color:var(--blue);'
                    f'text-decoration:underline">DE genes (TSV)</a>'
                )
            if pw_tsv:
                rel = f"tables/{Path(pw_tsv).parent.parent.name}_pathways.tsv"
                link_parts.append(
                    f'<a href="{rel}" style="color:var(--blue);'
                    f'text-decoration:underline">Pathways (TSV)</a>'
                )
            if link_parts:
                html_parts.append(
                    f'<p style="text-align:center;font-size:0.85rem;'
                    f'color:var(--muted);margin:0.5rem 0">'
                    f'Supplementary tables: ' + " &middot; ".join(link_parts)
                    + '</p>'
                )

        return "\n".join(html_parts)

    def _stage_artifacts(self, agent_results: dict, report_dir: Path):
        """
        Copy figures and tables from contrast working dirs into the
        report directory tree, so the report folder is self-contained.

        Source layout:
          <output_dir>/<contrast_slug>/figures/*.svg|*.png
          <output_dir>/<contrast_slug>/tables/*.tsv

        Destination layout:
          report_dir/figures/<contrast_slug>/*.svg|*.png
          report_dir/tables/<contrast_slug>_*.tsv
        """
        import shutil

        bulk = agent_results.get("bulk_rna_agent", {})

        contrasts = []
        if bulk and bulk.get("status") == "done":
            contrasts = bulk.get("findings", {}).get("contrasts", [])
        for c in contrasts:
            if c.get("status") != "success":
                continue

            contrast_dir = c.get("contrast_dir")
            if not contrast_dir:
                continue

            slug = Path(contrast_dir).name

            # Stage figures: report_dir/figures/<slug>/
            src_fig = Path(contrast_dir) / "figures"
            if src_fig.exists():
                dst_fig = report_dir / "figures" / slug
                dst_fig.mkdir(parents=True, exist_ok=True)
                for f in src_fig.iterdir():
                    if f.is_file():
                        try:
                            shutil.copy2(f, dst_fig / f.name)
                        except Exception as e:
                            log.warning(f"Stage figure {f.name}: {e}")

            # Stage tables: report_dir/tables/<slug>_*.tsv
            src_tbl = Path(contrast_dir) / "tables"
            if src_tbl.exists():
                for f in src_tbl.iterdir():
                    if f.is_file():
                        try:
                            shutil.copy2(
                                f,
                                report_dir / "tables" / f"{slug}_{f.name}",
                            )
                        except Exception as e:
                            log.warning(f"Stage table {f.name}: {e}")

        # Also stage shared figures + tables at report_dir level (v3.8)
        import shutil as _sh
        sqc = bulk.get("findings", {}).get("sample_qc", {}) if bulk else {}

        # PCA (existing)
        pca = sqc.get("pca_plot")
        if pca and Path(pca).exists():
            try:
                _sh.copy2(pca, report_dir / "figures" / "pca_all_samples.svg")
            except Exception as e:
                log.warning(f"Stage PCA: {e}")

        # MDS (new v3.8)
        mds = sqc.get("mds_plot")
        if mds and Path(mds).exists():
            try:
                _sh.copy2(mds, report_dir / "figures" / "mds_all_samples.svg")
            except Exception as e:
                log.warning(f"Stage MDS: {e}")

        # Individual PCA + MDS SVGs if they were saved separately
        for fname in ("pca.svg", "mds.svg", "pca_mds.svg"):
            for c in contrasts:
                cdir = c.get("contrast_dir")
                if cdir:
                    candidate = Path(cdir).parent / fname
                    if candidate.exists():
                        try:
                            _sh.copy2(candidate,
                                       report_dir / "figures" / fname)
                        except Exception:
                            pass
                    break

        # TPM supplementary table (new v3.8)
        # rna_bulk_de writes it into the output_dir (parent of contrast_dirs)
        if contrasts:
            cdir = contrasts[0].get("contrast_dir")
            if cdir:
                tpm_src = Path(cdir).parent / "counts_tpm.tsv"
                if tpm_src.exists():
                    try:
                        _sh.copy2(tpm_src,
                                   report_dir / "tables" / "counts_tpm.tsv")
                    except Exception as e:
                        log.warning(f"Stage TPM: {e}")

        # scRNA reports usually carry rich result dictionaries rather than
        # contrast_dir/tables folders. Export those dictionaries here.
        sc = (agent_results or {}).get("scrna_agent") or \
             (agent_results or {}).get("rna_agent") or {}
        if sc.get("status") == "done":
            try:
                from aria.agents import _narrative_scrna
                sc_findings = _narrative_scrna.unwrap_scrna_findings(sc)
                _narrative_scrna.export_supplementary_tables(
                    sc_findings,
                    report_dir / "tables",
                )
            except Exception as e:
                log.warning(f"Stage scRNA supplementary tables: {e}",
                            exc_info=True)


    def _build_findings_table(self, grouped: dict) -> str:
        rows = []
        conf_order = [
            ("high", "HIGH", "high"),
            ("medium", "MEDIUM", "medium"),
            ("low", "LOW", "low"),
            ("insufficient", "INSUFFICIENT", "insuff"),
        ]
        for key, label, css in conf_order:
            for f in grouped[key]:
                summary = self._format_finding_summary(f)
                agent   = _html.escape(str(f.get("agent", "")))
                rows.append(
                    f'<tr><td><span class="badge {css}">{label}</span></td>'
                    f'<td>{summary}</td>'
                    f'<td style="color:var(--muted)">{agent}</td></tr>'
                )
        return "\n".join(rows) if rows else \
               '<tr><td colspan="3">No findings recorded.</td></tr>'

    @staticmethod
    def _format_finding_summary(finding: dict) -> str:
        """Human-readable finding text; no raw Python dicts in reports."""
        summary = finding.get("summary")
        if isinstance(summary, dict):
            if summary.get("error"):
                text = f"Skipped: {summary['error']}"
            else:
                text = "; ".join(
                    f"{k}={v}" for k, v in summary.items()
                    if v is not None
                )
        elif summary:
            text = str(summary)
        elif finding.get("error"):
            text = f"Skipped: {finding['error']}"
        else:
            text = str(finding)

        if "Integration requires at least 2 modalities" in text:
            text = (
                "Integration was skipped because fewer than two modalities "
                "were available for this run."
            )
        text = re.sub(r"\s+", " ", text).strip()
        if len(text) > 200:
            text = text[:197].rstrip() + "..."
        return _html.escape(text)

    @staticmethod
    def _guard_bulk_interpretation(text: str) -> str:
        """Downgrade unsupported causal/mechanistic wording in bulk reports."""
        replacements = {
            r"\bacts as an? ([a-z -]+)\b":
                r"is consistent with a \1-associated expression pattern",
            r"\bsuppressing\b":
                "associated with lower expression of",
            r"\bposition ([^.]+?) as hierarchical gatekeepers\b":
                r"support the hypothesis that \1 influence",
            r"\bhierarchical gatekeepers\b":
                "candidate upstream modulators",
            r"\b([A-Za-z0-9_.-]+) enforces ([A-Za-z0-9_ -]+)\b":
                r"\1 perturbation is associated with \2-related changes",
            r"\b([A-Za-z0-9_.-]+) drives ([A-Za-z0-9_ -]+)\b":
                r"\1 perturbation is associated with \2",
        }
        guarded = text
        for pattern, repl in replacements.items():
            guarded = re.sub(pattern, repl, guarded, flags=re.IGNORECASE)

        causal_terms = (
            "These interpretations are hypotheses from differential "
            "expression and pathway enrichment, not direct evidence of "
            "transcriptional causality; orthogonal perturbation, chromatin, "
            "or binding assays are required to distinguish direct regulation "
            "from secondary state changes."
        )
        if "not direct evidence of transcriptional causality" not in guarded:
            guarded = f"{guarded}\n\n{causal_terms}"
        return guarded

    def _build_decisions_table(self, decisions: list) -> str:
        if not decisions:
            return '<tr><td colspan="4">No decisions recorded.</td></tr>'
        rows = []
        for d in decisions[:30]:
            cp        = _html.escape(str(d.get("checkpoint", "")))
            question  = _html.escape(str(d.get("question", ""))[:120])
            decision  = _html.escape(str(d.get("decision", ""))[:160])
            rationale = _html.escape(str(d.get("rationale", ""))[:300])
            made_by   = _html.escape(str(d.get("made_by", "")))
            made_tag  = (f' <span style="color:var(--dim);font-size:0.85em">'
                          f'({made_by})</span>' if made_by else '')
            rows.append(
                f'<tr>'
                f'<td style="color:var(--muted);white-space:nowrap">CP {cp}</td>'
                f'<td style="color:var(--text)">{question}{made_tag}</td>'
                f'<td style="color:var(--navy);font-weight:600">{decision}</td>'
                f'<td style="color:var(--dim);font-size:0.88em;font-style:italic">'
                f'{rationale}</td>'
                f'</tr>'
            )
        return "\n".join(rows)

    # ── Data helpers ──────────────────────────────────────────────────────

    def _group_findings(self, raw_findings: list) -> dict:
        """Group findings by confidence level."""
        groups = {
            "high":         [],
            "medium":       [],
            "low":          [],
            "insufficient": [],
        }
        for f in raw_findings:
            # findings arrive as dicts from MessageBus payloads
            conf = f.get("confidence", "medium")
            if hasattr(conf, "value"):
                conf = conf.value
            conf = str(conf).lower()
            if conf in groups:
                groups[conf].append(f)
            else:
                groups["medium"].append(f)
        return groups

    def _get_decisions_log(self, experiment_id: str) -> list:
        """Retrieve all parameter decisions for this experiment."""
        try:
            return self.memory.get_decisions(experiment_id)
        except Exception:
            return []

    def _fallback_executive_summary(self, grouped: dict,
                                     intent: dict) -> str:
        """Plain-text fallback when LLM is unavailable."""
        n_high = len(grouped["high"])
        n_total = sum(len(v) for v in grouped.values())
        question = intent.get("summary", "the submitted biological question")
        return (
            f"ARIA completed the analysis for {question}. "
            f"{n_total} findings were recorded across all modalities, "
            f"of which {n_high} reached HIGH confidence. "
            f"See the findings table below for the complete results."
        )

    def receive(self, message):
        pass
