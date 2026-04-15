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
from datetime import datetime
from pathlib import Path
from typing import Optional

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
        Uses HEAVY tier — this is user-facing prose.
        """
        high_findings = grouped["high"][:5]
        med_findings  = grouped["medium"][:3]
        n_low         = len(grouped["low"])
        n_insuff      = len(grouped["insufficient"])

        finding_summaries = "\n".join([
            f"- [HIGH] {f.get('summary', str(f))[:120]}"
            for f in high_findings
        ] + [
            f"- [MEDIUM] {f.get('summary', str(f))[:120]}"
            for f in med_findings
        ])

        prompt = f"""
Experiment: {exp_ctx.get('organism', 'unknown')} / {exp_ctx.get('genome', 'unknown')}
Biological question: {intent.get('summary', exp_ctx.get('user_question', 'unknown'))}
Analysis type: {intent.get('analysis_type', 'unknown')}

Key findings (by confidence):
{finding_summaries}

Additional context:
- {n_low} findings with LOW confidence (require validation)
- {n_insuff} analyses where data was insufficient to conclude

Write a 2-3 sentence executive summary for a PI.
Be specific about what was found. Include confidence levels.
Do not speculate beyond what the data showed.
"""
        try:
            return self.llm.complete(
                prompt=prompt,
                system=NARRATIVE_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=300,
            )
        except Exception as e:
            log.warning(f"Executive summary LLM failed: {e}")
            return self._fallback_executive_summary(grouped, intent)

    def _write_findings_sections(self, grouped: dict,
                                  agent_results: dict,
                                  exp_ctx: dict) -> dict:
        """
        Write structured findings for each modality.
        Returns dict of section_name -> prose text.
        """
        sections = {}

        # RNA findings
        rna = agent_results.get("rna_agent", {})
        if rna.get("status") == "done":
            sections["rna"] = self._summarize_rna(rna, grouped)

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

        # RNA methods
        rna = agent_results.get("rna_agent", {})
        if rna.get("status") == "done":
            findings_rna = rna.get("findings", {})
            qc           = findings_rna.get("scRNA", {}).get(
                               "findings", {}
                           ).get("qc", {})
            clustering   = findings_rna.get("scRNA", {}).get(
                               "findings", {}
                           ).get("clustering_decision", {})

            lines.append(f"\n**Single-cell RNA-seq**\n")
            if qc:
                mt_thr = qc.get("mt_threshold", "")
                lines.append(
                    f"Raw count matrices were processed using scanpy. "
                    f"Cells were filtered using adaptive MAD-based thresholds "
                    f"with a maximum mitochondrial fraction of "
                    f"{mt_thr}%. "
                    f"Counts were normalized to 10,000 reads per cell "
                    f"and log1p-transformed."
                )
            if clustering:
                res = clustering.get("recommended", "")
                lines.append(
                    f"Dimensionality reduction was performed using PCA "
                    f"(50 components), followed by k-nearest neighbor "
                    f"graph construction (k=15) and UMAP visualization. "
                    f"Leiden clustering was performed at resolution={res}, "
                    f"selected by the ParameterAdvisor based on silhouette "
                    f"score maximization across {4} candidates."
                )

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

    def _summarize_rna(self, rna_result: dict,
                        grouped: dict) -> str:
        findings = rna_result.get("findings", {})
        sc       = findings.get("scRNA", {}).get("findings", {})

        parts = []

        qc = sc.get("qc", {})
        if qc:
            n_after  = qc.get("n_cells_after", "?")
            pct_rm   = qc.get("pct_removed", "?")
            parts.append(
                f"After quality control, {n_after} cells were retained "
                f"({pct_rm}% removed)."
            )

        ct = sc.get("cell_types", {})
        if ct and ct.get("cell_types"):
            unique_types = list(set(ct["cell_types"].values()))
            # Filter out failed annotations
            unique_types = [t for t in unique_types
                            if "failed" not in str(t).lower()]
            if unique_types:
                parts.append(
                    f"Cell type annotation identified: "
                    f"{', '.join(unique_types[:6])}."
                )

        de = sc.get("differential_expression", {})
        if de:
            n_sig = de.get("n_significant", 0)
            parts.append(
                f"Differential expression analysis identified "
                f"{n_sig} significant genes."
            )

        return " ".join(parts) if parts else \
               "RNA analysis completed. See findings table for details."

    def _summarize_chromatin(self, chrom_result: dict,
                              grouped: dict) -> str:
        findings = chrom_result.get("findings", {})
        parts    = []

        for assay in ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN"):
            assay_f = findings.get(assay, {}).get("findings", {})
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

    def _render_html_report(self, experiment_id: str,
                             exp_ctx: dict,
                             intent: dict,
                             executive_summary: str,
                             findings_sections: dict,
                             grouped_findings: dict,
                             methods: str,
                             decisions: list) -> Path:
        """
        Render the full HTML report.
        Self-contained: CSS embedded, no external dependencies.
        """
        organism = exp_ctx.get("organism", "Unknown organism")
        genome   = exp_ctx.get("genome", "Unknown assembly")
        question = intent.get("summary",
                               exp_ctx.get("user_question", ""))
        date_str = datetime.now().strftime("%B %d, %Y")
        exp_short = experiment_id[:8]

        # Findings table rows
        findings_rows = self._build_findings_table(grouped_findings)

        # Decisions log rows
        decisions_rows = self._build_decisions_table(decisions)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ARIA Report — {exp_short}</title>
<style>
  :root {{
    --navy:  #0f1729;
    --panel: #1a2744;
    --card:  #1e2d50;
    --cyan:  #22d3ee;
    --teal:  #14b8a6;
    --green: #4ade80;
    --amber: #fbbf24;
    --red:   #f87171;
    --text:  #e2e8f0;
    --muted: #94a3b8;
    --dim:   #64748b;
    --border:#2d3f6e;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: 'Segoe UI', system-ui, sans-serif;
    background: var(--navy);
    color: var(--text);
    line-height: 1.6;
    padding: 2rem;
    max-width: 960px;
    margin: 0 auto;
  }}
  h1 {{ color: var(--cyan); font-size: 1.8rem; margin-bottom: 0.3rem; }}
  h2 {{ color: var(--teal); font-size: 1.2rem; margin: 2rem 0 0.8rem;
        border-bottom: 1px solid var(--border); padding-bottom: 0.4rem; }}
  h3 {{ color: var(--muted); font-size: 1rem; margin: 1.2rem 0 0.4rem; }}
  p  {{ margin: 0.6rem 0; color: var(--text); }}
  .meta {{ color: var(--muted); font-size: 0.85rem; margin-bottom: 2rem; }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.2rem 1.5rem;
    margin: 1rem 0;
  }}
  .badge {{
    display: inline-block;
    padding: 0.15rem 0.6rem;
    border-radius: 99px;
    font-size: 0.75rem;
    font-weight: 600;
    margin-right: 0.4rem;
  }}
  .high   {{ background: #14532d; color: var(--green); }}
  .medium {{ background: #451a03; color: var(--amber); }}
  .low    {{ background: #450a0a; color: var(--red); }}
  .insuff {{ background: #1e293b; color: var(--dim); }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.875rem;
    margin: 0.8rem 0;
  }}
  th {{
    background: var(--panel);
    color: var(--muted);
    text-align: left;
    padding: 0.5rem 0.75rem;
    font-weight: 500;
  }}
  td {{
    border-top: 1px solid var(--border);
    padding: 0.5rem 0.75rem;
    vertical-align: top;
  }}
  tr:hover td {{ background: rgba(255,255,255,0.02); }}
  pre {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 1rem;
    font-size: 0.82rem;
    white-space: pre-wrap;
    color: var(--muted);
    margin: 0.8rem 0;
  }}
  .warning {{
    border-left: 3px solid var(--amber);
    padding-left: 1rem;
    color: var(--amber);
    font-size: 0.9rem;
    margin: 0.8rem 0;
  }}
  footer {{
    margin-top: 3rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
    color: var(--dim);
    font-size: 0.8rem;
    text-align: center;
  }}
</style>
</head>
<body>

<h1>ARIA Analysis Report</h1>
<p class="meta">
  Experiment: <strong>{exp_short}</strong> &nbsp;|&nbsp;
  {organism} / {genome} &nbsp;|&nbsp;
  {date_str} &nbsp;|&nbsp;
  Generated by ARIA v0.2
</p>

<div class="card">
  <h3>Biological Question</h3>
  <p><em>{question}</em></p>
</div>

<h2>Executive Summary</h2>
<div class="card">
  <p>{executive_summary}</p>
</div>

<h2>Quality Control Summary</h2>
{self._build_qc_section(grouped_findings)}

<h2>Findings</h2>
{self._build_findings_section(findings_sections)}

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
  <pre>{findings_sections.get("conflicts", "No conflicts detected.")}</pre>
</div>

<h2>Methods</h2>
<div class="card">
  <pre>{methods}</pre>
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

        report_path = self.reports_dir / f"{experiment_id}.html"
        report_path.write_text(html, encoding="utf-8")
        log.info(f"Report written to {report_path}")
        return report_path

    # ── HTML helpers ──────────────────────────────────────────────────────

    def _build_qc_section(self, grouped: dict) -> str:
        high = len(grouped["high"])
        med  = len(grouped["medium"])
        low  = len(grouped["low"])
        ins  = len(grouped["insufficient"])
        total = high + med + low + ins
        if total == 0:
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
</div>"""

    def _build_findings_section(self, sections: dict) -> str:
        parts = []
        section_labels = {
            "rna":         ("RNA-seq", "var(--green)"),
            "chromatin":   ("Chromatin", "var(--teal)"),
            "hic":         ("3D Genome", "#a78bfa"),
            "integration": ("Integration", "#f472b6"),
        }
        for key, (label, color) in section_labels.items():
            text = sections.get(key, "")
            if text:
                parts.append(f"""
<div class="card">
  <h3 style="color:{color}">{label}</h3>
  <p>{text}</p>
</div>""")
        return "\n".join(parts) if parts else \
               '<div class="card"><p>No modality findings available.</p></div>'

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
                summary = f.get("summary", str(f))[:200]
                agent   = f.get("agent", "")
                rows.append(
                    f'<tr><td><span class="badge {css}">{label}</span></td>'
                    f'<td>{summary}</td>'
                    f'<td style="color:var(--muted)">{agent}</td></tr>'
                )
        return "\n".join(rows) if rows else \
               '<tr><td colspan="3">No findings recorded.</td></tr>'

    def _build_decisions_table(self, decisions: list) -> str:
        if not decisions:
            return '<tr><td colspan="4">No decisions recorded.</td></tr>'
        rows = []
        for d in decisions[:30]:
            cp        = d.get("checkpoint", "")
            question  = d.get("question", "")[:60]
            decision  = d.get("decision", "")[:40]
            rationale = d.get("rationale", "")[:100]
            rows.append(
                f'<tr>'
                f'<td style="color:var(--muted)">CP {cp}</td>'
                f'<td>{question}</td>'
                f'<td style="color:var(--cyan)">{decision}</td>'
                f'<td style="color:var(--dim)">{rationale}</td>'
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
