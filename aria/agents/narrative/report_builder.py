"""HTML report rendering / building helpers for :class:`NarrativeAgent`.

Extracted from ``narrative_agent.py`` (P2-8 god-file split, increment 3) as a
behavior-preserving mixin: the methods keep ``self`` and are mixed into
``NarrativeAgent`` through the MRO, so every external call site and every
internal ``self._...`` cross-call is unchanged. This module must not import
``narrative_agent`` (keeps the dependency one-way, no import cycle)."""

from __future__ import annotations

import json
import logging
import html as _html
from datetime import datetime
from pathlib import Path
from typing import Optional

from aria import __version__ as ARIA_VERSION
from aria.utils.provenance import collect_llm_usage, collect_provenance

log = logging.getLogger("aria.narrative")


class ReportBuilderMixin:
    """Report-rendering/build/staging methods mixed into ``NarrativeAgent``."""

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
        narrative_blocks = self._collect_narrative_blocks(agent_results, exp_ctx)

        # P-DEVIL: run the deterministic devil's advocate before rendering so its
        # info caveats appear in the HTML blocks (claim tiers were annotated when
        # the blocks were collected). P-LEDGER: build the planned-vs-run manifest
        # for the provenance section.
        try:
            from aria.agents.narrative.devils_advocate import build_devils_advocate
            devils_advocate = build_devils_advocate(
                narrative_blocks, agent_results, exp_ctx
            )
        except Exception as exc:
            log.warning(f"Devil's-advocate pass failed: {exc}", exc_info=True)
            devils_advocate = []
        try:
            from aria.agents.narrative.run_ledger import build_run_ledger
            run_ledger = build_run_ledger(exp_ctx, agent_results)
        except Exception as exc:
            log.warning(f"Run-ledger build failed: {exc}", exc_info=True)
            run_ledger = {"entries": [], "divergences": [], "n_divergences": 0}

        # W-LEDGER: verify that no associative-or-stronger claim cites a ledger
        # node the run marked not-run/skipped/error. This cross-references TWO
        # structures whose "ran" semantics can legitimately differ (the ledger's
        # finding-based detection vs a narrator's block-creation condition), so it
        # is RECORD-ONLY (fail-safe): a mismatch is recorded in methodology.json
        # and surfaced as a loud caveat in the report, but it never aborts a real
        # report. (Contrast W-CLAIM, which checks a block against its OWN evidence
        # card and is safe to hard-fail in render_blocks.)
        try:
            from aria.agents.narrative.run_ledger import verify_blocks_against_ledger
            ledger_verification = verify_blocks_against_ledger(
                narrative_blocks, run_ledger, strict=False
            )
            if isinstance(run_ledger, dict):
                run_ledger["claim_ledger_verification"] = ledger_verification
            if ledger_verification.get("n_violations"):
                log.warning(
                    "W-LEDGER: %d claim(s) cite a ledger node the run did not "
                    "execute: %s",
                    ledger_verification["n_violations"],
                    [v.get("claim_id") for v in ledger_verification["violations"]],
                )
        except Exception as exc:
            log.warning(f"Ledger claim verification failed: {exc}", exc_info=True)

        # Findings table rows
        findings_rows = self._build_findings_table(grouped_findings)

        # Decisions log rows
        decisions_rows = self._build_decisions_table(decisions)

        version = _html.escape(str(ARIA_VERSION))
        conflicts_html = self._plain_text_to_html(
            findings_sections.get("conflicts", "No conflicts detected.")
        )
        methods_html = self._plain_text_to_html(methods)
        provenance = collect_provenance()
        if isinstance(exp_ctx.get("provenance"), dict):
            provenance.update(exp_ctx["provenance"])
        llm_usage_since = provenance.get("timestamp_utc")
        llm_usage = collect_llm_usage(llm_usage_since)
        if reproducible:
            provenance = dict(provenance)
            provenance["timestamp_utc"] = "<timestamp redacted for byte-identity>"
        git_sha = str(provenance.get("git_sha") or provenance.get("git_commit") or "")
        git_short = _html.escape(git_sha[:12] if git_sha != "unknown" else git_sha)
        git_dirty = _html.escape(str(provenance.get("git_dirty", "")))
        workflow_hash = str(provenance.get("workflow_hash", ""))
        workflow_short = _html.escape(
            workflow_hash[:16] if workflow_hash != "unknown" else workflow_hash
        )
        provenance_html = self._build_provenance_section(
            provenance=provenance,
            input_files=exp_ctx.get("input_files", []),
            agent_results=agent_results,
            llm_usage=llm_usage,
            run_ledger=run_ledger,
        )
        findings_html = self._build_findings_section(
            findings_sections,
            agent_results,
            narrative_blocks=narrative_blocks,
            report_dir=report_dir,
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
  <br>
  Commit: <code>{git_short}</code> &nbsp;|&nbsp;
  Dirty: <code>{git_dirty}</code> &nbsp;|&nbsp;
  Workflow: <code>{workflow_short}</code>
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
{findings_html}

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
                    narrative_blocks=narrative_blocks,
                    run_ledger=run_ledger,
                    devils_advocate=devils_advocate,
                ),
                indent=2,
                sort_keys=True,
                default=str,
            ),
            encoding="utf-8",
        )
        # W-LEDGER: emit a machine-readable RO-Crate (W3C-PROV JSON-LD) next to
        # methodology.json so the run evidence graph is DOI/repository-ready. Pure
        # serialization of what was just written; never blocks report generation.
        try:
            from aria.agents.narrative.ledger_export import write_ro_crate
            write_ro_crate(report_dir)
        except Exception as exc:
            log.warning(f"RO-Crate export failed: {exc}", exc_info=True)
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
                                llm_usage: dict | None = None,
                                narrative_blocks: list | None = None,
                                run_ledger: dict | None = None,
                                devils_advocate: list | None = None) -> dict:
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
        tools = self._collect_tool_versions(
            (
                "scanpy", "anndata", "pydeseq2", "gseapy", "numpy",
                "pandas", "kb-python",
            )
        )
        if narrative_blocks is None:
            narrative_blocks = self._collect_narrative_blocks(agent_results, exp_ctx)
        # X14: per-claim evidence-tier manifests (claim_id -> tier, evidence,
        # limitations, confidence). Best-effort; never block methodology output.
        try:
            from aria.agents.narrative.claim_compiler import compile_claims
            claims = compile_claims(list(narrative_blocks or []), exp_ctx)
        except Exception as exc:
            log.warning(f"Claim manifest compilation failed: {exc}", exc_info=True)
            claims = []
        # P-DEVIL: deterministic devil's advocate over associative+ claims
        # (annotate_claim_tiers ran inside compile_claims, so tiers exist). The
        # build pass is idempotent, so recomputing here when a caller did not
        # pass it is safe.
        if devils_advocate is None:
            try:
                from aria.agents.narrative.devils_advocate import build_devils_advocate
                devils_advocate = build_devils_advocate(
                    list(narrative_blocks or []), agent_results, exp_ctx
                )
            except Exception as exc:
                log.warning(f"Devil's-advocate pass failed: {exc}", exc_info=True)
                devils_advocate = []
        # P-LEDGER: deterministic planned-vs-run manifest.
        if run_ledger is None:
            try:
                from aria.agents.narrative.run_ledger import build_run_ledger
                run_ledger = build_run_ledger(exp_ctx, agent_results)
            except Exception as exc:
                log.warning(f"Run-ledger build failed: {exc}", exc_info=True)
                run_ledger = {"entries": [], "divergences": [],
                              "n_divergences": 0}
        # W-LEDGER: link every compiled claim to its ledger node so each claim is
        # traceable to both an evidence card (W-CLAIM) and a run-ledger node. The
        # per-claim ledger_node_id/ledger_status are written in place; the summary
        # (incl. any contradiction) is recorded on the ledger manifest.
        try:
            from aria.agents.narrative.run_ledger import link_claims_to_ledger
            run_ledger["claim_linkage"] = link_claims_to_ledger(claims, run_ledger)
        except Exception as exc:
            log.warning(f"Claim-ledger linkage failed: {exc}", exc_info=True)
        try:
            from aria.agents.narrative.robustness import build_robustness_multiverse
            robustness_multiverse = build_robustness_multiverse(agent_results)
        except Exception as exc:
            log.warning(f"Robustness multiverse build failed: {exc}", exc_info=True)
            robustness_multiverse = {"status": "error", "details": str(exc)}
        return {
            "provenance": provenance,
            "inputs": exp_ctx.get("input_files", []),
            "raw_ingestion": exp_ctx.get("raw_ingestion", []),
            "narrative_blocks": [
                block.to_dict() for block in narrative_blocks or []
            ],
            "claims": claims,
            "devils_advocate": devils_advocate,
            "run_ledger": run_ledger,
            "robustness_multiverse": robustness_multiverse,
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
                                  llm_usage: dict | None = None,
                                  run_ledger: dict | None = None) -> str:
        rows = []
        for key in [
            "aria_version", "version_source", "git_sha", "git_commit",
            "git_dirty", "git_tree_sha", "git_describe", "workflow_hash",
            "workflow_hash_algorithm", "python_version", "platform",
            "conda_env", "timestamp_utc",
        ]:
            rows.append(
                "<tr>"
                f"<td>{_html.escape(key)}</td>"
                f"<td><code>{_html.escape(str(provenance.get(key, '')))}</code></td>"
                "</tr>"
            )
        # P2-2: cite the container image identity (digest) the report ran in.
        # `image` is a nested dict, so render it explicitly; when ARIA is not
        # running in a pinned image, say so honestly rather than omit it.
        image = provenance.get("image") or {}
        if isinstance(image, dict) and image.get("containerized"):
            for ikey, label in (
                ("kind", "image_kind"),
                ("digest", "image_digest"),
                ("reference", "image_reference"),
                ("revision", "image_revision"),
                ("env_lock_sha256", "image_env_lock_sha256"),
                ("validation", "image_validation"),
            ):
                val = image.get(ikey)
                if val:
                    rows.append(
                        "<tr>"
                        f"<td>{_html.escape(label)}</td>"
                        f"<td><code>{_html.escape(str(val))}</code></td>"
                        "</tr>"
                    )
        else:
            rows.append(
                "<tr><td>image</td>"
                "<td><code>not containerized</code></td></tr>"
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
            "total_tokens", "estimated_cost_usd", "deterministic",
            "degraded", "fallback_calls",
            "temperature", "seed", "models", "tiers",
        ):
            llm_rows.append(
                "<tr>"
                f"<td>{_html.escape(key)}</td>"
                f"<td><code>{_html.escape(str(llm_usage.get(key, 0)))}</code></td>"
                "</tr>"
            )
        ingestion_html = self._build_raw_ingestion_section(
            agent_results=agent_results,
            exp_ctx_records=[]
        )
        ledger_html = self._build_run_ledger_section(run_ledger)
        # W-CALIB: numerical-calibration badge. `provenance["calibration"]` is a
        # build property, present only when a real calibration run was attached;
        # absent on a normal report -> the badge says "not measured" honestly.
        calibration_html = self._build_calibration_badge(provenance.get("calibration"))
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
            + ingestion_html
            + "<h3>LLM Usage</h3>"
            + "<table><tr><th>Field</th><th>Value</th></tr>"
            + "".join(llm_rows)
            + "</table>"
            + ledger_html
            + calibration_html
            + "<h3>Conda Lockfiles</h3>"
            + self._build_lockfile_section()
            + "</div>"
        )






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
                out = sqc.get("candidate_outliers", sqc.get("outliers", []))
                sens_out = sqc.get("sensitivity_outliers_removed", [])
                ratio = sqc.get("size_ratio", "")
                ratio_str = (f" · library-size range {ratio:.1f}×"
                             if isinstance(ratio, float) else "")
                out_str = (f" · <span style='color:var(--amber)'>"
                           f"{len(out)} outlier(s) flagged; "
                           f"{len(sens_out)} removed in sensitivity</span>"
                           if out else " · no outliers flagged")
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
                                  agent_results: dict = None,
                                  narrative_blocks: list | None = None,
                                  report_dir: Path | None = None) -> str:
        """Build the findings cards. When agent_results provided, embed plots."""
        parts = []
        block_groups = {}
        if narrative_blocks:
            from aria.agents.narrative.render_blocks import (
                group_blocks_by_prefix,
                render_blocks,
            )
            block_groups = group_blocks_by_prefix(narrative_blocks)
        section_labels = {
            "bulk_rna":    ("Bulk RNA-seq", "var(--green)"),
            "scrna":       ("Single-cell RNA-seq", "var(--green)"),
            "rna":         ("RNA-seq", "var(--green)"),
            "chromatin":   ("Chromatin", "var(--teal)"),
            "hic":         ("3D Genome", "#a78bfa"),
            "integration": ("Integration", "#f472b6"),
            "synthesis":   ("Integrated Interpretation", "var(--navy)"),
        }
        section_prefix = {
            "bulk_rna": "bulk",
            "scrna": "scrna",
        }
        for key, (label, color) in section_labels.items():
            text = sections.get(key, "")
            blocks = block_groups.get(section_prefix.get(key, ""), [])
            if not text and not blocks:
                continue

            # Build plot embeds if we have bulk RNA results
            plot_html = ""
            if blocks:
                plot_html = render_blocks(blocks, report_dir=report_dir)
                body_html = ""
            elif key == "bulk_rna" and agent_results:
                plot_html = self._build_bulk_rna_plots(
                    agent_results.get("bulk_rna_agent", {})
                )
                body_html = self._plain_text_to_html(text)
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
                body_html = self._plain_text_to_html(text)
            else:
                body_html = self._plain_text_to_html(text)

            # Escape HTML-breaking newlines into <br> for readability
            parts.append(f"""
<div class="card">
  <h3 style="color:{color}">{label}</h3>
  {f'<p>{body_html}</p>' if body_html else ''}
  {plot_html}
</div>""")
        return "\n".join(parts) if parts else \
               '<div class="card"><p>No modality findings available.</p></div>'


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
        import html as _html
        rows = []
        conf_order = [
            ("high", "HIGH", "high"),
            ("medium", "MEDIUM", "medium"),
            ("low", "LOW", "low"),
            ("insufficient", "INSUFFICIENT", "insuff"),
        ]
        for key, label, css in conf_order:
            for f in grouped[key]:
                summary = _html.escape(self._format_finding_summary(f))
                agent   = _html.escape(str(f.get("agent", "")))
                rows.append(
                    f'<tr><td><span class="badge {css}">{label}</span></td>'
                    f'<td>{summary}</td>'
                    f'<td style="color:var(--muted)">{agent}</td></tr>'
                )
        return "\n".join(rows) if rows else \
               '<tr><td colspan="3">No findings recorded.</td></tr>'



    def _build_decisions_table(self, decisions: list) -> str:
        import html as _html
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
