"""Report rendering, figure orchestration, and speculative-section wiring."""

from __future__ import annotations

import math

from aria.agents.narrative.reporting._base import *


def _sanitize_json_value(obj):
    """Recursively replace non-finite floats with None for strict RFC-8259 JSON."""
    if isinstance(obj, bool) or obj is None:
        return obj
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _sanitize_json_value(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize_json_value(v) for v in obj]
    if isinstance(obj, tuple):
        return [_sanitize_json_value(v) for v in obj]
    if isinstance(obj, set):
        return [_sanitize_json_value(v) for v in sorted(obj, key=str)]
    # Handle NumPy scalar types cleanly without requiring numpy import.
    item = getattr(obj, "item", None)
    if callable(item) and type(obj).__module__.startswith("numpy"):
        try:
            return _sanitize_json_value(item())
        except Exception:  # aria: allow-broad-except cosmetic numpy scalar conversion fallback
            pass
    return obj


class ReportRenderMixin:
    def _collect_execution_llm_usage(self, since_utc: str | None = None) -> dict:
        provider = getattr(self, "llm", None)
        return collect_llm_usage(
            since_utc,
            experiment_id=getattr(provider, "experiment_id", None),
            usage_log=getattr(provider, "usage_log", None),
        )

    def _build_report_dir(self, experiment_id: str,
                           intent: dict, exp_ctx: dict) -> Path:
        """
        Construct (and mkdir) ~/.aria/reports/aria_<ts>_<slug>_<run>/
        with figures/ + tables/ subdirs. Call exactly once early in run() so
        figure subprocesses have a stable, execution-owned target.
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
            content_suffix = first_hash or "nohash"
        else:
            content_suffix = ""

        # A6: experiment_id is the execution identity (entrypoints allocate a
        # fresh UUID-derived id per run).  Keep it in every directory name so
        # two runs over identical inputs cannot resolve to the same artifact
        # bundle.  The short hash protects uniqueness if sanitization collapses
        # unusual caller-supplied ids.  Keep the legacy final-id suffix last for
        # the history/headless fallback lookup contract.
        raw_run_id = str(experiment_id or "run")
        safe_run_id = re.sub(r"[^A-Za-z0-9_-]+", "-", raw_run_id).strip("-_")
        safe_run_id = (safe_run_id or "run")[-48:]
        lookup_suffix = re.sub(
            r"[^A-Za-z0-9_-]+", "-", raw_run_id[-4:]
        ) or "run"
        run_token = (
            f"{hashlib.sha256(raw_run_id.encode('utf-8')).hexdigest()[:8]}_"
            f"{safe_run_id}_{lookup_suffix}"
        )
        suffix = f"{content_suffix}_{run_token}" if content_suffix else run_token
        report_name = (f"aria_{ts}_{slug}_{suffix}" if slug
                       else f"aria_{ts}_{suffix}")
        report_dir = self.reports_dir / report_name
        # Publish a new run identity exactly once.  A collision fails closed;
        # never delete or reuse a previous run's evidence bundle.
        report_dir.mkdir(parents=True, exist_ok=False)
        (report_dir / "figures").mkdir()
        (report_dir / "tables").mkdir()
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
        except Exception as e:  # aria: allow-broad-except record scRNA figure import error for W-WARN
            log.warning(f"Cannot import figure helpers: {e}")
            from aria.agents.narrative.protocols import record_figure_error
            if isinstance(sc.get("findings"), dict):
                record_figure_error(sc["findings"], "scrna_figures_import", e)
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
        except Exception as e:  # aria: allow-broad-except record scRNA figure generation error for W-WARN
            log.warning(f"scRNA figure generation failed: {e}", exc_info=True)
            from aria.agents.narrative.protocols import record_figure_error
            record_figure_error(findings, "scrna_figures", e)

    def _generate_chromatin_figures(self, agent_results: dict,
                                    report_dir: Path) -> None:
        """W0.1 (scATAC P0): render scATAC figures in the chromatin stack and
        attach their paths so ``ChromatinNarrator.figures()`` surfaces them.

        Mirrors ``_generate_scrna_figures``. Mutates the LIVE per-modality findings
        dict (``_narrative_chromatin.live_findings``) so the narrator sees the
        figures on its re-unwrap. No-op when there is no chromatin run, no
        clustered ``.h5ad``, or the env manager is unavailable (honest absence).
        """
        ch = (agent_results or {}).get("chromatin_agent") or {}
        try:
            from aria.agents import _narrative_chromatin
            from aria.agents.narrative.narrators.chromatin import (
                unwrap_chromatin_findings,
            )
            from aria.utils.environment_manager import env_manager
        except Exception as e:  # aria: allow-broad-except record chromatin figure import error for W-WARN
            log.warning(f"Cannot import chromatin figure helpers: {e}")
            from aria.agents.narrative.protocols import record_figure_error
            if isinstance(ch.get("findings"), dict):
                record_figure_error(ch["findings"], "chromatin_figures_import", e)
            return
        view = unwrap_chromatin_findings(ch)
        if not view:
            return
        lsi = view.get("lsi") or view.get("lsi_clustering") or {}
        h5ad_path = lsi.get("output_path")
        live = _narrative_chromatin.live_findings(ch)
        if live is None:
            return
        try:
            _narrative_chromatin.generate_figures(
                live,
                h5ad_path,
                Path(report_dir) / "figures",
                env_manager=env_manager,
            )
        except Exception as e:  # aria: allow-broad-except record chromatin figure generation error for W-WARN
            log.warning(f"chromatin figure generation failed: {e}", exc_info=True)
            from aria.agents.narrative.protocols import record_figure_error
            record_figure_error(live, "chromatin_figures", e)

    def _generate_integration_figures(self, agent_results: dict,
                                      report_dir: Path) -> None:
        """V48 Slice I7: render dual PNG+SVG integration figures inline and attach
        their paths onto ``integration_agent.findings['figures']`` so
        ``IntegrationNarrator`` embeds them into governed narrative blocks."""
        integ = (agent_results or {}).get("integration_agent") or {}
        if integ.get("status") != "done":
            return
        try:
            from aria.agents import _narrative_integration
            from aria.agents.narrative.narrators.integration import (
                unwrap_integration_findings,
            )
        except Exception as e:  # aria: allow-broad-except record integration figure import error for W-WARN
            log.warning(f"Cannot import integration figure helpers: {e}")
            from aria.agents.narrative.protocols import record_figure_error
            if isinstance(integ.get("findings"), dict):
                record_figure_error(integ["findings"], "integration_figures_import", e)
            return
        findings = unwrap_integration_findings(integ)
        if not findings:
            return
        try:
            _narrative_integration.generate_figures(
                findings=findings,
                output_dir=Path(report_dir) / "figures",
            )
        except Exception as e:  # aria: allow-broad-except record integration figure generation error for W-WARN
            log.warning(
                f"integration figure generation failed: {e}", exc_info=True
            )
            from aria.agents.narrative.protocols import record_figure_error
            record_figure_error(findings, "integration_figures", e)

    def _speculative_verification_state(
        self,
        run_ledger: dict | None,
        narrative_blocks: list | None,
        n_withheld: int = 0,
    ) -> "VerificationReceipt":
        """Resolve the REAL W-CLAIM / W-LEDGER state for the causal gate (rail #1).

        ADR-057 rail #1 says the agent only speculates downstream of W-CLAIM +
        W-LEDGER PASSING. The old tuple form synthesised ``(True, True)`` from the
        ABSENCE of the verification artifacts (``not any([])`` is True; a missing
        ledger record defaults ``n_violations`` to 0) — i.e. absence was read as
        approval, making the guarantee decorative exactly in the dangerous case.

        Round-3 H14: this returns a fail-closed ``VerificationReceipt`` that
        requires POSITIVE evidence on each side; absence yields ``complete=False``
        and the gate stays shut.

        - W-CLAIM: requires rendered blocks to exist. With blocks, it passes iff
          none recorded an ``unsupported`` W-CLAIM verification (strict W-CLAIM
          hard-fails at render, so a surviving block is normally supported). With
          ``narrative_blocks is None`` there is no evidence — incomplete, shut.
        - W-LEDGER: requires the ``claim_ledger_verification`` record (written
          unconditionally on a normal run, see the W-LEDGER block above). Present,
          it passes iff ``n_violations == 0``. Absent, W-LEDGER did not complete —
          incomplete, shut. W-LEDGER is record-only (never aborts the report), so
          this is the only place its finding actually gates anything.
        """
        from aria.agents.narrative.hypothesis import VerificationReceipt

        withheld_count = max(0, int(n_withheld or 0))
        if narrative_blocks is None:
            w_claim_passed = False
            w_claim_complete = False
            w_claim_evidence = "absent: no rendered narrative blocks"
        else:
            n_unsupported = sum(
                1
                for b in narrative_blocks
                if (getattr(b, "metadata", {}) or {})
                .get("claim_verification", {})
                .get("status")
                == "unsupported"
            )
            w_claim_passed = n_unsupported == 0 and withheld_count == 0
            w_claim_complete = True
            w_claim_evidence = (
                f"{len(narrative_blocks)} block(s), {n_unsupported} unsupported, "
                f"{withheld_count} withheld"
            )

        has_ledger_record = (
            isinstance(run_ledger, dict)
            and "claim_ledger_verification" in run_ledger
        )
        if not has_ledger_record:
            w_ledger_passed = False
            w_ledger_complete = False
            w_ledger_evidence = "absent: no claim_ledger_verification record"
        else:
            n_violations = int(
                (run_ledger["claim_ledger_verification"] or {}).get(
                    "n_violations", 0
                )
                or 0
            )
            w_ledger_passed = n_violations == 0
            w_ledger_complete = True
            w_ledger_evidence = f"n_violations={n_violations}"

        return VerificationReceipt(
            w_claim_passed=w_claim_passed,
            w_ledger_passed=w_ledger_passed,
            complete=w_claim_complete and w_ledger_complete,
            w_claim_evidence=w_claim_evidence,
            w_ledger_evidence=w_ledger_evidence,
        )

    def _build_speculative_section_html(self, agent_results: dict,
                                        run_ledger: dict | None,
                                        exp_ctx: dict,
                                        narrative_blocks: list | None = None,
                                        report_dir: Optional[Path] = None,
                                        n_withheld: int = 0,
                                        ) -> str:
        """ADR-057 S9: render the opt-in SPECULATIVE hypotheses section, or ''.

        Off unless ``exp_ctx['enable_hypotheses']`` is set (opt-in, never
        inferred). Uses the agent's LLM as the proposer; the HypothesisAgent's
        gates + quarantine wall the output. The causal gate (rail #1) is fed the
        run's REAL W-CLAIM/W-LEDGER state, not an unconditional pass. The
        structured result is persisted to a SEPARATE, non-promotable
        ``speculative_hypotheses.json`` so the speculative layer is auditable
        outside the HTML. Fully guarded — a failure here never breaks report
        generation.
        """
        try:
            if not (exp_ctx or {}).get("enable_hypotheses"):
                return ""
            from aria.agents.narrative.hypothesis import (
                LLMProposer,
                build_speculative_section,
                persist_speculative_manifest,
                render_speculative_section_html,
            )
            proposer = None
            llm = getattr(self, "llm", None)
            if llm is not None:
                proposer = LLMProposer.from_provider(llm)
            verification = self._speculative_verification_state(
                run_ledger, narrative_blocks, n_withheld=n_withheld
            )
            section = build_speculative_section(
                agent_results, run_ledger, exp_ctx, proposer=proposer,
                verification=verification,
            )
            manifest_warn = ""
            if report_dir is not None and section is not None:
                # Auditable, non-promotable manifest — separate from the audited
                # claim manifest. A persist failure must not lose the HTML.
                try:
                    persist_speculative_manifest(
                        section, report_dir,
                        reproducible=bool(exp_ctx.get("reproducible_mode")),
                    )
                except Exception as exc:  # aria: allow-broad-except surface speculative manifest persistence failure
                    log.warning(
                        "Persisting speculative manifest failed: %s",
                        exc, exc_info=True,
                    )
                    manifest_warn = (
                        "<div class='warning'><strong>Speculative manifest persistence failed:</strong> "
                        f"{_html.escape(f'{type(exc).__name__}: {exc}')}</div>"
                    )
            return render_speculative_section_html(section) + manifest_warn
        except Exception as exc:  # aria: allow-broad-except surface speculative hypotheses error in HTML and manifest
            log.warning(
                f"Speculative hypotheses section failed: {exc}", exc_info=True
            )
            err_msg = f"{type(exc).__name__}: {exc}"
            if report_dir is not None:
                try:
                    err_manifest = {
                        "schema_version": "aria.speculative_hypotheses.v1",
                        "status": "error",
                        "promotable": False,
                        "error": err_msg,
                        "hypotheses": [],
                    }
                    (Path(report_dir) / "speculative_hypotheses.json").write_text(
                        json.dumps(err_manifest, indent=2, sort_keys=True),
                        encoding="utf-8",
                    )
                except Exception as write_exc:  # aria: allow-broad-except best-effort speculative error manifest write
                    log.warning(
                        "Could not write speculative error manifest: %s",
                        write_exc,
                    )
            return (
                "<div class='card warning'>"
                "<h3>Speculative Hypotheses (Error)</h3>"
                f"<p>Speculative hypothesis generation failed "
                f"(<code>{_html.escape(err_msg)}</code>); "
                "no hypotheses were generated or promoted.</p>"
                "</div>"
            )

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
        organism = _html.escape(str(exp_ctx.get("organism", "Unknown organism")))
        genome = _html.escape(str(exp_ctx.get("genome", "Unknown assembly")))
        question = exp_ctx.get("user_question") or "Submitted question unavailable"
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
        (report_dir / "figures").mkdir(exist_ok=True)
        (report_dir / "tables").mkdir(exist_ok=True)
        self._stage_artifacts(agent_results, report_dir)
        integ_env = (agent_results or {}).get("integration_agent") or {}
        if integ_env.get("status") == "done":
            integ_f = integ_env.get("findings") if isinstance(integ_env.get("findings"), dict) else integ_env
            if isinstance(integ_f, dict) and not integ_f.get("figures"):
                self._generate_integration_figures(agent_results, report_dir)
        narrative_blocks = self._collect_narrative_blocks(agent_results, exp_ctx)
        executive_summary, executive_summary_warning = (
            self._govern_executive_summary(
                executive_summary=executive_summary,
                grouped_findings=grouped_findings,
                intent=intent,
                agent_results=agent_results,
                narrative_blocks=narrative_blocks,
            )
        )
        executive_summary, executive_summary_warning, executive_summary_block = (
            self._build_executive_summary_block(
                executive_summary=executive_summary,
                executive_summary_warning=executive_summary_warning,
                grouped_findings=grouped_findings,
                intent=intent,
                exp_ctx=exp_ctx,
                agent_results=agent_results,
                narrative_blocks=narrative_blocks,
            )
        )
        narrative_blocks = [executive_summary_block, *list(narrative_blocks or [])]

        # Devil's advocate consumes deterministic evidence tiers. This is tier
        # annotation only; the single public compilation happens after its
        # caveats are attached so HTML and methodology see identical manifests.
        from aria.agents.narrative.claim_compiler import annotate_claim_tiers
        annotate_claim_tiers(narrative_blocks, exp_ctx)

        # P-DEVIL: run the deterministic devil's advocate before rendering so its
        # info caveats appear in the HTML blocks (claim tiers were annotated when
        # the blocks were collected). P-LEDGER: build the planned-vs-run manifest
        # for the provenance section.
        governance_alerts: list[str] = []
        try:
            from aria.agents.narrative.devils_advocate import build_devils_advocate
            devils_advocate = build_devils_advocate(
                narrative_blocks, agent_results, exp_ctx
            )
        except Exception as exc:  # aria: allow-broad-except surface failed devil's-advocate governance pass
            log.warning(f"Devil's-advocate pass failed: {exc}", exc_info=True)
            da_err = f"{type(exc).__name__}: {exc}"
            devils_advocate = {
                "status": "error",
                "error": da_err,
                "entries": [],
            }
            governance_alerts.append(
                f"Devil's-advocate confound review pass failed ({da_err}); "
                "associative claims carry a blocking unverified-confound caveat."
            )
            from aria.agents.narrative.types import Caveat
            for blk in narrative_blocks or []:
                tier = ((getattr(blk, "metadata", None) or {}).get("claim") or {}).get("tier")
                if getattr(blk, "status", None) == "success" and tier in {
                    "associative", "weak_mechanistic", "strong_mechanistic",
                }:
                    blk.caveats.append(Caveat(
                        text=(
                            "Devil's-advocate confound review failed "
                            f"({da_err}); technical alternatives could not be evaluated."
                        ),
                        severity="blocking",
                    ))
        try:
            from aria.agents.narrative.run_ledger import build_run_ledger
            run_ledger = build_run_ledger(exp_ctx, agent_results)
            from aria.agents.narrative.run_ledger import (
                ensure_report_ledger_nodes,
            )
            ensure_report_ledger_nodes(run_ledger, narrative_blocks)
        except Exception as exc:  # aria: allow-broad-except surface failed run-ledger governance pass
            log.warning(f"Run-ledger build failed: {exc}", exc_info=True)
            run_ledger = {
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "entries": [],
                "divergences": [],
                "n_divergences": 0,
                "modalities": [],
            }
            from aria.agents.narrative.run_ledger import (
                ensure_report_ledger_nodes,
            )
            ensure_report_ledger_nodes(run_ledger, narrative_blocks)

        # C1: the single claim-compilation boundary for every public surface.
        # Raw MessageBus findings and legacy prose are not inputs: only typed,
        # evidence-bearing NarrativeBlocks can enter. Unsupported blocks are
        # withheld without exposing their text.
        from aria.agents.narrative.claim_compiler import (
            PublicClaimCompilation,
            compile_public_claims,
        )
        try:
            public_compilation = compile_public_claims(
                narrative_blocks, exp_ctx, run_ledger=run_ledger,
                report_dir=report_dir,
            )
        except Exception as exc:  # aria: allow-broad-except surface failed claim-compiler governance pass with withheld count
            log.warning(f"Public claim compilation failed: {exc}", exc_info=True)
            cc_err = f"{type(exc).__name__}: {exc}"
            candidates = list(narrative_blocks or [])
            public_compilation = PublicClaimCompilation(
                blocks=[],
                claims=[],
                withheld=[
                    {
                        "claim_id": str(getattr(b, "id", "unknown")),
                        "status": "withheld",
                        "reason": f"claim compiler governance pass failed: {cc_err}",
                    }
                    for b in candidates
                ],
                error=cc_err,
            )
            governance_alerts.append(
                f"Claim Compiler governance pass failed ({cc_err}); "
                f"{len(candidates)} candidate claim(s) withheld."
            )
        narrative_blocks = public_compilation.blocks
        public_claims = public_compilation.claims
        findings_html = self._build_findings_section(
            findings_sections,
            agent_results,
            narrative_blocks=narrative_blocks,
            report_dir=report_dir,
        )
        # Reconcile any block withheld during HTML Findings rendering with the
        # compiled claims and methodology manifest so HTML and JSON never drift.
        render_withheld_ids = {
            _html.unescape(match)
            for match in re.findall(
                r'data-block-id="([^"]+)"\s+data-withheld="true"',
                findings_html or "",
            )
        }
        if render_withheld_ids:
            narrative_blocks = [
                b for b in narrative_blocks
                if str(getattr(b, "id", "")) not in render_withheld_ids
            ]
            public_claims = [
                c for c in public_claims
                if str(c.get("claim_id") or "") not in render_withheld_ids
            ]
            already_withheld = {
                str(w.get("claim_id") or "") for w in public_compilation.withheld
            }
            for bid in sorted(render_withheld_ids - already_withheld):
                public_compilation.withheld.append({
                    "claim_id": bid,
                    "status": "withheld",
                    "reason": "withheld during HTML narrative rendering",
                })
        # Filter devil's advocate entries to surviving compiled claims so
        # withheld claim_ids never leak into methodology.json["devils_advocate"].
        # Preserve dict error payloads when the devil's-advocate pass failed.
        surviving_claim_ids = {
            str(c.get("claim_id"))
            for c in public_claims
            if isinstance(c, dict) and c.get("claim_id")
        }
        if isinstance(devils_advocate, list):
            devils_advocate = [
                entry for entry in devils_advocate
                if isinstance(entry, dict)
                and str(entry.get("claim_id") or "") in surviving_claim_ids
            ]
        executive_summary_block = next(
            (
                block for block in narrative_blocks
                if block.id == "executive_summary"
            ),
            None,
        )
        if executive_summary_block is not None:
            executive_summary = executive_summary_block.claim
        else:
            executive_summary = (
                "The executive summary was withheld because it could not be "
                "verified against structured evidence."
            )
        self._last_governed_executive_summary = executive_summary

        try:
            from aria.agents.narrative.run_ledger import link_claims_to_ledger
            run_ledger["claim_linkage"] = link_claims_to_ledger(
                public_claims, run_ledger
            )
        except Exception as exc:  # aria: allow-broad-except surface failed claim-ledger linkage in run ledger and report
            log.warning(f"Claim-ledger linkage failed: {exc}", exc_info=True)
            err_msg = f"{type(exc).__name__}: {exc}"
            if isinstance(run_ledger, dict):
                run_ledger["claim_linkage"] = {
                    "status": "error",
                    "error": err_msg,
                    "linked": 0,
                    "unlinked": len(public_claims),
                    "n_violations": len(public_claims),
                }
            governance_alerts.append(
                f"Claim-ledger linkage pass failed ({err_msg})."
            )

        # C3: the compiler already withheld any result claim without typed
        # evidence/ledger nodes. Re-verify the accepted set as a hard invariant;
        # an unexpected mismatch must stop publication, never degrade to a loud
        # but still-public caveat.
        from aria.agents.narrative.run_ledger import verify_blocks_against_ledger
        ledger_verification = verify_blocks_against_ledger(
            narrative_blocks, run_ledger, strict=True
        )
        if isinstance(run_ledger, dict):
            run_ledger["claim_ledger_verification"] = ledger_verification

        public_claim_rows = self._build_public_claims_table(public_claims)
        withheld_entries = list(public_compilation.withheld or [])
        withheld_suffix = (
            f", {len(withheld_entries)} withheld" if withheld_entries else ""
        )
        withheld_claims_html = ""
        if withheld_entries:
            reason_counts: dict[str, int] = {}
            for w in withheld_entries:
                rsn = str(w.get("reason") or "withheld")
                reason_counts[rsn] = reason_counts.get(rsn, 0) + 1
            w_rows = "".join(
                "<tr>"
                f"<td>{cnt}</td>"
                f"<td>{_html.escape(rsn)}</td>"
                "</tr>"
                for rsn, cnt in sorted(reason_counts.items())
            )
            withheld_claims_html = (
                "<div class='warning'>"
                f"<strong>Withheld Claims ({len(withheld_entries)}):</strong> "
                f"{len(withheld_entries)} claim(s) were withheld before publication "
                "(unverified text and claim IDs omitted from HTML; see "
                "<code>methodology.json</code> &rarr; <code>claim_compilation.withheld</code>)."
                "<table><tr><th>Count</th><th>Reason</th></tr>"
                f"{w_rows}</table></div>"
            )

        # Decisions log rows
        decisions_rows = self._build_decisions_table(decisions)

        version = _html.escape(str(ARIA_VERSION))
        conflicts_html = self._plain_text_to_html(
            self._public_conflict_notice(agent_results)
        )
        methods_html = self._plain_text_to_html(methods)
        provenance = collect_provenance()
        if isinstance(exp_ctx.get("provenance"), dict):
            provenance.update(exp_ctx["provenance"])
        llm_usage_since = provenance.get("timestamp_utc")
        llm_usage = self._collect_execution_llm_usage(llm_usage_since)
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
        from aria.agents.narrative.falsification_certificates import (
            build_falsification_certificates,
        )
        from aria.agents.narrative.robustness import build_robustness_multiverse
        try:
            _prov_rm = build_robustness_multiverse(agent_results)
        except Exception as exc:  # aria: allow-broad-except surface robustness multiverse error independently
            log.warning(
                f"Robustness multiverse build failed: {exc}", exc_info=True
            )
            _prov_rm = {"status": "error", "details": f"{type(exc).__name__}: {exc}"}
        try:
            _prov_fc = build_falsification_certificates(
                public_claims,
                agent_results,
                exp_ctx,
                robustness_multiverse=_prov_rm,
            )
        except Exception as exc:  # aria: allow-broad-except surface failed falsification certificates governance pass
            log.warning(
                f"Falsification certificates build failed: {exc}",
                exc_info=True,
            )
            fc_err = f"{type(exc).__name__}: {exc}"
            _prov_fc = {
                "schema_version": "aria.falsification_certificates.v1",
                "feeds_into_parameter_selection": False,
                "status": "error",
                "details": fc_err,
                "certificates": [],
            }
            governance_alerts.append(
                f"Falsification Certificates governance pass failed ({fc_err})."
            )
        provenance_html = self._build_provenance_section(
            provenance=provenance,
            input_files=exp_ctx.get("input_files", []),
            agent_results=agent_results,
            llm_usage=llm_usage,
            run_ledger=run_ledger,
            decision_manifest=exp_ctx.get("decision_manifest"),
            fork_ledger=exp_ctx.get("analytical_fork_ledger"),
            falsification_certificates=_prov_fc,
            raw_ingestion=exp_ctx.get("raw_ingestion", []),
        )
        # ADR-057 S9: opt-in SPECULATIVE hypotheses section. Off by default, so
        # existing reports are unchanged. Only built downstream of a non-aborted
        # report (W-CLAIM/W-LEDGER passed), over audited evidence only, and fully
        # guarded so it can never break report generation.
        n_withheld = int(
            (public_compilation.summary() or {}).get("n_withheld", 0) or 0
        )
        speculative_html = self._build_speculative_section_html(
            agent_results, run_ledger, exp_ctx,
            narrative_blocks=narrative_blocks, report_dir=report_dir,
            n_withheld=n_withheld,
        )
        governance_alerts_html = ""
        if governance_alerts:
            alert_items = "".join(
                f"<li>{_html.escape(msg)}</li>" for msg in governance_alerts
            )
            governance_alerts_html = (
                "<div class='warning'><strong>Governance Pass Alert:</strong>"
                f"<ul style='margin:0.4rem 0 0 1.2rem'>{alert_items}</ul></div>"
            )
        executive_summary_warning_html = ""
        if executive_summary_warning:
            executive_summary_warning_html = (
                "<p style='color:var(--amber);font-size:0.9rem;margin-top:0.75rem'>"
                f"{_html.escape(executive_summary_warning)}"
                "</p>"
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
  <h3>Submitted Biological Question (untrusted input)</h3>
  <p><em>{_html.escape(str(question))}</em></p>
</div>

<h2>Executive Summary</h2>
{governance_alerts_html}
<div class="card">
  <p>{self._plain_text_to_html(executive_summary)}</p>
  {executive_summary_warning_html}
</div>

<h2>Provenance</h2>
{provenance_html}

<h2>Quality Control Summary</h2>
{self._build_qc_section(grouped_findings, agent_results, exp_ctx)}

<h2>Findings</h2>
{findings_html}
{speculative_html}

<h2>Governed Claim Ledger ({len(public_claims)} published{withheld_suffix})</h2>
<p style="color: var(--muted); font-size: 0.875rem;">
  Every row below was compiled from an evidence-bearing NarrativeBlock and
  verified before rendering. MessageBus events are operational records, not
  public claims.
</p>
<table>
  <tr>
    <th>Claim ID</th>
    <th>Claim</th>
    <th>Evidence scope</th>
  </tr>
  {public_claim_rows}
</table>
{withheld_claims_html}

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
        methodology_obj = self._build_methodology_json(
            provenance=provenance,
            exp_ctx=exp_ctx,
            agent_results=agent_results,
            decisions=decisions,
            llm_usage=llm_usage,
            narrative_blocks=narrative_blocks,
            run_ledger=run_ledger,
            devils_advocate=devils_advocate,
            compiled_claims=public_claims,
            claim_compilation=public_compilation.summary(),
            robustness_multiverse=_prov_rm,
            falsification_certificates=_prov_fc,
        )
        methodology_path = report_dir / "methodology.json"
        methodology_clean = _sanitize_json_value(methodology_obj)
        methodology_path.write_text(
            json.dumps(
                methodology_clean,
                indent=2,
                sort_keys=True,
                default=str,
                allow_nan=False,
            ),
            encoding="utf-8",
        )
        if isinstance(methodology_clean.get("decision_manifest"), dict):
            dm_path = report_dir / "decision_manifest.json"
            dm_path.write_text(
                json.dumps(
                    methodology_clean["decision_manifest"],
                    indent=2,
                    sort_keys=True,
                    default=str,
                    allow_nan=False,
                ),
                encoding="utf-8",
            )
        # A2: publish the minimized experiment snapshot before the RO-Crate so the
        # crate can describe it as part of this run's output graph.
        if reproducible:
            self._write_memory_snapshot(report_dir, experiment_id)
        # W-LEDGER: emit a machine-readable RO-Crate (W3C-PROV JSON-LD) next to
        # methodology.json so the run evidence graph is DOI/repository-ready. Pure
        # serialization of what was just written; never blocks report generation.
        try:
            from aria.agents.narrative.ledger_export import write_ro_crate
            write_ro_crate(report_dir)
        except Exception as exc:  # aria: allow-broad-except record non-fatal RO-Crate export failure
            log.warning(f"RO-Crate export failed: {exc}", exc_info=True)
        log.info(f"Report written to {report_path}")
        return report_path

