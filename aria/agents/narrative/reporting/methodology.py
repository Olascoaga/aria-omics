"""Methodology, provenance, and scoped-memory report helpers."""

from __future__ import annotations

from aria.agents.narrative.reporting._base import *


class MethodologyMixin:
    def _write_memory_snapshot(self, report_dir: Path,
                               experiment_id: str | None = None) -> None:
        """Write a MINIMIZED, experiment-scoped memory snapshot next to the report.

        Preprint audit A2 (interim): the lab keeps one SQLite for ALL experiments
        (wings), so the previous ``shutil.copy2`` of the whole DB leaked every other
        experiment's state into the report/capsule. Export only this experiment's
        wing subtree via ``ARIAMemory.export_experiment_snapshot`` (consistent read
        under the shared lock). If scoping is impossible (no experiment_id or the
        export fails), write an EMPTY placeholder — never fall back to the full DB.
        """
        dest = report_dir / "memory_snapshot.sqlite"
        try:
            export = getattr(self.memory, "export_experiment_snapshot", None)
            db_path = getattr(self.memory, "db_path", "")
            if (experiment_id and callable(export)
                    and db_path and db_path != ":memory:"):
                result = export(experiment_id, dest)
                log.info(f"Scoped memory snapshot for {experiment_id}: "
                         f"{result.get('tables')}")
            else:
                dest.write_bytes(b"")
        except Exception as exc:
            log.warning(f"Could not write scoped memory snapshot: {exc}")
            # Fail closed: never leave (or fall back to) the full lab DB.
            try:
                dest.write_bytes(b"")
            except Exception:
                pass

    def _build_methodology_json(self, provenance: dict, exp_ctx: dict,
                                agent_results: dict, decisions: list,
                                llm_usage: dict | None = None,
                                narrative_blocks: list | None = None,
                                run_ledger: dict | None = None,
                                devils_advocate: list | None = None,
                                compiled_claims: list | None = None,
                                claim_compilation: dict | None = None) -> dict:
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
        # C3: public compilation depends on an existing typed run ledger, so the
        # standalone methodology path must build it before compiling claims.
        if run_ledger is None:
            try:
                from aria.agents.narrative.run_ledger import build_run_ledger
                run_ledger = build_run_ledger(exp_ctx, agent_results)
            except Exception as exc:
                log.warning(f"Run-ledger build failed: {exc}", exc_info=True)
                run_ledger = {"entries": [], "divergences": [],
                              "n_divergences": 0}
        from aria.agents.narrative.run_ledger import ensure_report_ledger_nodes
        ensure_report_ledger_nodes(run_ledger, list(narrative_blocks or []))
        # C1: the HTML and methodology consume the SAME pre-render compilation.
        # Standalone callers still use the public compiler, never the older
        # manifest-only path.
        if compiled_claims is None:
            try:
                from aria.agents.narrative.claim_compiler import (
                    compile_public_claims,
                )
                compilation = compile_public_claims(
                    list(narrative_blocks or []), exp_ctx,
                    run_ledger=run_ledger,
                )
                claims = compilation.claims
                claim_compilation = compilation.summary()
                narrative_blocks = compilation.blocks
            except Exception as exc:
                log.warning(
                    f"Public claim compilation failed: {exc}", exc_info=True
                )
                claims = []
                claim_compilation = {
                    "compiler": "compile_public_claims",
                    "n_published": 0,
                    "n_withheld": len(narrative_blocks or []),
                    "withheld": [],
                }
        else:
            claims = compiled_claims
        # ADR-057 rail #6: the active non-promotion wall. A machine hypothesis
        # (SPECULATIVE tier / hypothesis:// node) must never enter the audited
        # claim manifest; this raises loudly if one ever leaks here. It is a
        # no-op on clean claims (a speculation cannot become a premise).
        from aria.agents.narrative.hypothesis import assert_no_speculative_promotion
        assert_no_speculative_promotion(claims)
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
            "claim_compilation": claim_compilation or {
                "compiler": "compile_public_claims",
                "n_published": len(claims),
                "n_withheld": 0,
                "withheld": [],
            },
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
            "llm_usage": llm_usage or self._collect_execution_llm_usage(
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
            "temperature", "temperature_controlled", "seed", "seed_applied",
            "seed_deterministic", "models", "tiers",
        ):
            llm_rows.append(
                "<tr>"
                f"<td>{_html.escape(key)}</td>"
                f"<td><code>{_html.escape(str(llm_usage.get(key, 0)))}</code></td>"
                "</tr>"
            )
        backend_determinism = json.dumps(
            llm_usage.get("by_model", {}), sort_keys=True, default=str
        )
        llm_rows.append(
            "<tr>"
            "<td>backend_determinism</td>"
            f"<td><code>{_html.escape(backend_determinism)}</code></td>"
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

