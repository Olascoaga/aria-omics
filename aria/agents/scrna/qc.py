"""scRNAAgent QC / integration / clustering mixin (A7 extraction from scrna_agent.py; bodies verbatim)."""
from __future__ import annotations

from aria.agents.scrna._base import *  # noqa: F401,F403


class QCClusteringMixin:
    def _run_qc(self, experiment_id: str, files: list,
                exp_ctx: dict, intent: dict) -> dict:
        """
        Single-sample → call rna_qc directly.
        Multi-sample  → QC each sample, then concatenate via rna_concat so
                        downstream Harmony has a populated obs["batch"].
        """
        if not files:
            return {"status": "error",
                    "error_type": "NoInputs",
                    "details":    "scRNA modality is empty."}

        organism = exp_ctx.get("organism", "Homo sapiens")

        # Single-sample fast-path — backwards compatible with prior runs.
        if len(files) == 1:
            return self._qc_single(experiment_id, files[0], organism, intent)

        # Multi-sample: per-sample QC followed by concat. Each rna_qc call
        # gets its own sample_id so the script writes qc_filtered_{sid}.h5ad
        # without overwriting siblings.
        workspace = self._workspace(experiment_id, "qc")

        per_sample = []
        manifest   = []
        for f in files:
            sid = self._sample_id_from_path(f)
            sample_result = self.env.run_in_stack(
                stack="rna",
                script_path="aria/scripts/rna_qc.py",
                params={
                    "data_path":          f,
                    "organism":           organism,
                    "biological_context": intent,
                    "sample_id":          sid,
                    "output_dir":         str(workspace),
                },
            )
            if sample_result.get("status") != "success":
                log.warning(f"rna_qc.py failed on sample {sid}: "
                            f"{sample_result.get('error_type', '?')}")
                return {"status":     "error",
                        "error_type": "PerSampleQCFailed",
                        "details":    (f"Sample {sid}: "
                                       f"{sample_result.get('error_type', '?')} "
                                       f"— {sample_result.get('details', '')[:200]}"),
                        "failed_sample": sid}
            per_sample.append({
                "sample_id":      sid,
                "n_cells_before": sample_result.get("n_cells_before", 0),
                "n_cells_after":  sample_result.get("n_cells_after", 0),
                "pct_removed":    sample_result.get("pct_removed", 0),
                "scrublet":       sample_result.get("scrublet", {}),
            })
            manifest.append({
                "path":      sample_result["output_path"],
                "sample_id": sid,
            })

        # Concatenate the per-sample QC outputs into one .h5ad.
        concat_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_concat.py",
            params={
                "samples":    manifest,
                "output_dir": str(workspace),
                "join":       "inner",
            },
        )
        if concat_result.get("status") != "success":
            log.warning(f"rna_concat.py failed: "
                        f"{concat_result.get('error_type', '?')}")
            return {"status":     "error",
                    "error_type": "ConcatFailed",
                    "details":    concat_result.get("details", "")[:300],
                    "per_sample": per_sample}

        n_total = concat_result.get("n_cells_total", 0)
        self.publish_finding(
            experiment_id,
            {"summary": (f"Multi-sample QC + concat: "
                         f"{len(files)} samples → {n_total} cells "
                         f"({concat_result.get('n_genes_shared', 0)} shared genes)."),
             "per_sample": per_sample},
            Confidence.HIGH,
        )

        # Aggregate per-sample counts so downstream report consumers
        # (narrative_agent QC table, _narrative_scrna.summarize_scrna_text)
        # can show "before → after" without special-casing the multi-sample
        # shape.
        n_before_total = sum(int(p.get("n_cells_before", 0)) for p in per_sample)
        n_after_total  = sum(int(p.get("n_cells_after",  0)) for p in per_sample)
        pct_removed    = (
            round(100.0 * (n_before_total - n_after_total) / n_before_total, 2)
            if n_before_total else 0.0
        )
        return {
            "status":         "success",
            "output_path":    concat_result["output_path"],
            "n_samples":      len(files),
            "n_cells_before": n_before_total,
            "n_cells_after":  n_after_total,
            "n_cells_total":  n_total,
            "pct_removed":    pct_removed,
            "n_genes_shared": concat_result.get("n_genes_shared"),
            "per_sample":     per_sample,
            "batch_col":      concat_result.get("batch_col", "batch"),
        }

    def _qc_single(self, experiment_id: str, path: str,
                   organism: str, intent: dict) -> dict:
        workspace = self._workspace(experiment_id, "qc")
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_qc.py",
            params={
                "data_path":          path,
                "organism":           organism,
                "biological_context": intent,
                "output_dir":         str(workspace),
            },
        )

        if result.get("status") != "success":
            log.warning(f"rna_qc.py failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
            return result

        n_before = result.get("n_cells_before", 0)
        n_after  = result.get("n_cells_after", n_before)
        pct_rm   = result.get("pct_removed", 0)
        mt_thr   = result.get("mt_threshold_used", "?")

        if n_after < 100:
            self.publish_finding(
                experiment_id,
                {"summary": f"Only {n_after} cells passed QC — analysis unreliable"},
                Confidence.INSUFFICIENT,
            )
            return {**result, "status": "error",
                    "error_type": "InsufficientCells"}

        self.publish_finding(
            experiment_id,
            {"summary": f"scRNA QC: {n_before} → {n_after} cells "
                        f"({pct_rm:.1f}% removed). MT ≤ {mt_thr}%",
             "warnings": result.get("warnings", [])},
            Confidence.HIGH if pct_rm < 30 else Confidence.MEDIUM,
        )
        return result

    # ── Integration ──────────────────────────────────────────────────────

    def _resolve_batch_column(self, design: dict) -> str | None:
        """
        Pick a batch column from the user-confirmed design.
        Returns None if no batch correction should run.
        """
        batch = design.get("batch_factor") or design.get("batch_covariate")
        if batch:
            return str(batch)
        return None

    def _run_integration(self, experiment_id: str,
                          input_h5ad: str, batch_col: str) -> dict:
        workspace = self._workspace(experiment_id, "integration")
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_integration.py",
            params={
                "data_path": input_h5ad,
                "batch_col": batch_col,
                "output_dir": str(workspace),
            },
        )

        if result.get("status") == "success":
            delta = result.get("batch_correction_delta", 0)
            conf  = Confidence.HIGH if delta > 0.05 else Confidence.MEDIUM
            self.publish_finding(
                experiment_id,
                {"summary": f"Harmony batch correction over '{batch_col}': "
                            f"{result.get('n_batches', '?')} batches, "
                            f"batch silhouette "
                            f"{result.get('silhouette_before', 0):.3f} → "
                            f"{result.get('silhouette_after', 0):.3f} "
                            f"(Δ={delta:+.3f}; lower = better)."},
                conf,
            )
        elif result.get("status") == "skipped":
            log.info(f"Integration skipped: {result.get('reason', '?')}")
            if result.get("n_cells"):
                self.publish_finding(
                    experiment_id,
                    {"summary": ("Harmony batch correction skipped for "
                                 f"{result.get('n_cells')} cells: "
                                 f"{result.get('reason', '')}")},
                    Confidence.MEDIUM,
                )
        else:
            log.warning(f"Integration failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
        return result

    # ── Clustering ────────────────────────────────────────────────────────

    def _run_clustering(self, experiment_id: str,
                         input_h5ad: str, intent: dict,
                         cluster_col: str | None = None
                         ) -> tuple[dict, "ParameterDecision"]:
        # ParameterAdvisor evaluates candidates in aria-rna-env via subprocess.
        # When the caller passes a pre-existing cluster_col, ARIA skips Leiden
        # entirely; we still surface a decision record so the report can show
        # provenance.
        if cluster_col:
            from aria.llm.parameter_advisor import ParameterDecision
            decision = ParameterDecision(
                decision_id=f"clustering_{experiment_id}",
                experiment_id=experiment_id,
                analysis_type="cell_type_grouping",
                parameter_name="cluster_col",
                candidates=[],
                chosen_value=cluster_col,
                chosen_by="input_obs",
                biological_context=intent or {},
                justification=(
                    f"Using pre-existing obs column '{cluster_col}' as cluster "
                    "grouping (skipping Leiden)."
                ),
                warnings=[],
            )
            self._log_decision(
                experiment_id,
                checkpoint="scRNA",
                question="Cell grouping",
                decision=f"reuse obs['{cluster_col}']; skip Leiden",
                rationale=(
                    "A trusted h5ad obs cell-type column was available, so "
                    "ARIA reused it instead of reclustering or inventing labels."
                ),
                made_by="scrna_agent",
            )
        else:
            decision = self.advisor.advise_leiden_resolution(
                data_path=input_h5ad,
                experiment_id=experiment_id,
                biological_context=intent,
            )

            _, user_decision = self.publish_blocking_escalation(
                experiment_id=experiment_id,
                checkpoint=3,
                question=self.advisor.format_for_checkpoint(decision),
                options=[
                    f"Use recommended (resolution={decision.chosen_value})",
                    "Other: Enter custom resolution",
                    "Skip clustering",
                ],
                context={
                    "decision": decision.decision_id,
                    "analysis_type": "leiden_clustering",
                    "parameter_name": "resolution",
                },
            )
            choice = self.checkpoint_choice_text(user_decision)
            if not choice:
                return {
                    "status": "error",
                    "error_type": "CheckpointUnresolved",
                    "details": "Leiden resolution checkpoint was not resolved.",
                }, decision
            if "skip" in choice.lower():
                decision.chosen_by = "user"
                decision.approved_by_user = True
                return {
                    "status": "skipped",
                    "reason": "user_skipped_clustering",
                    "details": "User skipped Leiden clustering at checkpoint 3.",
                }, decision
            override = self.numeric_checkpoint_override(user_decision, float)
            if hasattr(self.advisor, "approve_decision"):
                decision = self.advisor.approve_decision(
                    decision,
                    user_override=override,
                )
            else:
                if override is not None:
                    decision.chosen_value = override
                    decision.chosen_by = "user"
                decision.approved_by_user = True

        # Run clustering with the chosen resolution (or skip Leiden when a
        # pre-existing cluster_col is provided — rna_clustering accepts it).
        params = {
            "data_path":  input_h5ad,
            "resolution": float(decision.chosen_value) if not cluster_col else 0.5,
            "max_cells":  100_000,
            "output_dir": str(self._workspace(experiment_id, "clustering")),
        }
        if cluster_col:
            params["cluster_col"] = cluster_col
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_clustering.py",
            params=params,
        )

        if result.get("status") == "success":
            n_clusters = result.get("n_clusters", 0)
            sketch_note = ""
            if result.get("sketch_used"):
                sketch_note = (f" on a {result.get('n_cells_used')} cell "
                               f"sketch from {result.get('n_cells_total')} "
                               "QC-passed cells")
            if result.get("predef_clusters"):
                summary = (
                    f"{n_clusters} groups from obs['{result.get('groupby')}'] "
                    "(Leiden skipped; existing annotation reused)"
                )
            else:
                summary = (
                    f"{n_clusters} clusters at "
                    f"resolution={decision.chosen_value} "
                    f"(rep={result.get('rep_used', 'X_pca')})"
                    f"{sketch_note}"
                )
            self.publish_finding(
                experiment_id,
                {"summary": summary,
                 "resolution":    decision.chosen_value,
                 "cluster_sizes": result.get("cluster_sizes", {}),
                 "groupby":       result.get("groupby", "leiden")},
                Confidence.MEDIUM if result.get("sketch_used") else Confidence.HIGH,
            )
        else:
            log.warning(f"Clustering failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")

        return result, decision

    # ── Cell type annotation: CellTypist anchors, LLM reinterprets ───────

