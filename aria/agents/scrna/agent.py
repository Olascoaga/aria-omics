"""scRNAAgent core (A7 extraction): __init__, run loop, receive, focus routing.

Composes the concern mixins (QC/clustering, annotation, DE/pathway, pseudobulk,
trajectory/cell-communication). The class is re-exported via aria.agents.scrna
and the aria/agents/scrna_agent.py compatibility facade, so no consumer import
changed. Bodies are verbatim; behavior is pinned by
tests/test_scrna_agent_contract.py.
"""
from __future__ import annotations

from aria.agents.scrna._base import *  # noqa: F401,F403
from aria.agents.scrna.qc import QCClusteringMixin
from aria.agents.scrna.annotation import AnnotationMixin
from aria.agents.scrna.de import DEMixin
from aria.agents.scrna.pseudobulk import PseudobulkMixin
from aria.agents.scrna.advanced import AdvancedMixin


class scRNAAgent(
    QCClusteringMixin,
    AnnotationMixin,
    DEMixin,
    PseudobulkMixin,
    AdvancedMixin,
    BaseAgent,
):
    name        = "scrna_agent"
    description = ("Single-cell RNA-seq: QC, integration, clustering, "
                   "annotation, trajectory, cell-comm, DE.")

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm=llm, api_key=api_key)
        self.advisor = ParameterAdvisor(memory, llm)
        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    @staticmethod
    def _workspace(experiment_id: str, *parts: str) -> Path:
        root = Path("~/.aria/workspace").expanduser() / str(experiment_id) / "scrna"
        path = root.joinpath(*parts) if parts else root
        try:
            path.mkdir(parents=True, exist_ok=True)
            return path
        except OSError:
            fallback_root = Path("/tmp/aria_workspace") / str(experiment_id) / "scrna"
            fallback = fallback_root.joinpath(*parts) if parts else fallback_root
            fallback.mkdir(parents=True, exist_ok=True)
            return fallback

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx    = context.get("exp_context", {})
        intent     = context.get("biological_intent", {})
        modalities = exp_ctx.get("modalities", {})
        files      = modalities.get("scRNA", [])

        if not files:
            return {"status": "failed", "reason": "no_scrna_files"}

        self.publish_status(experiment_id, "scRNAAgent starting...", 0.0)
        findings: dict = {}

        focus = self._prepare_focused_h5ads(experiment_id, files, exp_ctx, intent)
        if focus.get("status") == "success":
            files = focus.get("files", files)
            findings["cell_focus"] = focus
            self._log_decision(
                experiment_id,
                checkpoint="scRNA",
                question="Focused scRNA input",
                decision=(
                    f"{focus.get('groupby')} in "
                    f"{', '.join(focus.get('values', []))}"
                ),
                rationale=(
                    f"User requested a focused cell population; subsetting "
                    f"before QC reduced input cells from "
                    f"{focus.get('n_cells_before')} to {focus.get('n_cells_after')}."
                ),
                made_by="scrna_agent",
            )
        elif focus.get("status") == "error":
            findings["cell_focus"] = focus
            return {"status": "failed", "reason": "cell_focus_failed",
                    "findings": findings}

        # 1. QC ───────────────────────────────────────────────────────────
        qc = self._run_qc(experiment_id, files, exp_ctx, intent)
        findings["qc"] = qc
        if qc.get("status") != "success":
            return {"status": "failed", "reason": "qc_failed",
                    "findings": findings}

        current_h5ad = qc.get("output_path")
        if not current_h5ad or not Path(current_h5ad).exists():
            return {"status": "failed", "reason": "qc_no_output",
                    "findings": findings}
        raw_counts_h5ad = current_h5ad

        # 2. Integration (if batch present) ───────────────────────────────
        # When QC concatenated multiple samples, it returns the batch_col
        # it populated — prefer that over the design-declared one so we
        # never silently skip Harmony on a multi-sample run.
        design = exp_ctx.get("design", {})
        batch_col = qc.get("batch_col") or self._resolve_batch_column(design)
        if batch_col:
            self.publish_status(experiment_id,
                                "Batch correction (Harmony)...", 0.30)
            integration_result = self._run_integration(
                experiment_id, current_h5ad, batch_col
            )
            findings["integration"] = integration_result
            if integration_result.get("status") == "success":
                current_h5ad = integration_result["output_path"]

        # 3. Clustering with ParameterAdvisor (CP3) ───────────────────────
        # If DataAuditAgent inferred a usable cell-type column from h5ad obs
        # (Seurat-exported h5ads with `subclass`, `cell_type`, etc.), reuse it
        # as the cluster grouping and skip Leiden entirely. This preserves the
        # user's existing annotation and avoids re-clustering work that
        # depends on a clean log-normalised X — which Seurat-scaled inputs do
        # not provide.
        predef_celltype_col = self._predefined_celltype_col(current_h5ad,
                                                            exp_ctx)
        self.publish_status(experiment_id, "Clustering...", 0.50)
        cluster_result, cluster_decision = self._run_clustering(
            experiment_id, current_h5ad, intent,
            cluster_col=predef_celltype_col,
        )
        findings["clustering"] = cluster_result
        findings["clustering_decision"] = {
            "resolution":    cluster_decision.chosen_value,
            "justification": cluster_decision.justification,
            "n_clusters":    cluster_result.get("n_clusters", 0),
            "groupby":       cluster_result.get("groupby", "leiden"),
            "predef_clusters": cluster_result.get("predef_clusters", False),
        }
        if cluster_result.get("status") != "success":
            return {"status": "failed", "reason": "clustering_failed",
                    "findings": findings}
        current_h5ad = cluster_result["output_path"]

        # X8: integration QC red-flags (overcorrection / residual batch).
        # Surfaced as a structured finding instead of a passive silhouette
        # number, only when integration actually ran.
        integ = findings.get("integration") or {}
        if integ.get("status") == "success":
            try:
                from aria.utils.integration_qc import assess_integration_quality
                cluster_sil = (
                    cluster_result.get("silhouette")
                    if cluster_result.get("silhouette") is not None
                    else getattr(cluster_decision, "metadata", {}).get("silhouette")
                    if hasattr(cluster_decision, "metadata") else None
                )
                findings["integration_qc"] = assess_integration_quality(
                    integ.get("silhouette_before"),
                    integ.get("silhouette_after"),
                    cluster_sil,
                )
                if self._integration_qc_has_blocking(findings["integration_qc"]):
                    exp_ctx = dict(exp_ctx or {})
                    exp_ctx["integration_qc"] = findings["integration_qc"]
            except Exception as exc:
                log.warning(f"Integration QC assessment failed: {exc}")

        # P1-4: hidden (unmodeled) batch red-flags. Warns when a technical/batch
        # obs column is present but was neither declared, corrected, nor modeled
        # — the case integration_qc cannot see because integration never ran on
        # it. Name + design based, no cell-level values, no correction.
        try:
            from aria.utils.batch_qc import assess_hidden_batch
            inferred = (exp_ctx or {}).get("inferred_design", {}) or {}
            obs_cols_map = inferred.get("obs_columns") or {}
            obs_columns: list[str] = []
            if isinstance(obs_cols_map, dict):
                for cols in obs_cols_map.values():
                    obs_columns.extend(str(c) for c in (cols or []))
            elif isinstance(obs_cols_map, list):
                obs_columns = [str(c) for c in obs_cols_map]
            pb_cfg = design.get("pseudobulk", {}) or {}
            findings["batch_qc"] = assess_hidden_batch(
                sorted(set(obs_columns)),
                condition_col=(pb_cfg.get("condition_col")
                               or design.get("main_factor")),
                replicate_col=pb_cfg.get("replicate_col"),
                declared_batch=(self._resolve_batch_column(design)
                                or qc.get("batch_col")),
                integration_ran=(integ.get("status") == "success"),
            )
        except Exception as exc:
            log.warning(f"Hidden-batch QC assessment failed: {exc}")

        # 4. Annotation ───────────────────────────────────────────────────
        # When clustering used a pre-existing cell-type col, the annotation
        # IS that column — we skip CellTypist and synthesise the findings.
        self.publish_status(experiment_id, "Annotating cell types...", 0.65)
        if cluster_result.get("predef_clusters"):
            annotation = self._annotation_from_obs(
                experiment_id, current_h5ad,
                cell_type_col=cluster_result.get("groupby"),
                cluster_sizes=cluster_result.get("cluster_sizes", {}),
                top_markers=cluster_result.get("top_markers", {}),
            )
        else:
            annotation = self._annotate_cell_types(
                experiment_id,
                clustered_h5ad=current_h5ad,
                top_markers=cluster_result.get("top_markers", {}),
                exp_ctx=exp_ctx,
                intent=intent,
            )
        findings["cell_types"] = annotation
        # If CellTypist annotated successfully, downstream scripts can use
        # the cell_type_celltypist column on the new annotated.h5ad.
        if annotation.get("annotated_h5ad"):
            current_h5ad = annotation["annotated_h5ad"]

        # X9: annotation-coherence check. Reused obs labels are trusted with no
        # marker verification (the fast path leaves top_markers empty), so flag
        # them as unverified; for computed clusters, flag labels lacking a
        # distinct marker signature. Data-driven, no hardcoded marker map.
        try:
            from aria.utils.annotation_qc import assess_annotation_coherence
            top_markers = cluster_result.get("top_markers", {}) or {}
            reused = bool(cluster_result.get("predef_clusters"))
            markers_present = any(
                len(m or []) > 0 for m in top_markers.values()
            )
            findings["annotation_qc"] = assess_annotation_coherence(
                top_markers,
                cluster_result.get("cluster_sizes", {}),
                reused=reused,
                markers_verified=markers_present,
            )
        except Exception as exc:
            log.warning(f"Annotation QC assessment failed: {exc}")

        # P1-4: ambient-RNA contamination red-flag. Data-driven cross-cluster
        # top-marker ubiquity (no hardcoded genes, no correction); only
        # meaningful with computed per-cluster markers.
        try:
            from aria.utils.ambient_qc import assess_ambient_contamination
            findings["ambient_qc"] = assess_ambient_contamination(
                cluster_result.get("top_markers", {}) or {},
            )
        except Exception as exc:
            log.warning(f"Ambient QC assessment failed: {exc}")

        # 5. DE per cluster (always when ≥2 clusters) ─────────────────────
        de_result = None
        if cluster_result.get("n_clusters", 0) >= 2:
            self.publish_status(experiment_id,
                                "Differential expression per cluster...", 0.75)
            de_result = self._differential_expression(
                experiment_id, current_h5ad, intent, exp_ctx,
                groupby=cluster_result.get("groupby", "leiden")
            )
            findings["differential_expression"] = de_result

            # 5b. Pathway enrichment per cluster (depends on DE) ──────────
            if de_result.get("status") == "success" and \
                    de_result.get("n_significant_genes", 0) > 0:
                self.publish_status(experiment_id,
                                    "Pathway enrichment per cluster...", 0.82)
                findings["pathways"] = self._run_pathway_per_cluster(
                    experiment_id, de_result, exp_ctx, data_path=current_h5ad
                )

        # 5c. Pseudobulk DE between conditions ────────────────────────────
        # Triggered when DesignAgent identified ≥2 biological groups AND the
        # user's question carries comparison intent (aging, treatment vs
        # control, etc.). This makes the TUI path emit the same kind of
        # per-cell-type between-condition DE that v4.3.4's harness produces.
        if self._needs_pseudobulk(intent, exp_ctx):
            self.publish_status(experiment_id,
                                "Pseudobulk DE between conditions...", 0.84)
            pb_result = self._run_pseudobulk(
                experiment_id, current_h5ad, exp_ctx, intent, annotation,
                raw_counts_h5ad=raw_counts_h5ad,
            )
            if pb_result.get("status") == "success":
                findings["pseudobulk_de"] = pb_result.get("pseudobulk_de")
                da = pb_result.get("differential_abundance")
                if da:
                    findings["differential_abundance"] = da
                pwp = pb_result.get("pseudobulk_pathways")
                if pwp:
                    findings["pseudobulk_pathways"] = pwp
            elif pb_result.get("status") == "skipped":
                log.info(
                    f"Pseudobulk skipped: {pb_result.get('reason', '?')}"
                )

        # 6. Trajectory (developmental / time-course intent) ──────────────
        if ((self._needs_trajectory(intent)
                or self._design_intelligence_optional_selected(exp_ctx, "PAGA"))
                and not self._design_intelligence_blocks(exp_ctx, "PAGA/DPT")):
            self.publish_status(experiment_id,
                                "Trajectory analysis (PAGA + DPT)...", 0.85)
            findings["trajectory"] = self._run_trajectory(
                experiment_id, current_h5ad, annotation, intent
            )

        # 7. Cell-cell communication (tissue / signaling intent) ──────────
        if ((self._needs_cell_communication(intent)
                or self._design_intelligence_optional_selected(exp_ctx, "LIANA"))
                and not self._design_intelligence_blocks(exp_ctx, "LIANA")):
            self.publish_status(experiment_id,
                                "Cell-cell communication (LIANA)...", 0.92)
            findings["cell_communication"] = self._run_cell_communication(
                experiment_id, current_h5ad, exp_ctx, annotation=annotation
            )

        self.publish_status(experiment_id, "scRNAAgent complete.", 1.0)
        return {"status": "done", "findings": findings,
                "output_h5ad": current_h5ad}

    # ── QC ────────────────────────────────────────────────────────────────

    # ARIA's own intermediate outputs that should NEVER be promoted to a
    # sample_id — passing one of these back into rna_qc would create files
    # like `qc_filtered_annotated.h5ad` or `qc_filtered_qc_filtered.h5ad`.
    _ARIA_INTERMEDIATE_STEMS = {
        "qc_filtered",
        "concatenated",
        "integrated",
        "annotated",
        "clustered",
        "clustered_sketch",
        "trajectory",
        "with_condition",
    }

    def _prepare_focused_h5ads(self, experiment_id: str, files: list,
                               exp_ctx: dict, intent: dict) -> dict:
        """
        If the user asks to focus on specific obs cell types, materialize a
        focused h5ad before QC so every downstream stage avoids unrelated cells.
        """
        groupby = self._design_groupby_col(exp_ctx)
        if not groupby:
            return {"status": "skipped", "reason": "no_groupby_column"}
        focus_values = self._infer_cell_focus_values(files, groupby, exp_ctx, intent)
        if not focus_values:
            return {"status": "skipped", "reason": "no_cell_focus_requested"}

        focused_files = []
        summaries = []
        workspace = self._scrna_focus_workspace()

        try:
            import anndata as ad
        except ImportError:
            return {
                "status": "error",
                "error_type": "MissingDependency",
                "details": "anndata is required to subset focused h5ad inputs.",
            }

        for path in files:
            if not str(path).lower().endswith(".h5ad"):
                focused_files.append(path)
                continue
            try:
                adata = ad.read_h5ad(path)
            except Exception as e:
                return {
                    "status": "error",
                    "error_type": "FocusReadFailed",
                    "details": f"{Path(path).name}: {e}",
                }
            if groupby not in adata.obs:
                focused_files.append(path)
                continue

            labels = adata.obs[groupby].astype(str)
            keep = labels.isin(focus_values).to_numpy()
            n_before = int(adata.n_obs)
            n_after = int(keep.sum())
            if n_after == 0:
                return {
                    "status": "error",
                    "error_type": "EmptyCellFocus",
                    "details": (
                        f"No cells matched {groupby} in {sorted(focus_values)} "
                        f"for {Path(path).name}."
                    ),
                }
            if n_after == n_before:
                focused_files.append(path)
            else:
                out_path = workspace / f"focused_{Path(path).stem}_{uuid.uuid4().hex[:8]}.h5ad"
                adata[keep].copy().write_h5ad(out_path)
                focused_files.append(str(out_path))
            summaries.append({
                "input_path": str(path),
                "groupby": groupby,
                "values": sorted(focus_values),
                "n_cells_before": n_before,
                "n_cells_after": n_after,
            })

        if not summaries:
            return {"status": "skipped", "reason": "no_h5ad_inputs"}

        total_before = sum(s["n_cells_before"] for s in summaries)
        total_after = sum(s["n_cells_after"] for s in summaries)
        label = ", ".join(sorted(focus_values))
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Focused scRNA input on obs['{groupby}'] in {{{label}}}: "
                f"{total_before} → {total_after} cells before QC."
            ),
             "cell_focus": summaries},
            Confidence.HIGH,
        )
        return {
            "status": "success",
            "files": focused_files,
            "groupby": groupby,
            "values": sorted(focus_values),
            "n_cells_before": total_before,
            "n_cells_after": total_after,
            "per_file": summaries,
        }

    @staticmethod
    def _scrna_focus_workspace() -> Path:
        workspace = Path("~/.aria/workspace/scrna_focus").expanduser()
        try:
            workspace.mkdir(parents=True, exist_ok=True)
            probe = workspace / ".write_test"
            probe.write_text("ok")
            probe.unlink(missing_ok=True)
            return workspace
        except OSError:
            fallback = Path("/tmp/aria_workspace/scrna_focus")
            fallback.mkdir(parents=True, exist_ok=True)
            return fallback

    @staticmethod
    def _design_groupby_col(exp_ctx: dict) -> str | None:
        design = (exp_ctx or {}).get("design", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        col = pb_cfg.get("groupby_col")
        if not col:
            inferred = (exp_ctx or {}).get("inferred_design", {}) or {}
            col = inferred.get("groupby_col")
        return col if col and col != "leiden" else None

    @classmethod
    def _infer_cell_focus_values(cls, files: list, groupby: str,
                                 exp_ctx: dict, intent: dict) -> set[str]:
        available = cls._available_groupby_values(files, groupby)
        if not available:
            return set()
        text = cls._cell_focus_text(exp_ctx, intent)
        if not text:
            return set()
        focus: set[str] = set()
        for value in available:
            if re.search(rf"\b{re.escape(value.lower())}\b", text):
                focus.add(value)
        return focus if 0 < len(focus) < len(available) else set()

    @staticmethod
    def _cell_focus_text(exp_ctx: dict, intent: dict) -> str:
        raw = str((exp_ctx or {}).get("user_question", "") or "")
        if not raw:
            raw = str((intent or {}).get("user_question", "") or "")
        if not raw:
            raw = str((intent or {}).get("summary", "") or "")
        clauses = [
            c.strip() for c in re.split(r"[\n.;]+", raw)
            if c and c.strip()
        ]
        focus_markers = (
            "focus", "focused", "focusing", "restrict", "restricted",
            "subset", "only", "exclusively", "obs[", "==",
            "solo", "sólo", "unicamente", "únicamente", "enfoc",
            "centr", "limita", "limitar",
        )
        selected = [
            c for c in clauses
            if any(marker in c.lower() for marker in focus_markers)
        ]
        return " ".join(selected).lower()

    @staticmethod
    def _available_groupby_values(files: list, groupby: str) -> set[str]:
        values: set[str] = set()
        try:
            import anndata as ad
        except ImportError:
            return values
        for path in files[:3]:
            if not str(path).lower().endswith(".h5ad"):
                continue
            try:
                adata = ad.read_h5ad(path, backed="r")
                if groupby in adata.obs:
                    values.update(
                        str(v) for v in adata.obs[groupby].dropna().unique()
                        if str(v) and str(v).lower() != "nan"
                    )
                backing_file = getattr(adata, "file", None)
                if backing_file is not None:
                    backing_file.close()
            except Exception:
                continue
        return values

    @staticmethod
    def _sample_id_from_path(path: str) -> str:
        """
        Derive a stable per-sample label from a 10x .h5 / MEX / .h5ad path.
        Strips well-known 10x suffixes so sample accessions stay readable
        (e.g. sample01_raw_feature_bc_matrix.h5 -> sample01).

        Also strips a leading `qc_filtered_` and rejects ARIA intermediate
        stems so a user who accidentally re-feeds a workspace output (e.g.
        `clustered.h5ad`) does not produce recursive `qc_filtered_clustered`
        artefacts.
        """
        stem = Path(path).stem
        for suffix in ("_raw_feature_bc_matrix",
                       "_filtered_feature_bc_matrix",
                       "_feature_bc_matrix",
                       "_matrix"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        if stem.startswith("qc_filtered_"):
            stem = stem[len("qc_filtered_"):]
        if stem in scRNAAgent._ARIA_INTERMEDIATE_STEMS:
            # Fall back to a hash of the full path so reruns over the same
            # intermediate stay stable but never collide with the canonical
            # intermediate name.
            digest = hashlib.sha1(str(path).encode()).hexdigest()[:8]
            stem = f"sample_{digest}"
        return stem or "sample"

    def receive(self, message):
        pass

    def _log_decision(self, experiment_id: str, checkpoint: str,
                      question: str, decision: str, rationale: str,
                      made_by: str = "scrna_agent") -> None:
        try:
            digest = hashlib.sha1(
                f"{checkpoint}|{question}|{decision}".encode()
            ).hexdigest()[:10]
            self.memory.store_decision(
                decision_id=f"{experiment_id}_scrna_{digest}",
                wing_id=experiment_id,
                checkpoint=checkpoint,
                question=question,
                decision=decision,
                rationale=rationale,
                made_by=made_by,
            )
        except Exception as e:
            log.warning(f"scRNA decision logging failed: {e}")

    @staticmethod
    def _design_intelligence_blocks(exp_ctx: dict, token: str) -> bool:
        di = (exp_ctx or {}).get("design_intelligence", {}) or {}
        token_l = token.lower()
        return any(token_l in str(item).lower()
                   for item in di.get("unsupported", []) or [])

    @staticmethod
    def _design_intelligence_optional_selected(exp_ctx: dict, token: str) -> bool:
        if not (exp_ctx or {}).get("run_optional_supported"):
            return False
        di = (exp_ctx or {}).get("design_intelligence", {}) or {}
        token_l = token.lower()
        return any(token_l in str(item).lower()
                   for item in di.get("optional", []) or [])
