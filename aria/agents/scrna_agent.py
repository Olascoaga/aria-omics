"""
ARIA scRNAAgent (v4.3.1 — subprocess-only)
------------------------------------------
Single-cell RNA-seq orchestrator. The agent itself never imports scanpy.
All heavy computation runs in aria-rna-env via env_manager.run_in_stack().

Pipeline:
  1. QC              — rna_qc.py            → qc_filtered.h5ad
  2. Integration     — rna_integration.py   → integrated.h5ad  (if batch present)
  3. Resolution adv. — rna_advise_resolution.py (candidates for CP3)
  4. Clustering      — rna_clustering.py    → clustered.h5ad
  5. Annotation      — LLM proposes cell types from per-cluster top markers
  6. DE per cluster  — rna_de_per_cluster.py (respects global_padj/global_lfc)
  7. Trajectory      — rna_trajectory.py    (developmental intent only)
  8. Cell-comm       — rna_cellcomm.py      (signaling intent only)

Design principle: data flows as .h5ad paths between subprocesses. The agent
process keeps no AnnData in memory, so it cannot accidentally break the
conda-env isolation that the rest of ARIA relies on.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import uuid
from pathlib import Path

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence
from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.scrna")


SCRNA_SYSTEM = """
You are ARIA's scRNAAgent — a specialist in single-cell RNA-seq analysis.

Your expertise:
- scRNA-seq QC: doublet detection, mitochondrial filtering, MAD thresholds
- Normalization: scran, log1p; feature selection: HVGs
- Dimensionality reduction: PCA, UMAP, tSNE
- Clustering: Leiden algorithm, resolution selection
- Cell type annotation: marker-based, known cell type databases
- Differential expression: Wilcoxon per cluster, pseudo-bulk per condition
- Trajectory analysis: PAGA, DPT pseudotime, RNA velocity (scVelo)
- Cell-cell communication: LIANA, CellChat, NicheNet
- Batch correction: Harmony, scVI

Critical knowledge:
- Marker gene extraction: use names[cluster] NOT names[0]
- Single cells are not replicates: use pseudo-bulk for condition comparisons
- Mitochondrial % cutoffs must be context-aware (stressed cells have high MT%)
- Leiden resolution is data-dependent: always use ParameterAdvisor
- Harmony works on PCA embeddings: recompute neighbors after correction
- PAGA needs a root cell for meaningful pseudotime direction
""".strip()


class scRNAAgent(BaseAgent):

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
                experiment_id, current_h5ad, exp_ctx, intent, annotation
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

    _FOCUS_ALIASES = {
        "oligodendrocyte": {"Oligo"},
        "oligodendrocytes": {"Oligo"},
        "oligodendroglial": {"OPC", "Oligo"},
        "oligodendroglia": {"OPC", "Oligo"},
        "oligo": {"Oligo"},
        "oligos": {"Oligo"},
        "opc": {"OPC"},
        "opcs": {"OPC"},
        "microglia": {"Microglia"},
        "astrocyte": {"Astro"},
        "astrocytes": {"Astro"},
        "astrocito": {"Astro"},
        "astrocitos": {"Astro"},
        "microglía": {"Microglia"},
        "oligodendrocito": {"Oligo"},
        "oligodendrocitos": {"Oligo"},
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
        for token, values in cls._FOCUS_ALIASES.items():
            if re.search(rf"\b{re.escape(token)}\b", text):
                focus.update(v for v in values if v in available)
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
        Strips well-known 10x suffixes so accessions stay readable
        (e.g. GSE278576_hc11_raw_feature_bc_matrix.h5 → GSE278576_hc11).

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
        workspace = Path("~/.aria/workspace/scrna_multi").expanduser()
        workspace.mkdir(parents=True, exist_ok=True)

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
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_qc.py",
            params={
                "data_path":          path,
                "organism":           organism,
                "biological_context": intent,
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
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_integration.py",
            params={
                "data_path": input_h5ad,
                "batch_col": batch_col,
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

            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint=3,
                question=self.advisor.format_for_checkpoint(decision),
                options=[
                    f"Use recommended (resolution={decision.chosen_value})",
                    "Enter custom resolution",
                    "Skip clustering",
                ],
                context={"decision": decision.decision_id},
            )

        # Run clustering with the chosen resolution (or skip Leiden when a
        # pre-existing cluster_col is provided — rna_clustering accepts it).
        params = {
            "data_path":  input_h5ad,
            "resolution": float(decision.chosen_value) if not cluster_col else 0.5,
            "max_cells":  100_000,
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

    # Tissue keywords → CellTypist tissue_hint. Order matters: most-specific
    # first so "fetal brain" routes to fetal, not brain.
    _TISSUE_KEYWORDS = [
        ("fetal",     ["fetal", "fetus", "prenatal", "embryon"]),
        ("brain",     ["brain", "cortex", "neuron", "hippocamp", "cerebr"]),
        ("kidney",    ["kidney", "renal", "nephron"]),
        ("lung",      ["lung", "pulmonary", "alveol", "bronch"]),
        ("intestine", ["intestin", "gut", "colon", "ileum"]),
        ("skin",      ["skin", "epidermis", "dermal"]),
        ("pbmc",      ["pbmc", "peripheral blood", "blood mononuclear"]),
        ("immune",    ["immune", "t cell", "b cell", "monocyt", "lymph",
                       "macrophag", "dendritic", "nk cell"]),
    ]

    @classmethod
    def _infer_tissue_hint(cls, exp_ctx: dict, intent: dict) -> str:
        text = " ".join(filter(None, [
            intent.get("summary", ""),
            intent.get("user_question", ""),
            exp_ctx.get("user_question", ""),
            " ".join(intent.get("biological_entities", []) or []),
        ])).lower()
        for hint, keywords in cls._TISSUE_KEYWORDS:
            if any(kw in text for kw in keywords):
                return hint
        return "immune"  # CellTypist default — Immune_All_Low covers PBMC well

    def _predefined_celltype_col(self, h5ad_path: str,
                                  exp_ctx: dict) -> str | None:
        """
        Return a cell-type column name to reuse from obs, or None.

        DataAuditAgent populates exp_ctx['design']['pseudobulk']['groupby_col']
        when the input h5ad shipped with usable cell-type annotation (Seurat
        `subclass`, `cell_type_celltypist`, etc.). We trust that audit decision
        but also verify the column survives the current pipeline state (post-
        QC, post-concat) before promoting it to the canonical grouping.
        """
        design = (exp_ctx or {}).get("design", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        col = pb_cfg.get("groupby_col")
        if not col:
            inferred = (exp_ctx or {}).get("inferred_design", {}) or {}
            col = inferred.get("groupby_col")
        if not col or col == "leiden":
            return None
        # Verify the column survives in the working h5ad.
        try:
            import anndata as ad
            adata = ad.read_h5ad(h5ad_path, backed="r")
            try:
                if col not in adata.obs.columns:
                    return None
                vals = adata.obs[col].astype(str)
                levels = [v for v in vals.unique() if v and v.lower() != "nan"]
                if len(levels) < 2:
                    return None
            finally:
                backing = getattr(adata, "file", None)
                if backing is not None:
                    backing.close()
        except Exception as e:
            log.debug(f"Predefined celltype check failed for {col}: {e}")
            return None
        return col

    def _annotation_from_obs(self, experiment_id: str,
                              clustered_h5ad: str,
                              cell_type_col: str,
                              cluster_sizes: dict,
                              top_markers: dict) -> dict:
        """
        Build a findings.cell_types payload from a pre-existing obs column,
        bypassing CellTypist and the LLM reinterpretation. Each unique label
        becomes one entry whose `cluster_id` is the label itself (downstream
        groupby is the same column, so cluster IDs and labels coincide).
        """
        cell_types: dict = {}
        for label, n in cluster_sizes.items():
            cell_types[str(label)] = {
                "cell_type":         str(label),
                "celltypist_label":  None,
                "agrees_with_celltypist": None,
                "confidence":        "high",
                "rationale":         (
                    "Pre-existing annotation from input obs column "
                    f"'{cell_type_col}'. ARIA reused the user-supplied "
                    "labels and did not re-cluster or re-annotate."
                ),
                "key_markers":       top_markers.get(str(label), [])[:5],
                "n_cells":           int(n),
                "annotation_source": "input_obs",
            }

        labels_preview = [v["cell_type"] for v in list(cell_types.values())[:5]]
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Reused {len(cell_types)} cell types from input "
                f"obs['{cell_type_col}']: {labels_preview}"
             ),
             "cell_types":   cell_types,
             "celltypist":   {"ran": False, "reason": "predefined_obs"},
             "markers_used": {k: v[:5] for k, v in (top_markers or {}).items()}},
            Confidence.HIGH,
        )
        return {
            "cell_types":        cell_types,
            "markers_used":      top_markers or {},
            "n_clusters":        len(cell_types),
            "celltypist":        {"ran": False, "reason": "predefined_obs"},
            "tissue_hint":       None,
            "label_col":         cell_type_col,
            "annotated_h5ad":    clustered_h5ad,
            "annotation_source": "input_obs",
        }

    def _annotate_cell_types(self, experiment_id: str,
                              clustered_h5ad: str,
                              top_markers: dict,
                              exp_ctx: dict, intent: dict) -> dict:
        if not top_markers:
            return {"cell_types": {}, "markers_used": {}, "n_clusters": 0,
                    "celltypist": {"ran": False}}

        # ── Layer 1: CellTypist (database-backed, code-guarantee) ────────
        tissue_hint = self._infer_tissue_hint(exp_ctx, intent)
        celltypist_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_celltypist.py",
            params={
                "data_path":    clustered_h5ad,
                "organism":     exp_ctx.get("organism", "Homo sapiens"),
                "tissue_hint":  tissue_hint,
                "cluster_col":  "leiden",
                "majority_voting": True,
            },
        )

        # Trim markers for the LLM prompt (independent of CellTypist status).
        markers_for_prompt = {
            str(k): [g for g in (v or [])[:20] if g and str(g) != "nan"]
            for k, v in top_markers.items()
        }

        # ── Hard stop: when CellTypist rejects the matrix because it is
        # scaled (Seurat-style export with no recoverable log1p-CPM and no
        # pre-existing cell-type labels in obs), do NOT fall back to LLM-only
        # annotation on noisy markers. Surface the failure as an explicit
        # warning so the report does not silently report `annotation_failed`.
        unrecoverable = (
            celltypist_result.get("status") == "error"
            and celltypist_result.get("error_type") in {
                "InvalidExpressionMatrix",
                "NoUsableCountsOrLogNorm",
            }
        )
        if unrecoverable:
            reason = celltypist_result.get("details", "")[:300]
            log.warning(
                f"CellTypist rejected matrix: {celltypist_result.get('error_type')}. "
                "ARIA will not synthesise cell-type labels from noisy markers."
            )
            self.publish_finding(
                experiment_id,
                {"summary": (
                    "Cell-type annotation skipped: input expression matrix is "
                    "scaled/normalised in a way CellTypist cannot validate, "
                    "and no recoverable log1p-CPM or counts were found. "
                    "Provide raw counts or include cell-type labels in obs."
                 ),
                 "cell_types":   {},
                 "celltypist":   {
                     "ran":         False,
                     "error_type":  celltypist_result.get("error_type"),
                     "details":     reason,
                     "tissue_hint": tissue_hint,
                 },
                 "markers_used": {k: v[:5] for k, v in markers_for_prompt.items()}},
                Confidence.LOW,
            )
            return {
                "cell_types":     {},
                "markers_used":   markers_for_prompt,
                "n_clusters":     len(markers_for_prompt),
                "celltypist":     celltypist_result,
                "tissue_hint":    tissue_hint,
                "label_col":      None,
                "annotated_h5ad": None,
                "annotation_source": "unrecoverable_matrix",
                "warnings":       [
                    f"Cell-type annotation unavailable: "
                    f"{celltypist_result.get('error_type')}. {reason}"
                ],
            }

        # ── Layer 2: LLM reinterprets CellTypist results in biological
        #            context — does NOT invent labels from markers alone ──
        if celltypist_result.get("status") == "success":
            per_cluster = celltypist_result.get("per_cluster", {})
            celltypist_evidence = json.dumps(per_cluster, indent=2)
            prompt = f"""
You are reinterpreting database-backed cell type calls in their biological
context. CellTypist (model: {celltypist_result.get("model_used")}) produced
these per-cluster labels for {exp_ctx.get("organism", "?")} data.

Biological question: {intent.get("summary", exp_ctx.get("user_question", "?"))}
Tissue hint: {tissue_hint}

CellTypist per-cluster results (label = majority-voted, frequency = fraction
of cluster carrying that label, alt_labels = runner-up labels):
{celltypist_evidence}

Top marker genes per cluster (for cross-validation):
{json.dumps(markers_for_prompt, indent=2)}

Return JSON ONLY (no markdown fences):
{{
  "cluster_id": {{
    "cell_type":         "<chosen final label>",
    "celltypist_label":  "<what celltypist said>",
    "agrees_with_celltypist": true | false,
    "confidence":        "high" | "medium" | "low",
    "rationale":         "<one sentence: why this label, in this tissue>",
    "key_markers":       ["gene1", "gene2"]
  }},
  ...
}}

Rules:
- Default to the CellTypist label. Only override if the markers contradict
  it AND you can justify the override with a specific marker mismatch.
- HIGH confidence: CellTypist frequency >= 0.85 AND markers consistent.
- MEDIUM: frequency 0.5-0.85 OR markers partially support.
- LOW: frequency < 0.5 OR markers contradict — say so explicitly.
- If you override, set agrees_with_celltypist=false and explain in rationale.
- Do NOT invent labels not supported by either source.
"""
        else:
            log.warning(
                f"CellTypist failed ({celltypist_result.get('error_type', '?')}); "
                f"falling back to LLM-only annotation."
            )
            celltypist_evidence = None
            prompt = f"""
CellTypist annotation was not available for this run. Annotate clusters
from marker genes alone — be conservative about confidence.

Organism: {exp_ctx.get("organism", "unknown")}
Biological question: {intent.get("summary", "unknown")}
Cluster marker genes (top per cluster):
{json.dumps(markers_for_prompt, indent=2)}

Return JSON ONLY:
{{"cluster_id": {{"cell_type": "name", "confidence": "high|medium|low",
  "key_markers": ["gene1", "gene2"], "rationale": "<one sentence>"}}}}

Rules:
- If markers are ambiguous, say "ambiguous — possible types: X, Y"
- Without CellTypist evidence, max confidence is MEDIUM.
"""

        cell_types: dict = {}
        try:
            raw = self.llm.complete(
                prompt=prompt, system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY, max_tokens=1500,
            )
            cell_types = self._parse_annotation_json(raw)
        except Exception as e:
            log.warning(f"LLM annotation failed: {e}")

        if not cell_types:
            # Fall back to raw CellTypist labels if we have them; otherwise mark
            # clusters conservatively from canonical marker panels. The marker
            # fallback is intentionally low/medium confidence and avoids the
            # unhelpful all-"annotation_failed" state that breaks downstream
            # UMAP, trajectory, and cell-communication labels.
            if celltypist_result.get("status") == "success":
                cell_types = {
                    cl: {
                        "cell_type":               info["label"],
                        "celltypist_label":        info["label"],
                        "agrees_with_celltypist":  True,
                        "confidence":              "medium",
                        "rationale": (
                            f"LLM unavailable; using CellTypist label directly "
                            f"({info['frequency']*100:.0f}% of cluster)."
                        ),
                        "key_markers": markers_for_prompt.get(cl, [])[:5],
                    }
                    for cl, info in celltypist_result.get("per_cluster", {}).items()
                }
            else:
                cell_types = self._marker_based_annotation(markers_for_prompt)

        # Confidence summary for the bus message.
        agree_count = sum(
            1 for v in cell_types.values()
            if isinstance(v, dict) and v.get("agrees_with_celltypist") is True
        )
        labels_preview = [
            v.get("cell_type", "?") if isinstance(v, dict) else str(v)
            for v in list(cell_types.values())[:5]
        ]

        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Annotated {len(markers_for_prompt)} clusters "
                f"(CellTypist: {celltypist_result.get('model_used', 'N/A')}, "
                f"{agree_count}/{len(cell_types)} agree with LLM): "
                f"{labels_preview}"
             ),
             "cell_types":   cell_types,
             "celltypist":   {
                 "ran":         celltypist_result.get("status") == "success",
                 "model_used":  celltypist_result.get("model_used"),
                 "tissue_hint": tissue_hint,
                 "per_cluster": celltypist_result.get("per_cluster", {}),
             },
             "markers_used": {k: v[:5] for k, v in markers_for_prompt.items()}},
            Confidence.HIGH if (
                celltypist_result.get("status") == "success" and
                agree_count == len(cell_types)
            ) else Confidence.MEDIUM,
        )
        annotated_h5ad = celltypist_result.get("output_path")
        label_col = celltypist_result.get("label_col")
        if not annotated_h5ad and cell_types:
            applied = self.env.run_in_stack(
                stack="rna",
                script_path="aria/scripts/rna_apply_cluster_labels.py",
                params={
                    "data_path": clustered_h5ad,
                    "labels": cell_types,
                    "cluster_col": "leiden",
                    "label_col": "cell_type_marker",
                },
            )
            if applied.get("status") == "success":
                annotated_h5ad = applied.get("output_path")
                label_col = applied.get("label_col", "cell_type_marker")

        return {
            "cell_types":     cell_types,
            "markers_used":   markers_for_prompt,
            "n_clusters":     len(markers_for_prompt),
            "celltypist":     celltypist_result,
            "tissue_hint":    tissue_hint,
            "label_col":      label_col,
            "annotated_h5ad": annotated_h5ad,
        }

    @staticmethod
    def _marker_based_annotation(markers_by_cluster: dict) -> dict:
        panels = {
            "OPC": {"PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2", "SOX10"},
            "Oligodendrocyte": {"MBP", "PLP1", "MOG", "MOBP", "MAG", "TF"},
            "Microglia": {"P2RY12", "CX3CR1", "AIF1", "TYROBP", "C1QA",
                          "C1QB", "CSF1R", "AOAH", "HLA-A", "B2M"},
            "Astrocyte": {"AQP4", "GFAP", "ALDH1L1", "SLC1A3", "SLC1A2",
                          "CLU", "APOE"},
            "Excitatory neuron": {"SLC17A7", "SLC17A6", "CAMK2A", "SATB2",
                                  "RBFOX3", "SNAP25", "SYT1"},
            "Inhibitory neuron": {"GAD1", "GAD2", "SLC6A1", "DLX1", "DLX2",
                                  "RBFOX3", "SNAP25"},
            "Endothelial": {"CLDN5", "FLT1", "PECAM1", "VWF", "KDR"},
            "Pericyte / vascular smooth muscle": {"PDGFRB", "RGS5", "ACTA2",
                                                   "TAGLN", "MYH11"},
            "Ependymal": {"FOXJ1", "TTR", "PIFO", "DNAH5"},
        }
        out = {}
        for cluster, markers in markers_by_cluster.items():
            marker_set = {str(g).upper() for g in (markers or [])[:30]}
            hits = []
            for label, genes in panels.items():
                overlap = sorted(marker_set & genes)
                if overlap:
                    hits.append((label, overlap))
            hits.sort(key=lambda x: len(x[1]), reverse=True)
            if hits:
                label, overlap = hits[0]
                conf = "medium" if len(overlap) >= 2 else "low"
                out[str(cluster)] = {
                    "cell_type": label,
                    "confidence": conf,
                    "key_markers": overlap[:5],
                    "rationale": (
                        "Conservative marker-panel fallback; database-backed "
                        f"annotation unavailable. Matched: {', '.join(overlap[:5])}."
                    ),
                    "annotation_source": "marker_fallback",
                }
            else:
                out[str(cluster)] = {
                    "cell_type": f"Unresolved cluster {cluster}",
                    "confidence": "low",
                    "key_markers": list(markers or [])[:5],
                    "rationale": (
                        "No canonical brain/glia marker panel matched the top "
                        "cluster markers."
                    ),
                    "annotation_source": "marker_fallback",
                }
        return out

    @staticmethod
    def _parse_annotation_json(raw: str) -> dict:
        """
        Robust JSON extraction: handles ```json fences, leading prose, and
        bare JSON. Returns {} if nothing parseable is found.
        """
        if not raw:
            return {}
        stripped = raw.strip()
        # Markdown fence
        fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, re.DOTALL)
        candidate = fence.group(1) if fence else None
        if candidate is None:
            # Greedy first balanced object
            m = re.search(r"\{.*\}", stripped, re.DOTALL)
            candidate = m.group(0) if m else None
        if candidate is None:
            return {}
        try:
            data = json.loads(candidate)
            return data if isinstance(data, dict) else {}
        except (json.JSONDecodeError, TypeError):
            return {}

    # ── Differential expression per cluster ──────────────────────────────

    def _differential_expression(self, experiment_id: str,
                                   clustered_h5ad: str,
                                   intent: dict, exp_ctx: dict,
                                   groupby: str = "leiden") -> dict:
        # Respect user-confirmed thresholds from CP3 (Strict/Standard/etc.).
        padj_max = float(exp_ctx.get("global_padj", 0.05))
        lfc_min  = float(exp_ctx.get("global_lfc", 0.5))

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_de_per_cluster.py",
            params={
                "data_path": clustered_h5ad,
                "groupby":   groupby,
                "padj_max":  padj_max,
                "lfc_min":   lfc_min,
                "top_n":     20,
            },
        )

        if result.get("status") != "success":
            log.warning(f"DE per cluster failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
            return result

        n_sig = result.get("n_significant_total", 0)
        conf  = (Confidence.HIGH   if n_sig > 50 else
                 Confidence.MEDIUM if n_sig > 10 else
                 Confidence.LOW    if n_sig > 0  else
                 Confidence.INSUFFICIENT)

        # LLM interpretation grounded in actual numbers.
        de_for_llm = {
            k: v[:5] for k, v in result.get("de_genes_by_cluster", {}).items()
        }
        try:
            interpretation = self.llm.complete(
                prompt=(
                    f'Biological question: {intent.get("summary", "")}\n'
                    f'Organism: {exp_ctx.get("organism", "")}\n'
                    f'Per-cluster DE (top hits, padj<{padj_max}, '
                    f'|log2FC|>{lfc_min}):\n'
                    f'{json.dumps(de_for_llm, indent=2)}\n'
                    f'Total significant genes: {n_sig}\n\n'
                    f'Interpret in 3-4 sentences. Focus on biology, not stats.'
                ),
                system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY, max_tokens=400,
            )
        except Exception:
            interpretation = f"{n_sig} significant DE genes across clusters."

        self.publish_finding(
            experiment_id,
            {"summary":             interpretation,
             "n_significant_genes": n_sig,
             "thresholds":          {"padj_max": padj_max, "lfc_min": lfc_min}},
            conf,
        )

        return {
            "status":               result["status"],
            "padj_max":             padj_max,
            "lfc_min":              lfc_min,
            "n_significant_genes":  n_sig,
            "de_genes_by_cluster":  result.get("de_genes_by_cluster", {}),
            "n_sig_by_cluster":     result.get("n_sig_by_cluster", {}),
            "interpretation":       interpretation,
            "output_csv":           result.get("output_csv"),
        }

    # ── Pathway enrichment per cluster (ORA) ─────────────────────────────

    def _run_pathway_per_cluster(self, experiment_id: str,
                                   de_result: dict, exp_ctx: dict,
                                   data_path: str | None = None) -> dict:
        """
        Run ORA per Leiden cluster against GO_BP / KEGG / Reactome via the
        rna_pathway_per_cluster.py subprocess. Anchors cell-type biology in
        actual pathway hits instead of "what does the LLM guess about this
        marker gene".
        """
        de_by_cluster = de_result.get("de_genes_by_cluster", {})
        if not de_by_cluster:
            return {"status": "skipped", "reason": "no DE genes to enrich"}
        background_genes = self._expressed_gene_background(data_path)

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_pathway_per_cluster.py",
            params={
                "de_genes_by_cluster":   de_by_cluster,
                "organism":              exp_ctx.get("organism", "Homo sapiens"),
                "top_genes_per_cluster": 200,
                "background_genes":       background_genes,
                "padj_db_max":           0.05,
            },
            # Pathway enrichment hits Enrichr with rate limits; for a 10-cluster
            # × 3-database dataset that's 30 calls × 8s sleep = ~4 min minimum.
            # Allow up to 30 min so very dense datasets don't time out.
            timeout=1800,
        )

        if result.get("status") != "success":
            log.warning(
                f"Pathway per cluster failed: "
                f"{result.get('error_type', '?')} — "
                f"{result.get('details', '')[:200]}"
            )
            return result

        per_cluster = result.get("per_cluster", {})
        total_sig   = sum(c.get("n_significant", 0) for c in per_cluster.values())
        top_terms_preview = []
        for cl, info in list(per_cluster.items())[:4]:
            for db, terms in info.get("results", {}).items():
                if terms:
                    top_terms_preview.append(
                        f"cluster {cl}/{db}: {terms[0]['term']}"
                    )
                    break

        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Per-cluster ORA: {total_sig} significant pathway hits "
                f"across {len(per_cluster)} clusters "
                f"({', '.join(result.get('databases', {}).keys())}). "
                f"Top: {' | '.join(top_terms_preview[:3])}."
             ),
             "n_significant": total_sig,
             "databases":     result.get("databases", {}),
             "output_csv":    result.get("output_csv")},
            Confidence.HIGH if total_sig > 20 else
            Confidence.MEDIUM if total_sig > 0 else
            Confidence.INSUFFICIENT,
        )
        return result

    @staticmethod
    def _expressed_gene_background(data_path: str | None) -> list[str]:
        """Return genes detected at least once in the retained h5ad."""
        if not data_path:
            return []
        try:
            import numpy as np
            import anndata as ad
            from scipy import sparse
            adata = ad.read_h5ad(data_path)
            if adata.raw is not None:
                mat = adata.raw.X
                genes = list(adata.raw.var_names)
            else:
                mat = adata.X
                genes = list(adata.var_names)
            if sparse.issparse(mat):
                mask = np.asarray((mat > 0).sum(axis=0)).ravel() > 0
            else:
                mask = (np.asarray(mat) > 0).sum(axis=0) > 0
            return [str(g) for g, keep in zip(genes, mask) if keep]
        except Exception as exc:
            log.debug(f"Could not compute expressed-gene background: {exc}")
            return []

    # ── Pseudobulk DE between conditions ─────────────────────────────────

    PSEUDOBULK_KEYWORDS = (
        # Aging / time
        "aging", "age", "young", "old", "lifespan", "senescen",
        # Comparison / contrast verbs
        "compare", "comparison", "contrast", "between", "differential",
        "differentially expressed", "deg", "versus", " vs ", "vs.",
        # Disease vs healthy
        "disease", "disorder", "healthy", "control", "wildtype", "wild-type",
        # Treatment / perturbation
        "treatment", "treated", "untreated", "perturb", "stimulus",
        "knockout", "knock-out", "knockdown", "ko ", "wt ",
        # Common biological contrasts
        "tumor", "normal", "responder", "non-responder",
    )

    def _needs_pseudobulk(self, intent: dict, exp_ctx: dict) -> bool:
        """
        True when DesignAgent identified enough biological groups (the
        prerequisite for any between-condition DE) AND the user's intent
        carries comparison semantics. Conservative: returns False if either
        side is missing — we never run a between-condition test against a
        single group, and we don't surprise the user with extra compute
        when the question is purely descriptive (cell-type characterisation).

        v4.4 publication-readiness rule: n>=3 per group is recommended.
        n=2 remains supported with a low-power warning, but only when CP2
        selected optional supported analyses.
        """
        design = (exp_ctx or {}).get("design", {}) or {}
        groups = design.get("groups", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        has_obs_design = bool(
            pb_cfg.get("condition_col") and pb_cfg.get("replicate_col")
            and (len(groups) >= 2 or pb_cfg.get("comparisons"))
        )
        if len(groups) < 2:
            if not has_obs_design:
                return False
        if groups:
            if all(len(s or []) >= 3 for s in groups.values()):
                pass
            elif all(len(s or []) >= 2 for s in groups.values()):
                if not (exp_ctx or {}).get("run_optional_supported"):
                    return False
            else:
                return False
        text = (
            (intent.get("summary", "") or "") + " "
            + " ".join(intent.get("biological_entities", []) or [])
            + " " + (exp_ctx.get("user_question", "") or "")
        ).lower()
        return any(kw in text for kw in self.PSEUDOBULK_KEYWORDS)

    def _inject_condition_obs(self, h5ad_in: str, design: dict,
                                output_path: str) -> dict:
        """
        Write a copy of `h5ad_in` with `obs[main_factor]` populated from
        the sample_id → group mapping in `design.groups`. The cell-level
        sample column is preferred from obs in this priority:
        `sample_id` (set by rna_concat) → `batch` → `orig.ident`.

        Returns:
            {"status": "success" | "skipped" | "error",
             "output_path": str | None,
             "condition_col": str | None,
             "replicate_col": str | None,
             "matched_cells": int, "unmatched_cells": int}
        """
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_inject_condition.py",
            params={
                "data_path":    h5ad_in,
                "groups":       design.get("groups", {}),
                "factor":       design.get("main_factor", "condition"),
                "batch_col":    design.get("batch_covariate"),
                "output_path":  output_path,
            },
        )
        return result

    def _run_pseudobulk(self, experiment_id: str,
                         current_h5ad: str,
                         exp_ctx: dict, intent: dict,
                         annotation: dict) -> dict:
        """
        Pseudobulk DE between condition groups identified by DesignAgent.

        Aggregates counts per (cell_type × replicate), fits a pyDESeq2
        model with design ~ condition [+ batch], and runs ORA per
        (cell_type × comparison) on the top-N DE genes.
        """
        from pathlib import Path
        design = (exp_ctx or {}).get("design", {}) or {}
        groups = design.get("groups", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        factor = pb_cfg.get("condition_col") or design.get("main_factor", "condition")
        batch_cov = design.get("batch_covariate")

        comparisons = self._normalise_pseudobulk_comparisons(
            pb_cfg.get("comparisons")
        )
        if not comparisons:
            # Pairwise comparisons (alphabetical pairs, smaller=ref by default)
            group_names = sorted(groups.keys())
            for i, a in enumerate(group_names):
                for b in group_names[i + 1:]:
                    comparisons.append([b, a])  # test=b, ref=a
        if not comparisons:
            return {"status": "skipped", "reason": "no_comparisons"}

        # 1. Prefer h5ad-native obs design when CP1 inferred one. Otherwise
        # inject condition obs from sample → group mapping.
        workspace = Path(current_h5ad).parent / "pseudobulk"
        workspace.mkdir(parents=True, exist_ok=True)
        use_obs_design = bool(
            pb_cfg.get("from_obs") and factor and pb_cfg.get("replicate_col")
        )
        if use_obs_design:
            pb_input = current_h5ad
            replicate_col = pb_cfg.get("replicate_col")
        else:
            injected = workspace / "with_condition.h5ad"
            inj = self._inject_condition_obs(current_h5ad, design, str(injected))
            if inj.get("status") != "success":
                return {"status": "skipped",
                        "reason": f"inject_failed: {inj.get('reason', '?')}"}
            pb_input = str(injected)
            replicate_col = inj.get("replicate_col") or "sample_id"

        # 2. Choose groupby column: CellTypist labels > leiden
        cell_type_col = pb_cfg.get("groupby_col")
        if (cell_type_col == "cell_type_celltypist"
                and (annotation or {}).get("celltypist", {}).get("status") != "success"
                and (annotation or {}).get("label_col")):
            cell_type_col = (annotation or {}).get("label_col")
        if not cell_type_col:
            cell_type_col = (
                "cell_type_celltypist"
                if (annotation or {}).get("celltypist", {}).get("status") == "success"
                else (annotation or {}).get("label_col") or "leiden"
            )

        # Hard skip: when annotation failed for an unrecoverable matrix and
        # we would be left grouping by numeric Leiden IDs, pseudobulk between
        # conditions is not interpretable. Bail with a clear reason instead.
        if cell_type_col == "leiden" and (annotation or {}).get(
                "annotation_source") == "unrecoverable_matrix":
            return {"status": "skipped",
                    "reason": "no_celltype_groupby_after_annotation_failure"}
        covariates = list(pb_cfg.get("covariates") or [])
        if not covariates and batch_cov:
            covariates = [batch_cov]

        self._log_decision(
            experiment_id,
            checkpoint="scRNA",
            question="Pseudobulk design",
            decision=(
                f"groupby={cell_type_col}; condition={factor}; "
                f"replicate={replicate_col}; covariates="
                f"{', '.join(covariates) if covariates else 'none'}"
            ),
            rationale=(
                "Single-cell condition contrasts use donor/sample-level "
                "pseudobulk rather than treating cells as independent replicates."
            ),
            made_by="scrna_agent",
        )

        # 2b. (T1.1) Differential abundance BEFORE pseudobulk DE.
        # Two purposes: (a) report shifts in cell-type proportions as a
        # primary observation; (b) decide whether the pseudobulk design
        # should include a composition covariate to control for those
        # shifts when comparing within-cell-type expression.
        da_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_diff_abundance.py",
            params={
                "data_path":          str(pb_input),
                "groupby":            cell_type_col,
                "condition_col":      factor,
                "replicate_col":      replicate_col,
                "comparisons":        comparisons,
                "covariates":         covariates,
                "significance_alpha": 0.10,
                "output_dir":         str(workspace),
            },
        )
        da_significant = bool(
            da_result.get("status") == "success"
            and da_result.get("any_significant")
        )
        if da_result.get("status") == "success":
            self._log_decision(
                experiment_id,
                checkpoint="scRNA",
                question="Compositional shift before pseudobulk DE",
                decision=(
                    "composition_covariate=ON"
                    if da_significant else "composition_covariate=OFF"
                ),
                rationale=(
                    f"rna_diff_abundance flagged at least one cell type as "
                    f"significantly shifting (alpha=0.10); adding "
                    f"log(cells_in_group/total_cells) as a covariate in the "
                    f"DESeq2 design so per-cell-type DE is not confounded by "
                    f"mixture changes."
                    if da_significant else
                    "No cell type shifted significantly (alpha=0.10); the "
                    "DESeq2 design stays free of a composition covariate."
                ),
                made_by="scrna_agent",
            )

        group_sizes = [
            len(samples or []) for samples in (groups or {}).values()
        ]
        low_power_optional = (
            bool(group_sizes)
            and min(group_sizes) == 2
            and bool((exp_ctx or {}).get("run_optional_supported"))
        )
        min_replicates = 2 if low_power_optional else 3

        # 3. Run pseudobulk DE
        pb = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_pseudobulk_de.py",
            params={
                "data_path":     str(pb_input),
                "groupby":       cell_type_col,
                "condition_col": factor,
                "replicate_col": replicate_col,
                "comparisons":   comparisons,
                "covariates":    covariates,
                "composition_covariate": da_significant,
                "min_cells_per_pseudosample":   10,
                "min_replicates_per_condition": min_replicates,
                "padj_max":      0.05,
                "lfc_min":       0.5,
                "top_n":         50,
                "output_dir":    str(workspace),
            },
        )
        if pb.get("status") != "success":
            return {"status": "error",
                    "reason": f"pseudobulk_de_failed: {pb.get('error_type', '?')}",
                    "details": pb.get("details", "")[:300]}

        # 4. Pathway ORA per (group × comparison) — reuse the per-cluster
        # ORA script with synthetic cluster IDs like "Astro::old_vs_young".
        de_for_ora: dict = {}
        for group, info in (pb.get("per_group") or {}).items():
            if info.get("status") == "skipped":
                continue
            for comp_key, comp in (info.get("per_comparison") or {}).items():
                if comp.get("status") != "success" or not comp.get("all_sig"):
                    continue
                de_for_ora[f"{group}::{comp_key}"] = [
                    {"gene": r["gene"], "log2fc": r["log2fc"], "padj": r["padj"]}
                    for r in comp["all_sig"]
                ]

        pw_findings = None
        if de_for_ora:
            pw = self.env.run_in_stack(
                stack="rna",
                script_path="aria/scripts/rna_pathway_per_cluster.py",
                params={
                    "de_genes_by_cluster":   de_for_ora,
                    "organism":              exp_ctx.get("organism",
                                                          "Homo sapiens"),
                    "top_genes_per_cluster": 200,
                    "background_genes":       pb.get("background_genes", []),
                    "padj_db_max":           0.05,
                    "output_dir":            str(workspace),
                },
            )
            if pw.get("status") == "success":
                pw_findings = {
                    "organism":    pw.get("organism"),
                    "databases":   pw.get("databases", {}),
                    "background_size": pw.get("background_size"),
                    "background_source": pw.get("background_source"),
                    "per_cluster": pw.get("per_cluster", {}),
                }

        # 5. Build the finding payload + announce on the bus
        pb_payload = {
            "groupby":       pb.get("groupby"),
            "condition_col": pb.get("condition_col"),
            "replicate_col": pb.get("replicate_col"),
            "covariates":    pb.get("covariates", []),
            "from_obs":      use_obs_design,
            "thresholds":    pb.get("thresholds", {}),
            "multiple_testing": pb.get("multiple_testing", {}),
            "background_size": pb.get("background_size"),
            "background_source": pb.get("background_source"),
            "n_groups":      pb.get("n_groups"),
            "per_group":     pb.get("per_group", {}),
        }
        # Summary line
        per_group = pb.get("per_group", {}) or {}
        n_with_de = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success" and c.get("n_significant", 0) > 0
        )
        comp_str = ", ".join(f"{t}_vs_{r}" for t, r in comparisons)
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Pseudobulk DE (DESeq2) on {pb.get('n_groups', 0)} cell "
                f"types across {comp_str}: {n_with_de} (group × comparison) "
                f"blocks with significant DE."),
             "comparisons": comparisons,
             "n_groups":    pb.get("n_groups"),
             "output_csv":  pb.get("output_csv")},
            Confidence.HIGH if n_with_de > 0 else Confidence.INSUFFICIENT,
        )

        da_payload = None
        if da_result.get("status") == "success":
            da_payload = {
                "method":                 da_result.get("method"),
                "groupby":                da_result.get("groupby"),
                "condition_col":          da_result.get("condition_col"),
                "replicate_col":          da_result.get("replicate_col"),
                "covariates":             da_result.get("covariates", []),
                "significance_alpha":     da_result.get("significance_alpha"),
                "any_significant":        da_result.get("any_significant"),
                "n_replicates_per_group": da_result.get("n_replicates_per_group"),
                "per_comparison":         da_result.get("per_comparison", {}),
                "output_path":            da_result.get("output_path"),
                "warnings":               da_result.get("warnings", []),
            }
        elif da_result.get("status"):
            da_payload = {
                "status":     da_result.get("status"),
                "error_type": da_result.get("error_type"),
                "details":    (da_result.get("details") or "")[:300],
            }

        return {
            "status":                 "success",
            "differential_abundance": da_payload,
            "pseudobulk_de":          pb_payload,
            "pseudobulk_pathways":    pw_findings,
        }

    @staticmethod
    def _normalise_pseudobulk_comparisons(raw) -> list[list[str]]:
        if not raw:
            return []
        comparisons = []
        for comp in raw:
            if isinstance(comp, dict):
                test = comp.get("test") or comp.get("case") or comp.get("contrast")
                ref = comp.get("ref") or comp.get("reference") or comp.get("control")
                if test and ref:
                    comparisons.append([str(test), str(ref)])
            elif isinstance(comp, (list, tuple)) and len(comp) >= 2:
                comparisons.append([str(comp[0]), str(comp[1])])
        return comparisons

    # ── Trajectory ────────────────────────────────────────────────────────

    def _needs_trajectory(self, intent: dict) -> bool:
        keywords = ["differentiat", "develop", "pseudotime", "trajectory",
                    "progenitor", "stem", "lineage", "time course",
                    "progression", "maturation", "hematopoiesis"]
        text = (intent.get("summary", "") + " " +
                " ".join(intent.get("biological_entities", []))).lower()
        return any(kw in text for kw in keywords)

    def _run_trajectory(self, experiment_id: str,
                         clustered_h5ad: str,
                         annotation: dict, intent: dict) -> dict:
        # Prefer CellTypist labels for group naming if available — the
        # pseudotime-by-group output is far more interpretable with real
        # cell type names than with leiden numbers.
        cell_type_col = (
            "cell_type_celltypist"
            if annotation.get("celltypist", {}).get("status") == "success"
            else annotation.get("label_col") or "leiden"
        )
        self._log_decision(
            experiment_id,
            checkpoint="scRNA",
            question="Trajectory grouping",
            decision=f"PAGA/DPT grouped by {cell_type_col}",
            rationale=(
                "Trajectory analysis uses the same trusted cell grouping used "
                "for annotation and reporting."
            ),
            made_by="scrna_agent",
        )
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_trajectory.py",
            params={
                "data_path":      clustered_h5ad,
                "root_cell_type": intent.get("root_cell_type"),
                "cell_type_col":  cell_type_col,
            },
        )

        if result.get("status") == "success":
            paga       = result.get("paga", {})
            top_conn   = paga.get("top_connections", {})
            pseudotime = result.get("pseudotime", {})
            velocity   = result.get("velocity", {})

            self.publish_finding(
                experiment_id,
                {"summary": f"Trajectory: PAGA {paga.get('n_connections', 0)} transitions, "
                            f"DPT computed={pseudotime.get('computed', False)}, "
                            f"RNA velocity computed={velocity.get('computed', False)}.",
                 "paga_top":   top_conn,
                 "pseudotime": pseudotime},
                Confidence.MEDIUM,
            )
        else:
            log.warning(f"Trajectory failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
        return result

    # ── Cell-cell communication ───────────────────────────────────────────

    def _needs_cell_communication(self, intent: dict) -> bool:
        keywords = ["signal", "interact", "ligand", "receptor", "crosstalk",
                    "communication", "niche", "paracrine", "secreted",
                    "co-culture", "coculture", "microenvironment"]
        text = (intent.get("summary", "") + " " +
                " ".join(intent.get("biological_entities", []))).lower()
        return any(kw in text for kw in keywords)

    def _run_cell_communication(self, experiment_id: str,
                                 clustered_h5ad: str,
                                 exp_ctx: dict,
                                 annotation: dict | None = None) -> dict:
        # Use CellTypist labels as cell type groups when available; falls
        # back to leiden otherwise. rna_cellcomm.py already auto-falls to
        # leiden if the requested column is missing.
        cell_type_col = (
            "cell_type_celltypist"
            if (annotation or {}).get("celltypist", {}).get("status") == "success"
            else (annotation or {}).get("label_col") or "leiden"
        )
        self._log_decision(
            experiment_id,
            checkpoint="scRNA",
            question="Cell-cell communication grouping",
            decision=f"LIANA grouped by {cell_type_col}; n_perms=100",
            rationale=(
                "Ligand-receptor analysis requires annotated sender and "
                "receiver groups; ARIA reused the trusted grouping column."
            ),
            made_by="scrna_agent",
        )
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_cellcomm.py",
            params={
                "data_path":     clustered_h5ad,
                "cell_type_col": cell_type_col,
                "organism":      exp_ctx.get("organism", "Homo sapiens"),
                "n_perms":       100,
            },
        )

        if result.get("status") == "success":
            n_ia      = result.get("n_interactions", 0)
            n_types   = result.get("n_cell_types", 0)
            method    = result.get("method", "?")
            top_pairs = result.get("top_pairs", [])
            self.publish_finding(
                experiment_id,
                {"summary": f"Cell-cell communication ({method}): "
                            f"{n_ia} interactions across {n_types} clusters. "
                            f"Top pairs: {', '.join(top_pairs[:3])}."},
                Confidence.MEDIUM,
            )
        elif result.get("status") == "skipped":
            log.info(f"Cell-comm skipped: {result.get('reason', '?')}")
        else:
            log.warning(f"Cell-comm failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
        return result

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
