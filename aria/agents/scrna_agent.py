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
        design = exp_ctx.get("design", {})
        batch_col = self._resolve_batch_column(design)
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
        self.publish_status(experiment_id, "Clustering...", 0.50)
        cluster_result, cluster_decision = self._run_clustering(
            experiment_id, current_h5ad, intent
        )
        findings["clustering"] = cluster_result
        findings["clustering_decision"] = {
            "resolution":    cluster_decision.chosen_value,
            "justification": cluster_decision.justification,
            "n_clusters":    cluster_result.get("n_clusters", 0),
        }
        if cluster_result.get("status") != "success":
            return {"status": "failed", "reason": "clustering_failed",
                    "findings": findings}
        current_h5ad = cluster_result["output_path"]

        # 4. Annotation (LLM proposes, code stays in charge of mapping) ──
        self.publish_status(experiment_id, "Annotating cell types...", 0.65)
        annotation = self._annotate_cell_types(
            experiment_id,
            top_markers=cluster_result.get("top_markers", {}),
            exp_ctx=exp_ctx,
            intent=intent,
        )
        findings["cell_types"] = annotation

        # 5. DE per cluster (always when ≥2 clusters) ─────────────────────
        if cluster_result.get("n_clusters", 0) >= 2:
            self.publish_status(experiment_id,
                                "Differential expression per cluster...", 0.78)
            findings["differential_expression"] = self._differential_expression(
                experiment_id, current_h5ad, intent, exp_ctx
            )

        # 6. Trajectory (developmental / time-course intent) ──────────────
        if self._needs_trajectory(intent):
            self.publish_status(experiment_id,
                                "Trajectory analysis (PAGA + DPT)...", 0.85)
            findings["trajectory"] = self._run_trajectory(
                experiment_id, current_h5ad, annotation, intent
            )

        # 7. Cell-cell communication (tissue / signaling intent) ──────────
        if self._needs_cell_communication(intent):
            self.publish_status(experiment_id,
                                "Cell-cell communication (LIANA)...", 0.92)
            findings["cell_communication"] = self._run_cell_communication(
                experiment_id, current_h5ad, exp_ctx
            )

        self.publish_status(experiment_id, "scRNAAgent complete.", 1.0)
        return {"status": "done", "findings": findings,
                "output_h5ad": current_h5ad}

    # ── QC ────────────────────────────────────────────────────────────────

    def _run_qc(self, experiment_id: str, files: list,
                exp_ctx: dict, intent: dict) -> dict:
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_qc.py",
            params={
                "data_path":          files[0] if files else "",
                "organism":           exp_ctx.get("organism", "Homo sapiens"),
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
        else:
            log.warning(f"Integration failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
        return result

    # ── Clustering ────────────────────────────────────────────────────────

    def _run_clustering(self, experiment_id: str,
                         input_h5ad: str, intent: dict) -> tuple[dict, "ParameterDecision"]:
        # ParameterAdvisor evaluates candidates in aria-rna-env via subprocess.
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

        # Run clustering with the chosen resolution.
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_clustering.py",
            params={
                "data_path":  input_h5ad,
                "resolution": float(decision.chosen_value),
            },
        )

        if result.get("status") == "success":
            n_clusters = result.get("n_clusters", 0)
            self.publish_finding(
                experiment_id,
                {"summary": f"{n_clusters} clusters at "
                            f"resolution={decision.chosen_value} "
                            f"(rep={result.get('rep_used', 'X_pca')})",
                 "resolution":    decision.chosen_value,
                 "cluster_sizes": result.get("cluster_sizes", {})},
                Confidence.HIGH,
            )
        else:
            log.warning(f"Clustering failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")

        return result, decision

    # ── Cell type annotation (LLM proposes, agent maps) ──────────────────

    def _annotate_cell_types(self, experiment_id: str,
                              top_markers: dict,
                              exp_ctx: dict, intent: dict) -> dict:
        if not top_markers:
            return {"cell_types": {}, "markers_used": {}, "n_clusters": 0}

        # Trim to top-20 per cluster for the LLM prompt.
        markers_for_prompt = {
            str(k): [g for g in (v or [])[:20] if g and str(g) != "nan"]
            for k, v in top_markers.items()
        }

        prompt = f"""
Organism: {exp_ctx.get("organism", "unknown")}
Tissue/context: {intent.get("summary", "unknown")}
Cluster marker genes (top per cluster):
{json.dumps(markers_for_prompt, indent=2)}

Annotate each cluster. Return JSON ONLY (no markdown fences, no preamble):
{{"cluster_id": {{"cell_type": "name", "confidence": "high|medium|low",
  "key_markers": ["gene1", "gene2"]}}}}

Rules:
- If markers are ambiguous, say "ambiguous — possible types: X, Y"
- Use established cell type names (T cell, not "T lymphocyte")
- HIGH confidence: 3+ unambiguous markers
- MEDIUM: 1-2 markers or tissue context required
- LOW: no clear markers
"""
        cell_types: dict = {}
        try:
            raw = self.llm.complete(
                prompt=prompt, system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY, max_tokens=1000,
            )
            cell_types = self._parse_annotation_json(raw)
        except Exception as e:
            log.warning(f"LLM annotation failed: {e}")

        if not cell_types:
            cell_types = {
                k: {"cell_type": "annotation_failed", "confidence": "low"}
                for k in markers_for_prompt
            }

        self.publish_finding(
            experiment_id,
            {"summary": f"Annotated {len(markers_for_prompt)} clusters: "
                        f"{[v.get('cell_type', '?') if isinstance(v, dict) else v for v in list(cell_types.values())[:5]]}",
             "cell_types":   cell_types,
             "markers_used": {k: v[:5] for k, v in markers_for_prompt.items()}},
            Confidence.MEDIUM,
        )
        return {
            "cell_types":   cell_types,
            "markers_used": markers_for_prompt,
            "n_clusters":   len(markers_for_prompt),
        }

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
                                   intent: dict, exp_ctx: dict) -> dict:
        # Respect user-confirmed thresholds from CP3 (Strict/Standard/etc.).
        padj_max = float(exp_ctx.get("global_padj", 0.05))
        lfc_min  = float(exp_ctx.get("global_lfc", 0.5))

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_de_per_cluster.py",
            params={
                "data_path": clustered_h5ad,
                "groupby":   "leiden",
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
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_trajectory.py",
            params={
                "data_path":      clustered_h5ad,
                "root_cell_type": intent.get("root_cell_type"),
                "cell_type_col":  "leiden",  # use cluster IDs until CellTypist lands
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
                                 exp_ctx: dict) -> dict:
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_cellcomm.py",
            params={
                "data_path":     clustered_h5ad,
                "cell_type_col": "leiden",  # until CellTypist writes cell_type
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
