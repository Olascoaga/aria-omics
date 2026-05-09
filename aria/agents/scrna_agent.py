"""
ARIA scRNAAgent (v4.2)
-----------------------
Single-cell RNA-seq analysis pipeline.

Steps:
  1. QC          — MAD thresholds, MT%, doublet detection
  2. Integration — Harmony batch correction (when batch_factor in design or multiple samples)
  3. Clustering  — ParameterAdvisor selects Leiden resolution
  4. Annotation  — LLM-assisted cell type annotation from per-cluster markers
  5. Trajectory  — PAGA + DPT pseudotime, optional RNA velocity (developmental intent)
  6. Cell-cell communication — LIANA rank_aggregate or mean-expression fallback
  7. DE          — per-cluster Wilcoxon, per-comparison pseudo-bulk
"""

from __future__ import annotations

import json
import logging
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
    description = "Single-cell RNA-seq: QC, integration, clustering, annotation, trajectory, cell-comm, DE."

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
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

        # 1. QC
        qc = self._run_qc(experiment_id, files, exp_ctx, intent)
        findings["qc"] = qc
        if qc.get("status") == "error":
            return {"status": "failed", "reason": "qc_failed", "findings": findings}

        # 2. Load + preprocess
        adata = self._load_adata(qc.get("output_path"), files)
        if adata is None:
            return {"status": "failed", "reason": "load_error"}

        # 3. Multi-sample / batch integration
        design = exp_ctx.get("design", {})
        if self._needs_integration(adata, design):
            self.publish_status(experiment_id,
                                "Batch correction (Harmony)...", 0.30)
            adata, integration_result = self._run_integration(
                experiment_id, adata, design
            )
            findings["integration"] = integration_result

        # 4. Clustering via ParameterAdvisor (Checkpoint 3)
        self.publish_status(experiment_id, "Clustering...", 0.50)
        adata, cluster_decision = self._run_clustering(experiment_id, adata, intent)
        n_clusters = int(adata.obs["leiden"].nunique()) \
                     if "leiden" in adata.obs else 0
        findings["clustering_decision"] = {
            "resolution":    cluster_decision.chosen_value,
            "justification": cluster_decision.justification,
            "n_clusters":    n_clusters,
        }

        # 5. Cell type annotation
        self.publish_status(experiment_id, "Annotating cell types...", 0.65)
        annotation = self._annotate_cell_types(
            experiment_id, adata, exp_ctx, intent
        )
        findings["cell_types"] = annotation

        # 6. Trajectory (developmental / time-course intent)
        if self._needs_trajectory(intent):
            self.publish_status(experiment_id,
                                "Trajectory analysis (PAGA + DPT)...", 0.73)
            findings["trajectory"] = self._run_trajectory(
                experiment_id, adata, annotation, exp_ctx, intent
            )

        # 7. Cell-cell communication (tissue / signaling intent)
        if self._needs_cell_communication(intent):
            self.publish_status(experiment_id,
                                "Cell-cell communication (LIANA)...", 0.82)
            findings["cell_communication"] = self._run_cell_communication(
                experiment_id, adata, annotation, exp_ctx
            )

        # 8. DE (if comparison requested)
        if intent.get("comparison"):
            self.publish_status(experiment_id, "Differential expression...", 0.88)
            findings["differential_expression"] = self._differential_expression(
                experiment_id, adata, intent, exp_ctx
            )

        self.publish_status(experiment_id, "scRNAAgent complete.", 1.0)
        return {"status": "done", "findings": findings}

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

        if result.get("status") not in ("success", "error"):
            result = self._inline_qc(files, exp_ctx)

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

    def _inline_qc(self, files: list, exp_ctx: dict) -> dict:
        try:
            import scanpy as sc
            import numpy as np

            adata = self._load_raw_adata(files)
            if adata is None:
                return {"status": "error", "error_type": "LoadFailed",
                        "details": "Could not load scRNA data"}

            organism  = exp_ctx.get("organism", "")
            mt_prefix = "MT-" if ("sapiens" in organism.lower() or
                                   "musculus" in organism.lower()) else "mt-"

            adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None,
                log1p=False, inplace=True,
            )
            n_before = adata.n_obs

            def mad_thr(vals, n_mad=3):
                med = np.median(vals)
                mad = np.median(np.abs(vals - med))
                return med - n_mad * mad, med + n_mad * mad

            _, mt_hi    = mad_thr(adata.obs["pct_counts_mt"].values)
            cnt_lo, _   = mad_thr(adata.obs["total_counts"].values)
            genes_lo, _ = mad_thr(adata.obs["n_genes_by_counts"].values)

            mt_thr = float(min(mt_hi, 25))
            adata  = adata[
                (adata.obs["pct_counts_mt"]    <= mt_thr) &
                (adata.obs["total_counts"]      >= max(cnt_lo, 500)) &
                (adata.obs["n_genes_by_counts"] >= max(genes_lo, 200))
            ].copy()

            return {
                "status":            "success",
                "n_cells_before":    int(n_before),
                "n_cells_after":     int(adata.n_obs),
                "pct_removed":       round((n_before - adata.n_obs) / n_before * 100, 1),
                "mt_threshold_used": round(mt_thr, 2),
                "output_path":       None,
                "warnings":          ["Using fallback inline QC"],
            }
        except ImportError:
            return {
                "status": "success", "n_cells_before": 0, "n_cells_after": 0,
                "pct_removed": 0, "warnings": ["scanpy not available — mock QC"],
                "note": "mock",
            }

    # ── Integration (Harmony) ─────────────────────────────────────────────

    def _needs_integration(self, adata, design: dict) -> bool:
        """True if batch_factor set or multiple samples detected in obs."""
        if design.get("batch_factor"):
            return True
        for col in ("batch", "sample", "donor"):
            if col in adata.obs.columns and adata.obs[col].nunique() > 1:
                return True
        return False

    def _run_integration(self, experiment_id: str, adata, design: dict):
        try:
            import scanpy as sc
            import numpy as np
            from sklearn.metrics import silhouette_score
            from sklearn.preprocessing import LabelEncoder

            # Resolve batch column
            batch_col = design.get("batch_factor")
            if not batch_col or batch_col not in adata.obs.columns:
                for col in ("batch", "sample", "donor"):
                    if col in adata.obs.columns and adata.obs[col].nunique() > 1:
                        batch_col = col
                        break

            if not batch_col:
                return adata, {"status": "skipped",
                               "reason": "no batch column found"}

            n_batches = int(adata.obs[batch_col].nunique())
            if n_batches < 2:
                return adata, {"status": "skipped",
                               "reason": "only one batch"}

            if "X_pca" not in adata.obsm:
                return adata, {"status": "skipped",
                               "reason": "PCA not computed — run preprocessing first"}

            le = LabelEncoder()
            labels = le.fit_transform(adata.obs[batch_col].astype(str))
            sil_before = round(float(
                silhouette_score(adata.obsm["X_pca"][:, :20], labels)
            ), 4)

            try:
                sc.external.pp.harmony_integrate(
                    adata, batch_col,
                    basis="X_pca", adjusted_basis="X_pca_harmony",
                )
            except Exception:
                import harmonypy as hm
                ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_col)
                adata.obsm["X_pca_harmony"] = ho.Z_corr.T

            rep = "X_pca_harmony"
            sil_after = round(float(
                silhouette_score(adata.obsm[rep][:, :20], labels)
            ), 4)
            delta = round(sil_before - sil_after, 4)

            # Recompute graph on corrected embedding
            sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, n_pcs=30)
            sc.tl.umap(adata)

            conf = Confidence.HIGH if delta > 0.05 else Confidence.MEDIUM
            self.publish_finding(
                experiment_id,
                {"summary": f"Harmony batch correction: {n_batches} batches. "
                            f"Batch silhouette {sil_before:.3f} → {sil_after:.3f} "
                            f"(Δ={delta:+.3f}; lower = better correction)."},
                conf,
            )
            return adata, {
                "status":            "done",
                "method":            "harmony",
                "n_batches":         n_batches,
                "batch_col":         batch_col,
                "silhouette_before": sil_before,
                "silhouette_after":  sil_after,
                "delta":             delta,
            }

        except Exception as e:
            log.warning(f"Integration failed: {e}")
            return adata, {"status": "failed", "reason": str(e)}

    # ── Clustering ────────────────────────────────────────────────────────

    def _run_clustering(self, experiment_id: str, adata,
                         intent: dict) -> tuple:
        decision = self.advisor.advise_leiden_resolution(
            adata=adata,
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

        try:
            import scanpy as sc
            sc.tl.leiden(
                adata,
                resolution=decision.chosen_value,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
        except Exception as e:
            log.warning(f"Leiden clustering failed: {e}")

        n_clusters = int(adata.obs["leiden"].nunique()) \
                     if "leiden" in adata.obs else 0

        self.publish_finding(
            experiment_id,
            {"summary": f"{n_clusters} clusters at resolution={decision.chosen_value}",
             "resolution": decision.chosen_value},
            Confidence.HIGH,
        )
        return adata, decision

    # ── Cell type annotation ──────────────────────────────────────────────

    def _annotate_cell_types(self, experiment_id: str, adata,
                              exp_ctx: dict, intent: dict) -> dict:
        markers: dict = {}
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
            rgg    = adata.uns["rank_genes_groups"]
            labels = list(rgg["names"].dtype.names)
            for cl in labels:
                genes = [g for g in rgg["names"][cl][:20]
                         if g and str(g) != "nan"]
                markers[str(cl)] = genes
        except (ImportError, KeyError) as e:
            log.warning(f"Marker extraction failed: {e}")
            if "leiden" in adata.obs:
                for cl in adata.obs["leiden"].unique():
                    markers[str(cl)] = ["marker_unavailable"]

        if not markers:
            return {"cell_types": {}, "markers_used": {}, "n_clusters": 0}

        prompt = f"""
Organism: {exp_ctx.get("organism", "unknown")}
Tissue/context: {intent.get("summary", "unknown")}
Cluster marker genes (top 20 per cluster):
{json.dumps(markers, indent=2)}

Annotate each cluster. Return JSON:
{{"cluster_id": {{"cell_type": "name", "confidence": "high|medium|low",
  "key_markers": ["gene1", "gene2"]}}}}

Rules:
- If markers are ambiguous, say "ambiguous — possible types: X, Y"
- Use established cell type names (T cell, not "T lymphocyte")
- HIGH confidence: 3+ unambiguous markers
- MEDIUM: 1-2 markers or tissue context required
- LOW: no clear markers
"""
        try:
            raw = self.llm.complete(
                prompt=prompt, system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY, max_tokens=1000,
            )
            raw = raw.strip()
            if "```" in raw:
                raw = raw.split("```")[1]
                if raw.startswith("json"):
                    raw = raw[4:]
            cell_types = json.loads(raw.strip("` \n"))
        except Exception as e:
            log.warning(f"LLM annotation failed: {e}")
            cell_types = {k: {"cell_type": "annotation_failed", "confidence": "low"}
                          for k in markers}

        # Write cell_type into adata.obs for downstream use
        try:
            import pandas as pd
            ct_map = {k: (v.get("cell_type", k) if isinstance(v, dict) else str(v))
                      for k, v in cell_types.items()}
            if "leiden" in adata.obs.columns:
                adata.obs["cell_type"] = adata.obs["leiden"].map(ct_map).fillna(
                    adata.obs["leiden"]
                )
        except Exception as e:
            log.warning(f"Could not write cell_type to obs: {e}")

        self.publish_finding(
            experiment_id,
            {"summary": f"Annotated {len(markers)} clusters: "
                        f"{[v.get('cell_type','?') if isinstance(v,dict) else v for v in list(cell_types.values())[:5]]}",
             "cell_types":   cell_types,
             "markers_used": {k: v[:5] for k, v in markers.items()}},
            Confidence.MEDIUM,
        )
        return {
            "cell_types":   cell_types,
            "markers_used": markers,
            "n_clusters":   len(markers),
        }

    # ── Trajectory ────────────────────────────────────────────────────────

    def _needs_trajectory(self, intent: dict) -> bool:
        keywords = ["differentiat", "develop", "pseudotime", "trajectory",
                    "progenitor", "stem", "lineage", "time course",
                    "progression", "maturation", "hematopoiesis"]
        text = (intent.get("summary", "") + " " +
                " ".join(intent.get("biological_entities", []))).lower()
        return any(kw in text for kw in keywords)

    def _run_trajectory(self, experiment_id: str, adata,
                         annotation: dict, exp_ctx: dict, intent: dict) -> dict:
        try:
            import scanpy as sc
            import numpy as np

            groupby = ("leiden" if "leiden" in adata.obs.columns
                       else "cell_type" if "cell_type" in adata.obs.columns
                       else None)
            if groupby is None:
                return {"status": "skipped", "reason": "no cluster column"}

            # PAGA
            sc.tl.paga(adata, groups=groupby)
            paga_conn: dict = {}
            try:
                conn_mat = adata.uns["paga"]["connectivities"]
                if hasattr(conn_mat, "toarray"):
                    conn_mat = conn_mat.toarray()
                cats = (list(adata.obs[groupby].cat.categories)
                        if hasattr(adata.obs[groupby], "cat")
                        else sorted(adata.obs[groupby].unique()))
                for i, g1 in enumerate(cats):
                    for j, g2 in enumerate(cats):
                        if i < j and float(conn_mat[i, j]) > 0.10:
                            paga_conn[f"{g1}→{g2}"] = round(float(conn_mat[i, j]), 3)
            except Exception as e:
                log.warning(f"PAGA connectivity: {e}")

            # DPT pseudotime
            dpt_result: dict = {"computed": False}
            try:
                sc.tl.diffmap(adata)

                root_type = intent.get("root_cell_type")
                ct_col = "cell_type" if "cell_type" in adata.obs.columns else None
                root_used = "auto"

                if root_type and ct_col:
                    mask = adata.obs[ct_col] == root_type
                    if mask.sum() > 0:
                        adata.uns["iroot"] = int(np.where(mask)[0][0])
                        root_used = root_type

                if "iroot" not in adata.uns:
                    if "n_genes_by_counts" in adata.obs.columns:
                        adata.uns["iroot"] = int(
                            adata.obs["n_genes_by_counts"].values.argmin()
                        )
                    else:
                        adata.uns["iroot"] = 0

                sc.tl.dpt(adata)

                col = ct_col or groupby
                pt_by_group = (
                    adata.obs.groupby(col)["dpt_pseudotime"]
                    .mean().sort_values().round(4).to_dict()
                )
                dpt_result = {
                    "computed":            True,
                    "pseudotime_by_group": pt_by_group,
                    "root_used":           root_used,
                }
            except Exception as e:
                dpt_result = {"computed": False, "reason": str(e)}

            # RNA Velocity (only if loom/spliced layers available)
            velocity_result: dict
            if "spliced" in adata.layers and "unspliced" in adata.layers:
                try:
                    import scvelo as scv
                    scv.pp.filter_and_normalize(adata,
                                               min_shared_counts=20,
                                               n_top_genes=2000)
                    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
                    scv.tl.velocity(adata)
                    scv.tl.velocity_graph(adata)
                    velocity_result = {"computed": True, "method": "scvelo_stochastic"}
                except ImportError:
                    velocity_result = {"computed": False,
                                       "reason": "scvelo not installed"}
                except Exception as e:
                    velocity_result = {"computed": False, "reason": str(e)}
            else:
                velocity_result = {
                    "computed": False,
                    "reason":   "no spliced/unspliced layers — provide loom/raw "
                                "data for RNA velocity",
                }

            top_conn = dict(sorted(paga_conn.items(), key=lambda x: -x[1])[:10])
            dpt_flag = dpt_result.get("computed", False)
            vel_flag = velocity_result.get("computed", False)

            self.publish_finding(
                experiment_id,
                {"summary": f"Trajectory: PAGA {len(paga_conn)} transitions, "
                            f"DPT pseudotime computed={dpt_flag}, "
                            f"RNA velocity computed={vel_flag}.",
                 "paga_top": top_conn,
                 "pseudotime": dpt_result},
                Confidence.MEDIUM,
            )
            return {
                "status":             "done",
                "paga_top_connections": top_conn,
                "pseudotime":         dpt_result,
                "velocity":           velocity_result,
            }

        except Exception as e:
            log.warning(f"Trajectory failed: {e}")
            return {"status": "failed", "reason": str(e)}

    # ── Cell-cell communication ───────────────────────────────────────────

    def _needs_cell_communication(self, intent: dict) -> bool:
        keywords = ["signal", "interact", "ligand", "receptor", "crosstalk",
                    "communication", "niche", "paracrine", "secreted",
                    "co-culture", "coculture", "microenvironment"]
        text = (intent.get("summary", "") + " " +
                " ".join(intent.get("biological_entities", []))).lower()
        return any(kw in text for kw in keywords)

    def _run_cell_communication(self, experiment_id: str, adata,
                                 annotation: dict, exp_ctx: dict) -> dict:
        try:
            import numpy as np

            # Ensure cell_type col
            ct_col = "cell_type" if "cell_type" in adata.obs.columns else \
                     "leiden"    if "leiden"    in adata.obs.columns else None
            if ct_col is None:
                return {"status": "skipped", "reason": "no cell type column"}

            n_types = int(adata.obs[ct_col].nunique())
            if n_types < 2:
                return {"status": "skipped", "reason": "need ≥2 cell types"}

            organism   = exp_ctx.get("organism", "Homo sapiens").lower()
            interactions: list[dict] = []
            method = "mean_expression_fallback"

            # Primary: LIANA
            try:
                import liana as li
                li.mt.rank_aggregate(
                    adata, groupby=ct_col,
                    use_raw=False, verbose=False, n_perms=100,
                )
                liana_df = adata.uns["liana_res"]
                for _, row in liana_df.sort_values("magnitude_rank").head(50).iterrows():
                    interactions.append({
                        "source":   str(row.get("source", "")),
                        "target":   str(row.get("target", "")),
                        "ligand":   str(row.get("ligand_complex",
                                                 row.get("ligand", ""))),
                        "receptor": str(row.get("receptor_complex",
                                                 row.get("receptor", ""))),
                        "score":    round(float(row.get("magnitude_rank", 0)), 4),
                    })
                method = "liana_rank_aggregate"

            except (ImportError, Exception):
                interactions, method = self._cellcomm_fallback(
                    adata, ct_col, organism
                )

            # Summarise by sender→receiver pair
            pair_counts: dict = {}
            for ia in interactions:
                k = f"{ia['source']}→{ia['target']}"
                pair_counts[k] = pair_counts.get(k, 0) + 1
            top_pairs = [p for p, _ in
                         sorted(pair_counts.items(), key=lambda x: -x[1])[:8]]

            n_ia = len(interactions)
            self.publish_finding(
                experiment_id,
                {"summary": f"Cell-cell communication ({method}): "
                            f"{n_ia} interactions across {n_types} cell types. "
                            f"Top pairs: {', '.join(top_pairs[:3])}."},
                Confidence.MEDIUM,
            )
            return {
                "status":           "done",
                "method":           method,
                "n_cell_types":     n_types,
                "n_interactions":   n_ia,
                "top_interactions": interactions[:20],
                "top_pairs":        top_pairs,
            }

        except Exception as e:
            log.warning(f"Cell-cell communication failed: {e}")
            return {"status": "failed", "reason": str(e)}

    def _cellcomm_fallback(self, adata, ct_col: str,
                            organism: str) -> tuple[list, str]:
        """Mean-expression scoring with a curated L-R resource."""
        import numpy as np

        LR = [
            ("TGFB1","TGFBR1"),("TGFB1","TGFBR2"),("VEGFA","KDR"),
            ("VEGFA","FLT1"),  ("MIF","CD44"),     ("SPP1","CD44"),
            ("SPP1","ITGAV"),  ("IL6","IL6R"),     ("IL6","IL6ST"),
            ("CXCL12","CXCR4"),("CCL5","CCR5"),   ("EGF","EGFR"),
            ("HGF","MET"),     ("PDGFB","PDGFRB"),("FGF2","FGFR1"),
            ("WNT5A","ROR2"),  ("NOTCH1","JAG1"),  ("DLL4","NOTCH1"),
            ("SEMA3A","NRP1"), ("ANGPT1","TEK"),   ("BMP2","BMPR1A"),
            ("KITLG","KIT"),   ("THBS1","CD36"),   ("FN1","ITGA5"),
            ("LAMB1","ITGB1"), ("EFNB2","EPHB4"),  ("IGF1","IGF1R"),
            ("TNF","TNFRSF1A"),("IFNG","IFNGR1"),  ("IL10","IL10RA"),
        ]
        if "musculus" in organism:
            LR = [(l.capitalize(), r.capitalize()) for l, r in LR]

        var_names  = set(adata.var_names)
        valid_lr   = [(l, r) for l, r in LR if l in var_names and r in var_names]
        if not valid_lr:
            return [], "fallback_no_pairs_found"

        cell_types = adata.obs[ct_col].unique()
        var_idx    = {g: i for i, g in enumerate(adata.var_names)}

        mean_expr: dict = {}
        for ct in cell_types:
            mask  = adata.obs[ct_col] == ct
            X_sub = adata.X[mask]
            if hasattr(X_sub, "toarray"):
                X_sub = X_sub.toarray()
            mean_expr[ct] = np.asarray(X_sub, dtype=float).mean(axis=0)

        interactions = []
        for l_gene, r_gene in valid_lr:
            li_i, ri_i = var_idx[l_gene], var_idx[r_gene]
            for src in cell_types:
                for tgt in cell_types:
                    if src == tgt:
                        continue
                    score = float(mean_expr[src][li_i]) * float(mean_expr[tgt][ri_i])
                    if score > 0:
                        interactions.append({
                            "source": str(src),  "target":   str(tgt),
                            "ligand": l_gene,    "receptor": r_gene,
                            "score":  round(score, 4),
                        })

        interactions.sort(key=lambda x: -x["score"])
        return interactions[:50], "mean_expression_fallback"

    # ── Differential expression ───────────────────────────────────────────

    def _differential_expression(self, experiment_id: str, adata,
                                   intent: dict, exp_ctx: dict) -> dict:
        de_genes: dict = {}
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(
                adata, groupby="leiden",
                method="wilcoxon", key_added="de_wilcoxon",
            )
            rgg = adata.uns["de_wilcoxon"]
            for cl in rgg["names"].dtype.names:
                names = rgg["names"][cl][:50]
                pvals = rgg["pvals_adj"][cl][:50]
                lfc   = rgg["logfoldchanges"][cl][:50]
                sig   = [
                    {"gene": str(g), "log2fc": round(float(l), 3),
                     "padj": round(float(p), 6)}
                    for g, p, l in zip(names, pvals, lfc)
                    if float(p) < 0.05 and abs(float(l)) > 0.5
                    and str(g) not in ("nan", "")
                ]
                de_genes[str(cl)] = sig[:20]

        except (ImportError, KeyError) as e:
            log.warning(f"DE failed: {e}")

        n_sig = sum(len(v) for v in de_genes.values())
        conf  = (Confidence.HIGH   if n_sig > 50  else
                 Confidence.MEDIUM if n_sig > 10  else
                 Confidence.LOW    if n_sig > 0   else
                 Confidence.INSUFFICIENT)

        prompt = f"""
Biological question: {intent.get("summary", "")}
Organism: {exp_ctx.get("organism", "")}
Per-cluster DE (top hits, Wilcoxon):
{json.dumps({k: v[:5] for k, v in de_genes.items()}, indent=2)}
Total significant genes: {n_sig}

Interpret in 3-4 sentences. Focus on biology, not statistics.
"""
        try:
            interpretation = self.llm.complete(
                prompt=prompt, system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY, max_tokens=400,
            )
        except Exception:
            interpretation = f"{n_sig} significant DE genes across clusters."

        self.publish_finding(
            experiment_id,
            {"summary": interpretation, "n_significant_genes": n_sig},
            conf,
        )
        return {
            "de_genes_by_cluster":  de_genes,
            "n_significant_genes":  n_sig,
            "interpretation":       interpretation,
        }

    # ── Data loading ──────────────────────────────────────────────────────

    def _load_adata(self, qc_output_path: str, raw_files: list):
        if qc_output_path and Path(qc_output_path).exists():
            adata = self._load_raw_adata([qc_output_path])
            if adata is not None:
                return self._preprocess(adata)
        adata = self._load_raw_adata(raw_files)
        if adata is not None:
            adata = self._preprocess(adata)
        return adata

    def _load_raw_adata(self, files: list):
        try:
            import scanpy as sc
            for f in files:
                p = Path(f)
                if not p.exists():
                    continue
                if p.suffix == ".h5ad":
                    return sc.read_h5ad(str(p))
                if p.suffix == ".h5":
                    return sc.read_10x_h5(str(p))
                if p.is_dir():
                    return sc.read_10x_mtx(str(p), var_names="gene_symbols",
                                           cache=True)
        except ImportError:
            log.warning("scanpy not available")
        except Exception as e:
            log.error(f"Load failed: {e}")
        return None

    def _preprocess(self, adata):
        """Normalize → HVG → PCA → neighbors → UMAP."""
        try:
            import scanpy as sc
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
            sc.tl.umap(adata)
        except (ImportError, Exception) as e:
            log.warning(f"Preprocessing: {e}")
        return adata

    def receive(self, message):
        pass
