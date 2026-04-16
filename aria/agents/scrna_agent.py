"""
ARIA scRNAAgent
---------------
Single-cell RNA-seq analysis only.

Pipeline:
  1. QC          — MAD thresholds, MT%, doublet detection (scrublet)
                   Delegated to aria/scripts/rna_qc.py via EnvironmentManager
  2. Clustering  — ParameterAdvisor selects Leiden resolution
                   Delegated to aria/scripts/rna_clustering.py
  3. Annotation  — LLM-assisted cell type annotation from marker genes
                   Fix: each cluster gets ITS OWN markers (names[cluster])
  4. DE          — per-cluster Wilcoxon, per-comparison pseudo-bulk
                   Fix: per-cluster extraction (not names[0] for all)

Does NOT handle: bulk RNA-seq, pathway enrichment (→ BulkRNAAgent),
                 spatial transcriptomics (→ future SpatialAgent)
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
- Trajectory analysis: RNA velocity, diffusion pseudotime

Critical knowledge:
- Marker gene extraction: use names[cluster] NOT names[0]
  (names[0] gives the same genes for every cluster — a common bug)
- Single cells are not replicates for statistical testing
  Use pseudo-bulk for condition comparisons, not per-cell DE
- Mitochondrial % cutoffs must be context-aware:
  Stressed/activated cells have legitimately high MT%
- Leiden resolution is data-dependent: always use ParameterAdvisor
  Do not hardcode resolution values

Always distinguish biological signal from technical variation.
""".strip()


class scRNAAgent(BaseAgent):

    name        = "scrna_agent"
    description = "Single-cell RNA-seq: QC, clustering, annotation, DE."

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

        findings = {}

        # 1. QC via isolated environment
        qc = self._run_qc(experiment_id, files, exp_ctx, intent)
        findings["qc"] = qc
        if qc.get("status") == "error":
            return {"status": "failed", "reason": "qc_failed",
                    "findings": findings}

        # Load filtered AnnData for ParameterAdvisor
        adata = self._load_adata(qc.get("output_path"), files)
        if adata is None:
            return {"status": "failed", "reason": "load_error"}

        # 2. Clustering via ParameterAdvisor (Checkpoint 3)
        adata, cluster_decision = self._run_clustering(
            experiment_id, adata, intent
        )
        findings["clustering"] = {
            "resolution":    cluster_decision.chosen_value,
            "justification": cluster_decision.justification,
        }

        # 3. Cell type annotation
        self.publish_status(experiment_id, "Annotating cell types...", 0.7)
        annotation = self._annotate_cell_types(
            experiment_id, adata, exp_ctx, intent
        )
        findings["cell_types"] = annotation

        # 4. Differential expression (if comparison requested)
        if intent.get("comparison"):
            self.publish_status(experiment_id, "Running DE analysis...", 0.85)
            de = self._differential_expression(
                experiment_id, adata, intent, exp_ctx
            )
            findings["differential_expression"] = de

        self.publish_status(experiment_id, "scRNAAgent complete.", 1.0)
        return {"status": "done", "findings": findings}

    # ── QC ────────────────────────────────────────────────────────────────

    def _run_qc(self, experiment_id: str, files: list,
                exp_ctx: dict, intent: dict) -> dict:
        """QC via rna_qc.py in aria-rna-env."""
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
            # Fallback to inline QC
            result = self._inline_qc(files, exp_ctx)

        n_before = result.get("n_cells_before", 0)
        n_after  = result.get("n_cells_after",  result.get("n_cells_before", 0))
        pct_rm   = result.get("pct_removed",    0)
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
        """Fallback inline QC when aria-rna-env not available."""
        try:
            import scanpy as sc
            import numpy as np

            adata = self._load_raw_adata(files)
            if adata is None:
                return {"status": "error", "error_type": "LoadFailed",
                        "details": "Could not load scRNA data"}

            organism  = exp_ctx.get("organism", "")
            mt_prefix = "MT-" if "sapiens" in organism.lower() or \
                                  "musculus" in organism.lower() else "mt-"

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

            _, mt_hi   = mad_thr(adata.obs["pct_counts_mt"].values)
            cnt_lo, _  = mad_thr(adata.obs["total_counts"].values)
            genes_lo, _ = mad_thr(adata.obs["n_genes_by_counts"].values)

            mt_thr = float(min(mt_hi, 25))
            adata  = adata[
                (adata.obs["pct_counts_mt"]      <= mt_thr) &
                (adata.obs["total_counts"]        >= max(cnt_lo, 500)) &
                (adata.obs["n_genes_by_counts"]   >= max(genes_lo, 200))
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
                "status":         "success",
                "n_cells_before": 0,
                "n_cells_after":  0,
                "pct_removed":    0,
                "warnings":       ["scanpy not available — mock QC"],
                "note":           "mock",
            }

    # ── Clustering ────────────────────────────────────────────────────────

    def _run_clustering(self, experiment_id: str, adata,
                         intent: dict) -> tuple:
        """
        Cluster via ParameterAdvisor → Checkpoint 3 → Leiden.
        """
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

        # Apply clustering
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
            {"summary": f"{n_clusters} clusters found at "
                        f"resolution={decision.chosen_value}",
             "resolution": decision.chosen_value},
            Confidence.HIGH,
        )

        return adata, decision

    # ── Cell type annotation ──────────────────────────────────────────────

    def _annotate_cell_types(self, experiment_id: str, adata,
                              exp_ctx: dict, intent: dict) -> dict:
        """
        LLM annotation using per-cluster marker genes.
        FIX: uses names[cluster] not names[0].
        """
        markers = {}
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(
                adata, groupby="leiden", method="wilcoxon"
            )
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
                prompt=prompt,
                system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=1000,
            )
            raw = raw.strip()
            if "```" in raw:
                raw = raw.split("```")[1]
                if raw.startswith("json"):
                    raw = raw[4:]
            cell_types = json.loads(raw.strip("` \n"))
        except Exception as e:
            log.warning(f"LLM annotation failed: {e}")
            cell_types = {k: {"cell_type": "annotation_failed",
                               "confidence": "low"}
                          for k in markers}

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

    # ── Differential expression ───────────────────────────────────────────

    def _differential_expression(self, experiment_id: str, adata,
                                   intent: dict, exp_ctx: dict) -> dict:
        """
        Per-cluster DE using Wilcoxon.
        FIX: iterates names[cluster] not names[0].
        """
        de_genes = {}
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(
                adata, groupby="leiden",
                method="wilcoxon",
                key_added="de_wilcoxon",
            )
            rgg = adata.uns["de_wilcoxon"]
            for cl in rgg["names"].dtype.names:
                names = rgg["names"][cl][:50]
                pvals = rgg["pvals_adj"][cl][:50]
                lfc   = rgg["logfoldchanges"][cl][:50]
                sig   = [
                    {"gene": str(g),
                     "log2fc": round(float(l), 3),
                     "padj":   round(float(p), 6)}
                    for g, p, l in zip(names, pvals, lfc)
                    if float(p) < 0.05 and abs(float(l)) > 0.5
                    and str(g) not in ("nan", "")
                ]
                de_genes[str(cl)] = sig[:20]

        except (ImportError, KeyError) as e:
            log.warning(f"DE failed: {e}")
            de_genes = {}

        n_sig = sum(len(v) for v in de_genes.values())
        conf  = (Confidence.HIGH   if n_sig > 50  else
                 Confidence.MEDIUM if n_sig > 10  else
                 Confidence.LOW    if n_sig > 0   else
                 Confidence.INSUFFICIENT)

        # LLM interpretation
        prompt = f"""
Biological question: {intent.get("summary", "")}
Organism: {exp_ctx.get("organism", "")}
Per-cluster DE results (top hits, Wilcoxon):
{json.dumps({k: v[:5] for k, v in de_genes.items()}, indent=2)}
Total significant genes: {n_sig}

Interpret in 3-4 sentences. Focus on biology, not statistics.
What processes or cell states do these genes suggest?
"""
        try:
            interpretation = self.llm.complete(
                prompt=prompt,
                system=SCRNA_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=400,
            )
        except Exception:
            interpretation = f"{n_sig} significant DE genes across clusters."

        self.publish_finding(
            experiment_id,
            {"summary":            interpretation,
             "n_significant_genes": n_sig},
            conf,
        )

        return {
            "de_genes_by_cluster": de_genes,
            "n_significant_genes": n_sig,
            "interpretation":      interpretation,
        }

    # ── Data loading ──────────────────────────────────────────────────────

    def _load_adata(self, qc_output_path: str, raw_files: list):
        """Load AnnData — prefer QC-filtered output, fall back to raw."""
        if qc_output_path and Path(qc_output_path).exists():
            adata = self._load_raw_adata([qc_output_path])
            if adata is not None:
                adata = self._preprocess(adata)
                return adata
        adata = self._load_raw_adata(raw_files)
        if adata is not None:
            adata = self._preprocess(adata)
        return adata

    def _load_raw_adata(self, files: list):
        """Load from h5ad, 10x h5, or MEX directory."""
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
                    return sc.read_10x_mtx(
                        str(p), var_names="gene_symbols", cache=True
                    )
        except ImportError:
            log.warning("scanpy not available")
        except Exception as e:
            log.error(f"Load failed: {e}")
        return None

    def _preprocess(self, adata):
        """Normalize, HVG, PCA, neighbors for ParameterAdvisor."""
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
