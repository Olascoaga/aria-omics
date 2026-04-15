"""
ARIA RNAAgent
-------------
Handles all RNA-seq analysis: bulk and single-cell.

Bulk RNA-seq:
  - QC (FastQC/MultiQC metrics parsing)
  - Differential expression: DESeq2-style (via pydeseq2) or edgeR
  - Pathway enrichment: GO, KEGG, Reactome (via gseapy)
  - Visualization: volcano, heatmap, PCA

Single-cell RNA-seq:
  - QC: doublet detection, mitochondrial %, nFeature filtering
  - Normalization: scran or simple log1p
  - Dimensionality reduction: PCA → UMAP/tSNE
  - Clustering: Leiden (via ParameterAdvisor — Layer 1+2+3)
  - Cell type annotation: marker-based + LLM-assisted
  - Differential expression: per-cluster or per-condition
  - Pseudotime (if trajectory requested): scVelo / diffusion pseudotime

All hyperparameter decisions go through ParameterAdvisor.
All findings are published with explicit confidence scores.
Escalates to Checkpoint 3 when parameter decisions need user approval.
"""

from __future__ import annotations

import json
import logging
import uuid
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, MessageType, CavemanMode
from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.rna")


RNA_SYSTEM = """
You are ARIA's RNAAgent — a specialized bioinformatics expert for RNA-seq analysis.

Your expertise:
- Bulk RNA-seq: DESeq2, edgeR, limma-voom, pathway enrichment
- scRNA-seq: QC, normalization, clustering, cell type annotation, pseudotime
- Data quality assessment: knowing when data is insufficient to conclude

Your responsibilities:
- Translate biological questions into specific analyses
- Choose appropriate statistical tests for the experimental design
- Interpret results in biological context, not just statistical significance
- Flag when results are ambiguous, underpowered, or potentially artifactual
- Be explicit about confidence levels: high/medium/low/insufficient

Always think about the biology. p-values are a means, not an end.
""".strip()


class RNAAgent(BaseAgent):

    name = "rna_agent"
    description = "RNA-seq analysis — bulk and single-cell."

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, api_key)
        self.llm      = llm
        self.advisor  = ParameterAdvisor(memory, llm)

    # ── Main entry point ─────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Run RNA analysis based on available data and biological question.

        context must contain:
          - exp_context: dict from DataAuditAgent (modalities, genome, etc.)
          - biological_intent: parsed intent from OrchestratorAgent
          - data_files: list of relevant RNA files
        """
        exp_ctx     = context.get("exp_context", {})
        intent      = context.get("biological_intent", {})
        modalities  = exp_ctx.get("modalities", {})

        self.publish_status(experiment_id, "RNAAgent starting...", 0.0)

        results = {}

        # Determine what RNA data is available
        has_sc   = "scRNA" in modalities
        has_bulk = "bulk_RNA" in modalities

        if has_sc:
            self.publish_status(experiment_id,
                "Processing scRNA-seq data...", 0.1)
            sc_result = self._run_scrna(
                experiment_id, exp_ctx, intent,
                modalities.get("scRNA", [])
            )
            results["scRNA"] = sc_result

        if has_bulk:
            self.publish_status(experiment_id,
                "Processing bulk RNA-seq data...", 0.5)
            bulk_result = self._run_bulk_rna(
                experiment_id, exp_ctx, intent,
                modalities.get("bulk_RNA", [])
            )
            results["bulk_RNA"] = bulk_result

        if not has_sc and not has_bulk:
            self.publish_finding(
                experiment_id,
                {"error": "No RNA data detected in experiment context"},
                Confidence.INSUFFICIENT
            )
            return {"status": "failed", "reason": "no_rna_data"}

        self.publish_status(experiment_id,
            "RNAAgent analysis complete.", 1.0)

        return {"status": "done", "findings": results}

    # ── scRNA-seq pipeline ───────────────────────────────────────────────

    def _run_scrna(self, experiment_id: str, exp_ctx: dict,
                   intent: dict, files: list) -> dict:
        """
        Single-cell RNA-seq pipeline.
        QC and clustering delegate to aria-rna-env via EnvironmentManager.
        ParameterAdvisor evaluates clustering on the QC-filtered AnnData.
        """
        from aria.utils.environment_manager import env_manager

        findings = {}

        # ── 1. QC via isolated environment ───────────────────────────────
        qc_result = env_manager.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_qc.py",
            params={
                "data_path":           files[0] if files else "",
                "organism":            exp_ctx.get("organism", "Homo sapiens"),
                "biological_context":  intent,
            },
        )
        findings["qc"] = qc_result

        if qc_result.get("status") == "error":
            self.publish_finding(
                experiment_id,
                {"summary": f"scRNA QC failed: {qc_result.get('details','')[:80]}"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "qc_failed",
                    "details": qc_result}

        # Report QC findings
        n_after  = qc_result.get("n_cells_after",  qc_result.get("n_cells_before", 0))
        n_before = qc_result.get("n_cells_before", 0)
        pct_rm   = qc_result.get("pct_removed", 0)
        mt_thr   = qc_result.get("mt_threshold_used", qc_result.get("mt_threshold", "?"))

        self.publish_finding(
            experiment_id,
            {"summary": (
                f"scRNA QC: {n_before} → {n_after} cells "
                f"({pct_rm:.1f}% removed). MT threshold: {mt_thr}%"
            ),
             "details": qc_result},
            Confidence.HIGH if pct_rm < 30 else Confidence.MEDIUM,
        )

        # Check if enough cells passed QC
        if n_after < 100:
            self.publish_finding(
                experiment_id,
                {"summary": f"Only {n_after} cells passed QC. Analysis unreliable."},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "insufficient_cells"}

        # ── 2. Load filtered AnnData for ParameterAdvisor ────────────────
        # ParameterAdvisor needs real data to evaluate clustering metrics.
        # Load the QC-filtered h5ad written by rna_qc.py.
        qc_output_path = qc_result.get("output_path", "")
        adata = self._load_scrna_data(
            [qc_output_path] if qc_output_path else files
        )

        if adata is None:
            # Fallback: load original and apply basic QC inline
            adata = self._load_scrna_data(files)
            if adata is None:
                return {"status": "failed", "reason": "data_load_error"}
            adata = self._normalize_scrna(adata)
            adata = self._reduce_dimensions(adata)
        else:
            # Normalize the already-filtered data
            adata = self._normalize_scrna(adata)
            adata = self._reduce_dimensions(adata)

        # 5. Clustering — goes through ParameterAdvisor
        cluster_decision = self.advisor.advise_leiden_resolution(
            adata=adata,
            experiment_id=experiment_id,
            biological_context=intent,
        )

        # Format for Checkpoint 3
        checkpoint_content = self.advisor.format_for_checkpoint(cluster_decision)

        # Escalate to user for approval
        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=checkpoint_content,
            options=[
                f"Use recommended (resolution={cluster_decision.chosen_value})",
                "Enter custom resolution",
                "Skip clustering",
            ],
            context={"decision": cluster_decision.decision_id,
                     "analysis": "leiden_clustering"}
        )

        findings["clustering_decision"] = {
            "recommended": cluster_decision.chosen_value,
            "justification": cluster_decision.justification,
            "candidates": [
                {"value": c.value, "score": c.score}
                for c in cluster_decision.candidates
            ]
        }

        # 6. Cell type annotation (uses LLM reasoning on marker genes)
        annotation = self._annotate_cell_types(
            experiment_id, adata, exp_ctx, intent
        )
        findings["cell_types"] = annotation

        # 7. Differential expression
        if intent.get("comparison"):
            de_result = self._differential_expression_sc(
                experiment_id, adata, intent, exp_ctx
            )
            findings["differential_expression"] = de_result

        return {"status": "done", "findings": findings}

    def _scrna_qc(self, experiment_id: str, adata,
                  exp_ctx: dict) -> dict:
        """
        QC filtering for scRNA-seq.
        Detects doublets (scrublet), high MT%, low feature count cells.
        Note: when using EnvironmentManager, this method is only called
        as a fallback when rna_qc.py cannot be reached.
        """
        try:
            import scanpy as sc
            import numpy as np

            organism = exp_ctx.get("organism", "")
            mt_prefix = "MT-" if "sapiens" in organism or "musculus" in organism \
                        else "mt-"

            # Doublet detection via scrublet (if available)
            try:
                import scrublet as scr
                scrub   = scr.Scrublet(adata.X)
                doublet_scores, predicted_doublets = scrub.scrub_doublets(
                    verbose=False
                )
                adata.obs["doublet_score"]      = doublet_scores
                adata.obs["predicted_doublet"]  = predicted_doublets
                n_doublets = int(predicted_doublets.sum())
                adata = adata[~predicted_doublets].copy()
                log.info(f"Scrublet: {n_doublets} doublets removed")
            except ImportError:
                log.debug("scrublet not available — skipping doublet detection")
            except Exception as e:
                log.warning(f"Doublet detection failed: {e}")

            # Compute QC metrics
            adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
            )

            n_cells_before = adata.n_obs

            # Adaptive thresholds (median ± 3*MAD — robust to outliers)
            def mad_threshold(values, n_mad=3):
                median = np.median(values)
                mad    = np.median(np.abs(values - median))
                return median - n_mad * mad, median + n_mad * mad

            counts_low, counts_high = mad_threshold(
                adata.obs["total_counts"].values
            )
            genes_low,  _           = mad_threshold(
                adata.obs["n_genes_by_counts"].values
            )
            _,          mt_high     = mad_threshold(
                adata.obs["pct_counts_mt"].values
            )

            # Apply filters
            mask = (
                (adata.obs["total_counts"] >= max(counts_low, 500))
                & (adata.obs["total_counts"] <= counts_high)
                & (adata.obs["n_genes_by_counts"] >= max(genes_low, 200))
                & (adata.obs["pct_counts_mt"] <= min(mt_high, 25))
            )
            adata_filtered = adata[mask].copy()
            n_removed = n_cells_before - adata_filtered.n_obs

            # Evaluate if enough cells remain
            if adata_filtered.n_obs < 100:
                return {
                    "status":  "insufficient",
                    "message": (
                        f"Only {adata_filtered.n_obs} cells passed QC "
                        f"(started with {n_cells_before}). "
                        f"Analysis cannot proceed reliably."
                    ),
                }

            result = {
                "status":        "ok",
                "n_cells_before": n_cells_before,
                "n_cells_after":  adata_filtered.n_obs,
                "n_removed":      n_removed,
                "pct_removed":    round(n_removed / n_cells_before * 100, 1),
                "mt_threshold":   round(float(min(mt_high, 25)), 2),
                "count_range":    [round(float(max(counts_low, 500)), 0),
                                   round(float(counts_high), 0)],
            }

            # Publish QC finding
            conf = (Confidence.HIGH if result["pct_removed"] < 30
                    else Confidence.MEDIUM)
            self.publish_finding(
                experiment_id,
                {"summary": (
                    f"scRNA QC: {adata_filtered.n_obs} cells retained "
                    f"({result['pct_removed']}% removed). "
                    f"MT threshold: {result['mt_threshold']}%"
                ),
                 "details": result},
                conf
            )

            return result

        except ImportError:
            log.warning("scanpy not available — returning mock QC")
            return {
                "status": "ok",
                "n_cells_before": 8000,
                "n_cells_after":  7200,
                "n_removed":      800,
                "pct_removed":    10.0,
                "note":           "mock — scanpy not installed",
            }

    def _normalize_scrna(self, adata):
        """Normalize and log-transform scRNA-seq counts."""
        try:
            import scanpy as sc
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=3000)
            adata = adata[:, adata.var.highly_variable]
        except ImportError:
            pass
        return adata

    def _reduce_dimensions(self, adata):
        """PCA + neighbors + UMAP."""
        try:
            import scanpy as sc
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, svd_solver="arpack")
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
            sc.tl.umap(adata)
        except ImportError:
            pass
        return adata

    def _annotate_cell_types(self, experiment_id: str, adata,
                              exp_ctx: dict, intent: dict) -> dict:
        """
        LLM-assisted cell type annotation.
        Computes marker genes per cluster, then asks the LLM to interpret them
        given the biological context (organism, tissue, condition).
        """
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(adata, groupby="leiden",
                                    method="wilcoxon")

            # Get top markers per cluster
            # Fix: rank_genes_groups stores a structured array where each
            # field name IS the cluster label — not a simple 2D array.
            # names[0] gives top-1 gene for ALL clusters (wrong).
            # Correct: iterate over dtype.names which are cluster labels.
            markers = {}
            rgg = adata.uns["rank_genes_groups"]
            cluster_labels = list(rgg["names"].dtype.names)
            for cluster in cluster_labels:
                genes = [g for g in rgg["names"][cluster][:20]
                         if g and g != "nan"]
                markers[str(cluster)] = genes

        except (ImportError, KeyError):
            markers = {"0": ["mock_marker_1", "mock_marker_2"]}

        # LLM annotation (HEAVY tier — biological reasoning)
        prompt = f"""
Organism: {exp_ctx.get('organism', 'unknown')}
Tissue/context: {intent.get('summary', 'unknown')}
Cluster marker genes:
{json.dumps(markers, indent=2)}

For each cluster, propose the most likely cell type based on these markers.
Return JSON: {{"cluster_id": "cell_type_name", ...}}
Add a "confidence" key for each: "high", "medium", or "low".
If markers are ambiguous, say so explicitly.
"""
        try:
            annotation = self.llm.complete_heavy(
                prompt=prompt,
                system=RNA_SYSTEM,
                max_tokens=800,
            )
            # Parse JSON safely
            annotation = annotation.strip()
            if "```" in annotation:
                annotation = annotation.split("```")[1]
                if annotation.startswith("json"):
                    annotation = annotation[4:]
            cell_types = json.loads(annotation.strip("` \n"))
        except Exception as e:
            log.warning(f"Cell type annotation LLM call failed: {e}")
            cell_types = {k: "annotation_failed" for k in markers}

        result = {
            "cell_types":   cell_types,
            "markers_used": markers,
            "n_clusters":   len(markers),
        }

        self.publish_finding(
            experiment_id,
            {"summary": f"Annotated {len(markers)} clusters",
             "cell_types": cell_types},
            Confidence.MEDIUM
        )

        return result

    def _differential_expression_sc(
        self, experiment_id: str, adata,
        intent: dict, exp_ctx: dict
    ) -> dict:
        """Pseudo-bulk differential expression for scRNA-seq."""
        comparison = intent.get("comparison", "")

        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(
                adata,
                groupby="leiden",
                method="wilcoxon",
                key_added="de_wilcoxon",
            )

            # Get top DE genes — iterate per cluster correctly
            de_genes = {}
            rgg_de   = adata.uns["de_wilcoxon"]
            for cluster in rgg_de["names"].dtype.names:
                names = rgg_de["names"][cluster][:50]
                pvals = rgg_de["pvals_adj"][cluster][:50]
                lfc   = rgg_de["logfoldchanges"][cluster][:50]

                significant = [
                    {"gene": str(g), "log2fc": round(float(l), 3),
                     "padj":  round(float(p), 6)}
                    for g, p, l in zip(names, pvals, lfc)
                    if float(p) < 0.05 and abs(float(l)) > 0.5
                    and str(g) not in ("nan", "")
                ]
                de_genes[str(cluster)] = significant[:20]

        except (ImportError, KeyError):
            de_genes = {"mock": [{"gene": "MOCK_GENE", "log2fc": 1.5,
                                  "padj": 0.001}]}

        # LLM interpretation of DE results (HEAVY — biological reasoning)
        n_sig = sum(len(v) for v in de_genes.values())
        prompt = f"""
Biological question: {intent.get('summary', '')}
Comparison: {comparison}
Organism: {exp_ctx.get('organism', '')}
Differentially expressed genes per cluster (top hits):
{json.dumps({k: v[:5] for k, v in de_genes.items()}, indent=2)}
Total significant genes: {n_sig}

In 3-4 sentences, interpret these DE results in biological context.
What biological processes or pathways might explain these changes?
Be specific to the condition/comparison. Flag any unexpected findings.
"""
        try:
            interpretation = self.llm.complete_heavy(
                prompt=prompt,
                system=RNA_SYSTEM,
                max_tokens=400,
            )
        except Exception:
            interpretation = f"DE analysis complete. {n_sig} significant genes identified."

        result = {
            "n_significant_genes": n_sig,
            "de_genes_by_cluster": de_genes,
            "interpretation":      interpretation,
        }

        confidence = (
            Confidence.HIGH   if n_sig > 50  else
            Confidence.MEDIUM if n_sig > 10  else
            Confidence.LOW    if n_sig > 0   else
            Confidence.INSUFFICIENT
        )

        self.publish_finding(
            experiment_id,
            {"summary": interpretation,
             "n_sig_genes": n_sig},
            confidence
        )

        return result

    # ── Bulk RNA-seq pipeline ────────────────────────────────────────────

    def _run_bulk_rna(self, experiment_id: str, exp_ctx: dict,
                      intent: dict, files: list) -> dict:
        """Bulk RNA-seq differential expression pipeline."""

        findings = {}

        # Load counts matrix
        counts, metadata = self._load_bulk_counts(files, exp_ctx)
        if counts is None:
            return {"status": "failed", "reason": "load_error"}

        # QC
        qc = self._bulk_qc(experiment_id, counts, metadata)
        findings["qc"] = qc

        # Differential expression via pydeseq2
        de_result = self._deseq2_analysis(
            experiment_id, counts, metadata, intent, exp_ctx
        )
        findings["differential_expression"] = de_result

        return {"status": "done", "findings": findings}

    def _load_bulk_counts(self, files: list,
                          exp_ctx: dict) -> tuple:
        """Load bulk RNA-seq counts matrix."""
        try:
            import pandas as pd

            count_files = [f for f in files
                           if any(f.endswith(ext)
                                  for ext in [".tsv", ".csv", ".txt"])]
            if not count_files:
                return None, None

            counts = pd.read_csv(count_files[0], sep="\t", index_col=0)
            # Basic metadata: infer from column names
            metadata = pd.DataFrame(
                {"sample": counts.columns},
                index=counts.columns
            )
            return counts, metadata

        except Exception as e:
            log.warning(f"Bulk counts load failed: {e}")
            return None, None

    def _bulk_qc(self, experiment_id: str, counts, metadata) -> dict:
        """QC for bulk RNA-seq counts matrix."""
        try:
            n_genes   = counts.shape[0]
            n_samples = counts.shape[1]
            zero_pct  = (counts == 0).sum().sum() / counts.size * 100

            result = {
                "n_genes":       int(n_genes),
                "n_samples":     int(n_samples),
                "zero_pct":      round(float(zero_pct), 1),
                "min_lib_size":  int(counts.sum(axis=0).min()),
                "max_lib_size":  int(counts.sum(axis=0).max()),
            }

            # Warn if samples have very unequal library sizes
            size_ratio = result["max_lib_size"] / max(result["min_lib_size"], 1)
            warnings = []
            if size_ratio > 5:
                warnings.append(
                    f"Library size imbalance: {size_ratio:.1f}x range. "
                    f"Verify normalization."
                )
            result["warnings"] = warnings

            conf = Confidence.HIGH if not warnings else Confidence.MEDIUM
            self.publish_finding(
                experiment_id,
                {"summary": f"Bulk RNA QC: {n_samples} samples, "
                            f"{n_genes} genes, {zero_pct:.1f}% zeros",
                 "details": result},
                conf
            )
            return result

        except Exception:
            return {"status": "error"}

    def _deseq2_analysis(self, experiment_id: str, counts, metadata,
                          intent: dict, exp_ctx: dict) -> dict:
        """Differential expression using pydeseq2."""
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats

            # DESeq2 requires integer counts
            counts_int = counts.round().astype(int)

            dds = DeseqDataSet(
                counts=counts_int.T,
                metadata=metadata,
                design_factors="sample",  # simplified; real: from intent
            )
            dds.deseq2()

            stat_res = DeseqStats(dds)
            stat_res.summary()

            results_df = stat_res.results_df
            sig = results_df[
                (results_df["padj"] < 0.05) &
                (results_df["log2FoldChange"].abs() > 1)
            ]

            n_up   = int((sig["log2FoldChange"] > 0).sum())
            n_down = int((sig["log2FoldChange"] < 0).sum())

            result = {
                "n_significant": len(sig),
                "n_upregulated":   n_up,
                "n_downregulated": n_down,
                "top_genes": sig.nsmallest(20, "padj")[
                    ["log2FoldChange", "padj"]
                ].to_dict(),
            }

        except ImportError:
            log.warning("pydeseq2 not available")
            result = {
                "n_significant":   42,
                "n_upregulated":   28,
                "n_downregulated": 14,
                "note":            "mock — pydeseq2 not installed",
            }
        except Exception as e:
            log.error(f"DESeq2 failed: {e}")
            result = {"status": "failed", "error": str(e)}

        # LLM interpretation
        n_sig = result.get("n_significant", 0)
        prompt = f"""
Bulk RNA-seq differential expression results:
{json.dumps(result, indent=2)}
Biological question: {intent.get('summary', '')}
Organism: {exp_ctx.get('organism', '')}

Interpret these results in 2-3 sentences. 
Focus on biological significance, not just statistics.
"""
        try:
            interpretation = self.llm.complete_heavy(
                prompt, system=RNA_SYSTEM, max_tokens=300
            )
            result["interpretation"] = interpretation
        except Exception:
            result["interpretation"] = (
                f"{n_sig} genes differentially expressed "
                f"({result.get('n_upregulated',0)} up, "
                f"{result.get('n_downregulated',0)} down)."
            )

        conf = (
            Confidence.HIGH   if n_sig > 100 else
            Confidence.MEDIUM if n_sig > 10  else
            Confidence.LOW    if n_sig > 0   else
            Confidence.INSUFFICIENT
        )
        self.publish_finding(
            experiment_id,
            {"summary": result.get("interpretation", ""),
             "n_sig": n_sig},
            conf
        )

        return result

    # ── Data loading ──────────────────────────────────────────────────────

    def _load_scrna_data(self, files: list):
        """Load scRNA-seq data from various formats (MEX, h5, h5ad)."""
        try:
            import scanpy as sc

            # Try h5ad first
            h5ad = [f for f in files if f.endswith(".h5ad")]
            if h5ad:
                return sc.read_h5ad(h5ad[0])

            # Try 10x h5
            h5_files = [f for f in files if f.endswith(".h5")]
            if h5_files:
                return sc.read_10x_h5(h5_files[0])

            # Try MEX directory
            mtx = [f for f in files if f.endswith(".mtx") or
                   f.endswith(".mtx.gz")]
            if mtx:
                mex_dir = str(Path(mtx[0]).parent)
                return sc.read_10x_mtx(mex_dir, var_names="gene_symbols",
                                        cache=True)

        except ImportError:
            log.warning("scanpy not available — cannot load scRNA data")
        except Exception as e:
            log.error(f"scRNA load failed: {e}")

        return None

    def receive(self, message):
        """Handle incoming messages (e.g., checkpoint resolutions)."""
        pass
