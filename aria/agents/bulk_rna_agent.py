"""
ARIA BulkRNAAgent
-----------------
Bulk RNA-seq differential expression and pathway analysis only.

Delegates all computation to aria/scripts/rna_bulk_de.py
running inside aria-rna-env via EnvironmentManager.

Does NOT handle: single-cell, clustering, cell type annotation
                 (→ scRNAAgent)
"""

from __future__ import annotations

import logging
from pathlib import Path

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence
from aria.llm.provider import LLMProvider, TaskTier
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.bulk_rna")


BULK_RNA_SYSTEM = """
You are ARIA's BulkRNAAgent — a specialist in bulk RNA-seq analysis.

Your expertise:
- Experimental design: balanced replicates, confounders, batch effects
- Differential expression: DESeq2 (primary), edgeR, limma-voom
- Pathway enrichment: GO BP, KEGG, Reactome, GSEA
- Quality control: library size normalization, PCA outlier detection
- Interpretation: biological context over statistical significance

Critical knowledge:
- DESeq2 requires integer counts — round if needed
- Design factor must match the biological comparison, not "sample"
- Pseudo-replication invalidates DE results
- Low-count genes must be filtered before DE
- Pathway enrichment requires background gene set matching analysis
- Multiple testing correction: always use adjusted p-values
""".strip()


class BulkRNAAgent(BaseAgent):

    name        = "bulk_rna_agent"
    description = "Bulk RNA-seq: QC, DESeq2, pathway enrichment, plots."

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx    = context.get("exp_context", {})
        intent     = context.get("biological_intent", {})
        modalities = exp_ctx.get("modalities", {})
        files      = modalities.get("bulk_RNA", [])

        if not files:
            return {"status": "failed", "reason": "no_bulk_rna_files"}

        self.publish_status(experiment_id,
                            "BulkRNAAgent starting...", 0.0)

        # Extract design from intent
        design_factor, comparison = self._extract_design(intent, files)

        self.publish_status(experiment_id,
                            f"Running DESeq2 ({design_factor})...", 0.2)

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_bulk_de.py",
            params={
                "files":          files,
                "design_factor":  design_factor,
                "comparison":     comparison,
                "organism":       exp_ctx.get("organism", "Homo sapiens"),
                "genome":         exp_ctx.get("genome", "hg38"),
                "output_dir":     self._output_dir(files),
                "run_pathways":   True,
                "padj_threshold": 0.05,
                "lfc_threshold":  1.0,
            },
        )

        if result.get("status") == "error":
            self.publish_finding(
                experiment_id,
                {"summary": f"Bulk DE failed: "
                            f"{result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "findings": result}

        # Publish findings
        self._publish_findings(experiment_id, result)

        # LLM interpretation
        result["interpretation"] = self._interpret(result, intent, exp_ctx)

        self.publish_status(experiment_id,
                            "BulkRNAAgent complete.", 1.0)
        return {"status": "done", "findings": result}

    # ── Design extraction ─────────────────────────────────────────────────

    def _extract_design(self, intent: dict,
                         files: list) -> tuple[str, dict]:
        """Extract DESeq2 design factor and comparison from intent."""
        import re

        design_factor   = "condition"
        comparison_dict = {}
        comparison_str  = str(intent.get("comparison", "")).lower()

        if any(k in comparison_str for k in
               ["genotype", "knockout", "ko", "wt"]):
            design_factor = "genotype"
        elif any(k in comparison_str for k in
                 ["treatment", "treated", "drug"]):
            design_factor = "treatment"
        elif any(k in comparison_str for k in
                 ["time", "timepoint", "hour", "day",
                  "h vs", "min vs"]):
            design_factor = "timepoint"

        vs = re.search(
            r'([\w\-]+)\s+(?:vs\.?|versus)\s+([\w\-]+)',
            comparison_str
        )
        if vs:
            comparison_dict = {
                "numerator":   vs.group(1),
                "denominator": vs.group(2),
            }

        return design_factor, comparison_dict

    # ── Findings publisher ────────────────────────────────────────────────

    def _publish_findings(self, experiment_id: str, result: dict):
        n_sig  = result.get("n_significant",   0)
        n_up   = result.get("n_upregulated",   0)
        n_down = result.get("n_downregulated", 0)
        comp   = result.get("comparison_used", {})

        conf = (Confidence.HIGH   if n_sig > 100 else
                Confidence.MEDIUM if n_sig > 10  else
                Confidence.LOW    if n_sig > 0   else
                Confidence.INSUFFICIENT)

        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Bulk DE ({comp.get('numerator','?')} vs "
                f"{comp.get('denominator','?')}): "
                f"{n_sig} genes ({n_up} up, {n_down} down)"
            ),
             "n_significant":   n_sig,
             "n_upregulated":   n_up,
             "n_downregulated": n_down,
             "top_genes":       result.get("top_genes", [])[:10],
            },
            conf,
        )

        # QC
        qc = result.get("sample_qc", {})
        if qc:
            outliers = qc.get("outliers", [])
            self.publish_finding(
                experiment_id,
                {"summary": (
                    f"Bulk QC: {qc.get('n_samples','?')} samples. "
                    f"Size range: {qc.get('size_ratio',1):.1f}x."
                    + (f" Outliers removed: {outliers}."
                       if outliers else "")
                )},
                Confidence.HIGH if not outliers else Confidence.MEDIUM,
            )

        # Pathways
        for db, terms in result.get("pathways", {}).items():
            if isinstance(terms, list) and terms:
                self.publish_finding(
                    experiment_id,
                    {"summary": (
                        f"{db}: {len(terms)} pathways. "
                        f"Top: {', '.join(t['term'] for t in terms[:3])}"
                    ),
                     "pathways": terms[:10],
                     "database": db},
                    Confidence.MEDIUM,
                )

    # ── LLM interpretation ────────────────────────────────────────────────

    def _interpret(self, result: dict, intent: dict,
                    exp_ctx: dict) -> str:
        n_sig  = result.get("n_significant", 0)
        tops   = [g.get("gene", g) for g in result.get("top_genes", [])[:8]]
        top_pw = [t["term"]
                  for db, terms in result.get("pathways", {}).items()
                  if isinstance(terms, list)
                  for t in terms[:2]]

        prompt = f"""
Bulk RNA-seq DE: {result.get("comparison_used", {})}
Design: {result.get("design_used", "")}
Organism: {exp_ctx.get("organism", "")}
Question: {intent.get("summary", "")}

Results:
  {n_sig} significant genes (padj < 0.05, |log2FC| > 1)
  {result.get("n_upregulated", 0)} up / {result.get("n_downregulated", 0)} down
  Top genes: {tops}
  Top pathways: {top_pw[:5]}
  Warnings: {result.get("warnings", [])[:3]}

Write a 3-4 sentence biological interpretation.
Include pathway context. Flag data quality concerns.
"""
        try:
            return self.llm.complete(
                prompt=prompt,
                system=BULK_RNA_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=350,
            )
        except Exception:
            return (
                f"{n_sig} DE genes identified. "
                f"Top pathways: {', '.join(top_pw[:3]) if top_pw else 'none'}"
            )

    def _output_dir(self, files: list) -> str:
        if files:
            return str(Path(files[0]).parent / "aria_bulk_de")
        return "/tmp/aria_bulk_de"

    def receive(self, message):
        pass
