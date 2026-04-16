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



def _is_fastq(files: list) -> bool:
    """Detect if files are FASTQs (raw reads) vs count matrices."""
    if not files:
        return False
    ext = str(files[0]).lower()
    return any(ext.endswith(s) for s in
               [".fastq.gz", ".fq.gz", ".fastq", ".fq"])

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

        self.publish_status(experiment_id, "BulkRNAAgent starting...", 0.0)

        # ── Detect if files are FASTQs (raw) or count matrices ───────────
        is_raw = _is_fastq(files)

        if is_raw:
            # Full pipeline: FASTQ → trim → align → quantify → DE
            counts_files, preprocessing = self._run_preprocessing(
                experiment_id, files, exp_ctx, intent
            )
            if counts_files is None:
                return {"status": "failed",
                        "findings": preprocessing,
                        "reason": "preprocessing_failed"}
            # Hand off to DE with the generated counts matrix
            files = counts_files
        else:
            preprocessing = None

        # ── Differential expression ───────────────────────────────────────
        design_factor, comparison = self._extract_design(intent, files)
        self.publish_status(experiment_id,
                            f"Running DESeq2 ({design_factor})...", 0.7)

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
                {"summary": f"Bulk DE failed: {result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "findings": result}

        if preprocessing:
            result["preprocessing"] = preprocessing

        self._publish_findings(experiment_id, result)
        result["interpretation"] = self._interpret(result, intent, exp_ctx)

        self.publish_status(experiment_id, "BulkRNAAgent complete.", 1.0)
        return {"status": "done", "findings": result}

    # ── FASTQ preprocessing pipeline ──────────────────────────────────────

    def _run_preprocessing(self, experiment_id: str, fastq_files: list,
                            exp_ctx: dict, intent: dict) -> tuple:
        """
        Full pipeline: FASTQ → fastp → STAR → featureCounts → counts matrix.
        Returns (counts_file_list, preprocessing_summary) or (None, error).
        """
        from pathlib import Path

        fastq_dir  = str(Path(fastq_files[0]).parent)
        output_dir = str(Path(fastq_files[0]).parent.parent / "aria_processing")
        organism   = exp_ctx.get("organism", "Homo sapiens")
        genome_cfg = exp_ctx.get("genome_config", {})

        # Step 1: FASTQ QC + trimming
        self.publish_status(experiment_id, "Trimming reads (fastp)...", 0.05)
        qc_result = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_fastq_qc.py",
            params={
                "fastq_dir":  fastq_dir,
                "output_dir": str(Path(output_dir) / "qc"),
                "threads":    8,
            },
        )

        if qc_result.get("status") == "error":
            self.publish_finding(
                experiment_id,
                {"summary": f"FASTQ QC failed: {qc_result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return None, qc_result

        # Publish QC findings
        self._publish_fastq_qc_findings(experiment_id, qc_result)
        self.publish_status(
            experiment_id,
            f"QC complete: {qc_result.get('n_samples',0)} samples trimmed",
            0.20
        )

        # Step 2: Alignment with STAR
        self.publish_status(experiment_id, "Aligning to genome (STAR)...", 0.25)
        align_result = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_align.py",
            params={
                "samples":       qc_result.get("samples", []),
                "genome_dir":    genome_cfg.get("star_index", ""),
                "genome_fasta":  genome_cfg.get("fasta", ""),
                "gtf_file":      genome_cfg.get("gtf", ""),
                "output_dir":    str(Path(output_dir) / "aligned"),
                "threads":       8,
                "two_pass":      True,
            },
        )

        if align_result.get("status") == "error":
            self.publish_finding(
                experiment_id,
                {"summary": f"Alignment failed: {align_result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return None, align_result

        # Publish alignment findings
        self._publish_alignment_findings(experiment_id, align_result)
        n_aligned = align_result.get("n_aligned", 0)
        self.publish_status(
            experiment_id,
            f"Alignment complete: {n_aligned} samples mapped",
            0.55
        )

        # Step 3: Quantification with featureCounts
        self.publish_status(experiment_id, "Counting reads (featureCounts)...", 0.60)
        quant_result = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_quantify.py",
            params={
                "bam_files":  align_result.get("bam_files", []),
                "gtf_file":   genome_cfg.get("gtf", ""),
                "output_dir": str(Path(output_dir) / "counts"),
                "threads":    8,
                "paired":     True,
                "strand":     genome_cfg.get("strand", 0),
            },
        )

        if quant_result.get("status") == "error":
            self.publish_finding(
                experiment_id,
                {"summary": f"Quantification failed: {quant_result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return None, quant_result

        counts_path = quant_result.get("counts_matrix")
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Quantification complete: {quant_result.get('n_genes',0):,} genes × "
                f"{quant_result.get('n_samples',0)} samples"
            )},
            Confidence.HIGH,
        )
        self.publish_status(experiment_id, "Preprocessing complete.", 0.65)

        preprocessing = {
            "qc":         qc_result,
            "alignment":  align_result,
            "quantification": quant_result,
        }

        return [counts_path], preprocessing

    def _publish_fastq_qc_findings(self, experiment_id: str, qc: dict):
        samples = qc.get("samples", [])
        if not samples:
            return
        avg_pass = sum(s.get("pct_passed", 0) for s in samples) / len(samples)
        low_qual = [s["name"] for s in samples
                    if s.get("pct_passed", 100) < 80]
        conf = Confidence.HIGH if not low_qual else Confidence.MEDIUM
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"FASTQ QC: {len(samples)} samples trimmed. "
                f"Avg {avg_pass:.1f}% reads passed."
                + (f" Low quality: {low_qual}" if low_qual else "")
            ),
             "multiqc": qc.get("multiqc_report")},
            conf,
        )

    def _publish_alignment_findings(self, experiment_id: str, align: dict):
        bams = align.get("bam_files", [])
        if not bams:
            return
        ok_bams   = [b for b in bams if b.get("status") == "success"]
        avg_map   = sum(b.get("pct_unique", 0) for b in ok_bams) / max(len(ok_bams), 1)
        low_map   = [b["name"] for b in ok_bams if b.get("pct_unique", 100) < 70]
        conf = Confidence.HIGH if avg_map > 75 and not low_map else Confidence.MEDIUM
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"STAR alignment: {len(ok_bams)}/{len(bams)} samples mapped. "
                f"Avg unique mapping: {avg_map:.1f}%."
                + (f" Low mapping: {low_map}" if low_map else "")
            )},
            conf,
        )

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
