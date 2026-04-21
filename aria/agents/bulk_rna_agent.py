"""
ARIA BulkRNAAgent (v3.1)
----------------------
Bulk RNA-seq differential expression and pathway analysis.

v3.1 CHANGES (critical parameter injection fix):
  - Agent now properly pulls `global_padj` and `global_lfc` from exp_context
    (injected via Orchestrator Checkpoint 3) rather than hardcoding 0.05.
"""

from __future__ import annotations

import gzip
import logging
import re
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
- TF knockouts produce modest direct effects (log2FC 0.5-1); use lower
  LFC thresholds for direct target discovery
- Multiple testing correction: always use adjusted p-values
""".strip()

KNOWN_TFS = {
    "bmal1", "arntl", "clock", "per1", "per2", "per3", "cry1", "cry2", "nr1d1", "reverba", "nr1d2",
    "rora", "rorb", "rorc", "oct4", "pou5f1", "sox2", "nanog", "klf4", "lin28", "gata1", "gata2", 
    "gata3", "gata4", "gata6", "tbx5", "runx1", "runx2", "runx3", "foxa1", "foxa2", "foxp3", 
    "foxo1", "foxo3", "cebpa", "cebpb", "pparg", "ppara", "rela", "relb", "stat1", "stat3", 
    "stat5", "tp53", "trp53", "myc", "max", "sp1", "e2f1", "hif1a", "arnt", "yap1", "wwtr1",
}

def _is_fastq(files: list) -> bool:
    if not files: return False
    return any(str(files[0]).lower().endswith(s) for s in [".fastq.gz", ".fq.gz", ".fastq", ".fq"])

def _infer_lfc_threshold(intent: dict) -> float:
    entities = [str(e).lower() for e in intent.get("biological_entities", [])]
    text = re.sub(r'[\-\s]', '', " ".join([str(intent.get("summary", "")), str(intent.get("comparison", "")), *entities]).lower())
    if any(tf in text for tf in KNOWN_TFS) or re.search(r"\b(knockout|knockdown|ko|kd|overexpression|oe)\b", text) or "transcription factor" in text:
        return 0.58
    return 1.0


class BulkRNAAgent(BaseAgent):

    name        = "bulk_rna_agent"
    description = "Bulk RNA-seq: QC, DESeq2, pathway enrichment, plots."

    def __init__(self, memory: ARIAMemory, llm: LLMProvider, api_key: str = None):
        super().__init__(memory, llm, api_key)
        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx    = context.get("exp_context", {})
        intent     = context.get("biological_intent", {})
        modalities = exp_ctx.get("modalities", {})
        files      = modalities.get("bulk_RNA", [])

        if not files:
            return {"status": "failed", "reason": "no_bulk_rna_files"}

        self.publish_status(experiment_id, "BulkRNAAgent starting...", 0.0)

        if _is_fastq(files):
            counts_files, preprocessing = self._run_preprocessing(experiment_id, files, exp_ctx, intent)
            if counts_files is None:
                return {"status": "failed", "findings": preprocessing, "reason": "preprocessing_failed"}
            files = counts_files
        else:
            preprocessing = None

        sample_names, group_labels = self._discover_groups(files)

        if not group_labels:
            self.publish_finding(experiment_id, {"summary": "Could not infer experimental groups from sample names."}, Confidence.INSUFFICIENT)
            return {"status": "failed", "reason": "group_inference_failed"}

        self.publish_status(experiment_id, f"Detected groups: {list(group_labels.keys())}", 0.68)
        design_factor, contrasts = self._build_contrasts(intent, group_labels, experiment_id)
        
        # BUG FIX: Properly ingest modified thresholds from Orchestrator Checkpoint 3
        padj_thr = exp_ctx.get("global_padj", 0.05)
        lfc_thr  = exp_ctx.get("global_lfc", _infer_lfc_threshold(intent))

        self.publish_status(
            experiment_id,
            f"Running {len(contrasts)} contrast(s), padj < {padj_thr}, |log2FC| > {lfc_thr}...",
            0.70,
        )

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_bulk_de.py",
            params={
                "files":          files,
                "design_factor":  design_factor,
                "contrasts":      contrasts,
                "organism":       exp_ctx.get("organism", "Homo sapiens"),
                "genome":         exp_ctx.get("genome", "hg38"),
                "output_dir":     self._output_dir(files),
                "run_pathways":   True,
                "padj_threshold": padj_thr,
                "lfc_threshold":  lfc_thr,
            },
        )

        if result.get("status") == "error":
            self.publish_finding(experiment_id, {"summary": f"Bulk DE failed: {result.get('details','')[:100]}"}, Confidence.INSUFFICIENT)
            return {"status": "failed", "findings": result}

        if preprocessing:
            result["preprocessing"] = preprocessing
            
        # Ensure thresholds are passed to interpretation and memory 
        result["padj_threshold"] = padj_thr
        result["lfc_threshold"]  = lfc_thr

        self._publish_findings(experiment_id, result)
        result["interpretation"] = self._interpret(result, intent, exp_ctx)
        self._record_methodology_decisions(experiment_id, result)

        self.publish_status(experiment_id, "BulkRNAAgent complete.", 1.0)
        return {"status": "done", "findings": result}

    def _record_methodology_decisions(self, experiment_id: str, result: dict) -> None:
        try:
            methodology = (result.get("methodology") or {})
            decisions   = methodology.get("decisions", []) or []
            stage_to_cp = {"Differential expression (DESeq2)": 1, "PCA + MDS (sample-level structure)": 2, "Heatmap (padj top 50)": 3, "Heatmap (|log2FC| top 50)": 4, "Pathway enrichment (ORA)": 5, "GSEA (pre-ranked)": 6, "TPM (supplementary export)": 7}

            for d in decisions:
                step, cp = d.get("step", ""), stage_to_cp.get(d.get("step", ""), 0)
                decision_summary = f"{d.get('input','?')} | {d.get('normalization','?')} | {d.get('gene_filter','?')}"
                try:
                    self.memory.store_decision(
                        decision_id=f"{experiment_id[:8]}-auto-{cp:02d}", wing_id=experiment_id, checkpoint=cp,
                        question=step, decision=decision_summary, rationale=d.get("justification", "")[:500], made_by="bulk_rna_agent (auto)"
                    )
                except Exception as e: log.debug(f"Failed to store decision '{step}': {e}")

            try:
                self.memory.store_decision(
                    decision_id=f"{experiment_id[:8]}-auto-00-thr", wing_id=experiment_id, checkpoint=0,
                    question="Statistical thresholds for DE significance",
                    decision=f"padj < {result.get('padj_threshold')}, |log2FC| > {result.get('lfc_threshold')}",
                    rationale="Thresholds enforced via User Checkpoint 3 profile selection or inferred by Agent based on TF/KO targets.",
                    made_by="User / bulk_rna_agent"
                )
            except Exception as e: log.debug(f"Failed to store threshold decision: {e}")

            try:
                self.memory.store_decision(
                    decision_id=f"{experiment_id[:8]}-auto-00-design", wing_id=experiment_id, checkpoint=0,
                    question="DESeq2 design formula", decision=result.get("design_used", "~condition"),
                    rationale="Single-factor design inferred from sample labels. No batch or covariate adjustment applied.",
                    made_by="bulk_rna_agent (auto)"
                )
            except Exception as e: log.debug(f"Failed to store design decision: {e}")
        except Exception as e:
            log.warning(f"Decision logging failed (non-fatal): {e}")

    def _run_preprocessing(self, experiment_id: str, fastq_files: list, exp_ctx: dict, intent: dict) -> tuple:
        fastq_dir, output_dir, genome_cfg = str(Path(fastq_files[0]).parent), str(Path(fastq_files[0]).parent.parent / "aria_processing"), exp_ctx.get("genome_config", {})

        self.publish_status(experiment_id, "Trimming reads (fastp)...", 0.05)
        qc_result = self.env.run_in_stack("rnaseq", "aria/scripts/rna_fastq_qc.py", {"fastq_dir": fastq_dir, "output_dir": str(Path(output_dir) / "qc"), "threads": 8})
        if qc_result.get("status") == "error": return None, qc_result

        self._publish_fastq_qc_findings(experiment_id, qc_result)
        self.publish_status(experiment_id, f"QC complete: {qc_result.get('n_samples',0)} samples trimmed", 0.20)

        self.publish_status(experiment_id, "Aligning to genome (STAR)...", 0.25)
        align_result = self.env.run_in_stack("rnaseq", "aria/scripts/rna_align.py", {"samples": qc_result.get("samples", []), "genome_dir": genome_cfg.get("star_index", ""), "genome_fasta": genome_cfg.get("fasta", ""), "gtf_file": genome_cfg.get("gtf", ""), "output_dir": str(Path(output_dir) / "aligned"), "threads": 8, "two_pass": True})
        if align_result.get("status") == "error": return None, align_result

        self._publish_alignment_findings(experiment_id, align_result)
        self.publish_status(experiment_id, f"Alignment complete: {align_result.get('n_aligned', 0)} samples mapped", 0.55)

        self.publish_status(experiment_id, "Counting reads (featureCounts)...", 0.60)
        quant_result = self.env.run_in_stack("rnaseq", "aria/scripts/rna_quantify.py", {"bam_files": align_result.get("bam_files", []), "gtf_file": genome_cfg.get("gtf", ""), "output_dir": str(Path(output_dir) / "counts"), "threads": 8, "paired": True, "strand": genome_cfg.get("strand", "auto")})
        if quant_result.get("status") == "error": return None, quant_result

        self.publish_finding(experiment_id, {"summary": f"Quantification complete: {quant_result.get('n_genes',0):,} genes × {quant_result.get('n_samples',0)} samples"}, Confidence.HIGH)
        return [quant_result.get("counts_matrix")], {"qc": qc_result, "alignment": align_result, "quantification": quant_result}

    def _discover_groups(self, files: list) -> tuple[list, dict]:
        if not files: return [], {}
        sample_names = self._read_sample_names(files[0])
        if not sample_names: return [], {}
        groups = self._infer_groups_local(sample_names)
        if not groups: return sample_names, {}
        by_group = {}
        for sample, label in groups.items(): by_group.setdefault(label, []).append(sample)
        return sample_names, by_group

    @staticmethod
    def _read_sample_names(path: str) -> list[str]:
        try:
            p = Path(path)
            if not p.exists(): return []
            with (gzip.open if str(p).endswith(".gz") else open)(p, "rt") as f: header = f.readline().rstrip("\n")
            return [c for c in header.split("\t" if "\t" in header else ",") if c.lower().strip() not in {"gene_id", "geneid", "", "chr", "start", "end", "strand", "length", "gene_name", "symbol", "feature", "ensembl", "ensembl_id", "entrez", "entrez_id"}]
        except Exception as e: return []

    @staticmethod
    def _infer_groups_local(samples: list[str]) -> dict:
        for pattern in [r'^([A-Za-z][A-Za-z0-9]+)[_\-](\d+)$', r'^([A-Za-z][A-Za-z0-9]+)[_\-]([Rr]ep\d+|[A-Za-z]\d*)$', r'^([A-Za-z][A-Za-z0-9\-]*?)(\d.*)$']:
            m = {s: match.group(1).rstrip("_-") for s in samples if (match := re.compile(pattern).match(s))}
            if len(m) == len(samples) and len(set(m.values())) >= 2: return m
        if all("_" in s for s in samples):
            gr = {s: "_".join(s.split("_")[:-1]) for s in samples}
            if len(set(gr.values())) >= 2: return gr
        return {}

    def _build_contrasts(self, intent: dict, group_labels: dict, experiment_id: str) -> tuple[str, list]:
        design_factor, group_names = self._infer_design_factor(intent), list(group_labels.keys())
        entity_to_label = self._map_entities_to_labels(intent.get("biological_entities", []), group_names, intent)
        control, contrasts = self._identify_control(group_names), []
        
        if control:
            contrasts = [{"numerator": l, "denominator": control, "name": self._humanize_contrast(l, control, entity_to_label)} for l in group_names if l != control]
        else:
            sorted_g = sorted(group_names)
            contrasts = [{"numerator": sorted_g[i], "denominator": b, "name": f"{sorted_g[i]} vs {b}"} for i in range(len(sorted_g)) for b in sorted_g[i+1:]]

        if entity_to_label: self.publish_finding(experiment_id, {"summary": f"Entity-to-label mapping: {', '.join(f'{e}→{l}' for e,l in entity_to_label.items())}. Contrasts: {[c['name'] for c in contrasts]}"}, Confidence.HIGH)
        return design_factor, contrasts

    @staticmethod
    def _infer_design_factor(intent: dict) -> str:
        text = f"{intent.get('comparison', '')} {intent.get('summary', '')}".lower()
        if any(k in text for k in ["knockout", "ko ", "knockdown", "kd ", "genotype", "mutant", "wt ", "wildtype"]): return "genotype"
        if any(k in text for k in ["treat", "drug", "vehicle", "dmso"]): return "treatment"
        if any(k in text for k in ["time", "timepoint", "hour", "day", "min"]): return "timepoint"
        return "condition"

    def _map_entities_to_labels(self, entities: list, group_names: list, intent: dict) -> dict:
        mapping, used_labels, ent_clean = {}, set(), [(str(e), re.sub(r'[^a-z0-9]', '', str(e).lower())) for e in entities if re.sub(r'[^a-z0-9]', '', str(e).lower()) not in ("cells", "cell", "h9", "h1", "hesc")]
        for label in sorted(group_names, key=len, reverse=True):
            if label in used_labels: continue
            for original, norm in ent_clean:
                if original not in mapping and norm.startswith(label.lower()) and len(label.lower()) >= 1:
                    mapping[original] = label
                    used_labels.add(label)
                    break
        for original, norm in ent_clean:
            if original not in mapping and any(w in norm for w in {"wt", "wildtype", "control", "ctrl", "untreated"}):
                for label in group_names:
                    if label not in used_labels and label.lower() in {"wt", "wildtype", "control", "ctrl", "untreated"}:
                        mapping[original] = label
                        used_labels.add(label)
                        break
        if len(mapping) >= max(1, len(ent_clean) // 2): return mapping
        try: return self._llm_match_labels(entities, group_names, intent) or mapping
        except Exception: return mapping

    def _llm_match_labels(self, entities: list, group_names: list, intent: dict) -> dict:
        prompt = f"Entities: {entities}\nQuestion: {intent.get('summary', '')}\nLabels: {group_names}\nMatch entities to labels. Return JSON like: {{\"BMAL1\": \"B\", \"wildtype\": \"WT\"}}"
        result = self.think_structured(prompt, "Match biological names to data labels.", "Return JSON mapping entity to label.")
        return {k: v for k, v in result.items() if v in group_names and isinstance(v, str)} if isinstance(result, dict) else {}

    @staticmethod
    def _identify_control(group_names: list) -> str | None:
        for keyword in ["wt", "wildtype", "control", "ctrl", "vehicle", "dmso", "untreated", "scramble", "mock", "normal", "healthy", "baseline"]:
            for label in group_names:
                if label.lower() == keyword: return label
        for keyword in ["wt", "wildtype", "control", "ctrl", "vehicle", "dmso", "untreated", "scramble", "mock", "normal", "healthy", "baseline"]:
            for label in group_names:
                if keyword in label.lower(): return label
        return None

    @staticmethod
    def _humanize_contrast(num_label: str, den_label: str, entity_to_label: dict) -> str:
        label_to_entity = {v: k for k, v in entity_to_label.items()}
        return f"{label_to_entity.get(num_label, num_label)} vs {label_to_entity.get(den_label, den_label)}"

    def _publish_findings(self, experiment_id: str, result: dict):
        if not result.get("contrasts", []):
            self._publish_single_contrast_findings(experiment_id, result)
            return

        for c_result in result.get("contrasts", []):
            name, n_sig, n_up, n_down = c_result.get("name", "unknown"), c_result.get("n_significant", 0), c_result.get("n_upregulated", 0), c_result.get("n_downregulated", 0)
            conf = Confidence.HIGH if n_sig > 100 else Confidence.MEDIUM if n_sig > 10 else Confidence.LOW if n_sig > 0 else Confidence.INSUFFICIENT
            self.publish_finding(experiment_id, {"summary": f"[{name}] {n_sig} DE genes ({n_up} up, {n_down} down) at padj<{result.get('padj_threshold', 0.05)}, |log2FC|>{result.get('lfc_threshold',1.0)}", "contrast": name, "n_significant": n_sig, "n_upregulated": n_up, "n_downregulated": n_down, "top_genes": c_result.get("top_genes", [])[:10]}, conf)
            for db, terms in c_result.get("pathways", {}).items():
                if isinstance(terms, list) and terms:
                    self.publish_finding(experiment_id, {"summary": f"[{name}] {db}: {len(terms)} pathways. Top: {', '.join(t['term'] for t in terms[:3])}", "contrast": name, "pathways": terms[:10], "database": db}, Confidence.MEDIUM)

        if qc := result.get("sample_qc", {}):
            outliers = qc.get("outliers", [])
            self.publish_finding(experiment_id, {"summary": f"Bulk QC: {qc.get('n_samples','?')} samples. Lib size range: {qc.get('size_ratio',1):.1f}x." + (f" Outliers: {outliers}." if outliers else "")}, Confidence.HIGH if not outliers else Confidence.MEDIUM)

    def _publish_single_contrast_findings(self, experiment_id: str, result: dict):
        n_sig, n_up, n_down, comp = result.get("n_significant", 0), result.get("n_upregulated", 0), result.get("n_downregulated", 0), result.get("comparison_used", {})
        conf = Confidence.HIGH if n_sig > 100 else Confidence.MEDIUM if n_sig > 10 else Confidence.LOW if n_sig > 0 else Confidence.INSUFFICIENT
        self.publish_finding(experiment_id, {"summary": f"Bulk DE ({comp.get('numerator','?')} vs {comp.get('denominator','?')}): {n_sig} genes ({n_up} up, {n_down} down)", "n_significant": n_sig, "n_upregulated": n_up, "n_downregulated": n_down, "top_genes": result.get("top_genes", [])[:10]}, conf)

    def _publish_fastq_qc_findings(self, experiment_id: str, qc: dict):
        if samples := qc.get("samples", []):
            low_qual = [s["name"] for s in samples if s.get("pct_passed", 100) < 80]
            self.publish_finding(experiment_id, {"summary": f"FASTQ QC: {len(samples)} samples trimmed. Avg {sum(s.get('pct_passed', 0) for s in samples) / len(samples):.1f}% reads passed." + (f" Low quality: {low_qual}" if low_qual else ""), "multiqc": qc.get("multiqc_report")}, Confidence.MEDIUM if low_qual else Confidence.HIGH)

    def _publish_alignment_findings(self, experiment_id: str, align: dict):
        if bams := align.get("bam_files", []):
            ok_bams = [b for b in bams if b.get("status") == "success"]
            avg_map = sum(b.get("pct_unique", 0) for b in ok_bams) / max(len(ok_bams), 1)
            low_map = [b["name"] for b in ok_bams if b.get("pct_unique", 100) < 70]
            self.publish_finding(experiment_id, {"summary": f"STAR alignment: {len(ok_bams)}/{len(bams)} samples mapped. Avg unique mapping: {avg_map:.1f}%." + (f" Low mapping: {low_map}" if low_map else "")}, Confidence.MEDIUM if avg_map <= 75 or low_map else Confidence.HIGH)

    def _interpret(self, result: dict, intent: dict, exp_ctx: dict) -> str:
        if not (contrasts := result.get("contrasts", [])): return self._interpret_single(result, intent, exp_ctx)
        summaries = [f"  {c.get('name','?')}: {c.get('n_significant', 0)} DE genes ({c.get('n_upregulated', 0)} up, {c.get('n_downregulated', 0)} down). Top: {[g.get('symbol') or g.get('gene', g) for g in c.get('top_genes', [])[:5]]}. Pathways: {[t['term'] for db, terms in c.get('pathways', {}).items() if isinstance(terms, list) for t in terms[:2]][:3]}" for c in contrasts]
        pathway_warnings = [w for w in result.get("warnings", []) or [] if any(k in w.lower() for k in ("pathway", "enrichment", "enrichr", "gseapy", "go_bp", "kegg", "reactome", "symbol", "gtf"))]
        pathway_status_block = f"\nPATHWAY ENRICHMENT STATUS: NO RESULTS RETURNED.\nVerbatim warnings from the pipeline:\n{chr(10).join(f'  - {w}' for w in pathway_warnings) if pathway_warnings else '  - (no pathway-related warnings recorded)'}\nState the empty-result fact directly and quote ONE verbatim warning." if not any(c.get("pathways") for c in contrasts if c.get("status") == "success") else ""
        prompt = f"Bulk RNA-seq, multiple contrasts:\nOrganism: {exp_ctx.get('organism', '')}\nQuestion: {intent.get('summary', '')}\nThresholds: padj<{result.get('padj_threshold', 0.05)}, |log2FC|>{result.get('lfc_threshold', 1.0)}\n\nResults:\n{chr(10).join(summaries)}\n{pathway_status_block}\nWrite a 4-6 sentence biological synthesis comparing gene overlap and pathways."
        try: return self.llm.complete(prompt=prompt, system=BULK_RNA_SYSTEM, tier=TaskTier.HEAVY, max_tokens=500)
        except Exception: return f"{len(contrasts)} contrasts analyzed, {sum(c.get('n_significant', 0) for c in contrasts)} total DE genes."

    def _interpret_single(self, result: dict, intent: dict, exp_ctx: dict) -> str:
        prompt = f"Bulk RNA-seq DE: {result.get('comparison_used', {})}\nOrganism: {exp_ctx.get('organism', '')}\nQuestion: {intent.get('summary', '')}\n\nResults:\n  {result.get('n_significant', 0)} significant genes\n  Top genes: {[g.get('symbol') or g.get('gene', g) for g in result.get('top_genes', [])[:8]]}\nWrite a 3-4 sentence biological interpretation."
        try: return self.llm.complete(prompt=prompt, system=BULK_RNA_SYSTEM, tier=TaskTier.HEAVY, max_tokens=350)
        except Exception: return f"{result.get('n_significant', 0)} DE genes identified."

    def _output_dir(self, files: list) -> str: return str(Path(files[0]).parent / "aria_bulk_de") if files else "/tmp/aria_bulk_de"

    def receive(self, message): pass
