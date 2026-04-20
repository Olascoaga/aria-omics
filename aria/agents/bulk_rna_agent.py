"""
ARIA BulkRNAAgent (v3)
----------------------
Bulk RNA-seq differential expression and pathway analysis.

v3 CHANGES (critical biology fixes):
  1. Label-aware intent parsing — matches biological entities in the
     user's question (BMAL1, REV-ERBα) to actual sample label prefixes
     (B, R, WT). Uses heuristics + LLM fallback, not regex on free text.
  2. Multiple contrasts — when 3+ groups exist, runs all biologically
     meaningful pairwise comparisons (B vs WT, R vs WT) instead of one
     arbitrary contrast. Auto-identifies the control group.
  3. Context-aware LFC threshold — TF knockouts (BMAL1, etc.) use 0.58
     (1.5x) default because direct targets often have modest effect sizes.
     Non-TF perturbations keep the 1.0 (2x) default.
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


# ── Transcription factors — use lower LFC threshold for direct targets ────
KNOWN_TFS = {
    # Circadian
    "bmal1", "arntl", "clock", "per1", "per2", "per3",
    "cry1", "cry2", "nr1d1", "reverba", "nr1d2",
    "rora", "rorb", "rorc",
    # Pluripotency
    "oct4", "pou5f1", "sox2", "nanog", "klf4", "lin28",
    # Lineage TFs
    "gata1", "gata2", "gata3", "gata4", "gata6",
    "tbx5", "runx1", "runx2", "runx3",
    "foxa1", "foxa2", "foxp3", "foxo1", "foxo3",
    "cebpa", "cebpb", "pparg", "ppara",
    "rela", "relb", "stat1", "stat3", "stat5",
    "tp53", "trp53", "myc", "max", "sp1", "e2f1",
    "hif1a", "arnt", "yap1", "wwtr1",
}


def _is_fastq(files: list) -> bool:
    if not files:
        return False
    ext = str(files[0]).lower()
    return any(ext.endswith(s) for s in
               [".fastq.gz", ".fq.gz", ".fastq", ".fq"])


def _infer_lfc_threshold(intent: dict) -> float:
    """TF knockouts → 0.58 (1.5x). Others → 1.0 (2x)."""
    entities = [str(e).lower() for e in intent.get("biological_entities", [])]
    text     = " ".join([
        str(intent.get("summary",    "")),
        str(intent.get("comparison", "")),
        *entities,
    ]).lower()
    # Normalize hyphens (rev-erba → reverba)
    text_norm = re.sub(r'[\-\s]', '', text)

    if any(tf in text_norm for tf in KNOWN_TFS):
        return 0.58
    if re.search(r"\b(knockout|knockdown|ko|kd|overexpression|oe)\b", text):
        return 0.58
    if "transcription factor" in text or "regulator" in text:
        return 0.58
    return 1.0


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

        # ── Preprocessing if raw FASTQs ──────────────────────────────────
        if _is_fastq(files):
            counts_files, preprocessing = self._run_preprocessing(
                experiment_id, files, exp_ctx, intent
            )
            if counts_files is None:
                return {"status": "failed",
                        "findings": preprocessing,
                        "reason": "preprocessing_failed"}
            files = counts_files
        else:
            preprocessing = None

        # ── Discover actual group labels from counts matrix ──────────────
        sample_names, group_labels = self._discover_groups(files)

        if not group_labels:
            self.publish_finding(
                experiment_id,
                {"summary": "Could not infer experimental groups "
                            "from sample names."},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "group_inference_failed"}

        self.publish_status(
            experiment_id,
            f"Detected groups: {list(group_labels.keys())}",
            0.68,
        )

        # ── Build contrasts from intent + real labels ────────────────────
        design_factor, contrasts = self._build_contrasts(
            intent, group_labels, experiment_id
        )
        lfc_thr = _infer_lfc_threshold(intent)

        self.publish_status(
            experiment_id,
            f"Running {len(contrasts)} contrast(s), |log2FC|>{lfc_thr}...",
            0.70,
        )

        # ── Run all contrasts in one script invocation ───────────────────
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
                "padj_threshold": 0.05,
                "lfc_threshold":  lfc_thr,
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

        if preprocessing:
            result["preprocessing"] = preprocessing

        self._publish_findings(experiment_id, result)
        result["interpretation"] = self._interpret(result, intent, exp_ctx)

        # ── Auto-record methodology decisions to memory (NEW v3.9) ────
        # Each analytical choice becomes a row in the decisions log so
        # the narrative report can render them in the Decision Log table.
        # These are "auto" decisions (ARIA made them without user input).
        # In v4.0, DesignAgent will add interactive user decisions on top.
        self._record_methodology_decisions(experiment_id, result)

        self.publish_status(experiment_id, "BulkRNAAgent complete.", 1.0)
        return {"status": "done", "findings": result}

    def _record_methodology_decisions(self, experiment_id: str,
                                        result: dict) -> None:
        """
        Persist the methodology choices this agent made into the memory's
        decisions table. Each call is idempotent (INSERT OR REPLACE keyed
        on decision_id). Safe to call multiple times across reruns.
        """
        try:
            import uuid

            methodology = (result.get("methodology") or {})
            decisions   = methodology.get("decisions", []) or []
            if not decisions:
                return

            # Map methodology steps into the decision log schema.
            # checkpoint number is synthetic — reflects pipeline stage order.
            stage_to_cp = {
                "Differential expression (DESeq2)":        1,
                "PCA + MDS (sample-level structure)":      2,
                "Heatmap (padj top 50)":                   3,
                "Heatmap (|log2FC| top 50)":               4,
                "Pathway enrichment (ORA)":                5,
                "GSEA (pre-ranked)":                       6,
                "TPM (supplementary export)":              7,
            }

            for d in decisions:
                step = d.get("step", "")
                cp   = stage_to_cp.get(step, 0)
                # Build a compact decision summary: what tool / input / filter
                decision_summary = (
                    f"{d.get('input','?')} | "
                    f"{d.get('normalization','?')} | "
                    f"{d.get('gene_filter','?')}"
                )
                # Use deterministic ID so reruns overwrite (INSERT OR REPLACE)
                deterministic_id = f"{experiment_id[:8]}-auto-{cp:02d}"

                try:
                    self.memory.store_decision(
                        decision_id=deterministic_id,
                        wing_id=experiment_id,
                        checkpoint=cp,
                        question=step,
                        decision=decision_summary,
                        rationale=d.get("justification", "")[:500],
                        made_by="bulk_rna_agent (auto)",
                    )
                except Exception as e:
                    log.debug(f"Failed to store decision '{step}': {e}")

            # Also record two pipeline-level decisions (thresholds)
            try:
                self.memory.store_decision(
                    decision_id=f"{experiment_id[:8]}-auto-00-thr",
                    wing_id=experiment_id,
                    checkpoint=0,
                    question="Statistical thresholds for DE significance",
                    decision=(
                        f"padj < {result.get('padj_threshold', 0.05)}, "
                        f"|log2FC| > {result.get('lfc_threshold', 1.0)}"
                    ),
                    rationale=(
                        "padj < 0.05 is the community-standard FDR cutoff. "
                        "|log2FC| threshold is lower (0.58 ≈ 1.5-fold) for TF "
                        "knockouts because TFs mediate most effects indirectly "
                        "at modest magnitudes."
                    ),
                    made_by="bulk_rna_agent (auto)",
                )
            except Exception as e:
                log.debug(f"Failed to store threshold decision: {e}")

            # Design formula
            try:
                self.memory.store_decision(
                    decision_id=f"{experiment_id[:8]}-auto-00-design",
                    wing_id=experiment_id,
                    checkpoint=0,
                    question="DESeq2 design formula",
                    decision=result.get("design_used", "~condition"),
                    rationale=(
                        "Single-factor design inferred from sample labels. "
                        "No batch or covariate adjustment applied. "
                        "v4.0 DesignAgent will replace this with an "
                        "interactive design-confirmation checkpoint."
                    ),
                    made_by="bulk_rna_agent (auto)",
                )
            except Exception as e:
                log.debug(f"Failed to store design decision: {e}")

            # Report how many decisions were recorded (informational)
            log.info(
                f"Recorded {len(decisions) + 2} methodology decisions "
                f"to memory for experiment {experiment_id[:8]}."
            )
        except Exception as e:
            # Non-fatal: if decision logging fails, the pipeline still runs.
            log.warning(f"Decision logging failed (non-fatal): {e}")

    # ── FASTQ preprocessing pipeline ──────────────────────────────────────

    def _run_preprocessing(self, experiment_id: str, fastq_files: list,
                            exp_ctx: dict, intent: dict) -> tuple:
        fastq_dir  = str(Path(fastq_files[0]).parent)
        output_dir = str(Path(fastq_files[0]).parent.parent / "aria_processing")
        genome_cfg = exp_ctx.get("genome_config", {})

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
                {"summary": f"FASTQ QC failed: "
                            f"{qc_result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return None, qc_result

        self._publish_fastq_qc_findings(experiment_id, qc_result)
        self.publish_status(
            experiment_id,
            f"QC complete: {qc_result.get('n_samples',0)} samples trimmed",
            0.20,
        )

        self.publish_status(experiment_id,
                            "Aligning to genome (STAR)...", 0.25)
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
                {"summary": f"Alignment failed: "
                            f"{align_result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return None, align_result

        self._publish_alignment_findings(experiment_id, align_result)
        self.publish_status(
            experiment_id,
            f"Alignment complete: "
            f"{align_result.get('n_aligned', 0)} samples mapped",
            0.55,
        )

        self.publish_status(experiment_id,
                            "Counting reads (featureCounts)...", 0.60)
        quant_result = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_quantify.py",
            params={
                "bam_files":  align_result.get("bam_files", []),
                "gtf_file":   genome_cfg.get("gtf", ""),
                "output_dir": str(Path(output_dir) / "counts"),
                "threads":    8,
                "paired":     True,
                # "auto" triggers _detect_strandedness on the first BAM.
                # Override only if user/genome_config specifies an integer.
                "strand":     genome_cfg.get("strand", "auto"),
            },
        )
        if quant_result.get("status") == "error":
            self.publish_finding(
                experiment_id,
                {"summary": f"Quantification failed: "
                            f"{quant_result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT,
            )
            return None, quant_result

        counts_path = quant_result.get("counts_matrix")
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Quantification complete: "
                f"{quant_result.get('n_genes',0):,} genes × "
                f"{quant_result.get('n_samples',0)} samples"
            )},
            Confidence.HIGH,
        )
        self.publish_status(experiment_id, "Preprocessing complete.", 0.65)

        return [counts_path], {
            "qc":         qc_result,
            "alignment":  align_result,
            "quantification": quant_result,
        }

    # ── Group discovery ──────────────────────────────────────────────────

    def _discover_groups(self, files: list) -> tuple[list, dict]:
        """Peek at counts header, infer groups. Returns (samples, {label:[samples]})."""
        if not files:
            return [], {}

        sample_names = self._read_sample_names(files[0])
        if not sample_names:
            return [], {}

        groups = self._infer_groups_local(sample_names)
        if not groups:
            return sample_names, {}

        by_group: dict = {}
        for sample, label in groups.items():
            by_group.setdefault(label, []).append(sample)
        return sample_names, by_group

    @staticmethod
    def _read_sample_names(path: str) -> list[str]:
        try:
            p = Path(path)
            if not p.exists():
                return []
            opener = gzip.open if str(p).endswith(".gz") else open
            with opener(p, "rt") as f:
                header = f.readline().rstrip("\n")
            sep = "\t" if "\t" in header else ","
            cols = header.split(sep)
            skip = {
                "gene_id", "geneid", "", "chr", "start", "end",
                "strand", "length", "gene_name", "symbol",
                "feature", "ensembl", "ensembl_id", "entrez", "entrez_id",
            }
            return [c for c in cols if c.lower().strip() not in skip]
        except Exception as e:
            log.warning(f"Failed to read samples from {path}: {e}")
            return []

    @staticmethod
    def _infer_groups_local(samples: list[str]) -> dict:
        p1 = re.compile(r'^([A-Za-z][A-Za-z0-9]+)[_\-](\d+)$')
        m = {s: match.group(1) for s in samples if (match := p1.match(s))}
        if len(m) == len(samples) and len(set(m.values())) >= 2:
            return m

        p2 = re.compile(
            r'^([A-Za-z][A-Za-z0-9]+)[_\-]([Rr]ep\d+|[A-Za-z]\d*)$'
        )
        m = {s: match.group(1) for s in samples if (match := p2.match(s))}
        if len(m) == len(samples) and len(set(m.values())) >= 2:
            return m

        if all("_" in s for s in samples):
            gr = {s: "_".join(s.split("_")[:-1]) for s in samples}
            if len(set(gr.values())) >= 2:
                return gr

        p4 = re.compile(r'^([A-Za-z][A-Za-z0-9\-]*?)(\d.*)$')
        m = {s: match.group(1).rstrip("_-") for s in samples
             if (match := p4.match(s))}
        if len(m) == len(samples) and len(set(m.values())) >= 2:
            return m

        return {}

    # ── Intent ↔ label matching ───────────────────────────────────────────

    def _build_contrasts(self, intent: dict,
                          group_labels: dict,
                          experiment_id: str) -> tuple[str, list]:
        design_factor = self._infer_design_factor(intent)
        group_names   = list(group_labels.keys())
        entities      = intent.get("biological_entities", [])

        entity_to_label = self._map_entities_to_labels(
            entities, group_names, intent
        )

        control = self._identify_control(group_names)

        contrasts = []
        if control:
            for label in group_names:
                if label == control:
                    continue
                name = self._humanize_contrast(label, control, entity_to_label)
                contrasts.append({
                    "numerator":   label,
                    "denominator": control,
                    "name":        name,
                })
        else:
            sorted_g = sorted(group_names)
            for i, a in enumerate(sorted_g):
                for b in sorted_g[i+1:]:
                    contrasts.append({
                        "numerator":   a,
                        "denominator": b,
                        "name":        f"{a} vs {b}",
                    })

        if entity_to_label:
            mapping_str = ", ".join(
                f"{ent}→{lab}" for ent, lab in entity_to_label.items()
            )
            self.publish_finding(
                experiment_id,
                {"summary": f"Entity-to-label mapping: {mapping_str}. "
                            f"Contrasts: {[c['name'] for c in contrasts]}"},
                Confidence.HIGH,
            )

        return design_factor, contrasts

    @staticmethod
    def _infer_design_factor(intent: dict) -> str:
        text = (str(intent.get("comparison", "")).lower() + " " +
                str(intent.get("summary",    "")).lower())
        if any(k in text for k in ["knockout", "ko ", "knockdown", "kd ",
                                     "genotype", "mutant", "wt ", "wildtype"]):
            return "genotype"
        if any(k in text for k in ["treat", "drug", "vehicle", "dmso"]):
            return "treatment"
        if any(k in text for k in ["time", "timepoint", "hour", "day", "min"]):
            return "timepoint"
        return "condition"

    def _map_entities_to_labels(self, entities: list,
                                  group_names: list,
                                  intent: dict) -> dict:
        """Map biological entities → group labels. Heuristics first, LLM fallback."""
        mapping = {}
        used_labels = set()

        ent_clean = []
        for e in entities:
            s = re.sub(r'[^a-z0-9]', '', str(e).lower())
            if s and s not in ("cells", "cell", "h9", "h1", "hesc"):
                ent_clean.append((str(e), s))

        # Heuristic 1: label is prefix of (or equals) entity name
        #   e.g. label "B" matches entity "BMAL1" (B is prefix of bmal1)
        #   e.g. label "WT" matches entity "wildtype" (wt is prefix of wildtype)
        #   We iterate labels by descending length to prefer longer matches.
        sorted_labels = sorted(group_names, key=len, reverse=True)
        for label in sorted_labels:
            if label in used_labels:
                continue
            lbl_norm = label.lower()
            for original, norm in ent_clean:
                if original in mapping:
                    continue
                if norm.startswith(lbl_norm) and len(lbl_norm) >= 1:
                    mapping[original] = label
                    used_labels.add(label)
                    break

        # Heuristic 2: WT/wildtype/control entities → WT-like label
        wt_keywords = {"wt", "wildtype", "control", "ctrl", "untreated"}
        for original, norm in ent_clean:
            if original in mapping:
                continue
            if any(w in norm for w in wt_keywords):
                for label in group_names:
                    if label in used_labels:
                        continue
                    if label.lower() in wt_keywords:
                        mapping[original] = label
                        used_labels.add(label)
                        break

        # If we matched at least half the entities, trust heuristics
        if len(mapping) >= max(1, len(ent_clean) // 2):
            return mapping

        # LLM fallback
        try:
            llm_mapping = self._llm_match_labels(
                entities, group_names, intent
            )
            if llm_mapping:
                return llm_mapping
        except Exception as e:
            log.warning(f"LLM label matching failed: {e}")

        return mapping

    def _llm_match_labels(self, entities: list,
                           group_names: list,
                           intent: dict) -> dict:
        prompt = f"""
User asked about these biological entities: {entities}
Biological question: {intent.get('summary', '')}
Comparison described: {intent.get('comparison', '')}
Actual sample group labels found in the data: {group_names}

Match each biological entity to its corresponding group label.
Return JSON like: {{"BMAL1": "B", "REV-ERBa": "R", "wildtype": "WT"}}
Return only entities that have a clear match. If uncertain, omit.
"""
        result = self.think_structured(
            prompt=prompt,
            system="You match biological names to data labels. Be conservative.",
            schema_hint="Return a JSON object mapping entity names to label strings.",
        )
        if isinstance(result, dict):
            return {k: v for k, v in result.items()
                    if v in group_names and isinstance(v, str)}
        return {}

    @staticmethod
    def _identify_control(group_names: list) -> str | None:
        priority = [
            "wt", "wildtype", "control", "ctrl",
            "vehicle", "dmso", "untreated", "scramble",
            "mock", "normal", "healthy", "baseline",
        ]
        for keyword in priority:
            for label in group_names:
                if label.lower() == keyword:
                    return label
        for keyword in priority:
            for label in group_names:
                if keyword in label.lower():
                    return label
        return None

    @staticmethod
    def _humanize_contrast(num_label: str, den_label: str,
                            entity_to_label: dict) -> str:
        label_to_entity = {v: k for k, v in entity_to_label.items()}
        num_human = label_to_entity.get(num_label, num_label)
        den_human = label_to_entity.get(den_label, den_label)
        return f"{num_human} vs {den_human}"

    # ── Findings publishing ──────────────────────────────────────────────

    def _publish_findings(self, experiment_id: str, result: dict):
        contrast_results = result.get("contrasts", [])

        if not contrast_results:
            self._publish_single_contrast_findings(experiment_id, result)
            return

        for c_result in contrast_results:
            name   = c_result.get("name", "unknown")
            n_sig  = c_result.get("n_significant",   0)
            n_up   = c_result.get("n_upregulated",   0)
            n_down = c_result.get("n_downregulated", 0)

            conf = (Confidence.HIGH   if n_sig > 100 else
                    Confidence.MEDIUM if n_sig > 10  else
                    Confidence.LOW    if n_sig > 0   else
                    Confidence.INSUFFICIENT)

            self.publish_finding(
                experiment_id,
                {"summary": (
                    f"[{name}] {n_sig} DE genes "
                    f"({n_up} up, {n_down} down) "
                    f"at padj<0.05, |log2FC|>{result.get('lfc_threshold',1.0)}"
                ),
                 "contrast":        name,
                 "n_significant":   n_sig,
                 "n_upregulated":   n_up,
                 "n_downregulated": n_down,
                 "top_genes":       c_result.get("top_genes", [])[:10]},
                conf,
            )

            for db, terms in c_result.get("pathways", {}).items():
                if isinstance(terms, list) and terms:
                    self.publish_finding(
                        experiment_id,
                        {"summary": (
                            f"[{name}] {db}: {len(terms)} pathways. "
                            f"Top: {', '.join(t['term'] for t in terms[:3])}"
                        ),
                         "contrast": name,
                         "pathways": terms[:10],
                         "database": db},
                        Confidence.MEDIUM,
                    )

        qc = result.get("sample_qc", {})
        if qc:
            outliers = qc.get("outliers", [])
            self.publish_finding(
                experiment_id,
                {"summary": (
                    f"Bulk QC: {qc.get('n_samples','?')} samples. "
                    f"Lib size range: {qc.get('size_ratio',1):.1f}x."
                    + (f" Outliers: {outliers}." if outliers else "")
                )},
                Confidence.HIGH if not outliers else Confidence.MEDIUM,
            )

    def _publish_single_contrast_findings(self, experiment_id: str,
                                             result: dict):
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
             "top_genes":       result.get("top_genes", [])[:10]},
            conf,
        )

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
        ok_bams = [b for b in bams if b.get("status") == "success"]
        avg_map = sum(b.get("pct_unique", 0) for b in ok_bams) / \
                  max(len(ok_bams), 1)
        low_map = [b["name"] for b in ok_bams
                   if b.get("pct_unique", 100) < 70]
        conf = Confidence.HIGH if avg_map > 75 and not low_map \
               else Confidence.MEDIUM
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"STAR alignment: {len(ok_bams)}/{len(bams)} samples mapped. "
                f"Avg unique mapping: {avg_map:.1f}%."
                + (f" Low mapping: {low_map}" if low_map else "")
            )},
            conf,
        )

    # ── LLM interpretation ────────────────────────────────────────────────

    def _interpret(self, result: dict, intent: dict,
                    exp_ctx: dict) -> str:
        contrasts = result.get("contrasts", [])
        if not contrasts:
            return self._interpret_single(result, intent, exp_ctx)

        summaries = []
        for c in contrasts:
            tops = [g.get("symbol") or g.get("gene", g)
                    for g in c.get("top_genes", [])[:5]]
            top_pw = [
                t["term"]
                for db, terms in c.get("pathways", {}).items()
                if isinstance(terms, list)
                for t in terms[:2]
            ][:3]
            summaries.append(
                f"  {c.get('name','?')}: "
                f"{c.get('n_significant', 0)} DE genes "
                f"({c.get('n_upregulated', 0)} up, "
                f"{c.get('n_downregulated', 0)} down). "
                f"Top: {tops}. Pathways: {top_pw}"
            )

        # Collect ALL pathway-related warnings to feed verbatim to LLM.
        # This is the anti-hallucination guard: if pathways are empty,
        # the LLM should report the EXACT error from the script, not
        # invent a plausible-sounding cause.
        all_warnings = result.get("warnings", []) or []
        pathway_warnings = [
            w for w in all_warnings
            if any(k in w.lower() for k in
                   ("pathway", "enrichment", "enrichr", "gseapy",
                    "go_bp", "kegg", "reactome", "symbol", "gtf"))
        ]

        # Detect the empty-pathways case explicitly
        has_pathways = any(c.get("pathways") for c in contrasts
                            if c.get("status") == "success")
        pathway_status_block = ""
        if not has_pathways:
            pathway_status_block = f"""

PATHWAY ENRICHMENT STATUS: NO RESULTS RETURNED.
Verbatim warnings from the pipeline (use these EXACT facts —
do NOT invent or speculate about other causes):
{chr(10).join(f"  - {w}" for w in pathway_warnings) if pathway_warnings
              else "  - (no pathway-related warnings recorded)"}

When discussing this in your synthesis:
- State the empty-result fact directly
- Quote ONE of the verbatim warnings above as the cause if any exists
- Do NOT invent organism naming issues, annotation problems, or other
  causes that aren't in the warnings above
- If no warnings exist, just say "pathway enrichment did not return
  results; cause unclear from the available logs"
"""

        prompt = f"""
Bulk RNA-seq, multiple contrasts:
Organism: {exp_ctx.get("organism", "")}
Question: {intent.get("summary", "")}
LFC threshold used: {result.get('lfc_threshold', 1.0)}

Per-contrast results:
{chr(10).join(summaries)}
{pathway_status_block}
Write a 4-6 sentence biological synthesis:
- Compare gene overlap and pathway convergence across contrasts
- Identify shared vs unique biology between the perturbations
- Ground interpretation in the specific biology (no generic statements)
- If pathway enrichment is empty, state this honestly using the
  verbatim warning above — do NOT speculate about causes
"""
        try:
            return self.llm.complete(
                prompt=prompt,
                system=BULK_RNA_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=500,
            )
        except Exception:
            total_sig = sum(c.get("n_significant", 0) for c in contrasts)
            return (f"{len(contrasts)} contrasts analyzed, "
                    f"{total_sig} total DE genes across them.")

    def _interpret_single(self, result: dict, intent: dict,
                            exp_ctx: dict) -> str:
        n_sig  = result.get("n_significant", 0)
        tops   = [g.get("symbol") or g.get("gene", g)
                  for g in result.get("top_genes", [])[:8]]
        top_pw = [t["term"]
                  for db, terms in result.get("pathways", {}).items()
                  if isinstance(terms, list)
                  for t in terms[:2]]

        prompt = f"""
Bulk RNA-seq DE: {result.get("comparison_used", {})}
Organism: {exp_ctx.get("organism", "")}
Question: {intent.get("summary", "")}

Results:
  {n_sig} significant genes
  {result.get("n_upregulated", 0)} up / {result.get("n_downregulated", 0)} down
  Top genes: {tops}
  Top pathways: {top_pw[:5]}

Write a 3-4 sentence biological interpretation.
"""
        try:
            return self.llm.complete(
                prompt=prompt,
                system=BULK_RNA_SYSTEM,
                tier=TaskTier.HEAVY,
                max_tokens=350,
            )
        except Exception:
            return f"{n_sig} DE genes identified."

    def _output_dir(self, files: list) -> str:
        if files:
            return str(Path(files[0]).parent / "aria_bulk_de")
        return "/tmp/aria_bulk_de"

    def receive(self, message):
        pass
