"""
ARIA BulkRNAAgent (v3.10 → v4.0)
---------------------------------
Bulk RNA-seq differential expression and pathway analysis.

v4.0 CHANGE: if exp_context contains a 'design' key (from DesignAgent),
agent attempts to use it. If application fails, falls back to automatic
inference from file names (backward compatible).
"""

from __future__ import annotations

import csv
import gzip
import logging
import os
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
    if any(tf in text for tf in KNOWN_TFS) or re.search(r"\b(knockout|knockdown|ko|kd|overexpression|oe)\b", text) or "transcriptionfactor" in text:
        return 0.58
    return 1.0


def _normalise_sample_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


class BulkRNAAgent(BaseAgent):

    name        = "bulk_rna_agent"
    description = "Bulk RNA-seq: QC, DESeq2, pathway enrichment, plots."

    def __init__(self, memory: ARIAMemory, llm: LLMProvider, api_key: str = None):
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

        if _is_fastq(files):
            counts_files, preprocessing = self._run_preprocessing(experiment_id, files, exp_ctx, intent)
            if counts_files is None:
                return {"status": "failed", "findings": preprocessing, "reason": "preprocessing_failed"}
            files = counts_files
        else:
            preprocessing = None

        # ── v4.0: Attempt to apply DesignAgent output; fallback to legacy if fails ──
        design = exp_ctx.get("design")
        design_application_failed = False
        if design and design.get("groups"):
            self.publish_status(experiment_id, "Applying confirmed experimental design...", 0.65)
            try:
                sample_names, group_labels, design_factor, contrasts = \
                    self._apply_design(design, files, experiment_id)
                # Record the design fact
                self.memory.store_decision(
                    decision_id=f"{experiment_id[:8]}-design-01",
                    wing_id=experiment_id,
                    checkpoint=2,
                    question="Experimental design (source: DesignAgent)",
                    decision=design.get("design_formula", ""),
                    rationale="User confirmed design interactively.",
                    made_by="user"
                )
            except Exception as e:
                log.warning(f"Design application failed ({e})")
                design_application_failed = True
                design = None
        if not design or not design.get("groups"):
            # P0-6: filename/column-name group inference is a GUESS, not a
            # confirmed design. In production it is gated: a confirmed design
            # that failed to apply must STOP (never silently degrade to
            # guessing), and a run with no confirmed design at all also stops.
            # Opt in explicitly with ARIA_ALLOW_FILENAME_FALLBACK=1.
            if not self._filename_fallback_allowed():
                reason = ("design_application_failed" if design_application_failed
                          else "no_confirmed_design")
                detail = (
                    "the confirmed experimental design could not be applied to "
                    "the count matrix" if design_application_failed
                    else "no confirmed experimental design was provided"
                )
                self.publish_finding(experiment_id,
                    {"summary": (
                        f"Bulk DE stopped: {detail}. ARIA refuses to infer the "
                        "design from file/column names in production (it would "
                        "run DE on a guessed design). Provide a confirmed design, "
                        "or set ARIA_ALLOW_FILENAME_FALLBACK=1 to allow "
                        "name-based inference (NOT publication-grade)."
                    ), "reason": reason},
                    Confidence.INSUFFICIENT)
                return {"status": "failed", "reason": reason}
            # ── Legacy inference (explicit opt-in only) ────────────────
            self.publish_finding(experiment_id,
                {"summary": (
                    "ARIA_ALLOW_FILENAME_FALLBACK=1: inferring experimental "
                    "groups from file/column names. This is a GUESS, not a "
                    "confirmed design — results are not publication-grade."
                )},
                Confidence.LOW)
            sample_names, group_labels = self._discover_groups(files)
            if not group_labels:
                self.publish_finding(experiment_id,
                    {"summary": "Could not infer experimental groups from sample names."},
                    Confidence.INSUFFICIENT)
                return {"status": "failed", "reason": "group_inference_failed"}
            design_factor = self._infer_design_factor(intent)
            contrasts = []
        # ────────────────────────────────────────────────────────────────────

        self.publish_status(
            experiment_id,
            f"Detected groups: {sorted(set(group_labels.values()))}",
            0.68,
        )
        if not contrasts:
            return self._block_for_explicit_contrast(
                experiment_id, design_factor, group_labels
            )

        padj_thr = exp_ctx.get("global_padj", 0.05)
        lfc_thr  = exp_ctx.get("global_lfc", _infer_lfc_threshold(intent))
        output_dir = self._output_dir(files)
        metadata_file = (
            str(self._write_design_metadata(group_labels, design_factor, output_dir))
            if design and group_labels else ""
        )

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
                "metadata_file":   metadata_file,
                "contrasts":      contrasts,
                # P0-4: forward confirmed covariates (e.g. batch) so DESeq2 fits
                # `~ batch + condition`, not a bare `~ condition`.
                "covariates":     self._design_covariates(design) if design else [],
                "organism":       exp_ctx.get("organism", "Homo sapiens"),
                "genome":         exp_ctx.get("genome", "hg38"),
                "output_dir":     output_dir,
                "run_pathways":   True,
                "padj_threshold": padj_thr,
                "lfc_threshold":  lfc_thr,
            },
        )

        if result.get("status") == "error":
            self.publish_finding(experiment_id,
                {"summary": f"Bulk DE failed: {result.get('details','')[:100]}"},
                Confidence.INSUFFICIENT)
            return {"status": "failed", "findings": result}

        if preprocessing:
            result["preprocessing"] = preprocessing

        result["padj_threshold"] = padj_thr
        result["lfc_threshold"]  = lfc_thr

        self._publish_findings(experiment_id, result)
        result["interpretation"] = self._interpret(result, intent, exp_ctx)
        self._record_methodology_decisions(experiment_id, result)

        self.publish_status(experiment_id, "BulkRNAAgent complete.", 1.0)
        return {"status": "done", "findings": result}

    #  v4.0 Design integration ─────────────────────────────────────────
    def _apply_design(self, design: dict, files: list, experiment_id: str):
        """
        Use the DesignAgent-confirmed groups and factor to build
        sample mapping, contrasts, and design factor.
        Raises ValueError if mapping fails.
        """
        groups = design["groups"]
        factor = design.get("main_factor", "condition")
        sample_stems = {sample: group for group, samples in groups.items() for sample in samples}
        sample_aliases = design.get("sample_aliases", {}) or {}
        alias_to_group = {}
        for stem, grp in sample_stems.items():
            aliases = sample_aliases.get(stem, [stem])
            for alias in aliases:
                if alias:
                    alias_to_group[str(alias)] = grp

        # Read actual column names from the first count file (gzip-aware)
        try:
            p = Path(files[0])
            opener = gzip.open if str(p).endswith(".gz") else open
            with opener(p, "rt") as f:
                header = f.readline().rstrip("\n")
            sep      = "\t" if "\t" in header else ","
            colnames = [c.strip().strip('"').strip("'") for c in header.split(sep)]
        except Exception as e:
            raise ValueError(f"Cannot read header from counts file {files[0]}: {e}")

        _GENE_COLS = {"geneid", "gene_id", "ensembl_id", "", "gene",
                      "gene_name", "symbol", "feature", "name"}

        # Map stem→group to column→group using fuzzy match
        group_labels = {}
        for col in colnames:
            if col.lower() in _GENE_COLS:
                continue
            best_match = None
            col_norm = _normalise_sample_token(col)
            for stem, grp in {**sample_stems, **alias_to_group}.items():
                stem_norm = _normalise_sample_token(stem)
                if stem in col or stem_norm in col_norm:
                    if best_match is not None:
                        log.warning(
                            f"Ambiguous column '{col}' matches both group "
                            f"'{best_match}' and '{stem}' — skipping, will "
                            f"retry with suffix-trimming or positional fallback."
                        )
                        best_match = None
                        break
                    best_match = grp
            if best_match:
                group_labels[col] = best_match
                

                # Attempt to match by trimming potential technical replicate suffixes (_1, _2, etc.)
        ordered_cols = [c for c in colnames if c.lower() not in _GENE_COLS]
        if not group_labels:
            for col in ordered_cols:
                best_match = None
                col_norm = _normalise_sample_token(col)
                for stem, grp in {**sample_stems, **alias_to_group}.items():
                    # Trim common suffixes from col
                    col_base = re.sub(r"[_\-]?[12]$", "", col)
                    col_base_norm = _normalise_sample_token(col_base)
                    stem_norm = _normalise_sample_token(stem)
                    if (
                        col_base == stem
                        or stem.startswith(col_base)
                        or col_base.startswith(stem)
                        or col_base_norm == stem_norm
                        or stem_norm.startswith(col_base_norm)
                        or col_base_norm.startswith(stem_norm)
                        or stem_norm in col_norm
                    ):
                        if best_match is not None:
                            best_match = None
                            break
                        best_match = grp
                if best_match:
                    group_labels[col] = best_match

        # Fallback: if no match, assign from design directly (assumes column order)
        if not group_labels:
            ordered_samples = list(sample_stems.keys())
            if len(ordered_samples) == len(ordered_cols):
                for i, col in enumerate(ordered_cols):
                    group_labels[col] = sample_stems[ordered_samples[i]]

        # Step 4: column-name inference when GEO sample count mismatches matrix columns
        # (e.g. spike-in experiments where only a subset of GEO samples is in the matrix)
        if not group_labels or len(set(group_labels.values())) < 2:
            col_groups = self._infer_col_groups(ordered_cols)
            if len(set(col_groups.values())) >= 2:
                log.warning(
                    f"GEO design has {len(sample_stems)} samples but matrix has "
                    f"{len(ordered_cols)} columns — using column-name group inference."
                )
                group_labels = col_groups
                groups = {}
                for col, grp in col_groups.items():
                    groups.setdefault(grp, []).append(col)

        if not group_labels or len(set(group_labels.values())) < 2:
            raise ValueError(
                f"Could not map design groups to count matrix columns. "
                f"Design stems: {list(sample_stems.keys())}, Columns: {colnames}"
            )

        sample_names = list(group_labels.keys())
        contrasts = self._normalise_explicit_contrasts(
            design.get("plan_contrasts")
            or design.get("contrasts")
            or design.get("comparisons"),
            groups.keys(),
        )
        return sample_names, group_labels, factor, contrasts

    # ── Legacy inference (unchanged) ────────────────────────────────────
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
    def _filename_fallback_allowed() -> bool:
        """P0-6: inferring the experimental design from file/column names is a
        guess and is OFF in production. It is allowed only when
        ARIA_ALLOW_FILENAME_FALLBACK is explicitly set to a truthy value."""
        return os.environ.get(
            "ARIA_ALLOW_FILENAME_FALLBACK", ""
        ).strip().lower() in {"1", "true", "yes", "on"}

    @staticmethod
    def _design_covariates(design: dict) -> list[str]:
        """Covariates confirmed at DesignAgent CHECKPOINT 2.4, forwarded to the
        DESeq2 design (P0-4). The confirmed batch covariate plus any explicit
        covariates list, deduped and order-preserving."""
        covariates: list[str] = []
        batch = (design or {}).get("batch_covariate")
        if batch:
            covariates.append(batch)
        for cov in (design or {}).get("covariates", []) or []:
            if cov:
                covariates.append(cov)
        return list(dict.fromkeys(covariates))

    @staticmethod
    def _normalise_explicit_contrasts(raw, available_groups) -> list[dict]:
        """Return only caller-confirmed contrasts with valid test/ref levels."""
        if not raw:
            return []
        available = {str(g) for g in available_groups}
        by_token: dict[str, str | None] = {}
        for group in available:
            token = _normalise_sample_token(group)
            by_token[token] = group if token not in by_token else None

        def resolve_level(value) -> str:
            level = str(value or "").strip()
            if level in available:
                return level
            token = _normalise_sample_token(level)
            resolved = by_token.get(token)
            return resolved if resolved is not None else ""

        contrasts: list[dict] = []
        for comp in raw:
            if isinstance(comp, dict):
                num = comp.get("numerator") or comp.get("test") or comp.get("case")
                den = (
                    comp.get("denominator")
                    or comp.get("reference")
                    or comp.get("ref")
                    or comp.get("control")
                )
                name = comp.get("name")
            elif isinstance(comp, (list, tuple)) and len(comp) >= 2:
                num, den = comp[0], comp[1]
                name = None
            else:
                continue
            num = resolve_level(num)
            den = resolve_level(den)
            if not num or not den:
                continue
            contrasts.append({
                "numerator": num,
                "denominator": den,
                "name": str(name or f"{num} vs {den}"),
            })
        return contrasts

    @staticmethod
    def _write_design_metadata(
        group_labels: dict[str, str], design_factor: str, output_dir: str | Path
    ) -> Path:
        """Persist the confirmed sample-to-group mapping for rna_bulk_de.py."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "confirmed_design_metadata.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", design_factor])
            for sample, group in group_labels.items():
                writer.writerow([sample, group])
        return path

    @staticmethod
    def _suggest_contrasts_from_groups(group_labels: dict) -> list[dict]:
        """Suggest candidate contrasts; suggestions do not authorize DE."""
        if not group_labels:
            return []
        values = list(group_labels.values())
        if values and all(isinstance(v, (list, tuple, set)) for v in values):
            groups = sorted(str(g) for g in group_labels.keys())
        else:
            groups = sorted({str(g) for g in values})
        if len(groups) < 2:
            return []
        control = BulkRNAAgent._identify_control(groups)
        contrasts = []
        if control:
            for grp in groups:
                if grp != control:
                    contrasts.append({
                        "numerator": grp,
                        "denominator": control,
                        "name": f"{grp} vs {control}",
                    })
        else:
            for i, ref in enumerate(groups):
                for test in groups[i + 1:]:
                    contrasts.append({
                        "numerator": test,
                        "denominator": ref,
                        "name": f"{test} vs {ref}",
                    })
        return contrasts

    def _block_for_explicit_contrast(
        self, experiment_id: str, design_factor: str, group_labels: dict
    ) -> dict:
        suggestions = self._suggest_contrasts_from_groups(group_labels)
        self.publish_finding(
            experiment_id,
            {
                "summary": (
                    "Bulk DE was not run because no explicit numerator/reference "
                    "contrast was confirmed."
                ),
                "design_factor": design_factor,
                "suggested_contrasts": suggestions,
            },
            Confidence.INSUFFICIENT,
        )
        options = ["Skip bulk DE until a contrast is confirmed"]
        options.extend(
            f"Run {c['numerator']} vs {c['denominator']} "
            f"(reference={c['denominator']})"
            for c in suggestions[:6]
        )
        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint="bulk.contrast",
            question=(
                "Bulk differential expression requires an explicit contrast. "
                "Choose the numerator/test level and denominator/reference "
                "level before ARIA runs DE."
            ),
            options=options,
            context={
                "analysis_type": "bulk_differential_expression",
                "parameter_name": "contrast",
                "suggested_contrasts": suggestions,
                "design_factor": design_factor,
            },
        )
        return {
            "status": "failed",
            "reason": "explicit_contrast_required",
            "findings": {
                "status": "error",
                "error_type": "ExplicitContrastRequired",
                "details": (
                    "Confirm a contrast with explicit numerator/test and "
                    "denominator/reference levels before running DE."
                ),
                "design_factor": design_factor,
                "suggested_contrasts": suggestions,
            },
        }

    @staticmethod
    def _read_sample_names(path: str) -> list[str]:
        try:
            p = Path(path)
            if not p.exists(): return []
            with (gzip.open if str(p).endswith(".gz") else open)(p, "rt") as f: header = f.readline().rstrip("\n")
            return [c for c in header.split("\t" if "\t" in header else ",") if c.lower().strip() not in {"gene_id", "geneid", "gene", "", "chr", "start", "end", "strand", "length", "gene_name", "symbol", "feature", "ensembl", "ensembl_id", "entrez", "entrez_id"}]
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

    @staticmethod
    def _infer_col_groups(cols: list[str]) -> dict:
        """Infer group labels from column names using the same patterns as rna_bulk_de.py."""
        for pattern in [
            r'^([A-Za-z][A-Za-z0-9]+)[_\-](\d+)$',
            r'^([A-Za-z][A-Za-z0-9]+)[_\-]([Rr]ep\d+|[A-Za-z]\d*)$',
            r'^([A-Za-z][A-Za-z0-9\-]*?)(\d.*)$',
        ]:
            m = {c: match.group(1).rstrip("_-") for c in cols if (match := re.compile(pattern).match(c))}
            if len(m) == len(cols) and len(set(m.values())) >= 2:
                return m
        if all("_" in c for c in cols):
            gr = {c: "_".join(c.split("_")[:-1]) for c in cols}
            if len(set(gr.values())) >= 2:
                return gr
        return {}

    def _build_contrasts(self, intent: dict, group_labels: dict, experiment_id: str) -> tuple[str, list]:
        design_factor = self._infer_design_factor(intent)
        # P0-5: this legacy helper no longer authorizes executable contrasts.
        # The run path emits a confirmation checkpoint via
        # _block_for_explicit_contrast when no plan/design contrast is present.
        return design_factor, []

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
        prompt = (
            f"Entities: {entities}\n"
            f"Question: {intent.get('summary', '')}\n"
            f"Labels: {group_names}\n"
            "Match biological names to data labels. Return only JSON, "
            "for example: {\"entity_name\": \"label_name\"}."
        )
        result = self.think_structured(prompt, "Match biological names to data labels.", "Return JSON mapping entity to label.")
        return {k: v for k, v in result.items() if v in group_names and isinstance(v, str)} if isinstance(result, dict) else {}

    @staticmethod
    def _identify_control(group_names: list) -> str | None:
        for keyword in ["wt", "wildtype", "control", "ctrl", "ctr", "vehicle", "dmso", "untreated", "scramble", "mock", "normal", "healthy", "baseline"]:
            for label in group_names:
                if label.lower() == keyword: return label
        for keyword in ["wt", "wildtype", "control", "ctrl", "ctr", "vehicle", "dmso", "untreated", "scramble", "mock", "normal", "healthy", "baseline"]:
            for label in group_names:
                if keyword in label.lower(): return label
        return None

    @staticmethod
    def _humanize_contrast(num_label: str, den_label: str, entity_to_label: dict) -> str:
        label_to_entity = {v: k for k, v in entity_to_label.items()}
        return f"{label_to_entity.get(num_label, num_label)} vs {label_to_entity.get(den_label, den_label)}"

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
