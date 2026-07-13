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
from aria.llm.provider import LLMProvider
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.bulk_rna")


def _is_fastq(files: list) -> bool:
    if not files: return False
    return any(str(files[0]).lower().endswith(s) for s in [".fastq.gz", ".fq.gz", ".fastq", ".fq"])

# F1 (preprint audit 2026-06-19 / ADR-055): the |log2FC| DE-significance threshold
# is DATA- AND PROMPT-INDEPENDENT by default. It must never be derived from the
# question text or a hardcoded gene list — otherwise the inferential cutoff (and
# thus which genes are called DE) depends on how the question is worded, which is
# not reproducible and violates "LLM proposes / code guarantees". The cutoff
# deviates from the default ONLY via an explicit, user-confirmed CP3 threshold
# profile or a versioned `global_lfc` override.
DEFAULT_LFC_THRESHOLD = 1.0


def _default_lfc_threshold() -> float:
    """The fixed, reproducible default |log2FC| significance threshold."""
    return DEFAULT_LFC_THRESHOLD


def suggest_lfc_profile(intent: dict) -> str | None:
    """ADVISORY ONLY — never applied to the threshold.

    Returns a non-binding hint (``"exploratory_tf"`` or ``None``) that a
    perturbation study (knockout / knockdown / overexpression) often warrants the
    user-selectable Exploratory/TF CP3 profile (LFC 0.58), so the checkpoint can
    surface it. ARIA guides from the biological question, but the actual cutoff is
    only ever changed by an explicit user choice — code does NOT silently move the
    statistical threshold based on prompt text.

    The hint is derived ONLY from generic experimental-design keywords, never from
    a hardcoded gene list (ADR-055): the advisory must not depend on which
    specific genes ARIA happens to know about, otherwise it silently privileges a
    curated, incomplete set of genes and bakes biology into the code.
    """
    text = " ".join([
        str(intent.get("summary", "")),
        str(intent.get("comparison", "")),
        *[str(e) for e in intent.get("biological_entities", [])],
    ]).lower()
    perturbation = re.search(
        r"\b(knock[\s-]?outs?|knock[\s-]?downs?|ko|kd|over[\s-]?expression|oe|"
        r"transcription[\s-]?factors?)\b",
        text,
    )
    return "exploratory_tf" if perturbation else None


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
        replicate_units: dict[str, str] = {}
        if design and design.get("groups"):
            self.publish_status(experiment_id, "Applying confirmed experimental design...", 0.65)
            try:
                sample_names, group_labels, design_factor, contrasts = \
                    self._apply_design(design, files, experiment_id)
                replicate_units = self._resolve_replicate_units(
                    design, sample_names
                )
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
        lfc_thr  = exp_ctx.get("global_lfc", _default_lfc_threshold())
        # F1: record how the cutoff was chosen so the report can state it. An
        # explicit CP3 profile stamps its own provenance; a bare global_lfc
        # override is "explicit_global_lfc"; otherwise it is the fixed default.
        lfc_provenance = (
            exp_ctx.get("lfc_threshold_provenance")
            or ("explicit_global_lfc" if "global_lfc" in exp_ctx else "default"))
        output_dir = self._output_dir(files)
        metadata_file = (
            str(self._write_design_metadata(
                group_labels,
                design_factor,
                output_dir,
                replicate_units=replicate_units,
            ))
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
                "technical_replicate_col": (
                    "biological_unit" if replicate_units else ""
                ),
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
        result["lfc_threshold_provenance"] = lfc_provenance

        self._publish_findings(experiment_id, result)
        result["interpretation_status"] = {
            "ran": False,
            "reason": (
                "free_text_llm_interpretation_disabled_by_F3_governance; "
                "bulk RNA Results are generated from structured DE/pathway "
                "outputs and NarrativeBlock evidence cards"
            ),
            "governance": "F3_preprint_audit",
        }
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

        # B4 preprint-readiness: NO column-name inference rescue here. When the
        # confirmed design cannot be mapped to the matrix columns (e.g. a
        # sample-count mismatch from incomplete metadata), fail closed so the
        # caller's P0-6 gate handles it — blocking, or inferring only under the
        # explicit ARIA_ALLOW_FILENAME_FALLBACK opt-in via the sample-name path.
        # Rebuilding the design from ctrl_*/treat_* column names would silently
        # discard the confirmed design and run DE on a guess (P0-5/F5).
        if not group_labels or len(set(group_labels.values())) < 2:
            raise ValueError(
                f"Could not map the confirmed design groups to the count matrix "
                f"columns. Design stems: {list(sample_stems.keys())}, "
                f"Columns: {colnames}. ARIA refuses to infer the design from "
                f"column names; provide metadata that maps the confirmed design "
                f"to these columns."
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
        covariates list, deduped and order-preserving. E7: this delegates to the
        single shared RNA/ATAC design contract so RNA and both ATAC lanes extract
        covariates identically."""
        from aria.utils.design_matrix import resolve_design_contract

        return resolve_design_contract(design, {})["covariates"]

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
        group_labels: dict[str, str],
        design_factor: str,
        output_dir: str | Path,
        *,
        replicate_units: dict[str, str] | None = None,
    ) -> Path:
        """Persist the confirmed sample-to-group mapping for rna_bulk_de.py."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "confirmed_design_metadata.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            header = ["sample", design_factor]
            if replicate_units:
                header.append("biological_unit")
            writer.writerow(header)
            for sample, group in group_labels.items():
                row = [sample, group]
                if replicate_units:
                    row.append(replicate_units[sample])
                writer.writerow(row)
        return path

    @staticmethod
    def _resolve_replicate_units(
        design: dict, sample_names: list[str]
    ) -> dict[str, str]:
        """Map actual matrix columns to CP2.5 biological units, fail-closed."""
        handling = (design or {}).get("replicate_handling", {}) or {}
        if handling.get("mode") != "technical_aggregate":
            return {}
        declared = handling.get("sample_to_unit", {}) or {}
        if not declared:
            raise ValueError("technical replicate design has no sample-to-unit map")

        resolved: dict[str, str] = {}
        for sample in sample_names:
            token = _normalise_sample_token(sample)
            exact = [
                unit for declared_sample, unit in declared.items()
                if _normalise_sample_token(declared_sample) == token
            ]
            candidates = exact or [
                unit for declared_sample, unit in declared.items()
                if _normalise_sample_token(declared_sample) in token
                or token in _normalise_sample_token(declared_sample)
            ]
            candidates = list(dict.fromkeys(candidates))
            if len(candidates) != 1:
                raise ValueError(
                    f"sample '{sample}' maps to {len(candidates)} biological units"
                )
            resolved[sample] = candidates[0]
        return resolved

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
                    rationale=(
                        f"|log2FC| threshold provenance: "
                        f"{result.get('lfc_threshold_provenance', 'default')}. "
                        f"The default ({DEFAULT_LFC_THRESHOLD}) is fixed and "
                        f"prompt-independent; it changes only via an explicit User "
                        f"Checkpoint 3 profile or a versioned global_lfc override "
                        f"(never inferred from the question text)."),
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
        return (
            "Free-text LLM bulk RNA interpretation is disabled by F3 governance; "
            "use structured NarrativeBlock results and evidence cards."
        )

    def _interpret_single(self, result: dict, intent: dict, exp_ctx: dict) -> str:
        return (
            "Free-text LLM bulk RNA interpretation is disabled by F3 governance; "
            "use structured NarrativeBlock results and evidence cards."
        )

    def _output_dir(self, files: list) -> str: return str(Path(files[0]).parent / "aria_bulk_de") if files else "/tmp/aria_bulk_de"

    def receive(self, message): pass
