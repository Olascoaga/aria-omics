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
from aria.llm.prompt_boundary import (
    PromptDataField,
    build_untrusted_prompt,
    system_with_untrusted_boundary,
)
from aria.memory.memory import ARIAMemory
from aria.utils.bulk_qc_contract import sample_qc_failed, sample_qc_size_ratio

log = logging.getLogger("aria.bulk_rna")


def _is_fastq(files: list) -> bool:
    if not files: return False
    return any(str(files[0]).lower().endswith(s) for s in [".fastq.gz", ".fq.gz", ".fastq", ".fq"])


def _preprocessing_threads(exp_ctx: dict) -> int:
    """Resolve an explicit per-run CPU request without oversubscribing the host."""
    try:
        requested = int(exp_ctx.get("n_cpus", 8))
    except (TypeError, ValueError):
        requested = 8
    return max(1, min(requested, os.cpu_count() or requested))


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
        files      = modalities.get("bulk_RNA", []) or modalities.get("bulk_RNA_raw", [])

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

        _group_vals = list(group_labels.values())
        if _group_vals and all(isinstance(v, (list, tuple, set)) for v in _group_vals):
            _detected_groups = sorted(str(g) for g in group_labels.keys())
        else:
            _detected_groups = sorted({str(v) for v in _group_vals})
        self.publish_status(
            experiment_id,
            f"Detected groups: {_detected_groups}",
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
        confirmed_covariates = self._design_covariates(design) if design else []
        metadata_file = (
            str(self._write_design_metadata(
                group_labels,
                design_factor,
                output_dir,
                replicate_units=replicate_units,
                covariates=confirmed_covariates,
                exp_ctx=exp_ctx,
                design=design,
            ))
            if design and group_labels else ""
        )

        self.publish_status(
            experiment_id,
            f"Running {len(contrasts)} contrast(s), padj < {padj_thr}, |log2FC| > {lfc_thr}...",
            0.70,
        )

        from aria.utils.rna_fastq_protocol import resolve_bulk_rna_protocol

        proto_for_de = resolve_bulk_rna_protocol(exp_ctx)
        de_proto = (
            (preprocessing or {}).get("protocol")
            or proto_for_de.get("protocol")
            or "full_length"
        )
        de_suppress_tpm = bool(
            (preprocessing or {}).get(
                "suppress_tpm",
                proto_for_de.get("suppress_tpm", False),
            )
        )
        de_extra_params = {}
        if de_proto != "full_length":
            de_extra_params["library_protocol"] = de_proto
        if de_suppress_tpm:
            de_extra_params["suppress_tpm"] = True
        if proto_for_de.get("umi_config"):
            de_extra_params["umi_config"] = proto_for_de["umi_config"]
        if exp_ctx.get("count_union_policy"):
            de_extra_params["count_union_policy"] = str(exp_ctx["count_union_policy"])

        # D1-H3: never fabricate conflicting human/hg38 defaults when only one
        # of organism/genome is specified; read from exp_ctx, design, or
        # genome_config, and forward an explicit GTF when available.
        genome_cfg = exp_ctx.get("genome_config", {}) or {}
        resolved_organism = (
            exp_ctx.get("organism")
            or (design or {}).get("organism")
            or genome_cfg.get("organism")
            or ""
        )
        resolved_genome = (
            exp_ctx.get("genome")
            or (design or {}).get("genome")
            or genome_cfg.get("genome")
            or ""
        )
        resolved_gtf = (
            exp_ctx.get("gtf")
            or exp_ctx.get("gtf_file")
            or genome_cfg.get("gtf")
            or genome_cfg.get("gtf_file")
            or ""
        )
        species_params: dict = {}
        if resolved_organism:
            species_params["organism"] = resolved_organism
        if resolved_genome:
            species_params["genome"] = resolved_genome
        if resolved_gtf:
            species_params["gtf"] = str(resolved_gtf)

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
                "covariates":     confirmed_covariates,
                "technical_replicate_col": (
                    "biological_unit" if replicate_units else ""
                ),
                **species_params,
                "output_dir":     output_dir,
                "run_pathways":   True,
                "padj_threshold": padj_thr,
                "lfc_threshold":  lfc_thr,
                **({"n_cpus": _preprocessing_threads(exp_ctx)}
                   if "n_cpus" in exp_ctx else {}),
                **de_extra_params,
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

    _GENE_COLS = {
        "geneid", "gene_id", "ensembl", "ensembl_id", "entrez", "entrez_id",
        "", "gene", "gene_name", "symbol", "gene_symbol", "feature", "name",
        "chr", "start", "end", "strand", "length",
    }

    @staticmethod
    def _read_counts_header(path: str | Path) -> list[str]:
        """Read the first non-empty, non-comment header row of a counts file."""
        p = Path(path)
        opener = gzip.open if str(p).endswith(".gz") else open
        header = ""
        with opener(p, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                header = line.rstrip("\r\n")
                break
        if not header:
            return []
        sep = "\t" if "\t" in header else ","
        return [c.strip().strip('"').strip("'") for c in header.split(sep)]

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
        stem_aliases: dict[str, list[str]] = {stem: [stem] for stem in sample_stems}

        def _register_alias(stem: str, raw_val) -> None:
            if stem not in stem_aliases or raw_val in (None, ""):
                return
            vals = raw_val if isinstance(raw_val, (list, tuple, set)) else [raw_val]
            for v in vals:
                tok = str(v).strip()
                if tok and tok not in stem_aliases[stem]:
                    stem_aliases[stem].append(tok)

        for stem in sample_stems:
            _register_alias(stem, sample_aliases.get(stem))

        # Explicit ID maps (e.g. GSM<->SRR) carried by the design layer.
        for map_key in ("sample_id_map", "gsm_to_srr", "srr_to_gsm"):
            id_map = design.get(map_key)
            if not isinstance(id_map, dict):
                continue
            for k, v in id_map.items():
                ks = str(k).strip()
                if ks in sample_stems:
                    _register_alias(ks, v)
                vals = v if isinstance(v, (list, tuple, set)) else [v]
                for item in vals:
                    vs = str(item).strip()
                    if vs in sample_stems:
                        _register_alias(vs, ks)

        for row in design.get("sample_sheet") or []:
            if not isinstance(row, dict):
                continue
            sid = str(row.get("sample") or row.get("id") or "").strip()
            if sid in sample_stems:
                for col_key in (
                    "run_accession", "sra_run", "sra_runs", "sra_accessions",
                    "srr", "gsm", "alias", "aliases",
                ):
                    if row.get(col_key):
                        _register_alias(sid, row[col_key])

        alias_to_group = {}
        alias_to_stems: dict[str, set[str]] = {}
        for stem, grp in sample_stems.items():
            for alias in stem_aliases[stem]:
                alias_to_group[alias] = grp
                alias_to_stems.setdefault(alias, set()).add(stem)

        # Read actual column names from the first count file (gzip + comment-aware)
        try:
            colnames = self._read_counts_header(files[0])
            if not colnames:
                raise ValueError("empty header")
        except Exception as e:
            raise ValueError(f"Cannot read header from counts file {files[0]}: {e}")

        all_candidates = {**sample_stems, **alias_to_group}
        ordered_cols = [c for c in colnames if c.lower() not in self._GENE_COLS]

        def _boundary_matches(token: str, text: str) -> bool:
            if not token:
                return False
            return bool(re.search(
                rf"(?<![A-Za-z0-9]){re.escape(token)}(?![A-Za-z0-9])",
                text,
                flags=re.IGNORECASE,
            ))

        # Map stem→group to column→group: exact (raw or normalized) first so
        # e.g. S1 never collides with S10, then boundary-delimited match, then
        # substring fallback. Multiple aliases resolving to the SAME group are
        # unambiguous.
        group_labels = {}
        matched_stems: set[str] = set()
        for col in ordered_cols:
            col_norm = _normalise_sample_token(col)
            col_leaf = Path(col).name
            col_leaf_stem = re.sub(
                r"\.(?:sorted\.)?(?:bam|sam|cram|fastq\.gz|fq\.gz|fastq|fq|counts|tsv|csv|txt)$",
                "",
                col_leaf,
                flags=re.IGNORECASE,
            )
            col_leaf_norm = _normalise_sample_token(col_leaf_stem)

            exact_hits = [
                (stem, grp)
                for stem, grp in all_candidates.items()
                if stem == col
                or stem == col_leaf
                or stem == col_leaf_stem
                or (
                    _normalise_sample_token(stem)
                    and _normalise_sample_token(stem) in {col_norm, col_leaf_norm}
                )
            ]
            exact_groups = {grp for _, grp in exact_hits}
            if len(exact_groups) == 1:
                group_labels[col] = next(iter(exact_groups))
                for tok, _ in exact_hits:
                    matched_stems.update(alias_to_stems.get(tok, {tok}))
                continue
            if len(exact_groups) > 1:
                log.warning(
                    f"Ambiguous exact column '{col}' matches multiple groups "
                    f"{sorted(exact_groups)} — skipping."
                )
                continue

            boundary_hits = [
                (stem, grp)
                for stem, grp in all_candidates.items()
                if _boundary_matches(str(stem), col)
                or _boundary_matches(str(stem), col_leaf_stem)
            ]
            boundary_groups = {grp for _, grp in boundary_hits}
            if len(boundary_groups) == 1:
                group_labels[col] = next(iter(boundary_groups))
                for tok, _ in boundary_hits:
                    matched_stems.update(alias_to_stems.get(tok, {tok}))
                continue
            if len(boundary_groups) > 1:
                log.warning(
                    f"Ambiguous column '{col}' matches multiple groups "
                    f"{sorted(boundary_groups)} — skipping, will retry with "
                    f"suffix-trimming."
                )
                continue

            sub_hits = [
                (stem, grp)
                for stem, grp in all_candidates.items()
                if (stem and str(stem) in col)
                or (_normalise_sample_token(stem) and _normalise_sample_token(stem) in col_norm)
            ]
            sub_groups = {grp for _, grp in sub_hits}
            if len(sub_groups) == 1:
                group_labels[col] = next(iter(sub_groups))
                for tok, _ in sub_hits:
                    matched_stems.update(alias_to_stems.get(tok, {tok}))
            elif len(sub_groups) > 1:
                log.warning(
                    f"Ambiguous column '{col}' matches multiple groups "
                    f"{sorted(sub_groups)} — skipping, will retry with "
                    f"suffix-trimming."
                )

        # Attempt to match by trimming potential technical replicate suffixes (_1, _2, etc.)
        if not group_labels:
            for col in ordered_cols:
                col_norm = _normalise_sample_token(col)
                col_base = re.sub(r"[_\-]?[12]$", "", col)
                col_base_norm = _normalise_sample_token(col_base)
                exact_trim_hits = [
                    (stem, grp)
                    for stem, grp in all_candidates.items()
                    if col_base == stem
                    or (col_base_norm and col_base_norm == _normalise_sample_token(stem))
                ]
                exact_trim = {grp for _, grp in exact_trim_hits}
                if len(exact_trim) == 1:
                    group_labels[col] = next(iter(exact_trim))
                    for tok, _ in exact_trim_hits:
                        matched_stems.update(alias_to_stems.get(tok, {tok}))
                    continue
                if len(exact_trim) > 1:
                    continue
                fuzzy_trim_hits = [
                    (stem, grp)
                    for stem, grp in all_candidates.items()
                    if (
                        (col_base and (stem.startswith(col_base) or col_base.startswith(stem)))
                        or (
                            col_base_norm
                            and _normalise_sample_token(stem)
                            and (
                                _normalise_sample_token(stem).startswith(col_base_norm)
                                or col_base_norm.startswith(_normalise_sample_token(stem))
                                or _normalise_sample_token(stem) in col_norm
                            )
                        )
                    )
                ]
                fuzzy_trim = {grp for _, grp in fuzzy_trim_hits}
                if len(fuzzy_trim) == 1:
                    group_labels[col] = next(iter(fuzzy_trim))
                    for tok, _ in fuzzy_trim_hits:
                        matched_stems.update(alias_to_stems.get(tok, {tok}))

        # B4 / D7-001 preprint-readiness: NO positional fallback and NO
        # column-name inference rescue here. When the confirmed design cannot
        # be mapped to the matrix columns (including same-length GSM-vs-SRR
        # mismatches without an explicit alias/ID map), fail closed naming the
        # unmatched identifiers on both sides so the caller's P0-6 gate handles
        # it instead of silently flipping contrast signs.
        unmatched_design = [s for s in sample_stems if s not in matched_stems]
        unmatched_cols = [c for c in ordered_cols if c not in group_labels]
        if not group_labels or len(set(group_labels.values())) < 2:
            raise ValueError(
                f"Could not map the confirmed design groups to the count matrix "
                f"columns. Unmatched design IDs: {unmatched_design}; "
                f"Unmatched matrix columns: {unmatched_cols}. "
                f"Design stems: {list(sample_stems.keys())}, "
                f"Columns: {colnames}. ARIA refuses positional fallback or "
                f"column-name inference; provide explicit sample_aliases / "
                f"sample_id_map (e.g. GSM<->SRR) or metadata that maps the "
                f"confirmed design to these columns."
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

        from aria.utils.design_matrix import normalize_comparisons

        contrasts: list[dict] = []
        for comp in normalize_comparisons(raw, keep_name=True):
            num = resolve_level(comp["test"])
            den = resolve_level(comp["reference"])
            if not num or not den:
                continue
            contrasts.append({
                "numerator": num,
                "denominator": den,
                "name": str(comp.get("name") or f"{num} vs {den}"),
            })
        return contrasts

    @staticmethod
    def _load_sample_metadata_rows(source) -> dict[str, dict]:
        """Normalise sample metadata from a dict, list-of-dicts, or TSV/CSV path
        into ``{sample_id: {col: value, ...}}``."""
        if not source:
            return {}
        id_keys = (
            "sample", "sample_id", "sample_name", "library", "library_id",
            "id", "name", "run", "accession", "geo_accession", "",
        )
        rows_by_id: dict[str, dict] = {}
        if isinstance(source, dict):
            # Could be {sample_id: {col: val}}
            for sid, meta in source.items():
                if isinstance(meta, dict):
                    rows_by_id[str(sid).strip()] = dict(meta)
            return rows_by_id
        if isinstance(source, (str, Path)):
            p = Path(source)
            if not p.is_file():
                return {}
            try:
                opener = gzip.open if str(p).endswith(".gz") else open
                with opener(p, "rt", encoding="utf-8-sig", errors="replace") as fh:
                    sample_lines = []
                    for line in fh:
                        if line.strip() and not line.lstrip().startswith("#"):
                            sample_lines.append(line)
                    if not sample_lines:
                        return {}
                sep = "\t" if "\t" in sample_lines[0] else ","
                reader = csv.DictReader(sample_lines, delimiter=sep)
                source = list(reader)
            except (OSError, UnicodeDecodeError, csv.Error) as exc:
                raise ValueError(
                    f"Failed to read declared sample metadata file '{p}': {exc}"
                ) from exc
        if isinstance(source, list):
            for row in source:
                if not isinstance(row, dict):
                    continue
                lower_keys = {str(k or "").strip().lower(): k for k in row.keys()}
                sid_val = None
                for key in id_keys:
                    if key in lower_keys:
                        raw = row.get(lower_keys[key])
                        if raw is not None and str(raw).strip():
                            sid_val = str(raw).strip()
                            break
                if sid_val is None and row:
                    first_k = next(iter(row.keys()))
                    raw = row.get(first_k)
                    if raw is not None and str(raw).strip():
                        sid_val = str(raw).strip()
                if sid_val is not None:
                    rows_by_id[sid_val] = {
                        str(k).strip(): v for k, v in row.items() if k is not None
                    }
        return rows_by_id

    @classmethod
    def _resolve_covariate_values(
        cls,
        sample_names: list[str],
        covariates: list[str],
        exp_ctx: dict | None = None,
        design: dict | None = None,
    ) -> dict[str, dict[str, object]]:
        """Resolve confirmed covariate values per count-matrix sample ID from
        exp_ctx/design metadata sources (BI-02)."""
        if not sample_names or not covariates:
            return {}
        exp_ctx = exp_ctx or {}
        design = design or {}
        inferred = exp_ctx.get("inferred_design") or {}

        raw_sources = [
            exp_ctx.get("sample_metadata"),
            design.get("sample_metadata"),
            inferred.get("sample_sheet"),
            design.get("sample_sheet"),
            exp_ctx.get("metadata_file"),
            exp_ctx.get("sample_metadata_path"),
            design.get("metadata_file"),
        ]
        # CP2.5 paired biological block values confirmed in DesignAgent (e.g. ~ donor + condition)
        paired_covs = design.get("paired_covariate_values")
        if isinstance(paired_covs, dict) and paired_covs:
            synthesized_rows: dict[str, dict[str, object]] = {}
            for cov_name, s_map in paired_covs.items():
                if not cov_name or not isinstance(s_map, dict):
                    continue
                for s_id, val in s_map.items():
                    synthesized_rows.setdefault(str(s_id), {"sample": str(s_id)})[str(cov_name)] = val
            if synthesized_rows:
                raw_sources.insert(0, list(synthesized_rows.values()))
        parsed_sources = [
            rows for src in raw_sources
            if (rows := cls._load_sample_metadata_rows(src))
        ]
        if not parsed_sources:
            return {}

        # Build alias equivalence groups so a count column matching an alias
        # can still look up covariate metadata keyed by the canonical stem or
        # another alias.
        alias_groups: list[set[str]] = []
        sample_aliases = design.get("sample_aliases") or {}
        if isinstance(sample_aliases, dict):
            for stem, aliases in sample_aliases.items():
                grp = {str(stem)}
                if isinstance(aliases, (list, tuple, set)):
                    grp.update(str(a) for a in aliases if a)
                elif aliases:
                    grp.add(str(aliases))
                alias_groups.append(grp)

        def _match_row(sample: str, source_map: dict[str, dict]) -> dict | None:
            if sample in source_map:
                return source_map[sample]
            leaf = Path(sample).name
            leaf_stem = re.sub(
                r"\.(?:sorted\.)?(?:bam|sam|cram|fastq\.gz|fq\.gz|fastq|fq|counts|tsv|csv|txt)$",
                "",
                leaf,
                flags=re.IGNORECASE,
            )
            for cand in (leaf, leaf_stem):
                if cand in source_map:
                    return source_map[cand]
            s_norm = _normalise_sample_token(sample)
            leaf_norm = _normalise_sample_token(leaf_stem)
            for grp in alias_groups:
                grp_norms = {_normalise_sample_token(a) for a in grp}
                if sample in grp or leaf_stem in grp or s_norm in grp_norms or leaf_norm in grp_norms:
                    for alias in grp:
                        if alias in source_map:
                            return source_map[alias]
                    norm_hits = [
                        k for k in source_map
                        if _normalise_sample_token(k) in grp_norms
                    ]
                    if len(norm_hits) == 1:
                        return source_map[norm_hits[0]]
            norm_matches = [
                k for k in source_map
                if _normalise_sample_token(k) in {s_norm, leaf_norm}
                and _normalise_sample_token(k)
            ]
            if len(norm_matches) == 1:
                return source_map[norm_matches[0]]
            boundary_matches = [
                k for k in source_map
                if k and re.search(
                    rf"(?<![A-Za-z0-9]){re.escape(k)}(?![A-Za-z0-9])",
                    sample,
                    flags=re.IGNORECASE,
                )
            ]
            if len(boundary_matches) == 1:
                return source_map[boundary_matches[0]]
            return None

        def _extract_cov(row: dict, cov: str):
            if cov in row:
                return row[cov]
            cov_l = cov.lower()
            for k, v in row.items():
                if str(k).strip().lower() == cov_l:
                    return v
            return None

        def _is_missing(val) -> bool:
            if val is None:
                return True
            if isinstance(val, float):
                import math
                return not math.isfinite(val)
            return str(val).strip().lower() in ("", "nan", "none")

        resolved: dict[str, dict[str, object]] = {}
        for cov in dict.fromkeys(covariates):
            if not cov:
                continue
            per_sample: dict[str, object] = {}
            for sample in sample_names:
                for source_map in parsed_sources:
                    row = _match_row(sample, source_map)
                    if row is None:
                        continue
                    val = _extract_cov(row, cov)
                    if not _is_missing(val):
                        per_sample[sample] = val
                        break
            if per_sample:
                resolved[cov] = per_sample
        return resolved

    @classmethod
    def _write_design_metadata(
        cls,
        group_labels: dict[str, str],
        design_factor: str,
        output_dir: str | Path,
        *,
        replicate_units: dict[str, str] | None = None,
        covariates: list[str] | None = None,
        exp_ctx: dict | None = None,
        design: dict | None = None,
        covariate_values: dict[str, dict[str, object]] | None = None,
    ) -> Path:
        """Persist the confirmed sample-to-group mapping (and any resolved
        confirmed covariates) for rna_bulk_de.py."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "confirmed_design_metadata.tsv"

        reserved = {"sample", design_factor}
        if replicate_units:
            reserved.add("biological_unit")

        cov_list = [
            c for c in dict.fromkeys(covariates or [])
            if c and c not in reserved
        ]
        cov_data: dict[str, dict[str, object]] = {}
        if cov_list and (exp_ctx is not None or design is not None):
            cov_data.update(
                cls._resolve_covariate_values(
                    list(group_labels.keys()),
                    cov_list,
                    exp_ctx=exp_ctx,
                    design=design,
                )
            )
        if covariate_values:
            for cov, mapping in covariate_values.items():
                if cov and cov not in reserved and isinstance(mapping, dict):
                    cov_data[cov] = dict(mapping)
                    if cov not in cov_list:
                        cov_list.append(cov)

        active_covs = [c for c in cov_list if c in cov_data]

        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            header = ["sample", design_factor]
            if replicate_units:
                header.append("biological_unit")
            header.extend(active_covs)
            writer.writerow(header)
            for sample, group in group_labels.items():
                row = [sample, group]
                if replicate_units:
                    row.append(replicate_units[sample])
                for cov in active_covs:
                    val = cov_data[cov].get(sample, "")
                    row.append("" if val is None else val)
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

    @classmethod
    def _read_sample_names(cls, path: str) -> list[str]:
        try:
            p = Path(path)
            if not p.exists():
                return []
            cols = cls._read_counts_header(p)
            return [c for c in cols if c.lower().strip() not in cls._GENE_COLS]
        except Exception:
            return []

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
        prompt = build_untrusted_prompt(
            task=(
                "Match biological entity names to the available data labels. "
                "Use only labels present in the data."
            ),
            fields=[
                PromptDataField("entities", entities, "identifier_list",
                                "parsed_intent"),
                PromptDataField("biological_question", intent.get(
                    "summary", ""
                ), "user_text", "user"),
                PromptDataField("labels", group_names, "identifier_list",
                                "sample_metadata"),
            ],
            response_contract=(
                "Return only JSON mapping entity to label, for example: "
                '{"entity_name": "label_name"}.'
            ),
        )
        result = self.think_structured(
            prompt,
            system=system_with_untrusted_boundary(
                "Match biological names to data labels."
            ),
            schema_hint="Return JSON mapping entity to label.",
        )
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
            stage_to_cp = {
                "Differential expression (DESeq2)": 1,
                "PCA + MDS (sample-level structure)": 2,
                "Heatmap (padj top 50)": 3,
                "Heatmap (|log2FC| top 50)": 4,
                "Pathway enrichment (ORA)": 5,
                "GSEA (pre-ranked)": 6,
                "TPM (supplementary export)": 7,
                "CPM (supplementary export; TPM suppressed)": 7,
            }

            for d in decisions:
                step, cp = d.get("step", ""), stage_to_cp.get(d.get("step", ""), 0)
                decision_summary = f"{d.get('input','?')} | {d.get('normalization','?')} | {d.get('gene_filter','?')}"
                try:
                    self.memory.store_decision(
                        decision_id=f"{experiment_id[:8]}-auto-{cp:02d}", wing_id=experiment_id, checkpoint=cp,
                        question=step, decision=decision_summary, rationale=d.get("justification", "")[:500], made_by="bulk_rna_agent (auto)"
                    )
                except Exception as e:  # aria: allow-broad-except surface non-fatal decision persistence error
                    log.warning(f"Failed to store decision '{step}': {e}")
                    if isinstance(result, dict):
                        result.setdefault("warnings", []).append(
                            f"Decision persistence failed for '{step}': {type(e).__name__}: {e}"
                        )

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
            except Exception as e:  # aria: allow-broad-except surface non-fatal threshold decision persistence error
                log.warning(f"Failed to store threshold decision: {e}")
                if isinstance(result, dict):
                    result.setdefault("warnings", []).append(
                        f"Threshold decision persistence failed: {type(e).__name__}: {e}"
                    )

            try:
                contrasts_list = [
                    c for c in (result.get("contrasts") or []) if isinstance(c, dict)
                ]
                fitted_formulas = list(
                    dict.fromkeys(
                        str(c["fitted_design_formula"]).strip()
                        for c in contrasts_list
                        if str(c.get("fitted_design_formula") or "").strip()
                    )
                )
                design_decision = (
                    fitted_formulas[0]
                    if len(fitted_formulas) == 1
                    else (
                        "; ".join(fitted_formulas)
                        if fitted_formulas
                        else str(result.get("design_used") or "not reported")
                    )
                )
                covs_adjusted = list(
                    dict.fromkeys(
                        str(cov)
                        for c in contrasts_list
                        for cov in (c.get("covariates_adjusted") or [])
                        if str(cov or "").strip()
                    )
                )
                covs_dropped = list(
                    dict.fromkeys(
                        str(cov)
                        for c in contrasts_list
                        for cov in (c.get("covariates_dropped") or [])
                        if str(cov or "").strip()
                    )
                )
                if covs_adjusted:
                    design_rationale = (
                        f"Multi-factor DESeq2 design adjusting for confirmed covariate(s): "
                        f"{', '.join(covs_adjusted)}."
                    )
                    if covs_dropped:
                        design_rationale += (
                            f" Dropped unestimable/confounded covariate(s): "
                            f"{', '.join(covs_dropped)}."
                        )
                elif covs_dropped:
                    design_rationale = (
                        f"Single-factor DESeq2 design; candidate covariate(s) dropped "
                        f"before fitting: {', '.join(covs_dropped)}."
                    )
                else:
                    design_rationale = (
                        "Single-factor DESeq2 design with no additional batch or "
                        "covariate terms fitted."
                    )
                self.memory.store_decision(
                    decision_id=f"{experiment_id[:8]}-auto-00-design",
                    wing_id=experiment_id,
                    checkpoint=0,
                    question="DESeq2 design formula",
                    decision=design_decision,
                    rationale=design_rationale,
                    made_by="bulk_rna_agent (auto)",
                )
            except Exception as e:  # aria: allow-broad-except surface non-fatal design decision persistence error
                log.warning(f"Failed to store design decision: {e}")
                if isinstance(result, dict):
                    result.setdefault("warnings", []).append(
                        f"Design decision persistence failed: {type(e).__name__}: {e}"
                    )
        except Exception as e:  # aria: allow-broad-except surface non-fatal bulk decision logging failure
            log.warning(f"Decision logging failed (non-fatal): {e}")
            if isinstance(result, dict):
                result.setdefault("warnings", []).append(
                    f"Decision logging failed: {type(e).__name__}: {e}"
                )

    def _run_preprocessing(self, experiment_id: str, fastq_files: list, exp_ctx: dict, intent: dict) -> tuple:
        from aria.utils.rna_fastq_protocol import resolve_bulk_rna_protocol

        fastq_dir = str(Path(fastq_files[0]).parent)
        output_dir = str(Path(fastq_files[0]).parent.parent / "aria_processing")
        genome_cfg = exp_ctx.get("genome_config", {}) or {}
        threads = _preprocessing_threads(exp_ctx)

        proto = resolve_bulk_rna_protocol(exp_ctx)
        if proto["status"] == "invalid":
            err_msg = "; ".join(proto["errors"])
            self.publish_finding(
                experiment_id,
                {"summary": f"FASTQ protocol validation failed: {err_msg}"},
                Confidence.INSUFFICIENT,
            )
            return None, {
                "status": "error",
                "error_type": "InvalidRnaFastqProtocol",
                "details": err_msg,
            }

        canonical_proto = proto["protocol"]
        star_index_dir = genome_cfg.get("star_index", "")
        if not star_index_dir and genome_cfg.get("fasta"):
            star_index_dir = str(Path(output_dir) / "star_index")

        # ── Family 4: Multiplexed well-barcode + UMI plate bulk (STARsolo) ───
        if canonical_proto == "pooled_plate_bulk":
            self.publish_status(
                experiment_id,
                "Demultiplexing & UMI-quantifying pooled plate FASTQs (STARsolo)...",
                0.15,
            )
            solo_result = self.env.run_in_stack(
                "rnaseq",
                "aria/scripts/rna_starsolo_bulk.py",
                {
                    "fastq_files": fastq_files,
                    "genome_dir": star_index_dir,
                    "genome_fasta": genome_cfg.get("fasta", ""),
                    "gtf_file": genome_cfg.get("gtf", ""),
                    "pooled_plate_manifest": proto["pooled_plate_manifest"],
                    "preset_hint": proto.get("raw_token"),
                    "output_dir": str(Path(output_dir) / "starsolo_bulk"),
                    "threads": threads,
                },
            )
            if solo_result.get("status") == "error":
                return None, solo_result

            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        f"STARsolo pooled bulk quantification complete: "
                        f"{solo_result.get('n_genes', 0):,} genes × "
                        f"{solo_result.get('n_samples', 0)} barcoded samples "
                        f"({solo_result.get('pct_unique', 0.0):.1f}% unique mapping, "
                        "1MM_Directional UMI deduplicated)."
                    )
                },
                Confidence.HIGH,
            )
            return [solo_result.get("counts_matrix")], {
                "protocol": "pooled_plate_bulk",
                "suppress_tpm": True,
                "starsolo": solo_result,
                "quantification": solo_result,
            }

        # ── Families 1, 2, 3: Per-sample full_length / tag_3prime / tag_5prime / umi_bulk
        qc_params = {
            "fastq_dir": fastq_dir,
            "fastq_files": fastq_files,
            "output_dir": str(Path(output_dir) / "qc"),
            "threads": threads,
        }
        design_groups = ((exp_ctx.get("design") or {}).get("groups") or {})
        if isinstance(design_groups, dict) and design_groups:
            known_samples = [
                str(s)
                for samples in design_groups.values()
                for s in (samples if isinstance(samples, (list, tuple, set)) else [samples])
                if s
            ]
            if known_samples:
                qc_params["known_samples"] = known_samples
        if canonical_proto != "full_length":
            qc_params["library_protocol"] = canonical_proto
        if proto.get("umi_config"):
            qc_params["umi_config"] = proto["umi_config"]
        for opt_key in ("trim_poly_x", "trim_poly_g", "trim_front1", "r1_only", "min_length", "quality"):
            if opt_key in genome_cfg:
                qc_params[opt_key] = genome_cfg[opt_key]

        self.publish_status(experiment_id, "Trimming reads (fastp)...", 0.05)
        qc_result = self.env.run_in_stack(
            "rnaseq",
            "aria/scripts/rna_fastq_qc.py",
            qc_params,
        )
        if qc_result.get("status") == "error":
            return None, qc_result

        self._publish_fastq_qc_findings(experiment_id, qc_result)
        self.publish_status(experiment_id, f"QC complete: {qc_result.get('n_samples',0)} samples trimmed", 0.20)

        align_params = {
            "samples": qc_result.get("samples", []),
            "genome_dir": star_index_dir,
            "genome_fasta": genome_cfg.get("fasta", ""),
            "gtf_file": genome_cfg.get("gtf", ""),
            "output_dir": str(Path(output_dir) / "aligned"),
            "threads": threads,
            "two_pass": (
                bool(genome_cfg["two_pass"])
                if "two_pass" in genome_cfg
                else (canonical_proto == "full_length")
            ),
        }
        if canonical_proto != "full_length":
            align_params["library_protocol"] = canonical_proto

        self.publish_status(experiment_id, "Aligning to genome (STAR)...", 0.25)
        align_result = self.env.run_in_stack(
            "rnaseq",
            "aria/scripts/rna_align.py",
            align_params,
        )
        if align_result.get("status") == "error":
            return None, align_result

        self._publish_alignment_findings(experiment_id, align_result)
        self.publish_status(experiment_id, f"Alignment complete: {align_result.get('n_aligned', 0)} samples mapped", 0.50)

        bams_for_quant = align_result.get("bam_files", [])
        dedup_result = None
        if canonical_proto == "umi_bulk":
            self.publish_status(
                experiment_id,
                "Deduplicating inline UMIs (directional 1-MM)...",
                0.55,
            )
            dedup_result = self.env.run_in_stack(
                "rnaseq",
                "aria/scripts/rna_umi_dedup.py",
                {
                    "bam_files": bams_for_quant,
                    "output_dir": str(Path(output_dir) / "umi_dedup"),
                },
            )
            if dedup_result.get("status") == "error":
                return None, dedup_result
            bams_for_quant = dedup_result.get("bam_files", [])
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        f"UMI deduplication complete: "
                        f"{dedup_result.get('n_deduplicated', 0)} samples "
                        f"(mean PCR duplicate collapse "
                        f"{dedup_result.get('mean_umi_dedup_rate_pct', 0.0):.1f}%)."
                    )
                },
                Confidence.HIGH,
            )

        self.publish_status(experiment_id, "Counting reads (featureCounts)...", 0.60)
        quant_result = self.env.run_in_stack(
            "rnaseq",
            "aria/scripts/rna_quantify.py",
            {
                "bam_files": bams_for_quant,
                "gtf_file": genome_cfg.get("gtf", ""),
                "output_dir": str(Path(output_dir) / "counts"),
                "threads": threads,
                "strand": genome_cfg.get("strand", "auto"),
                **(
                    {"min_strand_reads": int(genome_cfg["min_strand_reads"])}
                    if "min_strand_reads" in genome_cfg
                    else {}
                ),
            },
        )
        if quant_result.get("status") == "error":
            return None, quant_result

        self.publish_finding(experiment_id, {"summary": f"Quantification complete: {quant_result.get('n_genes',0):,} genes × {quant_result.get('n_samples',0)} samples"}, Confidence.HIGH)
        preproc_payload = {
            "protocol": canonical_proto,
            "suppress_tpm": proto["suppress_tpm"],
            "qc": qc_result,
            "alignment": align_result,
            "quantification": quant_result,
        }
        if dedup_result is not None:
            preproc_payload["umi_dedup"] = dedup_result
        return [quant_result.get("counts_matrix")], preproc_payload

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
            ratio = sample_qc_size_ratio(qc)
            ratio_txt = f"{ratio:.1f}x" if ratio is not None else "not available"
            if sample_qc_failed(qc):
                err_txt = str(qc.get("error") or "QC failed")[:120]
                self.publish_finding(
                    experiment_id,
                    {
                        "summary": (
                            f"Bulk QC failed ({err_txt}): {qc.get('n_samples', '?')} samples. "
                            f"Lib size range: {ratio_txt}. Outlier screening not available."
                        ),
                        "status": "error",
                        "error": err_txt,
                    },
                    Confidence.INSUFFICIENT,
                )
            else:
                outliers = qc.get("candidate_outliers", qc.get("outliers", [])) or []
                self.publish_finding(
                    experiment_id,
                    {
                        "summary": (
                            f"Bulk QC: {qc.get('n_samples', '?')} samples. "
                            f"Lib size range: {ratio_txt}."
                            + (f" Outliers: {outliers}." if outliers else "")
                        )
                    },
                    Confidence.HIGH if not outliers else Confidence.MEDIUM,
                )

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
