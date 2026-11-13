"""
ARIA IntegrationAgent
---------------------
Synthesizes findings across omics modalities into a coherent biological picture.

Three integration strategies (in order of priority):

  1. WNN (Weighted Nearest Neighbors)
     - For: 2 modalities from the SAME cells (multiome: RNA + ATAC)
     - What it does: builds a joint embedding weighting each modality
       by its information content per cell
     - Output: shared cell embedding, modality weights per cell,
       joint clusters that reflect both transcriptome and chromatin
     - ParameterAdvisor: advises on k (neighbors)
     - DebateCouncil: invoked when modality weights are imbalanced
       (e.g. RNA=0.91, ATAC=0.09 suggests ATAC QC failure)

  2. Peak-to-gene linking
     - For: any experiment with both RNA and ATAC data
     - What it does: correlates peak accessibility with nearby gene
       expression across cells, identifying regulatory relationships
     - Output: peak-gene links stored as Tunnels in ARIAMemory
       (this is the answer to "what regulates this gene?")
     - DebateCouncil: invoked for high-confidence links that lack
       supporting CTCF or footprinting evidence

  3. MOFA+ (Multi-Omics Factor Analysis)
     - For: 3+ modalities OR when user wants latent factor decomposition
     - What it does: decomposes variance into shared and modality-specific
       factors, identifies which biological programs drive each factor
     - Output: factor scores per cell, factor loadings per feature,
       variance explained per modality
     - ParameterAdvisor: advises on number of factors
     - DebateCouncil: always invoked for factor interpretation claims

Integration order:
  Checkpoint 2 reveals available modalities
  OrchestratorAgent selects strategy based on modalities present
  IntegrationAgent executes and stores cross-modal links as Tunnels

Cross-modal conflicts:
  When RNA and ATAC disagree (open chromatin but no expression,
  or expressed gene but closed chromatin), the DebateCouncil is
  called to reason about the biology before reporting.
  These conflicts are often the most biologically interesting findings.
"""

from __future__ import annotations

import json
import logging
import uuid
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, MessageType
from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.integration")


INTEGRATION_SYSTEM = """
You are ARIA's IntegrationAgent — a specialist in multi-omics data integration.

Your expertise:
- WNN (Weighted Nearest Neighbors): joint embedding of RNA + ATAC
- Peak-to-gene linking: correlating chromatin accessibility with expression
- Candidate eRegulons: assembling associative TF -> cis peak -> target gene modules
- MOFA+: latent factor decomposition across 3+ modalities
- Cross-modal conflict resolution: reasoning about discordant signals

Critical knowledge:

WNN:
  - Modality weights reflect information content per cell, not quality alone
  - RNA weight > 0.85 often means ATAC data is low quality OR the biological
    question is primarily transcriptional (both are valid — distinguish them)
  - Cells near cluster boundaries often switch modality dominance
  - WNN clusters can differ from RNA-only clusters — this is expected and
    biologically meaningful, not a bug

PEAK-TO-GENE LINKING:
  - Correlation-based linking has high false positive rate
  - Distance cutoff (typically 500kb) must be justified, not arbitrary
  - Positive correlation and negative correlation are both descriptive
    observations until independent evidence supports interpretation
  - Always require at least 2 independent lines of evidence before
    claiming a regulatory relationship

CANDIDATE EREGULONS:
  - Candidate eRegulons combine PWM motif matches, cis peak-to-gene links,
    and TF-to-target expression co-variation across paired cells
  - Report descriptive sign classes (pos_tf_pos_p2g, neg_tf_pos_p2g) only;
    never label a TF as a proven direct activator or repressor without
    orthogonal perturbation or chromatin-binding assays
  - Flag paralogous TFs sharing highly overlapping motif cistromes and
    modules whose correlation collapses after library-depth conditioning

MOFA+:
  - Number of factors is the critical hyperparameter
  - Too few: biological programs are confounded
  - Too many: factors become noise-driven
  - Factor interpretation requires explicit evidence; do not infer biology
    from feature names alone
  - Variance explained per modality tells you which omics drives the biology
  - If one modality explains <10% variance across all factors, question its quality

CROSS-MODAL CONFLICTS:
  Gene expressed but chromatin closed:
    → Discordant evidence; flag it and avoid causal interpretation
  Chromatin open but gene not expressed:
    → Discordant evidence; flag it and avoid causal interpretation
  HiC loop present but no expression link:
    → Structural contact without expression support
    → Requires additional evidence before claiming enhancer-promoter interaction
""".strip()


class IntegrationAgent(BaseAgent):

    name        = "integration_agent"
    description = (
        "Multi-omics integration: WNN (RNA+ATAC), "
        "peak-to-gene linking, MOFA+ for 3+ modalities."
    )

    # Minimum cells required per modality for reliable WNN
    MIN_CELLS_WNN = 200

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.advisor = ParameterAdvisor(memory, llm)

        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Select and execute the appropriate integration strategy.

        context must contain:
          - exp_context:        from DataAuditAgent
          - biological_intent:  from OrchestratorAgent
          - rna_results:        from RNAAgent (optional)
          - chromatin_results:  from ChromatinAgent (optional)
          - hic_results:        from GenomeArchAgent (optional)
        """
        exp_ctx  = context.get("exp_context", {})
        intent   = context.get("biological_intent", {})
        mods     = dict(exp_ctx.get("modalities", {}) or {})

        self.publish_status(experiment_id,
                            "IntegrationAgent: selecting strategy...", 0.0)

        from aria.utils.multiome_contracts import (
            enforce_same_cell_integration_gate,
            infer_multiome_contract,
        )

        multiome_contract = infer_multiome_contract(exp_ctx)
        is_paired_mudata = (
            multiome_contract.get("object_type")
            in {"paired_mudata", "paired_cohort_mudata", "raw_multiome_fastq"}
            and multiome_contract.get("same_cell") is True
        )

        # Determine which modalities are available (a paired .h5mu / ARC /
        # raw Multiome FASTQ contract carries both scRNA and scATAC in a shared
        # cell namespace).
        has_scrna  = ("scRNA" in mods) or is_paired_mudata
        has_scatac = ("scATAC" in mods) or is_paired_mudata
        has_rna    = has_scrna or "bulk_RNA" in mods
        has_atac   = has_scatac or "bulk_ATAC" in mods
        has_hic    = "HiC"      in mods
        has_chip   = any(m in mods for m in ("ChIP", "CUT_AND_RUN", "CUT_AND_TAG"))

        n_modalities = sum([has_rna, has_atac, has_hic, has_chip])

        if n_modalities < 2:
            self.publish_finding(
                experiment_id,
                {"error": "Integration requires at least 2 modalities. "
                          f"Only found: {list(mods.keys())}"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "insufficient_modalities"}

        strategy  = self._select_strategy(
            has_scrna, has_scatac, has_hic, n_modalities, intent,
            exp_ctx=exp_ctx,
        )
        results   = {}
        stage_ctx = dict(exp_ctx)

        # V48 I1 + I8 hard gate (ADR-064 / ADR-071): ALL cross-modal integration
        # strategies MUST satisfy confirmed same-cell RNA+ATAC pairing before
        # touching downstream algorithms (closing any non-WNN bypass).
        gate = enforce_same_cell_integration_gate(exp_ctx)
        results["multiome_gate"] = gate
        if not gate["allowed"]:
            blocker_msgs = "; ".join(
                str(b.get("message") or b.get("check"))
                for b in gate.get("blockers", [])
            )
            self.publish_finding(
                experiment_id,
                {
                    "error": (
                        "Same-cell multi-omics integration blocked by "
                        f"multiome contract gate: {blocker_msgs}"
                    ),
                    "blockers": gate.get("blockers", []),
                },
                Confidence.INSUFFICIENT,
            )
            return {
                "status": "blocked",
                "reason": "multiome_same_cell_gate_blocked",
                "strategy": strategy,
                "multiome_gate": gate,
            }

        # ADR-079: Universal same-cell multiome assembly & multi-donor cohort
        # merger before joint QC (Cell Ranger ARC .h5 / MEX, paired matrices
        # with barcode translation or fragments+peaks, raw Multiome FASTQs, and
        # multi-container cohorts).
        effective_gate = dict(gate)
        pre_res = self._run_multiome_preprocessing(
            experiment_id, exp_ctx, intent, effective_gate, mods
        )
        if pre_res is not None:
            results["preprocessing"] = pre_res
            if pre_res.get("status") != "success":
                self.publish_finding(
                    experiment_id,
                    {
                        "error": (
                            "Multiome preprocessing / cohort assembly failed: "
                            f"{pre_res.get('details', pre_res.get('error_type', 'unknown'))}"
                        ),
                        "preprocessing": pre_res,
                    },
                    Confidence.INSUFFICIENT,
                )
                return {
                    "status": "failed",
                    "reason": "multiome_preprocessing_failed",
                    "strategy": strategy,
                    "findings": results,
                }
            assembled_out = pre_res.get("output_path")
            if assembled_out:
                effective_gate = {
                    **effective_gate,
                    "evidence_paths": [str(assembled_out)],
                    "contract": {
                        **dict(effective_gate.get("contract") or {}),
                        "object_type": "paired_mudata",
                        "evidence_paths": [str(assembled_out)],
                    },
                }
                is_paired_mudata = True
                if pre_res.get("inferred_batch_key") and not stage_ctx.get("batch_key"):
                    stage_ctx["batch_key"] = str(pre_res["inferred_batch_key"])

        if effective_gate["contract"].get("object_type") in {
            "paired_mudata",
            "paired_cohort_mudata",
        }:
            self.publish_status(
                experiment_id,
                "Running paired .h5mu contract & joint RNA+ATAC QC...",
                0.05,
            )
            paired_qc_res = self._run_paired_qc(
                experiment_id, stage_ctx, intent, effective_gate
            )
            results["paired_qc"] = paired_qc_res
            if paired_qc_res.get("status") != "success":
                self.publish_finding(
                    experiment_id,
                    {
                        "error": (
                            "Paired .h5mu joint QC failed: "
                            f"{paired_qc_res.get('details', paired_qc_res.get('error_type', 'unknown'))}"
                        ),
                        "paired_qc": paired_qc_res,
                    },
                    Confidence.INSUFFICIENT,
                )
                return {
                    "status": "failed",
                    "reason": "paired_qc_failed",
                    "strategy": strategy,
                    "findings": results,
                }
            validated_h5mu = paired_qc_res.get("output_path")
            evidence_paths = list(effective_gate.get("evidence_paths") or [])
            raw_exists = bool(
                evidence_paths and Path(str(evidence_paths[0])).exists()
            )
            if not validated_h5mu or (
                raw_exists and not Path(str(validated_h5mu)).is_file()
            ):
                self.publish_finding(
                    experiment_id,
                    {
                        "error": (
                            "Paired .h5mu joint QC succeeded without a valid "
                            f"post-QC output container: {validated_h5mu!r}"
                        ),
                        "paired_qc": paired_qc_res,
                    },
                    Confidence.INSUFFICIENT,
                )
                return {
                    "status": "failed",
                    "reason": "paired_qc_missing_output",
                    "strategy": strategy,
                    "findings": results,
                }
            validated_str = str(validated_h5mu)
            mods = {
                **mods,
                "scRNA": [validated_str],
                "scATAC": [validated_str],
            }
            stage_ctx = {
                **stage_ctx,
                "_validated_paired_h5mu": validated_str,
                "h5mu_path": validated_str,
                "paired_h5mu": validated_str,
            }

        self.publish_status(
            experiment_id,
            f"Integration strategy: {strategy}", 0.1
        )

        # ── Execute selected strategy ─────────────────────────────────────

        if strategy in ("wnn", "wnn_then_mofa"):
            if has_scrna and has_scatac:
                self.publish_status(experiment_id,
                                    "Running WNN (RNA + ATAC)...", 0.2)
                wnn_result = self._run_wnn(
                    experiment_id, stage_ctx, intent, mods
                )
                results["wnn"] = wnn_result

        if has_rna and has_atac:
            self.publish_status(experiment_id,
                                "Linking peaks to genes...", 0.5)
            p2g_result = self._run_peak_to_gene(
                experiment_id, stage_ctx, intent, mods,
                context.get("rna_results", {}),
                context.get("chromatin_results", {}),
            )
            results["peak_to_gene"] = p2g_result

        if has_scrna and has_scatac and self._should_run_eregulons(
            is_paired_mudata, stage_ctx, intent
        ):
            self.publish_status(
                experiment_id,
                "Assembling candidate eRegulons (TF -> cis peak -> gene)...",
                0.6,
            )
            ereg_result = self._run_eregulons(
                experiment_id,
                stage_ctx,
                intent,
                mods,
                results.get("peak_to_gene", {}),
            )
            results["eregulons"] = ereg_result

        if strategy in ("mofa", "wnn_then_mofa") and n_modalities >= 2:
            self.publish_status(experiment_id,
                                "Running MOFA+ factor decomposition...", 0.7)
            mofa_result = self._run_mofa(
                experiment_id, stage_ctx, intent, mods
            )
            results["mofa"] = mofa_result

        # ── Resolve cross-modal conflicts ─────────────────────────────────
        if results:
            self.publish_status(experiment_id,
                                "Resolving cross-modal conflicts...", 0.9)
            conflicts = self._resolve_conflicts(
                experiment_id, stage_ctx, intent, results
            )
            results["conflicts"] = conflicts

        self.publish_status(experiment_id,
                            "IntegrationAgent complete.", 1.0)
        return {"status": "done", "strategy": strategy, "findings": results}

    # ── Strategy selection ────────────────────────────────────────────────

    def _select_strategy(self, has_scrna: bool, has_scatac: bool,
                          has_hic: bool, n_modalities: int,
                          intent: dict,
                          exp_ctx: dict | None = None) -> str:
        """
        Select the optimal integration strategy.

        Decision tree:
          same-cell RNA + ATAC (multiome) → WNN (primary) + peak-to-gene
          3+ modalities OR factor request / enable_mofa → MOFA+
          enable_mofa=False → explicitly disables MOFA+
          RNA + ATAC (non-WNN fallback) → peak-to-gene only
        """
        question = (intent or {}).get("user_question", "").lower()
        ctx = exp_ctx or {}

        FACTOR_KEYWORDS = [
            "factor", "latent", "mofa", "decomposition",
            "shared variation", "modality-specific",
        ]

        enable_mofa = ctx.get("enable_mofa")
        if enable_mofa is False:
            return "wnn" if (has_scrna and has_scatac) else "peak_to_gene_only"

        user_wants_mofa = (
            (enable_mofa is True)
            or bool(ctx.get("mofa_params"))
            or any(k in question for k in FACTOR_KEYWORDS)
        )

        if has_scrna and has_scatac:
            if n_modalities >= 3 or user_wants_mofa:
                return "wnn_then_mofa"
            return "wnn"
        elif n_modalities >= 3 or user_wants_mofa:
            return "mofa"
        else:
            return "peak_to_gene_only"

    # ── Universal Multiome Intake & Cohort Merger (ADR-079) ───────────────

    def _run_multiome_preprocessing(
        self,
        experiment_id: str,
        exp_ctx: dict,
        intent: dict,
        gate: dict,
        mods: dict,
    ) -> dict | None:
        """Materialize raw FASTQs / ARC matrices / multi-donor cohorts into `.h5mu`."""
        contract = dict(gate.get("contract") or {})
        obj_type = str(contract.get("object_type") or "")
        evidence_paths = [str(p) for p in (gate.get("evidence_paths") or []) if p]
        assemble_overrides = dict(exp_ctx.get("assemble_params") or {})
        merge_overrides = dict(exp_ctx.get("cohort_merge_params") or {})
        genome_cfg = dict(exp_ctx.get("genome_config") or {})

        # 1. Raw Multiome FASTQ manifest (single-sample or multi-sample cohort)
        if obj_type == "raw_multiome_fastq":
            self.publish_status(
                experiment_id,
                "Assembling raw Multiome FASTQs into paired .h5mu...",
                0.03,
            )
            manifest = exp_ctx.get("multiome_fastq_manifest") or {}
            samples = list(
                (manifest.get("samples") if isinstance(manifest, dict) else None)
                or (manifest.get("libraries") if isinstance(manifest, dict) else None)
                or []
            )
            if not samples:
                samples = [
                    {
                        "sample_id": (contract.get("sample_ids") or ["sample_1"])[0],
                        "rna_fastqs": list(mods.get("scRNA") or []),
                        "atac_fastqs": list(mods.get("scATAC") or []),
                        "platform": contract.get("platform") or "10x_multiome_arc_v1",
                    }
                ]
            base_out = Path(
                str(
                    assemble_overrides.get("output_dir")
                    or self._get_output_dir(evidence_paths or ["multiome"])
                )
            )
            assembled_paths: list[str] = []
            sample_ids: list[str] = []
            sample_meta: dict[str, dict] = dict(exp_ctx.get("sample_metadata") or {})
            last_asm: dict = {}
            all_rna_fqs: list[str] = []
            all_atac_fqs: list[str] = []

            for idx, row in enumerate(samples):
                if not isinstance(row, dict):
                    continue
                sid = str(row.get("sample_id") or f"sample_{idx + 1}")
                sample_ids.append(sid)
                rna_fqs = [str(p) for p in (row.get("rna_fastqs") or []) if p]
                atac_fqs = [str(p) for p in (row.get("atac_fastqs") or []) if p]
                all_rna_fqs.extend(rna_fqs)
                all_atac_fqs.extend(atac_fqs)
                if isinstance(row.get("metadata"), dict):
                    sample_meta.setdefault(sid, {}).update(row["metadata"])
                if row.get("donor_id"):
                    sample_meta.setdefault(sid, {})["donor_id"] = row["donor_id"]
                if row.get("condition"):
                    sample_meta.setdefault(sid, {})["condition"] = row["condition"]

                asm_params: dict = {
                    "mode": "raw_multiome_fastq",
                    "platform": (
                        row.get("platform")
                        or (manifest.get("platform") if isinstance(manifest, dict) else None)
                        or contract.get("platform")
                        or "10x_multiome_arc_v1"
                    ),
                    "rna_fastqs": rna_fqs,
                    "atac_fastqs": atac_fqs,
                    "output_dir": str(base_out / f"assembled_{sid}"),
                }
                for key in (
                    "gex_whitelist_path",
                    "atac_whitelist_path",
                    "barcode_map_path",
                    "rna_counts_out",
                    "atac_fragments_out",
                    "peaks_path",
                    "kb_index",
                    "t2g_path",
                    "genome_fasta",
                    "chromap_index",
                ):
                    val = (
                        row.get(key)
                        or (manifest.get(key) if isinstance(manifest, dict) else None)
                        or contract.get(key)
                    )
                    if val:
                        asm_params[key] = val
                if not asm_params.get("genome_fasta") and genome_cfg.get("fasta"):
                    asm_params["genome_fasta"] = genome_cfg["fasta"]
                asm_params.update(
                    {k: v for k, v in assemble_overrides.items() if k != "output_dir"}
                )
                asm_res = self.env.run_in_stack(
                    stack="chromatin",
                    script_path="aria/scripts/integration_assemble_multiome.py",
                    params=asm_params,
                )
                last_asm = asm_res
                if asm_res.get("status") != "success" or not asm_res.get("output_path"):
                    return {
                        **asm_res,
                        "mode": "raw_multiome_fastq",
                        "rna_fastqs": all_rna_fqs,
                        "atac_fastqs": all_atac_fqs,
                        "environment_stack": "chromatin",
                    }
                assembled_paths.append(str(asm_res["output_path"]))

            if len(assembled_paths) == 1:
                return {
                    **last_asm,
                    "n_samples": 1,
                    "sample_ids": sample_ids,
                    "rna_fastqs": all_rna_fqs,
                    "atac_fastqs": all_atac_fqs,
                    "environment_stack": "chromatin",
                }

            merge_params = {
                "inputs": assembled_paths,
                "sample_ids": sample_ids,
                "sample_metadata": sample_meta,
                "output_dir": str(base_out / "cohort_merge"),
                **merge_overrides,
            }
            merge_res = self.env.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/integration_merge_cohort.py",
                params=merge_params,
            )
            if merge_res.get("status") != "success":
                return {
                    **merge_res,
                    "mode": "raw_multiome_fastq_cohort",
                    "assembly": last_asm,
                    "rna_fastqs": all_rna_fqs,
                    "atac_fastqs": all_atac_fqs,
                    "environment_stack": "chromatin",
                }
            return {
                **merge_res,
                "mode": "raw_multiome_fastq_cohort",
                "platform": last_asm.get("platform"),
                "assembly": {
                    **last_asm,
                    "artifacts": assembled_paths,
                    "rna_fastqs": all_rna_fqs,
                    "atac_fastqs": all_atac_fqs,
                },
                "cohort_merge": merge_res,
                "rna_fastqs": all_rna_fqs,
                "atac_fastqs": all_atac_fqs,
                "inferred_batch_key": "donor_id",
                "environment_stack": "chromatin",
            }

        # 2. Multi-donor / multi-sample paired containers -> cohort union
        if obj_type in {"paired_mudata", "paired_cohort_mudata"} and (
            len(evidence_paths) >= 2 or contract.get("cohort_merge_required") is True
        ):
            self.publish_status(
                experiment_id,
                f"Merging {len(evidence_paths)} paired multiome samples onto consensus peaks...",
                0.03,
            )
            base_out = Path(self._get_output_dir(evidence_paths)) / "cohort_merge"
            merge_params = {
                "inputs": evidence_paths,
                "output_dir": str(base_out),
            }
            sids = exp_ctx.get("sample_ids") or contract.get("sample_ids")
            if sids:
                merge_params["sample_ids"] = list(sids)
            if exp_ctx.get("sample_metadata"):
                merge_params["sample_metadata"] = dict(exp_ctx["sample_metadata"])
            merge_params.update(merge_overrides)
            merge_res = self.env.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/integration_merge_cohort.py",
                params=merge_params,
            )
            return {
                **merge_res,
                "mode": "paired_cohort_merge",
                "platform": contract.get("platform") or "paired_multiome",
                "barcode_translation": "shared_mudata_obs",
                "inputs": evidence_paths,
                "cohort_merge": merge_res,
                "inferred_batch_key": "donor_id",
                "environment_stack": "chromatin",
            }

        # 3. Single Cell Ranger ARC `.h5` or ARC dual MEX directory -> `.h5mu`
        if (
            obj_type == "paired_mudata"
            and len(evidence_paths) == 1
            and not evidence_paths[0].lower().endswith(".h5mu")
            and Path(evidence_paths[0]).exists()
        ):
            self.publish_status(
                experiment_id,
                "Assembling Cell Ranger ARC container into verified .h5mu...",
                0.03,
            )
            base_out = Path(self._get_output_dir(evidence_paths)) / "assembled_arc"
            asm_params = {
                "mode": "cellranger_arc",
                "arc_path": evidence_paths[0],
                "platform": contract.get("platform") or "10x_multiome_arc_v1",
                "output_dir": str(base_out),
                **assemble_overrides,
            }
            asm_res = self.env.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/integration_assemble_multiome.py",
                params=asm_params,
            )
            return {
                **asm_res,
                "arc_path": evidence_paths[0],
                "n_samples": 1,
                "environment_stack": "chromatin",
            }

        # 4. Explicit same-cell paired matrices requiring barcode translation
        #    or fragments+peaks interval quantification
        if obj_type in {"explicit_same_cell", "explicit_pairing"}:
            ns = str(
                contract.get("cell_namespace")
                or contract.get("barcode_namespace")
                or ""
            )
            rna_files = [str(f) for f in (mods.get("scRNA") or []) if f]
            atac_files = [str(f) for f in (mods.get("scATAC") or []) if f]
            has_frag_peaks = any(
                "fragment" in Path(f).name.lower() or f.lower().endswith(".tsv.gz")
                for f in atac_files
            )
            needs_assembly = bool(
                assemble_overrides
                or ns in {"10x_arc_whitelist_pair", "barcode_map_tsv"}
                or contract.get("barcode_map")
                or contract.get("gex_whitelist_path")
                or has_frag_peaks
            )
            if needs_assembly and rna_files and atac_files:
                self.publish_status(
                    experiment_id,
                    "Assembling paired RNA + ATAC matrices over verified cell namespace...",
                    0.03,
                )
                base_out = (
                    Path(self._get_output_dir(rna_files or atac_files))
                    / "assembled_paired"
                )
                asm_params = {
                    "mode": "paired_matrices",
                    "rna_path": rna_files[0],
                    "platform": contract.get("platform") or "paired_multiome",
                    "output_dir": str(base_out),
                }
                if has_frag_peaks:
                    frags = [
                        f
                        for f in atac_files
                        if "fragment" in Path(f).name.lower()
                        or f.lower().endswith(".tsv.gz")
                    ]
                    peaks = [
                        f
                        for f in atac_files
                        if f not in frags
                        and Path(f).name.lower().endswith(
                            (".bed", ".bed.gz", ".narrowpeak", ".broadpeak")
                        )
                    ]
                    if frags:
                        asm_params["fragments_path"] = frags[0]
                    if peaks:
                        asm_params["peaks_path"] = peaks[0]
                else:
                    asm_params["atac_path"] = atac_files[0]

                if ns == "10x_arc_whitelist_pair" or contract.get("gex_whitelist_path"):
                    asm_params["barcode_translation"] = "10x_arc_whitelist_pair"
                    if contract.get("gex_whitelist_path"):
                        asm_params["gex_whitelist_path"] = contract["gex_whitelist_path"]
                    if contract.get("atac_whitelist_path"):
                        asm_params["atac_whitelist_path"] = contract["atac_whitelist_path"]
                elif ns == "barcode_map_tsv" or contract.get("barcode_map"):
                    asm_params["barcode_translation"] = "barcode_map_tsv"
                    if contract.get("barcode_map"):
                        asm_params["barcode_map_path"] = contract["barcode_map"]

                asm_params.update(assemble_overrides)
                asm_res = self.env.run_in_stack(
                    stack="chromatin",
                    script_path="aria/scripts/integration_assemble_multiome.py",
                    params=asm_params,
                )
                return {
                    **asm_res,
                    "inputs": [rna_files[0], *atac_files],
                    "n_samples": 1,
                    "environment_stack": "chromatin",
                }

        return None

    # ── Paired .h5mu contract & joint QC (V48 I1) ─────────────────────────

    def _run_paired_qc(
        self,
        experiment_id: str,
        exp_ctx: dict,
        intent: dict,
        gate: dict,
    ) -> dict:
        """Validate paired `.h5mu` and run joint RNA+ATAC QC over shared cells."""
        evidence_paths = list(gate.get("evidence_paths") or [])
        if not evidence_paths:
            return {
                "status": "error",
                "error_type": "MissingPairedH5mu",
                "details": "No paired .h5mu path resolved from multiome contract.",
            }
        if len(evidence_paths) > 1:
            return {
                "status": "error",
                "error_type": "MultiplePairedMuDataUnsupported",
                "details": (
                    f"Found {len(evidence_paths)} paired .h5mu containers; "
                    "merge shared-donor containers first or select a single .h5mu."
                ),
            }

        h5mu_path = evidence_paths[0]
        qc_overrides = dict(exp_ctx.get("paired_qc_params") or {})
        params = {
            "h5mu_path": h5mu_path,
            "output_dir": qc_overrides.pop(
                "output_dir",
                self._get_output_dir(evidence_paths),
            ),
            "organism": exp_ctx.get("organism", "Homo sapiens"),
            "genome": exp_ctx.get("genome"),
        }
        for opt_key in (
            "sample_id",
            "reference_h5ad",
            "reference_label_col",
            "target_label_col",
            "reference_sample_col",
        ):
            if opt_key in exp_ctx and exp_ctx[opt_key] is not None:
                params[opt_key] = exp_ctx[opt_key]
        params.update(qc_overrides)

        qc_res = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/integration_paired_qc.py",
            params=params,
        )
        if qc_res.get("status") == "success":
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        f"Paired multiome joint QC retained "
                        f"{qc_res.get('n_cells_after')}/{qc_res.get('n_cells_before')} "
                        f"shared cells ({qc_res.get('n_genes')} genes, "
                        f"{qc_res.get('n_peaks')} peaks)."
                    ),
                    "n_cells_before": qc_res.get("n_cells_before"),
                    "n_cells_after": qc_res.get("n_cells_after"),
                    "label_join": qc_res.get("label_join"),
                },
                Confidence.HIGH,
            )
        return qc_res

    # ── WNN integration ───────────────────────────────────────────────────

    def _run_wnn(self, experiment_id: str, exp_ctx: dict,
                  intent: dict, mods: dict) -> dict:
        """
        WNN integration of scRNA-seq + scATAC-seq from the same cells.
        Uses ParameterAdvisor to select k (number of neighbors).
        Invokes DebateCouncil when modality weights are imbalanced.
        """
        # ParameterAdvisor for WNN k
        k_decision = self.advisor.advise_wnn_k(
            multimodal_adata=None,  # real adata loaded in script
            experiment_id=experiment_id,
            biological_context=intent,
        )

        # Escalate to user (Checkpoint 3) and block until the dispatch-loop UI
        # resolves it. Otherwise the advisor value would run before user input.
        _, user_decision = self.publish_blocking_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=self._format_wnn_checkpoint(k_decision),
            options=[
                f"Use recommended k={k_decision.chosen_value}",
                "Other: Enter custom k value",
                "Skip WNN integration",
            ],
            context={
                "wnn_k_decision": k_decision.decision_id,
                "analysis_type": "wnn_integration",
                "parameter_name": "k",
            },
        )
        choice = self.checkpoint_choice_text(user_decision)
        if not choice:
            return {
                "status": "error",
                "error_type": "CheckpointUnresolved",
                "details": "WNN k checkpoint was not resolved.",
            }
        if "skip" in choice.lower():
            k_decision.chosen_by = "user"
            k_decision.approved_by_user = True
            return {
                "status": "skipped",
                "reason": "user_skipped_wnn",
                "details": "User skipped WNN integration at checkpoint 3.",
            }
        override = self.numeric_checkpoint_override(user_decision, int)
        if hasattr(self.advisor, "approve_decision"):
            k_decision = self.advisor.approve_decision(
                k_decision,
                user_override=override,
            )
        else:
            if override is not None:
                k_decision.chosen_value = override
                k_decision.chosen_by = "user"
            k_decision.approved_by_user = True

        # Execute WNN in isolated environment (aria-chromatin-env provides the
        # locked muon + mudata + scanpy + leidenalg + umap-learn stack).
        rna_files  = list(mods.get("scRNA", []) or [])
        atac_files = list(mods.get("scATAC", []) or [])
        wnn_overrides = dict(exp_ctx.get("wnn_params") or {})

        h5mu_path = wnn_overrides.pop("h5mu_path", None)
        if not h5mu_path:
            for cand in [*rna_files, *atac_files]:
                if str(cand).lower().endswith(".h5mu"):
                    h5mu_path = str(cand)
                    break

        ref_paths = [h5mu_path] if h5mu_path else (rna_files or atac_files)
        wnn_params = {
            "rna_files":   rna_files,
            "atac_files":  atac_files,
            "genome":      exp_ctx.get("genome", "hg38"),
            "organism":    exp_ctx.get("organism", "Homo sapiens"),
            "k":           k_decision.chosen_value,
            "seed":        int(exp_ctx.get("seed", 0)),
            "output_dir":  wnn_overrides.pop(
                "output_dir",
                self._get_output_dir(ref_paths),
            ),
        }
        if h5mu_path:
            wnn_params["h5mu_path"] = h5mu_path
        if exp_ctx.get("batch_key"):
            wnn_params["batch_key"] = str(exp_ctx["batch_key"])
        wnn_params.update(wnn_overrides)
        validated_paired = exp_ctx.get("_validated_paired_h5mu")
        if validated_paired:
            wnn_params["h5mu_path"] = str(validated_paired)
            wnn_params["rna_files"] = [str(validated_paired)]
            wnn_params["atac_files"] = [str(validated_paired)]

        wnn_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/integration_wnn.py",
            params=wnn_params,
        )

        if wnn_result.get("status") == "success":
            self._evaluate_wnn_weights(
                experiment_id, wnn_result, intent
            )

        return wnn_result

    def _run_wnn_diagnostics(
        self,
        h5mu_path: str,
        *,
        seed: int = 0,
        audit_group_col: str | None = None,
        update_h5mu: bool = False,
        output_json: str | None = None,
    ) -> dict:
        """Execute standalone V48 Slice I3 WNN anti-degradation diagnostics."""
        diag_params: dict = {
            "h5mu_path": str(h5mu_path),
            "seed": int(seed),
            "update_h5mu": bool(update_h5mu),
        }
        if audit_group_col:
            diag_params["audit_group_col"] = str(audit_group_col)
        if output_json:
            diag_params["output_json"] = str(output_json)
        return self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/integration_wnn_diagnostics.py",
            params=diag_params,
        )

    def _evaluate_wnn_weights(self, experiment_id: str,
                               wnn_result: dict, intent: dict):
        """
        Evaluate WNN modality weights and V48 Slice I3 anti-degradation
        diagnostics, invoking DebateCouncil when imbalanced weights suggest
        data quality issues.
        """
        rna_weight  = wnn_result.get("mean_rna_weight")
        atac_weight = wnn_result.get("mean_atac_weight")
        anti_deg    = wnn_result.get("anti_degradation")

        if rna_weight is None or atac_weight is None:
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        "WNN integration completed without measured modality "
                        "weights (reported as null; not estimated)."
                    ),
                    "rna_weight": rna_weight,
                    "atac_weight": atac_weight,
                },
                Confidence.LOW,
            )
            return

        summary = (
            f"WNN integration: RNA weight={rna_weight:.2f}, "
            f"ATAC weight={atac_weight:.2f}"
        )

        if rna_weight > 0.85:
            # Imbalanced — debate required
            try:
                from aria.agents.debate_council import DebateCouncil
                council = DebateCouncil(llm=self.llm, max_rounds=2)

                claim = (
                    f"WNN modality weights are severely imbalanced: "
                    f"RNA={rna_weight:.2f}, ATAC={atac_weight:.2f}. "
                    f"The joint embedding is dominated by RNA signal."
                )
                debate = council.resolve(
                    topic="wnn_modality_weight_imbalance",
                    initial_claim=claim,
                    evidence={
                        "mean_rna_weight":  rna_weight,
                        "mean_atac_weight": atac_weight,
                        "n_cells":          wnn_result.get("n_cells", 0),
                        "atac_frip":        wnn_result.get("atac_frip", None),
                    },
                    biological_context=intent,
                )
                wnn_result["weight_debate"] = debate.consensus
                wnn_result["weight_verdict"] = debate.verdict.value

                conf = Confidence.LOW if rna_weight > 0.90 else Confidence.MEDIUM
            except Exception as e:
                log.warning(f"DebateCouncil unavailable: {e}")
                conf = Confidence.MEDIUM
        else:
            conf = Confidence.HIGH

        payload = {
            "summary": summary,
            "rna_weight": rna_weight,
            "atac_weight": atac_weight,
        }
        if isinstance(anti_deg, dict):
            g_status = str(anti_deg.get("guardrail_status", "pass"))
            g_codes = list(anti_deg.get("issue_codes") or [])
            payload["guardrail_status"] = g_status
            payload["guardrail_issue_codes"] = g_codes
            if g_status == "degraded":
                conf = Confidence.LOW
                payload["summary"] = (
                    f"{summary} (anti-degradation guardrail=degraded: "
                    f"{', '.join(g_codes) if g_codes else 'check per-modality structure'})"
                )
            elif g_status == "warn" and conf == Confidence.HIGH:
                conf = Confidence.MEDIUM

        self.publish_finding(
            experiment_id,
            payload,
            conf,
        )

    # ── Peak-to-gene linking ──────────────────────────────────────────────

    def _run_peak_to_gene(self, experiment_id: str, exp_ctx: dict,
                           intent: dict, mods: dict,
                           rna_results: dict,
                           chromatin_results: dict) -> dict:
        """
        Link ATAC peaks to nearby genes using paired cis correlation + background null.
        Stores strong associative links as Tunnels in ARIAMemory.
        """
        rna_files  = list(mods.get("scRNA", mods.get("bulk_RNA", [])) or [])
        atac_files = list(mods.get("scATAC", mods.get("bulk_ATAC", [])) or [])
        peaks_path = (chromatin_results.get("scATAC", {})
                      .get("findings", {})
                      .get("peaks", {})
                      .get("peaks_path", ""))
        p2g_overrides = dict(exp_ctx.get("peak2gene_params") or {})

        h5mu_path = (
            p2g_overrides.pop("h5mu_path", None)
            or exp_ctx.get("h5mu_path")
            or exp_ctx.get("paired_h5mu")
        )
        if not h5mu_path:
            for cand in [*rna_files, *atac_files]:
                if str(cand).lower().endswith(".h5mu"):
                    h5mu_path = str(cand)
                    break

        ref_paths = [h5mu_path] if h5mu_path else (rna_files or atac_files)
        p2g_params: dict = {
            "rna_files":    rna_files,
            "atac_files":   atac_files,
            "peaks_path":   peaks_path,
            "genome":       exp_ctx.get("genome", "hg38"),
            "organism":     exp_ctx.get("organism", "Homo sapiens"),
            "distance_kb":  int(exp_ctx.get("peak2gene_distance_kb", 500)),
            "min_corr":     float(exp_ctx.get("peak2gene_min_corr", 0.3)),
            "seed":         int(exp_ctx.get("seed", 0)),
            "output_dir":   p2g_overrides.pop(
                "output_dir",
                self._get_output_dir(ref_paths),
            ),
        }
        genome_cfg = exp_ctx.get("genome_config") or {}
        gtf_val = exp_ctx.get("gtf_path") or genome_cfg.get("gtf")
        fasta_val = exp_ctx.get("genome_fasta") or genome_cfg.get("fasta")
        if h5mu_path:
            p2g_params["h5mu_path"] = str(h5mu_path)
        if gtf_val:
            p2g_params["gtf_path"] = str(gtf_val)
        if fasta_val:
            p2g_params["genome_fasta"] = str(fasta_val)
        p2g_params.update(p2g_overrides)
        validated_paired = exp_ctx.get("_validated_paired_h5mu")
        if validated_paired:
            p2g_params["h5mu_path"] = str(validated_paired)
            p2g_params["rna_files"] = [str(validated_paired)]
            p2g_params["atac_files"] = [str(validated_paired)]

        p2g_result = self.env.run_in_stack(
            stack="integration",
            script_path="aria/scripts/integration_peak2gene.py",
            params=p2g_params,
        )

        if p2g_result.get("status") == "success":
            self._store_peak_gene_tunnels(
                experiment_id, p2g_result
            )
            self._debate_peak_gene_links(
                experiment_id, p2g_result, intent
            )

        return p2g_result

    def _store_peak_gene_tunnels(self, experiment_id: str,
                                  p2g_result: dict):
        """
        Store strong associative peak-to-gene links as Tunnels in ARIAMemory.
        Tunnels represent cross-modal connections (ATAC <-> RNA).
        """
        top_links = p2g_result.get("top_links", [])[:50]
        n_stored = 0

        for link in top_links:
            corr_val = float(link.get("correlation", 0.0) or 0.0)
            if abs(corr_val) >= 0.4:
                try:
                    self.memory.create_tunnel(
                        tunnel_id=str(uuid.uuid4())[:8],
                        wing_id=experiment_id,
                        from_hall=f"{experiment_id}_scATAC",
                        to_hall=f"{experiment_id}_scRNA",
                        entity=link.get("gene", "unknown"),
                        description=(
                            f"Peak {link.get('peak', '?')} <-> "
                            f"Gene {link.get('gene', '?')} (associative): "
                            f"r={corr_val:.3f}, "
                            f"dist={float(link.get('distance_kb', 0.0) or 0.0):.1f}kb"
                        ),
                        confidence=(
                            "high" if abs(corr_val) >= 0.6 else "medium"
                        ),
                    )
                    n_stored += 1
                except Exception as e:
                    log.warning(f"Could not store tunnel: {e}")

        log.info(f"Stored {n_stored} peak-gene links as Tunnels in ARIAMemory")

    def _debate_peak_gene_links(self, experiment_id: str,
                                 p2g_result: dict, intent: dict):
        """
        Invoke DebateCouncil for associative peak-gene links.
        Requires independent orthogonal evidence before claiming mechanism.
        """
        n_links   = p2g_result.get("n_links", 0)
        top_links = p2g_result.get("top_links", [])[:5]

        if n_links == 0:
            self.publish_finding(
                experiment_id,
                {"summary": "No significant peak-gene links found",
                 "n_links": 0},
                Confidence.INSUFFICIENT,
            )
            return

        try:
            from aria.agents.debate_council import DebateCouncil
            council = DebateCouncil(llm=self.llm, max_rounds=2)

            top_genes = [l.get("gene", "?") for l in top_links]
            claim = (
                f"{n_links} cis peak-gene accessibility-expression links identified. "
                f"Top linked genes: {', '.join(top_genes)}"
            )

            debate = council.resolve(
                topic="peak_to_gene_regulatory_links",
                initial_claim=claim,
                evidence={
                    "n_links":                 n_links,
                    "n_positive_correlations": p2g_result.get("n_positive_correlations", 0),
                    "n_negative_correlations": p2g_result.get("n_negative_correlations", 0),
                    "top_links":               top_links,
                    "distance_cutoff_kb":      p2g_result.get("distance_kb", 500),
                    "min_correlation":         p2g_result.get("min_corr", 0.3),
                    "max_padj":                p2g_result.get("max_padj"),
                    "background_method":       p2g_result.get("background_method"),
                    "ctcf_validated":          p2g_result.get("ctcf_validated", False),
                    "hic_corroborated":        p2g_result.get("hic_corroborated", False),
                },
                biological_context=intent,
            )

            conf_map = {
                "high":         Confidence.HIGH,
                "medium":       Confidence.MEDIUM,
                "low":          Confidence.LOW,
                "insufficient": Confidence.INSUFFICIENT,
            }
            self.publish_finding(
                experiment_id,
                {"summary":     debate.consensus,
                 "n_links":     n_links,
                 "top_genes":   top_genes,
                 "verdict":     debate.verdict.value,
                 "limitations": debate.limitations},
                conf_map.get(debate.confidence, Confidence.MEDIUM),
            )

        except Exception as e:
            log.warning(f"DebateCouncil unavailable for p2g: {e}")
            conf = (Confidence.HIGH   if n_links > 1000 else
                    Confidence.MEDIUM if n_links > 100  else
                    Confidence.LOW)
            self.publish_finding(
                experiment_id,
                {"summary": f"{n_links} peak-gene links identified",
                 "top_genes": [l.get("gene") for l in top_links]},
                conf,
            )

    # ── Candidate eRegulons (V48 Slice I6 — Tier 2b) ─────────────────────

    @staticmethod
    def _should_run_eregulons(
        is_paired_mudata: bool,
        exp_ctx: dict,
        intent: dict,
    ) -> bool:
        """Determine whether paired candidate eRegulon assembly should run."""
        if exp_ctx.get("enable_eregulons") is False:
            return False
        if exp_ctx.get("enable_eregulons") is True:
            return True
        if exp_ctx.get("eregulon_params"):
            return True
        if is_paired_mudata:
            return True
        question = str(intent.get("user_question") or "").lower()
        return any(
            kw in question
            for kw in ("eregulon", "regulon", "grn", "gene regulatory", "scenic")
        )

    def _run_eregulons(
        self,
        experiment_id: str,
        exp_ctx: dict,
        intent: dict,
        mods: dict,
        p2g_result: dict | None = None,
    ) -> dict:
        """
        Assemble associative candidate eRegulons (`TF -> cis peak -> gene`)
        across shared same-cell RNA + ATAC profiles inside `aria-integration-env`.
        """
        rna_files = list(mods.get("scRNA", []) or [])
        atac_files = list(mods.get("scATAC", []) or [])
        ereg_overrides = dict(exp_ctx.get("eregulon_params") or {})

        h5mu_path = (
            ereg_overrides.pop("h5mu_path", None)
            or exp_ctx.get("h5mu_path")
            or exp_ctx.get("paired_h5mu")
        )
        if not h5mu_path:
            for cand in [*rna_files, *atac_files]:
                if str(cand).lower().endswith(".h5mu"):
                    h5mu_path = str(cand)
                    break

        ref_paths = [h5mu_path] if h5mu_path else (rna_files or atac_files)
        ereg_params: dict = {
            "rna_files": rna_files,
            "atac_files": atac_files,
            "genome": exp_ctx.get("genome", "hg38"),
            "organism": exp_ctx.get("organism", "Homo sapiens"),
            "distance_kb": int(exp_ctx.get("peak2gene_distance_kb", 500)),
            "p2g_min_corr": float(exp_ctx.get("eregulon_p2g_min_corr", 0.05)),
            "p2g_max_padj": float(exp_ctx.get("eregulon_p2g_max_padj", 0.05)),
            "motif_pvalue": float(exp_ctx.get("eregulon_motif_pvalue", 5e-5)),
            "min_tf_gene_corr": float(exp_ctx.get("eregulon_min_tf_gene_corr", 0.05)),
            "min_target_genes": int(exp_ctx.get("eregulon_min_target_genes", 5)),
            "max_regulon_padj": float(exp_ctx.get("eregulon_max_padj", 0.05)),
            "seed": int(exp_ctx.get("seed", 0)),
            "output_dir": ereg_overrides.pop(
                "output_dir",
                self._get_output_dir(ref_paths),
            ),
        }
        genome_cfg = exp_ctx.get("genome_config") or {}
        if h5mu_path:
            ereg_params["h5mu_path"] = str(h5mu_path)
        gtf_val = exp_ctx.get("gtf_path") or genome_cfg.get("gtf")
        fasta_val = exp_ctx.get("genome_fasta") or genome_cfg.get("fasta")
        if gtf_val:
            ereg_params["gtf_path"] = str(gtf_val)
        if fasta_val:
            ereg_params["genome_fasta"] = str(fasta_val)
        for opt_key in ("motif_collection", "motif_file"):
            if exp_ctx.get(opt_key):
                ereg_params[opt_key] = str(exp_ctx[opt_key])
        ereg_params.update(ereg_overrides)
        validated_paired = exp_ctx.get("_validated_paired_h5mu")
        if validated_paired:
            ereg_params["h5mu_path"] = str(validated_paired)
            ereg_params["rna_files"] = [str(validated_paired)]
            ereg_params["atac_files"] = [str(validated_paired)]

        ereg_result = self.env.run_in_stack(
            stack="integration",
            script_path="aria/scripts/integration_eregulons.py",
            params=ereg_params,
        )

        if ereg_result.get("status") == "success":
            self._store_eregulon_tunnels(experiment_id, ereg_result)
            self._debate_eregulons(experiment_id, ereg_result, intent)

        return ereg_result

    def _store_eregulon_tunnels(
        self,
        experiment_id: str,
        ereg_result: dict,
    ) -> None:
        """
        Store top associative candidate eRegulons as cross-modal Tunnels
        in ARIAMemory without assigning causal activation/repression claims.
        """
        top_regs = list(ereg_result.get("top_eregulons") or [])[:30]
        n_stored = 0
        for reg in top_regs:
            if reg.get("depth_confounded"):
                continue
            padj = float(reg.get("regulon_padj", 1.0) or 1.0)
            if padj > 0.05:
                continue
            amb_flag = bool(reg.get("motif_family_ambiguity"))
            try:
                self.memory.create_tunnel(
                    tunnel_id=str(uuid.uuid4())[:8],
                    wing_id=experiment_id,
                    from_hall=f"{experiment_id}_scATAC",
                    to_hall=f"{experiment_id}_scRNA",
                    entity=str(reg.get("tf", "unknown")),
                    description=(
                        f"Candidate eRegulon {reg.get('regulon_id', '?')} "
                        f"({reg.get('sign_class', '?')}, associative): "
                        f"{reg.get('n_enhancer_peaks', 0)} cis peaks -> "
                        f"{reg.get('n_target_genes', 0)} genes, "
                        f"padj={padj:.3g}, motif_ambiguity={amb_flag}"
                    ),
                    confidence=(
                        "low"
                        if amb_flag
                        else ("high" if padj <= 0.005 else "medium")
                    ),
                )
                n_stored += 1
            except Exception as e:
                log.warning(f"Could not store eRegulon tunnel: {e}")

        log.info(f"Stored {n_stored} candidate eRegulons as Tunnels in ARIAMemory")

    def _debate_eregulons(
        self,
        experiment_id: str,
        ereg_result: dict,
        intent: dict,
    ) -> None:
        """Invoke DebateCouncil on associative candidate eRegulon modules."""
        n_ereg = int(ereg_result.get("n_eregulons", 0) or 0)
        top_regs = list(ereg_result.get("top_eregulons") or [])[:5]
        if n_ereg == 0:
            self.publish_finding(
                experiment_id,
                {
                    "summary": "No significant candidate eRegulons met module FDR criteria.",
                    "n_eregulons": 0,
                },
                Confidence.INSUFFICIENT,
            )
            return

        top_tfs = [str(r.get("tf", "?")) for r in top_regs]
        try:
            from aria.agents.debate_council import DebateCouncil
            council = DebateCouncil(llm=self.llm, max_rounds=2)

            claim = (
                f"{n_ereg} associative candidate eRegulon module(s) identified "
                f"across {ereg_result.get('n_unique_tfs', 0)} TF(s). "
                f"Top candidate TFs: {', '.join(top_tfs)}"
            )
            debate = council.resolve(
                topic="candidate_eregulon_modules",
                initial_claim=claim,
                evidence={
                    "n_eregulons": n_ereg,
                    "n_triplets": ereg_result.get("n_triplets", 0),
                    "n_unique_tfs": ereg_result.get("n_unique_tfs", 0),
                    "n_pos_tf_pos_p2g": ereg_result.get("n_pos_tf_pos_p2g", 0),
                    "n_neg_tf_pos_p2g": ereg_result.get("n_neg_tf_pos_p2g", 0),
                    "n_paralog_ambiguous": ereg_result.get("n_paralog_ambiguous", 0),
                    "n_depth_confounded": ereg_result.get("n_depth_confounded", 0),
                    "top_eregulons": top_regs,
                    "associative_only": ereg_result.get("associative_only", True),
                    "perturbation_validated": ereg_result.get(
                        "perturbation_validated", False
                    ),
                },
                biological_context=intent,
            )
            conf_map = {
                "high": Confidence.HIGH,
                "medium": Confidence.MEDIUM,
                "low": Confidence.LOW,
                "insufficient": Confidence.INSUFFICIENT,
            }
            self.publish_finding(
                experiment_id,
                {
                    "summary": debate.consensus,
                    "n_eregulons": n_ereg,
                    "top_tfs": top_tfs,
                    "verdict": debate.verdict.value,
                    "limitations": debate.limitations,
                },
                conf_map.get(debate.confidence, Confidence.MEDIUM),
            )
        except Exception as e:
            log.warning(f"DebateCouncil unavailable for eRegulons: {e}")
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        f"{n_ereg} associative candidate eRegulon module(s) identified"
                    ),
                    "top_tfs": top_tfs,
                },
                Confidence.MEDIUM if n_ereg >= 5 else Confidence.LOW,
            )

    # ── MOFA+ integration ─────────────────────────────────────────────────

    def _run_mofa(self, experiment_id: str, exp_ctx: dict,
                   intent: dict, mods: dict) -> dict:
        """
        MOFA+ factor decomposition for 2+ modalities.
        Identifies shared and modality-specific sources of variation.
        Always uses DebateCouncil for factor interpretation.
        """
        n_factors_decision = self._advise_mofa_factors(
            experiment_id, intent, mods
        )
        if n_factors_decision is None:
            return {
                "status": "skipped",
                "reason": "user_skipped_mofa",
                "details": "User skipped MOFA+ at checkpoint 3.",
            }

        mofa_overrides = dict(exp_ctx.get("mofa_params") or {})
        mofa_mods = {
            k: list(v) for k, v in mods.items()
            if k in ("scRNA", "bulk_RNA", "scATAC", "bulk_ATAC", "HiC")
        }
        first_files = next(iter(mofa_mods.values()), [])
        h5mu_cand = (
            mofa_overrides.pop("h5mu_path", None)
            or exp_ctx.get("h5mu_path")
            or exp_ctx.get("paired_h5mu")
        )
        if not h5mu_cand:
            for cand in [f for files in mofa_mods.values() for f in (files or [])]:
                if str(cand).lower().endswith(".h5mu"):
                    h5mu_cand = str(cand)
                    break
        ref_paths = [str(h5mu_cand)] if h5mu_cand else first_files
        mofa_params: dict = {
            "modalities": mofa_mods,
            "genome": exp_ctx.get("genome", "hg38"),
            "organism": exp_ctx.get("organism", "Homo sapiens"),
            "n_factors": n_factors_decision,
            "seed": int(exp_ctx.get("seed", 0)),
            "max_cells": int(exp_ctx.get("max_cells_mofa", 5000)),
            "output_dir": mofa_overrides.pop(
                "output_dir",
                exp_ctx.get("mofa_output_dir")
                or (self._get_output_dir(ref_paths) if ref_paths else "/tmp/aria_mofa"),
            ),
        }
        if h5mu_cand:
            mofa_params["h5mu_path"] = str(h5mu_cand)
        if isinstance(exp_ctx.get("cell_cycle_genes"), dict):
            mofa_params["cell_cycle_genes"] = exp_ctx["cell_cycle_genes"]
        if exp_ctx.get("batch_key"):
            mofa_params["batch_key"] = str(exp_ctx["batch_key"])
        cond_col = exp_ctx.get("condition_col") or (
            (exp_ctx.get("design") or {}).get("condition_col")
        )
        if cond_col:
            mofa_params["condition_col"] = str(cond_col)
        mofa_params.update(mofa_overrides)
        validated_paired = exp_ctx.get("_validated_paired_h5mu")
        if validated_paired:
            mofa_params["h5mu_path"] = str(validated_paired)
            mofa_params["modalities"] = {
                **dict(mofa_params.get("modalities") or {}),
                "scRNA": [str(validated_paired)],
                "scATAC": [str(validated_paired)],
            }

        mofa_result = self.env.run_in_stack(
            stack="integration",
            script_path="aria/scripts/integration_mofa.py",
            params=mofa_params,
        )

        if mofa_result.get("status") == "success":
            self._interpret_mofa_factors(
                experiment_id, mofa_result, intent
            )

        return mofa_result

    def _advise_mofa_factors(self, experiment_id: str,
                              intent: dict, mods: dict) -> int | None:
        """
        Advise on number of MOFA+ factors using biological context.
        Stores decision in ARIAMemory for reproducibility.
        """
        question    = intent.get("user_question", "").lower()
        n_modalities = len([m for m in mods if m != "unknown"])

        COMPLEX_KEYWORDS = ["atlas", "trajectory", "heterogeneous",
                             "many cell types", "complex"]
        complexity = "complex" if any(k in question for k in COMPLEX_KEYWORDS) \
                     else "moderate"

        # Layer 1: intent-based range
        if complexity == "complex":
            n_factors = 20
        elif n_modalities >= 3:
            n_factors = 15
        else:
            n_factors = 10

        question = (
            "MOFA+ factor count\n\n"
            f"Recommended n_factors={n_factors} based on {n_modalities} "
            f"modalities and {complexity} experimental complexity.\n\n"
            "Approve the factor count to proceed with MOFA+."
        )
        _, user_decision = self.publish_blocking_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=question,
            options=[
                f"Use recommended n_factors={n_factors}",
                "Other: Enter custom n_factors",
                "Skip MOFA+",
            ],
            context={
                "analysis_type": "mofa_integration",
                "parameter_name": "n_factors",
                "recommended_value": n_factors,
            },
        )
        choice = self.checkpoint_choice_text(user_decision)
        if not choice:
            return None
        if "skip" in choice.lower():
            return None
        override = self.numeric_checkpoint_override(user_decision, int)
        chosen_factors = override if override is not None else n_factors

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=3,
            question="MOFA+ number of factors",
            decision=str(chosen_factors),
            rationale=(
                f"n_factors={chosen_factors} based on user approval of the "
                f"checkpoint recommendation for {n_modalities} modalities and "
                f"{complexity} experimental complexity."
            ),
            made_by="user" if override is not None else "advisor",
        )

        return chosen_factors

    def _interpret_mofa_factors(self, experiment_id: str,
                                 mofa_result: dict, intent: dict):
        """
        Interpret MOFA+ factors using DebateCouncil.
        Always check if Factor 1 captures technical variation.
        """
        try:
            from aria.agents.debate_council import DebateCouncil
            council = DebateCouncil(llm=self.llm, max_rounds=2)

            factors        = mofa_result.get("top_factors", [])
            variance_exp   = mofa_result.get("variance_explained", {})
            factor1_top    = mofa_result.get("factor1_top_features", [])
            shared_factors = mofa_result.get("shared_factors", [])
            excluded_list  = mofa_result.get("excluded_factors", [])
            excluded_ids   = [
                f.get("factor_id") for f in excluded_list
                if isinstance(f, dict) and f.get("factor_id") is not None
            ]
            if not excluded_ids and mofa_result.get("technical_factor_flag"):
                tfid = mofa_result.get("technical_factor_id")
                if tfid is not None:
                    excluded_ids.append(tfid)

            claim = (
                f"MOFA+ identified {len(factors)} factors "
                f"({len(shared_factors)} shared cross-modal, "
                f"{len(excluded_ids)} excluded confounder). "
                f"Variance explained: {list(variance_exp.values())[:3]}. "
                f"Factor 1 top features: {factor1_top[:5]}"
            )

            debate = council.resolve(
                topic="mofa_factor_interpretation",
                initial_claim=claim,
                evidence={
                    "n_factors":        len(factors),
                    "n_shared_factors": len(shared_factors),
                    "variance_explained": variance_exp,
                    "factor1_features": factor1_top,
                    "technical_factor_flag": mofa_result.get(
                        "technical_factor_flag", False
                    ),
                    "technical_factor_id": mofa_result.get("technical_factor_id"),
                    "cell_cycle_factor_flag": mofa_result.get(
                        "cell_cycle_factor_flag", False
                    ),
                    "excluded_factor_ids": excluded_ids,
                },
                biological_context=intent,
            )

            conf_map = {
                "high": Confidence.HIGH, "medium": Confidence.MEDIUM,
                "low": Confidence.LOW,  "insufficient": Confidence.INSUFFICIENT,
            }
            self.publish_finding(
                experiment_id,
                {"summary":             debate.consensus,
                 "n_factors":           len(factors),
                 "n_shared_factors":    len(shared_factors),
                 "excluded_factor_ids": excluded_ids,
                 "verdict":             debate.verdict.value,
                 "limitations":         debate.limitations},
                conf_map.get(debate.confidence, Confidence.MEDIUM),
            )

        except Exception as e:
            log.warning(f"MOFA factor debate failed: {e}")
            self.publish_finding(
                experiment_id,
                {"summary": f"MOFA+ factors: {mofa_result.get('n_factors', 0)}"},
                Confidence.MEDIUM,
            )

    # ── Cross-modal conflict resolution ───────────────────────────────────

    def _resolve_conflicts(self, experiment_id: str, exp_ctx: dict,
                            intent: dict, results: dict) -> dict:
        """
        Identify and reason about cross-modal discordances.

        Cases are descriptive discordances only; ARIA stores them with explicit
        confidence and limitations, without assigning a mechanism from pattern
        alone.
        """
        conflicts = []

        # Check WNN vs RNA cluster discordance
        wnn    = results.get("wnn", {})
        p2g    = results.get("peak_to_gene", {})
        ereg   = results.get("eregulons", {})
        mofa   = results.get("mofa", {})

        if wnn.get("status") == "success":
            n_discordant = wnn.get("n_discordant_clusters", 0)
            if n_discordant > 0:
                conflict = {
                    "type":        "wnn_rna_cluster_discordance",
                    "description": (
                        f"{n_discordant} cell clusters differ between "
                        f"RNA-only and WNN joint embedding. "
                        f"Chromatin accessibility provides additional "
                        f"cell identity information beyond transcriptomics."
                    ),
                    "confidence":  "medium",
                    "action":      "Review discordant clusters manually",
                }
                conflicts.append(conflict)

            anti_deg = wnn.get("anti_degradation")
            if isinstance(anti_deg, dict):
                issue_codes = set(anti_deg.get("issue_codes") or [])
                if (
                    "rare_population_degraded" in issue_codes
                    or "reference_rare_population_degraded" in issue_codes
                ):
                    conflicts.append(
                        {
                            "type": "wnn_rare_population_degradation",
                            "description": (
                                "One or more rare per-modality populations lost "
                                "local neighborhood cohesion in the joint WNN "
                                "embedding. Retain per-modality clustering as "
                                "primary reference ('lens, not replacement')."
                            ),
                            "confidence": "low",
                            "action": "Inspect rare per-modality clusters directly",
                        }
                    )
                if "spurious_bridge_candidate" in issue_codes:
                    flagged_cl = (
                        anti_deg.get("spurious_bridge_screen", {}).get(
                            "flagged_clusters", []
                        )
                    )
                    conflicts.append(
                        {
                            "type": "wnn_spurious_bridge_candidate",
                            "description": (
                                f"Joint cluster(s) {flagged_cl} mix distinct "
                                f"single-modality populations without native kNN "
                                f"support. Do not interpret as transitional states."
                            ),
                            "confidence": "low",
                            "action": "Screen flagged clusters for doublets or graph artifacts",
                        }
                    )
                if issue_codes & {
                    "rna_structure_degraded",
                    "atac_structure_degraded",
                    "held_out_joint_degraded",
                }:
                    conflicts.append(
                        {
                            "type": "wnn_modality_structure_degradation",
                            "description": (
                                "Joint WNN graph degraded informative per-modality "
                                "neighborhood or held-out latent predictability."
                            ),
                            "confidence": "low",
                            "action": "Use single-modality embeddings as primary reference",
                        }
                    )

            batch_int = wnn.get("batch_integration")
            if isinstance(batch_int, dict):
                qa = batch_int.get("quality_assessment") or {}
                if qa.get("possible_overcorrection"):
                    conflicts.append(
                        {
                            "type": "wnn_harmony_possible_overcorrection",
                            "description": (
                                "Batch/donor Harmony harmonization aggressively reduced "
                                "batch separation while joint WNN cluster cohesion was "
                                "weak; biological subpopulation variation may be compressed."
                            ),
                            "confidence": "low",
                            "action": "Compare harmonized WNN clusters against unintegrated per-modality clusters",
                        }
                    )

        if mofa.get("status") == "success" and mofa.get("batch_factor_flag"):
            b_fids = mofa.get("batch_factor_ids") or []
            conflicts.append(
                {
                    "type": "mofa_batch_confounded_factors",
                    "description": (
                        f"MOFA+ factor(s) {b_fids} strongly associate with "
                        f"batch/donor '{mofa.get('batch_key')}' and were flagged "
                        "for exclusion from biological interpretation."
                    ),
                    "confidence": "medium",
                    "action": "Restrict downstream biological claims to non-confounded factors",
                }
            )

        # Check peak-gene negative correlations (discordant evidence)
        if p2g.get("status") == "success":
            n_negative = p2g.get("n_negative_correlations", 0)
            if n_negative > 0:
                conflict = {
                    "type":        "accessibility_expression_discordance",
                    "description": (
                        f"{n_negative} peak-gene pair(s) show negative "
                        f"accessibility-expression correlation. "
                        f"Treat these as discordant evidence until an "
                        f"independent validation assay supports mechanism."
                    ),
                    "n_events":    n_negative,
                    "confidence":  "medium",
                    "action":      "Prioritize for independent validation",
                }
                conflicts.append(conflict)
                # Store as findings for NarrativeAgent
                self.publish_finding(
                    experiment_id,
                    {"summary": conflict["description"],
                     "n_events": n_negative,
                     "type": "accessibility_expression_discordance"},
                    Confidence.MEDIUM,
                )

        # Check candidate eRegulon discordances, paralog ambiguity & depth confounding
        if ereg.get("status") == "success":
            n_neg_tf = int(ereg.get("n_neg_tf_pos_p2g", 0) or 0)
            if n_neg_tf > 0:
                conflicts.append(
                    {
                        "type": "eregulon_negative_tf_gene_discordance",
                        "description": (
                            f"{n_neg_tf} candidate eRegulon(s) show inverse TF-to-target "
                            "expression correlation across positively linked cis peaks "
                            "(neg_tf_pos_p2g). Treat as associative discordance rather "
                            "than direct transcriptional repression without orthogonal "
                            "perturbation or chromatin-binding assays."
                        ),
                        "n_events": n_neg_tf,
                        "confidence": "medium",
                        "action": "Validate inverse TF-target candidates by perturbation or CUT&RUN",
                    }
                )
            n_ambig = int(ereg.get("n_paralog_ambiguous", 0) or 0)
            if n_ambig > 0:
                conflicts.append(
                    {
                        "type": "eregulon_paralog_motif_ambiguity",
                        "description": (
                            f"{n_ambig} candidate eRegulon(s) share substantially "
                            "overlapping PWM enhancer cistromes with co-expressed "
                            "paralogous TFs; motif enrichment alone cannot disambiguate "
                            "which family member binds in vivo."
                        ),
                        "n_events": n_ambig,
                        "confidence": "low",
                        "action": "Disambiguate paralogs via TF ChIP/CUT&RUN or Perturb-seq",
                    }
                )
            n_depth = int(ereg.get("n_depth_confounded", 0) or 0)
            if n_depth > 0:
                conflicts.append(
                    {
                        "type": "eregulon_depth_confounded",
                        "description": (
                            f"{n_depth} candidate eRegulon(s) lose substantial "
                            "TF-target co-variation after conditioning on RNA/ATAC "
                            "library depth and are flagged as depth-confounded."
                        ),
                        "n_events": n_depth,
                        "confidence": "low",
                        "action": "Exclude depth-confounded modules from biological claims",
                    }
                )

        # Check MOFA+ confounder-associated factors (depth / phase covariates)
        if mofa.get("status") == "success":
            excl_ids = [
                f.get("factor_id")
                for f in (mofa.get("excluded_factors") or [])
                if isinstance(f, dict) and f.get("factor_id") is not None
            ]
            if not excl_ids and mofa.get("technical_factor_flag"):
                tfid = mofa.get("technical_factor_id")
                if tfid is not None:
                    excl_ids.append(tfid)
            if not excl_ids and mofa.get("cell_cycle_factor_flag"):
                excl_ids.extend(mofa.get("cell_cycle_factor_ids") or [])
            if excl_ids:
                conflicts.append(
                    {
                        "type": "mofa_confounder_factor_excluded",
                        "description": (
                            f"MOFA+ factor(s) {excl_ids} strongly correlate with "
                            "measured technical sequencing depth or phase covariates "
                            "and are excluded from biological lineage interpretation."
                        ),
                        "excluded_factor_ids": excl_ids,
                        "confidence": "medium",
                        "action": "Exclude flagged confounder factors from biological claims",
                    }
                )

        if conflicts:
            log.info(f"Identified {len(conflicts)} cross-modal conflicts")

        return {
            "n_conflicts": len(conflicts),
            "conflicts":   conflicts,
        }

    # ── Checkpoint formatters ─────────────────────────────────────────────

    def _format_wnn_checkpoint(self, k_decision) -> str:
        lines = ["WNN Integration — k Parameter\n"]
        lines.append("Candidates evaluated:\n")
        for c in k_decision.candidates:
            marker = " [RECOMMENDED]" if c.recommended else ""
            if c.metrics.get("measured") is False:
                basis = c.metrics.get("basis", "prior")
                lines.append(
                    f"  k={c.value:2d}  basis={basis}  "
                    f"pre_run_metrics=not_computed{marker}"
                )
            else:
                rna_w = c.metrics.get("rna_weight")
                atac_w = c.metrics.get("atac_weight")
                if rna_w is not None and atac_w is not None:
                    lines.append(
                        f"  k={c.value:2d}  RNA_weight={rna_w:.2f}  "
                        f"ATAC_weight={atac_w:.2f}{marker}"
                    )
                else:
                    lines.append(f"  k={c.value:2d}{marker}")
        lines.append(f"\n  {k_decision.justification}")
        lines.append(
            "\n  Note: RNA/ATAC modality weights are reported only after the "
            "WNN script executes; this checkpoint does not display estimated "
            "weights as measured metrics."
        )
        lines.append("\nApprove k to proceed with WNN integration:")
        return "\n".join(lines)

    # ── Helpers ───────────────────────────────────────────────────────────

    def _get_output_dir(self, files: list) -> str:
        if files:
            parent = Path(files[0]).parent
            if parent.name == "aria_integration":
                return str(parent)
            return str(parent / "aria_integration")
        return "/tmp/aria_integration"

    def receive(self, message):
        pass
