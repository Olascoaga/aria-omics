"""
ARIA GenomeArchAgent (ADR-082 — Production 3D Genome Architecture Engine)
-------------------------------------------------------------------------
3D genome organization analysis: Hi-C / Micro-C preprocessing, QC, ICE
balancing, P(s) polymer decay, HiCRep SCC reproducibility, GC/GTF-phased A/B
compartments, Crane insulation TADs, expected_cis-calibrated chromatin loops,
and replicate-gated differential 3D architecture (ΔE1 & ΔInsulation).

Supported formats:
  .cool / .mcool        — Cooler format (single or multi-resolution)
  .hic                  — Juicer / 4D Nucleome format (converted out-of-core)
  .pairs / .pairs.gz    — 4DN contact pairs (cload + zoomify)
  .bam / .sam           — Chimeric Hi-C alignments (pairtools parse/sort/dedup)
  .fastq / .fastq.gz    — Raw paired Hi-C reads (bwa-mem2 -SP5M -> pairtools -> cooler)

Memory strategy (critical for Hi-C):
  - NEVER load the full genome-wide matrix into RAM
  - All analysis uses out-of-core operations (chunked pixels & chromosome views)
"""

from __future__ import annotations

import logging
import uuid
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.agents.chromatin_helpers import is_fastq as _is_fastq
from aria.bus.message_bus import Confidence
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.llm.provider import LLMProvider
from aria.memory.memory import ARIAMemory
from aria.utils.hic_constants import RAM_ESTIMATES_GB

log = logging.getLogger("aria.hic")


GENOME_ARCH_SYSTEM = """
You are ARIA's GenomeArchAgent — a specialist in 3D genome organization.

Your expertise:
- Hi-C / Micro-C data processing: ICE/KR balancing, P(s) polymer decay, HiCRep SCC
- Compartments A/B: eigenvector decomposition phased by GC or gene density
- TAD calling: Crane Insulation Score (primary), Li prominence boundary strength
- Chromatin loops: expected_cis-calibrated donut dot calling, promoter-enhancer anchors
- Multi-replicate differential 3D organization: moderated linear ΔE1 and ΔInsulation

Critical knowledge:

COMPARTMENTS:
  - Compartment A = active chromatin (high GC / gene density, euchromatin)
  - Compartment B = inactive / heterochromatic chromatin
  - Compartment B is NOT always transcriptionally silent. Never claim "B = silent".
  - PC1 sign is arbitrary unless phased against GC content or gene density; always
    check pc1_validated before interpreting A vs B polarity.

TADs:
  - Insulation Score (Crane 2015) with window_size calibrated on chr1.
  - ~30% of TAD boundaries vary across cell types; never claim cell-type specificity
    without replicate-backed differential insulation testing.

LOOPS:
  - Physical proximity (3D contact enrichment) is associative and does NOT by itself
    prove direct transcriptional activation without perturbation evidence.
""".strip()


GENOME_SCALE = {
    "Homo sapiens": 1.0,
    "Mus musculus": 0.9,
    "Drosophila melanogaster": 0.02,
    "C. elegans": 0.003,
    "Danio rerio": 0.5,
    "S. cerevisiae": 0.001,
}


class GenomeArchAgent(BaseAgent):

    name = "genome_arch_agent"
    description = (
        "3D genome organization: Hi-C/Micro-C preprocessing, QC, ICE balancing, "
        "compartments A/B, TAD calling, loop calling, and differential 3D topology."
    )
    validation_level = "production"  # ADR-082
    dispatch_enabled = True
    REQUIRED_SCRIPTS = (
        "aria/scripts/hic_preprocess.py",
        "aria/scripts/hic_inspect.py",
        "aria/scripts/hic_qc_and_balance.py",
        "aria/scripts/hic_topology.py",
        "aria/scripts/hic_figure.py",
    )
    PLANNED_SCRIPTS = ()

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.advisor = ParameterAdvisor(memory, llm)

        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    # ── Main entry point ──────────────────────────────────────────────────

    def _resolve_reference_assets(self, exp_ctx: dict) -> tuple[Optional[str], Optional[str]]:
        genome_fasta = exp_ctx.get("genome_fasta") or exp_ctx.get("reference_fasta")
        if not genome_fasta:
            try:
                from aria.utils import genomes
                local, _src = genomes.resolve_local_genome_fasta(
                    exp_ctx.get("genome", "hg38")
                )
                genome_fasta = local
            except Exception:
                genome_fasta = None
        gtf_path = (
            exp_ctx.get("gtf")
            or exp_ctx.get("gtf_path")
            or exp_ctx.get("gene_annotation_gtf")
        )
        return (
            str(genome_fasta) if genome_fasta else None,
            str(gtf_path) if gtf_path else None,
        )

    def _align_hic_fastqs(self, experiment_id: str, exp_ctx: dict,
                          fastq_files: list, genome_fasta: Optional[str],
                          out_dir: Path) -> tuple[Optional[list[str]], dict]:
        """Align raw paired Hi-C FASTQs with ``bwa-mem2 mem -SP5M`` (rnaseq+atacseq
        stacks) to produce chimeric unsorted BAMs for ``pairtools parse``."""
        assembly = exp_ctx.get("genome", "hg38")
        if not genome_fasta:
            return None, {
                "status": "skipped",
                "ran": False,
                "analysis": "hic_fastq_alignment",
                "data_type": "HiC",
                "validation_level": "production",
                "reason": "no_reference_genome_fasta",
                "message": (
                    f"Raw Hi-C FASTQ alignment requires a local reference FASTA "
                    f"for assembly '{assembly}'."
                ),
            }

        self.publish_status(experiment_id, "Hi-C FASTQ QC (fastp)...", 0.05)
        qc_params = {
            "fastq_dir": str(Path(fastq_files[0]).parent),
            "fastq_files": [str(f) for f in fastq_files],
            "output_dir": str(out_dir / "fastq_qc"),
            "threads": int(exp_ctx.get("threads", 8)),
        }
        for key in ("bulk_fastq_manifest", "sample_ids"):
            if exp_ctx.get(key):
                qc_params[key] = exp_ctx[key]
        fq_qc = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_fastq_qc.py",
            params=qc_params,
        )
        if fq_qc.get("status") != "success" or not fq_qc.get("samples"):
            return None, {
                "status": "failed",
                "ran": False,
                "analysis": "hic_fastq_alignment",
                "data_type": "HiC",
                "stage": "fastq_qc",
                "details": fq_qc,
            }

        self.publish_status(experiment_id, "Hi-C chimeric alignment (bwa-mem2 -SP5M)...", 0.10)
        align = self.env.run_in_stack(
            stack="atacseq",
            script_path="aria/scripts/atac_align.py",
            params={
                "samples": fq_qc.get("samples", []),
                "genome_fasta": genome_fasta,
                "bwa_index": exp_ctx.get("bwa_index"),
                "output_dir": str(out_dir / "aligned_chimeric"),
                "threads": int(exp_ctx.get("threads", 8)),
                "assay_type": "HiC",
                "hic_chimeric": True,
                "remove_duplicates": False,
                "mapq": 0,
                "remove_mito": False,
            },
        )
        bams = [
            b.get("bam") for b in (align.get("bam_files") or [])
            if isinstance(b, dict) and b.get("status") == "success" and b.get("bam")
        ]
        if align.get("status") != "success" or not bams:
            return None, {
                "status": "failed",
                "ran": False,
                "analysis": "hic_fastq_alignment",
                "data_type": "HiC",
                "stage": "align",
                "details": align,
            }
        return bams, {
            "status": "success",
            "ran": True,
            "analysis": "hic_fastq_alignment",
            "data_type": "HiC",
            "validation_level": "production",
            "bam_files": align.get("bam_files"),
        }

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx = context.get("exp_context", {})
        intent = context.get("biological_intent", {})
        mods = exp_ctx.get("modalities", {})
        hic_files = [str(f) for f in (mods.get("HiC") or []) if f]

        if not hic_files:
            self.publish_finding(
                experiment_id,
                {"error": "No Hi-C/Micro-C files detected"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "no_hic_data"}

        self.publish_status(experiment_id, "GenomeArchAgent starting...", 0.0)
        work_root = Path(
            str(exp_ctx.get("output_dir") or (Path(hic_files[0]).parent / "aria_hic"))
        )
        work_root.mkdir(parents=True, exist_ok=True)

        genome_fasta, gtf_path = self._resolve_reference_assets(exp_ctx)
        findings: dict = {}
        active_files = list(hic_files)

        # Phase 0a: Raw FASTQ chimeric alignment if raw reads were supplied
        if _is_fastq(active_files):
            aligned_bams, fq_prep = self._align_hic_fastqs(
                experiment_id, exp_ctx, active_files, genome_fasta, work_root / "fastq_align"
            )
            findings["fastq_alignment"] = fq_prep
            if not aligned_bams:
                return {
                    "status": "failed",
                    "findings": findings,
                    "reason": fq_prep.get("reason", "hic_fastq_alignment_failed"),
                }
            active_files = aligned_bams

        # Phase 0b: Preprocess .pairs / .hic / .bam (or single-res .cool needing zoomify)
        needs_preprocess = any(
            Path(f).name.lower().endswith((".pairs", ".pairs.gz", ".hic", ".bam", ".sam"))
            for f in active_files
        ) or bool(exp_ctx.get("zoomify_cool", False))

        if needs_preprocess:
            self.publish_status(experiment_id, "Preprocessing Hi-C contact files...", 0.12)
            prep_params: dict = {
                "files": active_files,
                "genome": exp_ctx.get("genome", "hg38"),
                "resolution": int(exp_ctx.get("base_resolution", 10_000)),
                "zoomify": bool(exp_ctx.get("zoomify", True)),
                "output_dir": str(work_root / "preprocessed"),
                "threads": int(exp_ctx.get("threads", 4)),
            }
            for opt_key, val in (
                ("sample_ids", exp_ctx.get("sample_ids")),
                ("chrom_sizes", exp_ctx.get("chrom_sizes")),
                ("genome_fasta", genome_fasta),
                ("genome_fai", exp_ctx.get("genome_fai")),
                ("zoom_resolutions", exp_ctx.get("zoom_resolutions")),
            ):
                if val is not None:
                    prep_params[opt_key] = val
            prep_res = self.env.run_in_stack(
                stack="hic",
                script_path="aria/scripts/hic_preprocess.py",
                params=prep_params,
            )
            findings["preprocessing"] = prep_res
            if prep_res.get("status") != "success" or not prep_res.get("cooler_files"):
                return {
                    "status": "failed",
                    "findings": findings,
                    "reason": prep_res.get("error_type", "hic_preprocess_failed"),
                }
            active_files = list(prep_res["cooler_files"])
            if prep_res.get("sample_ids") and not exp_ctx.get("sample_ids"):
                exp_ctx = {**exp_ctx, "sample_ids": prep_res["sample_ids"]}

        # Phase 1: Inspect files and available resolutions
        file_info = self._inspect_hic_files(active_files, exp_ctx)
        findings["file_info"] = file_info

        # Checkpoint 2: Resolution selection
        resolution_decision = self._advise_resolution(
            experiment_id, intent, file_info, exp_ctx
        )
        findings["resolution_decision"] = resolution_decision

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=2,
            question=self._format_resolution_checkpoint(
                file_info, resolution_decision
            ),
            options=[
                f"Use recommended resolution ({resolution_decision['recommended_resolution']:,}bp)",
                "Select different resolution",
                "Skip 3D genome analysis",
            ],
            context={
                "file_info": file_info,
                "resolution_decision": resolution_decision,
            },
        )

        # Phase 2: QC, ICE balancing, P(s) contact decay, and HiCRep SCC
        self.publish_status(experiment_id, "Hi-C QC and ICE balancing...", 0.25)
        qc_params: dict = {
            "files": active_files,
            "genome": exp_ctx.get("genome", "hg38"),
            "resolution": resolution_decision["recommended_resolution"],
            "output_dir": str(work_root / "qc_balance"),
        }
        for opt_key in ("sample_ids", "min_cis_trans_ratio", "min_valid_pairs", "force_rebalance"):
            if exp_ctx.get(opt_key) is not None:
                qc_params[opt_key] = exp_ctx[opt_key]
        qc_result = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_qc_and_balance.py",
            params=qc_params,
        )
        findings["qc"] = qc_result

        if qc_result.get("status") == "error":
            log.warning("Hi-C QC failed: %s", qc_result.get("details", ""))

        if qc_result.get("warnings"):
            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint=3,
                question=self._format_qc_checkpoint(qc_result),
                options=[
                    "Continue — warnings noted in report",
                    "Abort — data quality insufficient",
                ],
                context={"qc": qc_result},
            )

        topology_files = (
            qc_result.get("balanced_files")
            if isinstance(qc_result, dict) and qc_result.get("balanced_files")
            else active_files
        )

        # Phase 3: Topology (compartments, TADs, loops, multi-replicate differential 3D)
        self.publish_status(experiment_id, "Computing 3D topology...", 0.55)
        topology_result = self._run_topology(
            experiment_id, exp_ctx, intent,
            topology_files, resolution_decision, qc_result,
            genome_fasta=genome_fasta, gtf_path=gtf_path,
            work_root=work_root,
        )
        findings["topology"] = topology_result

        # Phase 4: Dual PNG + SVG publication figures
        self.publish_status(experiment_id, "Rendering 3D genome figures...", 0.90)
        fig_params: dict = {
            "output_dir": str(work_root / "figures"),
            "cooler_path": topology_files[0] if topology_files else None,
            "resolution": resolution_decision["recommended_resolution"],
            "ps_curve_tsv": qc_result.get("ps_curve_tsv") if isinstance(qc_result, dict) else None,
            "per_file": qc_result.get("per_file") if isinstance(qc_result, dict) else None,
            "compartments_tsv": (topology_result.get("compartments") or {}).get("compartments_tsv"),
            "insulation_tsv": (topology_result.get("tads") or {}).get("insulation_tsv"),
            "tads_bed": (topology_result.get("tads") or {}).get("tads_bed"),
            "differential_compartments_csv": (
                (topology_result.get("differential_3d") or {}).get("differential_compartments_csv")
            ),
        }
        fig_res = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_figure.py",
            params={k: v for k, v in fig_params.items() if v is not None},
        )
        if isinstance(fig_res, dict) and fig_res.get("status") == "success":
            findings["figures"] = dict(fig_res.get("figures") or {})
            findings["figure_bundles"] = dict(fig_res.get("figure_bundles") or {})

        self.publish_status(experiment_id, "GenomeArchAgent complete.", 1.0)
        return {"status": "done", "findings": findings}

    # ── Phase 1: File inspection ──────────────────────────────────────────

    def _inspect_hic_files(self, files: list, exp_ctx: dict) -> dict:
        result = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_inspect.py",
            params={
                "files": files,
                "genome": exp_ctx.get("genome", "hg38"),
            },
        )
        if result.get("status") == "error":
            log.warning("File inspection failed — using filename inference")
            result = self._infer_file_info(files)
        return result

    def _infer_file_info(self, files: list) -> dict:
        formats = []
        for f in files:
            s = str(f)
            if s.endswith(".cool") or s.endswith(".mcool"):
                formats.append("cooler")
            elif s.endswith(".hic"):
                formats.append("hic")
            elif s.endswith(".pairs") or s.endswith(".pairs.gz"):
                formats.append("pairs")
        return {
            "status": "partial",
            "files": files,
            "formats": list(set(formats)),
            "available_resolutions": [1_000_000, 500_000, 100_000, 40_000, 10_000, 5_000],
            "note": "Full inspection requires aria-hic-env",
        }

    # ── Resolution ParameterAdvisor ───────────────────────────────────────

    def _advise_resolution(self, experiment_id: str, intent: dict,
                           file_info: dict, exp_ctx: dict) -> dict:
        question = intent.get("user_question", "").lower()
        organism = exp_ctx.get("organism", "Homo sapiens")
        scale = GENOME_SCALE.get(organism, 1.0)
        available = file_info.get("available_resolutions") or [
            1_000_000, 100_000, 40_000, 10_000
        ]

        explicit_res = exp_ctx.get("resolution")
        LOOP_KEYWORDS = ["loop", "enhancer", "promoter", "ctcf anchor",
                         "hi-c loop", "chromatin loop", "dot"]
        TAD_KEYWORDS = ["tad", "domain", "topological", "boundary",
                        "insulation", "compartment switch"]
        COMP_KEYWORDS = ["compartment", "a/b", "heterochromatin",
                         "euchromatin", "a-b", "lamina"]

        needs_loops = any(k in question for k in LOOP_KEYWORDS)
        needs_tads = any(k in question for k in TAD_KEYWORDS)
        needs_compartments = any(k in question for k in COMP_KEYWORDS)

        if explicit_res is not None:
            target_res = int(explicit_res)
            analysis = (
                "loop_calling" if target_res <= 10_000
                else "tad_calling" if target_res <= 50_000
                else "compartments"
            )
        elif needs_loops:
            target_res = 10_000 if 5_000 not in available and 10_000 in available else 5_000
            analysis = "loop_calling"
        elif needs_tads:
            target_res = 40_000
            analysis = "tad_calling"
        elif needs_compartments:
            target_res = 100_000
            analysis = "compartments"
        else:
            target_res = 40_000
            analysis = "tad_calling"

        if target_res not in available:
            coarser = [r for r in sorted(available) if r >= target_res]
            target_res = coarser[0] if coarser else sorted(available)[0]

        base_ram = RAM_ESTIMATES_GB.get(target_res, 1.0)
        req_ram = base_ram * scale

        sorted_avail = sorted(available)
        idx = sorted_avail.index(target_res) if target_res in sorted_avail else 0

        candidates = []
        seen_idx: set[int] = set()
        for i in [max(0, idx - 1), idx, min(len(sorted_avail) - 1, idx + 1)]:
            if i in seen_idx:
                continue
            seen_idx.add(i)
            res = sorted_avail[i]
            ram_need = RAM_ESTIMATES_GB.get(res, 1.0) * scale
            candidates.append({
                "resolution": res,
                "ram_required_gb": round(ram_need, 1),
                "analysis_depth": self._resolution_to_depth(res),
                "recommended": (res == target_res),
            })

        decision = {
            "recommended_resolution": target_res,
            "analysis_type": analysis,
            "ram_required_gb": round(req_ram, 1),
            "available_resolutions": sorted_avail,
            "candidates": candidates,
            "rationale": (
                f"Resolution {target_res:,}bp selected for {analysis}. "
                f"Requires ~{req_ram:.1f}GB RAM for {organism}. "
                f"Intent: '{question[:80]}'"
            ),
        }

        if self.advisor is not None:
            try:
                hist = self.advisor._recall_similar_decisions(
                    experiment_id, "hic_resolution", intent
                )
                if hist:
                    decision["historical_precedent"] = hist[:2]
            except Exception:
                pass

        try:
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=2,
                question=f"Hi-C resolution for {analysis}",
                decision=str(target_res),
                rationale=decision["rationale"],
                made_by="advisor",
            )
        except Exception as e:
            log.warning("Could not store resolution decision: %s", e)

        return decision

    # ── Phase 3: Topology ─────────────────────────────────────────────────

    def _run_topology(self, experiment_id: str, exp_ctx: dict,
                      intent: dict, files: list,
                      resolution_decision: dict,
                      qc_result: dict,
                      genome_fasta: Optional[str] = None,
                      gtf_path: Optional[str] = None,
                      work_root: Optional[Path] = None) -> dict:
        resolution = int(resolution_decision["recommended_resolution"])
        analysis = resolution_decision["analysis_type"]
        organism = exp_ctx.get("organism", "Homo sapiens")
        avail_res = resolution_decision.get("available_resolutions") or [resolution]
        topo_root = (work_root / "topology") if work_root else (Path(files[0]).parent / "aria_hic_topology")

        # 1. Compartments A/B (prefer explicit comp_resolution, else >= 100kb when available)
        if exp_ctx.get("comp_resolution") is not None:
            comp_resolution = int(exp_ctx["comp_resolution"])
        else:
            comp_candidates = [int(r) for r in avail_res if 100_000 <= int(r) <= 250_000]
            if comp_candidates:
                comp_resolution = min(comp_candidates)
            else:
                ge_100 = [int(r) for r in avail_res if int(r) >= 100_000]
                comp_resolution = min(ge_100) if ge_100 else max(int(r) for r in avail_res)

        comp_params: dict = {
            "files": files,
            "genome": exp_ctx.get("genome", "hg38"),
            "organism": organism,
            "analysis": "compartments",
            "resolution": comp_resolution,
            "chromosomes": exp_ctx.get("chromosomes", "all"),
            "out_of_core": True,
            "output_dir": str(topo_root / "compartments"),
        }
        if genome_fasta:
            comp_params["genome_fasta"] = genome_fasta
        if gtf_path:
            comp_params["gtf"] = gtf_path

        compartments = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_topology.py",
            params=comp_params,
        )

        # 2. TAD calling via Crane Insulation Score
        tads: dict = {}
        window_size = exp_ctx.get("window_size")
        if resolution <= 100_000:
            if window_size is None:
                window_decision = self._advise_insulation_window(
                    experiment_id, intent, files,
                    resolution, exp_ctx.get("genome", "hg38"),
                    out_dir=topo_root / "calibration",
                )
                window_size = window_decision["recommended_window"]
            tads = self.env.run_in_stack(
                stack="hic",
                script_path="aria/scripts/hic_topology.py",
                params={
                    "files": files,
                    "genome": exp_ctx.get("genome", "hg38"),
                    "organism": organism,
                    "analysis": "tads",
                    "resolution": resolution,
                    "window_size": int(window_size),
                    "chromosomes": exp_ctx.get("chromosomes", "all"),
                    "out_of_core": True,
                    "output_dir": str(topo_root / "tads"),
                },
            )
            if tads.get("status") == "success":
                self._publish_tads_finding(experiment_id, tads)

        # 3. Loop calling (whenever resolution <= 25kb and not explicitly disabled)
        loops: dict = {}
        finest_avail = min(int(r) for r in avail_res) if avail_res else resolution
        loop_res = resolution if resolution <= 25_000 else finest_avail
        run_loops = (
            exp_ctx.get("enable_loops") is not False
            and (loop_res <= 25_000 or analysis == "loop_calling")
        )
        if run_loops:
            loop_params: dict = {
                "files": files,
                "genome": exp_ctx.get("genome", "hg38"),
                "organism": organism,
                "analysis": "loops",
                "resolution": loop_res,
                "chromosomes": exp_ctx.get("chromosomes", "all"),
                "out_of_core": True,
                "output_dir": str(topo_root / "loops"),
            }
            if gtf_path:
                loop_params["gtf"] = gtf_path
            if exp_ctx.get("loop_fdr_threshold") is not None:
                loop_params["fdr_threshold"] = float(exp_ctx["loop_fdr_threshold"])
            loops = self.env.run_in_stack(
                stack="hic",
                script_path="aria/scripts/hic_topology.py",
                params=loop_params,
            )
            if loops.get("status") == "success":
                self._publish_loops_finding(experiment_id, loops, intent)

        # 4. Multi-replicate differential 3D organization (ΔE1 & ΔInsulation)
        differential_3d: dict = {}
        from aria.utils.design_matrix import resolve_design_contract

        contract = resolve_design_contract(exp_ctx.get("design"), exp_ctx)
        comparisons = (
            contract.get("comparisons")
            or exp_ctx.get("comparisons")
            or intent.get("comparisons")
            or intent.get("comparison")
        )
        sample_meta = (
            exp_ctx.get("sample_metadata")
            or exp_ctx.get("sample_metadata_path")
            or (exp_ctx.get("design") or {}).get("sample_sheet")
        )
        if len(files) >= 2 and sample_meta and comparisons:
            diff_params: dict = {
                "files": files,
                "sample_ids": exp_ctx.get("sample_ids"),
                "sample_metadata": sample_meta,
                "condition_col": contract.get("condition_col") or "condition",
                "replicate_col": contract.get("replicate_col") or "replicate",
                "covariates": contract.get("covariates") or [],
                "comparisons": comparisons,
                "analysis": "differential_3d",
                "genome": exp_ctx.get("genome", "hg38"),
                "organism": organism,
                "comp_resolution": comp_resolution,
                "tad_resolution": resolution,
                "window_size": int(window_size or (resolution * 5)),
                "chromosomes": exp_ctx.get("chromosomes", "all"),
                "output_dir": str(topo_root / "differential_3d"),
            }
            if genome_fasta:
                diff_params["genome_fasta"] = genome_fasta
            if gtf_path:
                diff_params["gtf"] = gtf_path
            if exp_ctx.get("min_replicates_per_condition") is not None:
                diff_params["min_replicates_per_condition"] = int(
                    exp_ctx["min_replicates_per_condition"]
                )
            differential_3d = self.env.run_in_stack(
                stack="hic",
                script_path="aria/scripts/hic_topology.py",
                params={k: v for k, v in diff_params.items() if v is not None},
            )
            if (isinstance(differential_3d, dict)
                    and differential_3d.get("ran")
                    and isinstance(compartments, dict)
                    and compartments.get("status") == "success"):
                compartments["n_ab_switches"] = int(
                    differential_3d.get("n_ab_switches", 0)
                )

        if isinstance(compartments, dict) and compartments.get("status") == "success":
            self._publish_compartments_finding(experiment_id, compartments, intent)

        out = {
            "compartments": compartments,
            "tads": tads,
            "loops": loops,
        }
        if differential_3d:
            out["differential_3d"] = differential_3d
        return out

    # ── Insulation Score window_size ParameterAdvisor ────────────────────

    def _advise_insulation_window(self, experiment_id: str,
                                  intent: dict, files: list,
                                  resolution: int,
                                  genome: str,
                                  out_dir: Optional[Path] = None) -> dict:
        question = intent.get("user_question", "").lower()

        MEGADOMAIN_KEYWORDS = ["megadomain", "large domain", "broad",
                               "topological domain large"]
        SUBDOMAIN_KEYWORDS = ["sub-tad", "subdomain", "fine", "small", "local"]

        if any(k in question for k in MEGADOMAIN_KEYWORDS):
            windows = [300_000, 500_000, 1_000_000]
        elif any(k in question for k in SUBDOMAIN_KEYWORDS):
            windows = [max(resolution * 3, 50_000), max(resolution * 5, 100_000), max(resolution * 8, 200_000)]
        else:
            windows = [resolution * 3, resolution * 5, resolution * 10]

        cal_params: dict = {
            "files": files,
            "genome": genome,
            "analysis": "insulation_calibration",
            "resolution": resolution,
            "windows": windows,
            "chromosomes": ["chr1"],
            "out_of_core": True,
        }
        if out_dir is not None:
            cal_params["output_dir"] = str(out_dir)

        calibration = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_topology.py",
            params=cal_params,
        )

        best_window = windows[1]
        best_strength = 0.0

        if calibration.get("status") == "success":
            for w, strength in calibration.get("boundary_strengths", {}).items():
                if float(strength) > best_strength:
                    best_strength = float(strength)
                    best_window = int(w)

        decision = {
            "recommended_window": best_window,
            "calibration_results": calibration.get("boundary_strengths", {}),
            "windows_tested": windows,
            "rationale": (
                f"window_size={best_window:,}bp selected based on highest "
                f"mean boundary strength ({best_strength:.3f}) across "
                f"{calibration.get('calibration_chrom', 'chr1')}."
            ),
        }

        try:
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=3,
                question=f"Insulation Score window_size at {resolution:,}bp",
                decision=str(best_window),
                rationale=decision["rationale"],
                made_by="advisor",
            )
        except Exception as e:
            log.warning("Could not store window decision: %s", e)

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=self._format_window_checkpoint(decision, resolution),
            options=[
                f"Use recommended window ({best_window:,}bp)",
                "Enter custom window size",
                "Skip TAD calling",
            ],
            context={"window_decision": decision},
        )
        return decision

    # ── Findings publishers ───────────────────────────────────────────────

    def _publish_compartments_finding(self, experiment_id: str,
                                      result: dict, intent: dict):
        n_switches = result.get("n_ab_switches", 0)
        summary = (
            f"Compartment A/B: {result.get('pct_A', 0):.1f}% A-compartment, "
            f"{result.get('pct_B', 0):.1f}% B-compartment"
        )
        if n_switches > 0:
            summary += f", {n_switches} A/B switches detected"
            try:
                from aria.agents.debate_council import DebateCouncil
                council = DebateCouncil(llm=self.llm, max_rounds=2)
                claim = (
                    f"{n_switches} genomic regions show compartment switching "
                    f"(A-to-B or B-to-A) compared to the reference condition."
                )
                debate = council.resolve(
                    topic="compartment_switching",
                    initial_claim=claim,
                    evidence={
                        "n_switches": n_switches,
                        "switch_regions": result.get("switch_regions", [])[:5],
                        "pct_A": result.get("pct_A", 0),
                        "pct_B": result.get("pct_B", 0),
                        "pc1_validated": result.get("pc1_validated", False),
                    },
                    biological_context=intent,
                )
                result["debate_verdict"] = debate.verdict.value
                result["debate_consensus"] = debate.consensus
                result["limitations"] = debate.limitations
            except Exception as e:
                log.warning("DebateCouncil unavailable: %s", e)

        conf = Confidence.HIGH if result.get("pc1_validated") else Confidence.MEDIUM
        self.publish_finding(
            experiment_id,
            {"summary": summary, "n_ab_switches": n_switches},
            conf,
        )

    def _publish_tads_finding(self, experiment_id: str, result: dict):
        n_tads = result.get("n_tads", 0)
        conf = (Confidence.HIGH if n_tads > 100 else
                Confidence.MEDIUM if n_tads >= 5 else
                Confidence.LOW)
        self.publish_finding(
            experiment_id,
            {"summary": f"TADs: {n_tads:,} domains called",
             "n_tads": n_tads,
             "median_tad_size_kb": result.get("median_size_kb", 0)},
            conf,
        )

    def _publish_loops_finding(self, experiment_id: str,
                               result: dict, intent: dict):
        n_loops = result.get("n_loops", 0)
        if n_loops > 0:
            try:
                from aria.agents.debate_council import DebateCouncil
                council = DebateCouncil(llm=self.llm, max_rounds=2)
                top_loops = result.get("top_loops", [])[:5]
                claim = (
                    f"{n_loops:,} chromatin loops identified. "
                    f"Top loops span loci: {', '.join(str(l) for l in top_loops)}"
                )
                debate = council.resolve(
                    topic="chromatin_loop_regulatory_significance",
                    initial_claim=claim,
                    evidence={
                        "n_loops": n_loops,
                        "top_loops": top_loops,
                        "resolution": result.get("resolution", 0),
                        "ctcf_validated": result.get("ctcf_validated", False),
                        "atac_corroborated": result.get("atac_corroborated", False),
                    },
                    biological_context=intent,
                )
                result["debate_verdict"] = debate.verdict.value
                result["debate_consensus"] = debate.consensus
            except Exception as e:
                log.warning("DebateCouncil unavailable: %s", e)

        conf = (Confidence.HIGH if n_loops > 100 else
                Confidence.MEDIUM if n_loops >= 1 else
                Confidence.LOW)
        self.publish_finding(
            experiment_id,
            {"summary": f"Loops: {n_loops:,} chromatin loops called",
             "n_loops": n_loops},
            conf,
        )

    # ── Checkpoint formatters ─────────────────────────────────────────────

    def _format_resolution_checkpoint(self, file_info: dict,
                                      decision: dict) -> str:
        lines = ["Hi-C Resolution Selection\n"]
        lines.append(f"Files: {len(file_info.get('files', []))} Hi-C file(s)")
        lines.append(
            f"Available resolutions: "
            f"{[f'{r:,}bp' for r in file_info.get('available_resolutions', [])]}"
        )
        lines.append("")
        lines.append("Resolution options:")
        for c in decision.get("candidates", []):
            marker = " [RECOMMENDED]" if c.get("recommended") else ""
            lines.append(
                f"  {c['resolution']:>10,}bp  "
                f"RAM ~{c['ram_required_gb']:.1f}GB  "
                f"{c['analysis_depth']}{marker}"
            )
        lines.append("")
        lines.append(f"Rationale: {decision.get('rationale', '')}")
        lines.append("\nSelect resolution to proceed:")
        return "\n".join(lines)

    def _format_qc_checkpoint(self, qc_result: dict) -> str:
        lines = ["Hi-C QC Warnings:\n"]
        for w in qc_result.get("warnings", []):
            lines.append(f"  * {w}")
        lines.append(
            f"\n  Cis/trans ratio: {qc_result.get('cis_trans_ratio', 'N/A')}"
        )
        lines.append(
            f"  Valid pairs: {qc_result.get('n_valid_pairs', 'N/A'):,}"
            if isinstance(qc_result.get("n_valid_pairs"), int) else
            f"  Valid pairs: {qc_result.get('n_valid_pairs', 'N/A')}"
        )
        lines.append("\nHow would you like to proceed?")
        return "\n".join(lines)

    def _format_window_checkpoint(self, decision: dict,
                                  resolution: int) -> str:
        lines = [
            f"Insulation Score Window Size (resolution: {resolution:,}bp)\n"
        ]
        lines.append("Calibration results (chr1 boundary strength):")
        for w, strength in decision.get("calibration_results", {}).items():
            marker = " [RECOMMENDED]" if int(w) == decision["recommended_window"] else ""
            lines.append(
                f"  window={int(w):>10,}bp  "
                f"boundary_strength={float(strength):.4f}{marker}"
            )
        lines.append(f"\n  {decision.get('rationale', '')}")
        lines.append("\nApprove window size to proceed with full genome TAD calling:")
        return "\n".join(lines)

    def _resolution_to_depth(self, resolution: int) -> str:
        DEPTH = {
            1_000_000: "compartments only",
            500_000: "compartments",
            250_000: "compartments",
            100_000: "compartments + coarse TADs",
            50_000: "TADs",
            40_000: "TADs",
            25_000: "TADs + coarse loops",
            10_000: "TADs + loops",
            5_000: "fine loops (HPC recommended)",
        }
        return DEPTH.get(resolution, "custom resolution")

    def receive(self, message):
        pass
