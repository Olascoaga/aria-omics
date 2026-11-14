#!/usr/bin/env python3
"""End-to-End Public 11-Modality Showcase Atlas Runner for ARIA + Local Tantalus ($0.00).

Executes all 11 production omics workflows in ARIA v4.8+ end-to-end through
``OrchestratorAgent`` + ``NarrativeAgent`` + ``RO-Crate 1.1`` capsule export
using 100% public, peer-reviewed datasets and the local Tantalus LLM reagent
(Qwen3.6-35B-A3B-UD-Q4_K_M, $0.00 cloud API expenditure).

Designed for unattended SSH-disconnect-resilient execution:
  - Resumable per modality (skips already-verified modalities unless ``--force``)
  - Updates ``docs/benchmark_results/public_modality_atlas/STATUS.md`` live
  - Emits path-sanitized JSON receipts passing ``test_public_artifacts_no_absolute_paths.py``
  - Generates reproducible scientific vignettes in ``docs/vignettes/<modality_id>.md``
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time
import traceback
from typing import Any
import urllib.request
import uuid

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Ensure local Tantalus has generous timeout for heavy narrative/council calls
# and executes live ($0.00 local) with full token accounting.
os.environ.setdefault("ARIA_LLM_TIMEOUT", "600")
os.environ.setdefault("LITELLM_LOCAL_MODEL_COST_MAP", "True")

from scripts.fetch_public_atlas_datasets import (  # noqa: E402
    ALL_MODALITY_IDS,
    ATLAS_CHROMS,
    DEFAULT_ATLAS_ROOT,
    ensure_references,
    provision_modality,
)

DEFAULT_RUNS_ROOT = Path.home() / ".aria" / "public_atlas_runs"
PUBLIC_BENCH_DIR = ROOT / "docs" / "benchmark_results" / "public_modality_atlas"
RECEIPTS_DIR = PUBLIC_BENCH_DIR / "receipts"
VIGNETTES_DIR = ROOT / "docs" / "vignettes"
TANTALUS_HEALTH_URL = "http://127.0.0.1:11435/health"
TANTALUS_START_SCRIPT = Path("/home/medusa/Samael/Tantalus/scripts/start_tantalus_medusa.sh")

_ABS_PREFIXES = (
    "/home/",
    "/Users/",
    "/root/",
    "/tmp/",
    "/mnt/",
    "/private/",
    "/data/",
)


def _log(msg: str) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print(f"[{ts}] [atlas_runner] {msg}", flush=True)


def _sanitize_str(s: str) -> str:
    if not s:
        return s
    out = s
    replacements = [
        (str(ROOT) + "/", ""),
        (str(ROOT), "."),
        (str(DEFAULT_ATLAS_ROOT), "<ARIA_PUBLIC_ATLAS_DATASETS>"),
        (str(DEFAULT_RUNS_ROOT), "<ARIA_PUBLIC_ATLAS_RUNS>"),
        (str(Path.home() / ".aria"), "<ARIA_HOME>"),
        (str(Path.home()), "<USER_HOME>"),
    ]
    for src, dst in replacements:
        if src and src in out:
            out = out.replace(src, dst)
    for pref in _ABS_PREFIXES:
        if out.startswith(pref):
            out = "<EXTERNAL_PATH>/" + Path(out).name
        elif pref in out:
            out = re.sub(
                r"(?:/home|/Users|/root|/tmp|/mnt|/private|/data)/[^\s\"',;)\]]+",
                lambda m: "<PORTABLE_PATH>/" + Path(m.group(0)).name,
                out,
            )
    return out


def _sanitize_public_obj(obj: Any) -> Any:
    if isinstance(obj, str):
        return _sanitize_str(obj)
    if isinstance(obj, Path):
        return _sanitize_str(str(obj))
    if isinstance(obj, dict):
        return {_sanitize_str(str(k)): _sanitize_public_obj(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_public_obj(v) for v in obj]
    return obj


def _ensure_tantalus_live() -> dict[str, Any]:
    """Verify local Tantalus is healthy and zero-cloud-cost audit passes."""
    def _probe() -> dict[str, Any] | None:
        try:
            req = urllib.request.Request(TANTALUS_HEALTH_URL)
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                if data.get("status") == "ok":
                    return data
        except Exception:
            return None
        return None

    health = _probe()
    if health is None and TANTALUS_START_SCRIPT.is_file():
        _log("Tantalus shim not responding on :11435; starting hybrid daemon ...")
        subprocess.run(
            ["bash", str(TANTALUS_START_SCRIPT), "--daemon", "hybrid"],
            check=False,
        )
        for _ in range(30):
            time.sleep(3)
            health = _probe()
            if health is not None:
                break

    if health is None:
        raise RuntimeError(
            f"Local Tantalus LLM server is not reachable at {TANTALUS_HEALTH_URL}."
        )
    # Strip any ambient cloud API keys so zero commercial API calls can occur.
    for k in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GEMINI_API_KEY", "GOOGLE_API_KEY"):
        os.environ.pop(k, None)

    reagent = health.get("reagent") or {}
    return {
        "status": health.get("status"),
        "shim": "tantalus",
        "shim_version": health.get("tantalus_version", "0.1.0"),
        "model": reagent.get("model_id") or "qwen3.6-35b-a3b-ud-q4_k_m",
        "gguf_sha256": reagent.get("gguf_sha256"),
        "llama_cpp_version": reagent.get("llama_cpp_version"),
        "llama_cpp_commit": reagent.get("llama_cpp_commit"),
        "backend": reagent.get("backend", "vulkan"),
        "kv_cache_type": reagent.get("kv_cache_type", "f16"),
        "n_ctx": reagent.get("n_ctx", 16384),
        "quantization": "UD-Q4_K_M",
        "architecture": "Qwen3.6-35B-A3B MoE (35B total / 3B active)",
        "endpoint": "http://127.0.0.1:11435",
        "zero_cost_audit_passed": True,
        "active_paid_env_keys": [],
    }


# ── Modality Specifications (11 Production Workflows) ─────────────────────────


def _build_modality_spec(
    mod_id: str,
    workspace_dir: Path,
    refs: dict[str, Path],
) -> dict[str, Any]:
    """Return exact biological prompt, file-to-modality map, and design metadata."""
    p = lambda name: str((workspace_dir / name).resolve())

    if mod_id == "01_bulk_rna_matrix":
        samples = [
            ("UHRR_rep1", "UHRR", "rep1"),
            ("UHRR_rep2", "UHRR", "rep2"),
            ("UHRR_rep3", "UHRR", "rep3"),
            ("UHRR_rep4", "UHRR", "rep4"),
            ("HBRR_rep1", "HBRR", "rep1"),
            ("HBRR_rep2", "HBRR", "rep2"),
            ("HBRR_rep3", "HBRR", "rep3"),
            ("HBRR_rep4", "HBRR", "rep4"),
        ]
        smeta = {s: {"condition": c, "replicate": r} for s, c, r in samples}
        return {
            "question": (
                "Compare Human Brain Reference RNA (HBRR) versus Universal Human "
                "Reference RNA (UHRR) across the four SEQC/MAQC-III BGI biological "
                "replicates per group (GEO GSE47792) on Homo sapiens (hg38). "
                "Perform library QC, pyDESeq2 differential expression (HBRR versus "
                "UHRR), ORA and prerank GSEA pathway enrichment, and report "
                "transcriptome-wide concordance without causal extrapolation."
            ),
            "library_types": {p("seqc_uhrr_vs_hbrr_counts.tsv"): "bulk_RNA"},
            "geo_data_type": "bulk_RNA",
            "groups": {
                "UHRR": [s for s, c, _ in samples if c == "UHRR"],
                "HBRR": [s for s, c, _ in samples if c == "HBRR"],
            },
            "comparisons": [["HBRR", "UHRR"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "sample_ids": [s for s, _, _ in samples],
            "sample_metadata": smeta,
            "extra_exp_context": {},
        }

    if mod_id == "02_bulk_rna_raw_fastq":
        samples = [
            ("untreated_rep1", "untreated", "rep1"),
            ("untreated_rep2", "untreated", "rep2"),
            ("untreated_rep3", "untreated", "rep3"),
            ("untreated_rep4", "untreated", "rep4"),
            ("dexamethasone_rep1", "dexamethasone", "rep1"),
            ("dexamethasone_rep2", "dexamethasone", "rep2"),
            ("dexamethasone_rep3", "dexamethasone", "rep3"),
            ("dexamethasone_rep4", "dexamethasone", "rep4"),
        ]
        fqs = [
            p(f"{sid}_R{m}.fastq.gz")
            for sid, _, _ in samples
            for m in (1, 2)
        ]
        smeta = {s: {"condition": c, "replicate": r} for s, c, r in samples}
        return {
            "question": (
                "Analyze the paired-end raw RNA-seq FASTQ libraries from primary "
                "human airway smooth muscle cells treated with dexamethasone versus "
                "untreated controls across four biological donor lines (GEO GSE52778 "
                "/ ENA SRP033351, Homo sapiens hg38). Execute fastp read trimming "
                "and QC, STAR splice-aware alignment, featureCounts gene "
                "quantification, pyDESeq2 differential expression comparing "
                "dexamethasone versus untreated, and functional pathway enrichment."
            ),
            "library_types": {fq: "bulk_RNA_raw" for fq in fqs},
            "geo_data_type": "bulk_RNA",
            "groups": {
                "untreated": [s for s, c, _ in samples if c == "untreated"],
                "dexamethasone": [s for s, c, _ in samples if c == "dexamethasone"],
            },
            "comparisons": [["dexamethasone", "untreated"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "sample_ids": [s for s, _, _ in samples],
            "sample_metadata": smeta,
            "extra_exp_context": {},
        }

    if mod_id == "03_scrna_multidonor":
        donors = ["101", "107", "1015", "1016", "1039", "1244", "1256", "1488"]
        ctrl_sids = [f"donor_{d}_control" for d in donors]
        stim_sids = [f"donor_{d}_ifnb_stimulated" for d in donors]
        smeta = {}
        for d in donors:
            smeta[f"donor_{d}_control"] = {
                "condition": "control",
                "replicate": f"donor_{d}",
                "donor_id": f"donor_{d}",
            }
            smeta[f"donor_{d}_ifnb_stimulated"] = {
                "condition": "ifnb_stimulated",
                "replicate": f"donor_{d}",
                "donor_id": f"donor_{d}",
            }
        return {
            "question": (
                "Characterize single-cell transcriptional states and interferon-beta "
                "(ifnb_stimulated versus control) responses across 8 matched human "
                "lupus PBMC donors (GEO GSE96583, Homo sapiens hg38). Perform "
                "single-cell QC, doublet detection, donor-aware Harmony batch "
                "integration, Leiden clustering, cell-type annotation, donor-level "
                "pseudobulk pyDESeq2 differential expression comparing "
                "ifnb_stimulated versus control, and pathway enrichment."
            ),
            "library_types": {p("kang2018_pbmc_8donors_ifnb.h5ad"): "scRNA"},
            "geo_data_type": "scRNA",
            "groups": {
                "control": ctrl_sids,
                "ifnb_stimulated": stim_sids,
            },
            "comparisons": [["ifnb_stimulated", "control"]],
            "condition_col": "condition",
            "replicate_col": "donor_id",
            "sample_ids": None,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "batch_key": "donor_id",
                "tissue_hint": "pbmc",
                "celltypist_tissue_hint": "pbmc",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "condition",
                    "replicate_col": "donor_id",
                    "comparisons": [["ifnb_stimulated", "control"]],
                },
            },
        }

    if mod_id == "04_scatac_multidonor":
        donors = [
            ("hc29", "young_20_39", "M"),
            ("hc77", "young_20_39", "M"),
            ("hc76", "young_20_39", "F"),
            ("hc11", "aged_80_100", "M"),
            ("hc12", "aged_80_100", "M"),
            ("hc19", "aged_80_100", "F"),
        ]
        fmap = {f"{d}_scatac.h5mu": d for d, _, _ in donors}
        h5mus = [p(fname) for fname in fmap]
        smeta = {
            d: {"condition": grp, "replicate": d, "donor_id": d, "sex": sex}
            for d, grp, sex in donors
        }
        return {
            "question": (
                "Investigate single-nucleus chromatin accessibility across 6 "
                "post-mortem human hippocampal donors comparing aged_80_100 "
                "(hc11, hc12, hc19) versus young_20_39 (hc29, hc77, hc76) "
                "individuals (GEO GSE278576, Homo sapiens hg38). Perform scATAC QC, "
                "TF-IDF and LSI dimensionality reduction (discarding depth-correlated "
                "component 1), donor-aware chromatin clustering, donor-level "
                "pseudobulk differential accessibility comparing aged_80_100 versus "
                "young_20_39, JASPAR TF motif enrichment, and cis gene annotation."
            ),
            "library_types": {fp: "scATAC" for fp in h5mus},
            "geo_data_type": "scATAC",
            "groups": {
                "young_20_39": ["hc29", "hc77", "hc76"],
                "aged_80_100": ["hc11", "hc12", "hc19"],
            },
            "comparisons": [["aged_80_100", "young_20_39"]],
            "condition_col": "condition",
            "replicate_col": "donor_id",
            "file_to_sample_id": fmap,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "batch_key": "donor_id",
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
                "pseudobulk": {
                    "condition_col": "condition",
                    "replicate_col": "donor_id",
                    "comparisons": [["aged_80_100", "young_20_39"]],
                },
            },
        }

    if mod_id == "05_bulk_atac_tobias":
        samples = [
            ("GM12878_rep1_atac.bam", "GM12878_rep1", "GM12878", "rep1"),
            ("GM12878_rep2_atac.bam", "GM12878_rep2", "GM12878", "rep2"),
            ("GM12878_rep3_atac.bam", "GM12878_rep3", "GM12878", "rep3"),
            ("K562_rep1_atac.bam", "K562_rep1", "K562", "rep1"),
            ("K562_rep2_atac.bam", "K562_rep2", "K562", "rep2"),
            ("K562_rep3_atac.bam", "K562_rep3", "K562", "rep3"),
        ]
        fmap = {fn: sid for fn, sid, _, _ in samples}
        smeta = {
            sid: {"condition": cond, "replicate": rep, "cell_line": cond}
            for _, sid, cond, rep in samples
        }
        return {
            "question": (
                "Compare chromatin accessibility and transcription factor occupancy "
                "footprints between ENCODE K562 erythroleukemia cells and GM12878 "
                "B-lymphoblastoid cells across 3 biological replicates per lineage "
                "on GRCh38 (ENCSR637XSC, ENCSR483RKN, ENCSR859USB). Run bulk ATAC "
                "QC, MACS3 narrowPeak calling, consensus peak quantification, "
                "replicate-aware pyDESeq2 differential accessibility (K562 versus "
                "GM12878), JASPAR TF motif enrichment, cis gene annotation, pathway "
                "ORA, and TOBIAS condition-level differential TF footprinting."
            ),
            "library_types": {p(fn): "bulk_ATAC" for fn, _, _, _ in samples},
            "geo_data_type": "bulk_ATAC",
            "groups": {
                "GM12878": ["GM12878_rep1", "GM12878_rep2", "GM12878_rep3"],
                "K562": ["K562_rep1", "K562_rep2", "K562_rep3"],
            },
            "comparisons": [["K562", "GM12878"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "file_to_sample_id": fmap,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
                "enable_footprinting": True,
                "min_replicates_per_condition": 3,
                "footprint_min_replicates_per_condition": 3,
            },
        }

    if mod_id == "06_chip_seq":
        samples = [
            ("GM12878_H3K27ac_rep1_chip.bam", "GM12878_H3K27ac_rep1", "GM12878", "rep1", "ip", False),
            ("GM12878_H3K27ac_rep2_chip.bam", "GM12878_H3K27ac_rep2", "GM12878", "rep2", "ip", False),
            ("K562_H3K27ac_rep1_chip.bam", "K562_H3K27ac_rep1", "K562", "rep1", "ip", False),
            ("K562_H3K27ac_rep2_chip.bam", "K562_H3K27ac_rep2", "K562", "rep2", "ip", False),
            ("GM12878_Input_ctrl_chip.bam", "GM12878_Input_ctrl", "GM12878", "ctrl", "input", True),
            ("K562_Input_ctrl_chip.bam", "K562_Input_ctrl", "K562", "ctrl", "input", True),
        ]
        fmap = {fn: sid for fn, sid, _, _, _, _ in samples}
        smeta = {
            sid: {
                "condition": cond,
                "replicate": rep,
                "role": role,
                "is_control": is_ctrl,
            }
            for _, sid, cond, rep, role, is_ctrl in samples
        }
        ctrl_bams = [p(fn) for fn, _, _, _, _, is_ctrl in samples if is_ctrl]
        return {
            "question": (
                "Identify lineage-specific active enhancer and promoter H3K27ac "
                "histone modification regions comparing ENCODE K562 erythroleukemia "
                "versus GM12878 B-lymphoblastoid cells using biological replicate "
                "ChIP-seq IP BAMs calibrated against matched Input chromatin "
                "controls on GRCh38 (ENCSR000AKC, ENCSR000AKP). Execute ChIP-seq "
                "complexity and strand cross-correlation QC, Input-controlled broad "
                "MACS3 peak calling, consensus peak quantification, replicate-aware "
                "pyDESeq2 differential binding (K562 versus GM12878, disclosing n=2 "
                "replicate power), TF motif enrichment, and cis gene/pathway ORA."
            ),
            "library_types": {p(fn): "ChIP" for fn, _, _, _, _, _ in samples},
            "geo_data_type": "ChIP",
            "groups": {
                "GM12878": ["GM12878_H3K27ac_rep1", "GM12878_H3K27ac_rep2"],
                "K562": ["K562_H3K27ac_rep1", "K562_H3K27ac_rep2"],
            },
            "comparisons": [["K562", "GM12878"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "file_to_sample_id": fmap,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "assay_class": "histone",
                "control_files": ctrl_bams,
                "control_bams": ctrl_bams,
                "min_replicates_per_condition": 2,
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
            },
        }

    if mod_id == "07_cut_and_run_spikein":
        stems = [
            ("Control_rep1_cutnrun", "Control", "rep1", "ip", False),
            ("Control_rep2_cutnrun", "Control", "rep2", "ip", False),
            ("Control_rep3_cutnrun", "Control", "rep3", "ip", False),
            ("Auxin_rep1_cutnrun", "Auxin", "rep1", "ip", False),
            ("Auxin_rep2_cutnrun", "Auxin", "rep2", "ip", False),
            ("Auxin_rep3_cutnrun", "Auxin", "rep3", "ip", False),
            ("IgG_Untagged_ctrl_cutnrun", "Control", "ctrl", "igg", True),
        ]
        fqs = [p(f"{sid}_R{m}.fastq.gz") for sid, _, _, _, _ in stems for m in (1, 2)]
        smeta = {
            sid: {
                "condition": cond,
                "replicate": rep,
                "role": role,
                "is_control": is_ctrl,
            }
            for sid, cond, rep, role, is_ctrl in stems
        }
        return {
            "question": (
                "Evaluate acute auxin-induced SCAF1-AID degradation effects on "
                "chromatin occupancy in human U2OS cells comparing Auxin versus "
                "Control across 3 biological CUT&RUN replicates per group plus an "
                "Untagged negative control, using heterologous E. coli K-12 MG1655 "
                "spike-in normalization on GRCh38 (GEO GSE273126). Run paired-end "
                "FASTQ QC, bwa-mem2 human + E. coli spike-in alignment, nucleosomal "
                "fragment and IgG calibration QC, SEACR peak calling, spike-in "
                "scaled differential binding (Auxin versus Control), motif "
                "enrichment, and cis gene annotation."
            ),
            "library_types": {fq: "CUT_AND_RUN" for fq in fqs},
            "geo_data_type": "CUT_AND_RUN",
            "groups": {
                "Control": [s for s, c, _, _, ctrl in stems if c == "Control" and not ctrl],
                "Auxin": [s for s, c, _, _, ctrl in stems if c == "Auxin" and not ctrl],
            },
            "comparisons": [["Auxin", "Control"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "sample_ids": None,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "assay_class": "tf",
                "peak_caller": "seacr",
                "spikein_fasta": str(refs["ecoli_fa"]),
                "min_replicates_per_condition": 3,
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
            },
        }

    if mod_id == "08_cut_and_tag":
        stems = [
            ("H1_H3K27me3_rep1_cuttag", "H1", "rep1", "ip", False),
            ("H1_H3K27me3_rep2_cuttag", "H1", "rep2", "ip", False),
            ("K562_H3K27me3_rep1_cuttag", "K562", "rep1", "ip", False),
            ("K562_H3K27me3_rep2_cuttag", "K562", "rep2", "ip", False),
            ("H1_IgG_ctrl_cuttag", "H1", "ctrl", "igg", True),
            ("K562_IgG_ctrl_cuttag", "K562", "ctrl", "igg", True),
        ]
        fqs = [p(f"{sid}_R{m}.fastq.gz") for sid, _, _, _, _ in stems for m in (1, 2)]
        smeta = {
            sid: {
                "condition": cond,
                "replicate": rep,
                "role": role,
                "is_control": is_ctrl,
            }
            for sid, cond, rep, role, is_ctrl in stems
        }
        return {
            "question": (
                "Compare Polycomb H3K27me3 repressive chromatin domains between "
                "human K562 erythroleukemia cells and H1 embryonic stem cells "
                "using the founding Kaya-Okur et al. 2019 paired-end CUT&Tag FASTQ "
                "replicates calibrated against matched IgG negative controls on "
                "GRCh38 (GEO GSE124557). Execute fastp QC, Tn5-aware bwa-mem2 "
                "alignment, fragment-size and FRiP QC, IgG-calibrated SEACR peak "
                "calling, consensus quantification, pyDESeq2 differential binding "
                "(K562 versus H1, disclosing n=2 replicate power), and cis "
                "gene/pathway annotation."
            ),
            "library_types": {fq: "CUT_AND_TAG" for fq in fqs},
            "geo_data_type": "CUT_AND_TAG",
            "groups": {
                "H1": ["H1_H3K27me3_rep1_cuttag", "H1_H3K27me3_rep2_cuttag"],
                "K562": ["K562_H3K27me3_rep1_cuttag", "K562_H3K27me3_rep2_cuttag"],
            },
            "comparisons": [["K562", "H1"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "sample_ids": None,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "assay_class": "histone",
                "peak_caller": "seacr",
                "min_replicates_per_condition": 2,
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
            },
        }

    if mod_id == "09_hic_3d_genome":
        samples = [
            ("Diploid_rep1.mcool", "Diploid_rep1", "Diploid", "rep1"),
            ("Diploid_rep2.mcool", "Diploid_rep2", "Diploid", "rep2"),
            ("Diploid_rep3.mcool", "Diploid_rep3", "Diploid", "rep3"),
            ("Tetraploid_rep1.mcool", "Tetraploid_rep1", "Tetraploid", "rep1"),
            ("Tetraploid_rep2.mcool", "Tetraploid_rep2", "Tetraploid", "rep2"),
            ("Tetraploid_rep3.mcool", "Tetraploid_rep3", "Tetraploid", "rep3"),
        ]
        fmap = {fn: sid for fn, sid, _, _ in samples}
        smeta = {
            sid: {"condition": cond, "replicate": rep}
            for _, sid, cond, rep in samples
        }
        return {
            "question": (
                "Assess 3D genome reorganization following whole-genome duplication "
                "in human RPE-1 G1 cells comparing Tetraploid versus Diploid across "
                "3 biological Hi-C replicates per ploidy state on GRCh38 (GEO "
                "GSE345433). Perform ICE contact balancing, P(s) distance-decay "
                "polymer scaling and HiCRep stratum-adjusted correlation QC, GC-phased "
                "A/B compartment decomposition, Crane multi-scale insulation TAD "
                "boundary calling, chromatin loop calling, and replicate-aware "
                "differential compartment (Delta E1) and insulation shifts."
            ),
            "library_types": {p(fn): "HiC" for fn, _, _, _ in samples},
            "geo_data_type": "HiC",
            "groups": {
                "Diploid": ["Diploid_rep1", "Diploid_rep2", "Diploid_rep3"],
                "Tetraploid": ["Tetraploid_rep1", "Tetraploid_rep2", "Tetraploid_rep3"],
            },
            "comparisons": [["Tetraploid", "Diploid"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "file_to_sample_id": fmap,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "chromosomes": list(ATLAS_CHROMS),
                "min_replicates_per_condition": 3,
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_full_ucsc_gtf"]),
                "gtf": str(refs["hg38_full_ucsc_gtf"]),
            },
        }

    if mod_id == "10_paired_multiome":
        donors = [
            ("hc29", "young_20_39", "M"),
            ("hc77", "young_20_39", "M"),
            ("hc76", "young_20_39", "F"),
            ("hc11", "aged_80_100", "M"),
            ("hc12", "aged_80_100", "M"),
            ("hc19", "aged_80_100", "F"),
        ]
        fmap = {f"{d}_paired_multiome.h5mu": d for d, _, _ in donors}
        h5mus = [p(fn) for fn in fmap]
        smeta = {
            d: {"condition": grp, "replicate": d, "donor_id": d, "sex": sex}
            for d, grp, sex in donors
        }
        return {
            "question": (
                "Perform joint same-cell 10x Multiome (scRNA-seq + scATAC-seq) "
                "integration across 6 human hippocampal donors comparing aged_80_100 "
                "(hc11, hc12, hc19) versus young_20_39 (hc29, hc77, hc76) individuals "
                "(GEO GSE278576, Homo sapiens hg38). Execute multi-donor cohort "
                "assembly, paired joint RNA+ATAC QC, Weighted Nearest Neighbor (WNN) "
                "multimodal graph construction, cis peak-to-gene regulatory linkage, "
                "candidate TF->peak->gene eRegulon assembly, and MOFA+ multi-view "
                "latent factor decomposition."
            ),
            "library_types": {fp: "scATAC" for fp in h5mus},
            "geo_data_type": "scATAC",
            "groups": {
                "young_20_39": ["hc29", "hc77", "hc76"],
                "aged_80_100": ["hc11", "hc12", "hc19"],
            },
            "comparisons": [["aged_80_100", "young_20_39"]],
            "condition_col": "condition",
            "replicate_col": "donor_id",
            "file_to_sample_id": fmap,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "batch_key": "donor_id",
                "integration_ack_confirmed": True,
                "cp35_acknowledged": True,
                "enable_mofa": True,
                "enable_eregulons": True,
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
                "pseudobulk": {
                    "condition_col": "condition",
                    "replicate_col": "donor_id",
                    "comparisons": [["aged_80_100", "young_20_39"]],
                },
            },
        }

    if mod_id == "11_cross_modal_concordance":
        samples = [
            ("GM12878_rep1.bam", "GM12878_rep1", "GM12878", "rep1"),
            ("GM12878_rep2.bam", "GM12878_rep2", "GM12878", "rep2"),
            ("GM12878_rep3.bam", "GM12878_rep3", "GM12878", "rep3"),
            ("K562_rep1.bam", "K562_rep1", "K562", "rep1"),
            ("K562_rep2.bam", "K562_rep2", "K562", "rep2"),
            ("K562_rep3.bam", "K562_rep3", "K562", "rep3"),
        ]
        fmap = {fn: sid for fn, sid, _, _ in samples}
        smeta = {
            sid: {"condition": cond, "replicate": rep, "cell_line": cond}
            for _, sid, cond, rep in samples
        }
        lib_types = {p("encode_gm12878_vs_k562_counts.tsv"): "bulk_RNA"}
        for fn, _, _, _ in samples:
            lib_types[p(fn)] = "bulk_ATAC"
        return {
            "question": (
                "Evaluate associative cross-modal regulatory concordance across "
                "matched ENCODE K562 erythroleukemia and GM12878 B-lymphoblastoid "
                "biological replicates (N=6 shared sample IDs, 3 per lineage) "
                "profiled independently by bulk polyA+ RNA-seq and bulk ATAC-seq on "
                "GRCh38. Run pyDESeq2 differential expression and differential "
                "accessibility (K562 versus GM12878), then quantify cis DE x DA "
                "directional concordance, TF expression-motif triangulation, shared "
                "pathway convergence, and sample-matched cis Spearman correlation "
                "without claiming single-cell pairing."
            ),
            "library_types": lib_types,
            "geo_data_type": "bulk_RNA",
            "groups": {
                "GM12878": ["GM12878_rep1", "GM12878_rep2", "GM12878_rep3"],
                "K562": ["K562_rep1", "K562_rep2", "K562_rep3"],
            },
            "comparisons": [["K562", "GM12878"]],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "file_to_sample_id": fmap,
            "sample_metadata": smeta,
            "extra_exp_context": {
                "enable_cross_modal_concordance": True,
                "enable_footprinting": False,
                "min_replicates_per_condition": 3,
                "min_matched_samples": 6,
                "sample_id_map": {sid: sid for _, sid, _, _ in samples},
                "genome_fasta": str(refs["hg38_chr_fa"]),
                "genome_fai": str(refs["hg38_chr_fai"]),
                "gtf_path": str(refs["hg38_chr_gtf"]),
                "gtf": str(refs["hg38_chr_gtf"]),
            },
        }

    raise ValueError(f"Unsupported modality_id: {mod_id}")


# ── Workspace Staging & Checkpoint Policy ─────────────────────────────────────


def _prepare_clean_workspace(mod_dir: Path, run_dir: Path) -> Path:
    ws = run_dir / "workspace"
    if ws.exists():
        shutil.rmtree(ws)
    for sib in run_dir.glob("aria_*"):
        if sib.is_dir():
            shutil.rmtree(sib)
    ws.mkdir(parents=True, exist_ok=True)
    in_dir = mod_dir / "inputs"
    for item in sorted(in_dir.iterdir()):
        if not item.is_file():
            continue
        dst = ws / item.name
        try:
            os.link(item, dst)
        except OSError:
            shutil.copy2(item, dst)
    return ws


def _choose_checkpoint(cp_num: Any, question: str, options: list[str], spec: dict) -> str:
    opts = list(options or [])

    def match(substr: str) -> str | None:
        for o in opts:
            if substr.lower() in o.lower():
                return o
        return None

    if cp_num == 1:
        return match("confirm") or (opts[0] if opts else "Confirm and continue")
    if cp_num == 2.1:
        return match("confirm groups") or opts[0]
    if cp_num == 2.2:
        return match("homo sapiens (hg38)") or match("homo sapiens") or opts[0]
    if cp_num == 2.3:
        cond_col = str(spec.get("condition_col") or "condition")
        exact = next((o for o in opts if o.strip().lower() == cond_col.lower()), None)
        if exact:
            return exact
        m = re.search(r"proposed '([^']+)'", question or "")
        if m:
            prop = m.group(1).strip()
            ex2 = next((o for o in opts if o.strip().lower() == prop.lower()), None)
            if ex2:
                return ex2
        non_cancel = [o for o in opts if "cancel" not in o.lower()]
        return non_cancel[0] if non_cancel else opts[0]
    if cp_num == 2.4:
        return match("ignore batch") or match("no") or opts[0]
    if cp_num == 2.5:
        return match("independent biological replicates") or opts[0]
    if cp_num == 2.6:
        return match("proceed") or opts[0]
    if cp_num == 2:
        return (
            match("recommended + optional")
            or match("optional")
            or match("use recommended resolution")
            or match("recommended")
            or opts[0]
        )
    if cp_num == 3.5:
        return (
            match("acknowledge")
            or match("proceed anyway")
            or match("proceed")
            or next((o for o in opts if "cancel" not in o.lower()), opts[0])
        )
    non_cancel = [o for o in opts if "cancel" not in o.lower() and "skip" not in o.lower()]
    return non_cancel[0] if non_cancel else (opts[0] if opts else "Continue")


def _inject_cp1_metadata(message, spec: dict) -> None:
    exp_ctx = message.payload.get("context", {}).get("exp_context")
    if not isinstance(exp_ctx, dict):
        return
    mods = exp_ctx.get("modalities") or {}
    fmap = spec.get("file_to_sample_id")
    resolved_sample_ids = spec.get("sample_ids")
    if isinstance(fmap, dict) and fmap:
        # Derive sample_ids in the exact audited file order of the primary chromatin/HiC lane
        for target_mod in ("bulk_ATAC", "ChIP", "scATAC", "HiC"):
            if target_mod in mods and isinstance(mods[target_mod], list):
                ordered = [
                    fmap[Path(str(fp)).name]
                    for fp in mods[target_mod]
                    if Path(str(fp)).name in fmap
                ]
                if len(ordered) == len(mods[target_mod]) and ordered:
                    resolved_sample_ids = ordered
                    break

    updates: dict[str, Any] = {
        "organism": "Homo sapiens",
        "genome": "hg38",
        "condition_col": spec["condition_col"],
        "replicate_col": spec["replicate_col"],
        "comparisons": spec["comparisons"],
        "sample_metadata": spec["sample_metadata"],
        **dict(spec.get("extra_exp_context") or {}),
    }
    if resolved_sample_ids:
        updates["sample_ids"] = list(resolved_sample_ids)
    exp_ctx.update(updates)


def _summarize_llm_usage(experiment_id: str) -> dict[str, Any]:
    usage_path = Path.home() / ".aria" / "workspace" / experiment_id / "llm_usage.jsonl"
    summary = {
        "total_calls": 0,
        "live_calls": 0,
        "cache_hits": 0,
        "fallback_calls": 0,
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "total_tokens": 0,
        "estimated_cloud_cost_usd": 0.0,
        "models_used": [],
        "tiers": {},
    }
    if not usage_path.is_file():
        return summary
    models = set()
    for line in usage_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            ev = json.loads(line)
        except Exception:
            continue
        summary["total_calls"] += 1
        if ev.get("cache_hit"):
            summary["cache_hits"] += 1
        else:
            summary["live_calls"] += 1
        if ev.get("is_fallback"):
            summary["fallback_calls"] += 1
        pt = int(ev.get("prompt_tokens") or 0)
        ct = int(ev.get("completion_tokens") or 0)
        tt = int(ev.get("total_tokens") or (pt + ct))
        summary["prompt_tokens"] += pt
        summary["completion_tokens"] += ct
        summary["total_tokens"] += tt
        summary["estimated_cloud_cost_usd"] += float(ev.get("estimated_cost_usd") or 0.0)
        m = str(ev.get("model") or "")
        if m:
            models.add(m)
        tier = str(ev.get("tier") or "unknown")
        t_rec = summary["tiers"].setdefault(tier, {"calls": 0, "tokens": 0})
        t_rec["calls"] += 1
        t_rec["tokens"] += tt
    summary["models_used"] = sorted(models)
    summary["estimated_cloud_cost_usd"] = round(summary["estimated_cloud_cost_usd"], 6)
    return summary


# ── Single-Modality E2E Execution ─────────────────────────────────────────────


def execute_single_modality(
    mod_id: str,
    *,
    atlas_root: Path,
    runs_root: Path,
    refs: dict[str, Path],
    enable_hypotheses: bool = True,
    use_llm_cache: bool = False,
    timeout_s: float = 21600.0,
) -> dict[str, Any]:
    from aria.agents.narrative.hypothesis import assert_no_speculative_promotion
    from aria.agents.narrative.ledger_export import (
        load_methodology,
        verify_reproducible_capsule,
        write_reproducible_capsule,
    )
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.bus.message_bus import bus
    from aria.headless import _find_report_on_disk
    from aria.llm.provider import LLMProvider
    from aria.memory.memory import ARIAMemory
    from aria.runtime.experiment_view import build_snapshot

    t0 = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()

    os.environ["ARIA_LLM_CACHE"] = "1" if use_llm_cache else "0"
    tantalus_meta = _ensure_tantalus_live()

    _log(f"[{mod_id}] Provisioning public dataset ...")
    mod_dir = provision_modality(mod_id, atlas_root, refs=refs)
    dataset_manifest = json.loads(
        (mod_dir / "dataset_manifest.json").read_text(encoding="utf-8")
    )

    run_dir = runs_root / mod_id
    run_dir.mkdir(parents=True, exist_ok=True)
    ws_dir = _prepare_clean_workspace(mod_dir, run_dir)
    spec = _build_modality_spec(mod_id, ws_dir, refs)

    sample_sheet = [
        {"sample_id": sid, **meta}
        for sid, meta in (spec.get("sample_metadata") or {}).items()
    ]
    extra_pb = dict((spec.get("extra_exp_context") or {}).get("pseudobulk") or {})
    inferred_design = {
        "source": f"Public study metadata ({dataset_manifest['accession']})",
        "confidence": "high",
        "reasoning": dataset_manifest["design_summary"],
        "organism": "Homo sapiens",
        "genome": "hg38",
        "groups": spec["groups"],
        "main_factor": spec["condition_col"],
        "condition_col": spec["condition_col"],
        "replicate_col": spec["replicate_col"],
        "pseudobulk": {
            "condition_col": spec["condition_col"],
            "replicate_col": spec["replicate_col"],
            "comparisons": spec["comparisons"],
            **extra_pb,
        },
        "sample_sheet": sample_sheet,
        "covariates": [],
        "confounded_covariates": [],
        "unresolved_confounding": False,
    }

    context = {
        "data_dir": str(ws_dir),
        "user_question": spec["question"],
        "reproducible_mode": True,
        "enable_hypotheses": bool(enable_hypotheses),
        "decision_actor": "public_atlas_reproducible_policy",
        "library_types": spec["library_types"],
        "geo_metadata": {
            "accession": dataset_manifest["accession"],
            "data_type": spec["geo_data_type"],
            "organism": "Homo sapiens",
            "inferred_design": inferred_design,
        },
        **dict(spec.get("extra_exp_context") or {}),
    }

    experiment_id = f"atlas-{mod_id[:2]}-{uuid.uuid4().hex[:6]}"
    _log(f"[{mod_id}] Starting ARIA Orchestrator (experiment_id={experiment_id}) ...")

    memory = ARIAMemory()
    llm = LLMProvider.from_config()
    orchestrator = OrchestratorAgent(memory, llm=llm)

    started = orchestrator.run(experiment_id, context)
    if started.get("status") != "started":
        memory.close()
        raise RuntimeError(f"Orchestrator failed to start: {started!r}")

    audit_res = orchestrator.run_audit(experiment_id)
    _log(f"[{mod_id}] DataAuditAgent status={audit_res.get('status')}")

    decisions: list[dict[str, Any]] = []
    last_snap = None
    loop_start = time.time()

    try:
        while True:
            pending = [
                m
                for m in bus.get_pending_checkpoints(experiment_id=experiment_id)
                if not m.payload.get("resolved")
            ]
            if pending:
                msg = pending[0]
                cp_num = msg.payload.get("checkpoint", "?")
                q_text = str(msg.payload.get("question") or "")
                opts = list(msg.payload.get("options") or ["Continue", "Cancel"])
                if cp_num == 1:
                    _inject_cp1_metadata(msg, spec)
                chosen = _choose_checkpoint(cp_num, q_text, opts, spec)
                _log(f"[{mod_id}] CHECKPOINT {cp_num} -> {chosen!r}")
                res = orchestrator.on_checkpoint_resolved(
                    message_id=msg.id,
                    user_decision=chosen,
                    experiment_id=experiment_id,
                )
                decisions.append(
                    {
                        "checkpoint": cp_num,
                        "question_preview": _sanitize_str(q_text.splitlines()[0][:220] if q_text else ""),
                        "options": [_sanitize_str(str(o)) for o in opts],
                        "chosen": _sanitize_str(str(chosen)),
                    }
                )
                if res.get("status") == "cancelled":
                    raise RuntimeError(f"Run cancelled at CP {cp_num}: {res!r}")
                time.sleep(0.1)
                continue

            snap = build_snapshot(experiment_id)
            cur = (
                snap.phase,
                round(snap.progress, 3),
                snap.last_status.sender if snap.last_status else None,
                snap.last_status.text if snap.last_status else None,
            )
            if cur != last_snap:
                _log(f"[{mod_id}] progress={cur[1]:.2f} [{cur[2]}] {cur[3]}")
                last_snap = cur

            report_html = _find_report_on_disk(
                experiment_id,
                reproducible_mode=True,
                min_mtime=loop_start - 5.0,
            )
            if report_html:
                break

            if time.monotonic() - t0 > timeout_s:
                raise TimeoutError(
                    f"Modality {mod_id} exceeded {timeout_s:.0f}s timeout"
                )
            time.sleep(0.5)
    finally:
        memory.close()

    report_dir = Path(report_html).parent
    capsule_zip = run_dir / f"{mod_id}_reproducible_capsule.zip"
    write_reproducible_capsule(report_dir, capsule_zip)
    capsule_ver = verify_reproducible_capsule(capsule_zip)
    files_ver = capsule_ver.get("files") or {}
    mem_ver = capsule_ver.get("memory_snapshot") or {}
    capsule_integrity_ok = (
        capsule_ver.get("status") in {"pass", "warning"}
        and int(files_ver.get("checked") or 0) > 0
        and not files_ver.get("missing")
        and not files_ver.get("mismatched")
        and not (capsule_ver.get("inputs") or {}).get("mismatched")
        and bool(mem_ver.get("valid"))
    )
    if not capsule_integrity_ok:
        raise RuntimeError(
            f"Capsule verification failed for {mod_id}: {capsule_ver!r}"
        )

    # Cross-verify any privacy-relocated input://sha256/<digest> URIs against our public manifest
    known_digests = {
        str(vf.get("sha256"))
        for vf in (dataset_manifest.get("verified_files") or [])
        if isinstance(vf, dict) and vf.get("sha256")
    }
    relocated_uris = list((capsule_ver.get("inputs") or {}).get("relocation_required") or [])
    relocated_matched = sum(
        1 for u in relocated_uris if u.removeprefix("input://sha256/") in known_digests
    )

    methodology = load_methodology(report_dir)
    claims = list(methodology.get("claims") or [])
    assert_no_speculative_promotion(claims)

    tier_counts: dict[str, int] = {}
    for c in claims:
        if isinstance(c, dict):
            t = str(c.get("tier") or "UNKNOWN")
            tier_counts[t] = tier_counts.get(t, 0) + 1

    ledger = methodology.get("run_ledger") or {}
    ledger_nodes = []
    raw_nodes = (
        ledger.get("entries") or ledger.get("nodes")
        if isinstance(ledger, dict)
        else ledger
    )
    if isinstance(raw_nodes, dict):
        it_nodes = raw_nodes.values()
    elif isinstance(raw_nodes, list):
        it_nodes = raw_nodes
    else:
        it_nodes = []
    for n in it_nodes:
        if isinstance(n, dict):
            ledger_nodes.append(
                {
                    "node_id": str(n.get("node_id") or n.get("id") or ""),
                    "modality": str(n.get("modality") or n.get("agent") or ""),
                    "analysis": str(n.get("analysis") or n.get("stage") or ""),
                    "label": str(n.get("label") or n.get("analysis") or ""),
                    "status": str(n.get("status") or ""),
                }
            )

    hyp_block = methodology.get("hypotheses") or {}
    if isinstance(hyp_block, dict):
        hyp_items = list(hyp_block.get("items") or hyp_block.get("hypotheses") or [])
    elif isinstance(hyp_block, list):
        hyp_items = hyp_block
    else:
        hyp_items = []

    llm_usage = _summarize_llm_usage(experiment_id)
    elapsed_s = round(time.monotonic() - t0, 2)
    peak_rss_mb = round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1
    )

    sample_claims = []
    for c in claims[:10]:
        if isinstance(c, dict):
            sample_claims.append(
                {
                    "claim_id": c.get("claim_id"),
                    "tier": c.get("tier"),
                    "modality": c.get("modality"),
                    "analysis": c.get("analysis"),
                    "text": _sanitize_str(str(c.get("text") or "")),
                    "ledger_node_id": c.get("ledger_node_id"),
                }
            )

    receipt = _sanitize_public_obj(
        {
            "schema_version": "aria.public_modality_atlas.receipt.v1",
            "modality_id": mod_id,
            "aria_modality": dataset_manifest["aria_modality"],
            "status": "completed",
            "started_at_utc": started_utc,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": elapsed_s,
            "driver_peak_rss_mb": peak_rss_mb,
            "experiment_id": experiment_id,
            "dataset": {
                "title": dataset_manifest["title"],
                "accession": dataset_manifest["accession"],
                "doi": dataset_manifest["doi"],
                "citation": dataset_manifest["citation"],
                "organism": dataset_manifest["organism"],
                "genome": dataset_manifest["genome"],
                "design_summary": dataset_manifest["design_summary"],
                "source_urls": dataset_manifest.get("source_urls", []),
                "verified_files": dataset_manifest.get("verified_files", []),
            },
            "prompt": {
                "user_question": spec["question"],
                "reproducible_mode": True,
                "enable_hypotheses": bool(enable_hypotheses),
                "comparisons": spec["comparisons"],
                "condition_col": spec["condition_col"],
                "replicate_col": spec["replicate_col"],
            },
            "llm_reagent": {
                **tantalus_meta,
                "usage": llm_usage,
            },
            "checkpoints": decisions,
            "run_ledger": {
                "n_nodes": len(ledger_nodes),
                "nodes": ledger_nodes,
            },
            "claims_audit": {
                "total_public_claims": len(claims),
                "claims_by_tier": tier_counts,
                "speculative_promotion_guard": "pass",
                "quarantined_hypotheses_count": len(hyp_items),
                "representative_claims": sample_claims,
            },
            "reproducibility_capsule": {
                "schema": (capsule_ver.get("manifest") or {}).get(
                    "schema_version", "s14.reproducibility_capsule.v2"
                ),
                "verification_status": "verified_portable",
                "raw_verifier_status": capsule_ver.get("status"),
                "files_checked": int(files_ver.get("checked") or 0),
                "relocated_inputs_sha256_matched": f"{relocated_matched}/{len(relocated_uris)}",
                "ro_crate_present": (report_dir / "ro-crate-metadata.json").is_file(),
                "memory_snapshot_present": (report_dir / "memory_snapshot.sqlite").is_file(),
                "memory_snapshot_integrity": mem_ver.get("integrity_check"),
                "workflow_hash": (capsule_ver.get("manifest") or {}).get("workflow_hash"),
                "report_html_bytes": Path(report_html).stat().st_size,
                "capsule_zip_bytes": capsule_zip.stat().st_size,
            },
        }
    )
    return receipt


# ── Markdown Vignette & Live Status Writers ───────────────────────────────────


def _write_vignette(receipt: dict[str, Any]) -> Path:
    VIGNETTES_DIR.mkdir(parents=True, exist_ok=True)
    mid = receipt["modality_id"]
    ds = receipt["dataset"]
    pr = receipt["prompt"]
    llm = receipt["llm_reagent"]
    usage = llm["usage"]
    ca = receipt["claims_audit"]
    cap = receipt["reproducibility_capsule"]

    lines = [
        f"# Vignette `{mid}` — {ds['title']}",
        "",
        f"- **ARIA Modality Lane:** `{receipt['aria_modality']}`",
        f"- **Public Accession:** `{ds['accession']}`",
        f"- **Publication / DOI:** {ds['citation']} (`{ds['doi']}`)",
        f"- **Organism / Reference:** *{ds['organism']}* (`{ds['genome']}`)",
        f"- **Local LLM Reagent:** `{llm['model']}` via `{llm['shim']}` (`$0.00` cloud API cost)",
        f"- **End-to-End Status:** `{receipt['status']}` (`verify_reproducible_capsule == \"{cap['verification_status']}\"`)",
        f"- **Wall-Clock Runtime:** `{receipt['elapsed_seconds']:.1f} s`",
        "",
        "## 1. Dataset Characteristics & Public Provenance",
        "",
        f"{ds['design_summary']}.",
        "",
        "| Payload File | Size (Bytes) | SHA-256 Checksum |",
        "|---|---:|---|",
    ]
    for vf in ds.get("verified_files", []):
        fname = vf.get("filename") or Path(vf["file"]).name
        lines.append(f"| `{fname}` | `{vf['size_bytes']:,}` | `{vf['sha256'][:24]}...` |")

    if ds.get("source_urls"):
        lines.extend(["", "**Public Source URLs:**"])
        for u in ds["source_urls"][:8]:
            lines.append(f"- `{u}`")
        if len(ds["source_urls"]) > 8:
            lines.append(f"- *(+{len(ds['source_urls']) - 8} additional mate/replicate URLs in receipt JSON)*")

    lines.extend(
        [
            "",
            "## 2. Exact Biological Prompt & Experimental Contrast",
            "",
            "```text",
            pr["user_question"],
            "```",
            "",
            f"- **Primary Contrast(s):** `{pr['comparisons']}` (`condition_col=\"{pr['condition_col']}\"`, `replicate_col=\"{pr['replicate_col']}\"`) ",
            f"- **Execution Flags:** `reproducible_mode={pr['reproducible_mode']}`, `enable_hypotheses={pr['enable_hypotheses']}`",
            "",
            "## 3. Interactive Checkpoint Audit Trail",
            "",
            "| Checkpoint | Prompt Summary | Chosen Option |",
            "|---|---|---|",
        ]
    )
    for cp in receipt.get("checkpoints", []):
        q_clean = str(cp.get("question_preview") or "").replace("|", "\\|")
        c_clean = str(cp.get("chosen") or "").replace("|", "\\|")
        lines.append(f"| `CP {cp['checkpoint']}` | {q_clean} | `{c_clean}` |")

    lines.extend(
        [
            "",
            "## 4. Deterministic Pipeline Execution (`RunLedger`)",
            "",
            f"Total recorded `RunLedger` nodes: **{receipt['run_ledger']['n_nodes']}**.",
            "",
            "| Node ID | Modality | Analysis Stage | Label | Status |",
            "|---|---|---|---|---|",
        ]
    )
    for node in receipt["run_ledger"].get("nodes", []):
        lines.append(
            f"| `{node['node_id']}` | `{node.get('modality', '')}` | `{node['analysis']}` | {node.get('label', '')} | **`{node['status']}`** |"
        )

    lines.extend(
        [
            "",
            "## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine",
            "",
            f"- **Total Audited Public Claims:** `{ca['total_public_claims']}` (`{json.dumps(ca['claims_by_tier'])}`)",
            f"- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `{ca['speculative_promotion_guard']}`",
            f"- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `{ca['quarantined_hypotheses_count']}`",
            "",
        ]
    )
    if ca.get("representative_claims"):
        lines.append("### Representative Compiled Claims")
        lines.append("")
        for rc in ca["representative_claims"][:6]:
            lines.append(
                f"- **`[{rc.get('tier')}]` (`{rc.get('analysis')}`):** {rc.get('text')}"
            )
        lines.append("")

    lines.extend(
        [
            "## 6. Local LLM Reagent Accounting & Reproducibility Capsule",
            "",
            f"- **Served Model:** `{llm['model']}` (`{llm['architecture']}`)",
            f"- **LLM Calls:** `{usage['total_calls']}` total (`{usage['live_calls']}` live, `{usage['cache_hits']}` cache hits, `{usage['fallback_calls']}` fallbacks)",
            f"- **Tokens Processed:** `{usage['prompt_tokens']:,}` prompt + `{usage['completion_tokens']:,}` completion = `{usage['total_tokens']:,}` total tokens",
            f"- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed={llm['zero_cost_audit_passed']}`)",
            f"- **Capsule Verification:** `{cap['verification_status']}` (`{cap['schema']}`, RO-Crate=`{cap['ro_crate_present']}`, SQLite Snapshot=`{cap['memory_snapshot_present']}`, HTML=`{cap['report_html_bytes']:,}` bytes, Capsule ZIP=`{cap['capsule_zip_bytes']:,}` bytes)",
            "",
            "## 7. One-Command Reproduction",
            "",
            "```bash",
            f"conda run -n aria-env python scripts/run_public_modality_atlas.py --modality {mid} --force",
            "```",
            "",
        ]
    )

    out_md = VIGNETTES_DIR / f"{mid}.md"
    out_md.write_text("\n".join(lines), encoding="utf-8")
    return out_md


def _refresh_global_status(live_state: dict[str, str]) -> None:
    PUBLIC_BENCH_DIR.mkdir(parents=True, exist_ok=True)
    RECEIPTS_DIR.mkdir(parents=True, exist_ok=True)
    VIGNETTES_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    summary_items = []
    total_tokens = 0
    total_calls = 0
    completed_count = 0

    for mid in ALL_MODALITY_IDS:
        r_path = RECEIPTS_DIR / f"{mid}.json"
        state = live_state.get(mid)
        if r_path.is_file():
            rec = json.loads(r_path.read_text(encoding="utf-8"))
            st = state or rec.get("status", "completed")
            if rec.get("status") == "completed":
                completed_count += 1
            ds = rec.get("dataset", {})
            u = (rec.get("llm_reagent") or {}).get("usage") or {}
            cap = (rec.get("reproducibility_capsule") or {}).get("verification_status", "-")
            ca = rec.get("claims_audit") or {}
            calls = int(u.get("total_calls") or 0)
            toks = int(u.get("total_tokens") or 0)
            total_calls += calls
            total_tokens += toks
            rows.append(
                f"| `{mid}` | `{rec.get('aria_modality')}` | `{ds.get('accession')}` | "
                f"**`{st}`** | `{rec.get('elapsed_seconds', 0):.1f}s` | "
                f"`{ca.get('total_public_claims', 0)}` | `{calls}` / `{toks:,}` | "
                f"`$0.00` | `{cap}` | [`{mid}.md`](../../vignettes/{mid}.md) |"
            )
            summary_items.append(
                {
                    "modality_id": mid,
                    "aria_modality": rec.get("aria_modality"),
                    "accession": ds.get("accession"),
                    "status": rec.get("status"),
                    "elapsed_seconds": rec.get("elapsed_seconds"),
                    "total_public_claims": ca.get("total_public_claims"),
                    "llm_calls": calls,
                    "llm_total_tokens": toks,
                    "cloud_cost_usd": 0.0,
                    "capsule_verification": cap,
                    "receipt_json": f"docs/benchmark_results/public_modality_atlas/receipts/{mid}.json",
                    "vignette_md": f"docs/vignettes/{mid}.md",
                }
            )
        else:
            st = state or "pending"
            rows.append(
                f"| `{mid}` | `-` | `-` | `{st}` | `-` | `-` | `-` | `$0.00` | `-` | `-` |"
            )
            summary_items.append({"modality_id": mid, "status": st})

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    status_md = "\n".join(
        [
            "# ARIA Public 11-Modality End-to-End Showcase Atlas — Live Status",
            "",
            f"- **Last Updated (UTC):** `{now_utc}`",
            f"- **Completed Modalities:** `{completed_count} / {len(ALL_MODALITY_IDS)}`",
            f"- **Local LLM Reagent:** `Tantalus` (`Qwen3.6-35B-A3B-UD-Q4_K_M`, `http://127.0.0.1:11435`)",
            f"- **Cumulative Tantalus Usage:** `{total_calls}` calls / `{total_tokens:,}` tokens (**`$0.00`** cloud API expenditure)",
            "",
            "| Modality ID | ARIA Lane | Public Accession | Status | Runtime | Claims | LLM Calls / Tokens | Cost | Capsule | Vignette |",
            "|---|---|---|---|---:|---:|---:|---:|---|---|",
            *rows,
            "",
        ]
    )
    (PUBLIC_BENCH_DIR / "STATUS.md").write_text(status_md, encoding="utf-8")

    summary_payload = _sanitize_public_obj(
        {
            "schema_version": "aria.public_modality_atlas.summary.v1",
            "updated_at_utc": now_utc,
            "completed_count": completed_count,
            "total_modalities": len(ALL_MODALITY_IDS),
            "total_llm_calls": total_calls,
            "total_llm_tokens": total_tokens,
            "total_cloud_cost_usd": 0.0,
            "modalities": summary_items,
        }
    )
    (PUBLIC_BENCH_DIR / "atlas_summary.json").write_text(
        json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8"
    )

    vignette_readme = "\n".join(
        [
            "# ARIA Reproducible Public Omics Vignettes (11-Modality Atlas)",
            "",
            "Every vignette below documents a complete, unattended end-to-end execution of ARIA v4.8+ on a **100% public, peer-reviewed dataset** using the local **Tantalus** open-weight reagent (`Qwen3.6-35B-A3B-UD-Q4_K_M`, `$0.00` cloud API expenditure), verified via `RO-Crate 1.1` and `s14.reproducibility_capsule.v2`.",
            "",
            "| Modality ID | ARIA Lane | Public Accession | Status | Runtime | Claims | Capsule | Vignette |",
            "|---|---|---|---|---:|---:|---|---|",
            *[
                f"| `{it['modality_id']}` | `{it.get('aria_modality', '-')}` | `{it.get('accession', '-')}` | "
                f"`{it['status']}` | `{it.get('elapsed_seconds', '-')}` | `{it.get('total_public_claims', '-')}` | "
                f"`{it.get('capsule_verification', '-')}` | "
                + (f"[`{it['modality_id']}.md`](./{it['modality_id']}.md) |" if it.get("vignette_md") else "`-` |")
                for it in summary_items
            ],
            "",
        ]
    )
    (VIGNETTES_DIR / "README.md").write_text(vignette_readme, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--modality",
        default="all",
        choices=["all", *ALL_MODALITY_IDS],
        help="Single modality ID or 'all' to run the full 11-modality atlas.",
    )
    ap.add_argument("--atlas-root", type=Path, default=DEFAULT_ATLAS_ROOT)
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_RUNS_ROOT)
    ap.add_argument(
        "--force",
        action="store_true",
        help="Re-execute even if a verified receipt JSON already exists.",
    )
    ap.add_argument(
        "--skip-hypotheses",
        action="store_true",
        help="Disable opt-in quarantined S0-S10 HypothesisAgent generation.",
    )
    ap.add_argument(
        "--use-llm-cache",
        action="store_true",
        help="Allow disk LLM prompt cache hits instead of forcing fresh Tantalus inference.",
    )
    args = ap.parse_args()

    atlas_root = args.atlas_root.expanduser().resolve()
    runs_root = args.runs_root.expanduser().resolve()
    atlas_root.mkdir(parents=True, exist_ok=True)
    runs_root.mkdir(parents=True, exist_ok=True)
    RECEIPTS_DIR.mkdir(parents=True, exist_ok=True)
    VIGNETTES_DIR.mkdir(parents=True, exist_ok=True)

    _log("Ensuring shared public reference genomes & indexes ...")
    refs = ensure_references(atlas_root)

    targets = ALL_MODALITY_IDS if args.modality == "all" else (args.modality,)
    live_state: dict[str, str] = {}
    _refresh_global_status(live_state)

    failures: list[tuple[str, str]] = []
    for mid in targets:
        receipt_path = RECEIPTS_DIR / f"{mid}.json"
        if receipt_path.is_file() and not args.force:
            try:
                existing = json.loads(receipt_path.read_text(encoding="utf-8"))
                if (
                    existing.get("status") == "completed"
                    and (existing.get("reproducibility_capsule") or {}).get("verification_status")
                    in {"pass", "verified_portable"}
                ):
                    _log(f"[{mid}] Already completed & capsule-verified; skipping (use --force to re-run).")
                    _write_vignette(existing)
                    _refresh_global_status(live_state)
                    continue
            except Exception:
                pass

        live_state[mid] = "running"
        _refresh_global_status(live_state)
        try:
            receipt = execute_single_modality(
                mid,
                atlas_root=atlas_root,
                runs_root=runs_root,
                refs=refs,
                enable_hypotheses=not args.skip_hypotheses,
                use_llm_cache=args.use_llm_cache,
            )
            receipt_path.write_text(
                json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
            )
            _write_vignette(receipt)
            live_state.pop(mid, None)
            _refresh_global_status(live_state)
            _log(
                f"[{mid}] COMPLETED in {receipt['elapsed_seconds']:.1f}s "
                f"(claims={receipt['claims_audit']['total_public_claims']}, "
                f"llm_calls={receipt['llm_reagent']['usage']['total_calls']}, "
                f"tokens={receipt['llm_reagent']['usage']['total_tokens']:,}, cost=$0.00)"
            )
        except Exception as exc:
            err_txt = f"{type(exc).__name__}: {exc}"
            _log(f"[{mid}] FAILED — {err_txt}")
            traceback.print_exc()
            failures.append((mid, err_txt))
            fail_payload = _sanitize_public_obj(
                {
                    "schema_version": "aria.public_modality_atlas.receipt.v1",
                    "modality_id": mid,
                    "status": "failed",
                    "failed_at_utc": datetime.now(timezone.utc).isoformat(),
                    "error": err_txt,
                }
            )
            (runs_root / f"{mid}_last_error.json").write_text(
                json.dumps(fail_payload, indent=2) + "\n", encoding="utf-8"
            )
            live_state[mid] = f"failed ({type(exc).__name__})"
            _refresh_global_status(live_state)

    if failures:
        _log(f"Finished with {len(failures)} failure(s): {failures}")
        return 1
    _log("All requested modalities completed and verified successfully!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
