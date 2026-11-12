"""ADR-077: Universal Bulk ATAC-seq Engine & Production Promotion E2E Suite.

Covers:
1. Pure-Python governance, ``BulkATACAuditAgent`` green/yellow/red gating,
   single-cell barcoded FASTQ misrouting rejection, ``ScriptContract`` &
   registry integrity, ``RunLedger`` ``ledger://chromatin/bulk_atac_preprocessing``
   provenance, and strict ``W-CLAIM`` ``ChromatinNarrator`` compilation (runs in
   every CI lane).
2. Live multi-stack binary end-to-end verification across ``aria-rnaseq-env``
   (``fastp``), ``aria-atacseq-env`` (``bwa-mem2`` + ``samtools``),
   ``aria-chromatin-env`` (``pysam`` + ``bedtools`` + ``macs3`` + ``snapatac2``),
   and ``aria-rna-env`` (``pyDESeq2`` + ``matplotlib``):
   - Trajectory 1: 6-replicate paired-end raw FASTQ cohort (3 ctrl vs 3 treat)
                   with nucleosome-free & mononucleosomal insert sizes,
                   mitochondrial ``chrM`` contaminant reads, and blacklist
                   artifact reads -> ``fastp`` -> ``atac_align`` (chrM stripped)
                   -> local-GTF ENCODE QC (FRiP + TSSe + dup_rate + nucleosome
                   banding) -> ``macs3`` consensus peak calling + blacklist
                   filtering -> raw + CPM peak counts -> replicate-gated
                   ``pyDESeq2`` DA (driven by ``design.contrasts`` without
                   ``intent["comparison"]``) -> genomic peak annotation + gene
                   ORA + motifs -> dual PNG+SVG publication figures
                   (``bulk_atac_pca``, ``bulk_atac_qc_summary``, ``da_volcano``,
                   ``fragment_size``, ``sample_correlation``,
                   ``peak_annotation``) -> ``NarrativeAgent`` ``C3`` HTML report
                   with ``0`` withheld blocks!
   - Trajectory 2: Single-end raw FASTQ -> ``atac_align`` SE mode (``paired=False``,
                   no ``fixmate`` / ``-f 2``) -> ``chromatin_peaks`` automatic
                   SE Tn5-shift detection (``-f BAM --nomodel --shift -100
                   --extsize 200``) -> ``A4`` content-addressed stage-manifest
                   cache hit on second execution!
   - Trajectory 3: Direct pre-computed ``peak_id x sample`` TSV count matrix
                   intake in ``ChromatinAgent._run_bulk_atac`` -> ``pyDESeq2``
                   differential accessibility recovery!
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def _write_fastq_gz(path: Path, records: list[tuple[str, str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for qname, seq, qual in records:
            fh.write(f"@{qname}\n{seq}\n+\n{qual}\n")
    return path


def _bulk_atac_binaries_available() -> bool:
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    return bool(
        envs.get("rnaseq")
        and envs.get("atacseq")
        and envs.get("chromatin")
        and envs.get("rna")
    )


# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Pure-Python Contracts, BulkATACAuditAgent, Misrouting & C3 Ledger
# ══════════════════════════════════════════════════════════════════════════════


def test_bulk_atac_production_promotion_and_audit_governance(tmp_path):
    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.orchestrator_agent import (
        INTEGRATION_VALIDATION,
        MODALITY_VALIDATION,
    )
    from aria.utils.registry_integrity import check_registry_integrity
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    assert MODALITY_VALIDATION["bulk_ATAC"] == {
        "level": "production",
        "dispatch_enabled": True,
        "promotion_adr": "ADR-077",
    }
    assert INTEGRATION_VALIDATION["production_promoted"] is False
    assert check_registry_integrity() == []

    for script in (
        "aria/scripts/rna_fastq_qc.py",
        "aria/scripts/atac_align.py",
        "aria/scripts/chromatin_qc.py",
        "aria/scripts/chromatin_peaks.py",
        "aria/scripts/chromatin_peak_counts.py",
        "aria/scripts/chromatin_bulk_diffacc.py",
        "aria/scripts/chromatin_peak_annotation.py",
        "aria/scripts/chromatin_peak_ora.py",
        "aria/scripts/chromatin_motifs.py",
        "aria/scripts/chromatin_figure_bulk.py",
    ):
        assert SCRIPT_CONTRACTS[script].validation_level == "production", script

    assert (
        SCRIPT_CONTRACTS["aria/scripts/chromatin_footprint_tobias.py"].validation_level
        == "beta"
    )

    bam = tmp_path / "cohort_ctrl_1.bam"
    bam.write_bytes(b"BAM\x01")

    # 1) Non-contrasted exploratory run -> green / auto
    m_exp = build_capability_matrix(
        {"genome": "hg38", "modalities": {"bulk_ATAC": [str(bam)]}},
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_exp["cards"]["bulk_ATAC"]["status"] == "green"
    assert m_exp["cards"]["bulk_ATAC"]["dispatch_policy"] == "auto"
    assert m_exp["cards"]["bulk_ATAC"]["validation_level"] == "production"

    # 2) Replicated n=3 vs n=3 contrast -> green / auto
    m_green = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(bam)]},
            "design": {
                "main_factor": "condition",
                "replicate_col": "replicate",
                "groups": {
                    "ctrl": ["c1", "c2", "c3"],
                    "treat": ["t1", "t2", "t3"],
                },
                "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
            },
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_green["cards"]["bulk_ATAC"]["status"] == "green"
    assert m_green["cards"]["bulk_ATAC"]["dispatch_policy"] == "auto"

    # 3) Low-power n=2 vs n=2 contrast -> yellow / requires_ack
    m_yellow = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(bam)]},
            "design": {
                "main_factor": "condition",
                "replicate_col": "replicate",
                "groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]},
                "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
            },
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_yellow["cards"]["bulk_ATAC"]["status"] == "yellow"
    assert m_yellow["cards"]["bulk_ATAC"]["dispatch_policy"] == "requires_ack"
    assert m_yellow["cards"]["bulk_ATAC"]["findings"]

    # 4) Unreplicated n=1 vs n=1 contrast -> red / blocked
    m_red = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(bam)]},
            "design": {
                "main_factor": "condition",
                "replicate_col": "replicate",
                "groups": {"ctrl": ["c1"], "treat": ["t1"]},
                "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
            },
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_red["cards"]["bulk_ATAC"]["status"] == "red"
    assert m_red["cards"]["bulk_ATAC"]["dispatch_policy"] == "blocked"

    # 5) Batch-confounded n=3 vs n=3 contrast -> red / blocked
    m_conf = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(bam)]},
            "sample_metadata": {
                "c1": {"condition": "ctrl", "batch": "B1"},
                "c2": {"condition": "ctrl", "batch": "B1"},
                "c3": {"condition": "ctrl", "batch": "B1"},
                "t1": {"condition": "treat", "batch": "B2"},
                "t2": {"condition": "treat", "batch": "B2"},
                "t3": {"condition": "treat", "batch": "B2"},
            },
            "design": {
                "main_factor": "condition",
                "batch_factors": ["batch"],
                "groups": {
                    "ctrl": ["c1", "c2", "c3"],
                    "treat": ["t1", "t2", "t3"],
                },
                "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
            },
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_conf["cards"]["bulk_ATAC"]["status"] == "red"
    assert m_conf["cards"]["bulk_ATAC"]["dispatch_policy"] == "blocked"

    # 6) Barcoded 3-read scATAC FASTQ misrouted to bulk_ATAC -> red / blocked
    sc_r1 = _write_fastq_gz(
        tmp_path / "sc1_S1_L001_R1_001.fastq.gz", [("q1", "ACGT" * 13, "I" * 52)]
    )
    sc_r2 = _write_fastq_gz(
        tmp_path / "sc1_S1_L001_R2_001.fastq.gz", [("q1", "ACGTACGTACGTACGT", "I" * 16)]
    )
    sc_r3 = _write_fastq_gz(
        tmp_path / "sc1_S1_L001_R3_001.fastq.gz", [("q1", "TGCA" * 13, "I" * 52)]
    )
    m_misrouted = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(sc_r1), str(sc_r2), str(sc_r3)]},
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_misrouted["cards"]["bulk_ATAC"]["status"] == "red"
    assert m_misrouted["cards"]["bulk_ATAC"]["dispatch_policy"] == "blocked"
    assert any(
        "single_cell" in str(f).lower() or "scatac" in str(f).lower()
        for f in m_misrouted["cards"]["bulk_ATAC"]["findings"]
    )


def test_bulk_atac_c3_run_ledger_and_w_claim_strict_narrative(tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.narrative.evidence_verifier import verify_block_claim_support
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.memory.memory import ARIAMemory

    counts_tsv = tmp_path / "bulk_atac_peak_counts.tsv"
    cpm_tsv = tmp_path / "bulk_atac_peak_cpm.tsv"
    da_csv = tmp_path / "da_treat_vs_ctrl.csv"
    counts_tsv.write_text("peak_id\tc1\tc2\tc3\tt1\tt2\tt3\nchr1:100-600\t8\t9\t7\t80\t85\t78\n")
    cpm_tsv.write_text("peak_id\tc1\tc2\tc3\tt1\tt2\tt3\nchr1:100-600\t100.0\t100.0\t100.0\t900.0\t900.0\t900.0\n")
    da_csv.write_text(
        "peak_id,base_mean,log2_fc,lfc_se,stat,pvalue,padj\n"
        "chr1:100-600,44.5,3.32,0.31,10.7,1e-18,1e-15\n"
    )
    bam_1 = tmp_path / "c1.filtered.bam"
    bam_1.write_bytes(b"BAM\x01")
    fq_1 = tmp_path / "c1_R1.fastq.gz"
    fq_1.write_bytes(b"\x1f\x8b")

    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "bulk_ATAC": {
                    "status": "done",
                    "findings": {
                        "preprocessing": {
                            "status": "success",
                            "ran": True,
                            "data_type": "bulk_ATAC",
                            "validation_level": "production",
                            "aligner": "bwa-mem2",
                            "read_layout": "paired",
                            "n_samples_aligned": 6,
                            "fastq_files": [str(fq_1)],
                            "bam_files": [
                                {"name": "c1", "sample_id": "c1", "status": "success", "bam": str(bam_1)}
                            ],
                        },
                        "qc": {
                            "status": "success",
                            "data_type": "bulk_ATAC",
                            "validation_level": "production",
                            "n_samples": 6,
                            "mapped_reads": 120000,
                            "frip": 0.42,
                            "tss_enrichment": 8.4,
                            "dup_rate": 0.08,
                            "mito_fraction": 0.0,
                            "qc_complete": True,
                            "pass_qc": True,
                        },
                        "peaks": {
                            "status": "success",
                            "data_type": "bulk_ATAC",
                            "validation_level": "production",
                            "n_peaks": 240,
                            "n_consensus_peaks": 210,
                            "frip": 0.42,
                            "peaks_path": str(tmp_path / "pooled.narrowPeak"),
                            "consensus_peaks_path": str(tmp_path / "consensus.bed"),
                        },
                        "peak_counts": {
                            "status": "success",
                            "data_type": "bulk_ATAC",
                            "validation_level": "production",
                            "n_peaks": 210,
                            "n_samples": 6,
                            "counting_method": "bedtools coverage -counts",
                            "counts_matrix_path": str(counts_tsv),
                            "counts_cpm_path": str(cpm_tsv),
                            "sample_metadata_path": str(tmp_path / "meta.tsv"),
                        },
                        "differential_accessibility": {
                            "status": "success",
                            "ran": True,
                            "data_type": "bulk_ATAC",
                            "validation_level": "production",
                            "method": "DESeq2",
                            "design_formula": "~ condition",
                            "n_peaks_tested": 210,
                            "comparisons": [
                                {
                                    "status": "success",
                                    "comparison": "treat_vs_ctrl",
                                    "test": "treat",
                                    "reference": "ctrl",
                                    "n_test": 3,
                                    "n_reference": 3,
                                    "n_tested": 210,
                                    "n_sig": 12,
                                    "n_up": 8,
                                    "n_down": 4,
                                    "output_csv": str(da_csv),
                                    "top_up": [
                                        {"peak_id": "chr1:100-600", "log2_fc": 3.32, "padj": 1e-15}
                                    ],
                                    "top_down": [],
                                }
                            ],
                        },
                    },
                }
            },
        }
    }
    exp_ctx = {
        "genome": "hg38",
        "modalities": {"bulk_ATAC": [str(fq_1)]},
        "design": {
            "main_factor": "condition",
            "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
        },
    }

    ledger = build_run_ledger(exp_ctx, agent_results)
    pre_entry = next(
        e for e in ledger["entries"] if e["analysis"] == "bulk_atac_preprocessing"
    )
    assert pre_entry["node_id"] == "ledger://chromatin/bulk_atac_preprocessing"
    assert pre_entry["status"] == "executed"
    assert pre_entry["provenance"]["script"] == "aria/scripts/atac_align.py"
    assert str(bam_1) in pre_entry["provenance"]["artifacts"]

    narrator = ChromatinNarrator()
    blocks = narrator.collect("chromatin_agent", agent_results["chromatin_agent"], exp_ctx)
    block_ids = [b.id for b in blocks]
    assert "chromatin.bulk_atac_preprocessing" in block_ids

    for blk in blocks:
        ver = verify_block_claim_support(blk, strict=True)
        assert ver["status"] == "supported", (blk.id, ver)
        assert blk.metadata.get("validation_level") == "production", (blk.id, blk.metadata)

    tables = narrator.tables("chromatin_agent", agent_results["chromatin_agent"])
    table_paths = {t.get("path") for t in tables if isinstance(t, dict)}
    assert str(counts_tsv) in table_paths
    assert str(cpm_tsv) in table_paths
    assert str(da_csv) in table_paths

    import json

    mem = ARIAMemory(db_path=str(tmp_path / "mem_narr.sqlite"))
    na = NarrativeAgent(memory=mem, llm=None)
    report_path = na._render_html_report(
        experiment_id="exp_bulk_atac_c3",
        exp_ctx=exp_ctx,
        intent={"summary": "Compare open chromatin between treat and ctrl"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Bulk ATAC workflow.",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report_out",
    )
    assert Path(report_path).is_file()
    methodology = json.loads((Path(report_path).parent / "methodology.json").read_text())
    assert methodology["claim_compilation"]["n_withheld"] == 0, methodology["claim_compilation"]
    assert methodology["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0


# ══════════════════════════════════════════════════════════════════════════════
# Part 2: Live Multi-Stack Binary Verification (fastp + bwa-mem2 + samtools +
#         pysam + bedtools + macs3 + pyDESeq2 + matplotlib)
# ══════════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def mini_bulk_atac_reference(tmp_path_factory):
    """Build a deterministic 3-contig (chr1, chr2, chrM) FASTA + .fai + GTF +
    blacklist BED with real human p53-target gene symbols at treat_up promoters."""
    if not _bulk_atac_binaries_available():
        pytest.skip(
            "Live bulk ATAC binaries (rnaseq/atacseq/chromatin/rna stacks) unavailable"
        )

    ref_dir = tmp_path_factory.mktemp("mini_bulk_atac_ref")
    rng = random.Random(20261112)
    alphabet = "ACGT"
    contigs = {
        "chr1": "".join(rng.choice(alphabet) for _ in range(100_000)),
        "chr2": "".join(rng.choice(alphabet) for _ in range(100_000)),
        "chrM": "".join(rng.choice(alphabet) for _ in range(16_000)),
    }

    fasta = ref_dir / "genome.fa"
    fai = ref_dir / "genome.fa.fai"
    offset = 0
    fai_lines = []
    with open(fasta, "w", encoding="utf-8") as fh:
        for cname, seq in contigs.items():
            header = f">{cname}\n"
            fh.write(header)
            offset += len(header)
            line_w = 80
            for pos in range(0, len(seq), line_w):
                fh.write(seq[pos : pos + line_w] + "\n")
            fai_lines.append(f"{cname}\t{len(seq)}\t{offset}\t{line_w}\t{line_w + 1}")
            n_lines = (len(seq) + line_w - 1) // line_w
            offset += len(seq) + n_lines
    fai.write_text("\n".join(fai_lines) + "\n", encoding="utf-8")

    p53_genes = ["TP53", "CDKN1A", "MDM2", "BAX", "BBC3", "GADD45A", "FAS", "DDB2"]
    ctrl_genes = ["ACTB", "GAPDH", "RPLP0", "TBP", "HPRT1", "B2M", "PPIA", "EEF1A1"]
    bg_genes = [
        "RPL13A", "RPS18", "PGK1", "LDHA", "ALDOA", "ENO1", "PKM", "TPI1",
        "YWHAZ", "SDHA", "UBC", "HMBS", "GUSB", "TFRC", "NONO", "VCP",
        "CANX", "CALR", "HSPA8", "HSP90AB1", "ATP5F1B", "COX4I1", "NDUFA1", "UQCRC1",
    ]

    peaks = []
    gtf_lines = []
    bg_idx = 0
    peak_idx = 0
    for cname in ("chr1", "chr2"):
        for k in range(20):
            peak_idx += 1
            center = 6_000 + k * 4_000
            p_start = center - 200
            p_end = center + 200
            if cname == "chr1" and k < 8:
                regime = "treat_up"
                gname = p53_genes[k]
            elif cname == "chr2" and k < 8:
                regime = "ctrl_up"
                gname = ctrl_genes[k]
            else:
                regime = "bg"
                gname = bg_genes[bg_idx % len(bg_genes)]
                bg_idx += 1
            peaks.append(
                {
                    "peak_idx": peak_idx,
                    "chrom": cname,
                    "start": p_start,
                    "end": p_end,
                    "center": center,
                    "gene_name": gname,
                    "regime": regime,
                }
            )
            # Place gene TSS right at `center` so peak reads fall in [-50, +50] of TSS
            # and genomic peak annotation classifies them as Promoter.
            tss_1based = center + 1
            g_end = tss_1based + 1_500
            gtf_lines.append(
                f'{cname}\tsynthetic\tgene\t{tss_1based}\t{g_end}\t.\t+\t.\tgene_id "ENSG_{peak_idx:04d}"; gene_name "{gname}";\n'
                f'{cname}\tsynthetic\texon\t{tss_1based}\t{tss_1based + 400}\t.\t+\t.\tgene_id "ENSG_{peak_idx:04d}"; gene_name "{gname}";\n'
            )

    gtf = ref_dir / "annotation.gtf"
    gtf.write_text("".join(gtf_lines), encoding="utf-8")

    # Blacklist locus at chr2:90000-91000
    blacklist = ref_dir / "blacklist.bed"
    blacklist.write_text("chr2\t90000\t91000\tartefact_locus\n", encoding="utf-8")

    return {
        "ref_dir": ref_dir,
        "fasta": str(fasta),
        "fai": str(fai),
        "gtf": str(gtf),
        "blacklist": str(blacklist),
        "contigs": contigs,
        "peaks": peaks,
    }


def _emit_pe_pair(
    contigs: dict[str, str],
    chrom: str,
    frag_start: int,
    frag_len: int,
    read_len: int = 50,
) -> tuple[str, str]:
    seq = contigs[chrom]
    r1 = seq[frag_start : frag_start + read_len]
    r2 = _revcomp(seq[frag_start + frag_len - read_len : frag_start + frag_len])
    return r1, r2


def test_live_paired_end_fastq_to_bam_qc_macs3_da_ora_figures_c3_e2e(
    tmp_path, mini_bulk_atac_reference
):
    """Full live ADR-077 E2E: 6 paired-end FASTQ samples -> fastp -> bwa-mem2
    (chrM stripped) -> local-GTF ENCODE QC -> MACS3 consensus peaks + blacklist
    filtering -> bedtools peak counts + CPM -> pyDESeq2 DA -> annotation + ORA
    -> dual PNG+SVG publication figures -> NarrativeAgent C3 HTML report."""
    pytest.importorskip("litellm")
    from aria.agents import _narrative_chromatin
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.privacy import EgressPolicy, use_egress_policy

    contigs = mini_bulk_atac_reference["contigs"]
    peaks = mini_bulk_atac_reference["peaks"]
    rng = random.Random(77001)

    samples = [
        ("ctrl_1", "ctrl"),
        ("ctrl_2", "ctrl"),
        ("ctrl_3", "ctrl"),
        ("treat_1", "treat"),
        ("treat_2", "treat"),
        ("treat_3", "treat"),
    ]
    raw_dir = tmp_path / "raw_fastq"
    fastq_files: list[str] = []
    sample_metadata: dict[str, dict[str, str]] = {}

    for s_idx, (sid, cond) in enumerate(samples):
        r1_recs: list[tuple[str, str, str]] = []
        r2_recs: list[tuple[str, str, str]] = []
        q_id = 0

        for pk in peaks:
            if pk["regime"] == "treat_up":
                n_frags = rng.randint(48, 58) if cond == "treat" else rng.randint(5, 8)
            elif pk["regime"] == "ctrl_up":
                n_frags = rng.randint(48, 58) if cond == "ctrl" else rng.randint(5, 8)
            else:
                n_frags = rng.randint(18, 24)

            for f_i in range(n_frags):
                q_id += 1
                # Alternate NFR (~82 bp) and mononucleosomal (~205 bp) fragments
                # with unique start/end offsets so samtools markdup -r retains them!
                if f_i % 2 == 0:
                    frag_len = 76 + (f_i % 20)  # 76..95 bp (<100 NFR)
                    offset = ((f_i * 7 + s_idx * 5) % 70) - 35
                    frag_start = pk["center"] + offset
                else:
                    frag_len = 190 + (f_i % 45)  # 190..234 bp (mononucleosomal)
                    offset = ((f_i * 11 + s_idx * 3) % 80) - 40
                    frag_start = pk["center"] + offset - 90
                r1_s, r2_s = _emit_pe_pair(contigs, pk["chrom"], frag_start, frag_len, 50)
                r1_recs.append((f"{sid}_q{q_id}", r1_s, "I" * 50))
                r2_recs.append((f"{sid}_q{q_id}", r2_s, "I" * 50))

        # Sparse TSS distal-flank reads in [-1980, -1920] so TSSe background flank > 0
        for flank_k in range(6):
            q_id += 1
            flank_start = peaks[flank_k]["center"] - 1965 + flank_k * 5 + s_idx
            r1_s, r2_s = _emit_pe_pair(contigs, peaks[flank_k]["chrom"], flank_start, 88, 50)
            r1_recs.append((f"{sid}_flank{q_id}", r1_s, "I" * 50))
            r2_recs.append((f"{sid}_flank{q_id}", r2_s, "I" * 50))

        # Blacklist artifact locus at chr2:90200-90600 (to prove blacklist_bed strips it)
        for bl_i in range(45):
            q_id += 1
            bl_start = 90_300 + (bl_i * 3) + s_idx
            r1_s, r2_s = _emit_pe_pair(contigs, "chr2", bl_start, 195 + (bl_i % 25), 50)
            r1_recs.append((f"{sid}_bl{q_id}", r1_s, "I" * 50))
            r2_recs.append((f"{sid}_bl{q_id}", r2_s, "I" * 50))

        # Mitochondrial contaminant reads on chrM (to prove atac_align strips chrM)
        for mt_i in range(25):
            q_id += 1
            mt_start = 1_000 + mt_i * 120 + s_idx * 7
            r1_s, r2_s = _emit_pe_pair(contigs, "chrM", mt_start, 180, 50)
            r1_recs.append((f"{sid}_mt{q_id}", r1_s, "I" * 50))
            r2_recs.append((f"{sid}_mt{q_id}", r2_s, "I" * 50))

        p1 = _write_fastq_gz(raw_dir / f"{sid}_R1.fastq.gz", r1_recs)
        p2 = _write_fastq_gz(raw_dir / f"{sid}_R2.fastq.gz", r2_recs)
        fastq_files.extend([str(p1), str(p2)])
        sample_metadata[sid] = {"condition": cond, "replicate": sid}

    mem = ARIAMemory(db_path=str(tmp_path / "mem_bulk_atac.sqlite"))
    agent = ChromatinAgent(memory=mem, llm=None)

    exp_ctx = {
        "genome": "sacCer3",  # compact effective genome size (1.2e7) while Homo sapiens organism drives human GMT ORA
        "organism": "Homo sapiens",
        "threads": 2,
        "genome_fasta": mini_bulk_atac_reference["fasta"],
        "genome_fai": mini_bulk_atac_reference["fai"],
        "gtf": mini_bulk_atac_reference["gtf"],
        "blacklist_bed": mini_bulk_atac_reference["blacklist"],
        "modalities": {"bulk_ATAC": fastq_files},
        "sample_metadata": sample_metadata,
        "design": {
            "main_factor": "condition",
            "replicate_col": "replicate",
            "design_formula": "~ condition",
            "groups": {
                "ctrl": ["ctrl_1", "ctrl_2", "ctrl_3"],
                "treat": ["treat_1", "treat_2", "treat_3"],
            },
            "contrasts": [
                {"numerator": "treat", "denominator": "ctrl", "name": "treat_vs_ctrl"}
            ],
        },
    }

    # Execute 100% air-gapped — notice biological_intent intentionally omits `comparison`
    # so `design.contrasts` alone must trigger peak_counts -> DA -> annotation -> ORA!
    with use_egress_policy(EgressPolicy(air_gapped=True, reason="test_adr077")):
        run_out = agent.run(
            "exp_bulk_atac_adr077",
            {
                "exp_context": exp_ctx,
                "biological_intent": {
                    "summary": "Evaluate chromatin accessibility remodeling in treat vs ctrl"
                },
            },
        )

    assert run_out["status"] == "done", run_out
    bk = run_out["findings"]["bulk_ATAC"]["findings"]

    # 1. Preprocessing: fastp + bwa-mem2 + chrM stripping + A4 stage manifests
    pre = bk["preprocessing"]
    assert pre["status"] == "success", pre
    assert pre["validation_level"] == "production"
    assert pre["data_type"] == "bulk_ATAC"
    assert pre["read_layout"] == "paired"
    assert pre["n_samples_aligned"] == 6
    assert all(b["pct_mito"] > 0 for b in pre["bam_files"])

    # 2. Local ENCODE QC: FRiP merged post-peaks, TSSe from local GTF, mito_fraction == 0.0
    qc = bk["qc"]
    assert qc["status"] == "success", qc
    assert qc["validation_level"] == "production"
    assert qc["n_samples"] == 6
    assert qc["mito_fraction"] == 0.0
    assert qc["frip"] is not None and qc["frip"] > 0.20, qc
    assert qc["tss_enrichment"] is not None and qc["tss_enrichment"] > 5.0, qc
    assert qc["qc_complete"] is True, qc
    assert qc["pass_qc"] is True, qc
    assert qc["fragment_sizes"]["nucleosome_banding_detected"] is True, qc["fragment_sizes"]

    # 3. MACS3 peak calling + consensus reproducibility + blacklist filtering
    pks = bk["peaks"]
    assert pks["status"] == "success", pks
    assert pks["validation_level"] == "production"
    assert pks["n_peaks"] >= 20, pks
    assert pks["n_blacklist_filtered"] >= 1, pks
    assert pks["consensus_peaks_path"] and Path(pks["consensus_peaks_path"]).is_file()

    # 4. Peak counts + CPM matrix (sample IDs cleanly stripped from *.filtered.bam!)
    pc = bk["peak_counts"]
    assert pc["status"] == "success", pc
    assert pc["validation_level"] == "production"
    assert pc["n_samples"] == 6
    assert Path(pc["counts_matrix_path"]).is_file()
    assert Path(pc["counts_cpm_path"]).is_file()
    cpm_df = pd.read_csv(pc["counts_cpm_path"], sep="\t")
    assert set(cpm_df.columns[1:]) == {
        "ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"
    }

    # 5. Replicate-gated pyDESeq2 differential accessibility
    da = bk["differential_accessibility"]
    assert da["status"] == "success" and da["ran"] is True, da
    assert da["validation_level"] == "production"
    comp = da["comparisons"][0]
    assert comp["status"] == "success", comp
    assert comp["n_up"] >= 4 and comp["n_down"] >= 4, comp
    da_df = pd.read_csv(comp["full_results_csv"])
    peak_col = "peak" if "peak" in da_df.columns else da_df.columns[0]
    lfc_raw_col = "log2FoldChange_raw" if "log2FoldChange_raw" in da_df.columns else "log2FoldChange"
    sig_mask = (
        da_df["significant"].fillna(False).astype(bool)
        if "significant" in da_df.columns
        else (da_df["padj"] < 0.05)
    )
    sig_up_chr1 = da_df[
        sig_mask
        & (da_df[lfc_raw_col] > 1.0)
        & da_df[peak_col].astype(str).str.startswith("chr1:")
    ]
    sig_dn_chr2 = da_df[
        sig_mask
        & (da_df[lfc_raw_col] < -1.0)
        & da_df[peak_col].astype(str).str.startswith("chr2:")
    ]
    assert len(sig_up_chr1) >= 4
    assert len(sig_dn_chr2) >= 4

    # 6. Genomic peak annotation + functional ORA
    ann = bk["peak_annotation"]
    assert ann["status"] == "success" and ann["ran"] is True, ann
    assert ann["validation_level"] == "production"
    assert ann["feature_distribution_overall"].get("Promoter", 0) >= 8

    ora = bk["peak_ora"]
    assert ora["validation_level"] == "production"
    if ora.get("ran"):
        assert ora["status"] == "success"

    # 7. Dual PNG + SVG publication figures (including bulk_atac_pca & bulk_atac_qc_summary)
    fig_dir = tmp_path / "figures_out"
    _narrative_chromatin.generate_figures(
        bk, None, fig_dir, env_manager=agent.env
    )
    figs = bk["figures"]
    for fig_key in (
        "bulk_atac_pca",
        "bulk_atac_qc_summary",
        "fragment_size",
        "bulk_atac_sample_correlation",
        "bulk_atac_peak_annotation",
    ):
        assert fig_key in figs, (fig_key, list(figs))
        png_p = Path(figs[fig_key])
        assert png_p.is_file() and png_p.with_suffix(".svg").is_file(), png_p
    assert any(k.startswith("da_volcano_") for k in figs), list(figs)

    # 8. C3 RunLedger + NarrativeAgent strict compilation (0 withheld blocks)
    import json

    agent_results = {"chromatin_agent": run_out}
    ledger = build_run_ledger(exp_ctx, agent_results)
    pre_entry = next(
        e for e in ledger["entries"] if e["analysis"] == "bulk_atac_preprocessing"
    )
    assert pre_entry["status"] == "executed", pre_entry

    na = NarrativeAgent(memory=mem, llm=None)
    html_path = na._render_html_report(
        experiment_id="exp_bulk_atac_adr077",
        exp_ctx=exp_ctx,
        intent={"summary": "Evaluate chromatin accessibility remodeling in treat vs ctrl"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Bulk ATAC workflow.",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report_live",
    )
    assert Path(html_path).is_file()
    methodology_live = json.loads((Path(html_path).parent / "methodology.json").read_text())
    assert methodology_live["claim_compilation"]["n_withheld"] == 0, methodology_live["claim_compilation"]
    assert methodology_live["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0


def test_live_single_end_fastq_tn5_shift_and_a4_stage_manifest_resume(
    tmp_path, mini_bulk_atac_reference
):
    """Live single-end ATAC FASTQ -> atac_align (paired=False) -> chromatin_peaks
    auto-detecting BAM_SE_TN5_SHIFT (--nomodel --shift -100 --extsize 200) ->
    A4 content-addressed stage-manifest cache reuse on rerun."""
    from aria.utils.environment_manager import env_manager

    contigs = mini_bulk_atac_reference["contigs"]
    peaks = mini_bulk_atac_reference["peaks"][:10]

    se_recs: list[tuple[str, str, str]] = []
    q_id = 0
    for pk in peaks:
        for f_i in range(40):
            q_id += 1
            pos = pk["center"] - 80 + (f_i * 4)
            seq = contigs[pk["chrom"]][pos : pos + 50]
            if f_i % 2 == 1:
                seq = _revcomp(seq)
            se_recs.append((f"se_q{q_id}", seq, "I" * 50))

    se_fq = _write_fastq_gz(tmp_path / "se_sample_R1.fastq.gz", se_recs)
    align_dir = tmp_path / "se_aligned"

    align_1 = env_manager.run_in_stack(
        stack="atacseq",
        script_path="aria/scripts/atac_align.py",
        params={
            "samples": [{"name": "se_sample", "r1": str(se_fq), "paired": False}],
            "genome_fasta": mini_bulk_atac_reference["fasta"],
            "output_dir": str(align_dir),
            "threads": 2,
        },
    )
    assert align_1["status"] == "success", align_1
    assert align_1["validation_level"] == "production"
    assert align_1["read_layout"] == "single"
    assert align_1["filters"]["samtools_f"] == 0
    assert align_1["bam_files"][0]["reused_from_stage_manifest"] is False
    se_bam = align_1["bam_files"][0]["bam"]
    assert Path(se_bam).is_file()

    # Second run -> A4 stage manifest cache hit!
    align_2 = env_manager.run_in_stack(
        stack="atacseq",
        script_path="aria/scripts/atac_align.py",
        params={
            "samples": [{"name": "se_sample", "r1": str(se_fq), "paired": False}],
            "genome_fasta": mini_bulk_atac_reference["fasta"],
            "output_dir": str(align_dir),
            "threads": 2,
        },
    )
    assert align_2["status"] == "success", align_2
    assert align_2["bam_files"][0]["reused_from_stage_manifest"] is True

    peaks_dir = tmp_path / "se_peaks"
    pk_1 = env_manager.run_in_stack(
        stack="chromatin",
        script_path="aria/scripts/chromatin_peaks.py",
        params={
            "data_type": "bulk_ATAC",
            "files": [se_bam],
            "genome": "sacCer3",
            "macs3_params": {"format": "BAMPE", "qvalue": 0.05},
            "output_dir": str(peaks_dir),
        },
    )
    assert pk_1["status"] == "success", pk_1
    assert pk_1["effective_format"] == "BAM_SE_TN5_SHIFT"
    assert "--shift -100" in pk_1["macs3_cmd"]
    assert pk_1["n_peaks"] >= 5
    assert pk_1["reused_from_stage_manifest"] is False

    # Second peak-calling run -> A4 stage manifest cache hit!
    pk_2 = env_manager.run_in_stack(
        stack="chromatin",
        script_path="aria/scripts/chromatin_peaks.py",
        params={
            "data_type": "bulk_ATAC",
            "files": [se_bam],
            "genome": "sacCer3",
            "macs3_params": {"format": "BAMPE", "qvalue": 0.05},
            "output_dir": str(peaks_dir),
        },
    )
    assert pk_2["status"] == "success", pk_2
    assert pk_2["reused_from_stage_manifest"] is True


def test_live_precomputed_peak_counts_matrix_intake_to_pydeseq2_da(
    tmp_path, mini_bulk_atac_reference
):
    """Direct pre-computed peak_id x sample TSV count matrix routed through
    ChromatinAgent._run_bulk_atac -> pyDESeq2 DA recovery without BAMs."""
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.memory.memory import ARIAMemory

    rng = np.random.default_rng(42077)
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    rows = []
    for k in range(120):
        pid = f"chr1:{1000 + k * 600}-{1400 + k * 600}"
        if k < 15:
            ctrl_c = rng.poisson(18, size=3)
            treat_c = rng.poisson(160, size=3)
        elif k < 30:
            ctrl_c = rng.poisson(160, size=3)
            treat_c = rng.poisson(18, size=3)
        else:
            ctrl_c = rng.poisson(55, size=3)
            treat_c = rng.poisson(55, size=3)
        rows.append([pid, *ctrl_c.tolist(), *treat_c.tolist()])

    mat_tsv = tmp_path / "precomputed_peak_counts.tsv"
    pd.DataFrame(rows, columns=["peak_id", *samples]).to_csv(
        mat_tsv, sep="\t", index=False
    )

    mem = ARIAMemory(db_path=str(tmp_path / "mem_precomp.sqlite"))
    agent = ChromatinAgent(memory=mem, llm=None)
    run_out = agent.run(
        "exp_bulk_atac_precomputed",
        {
            "exp_context": {
                "genome": "hg38",
                "gtf": mini_bulk_atac_reference["gtf"],
                "modalities": {"bulk_ATAC": [str(mat_tsv)]},
                "sample_metadata": {
                    s: {"condition": "ctrl" if s.startswith("ctrl") else "treat", "replicate": s}
                    for s in samples
                },
                "design": {
                    "main_factor": "condition",
                    "replicate_col": "replicate",
                    "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
                },
            },
            "biological_intent": {"summary": "Test precomputed peak matrix DA"},
        },
    )
    assert run_out["status"] == "done", run_out
    bk = run_out["findings"]["bulk_ATAC"]["findings"]
    assert bk["peak_counts"]["source"] == "precomputed_matrix"
    assert bk["peak_counts"]["validation_level"] == "production"
    da = bk["differential_accessibility"]
    assert da["status"] == "success" and da["ran"] is True, da
    assert da["validation_level"] == "production"
    comp = da["comparisons"][0]
    assert comp["n_up"] >= 10 and comp["n_down"] >= 10, comp
