"""ADR-076: Universal Multi-Platform scATAC-seq Engine & Production Promotion Suite.

Covers:
1. Pure-Python readiness, contract, manifest geometry, auto-RC orientation,
   ChromatinNarrator, and C3 RunLedger claim compilation (runs in every CI lane).
2. Live binary multi-platform end-to-end verification across ``aria-atacseq-env``,
   ``aria-chromatin-env``, and ``aria-rna-env`` (real ``chromap``, ``macs3``,
   ``snapatac2``, ``pyDESeq2``, ``chromatin_peak_annotation``,
   ``chromatin_peak_ora``, ``chromatin_merge_matrices``, ``chromatin_regulatory``,
   and ``NarrativeAgent`` C3 compilation):
   - Trajectory 1: 6-donor (3 ctrl vs 3 treat), multi-lane (L001+L002) raw FASTQ
                   cohort with reverse-complemented cell-barcode reads on half
                   the libraries -> auto-RC detection + A4 stage resume ->
                   100% air-gapped snapatac2 fragments-to-matrix -> LSI ->
                   donor-aware pyDESeq2 DA -> peak annotation & ORA -> C3 report!
   - Trajectory 2: Non-10x barcode geometries (``sci_atac`` & ``custom``
                   ``--read-format``) on live ``chromap`` + fail-closed
                   ``barcode_whitelist_mismatch`` rejection.
   - Trajectory 3: Multi-sample ``fragments.tsv.gz`` cohort routed through
                   ``ChromatinAgent`` with zero dropped samples.
   - Trajectory 4: Multi-file ``.h5ad`` cohort merged via genomic-overlap
                   consensus peaks (``chromatin_merge_matrices.py``, ADR-043) ->
                   donor-aware pyDESeq2 DA + Signac-calibrated ``peak2gene`` null.
"""

from __future__ import annotations

import gzip
import json
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


def _scatac_binaries_available() -> bool:
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    return bool(
        envs.get("atacseq") and envs.get("chromatin") and envs.get("rna")
    )


def _generate_barcodes(n: int, length: int = 16, seed: int = 101) -> list[str]:
    """Generate distinct non-homopolymer barcodes at >= 3 Hamming distance."""
    rng = random.Random(seed)
    alphabet = "ACGT"
    chosen: list[str] = []
    while len(chosen) < n:
        cand = "".join(rng.choice(alphabet) for _ in range(length))
        if any(b * 5 in cand for b in alphabet):
            continue
        rc = _revcomp(cand)
        if cand == rc:
            continue
        if all(
            sum(x != y for x, y in zip(cand, prev)) >= 4
            and sum(x != y for x, y in zip(rc, prev)) >= 4
            for prev in chosen
        ):
            chosen.append(cand)
    return chosen


# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Pure-Python Contracts, ScATACAuditAgent, Geometries & C3 Ledger
# ══════════════════════════════════════════════════════════════════════════════


def test_scatac_production_promotion_and_scatac_audit_governance(tmp_path):
    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.orchestrator_agent import (
        INTEGRATION_VALIDATION,
        MODALITY_VALIDATION,
    )
    from aria.utils.registry_integrity import check_registry_integrity
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    assert MODALITY_VALIDATION["scATAC"] == {
        "level": "production",
        "dispatch_enabled": True,
        "promotion_adr": "ADR-076",
    }
    assert MODALITY_VALIDATION["bulk_ATAC"]["level"] == "production"
    assert MODALITY_VALIDATION["bulk_ATAC"]["promotion_adr"] == "ADR-077"
    assert INTEGRATION_VALIDATION["production_promoted"] is False
    assert check_registry_integrity() == []

    for script in (
        "aria/scripts/chromatin_scatac_align.py",
        "aria/scripts/chromatin_fragments_to_matrix.py",
        "aria/scripts/chromatin_merge_matrices.py",
        "aria/scripts/chromatin_qc.py",
        "aria/scripts/chromatin_peaks.py",
        "aria/scripts/chromatin_lsi_clustering.py",
        "aria/scripts/chromatin_diffacc.py",
        "aria/scripts/chromatin_peak_annotation.py",
        "aria/scripts/chromatin_peak_ora.py",
        "aria/scripts/chromatin_motifs.py",
        "aria/scripts/chromatin_figure_clusters.py",
    ):
        assert SCRIPT_CONTRACTS[script].validation_level == "production", script

    assert SCRIPT_CONTRACTS["aria/scripts/chromatin_footprint_tobias.py"].validation_level == "production"
    assert SCRIPT_CONTRACTS["aria/scripts/chromatin_peak_counts.py"].validation_level == "production"

    h5mu = tmp_path / "cohort.h5mu"
    h5mu.write_bytes(b"stub")

    # 1) Non-contrasted exploratory run -> green / auto
    m_exp = build_capability_matrix(
        {"genome": "hg38", "modalities": {"scATAC": [str(h5mu)]}},
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_exp["cards"]["scATAC"]["status"] == "green"
    assert m_exp["cards"]["scATAC"]["dispatch_policy"] == "auto"

    # 2) Donor-aware n=3 vs n=3 contrast -> green / auto
    m_green = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"scATAC": [str(h5mu)]},
            "design": {
                "main_factor": "condition",
                "replicate_col": "donor_id",
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
    assert m_green["cards"]["scATAC"]["status"] == "green"
    assert m_green["cards"]["scATAC"]["dispatch_policy"] == "auto"

    # 3) Low-power n=2 vs n=2 contrast -> yellow / requires_ack
    m_yellow = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"scATAC": [str(h5mu)]},
            "design": {
                "main_factor": "condition",
                "replicate_col": "donor_id",
                "groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]},
                "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
            },
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_yellow["cards"]["scATAC"]["status"] == "yellow"
    assert m_yellow["cards"]["scATAC"]["dispatch_policy"] == "requires_ack"

    # 4) Unreplicated n=1 vs n=1 contrast -> red / blocked
    m_red = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"scATAC": [str(h5mu)]},
            "design": {
                "main_factor": "condition",
                "replicate_col": "donor_id",
                "groups": {"ctrl": ["c1"], "treat": ["t1"]},
                "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
            },
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert m_red["cards"]["scATAC"]["status"] == "red"
    assert m_red["cards"]["scATAC"]["dispatch_policy"] == "blocked"


def test_scatac_manifest_multiplatform_multilane_and_auto_rc(tmp_path):
    from aria.scripts.chromatin_scatac_align import (
        detect_barcode_whitelist_orientation,
    )
    from aria.utils.scatac_fastq_manifest import (
        PLATFORM_ROLES,
        resolve_scatac_fastq_manifest,
    )

    wl_barcodes = _generate_barcodes(25, length=16, seed=7)
    wl = tmp_path / "whitelist.txt"
    wl.write_text("\n".join(wl_barcodes) + "\n", encoding="utf-8")

    r1_l1 = _write_fastq_gz(tmp_path / "s1_L001_R1.fastq.gz", [("q1", "ACGT" * 12, "I" * 48)])
    r1_l2 = _write_fastq_gz(tmp_path / "s1_L002_R1.fastq.gz", [("q2", "TGCA" * 12, "I" * 48)])
    r3_l1 = _write_fastq_gz(tmp_path / "s1_L001_R3.fastq.gz", [("q1", "ACGT" * 12, "I" * 48)])
    r3_l2 = _write_fastq_gz(tmp_path / "s1_L002_R3.fastq.gz", [("q2", "TGCA" * 12, "I" * 48)])
    bc_rc_l1 = _write_fastq_gz(
        tmp_path / "s1_L001_R2.fastq.gz",
        [(f"q{i}", _revcomp(wl_barcodes[i % 25]), "I" * 16) for i in range(30)],
    )
    bc_rc_l2 = _write_fastq_gz(
        tmp_path / "s1_L002_R2.fastq.gz",
        [(f"q{i+30}", _revcomp(wl_barcodes[i % 25]), "I" * 16) for i in range(30)],
    )

    for platform, read_fmt in (
        ("10x_atac", None),
        ("10x_multiome_atac", None),
        ("sci_atac", None),
        ("share_seq", None),
        ("surecell", None),
        ("custom", "bc:0:15,r1:0:-1,r2:0:-1"),
    ):
        fastqs = {"R1": [str(r1_l1), str(r1_l2)]}
        if "R2" in PLATFORM_ROLES[platform] or platform == "custom":
            fastqs["R2"] = [str(bc_rc_l1), str(bc_rc_l2)]
        if "R3" in PLATFORM_ROLES[platform] or platform == "custom":
            fastqs["R3"] = [str(r3_l1), str(r3_l2)]
        lib = {
            "library_id": f"lib_{platform}",
            "sample_id": "s1",
            "donor_id": "d1",
            "condition": "ctrl",
            "platform": platform,
            "barcode_orientation": "auto",
            "fastqs": fastqs,
            "barcode_whitelist": str(wl),
        }
        if read_fmt:
            lib["read_format"] = read_fmt
        res = resolve_scatac_fastq_manifest(
            {
                "scatac_fastq_manifest": {
                    "schema_version": "1",
                    "libraries": [lib],
                }
            },
            require_paths=True,
        )
        assert res["status"] == "valid", (platform, res)
        norm_lib = res["manifest"]["libraries"][0]
        assert norm_lib["platform"] == platform
        assert norm_lib["fastq_lanes"]["R1"] == [
            str(r1_l1.resolve()),
            str(r1_l2.resolve()),
        ]
        if platform == "custom":
            assert norm_lib["read_format"] == read_fmt

    # Automatic RC detection identifies reverse-complemented cell barcodes
    orient = detect_barcode_whitelist_orientation(
        [str(bc_rc_l1), str(bc_rc_l2)],
        str(wl),
    )
    assert orient["status"] == "sampled"
    assert orient["orientation"] == "reverse_complement"
    assert orient["rc_rate"] == pytest.approx(1.0)

    # Incompatible barcodes report whitelist_mismatch so caller fails closed
    bad_bc = _write_fastq_gz(
        tmp_path / "bad_R2.fastq.gz",
        [(f"b{i}", "CCCCCCCCCCCCCCCC", "I" * 16) for i in range(30)],
    )
    bad = detect_barcode_whitelist_orientation(
        [str(bad_bc)],
        str(wl),
    )
    assert bad["status"] == "whitelist_mismatch"


def test_scatac_preprocessing_narrative_and_c3_ledger_linkage(tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.agents.narrative.validators import validate_blocks
    from aria.agents.narrative_agent import NarrativeAgent

    findings = {
        "scATAC": {
            "status": "done",
            "findings": {
                "fastq_to_fragments": {
                    "status": "success",
                    "aligner": "chromap",
                    "n_libraries": 6,
                    "n_fragments": 7200,
                    "validation_level": "production",
                    "library_manifest": [
                        {"barcode_orientation": "forward"},
                        {"barcode_orientation": "reverse_complement"},
                    ],
                },
                "fragments_to_matrix": {
                    "status": "success",
                    "n_cells": 96,
                    "n_peaks": 40,
                    "peak_mode": "provided",
                    "validation_level": "production",
                },
                "qc": {
                    "status": "success",
                    "data_type": "scATAC",
                    "n_cells": 96,
                    "n_peaks": 40,
                    "n_fragments": 7200,
                },
                "lsi": {
                    "status": "success",
                    "n_cells": 96,
                    "n_peaks": 40,
                    "n_clusters": 2,
                    "n_components_kept": 29,
                },
                "differential_accessibility": {
                    "status": "success",
                    "groupby": "leiden",
                    "n_clusters": 2,
                    "padj_max": 0.05,
                    "lfc_min": 0.5,
                    "per_cluster": {
                        "ran": True,
                        "n_da_total": 8,
                        "n_da_by_cluster": {"0": 4, "1": 4},
                        "da_peaks_by_cluster": {},
                    },
                    "pseudobulk": {
                        "ran": True,
                        "comparisons": [
                            {
                                "test": "treat",
                                "reference": "ctrl",
                                "status": "success",
                                "n_tested": 40,
                                "n_sig": 8,
                                "n_up": 4,
                                "n_down": 4,
                                "top_peaks": [],
                            }
                        ],
                    },
                },
                "peak_annotation": {
                    "status": "success",
                    "ran": True,
                    "analysis": "peak_annotation",
                    "data_type": "scATAC",
                    "validation_level": "production",
                    "gtf": "synthetic.gtf",
                    "promoter_upstream": 2000,
                    "promoter_downstream": 500,
                    "feature_distribution_overall": {"Promoter": 8},
                    "comparisons": [
                        {
                            "test": "treat",
                            "reference": "ctrl",
                            "n_annotated": 8,
                            "n_unparsed": 0,
                            "feature_distribution": {"Promoter": 8},
                        }
                    ],
                },
            },
        }
    }

    narrator = ChromatinNarrator()
    blocks = validate_blocks(
        narrator.collect("chromatin_agent", {"status": "done", "findings": findings})
    )
    by_id = {b.id: b for b in blocks}
    assert "chromatin.scatac_preprocessing" in by_id
    assert by_id["chromatin.scatac_preprocessing"].analysis == "scatac_preprocessing"
    assert "chromatin.peak_annotation" in by_id
    assert "scATAC" in by_id["chromatin.peak_annotation"].title

    ledger = build_run_ledger(
        {}, {"chromatin_agent": {"status": "done", "findings": findings}}
    )
    by_node = {e["node_id"]: e for e in ledger["entries"]}
    assert by_node["ledger://chromatin/scatac_preprocessing"]["status"] == "ran"
    assert by_node["ledger://chromatin/peak_annotation"]["status"] == "ran"

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    report = agent._render_html_report(
        experiment_id="exp_scatac_c3",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "Universal scATAC multi-donor DA verification"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="scATAC chromap + snapatac2 + pyDESeq2",
        decisions=[],
        agent_results={"chromatin_agent": {"status": "done", "findings": findings}},
        report_dir=tmp_path / "report",
    )
    assert report.exists()
    methodology = json.loads((report.parent / "methodology.json").read_text())
    assert methodology["claim_compilation"]["n_withheld"] == 0
    assert methodology["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0
    claim_analyses = {c["analysis"] for c in methodology["claims"]}
    assert "scatac_preprocessing" in claim_analyses
    assert "peak_annotation" in claim_analyses


# ══════════════════════════════════════════════════════════════════════════════
# Part 2: Live Binary Multi-Platform E2E Suite
#         (chromap + macs3 + snapatac2 + pyDESeq2 + peak_annotation + peak_ora)
# ══════════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def mini_chromatin_reference(tmp_path_factory):
    """Build a non-repetitive 2-contig genome (chr1, chr2; 100 kb each) with 40
    well-spaced open chromatin loci and matching neutral GTF gene models."""
    if not _scatac_binaries_available():
        pytest.skip("aria-atacseq-env / aria-chromatin-env / aria-rna-env not installed")

    ref_dir = tmp_path_factory.mktemp("mini_scatac_ref")
    rng = random.Random(20261112)
    bases = "ACGT"
    chrom_len = 100_000
    contigs = {
        "chr1": "".join(rng.choice(bases) for _ in range(chrom_len)),
        "chr2": "".join(rng.choice(bases) for _ in range(chrom_len)),
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
                fh.write(seq[pos: pos + line_w] + "\n")
            fai_lines.append(f"{cname}\t{len(seq)}\t{offset}\t{line_w}\t{line_w + 1}")
            n_lines = (len(seq) + line_w - 1) // line_w
            offset += len(seq) + n_lines
    fai.write_text("\n".join(fai_lines) + "\n", encoding="utf-8")

    # 40 peak centers (20 on chr1, 20 on chr2) spaced every 4,000 bp from 6,000..82,000.
    # Peaks 0..3 on chr1 are 8x higher in treat; peaks 0..3 on chr2 are 8x higher in ctrl;
    # remaining 32 peaks are constant background loci.
    peaks = []
    gtf_lines = []
    idx = 0
    for cname in ("chr1", "chr2"):
        for k in range(20):
            idx += 1
            center = 6_000 + k * 4_000
            p_start = center - 250
            p_end = center + 250
            gid = f"LOC_{idx:04d}"
            regime = "bg"
            if cname == "chr1" and k < 4:
                regime = "treat_up"
            elif cname == "chr2" and k < 4:
                regime = "ctrl_up"
            peaks.append(
                {
                    "peak_idx": idx,
                    "chrom": cname,
                    "start": p_start,
                    "end": p_end,
                    "center": center,
                    "gene_id": gid,
                    "regime": regime,
                }
            )
            # Place TSS right at p_end + 100 so the peak falls cleanly in Promoter!
            tss = p_end + 100
            g_end = tss + 1_500
            gtf_lines.append(
                f'{cname}\tsynthetic\tgene\t{tss}\t{g_end}\t.\t+\t.\tgene_id "{gid}"; gene_name "{gid}";\n'
                f'{cname}\tsynthetic\texon\t{tss}\t{tss + 400}\t.\t+\t.\tgene_id "{gid}"; gene_name "{gid}";\n'
            )

    gtf = ref_dir / "genes.gtf"
    gtf.write_text("".join(gtf_lines), encoding="utf-8")

    return {
        "ref_dir": ref_dir,
        "fasta": str(fasta),
        "fai": str(fai),
        "gtf": str(gtf),
        "contigs": contigs,
        "peaks": peaks,
    }


def _emit_fragment_pair(contigs: dict[str, str], chrom: str, frag_start: int, frag_len: int = 210, read_len: int = 50) -> tuple[str, str]:
    seq = contigs[chrom]
    r1 = seq[frag_start: frag_start + read_len]
    r3 = _revcomp(seq[frag_start + frag_len - read_len: frag_start + frag_len])
    return r1, r3


def test_live_multidonor_multilane_autorc_airgapped_raw_fastq_scatac_e2e(
    tmp_path, mini_chromatin_reference
):
    """Full live ADR-076 E2E: 6 donors (3 ctrl vs 3 treat), 2 lanes/donor,
    mixed forward + reverse-complement barcode FASTQs, air-gapped snapatac2
    matrix construction, LSI, donor-aware pyDESeq2 DA, peak annotation/ORA, and
    NarrativeAgent C3 verification."""
    pytest.importorskip("litellm")
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils import privacy
    from aria.utils.privacy import EgressPolicy, use_egress_policy

    contigs = mini_chromatin_reference["contigs"]
    peaks = mini_chromatin_reference["peaks"]
    rng = random.Random(9001)

    # Shared 96-barcode 16-bp whitelist (16 distinct cells per donor across 6 donors)
    all_barcodes = _generate_barcodes(96, length=16, seed=55)
    wl_path = tmp_path / "10x_atac_whitelist.txt"
    wl_path.write_text("\n".join(all_barcodes) + "\n", encoding="utf-8")

    donors = [
        ("ctrl_1", "ctrl", False),
        ("ctrl_2", "ctrl", True),   # Reverse-complemented cell barcodes!
        ("ctrl_3", "ctrl", False),
        ("treat_1", "treat", True),  # Reverse-complemented cell barcodes!
        ("treat_2", "treat", False),
        ("treat_3", "treat", True),  # Reverse-complemented cell barcodes!
    ]

    libraries = []
    all_fastq_files = []
    raw_dir = tmp_path / "raw_cohort"

    for d_idx, (donor_id, cond, use_rc) in enumerate(donors):
        donor_cells = all_barcodes[d_idx * 16: (d_idx + 1) * 16]
        lane_records = {1: ([], [], []), 2: ([], [], [])}
        q_counter = 0

        for c_idx, cell_bc in enumerate(donor_cells):
            emitted_bc = _revcomp(cell_bc) if use_rc else cell_bc
            for pk in peaks:
                if pk["regime"] == "treat_up":
                    n_frags = rng.randint(7, 10) if cond == "treat" else rng.randint(0, 1)
                elif pk["regime"] == "ctrl_up":
                    n_frags = rng.randint(7, 10) if cond == "ctrl" else rng.randint(0, 1)
                else:
                    n_frags = rng.randint(2, 4)

                for f_i in range(n_frags):
                    q_counter += 1
                    lane = 1 if (q_counter % 2 == 1) else 2
                    # Distinct start offsets so chromap retains unique fragments per cell!
                    jitter = ((c_idx * 11 + f_i * 19 + d_idx * 7) % 240) - 120
                    frag_start = pk["center"] + jitter - 105
                    frag_len = 200 + ((f_i * 7 + c_idx * 3) % 35)
                    r1_seq, r3_seq = _emit_fragment_pair(
                        contigs, pk["chrom"], frag_start, frag_len=frag_len, read_len=50
                    )
                    qname = f"{donor_id}_q{q_counter}"
                    lane_records[lane][0].append((qname, r1_seq, "I" * 50))
                    lane_records[lane][1].append((qname, emitted_bc, "I" * 16))
                    lane_records[lane][2].append((qname, r3_seq, "I" * 50))

        r1_lanes, r2_lanes, r3_lanes = [], [], []
        for lane in (1, 2):
            p1 = _write_fastq_gz(raw_dir / f"{donor_id}_L00{lane}_R1.fastq.gz", lane_records[lane][0])
            p2 = _write_fastq_gz(raw_dir / f"{donor_id}_L00{lane}_R2.fastq.gz", lane_records[lane][1])
            p3 = _write_fastq_gz(raw_dir / f"{donor_id}_L00{lane}_R3.fastq.gz", lane_records[lane][2])
            r1_lanes.append(str(p1))
            r2_lanes.append(str(p2))
            r3_lanes.append(str(p3))
            all_fastq_files.extend([str(p1), str(p2), str(p3)])

        libraries.append(
            {
                "library_id": donor_id,
                "sample_id": donor_id,
                "donor_id": donor_id,
                "condition": cond,
                "platform": "10x_atac",
                "barcode_orientation": "auto",
                "fastqs": {"R1": r1_lanes, "R2": r2_lanes, "R3": r3_lanes},
                "barcode_whitelist": str(wl_path),
            }
        )

    mem = ARIAMemory(db_path=str(tmp_path / "mem_scatac.sqlite"))
    agent = ChromatinAgent(memory=mem, llm=None)

    exp_ctx = {
        "genome": "hg38",
        "organism": "Homo sapiens",
        "threads": 2,
        "min_fragments": 20,
        "genome_fasta": mini_chromatin_reference["fasta"],
        "genome_fai": mini_chromatin_reference["fai"],
        "gtf": mini_chromatin_reference["gtf"],
        "enable_footprinting": False,
        "modalities": {"scATAC": all_fastq_files},
        "scatac_fastq_manifest": {"schema_version": "1", "libraries": libraries},
        "design": {
            "main_factor": "condition",
            "replicate_col": "donor_id",
            "design_formula": "~ condition",
            "groups": {
                "ctrl": ["ctrl_1", "ctrl_2", "ctrl_3"],
                "treat": ["treat_1", "treat_2", "treat_3"],
            },
            "contrasts": [
                {"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}
            ],
        },
    }

    # Execute 100% air-gapped (egress_allowed() == False)!
    with use_egress_policy(EgressPolicy(air_gapped=True, reason="test_adr076")):
        assert privacy.egress_allowed() is False
        run_out = agent.run(
            "exp_scatac_adr076",
            {
                "exp_context": exp_ctx,
                "biological_intent": {
                    "summary": "Compare open chromatin accessibility between treat and ctrl donors"
                },
            },
        )

    assert run_out["status"] == "done", run_out
    sc_findings = run_out["findings"]["scATAC"]["findings"]

    # 1. Verify chromap multi-lane + auto-RC alignment & A4 stage manifest
    align = sc_findings["fastq_to_fragments"]
    assert align["status"] == "success", align
    assert align["n_libraries"] == 6
    assert align["n_fragments"] > 3_000
    orientations = {
        row["library_id"]: row["barcode_orientation"]
        for row in align["library_manifest"]
    }
    assert orientations == {
        "ctrl_1": "direct",
        "ctrl_2": "reverse_complement",
        "ctrl_3": "direct",
        "treat_1": "reverse_complement",
        "treat_2": "direct",
        "treat_3": "reverse_complement",
    }
    align_dir = Path(align["fragments_file"]).parent
    assert (align_dir / "ctrl_1" / "chromatin_scatac_align.stage_manifest.json").is_file()
    assert (align_dir / "ctrl_2" / "barcode_whitelist.rc.txt").is_file()

    # Verify A4 stage manifest cache hit on second alignment invocation
    second_align = agent._run_scatac_fastq_to_fragments(
        "exp_scatac_adr076", exp_ctx, all_fastq_files
    )
    assert second_align["status"] == "success"
    assert second_align.get("resumed") is True

    # 2. Verify offline air-gapped snapatac2 fragments->matrix bridge
    f2m = sc_findings["fragments_to_matrix"]
    assert f2m["status"] == "success", f2m
    assert f2m["n_cells"] == 96
    assert f2m["n_peaks"] >= 30
    assert (
        Path(f2m["output_path"]).parent
        / "chromatin_fragments_to_matrix.stage_manifest.json"
    ).is_file()

    # 3. Verify LSI clustering + donor-aware pyDESeq2 pseudobulk DA
    assert sc_findings["lsi"]["status"] == "success"
    da = sc_findings["differential_accessibility"]
    assert da["status"] == "success", da
    pb = da["pseudobulk"]
    assert pb["ran"] is True, pb
    comp = pb["comparisons"][0]
    assert comp["status"] == "success", comp
    assert comp["n_sig"] >= 6, comp
    assert comp["n_up"] >= 3 and comp["n_down"] >= 3, comp
    assert comp.get("full_results_csv") and Path(comp["full_results_csv"]).is_file()

    # 4. Verify live peak_annotation + peak_ora subprocesses in aria-rna-env
    ann = sc_findings["peak_annotation"]
    assert ann["status"] == "success", ann
    assert ann["data_type"] == "scATAC"
    assert ann["validation_level"] == "production"
    assert ann["comparisons"][0]["n_annotated"] == comp["n_sig"]
    assert Path(ann["comparisons"][0]["annotation_csv"]).name.startswith("scatac_peak_annotation_")

    ora = sc_findings["peak_ora"]
    assert ora["status"] in ("success", "skipped"), ora
    assert ora["data_type"] == "scATAC"
    assert ora["validation_level"] == "production"

    # 5. Verify NarrativeAgent compiles claims with 0 C3 RunLedger violations!
    narr = NarrativeAgent.__new__(NarrativeAgent)
    narr.reports_dir = tmp_path
    narr.memory = type("M", (), {"db_path": ":memory:"})()
    report = narr._render_html_report(
        experiment_id="exp_scatac_live_report",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "Multi-donor air-gapped raw scATAC DA"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="chromap + snapatac2 + pyDESeq2",
        decisions=[],
        agent_results={"chromatin_agent": run_out},
        report_dir=tmp_path / "live_report",
    )
    methodology = json.loads((report.parent / "methodology.json").read_text())
    assert methodology["claim_compilation"]["n_withheld"] == 0, methodology["claim_compilation"]
    assert methodology["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0
    claim_analyses = {c["analysis"] for c in methodology["claims"]}
    assert "scatac_preprocessing" in claim_analyses
    assert "peak_annotation" in claim_analyses


def test_live_non10x_custom_geometry_and_fail_closed_whitelist_mismatch(
    tmp_path, mini_chromatin_reference
):
    """Live chromap verification on sci_atac (18-bp barcode) and custom spacer
    geometry (4-bp adapter + 16-bp barcode), plus fail-closed mismatch."""
    from aria.utils.environment_manager import env_manager

    contigs = mini_chromatin_reference["contigs"]
    peaks = mini_chromatin_reference["peaks"]

    # 1) sci_atac (18-bp combinatorial barcode)
    sci_bcs = _generate_barcodes(12, length=18, seed=303)
    sci_wl = tmp_path / "sci_wl.txt"
    sci_wl.write_text("\n".join(sci_bcs) + "\n", encoding="utf-8")

    r1_recs, r2_recs, r3_recs = [], [], []
    for i in range(60):
        pk = peaks[i % 10]
        r1_s, r3_s = _emit_fragment_pair(
            contigs, pk["chrom"], pk["center"] - 100 + (i % 25), frag_len=210
        )
        r1_recs.append((f"sci_{i}", r1_s, "I" * 50))
        r2_recs.append((f"sci_{i}", sci_bcs[i % 12], "I" * 18))
        r3_recs.append((f"sci_{i}", r3_s, "I" * 50))

    sci_r1 = _write_fastq_gz(tmp_path / "sci" / "R1.fastq.gz", r1_recs)
    sci_r2 = _write_fastq_gz(tmp_path / "sci" / "R2.fastq.gz", r2_recs)
    sci_r3 = _write_fastq_gz(tmp_path / "sci" / "R3.fastq.gz", r3_recs)

    sci_res = env_manager.run_in_stack(
        stack="atacseq",
        script_path="aria/scripts/chromatin_scatac_align.py",
        params={
            "r1_fastq": str(sci_r1),
            "r3_fastq": str(sci_r3),
            "barcode_fastq": str(sci_r2),
            "barcode_whitelist": str(sci_wl),
            "platform": "sci_atac",
            "genome_fasta": mini_chromatin_reference["fasta"],
            "output_dir": str(tmp_path / "sci_out"),
            "threads": 2,
        },
    )
    assert sci_res["status"] == "success", sci_res
    assert sci_res["n_fragments"] >= 30

    # 2) custom spacer geometry: 4-bp fixed adapter ("GCTA") + 16-bp barcode -> bc:4:19
    cust_bcs = _generate_barcodes(10, length=16, seed=404)
    cust_wl = tmp_path / "cust_wl.txt"
    cust_wl.write_text("\n".join(cust_bcs) + "\n", encoding="utf-8")
    c_r2 = _write_fastq_gz(
        tmp_path / "cust" / "R2.fastq.gz",
        [(f"sci_{i}", "GCTA" + cust_bcs[i % 10], "I" * 20) for i in range(60)],
    )
    cust_res = env_manager.run_in_stack(
        stack="atacseq",
        script_path="aria/scripts/chromatin_scatac_align.py",
        params={
            "r1_fastq": str(sci_r1),
            "r3_fastq": str(sci_r3),
            "barcode_fastq": str(c_r2),
            "barcode_whitelist": str(cust_wl),
            "platform": "custom",
            "read_format": "bc:4:19,r1:0:-1,r2:0:-1",
            "genome_fasta": mini_chromatin_reference["fasta"],
            "output_dir": str(tmp_path / "cust_out"),
            "threads": 2,
        },
    )
    assert cust_res["status"] == "success", cust_res
    assert cust_res["n_fragments"] >= 30

    # 3) Fail-closed barcode_whitelist_mismatch when wrong whitelist is supplied
    mismatch = env_manager.run_in_stack(
        stack="atacseq",
        script_path="aria/scripts/chromatin_scatac_align.py",
        params={
            "r1_fastq": str(sci_r1),
            "r3_fastq": str(sci_r3),
            "barcode_fastq": str(sci_r2),
            "barcode_whitelist": str(cust_wl),
            "genome_fasta": mini_chromatin_reference["fasta"],
            "output_dir": str(tmp_path / "mismatch_out"),
            "threads": 2,
        },
    )
    assert mismatch["status"] == "skipped"
    assert mismatch["reason"] == "barcode_whitelist_mismatch"


def test_live_multisample_fragments_cohort_preserves_all_samples(
    tmp_path, mini_chromatin_reference
):
    """Verify multi-file fragments.tsv.gz inputs passed to ChromatinAgent union
    every sample onto one matrix without dropping non-first fragment files."""
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.memory.memory import ARIAMemory

    peaks = mini_chromatin_reference["peaks"]
    bcs = _generate_barcodes(20, length=16, seed=707)
    frag_paths = []

    for s_idx, sample_id in enumerate(("sampleA", "sampleB")):
        fpath = tmp_path / f"{sample_id}.fragments.tsv.gz"
        rows = []
        for pk in peaks[:25]:
            for c_i in range(10):
                bc = bcs[s_idx * 10 + c_i]
                for rep in range(3):
                    s = pk["start"] + rep * 25 + c_i * 5
                    rows.append((pk["chrom"], s, s + 120, bc, 1))
        rows.sort(key=lambda r: (r[0], r[1], r[2]))
        with gzip.open(fpath, "wt", encoding="utf-8") as fh:
            for chrom, s, e, bc, cnt in rows:
                fh.write(f"{chrom}\t{s}\t{e}\t{bc}\t{cnt}\n")
        frag_paths.append(str(fpath))

    mem = ARIAMemory(db_path=str(tmp_path / "mem_frags.sqlite"))
    agent = ChromatinAgent(memory=mem, llm=None)
    res = agent.run(
        "exp_multi_frags",
        {
            "exp_context": {
                "genome": "hg38",
                "min_fragments": 15,
                "genome_fasta": mini_chromatin_reference["fasta"],
                "genome_fai": mini_chromatin_reference["fai"],
                "enable_footprinting": False,
                "sample_ids": ["sampleA", "sampleB"],
                "sample_metadata": {
                    "sampleA": {"condition": "ctrl", "donor_id": "sampleA"},
                    "sampleB": {"condition": "treat", "donor_id": "sampleB"},
                },
                "modalities": {"scATAC": frag_paths},
            },
            "biological_intent": {"summary": "Exploratory multi-sample fragments QC & clustering"},
        },
    )
    assert res["status"] == "done", res
    f2m = res["findings"]["scATAC"]["findings"]["fragments_to_matrix"]
    assert f2m["status"] == "success", f2m
    assert f2m["n_cells"] == 20
    assert f2m["n_libraries"] == 2


def test_live_multifile_matrix_cohort_consensus_union_and_signac_peak2gene_e2e(
    tmp_path, mini_chromatin_reference
):
    """Verify multi-file .h5ad cohort union (jittered peak coordinates unified by
    genomic overlap per ADR-043) -> donor-aware pyDESeq2 DA + Signac-calibrated
    peak2gene null in standalone chromatin_regulatory.py."""
    import anndata as ad
    from scipy import sparse

    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.environment_manager import env_manager

    peaks = mini_chromatin_reference["peaks"]
    rng = np.random.default_rng(2026)
    matrix_files = []
    sample_ids = []
    sample_metadata = {}

    all_rna_X = []
    all_rna_obs = []

    donors = [
        ("ctrl_1", "ctrl", -30),
        ("ctrl_2", "ctrl", 0),
        ("ctrl_3", "ctrl", +25),
        ("treat_1", "treat", -20),
        ("treat_2", "treat", +15),
        ("treat_3", "treat", +35),
    ]

    for donor_id, cond, coord_shift in donors:
        n_cells = 15
        X = np.zeros((n_cells, len(peaks)), dtype=np.int32)
        rna_X = np.zeros((n_cells, len(peaks)), dtype=np.int32)
        peak_ids = []
        for j, pk in enumerate(peaks):
            # Jitter peak intervals per donor so string intersection would be 0%,
            # while genomic interval overlap is >85%!
            p_s = pk["start"] + coord_shift
            p_e = pk["end"] + coord_shift
            peak_ids.append(f"{pk['chrom']}:{p_s}-{p_e}")
            if pk["regime"] == "treat_up":
                lam = 9.0 if cond == "treat" else 1.0
            elif pk["regime"] == "ctrl_up":
                lam = 9.0 if cond == "ctrl" else 1.0
            else:
                lam = 3.0
            counts = rng.poisson(lam, size=n_cells).astype(np.int32)
            X[:, j] = counts
            # Paired RNA expression correlated with open chromatin for first 8 loci
            if j < 8:
                rna_X[:, j] = (counts * 2 + rng.poisson(1.0, size=n_cells)).astype(np.int32)
            else:
                rna_X[:, j] = rng.poisson(3.0, size=n_cells).astype(np.int32)

        obs_names = [f"cell_{i:02d}" for i in range(n_cells)]
        adata = ad.AnnData(
            X=sparse.csr_matrix(X),
            obs=pd.DataFrame(
                {"condition": cond, "donor_id": donor_id},
                index=obs_names,
            ),
            var=pd.DataFrame(index=peak_ids),
        )
        mpath = tmp_path / "donor_matrices" / f"{donor_id}.h5ad"
        mpath.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(mpath)
        matrix_files.append(str(mpath))
        sample_ids.append(donor_id)
        sample_metadata[donor_id] = {"condition": cond, "donor_id": donor_id}

        all_rna_X.append(rna_X)
        all_rna_obs.extend([f"{donor_id}#{c}" for c in obs_names])

    # Write paired RNA h5ad matching the merged `{sample_id}#{barcode}` obs_names
    rna_adata = ad.AnnData(
        X=sparse.csr_matrix(np.vstack(all_rna_X)),
        obs=pd.DataFrame(index=all_rna_obs),
        var=pd.DataFrame(index=[pk["gene_id"] for pk in peaks]),
    )
    rna_path = tmp_path / "paired_rna.h5ad"
    rna_adata.write_h5ad(rna_path)

    mem = ARIAMemory(db_path=str(tmp_path / "mem_cohort.sqlite"))
    agent = ChromatinAgent(memory=mem, llm=None)
    res = agent.run(
        "exp_cohort_union",
        {
            "exp_context": {
                "genome": "hg38",
                "genome_fasta": mini_chromatin_reference["fasta"],
                "genome_fai": mini_chromatin_reference["fai"],
                "gtf": mini_chromatin_reference["gtf"],
                "enable_footprinting": False,
                "rna_data_path": str(rna_path),
                "sample_ids": sample_ids,
                "sample_metadata": sample_metadata,
                "modalities": {"scATAC": matrix_files},
                "design": {
                    "main_factor": "condition",
                    "replicate_col": "donor_id",
                    "design_formula": "~ condition",
                    "groups": {
                        "ctrl": ["ctrl_1", "ctrl_2", "ctrl_3"],
                        "treat": ["treat_1", "treat_2", "treat_3"],
                    },
                    "contrasts": [
                        {"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}
                    ],
                },
            },
            "biological_intent": {
                "summary": "Compare consensus chromatin accessibility and cis peak-to-gene links"
            },
        },
    )
    assert res["status"] == "done", res
    sc_f = res["findings"]["scATAC"]["findings"]

    # 1. Verify genomic-overlap consensus peak union preserved all 40 loci across 90 cells!
    union = sc_f["cohort_matrix_union"]
    assert union["status"] == "success", union
    assert union["n_cells"] == 90
    assert union["n_consensus_peaks"] == 40

    # 2. Verify pyDESeq2 recovered the 8 differential consensus peaks!
    pb = sc_f["differential_accessibility"]["pseudobulk"]
    assert pb["ran"] is True
    assert pb["comparisons"][0]["n_sig"] == 8

    # 3. Verify standalone chromatin_regulatory.py ran Signac-calibrated null peak2gene!
    reg = env_manager.run_in_stack(
        stack="chromatin",
        script_path="aria/scripts/chromatin_regulatory.py",
        params={
            "data_path": union["output_path"],
            "rna_data_path": str(rna_path),
            "gtf_path": mini_chromatin_reference["gtf"],
            "genome_fasta": mini_chromatin_reference["fasta"],
            "genome_fai": mini_chromatin_reference["fai"],
            "min_shared_cells": 50,
            "min_peak_gene_corr": 0.25,
            "output_dir": str(tmp_path / "reg_out"),
        },
    )
    assert reg["status"] == "success", reg
    p2g = reg["peak_to_gene"]
    assert p2g["ran"] is True, p2g
    assert p2g["null_calibrated"] is True, p2g
    assert p2g["n_links"] >= 4
    top0 = p2g["top_links"][0]
    assert top0["z_score"] is not None
    assert top0["pvalue"] is not None
    assert top0["fdr"] is not None
