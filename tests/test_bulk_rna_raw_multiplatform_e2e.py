"""ADR-075: Universal Raw FASTQ Multi-Platform RNA Engine & bulk_RNA_raw Production Suite.

Covers:
1. Pure-Python governance, contract, readiness, chemistry, CPM/TPM, narrative,
   C3 ledger linkage, and directional UMI graph unit checks (runs in every CI lane).
2. Live binary multi-platform end-to-end verification in ``aria-rnaseq-env`` +
   ``aria-rna-env`` (real ``fastp``, ``STAR``, ``STARsolo``, ``pysam``,
   ``samtools``, ``featureCounts``, and ``pyDESeq2``):
   - Family 1: ``full_length`` PE multi-lane (L001+L002), reverse-stranded (-s 2)
               + A4 ``counts_matrix.stage.json`` resume & param-invalidation
   - Family 2a: ``tag_3prime`` SE poly(A)-tailed 3' tags (-s 1) -> BulkRNAAgent E2E
                (fastp polyX trim -> STAR -> featureCounts -> pyDESeq2; CPM written,
                TPM suppressed)
   - Family 2b: ``tag_5prime`` SE 5'-TSS CAGE-style tags (-s 1)
   - Family 3: ``umi_bulk`` SE inline-UMI with 4x PCR amplification + 1-bp UMI
               sequencing errors -> exact molecule count recovery via rna_umi_dedup
   - Family 4: ``pooled_plate_bulk`` (BRB-seq 6-well multiplexed R1 barcode+UMI /
               R2 cDNA) -> STARsolo CB_UMI_Simple 1MM_Directional_UMItools ->
               BulkRNAAgent pyDESeq2 E2E
   - Fail-closed ``AmbiguousStrandedness`` guard on a real BAM with 30% sense reads
"""

from __future__ import annotations

from collections import Counter
import gzip
import json
import random
from pathlib import Path

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


def _rnaseq_binaries_available() -> bool:
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    return bool(envs.get("rnaseq") and envs.get("rna"))


# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Pure-Python Contracts, Readiness, Chemistry, Transforms & Narrative
# ══════════════════════════════════════════════════════════════════════════════


def test_bulk_rna_raw_production_promotion_and_fail_closed_contracts(tmp_path):
    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION
    from aria.utils.registry_integrity import check_registry_integrity
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    assert MODALITY_VALIDATION["bulk_RNA_raw"]["level"] == "production"
    assert MODALITY_VALIDATION["bulk_RNA_raw"]["dispatch_enabled"] is True

    for script in (
        "aria/scripts/rna_fastq_qc.py",
        "aria/scripts/rna_align.py",
        "aria/scripts/rna_umi_dedup.py",
        "aria/scripts/rna_quantify.py",
        "aria/scripts/rna_starsolo_bulk.py",
    ):
        assert SCRIPT_CONTRACTS[script].validation_level == "production", script

    assert check_registry_integrity() == []

    # 1) Valid n=3 vs n=3 bulk_RNA_raw dispatches automatically (green / auto)
    valid_fastqs = []
    for grp in ("ctrl", "treat"):
        for rep in (1, 2, 3):
            fq = _write_fastq_gz(
                tmp_path / f"{grp}_{rep}_R1.fastq.gz",
                [(f"r{i}", "ACGT" * 20, "I" * 80) for i in range(15)],
            )
            valid_fastqs.append(str(fq))

    green_matrix = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_RNA_raw": valid_fastqs},
            "design": {
                "groups": {
                    "ctrl": ["ctrl_1", "ctrl_2", "ctrl_3"],
                    "treat": ["treat_1", "treat_2", "treat_3"],
                }
            },
        },
        modality_validation=MODALITY_VALIDATION,
    )
    assert green_matrix["status"] == "green"
    assert green_matrix["cards"]["bulk_RNA_raw"]["dispatch_policy"] == "auto"

    # 2) Under-replicated n=1 bulk_RNA_raw blocks fail-closed
    red_n1 = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_RNA_raw": valid_fastqs[:2]},
            "design": {"groups": {"ctrl": ["ctrl_1"], "treat": ["treat_1"]}},
        },
        modality_validation=MODALITY_VALIDATION,
    )
    assert red_n1["status"] == "red"
    assert "bulk_RNA_raw" in red_n1["dispatch"]["blocked"]

    # 3) Unmanifested asymmetric barcode+UMI geometry (R1=28bp, R2=90bp) blocks
    asym_r1 = _write_fastq_gz(
        tmp_path / "asym" / "pool_R1.fastq.gz",
        [(f"r{i}", "A" * 14 + "C" * 14, "I" * 28) for i in range(20)],
    )
    asym_r2 = _write_fastq_gz(
        tmp_path / "asym" / "pool_R2.fastq.gz",
        [(f"r{i}", "ACGT" * 22 + "AC", "I" * 90) for i in range(20)],
    )
    red_asym = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_RNA_raw": [str(asym_r1), str(asym_r2)]},
            "design": {
                "groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]}
            },
        },
        modality_validation=MODALITY_VALIDATION,
    )
    assert red_asym["status"] == "red"
    checks_asym = {f["check"] for f in red_asym["findings"]}
    assert "assay_contract_asymmetric_barcode_umi_fastq" in checks_asym

    # 4) Long-read FASTQ (>=600 bp) blocks fail-closed
    lr_fq = _write_fastq_gz(
        tmp_path / "lr" / "nanopore_R1.fastq.gz",
        [(f"r{i}", "ACGT" * 200, "I" * 800) for i in range(15)],
    )
    red_lr = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"bulk_RNA_raw": [str(lr_fq)]},
            "design": {
                "groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]}
            },
        },
        modality_validation=MODALITY_VALIDATION,
    )
    assert red_lr["status"] == "red"
    assert any(
        f["check"] == "assay_contract_long_read_rna_unsupported"
        for f in red_lr["findings"]
    )

    # 5) V(D)J TCR/BCR repertoire FASTQ blocks fail-closed
    vdj_fq = _write_fastq_gz(
        tmp_path / "vdj" / "donor1_vdj_t_R1.fastq.gz",
        [(f"r{i}", "ACGT" * 25, "I" * 100) for i in range(15)],
    )
    red_vdj = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"scRNA": [str(vdj_fq)]},
        },
        modality_validation=MODALITY_VALIDATION,
    )
    assert red_vdj["status"] == "red"
    assert any(
        f["check"] == "assay_contract_vdj_repertoire_unsupported"
        for f in red_vdj["findings"]
    )


def test_scrna_kb_chemistry_universal_3p_5p_and_plate_resolution(tmp_path):
    from aria.utils.raw_ingestion import build_kb_count_command, resolve_kb_chemistry

    f1 = tmp_path / "s_R1.fastq.gz"
    f2 = tmp_path / "s_R2.fastq.gz"
    idx = tmp_path / "idx.idx"
    t2g = tmp_path / "t2g.txt"
    for p in (f1, f2, idx, t2g):
        p.write_bytes(b"stub")

    assert resolve_kb_chemistry("10xv3") == "10xv3"
    assert resolve_kb_chemistry("10X_5P") == "10XV2"
    assert resolve_kb_chemistry("10x_5p_pe") == "0,0,16:0,16,26:0,39,0,1,0,0"
    assert resolve_kb_chemistry("smartseq3") == "smartseq3"
    assert resolve_kb_chemistry("smart_seq3") == "SMARTSEQ3"
    assert resolve_kb_chemistry("bd_rhapsody") == "BDWTA"
    assert resolve_kb_chemistry("celseq2") == "celseq2"
    assert resolve_kb_chemistry("cel_seq2") == "CELSEQ2"
    assert resolve_kb_chemistry("splitseq") == "SPLIT-SEQ"
    assert resolve_kb_chemistry("visium") == "visium"
    assert resolve_kb_chemistry("10x_visium") == "VISIUM"
    assert resolve_kb_chemistry("0,0,12:0,12,20:1,0,0") == "0,0,12:0,12,20:1,0,0"

    cmd_5p_pe = build_kb_count_command(
        fastq_files=[f1, f2],
        index_path=idx,
        t2g_path=t2g,
        output_dir=tmp_path / "out",
        chemistry="10x_5p_pe",
    )
    x_idx = cmd_5p_pe.index("-x")
    assert cmd_5p_pe[x_idx + 1] == "0,0,16:0,16,26:0,39,0,1,0,0"

    with pytest.raises(ValueError, match="Unsupported scRNA kb chemistry"):
        resolve_kb_chemistry("totally_unknown_platform_xyz")


def test_directional_umi_graph_collapse_unit():
    from aria.scripts.rna_umi_dedup import _directional_winning_umis

    # 1 molecule amplified to 12 copies + two 1-bp sequencing error UMIs (counts 2 and 1)
    # plus a second genuine molecule at Hamming distance 1 with high count (count 10:
    # 12 < 2*10 - 1 = 19, so AAAAACGT must NOT be absorbed into AAAAACGA!).
    counts = Counter({
        "AAAAACGA": 12,
        "AAAAACGG": 2,   # 12 >= 2*2 - 1 = 3 -> absorbed
        "TAAAACGA": 1,   # 12 >= 2*1 - 1 = 1 -> absorbed
        "AAAAACGT": 10,  # 12 < 19 -> independent molecule retained!
        "CCCCCCCC": 4,   # distant UMI -> independent molecule retained!
    })
    winners = _directional_winning_umis(counts)
    assert winners == {"AAAAACGA", "AAAAACGT", "CCCCCCCC"}


def test_cpm_export_and_tpm_suppression_policy():
    from aria.scripts.rna_bulk.transforms import _compute_cpm
    from aria.utils.rna_fastq_protocol import is_end_tag_or_umi_bias_protocol

    assert is_end_tag_or_umi_bias_protocol("full_length") is False
    assert is_end_tag_or_umi_bias_protocol("tag_3prime") is True
    assert is_end_tag_or_umi_bias_protocol("quantseq") is True
    assert is_end_tag_or_umi_bias_protocol("tag_5prime") is True
    assert is_end_tag_or_umi_bias_protocol("cage") is True
    assert is_end_tag_or_umi_bias_protocol("pooled_plate_bulk") is True
    assert is_end_tag_or_umi_bias_protocol("brbseq") is True
    assert is_end_tag_or_umi_bias_protocol("umi_bulk") is True
    assert is_end_tag_or_umi_bias_protocol("umi_bulk", {"is_tag_assay": False}) is False

    df = pd.DataFrame(
        {"s1": [250, 750], "s2": [100, 100]},
        index=["G1", "G2"],
    )
    cpm = _compute_cpm(df, [])
    assert cpm.loc["G1", "s1"] == pytest.approx(250_000.0)
    assert cpm.loc["G2", "s1"] == pytest.approx(750_000.0)
    assert cpm.loc["G1", "s2"] == pytest.approx(500_000.0)


def test_bulk_preprocessing_narrative_and_c3_ledger_linkage(tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.agents.narrative.validators import validate_blocks
    from aria.agents.narrative_agent import NarrativeAgent

    findings = {
        "preprocessing": {
            "protocol": "tag_3prime",
            "suppress_tpm": True,
            "qc": {
                "n_samples": 6,
                "samples": [{"name": f"s{i}", "pct_passed": 96.5} for i in range(6)],
            },
            "alignment": {
                "n_aligned": 6,
                "bam_files": [
                    {"name": f"s{i}", "status": "success", "pct_unique": 91.2}
                    for i in range(6)
                ],
            },
            "quantification": {
                "n_genes": 12450,
                "n_samples": 6,
                "strand_used": 1,
                "read_layout": "single-end",
                "assignment_summary": {"pct_assigned": 84.3},
            },
        },
        "sample_qc": {"n_samples": 6, "size_ratio": 1.2, "outliers": []},
        "library_protocol": "tag_3prime",
        "tpm_status": {
            "computed": False,
            "reason": "suppressed_for_end_tag_or_umi_protocol_length_bias",
            "alternative": "counts_cpm.tsv",
            "library_protocol": "tag_3prime",
        },
        "contrasts": [
            {
                "name": "treat_vs_ctrl",
                "status": "success",
                "n_significant": 12,
                "n_upregulated": 8,
                "n_downregulated": 4,
            }
        ],
    }

    narrator = BulkRnaNarrator()
    blocks = validate_blocks(
        narrator.collect("bulk_rna_agent", {"status": "done", "findings": findings})
    )
    by_id = {b.id: b for b in blocks}
    assert "bulk.preprocessing" in by_id
    assert by_id["bulk.preprocessing"].analysis == "preprocessing"
    methods = narrator.methods("bulk_rna_agent", {"status": "done", "findings": findings})
    assert any("protocol=tag_3prime" in m for m in methods)
    assert any("TPM" in m and "counts_cpm.tsv" in m for m in methods)

    ledger = build_run_ledger({}, {"bulk_rna_agent": {"findings": findings}})
    by_node = {e["node_id"]: e for e in ledger["entries"]}
    assert by_node["ledger://bulk/preprocessing"]["status"] == "ran"
    assert ledger["n_divergences"] == 0

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    report = agent._render_html_report(
        experiment_id="exp_bulk_raw_c3",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "3prime tag raw FASTQ test"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="\n".join(methods),
        decisions=[],
        agent_results={"bulk_rna_agent": {"findings": findings}},
        report_dir=tmp_path / "report",
    )
    assert report.exists()
    methodology = json.loads((report.parent / "methodology.json").read_text())
    assert methodology["claim_compilation"]["n_withheld"] == 0
    assert methodology["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0
    claim_analyses = {c["analysis"] for c in methodology["claims"]}
    assert "preprocessing" in claim_analyses


# ══════════════════════════════════════════════════════════════════════════════
# Part 2: Live Binary Multi-Platform E2E Suite (fastp + STAR + STARsolo + pysam +
#         featureCounts + pyDESeq2)
# ══════════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def mini_spliced_reference(tmp_path_factory):
    """Build a 2-contig, 12-gene spliced reference FASTA + GTF + STAR index."""
    if not _rnaseq_binaries_available():
        pytest.skip("aria-rnaseq-env binaries (fastp/STAR/featureCounts/pysam) unavailable")

    from aria.utils.environment_manager import env_manager

    ref_dir = tmp_path_factory.mktemp("synth_rna_ref")
    rng = random.Random(20261112)
    alphabet = "ACGT"
    contigs = {
        "chrSynth1": [rng.choice(alphabet) for _ in range(30_000)],
        "chrSynth2": [rng.choice(alphabet) for _ in range(30_000)],
    }

    genes: list[dict] = []
    gtf_lines: list[str] = []

    for idx in range(12):
        gid = f"ENSG000009900{idx + 1:02d}"
        sym = f"SYNTH{idx + 1:02d}"
        chrom = "chrSynth1" if idx < 6 else "chrSynth2"
        local_i = idx % 6
        strand = "+" if idx % 2 == 0 else "-"
        g_start = 1500 + local_i * 4400  # 0-based
        # Two exons: [g_start, g_start+500) and [g_start+900, g_start+1600)
        ex1_s, ex1_e = g_start, g_start + 500
        in_s, in_e = g_start + 500, g_start + 900
        ex2_s, ex2_e = g_start + 900, g_start + 1600

        seq_arr = contigs[chrom]
        if strand == "+":
            seq_arr[in_s:in_s + 2] = list("GT")
            seq_arr[in_e - 2:in_e] = list("AG")
        else:
            # Reverse-complement of GT...AG on the + strand is CT...AC
            seq_arr[in_s:in_s + 2] = list("CT")
            seq_arr[in_e - 2:in_e] = list("AC")

        genes.append({
            "gene_id": gid,
            "symbol": sym,
            "chrom": chrom,
            "strand": strand,
            "exons_0based": [(ex1_s, ex1_e), (ex2_s, ex2_e)],
            "gene_span_0based": (ex1_s, ex2_e),
        })

        attrs = (
            f'gene_id "{gid}"; gene_name "{sym}"; '
            'gene_biotype "protein_coding";'
        )
        gtf_lines.append(
            f"{chrom}\tsynth\tgene\t{ex1_s + 1}\t{ex2_e}\t.\t{strand}\t.\t{attrs}"
        )
        t_attrs = f'{attrs} transcript_id "{gid}.t1";'
        gtf_lines.append(
            f"{chrom}\tsynth\ttranscript\t{ex1_s + 1}\t{ex2_e}\t.\t{strand}\t.\t{t_attrs}"
        )
        gtf_lines.append(
            f"{chrom}\tsynth\texon\t{ex1_s + 1}\t{ex1_e}\t.\t{strand}\t.\t{t_attrs} exon_number \"1\";"
        )
        gtf_lines.append(
            f"{chrom}\tsynth\texon\t{ex2_s + 1}\t{ex2_e}\t.\t{strand}\t.\t{t_attrs} exon_number \"2\";"
        )

    fasta_path = ref_dir / "genome.fa"
    with open(fasta_path, "w", encoding="utf-8") as fh:
        for chrom, arr in contigs.items():
            s = "".join(arr)
            fh.write(f">{chrom}\n")
            for pos in range(0, len(s), 80):
                fh.write(s[pos:pos + 80] + "\n")

    gtf_path = ref_dir / "genes.gtf"
    gtf_path.write_text("\n".join(gtf_lines) + "\n", encoding="utf-8")

    # Compute mature 5'->3' mRNA sequence (1200 bp) per gene
    contig_strs = {c: "".join(arr) for c, arr in contigs.items()}
    for g in genes:
        c_seq = contig_strs[g["chrom"]]
        (s1, e1), (s2, e2) = g["exons_0based"]
        spliced_plus = c_seq[s1:e1] + c_seq[s2:e2]
        g["mrna"] = spliced_plus if g["strand"] == "+" else _revcomp(spliced_plus)

    # Pre-build shared STAR index once via rna_align (also exercises dynamic SAindexNbases)
    star_index = ref_dir / "star_index"
    dummy_fq = _write_fastq_gz(
        ref_dir / "init" / "probe_R1.fastq.gz",
        [(f"p{i}", genes[0]["mrna"][100:175], "I" * 75) for i in range(10)],
    )
    init_res = env_manager.run_in_stack(
        "rnaseq",
        "aria/scripts/rna_align.py",
        {
            "samples": [{
                "name": "probe",
                "r1_trimmed": str(dummy_fq),
                "r2_trimmed": None,
                "paired": False,
                "read_layout": "single-end",
            }],
            "genome_dir": str(star_index),
            "genome_fasta": str(fasta_path),
            "gtf_file": str(gtf_path),
            "output_dir": str(ref_dir / "init_align"),
            "threads": 2,
            "two_pass": False,
        },
    )
    assert init_res["status"] == "success", init_res

    return {
        "fasta": str(fasta_path),
        "gtf": str(gtf_path),
        "star_index": str(star_index),
        "genes": genes,
    }


def test_live_family1_full_length_pe_multilane_autostrand_and_a4_resume(
    tmp_path, mini_spliced_reference
):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.environment_manager import env_manager

    genes = mini_spliced_reference["genes"]
    fq_dir = tmp_path / "raw_pe"
    fastq_files: list[str] = []

    # Generate reverse-stranded (-s 2) paired-end spliced reads across L001 + L002
    for sample_name in ("ctrl_1", "treat_1"):
        for lane in ("L001", "L002"):
            r1_recs = []
            r2_recs = []
            q_idx = 0
            for g_i, g in enumerate(genes[:6]):
                n_frags = 30 if (sample_name == "treat_1" and g_i < 2) else 12
                mrna = g["mrna"]
                for k in range(n_frags):
                    # Span the splice junction (400..620 crosses 500; 620..840 crosses 700)
                    start = 420 + (k % 5) * 8 if g["strand"] == "+" else 620 + (k % 5) * 8
                    frag = mrna[start:start + 210]
                    # TruSeq Stranded (-s 2): R1 is antisense, R2 is sense
                    r1_seq = _revcomp(frag)[:75]
                    r2_seq = frag[:75]
                    q_idx += 1
                    r1_recs.append((f"{sample_name}_{lane}_{q_idx}", r1_seq, "I" * 75))
                    r2_recs.append((f"{sample_name}_{lane}_{q_idx}", r2_seq, "I" * 75))
            p1 = _write_fastq_gz(fq_dir / f"{sample_name}_{lane}_R1_001.fastq.gz", r1_recs)
            p2 = _write_fastq_gz(fq_dir / f"{sample_name}_{lane}_R2_001.fastq.gz", r2_recs)
            fastq_files.extend([str(p1), str(p2)])

    mem = ARIAMemory(db_path=str(tmp_path / "mem.sqlite"))
    agent = BulkRNAAgent(memory=mem, llm=None)
    exp_ctx = {
        "n_cpus": 2,
        "genome_config": {
            "star_index": mini_spliced_reference["star_index"],
            "fasta": mini_spliced_reference["fasta"],
            "gtf": mini_spliced_reference["gtf"],
            "strand": "auto",
            "min_strand_reads": 40,
        },
    }
    counts_files, preproc = agent._run_preprocessing("exp_f1", fastq_files, exp_ctx, {})
    assert counts_files is not None, preproc
    quant = preproc["quantification"]
    assert quant["status"] == "success"
    assert quant["strand_used"] == 2  # Auto-detected reverse-stranded!
    assert quant["read_layout"] == "paired-end"
    assert quant["n_genes"] == 6
    assert quant["assignment_summary"]["pct_assigned"] > 80.0

    df = pd.read_csv(counts_files[0], sep="\t", index_col=0)
    assert list(df.columns) == ["ctrl_1", "treat_1"]
    # Each lane had 12 fragments -> 24 paired fragments in ctrl_1; 60 in treat_1 for GENE01
    assert df.loc[genes[0]["gene_id"], "ctrl_1"] == 24
    assert df.loc[genes[0]["gene_id"], "treat_1"] == 60

    # Verify A4 resume on rna_quantify, then verify invalidation when strand param changes
    counts_dir = Path(counts_files[0]).parent
    resumed = env_manager.run_in_stack(
        "rnaseq",
        "aria/scripts/rna_quantify.py",
        {
            "bam_files": preproc["alignment"]["bam_files"],
            "gtf_file": mini_spliced_reference["gtf"],
            "output_dir": str(counts_dir),
            "threads": 2,
            "strand": 2,
        },
    )
    assert resumed.get("resumed") is True

    # Switching strand to 1 (wrong strand for TruSeq -s 2) MUST invalidate A4 manifest
    rerun_wrong_strand = env_manager.run_in_stack(
        "rnaseq",
        "aria/scripts/rna_quantify.py",
        {
            "bam_files": preproc["alignment"]["bam_files"],
            "gtf_file": mini_spliced_reference["gtf"],
            "output_dir": str(counts_dir),
            "threads": 2,
            "strand": 1,
        },
    )
    assert not rerun_wrong_strand.get("resumed")
    assert rerun_wrong_strand["assignment_summary"]["pct_assigned"] < 5.0


def test_live_family2_tag_3prime_and_5prime_e2e_with_tpm_suppression(
    tmp_path, mini_spliced_reference
):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.memory.memory import ARIAMemory

    genes = mini_spliced_reference["genes"]
    rng = random.Random(77)

    # ── 2A: QuantSeq-style 3'-tag SE with 20-bp poly(A) tails -> full DE E2E ──
    fq3_dir = tmp_path / "tag3_study" / "raw"
    fq3_files: list[str] = []
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    for s_idx, s_name in enumerate(samples):
        is_treat = s_name.startswith("treat")
        recs = []
        r_i = 0
        for g_idx, g in enumerate(genes):
            base = 35 + ((s_idx + g_idx) % 5) * 2
            n_tags = base * 5 if (is_treat and g_idx < 3) else base
            mrna = g["mrna"]
            for k in range(n_tags):
                offset = 1090 + (k % 25)
                # 55-bp sense 3'-UTR sequence + 20-bp poly(A) tail = 75 bp
                seq = mrna[offset:offset + 55] + ("A" * 20)
                r_i += 1
                recs.append((f"{s_name}_r{r_i}", seq, "I" * len(seq)))
        rng.shuffle(recs)
        p = _write_fastq_gz(fq3_dir / f"{s_name}.fastq.gz", recs)
        fq3_files.append(str(p))

    mem = ARIAMemory(db_path=str(tmp_path / "mem3.sqlite"))
    agent = BulkRNAAgent(memory=mem, llm=None)
    run_out = agent.run(
        "exp_tag3",
        {
            "exp_context": {
                "n_cpus": 2,
                "rna_library_protocol": "quantseq_fwd",
                "modalities": {"bulk_RNA_raw": fq3_files},
                "genome_config": {
                    "star_index": mini_spliced_reference["star_index"],
                    "fasta": mini_spliced_reference["fasta"],
                    "gtf": mini_spliced_reference["gtf"],
                    "strand": "auto",
                    "min_strand_reads": 40,
                },
                "design": {
                    "main_factor": "condition",
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
            "biological_intent": {"summary": "QuantSeq 3prime tag DE"},
        },
    )
    assert run_out["status"] == "done", run_out
    findings = run_out["findings"]
    assert findings["library_protocol"] == "tag_3prime"
    assert findings["preprocessing"]["protocol"] == "tag_3prime"
    assert findings["preprocessing"]["quantification"]["strand_used"] == 1
    assert findings["tpm_status"]["computed"] is False

    de_dir = Path(findings["tpm_status"]["alternative"]).parent
    assert (de_dir / "counts_cpm.tsv").exists()
    assert not (de_dir / "counts_tpm.tsv").exists()

    sig_ids = set(findings["contrasts"][0]["up_gene_ids"])
    expected_up = {genes[0]["gene_id"], genes[1]["gene_id"], genes[2]["gene_id"]}
    assert expected_up <= sig_ids

    # ── 2B: 5'-CAGE style TSS sense tags ──
    fq5_dir = tmp_path / "tag5_study" / "raw"
    fq5_recs = []
    for g_idx, g in enumerate(genes[:6]):
        for k in range(20 + g_idx * 5):
            seq = g["mrna"][k % 10:(k % 10) + 65]
            fq5_recs.append((f"cage_{g_idx}_{k}", seq, "I" * 65))
    fq5_path = _write_fastq_gz(fq5_dir / "cage_rep1.fastq.gz", fq5_recs)
    counts5, pre5 = agent._run_preprocessing(
        "exp_tag5",
        [str(fq5_path)],
        {
            "n_cpus": 2,
            "rna_library_protocol": "cage",
            "genome_config": {
                "star_index": mini_spliced_reference["star_index"],
                "fasta": mini_spliced_reference["fasta"],
                "gtf": mini_spliced_reference["gtf"],
                "strand": "auto",
                "min_strand_reads": 40,
            },
        },
        {},
    )
    assert counts5 is not None, pre5
    assert pre5["protocol"] == "tag_5prime"
    assert pre5["suppress_tpm"] is True
    assert pre5["quantification"]["strand_used"] == 1
    df5 = pd.read_csv(counts5[0], sep="\t", index_col=0)
    assert df5.loc[genes[0]["gene_id"], "cage_rep1"] == 20
    assert df5.loc[genes[5]["gene_id"], "cage_rep1"] == 45


def test_live_family3_umi_bulk_collapses_pcr_duplicates_and_1bp_errors(
    tmp_path, mini_spliced_reference
):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.memory.memory import ARIAMemory

    genes = mini_spliced_reference["genes"]

    # Build 80 well-separated 8-bp UMIs (pairwise Hamming >= 3)
    umis: list[str] = []
    rng = random.Random(999)
    while len(umis) < 80:
        cand = "".join(rng.choice("ACGT") for _ in range(8))
        if all(sum(a != b for a, b in zip(cand, u)) >= 3 for u in umis):
            umis.append(cand)

    def _mutate_1bp(u: str) -> str:
        alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[u[-1]]
        return u[:-1] + alt

    # Sample has 25 true molecules on GENE01 and 60 true molecules on GENE02,
    # each sequenced 4x (3 exact copies + 1 copy with a 1-bp UMI error -> 340 raw reads!)
    recs = []
    r_id = 0
    for g, n_mol in ((genes[0], 25), (genes[1], 60)):
        tag_seq = g["mrna"][1100:1160]  # fixed 5' anchor per gene!
        for m in range(n_mol):
            true_umi = umis[m]
            err_umi = _mutate_1bp(true_umi)
            for copy_umi in (true_umi, true_umi, true_umi, err_umi):
                r_id += 1
                # 8-bp inline UMI + 4-bp spacer ("TATA") + 60-bp cDNA tag
                full_seq = copy_umi + "TATA" + tag_seq
                recs.append((f"umi_read_{r_id}", full_seq, "I" * len(full_seq)))

    fq_path = _write_fastq_gz(tmp_path / "umi_study" / "raw" / "umi_s1.fastq.gz", recs)
    mem = ARIAMemory(db_path=str(tmp_path / "mem_umi.sqlite"))
    agent = BulkRNAAgent(memory=mem, llm=None)

    counts_files, preproc = agent._run_preprocessing(
        "exp_umi",
        [str(fq_path)],
        {
            "n_cpus": 2,
            "rna_library_protocol": "umi_bulk",
            "umi_config": {
                "umi_loc": "read1",
                "umi_len": 8,
                "umi_skip": 4,
                "is_tag_assay": True,
            },
            "genome_config": {
                "star_index": mini_spliced_reference["star_index"],
                "fasta": mini_spliced_reference["fasta"],
                "gtf": mini_spliced_reference["gtf"],
                "strand": "auto",
                "min_strand_reads": 30,
            },
        },
        {},
    )
    assert counts_files is not None, preproc
    assert preproc["protocol"] == "umi_bulk"
    assert preproc["umi_dedup"]["mean_umi_dedup_rate_pct"] == pytest.approx(75.0, abs=0.5)

    df = pd.read_csv(counts_files[0], sep="\t", index_col=0)
    # Exact ground-truth molecule counts (25 and 60), NOT 100 and 240!
    assert int(df.loc[genes[0]["gene_id"], "umi_s1"]) == 25
    assert int(df.loc[genes[1]["gene_id"], "umi_s1"]) == 60


def test_live_family4_pooled_plate_brbseq_starsolo_to_deseq2_e2e(
    tmp_path, mini_spliced_reference
):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.memory.memory import ARIAMemory

    genes = mini_spliced_reference["genes"]
    rng = random.Random(2026)

    # 6 well barcodes of length 14 bp (BRB-seq preset: cb_len=14, umi_len=14 on R1; cDNA on R2)
    barcode_to_sample = {
        "AAACCCGGGTTTAA": "ctrl_1",
        "CCCGGGTTTAAACC": "ctrl_2",
        "GGGTTTAAACCCGG": "ctrl_3",
        "TTTAAACCCGGGTT": "treat_1",
        "ACGTACGTACGTAC": "treat_2",
        "TGCATGCATGCATG": "treat_3",
    }

    umi_pool: list[str] = []
    while len(umi_pool) < 300:
        u = "".join(rng.choice("ACGT") for _ in range(14))
        if all(sum(a != b for a, b in zip(u, prev)) >= 4 for prev in umi_pool):
            umi_pool.append(u)

    r1_recs = []
    r2_recs = []
    q_idx = 0
    for s_idx, (bc, s_name) in enumerate(barcode_to_sample.items()):
        is_treat = s_name.startswith("treat")
        for g_idx, g in enumerate(genes):
            base_mol = 28 + ((s_idx + g_idx) % 4) * 2
            n_mol = base_mol * 5 if (is_treat and g_idx < 3) else base_mol
            c_tag = g["mrna"][1080:1150]  # 70-bp sense 3' tag
            for m in range(n_mol):
                umi = umi_pool[m]
                # Emit 2 PCR copies per molecule to verify STARsolo UMI collapse!
                for _ in range(2):
                    q_idx += 1
                    r1_seq = bc + umi  # 28 bp barcode+UMI read
                    r1_recs.append((f"pool_{q_idx}", r1_seq, "I" * 28))
                    r2_recs.append((f"pool_{q_idx}", c_tag, "I" * 70))

    raw_dir = tmp_path / "brbseq_study" / "raw"
    r1_fq = _write_fastq_gz(raw_dir / "plate1_R1.fastq.gz", r1_recs)
    r2_fq = _write_fastq_gz(raw_dir / "plate1_R2.fastq.gz", r2_recs)

    mem = ARIAMemory(db_path=str(tmp_path / "mem_brb.sqlite"))
    agent = BulkRNAAgent(memory=mem, llm=None)

    run_out = agent.run(
        "exp_brbseq",
        {
            "exp_context": {
                "n_cpus": 2,
                "rna_library_protocol": "brbseq",
                "pooled_plate_manifest": {
                    "preset": "brbseq",
                    "barcode_to_sample": barcode_to_sample,
                },
                "modalities": {"bulk_RNA_raw": [str(r1_fq), str(r2_fq)]},
                "genome_config": {
                    "star_index": mini_spliced_reference["star_index"],
                    "fasta": mini_spliced_reference["fasta"],
                    "gtf": mini_spliced_reference["gtf"],
                },
                "design": {
                    "main_factor": "condition",
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
            "biological_intent": {"summary": "BRB-seq pooled 3prime bulk DE"},
        },
    )
    assert run_out["status"] == "done", run_out
    findings = run_out["findings"]
    solo = findings["preprocessing"]["starsolo"]
    assert solo["n_samples"] == 6
    assert solo["n_genes"] == 12
    assert solo["umi_deduplicated"] is True

    # Verify STARsolo collapsed the 2x PCR duplicates back to exact molecule counts!
    mat = pd.read_csv(solo["counts_matrix"], sep="\t", index_col=0)
    assert list(mat.columns) == sorted(barcode_to_sample, key=lambda k: k) or set(mat.columns) == set(barcode_to_sample.values())
    # For ctrl_1 (s_idx=0 where sorted keys put AAACCCGGGTTTAA first), g_idx=0 has 28 molecules (NOT 56!)
    assert int(mat.loc[genes[0]["gene_id"], "ctrl_1"]) == 28

    # Verify pyDESeq2 recovered the 3 upregulated genes and suppressed TPM!
    assert findings["tpm_status"]["computed"] is False
    de_dir = Path(findings["tpm_status"]["alternative"]).parent
    assert (de_dir / "counts_cpm.tsv").exists()
    assert not (de_dir / "counts_tpm.tsv").exists()
    up_ids = set(findings["contrasts"][0]["up_gene_ids"])
    assert {genes[0]["gene_id"], genes[1]["gene_id"], genes[2]["gene_id"]} <= up_ids


def test_live_fail_closed_ambiguous_strandedness_guard(
    tmp_path, mini_spliced_reference
):
    from aria.utils.environment_manager import env_manager

    genes = mini_spliced_reference["genes"]
    # Build a FASTQ with 30% sense reads and 70% antisense reads -> falls in [0.20, 0.40] ambiguous band!
    recs = []
    for k in range(30):
        seq = genes[0]["mrna"][100 + (k % 10):175 + (k % 10)]  # sense
        recs.append((f"sense_{k}", seq, "I" * 75))
    for k in range(70):
        seq = _revcomp(genes[1]["mrna"][100 + (k % 10):175 + (k % 10)])  # antisense
        recs.append((f"anti_{k}", seq, "I" * 75))

    fq = _write_fastq_gz(tmp_path / "ambig_R1.fastq.gz", recs)
    align_res = env_manager.run_in_stack(
        "rnaseq",
        "aria/scripts/rna_align.py",
        {
            "samples": [{
                "name": "ambig",
                "r1_trimmed": str(fq),
                "r2_trimmed": None,
                "paired": False,
                "read_layout": "single-end",
            }],
            "genome_dir": mini_spliced_reference["star_index"],
            "genome_fasta": mini_spliced_reference["fasta"],
            "gtf_file": mini_spliced_reference["gtf"],
            "output_dir": str(tmp_path / "aligned"),
            "threads": 2,
            "two_pass": False,
        },
    )
    assert align_res["status"] == "success", align_res

    quant_res = env_manager.run_in_stack(
        "rnaseq",
        "aria/scripts/rna_quantify.py",
        {
            "bam_files": align_res["bam_files"],
            "gtf_file": mini_spliced_reference["gtf"],
            "output_dir": str(tmp_path / "counts"),
            "threads": 2,
            "strand": "auto",
            "min_strand_reads": 50,
        },
    )
    assert quant_res["status"] == "error", quant_res
    assert quant_res["error_type"] == "AmbiguousStrandedness"
    assert 0.25 <= quant_res["sense_fraction"] <= 0.35
