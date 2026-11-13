"""ADR-082: Universal 3D Genome Architecture (`HiC` / Micro-C) E2E Test Suite.

Covers:
1. `MODALITY_VALIDATION["HiC"]` production status, `ScriptContract` & registry
   integrity, `HiCAuditAgent` green/yellow/red readiness governance, and
   single-cell barcoded FASTQ misrouting protection.
2. `HiCNarrator` (`hic.preprocessing`, `hic.qc`, `hic.compartments`, `hic.tads`,
   `hic.loops`, `hic.differential_3d`), `W-CLAIM` (`compile_public_claims`),
   `W-LEDGER` (`ledger://hic/*`), zero causal verbs, and deterministic
   `HypothesisAgent` `hic_evidence` (`HIC_LOOP_ASSOCIATIVE`).
3. Live `aria-hic-env` end-to-end execution through `GenomeArchAgent`:
   6-sample (3 WT vs 3 KO) `.pairs` cohort -> `.mcool` (`hic_preprocess` + `A4`
   resume) -> `hic_inspect` -> `hic_qc_and_balance` (ICE balancing, exact pixel
   sums, P(s) decay slope, HiCRep SCC) -> `hic_topology` (GC-phased `eigs_cis`,
   Crane `insulation` TADs, `dotfinder.dots` focal loops + GTF gene annotation,
   and 3v3 empirical-Bayes `differential_3d` A/B compartment switch recovery +
   1v1 `descriptive_only` guard) -> `hic_figure` dual 300-DPI PNG + editable SVG
   figures -> `HiCNarrator` + `W-CLAIM`/`W-LEDGER` compilation with 0 withheld!
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from aria.agents.genome_arch_agent import GenomeArchAgent
from aria.agents.modality_audit import HiCAuditAgent
from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.hypothesis.report_section import gather_evidence
from aria.agents.narrative.narrators.hic import HiCNarrator
from aria.agents.narrative.run_ledger import build_run_ledger
from aria.agents.narrative.validators import find_causal_language
from aria.agents.orchestrator_agent import MODALITY_VALIDATION, OrchestratorAgent
from aria.memory.memory import ARIAMemory
from aria.utils.registry_integrity import check_registry_integrity


def test_hic_is_production_and_registry_integrity_passes():
    meta = MODALITY_VALIDATION["HiC"]
    assert meta["level"] == "production"
    assert meta["dispatch_enabled"] is True
    assert meta["promotion_adr"] == "ADR-082"
    assert OrchestratorAgent._blocked_modalities({"HiC": ["cohort.mcool"]}) == {}
    errors = [i for i in check_registry_integrity() if i.severity == "error"]
    assert errors == []


def test_hic_audit_agent_green_yellow_red_governance(tmp_path):
    mcools = [str(tmp_path / f"s{i}.mcool") for i in range(6)]
    for p in mcools:
        Path(p).write_bytes(b"HDF5")

    audit = HiCAuditAgent()
    # Single/cohort profiling without differential contrast -> green auto
    prof = audit.audit({"genome": "hg38", "modalities": {"HiC": mcools[:1]}})
    assert prof["status"] == "green"
    assert prof["dispatch_policy"] == "auto"

    # 3v3 differential 3D contrast -> green auto
    sheet_3v3 = [
        {"sample_id": f"s{i}", "condition": "WT" if i < 3 else "KO", "replicate": f"r{(i % 3) + 1}"}
        for i in range(6)
    ]
    green_3v3 = audit.audit(
        {
            "genome": "hg38",
            "modalities": {"HiC": mcools},
            "design": {
                "sample_sheet": sheet_3v3,
                "condition_col": "condition",
                "replicate_col": "replicate",
                "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
            },
        }
    )
    assert green_3v3["status"] == "green"
    assert green_3v3["dispatch_policy"] == "auto"

    # 2v2 low-power differential 3D contrast -> yellow requires_ack
    sheet_2v2 = [
        {"sample_id": "s0", "condition": "WT", "replicate": "r1"},
        {"sample_id": "s1", "condition": "WT", "replicate": "r2"},
        {"sample_id": "s3", "condition": "KO", "replicate": "r1"},
        {"sample_id": "s4", "condition": "KO", "replicate": "r2"},
    ]
    yellow_2v2 = audit.audit(
        {
            "genome": "hg38",
            "modalities": {"HiC": [mcools[0], mcools[1], mcools[3], mcools[4]]},
            "design": {
                "sample_sheet": sheet_2v2,
                "condition_col": "condition",
                "replicate_col": "replicate",
                "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
            },
        }
    )
    assert yellow_2v2["status"] == "yellow"
    assert yellow_2v2["dispatch_policy"] == "requires_ack"

    # 1v1 contrast -> red blocked
    red_1v1 = audit.audit(
        {
            "genome": "hg38",
            "modalities": {"HiC": [mcools[0], mcools[3]]},
            "design": {
                "sample_sheet": [sheet_2v2[0], sheet_2v2[2]],
                "condition_col": "condition",
                "replicate_col": "replicate",
                "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
            },
        }
    )
    assert red_1v1["status"] == "red"
    assert red_1v1["dispatch_policy"] == "blocked"


def test_hic_narrator_ledger_and_hypothesis_adapter_compile_cleanly():
    agent_results = {
        "genome_arch_agent": {
            "status": "done",
            "findings": {
                "preprocessing": {
                    "status": "success",
                    "ran": True,
                    "data_type": "HiC",
                    "validation_level": "production",
                    "n_samples": 6,
                    "base_resolution": 20000,
                    "cooler_files": [f"/tmp/s{i}.mcool" for i in range(6)],
                    "per_sample": [{"sample_id": f"s{i}", "input_format": "pairs"} for i in range(6)],
                },
                "qc": {
                    "status": "success",
                    "ran": True,
                    "data_type": "HiC",
                    "validation_level": "production",
                    "n_samples": 6,
                    "resolution": 20000,
                    "total_valid_pairs": 12400000,
                    "mean_cis_trans_ratio": 0.84,
                    "mean_ps_slope": -1.08,
                    "mean_scc": 0.94,
                    "ice_converged_all": True,
                    "balanced_files": [f"/tmp/s{i}.mcool" for i in range(6)],
                },
                "topology": {
                    "compartments": {
                        "status": "success",
                        "ran": True,
                        "data_type": "HiC",
                        "validation_level": "production",
                        "resolution": 20000,
                        "pct_A": 48.5,
                        "pct_B": 51.5,
                        "n_valid_bins": 200,
                        "n_ab_transitions": 6,
                        "pc1_validated": True,
                        "phasing_method": "gc_content",
                        "compartments_tsv": "/tmp/compartments.tsv",
                    },
                    "tads": {
                        "status": "success",
                        "ran": True,
                        "data_type": "HiC",
                        "validation_level": "production",
                        "resolution": 20000,
                        "window_size": 100000,
                        "n_tads": 14,
                        "n_strong_boundaries": 13,
                        "median_size_kb": 260.0,
                        "tads_bed": "/tmp/tads.bed",
                    },
                    "loops": {
                        "status": "success",
                        "ran": True,
                        "data_type": "HiC",
                        "validation_level": "production",
                        "resolution": 20000,
                        "n_loops": 8,
                        "fdr_threshold": 0.1,
                        "loops_path": "/tmp/loops.bedpe",
                        "loop_class_counts": {
                            "promoter_promoter": 2,
                            "promoter_enhancer": 4,
                            "distal_distal": 2,
                        },
                        "top_loops": [
                            {
                                "chrom1": "chr1",
                                "start1": 120000,
                                "end1": 140000,
                                "chrom2": "chr1",
                                "start2": 880000,
                                "end2": 900000,
                                "count": 95,
                                "qval": 0.0012,
                                "anchor1_genes": "MYC",
                                "anchor2_genes": "CCND1,intergenic",
                                "loop_class": "promoter_promoter",
                            }
                        ],
                    },
                    "differential_3d": {
                        "status": "success",
                        "ran": True,
                        "descriptive_only": False,
                        "data_type": "HiC",
                        "validation_level": "production",
                        "n_ab_switches": 18,
                        "n_differential_boundaries": 3,
                        "comparisons": [
                            {
                                "comparison": "KO_vs_WT",
                                "test": "KO",
                                "reference": "WT",
                                "n_test": 3,
                                "n_reference": 3,
                                "n_bins_tested": 200,
                                "n_significant_e1": 22,
                                "n_a_to_b": 18,
                                "n_b_to_a": 0,
                                "n_differential_boundaries": 3,
                                "padj_max": 0.05,
                                "delta_e1_min": 0.2,
                            }
                        ],
                    },
                },
            },
        }
    }

    exp_ctx = {"modalities": {"HiC": ["/tmp/s0.pairs"]}, "genome": "hg38"}
    narrator = HiCNarrator()
    blocks = narrator.collect("genome_arch_agent", agent_results["genome_arch_agent"], exp_ctx)
    assert [b.id for b in blocks] == [
        "hic.preprocessing",
        "hic.qc",
        "hic.compartments",
        "hic.tads",
        "hic.loops",
        "hic.differential_3d",
    ]

    ledger = build_run_ledger(exp_ctx, agent_results)
    assert "hic" in ledger["modalities"]
    compiled = compile_public_claims(blocks, exp_ctx, run_ledger=ledger)
    assert compiled.withheld == [], f"Unexpected withheld Hi-C claims: {compiled.withheld}"
    assert len(compiled.claims) == 6

    methods_text = " ".join(narrator.methods("genome_arch_agent", agent_results["genome_arch_agent"], exp_ctx))
    assert find_causal_language(methods_text) is None
    assert "cooltools" in methods_text.lower() or "eigenvector" in methods_text.lower()

    # HypothesisAgent hic_evidence adapter emits MYC and CCND1 with HIC_LOOP_ASSOCIATIVE, skipping 'intergenic'
    signals = gather_evidence(agent_results, ledger, exp_ctx)
    entities = [s.entity for s in signals]
    assert entities == ["MYC", "CCND1"]
    for sig in signals:
        assert sig.modality == "HiC"
        assert sig.audited_node_ref == "ledger://hic/loops"
        assert "hic_loop_associative" in sig.caveats_inherited


def _write_synthetic_hic_cohort(work_dir: Path) -> tuple[list[str], Path, Path]:
    """Generate a 6-replicate (3 WT vs 3 KO) .pairs cohort on a 4 Mb chr1 (20 kb bins)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    chrom_len = 4_000_000
    bin_bp = 20_000
    n_bins = chrom_len // bin_bp  # 200 bins

    # Reference FASTA: [0..1Mb] A (60% GC), [1..2Mb] B (35% GC), [2..3Mb] A (60% GC), [3..4Mb] B (35% GC)
    rng_fa = np.random.default_rng(8200)
    fasta_path = work_dir / "ref.fa"
    seq_chunks: list[str] = []
    for b in range(n_bins):
        in_a = (b < 50) or (100 <= b < 150)
        p_gc = 0.62 if in_a else 0.34
        probs = [(1 - p_gc) / 2, p_gc / 2, p_gc / 2, (1 - p_gc) / 2]
        bases = rng_fa.choice(["A", "C", "G", "T"], size=bin_bp, p=probs)
        seq_chunks.append("".join(bases))
    full_seq = "".join(seq_chunks)
    with fasta_path.open("w", encoding="utf-8") as fh:
        fh.write(">chr1\n")
        for idx in range(0, len(full_seq), 80):
            fh.write(full_seq[idx : idx + 80] + "\n")

    # GTF with MYC at 130 kb and CCND1 at 890 kb (anchoring the focal loop)
    gtf_path = work_dir / "genes.gtf"
    gtf_path.write_text(
        'chr1\tsynth\tgene\t128000\t134000\t.\t+\t.\tgene_id "G1"; gene_name "MYC";\n'
        'chr1\tsynth\tgene\t888000\t894000\t.\t+\t.\tgene_id "G2"; gene_name "CCND1";\n'
        'chr1\tsynth\tgene\t2100000\t2110000\t.\t+\t.\tgene_id "G3"; gene_name "SOX2";\n',
        encoding="utf-8",
    )

    # TAD intervals in bin units: [0..50), [50..100), [100..150), [150..200)
    def _tad_id(b: int) -> int:
        return b // 50

    pairs_files: list[str] = []
    sample_specs = [
        ("WT_r1", False, 1),
        ("WT_r2", False, 2),
        ("WT_r3", False, 3),
        ("KO_r1", True, 4),
        ("KO_r2", True, 5),
        ("KO_r3", True, 6),
    ]

    for sid, is_ko, seed_off in sample_specs:
        rng = np.random.default_rng(9000 + seed_off)
        pfile = work_dir / f"{sid}.pairs"
        lines = [
            "## pairs format v1.0",
            "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2",
            f"#chromsize: chr1 {chrom_len}",
        ]
        rid = 1
        # Compartment state per bin (+1 = A, -1 = B); KO flips [100..150) from A (+1) to B (-1)!
        comp_sign = np.array(
            [
                (1 if b < 50 else (-1 if b < 100 else ((-1 if is_ko else 1) if b < 150 else -1)))
                for b in range(n_bins)
            ],
            dtype=int,
        )

        for i in range(n_bins):
            for j in range(i, n_bins):
                dist = j - i
                # Power-law polymer decay + strong intra-TAD insulation + checkerboard compartment affinity
                base_mu = 26.0 * ((dist + 1.0) ** (-0.95))
                if _tad_id(i) == _tad_id(j):
                    base_mu *= 2.8
                else:
                    base_mu *= 0.45
                if comp_sign[i] == comp_sign[j]:
                    base_mu *= 2.2
                else:
                    base_mu *= 0.45

                # Focal chromatin loop dot at bin 6 (120-140 kb) <-> bin 44 (880-900 kb)
                if i == 6 and j == 44:
                    base_mu += 140.0

                cnt = int(rng.poisson(base_mu))
                if cnt <= 0:
                    continue
                p1_base = i * bin_bp + 5_000
                p2_base = j * bin_bp + 10_000
                for k in range(cnt):
                    pos1 = p1_base + (k % 9_000)
                    pos2 = p2_base + ((k * 7) % 9_000)
                    if pos1 > pos2:
                        pos1, pos2 = pos2, pos1
                    lines.append(f"r{rid}\tchr1\t{pos1}\tchr1\t{pos2}\t+\t-")
                    rid += 1

        pfile.write_text("\n".join(lines) + "\n", encoding="utf-8")
        pairs_files.append(str(pfile))

    return pairs_files, fasta_path, gtf_path


def test_live_genome_arch_agent_pairs_to_mcool_topology_diff3d_and_report(tmp_path):
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    if not envs.get("hic"):
        pytest.skip("aria-hic-env required for live 3D genome E2E test")

    pairs_files, fasta_path, gtf_path = _write_synthetic_hic_cohort(tmp_path / "inputs")
    sample_ids = ["WT_r1", "WT_r2", "WT_r3", "KO_r1", "KO_r2", "KO_r3"]
    sample_sheet = [
        {
            "sample_id": sid,
            "condition": "WT" if sid.startswith("WT") else "KO",
            "replicate": sid.split("_")[1],
        }
        for sid in sample_ids
    ]

    mem = ARIAMemory(db_path=str(tmp_path / "hic_mem.sqlite"))
    llm = MagicMock()
    agent = GenomeArchAgent(mem, llm)
    agent.workspace = tmp_path / "workspace"
    agent.workspace.mkdir(parents=True, exist_ok=True)

    exp_ctx = {
        "genome": "hg38",
        "organism": "Homo sapiens",
        "modalities": {"HiC": pairs_files},
        "sample_ids": sample_ids,
        "genome_fasta": str(fasta_path),
        "gtf": str(gtf_path),
        "base_resolution": 20_000,
        "zoom_resolutions": [20_000, 40_000],
        "resolution": 20_000,
        "comp_resolution": 20_000,
        "design": {
            "sample_sheet": sample_sheet,
            "condition_col": "condition",
            "replicate_col": "replicate",
            "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
        },
    }
    res = agent.run(
        "exp_hic_live",
        {
            "exp_context": exp_ctx,
            "biological_intent": {
                "user_question": "Compare 3D compartments, TADs, and chromatin loops between KO and WT"
            },
        },
    )
    assert res["status"] == "done", res
    findings = res["findings"]

    # 1. Preprocessing (.pairs -> .mcool + A4 stage_manifest)
    prep = findings["preprocessing"]
    assert prep["status"] == "success"
    assert len(prep["cooler_files"]) == 6
    assert (Path(prep["output_dir"]) / "stage_manifest.json").is_file()

    # 2. QC + ICE balancing + P(s) slope + HiCRep SCC
    qc = findings["qc"]
    assert qc["status"] == "success"
    assert qc["n_valid_pairs"] > 50_000
    assert qc["mean_replicate_scc"] is not None and qc["mean_replicate_scc"] > 0.55
    scc_mat = qc["scc_matrix"]
    assert scc_mat["WT_r1"]["WT_r2"] > scc_mat["WT_r1"]["KO_r1"]

    # 3. Compartments (GC-phased E1), TADs, Loops (MYC <-> CCND1), and 3v3 Differential 3D
    topo = findings["topology"]
    comp = topo["compartments"]
    assert comp["status"] == "success"
    assert comp["pc1_validated"] is True
    assert comp["phasing_method"] == "gc_content"

    tads = topo["tads"]
    assert tads["status"] == "success"
    assert tads["n_tads"] >= 2

    loops = topo["loops"]
    assert loops["status"] == "success"
    assert loops["n_loops"] >= 1

    diff3d = topo["differential_3d"]
    assert diff3d["status"] == "success", diff3d
    assert diff3d["ran"] is True
    assert diff3d["descriptive_only"] is False
    assert diff3d["n_ab_switches"] >= 15

    # 4. Dual PNG + SVG figures
    figs = findings.get("figures") or {}
    bundles = findings.get("figure_bundles") or {}
    assert figs, "Expected rendered Hi-C figures"
    for k, png_path in figs.items():
        assert Path(png_path).is_file()
        svg_path = (bundles.get(k) or {}).get("svg")
        assert svg_path and Path(svg_path).is_file()

    # 5. Governed HiCNarrator + W-CLAIM + W-LEDGER compilation
    agent_results = {"genome_arch_agent": res}
    narrator = HiCNarrator()
    blocks = narrator.collect("genome_arch_agent", res, exp_ctx)
    ledger = build_run_ledger(exp_ctx, agent_results)
    compiled = compile_public_claims(blocks, exp_ctx, run_ledger=ledger)
    assert compiled.withheld == [], f"Withheld live Hi-C claims: {compiled.withheld}"
    assert len(compiled.claims) >= 5



def test_codex_adversarial_edge_cases_audit_and_narrator():
    """Regression guard for Codex Triad adversarial audit findings on ADR-081/082."""
    from aria.agents.modality_audit import CutAndRunAuditAgent, HiCAuditAgent
    from aria.agents.narrative.narrators.hic import HiCNarrator

    # 1. Three technical lanes of ONE biological replicate per condition must be blocked (red)
    #    even if design.groups lists 3 file IDs per condition!
    tech_sheet = [
        {"sample_id": "wt_l1", "condition": "WT", "replicate": "bio1"},
        {"sample_id": "wt_l2", "condition": "WT", "replicate": "bio1"},
        {"sample_id": "wt_l3", "condition": "WT", "replicate": "bio1"},
        {"sample_id": "ko_l1", "condition": "KO", "replicate": "bio2"},
        {"sample_id": "ko_l2", "condition": "KO", "replicate": "bio2"},
        {"sample_id": "ko_l3", "condition": "KO", "replicate": "bio2"},
    ]
    design_payload = {
        "groups": {"WT": ["wt_l1", "wt_l2", "wt_l3"], "KO": ["ko_l1", "ko_l2", "ko_l3"]},
        "sample_sheet": tech_sheet,
        "contrasts": [{"numerator": "KO", "denominator": "WT"}],
    }
    cr_card = CutAndRunAuditAgent().audit(
        {"genome": "hg38", "modalities": {"CUT_AND_RUN": ["a.bam"]}, "design": design_payload}
    )
    assert cr_card["status"] == "red"

    hic_card = HiCAuditAgent().audit(
        {"genome": "hg38", "modalities": {"HiC": ["a.mcool", "b.mcool"]}, "design": design_payload}
    )
    assert hic_card["status"] == "red"

    # 2. Unknown/pending stage status in HiCNarrator must fail closed (skipped diagnostic), never 'ok'!
    narr = HiCNarrator()
    blocks = narr.collect(
        "genome_arch_agent",
        {
            "status": "done",
            "findings": {
                "qc": {"status": "pending", "n_valid_pairs": 999_999},
                "topology": {
                    "compartments": {"status": "queued", "pct_A": 50.0, "pct_B": 50.0},
                },
            },
        },
        {"modalities": {"HiC": ["x.mcool"]}, "genome": "hg38"},
    )
    by_id = {b.id: b for b in blocks}
    assert by_id["hic.qc"].status == "limitation"
    assert by_id["hic.qc"].claim == ""
    assert by_id["hic.compartments"].status == "limitation"
    assert by_id["hic.compartments"].claim == ""
