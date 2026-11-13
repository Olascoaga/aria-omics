"""ADR-081: End-to-end governance, readiness, ChromatinAgent, ChromatinNarrator
(W-CLAIM / W-LEDGER), and HypothesisAgent tests for ChIP, CUT_AND_RUN, CUT_AND_TAG.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from aria.agents.modality_audit import (
    ChIPAuditAgent,
    CutAndRunAuditAgent,
    CutAndTagAuditAgent,
    build_capability_matrix,
)
from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.hypothesis.report_section import gather_evidence
from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
from aria.agents.narrative.run_ledger import build_run_ledger
from aria.agents.narrative.validators import find_causal_language
from aria.agents.orchestrator_agent import MODALITY_VALIDATION, OrchestratorAgent
from aria.utils.assay_contracts import validate_assay_contract

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("modality", ["ChIP", "CUT_AND_RUN", "CUT_AND_TAG"])
def test_epigenomic_modalities_are_production_and_dispatchable(modality):
    meta = MODALITY_VALIDATION[modality]
    assert meta["level"] == "production"
    assert meta["dispatch_enabled"] is True
    assert meta["promotion_adr"] == "ADR-081"
    assert OrchestratorAgent._blocked_modalities({modality: ["sample.bam"]}) == {}
    assert OrchestratorAgent._experimental_modalities({modality: ["sample.bam"]}) == {}


def test_chip_cutrun_cuttag_audit_agents_and_pooled_input_control_handling(tmp_path):
    ip_bams = [str(tmp_path / f"ip_{i}.bam") for i in range(6)]
    ctrl_bam = str(tmp_path / "pooled_input.bam")
    for p in ip_bams + [ctrl_bam]:
        Path(p).write_bytes(b"BAM\x01")

    # 3v3 IP contrast + 1 pooled Input control row -> must NOT trigger under_replicated!
    sample_sheet = [
        {"sample_id": "ip_0", "condition": "WT", "replicate": "r1", "role": "ip"},
        {"sample_id": "ip_1", "condition": "WT", "replicate": "r2", "role": "ip"},
        {"sample_id": "ip_2", "condition": "WT", "replicate": "r3", "role": "ip"},
        {"sample_id": "ip_3", "condition": "KO", "replicate": "r1", "role": "ip"},
        {"sample_id": "ip_4", "condition": "KO", "replicate": "r2", "role": "ip"},
        {"sample_id": "ip_5", "condition": "KO", "replicate": "r3", "role": "ip"},
        {"sample_id": "pooled_input", "condition": "Input", "replicate": "r1", "is_control": True},
    ]
    exp_ctx = {
        "genome": "mm39",
        "modalities": {"ChIP": ip_bams + [ctrl_bam]},
        "control_bams": [ctrl_bam],
        "design": {
            "sample_sheet": sample_sheet,
            "condition_col": "condition",
            "replicate_col": "replicate",
            "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
        },
    }
    chip_rec = ChIPAuditAgent().audit(exp_ctx)
    assert chip_rec["status"] == "green", chip_rec
    assert chip_rec["dispatch_policy"] == "auto"
    cap = build_capability_matrix(exp_ctx, modality_validation=MODALITY_VALIDATION)
    assert cap["cards"]["ChIP"]["status"] == "green", cap["cards"]["ChIP"]
    assert "ChIP" in cap["dispatch"]["allowed"]

    # Assay contract check also excludes pooled Input control from biological IP group size
    contract_res = validate_assay_contract(
        exp_ctx,
        "ChIP",
        ip_bams + [ctrl_bam],
    )
    assert [i for i in contract_res["findings"] if i["severity"] == "blocking"] == []

    # ChIP differential contrast without Input/IgG control -> yellow requires_ack
    no_ctrl_ctx = {
        "genome": "mm39",
        "modalities": {"ChIP": ip_bams},
        "design": {
            "sample_sheet": sample_sheet[:-1],
            "condition_col": "condition",
            "replicate_col": "replicate",
            "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
        },
    }
    chip_yellow = ChIPAuditAgent().audit(no_ctrl_ctx)
    assert chip_yellow["status"] == "yellow"
    assert chip_yellow["dispatch_policy"] == "requires_ack"

    # CUT&TAG control-free contrast with n>=3 IP replicates -> green auto
    ct_green = CutAndTagAuditAgent().audit(
        {
            "genome": "mm39",
            "modalities": {"CUT_AND_TAG": ip_bams},
            "design": no_ctrl_ctx["design"],
        }
    )
    assert ct_green["status"] == "green"
    assert ct_green["dispatch_policy"] == "auto"

    # Under-replicated n=1 contrast -> red blocked
    red_rec = CutAndRunAuditAgent().audit(
        {
            "genome": "mm39",
            "modalities": {"CUT_AND_RUN": ip_bams[:2]},
            "design": {
                "sample_sheet": [sample_sheet[0], sample_sheet[3]],
                "condition_col": "condition",
                "replicate_col": "replicate",
                "contrasts": [{"name": "KO_vs_WT", "numerator": "KO", "denominator": "WT"}],
            },
        }
    )
    assert red_rec["status"] == "red"
    assert red_rec["dispatch_policy"] == "blocked"


@pytest.mark.parametrize("modality", ["ChIP", "CUT_AND_RUN", "CUT_AND_TAG"])
def test_chromatin_narrator_and_ledger_compile_epigenomic_claims_with_zero_withheld(modality):
    caller = "macs3" if modality == "ChIP" else "seacr"
    remove_dup = False if modality == "CUT_AND_TAG" else True
    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                modality: {
                    "status": "done",
                    "findings": {
                        "preprocessing": {
                            "status": "success",
                            "data_type": modality,
                            "validation_level": "production",
                            "n_samples_aligned": 6,
                            "read_layout": "paired",
                            "filter": {
                                "mapq": 30,
                                "remove_mito": True,
                                "remove_duplicates": remove_dup,
                            },
                            "spikein_species": "sacCer3",
                            "bam_files": [
                                {
                                    "sample": f"s{i}",
                                    "status": "success",
                                    "bam": f"/tmp/s{i}.bam",
                                    "paired": True,
                                    "mapped_reads_post_filter": 1_200_000 + i * 10_000,
                                    "spikein_reads": 15_000 + i * 200,
                                }
                                for i in range(6)
                            ],
                        },
                        "qc": {
                            "status": "success",
                            "data_type": modality,
                            "validation_level": "production",
                            "n_samples": 6,
                            "mito_fraction": 0.012,
                            "dup_rate": 0.085,
                            "frip": 0.34,
                            "tss_enrichment": 9.8,
                            "nrf": 0.91,
                            "pbc1": 0.92,
                            "pbc2": 11.4,
                            "nsc": 1.18,
                            "rsc": 1.42,
                            "ip_control_jsd": 0.41,
                            "qc_complete": True,
                            "pass_qc": True,
                        },
                        "peaks": {
                            "status": "success",
                            "data_type": modality,
                            "validation_level": "production",
                            "peak_caller": caller,
                            "genome": "mm39",
                            "n_peaks": 14520,
                            "frip": 0.34,
                            "peaks_path": "/tmp/consensus.bed",
                            "peak_reproducibility": {
                                "strategy": "idr",
                                "status": "computed",
                                "idr": {
                                    "ran": True,
                                    "method": "gaussian_copula_em_idr",
                                    "idr_threshold": 0.05,
                                    "n_reproducible_peaks": 11240,
                                },
                            },
                        },
                        "peak_counts": {
                            "status": "success",
                            "ran": True,
                            "data_type": modality,
                            "validation_level": "production",
                            "n_peaks": 11240,
                            "n_samples": 6,
                            "matrix_path": "/tmp/peak_counts.tsv",
                        },
                        "differential_accessibility": {
                            "status": "success",
                            "ran": True,
                            "data_type": modality,
                            "validation_level": "production",
                            "method": "pydeseq2_spikein_size_factors",
                            "spikein_normalization": {
                                "used": True,
                                "method": "heterologous_spikein_size_factors",
                            },
                            "comparisons": [
                                {
                                    "comparison": "T3_vs_Control",
                                    "numerator": "T3",
                                    "denominator": "Control",
                                    "status": "success",
                                    "n_peaks_tested": 11240,
                                    "n_significant": 842,
                                    "n_up": 610,
                                    "n_down": 232,
                                    "padj_threshold": 0.05,
                                    "lfc_threshold": 1.0,
                                }
                            ],
                        },
                        "motifs": {
                            "status": "success",
                            "ran": True,
                            "data_type": modality,
                            "validation_level": "production",
                            "padj_max": 0.05,
                            "per_group": {
                                "T3_vs_Control::up_in_T3": {
                                    "n_peaks": 610,
                                    "n_enriched": 2,
                                    "top_motifs": [
                                        {
                                            "motif_id": "MA0139.1",
                                            "name": "CTCF",
                                            "log2_enrichment": 2.41,
                                            "padj": 1.2e-14,
                                        }
                                    ],
                                }
                            },
                        },
                    },
                }
            },
        }
    }

    exp_ctx = {"modalities": {modality: ["/tmp/s0.bam"]}, "genome": "mm39"}
    narrator = ChromatinNarrator()
    blocks = narrator.collect("chromatin_agent", agent_results["chromatin_agent"], exp_ctx)
    assert len(blocks) >= 5

    ledger = build_run_ledger(exp_ctx, agent_results)
    compiled = compile_public_claims(blocks, exp_ctx, run_ledger=ledger)
    assert compiled.withheld == [], f"Unexpected withheld claims for {modality}: {compiled.withheld}"
    assert len(compiled.claims) >= 5

    methods_text = " ".join(narrator.methods("chromatin_agent", agent_results["chromatin_agent"], exp_ctx))
    assert find_causal_language(methods_text) is None
    if modality != "ChIP":
        assert "SEACR" in methods_text
    assert "spike-in" in methods_text.lower() or "heterologous" in methods_text.lower()

    # HypothesisAgent adapter stamps actual epigenomic modality
    signals = gather_evidence(agent_results, ledger, exp_ctx)
    assert len(signals) == 1
    assert signals[0].entity == "CTCF"
    assert signals[0].modality == modality
    assert "motif_not_binding" in signals[0].caveats_inherited


def test_live_chromatin_agent_cut_and_run_seacr_idr_spikein_e2e(tmp_path):
    from unittest.mock import MagicMock

    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    if not (envs.get("chromatin") and envs.get("rna")):
        pytest.skip("aria-chromatin-env and aria-rna-env required for live CUT&RUN E2E")

    # Build 6 synthetic coordinate-sorted paired-end CUT&RUN BAMs + 1 IgG BAM via aria-chromatin-env
    gen_script = tmp_path / "make_cutrun_cohort.py"
    gen_script.write_text(
        """import sys
from pathlib import Path
import numpy as np
import pysam

out_dir = Path(sys.argv[1])
rng = np.random.default_rng(1081)
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 800_000}]}

peak_loci = [20_000 + i * 22_000 for i in range(30)]

def write_bam(path, is_treated, is_igg=False, seed_off=0):
    local_rng = np.random.default_rng(1081 + seed_off)
    unsorted = out_dir / (path.stem + ".u.bam")
    q = 0
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        # Uniform sparse background
        for pos in range(2_000, 790_000, 3_500):
            for _ in range(int(local_rng.integers(1, 3))):
                a = pysam.AlignedSegment(out.header)
                a.query_name = f"bg_{q}"; q += 1
                a.query_sequence = "A" * 50
                a.flag = 99; a.reference_id = 0; a.reference_start = pos
                a.mapping_quality = 42; a.cigar = [(0, 50)]
                a.next_reference_id = 0; a.next_reference_start = pos + 110; a.template_length = 160
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(a)
                b = pysam.AlignedSegment(out.header)
                b.query_name = a.query_name
                b.query_sequence = "T" * 50
                b.flag = 147; b.reference_id = 0; b.reference_start = pos + 110
                b.mapping_quality = 42; b.cigar = [(0, 50)]
                b.next_reference_id = 0; b.next_reference_start = pos; b.template_length = -160
                b.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(b)
        if not is_igg:
            for idx, pstart in enumerate(peak_loci):
                # First 18 peaks gain 4x in treated; remaining 12 invariant
                mult = 4.0 if (is_treated and idx < 18) else 1.0
                n_frags = int(local_rng.poisson(28.0 * mult))
                for _ in range(n_frags):
                    jit = int(local_rng.integers(0, 220))
                    s = pstart + jit
                    a = pysam.AlignedSegment(out.header)
                    a.query_name = f"pk_{q}"; q += 1
                    a.query_sequence = "C" * 50
                    a.flag = 99; a.reference_id = 0; a.reference_start = s
                    a.mapping_quality = 45; a.cigar = [(0, 50)]
                    a.next_reference_id = 0; a.next_reference_start = s + 100; a.template_length = 150
                    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                    out.write(a)
                    b = pysam.AlignedSegment(out.header)
                    b.query_name = a.query_name
                    b.query_sequence = "G" * 50
                    b.flag = 147; b.reference_id = 0; b.reference_start = s + 100
                    b.mapping_quality = 45; b.cigar = [(0, 50)]
                    b.next_reference_id = 0; b.next_reference_start = s; b.template_length = -150
                    b.query_qualities = pysam.qualitystring_to_array("I" * 50)
                    out.write(b)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    unsorted.unlink(missing_ok=True)

for r in range(1, 4):
    write_bam(out_dir / f"ctrl_r{r}.bam", False, False, r)
    write_bam(out_dir / f"treat_r{r}.bam", True, False, 10 + r)
write_bam(out_dir / "igg_ctrl.bam", False, True, 99)
""",
        encoding="utf-8",
    )
    subprocess.run(
        ["conda", "run", "-n", "aria-chromatin-env", "python", str(gen_script), str(tmp_path)],
        check=True,
    )

    ip_files = [str(tmp_path / f"ctrl_r{r}.bam") for r in range(1, 4)] + [
        str(tmp_path / f"treat_r{r}.bam") for r in range(1, 4)
    ]
    igg_file = str(tmp_path / "igg_ctrl.bam")

    meta_df = pd.DataFrame(
        {
            "sample_id": [f"ctrl_r{r}" for r in range(1, 4)] + [f"treat_r{r}" for r in range(1, 4)],
            "condition": ["ctrl"] * 3 + ["treat"] * 3,
            "replicate": ["r1", "r2", "r3", "r1", "r2", "r3"],
            "spikein_reads": [20_000, 20_500, 19_800, 20_100, 19_900, 20_200],
        }
    )
    meta_path = tmp_path / "sample_metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)

    mem = ARIAMemory(db_path=str(tmp_path / "mem.sqlite"))
    llm = MagicMock()
    agent = ChromatinAgent(mem, llm)
    agent.workspace = tmp_path / "workspace"
    agent.workspace.mkdir(parents=True, exist_ok=True)

    payload = {
        "exp_context": {
            "genome": "mm39",
            "organism": "Mus musculus",
            "modalities": {"CUT_AND_RUN": ip_files},
            "control_bams": [igg_file],
            "sample_metadata_path": str(meta_path),
            "design": {
                "condition_col": "condition",
                "replicate_col": "replicate",
                "contrasts": [{"name": "treat_vs_ctrl", "numerator": "treat", "denominator": "ctrl"}],
            },
        },
        "biological_intent": {"question": "CUT&RUN differential binding between treat and ctrl"},
    }
    res = agent.run("exp_cutrun_live", payload)
    assert res["status"] == "done", res
    cr_findings = res["findings"]["CUT_AND_RUN"]["findings"]
    assert cr_findings["peaks"]["status"] == "success"
    assert cr_findings["peaks"]["peak_caller"] == "seacr"
    assert cr_findings["peaks"]["n_peaks"] >= 18
    da = cr_findings["differential_accessibility"]
    assert da["status"] == "success"
    assert da["ran"] is True
    assert da["spikein_normalization"]["used"] is True
    comp = da["comparisons"][0]
    assert comp["n_up"] >= 12
    assert comp["n_down"] == 0
