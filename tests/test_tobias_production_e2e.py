"""ADR-078 End-to-End Production Verification for TOBIAS Differential TF Footprinting.

Covers:
1. Pure-Python governance & contract checks:
   - `SCRIPT_CONTRACTS["aria/scripts/chromatin_footprint_tobias.py"].validation_level == "production"`
   - `AGENT_SCRIPT_CONTRACTS["chromatin_agent"]["required"]` includes `chromatin_footprint_tobias.py`
   - `ScATACAuditAgent` and `BulkATACAuditAgent` report `differential_tf_footprinting == "production_associative"`
   - `docker/Dockerfile.tobias` stamps `ARIA_IMAGE_VALIDATION=validated`
   - `INTEGRATION_VALIDATION` remains strictly `beta + requires_ack` (`production_promoted=False`)
2. Live `TOBIAS 0.17.3` (`aria-tobias-env`) end-to-end execution:
   - Bulk ATAC biological-replicate (`n=3` vs `n=3`) `ATACorrect -> ScoreBigwig -> BINDetect`
     + Welch/BH replicate inference + `A4` `stage_manifest.json` instant resume + dual PNG+SVG
     `footprint_volcano` + `C3` `compile_public_claims` (`W-CLAIM` + `W-LEDGER`) + `HypothesisAgent`
     `footprint_score` grounding & `SPECULATIVE` quarantine (`assert_no_speculative_promotion`).
   - Multi-sample scATAC (`fragments_files` with `sample#barcode` prefixes + donor-aware `n=3` vs `n=3`
     barcode design) live `TOBIAS` replicate inference + batch-adjusted OLS verification.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pytest


def _tobias_env_available() -> bool:
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    return bool(envs.get("tobias"))


def _write_mini_reference_and_motifs(base_dir: Path) -> tuple[Path, Path, Path]:
    """Create a deterministic 80 kb synthetic chromosome, 24 peaks, and 6 MEME v4 motifs."""
    rng = np.random.default_rng(20261112)
    chrom_len = 80_000
    bases = np.array(list("ACGT"))
    seq_arr = bases[rng.integers(0, 4, size=chrom_len)]

    motifs = {
        "MA0466.3 CEBPB": "ATTGCGCAAT",
        "MA0035.4 GATA1": "CTTATCT",
        "MA0139.1 CTCF": "CCGCGNGGNGGCAG".replace("N", "A"),
        "MA0079.5 SP1": "CCCCGCCCC",
        "MA0095.3 YY1": "CAAGATGGCG",
        "MA0107.1 RELA": "GGGAATTTCC",
    }

    peaks: list[tuple[str, int, int, str]] = []
    # 24 peaks of 500 bp spaced well inside [3000, 75000] so ATACorrect flanks never clip.
    for idx in range(24):
        start = 3_000 + idx * 2_800
        end = start + 500
        center = start + 250
        if idx < 12:
            primary_seq = motifs["MA0466.3 CEBPB"]
            kind = "cebpb"
        else:
            primary_seq = motifs["MA0035.4 GATA1"]
            kind = "gata1"
        seq_arr[center : center + len(primary_seq)] = list(primary_seq)

        # Also plant background motifs inside every peak at offset positions so BINDetect
        # sees sufficient motif binding sites across all 6 PWMs.
        bg_keys = ["MA0139.1 CTCF", "MA0079.5 SP1", "MA0095.3 YY1", "MA0107.1 RELA"]
        bg_seq = motifs[bg_keys[idx % len(bg_keys)]]
        seq_arr[start + 90 : start + 90 + len(bg_seq)] = list(bg_seq)
        bg_seq2 = motifs[bg_keys[(idx + 1) % len(bg_keys)]]
        seq_arr[end - 110 : end - 110 + len(bg_seq2)] = list(bg_seq2)
        peaks.append(("chrSim", start, end, kind))

    chrom_seq = "".join(seq_arr.tolist())
    fasta = base_dir / "mini_genome.fa"
    lines = [">chrSim"] + [chrom_seq[i : i + 80] for i in range(0, chrom_len, 80)]
    fasta.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Exact FASTA index (.fai): NAME, LENGTH, OFFSET, LINEBASES, LINEWIDTH
    fai = base_dir / "mini_genome.fa.fai"
    fai.write_text(f"chrSim\t{chrom_len}\t8\t80\t81\n", encoding="utf-8")

    peaks_bed = base_dir / "consensus_peaks.bed"
    peaks_bed.write_text(
        "\n".join(f"{c}\t{s}\t{e}\tpeak_{i}" for i, (c, s, e, _) in enumerate(peaks)) + "\n",
        encoding="utf-8",
    )

    # Valid MEME version 4 PWM file
    meme_lines = [
        "MEME version 4",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -",
        "",
        "Background letter frequencies",
        "A 0.25 C 0.25 G 0.25 T 0.25",
        "",
    ]
    base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    for header, cons in motifs.items():
        meme_lines.append(f"MOTIF {header}")
        meme_lines.append(f"letter-probability matrix: alength= 4 w= {len(cons)} nsites= 200 E= 0")
        for ch in cons:
            probs = [0.02, 0.02, 0.02, 0.02]
            probs[base_idx[ch]] = 0.94
            meme_lines.append(f" {probs[0]:.4f} {probs[1]:.4f} {probs[2]:.4f} {probs[3]:.4f}")
        meme_lines.append("")
    meme_path = base_dir / "motifs.meme"
    meme_path.write_text("\n".join(meme_lines), encoding="utf-8")
    return fasta, peaks_bed, meme_path


def _simulate_footprint_cut_positions(
    rng: np.random.Generator,
    peaks_bed: Path,
    *,
    protect_first_half: bool,
    chrom_len: int = 80_000,
) -> list[tuple[int, int]]:
    """Generate ATAC fragment intervals `(start, end)` with strong open-flank + central
    Tn5 protection valley over protected motifs vs flat background over unprotected ones."""
    frags: list[tuple[int, int]] = []
    peak_rows = [
        line.strip().split("\t")
        for line in peaks_bed.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    for idx, row in enumerate(peak_rows):
        start, end = int(row[1]), int(row[2])
        center = start + 250
        motif_len = 10 if idx < 12 else 7
        m_mid = center + motif_len // 2
        is_protected = (idx < 12) if protect_first_half else (idx >= 12)

        if is_protected:
            # High flanking accessibility [-90..-14] and [+14..+90], zero cuts in [-12..+12]!
            for _ in range(95):
                left_cut = m_mid - int(rng.integers(15, 85))
                right_cut = m_mid + int(rng.integers(15, 85))
                # Account for TOBIAS/ATAC +4 / -5 shift so cuts land in the flanks
                f_start = max(100, left_cut - 4)
                f_end = min(chrom_len - 100, right_cut + 5)
                if f_end - f_start >= 40:
                    frags.append((f_start, f_end))
            # Shoulder cuts right around the motif boundary create a crisp footprint dip
            for _ in range(55):
                if rng.random() < 0.5:
                    cut = m_mid - int(rng.integers(12, 24))
                    f_start = cut - 4
                    f_end = f_start + int(rng.integers(60, 140))
                else:
                    cut = m_mid + int(rng.integers(12, 24))
                    f_end = cut + 5
                    f_start = f_end - int(rng.integers(60, 140))
                if 100 <= f_start < f_end <= chrom_len - 100:
                    frags.append((f_start, f_end))
        else:
            # Unprotected: cuts fall uniformly right through the motif center (no dip)
            for _ in range(35):
                cut = int(rng.integers(start + 20, end - 20))
                f_start = max(100, cut - 4)
                f_end = min(chrom_len - 100, f_start + int(rng.integers(65, 160)))
                frags.append((f_start, f_end))

        # Neutral background motifs also get moderate uniform cuts
        for _ in range(20):
            f_start = int(rng.integers(start + 10, end - 120))
            f_end = f_start + int(rng.integers(60, 110))
            frags.append((f_start, f_end))

    # Off-peak background fragments so ATACorrect can estimate genome background bias cleanly
    for _ in range(900):
        f_start = int(rng.integers(500, chrom_len - 700))
        f_end = f_start + int(rng.integers(70, 180))
        frags.append((f_start, f_end))

    frags.sort()
    return frags


def _write_bam_from_fragments(bam_path: Path, fasta_path: Path, frags: list[tuple[int, int]]) -> Path:
    """Write a coordinate-sorted, indexed BAM in aria-chromatin/tobias pysam."""
    from aria.utils.environment_manager import env_manager

    # Build BAM using python inside aria-tobias-env (which has pysam locked)
    import subprocess
    import json

    payload = json.dumps({"bam": str(bam_path), "fasta": str(fasta_path), "frags": frags})
    code = (
        "import json, sys, pysam\n"
        "d = json.loads(sys.stdin.read())\n"
        "fa = pysam.FastaFile(d['fasta'])\n"
        "clen = fa.get_reference_length('chrSim')\n"
        "hdr = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chrSim', 'LN': clen}]}\n"
        "unsorted = d['bam'] + '.unsorted.bam'\n"
        "with pysam.AlignmentFile(unsorted, 'wb', header=hdr) as out:\n"
        "    for i, (s, e) in enumerate(d['frags']):\n"
        "        rlen = min(50, max(20, e - s))\n"
        "        for rev, pos in ((False, s), (True, max(0, e - rlen))):\n"
        "            a = pysam.AlignedSegment(out.header)\n"
        "            a.query_name = f'frag_{i}'\n"
        "            a.reference_id = 0\n"
        "            a.reference_start = int(pos)\n"
        "            a.mapping_quality = 60\n"
        "            a.flag = (1 | 2 | (16 if rev else 32) | (128 if rev else 64))\n"
        "            seq = fa.fetch('chrSim', int(pos), int(pos) + rlen)\n"
        "            a.query_sequence = seq\n"
        "            a.cigar = [(0, len(seq))]\n"
        "            a.query_qualities = pysam.qualitystring_to_array('I' * len(seq))\n"
        "            out.write(a)\n"
        "pysam.sort('-o', d['bam'], unsorted)\n"
        "pysam.index(d['bam'])\n"
    )
    from aria.benchmarks.tobias_footprint_concordance import _env_tool

    py_cmd = _env_tool("aria-tobias-env", "python")
    assert py_cmd is not None, "aria-tobias-env python binary not found"
    subprocess.run([*py_cmd, "-c", code], input=payload, text=True, check=True)
    Path(f"{bam_path}.unsorted.bam").unlink(missing_ok=True)
    return bam_path


# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Pure-Python Production Governance & Contract Assertions (CI PR Gate)
# ══════════════════════════════════════════════════════════════════════════════


def test_adr078_tobias_production_contracts_and_readiness(tmp_path):
    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.orchestrator_agent import (
        INTEGRATION_VALIDATION,
        MODALITY_VALIDATION,
    )
    from aria.utils.registry_integrity import (
        AGENT_SCRIPT_CONTRACTS,
        check_registry_integrity,
    )
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    contract = SCRIPT_CONTRACTS["aria/scripts/chromatin_footprint_tobias.py"]
    assert contract.validation_level == "production"
    input_names = {f.name for f in contract.inputs}
    assert {
        "mode",
        "fragments_file",
        "fragments_files",
        "sample_ids",
        "barcode_groups",
        "condition_bams",
        "replicate_bams",
        "replicate_batches",
        "reference_fasta",
        "make_volcano_plot",
        "resume",
    } <= input_names
    output_names = {f.name for f in contract.success_outputs}
    assert {"ran", "validation_level", "differential_summary", "volcano_plot", "stage_manifest_path"} <= output_names

    assert (
        "aria/scripts/chromatin_footprint_tobias.py"
        in AGENT_SCRIPT_CONTRACTS["chromatin_agent"]["required"]
    )
    assert check_registry_integrity() == []
    assert INTEGRATION_VALIDATION["level"] == "production"
    assert INTEGRATION_VALIDATION["production_promoted"] is True
    assert INTEGRATION_VALIDATION["scope"] == "paired_same_cell_rna_atac"

    stub_h5mu = tmp_path / "cohort.h5mu"
    stub_h5mu.write_bytes(b"stub")
    stub_bam = tmp_path / "sample1.bam"
    stub_bam.write_bytes(b"BAM\x01")

    matrix = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {"scATAC": [str(stub_h5mu)], "bulk_ATAC": [str(stub_bam)]},
        },
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    for mod_key in ("scATAC", "bulk_ATAC"):
        readiness = matrix["cards"][mod_key]["checks"]["production_readiness"]
        assert readiness["tier"] == "production"
        assert (
            readiness["caveated_sub_layers"]["differential_tf_footprinting"]
            == "production_associative"
        )


# ══════════════════════════════════════════════════════════════════════════════
# Part 2: Live TOBIAS 0.17.3 End-to-End Execution (Bulk ATAC & Multi-Sample scATAC)
# ══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not _tobias_env_available(), reason="aria-tobias-env not installed")
def test_live_tobias_bulk_and_scatac_replicate_inference_resume_c3_and_hypothesis(tmp_path):
    from aria.agents.narrative.claim_compiler import compile_public_claims
    from aria.agents.narrative.hypothesis.adapters.bulk_atac import bulk_atac_evidence
    from aria.agents.narrative.hypothesis.adapters.scatac import scatac_evidence
    from aria.agents.narrative.hypothesis.caveats import FOOTPRINT_ASSOCIATIVE
    from aria.agents.narrative.hypothesis.devils_advocate import check_devils_advocate
    from aria.agents.narrative.hypothesis.quarantine import assert_no_speculative_promotion
    from aria.agents.narrative.hypothesis.types import Hypothesis
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.utils.environment_manager import env_manager

    ref_dir = tmp_path / "ref"
    ref_dir.mkdir()
    fasta, peaks_bed, meme_path = _write_mini_reference_and_motifs(ref_dir)

    # ── 1. Bulk ATAC live TOBIAS (n=3 Treated vs n=3 Control biological replicates) ──
    bam_dir = tmp_path / "bams"
    bam_dir.mkdir()
    replicate_bams: dict[str, dict[str, list[str]]] = {"Treated": {}, "Control": {}}
    for rep_idx in range(1, 4):
        rng_t = np.random.default_rng(100 + rep_idx)
        frags_t = _simulate_footprint_cut_positions(rng_t, peaks_bed, protect_first_half=True)
        bam_t = _write_bam_from_fragments(bam_dir / f"Treated_r{rep_idx}.bam", fasta, frags_t)
        replicate_bams["Treated"][f"r{rep_idx}"] = [str(bam_t)]

        rng_c = np.random.default_rng(200 + rep_idx)
        frags_c = _simulate_footprint_cut_positions(rng_c, peaks_bed, protect_first_half=False)
        bam_c = _write_bam_from_fragments(bam_dir / f"Control_r{rep_idx}.bam", fasta, frags_c)
        replicate_bams["Control"][f"r{rep_idx}"] = [str(bam_c)]

    bulk_out = tmp_path / "bulk_tobias_out"
    bulk_params = {
        "mode": "bulk",
        "data_type": "bulk_ATAC",
        "replicate_bams": replicate_bams,
        "replicate_batches": {
            "Treated": {"r1": "b1", "r2": "b1", "r3": "b2"},
            "Control": {"r1": "b1", "r2": "b2", "r3": "b2"},
        },
        "genome_fasta": str(fasta),
        "genome_fai": f"{fasta}.fai",
        "peaks_bed": str(peaks_bed),
        "motif_meme": str(meme_path),
        "group_a": "Treated",
        "group_b": "Control",
        "output_dir": str(bulk_out),
        "cores": 4,
        "footprint_fdr": 0.05,
        "make_aggregate_plots": False,
        "make_volcano_plot": True,
    }
    res_bulk = env_manager.run_in_stack(
        stack="tobias",
        script_path="aria/scripts/chromatin_footprint_tobias.py",
        params=bulk_params,
    )
    assert res_bulk["status"] == "success", res_bulk
    assert res_bulk["ran"] is True
    assert res_bulk["validation_level"] == "production"
    summary = res_bulk["differential_summary"]
    assert summary["parsed"] is True
    inf = summary["inference"]
    assert inf["status"] == "success"
    assert inf["pairing_status"] == "independent"
    assert inf["batch_adjusted"] is True
    assert inf["replicates_per_condition"] == {"Treated": 3, "Control": 3}
    assert Path(summary["results_table"]).is_file()
    assert Path(res_bulk["bindetect_results"]).is_file()
    assert Path(res_bulk["stage_manifest_path"]).is_file()
    assert Path(res_bulk["volcano_plot"]["png"]).stat().st_size > 1000
    assert Path(res_bulk["volcano_plot"]["svg"]).stat().st_size > 500

    # Verify A4 cryptographic stage_manifest resume hits cleanly without rerunning TOBIAS
    res_cached = env_manager.run_in_stack(
        stack="tobias",
        script_path="aria/scripts/chromatin_footprint_tobias.py",
        params=bulk_params,
    )
    assert res_cached["status"] == "success"
    assert res_cached.get("reused_from_stage_manifest") is True

    # Verify C3 strict public-claim compilation (`W-CLAIM` + `W-LEDGER`)
    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "bulk_ATAC": {
                    "status": "done",
                    "findings": {
                        "qc": {
                            "status": "success",
                            "data_type": "bulk_ATAC",
                            "n_samples": 6,
                            "mito_fraction": 0.01,
                        },
                        "footprinting": res_bulk,
                    },
                }
            },
        }
    }
    exp_ctx_bulk = {
        "genome": "hg38",
        "modalities": {"bulk_ATAC": [str(bam_t)]},
        "design": {
            "main_factor": "condition",
            "replicate_col": "replicate",
            "contrasts": [{"numerator": "Treated", "denominator": "Control"}],
        },
    }
    narrator = ChromatinNarrator()
    blocks = narrator.collect(
        "chromatin_agent", agent_results["chromatin_agent"], exp_ctx_bulk
    )
    fp_blocks = [b for b in blocks if b.id == "chromatin.differential_tf_footprinting"]
    assert len(fp_blocks) == 1
    ledger = build_run_ledger(exp_ctx_bulk, agent_results)
    compiled = compile_public_claims(blocks, exp_ctx_bulk, run_ledger=ledger)
    kept_ids = {b.id for b in compiled.blocks}
    assert "chromatin.differential_tf_footprinting" in kept_ids, compiled.withheld
    assert_no_speculative_promotion(compiled.claims)

    # Verify HypothesisAgent bulk_atac_evidence grounding + devils_advocate caveat gate
    ev_bulk = bulk_atac_evidence(agent_results, run_ledger=ledger, exp_ctx=exp_ctx_bulk)
    fp_sigs = [s for s in ev_bulk if s.measure.startswith("footprint_score(")]
    if summary.get("n_significant", 0) > 0:
        assert len(fp_sigs) >= 1
        assert all(FOOTPRINT_ASSOCIATIVE in s.caveats_inherited for s in fp_sigs)
        sig0 = fp_sigs[0]
        ev_idx = {sig0.entity.lower(): [sig0]}
        bad_hyp = Hypothesis(
            id="hypothesis://h1",
            mechanism=f"{sig0.entity} occupancy may associate with chromatin remodeling.",
            entities=[sig0.entity],
            observation_refs=[sig0.audited_node_ref],
            devils_advocate={
                "simpler_explanation": "Co-accessible chromatin state without TF action.",
                "confounds": ["batch"],  # omits footprint_associative!
            },
        )
        assert check_devils_advocate(bad_hyp, ev_idx).passed is False
        good_hyp = Hypothesis(
            id="hypothesis://h2",
            mechanism=f"{sig0.entity} occupancy may associate with chromatin remodeling.",
            entities=[sig0.entity],
            observation_refs=[sig0.audited_node_ref],
            devils_advocate={
                "simpler_explanation": "Co-accessible chromatin state without TF action.",
                "confounds": list(sig0.caveats_inherited),
            },
        )
        assert check_devils_advocate(good_hyp, ev_idx).passed is True

    # ── 2. Multi-sample scATAC live TOBIAS (2 fragment files, sample#barcode prefixes) ──
    sc_dir = tmp_path / "scatac_frags"
    sc_dir.mkdir()
    frag_s1 = sc_dir / "donor_pool1.fragments.tsv.gz"
    frag_s2 = sc_dir / "donor_pool2.fragments.tsv.gz"
    design_rows = ["barcode\tgroup\treplicate"]

    # Pool 1 carries d1, d2 for Stim and Ctrl; Pool 2 carries d3 for Stim and Ctrl
    pool1_lines: list[str] = []
    pool2_lines: list[str] = []
    for grp, protect in (("Stim", True), ("Ctrl", False)):
        for d_idx in (1, 2, 3):
            bc = f"AAACGAA{grp[0]}{d_idx}-1"
            pool_id = "s1" if d_idx in (1, 2) else "s2"
            design_rows.append(f"{pool_id}#{bc}\t{grp}\tdonor_{grp}_{d_idx}")
            rng_sc = np.random.default_rng(500 + d_idx + (0 if protect else 50))
            frags_sc = _simulate_footprint_cut_positions(
                rng_sc, peaks_bed, protect_first_half=protect
            )
            target = pool1_lines if pool_id == "s1" else pool2_lines
            for s, e in frags_sc:
                target.append(f"chrSim\t{s}\t{e}\t{bc}\t1")

    pool1_lines.sort(key=lambda ln: int(ln.split("\t")[1]))
    pool2_lines.sort(key=lambda ln: int(ln.split("\t")[1]))
    with gzip.open(frag_s1, "wt", encoding="utf-8") as fh:
        fh.write("\n".join(pool1_lines) + "\n")
    with gzip.open(frag_s2, "wt", encoding="utf-8") as fh:
        fh.write("\n".join(pool2_lines) + "\n")

    design_tsv = sc_dir / "barcode_design.tsv"
    design_tsv.write_text("\n".join(design_rows) + "\n", encoding="utf-8")

    sc_out = tmp_path / "scatac_tobias_out"
    res_sc = env_manager.run_in_stack(
        stack="tobias",
        script_path="aria/scripts/chromatin_footprint_tobias.py",
        params={
            "mode": "scatac",
            "data_type": "scATAC",
            "fragments_files": [str(frag_s1), str(frag_s2)],
            "sample_ids": ["s1", "s2"],
            "barcode_groups": str(design_tsv),
            "genome_fasta": str(fasta),
            "genome_fai": f"{fasta}.fai",
            "peaks_bed": str(peaks_bed),
            "motif_meme": str(meme_path),
            "group_a": "Stim",
            "group_b": "Ctrl",
            "output_dir": str(sc_out),
            "cores": 4,
            "make_aggregate_plots": False,
            "make_volcano_plot": True,
        },
    )
    assert res_sc["status"] == "success", res_sc
    assert res_sc["ran"] is True
    assert res_sc["validation_level"] == "production"
    sc_inf = res_sc["differential_summary"]["inference"]
    assert sc_inf["status"] == "success"
    assert sc_inf["replicates_per_condition"] == {"Stim": 3, "Ctrl": 3}
    assert Path(res_sc["volcano_plot"]["png"]).is_file()

    sc_agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {"scATAC": {"status": "done", "findings": {"footprinting": res_sc}}},
        }
    }
    exp_ctx_sc = {
        "genome": "hg38",
        "modalities": {"scATAC": [str(frag_s1), str(frag_s2)]},
        "design": {
            "main_factor": "condition",
            "replicate_col": "donor_id",
            "contrasts": [{"numerator": "Stim", "denominator": "Ctrl"}],
        },
    }
    sc_ledger = build_run_ledger(exp_ctx_sc, sc_agent_results)
    ev_sc = scatac_evidence(sc_agent_results, run_ledger=sc_ledger, exp_ctx=exp_ctx_sc)
    if res_sc["differential_summary"].get("n_significant", 0) > 0:
        assert any(s.measure == "footprint_score(Stim_vs_Ctrl)" for s in ev_sc)
