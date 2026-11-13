"""ADR-081: Unit & live statistical tests for SEACR, Copula IDR, Broad Domains,
and Heterologous Spike-In PyDESeq2 Calibration.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from aria.utils.epigenomic_spikein_de import resolve_replicate_spikein_size_factors
from aria.utils.seacr_idr import (
    call_seacr_peaks,
    compute_idr_reproducible_peaks,
    fit_copula_idr_pair,
)

ROOT = Path(__file__).resolve().parents[1]


def test_scientific_utils_respect_subprocess_import_boundary():
    forbidden = ("aria.agents", "litellm", "pydantic")
    for rel in ("aria/utils/seacr_idr.py", "aria/utils/epigenomic_spikein_de.py"):
        tree = ast.parse((ROOT / rel).read_text(encoding="utf-8"), filename=rel)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert not alias.name.startswith(forbidden), f"{rel} imports {alias.name}"
            elif isinstance(node, ast.ImportFrom) and node.module:
                assert not node.module.startswith(forbidden), f"{rel} imports {node.module}"


def test_seacr_empirical_igg_stringent_and_relaxed_plus_spikein_non_mode(tmp_path):
    rng = np.random.default_rng(20261113)
    target_bg = tmp_path / "target.bedgraph"
    igg_bg = tmp_path / "igg.bedgraph"

    t_lines: list[str] = []
    c_lines: list[str] = []
    pos = 10_000
    # 250 low-signal background blocks in both target and IgG
    for _ in range(250):
        w = int(rng.integers(120, 320))
        t_val = float(rng.uniform(1.0, 3.0))
        c_val = float(rng.uniform(1.0, 3.0))
        t_lines.append(f"chr1\t{pos}\t{pos + w}\t{t_val:.3f}")
        c_lines.append(f"chr1\t{pos}\t{pos + w}\t{c_val:.3f}")
        pos += 2_000

    # 30 true focal enriched peaks in target only, plus 10 moderate shoulder domains
    true_starts: set[int] = set()
    for i in range(30):
        true_starts.add(pos)
        t_lines.append(f"chr1\t{pos}\t{pos + 600}\t{float(25.0 + i):.3f}")
        pos += 3_000
    for i in range(10):
        t_lines.append(f"chr1\t{pos}\t{pos + 350}\t{float(6.5 + 0.3 * i):.3f}")
        pos += 3_000

    target_bg.write_text("\n".join(t_lines) + "\n", encoding="utf-8")
    igg_bg.write_text("\n".join(c_lines) + "\n", encoding="utf-8")

    str_bed = tmp_path / "seacr_stringent.bed"
    rel_bed = tmp_path / "seacr_relaxed.bed"
    res_str = call_seacr_peaks(
        str(target_bg),
        str_bed,
        control_inputs=str(igg_bg),
        mode="stringent",
        norm="norm",
    )
    res_rel = call_seacr_peaks(
        str(target_bg),
        rel_bed,
        control_inputs=str(igg_bg),
        mode="relaxed",
        norm="norm",
    )

    assert res_str["status"] == "success"
    assert res_rel["status"] == "success"
    assert res_str["norm"] == "norm"
    assert res_str["n_peaks"] >= 28
    assert res_rel["n_peaks"] >= res_str["n_peaks"]
    assert res_rel["threshold"] <= res_str["threshold"]

    called_starts = {
        int(line.split("\t")[1])
        for line in str_bed.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    assert len(called_starts & true_starts) >= 28

    # Spike-in pre-scaled tracks must force norm="non" so global scale is preserved
    res_spike = call_seacr_peaks(
        str(target_bg),
        tmp_path / "seacr_spike.bed",
        control_inputs=str(igg_bg),
        mode="stringent",
        norm="norm",
        already_spikein_scaled=True,
    )
    assert res_spike["norm"] == "non"
    assert res_spike["already_spikein_scaled"] is True

    # Control-free top-quantile fallback
    res_topq = call_seacr_peaks(
        str(target_bg),
        tmp_path / "seacr_top005.bed",
        control_inputs=None,
        top_quantile=0.05,
    )
    assert res_topq["status"] == "success"
    assert 14 <= res_topq["n_peaks"] <= 16


def test_copula_idr_separates_reproducible_peaks_from_rank_noise(tmp_path):
    rng = np.random.default_rng(42)
    # 80 reproducible high-signal peaks + 220 uncorrelated noise peaks
    latent = rng.normal(loc=3.2, scale=0.6, size=80)
    s1_true = np.exp(latent + rng.normal(0.0, 0.2, size=80))
    s2_true = np.exp(latent + rng.normal(0.0, 0.2, size=80))

    s1_noise = np.exp(rng.normal(loc=0.5, scale=0.5, size=220))
    s2_noise = np.exp(rng.normal(loc=0.5, scale=0.5, size=220))

    s1 = np.concatenate([s1_true, s1_noise])
    s2 = np.concatenate([s2_true, s2_noise])

    fit = fit_copula_idr_pair(s1, s2)
    assert fit["params"]["rho"] > 0.55
    gidr = fit["global_idr"]
    # Vast majority of true reproducible peaks pass IDR <= 0.05, noise rejected
    assert int(np.sum(gidr[:80] <= 0.05)) >= 55
    assert int(np.sum(gidr[80:] <= 0.05)) <= 5

    # Write 3 replicate narrowPeak files and verify consensus IDR output
    peak_files: list[str] = []
    for rep_idx in range(3):
        rep_true = np.exp(latent + rng.normal(0.0, 0.22, size=80))
        rep_noise = np.exp(rng.normal(0.5, 0.5, size=220))
        scores = np.concatenate([rep_true, rep_noise])
        pfile = tmp_path / f"rep{rep_idx + 1}.narrowPeak"
        lines = []
        for idx, sc in enumerate(scores):
            start = 10_000 + idx * 2_000
            end = start + 400
            lines.append(f"chr1\t{start}\t{end}\tpk_{idx}\t500\t.\t{sc:.4f}\t10.0\t8.0\t200")
        pfile.write_text("\n".join(lines) + "\n", encoding="utf-8")
        peak_files.append(str(pfile))

    out_idr = tmp_path / "idr_consensus.narrowPeak"
    idr_res = compute_idr_reproducible_peaks(peak_files, out_idr, idr_threshold=0.05)
    assert idr_res["ran"] is True
    assert idr_res["n_pairs_evaluated"] == 3
    assert idr_res["n_reproducible_peaks"] >= 50
    assert idr_res["self_consistency_ratio"] is not None
    assert idr_res["self_consistency_ratio"] < 2.0
    assert out_idr.is_file()


def test_spikein_technical_replicate_sum_first_invariant():
    raw_meta = pd.DataFrame(
        {
            "condition": ["ctrl", "ctrl", "ctrl", "treat", "treat"],
            "replicate": ["r1", "r1", "r2", "r1", "r2"],
        },
        index=["C1_t1", "C1_t2", "C2", "T1", "T2"],
    )
    rep_meta = pd.DataFrame(
        {"condition": ["ctrl", "ctrl", "treat", "treat"]},
        index=["ctrl::r1", "ctrl::r2", "treat::r1", "treat::r2"],
    )

    # C1_t1 (1000) + C1_t2 (3000) must sum to 4000 == C2 (4000)!
    spike_counts = {
        "C1_t1": 1000.0,
        "C1_t2": 3000.0,
        "C2": 4000.0,
        "T1": 2000.0,
        "T2": 2000.0,
    }
    sf_counts, info_counts = resolve_replicate_spikein_size_factors(
        raw_meta,
        rep_meta,
        "condition",
        "replicate",
        spikein_counts=spike_counts,
    )
    assert sf_counts is not None
    assert info_counts["used"] is True
    assert np.isclose(sf_counts["ctrl::r1"], sf_counts["ctrl::r2"], rtol=1e-9)
    assert np.isclose(sf_counts["ctrl::r1"] / sf_counts["treat::r1"], 2.0, rtol=1e-9)
    assert np.isclose(np.exp(np.mean(np.log(sf_counts.to_numpy()))), 1.0, rtol=1e-9)

    # Supplying inverse coverage scale factors C / N_spike must yield identical size_factors!
    spike_scales = {k: 10_000.0 / v for k, v in spike_counts.items()}
    sf_scales, _ = resolve_replicate_spikein_size_factors(
        raw_meta,
        rep_meta,
        "condition",
        "replicate",
        spikein_scale_factors=spike_scales,
    )
    assert sf_scales is not None
    np.testing.assert_allclose(sf_counts.to_numpy(), sf_scales.to_numpy(), rtol=1e-9)

    # Missing sample spike-in fails closed
    bad_counts = dict(spike_counts)
    bad_counts.pop("T2")
    with pytest.raises(ValueError, match="missing values"):
        resolve_replicate_spikein_size_factors(
            raw_meta,
            rep_meta,
            "condition",
            "replicate",
            spikein_counts=bad_counts,
        )


def test_live_spikein_deseq2_recovers_global_unidirectional_loss(tmp_path):
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    if not envs.get("rna"):
        pytest.skip("aria-rna-env unavailable")

    rng = np.random.default_rng(2026)
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "ezh2i_1", "ezh2i_2", "ezh2i_3"]
    n_peaks = 250
    rows: dict[str, list[int]] = {}
    # 90% of H3K27me3 peaks lose ~80% of signal under EZH2i; 10% invariant background
    for i in range(n_peaks):
        pid = f"chr1:{i * 2000 + 1000}-{i * 2000 + 1800}"
        base = float(rng.uniform(120.0, 450.0))
        is_lost = i < 225
        ctrl_mu = base
        treat_mu = base * 0.20 if is_lost else base
        vals = [
            int(rng.negative_binomial(n=25, p=25.0 / (25.0 + ctrl_mu)))
            for _ in range(3)
        ] + [
            int(rng.negative_binomial(n=25, p=25.0 / (25.0 + treat_mu)))
            for _ in range(3)
        ]
        rows[pid] = vals

    counts_df = pd.DataFrame.from_dict(rows, orient="index", columns=samples)
    counts_df.index.name = "peak_id"
    counts_path = tmp_path / "h3k27me3_counts.tsv"
    counts_df.to_csv(counts_path, sep="\t")

    meta_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition": ["ctrl", "ctrl", "ctrl", "ezh2i", "ezh2i", "ezh2i"],
            "replicate": ["r1", "r2", "r3", "r1", "r2", "r3"],
            "spikein_reads": [50_000, 51_000, 49_000, 50_500, 49_500, 50_000],
        }
    )
    meta_path = tmp_path / "metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)

    res_spike = env_manager.run_in_stack(
        "rna",
        str(ROOT / "aria" / "scripts" / "chromatin_bulk_diffacc.py"),
        {
            "data_type": "CUT_AND_RUN",
            "counts_matrix_path": str(counts_path),
            "sample_metadata_path": str(meta_path),
            "condition_col": "condition",
            "replicate_col": "replicate",
            "comparisons": [{"name": "ezh2i_vs_ctrl", "numerator": "ezh2i", "denominator": "ctrl"}],
            "output_dir": str(tmp_path / "spike_de"),
            "padj_threshold": 0.05,
            "lfc_threshold": 1.0,
        },
    )
    assert res_spike["status"] == "success", res_spike
    assert res_spike["ran"] is True
    assert res_spike["spikein_normalization"]["used"] is True
    comp = res_spike["comparisons"][0]
    # Spike-in calibration preserves unidirectional loss: strong down calls, zero false up calls!
    assert comp["n_down"] >= 180
    assert comp["n_up"] == 0


def test_live_bam_bedgraph_broad_domains_and_bigwig(tmp_path):
    from aria.utils.environment_manager import env_manager

    envs = env_manager.check_environments()
    if not envs.get("chromatin"):
        pytest.skip("aria-chromatin-env unavailable")

    helper = tmp_path / "run_chrom_seacr_check.py"
    helper.write_text(
        """import json, sys
from pathlib import Path
import pysam
from aria.utils.seacr_idr import (
    bam_to_fragment_bedgraph,
    call_broad_domains,
    write_normalized_bigwig,
)

work = Path(sys.argv[1])
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 300_000}, {"SN": "chrM", "LN": 16_569}]}

def _make_bam(path, enriched=True):
    unsorted = work / (path.stem + ".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        q = 0
        for pos in range(1000, 290_000, 2500):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"bg_{q}"; q += 1
            a.query_sequence = "A" * 50
            a.flag = 99; a.reference_id = 0; a.reference_start = pos
            a.mapping_quality = 40; a.cigar = [(0, 50)]
            a.next_reference_id = 0; a.next_reference_start = pos + 130; a.template_length = 180
            a.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(a)
        if enriched:
            for pos in range(80_000, 140_000, 80):
                a = pysam.AlignedSegment(out.header)
                a.query_name = f"dom_{q}"; q += 1
                a.query_sequence = "G" * 50
                a.flag = 99; a.reference_id = 0; a.reference_start = pos
                a.mapping_quality = 42; a.cigar = [(0, 50)]
                a.next_reference_id = 0; a.next_reference_start = pos + 120; a.template_length = 170
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(a)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))

ip_bam = work / "h3k9me3_ip.bam"
ctrl_bam = work / "input_ctrl.bam"
_make_bam(ip_bam, enriched=True)
_make_bam(ctrl_bam, enriched=False)

bg_res = bam_to_fragment_bedgraph(ip_bam, work / "ip.bedgraph", scale_factor=2.0)
dom_res = call_broad_domains([str(ip_bam)], work / "domains.broadPeak", control_bams=[str(ctrl_bam)], bin_size=10_000, min_domain_bp=30_000)
bw_res = write_normalized_bigwig(ip_bam, work / "ip.bw", spikein_scale_factor=2.0)
print(json.dumps({"bg": bg_res, "dom": dom_res, "bw": bw_res}))
""",
        encoding="utf-8",
    )
    import os
    import subprocess

    proc = subprocess.run(
        ["conda", "run", "-n", "aria-chromatin-env", "python", str(helper), str(tmp_path)],
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        capture_output=True,
        text=True,
        check=True,
    )
    out = json.loads(proc.stdout.strip().splitlines()[-1])
    assert out["bg"]["n_fragments"] > 500
    assert out["dom"]["status"] == "success"
    assert out["dom"]["n_peaks"] >= 1
    assert Path(out["bw"]["bigwig_path"]).is_file()
    assert out["bw"]["normalization"] == "SPIKEIN"
    assert (tmp_path / "ip.bw").stat().st_size > 0



def test_codex_adversarial_edge_cases_seacr_idr_and_spikein(tmp_path):
    """Regression guard for Codex Triad adversarial audit findings on ADR-081."""
    from aria.utils.epigenomic_spikein_de import resolve_replicate_spikein_size_factors
    from aria.utils.seacr_idr import _pair_two_replicate_peaks, call_seacr_peaks

    # 1. Bookended/adjacent 0-based half-open BED intervals [0, 10) and [10, 20)
    #    share 0 bp and must NEVER be paired as an overlapping IDR peak!
    p1 = [("chr1", 0, 10, 50.0), ("chr1", 100, 200, 80.0)]
    p2 = [("chr1", 10, 20, 45.0), ("chr1", 150, 250, 75.0)]
    paired = _pair_two_replicate_peaks(p1, p2)
    assert len(paired) == 1
    assert paired[0][:3] == ("chr1", 100, 250)

    # 2. SEACR bedGraph inputs with explicit target_scale_factors actually scale signal!
    bg = tmp_path / "raw.bedgraph"
    bg.write_text("chr1\t0\t100\t10.0\nchr1\t500\t600\t2.0\n", encoding="utf-8")
    out_bed = tmp_path / "scaled.bed"
    res = call_seacr_peaks([str(bg)], out_bed, top_quantile=0.5, target_scale_factors=[3.0], work_dir=tmp_path / "w")
    assert res["status"] == "success"
    assert res["already_spikein_scaled"] is True
    cols = out_bed.read_text(encoding="utf-8").strip().splitlines()[0].split("\t")
    # 10.0 * 3.0 * 100 bp = 3000.0 total_signal
    assert abs(float(cols[3]) - 3000.0) < 1e-3

    # 3. Partial/invalid metadata spike-in column fails closed instead of silent conventional fallback!
    raw_meta = pd.DataFrame(
        {
            "condition": ["WT", "WT", "KO", "KO"],
            "replicate": ["r1", "r2", "r1", "r2"],
            "spikein_reads": [50_000, -10, 48_000, None],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    rep_meta = pd.DataFrame(
        {"condition": ["WT", "WT", "KO", "KO"]},
        index=["WT::r1", "WT::r2", "KO::r1", "KO::r2"],
    )
    with pytest.raises(ValueError, match="refusing silent fallback"):
        resolve_replicate_spikein_size_factors(raw_meta, rep_meta, "condition", "replicate")
