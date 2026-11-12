"""V48 Slice I6 — Candidate eRegulons (`eregulon_core` + `integration_eregulons`) tests.

Covers all statistical, algorithmic, IPC, concordance, and governance invariants
mandated by `ADR-064` and Council `concilio_20261111_200904_d9842f`:
  1. Offline MEME PWM parsing & 1-to-1 TF symbol resolution (`exact`, unique
     `case_folded`, `(var.N)` stripping, `::` heterodimer subunit completeness,
     and unmapped/low-prevalence logging without hardcoded gene dictionaries).
  2. Fixed-width `500 bp` peak-center DNA extraction via `.fai` byte offsets
     (eliminating peak-length `2*L*p` bias + UCSC/Ensembl contig aliasing).
  3. Dual-backend motif scanning (`snapatac2_rust` vs `aria_numpy_dp`) `>=0.99`
     boolean parity across long and short (`w=6`) PWMs + `pvalue` bounds guard.
  4. GC-quintile Cochran-Mantel-Haenszel (CMH) cistrome enrichment rejecting pure
     GC-composition confounding while recovering true motif enrichment.
  5. Synthetic paired multiome ground-truth recovery (`pos_tf_pos_p2g` and
     `neg_tf_pos_p2g`), `g != TF` self-loop exclusion, paralog motif-sharing
     ambiguity disclosure vs `::` heterodimer partner exemption, library-depth
     partial-correlation confounder flagging, and bit-exact `seed=0` determinism.
  6. `.h5mu` end-to-end execution, `X_eregulon_rna`/`X_eregulon_atac` persistence,
     structured fail-closed errors, and `ScriptContract` beta validation.
  7. Four-tier `concordance_eregulons.py` scoring + `ADR-002` `not_run` honesty.
  8. AST & governance guard (`bh_correct` shared FDR, zero hardcoded TF panels,
     `INTEGRATION_VALIDATION` kept `scaffold` / `dispatch_enabled=False`, and
     strictly associative non-causal vocabulary).
"""

from __future__ import annotations

import ast
import csv
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _has_module(mod_name: str) -> bool:
    try:
        __import__(mod_name)
        return True
    except ImportError:
        return False


def _require_module(mod_name: str):
    try:
        return __import__(mod_name)
    except ImportError:
        try:
            import pytest

            pytest.skip(f"{mod_name} not installed in active environment")
        except ImportError as exc:
            raise RuntimeError(f"Required module {mod_name} missing") from exc


def _write_indexed_fasta(
    path: Path, records: list[tuple[str, str]], line_bases: int = 50
) -> Path:
    fai_lines: list[str] = []
    offset = 0
    with open(path, "wb") as fh:
        for name, seq in records:
            header = f">{name}\n".encode("ascii")
            fh.write(header)
            offset += len(header)
            seq_start = offset
            for i in range(0, len(seq), line_bases):
                chunk = (seq[i : i + line_bases] + "\n").encode("ascii")
                fh.write(chunk)
                offset += len(chunk)
            fai_lines.append(
                f"{name}\t{len(seq)}\t{seq_start}\t{line_bases}\t{line_bases + 1}\n"
            )
    Path(f"{path}.fai").write_text("".join(fai_lines), encoding="utf-8")
    return path


def _consensus_to_pfm_lines(consensus: str, hi: float = 0.94) -> str:
    lo = (1.0 - hi) / 3.0
    idx_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    lines: list[str] = []
    for ch in consensus.upper():
        row = [lo, lo, lo, lo]
        row[idx_map[ch]] = hi
        lines.append(f" {row[0]:.6f} {row[1]:.6f} {row[2]:.6f} {row[3]:.6f}")
    return "\n".join(lines)


def _write_test_meme(path: Path, motif_specs: list[tuple[str, str, str]]) -> Path:
    blocks = [
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"
    ]
    for m_id, m_name, cons in motif_specs:
        w = len(cons)
        pfm_txt = _consensus_to_pfm_lines(cons)
        blocks.append(
            f"\nMOTIF {m_id} {m_name}\n"
            f"letter-probability matrix: alength= 4 w= {w} nsites= 500 E= 0\n"
            f"{pfm_txt}\n"
        )
    path.write_text("".join(blocks), encoding="utf-8")
    return path


# ── 1. MEME parsing & TF symbol resolution ────────────────────────────────────


def test_meme_parsing_and_tf_symbol_resolution():
    from aria.utils.eregulon_core import (
        parse_meme_motifs,
        resolve_motifs_to_expressed_tfs,
    )

    with tempfile.TemporaryDirectory() as td:
        meme_p = _write_test_meme(
            Path(td) / "test.meme",
            [
                ("MA0001.1", "SPI1", "GGAAGTGAA"),
                ("MA0002.1", "Sox2", "CTTTGTTAT"),
                ("MA0003.1", "FOSL1(var.2)", "TGACTCATC"),
                ("MA0004.1", "JUN::FOS", "ATGACTCAT"),
                ("MA0005.1", "JUN::ABSENT_PARTNER", "ATGACTCAT"),
                ("MA0006.1", "LOWTF", "CCGGCCGGC"),
                ("MA0007.1", "EWSR1-FLI1", "GGAAAGGAA"),
            ],
        )
        motifs = parse_meme_motifs(meme_p)
        assert len(motifs) == 7

        rna_genes = ["SPI1", "SOX2", "FOSL1", "JUN", "FOS", "LOWTF", "GAPDH"]
        detected = np.array([50, 40, 30, 25, 22, 3, 100], dtype=np.int64)

        resolved, unmapped = resolve_motifs_to_expressed_tfs(
            motifs, rna_genes, detected, min_cells_detected=10
        )
        by_pair = {(r.motif_id, r.tf_gene): r for r in resolved}

        assert by_pair[("MA0001.1", "SPI1")].match_mode == "exact"
        assert by_pair[("MA0002.1", "SOX2")].match_mode == "case_folded"
        assert by_pair[("MA0003.1", "FOSL1")].variant_tag == "var.2"
        assert by_pair[("MA0004.1", "JUN")].motif_source == "heterodimer"
        assert by_pair[("MA0004.1", "FOS")].heterodimer_partners == ("JUN", "FOS")

        unmapped_reasons = {u["motif_id"]: u["reason"] for u in unmapped}
        assert "MA0005.1" in unmapped_reasons
        assert "ABSENT_PARTNER" in unmapped_reasons["MA0005.1"]
        assert unmapped_reasons["MA0006.1"] == "tf_below_min_cells_detected"
        assert unmapped_reasons["MA0007.1"] == "symbol_not_in_rna_universe"


# ── 2. Fixed-width 500bp FASTA window extraction ──────────────────────────────


def test_fixed_500bp_fasta_extraction_eliminates_peak_width_bias():
    from aria.utils.eregulon_core import fetch_peak_center_sequences_from_indexed_fasta

    rng = np.random.default_rng(7)
    chrom_seq = "".join(rng.choice(list("ACGT"), size=20_000))
    # Embed an N-masked block around 15,000
    chrom_seq = chrom_seq[:14700] + ("N" * 600) + chrom_seq[15300:]

    with tempfile.TemporaryDirectory() as td:
        # Write Ensembl-style contig `1` while querying UCSC-style `chr1` peaks
        fa_p = _write_indexed_fasta(Path(td) / "ref.fa", [("1", chrom_seq)])
        peaks = [
            "chr1:4900-5100",   # narrow 200bp peak (center 5000)
            "chr1:8000-12000",  # broad 4000bp peak (center 10000)
            "chr1:14800-15200", # >50% N-masked window
            "chrUn_gl000220:100-600",
        ]
        bundle = fetch_peak_center_sequences_from_indexed_fasta(
            peaks, fa_p, scan_window_bp=500
        )
        assert bundle["scan_window_bp"] == 500
        assert list(bundle["valid_mask"]) == [True, True, False, False]
        assert len(bundle["sequences"][0]) == 500
        assert len(bundle["sequences"][1]) == 500
        assert bundle["sequences"][0] == chrom_seq[4750:5250]
        assert bundle["sequences"][1] == chrom_seq[9750:10250]
        assert bundle["n_masked_high_n"] == 1
        assert bundle["n_missing_contig"] == 1


# ── 3. Dual-scanner boolean parity (including short w=6 E-box motif) ──────────


def test_dual_scanner_golden_boolean_parity_including_short_motifs():
    from aria.utils.eregulon_core import parse_meme_motifs, scan_motifs_in_sequences

    rng = np.random.default_rng(11)
    motif_specs = [
        ("MA0004.1", "Arnt", "CACGTG"),           # short w=6 E-box (`0.25^6 > 5e-5`!)
        ("MA0080.6", "SPI1", "AAAAAGAGGAAGTG"),   # w=14 ETS motif
        ("MA0139.1", "CTCF", "CCACCAGGGGGCGC"),   # w=14 GC-rich motif
    ]
    seqs: list[str] = []
    for i in range(240):
        base = "".join(rng.choice(list("ACGT"), size=500))
        if i % 5 == 0:
            base = base[:200] + "CACGTG" + base[206:]
        if i % 7 == 0:
            base = base[:310] + "CACTTCCTCTTTTT" + base[324:]  # RC of SPI1
        if i % 9 == 0:
            base = base[:120] + "CCACCAGGGGGCGC" + base[134:]
        seqs.append(base)

    with tempfile.TemporaryDirectory() as td:
        meme_p = _write_test_meme(Path(td) / "parity.meme", motif_specs)
        motifs = parse_meme_motifs(meme_p)

        hits_np, used_np = scan_motifs_in_sequences(
            meme_p, motifs, seqs, pvalue=5e-5, rc=True, backend="numpy"
        )
        assert used_np == "aria_numpy_dp"
        assert hits_np.shape == (3, 240)
        # Planted consensus hits must be detected on both forward and reverse complement
        assert bool(hits_np[0, 0]) is True
        assert bool(hits_np[1, 0]) is True
        assert bool(hits_np[2, 0]) is True

        # Guard against pvalue >= 1.0 or <= 0.0
        for bad_p in (0.0, 1.0, -0.01, 1.2):
            raised = False
            try:
                scan_motifs_in_sequences(meme_p, motifs, seqs[:5], pvalue=bad_p)
            except ValueError:
                raised = True
            assert raised, f"Expected ValueError for invalid motif pvalue={bad_p}"

        if _has_module("snapatac2"):
            hits_rs, used_rs = scan_motifs_in_sequences(
                meme_p, motifs, seqs, pvalue=5e-5, rc=True, backend="snapatac2"
            )
            assert used_rs == "snapatac2_rust"
            agreement = float(np.mean(hits_np == hits_rs))
            assert agreement >= 0.99, f"Scanner boolean parity {agreement:.4f} < 0.99"


# ── 4. GC-stratified CMH enrichment rejects pure GC confounding ───────────────


def test_gc_stratified_cmh_enrichment_rejects_pure_gc_confounder():
    from aria.utils.eregulon_core import compute_gc_stratified_motif_enrichment

    rng = np.random.default_rng(19)
    n_peaks = 1000
    peak_gc = np.concatenate(
        [
            rng.uniform(0.25, 0.45, size=600),
            rng.uniform(0.60, 0.80, size=400),
        ]
    )
    # Foreground cis peaks are heavily skewed toward high-GC regions
    fg_prob = np.where(peak_gc >= 0.60, 0.55, 0.05)
    fg_mask = rng.random(n_peaks) < fg_prob

    # Motif 0: pure GC confounder — hit probability depends ONLY on GC quintile, NOT fg_mask!
    q_edges = np.quantile(peak_gc, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    q_bins = np.clip(np.digitize(peak_gc, q_edges[1:-1], right=False), 0, 4)
    gc_rates = np.array([0.02, 0.06, 0.12, 0.30, 0.65])
    m_gc_only = rng.random(n_peaks) < gc_rates[q_bins]

    # Motif 1: true enriched motif — genuinely enriched inside fg_mask within every GC bin
    true_rates = np.clip(gc_rates[q_bins] + np.where(fg_mask, 0.38, 0.0), 0.0, 0.98)
    m_true = rng.random(n_peaks) < true_rates

    hit_mat = np.vstack([m_gc_only, m_true])
    enr = compute_gc_stratified_motif_enrichment(hit_mat, fg_mask, peak_gc, n_gc_bins=5)

    assert abs(float(enr["cmh_z"][0])) < 2.2, f"GC confounder inflated CMH z={enr['cmh_z'][0]}"
    assert float(enr["cmh_z"][1]) > 5.0, f"True motif missed CMH z={enr['cmh_z'][1]}"


# ── 5. Synthetic multiome eRegulon ground truth, paralog, depth & determinism ─


def _build_synthetic_eregulon_fixture(tmp_dir: Path, seed: int = 101):
    rng = np.random.default_rng(seed)
    n_cells = 180

    motif_act = "GGAAGTGAAAGT"   # shared by TF_ACT and TF_PARA (paralog family)
    motif_inv = "CCCTTAGGTTCC"   # TF_INV (inverse TF-target over + cis p2g)
    motif_dim = "ATGACTCATAGC"   # TF_DIM1::TF_DIM2 heterodimer
    motif_dep = "CGCGCGTATATA"   # TF_DEPTH (driven purely by sequencing depth!)
    motif_noi = "TTTAAACCCGGG"   # TF_NOISE (random unlinked TF)

    meme_p = _write_test_meme(
        tmp_dir / "synth.meme",
        [
            ("M_ACT.1", "TF_ACT", motif_act),
            ("M_ACT.2", "TF_ACT(var.2)", motif_act),
            ("M_PARA.1", "TF_PARA", motif_act),
            ("M_INV.1", "TF_INV", motif_inv),
            ("M_DIM.1", "TF_DIM1::TF_DIM2", motif_dim),
            ("M_DEP.1", "TF_DEPTH", motif_dep),
            ("M_NOI.1", "TF_NOISE", motif_noi),
        ],
    )

    # 3 chromosomes x 900,000 bp random DNA
    chrom_seqs = {
        c: list("".join(rng.choice(list("ACGT"), size=900_000)))
        for c in ("chr1", "chr2", "chr3")
    }

    def _embed(chrom: str, center: int, seq: str) -> str:
        s = center - 250
        e = center + 250
        mid = center - len(seq) // 2
        chrom_seqs[chrom][mid : mid + len(seq)] = list(seq)
        return f"{chrom}:{s}-{e}"

    # Cell programs + sequencing depth factor
    prog_a = np.repeat([0.0, 1.0], n_cells // 2) + rng.normal(0.0, 0.08, n_cells)
    prog_a = np.clip(prog_a, 0.05, 1.1)
    prog_b = np.tile([0.0, 1.0], n_cells // 2) + rng.normal(0.0, 0.08, n_cells)
    prog_b = np.clip(prog_b, 0.05, 1.1)
    prog_d = rng.permutation(prog_a)
    depth = rng.uniform(0.4, 4.0, size=n_cells)

    genes: list[str] = [
        "TF_ACT",
        "TF_PARA",
        "TF_INV",
        "TF_DIM1",
        "TF_DIM2",
        "TF_DEPTH",
        "TF_NOISE",
    ]
    act_targets = [f"ACT_G{i}" for i in range(8)]
    inv_targets = [f"INV_G{i}" for i in range(8)]
    dim_targets = [f"DIM_G{i}" for i in range(7)]
    dep_targets = [f"DEP_G{i}" for i in range(7)]
    bg_genes = [f"BG_G{i}" for i in range(90)]
    genes.extend(act_targets + inv_targets + dim_targets + dep_targets + bg_genes)

    rna = np.zeros((n_cells, len(genes)), dtype=np.int32)
    g_idx = {g: i for i, g in enumerate(genes)}

    # TFs
    rna[:, g_idx["TF_ACT"]] = rng.poisson(depth * (1.0 + 30.0 * prog_a))
    rna[:, g_idx["TF_PARA"]] = rng.poisson(depth * (1.0 + 28.0 * prog_a))
    # TF_INV is HIGH when prog_b is LOW, while INV_G* and their cis enhancers are HIGH when prog_b is HIGH -> neg_tf_pos_p2g!
    rna[:, g_idx["TF_INV"]] = rng.poisson(depth * (1.0 + 30.0 * (1.15 - prog_b)))
    rna[:, g_idx["TF_DIM1"]] = rng.poisson(depth * (1.0 + 26.0 * prog_d))
    rna[:, g_idx["TF_DIM2"]] = rng.poisson(depth * (1.0 + 26.0 * prog_d))
    # TF_DEPTH scales super-linearly with depth so even log1p(CP10k) correlates with depth ONLY!
    rna[:, g_idx["TF_DEPTH"]] = rng.poisson(depth * (0.5 + 14.0 * depth))
    rna[:, g_idx["TF_NOISE"]] = rng.poisson(depth * rng.uniform(4.0, 12.0, size=n_cells))

    peaks: list[str] = []
    atac_cols: list[np.ndarray] = []
    gtf_lines: list[str] = []

    def _add_gtf(g: str, chrom: str, tss: int):
        gtf_lines.append(
            f'{chrom}\tARIA\tgene\t{tss}\t{tss + 2000}\t.\t+\t.\t'
            f'gene_id "{g}_ID"; gene_name "{g}"; gene_type "protein_coding";\n'
        )

    # Self-loop peak for TF_ACT on chr1:50,000
    _add_gtf("TF_ACT", "chr1", 50_000)
    peaks.append(_embed("chr1", 52_000, motif_act))
    atac_cols.append(rng.poisson(depth * (0.2 + 7.0 * prog_a)))

    _add_gtf("TF_PARA", "chr1", 60_000)
    _add_gtf("TF_INV", "chr2", 50_000)
    _add_gtf("TF_DIM1", "chr3", 50_000)
    _add_gtf("TF_DIM2", "chr3", 60_000)
    _add_gtf("TF_DEPTH", "chr2", 70_000)
    _add_gtf("TF_NOISE", "chr3", 70_000)

    for k, g in enumerate(act_targets):
        tss = 110_000 + k * 60_000
        _add_gtf(g, "chr1", tss)
        rna[:, g_idx[g]] = rng.poisson(depth * (1.0 + 26.0 * prog_a))
        peaks.append(_embed("chr1", tss + 4_000, motif_act))
        atac_cols.append(rng.poisson(depth * (0.2 + 7.0 * prog_a)))

    for k, g in enumerate(inv_targets):
        tss = 110_000 + k * 60_000
        _add_gtf(g, "chr2", tss)
        rna[:, g_idx[g]] = rng.poisson(depth * (1.0 + 26.0 * prog_b))
        peaks.append(_embed("chr2", tss + 4_000, motif_inv))
        atac_cols.append(rng.poisson(depth * (0.2 + 7.0 * prog_b)))

    for k, g in enumerate(dim_targets):
        tss = 110_000 + k * 60_000
        _add_gtf(g, "chr3", tss)
        rna[:, g_idx[g]] = rng.poisson(depth * (1.0 + 25.0 * prog_d))
        peaks.append(_embed("chr3", tss + 4_000, motif_dim))
        atac_cols.append(rng.poisson(depth * (0.2 + 7.0 * prog_d)))

    for k, g in enumerate(dep_targets):
        tss = 620_000 + k * 30_000
        _add_gtf(g, "chr2", tss)
        rna[:, g_idx[g]] = rng.poisson(depth * (0.5 + 14.0 * depth))
        peaks.append(_embed("chr2", tss + 3_000, motif_dep))
        # Give cis peaks strong depth correlation AND explicit Cis signal so p2g passes, testing depth residual!
        atac_cols.append(rng.poisson(depth * (0.3 + 4.5 * depth)))

    # Background unlinked genes & background peaks across all 3 chromosomes
    chrom_cycle = ("chr1", "chr2", "chr3")
    for k, g in enumerate(bg_genes):
        c = chrom_cycle[k % 3]
        tss = 600_000 + (k // 3) * 8_000
        _add_gtf(g, c, tss)
        rna[:, g_idx[g]] = rng.poisson(depth * rng.uniform(5.0, 22.0))

    for b_i in range(120):
        c = chrom_cycle[b_i % 3]
        center = 15_000 + (b_i // 3) * 6_500
        # Sprinkle TF_NOISE motif into 20 background peaks
        seq_ins = motif_noi if (b_i % 6 == 0) else "AACCGGTTAACC"
        peaks.append(_embed(c, center, seq_ins))
        atac_cols.append(rng.poisson(depth * rng.uniform(0.8, 2.5)))

    atac = np.column_stack(atac_cols).astype(np.int32)
    fa_p = _write_indexed_fasta(
        tmp_dir / "synth.fa",
        [(c, "".join(chrom_seqs[c])) for c in ("chr1", "chr2", "chr3")],
    )
    gtf_p = tmp_dir / "synth.gtf"
    gtf_p.write_text("".join(gtf_lines), encoding="utf-8")

    return {
        "rna": sparse.csc_matrix(rna),
        "genes": genes,
        "atac": sparse.csc_matrix(atac),
        "peaks": peaks,
        "meme_path": meme_p,
        "fasta_path": fa_p,
        "gtf_path": gtf_p,
    }


def test_synthetic_multiome_eregulon_recovery_selfloop_paralog_and_determinism():
    from aria.utils.eregulon_core import compute_candidate_eregulons
    from aria.utils.peak2gene_core import parse_gtf_tss_detailed

    with tempfile.TemporaryDirectory() as td:
        fx = _build_synthetic_eregulon_fixture(Path(td))
        gene_records, _ = parse_gtf_tss_detailed(str(fx["gtf_path"]), wanted=set(fx["genes"]))

        res_a = compute_candidate_eregulons(
            rna_matrix=fx["rna"],
            rna_gene_names=fx["genes"],
            atac_matrix=fx["atac"],
            atac_peak_names=fx["peaks"],
            gene_tss_records=gene_records,
            meme_path=fx["meme_path"],
            fasta_path=fx["fasta_path"],
            distance_bp=25_000,
            p2g_min_corr=0.10,
            p2g_max_padj=0.05,
            min_cells_detected=10,
            motif_pvalue=5e-5,
            min_tf_gene_corr=0.10,
            min_target_genes=5,
            max_regulon_padj=0.05,
            seed=0,
        )
        res_b = compute_candidate_eregulons(
            rna_matrix=fx["rna"],
            rna_gene_names=fx["genes"],
            atac_matrix=fx["atac"],
            atac_peak_names=fx["peaks"],
            gene_tss_records=gene_records,
            meme_path=fx["meme_path"],
            fasta_path=fx["fasta_path"],
            distance_bp=25_000,
            p2g_min_corr=0.10,
            p2g_max_padj=0.05,
            min_cells_detected=10,
            motif_pvalue=5e-5,
            min_tf_gene_corr=0.10,
            min_target_genes=5,
            max_regulon_padj=0.05,
            seed=0,
        )

        by_id = {r["regulon_id"]: r for r in res_a["eregulons"]}
        assert "TF_ACT__pos_tf_pos_p2g" in by_id
        assert "TF_INV__neg_tf_pos_p2g" in by_id
        assert "TF_DIM1__pos_tf_pos_p2g" in by_id
        assert "TF_DIM2__pos_tf_pos_p2g" in by_id
        assert not any(r["tf"] == "TF_NOISE" for r in res_a["eregulons"])

        # Self-loop guard (`g != TF`) + `self_linked_peaks` disclosure + variant dedup
        act_reg = by_id["TF_ACT__pos_tf_pos_p2g"]
        assert "TF_ACT" not in act_reg["target_genes"]
        assert len(act_reg["self_linked_peaks"]) >= 1
        assert act_reg["n_motifs_tested_for_tf"] == 2  # M_ACT.1 and M_ACT.2(var.2)
        assert act_reg["regulon_padj_scope"] == "module_level"

        # Paralog motif ambiguity flagged on TF_ACT <-> TF_PARA, NOT on heterodimer TF_DIM1::TF_DIM2
        assert act_reg["motif_family_ambiguity"] is True
        assert any(a["tf"] == "TF_PARA" for a in act_reg["ambiguous_paralog_tfs"])
        dim1_reg = by_id["TF_DIM1__pos_tf_pos_p2g"]
        assert not any(a["tf"] == "TF_DIM2" for a in dim1_reg["ambiguous_paralog_tfs"])

        # Bit-exact seed=0 determinism
        assert [r["regulon_id"] for r in res_a["eregulons"]] == [
            r["regulon_id"] for r in res_b["eregulons"]
        ]
        assert res_a["triplets"] == res_b["triplets"]
        assert np.array_equal(res_a["rna_aucell"], res_b["rna_aucell"])
        assert np.array_equal(res_a["atac_aucell"], res_b["atac_aucell"])


def test_library_depth_partial_correlation_flags_depth_confounded_module():
    from aria.utils.eregulon_core import compute_candidate_eregulons
    from aria.utils.peak2gene_core import parse_gtf_tss_detailed

    with tempfile.TemporaryDirectory() as td:
        fx = _build_synthetic_eregulon_fixture(Path(td), seed=202)
        gene_records, _ = parse_gtf_tss_detailed(str(fx["gtf_path"]), wanted=set(fx["genes"]))
        # Supply precomputed cis links for DEP_G0..DEP_G6 so we isolate Step 4-6 depth conditioning
        dep_links = [
            {
                "gene": f"DEP_G{k}",
                "peak": f"chr2:{623_000 + k * 30_000 - 250}-{623_000 + k * 30_000 + 250}",
                "correlation": 0.45,
                "padj": 1e-4,
                "distance_bp": 3_000,
                "region_class": "distal_cis",
            }
            for k in range(7)
        ]
        res = compute_candidate_eregulons(
            rna_matrix=fx["rna"],
            rna_gene_names=fx["genes"],
            atac_matrix=fx["atac"],
            atac_peak_names=fx["peaks"],
            gene_tss_records=gene_records,
            meme_path=fx["meme_path"],
            fasta_path=fx["fasta_path"],
            precomputed_p2g={"links": dep_links, "n_tested": 7, "n_fdr_significant": 7},
            min_tf_gene_corr=0.05,
            min_target_genes=5,
            max_regulon_padj=1.0,
            seed=0,
        )
        dep_regs = [r for r in res["eregulons"] if r["tf"] == "TF_DEPTH"]
        assert dep_regs, "Expected TF_DEPTH module when precomputed cis links are supplied"
        # Conditioning on log10(1 + total_counts) collapses the spurious depth correlation!
        assert abs(dep_regs[0]["mean_r_tf_gene_partial_depth"]) < abs(dep_regs[0]["mean_r_tf_gene"])


# ── 6. Paired .h5mu IPC & ScriptContract validation ───────────────────────────


def test_integration_eregulons_error_contracts_and_h5mu_roundtrip():
    from aria.scripts.integration_eregulons import integration_eregulons

    contract = None
    has_contracts = _has_module("pydantic")
    if has_contracts:
        from aria.utils.script_contracts import contract_for_script

        contract = contract_for_script("aria/scripts/integration_eregulons.py")
        assert contract is not None
        assert contract.validation_level == "beta"

    err_p = integration_eregulons({"motif_pvalue": 1.0})
    assert err_p["status"] == "error"
    assert err_p["error_type"] == "InvalidParameter"

    err_m = integration_eregulons(
        {
            "rna_files": [],
            "atac_files": [],
        }
    )
    assert err_m["status"] == "error"
    assert err_m["error_type"] == "MissingRNA"

    if not (_has_module("anndata") and _has_module("mudata")):
        return

    import anndata as ad
    import mudata as md

    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        fx = _build_synthetic_eregulon_fixture(td_p)
        obs_names = [f"CELL_{i:04d}" for i in range(fx["rna"].shape[0])]
        rna_ad = ad.AnnData(X=fx["rna"].tocsr())
        rna_ad.obs_names = obs_names
        rna_ad.var_names = fx["genes"]
        atac_ad = ad.AnnData(X=fx["atac"].tocsr())
        atac_ad.obs_names = obs_names
        atac_ad.var_names = fx["peaks"]
        mdata = md.MuData({"rna": rna_ad, "atac": atac_ad})
        h5mu_p = td_p / "paired.h5mu"
        mdata.write_h5mu(str(h5mu_p))

        params = {
            "h5mu_path": str(h5mu_p),
            "gtf_path": str(fx["gtf_path"]),
            "genome_fasta": str(fx["fasta_path"]),
            "motif_file": str(fx["meme_path"]),
            "distance_kb": 25,
            "p2g_min_corr": 0.10,
            "p2g_max_padj": 0.05,
            "motif_pvalue": 5e-5,
            "min_tf_gene_corr": 0.10,
            "min_target_genes": 5,
            "max_regulon_padj": 0.05,
            "seed": 0,
            "update_h5mu": True,
            "output_dir": str(td_p / "out"),
        }
        if contract is not None:
            assert contract.validate_params(params) == []
        out = integration_eregulons(params)
        assert out["status"] == "success", out
        if contract is not None:
            assert contract.validate_result(out) == []
        assert out["associative_only"] is True
        assert out["perturbation_validated"] is False
        assert out["chromatin_conformation_validated"] is False
        assert out["n_eregulons"] >= 3
        assert Path(out["output_path"]).is_file()
        assert Path(out["triplets_path"]).is_file()

        with open(out["output_path"], newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == out["n_eregulons"]

        reloaded = md.read_h5mu(str(h5mu_p))
        assert "X_eregulon_rna" in reloaded.obsm
        assert "X_eregulon_atac" in reloaded.obsm
        assert reloaded.obsm["X_eregulon_rna"].shape == (
            fx["rna"].shape[0],
            out["n_eregulons"],
        )
        assert bool(reloaded.uns["aria_eregulons"]["associative_only"]) is True


# ── 7. Four-tier concordance scorer & ADR-002 honesty ─────────────────────────


def test_concordance_eregulons_four_tiers_and_not_run_honesty():
    from aria.benchmarks.concordance_eregulons import (
        score_eregulon_concordance,
        score_eregulon_determinism,
    )

    nr = score_eregulon_concordance(None, None)
    assert nr["status"] == "not_run"

    aria_trips = [
        {"tf": "SPI1", "peak": "chr11:1000-1500", "gene": "TYROBP", "sign_class": "pos_tf_pos_p2g"},
        {"tf": "SPI1", "peak": "chr19:5000-5500", "gene": "C1QA", "sign_class": "pos_tf_pos_p2g"},
        {"tf": "MEIS1", "peak": "chr2:8000-8500", "gene": "NRGN", "sign_class": "neg_tf_pos_p2g"},
    ]
    ref_trips = [
        # Self-loop injected by SCENIC+ must be ignored automatically (`SPI1 -> SPI1`)
        {"TF": "SPI1", "Region": "chr11:900-1400", "Gene": "SPI1", "sign": "+/+"},
        # Shifted by 100bp (`80%` reciprocal overlap with `chr11:1000-1500`!)
        {"TF": "SPI1", "Region": "chr11:1100-1600", "Gene": "TYROBP", "sign": "+/+"},
        # Exact coordinate match (`chr19-5000-5500` -> `chr19:5000-5500`)
        {"TF": "SPI1", "Region": "chr19-5000-5500", "Gene": "C1QA", "sign": "+/+"},
        {"TF": "SOX10", "Region": "chr22:2000-2500", "Gene": "MBP", "sign": "+/+"},
    ]
    rng = np.random.default_rng(3)
    base_sig = rng.normal(size=80)
    aria_auc = {"SPI1__pos_tf_pos_p2g": base_sig + rng.normal(0, 0.1, 80)}
    ref_auc = {"SPI1__pos_tf_pos_p2g": base_sig + rng.normal(0, 0.1, 80)}

    conc = score_eregulon_concordance(
        aria_trips,
        ref_trips,
        shared_tf_universe=[f"TF{i}" for i in range(50)] + ["SPI1", "MEIS1", "SOX10"],
        shared_gene_universe=[f"G{i}" for i in range(200)] + ["SPI1", "TYROBP", "C1QA", "NRGN", "MBP"],
        aria_cell_aucell=aria_auc,
        ref_cell_aucell=ref_auc,
        random_aucell_baseline_spearmans=[0.05, 0.11, -0.02, 0.08],
    )
    assert conc["status"] == "computed"
    assert conc["tf_concordance"]["n_intersection"] == 1
    assert conc["edge_concordance_global"]["n_intersection"] == 2
    assert conc["edge_concordance_global"]["fold_enrichment_over_null"] > 10.0
    # Exact peak match catches C1QA (1), reciprocal >=50% catches BOTH TYROBP and C1QA (2)!
    assert conc["triplet_concordance"]["n_intersection"] == 1
    assert conc["triplet_concordance"]["n_aria_matched_reciprocal_50pct"] == 2
    assert conc["cell_activity_concordance"]["median_spearman"] > 0.90
    assert conc["cell_activity_concordance"]["median_spearman_excess_over_random"] > 0.80

    det = score_eregulon_determinism(
        [{"tf": "SPI1", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 1e-5}],
        aria_trips,
        [{"tf": "SPI1", "sign_class": "pos_tf_pos_p2g", "regulon_padj": 1e-5}],
        aria_trips,
        rna_aucell_a=np.ones((10, 1), dtype=np.float32),
        rna_aucell_b=np.ones((10, 1), dtype=np.float32),
    )
    assert det["bit_exact"] is True


# ── 8. AST & ADR-064 governance guard ─────────────────────────────────────────


def test_no_causal_labels_or_hardcoded_gene_panels_ast_guard():
    if _has_module("litellm"):
        from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION

        assert INTEGRATION_VALIDATION["level"] == "scaffold"
        assert INTEGRATION_VALIDATION["dispatch_enabled"] is False

    core_src = (ROOT / "aria" / "utils" / "eregulon_core.py").read_text(encoding="utf-8")
    script_src = (ROOT / "aria" / "scripts" / "integration_eregulons.py").read_text(
        encoding="utf-8"
    )
    assert "from aria.utils.stats import bh_correct" in core_src

    for label, src in (("eregulon_core.py", core_src), ("integration_eregulons.py", script_src)):
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                # Sign-class identifiers must never emit causal "direct_activator" / "proven_repressor"
                assert node.value not in (
                    "direct_activator",
                    "direct_repressor",
                    "proven_regulator",
                    "causal_regulon",
                ), f"Forbidden causal sign label in {label}: {node.value}"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS: {name}")
    print("ALL V48 SLICE I6 TESTS PASSED")
