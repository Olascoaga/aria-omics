"""D8 convergence regression tests for genomic coordinates, contigs, and GTF parsers (D8-004..D8-007, D8-011)."""

from __future__ import annotations

from pathlib import Path

import pytest


def test_d8_006_contig_normalization_copies_converge():
    """D8-006: all 8 contig normalizers must reconcile UCSC/Ensembl aliases and preserve unplaced scaffolds."""
    from aria.scripts.chromatin_motifs import _resolve_contig_in_valid
    from aria.scripts.chromatin_peak_annotation import _norm_chrom as ann_norm
    from aria.scripts.chromatin_peak_counts import _norm_contig as cnt_norm
    from aria.scripts.chromatin_qc import _norm_contig as qc_norm
    from aria.scripts.hic_figure import _chrom_aliases as hic_fig_aliases
    from aria.scripts.hic_topology import _chrom_aliases as hic_aliases
    from aria.utils.cross_modal_concordance import _norm_chrom as xm_norm
    from aria.utils.peak2gene_core import canonical_contig as p2g_canon

    pairs = [("chrMt", "MT"), ("chrm", "chrM"), ("chrx", "X"), (" chr1 ", "1")]
    funcs = {
        "peak2gene_core.canonical_contig": p2g_canon,
        "chromatin_peak_annotation._norm_chrom": ann_norm,
        "cross_modal_concordance._norm_chrom": xm_norm,
        "chromatin_peak_counts._norm_contig": cnt_norm,
        "chromatin_qc._norm_contig": qc_norm,
    }
    failures: list[str] = []
    for a, b in pairs:
        assert p2g_canon(a) == p2g_canon(b)
        for name, fn in funcs.items():
            if fn(a) != fn(b):
                failures.append(f"{name}: {a!r}->{fn(a)!r} != {b!r}->{fn(b)!r}")

    scaff = "chrUn_KI270442v1"
    assert p2g_canon(scaff) == scaff
    for name, fn in funcs.items():
        if fn(scaff) != scaff:
            failures.append(f"{name} mangled {scaff!r} -> {fn(scaff)!r}")

    if _resolve_contig_in_valid("chrx", {"chrX"}) != "chrX":
        failures.append("chromatin_motifs._resolve_contig_in_valid('chrx', {'chrX'}) -> None")
    if "chrX" not in hic_aliases("chrx"):
        failures.append(f"hic_topology._chrom_aliases('chrx')={hic_aliases('chrx')}")
    if "chrX" not in hic_fig_aliases("chrx"):
        failures.append(f"hic_figure._chrom_aliases('chrx')={hic_fig_aliases('chrx')}")
    assert not failures, "\n".join(failures)

    from aria.utils.genomic_coords import canonical_contig, contig_aliases, ensembl_contig

    assert p2g_canon is canonical_contig
    assert hic_aliases is contig_aliases
    assert hic_fig_aliases is contig_aliases
    assert cnt_norm is ensembl_contig
    assert qc_norm is ensembl_contig
    assert xm_norm is ensembl_contig


def _write_bam(path: Path, chrom: str, chrom_len: int, positions: list[tuple[int, int]]):
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": chrom, "LN": chrom_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, length) in enumerate(positions):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag = f"r{i}", "A" * length, 0
            a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
            a.cigar = [(0, length)]
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(a)
    pysam.index(str(path))


def test_d8_007_chromatin_peaks_consensus_and_qc_jsd_contig_convergence(tmp_path: Path):
    """D8-007: consensus peak merging and IP-vs-Control JSD must normalize UCSC vs Ensembl contigs."""
    pytest.importorskip("pysam")
    from aria.scripts.chromatin_peaks import _filter_blacklist, _merge_intervals_with_support
    from aria.scripts.chromatin_qc import _compute_ip_control_jsd, _frip_from_intervals

    rep1 = tmp_path / "rep1.bed"
    rep2 = tmp_path / "rep2.bed"
    bl = tmp_path / "bl.bed"
    rep1_copy = tmp_path / "rep1_c.bed"
    rep1.write_text("chr1\t1000\t2000\tp1\n", encoding="utf-8")
    rep2.write_text("1\t1000\t2000\tp1\n", encoding="utf-8")
    bl.write_text("1\t1000\t2000\tbl1\n", encoding="utf-8")
    rep1_copy.write_text("chr1\t1000\t2000\tp1\n", encoding="utf-8")

    merged = _merge_intervals_with_support([str(rep1), str(rep2)])
    max_support = max(len(m[3]) for m in merged)
    removed = _filter_blacklist(str(rep1_copy), str(bl))

    pos_ip = [(b * 5000 + k * 20, 50) for b in range(25) for k in range(30 if b == 5 else 2)]
    pos_ct = [(b * 5000 + k * 20, 50) for b in range(25) for k in range(5)]
    ip_bam = tmp_path / "ip.bam"
    ct_ucsc = tmp_path / "ct_u.bam"
    ct_ens = tmp_path / "ct_e.bam"
    _write_bam(ip_bam, "chr1", 200_000, pos_ip)
    _write_bam(ct_ucsc, "chr1", 200_000, pos_ct)
    _write_bam(ct_ens, "1", 200_000, pos_ct)

    jsd_same = _compute_ip_control_jsd([str(ip_bam)], [str(ct_ucsc)])
    jsd_mixed = _compute_ip_control_jsd([str(ip_bam)], [str(ct_ens)])
    frip_mixed, _ = _frip_from_intervals([("1", 1000, 1100)], {"chr1": [(900, 1200)]})

    assert removed == 1 and jsd_same is not None and frip_mixed == 1.0
    assert max_support == 2
    assert jsd_mixed == jsd_same


def test_d8_004_peak_coordinate_parsers_converge(tmp_path: Path):
    """D8-004: all 5 peak coordinate parsers must accept Signac, Ensembl, and unplaced scaffold peaks."""
    from aria.benchmarks.bulk_atac_marker_mapping import parse_peak_interval
    from aria.benchmarks.concordance_atac import _parse_peak as conc_parse
    from aria.scripts.benchmark_peak2gene_export import _parse_peak as exp_parse
    from aria.scripts.chromatin_peak_annotation import annotate_peaks, parse_gtf_gene_models
    from aria.utils.mudata_io import _canonical_peak
    from aria.utils.peak2gene_core import genes_to_peak_indices, parse_gtf_tss, parse_peak_coord

    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        'chr1\tHAVANA\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n'
        'chr1\tHAVANA\texon\t1000\t1600\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n',
        encoding="utf-8",
    )
    signac_peak = "chr1-1100-1500"
    assert parse_peak_coord(signac_peak) == ("chr1", 1100, 1500)
    assert genes_to_peak_indices([signac_peak], parse_gtf_tss(str(gtf)), window_bp=5000) == {"GENE1": [0]}
    assert parse_peak_coord("1:1000-1500") == ("1", 1000, 1500)
    assert parse_peak_coord("chr1_KI270706v1_random:1000-1500") == ("chr1_KI270706v1_random", 1000, 1500)

    ann = annotate_peaks([signac_peak], parse_gtf_gene_models(str(gtf)))
    assert parse_peak_interval(signac_peak) == ("chr1", 1100, 1500)
    assert ann["n_annotated"] == 1
    assert _canonical_peak(signac_peak) == "chr1:1100-1500"
    assert conc_parse("1:1000-1500") == ("1", 1000, 1500)
    assert exp_parse("chr1_KI270706v1_random:1000-1500") == ("chr1_KI270706v1_random", 1000, 1500)

    from aria.utils.genomic_coords import parse_peak_coord as canonical_parse

    assert parse_peak_coord is canonical_parse
    assert conc_parse is canonical_parse
    assert exp_parse is canonical_parse
