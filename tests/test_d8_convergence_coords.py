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


def test_d8_005_gtf_tss_parsers_converge(tmp_path: Path):
    """D8-005: all GTF gene/TSS parsers emit 0-based TSS, support transcript-only GTFs, and drop _PAR_Y."""
    from aria.scripts.chromatin_peak_annotation import parse_gtf_gene_models
    from aria.scripts.chromatin_qc import _read_local_tss
    from aria.scripts.hic_topology import _parse_gtf_tss as hic_parse_tss
    from aria.utils.peak2gene_core import parse_gtf_tss_detailed

    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG000001"; gene_name "GENE_PLUS"; gene_type "protein_coding";\n'
        'chr1\tHAVANA\tgene\t5000\t6000\t.\t-\t.\tgene_id "ENSG000002"; gene_name "GENE_MINUS"; gene_type "protein_coding";\n'
        'chrX\tHAVANA\tgene\t10000\t20000\t.\t+\t.\tgene_id "ENSG000003"; gene_name "CSF2RA"; gene_type "protein_coding";\n'
        'chrY\tHAVANA\tgene\t10000\t20000\t.\t+\t.\tgene_id "ENSG000003_PAR_Y"; gene_name "CSF2RA"; tag "PAR";\n',
        encoding="utf-8",
    )
    p2g_recs, _ = parse_gtf_tss_detailed(str(gtf))
    ann_models = parse_gtf_gene_models(str(gtf))
    qc_tss = _read_local_tss(str(gtf))
    hic_tss = hic_parse_tss(str(gtf))

    tx_gtf = tmp_path / "tx_only.gtf"
    tx_gtf.write_text(
        'chr1\tRefSeq\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "TP53"; gene_name "TP53";\n'
        'chr1\tRefSeq\texon\t1000\t1200\t.\t+\t.\tgene_id "TP53"; gene_name "TP53";\n',
        encoding="utf-8",
    )
    tx_p2g, _ = parse_gtf_tss_detailed(str(tx_gtf))
    tx_ann = parse_gtf_gene_models(str(tx_gtf))
    tx_qc = _read_local_tss(str(tx_gtf))
    tx_hic = hic_parse_tss(str(tx_gtf))

    p2g_plus = p2g_recs["GENE_PLUS"]["tss"]
    p2g_minus = p2g_recs["GENE_MINUS"]["tss"]
    ann_plus = [g["tss"] for g in ann_models["1"] if g["gene"] == "GENE_PLUS"][0]
    ann_minus = [g["tss"] for g in ann_models["1"] if g["gene"] == "GENE_MINUS"][0]
    qc_plus = qc_tss["chr1"][1][0]
    qc_minus = qc_tss["chr1"][1][1]
    hic_plus = [t for _, t, g in hic_tss if g == "GENE_PLUS"][0]
    hic_minus = [t for _, t, g in hic_tss if g == "GENE_MINUS"][0]

    assert {p2g_plus, ann_plus, qc_plus, hic_plus} == {999}
    assert {p2g_minus, ann_minus, qc_minus, hic_minus} == {5999}
    assert len(tx_p2g) == 1 and len(tx_ann) == 1 and len(tx_qc) == 1 and len(tx_hic) == 1
    assert not any(r["canon_chrom"] == "chrY" for r in p2g_recs.values())
    assert "Y" not in ann_models
    assert not any(c == "chrY" for c, _, _ in hic_tss)

    from aria.utils.genomic_coords import (
        parse_gtf_gene_models as canonical_gene_models,
        parse_gtf_tss_records as canonical_tss_records,
    )

    assert parse_gtf_tss_detailed is canonical_tss_records
    assert parse_gtf_gene_models is canonical_gene_models


def test_d8_011_gtf_symbol_map_builders_converge(tmp_path: Path):
    """D8-011: all GTF gene-ID-to-symbol builders must prefer gene over exon, support RefSeq/GFF3, and keep dotted IDs."""
    from aria.scripts.rna_bulk.gtf_io import _gtf_to_symbol_map
    from aria.scripts.rna_quantify import _build_ensembl_to_symbol_map
    from aria.utils.cross_modal_concordance import gtf_symbols

    gtf = tmp_path / "mixed.gtf"
    gtf.write_text(
        'chr1\tRefSeq\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG00000141510.17"; gene_name "TP53";\n'
        'chr1\tHAVANA\texon\t600\t700\t.\t+\t.\tgene_id "ENSG00000000003.15"; gene_name "OLD_EXON_SYM";\n'
        'chr1\tHAVANA\tgene\t600\t900\t.\t+\t.\tgene_id "ENSG00000000003.15"; gene_name "TSPAN6";\n'
        'chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tgene_id "LOC100.1"; gene_name "DOTTED_GENE";\n'
        "chr1\tGFF3\tgene\t3000\t4000\t.\t+\t.\tgene_id=ENSG00000100000.2; gene_name=BRCA1;\n",
        encoding="utf-8",
    )
    m_bulk = _gtf_to_symbol_map(str(gtf))
    m_quant = _build_ensembl_to_symbol_map(str(gtf), [])
    m_xmodal = gtf_symbols(str(gtf))
    expected = {
        "ENSG00000141510": "TP53",
        "ENSG00000000003": "TSPAN6",
        "LOC100.1": "DOTTED_GENE",
        "ENSG00000100000": "BRCA1",
    }
    assert m_bulk == expected
    assert m_quant == expected
    assert m_xmodal == expected

    from aria.utils.genomic_coords import parse_gtf_symbol_map

    assert _gtf_to_symbol_map is parse_gtf_symbol_map
    assert gtf_symbols is parse_gtf_symbol_map
