"""B2 (bulk ATAC pre-print): genomic annotation of DA peaks.

Pure, deterministic peak annotation against a GTF gene model — nearest-TSS gene,
signed strand-aware distance-to-TSS, and a feature class (Promoter/Exonic/Intronic/
Distal Intergenic). The dispatchable entry auto-resolves a GTF and honest-skips when
none is found (annotation is a positional fact; without a gene model there is nothing
to assert)."""

from __future__ import annotations

from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

# A minimal two-gene GTF: GENEA on +, GENEB on -.
_GTF = "\n".join([
    # GENEA (+), TSS = 1000, body 1000-5000, exons 1000-1200/3000-3200/4800-5000
    'chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG_A"; gene_name "GENEA";',
    'chr1\ttest\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG_A"; gene_name "GENEA";',
    'chr1\ttest\texon\t3000\t3200\t.\t+\t.\tgene_id "ENSG_A"; gene_name "GENEA";',
    'chr1\ttest\texon\t4800\t5000\t.\t+\t.\tgene_id "ENSG_A"; gene_name "GENEA";',
    # GENEB (-), TSS = 25000, body 20000-25000, exons 20000-20200/24800-25000
    'chr1\ttest\tgene\t20000\t25000\t.\t-\t.\tgene_id "ENSG_B"; gene_name "GENEB";',
    'chr1\ttest\texon\t20000\t20200\t.\t-\t.\tgene_id "ENSG_B"; gene_name "GENEB";',
    'chr1\ttest\texon\t24800\t25000\t.\t-\t.\tgene_id "ENSG_B"; gene_name "GENEB";',
    "",
])


@pytest.fixture()
def gtf(tmp_path) -> Path:
    p = tmp_path / "annotation.gtf"
    p.write_text(_GTF)
    return p


def test_parse_gtf_gene_models(gtf):
    from aria.scripts.chromatin_peak_annotation import parse_gtf_gene_models

    models = parse_gtf_gene_models(str(gtf))
    # contigs are stored normalized (chr1 -> 1) so chr-prefixed peaks match.
    assert set(models) == {"1"}
    by_name = {g["gene"]: g for g in models["1"]}
    assert by_name["GENEA"]["tss"] == 1000 and by_name["GENEA"]["strand"] == "+"
    # minus-strand TSS is the gene END.
    assert by_name["GENEB"]["tss"] == 25000 and by_name["GENEB"]["strand"] == "-"
    assert (3000, 3200) in by_name["GENEA"]["exons"]
    # per-chrom list is TSS-sorted.
    assert [g["tss"] for g in models["1"]] == sorted(
        g["tss"] for g in models["1"])


def test_annotate_peaks_feature_classes_and_distance(gtf):
    from aria.scripts.chromatin_peak_annotation import (
        annotate_peaks, parse_gtf_gene_models)

    models = parse_gtf_gene_models(str(gtf))
    peaks = [
        "chr1:1050-1150",    # center 1100 -> Promoter of GENEA (+100 from TSS)
        "chr1:3050-3150",    # center 3100 -> within an exon -> Exonic
        "chr1:3500-3600",    # center 3550 -> in body, no exon -> Intronic
        "chr1:50000-50100",  # far -> Distal Intergenic
        "chr1:24900-25000",  # center 24950 -> Promoter of GENEB (- strand, +50)
        "not_a_peak",        # unparseable -> skipped + counted
    ]
    res = annotate_peaks(peaks, models)
    by_peak = {r["peak"]: r for r in res["rows"]}

    assert by_peak["chr1:1050-1150"]["feature"] == "Promoter"
    assert by_peak["chr1:1050-1150"]["nearest_gene"] == "GENEA"
    assert by_peak["chr1:1050-1150"]["distance_to_tss"] == 100
    assert by_peak["chr1:3050-3150"]["feature"] == "Exonic"
    assert by_peak["chr1:3500-3600"]["feature"] == "Intronic"
    assert by_peak["chr1:50000-50100"]["feature"] == "Distal Intergenic"
    # minus-strand promoter: signed distance is TSS-center = 25000-24950 = 50.
    assert by_peak["chr1:24900-25000"]["feature"] == "Promoter"
    assert by_peak["chr1:24900-25000"]["distance_to_tss"] == 50

    assert res["n_unparsed"] == 1
    assert res["n_annotated"] == 5
    assert res["feature_distribution"]["Promoter"] == 2


def test_chrom_naming_reconciled_ensembl_gtf_vs_chr_peaks(tmp_path):
    """Ensembl GTF contigs (`1`/`X`) must match UCSC chr-prefixed ATAC peaks
    (`chr1`/`chrX`) — otherwise real ENCODE peaks annotate as all-intergenic."""
    from aria.scripts.chromatin_peak_annotation import (
        annotate_peaks, parse_gtf_gene_models)

    ens = "\n".join([
        # Ensembl-style contig "7" (no chr prefix).
        '7\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG_X"; gene_name "GENEX";',
        '7\ttest\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG_X"; gene_name "GENEX";',
        "",
    ])
    p = tmp_path / "ens.gtf"
    p.write_text(ens)
    models = parse_gtf_gene_models(str(p))
    assert "7" in models  # stored normalized
    # chr-prefixed peak still resolves to the gene.
    res = annotate_peaks(["chr7:1050-1150"], models)
    assert res["rows"][0]["nearest_gene"] == "GENEX"
    assert res["rows"][0]["feature"] == "Promoter"


def test_resolve_gtf_prefers_experiment_assembly(tmp_path, monkeypatch):
    """Auto-resolution must prefer ~/.aria/genomes/<genome>/ for THIS experiment's
    assembly, not the first genome dir found (a human run must not get a mouse GTF)."""
    from aria.scripts import chromatin_peak_annotation as mod

    home = tmp_path / "home"
    hg = home / ".aria" / "genomes" / "hg38"
    hg.mkdir(parents=True)
    (hg / "annotation.gtf").write_bytes(b"x" * 200_000)
    monkeypatch.delenv("ARIA_GTF", raising=False)
    monkeypatch.setattr(mod.Path, "home", classmethod(lambda cls: home))

    resolved = mod._resolve_gtf({"genome": "hg38"}, tmp_path / "out")
    assert resolved == hg / "annotation.gtf"


def test_chromatin_peak_annotation_honest_skip_without_gtf(tmp_path, monkeypatch):
    from aria.scripts.chromatin_peak_annotation import chromatin_peak_annotation

    monkeypatch.delenv("ARIA_GTF", raising=False)
    # Force "no GTF locatable" regardless of the machine's ~/.aria/genomes
    # (auto-resolution finds a real GTF on Samael's box — that path is exercised
    # by the end-to-end test below).
    monkeypatch.setattr(
        "aria.scripts.chromatin_peak_annotation._locate_gtf",
        lambda *a, **k: None)
    # No gtf hint, no ARIA_GTF, and no locatable GTF -> honest skip.
    res = chromatin_peak_annotation({
        "data_type": "bulk_ATAC",
        "comparisons": [{"test": "A", "reference": "B",
                         "full_results_csv": "nope.csv"}],
        "output_dir": str(tmp_path / "ann"),
    })
    assert res["status"] == "skipped" and res["ran"] is False
    assert res["reason"] == "no_gtf_for_peak_annotation"


def test_chromatin_peak_annotation_end_to_end(tmp_path, gtf):
    from aria.scripts.chromatin_peak_annotation import chromatin_peak_annotation

    da_csv = tmp_path / "bulk_da_A_vs_B_full.csv"
    pd.DataFrame({
        "peak": ["chr1:1050-1150", "chr1:3500-3600", "chr1:50000-50100",
                 "chr1:9000-9100"],
        "log2FoldChange": [2.0, -1.5, 1.2, 0.1],
        "padj": [0.001, 0.002, 0.01, 0.9],
        "comparison": ["A_vs_B"] * 4,
        "significant": [True, True, True, False],  # last peak is NOT significant
    }).to_csv(da_csv, index=False)

    res = chromatin_peak_annotation({
        "data_type": "bulk_ATAC",
        "comparisons": [{"test": "A", "reference": "B",
                         "full_results_csv": str(da_csv)}],
        "gtf": str(gtf),
        "output_dir": str(tmp_path / "ann"),
    })
    assert res["status"] == "success" and res["ran"] is True
    dist = res["feature_distribution_overall"]
    # Only the 3 significant peaks are annotated (the non-sig one is excluded).
    assert sum(dist.values()) == 3
    assert dist.get("Promoter") == 1 and dist.get("Intronic") == 1
    assert dist.get("Distal Intergenic") == 1
    comp = res["comparisons"][0]
    assert comp["n_annotated"] == 3
    ann_csv = comp["annotation_csv"]
    assert ann_csv and Path(ann_csv).exists()
    written = pd.read_csv(ann_csv)
    assert list(written.columns) == ["peak", "nearest_gene",
                                     "distance_to_tss", "feature"]


def test_render_peak_annotation_dual_and_honest(tmp_path):
    from aria.agents._narrative_chromatin import render_peak_annotation

    ann = {"feature_distribution_overall": {
        "Promoter": 120, "Exonic": 30, "Intronic": 200, "Distal Intergenic": 90}}
    path = render_peak_annotation(ann, tmp_path / "ann.png")
    assert path and Path(path).exists()
    assert Path(path).with_suffix(".svg").exists()
    # No distribution / not-run -> honest None.
    assert render_peak_annotation({}, tmp_path / "none.png") is None
    assert render_peak_annotation(
        {"feature_distribution_overall": {}}, tmp_path / "none2.png") is None


def test_generate_figures_attaches_bulk_atac_peak_annotation(tmp_path):
    from aria.agents import _narrative_chromatin

    findings = {
        "differential_accessibility": {"data_type": "bulk_ATAC",
                                       "comparisons": []},
        "peak_annotation": {
            "ran": True,
            "feature_distribution_overall": {"Promoter": 10, "Intronic": 5}},
    }
    out = _narrative_chromatin.generate_figures(
        findings, None, tmp_path / "figures", env_manager=None)
    assert "bulk_atac_peak_annotation" in out["figures"]


def test_peak_annotation_narrator_block():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    ann = {
        "status": "success", "ran": True, "data_type": "bulk_ATAC",
        "analysis": "peak_annotation", "validation_level": "beta",
        "method": "nearest-TSS genomic annotation",
        "gtf": "annotation.gtf",
        "promoter_upstream": 2000, "promoter_downstream": 2000,
        "feature_distribution_overall": {
            "Promoter": 25, "Exonic": 10, "Intronic": 40,
            "Distal Intergenic": 25},
        "comparisons": [{"test": "A", "reference": "B",
                         "annotation_csv": "/tmp/ann.csv"}],
    }
    blocks = ChromatinNarrator().collect(
        "chromatin_agent",
        {"findings": {"peak_annotation": ann}})
    ids = {b.id for b in blocks}
    assert "chromatin.peak_annotation" in ids
    block = next(b for b in blocks if b.id == "chromatin.peak_annotation")
    assert block.status == "success"
    # Annotation must stay positional/associative, never a regulatory claim.
    assert "regulate" not in block.claim.lower()
    assert any("positional" in c.text.lower() for c in block.caveats)
