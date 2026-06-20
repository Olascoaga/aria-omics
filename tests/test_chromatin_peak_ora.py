"""B3 (bulk ATAC pre-print): functional ORA of genes near DA peaks.

Reuses the B2 peak→gene assignment + the bulk RNA ORA engine (local hypergeometric
over versioned GMTs). Foreground = genes near significant peaks; background = genes
near all tested peaks. Honest-skip without a GTF / significant genes; enrichment is
associative, never a regulatory claim."""

from __future__ import annotations

from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

_GTF = "\n".join([
    'chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG_A"; gene_name "GENEA";',
    'chr1\ttest\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG_A"; gene_name "GENEA";',
    'chr1\ttest\tgene\t20000\t25000\t.\t-\t.\tgene_id "ENSG_B"; gene_name "GENEB";',
    'chr1\ttest\texon\t24800\t25000\t.\t-\t.\tgene_id "ENSG_B"; gene_name "GENEB";',
    "",
])


@pytest.fixture()
def gtf(tmp_path) -> Path:
    p = tmp_path / "annotation.gtf"
    p.write_text(_GTF)
    return p


@pytest.fixture()
def da_csv(tmp_path) -> Path:
    p = tmp_path / "bulk_da_A_vs_B_full.csv"
    pd.DataFrame({
        "peak": ["chr1:1050-1150", "chr1:24900-25000", "chr1:99000-99100"],
        "log2FoldChange": [2.0, -1.5, 0.1],
        "padj": [0.001, 0.002, 0.9],
        "comparison": ["A_vs_B"] * 3,
        "significant": [True, True, False],
    }).to_csv(p, index=False)
    return p


def test_organism_from_genome():
    from aria.scripts.chromatin_peak_ora import _organism_from_genome

    assert _organism_from_genome("hg38") == "Homo sapiens"
    assert _organism_from_genome("GRCh38") == "Homo sapiens"
    assert _organism_from_genome("mm10") == "Mus musculus"
    assert _organism_from_genome("GRCm39") == "Mus musculus"


def test_genes_for_peaks_unique_nearest(gtf):
    from aria.scripts.chromatin_peak_annotation import parse_gtf_gene_models
    from aria.scripts.chromatin_peak_ora import _genes_for_peaks

    models = parse_gtf_gene_models(str(gtf))
    genes = _genes_for_peaks(
        ["chr1:1050-1150", "chr1:1060-1160", "chr1:24900-25000"], models)
    # Two distinct nearest genes, de-duplicated.
    assert set(genes) == {"GENEA", "GENEB"}


def test_peak_ora_honest_skip_without_gtf(tmp_path, monkeypatch):
    from aria.scripts import chromatin_peak_ora as mod

    monkeypatch.setattr(mod, "_resolve_gtf", lambda *a, **k: None)
    res = mod.chromatin_peak_ora({
        "data_type": "bulk_ATAC",
        "comparisons": [{"test": "A", "reference": "B",
                         "full_results_csv": "x.csv"}],
        "output_dir": str(tmp_path / "ora"),
    })
    assert res["status"] == "skipped" and res["ran"] is False
    assert res["reason"] == "no_gtf_for_peak_ora"


def test_peak_ora_end_to_end(tmp_path, gtf, da_csv, monkeypatch):
    """Foreground = sig-peak genes, background = all-tested-peak genes, directional
    up/down split by log2FC; ORA engine + dotplot reused (patched here so the test is
    deterministic and GMT/seaborn-independent)."""
    import aria.scripts.rna_bulk.ora as ora_mod
    import aria.scripts.rna_pathway_viz as viz_mod
    from aria.scripts.chromatin_peak_ora import chromatin_peak_ora

    captured = {}

    def fake_enrich(sig_genes, up_genes, down_genes, organism, output_dir,
                    symbol_map=None, background_genes=None, allow_mock=False):
        captured["sig"] = list(sig_genes)
        captured["up"] = list(up_genes)
        captured["down"] = list(down_genes)
        captured["background"] = list(background_genes or [])
        pathways = {"GO_BP": [
            {"term": "fake pathway", "padj": 1e-4, "overlap": "2/10",
             "odds_ratio": 5.0, "combined_score": 30.0, "genes": "GENEA;GENEB"}]}
        meta = {"method": "local_hypergeometric",
                "gene_set_versions": {"GO_BP": {"release": "2021"}},
                "databases_skipped": []}
        return pathways, [], meta

    def fake_dotplot(terms, db, contrast, output_path, top_n=15):
        Path(output_path).write_bytes(b"PNG")
        return output_path

    monkeypatch.setattr(ora_mod, "_run_pathway_enrichment", fake_enrich)
    monkeypatch.setattr(viz_mod, "make_ora_dotplot", fake_dotplot)

    res = chromatin_peak_ora({
        "data_type": "bulk_ATAC",
        "comparisons": [{"test": "A", "reference": "B",
                         "full_results_csv": str(da_csv)}],
        "gtf": str(gtf),
        "genome": "hg38",
        "output_dir": str(tmp_path / "ora"),
    })

    assert res["status"] == "success" and res["ran"] is True
    assert res["organism"] == "Homo sapiens"
    # Only the 2 significant peaks contribute foreground genes; the non-sig peak
    # is in the background only.
    assert set(captured["sig"]) == {"GENEA", "GENEB"}
    assert captured["up"] == ["GENEA"]      # log2FC > 0 peak near GENEA
    assert captured["down"] == ["GENEB"]    # log2FC < 0 peak near GENEB
    comp = res["comparisons"][0]
    assert comp["pathways_summary"]["GO_BP"] == 1
    assert comp["ora_method"] == "local_hypergeometric"
    assert f"peak_ora_GO_BP_A_vs_B" in res["figures"]
    assert Path(res["figures"]["peak_ora_GO_BP_A_vs_B"]).exists()


def test_peak_ora_narrator_block():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    ora = {
        "status": "success", "ran": True, "data_type": "bulk_ATAC",
        "analysis": "peak_ora", "validation_level": "beta",
        "method": "functional ORA (local hypergeometric)",
        "organism": "Homo sapiens", "gtf": "annotation.gtf",
        "comparisons": [{
            "test": "A", "reference": "B", "ora_method": "local_hypergeometric",
            "pathways_summary": {"GO_BP": 5, "KEGG": 2},
            "databases_skipped": ["Reactome"]}],
    }
    blocks = ChromatinNarrator().collect(
        "chromatin_agent", {"findings": {"peak_ora": ora}})
    block = next(b for b in blocks if b.id == "chromatin.peak_ora")
    assert block.status == "success"
    # Enrichment stays associative, never regulatory.
    assert "regulate" not in block.claim.lower()
    assert any("association" in c.text.lower() for c in block.caveats)
    # Provisioning gaps are disclosed, not hidden.
    assert any("Reactome" in c.text for c in block.caveats)


def test_generate_figures_surfaces_peak_ora_dotplots(tmp_path):
    from aria.agents import _narrative_chromatin

    findings = {
        "differential_accessibility": {"data_type": "bulk_ATAC",
                                       "comparisons": []},
        "peak_ora": {"ran": True, "figures": {
            "peak_ora_GO_BP_A_vs_B": "/tmp/x.png"}},
    }
    out = _narrative_chromatin.generate_figures(
        findings, None, tmp_path / "figures", env_manager=None)
    assert out["figures"].get("peak_ora_GO_BP_A_vs_B") == "/tmp/x.png"
