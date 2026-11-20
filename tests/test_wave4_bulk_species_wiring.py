"""Wave 4B: species, annotation, external ORA, and primary global-FDR wiring."""
from __future__ import annotations

import gzip
import importlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd
import pytest

from aria.scripts.rna_bulk import contrasts, gtf_io, ora
from aria.scripts import chromatin_peak_ora as peak
from aria.scripts import rna_pathway_per_cluster as cluster
from aria.utils import ora as local_ora


MOUSE = ["mouse", "mm", "mm9", "mm10", "mm39", "grcm38", "grcm39",
         "musculus", "mus musculus", "m. musculus"]
HUMAN = ["human", "hs", "hg19", "hg38", "grch37", "grch38", "sapiens",
         "homo sapiens", "h. sapiens"]


@pytest.mark.parametrize("alias,species,api,kegg", [
    *[(a, "Mus musculus", "mouse", "KEGG_2019_Mouse") for a in MOUSE],
    *[(a, "Homo sapiens", "human", "KEGG_2021_Human") for a in HUMAN],
])
def test_species_aliases_are_identical_across_modalities(alias, species, api, kegg):
    for value in (alias, f" {alias.upper()} "):
        assert gtf_io._normalize_species(value) == species
        assert peak._organism_from_genome(value) == species
        for mod in (ora, cluster):
            assert mod._gseapy_organism(value) == api
            assert mod._get_gene_sets(value) == {
                "GO_BP": "GO_Biological_Process_2021", "KEGG": kegg,
                "Reactome": "Reactome_2022"}


@pytest.mark.parametrize("organism,api", [
    ("Rattus norvegicus", "rat"), ("zebrafish", "fish"),
    ("Drosophila melanogaster", "fly"), ("worm", "worm"), ("yeast", "yeast"),
])
def test_other_organisms_keep_support(organism, api):
    for mod in (ora, cluster):
        assert mod._gseapy_organism(organism) == api
        assert "KEGG" not in mod._get_gene_sets(organism)


def _external_rows():
    return pd.DataFrame([
        {"Term": "missing", "Adjusted P-value": None},
        {"Term": "nan", "Adjusted P-value": np.nan},
        {"Term": "invalid", "Adjusted P-value": "bad"},
        {"Term": "negative", "Adjusted P-value": -1},
        {"Term": "inf", "Adjusted P-value": np.inf},
        {"Term": "zero", "Adjusted P-value": 0.0, "P-value": 0.0,
         "Odds Ratio": np.inf, "Combined Score": np.nan, "Genes": np.nan},
        {"Term": "tiny", "Adjusted P-value": 4.795123e-19,
         "P-value": 1.1e-20, "Odds Ratio": np.nan, "Combined Score": np.inf,
         "Genes": 42},
        {"Term": "finite", "Adjusted P-value": "0.001", "Odds Ratio": "4.21",
         "Combined Score": "12.3", "Genes": "Sox2;;nan;Nanog"},
    ])


def test_enrichr_sanitizes_nonfinite_fields_and_keeps_full_precision(monkeypatch):
    calls = []

    def enrichr(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(results=_external_rows())

    monkeypatch.setitem(sys.modules, "gseapy", SimpleNamespace(enrichr=enrichr))
    pathways, warnings = ora._enrichr_enrichment(
        ["Sox2"], ["KEGG"], ora._get_gene_sets("mm39"), "mm39",
        ["SOX2", "NANOG"], False, ["Sox2"])
    assert not warnings
    assert calls[0]["organism"] == "mouse"
    assert calls[0]["gene_sets"] == "KEGG_2019_Mouse"
    entries = pathways["KEGG"]
    assert [e["term"] for e in entries] == ["zero", "tiny", "finite"]
    assert entries[1]["padj"] == entries[1]["adjusted_p"] == 4.795123e-19
    assert entries[0]["odds_ratio"] is entries[0]["combined_score"] is None
    assert entries[0]["genes"] == entries[1]["genes"] == []
    assert entries[2]["genes"] == ["Sox2", "Nanog"]
    json.dumps(pathways, allow_nan=False)


def test_enrichr_sort_handles_missing_pvalue_and_score_ties():
    entries = ora._enrichr_entries(pd.DataFrame([
        {"Term": "missing", "Adjusted P-value": 0.01, "P-value": None,
         "Combined Score": np.nan},
        {"Term": "score", "Adjusted P-value": 0.01, "P-value": np.nan,
         "Combined Score": 2},
        {"Term": "pvalue", "Adjusted P-value": 0.01, "P-value": 0.001,
         "Combined Score": None},
    ]))
    assert [e["term"] for e in entries] == ["pvalue", "score", "missing"]


@pytest.mark.parametrize("direction", ["up", "down"])
def test_local_mouse_ora_runs_single_direction_and_matches_background_case(
        tmp_path, monkeypatch, direction):
    sets = {"stemness": ["SOX2", "NANOG", "POU5F1"]}
    loaded = []

    def load(name):
        loaded.append(name)
        return sets, {"library": name, "release": "test"}

    monkeypatch.setattr(local_ora, "load_local_library", load)
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    foreground = ["Sox2", "Nanog", "Pou5f1"]
    bg = [g.upper() for g in foreground] + [f"BG{i}" for i in range(1000)]
    pw, _, meta = ora._run_pathway_enrichment(
        foreground, foreground if direction == "up" else [],
        foreground if direction == "down" else [], "mm10", str(tmp_path),
        background_genes=bg)
    assert "KEGG_2019_Mouse" in loaded and "KEGG_2021_Human" not in loaded
    assert pw[f"KEGG_{direction}"][0]["genes"] == ["NANOG", "POU5F1", "SOX2"]
    assert f"KEGG_{'down' if direction == 'up' else 'up'}" not in pw
    entry = pw["KEGG"][0]
    assert 0 < entry["padj"] < 1e-6
    assert entry["padj"] == pytest.approx(
        1.0 / math.comb(len(bg), len(foreground)), rel=1e-9)
    if "adjusted_p" in entry:
        assert entry["adjusted_p"] == entry["padj"]
    assert meta["method"] == "local_hypergeometric"


def test_cluster_filters_background_before_cap_and_ranks_finite_significance(
        tmp_path, monkeypatch):
    submitted = []
    real_run = local_ora.run_ora
    sets = {"stemness": ["SOX2", "NANOG"]}
    monkeypatch.setattr(local_ora, "load_local_library",
                        lambda name: (sets, {"library": name}))

    def run(genes, *args, **kwargs):
        submitted.append(list(genes))
        return real_run(genes, *args, **kwargs)

    monkeypatch.setattr(local_ora, "run_ora", run)
    result = cluster.rna_pathway_per_cluster({
        "organism": "mm39", "top_genes_per_cluster": 2,
        "output_dir": str(tmp_path),
        "background_genes_by_cluster": {"c": ["SOX2", "NANOG", "BAD"] +
                                       [f"BG{i}" for i in range(100)]},
        "de_genes_by_cluster": {"c": [
            {"gene": "outside", "padj": 1e-100, "log2fc": 100},
            {"gene": "bad", "padj": None, "log2fc": None},
            {"gene": "Nanog", "padj": 0.002, "log2fc": np.inf},
            {"gene": "Sox2", "padj": 1e-10, "log2fc": 0.1},
            {"gene": "SOX2", "padj": np.nan, "log2fc": np.nan},
        ]}})
    assert result["status"] == "success"
    assert submitted == [["Sox2", "Nanog"]] * 3
    assert result["per_cluster"]["c"]["n_input_genes"] == 2
    assert result["databases"]["KEGG"] == "KEGG_2019_Mouse"
    assert len(pd.read_csv(result["output_csv"])) == 3
    assert not result["warnings"]
    json.dumps(result, allow_nan=False)


def test_cluster_enrichr_uses_shared_safe_rows(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(local_ora, "load_local_library", lambda name: None)
    monkeypatch.setattr(cluster, "_get_gene_sets", lambda organism: {"KEGG": "KEGG_2019_Mouse"})
    monkeypatch.setattr(local_ora, "enrichr_opt_in", lambda: True)
    monkeypatch.setattr("aria.utils.privacy.egress_allowed", lambda: True)

    def enrichr(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(results=_external_rows())

    monkeypatch.setitem(sys.modules, "gseapy", SimpleNamespace(enrichr=enrichr))
    result = cluster.rna_pathway_per_cluster({
        "organism": "mouse", "output_dir": str(tmp_path),
        "background_genes": ["SOX2"],
        "de_genes_by_cluster": {"0": [{"gene": "Sox2", "padj": None, "log2fc": None}]}})
    assert calls[0]["organism"] == "mouse"
    assert calls[0]["gene_list"] == ["Sox2"]
    assert result["per_cluster"]["0"]["n_significant"] == 3
    json.dumps(result, allow_nan=False)


def test_cluster_database_errors_are_warnings_not_csv_terms(tmp_path, monkeypatch):
    monkeypatch.setattr(local_ora, "load_local_library", lambda name: ({"term": ["Sox2"]}, {}))

    def fail(*args, **kwargs):
        raise ValueError("test database failure")

    monkeypatch.setattr(local_ora, "run_ora", fail)
    result = cluster.rna_pathway_per_cluster({
        "organism": "mouse", "output_dir": str(tmp_path),
        "de_genes_by_cluster": {"c": [{"gene": "Sox2"}]}})
    assert len(result["warnings"]) == 3
    assert "test database failure" in result["warnings"][0]
    assert "ERROR" not in Path(result["output_csv"]).read_text()


def _write_gtf(path, gene="ENSMUSG000001.7", symbol="Sox2"):
    text = (f'chr1\ttest\tgene\t1\t1000\t.\t+\t.\tgene_id "{gene}"; '
            f'gene_name "{symbol}"; gene_type "protein_coding";\n'
            f'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "{gene}";\n')
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("alias", MOUSE + HUMAN)
def test_gtf_cache_is_species_specific_and_deterministic(tmp_path, monkeypatch, alias):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    root = tmp_path / ".aria/genomes"
    human = _write_gtf(root / "hg38/annotation.gtf", "ENSG000001.3", "SOX2")
    mouse = _write_gtf(root / "mm39/annotation.gtf.gz")
    counts = tmp_path / "data/counts.csv"
    expected = mouse if alias in MOUSE else human
    assert gtf_io._locate_gtf(counts, organism=alias) == expected
    assert gtf_io._locate_gtf(counts, genome=alias) == expected
    if alias in MOUSE:
        assert gtf_io._load_gene_annotation([str(counts)], [], organism=alias) == {
            "biotype": {"ENSMUSG000001": "protein_coding"},
            "length": {"ENSMUSG000001": 100}}


def test_gtf_prefers_exact_assembly_and_never_uses_other_species(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    root = tmp_path / ".aria/genomes"
    _write_gtf(root / "hg38/annotation.gtf")
    counts = tmp_path / "data/counts.tsv"
    assert gtf_io._locate_gtf(counts, organism="mm39") is None
    _write_gtf(root / "mm10/annotation.gtf")
    exact = _write_gtf(root / "mm39/annotation.gtf.gz")
    assert gtf_io._locate_gtf(counts, genome="mm39") == exact
    assert gtf_io._locate_gtf(counts, genome="GRCm39") == exact
    sibling = _write_gtf(counts.parent / "annotation.gtf")
    assert gtf_io._locate_gtf(counts, genome="mm39") == sibling
    explicit = _write_gtf(tmp_path / "explicit.gtf.gz")
    assert gtf_io._locate_gtf(counts, str(explicit), genome="mm39") == explicit


def test_gtf_csv_symbol_map_persists_and_reload_preserves_samples(tmp_path):
    counts = tmp_path / "counts.csv"
    pd.DataFrame({"s1": [10, 20], "s2": [30, 40]},
                 index=["ENSMUSG000001.7", "ENSMUSG000002.1"]).to_csv(counts)
    gtf = _write_gtf(tmp_path / "tiny.gtf.gz")
    warnings = []
    mapping = gtf_io._load_symbol_map([str(counts)], warnings, str(gtf), organism="mm10")
    assert mapping == {"ENSMUSG000001": "Sox2"}
    table = pd.read_csv(tmp_path / "counts_with_symbols.tsv", sep="\t", index_col=0)
    assert list(table.columns) == ["gene_symbol", "s1", "s2"]
    assert table.loc["ENSMUSG000001.7", "gene_symbol"] == "Sox2"
    assert table["s2"].tolist() == [30, 40]
    # Unmapped Ensembl fallback entries (ENSMUSG000002.1) must not reload into
    # symbol_map and pollute downstream ORA (D1-H7).
    assert gtf_io._load_symbol_map([str(counts)], []) == {
        "ENSMUSG000001": "Sox2"}


@pytest.mark.parametrize("compressed", [False, True])
def test_gtf_edges_version_ids_exon_union_and_gene_length_fallback(tmp_path, compressed):
    text = '\n'.join([
        '# comment', 'malformed',
        'chr1\tt\tgene\t1\t1000\t.\t+\t.\tgene_id   "ENSMUSG000001.7"; gene_name "Sox2"; gene_type "protein_coding";',
        'chr1\tt\texon\t1\t10\t.\t+\t.\tgene_id "ENSMUSG000001.7";',
        'chr1\tt\texon\t5\t20\t.\t+\t.\tgene_id "ENSMUSG000001.7";',
        'chr1\tt\texon\t30\t40\t.\t+\t.\tgene_id "ENSMUSG000001.7";',
        'chr1\tt\texon\t50\t45\t.\t+\t.\tgene_id "ENSMUSG000001.7";',
        'chr1\tt\texon\tbad\t45\t.\t+\t.\tgene_id "ENSMUSG000001.7";',
        'chr1\tt\tgene\t100\t199\t.\t+\t.\tgene_id ENSG000002.3; gene_symbol HLA.DRA; gene_biotype lncRNA;',
        'chr1\tt\texon\t10\t20\t.\t+\t.\tgene_id=ENSG000003.2; gene_name=ExonOnly; gene_type=protein_coding;',
    ])
    path = tmp_path / ("edges.gtf.gz" if compressed else "edges.gtf")
    if compressed:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    assert gtf_io._gtf_to_symbol_map(str(path)) == {
        "ENSMUSG000001": "Sox2", "ENSG000002": "HLA.DRA", "ENSG000003": "ExonOnly"}
    biotype, length = gtf_io._parse_gtf_biotype_and_length(str(path))
    assert biotype == {"ENSMUSG000001": "protein_coding", "ENSG000002": "lncRNA",
                       "ENSG000003": "protein_coding"}
    assert length == {"ENSMUSG000001": 31, "ENSG000002": 100, "ENSG000003": 11}
    assert gtf_io._to_symbols(["ENSMUSG000001.7", "Sox2", "HLA.DRA", None,
                               np.nan, "nan", "", "ENSG_UNMAPPED.2"],
                              {"ENSMUSG000001": "Sox2"}) == ["Sox2", "HLA.DRA"]
    assert gtf_io._to_symbols(["Sox2", "Sox2", None, np.nan, "nan", ""], {}) == ["Sox2"]


def test_top_genes_skip_all_nonfinite_lfc_and_padj():
    frame = pd.DataFrame({"log2FoldChange": [2, -2, np.nan, np.inf, -np.inf, 2, 2, 2],
                          "padj": [1e-19, 0, 0.01, 0.01, 0.01, np.nan, np.inf, None]},
                         index=list("abcdefgh"))
    top = contrasts._format_top_genes({"results": frame, "sig_genes": list(frame.index)})
    assert [g["gene"] for g in top] == ["a", "b"]
    assert [g["direction"] for g in top] == ["up", "down"]
    assert top[0]["padj"] == 1e-19
    json.dumps(top, allow_nan=False)


def test_overlap_honors_empty_full_list_and_consistent_id_namespace():
    result = contrasts._contrast_overlap([
        {"status": "success", "name": "a", "all_sig_genes": [],
         "top_genes": [{"gene": "ENSG1"}]},
        {"status": "success", "name": "b", "all_sig_genes": ["ENSG1"]},
    ])["a ∩ b"]
    assert result["n_in_first"] == result["n_shared"] == 0
    result = contrasts._contrast_overlap([
        {"status": "success", "name": "a", "all_sig_genes": [], "all_sig_gene_ids": ["ENSG1"]},
        {"status": "success", "name": "b", "all_sig_genes": ["Sox2"], "all_sig_gene_ids": ["ENSG1"]},
    ])["a ∩ b"]
    assert result["n_shared"] == 1


def _bulk_case(tmp_path, monkeypatch, *, family, pathways, n_contrasts=2, lose_all=False,
               failed_second=False):
    mod = importlib.import_module("aria.scripts.rna_bulk_de")
    viz = importlib.import_module("aria.scripts.rna_pathway_viz")
    samples = [f"{group}{i}" for group in ("A", "B", "C") for i in range(3)]
    genes = ["ENSMUSG000001", "ENSMUSG000002", "ENSMUSG000003", "ENSMUSG000004"]
    counts_path, meta_path = tmp_path / "counts.tsv", tmp_path / "metadata.tsv"
    pd.concat([pd.DataFrame(30, index=genes, columns=samples),
               pd.DataFrame(0, index=[f"low{i}" for i in range(8)], columns=samples)]
              ).to_csv(counts_path, sep="\t")
    pd.DataFrame({"donor": range(9), "treatment": [s[0] for s in samples]},
                 index=samples).to_csv(meta_path, sep="\t")
    observed = {"fits": [], "plots": {}, "ora": [], "gsea": {}, "lobo": [], "annotation": []}
    frames = {}

    def fit(counts, metadata, design_factor, numerator, denominator, *args, **kwargs):
        observed["fits"].append(numerator)
        if failed_second and numerator == "C":
            return {"status": "error", "details": "test fit failure"}, []
        pvals = ([0.01, 0.02, 0.8, 0.9] if lose_all else [0.0001, 0.02, 0.8, 0.9])
        if numerator == "C":
            pvals = [0.8, 0.9, 0.9, 0.9]
        from aria.utils.stats import bh_correct
        frame = pd.DataFrame({"pvalue": pvals, "padj": bh_correct(pvals),
                              "log2FoldChange": [0.2, -2, 0.1, -0.1],
                              "baseMean": [30.0] * 4}, index=genes)
        frames[numerator] = frame.copy()
        sig = frame[frame["padj"] < 0.05]
        return {"status": "success", "results": frame,
                "sig_genes": list(sig.index), "n_sig": len(sig),
                "up_genes": list(sig.index[sig["log2FoldChange"] > 0]),
                "down_genes": list(sig.index[sig["log2FoldChange"] < 0]),
                "n_up": int((sig["log2FoldChange"] > 0).sum()),
                "n_down": int((sig["log2FoldChange"] < 0).sum())}, []

    def qc(*args, design_factor=None, **kwargs):
        observed["qc_factor"] = design_factor
        return {"outliers": [], "n_samples": 9}

    def annotation(*args, **kwargs):
        observed["annotation"].append(kwargs)
        return {}

    def symbols(*args, **kwargs):
        observed["symbol_kwargs"] = kwargs
        return dict(zip(genes, ["Sox2", "Nanog", "Other1", "Other2"]))

    def plots(**kwargs):
        if family == "global":
            assert len(observed["fits"]) == n_contrasts
        result = kwargs["de_result"]
        observed["plots"][kwargs["title_suffix"]] = (
            result["results"].copy(), list(result["sig_genes"]),
            result["n_sig"], result["n_up"], result["n_down"])
        return {}

    def enrich(**kwargs):
        observed["ora"].append(kwargs)
        return {"GO_BP": [{"term": "stemness", "padj": 1e-9, "genes": ["Sox2"],
                            "odds_ratio": 2, "overlap": "1/2"}]}, [], {"method": "test"}

    def gsea(**kwargs):
        assert pathways, "GSEA ran despite run_pathways=False"
        observed["gsea"][kwargs["contrast_name"]] = kwargs["de_results_df"].copy()
        return {}

    def lobo(*args, **kwargs):
        observed["lobo"].append(kwargs)
        return {"status": "success"}

    monkeypatch.setattr(mod, "_run_deseq2", fit)
    monkeypatch.setattr(mod, "_sample_qc", qc)
    monkeypatch.setattr(mod, "_load_gene_annotation", annotation)
    monkeypatch.setattr(mod, "_load_symbol_map", symbols)
    monkeypatch.setattr(mod, "_generate_plots", plots)
    monkeypatch.setattr(mod, "_run_pathway_enrichment", enrich)
    monkeypatch.setattr(mod, "compute_counts_lobo_jackknife", lobo)
    monkeypatch.setattr(viz, "make_gsea_running_sums", gsea)
    monkeypatch.setattr(viz, "make_ora_dotplot", lambda **kwargs: None)
    requested = [{"numerator": "A", "denominator": "B", "name": "A vs B"},
                 {"numerator": "C", "denominator": "B", "name": "C vs B"}][:n_contrasts]
    result = mod.bulk_rna_de({
        "files": [str(counts_path)], "metadata_file": str(meta_path),
        "design_factor": "treatment", "organism": "mouse", "genome": "mm39",
        "gtf": str(tmp_path / "explicit.gtf"), "contrasts": requested,
        "fdr_family": family, "run_pathways": pathways, "output_dir": str(tmp_path / "out"),
    })
    assert result["status"] == "success", result
    return result, observed, frames


@pytest.mark.parametrize("family,n_contrasts,lose_all", [
    ("global", 2, False), ("global", 2, True), ("global", 1, False),
    ("per_contrast", 2, False), ("per_contrast", 1, False),
])
def test_primary_fdr_reaches_every_artifact(tmp_path, monkeypatch, family, n_contrasts, lose_all):
    result, seen, local_frames = _bulk_case(
        tmp_path, monkeypatch, family=family, pathways=True,
        n_contrasts=n_contrasts, lose_all=lose_all)
    assert result["status"] == "success"
    assert seen["qc_factor"] == "treatment"
    assert seen["annotation"][0] == seen["symbol_kwargs"] == {
        "gtf_hint": str(tmp_path / "explicit.gtf"), "genome": "mm39", "organism": "mouse"}
    # Independent pooled BH oracle; no wrapper implementation or adjusted field
    # is used to construct expected primary significance.
    from aria.utils.stats import bh_correct
    flattened = np.concatenate([f["pvalue"].to_numpy() for f in local_frames.values()])
    pooled = bh_correct(flattened)
    expected_ora, total = [], 0
    for i, comp in enumerate(result["contrasts"]):
        source = local_frames[comp["numerator"]]
        expected_padj = pooled[i * 4:(i + 1) * 4] if family == "global" else source["padj"]
        expected_sig = list(source.index[np.asarray(expected_padj) < 0.05])
        expected_up = [g for g in expected_sig if source.loc[g, "log2FoldChange"] > 0]
        expected_down = [g for g in expected_sig if source.loc[g, "log2FoldChange"] < 0]
        total += len(expected_sig)
        assert comp["n_significant"] == len(expected_sig)
        assert comp["n_upregulated"] == len(expected_up)
        assert comp["n_downregulated"] == len(expected_down)
        assert comp["all_sig_gene_ids"] == expected_sig
        assert comp["up_gene_ids"] == expected_up
        assert comp["down_gene_ids"] == expected_down
        assert [g["gene"] for g in comp["top_genes"]] == expected_sig
        frame, sig, nsig, nup, ndown = seen["plots"][comp["name"]]
        np.testing.assert_allclose(frame["padj"], expected_padj)
        assert (sig, nsig, nup, ndown) == (
            expected_sig, len(expected_sig), len(expected_up), len(expected_down))
        np.testing.assert_allclose(seen["gsea"][comp["name"]]["padj"], expected_padj)
        table = pd.read_csv(comp["full_results_csv"], sep="\t")
        np.testing.assert_allclose(table["padj"], expected_padj)
        de_table = comp["plots"]["tables"].get("de_genes")
        if expected_sig:
            significant = pd.read_csv(de_table, sep="\t")
            assert significant["gene_id"].tolist() == expected_sig
            np.testing.assert_allclose(significant["padj"], np.asarray(expected_padj)[:len(expected_sig)])
            expected_ora.append(expected_sig)
        else:
            assert de_table is None
            assert not comp["pathways"] and not comp["top_genes"]
        assert seen["lobo"][i]["primary_significant_features"] == expected_sig
    assert [call["sig_genes"] for call in seen["ora"]] == expected_ora
    assert result["n_significant"] == total
    assert result["fdr_family"]["n_tests_family"] == 4 * n_contrasts
    assert not any("failed" in warning.lower() for warning in result["warnings"])


@pytest.mark.parametrize("family", ["global", "per_contrast"])
def test_run_pathways_false_disables_both_ora_and_gsea(tmp_path, monkeypatch, family):
    result, seen, _ = _bulk_case(tmp_path, monkeypatch, family=family, pathways=False)
    assert result["status"] == "success"
    assert not seen["ora"] and not seen["gsea"]
    assert all(c["full_results_csv"] for c in result["contrasts"])
    assert not any("visualization failed" in w for w in result["warnings"])


def test_global_family_excludes_failed_contrast_without_losing_success(tmp_path, monkeypatch):
    result, seen, frames = _bulk_case(
        tmp_path, monkeypatch, family="global", pathways=False, failed_second=True)
    assert [c["status"] for c in result["contrasts"]] == ["success", "error"]
    assert result["fdr_family"]["n_tests_family"] == 4
    np.testing.assert_allclose(seen["plots"]["A vs B"][0]["padj"], frames["A"]["padj"])


@pytest.mark.parametrize("data_type", ["bulk_ATAC", "scATAC"])
@pytest.mark.parametrize("svg_mode", ["paired", "returned", "absent"])
def test_peak_ora_canonical_mouse_and_svg_parity(tmp_path, monkeypatch, data_type, svg_mode):
    viz = importlib.import_module("aria.scripts.rna_pathway_viz")
    gtf = _write_gtf(tmp_path / "annotation.gtf")
    csv = tmp_path / "da.csv"
    pd.DataFrame({"peak": ["chr1:10-20"], "significant": [True],
                  "log2FoldChange": [2]}).to_csv(csv, index=False)
    captured = []

    def enrich(*args, **kwargs):
        captured.append(args[3])
        return {"KEGG": [{"term": "stemness", "padj": 1e-10}]}, [], {"method": "test"}

    def dotplot(terms, db, name, output):
        path = Path(output)
        path.write_bytes(b"PNG")
        if svg_mode != "absent":
            path.with_suffix(".svg").write_text("<svg/>")
        return str(path.with_suffix(".svg")) if svg_mode == "returned" else str(path)

    monkeypatch.setattr(ora, "_run_pathway_enrichment", enrich)
    monkeypatch.setattr(viz, "make_ora_dotplot", dotplot)
    result = peak.chromatin_peak_ora({
        "data_type": data_type, "gtf": str(gtf), "organism": "mm39",
        "comparisons": [{"test": "A", "reference": "B", "full_results_csv": str(csv)}],
        "output_dir": str(tmp_path / "out")})
    assert result["status"] == "success"
    assert captured == ["Mus musculus"] and result["organism"] == "Mus musculus"
    comp = result["comparisons"][0]
    prefix = "scatac_" if data_type == "scATAC" else ""
    key = f"{prefix}peak_ora_KEGG_A_vs_B"
    assert key in result["figures"] and result["figures"] == comp["figures"]
    assert result["figures_svg"] == comp["figures_svg"]
    if svg_mode == "absent":
        assert result["figures_svg"] == {}
    else:
        assert key in result["figures_svg"]
        assert Path(result["figures_svg"][key]).is_file()


def test_safe_enrichr_rows_render_wave2_dotplot_with_svg(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    pytest.importorskip("seaborn")
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mpl"))
    from aria.scripts.rna_pathway_viz import make_ora_dotplot

    rows = ora._enrichr_entries(_external_rows())
    path = tmp_path / "external_ora.png"
    assert make_ora_dotplot(rows, "KEGG", "mouse", str(path)) == str(path)
    assert path.stat().st_size > 0
    assert path.with_suffix(".svg").is_file()
