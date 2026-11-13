"""ADR-080: scientific direction, co-tested backgrounds, and honest non-runs."""
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact, pearsonr, spearmanr

from aria.scripts.cross_modal_concordance import cross_modal_concordance
from aria.utils.cross_modal_concordance import (
    CAUSAL_PATTERNS,
    bh_adjust,
    find_causal_language,
    read_cpm,
    tf_subunits,
)


@pytest.fixture(autouse=True)
def local_references_only(monkeypatch):
    # Tests never use machine-local reference installations.
    monkeypatch.setattr("aria.scripts.cross_modal_concordance._resolve_gtf",
                        lambda params, _: Path(params["gtf_path"]) if params.get("gtf_path") and Path(params["gtf_path"]).is_file() else None)


@pytest.fixture
def bulk_inputs(tmp_path):
    genes = ["FOS", "JUN", "GATA1", "GATA2", "STAT1", "STAT3", "SOX2", "SOX4", "KLF4", "KLF5", "MYC", "TP53"]
    positions = [10000 + i * 300000 for i in range(len(genes))]
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("".join(f'chr1\ttest\tgene\t{pos}\t{pos + 3000}\t.\t+\t.\tgene_id "ENSG{i:06}"; gene_name "{g}";\n'
                           for i, (g, pos) in enumerate(zip(genes, positions))))
    peaks = [f"chr1:{pos - 100}-{pos + 100}" for pos in positions]
    tables = tmp_path / "rna" / "treated_vs_control" / "tables"
    tables.mkdir(parents=True)
    rna_lfc = [4, 3, 2, -4, -3, -2, 1, -1, 2, .2, 3, .1]
    atac_lfc = [4, 3, 2, -4, -3, -2, -1, 1, .2, 2, .1, .1]
    rna_p = [.001] * 9 + [.8, .001, .8]
    atac_p = [.002] * 8 + [.8, .002, .8, .8]
    rna_table = tables / "de_genes.tsv"
    pd.DataFrame({"gene_id": [f"ENSG{i:06}.1" for i in range(12)], "log2FoldChange": rna_lfc, "padj": rna_p}).to_csv(rna_table, sep="\t", index=False)
    atac_table = tmp_path / "bulk_atac_da_treated_vs_control_full.csv"
    pd.DataFrame({"peak": peaks, "log2FoldChange": atac_lfc, "padj": atac_p,
                  "significant": np.array(atac_p) < .05}).to_csv(atac_table, index=False)
    annotation = tmp_path / "significant_peak_annotation.csv"
    pd.DataFrame({"peak": peaks[:8], "nearest_gene": genes[:8], "distance_to_tss": [0] * 8, "feature": ["Promoter"] * 8}).to_csv(annotation, index=False)
    pathways = {"GO": [{"term": "Cell cycle (GO:0007049)", "padj": .001},
                        {"term": "TP53 regulates transcription R-HSA-12345", "padj": .002},
                        {"term": "RNA processing", "padj": .003}]}
    motifs = {"ran": True, "per_group": {"treated_vs_control_up": {"top_motifs": [
        {"motif_id": "MA0099.3", "name": "FOS::JUN(var.2)", "family": "AP-1", "padj": .001, "log2_enrichment": 2},
        {"motif_id": "MA0035.4", "name": "GATA1", "family": "GATA", "padj": .002, "log2_enrichment": 1.5},
        {"motif_id": "MA0036.3", "name": "GATA2", "family": "GATA", "padj": .003, "log2_enrichment": 1.2},
        {"motif_id": "MA0493.1", "name": "KLF5", "family": "KLF", "padj": .004, "log2_enrichment": 1},
        {"motif_id": "MA0100.1", "name": "ABSENT", "family": "other", "padj": .005, "log2_enrichment": .8},
    ]}}}
    footprint = {"ran": True, "group_a": "treated", "group_b": "control", "differential_summary": {
        "parsed": True, "inference": {"status": "success"},
        "all_results": [{"tf": "FOS", "change": 2, "padj": .002}, {"tf": "STAT1", "change": -1, "padj": .003}]}}
    rc = {"name": "treated_vs_control", "numerator": "treated", "denominator": "control", "tables_dir": str(tables), "pathways": pathways}
    ac = {"test": "treated", "reference": "control", "full_results_csv": str(atac_table)}
    atac = {"differential_accessibility": {"comparisons": [ac]},
            "peak_annotation": {"comparisons": [{"test": "treated", "reference": "control", "annotation_csv": str(annotation)}]},
            "motifs": motifs, "footprinting": footprint,
            "peak_ora": {"comparisons": [{"test": "treated", "reference": "control", "pathways": {"Reactome": [
                {"term": "CELL CYCLE", "padj": .004}, {"term": "TP53 regulates transcription (R-HSA-12345)", "padj": .005},
                {"term": "Chromatin organization", "padj": .006}]}}]}}
    return {"rna_results": {"status": "done", "findings": {"contrasts": [rc]}},
            "chromatin_results": {"status": "done", "findings": {"bulk_ATAC": {"findings": atac}}},
            "gtf_path": str(gtf), "output_dir": str(tmp_path / "out")}


def _rna(params):
    return params["rna_results"]["findings"]["contrasts"][0]


def _atac(params):
    return params["chromatin_results"]["findings"]["bulk_ATAC"]["findings"]


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for k, v in value.items():
            yield k
            yield from _strings(v)
    elif isinstance(value, list):
        for v in value:
            yield from _strings(v)


def test_contract_production_and_optional_inputs(tmp_path):
    pytest.importorskip("pydantic")
    from aria.utils.script_contracts import contract_for_script
    try:
        from aria.agents.narrative.validators import CAUSAL_PATTERNS as NARR_CAUSAL
        assert tuple(CAUSAL_PATTERNS) == tuple(NARR_CAUSAL)
    except ImportError:
        pass
    contract = contract_for_script("aria/scripts/cross_modal_concordance.py")
    assert contract.validation_level == "production"
    assert contract.validate_params({}) == []
    assert {f.name for f in contract.inputs} == {
        "rna_results", "chromatin_results", "rna_cpm_path", "atac_cpm_path", "sample_id_map",
        "gtf_path", "organism", "genome", "output_dir", "padj_threshold", "cis_distance_bp",
        "promoter_bp", "min_matched_samples", "n_permutations", "seed", "resume"}
    assert contract.validate_params({"sample_id_map": "invalid"})
    result = cross_modal_concordance({"output_dir": str(tmp_path)})
    assert result["ran"] is False
    assert contract.validate_result(result) == []
    result["ran"] = True
    assert contract.validate_result(result) == []


def test_full_bulk_run_background_statistics_and_artifacts(bulk_inputs):
    result = cross_modal_concordance(bulk_inputs)
    assert result["status"] == "success" and result["ran"]
    assert result["rna_modality"] == "bulk_RNA"
    assert result["atac_modality"] == "bulk_ATAC"
    assert result["n_matched_contrasts"] == 1
    summary = result["matched_contrasts"][0]["de_da_summary"]
    assert summary["annotation_scope"] == "all_tested_peaks"
    assert summary["n_pairs_tested"] == summary["n_genes_tested"] == 12
    for key, count in (("n_concordant_up", 3), ("n_concordant_down", 3), ("n_discordant_total", 2),
                       ("n_rna_de_atac_unchanged", 2), ("n_atac_da_rna_unchanged", 1), ("n_both_unchanged", 1)):
        assert summary[key] == count
    expected_odds, expected_p = fisher_exact([[3, 1], [1, 3]], alternative="two-sided")
    assert summary["fisher_odds_ratio"] == pytest.approx(expected_odds)
    assert summary["fisher_pvalue"] == pytest.approx(expected_p)
    assert summary["concordance_fraction"] == .75
    assert summary["spearman_r"] > 0
    frame = pd.read_csv(result["tables"]["de_da_concordance_csv"])
    assert summary["spearman_r"] == pytest.approx(spearmanr(frame.rna_log2fc, frame.atac_log2fc).statistic)
    assert summary["pearson_r"] == pytest.approx(pearsonr(frame.rna_log2fc, frame.atac_log2fc).statistic)
    assert summary["promoter_spearman_r"] == summary["spearman_r"]
    assert summary["enrichment_pvalue"] == pytest.approx(fisher_exact([[8, 2], [1, 1]]).pvalue)
    assert frame.is_primary_peak.all()
    assert result["pathway_concordance"]["n_shared_terms"] == 2
    assert result["pathway_concordance"]["jaccard_index"] == .5
    tf = pd.read_csv(result["tables"]["tf_triangulation_csv"])
    assert set(tf.tf) >= {"FOS", "JUN", "GATA1", "GATA2", "KLF5", "STAT1", "ABSENT"}
    assert set(tf[tf.tf == "FOS"].source) == {"motif", "footprint"}
    assert tf.loc[tf.tf == "KLF5", "classification"].iloc[0] == "chromatin_shift_rna_unchanged"
    assert tf.loc[tf.tf == "ABSENT", "classification"].iloc[0] == "not_detected_in_rna"
    for key in ("de_da_concordance_csv", "tf_triangulation_csv", "shared_pathways_csv"):
        assert Path(result["tables"][key]).stat().st_size > 100
    for key in ("figures", "figures_svg"):
        assert Path(result[key]["cross_modal_starburst"]).stat().st_size > 5000
    try:
        from aria.utils.script_contracts import contract_for_script
        assert contract_for_script("aria/scripts/cross_modal_concordance.py").validate_result(result) == []
    except ImportError:
        pass
    json.dumps(result, allow_nan=False)


def test_reverse_sign_alignment_including_motifs_and_footprints(bulk_inputs):
    comp = _atac(bulk_inputs)["differential_accessibility"]["comparisons"][0]
    comp.update(test="control", reference="treated")
    frame = pd.read_csv(comp["full_results_csv"])
    frame.log2FoldChange *= -1
    frame.to_csv(comp["full_results_csv"], index=False)
    atac = _atac(bulk_inputs)
    # Use generic up/down groups in the reversed ATAC direction.
    atac["motifs"]["per_group"] = {"down": atac["motifs"]["per_group"].pop("treated_vs_control_up")}
    atac["footprinting"].update(group_a="control", group_b="treated")
    for row in atac["footprinting"]["differential_summary"]["all_results"]:
        row["change"] *= -1
    result = cross_modal_concordance(bulk_inputs)
    matched = result["matched_contrasts"][0]
    assert matched["atac_sign_flipped"] is True
    assert matched["de_da_summary"]["n_concordant_up"] == 3
    assert matched["de_da_summary"]["n_concordant_down"] == 3
    tf = pd.read_csv(result["tables"]["tf_triangulation_csv"])
    assert (tf.loc[tf.tf == "FOS", "classification"] == "concordant").all()
    assert (tf.loc[tf.tf == "FOS", "chromatin_direction"] == "up").all()


def test_incompatible_contrast_denominators(tmp_path):
    result = cross_modal_concordance({"rna_results": {"contrasts": [{"name": "A_vs_B"}]},
        "chromatin_results": {"differential_accessibility": {"comparisons": [{"test": "A", "reference": "C"}]}},
        "output_dir": str(tmp_path)})
    assert result["status"] == "success" and result["ran"] is False
    assert result["reason"] == "incompatible_contrast_denominators"
    assert result["n_matched_contrasts"] == 0
    assert result["warnings"] and not result["figures"]


def test_dimer_and_family_ambiguity(bulk_inputs):
    assert tf_subunits("FOS::JUN(var.2)") == ["FOS", "JUN"]
    path = Path(_rna(bulk_inputs)["tables_dir"]) / "de_genes.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame.loc[1, "log2FoldChange"] = -3
    frame.to_csv(path, sep="\t", index=False)
    result = cross_modal_concordance(bulk_inputs)
    tf = pd.read_csv(result["tables"]["tf_triangulation_csv"])
    assert tf.loc[tf.tf.isin(["FOS", "JUN", "GATA1", "GATA2"]), "paralog_motif_ambiguity"].all()
    assert tf.loc[(tf.tf == "JUN") & (tf.source == "motif"), "classification"].iloc[0] == "discordant"
    assert result["tf_triangulation"]["n_paralog_ambiguous"] >= 4


def _sample_inputs(tmp_path, n=8, disjoint=False, noise=False, n_genes=1):
    rng = np.random.default_rng(15)
    donors = [f"D{i}" for i in range(n)]
    gtf = tmp_path / "sample_genes.gtf"
    genes = [f"G{i}" for i in range(n_genes)]
    positions = [10000 + i * 300000 for i in range(n_genes)]
    gtf.write_text("".join(f'1\ttest\tgene\t{pos}\t{pos + 1000}\t.\t+\t.\tgene_id "id{i}"; gene_name "{g}";\n'
                           for i, (g, pos) in enumerate(zip(genes, positions))))
    x = rng.normal(5, 1.5, size=(n_genes, n))
    y = np.array([rng.permutation(row) for row in x]) if noise else x + rng.normal(0, .005, size=x.shape)
    r = pd.DataFrame(np.exp2(x), index=genes, columns=donors)
    a = pd.DataFrame(np.exp2(y), index=[f"chr1:{p - 10}-{p + 10}" for p in positions], columns=[f"ATAC_{i}" for i in range(n)] if disjoint else donors)
    rp, ap = tmp_path / "rna_cpm.tsv", tmp_path / "atac_cpm.tsv"
    r.to_csv(rp, sep="\t", index_label="gene")
    a.to_csv(ap, sep="\t", index_label="peak")
    return {"rna_cpm_path": str(rp), "atac_cpm_path": str(ap), "gtf_path": str(gtf),
            "n_permutations": 200, "seed": 11, "output_dir": str(tmp_path / "sample_out")}


def test_sample_correlations_empirical_recovery_and_heldout(tmp_path):
    result = cross_modal_concordance(_sample_inputs(tmp_path))
    layer = result["sample_matched_links"]
    assert layer["ran"] and result["ran"]
    assert layer["n_shared_samples"] == 8
    assert layer["n_significant_links"] == layer["n_positive_links"] == 1
    assert layer["top_links"][0]["padj"] < .05
    assert layer["top_links"][0]["pearson_r"] > .99
    assert layer["top_links"][0]["spearman_r"] > .99
    assert layer["permutation_null"]["n_null_statistics"] == 200
    assert layer["permutation_null"]["held_out_n_significant"] == 0
    assert Path(result["tables"]["sample_matched_links_csv"]).stat().st_size > 100


def test_sample_noise_calibration(tmp_path):
    result = cross_modal_concordance(_sample_inputs(tmp_path, noise=True, n_genes=80))
    layer = result["sample_matched_links"]
    assert layer["ran"] and layer["n_pairs_tested"] == 80
    assert layer["n_significant_links"] <= 2
    assert layer["permutation_null"]["held_out_n_significant"] <= 2


@pytest.mark.parametrize("n,disjoint,reason,count", [(4, False, "insufficient_matched_biological_units", 4),
                                                    (8, True, "unmatched_sample_namespaces", 0)])
def test_sample_honest_nonruns(tmp_path, n, disjoint, reason, count):
    layer = cross_modal_concordance(_sample_inputs(tmp_path, n, disjoint))["sample_matched_links"]
    assert layer["ran"] is False and layer["reason"] == reason
    assert layer["n_shared_samples"] == count and layer["min_required"] == 6


def test_explicit_sample_map_and_missing_gtf(tmp_path):
    params = _sample_inputs(tmp_path, disjoint=True)
    params["sample_id_map"] = {f"D{i}": f"ATAC_{i}" for i in range(8)}
    result = cross_modal_concordance(params)
    assert result["sample_matched_links"]["ran"]
    params.pop("gtf_path")
    result = cross_modal_concordance(params)
    assert result["sample_matched_links"]["reason"] == "gtf_gene_models_unavailable"
    params["sample_id_map"]["D1"] = "ATAC_0"
    assert cross_modal_concordance(params)["sample_matched_links"]["reason"] == "non_unique_sample_id_map"


@pytest.mark.parametrize("stratified", [False, True])
def test_scrna_scatac_cohort_matching_and_significant_csv_scope(bulk_inputs, stratified):
    rc = _rna(bulk_inputs)
    genes = [{"gene": "FOS", "log2fc": 4, "padj": .001}, {"gene": "JUN", "log2fc": 3, "padj": .001}]
    pb_path = Path(bulk_inputs["output_dir"]).parent / "pseudobulk_de.csv"
    pd.DataFrame([{**g, "group": group, "comparison": "treated_vs_control"} for group in ("T", "B") for g in genes]).to_csv(pb_path, index=False)
    pb = {"per_group": {g: {"per_comparison": {"treated_vs_control": {"status": "success", "all_sig": genes}}} for g in ("T", "B")}, "output_csv": str(pb_path)}
    atac = copy.deepcopy(_atac(bulk_inputs))
    comps = atac["differential_accessibility"]["comparisons"]
    if stratified:
        comps = [{**comps[0], "group": g} for g in ("T", "B")]
    atac["differential_accessibility"] = {"pseudobulk": {"ran": True, "comparisons": comps}}
    bulk_inputs["rna_results"] = {"findings": {"pseudobulk_de": pb, "pathways": {"per_cluster": {"T": {"results": rc["pathways"]}}}}}
    bulk_inputs["chromatin_results"] = {"findings": {"scATAC": {"findings": atac}}}
    result = cross_modal_concordance(bulk_inputs)
    assert result["ran"] and result["n_matched_contrasts"] == 2
    assert result["rna_modality"] == "scRNA" and result["atac_modality"] == "scATAC"
    for matched in result["matched_contrasts"]:
        assert matched["matching_granularity"] == ("cell_group_matched" if stratified else "condition_pseudobulk_vs_cell_groups")
        assert matched["de_da_summary"]["rna_scope"] == "reported_genes_only"
    table = pd.read_csv(result["tables"]["de_da_concordance_csv"])
    assert set(table.best_cell_group) == {"T", "B"}
    assert table.rna_log2fc.tolist() == [4, 3, 4, 3]


def test_manifest_resume_table_threshold_and_output_invalidation(bulk_inputs):
    first = cross_modal_concordance(bulk_inputs)
    assert first["resumed_from_manifest"] is False
    bulk_inputs["resume"] = True
    second = cross_modal_concordance(bulk_inputs)
    assert second["resumed_from_manifest"] is True
    path = Path(_rna(bulk_inputs)["tables_dir"]) / "de_genes.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame.loc[0, "log2FoldChange"] *= -1
    frame.to_csv(path, sep="\t", index=False)
    third = cross_modal_concordance(bulk_inputs)
    assert third["resumed_from_manifest"] is False
    assert third["de_da_concordance"]["n_concordant_up"] == 2
    assert cross_modal_concordance(bulk_inputs)["resumed_from_manifest"] is True
    bulk_inputs["padj_threshold"] = .0001
    assert cross_modal_concordance(bulk_inputs)["resumed_from_manifest"] is False
    cached = cross_modal_concordance(bulk_inputs)
    Path(cached["figures"]["cross_modal_starburst"]).unlink()
    assert cross_modal_concordance(bulk_inputs)["resumed_from_manifest"] is False
    manifest = json.loads(Path(first["stage_manifest_path"]).read_text())
    paths = {i["path"] for i in manifest["inputs"]}
    assert str(path) in paths and bulk_inputs["gtf_path"] in paths


def test_causal_guard_json_svg_and_raw_database_csv(bulk_inputs):
    result = cross_modal_concordance(bulk_inputs)
    assert all(find_causal_language(s) is None for s in _strings(result))
    svg = Path(result["figures_svg"]["cross_modal_starburst"]).read_text()
    assert find_causal_language(svg) is None
    assert "<text" in svg and "TP53 association" in svg and "transcription" in svg
    pathways = pd.read_csv(result["tables"]["shared_pathways_csv"])
    assert any("regulates" in x for x in pathways.term)
    assert all(find_causal_language(x) is None for x in pathways.safe_label)
    assert all(not Path(s).is_absolute() for s in _strings(result["warnings"]))


def test_annotation_only_discloses_selected_background(bulk_inputs):
    # Symbol mapping is supplied upstream when the GTF is absent.
    path = Path(_rna(bulk_inputs)["tables_dir"]) / "de_genes.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame["gene_symbol"] = ["FOS", "JUN", "GATA1", "GATA2", "STAT1", "STAT3", "SOX2", "SOX4", "KLF4", "KLF5", "MYC", "TP53"]
    frame.to_csv(path, sep="\t", index=False)
    bulk_inputs.pop("gtf_path")
    result = cross_modal_concordance(bulk_inputs)
    summary = result["matched_contrasts"][0]["de_da_summary"]
    assert summary["annotation_scope"] == "significant_peaks_only"
    assert summary["n_pairs_tested"] == 8
    assert summary["n_rna_de_atac_unchanged"] == 0
    assert any("significant_peaks_only" in w for w in result["warnings"])


def test_in_memory_effects_no_fabricated_magnitude(tmp_path):
    result = cross_modal_concordance({"rna_results": {"contrasts": [{"name": "A_vs_B", "up_gene_ids": ["X"], "down_gene_ids": ["Y"]}]},
        "chromatin_results": {"comparisons": [{"test": "A", "reference": "B", "top_peaks": [
            {"peak": "chr1:90-110", "nearest_gene": "X", "distance_to_tss": 0, "log2fc": 2, "padj": .001},
            {"peak": "chr1:290-310", "nearest_gene": "Y", "distance_to_tss": 0, "log2fc": -2, "padj": .001}]}]},
        "output_dir": str(tmp_path)})
    summary = result["de_da_concordance"]
    assert result["ran"] and summary["n_concordant_total"] == 2
    assert summary["spearman_r"] is None
    assert summary["fisher_odds_ratio"] is None  # Infinite odds, finite p-value.
    assert summary["fisher_pvalue"] is not None
    table = pd.read_csv(result["tables"]["de_da_concordance_csv"])
    assert table.rna_log2fc.isna().all()
    assert set(table.gene) == {"X", "Y"}


def test_descriptive_footprints_and_incompatible_motif_groups_are_skipped(bulk_inputs):
    atac = _atac(bulk_inputs)
    atac["motifs"]["per_group"] = {"other_vs_reference_up": atac["motifs"]["per_group"].pop("treated_vs_control_up")}
    atac["footprinting"]["differential_summary"]["inference"]["status"] = "descriptive_only"
    result = cross_modal_concordance(bulk_inputs)
    assert result["tf_triangulation"]["ran"] is False
    assert result["tables"]["tf_triangulation_csv"] is None
    assert any("descriptive_ranking" in w for w in result["warnings"])


def test_primary_peak_ties_and_unique_gene_enrichment(bulk_inputs):
    comp = _atac(bulk_inputs)["differential_accessibility"]["comparisons"][0]
    path = comp["full_results_csv"]
    frame = pd.read_csv(path)
    # Three FOS peaks share padj; the closest TSS wins, then the larger effect.
    frame.loc[0, "peak"] = "chr1:10400-10600"
    extra = pd.DataFrame({"peak": ["chr1:10000-10200", "chr1:9800-10000"],
                          "log2FoldChange": [2, -3], "padj": [.002, .002], "significant": [True, True]})
    pd.concat([frame, extra]).to_csv(path, index=False)
    result = cross_modal_concordance(bulk_inputs)
    summary = result["matched_contrasts"][0]["de_da_summary"]
    assert summary["n_pairs_tested"] == 14 and summary["n_genes_tested"] == 12
    table = pd.read_csv(result["tables"]["de_da_concordance_csv"])
    fos = table[(table.gene == "FOS") & table.is_primary_peak]
    assert fos.peak.tolist() == ["chr1:9800-10000"]
    assert summary["gene_counts"]["n_discordant_total"] == 3
    assert summary["gene_enrichment_contingency"] == [[8, 2], [1, 1]]


def test_missing_sc_pathways_are_not_a_completed_layer(bulk_inputs):
    atac = _atac(bulk_inputs)
    atac.pop("motifs")
    atac.pop("footprinting")
    bulk_inputs["rna_results"] = {"pseudobulk_de": {"per_group": {"T": {
        "per_comparison": {"treated_vs_control": {"all_sig": []}}}}}}
    result = cross_modal_concordance(bulk_inputs)
    assert not result["ran"]
    assert not result["pathway_concordance"]["ran"]
    assert result["tables"]["shared_pathways_csv"] is None


def test_upstream_blank_symbols_and_significant_only_scope(bulk_inputs):
    path = Path(_rna(bulk_inputs)["tables_dir"]) / "de_genes.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame[frame.padj < .05].copy()
    frame["gene_symbol"] = ""
    frame.to_csv(path, sep="\t", index=False)
    result = cross_modal_concordance(bulk_inputs)
    summary = result["matched_contrasts"][0]["de_da_summary"]
    assert summary["ran"] and summary["n_genes_tested"] == 10
    assert summary["rna_scope"] == "reported_genes_only"
    assert summary["enrichment_pvalue"] is None


def test_auto_cpm_discovery_and_input_signatures(bulk_inputs):
    root = Path(bulk_inputs["output_dir"]).parent
    params = _sample_inputs(root)
    rp = root / "rna" / "counts_cpm.tsv"
    ap = root / "bulk_atac_peak_cpm.tsv"
    Path(params["rna_cpm_path"]).rename(rp)
    Path(params["atac_cpm_path"]).rename(ap)
    # The CPM features share the same biological namespaces and explicit GTF.
    bulk_inputs["gtf_path"] = params["gtf_path"]
    _atac(bulk_inputs)["peak_counts"] = {"counts_cpm_path": str(ap)}
    result = cross_modal_concordance(bulk_inputs)
    assert result["sample_matched_links"]["ran"]
    manifest = json.loads(Path(result["stage_manifest_path"]).read_text())
    assert {str(rp), str(ap)}.issubset({x["path"] for x in manifest["inputs"]})
    bulk_inputs["resume"] = True
    assert cross_modal_concordance(bulk_inputs)["resumed_from_manifest"]
    table = pd.read_csv(rp, sep="\t")
    table.iloc[0, 1] *= 1.01
    table.to_csv(rp, sep="\t", index=False)
    assert not cross_modal_concordance(bulk_inputs)["resumed_from_manifest"]


def test_bh_and_long_donor_cpm(tmp_path):
    np.testing.assert_allclose(bh_adjust([.02, .001, .9]), [.03, .003, .9])
    path = tmp_path / "long.csv"
    pd.DataFrame({"gene": ["G", "G"], "sample_id": ["001", "002"], "cpm": [3, 7]}).to_csv(path, index=False)
    table = read_cpm(str(path), {}, [])
    assert table.loc["G"].to_dict() == {"001": 3, "002": 7}


@pytest.mark.parametrize("key,value", [("min_matched_samples", 2), ("n_permutations", 0), ("padj_threshold", -1), ("seed", -1)])
def test_invalid_parameters_are_honest_nonruns(tmp_path, key, value):
    result = cross_modal_concordance({key: value, "output_dir": str(tmp_path)})
    assert result["status"] == "success" and not result["ran"] and result["reason"] == "invalid_parameters"
