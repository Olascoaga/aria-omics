"""Wave 1: finite statistics, measured evidence and shared utility contracts."""
import json
import math

import numpy as np
import pandas as pd
import pytest

from aria.utils import cross_modal_concordance as cm
from aria.utils import lobo_jackknife as lobo
from aria.utils.batch_qc import assess_hidden_batch
from aria.utils.ora import run_ora
from aria.utils.power_estimation import _as_pair, bulk_power_estimate
from aria.utils.sensitivity import classify_sensitivity
from aria.utils.stats import bh_correct, contrast_family_significance, pooled_bh_across_groups


@pytest.mark.parametrize("correct", [bh_correct, cm.bh_adjust])
def test_bh_nonfinite_clipping_ties_and_empty(correct):
    values = [.01, np.nan, .04, np.inf, -np.inf, -.1, 1.5, .01]
    result = correct(values)
    np.testing.assert_allclose(result, [1 / 60, np.nan, .05, np.nan, np.nan, 0, 1, 1 / 60], equal_nan=True)
    assert correct([]).size == 0
    assert np.isnan(correct([np.nan, np.inf])).all()
    np.testing.assert_allclose(correct([.03, .01, .01]), [.03, .015, .015])


def test_pooled_bh_preserves_discoveries_across_filtered_contrasts():
    groups = {"a": {"hit": .001, "cooks": np.nan, "missing": None},
              "b": {"hit2": .01, "bad": np.inf}}
    adjusted = pooled_bh_across_groups(groups)
    assert adjusted["a"]["hit"] == .002
    assert adjusted["b"]["hit2"] == .01
    assert math.isnan(adjusted["a"]["cooks"])
    assert math.isnan(adjusted["b"]["bad"])
    assert "missing" not in adjusted["a"]
    stats = {c: {g: {"pvalue": p, "log2fc": 2} for g, p in genes.items()}
             for c, genes in groups.items()}
    assert [v["n_sig"] for v in contrast_family_significance(stats, padj_max=.05, lfc_min=0).values()] == [1, 1]


@pytest.mark.parametrize("value", [True, False, np.bool_(True), np.nan, np.inf, -np.inf])
def test_finite_rejects_boolean_and_nonfinite(value):
    assert cm.finite(value) is None


@pytest.mark.parametrize("missing", [None, "", "  ", np.nan, np.inf, -np.inf, "NaN", pd.NA, True])
def test_value_aliases_case_and_missing_fallback(missing):
    assert cm._value({"LoG2Fc": missing, "LOG2FOLDCHANGE": 2.5}, "log2fc", "log2FoldChange") == 2.5
    assert cm._value({"TERM": "Measured term"}, "term") == "Measured term"


@pytest.mark.parametrize("first", ["unreadable", "wrong_contrast", "wrong_group"])
def test_rna_loader_continues_after_empty_candidate(tmp_path, first):
    bad = tmp_path / "bad.csv"
    if first == "unreadable":
        bad.write_text('gene,padj\n"unterminated')
    else:
        pd.DataFrame({"gene": ["G0"], "padj": [.1], "comparison": ["X_vs_Y" if first == "wrong_contrast" else "A_vs_B"],
                      "group": ["other" if first == "wrong_group" else "T"]}).to_csv(bad, index=False)
    pd.DataFrame({"gene": ["G1"], "LOG2FC": [2], "PADJ": [.2]}).to_csv(tmp_path / "de_all_genes.tsv", sep="\t", index=False)
    frame, scope = cm.load_rna_table({"name": "A_vs_B", "group": "T", "full_results_csv": str(bad), "tables_dir": str(tmp_path)}, {}, [])
    assert frame.gene.tolist() == ["G1"]
    assert frame.rna_log2fc.tolist() == [2]
    assert scope == "all_tested_genes"


@pytest.mark.parametrize("conv", [False, 0, 0.0, "false", " False ", "0", " 0.0 "])
def test_atac_rejects_unconverged_peaks(conv):
    comp = {"top_peaks": [{"peak": "chr1:1-2", "nearest_gene": "G", "distance_to_tss": 0,
                            "log2fc": 2, "padj": .01, "lfc_converged": conv}]}
    frame, _ = cm.load_atac_table(comp, {}, {}, {"cis_distance_bp": 100}, False, [])
    assert frame.empty
    comp["top_peaks"][0]["lfc_converged"] = True
    assert len(cm.load_atac_table(comp, {}, {}, {"cis_distance_bp": 100}, False, [])[0]) == 1


def test_directionless_significant_pairs_have_no_directional_fisher():
    rna, _ = cm.load_rna_table({"all_sig": ["G1", "G2"]}, {}, [])
    atac = pd.DataFrame([{"gene_key": g, "peak": f"chr1:{i}-{i+1}", "atac_log2fc": 2,
                          "atac_padj": .01, "distance_to_tss": 0} for i, g in enumerate(["G1", "G2"])])
    _, summary = cm.de_da_layer(rna, atac, {"padj_threshold": .05, "promoter_bp": 2000}, "A_vs_B")
    assert summary["n_both_significant"] == 2
    assert summary["n_concordant_total"] == summary["n_discordant_total"] == 0
    assert summary["fisher_pvalue"] is summary["fisher_odds_ratio"] is None


@pytest.mark.parametrize("comp,foot", [({}, {"group_a": "A", "group_b": "B"}), ({"name": "A_vs_B"}, {}), ({}, {})])
def test_unoriented_footprints_skip(comp, foot):
    warnings = []
    rows, _ = cm._footprint_records({"ran": True, "all_results": [{"tf": "G", "padj": .01, "change": 2}], **foot}, comp, warnings)
    assert rows == [] and any("orientation_unavailable" in w for w in warnings)


@pytest.mark.parametrize("padj", [-.01, np.nan, np.inf, -.5])
def test_tf_and_pathways_reject_invalid_adjusted_pvalues(padj):
    atac = {"motifs": {"per_group": {"up": {"top_motifs": [{"name": "G", "padj": padj}]}}},
            "footprinting": {"ran": True, "group_a": "A", "group_b": "B",
                             "all_results": [{"tf": "G", "padj": padj, "change": 2}]}}
    assert cm.tf_layer(atac, {"name": "A_vs_B"}, pd.DataFrame(), False, .05, [])[0].empty
    assert cm.pathway_terms([{"term": "T", "padj": padj}], .05) == {}


def test_zero_motif_padj_does_not_invent_magnitude(tmp_path):
    rna, _ = cm.load_rna_table({"all_sig": [{"gene": "G", "log2fc": 2, "padj": .01}]}, {}, [])
    atac = {"motifs": {"per_group": {"up": {"top_motifs": [{"name": "G", "padj": 0}]}}}}
    frame, summary = cm.tf_layer(atac, {"name": "A_vs_B"}, rna, False, .05, [])
    assert pd.isna(frame.iloc[0].chromatin_score)
    assert frame.iloc[0].chromatin_padj == 0
    assert summary["top_concordant_tfs"][0]["chromatin_score"] is None
    json.dumps(summary, allow_nan=False)
    cm.render_starburst(pd.DataFrame(), frame, pd.DataFrame(), {}, tmp_path)
    atac["motifs"]["per_group"]["up"]["top_motifs"][0]["padj"] = 1e-5
    assert cm.tf_layer(atac, {"name": "A_vs_B"}, rna, False, .05, [])[0].iloc[0].chromatin_score == 5


def test_pathway_precision_zero_order_and_deterministic_ties():
    pairs = {"tiny": (1e-200, 1e-200), "balanced": (1e-50, 1e-50), "zero": (0, .049),
             "tie_b": (1e-10, 1e-8), "tie_a": (1e-10, 1e-8), "tie_wide": (1e-12, 1e-6)}
    rp = [{"term": term, "padj": p[0]} for term, p in pairs.items()]
    ap = [{"term": term, "padj": p[1]} for term, p in pairs.items()]
    frame, _, rt, at = cm.pathway_layer({}, {"peak_ora": {"pathways": ap}}, {"pathways": rp}, {}, .05)
    assert frame.term.tolist() == ["tiny", "balanced", "tie_a", "tie_b", "tie_wide", "zero"]
    assert frame.iloc[0].geometric_mean_padj == pytest.approx(1e-200, abs=0)
    assert frame.loc[frame.term == "zero", "geometric_mean_padj"].iloc[0] == 0
    aggregated = cm.aggregate_pathways(rt, at, frame.iloc[::-1], True)
    assert [r["term"] for r in aggregated["shared_terms"]] == frame.term.tolist()


@pytest.mark.parametrize("acomp,entry", [({}, {"name": "A_vs_B"}), ({"name": "A_vs_B"}, {}), ({}, {})])
def test_pathway_comparisons_require_orientation(acomp, entry):
    frame, summary, _, _ = cm.pathway_layer({}, {"peak_ora": {"comparisons": [{**entry, "pathways": [{"term": "T", "padj": .01}]}]}},
                                           {"pathways": [{"term": "T", "padj": .01}]}, acomp, .05)
    assert frame.empty and not summary["ran"]


@pytest.mark.parametrize("long_format", [False, True])
def test_cpm_gtf_alias_collisions_aggregate_but_raw_duplicates_fail(tmp_path, long_format):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "ID1"; gene_name "G";\n'
                   'chrY\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "ID2"; gene_name "G";\n')
    aliases = cm.gtf_symbols(gtf)
    path = tmp_path / "cpm.tsv"
    if long_format:
        data = pd.DataFrame({"gene_id": ["ID1", "ID2"], "sample_id": ["001", "001"], "cpm": [2, 3]})
    else:
        data = pd.DataFrame({"gene_id": ["ID1", "ID2"], "001": [2, 3]})
    data.to_csv(path, sep="\t", index=False)
    warnings = []
    result = cm.read_cpm(path, aliases, warnings)
    assert result.loc["G", "001"] == 5
    assert any("summed_mapped_loci" in w for w in warnings)
    data["gene_id"] = "ID1"
    data.to_csv(path, sep="\t", index=False)
    assert cm.read_cpm(path, {}, []).empty
    path.write_text("gene\tS\tS\nG\t1\t2\n")
    assert cm.read_cpm(path, aliases, []).empty


@pytest.mark.parametrize("n", [2, 3, 6, 8])
def test_sample_null_rejects_identity_and_chunks_heldout(monkeypatch, n):
    class Draws:
        def __init__(self):
            self.calls = 0
        def permutation(self, size):
            self.calls += 1
            return np.arange(size) if self.calls % 2 else np.roll(np.arange(size), 1)
    draws = Draws()
    monkeypatch.setattr(cm.np.random, "default_rng", lambda _: draws)
    real_einsum = np.einsum
    sizes = []
    def einsum(expr, left, right):
        sizes.append(len(left))
        return real_einsum(expr, left, right)
    monkeypatch.setattr(cm.np, "einsum", einsum)
    columns = [f"S{i}" for i in range(n)]
    rna = pd.DataFrame([np.arange(n) + 1], index=["G"], columns=columns)
    # Exercise two chunks without allocating a large expression matrix.
    atac = pd.DataFrame(np.tile(np.arange(n) + 1, (4097, 1)),
                        index=[f"chr1:{i}-{i+1}" for i in range(4097)], columns=columns)
    models = {"chr1": [{"gene": "G", "tss": 0, "strand": "+"}]}
    frame, summary = cm.sample_matched_layer(rna, atac, models, {"cis_distance_bp": 10000, "min_matched_samples": 2,
                           "seed": 0, "n_permutations": 3, "padj_threshold": .05}, [])
    assert len(frame) == 4097 and draws.calls == 8
    assert sizes[-2:] == [4096, 1]
    assert max(sizes[2:]) <= 4096
    assert "held_out_discovery_fraction" in summary["permutation_null"]


def _lobo_inputs():
    units = [f"S{i}" for i in range(6)]
    groups = dict(zip(units, ["A"] * 3 + ["B"] * 3))
    matrix = np.random.default_rng(3).normal(size=(25, 6))
    return units, groups, matrix


@pytest.mark.parametrize("missing", [None, pd.NA, np.nan, "", " ", "nAn", "None", "NA", "null"])
def test_lobo_missing_donors_never_pair(missing):
    units, groups, matrix = _lobo_inputs()
    donors = {u: missing for u in units}
    result = lobo.compute_lobo_jackknife([f"G{i}" for i in range(25)], matrix, units, groups,
                                       numerator="A", denominator="B", unit_to_donor=donors)
    assert result["status"] == "measured"
    assert result["design_mode"] == "unpaired_biological_replicate"
    counts = pd.DataFrame(np.exp2(matrix), index=[f"G{i}" for i in range(25)], columns=units)
    meta = pd.DataFrame({"condition": [groups[u] for u in units], "donor": [missing] * 6}, index=units)
    wrapped = lobo.compute_counts_lobo_jackknife(counts, meta, "condition", "A", "B", paired_unit_col="donor")
    assert wrapped["design_mode"] == "unpaired_biological_replicate"


@pytest.mark.parametrize("sf", [None, [1] * 6, [0] * 6])
def test_lobo_sanitizes_counts_before_normalizing(monkeypatch, sf):
    units, groups, matrix = _lobo_inputs()
    raw = np.exp2(matrix)
    raw[0, :4] = [np.nan, -2, np.inf, -np.inf]
    clean = np.where(np.isfinite(raw), np.maximum(raw, 0), 0)
    def capture(**kwargs):
        actual = kwargs["log_norm_matrix"]
        expected = np.log2(clean + 1) if sf and sf[0] else np.log2(clean / np.maximum(clean.sum(axis=0), 1) * 1e6 + 1)
        np.testing.assert_allclose(actual, expected)
        assert np.isfinite(actual).all()
        return {}
    monkeypatch.setattr(lobo, "compute_lobo_jackknife", capture)
    lobo.compute_counts_lobo_jackknife(pd.DataFrame(raw, columns=units),
                                     pd.DataFrame({"condition": list(groups.values())}, index=units),
                                     "condition", "A", "B", size_factors=sf)


def test_lobo_deduplicates_target_features():
    units, groups, matrix = _lobo_inputs()
    result = lobo.compute_lobo_jackknife([f"G{i}" for i in range(25)], matrix, units, groups,
                                       numerator="A", denominator="B", target_features=["G2", "G2", "G1", "missing"],
                                       max_reported_features=1)
    assert result["n_features_evaluated"] == 2
    assert result["n_features_lobo_stable"] + result["n_single_unit_driven"] == 2
    assert list(result["per_feature"]) == ["G2"] and result["per_feature_truncated"]


@pytest.mark.parametrize("broken", ["fold_lfc", "fold_padj", "full_lfc", "full_padj"])
@pytest.mark.parametrize("primary", [None, ["G"]])
def test_lobo_nonfinite_stats_are_unstable_and_strict_json(monkeypatch, broken, primary):
    units, groups, _ = _lobo_inputs()
    calls = []
    def stats(*args):
        calls.append(1)
        full = len(calls) == 1
        lfc, padj = np.array([2.]), np.array([.001])
        if broken == ("full_lfc" if full else "fold_lfc"):
            lfc[0] = np.nan
        if broken == ("full_padj" if full else "fold_padj"):
            padj[0] = np.inf
        return lfc, padj, padj
    monkeypatch.setattr(lobo, "_unpaired_fold_stats", stats)
    result = lobo.compute_lobo_jackknife(["G"], np.ones((1, 6)), units, groups, numerator="A", denominator="B",
                                       primary_significant_features=primary)
    assert result["per_feature"]["G"]["single_unit_driven"] is True
    assert result["n_single_unit_driven"] == 1 and result["n_features_lobo_stable"] == 0
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("value", [[4], (4,), {"A": 4}, 4])
def test_power_single_group_containers(value):
    assert _as_pair(value) == (4., 4.)
    assert bulk_power_estimate(value) == bulk_power_estimate(4)


@pytest.mark.parametrize("key", ["n_per_group", "target_log2fc", "alpha"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, True, False])
def test_power_invalid_inputs_return_zero(key, value):
    args = {"n_per_group": 4, key: value}
    assert bulk_power_estimate(**args) == 0


@pytest.mark.parametrize("key", ["dispersion", "mean_expression"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, True, False])
def test_power_nuisance_inputs_are_safe(key, value):
    result = bulk_power_estimate(4, **{key: value})
    assert math.isfinite(result) and 0 <= result <= 1
    if isinstance(value, bool):
        assert result == 0


def test_power_zero_alpha_zero_effect_and_none_defaults():
    assert bulk_power_estimate(4, alpha=0) == bulk_power_estimate(4, alpha=-.1) == 0
    assert bulk_power_estimate(4, target_log2fc=0) == .05
    assert bulk_power_estimate(4, alpha=None) == bulk_power_estimate(4)
    assert math.isfinite(bulk_power_estimate(4, mean_expression=None, dispersion=None, target_log2fc=None))


@pytest.mark.parametrize("name", ["mRNA", "pseudobulk", "neutrophil", "graphics", "stage", "passage", "coverage", "lineage", "macrophage"])
def test_sensitivity_avoids_substrings(name):
    result = classify_sensitivity("human", [name, "sex"])
    assert result["level"] == "elevated"
    assert result["phi_signals"] == result["quasi_identifier_signals"] == []


@pytest.mark.parametrize("name,token", [("patient_name", "patient"), ("mrn", "mrn"), ("mrn_id", "mrn"),
    ("dob", "dob"), ("patientDOB", "dob"), ("patient_dob", "dob"), ("ssn", "ssn"),
    ("email", "email"), ("phone", "phone"), ("ICD_code", "icd"), ("PHI", "phi"),
    ("medical_record", "medicalrecord"), ("dateOfBirth", "dateofbirth"), ("birthDate", "birthdate")])
def test_sensitivity_keeps_real_phi(name, token):
    result = classify_sensitivity("human", [name])
    assert result["level"] == "high" and token in result["phi_signals"]


def test_sensitivity_quasi_identifiers_and_path_boundaries():
    assert classify_sensitivity("human", ["subjectAge", "gender"])["quasi_identifier_signals"] == ["age+sex"]
    assert classify_sensitivity(path_hints=["/data/clinical/patientDOB.csv"])["level"] == "high"


@pytest.mark.parametrize("container", [pd.Index, pd.Series, np.array])
@pytest.mark.parametrize("background", [None, [], [" a ", "B", "C", "D", "E"]])
def test_ora_array_inputs_and_normalized_symbols(container, background):
    bg = None if background is None else container(background)
    actual = run_ora(container([" a ", "b"]), {"T": [" A", "b "]}, bg, padj_max=1.01)
    expected = run_ora(["A", "B"], {"T": ["A", "B"]}, background, padj_max=1.01)
    assert actual == expected and actual[0]["genes"] == ["A", "B"]


def test_ora_overlap_filter_does_not_shrink_bh_family():
    query = ["G0", "G1"]
    sets = {"both": query, "one": ["G0"], "none": ["G9"]}
    universe = [f"G{i}" for i in range(100)]
    all_terms = run_ora(query, sets, universe, padj_max=1.01, min_overlap=1)
    filtered = run_ora(query, sets, universe, padj_max=1.01, min_overlap=2)
    assert [r["term"] for r in filtered] == ["both"]
    assert filtered[0]["padj"] == next(r["padj"] for r in all_terms if r["term"] == "both")
    assert filtered[0]["padj"] == pytest.approx(filtered[0]["pvalue"] * 2)


@pytest.mark.parametrize("col", ["runx1_score", "cut_and_run", "pruned", "ion_channel_score"])
def test_batch_false_positives(col):
    assert assess_hidden_batch([col])["status"] == "clean"


@pytest.mark.parametrize("col", ["batch", "run_id", "lane", "seq_run", "channel_1", "donor_batch", "runID", "sequencingLane", "flowCell"])
def test_batch_true_tokens(col):
    assert assess_hidden_batch([col])["candidate_batch_columns"] == [col]


@pytest.mark.parametrize("modeled", ["condition_col", "replicate_col", "declared_batch"])
def test_batch_modeled_names_case_insensitive(modeled):
    assert assess_hidden_batch(["BATCH"], **{modeled: "batch"})["status"] == "clean"
