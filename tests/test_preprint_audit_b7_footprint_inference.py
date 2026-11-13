"""Preprint-readiness blocker B7: replicate-aware footprint inference.

The inferential unit must be an explicit biological replicate/donor. TOBIAS
site-level p-values are never promoted to significance. Replicate mean scores
are tested, adjusted across motifs with BH, and stress-tested by label
permutations. Under-replicated inputs remain descriptive.
"""
from __future__ import annotations

import math
import json
from pathlib import Path

import pytest


def _synthetic_scores(n_per_group: int = 6) -> tuple[dict, dict]:
    labels = {}
    for group in ("state_a", "state_b"):
        for idx in range(n_per_group):
            labels[f"{group}_r{idx + 1}"] = group

    table = {}
    for idx in range(40):
        # Twenty strong effects plus twenty deterministic null motifs. Keeping
        # multiple true effects also exercises BH rather than a single-test path.
        effect = 2.5 if idx < 20 else 0.0
        scores = {}
        for rep_idx in range(n_per_group):
            jitter = (rep_idx - 2.5) * 0.04
            scores[f"state_a_r{rep_idx + 1}"] = effect + jitter + idx * 0.001
            scores[f"state_b_r{rep_idx + 1}"] = -jitter + idx * 0.001
        table[f"motif_{idx:02d}"] = {"scores": scores, "n_sites": 100 + idx}
    return table, labels


def test_bh_adjustment_is_monotone_and_bounded():
    from aria.scripts.chromatin_footprint_tobias import _bh_adjust

    adjusted = _bh_adjust([0.01, 0.04, 0.03, math.nan, 1.0])
    assert adjusted[:3] == [0.04, 0.05333333333333334, 0.05333333333333334]
    assert math.isnan(adjusted[3])
    assert adjusted[4] == 1.0


def test_parse_bindetect_uses_replicate_mean_scores_not_site_pvalues(tmp_path: Path):
    from aria.scripts.chromatin_footprint_tobias import (
        parse_bindetect_replicate_scores,
    )

    path = tmp_path / "bindetect_results.txt"
    path.write_text(
        "name\ttotal_tfbs\ta1_mean_score\ta2_mean_score\tb1_mean_score\t"
        "b2_mean_score\ta1_b1_pvalue\n"
        "TF_A\t120\t2.0\t2.2\t0.1\t0.2\t1e-200\n",
        encoding="utf-8",
    )
    parsed = parse_bindetect_replicate_scores(
        str(path), {"a1": "state_a", "a2": "state_a", "b1": "state_b", "b2": "state_b"}
    )
    assert parsed["TF_A"]["scores"] == {
        "a1": 2.0, "a2": 2.2, "b1": 0.1, "b2": 0.2,
    }
    assert parsed["TF_A"]["tf"] == "TF_A"
    assert parsed["TF_A"]["motif_id"] == "TF_A"
    assert parsed["TF_A"]["n_sites"] == 120
    assert all("pvalue" not in key for key in parsed["TF_A"])


def test_parse_bindetect_preserves_duplicate_tf_names_as_distinct_motifs(tmp_path: Path):
    from aria.scripts.chromatin_footprint_tobias import (
        parse_bindetect_replicate_scores,
    )

    path = tmp_path / "bindetect_results.txt"
    path.write_text(
        "name\tmotif_id\ttotal_tfbs\ta1_mean_score\tb1_mean_score\n"
        "TF_FAMILY\tMA0001.1\t100\t1.0\t0.0\n"
        "TF_FAMILY\tMA0002.1\t120\t1.1\t0.1\n",
        encoding="utf-8",
    )
    parsed = parse_bindetect_replicate_scores(
        str(path), {"a1": "A", "b1": "B"})
    assert set(parsed) == {"MA0001.1", "MA0002.1"}
    assert {row["tf"] for row in parsed.values()} == {"TF_FAMILY"}


def test_replicate_inference_applies_bh_and_reports_null_controls():
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    table, labels = _synthetic_scores()
    result = infer_replicate_footprints(
        table, labels, "state_a", "state_b", alpha=0.05,
        min_replicates_per_condition=3, max_label_permutations=40, random_seed=7,
    )
    assert result["parsed"] is True
    assert result["inference"]["status"] == "success"
    assert result["inference"]["inferential_unit"] == "biological_replicate_or_donor"
    assert result["inference"]["test"] == "welch_t_test_on_replicate_mean_scores"
    assert result["inference"]["multiple_testing"] == "benjamini_hochberg_across_motifs"
    assert result["n_significant"] == 20
    assert result["n_motifs_tested"] == 40
    assert result["top_toward_state_a"][0]["padj"] <= 0.05
    null = result["inference"]["null_label_controls"]
    assert null["n_permutations"] > 0
    assert 0.0 <= null["mean_discovery_fraction"] <= 1.0
    assert null["observed_labels_excluded"] is True


def test_null_scores_have_calibrated_observed_fdr_and_permutation_diagnostic():
    import numpy as np

    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    labels = {f"a{i}": "A" for i in range(8)} | {f"b{i}": "B" for i in range(8)}
    table = {}
    reps = list(labels)
    null_scores = np.random.default_rng(0).normal(size=(120, len(reps)))
    for motif_idx in range(120):
        table[f"null_{motif_idx:03d}"] = {
            "scores": dict(zip(reps, null_scores[motif_idx])), "n_sites": 200,
        }
    result = infer_replicate_footprints(
        table, labels, "A", "B", alpha=0.05,
        min_replicates_per_condition=3, max_label_permutations=50, random_seed=11,
    )
    assert result["n_significant"] == 0
    assert result["inference"]["null_label_controls"]["mean_discovery_fraction"] <= 0.05


def test_insufficient_replicates_downgrade_without_significance_language():
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    table, labels = _synthetic_scores(n_per_group=2)
    result = infer_replicate_footprints(
        table, labels, "state_a", "state_b", min_replicates_per_condition=3,
    )
    assert result["inference"]["status"] == "descriptive_only"
    assert result["inference"]["reason"] == "insufficient_biological_replicates"
    assert "n_significant" not in result
    assert result["ranking_basis"]["fdr_controlled"] is False


def test_complete_donor_pairs_use_paired_test_and_within_donor_label_swaps():
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    groups = {}
    replicate_ids = {}
    scores = {}
    donor_baselines = [0.0, 5.0, -4.0, 8.0, -7.0, 3.0]
    for idx, baseline in enumerate(donor_baselines, start=1):
        a, b = f"A_signal_{idx}", f"B_signal_{idx}"
        groups[a], groups[b] = "A", "B"
        replicate_ids[a] = replicate_ids[b] = f"donor_{idx}"
        scores[a] = baseline + 1.0 + idx * 0.02
        scores[b] = baseline
    result = infer_replicate_footprints(
        {"TF_PAIRED": {"scores": scores, "n_sites": 500}},
        groups, "A", "B", min_replicates_per_condition=3,
        max_label_permutations=30, replicate_ids=replicate_ids,
    )
    inference = result["inference"]
    assert inference["pairing_status"] == "complete"
    assert inference["test"] == "paired_t_test_on_within_donor_mean_score_differences"
    assert inference["null_label_controls"]["method"] == (
        "within_donor_condition_label_swap_with_bh")
    assert result["n_significant"] == 1


def test_barcode_design_requires_explicit_replicate_column_for_inference(tmp_path: Path):
    from aria.scripts.chromatin_footprint_tobias import _load_barcode_design

    design = tmp_path / "barcode_design.tsv"
    design.write_text(
        "barcode\tgroup\treplicate\nBC1\tA\tdonor_1\nBC2\tB\tdonor_2\n",
        encoding="utf-8",
    )
    assert _load_barcode_design(str(design)) == {
        "BC1": {"group": "A", "replicate": "donor_1"},
        "BC2": {"group": "B", "replicate": "donor_2"},
    }


def test_narrator_distinguishes_bh_inference_from_descriptive_ranking():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    footprinting = {
        "ran": True,
        "group_a": "state_a",
        "group_b": "state_b",
        "group_label": "Conditions",
        "group_kind": "conditions",
        "differential_summary": {
            "parsed": True,
            "n_motifs_tested": 100,
            "n_significant": 4,
            "inference": {
                "status": "success",
                "inferential_unit": "biological_replicate_or_donor",
                "test": "welch_t_test_on_replicate_mean_scores",
                "multiple_testing": "benjamini_hochberg_across_motifs",
                "replicates_per_condition": {"state_a": 4, "state_b": 4},
                "null_label_controls": {"n_permutations": 20,
                    "mean_discovery_fraction": 0.01},
            },
        },
    }
    block = ChromatinNarrator()._footprint_block(footprinting)
    low = block.claim.lower()
    assert "4 of 100" in low
    assert "benjamini–hochberg" in low
    assert "biological replicate" in low
    assert "descriptively ranked" not in low
    assert any(e.label.startswith("BH-significant") for e in block.evidence)
    assert "associative" in " ".join(c.text.lower() for c in block.caveats)


def test_bulk_agent_groups_technical_bams_by_explicit_biological_replicate():
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    design = agent._bulk_replicate_bams(
        {
            "sample_ids": ["s1_lane1", "s1_lane2", "s2"],
            "sample_metadata": {
                "s1_lane1": {"condition": "A", "replicate": "donor_1"},
                "s1_lane2": {"condition": "A", "replicate": "donor_1"},
                "s2": {"condition": "B", "replicate": "donor_2"},
            },
        },
        ["lane1.bam", "lane2.bam", "s2.bam"],
    )
    assert design == {
        "A": {"donor_1": ["lane1.bam", "lane2.bam"]},
        "B": {"donor_2": ["s2.bam"]},
    }
    # Filename-derived sample identity is forbidden for the inferential design.
    assert agent._bulk_replicate_bams(
        {"sample_metadata": {"lane1": {"condition": "A", "replicate": "r1"}}},
        ["lane1.bam"],
    ) == {}


def test_bulk_driver_selects_replicate_route_and_writes_bh_table(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    genome = tmp_path / "genome.fa"
    peaks = tmp_path / "peaks.bed"
    motifs = tmp_path / "motifs.meme"
    for path in (genome, peaks, motifs, tmp_path / "genome.fa.fai"):
        path.write_text("x\n", encoding="utf-8")
    replicate_groups = {
        "A_rep1_r1": "A", "A_rep2_r2": "A",
        "B_rep1_r1": "B", "B_rep2_r2": "B",
    }
    replicate_info = {
        name: {"bam": str(tmp_path / f"{name}.bam"), "n_fragments": 100}
        for name in replicate_groups
    }
    monkeypatch.setattr(mod.shutil, "which", lambda _: "/usr/bin/TOBIAS")
    monkeypatch.setattr(
        mod, "_prepare_bulk_replicate_bams",
        lambda *args, **kwargs: (replicate_info, replicate_groups,
                       {name: name for name in replicate_groups}),
    )

    result_table = tmp_path / "bindetect_results.txt"
    result_table.write_text(
        "name\ttotal_tfbs\tA_rep1_r1_mean_score\tA_rep2_r2_mean_score\t"
        "B_rep1_r1_mean_score\tB_rep2_r2_mean_score\n"
        "TF_A\t100\t2.0\t2.1\t0.0\t0.1\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        mod, "_tobias_replicate_pipeline",
        lambda *args, **kwargs: {
            "bindetect_results": str(result_table),
            "bindetect_outdir": str(tmp_path / "bindetect"),
            "replicate_groups": replicate_groups,
            "corrected_signals": {},
            "footprint_signals": {},
        },
    )
    monkeypatch.setattr(mod, "_aggregate_plots", lambda *args, **kwargs: {})
    result = mod.chromatin_footprint_tobias_bulk({
        "replicate_bams": {
            "A": {"r1": ["a1.bam"], "r2": ["a2.bam"]},
            "B": {"r1": ["b1.bam"], "r2": ["b2.bam"]},
        },
        "genome_fasta": str(genome), "peaks_bed": str(peaks),
        "motif_meme": str(motifs), "group_a": "A", "group_b": "B",
        "output_dir": str(tmp_path), "min_replicates_per_condition": 2,
    })
    summary = result["differential_summary"]
    assert result["method"].endswith("welch_bh")
    assert summary["inference"]["status"] == "success"
    assert summary["inference"]["low_power_warning"] is True
    assert summary["n_significant"] == 1
    assert Path(summary["results_table"]).is_file()


def test_historical_pooled_artifacts_are_reclassified_as_descriptive():
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "docs/benchmark_results/bulk_atac_footprint/"
               "b4_bulk_footprint_tobias_bindetect.json",
        root / "docs/benchmark_results/scatac_footprint/"
               "p4_footprint_tobias_bindetect.json",
    ]
    for path in paths:
        summary = json.loads(path.read_text(encoding="utf-8"))["differential_summary"]
        assert "n_significant" not in summary
        assert summary["n_ranked_candidates"] > 0
        assert summary["ranking_basis"]["fdr_controlled"] is False
        assert summary["inference"]["status"] == "descriptive_only"


def test_standalone_bulk_wrapper_propagates_repo_pythonpath(tmp_path, monkeypatch):
    from scripts import run_bulk_atac_footprint_tobias as wrapper

    captured = []

    def fake_run(command, **kwargs):
        captured.append({"command": command, "env": kwargs.get("env") or {}})

    monkeypatch.setattr(wrapper.subprocess, "run", fake_run)
    result = wrapper.main([
        "--replicate-bams", str(tmp_path / "design.json"),
        "--genome-fasta", str(tmp_path / "genome.fa"),
        "--peaks-bed", str(tmp_path / "peaks.bed"),
        "--motif-meme", str(tmp_path / "motifs.meme"),
        "--group-a", "A", "--group-b", "B",
        "--work-dir", str(tmp_path / "work"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert result == 0
    child = next(call for call in captured
                 if "aria-tobias-env" in call["command"])
    assert str(wrapper.ROOT) in child["env"]["PYTHONPATH"].split(":")
    manifest = json.loads((tmp_path / "out" /
                           "b4_bulk_footprint_tobias_bindetect.json").read_text())
    assert manifest["provenance"]["environment"]["conda_prefix"] is None


def _batch_scores():
    import numpy as np

    groups = {f"a{i}": "A" for i in range(6)} | {f"b{i}": "B" for i in range(6)}
    # A is enriched in the low batch; the raw contrast hides a true positive effect.
    batches = dict(zip(groups, ["low"] * 4 + ["high"] * 2
                       + ["low"] * 2 + ["high"] * 4))
    rng = np.random.default_rng(42)
    table = {}
    for i, effect in enumerate([2.0, -2.0, 0.0]):
        table[f"m{i}"] = {"scores": {
            rep: (1000.0 if batches[rep] == "high" else 0.0)
            + (effect if group == "A" else 0.0) + rng.normal(0, 0.06)
            for rep, group in groups.items()
        }}
    return table, groups, batches


def test_batch_adjusted_ols_matches_numpy_scipy_and_recovers_masked_signal():
    import numpy as np
    from scipy.stats import false_discovery_control, t
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    table, groups, batches = _batch_scores()
    result = infer_replicate_footprints(table, groups, "A", "B", replicate_batches=batches)
    unadjusted = infer_replicate_footprints(table, groups, "A", "B")
    assert unadjusted["n_significant"] == 0
    x = np.array([[1, batches[r] == "low", groups[r] == "A"] for r in groups], float)
    covariance = np.linalg.inv(x.T @ x)
    expected_p = []
    for row in result["all_results"]:
        y = np.array([table[row["motif_id"]]["scores"][r] for r in groups])
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
        df = len(y) - x.shape[1]
        variance = np.sum((y - x @ beta) ** 2) / df
        pvalue = 2 * t.sf(abs(beta[-1] / np.sqrt(variance * covariance[-1, -1])), df)
        expected_p.append(pvalue)
        assert row["pvalue"] == pytest.approx(pvalue, rel=1e-8)
        assert row["change"] == round(beta[-1], 4)
        assert row["raw_change"] == round(y[:6].mean() - y[6:].mean(), 4)
        assert row["mean_score_a"] == pytest.approx(y[:6].mean())
        assert row["mean_score_b"] == pytest.approx(y[6:].mean())
    assert [r["padj"] for r in result["all_results"]] == pytest.approx(
        false_discovery_control(expected_p), rel=1e-8)
    assert result["all_results"][0]["change"] > 1.9
    assert result["all_results"][0]["raw_change"] < -300
    assert result["n_significant"] == 2
    inference = result["inference"]
    assert inference["test"] == "ols_t_test_on_replicate_mean_scores_batch_adjusted"
    assert inference["batch_adjusted"] is True
    assert inference["df_residual"] == 9
    assert inference["null_label_controls"]["method"] == (
        "within_batch_condition_label_permutation_with_bh")


@pytest.mark.parametrize("limit", [100, 5, 0])
def test_within_batch_permutations_preserve_counts_and_exclude_observed(monkeypatch, limit):
    import itertools
    import numpy as np
    from aria.scripts import chromatin_footprint_tobias as mod

    batches = ["x", "x", "y", "y", "x", "x", "y", "y"]
    matrix = np.random.default_rng(1).normal(size=(4, 8))
    designs, expected_fractions = [], []
    original = mod._ols_group_test

    def capture(data, x):
        designs.append(x.copy())
        result = original(data, x)
        expected_fractions.append(sum(p <= 0.05 for p in mod._bh_adjust(result[1])) / 4)
        return result

    monkeypatch.setattr(mod, "_ols_group_test", capture)
    result = mod._within_batch_null_controls(matrix, 4, batches, 0.05, limit, 7)
    assert result["n_permutations"] == min(limit, 34)  # 6 * 6 - observed - complement
    assignments = set()
    for x in designs:
        selected = tuple(np.flatnonzero(x[:, -1]))
        assert selected not in [(0, 1, 2, 3), (4, 5, 6, 7)]
        assert len(selected) == 4
        assert sum(batches[i] == "x" for i in selected) == 2
        assert sum(batches[i] == "y" for i in selected) == 2
        assert np.linalg.matrix_rank(x) == 3
        assignments.add(selected)
    assert len(assignments) == len(designs)
    if limit == 100:
        expected = {c for c in itertools.combinations(range(8), 4)
                    if sum(batches[i] == "x" for i in c) == 2}
        expected -= {(0, 1, 2, 3), (4, 5, 6, 7)}
        assert assignments == expected
    assert result["mean_discovery_fraction"] == pytest.approx(
        np.mean(expected_fractions) if expected_fractions else 0)
    assert result["max_discovery_fraction"] == max(expected_fractions, default=0)
    assert result["fraction_with_any_discoveries"] == pytest.approx(
        np.mean(np.array(expected_fractions) > 0) if expected_fractions else 0)
    assert result["observed_labels_excluded"] is True
    assert result["role"] == "diagnostic_only_not_a_post_hoc_significance_gate"


@pytest.mark.parametrize("case", ["disjoint", "no_residual_df"])
def test_batch_confounding_is_honestly_descriptive(case):
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    table, groups, _ = _batch_scores()
    if case == "disjoint":
        batches = {rep: group for rep, group in groups.items()}
    else:
        # One shared batch makes the contrast identifiable, but leaves no residual df.
        batches = {rep: rep for rep in groups}
        batches["a0"] = batches["b0"] = "shared"
    result = infer_replicate_footprints(table, groups, "A", "B", replicate_batches=batches)
    assert result["parsed"] is True
    assert result["inference"]["status"] == "descriptive_only"
    assert result["inference"]["reason"] == "batch_condition_confounded"
    assert result["ranking_basis"]["fdr_controlled"] is False
    assert result["n_ranked_candidates"] == len(table)
    assert "n_significant" not in result
    assert all("padj" not in row for row in result["top_toward_A"])


@pytest.mark.parametrize("batch_mode", ["none", "empty", "single", "missing", "blank", "paired"])
def test_batch_metadata_keeps_legacy_inference_when_inapplicable(batch_mode):
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    table, groups, batches = _batch_scores()
    identities = None
    if batch_mode == "none":
        batches = None
    elif batch_mode == "empty":
        batches = {}
    elif batch_mode == "single":
        batches = {rep: "one" for rep in groups}
    elif batch_mode == "missing":
        batches.pop("a0")
    elif batch_mode == "blank":
        batches["a0"] = " "
    elif batch_mode == "paired":
        identities = {rep: rep[1:] for rep in groups}
    baseline = infer_replicate_footprints(table, groups, "A", "B", replicate_ids=identities)
    result = infer_replicate_footprints(table, groups, "A", "B", replicate_ids=identities,
                                      replicate_batches=batches)
    assert result == baseline


def test_ols_zero_standard_error_has_finite_explicit_pvalues():
    import numpy as np
    from aria.scripts.chromatin_footprint_tobias import _batch_design, _ols_group_test

    x = _batch_design(["a", "b", "a", "b", "a", "b", "a", "b"], range(4))
    _, pvalues, _ = _ols_group_test(np.array([np.zeros(8), x[:, -1] * 2]), x)
    assert pvalues == [1.0, 0.0]


@pytest.mark.parametrize("header", ["batch", "batch_id", "pool", "lane", None])
def test_barcode_design_reads_batch_aliases_and_headerless_column(tmp_path, header):
    from aria.scripts.chromatin_footprint_tobias import _load_barcode_design

    path = tmp_path / "design.tsv"
    path.write_text((f"barcode\tgroup\tdonor\t{header}\n" if header else "")
                    + "bc1\tA\td1\tpool1\nbc2\tB\td2\t\n")
    assert _load_barcode_design(str(path)) == {
        "bc1": {"group": "A", "replicate": "d1", "batch": "pool1"},
        "bc2": {"group": "B", "replicate": "d2", "batch": None},
    }


@pytest.mark.parametrize("mode,batch_source", [
    ("scatac", "design"), ("scatac", "replicate"), ("scatac", "signal"),
    ("scatac", "inconsistent_design"), ("bulk", "replicate"), ("bulk", "nested"),
])
def test_drivers_plumb_batches_and_inferential_cache_requires_results_table(
    tmp_path, monkeypatch, mode, batch_source,
):
    from aria.scripts import chromatin_footprint_tobias as mod

    for name in ["genome.fa", "peaks.bed", "motifs.meme", "fragments.tsv"]:
        (tmp_path / name).write_text("input\n")
    (tmp_path / "genome.fa.fai").write_text("chr1\t1000\n")
    design = tmp_path / "design.tsv"
    design_rows = ["barcode\tgroup\treplicate\tbatch"]
    replicate_design, groups, ids, batches, by_replicate = {}, {}, {}, {}, {}
    for group in ["A", "B"]:
        replicate_design[group] = {}
        for i in range(3):
            replicate_id = f"{group}_donor{i}"
            batch = "lane1" if i % 2 == 0 else "lane2"
            signal = mod._safe_signal_label(group, replicate_id, i + 1)
            groups[signal], ids[signal], batches[signal] = group, replicate_id, batch
            by_replicate[replicate_id] = batch
            path = tmp_path / f"{replicate_id}.cram"
            path.write_bytes(b"cram")
            replicate_design[group][replicate_id] = [str(path)]
            design_rows.append(f"{replicate_id}#BC1\t{group}\t{replicate_id}\t{batch}")
    if batch_source == "inconsistent_design":
        design_rows.append("other_barcode\tA\tA_donor0\tother_batch")
    design.write_text("\n".join(design_rows) + "\n")
    params = {
        "mode": mode, "genome_fasta": str(tmp_path / "genome.fa"),
        "peaks_bed": str(tmp_path / "peaks.bed"),
        "motif_meme": str(tmp_path / "motifs.meme"),
        "group_a": "A", "group_b": "B", "output_dir": str(tmp_path / "out"),
        "fragments_file": str(tmp_path / "fragments.tsv"), "barcode_groups": str(design),
        "replicate_bams": replicate_design,
        "make_volcano_plot": False, "make_aggregate_plots": False,
    }
    if batch_source == "replicate":
        params["replicate_batches"] = by_replicate
    elif batch_source == "signal":
        params["replicate_batches"] = batches
    elif batch_source == "nested":
        params["replicate_batches"] = {
            group: {ids[s]: batches[s] for s in groups if groups[s] == group}
            for group in ["A", "B"]}
    monkeypatch.setattr(mod.shutil, "which", lambda _: "/bin/TOBIAS")

    def bam_info(names):
        return {name: {"bam": str(tmp_path / f"{name}.bam"), "n_fragments": 100}
                for name in names}

    def prepare_sc(files, barcode_signals, target_groups, chrom_sizes, out_dir, **kwargs):
        assert set(barcode_signals.values()) == set(groups)
        assert set(target_groups) == set(groups)
        assert kwargs["reuse"] is False
        return bam_info(target_groups)

    def merge_bulk(inputs, out_dir, **kwargs):
        assert kwargs == {"reference_fasta": params["genome_fasta"], "reuse": False}
        return bam_info(inputs)

    monkeypatch.setattr(mod, "_write_group_bams", prepare_sc)
    monkeypatch.setattr(mod, "_merge_condition_bams", merge_bulk)
    pipeline_calls = []

    def pipeline(bams, replicate_groups, *args, **kwargs):
        pipeline_calls.append(True)
        assert replicate_groups == groups
        path = tmp_path / "out" / "bindetect_results.txt"
        scores = [str((10 if batches[s] == "lane2" else 0)
                      + (2 if groups[s] == "A" else 0) + i * 0.01)
                  for i, s in enumerate(groups)]
        path.write_text("name\ttotal_tfbs\t" + "\t".join(f"{s}_mean_score" for s in groups)
                        + "\nTF\t100\t" + "\t".join(scores) + "\n")
        return {"replicate_groups": groups, "bindetect_results": str(path)}

    monkeypatch.setattr(mod, "_tobias_replicate_pipeline", pipeline)
    driver = mod.chromatin_footprint_tobias_dispatch
    result = driver(params)
    if batch_source == "inconsistent_design":
        assert "replicate_batches" not in result
        assert result["differential_summary"]["inference"]["test"] == (
            "welch_t_test_on_replicate_mean_scores")
    else:
        assert result["replicate_batches"] == batches
        assert result["differential_summary"]["inference"]["batch_adjusted"] is True
        assert "raw_change" in Path(result["differential_summary"]["results_table"]).read_text()
    assert result["validation_level"] == "production"
    assert driver(params)["reused_from_stage_manifest"] is True
    assert len(pipeline_calls) == 1
    Path(result["differential_summary"]["results_table"]).unlink()
    rerun = driver(params)
    assert len(pipeline_calls) == 2
    assert "reused_from_stage_manifest" not in rerun
    assert Path(rerun["differential_summary"]["results_table"]).is_file()


def test_bulk_driver_unpaired_by_default_for_shared_replicate_labels_and_opt_in_paired(
    tmp_path, monkeypatch,
):
    """In bulk ATAC, replicate_bams commonly uses generic labels like r1/r2/r3
    inside each condition. By default they must be modeled as independent
    condition::replicate units (allowing Welch or batch-adjusted OLS), while
    paired_replicates=True opts into within-donor paired t-tests."""
    from aria.scripts import chromatin_footprint_tobias as mod

    for name in ["genome.fa", "peaks.bed", "motifs.meme"]:
        (tmp_path / name).write_text("input\n", encoding="utf-8")
    (tmp_path / "genome.fa.fai").write_text("chr1\t1000\n", encoding="utf-8")

    replicate_design = {
        "Treated": {"r1": ["t1.bam"], "r2": ["t2.bam"], "r3": ["t3.bam"]},
        "Control": {"r1": ["c1.bam"], "r2": ["c2.bam"], "r3": ["c3.bam"]},
    }
    monkeypatch.setattr(mod.shutil, "which", lambda _: "/usr/bin/TOBIAS")
    monkeypatch.setattr(
        mod, "_merge_condition_bams",
        lambda inputs, *args, **kwargs: {
            k: {"bam": str(tmp_path / f"{k}.bam"), "n_fragments": 100, "n_bams": 1}
            for k in inputs
        },
    )

    def fake_pipeline(bams, replicate_groups, *args, **kwargs):
        path = tmp_path / "bindetect_results.txt"
        cols = list(replicate_groups)
        # Scores with identical delta across r1..r3 (zero paired variance if
        # falsely paired!) plus a batch shift on r3.
        vals = []
        for s in cols:
            is_t = replicate_groups[s] == "Treated"
            idx = 0 if "r1" in s else (1 if "r2" in s else 2)
            batch_shift = 5.0 if idx == 2 else 0.0
            vals.append(f"{(3.0 if is_t else 0.2) + idx * 0.15 + batch_shift:.4f}")
        path.write_text(
            "name\ttotal_tfbs\t" + "\t".join(f"{s}_mean_score" for s in cols)
            + "\nTF_X\t120\t" + "\t".join(vals) + "\n",
            encoding="utf-8",
        )
        return {
            "bindetect_results": str(path),
            "bindetect_outdir": str(tmp_path / "bd"),
            "replicate_groups": replicate_groups,
            "corrected_signals": {},
            "footprint_signals": {},
        }

    monkeypatch.setattr(mod, "_tobias_replicate_pipeline", fake_pipeline)
    monkeypatch.setattr(mod, "_aggregate_plots", lambda *args, **kwargs: {})

    base_params = {
        "replicate_bams": replicate_design,
        "replicate_batches": {
            "Treated": {"r1": "b1", "r2": "b1", "r3": "b2"},
            "Control": {"r1": "b1", "r2": "b1", "r3": "b2"},
        },
        "genome_fasta": str(tmp_path / "genome.fa"),
        "peaks_bed": str(tmp_path / "peaks.bed"),
        "motif_meme": str(tmp_path / "motifs.meme"),
        "group_a": "Treated",
        "group_b": "Control",
        "make_volcano_plot": False,
        "make_aggregate_plots": False,
    }
    unpaired = mod.chromatin_footprint_tobias_bulk(
        {**base_params, "output_dir": str(tmp_path / "unpaired")}
    )
    inf_u = unpaired["differential_summary"]["inference"]
    assert inf_u["pairing_status"] == "independent"
    assert inf_u["batch_adjusted"] is True
    assert "batch" in inf_u["test"]

    paired = mod.chromatin_footprint_tobias_bulk(
        {
            **base_params,
            "paired_replicates": True,
            "output_dir": str(tmp_path / "paired"),
        }
    )
    inf_p = paired["differential_summary"]["inference"]
    assert inf_p["pairing_status"] == "complete"
    assert inf_p["test"] == "paired_t_test_on_within_donor_mean_score_differences"
