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
        lambda *args: (replicate_info, replicate_groups,
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
