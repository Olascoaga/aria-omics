"""P3a — scATAC external-concordance harness (ARIA vs SnapATAC2/ArchR/Signac).

These guard the PURE scoring layer (cluster ARI/NMI, DA-peak genomic overlap,
motif concordance, seed stability) and the honest-gating contract of the
dispatchable runner: with no external-tool outputs the harness reports an honest
``not_run``/``MissingDependency`` instead of a fabricated concordance of 1.0.

Neutral synthetic labels only (ADR-011): integer clusters, ``PEAK_*`` /
coordinate peaks, ``MOTIF_*``. This is the buildable-now P3a scaffold; the actual
SnapATAC2 driver in ``aria-bench-env`` is the infra-gated P3b drop-in.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# --- cluster concordance: ARI / NMI -----------------------------------------

def test_ari_is_one_for_identical_partition_under_relabeling():
    from aria.benchmarks.concordance_atac import adjusted_rand_index

    a = [0, 0, 1, 1, 2, 2]
    b = ["x", "x", "y", "y", "z", "z"]  # same partition, different label names
    assert adjusted_rand_index(a, b) == pytest.approx(1.0)


def test_ari_near_zero_for_independent_partitions():
    from aria.benchmarks.concordance_atac import adjusted_rand_index

    a = [0, 1, 0, 1, 0, 1, 0, 1]
    b = [0, 0, 0, 0, 1, 1, 1, 1]
    assert adjusted_rand_index(a, b) == pytest.approx(0.0, abs=0.2)


def test_nmi_one_for_identical_and_zero_for_no_information():
    from aria.benchmarks.concordance_atac import normalized_mutual_info

    a = [0, 0, 1, 1, 2, 2]
    assert normalized_mutual_info(a, a) == pytest.approx(1.0)
    # one labeling carries no cluster structure -> no shared information
    constant = [0, 0, 0, 0, 0, 0]
    assert normalized_mutual_info(a, constant) == pytest.approx(0.0)


# --- DA-peak overlap (exact + genomic) --------------------------------------

def test_peak_overlap_exact_jaccard():
    from aria.benchmarks.concordance_atac import peak_set_overlap

    a = ["PEAK_1", "PEAK_2", "PEAK_3", "PEAK_4"]
    b = ["PEAK_3", "PEAK_4", "PEAK_5"]
    res = peak_set_overlap(a, b, mode="exact")
    # intersection 2, union 5
    assert res["jaccard"] == pytest.approx(2 / 5)
    assert res["n_a"] == 4 and res["n_b"] == 3


def test_peak_overlap_genomic_matches_boundary_shifted_peaks():
    from aria.benchmarks.concordance_atac import peak_set_overlap

    # boundary-shifted but overlapping intervals are matched in genomic mode,
    # which exact-string overlap would miss (the ADR-043 lesson).
    a = ["chr1:100-300", "chr1:1000-1200", "chr2:50-150"]
    b = ["chr1:200-400", "chr2:60-160"]  # overlap a[0] and a[2]; a[1] unmatched
    exact = peak_set_overlap(a, b, mode="exact")
    genomic = peak_set_overlap(a, b, mode="genomic")
    assert exact["jaccard"] == pytest.approx(0.0)
    assert genomic["matched_a"] == 2
    assert genomic["matched_b"] == 2
    assert genomic["overlap_fraction"] > exact["jaccard"]


# --- motif concordance ------------------------------------------------------

def test_motif_concordance_top_k_overlap_and_rank():
    from aria.benchmarks.concordance_atac import motif_concordance

    a = ["MOTIF_A", "MOTIF_B", "MOTIF_C", "MOTIF_D"]
    b = ["MOTIF_A", "MOTIF_B", "MOTIF_E"]
    res = motif_concordance(a, b, top_k=3)
    # top-3 of each: {A,B,C} vs {A,B,E} -> intersection 2, union 4
    assert res["top_k"] == 3
    assert res["jaccard_top_k"] == pytest.approx(2 / 4)
    assert 0.0 <= res["rank_biased_overlap"] <= 1.0


# --- seed stability ---------------------------------------------------------

def test_seed_stability_is_one_for_identical_label_sets():
    from aria.benchmarks.concordance_atac import seed_stability

    labels = [[0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2], [2, 2, 0, 0, 1, 1]]
    res = seed_stability(labels)
    assert res["mean_pairwise_ari"] == pytest.approx(1.0)
    assert res["n_seeds"] == 3


# --- aggregator: honest gating ----------------------------------------------

def test_score_concordance_is_not_run_without_external_outputs():
    from aria.benchmarks.concordance_atac import score_atac_concordance

    aria = {"clusters": [0, 0, 1, 1], "da_peaks": ["chr1:1-2"], "motifs": ["MOTIF_A"]}
    res = score_atac_concordance(aria, external=None, tool="SnapATAC2")
    assert res["status"] == "not_run"
    assert "reason" in res
    # never fabricate a perfect score when the comparator did not run
    assert "cluster_concordance" not in res or res.get("cluster_concordance") is None


def test_score_concordance_intrinsic_only_with_seeds_no_external():
    from aria.benchmarks.concordance_atac import score_atac_concordance

    aria = {
        "clusters": [0, 0, 1, 1, 2, 2],
        "cluster_seeds": [[0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2], [2, 2, 0, 0, 1, 1]],
        "da_peaks": ["chr1:1-2"],
        "motifs": ["MOTIF_A"],
    }
    res = score_atac_concordance(aria, external=None, tool="SnapATAC2")
    # ARIA-intrinsic seed stability does not need a comparator; cross-tool axes do.
    assert res["status"] == "intrinsic_only"
    assert res["seed_stability"]["mean_pairwise_ari"] == pytest.approx(1.0)
    assert res.get("cluster_concordance") is None
    assert res.get("da_peak_concordance") is None
    assert "reason" in res


def test_score_concordance_success_with_paired_outputs():
    from aria.benchmarks.concordance_atac import score_atac_concordance

    aria = {
        "clusters": [0, 0, 1, 1, 2, 2],
        "da_peaks": ["chr1:100-300", "chr1:1000-1200"],
        "motifs": ["MOTIF_A", "MOTIF_B", "MOTIF_C"],
    }
    external = {
        "clusters": [0, 0, 1, 1, 2, 2],
        "da_peaks": ["chr1:150-350", "chr2:5-9"],
        "motifs": ["MOTIF_A", "MOTIF_B", "MOTIF_D"],
    }
    res = score_atac_concordance(aria, external, tool="SnapATAC2")
    assert res["status"] == "success"
    assert res["tool"] == "SnapATAC2"
    assert res["cluster_concordance"]["ari"] == pytest.approx(1.0)
    assert res["da_peak_concordance"]["matched_a"] == 1
    assert 0.0 <= res["motif_concordance"]["jaccard_top_k"] <= 1.0


# --- dispatchable runner ----------------------------------------------------

def test_concordance_contract_registered_as_beta():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    key = "aria/scripts/benchmark_scatac_concordance.py"
    assert key in SCRIPT_CONTRACTS
    assert SCRIPT_CONTRACTS[key].validation_level == "beta"


def test_runner_honest_not_run_when_external_outputs_missing(tmp_path):
    from aria.scripts.benchmark_scatac_concordance import benchmark_scatac_concordance

    aria_json = tmp_path / "aria.json"
    aria_json.write_text(json.dumps({
        "clusters": [0, 0, 1, 1],
        "da_peaks": ["chr1:1-2"],
        "motifs": ["MOTIF_A"],
    }), encoding="utf-8")

    res = benchmark_scatac_concordance({
        "aria_outputs_json": str(aria_json),
        "external_outputs_json": str(tmp_path / "does_not_exist.json"),
        "tool": "SnapATAC2",
        "output_dir": str(tmp_path),
    })
    assert res["status"] == "not_run"
    assert res["error_type"] == "MissingDependency"


def test_runner_intrinsic_only_when_seeds_present_but_no_external(tmp_path):
    from aria.scripts.benchmark_scatac_concordance import benchmark_scatac_concordance

    aria_json = tmp_path / "aria.json"
    aria_json.write_text(json.dumps({
        "clusters": [0, 0, 1, 1, 2, 2],
        "cluster_seeds": [[0, 0, 1, 1, 2, 2], [2, 2, 0, 0, 1, 1]],
        "da_peaks": ["chr1:1-2"],
        "motifs": ["MOTIF_A"],
    }), encoding="utf-8")

    res = benchmark_scatac_concordance({
        "aria_outputs_json": str(aria_json),
        "external_outputs_json": str(tmp_path / "absent.json"),
        "tool": "SnapATAC2",
        "output_dir": str(tmp_path),
    })
    assert res["status"] == "intrinsic_only"
    assert res["seed_stability"]["mean_pairwise_ari"] == pytest.approx(1.0)
    assert res["cluster_concordance"] is None
    assert res["inputs"]["external_tool_present"] is False
    assert res["inputs"]["n_cluster_seeds"] == 2
    assert Path(tmp_path / Path(res["artifacts"]["manifest_json"]).name).exists()


def test_runner_scores_paired_outputs_and_writes_manifest(tmp_path):
    from aria.scripts.benchmark_scatac_concordance import benchmark_scatac_concordance

    aria_json = tmp_path / "aria.json"
    ext_json = tmp_path / "ext.json"
    aria_json.write_text(json.dumps({
        "clusters": [0, 0, 1, 1, 2, 2],
        "da_peaks": ["chr1:100-300", "chr1:1000-1200"],
        "motifs": ["MOTIF_A", "MOTIF_B"],
    }), encoding="utf-8")
    ext_json.write_text(json.dumps({
        "clusters": [0, 0, 1, 1, 2, 2],
        "da_peaks": ["chr1:150-350"],
        "motifs": ["MOTIF_A", "MOTIF_C"],
    }), encoding="utf-8")

    res = benchmark_scatac_concordance({
        "aria_outputs_json": str(aria_json),
        "external_outputs_json": str(ext_json),
        "tool": "SnapATAC2",
        "output_dir": str(tmp_path),
    })
    assert res["status"] == "success"
    manifest_json = res["artifacts"]["manifest_json"]
    # T10 hygiene: no machine-absolute path in the public manifest.
    assert not manifest_json.startswith("/")
    assert (tmp_path / Path(manifest_json).name).exists()
    assert res["cluster_concordance"]["ari"] == pytest.approx(1.0)


# --- P3b reference driver (infra-gated SnapATAC2 producer) -------------------

def test_reference_driver_contract_registered_as_scaffold():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    key = "aria/scripts/benchmark_scatac_reference_snapatac2.py"
    assert key in SCRIPT_CONTRACTS
    assert SCRIPT_CONTRACTS[key].validation_level == "scaffold"


def test_reference_driver_honest_not_run_without_input():
    from aria.scripts.benchmark_scatac_reference_snapatac2 import (
        benchmark_scatac_reference_snapatac2,
    )

    res = benchmark_scatac_reference_snapatac2({"input_path": "/nope/missing.h5ad"})
    assert res["status"] == "not_run"
    assert res["error_type"] == "MissingDependency"


def test_reference_driver_honest_not_run_when_snapatac2_absent(tmp_path):
    # aria-env has no SnapATAC2; the driver must refuse, not fabricate outputs.
    pytest.importorskip  # keep import side-effect-free
    try:
        import snapatac2  # noqa: F401
        pytest.skip("SnapATAC2 present; honest-skip path not exercised here")
    except ImportError:
        pass

    from aria.scripts.benchmark_scatac_reference_snapatac2 import (
        benchmark_scatac_reference_snapatac2,
    )

    h5ad = tmp_path / "matrix.h5ad"
    h5ad.write_bytes(b"not-a-real-h5ad")  # existence is enough to pass the path gate
    res = benchmark_scatac_reference_snapatac2({"input_path": str(h5ad)})
    assert res["status"] == "not_run"
    assert res["error_type"] == "MissingDependency"
    assert "aria-bench-atac-env" in res["reason"]


def test_ranked_motifs_from_enrichment_orders_and_dedups():
    # Pure helper: rank most-enriched first (smallest adjusted p, then largest
    # log2FC), upper-case TF names, keep each motif's strongest cross-group signal.
    from aria.scripts.benchmark_scatac_reference_snapatac2 import (
        _ranked_motifs_from_enrichment,
    )

    class _FakeDF:
        def __init__(self, rows):
            self._rows = rows
            self.columns = ["name", "adjusted p-value", "log2(fold change)"]

        def iter_rows(self, named=True):  # noqa: ARG002 - mirror polars API
            return iter(self._rows)

    enrichment = {
        "0": _FakeDF([
            {"name": "Klf5_1.00", "adjusted p-value": 0.10, "log2(fold change)": 1.0},
            {"name": "TEAD1", "adjusted p-value": 0.50, "log2(fold change)": 0.2},
        ]),
        "1": _FakeDF([
            # stronger signal for the same KLF5 motif -> must win the dedup
            {"name": "KLF5", "adjusted p-value": 0.001, "log2(fold change)": 3.0},
            {"name": "SP1", "adjusted p-value": 0.02, "log2(fold change)": 2.0},
        ]),
    }
    ranked = _ranked_motifs_from_enrichment(enrichment)
    assert ranked == ["KLF5", "SP1", "TEAD1"]  # by best adjusted p-value
    assert ranked.count("KLF5") == 1  # deduped across groups
