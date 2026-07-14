"""Plumbing guards for the c5_multimodal_null_permutations preprint-freeze lane.

These validate the lane registration and the null-permutation logic (label
permutation that breaks the true grouping, real-compiler classification, per
modality scoring) WITHOUT the aria-rna-env / aria-chromatin-env subprocesses:
fixture DESeq2 results stand in for the real runs. The lane is code-only (no
human, no external data resource) and reaches a real receipt at freeze time.
"""

from __future__ import annotations

from pathlib import Path

from aria.benchmarks import multimodal_null as mn
from aria.benchmarks.preprint_freeze import LANES


def _de(n_sig: int, n_up: int = 0, n_down: int = 0):
    genes = [f"F{i:03d}" for i in range(n_sig)]
    return {
        "status": "success", "n_tested": 120, "n_sig": n_sig,
        "n_up": n_up, "n_down": n_down, "sig_peaks": genes,
        "lfc_by_peak": {g: 2.0 for g in genes},
        "padj_by_peak": {g: 0.01 for g in genes},
    }


def _lane():
    return next(x for x in LANES if x["lane_id"] == "c5_multimodal_null_permutations")


def test_lane_is_registered_and_executable():
    lane = _lane()
    assert lane["claims"] == ["claim_5"]
    assert lane["implementation"] == "available"
    assert lane["evidence_kind"] == "multimodal_e2e"
    assert "run_c5_multimodal_null_permutations.py" in lane["command"]


def test_lane_binds_all_three_environments_and_has_no_human_gate():
    lane = _lane()
    assert lane["resources"] == [
        "env:aria-env", "env:aria-rna-env", "env:aria-chromatin-env"
    ]
    assert tuple(lane["expected_artifacts"]) == (
        "claim_5/multimodal_null/multimodal_null.json",
    )
    # Code-only: no human or external-data resource gates the receipt.
    assert not any(r.startswith("human:") or r.startswith("data:")
                   for r in lane["resources"])


def test_referenced_runner_scripts_exist():
    root = Path(__file__).resolve().parents[1]
    assert (root / "scripts" / "run_c5_multimodal_null_permutations.py").is_file()
    assert (root / "scripts" / "aria_atac_pseudobulk_matrix.py").is_file()
    assert (root / "scripts" / "aria_pseudobulk_da_from_tsv.py").is_file()


def test_permutation_breaks_the_true_grouping():
    import pandas as pd

    meta = pd.DataFrame({
        "sample": [f"A{i}" for i in range(6)] + [f"B{i}" for i in range(6)],
        "condition": ["COND_A"] * 6 + ["COND_B"] * 6,
    })
    true_group = mn._grouping(list(meta["sample"]), list(meta["condition"]))
    complement = frozenset(meta["sample"]) - true_group
    for seed in range(25):
        perm = mn.permute_conditions(meta, seed=seed)
        # Label counts are preserved (a permutation, not a resample).
        assert list(perm["condition"]).count("COND_B") == 6
        group = mn._grouping(list(perm["sample"]), list(perm["condition"]))
        # Never the observed grouping nor its complement (those retain signal).
        assert group != true_group and group != complement


def test_classify_run_uses_the_real_compiler_per_modality():
    for modality in ("rna", "atac"):
        signal = mn.classify_run(modality, _de(12, 6, 6), is_permuted=False)
        null = mn.classify_run(modality, _de(0), is_permuted=True)
        # Positive control: a real signal produces an emitted significant claim.
        assert signal["emitted"] and signal["asserts_significant"]
        assert signal["false_positive"] is False
        # Null: an honest "no features" claim is not a false-positive narrative.
        assert null["asserts_significant"] is False
        assert null["false_positive"] is False


def test_score_flags_a_null_false_positive():
    true_run = mn.classify_run("rna", _de(10, 5, 5), is_permuted=False)
    clean = [mn.classify_run("rna", _de(0), is_permuted=True) for _ in range(10)]
    s_clean = mn.score_modality("rna", true_run, clean)
    assert s_clean["false_positive_narrative_rate"] == 0.0
    assert all(s_clean["axis_pass"].values())

    # A permutation that emits a significant claim is a false-positive narrative;
    # a rate above the FDR tolerance fails the bounded-null axis.
    leaky = [mn.classify_run("rna", _de(0), is_permuted=True) for _ in range(3)]
    leaky.append(mn.classify_run("rna", _de(9, 4, 5), is_permuted=True))
    s_leaky = mn.score_modality("rna", true_run, leaky)
    assert s_leaky["n_false_positive_narratives"] == 1
    assert s_leaky["false_positive_narrative_rate"] == 0.25
    assert s_leaky["axis_pass"]["null_false_positive_bounded"] is False


def test_multimodal_manifest_passes_only_when_every_axis_passes():
    good = {
        "rna": {"axis_pass": {"positive_control": True,
                              "null_false_positive_bounded": True}},
        "atac": {"axis_pass": {"positive_control": True,
                               "null_false_positive_bounded": True}},
    }
    assert mn.score_multimodal_null(good)["status"] == "pass"
    bad = {"rna": good["rna"],
           "atac": {"axis_pass": {"positive_control": True,
                                  "null_false_positive_bounded": False}}}
    assert mn.score_multimodal_null(bad)["status"] == "fail"
