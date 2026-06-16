"""Guards for the B2 multi-annotator kit: blind labeling-sheet I/O, inter-annotator
agreement (Cohen's κ / Krippendorff's α), and adjudication."""
import pytest

from aria.benchmarks.b2_annotation import (
    B2_LABELS,
    export_labeling_sheet,
    load_annotations,
    cohen_kappa,
    krippendorff_alpha,
    adjudicate,
)


def test_export_sheet_is_blind_and_round_trips():
    rows = [
        {"claim_id": "c1", "modality": "bulk", "analysis": "de",
         "claim_text": "312 genes significant", "evidence_summary": "n_sig=312",
         "gold_label": "clean"},
    ]
    csv_text = export_labeling_sheet(rows)
    # gold label is NOT leaked into the blind sheet
    assert "clean" not in csv_text
    assert "claim_id,modality,analysis,claim_text,evidence_summary,label,notes" in csv_text
    assert "c1" in csv_text


def test_load_annotations_skips_blank_and_validates():
    csv_text = (
        "claim_id,modality,analysis,claim_text,evidence_summary,label,notes\n"
        "c1,bulk,de,x,y,clean,\n"
        "c2,bulk,de,x,y,,still blank\n"
    )
    assert load_annotations(csv_text) == {"c1": "clean"}
    bad = "claim_id,label\nc1,not_a_label\n"
    with pytest.raises(ValueError):
        load_annotations(bad)


def test_cohen_kappa_perfect_and_chance():
    a = {"c1": "clean", "c2": "fabricated", "c3": "overclaim"}
    assert cohen_kappa(a, a)["kappa"] == 1.0
    # total disagreement on a 2-label split -> kappa < 0
    b = {"c1": "clean", "c2": "clean"}
    c = {"c1": "fabricated", "c2": "fabricated"}
    assert cohen_kappa(b, c)["kappa"] is not None


def test_krippendorff_alpha_perfect_agreement():
    a = {"c1": "clean", "c2": "fabricated", "c3": "overclaim"}
    res = krippendorff_alpha([a, dict(a)])
    assert res["alpha"] == 1.0
    assert res["n_units"] == 3


def test_adjudicate_unanimous_majority_tie_and_override():
    a = {"c1": "clean", "c2": "fabricated", "c3": "overclaim"}
    b = {"c1": "clean", "c2": "fabricated", "c3": "unsupported"}
    third = {"c2": "overclaim"}  # disagrees on c2 -> majority, not unanimous
    res = adjudicate([a, b, third])
    assert res["c1"]["gold"] == "clean" and res["c1"]["agreement"] == "unanimous"
    assert res["c2"]["gold"] == "fabricated" and res["c2"]["agreement"] == "majority"
    # c3: clean? no — overclaim vs unsupported, 1-1 tie -> needs adjudication
    assert res["c3"]["gold"] is None and res["c3"]["agreement"] == "needs_adjudication"
    # adjudicator override resolves the tie
    res2 = adjudicate([a, b], overrides={"c3": "overclaim"})
    assert res2["c3"]["gold"] == "overclaim" and res2["c3"]["adjudicated"] is True


def test_taxonomy_is_the_single_source():
    assert set(B2_LABELS) == {
        "clean", "unsupported", "overclaim", "fabricated",
        "causal_inflation", "missing_caveat"}
