"""P1-4: hidden (unmodeled) batch detector — WARN only, never corrects.

`assess_hidden_batch` inspects obs column NAMES and the user-confirmed design
(condition / replicate / declared batch) plus whether integration ran, and
flags candidate technical/batch columns that are present in the data but were
neither declared, corrected, nor modeled. It is name+design based (generic
technical tokens, ADR-011 exception like sensitivity.py) and never reads
cell-level values nor applies any correction.
"""

from aria.utils.batch_qc import assess_hidden_batch


def test_clean_when_no_candidate_batch_columns():
    res = assess_hidden_batch(
        obs_columns=["cell_type", "total_counts", "pct_counts_mt"],
        condition_col="condition",
        replicate_col="donor_id",
        declared_batch=None,
        integration_ran=False,
    )
    assert res["status"] == "clean"
    assert res["issues"] == []
    assert res["candidate_batch_columns"] == []


def test_flags_unmodeled_technical_column():
    # 'sequencing_lane' is technical, undeclared, uncorrected -> warn.
    res = assess_hidden_batch(
        obs_columns=["cell_type", "condition", "donor_id", "sequencing_lane"],
        condition_col="condition",
        replicate_col="donor_id",
        declared_batch=None,
        integration_ran=False,
    )
    checks = {i["check"] for i in res["issues"]}
    assert "unmodeled_batch" in checks
    assert "sequencing_lane" in res["candidate_batch_columns"]
    # Pure detector: it must recommend, not claim it corrected anything.
    msg = next(i for i in res["issues"] if i["check"] == "unmodeled_batch")
    assert "correct" not in msg["message"].lower() or "not" in msg["message"].lower()


def test_declared_batch_is_not_flagged():
    res = assess_hidden_batch(
        obs_columns=["cell_type", "condition", "batch"],
        condition_col="condition",
        replicate_col=None,
        declared_batch="batch",
        integration_ran=True,
    )
    assert res["status"] == "clean"
    assert res["candidate_batch_columns"] == []


def test_replicate_and_condition_columns_are_not_treated_as_hidden_batch():
    # donor_id is the replicate (modeled at pseudobulk); condition is the factor.
    res = assess_hidden_batch(
        obs_columns=["condition", "donor_id"],
        condition_col="condition",
        replicate_col="donor_id",
        declared_batch=None,
        integration_ran=False,
    )
    assert res["status"] == "clean"
    assert res["issues"] == []


def test_candidate_column_corrected_by_integration_is_downgraded():
    # 'batch' present and integration ran on it -> not a hidden/unmodeled batch.
    res = assess_hidden_batch(
        obs_columns=["condition", "batch"],
        condition_col="condition",
        replicate_col=None,
        declared_batch="batch",
        integration_ran=True,
    )
    assert all(i["check"] != "unmodeled_batch" for i in res["issues"])


def test_graceful_on_empty_inputs():
    res = assess_hidden_batch(
        obs_columns=None,
        condition_col=None,
        replicate_col=None,
        declared_batch=None,
        integration_ran=False,
    )
    assert res["status"] == "clean"
    assert res["issues"] == []
