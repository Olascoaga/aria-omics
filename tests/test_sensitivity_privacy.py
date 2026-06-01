"""P1-8 (W-PRIV): input sensitivity checkpoint (a) + failed-run redaction (c).

ARIA never auto-disables egress. It classifies input sensitivity, surfaces it at
the first checkpoint, and lets the user opt into air-gapped mode (which blocks
ALL network egress). Failed-run archives are redacted by default.
"""

from __future__ import annotations

import json

import pytest

from aria.utils import privacy
from aria.utils import sensitivity as sens


@pytest.fixture(autouse=True)
def _restore_air_gapped_env():
    """`enable_air_gapped_runtime()` mutates os.environ on purpose (so dispatched
    subprocesses inherit it). Snapshot + restore around every test so the flag
    never leaks into other tests in the session."""
    import os
    prev = os.environ.get("ARIA_AIR_GAPPED")
    prev_reason = privacy._runtime_air_gapped_reason
    yield
    if prev is None:
        os.environ.pop("ARIA_AIR_GAPPED", None)
    else:
        os.environ["ARIA_AIR_GAPPED"] = prev
    privacy._runtime_air_gapped_reason = prev_reason


# ── (a) sensitivity classifier ───────────────────────────────────────────────

def test_classify_clinical_phi_is_high_and_recommends_air_gap():
    a = sens.classify_sensitivity(
        organism="Homo sapiens",
        field_names=["patient_id", "diagnosis", "condition", "cell_type"],
        path_hints=["/data/clinical_cohort"],
    )
    assert a["level"] == "high"
    assert a["is_human"] is True
    assert a["recommend_air_gapped"] is True
    assert "patient" in a["phi_signals"]


def test_classify_human_without_phi_is_elevated_no_recommend():
    a = sens.classify_sensitivity(
        organism="Homo sapiens",
        field_names=["condition", "Donor", "cell_type", "leiden"],
        path_hints=["pbmc.h5ad"],
    )
    assert a["level"] == "elevated"
    assert a["is_human"] is True
    assert a["recommend_air_gapped"] is False
    assert a["phi_signals"] == []


def test_classify_nonhuman_is_low():
    a = sens.classify_sensitivity(
        organism="Mus musculus",
        field_names=["condition", "replicate", "cell_type"],
        path_hints=["mouse_brain"],
    )
    assert a["level"] == "low"
    assert a["is_human"] is False
    assert a["recommend_air_gapped"] is False


def test_classify_human_age_sex_quasi_identifiers_is_high():
    a = sens.classify_sensitivity(
        organism="human",
        field_names=["age", "sex", "condition"],
        path_hints=[],
    )
    assert a["level"] == "high"
    assert a["quasi_identifier_signals"] == ["age+sex"]
    assert a["recommend_air_gapped"] is True


# ── (a) checkpoint contract ──────────────────────────────────────────────────

def test_checkpoint_options_keep_continue_first_and_offer_air_gap():
    opts = sens.checkpoint_options(sens.classify_sensitivity(organism="human"))
    assert opts[0] == "Confirm and continue"          # default never changes
    assert sens.AIR_GAPPED_OPTION in opts
    assert opts[-1] == "Cancel"


def test_decision_enables_air_gapped_matches_option():
    assert sens.decision_enables_air_gapped(sens.AIR_GAPPED_OPTION) is True
    assert sens.decision_enables_air_gapped("Confirm and continue") is False


def test_annotate_question_flags_recommendation_when_sensitive():
    high = sens.classify_sensitivity(organism="human", field_names=["patient_id"])
    q = sens.annotate_checkpoint_question("Detected data ...", high)
    assert "Data sensitivity (HIGH)" in q
    assert "(RECOMMENDED)" in q


# ── (a) runtime air-gapped enable ────────────────────────────────────────────

def test_enable_air_gapped_runtime_blocks_egress(monkeypatch):
    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)
    assert privacy.egress_allowed() is True
    privacy.enable_air_gapped_runtime(reason="sensitivity_checkpoint")
    assert privacy.air_gapped_enabled() is True
    assert privacy.egress_allowed() is False
    assert privacy.air_gapped_runtime_reason() == "sensitivity_checkpoint"


# ── (c) failed-run output/error redaction ────────────────────────────────────

def _make_env(tmp_path, monkeypatch):
    monkeypatch.delenv("ARIA_PRESERVE_FAILED_INPUTS", raising=False)
    from aria.utils.environment_manager import EnvironmentManager
    em = EnvironmentManager.__new__(EnvironmentManager)
    em.workspace = tmp_path / "ws"
    em.workspace.mkdir(parents=True)
    return em


def test_failed_run_redacts_output_and_error(tmp_path, monkeypatch):
    em = _make_env(tmp_path, monkeypatch)
    input_file = tmp_path / "input.json"
    output_file = tmp_path / "output.json"
    input_file.write_text(json.dumps({"data_path": "/secret/patient/x.h5ad"}))
    output_file.write_text(json.dumps({"output_path": "/secret/out/de.csv",
                                       "status": "error"}))
    result = {"status": "error", "script_path": "/secret/aria/rna_bulk_de.py",
              "api_key": "sk-must-not-leak"}

    em._archive_failed_run("r1", "rna", input_file, output_file, result)

    run_dir = em.workspace / "failed" / "rna_r1"
    # raw artifacts must NOT exist; only redacted ones
    assert not (run_dir / "output.json").exists()
    assert not (run_dir / "error.json").exists()
    out = json.loads((run_dir / "output.redacted.json").read_text())
    err = json.loads((run_dir / "error.redacted.json").read_text())
    inp = json.loads((run_dir / "input.redacted.json").read_text())
    assert out["output_path"] == "<path:de.csv>"
    assert err["script_path"] == "<path:rna_bulk_de.py>"
    assert err["api_key"] == "<redacted>"
    assert inp["data_path"] == "<path:x.h5ad>"


def test_failed_run_preserve_opt_in_keeps_raw(tmp_path, monkeypatch):
    monkeypatch.setenv("ARIA_PRESERVE_FAILED_INPUTS", "1")
    from aria.utils.environment_manager import EnvironmentManager
    em = EnvironmentManager.__new__(EnvironmentManager)
    em.workspace = tmp_path / "ws"
    em.workspace.mkdir(parents=True)
    input_file = tmp_path / "input.json"
    output_file = tmp_path / "output.json"
    input_file.write_text(json.dumps({"data_path": "/secret/x.h5ad"}))
    output_file.write_text(json.dumps({"output_path": "/secret/out.csv"}))

    em._archive_failed_run("r2", "rna", input_file, output_file,
                           {"status": "error"})

    run_dir = em.workspace / "failed" / "rna_r2"
    assert (run_dir / "input.json").exists()       # raw preserved under opt-in
    assert (run_dir / "output.json").exists()
    assert not (run_dir / "output.redacted.json").exists()
