"""Plumbing guards for the c1_h9_fastq_e2e preprint-freeze lane.

These validate the lane registration and the harness's capsule-assembly /
artifact contract cheaply, without the multi-hour real FASTQ-to-report run.
The real receipt is regenerated once, against the final clean source snapshot,
during the freeze regeneration step (see memory/NEXT_SESSION.md).
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from aria.benchmarks.preprint_freeze import LANES


REPO_ROOT = Path(__file__).resolve().parents[1]
HARNESS_PATH = REPO_ROOT / "scripts" / "run_c1_h9_fastq_e2e.py"


def _load_harness():
    spec = importlib.util.spec_from_file_location(
        "run_c1_h9_fastq_e2e", HARNESS_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _lane():
    return next(item for item in LANES if item["lane_id"] == "c1_h9_fastq_e2e")


def test_lane_is_registered_and_executable():
    lane = _lane()
    assert lane["claims"] == ["claim_1"]
    assert lane["implementation"] == "available"
    assert lane["evidence_kind"] == "e2e_capsule"
    assert lane["command"] and "run_c1_h9_fastq_e2e.py" in lane["command"]
    assert lane["required_for_freeze"] is True


def test_lane_declares_the_five_emitted_artifacts():
    lane = _lane()
    assert tuple(lane["expected_artifacts"]) == (
        "claim_1/h9_e2e/capsule.json",
        "claim_1/h9_e2e/report.html",
        "claim_1/h9_e2e/methodology.json",
        "claim_1/h9_e2e/de_results.tsv",
        "claim_1/h9_e2e/fig1_h9_bulk_de.svg",
    )


def test_lane_binds_orchestrator_and_scientific_environments():
    lane = _lane()
    assert lane["resources"] == [
        "env:aria-env", "env:aria-rnaseq-env", "env:aria-rna-env", "data:h9_fastq"
    ]
    # The command drives the orchestrator env; EnvironmentManager dispatches the
    # heavy science into the scientific envs.
    assert "conda run -n aria-env" in lane["command"]


def test_h9_fastq_resource_is_probed():
    from aria.benchmarks.preprint_freeze import probe_resources
    probed = probe_resources(REPO_ROOT / "docs/benchmark_results/preprint_v1")
    assert "data:h9_fastq" in probed


def test_frozen_llm_double_satisfies_provider_seam():
    harness = _load_harness()
    double = harness.FrozenLLMDouble()
    scoped = double.for_execution("exp", "/tmp/usage.jsonl", egress_policy=None)
    # A run-scoped view must still answer completion calls deterministically and
    # without contacting any provider.
    assert scoped.complete("anything") == double.marker
    assert scoped.complete_heavy("x") == double.marker
    assert scoped.complete_light("x") == double.marker
    assert double.get_active_model() is None


def _make_fake_report(tmp_path: Path) -> Path:
    report_dir = tmp_path / "report_exp1234"
    (report_dir / "figures").mkdir(parents=True)
    (report_dir / "tables").mkdir(parents=True)
    (report_dir / "report.html").write_text("<html>report</html>", encoding="utf-8")
    (report_dir / "methodology.json").write_text(
        json.dumps({
            "contrasts": [
                {"name": "BMAL1_KO vs WT", "n_significant": 42,
                 "n_upregulated": 20, "n_downregulated": 22},
                {"name": "REV_ERBa_KO vs WT", "n_significant": 17,
                 "n_upregulated": 9, "n_downregulated": 8},
            ]
        }),
        encoding="utf-8",
    )
    (report_dir / "tables" / "deseq2_BMAL1_KO_vs_WT.tsv").write_text(
        "gene\tlog2fc\tpadj\nGENE1\t2.0\t0.001\n", encoding="utf-8"
    )
    (report_dir / "figures" / "volcano_BMAL1_KO_vs_WT.svg").write_text(
        "<svg/>", encoding="utf-8"
    )
    return report_dir


def test_publish_assembles_capsule_and_copies_canonical_artifacts(tmp_path):
    harness = _load_harness()
    report_dir = _make_fake_report(tmp_path)
    output_dir = tmp_path / "out"

    capsule = harness.publish(
        REPO_ROOT, report_dir, output_dir,
        experiment_id="exp1234", decisions=[(1, "Confirm"), (2.3, "condition")],
    )

    for name in ("capsule.json", "report.html", "methodology.json",
                 "de_results.tsv", "fig1_h9_bulk_de.svg"):
        assert (output_dir / name).is_file(), name

    assert capsule["schema_version"] == harness.CAPSULE_SCHEMA
    assert capsule["lane_id"] == "c1_h9_fastq_e2e"
    assert {c["name"] for c in capsule["canonical_artifacts"]} == {
        "report.html", "methodology.json", "de_results.tsv", "fig1_h9_bulk_de.svg"
    }
    # DE summary is extracted schema-agnostically from methodology.json.
    names = {row.get("name") for row in capsule["de_summary"]}
    assert "BMAL1_KO vs WT" in names and "REV_ERBa_KO vs WT" in names
    # Every environment lock is enumerated by name (sha256 present when the lock
    # is committed).
    assert [e["env_name"] for e in capsule["environments"]] == [
        "aria-env", "aria-rnaseq-env", "aria-rna-env"
    ]


def test_publish_fails_loudly_without_a_de_table(tmp_path):
    harness = _load_harness()
    report_dir = _make_fake_report(tmp_path)
    # Remove the only DE table: a completed report without DE evidence must not
    # yield a partial capsule.
    (report_dir / "tables" / "deseq2_BMAL1_KO_vs_WT.tsv").unlink()

    with pytest.raises(RuntimeError, match="differential-expression table"):
        harness.publish(
            REPO_ROOT, report_dir, tmp_path / "out",
            experiment_id="exp1234", decisions=[],
        )


def test_capsule_carries_no_absolute_machine_paths(tmp_path):
    harness = _load_harness()
    report_dir = _make_fake_report(tmp_path)
    output_dir = tmp_path / "out"
    harness.publish(
        REPO_ROOT, report_dir, output_dir,
        experiment_id="exp1234", decisions=[],
    )
    text = (output_dir / "capsule.json").read_text(encoding="utf-8")
    assert "/home/" not in text
    assert str(tmp_path) not in text
