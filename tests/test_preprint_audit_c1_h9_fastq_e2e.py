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
    parsed = json.loads(scoped.complete(
        'Return JSON containing "analysis_type" and "key_modalities_needed".'
    ))
    assert parsed["analysis_type"] == "differential"
    assert parsed["key_modalities_needed"] == ["RNA"]
    plan = json.loads(scoped.complete(
        'Return JSON containing "steps" and "contrasts".'
    ))
    assert plan["contrasts"] == [
        {"numerator": "B", "denominator": "WT"},
        {"numerator": "R", "denominator": "WT"},
    ]
    assert scoped.complete_heavy("x") == double.marker
    assert scoped.complete_light("x") == double.marker
    assert double.get_active_model() is None


def test_frozen_h9_policy_drives_real_design_checkpoints(tmp_path):
    harness = _load_harness()
    from aria.agents.design_agent import DesignAgent
    from aria.bus.message_bus import bus
    from aria.memory.memory import ARIAMemory

    fastqs = []
    for samples in harness.H9_GROUPS.values():
        for sample in samples:
            for read in (1, 2):
                path = tmp_path / f"{sample}_{read}.fq.gz"
                path.touch()
                fastqs.append(str(path))

    experiment_id = "h9_frozen_design"
    memory = ARIAMemory(db_path=str(tmp_path / "memory.db"))
    agent = DesignAgent(memory=memory, llm=harness.FrozenLLMDouble())
    result = agent.start_design(
        experiment_id,
        {
            "modalities": {"bulk_RNA_raw": fastqs},
            "organism": "Homo sapiens",
            "genome": "hg38",
        },
        {"question": harness.DEFAULT_QUESTION},
    )
    decisions = []
    for _ in range(20):
        pending = [
            msg for msg in bus.get_pending_checkpoints(
                experiment_id=experiment_id
            )
            if not msg.payload.get("resolved")
        ]
        if not pending:
            break
        msg = pending[0]
        cp_num = msg.payload.get("checkpoint")
        choice = harness.frozen_h9_answer_policy(
            cp_num,
            msg.payload.get("question", ""),
            msg.payload.get("options", []),
        )
        bus.resolve_checkpoint(msg.id, {"choice": choice})
        decisions.append((cp_num, choice))
        result = agent.handle_user_response(
            experiment_id=experiment_id,
            checkpoint_num=cp_num,
            choice=choice,
        )
        if result.get("status") in {"done", "cancelled"}:
            break

    memory.close()
    assert result["status"] == "done", result
    design = result["design"]
    assert design["groups"] == {
        group: list(samples) for group, samples in harness.H9_GROUPS.items()
    }
    assert design["replicates"] == {"B": 3, "R": 3, "WT": 3}
    assert design["main_factor"] == "condition"
    assert design["design_formula"] == "~ condition"
    assert design["design_status"] == "ready"
    assert dict(decisions)[2.1] == harness.H9_MANUAL_GROUP_ASSIGNMENT
    assert dict(decisions)[2.6] == "Yes — proceed"


def test_frozen_h9_design_maps_to_bulk_matrix_and_explicit_contrasts(tmp_path):
    harness = _load_harness()
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    counts = tmp_path / "counts.tsv"
    sample_names = [
        sample for samples in harness.H9_GROUPS.values() for sample in samples
    ]
    counts.write_text(
        "gene\t" + "\t".join(sample_names) + "\n"
        "GENE_1\t" + "\t".join(["10"] * len(sample_names)) + "\n",
        encoding="utf-8",
    )
    design = {
        "groups": {
            group: list(samples) for group, samples in harness.H9_GROUPS.items()
        },
        "main_factor": "condition",
        "plan_contrasts": list(harness.H9_CONTRASTS),
    }

    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    mapped_samples, labels, factor, contrasts = agent._apply_design(
        design, [str(counts)], "h9_design_mapping"
    )

    assert mapped_samples == sample_names
    assert labels == {
        sample: group
        for group, samples in harness.H9_GROUPS.items()
        for sample in samples
    }
    assert factor == "condition"
    assert contrasts == [
        {"numerator": "B", "denominator": "WT", "name": "B vs WT"},
        {"numerator": "R", "denominator": "WT", "name": "R vs WT"},
    ]


def _make_fake_report(tmp_path: Path) -> Path:
    report_dir = tmp_path / "report_exp1234"
    (report_dir / "figures" / "b_vs_wt").mkdir(parents=True)
    (report_dir / "tables").mkdir(parents=True)
    (report_dir / "report.html").write_text(
        "<html><a href='/home/test/run/b_vs_wt_de_genes.tsv'>report</a></html>",
        encoding="utf-8",
    )
    (report_dir / "methodology.json").write_text(
        json.dumps({
            "contrasts": [
                {"name": "BMAL1_KO vs WT", "n_significant": 42,
                 "n_upregulated": 20, "n_downregulated": 22},
                {"name": "REV_ERBa_KO vs WT", "n_significant": 17,
                 "n_upregulated": 9, "n_downregulated": 8},
            ],
            "source_table": "/home/test/run/b_vs_wt_de_genes.tsv",
        }),
        encoding="utf-8",
    )
    (report_dir / "tables" / "b_vs_wt_de_genes.tsv").write_text(
        "gene\tlog2fc\tpadj\nGENE1\t2.0\t0.001\n", encoding="utf-8"
    )
    (report_dir / "figures" / "b_vs_wt" / "volcano.svg").write_text(
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
    assert capsule["requested_cpus"] == 30
    assert capsule["checkpoint_policy"].endswith("frozen_h9_answer_policy")
    assert capsule["design_contract"] == {
        "groups": {"B": ["B1", "B2", "B3"],
                   "R": ["R1", "R2", "R3"],
                   "WT": ["WT1", "WT2", "WT3"]},
        "contrasts": [
            {"numerator": "B", "denominator": "WT"},
            {"numerator": "R", "denominator": "WT"},
        ],
        "factor": "condition",
    }
    assert {c["name"] for c in capsule["canonical_artifacts"]} == {
        "report.html", "methodology.json", "de_results.tsv", "fig1_h9_bulk_de.svg"
    }
    # DE summary is extracted schema-agnostically from methodology.json.
    names = {row.get("name") for row in capsule["de_summary"]}
    assert "BMAL1_KO vs WT" in names and "REV_ERBa_KO vs WT" in names
    assert "/home/" not in (output_dir / "methodology.json").read_text()
    assert "/home/" not in (output_dir / "report.html").read_text()
    assert json.loads((output_dir / "methodology.json").read_text())[
        "source_table"
    ] == "b_vs_wt_de_genes.tsv"
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
    (report_dir / "tables" / "b_vs_wt_de_genes.tsv").unlink()

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
