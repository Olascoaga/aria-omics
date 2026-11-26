"""Plumbing guards for the c1_public_fastq_e2e (public paired GSE52778) preprint-freeze lane."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from aria.benchmarks.preprint_freeze import LANES, probe_resources


REPO_ROOT = Path(__file__).resolve().parents[1]
HARNESS_PATH = REPO_ROOT / "scripts" / "run_c1_public_fastq_e2e.py"


def _load_harness():
    spec = importlib.util.spec_from_file_location(
        "run_c1_public_fastq_e2e", HARNESS_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _lane():
    return next(item for item in LANES if item["lane_id"] == "c1_public_fastq_e2e")


def test_lane_is_registered_and_executable():
    lane = _lane()
    assert lane["claims"] == ["claim_1"]
    assert lane["implementation"] == "available"
    assert lane["evidence_kind"] == "e2e_capsule"
    assert lane["command"] and "run_c1_public_fastq_e2e.py" in lane["command"]
    assert lane["required_for_freeze"] is True


def test_lane_declares_the_five_emitted_artifacts():
    lane = _lane()
    assert tuple(lane["expected_artifacts"]) == (
        "claim_1/public_fastq_e2e/capsule.json",
        "claim_1/public_fastq_e2e/report.html",
        "claim_1/public_fastq_e2e/methodology.json",
        "claim_1/public_fastq_e2e/de_results.tsv",
        "claim_1/public_fastq_e2e/fig1_gse52778_bulk_de.svg",
    )


def test_lane_binds_orchestrator_and_scientific_environments():
    lane = _lane()
    assert lane["resources"] == [
        "env:aria-env", "env:aria-rnaseq-env", "env:aria-rna-env", "data:gse52778_fastq"
    ]
    assert "conda run -n aria-env" in lane["command"]


def test_public_gse52778_resource_is_probed_without_private_paths():
    probed = probe_resources(REPO_ROOT / "docs/benchmark_results/preprint_v1")
    assert "data:gse52778_fastq" in probed
    assert "Samael" not in json.dumps(probed)
    assert "Samael" not in json.dumps(LANES)


def test_frozen_llm_double_satisfies_provider_seam():
    harness = _load_harness()
    double = harness.FrozenLLMDouble()
    scoped = double.for_execution("exp", "/tmp/usage.jsonl", egress_policy=None)
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
        {"numerator": "dexamethasone", "denominator": "untreated"},
    ]
    assert double.get_active_model() is None


def test_frozen_gse52778_policy_selects_paired_donor_block_at_cp25(tmp_path):
    harness = _load_harness()
    from aria.agents.design_agent import DesignAgent
    from aria.bus.message_bus import bus
    from aria.memory.memory import ARIAMemory

    fastqs = []
    for samples in harness.GSE52778_GROUPS.values():
        for sample in samples:
            for read in (1, 2):
                path = tmp_path / f"{sample}_R{read}.fastq.gz"
                path.touch()
                fastqs.append(str(path))

    experiment_id = "gse52778_frozen_design"
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
        choice = harness.frozen_gse52778_answer_policy(
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
    assert design["replicates"] == {"dexamethasone": 4, "untreated": 4}
    assert design["main_factor"] == "condition"
    assert design["batch_covariate"] == "donor"
    assert design["design_formula"] == "~ donor + condition"
    assert design["replicate_handling"]["mode"] == "paired_biological_block"
    assert design["design_status"] == "ready"
    assert "Paired biological replicates" in dict(decisions)[2.5]
    assert dict(decisions)[2.6] == "Yes — proceed"


def _make_fake_report(tmp_path: Path) -> Path:
    report_dir = tmp_path / "report_exp1234"
    (report_dir / "figures" / "dex_vs_untrt").mkdir(parents=True)
    (report_dir / "tables").mkdir(parents=True)
    (report_dir / "report.html").write_text(
        "<html><a href='/home/test/run/dex_de_genes.tsv'>report</a></html>",
        encoding="utf-8",
    )
    (report_dir / "methodology.json").write_text(
        json.dumps({
            "design": {"design_formula": "~ donor + condition"},
            "contrasts": [
                {"name": "dexamethasone vs untreated", "n_significant": 24,
                 "n_upregulated": 14, "n_downregulated": 10},
            ],
            "source_table": "/home/test/run/dex_de_genes.tsv",
        }),
        encoding="utf-8",
    )
    (report_dir / "tables" / "dex_de_genes.tsv").write_text(
        "gene\tlog2fc\tpadj\nFKBP5\t3.2\t0.0001\n", encoding="utf-8"
    )
    (report_dir / "figures" / "dex_vs_untrt" / "volcano.svg").write_text(
        "<svg/>", encoding="utf-8"
    )
    return report_dir


def test_publish_assembles_portable_public_capsule(tmp_path):
    harness = _load_harness()
    report_dir = _make_fake_report(tmp_path)
    output_dir = tmp_path / "out"

    capsule = harness.publish(
        REPO_ROOT, report_dir, output_dir,
        experiment_id="exp1234", decisions=[(1, "Confirm"), (2.5, "Paired biological replicates")],
    )

    for name in ("capsule.json", "report.html", "methodology.json",
                 "de_results.tsv", "fig1_gse52778_bulk_de.svg"):
        assert (output_dir / name).is_file(), name

    assert capsule["schema_version"] == harness.CAPSULE_SCHEMA
    assert capsule["lane_id"] == "c1_public_fastq_e2e"
    assert capsule["design_formula_recorded"] == "~ donor + condition"
    assert "/home/" not in (output_dir / "report.html").read_text(encoding="utf-8")
    assert "/home/" not in (output_dir / "methodology.json").read_text(encoding="utf-8")
