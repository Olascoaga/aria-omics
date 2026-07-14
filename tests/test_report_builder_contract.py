"""Characterization contract for the A7 ``ReportBuilderMixin`` extraction."""

from __future__ import annotations

from pathlib import Path

from aria.agents.narrative.report_builder import ReportBuilderMixin


_REPORT_BUILDER_METHODS = {
    "_collect_execution_llm_usage",
    "_build_report_dir",
    "_generate_scrna_figures",
    "_generate_chromatin_figures",
    "_speculative_verification_state",
    "_build_speculative_section_html",
    "_render_html_report",
    "_govern_executive_summary",
    "_unsupported_executive_summary_numbers",
    "_build_executive_summary_block",
    "_make_executive_summary_block",
    "_write_memory_snapshot",
    "_build_methodology_json",
    "_build_provenance_section",
    "_build_qc_section",
    "_build_findings_section",
    "_public_conflict_notice",
    "_build_methodology_table",
    "_build_bulk_rna_plots",
    "_stage_artifacts",
    "_build_findings_table",
    "_build_public_claims_table",
    "_build_decisions_table",
}


class _Host(ReportBuilderMixin):
    def _collect_tool_versions(self, _packages):
        return {"python": "test"}

    def _collect_execution_llm_usage(self, _since_utc=None):
        return {"calls": 0}

    def _collect_narrative_blocks(self, _agent_results, _exp_ctx):
        return []


def test_report_builder_exposes_the_complete_bound_method_contract():
    missing = {
        name for name in _REPORT_BUILDER_METHODS
        if not callable(getattr(ReportBuilderMixin, name, None))
    }
    assert missing == set()


def test_methodology_json_preserves_report_provenance_contract():
    host = _Host()
    provenance = {"timestamp_utc": "2026-07-13T00:00:00Z", "git_sha": "abc"}
    exp_ctx = {
        "input_files": [{"path": "input.tsv", "sha256": "123"}],
        "raw_ingestion": [{"status": "not_required"}],
        "design": {"main_factor": "condition"},
        "design_intelligence": {"status": "confirmed"},
    }
    ledger = {
        "entries": [],
        "divergences": [],
        "n_divergences": 0,
        "claim_ledger_verification": {"n_violations": 0},
    }
    compilation = {
        "compiler": "compile_public_claims",
        "n_published": 0,
        "n_withheld": 0,
        "withheld": [],
    }

    methodology = host._build_methodology_json(
        provenance=provenance,
        exp_ctx=exp_ctx,
        agent_results={},
        decisions=[{"checkpoint": "3", "decision": "standard"}],
        llm_usage={"calls": 0, "seed_deterministic": False},
        narrative_blocks=[],
        run_ledger=ledger,
        devils_advocate=[],
        compiled_claims=[],
        claim_compilation=compilation,
    )

    assert set(methodology) == {
        "provenance",
        "inputs",
        "raw_ingestion",
        "narrative_blocks",
        "claims",
        "claim_compilation",
        "devils_advocate",
        "run_ledger",
        "robustness_multiverse",
        "design",
        "design_intelligence",
        "thresholds",
        "seeds",
        "tools",
        "llm_usage",
        "decisions",
    }
    assert methodology["provenance"] == provenance
    assert methodology["inputs"] == exp_ctx["input_files"]
    assert methodology["claims"] == []
    assert methodology["claim_compilation"] == compilation
    assert methodology["run_ledger"]["claim_linkage"] == {
        "linked": 0,
        "unlinked": 0,
        "violations": [],
        "n_violations": 0,
    }
    assert methodology["llm_usage"]["seed_deterministic"] is False
    assert methodology["seeds"] == {"global": 0, "scanpy": 0, "harmony": 0}


def test_stage_artifacts_keeps_bulk_destination_layout(tmp_path: Path):
    contrast_dir = tmp_path / "work" / "condition_a_vs_condition_b"
    figures = contrast_dir / "figures"
    tables = contrast_dir / "tables"
    figures.mkdir(parents=True)
    tables.mkdir()
    (figures / "volcano.svg").write_text("<svg/>", encoding="utf-8")
    (tables / "de_genes.tsv").write_text("gene\nGENE_1\n", encoding="utf-8")
    (contrast_dir.parent / "counts_tpm.tsv").write_text(
        "gene\tsample\nGENE_1\t1\n", encoding="utf-8"
    )
    report_dir = tmp_path / "report"
    (report_dir / "figures").mkdir(parents=True)
    (report_dir / "tables").mkdir()
    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "contrasts": [{
                    "status": "success",
                    "contrast_dir": str(contrast_dir),
                }],
                "sample_qc": {},
            },
        }
    }

    _Host()._stage_artifacts(agent_results, report_dir)

    assert (
        report_dir / "figures" / contrast_dir.name / "volcano.svg"
    ).read_text(encoding="utf-8") == "<svg/>"
    assert (
        report_dir / "tables" / f"{contrast_dir.name}_de_genes.tsv"
    ).is_file()
    assert (report_dir / "tables" / "counts_tpm.tsv").is_file()
