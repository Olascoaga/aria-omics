"""FAIL->PASS regression tests for D5 CRITICAL findings (D5-001, 002, 005, 006, 008)."""

from __future__ import annotations

import json
import pytest

from aria.agents.biological_synthesis_agent import BiologicalSynthesisAgent
from aria.agents.chromatin_agent import ChromatinAgent
from aria.agents.narrative_agent import NarrativeAgent
from aria.agents.orchestrator_agent import OrchestratorAgent
from aria.memory.memory import ARIAMemory


def _bulk_ok() -> dict:
    return {
        "status": "done",
        "findings": {
            "padj_threshold": 0.05,
            "lfc_threshold": 1.0,
            "sample_qc": {"n_samples": 6, "size_ratio": 1.2, "outliers": []},
            "contrasts": [
                {
                    "name": "treat_vs_ctrl",
                    "numerator": "treat",
                    "denominator": "ctrl",
                    "status": "success",
                    "n_significant": 12,
                    "n_upregulated": 7,
                    "n_downregulated": 5,
                    "top_genes": [
                        {"gene": "TP53", "symbol": "TP53", "log2fc": 2.3, "padj": 1e-5}
                    ],
                }
            ],
        },
    }


def _make_agent(tmp_path, memory=None) -> NarrativeAgent:
    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = memory or type(
        "M", (), {"db_path": ":memory:", "get_decisions": lambda self, eid: []}
    )()
    agent.llm = None
    return agent


def test_d5_001_narrator_exception_isolated_preserves_healthy_modalities_and_emits_error_block(
    tmp_path,
):
    agent = _make_agent(tmp_path)
    agent_results = {
        "bulk_rna_agent": _bulk_ok(),
        "scrna_agent": {
            "status": "done",
            "findings": {
                "qc": {"status": "success", "n_cells_before": 1500, "n_cells_after": "1200"},
            },
        },
    }
    report = agent._render_html_report(
        experiment_id="d5_001",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report",
    )
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    claim_ids = [c["claim_id"] for c in meth["claims"]]
    err_blocks = [b for b in meth["narrative_blocks"] if b.get("status") == "error"]
    assert "bulk.contrast.treat_vs_ctrl" in claim_ids, (
        f"Healthy bulk RNA claim disappeared when scRNA narrator raised: claims={claim_ids}"
    )
    assert any(b["id"].startswith("scrna.") for b in err_blocks), (
        f"Failing scRNA narrator did not emit an error block: blocks={meth['narrative_blocks']}"
    )
    assert any("ValueError" in str(b.get("error") or "") for b in err_blocks)


def test_d5_002_run_ledger_guards_non_dict_design_intelligence(tmp_path):
    agent = _make_agent(tmp_path)
    report = agent._render_html_report(
        experiment_id="d5_002",
        exp_ctx={
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_intelligence": ["differential expression"],
        },
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"bulk_rna_agent": _bulk_ok()},
        report_dir=tmp_path / "report",
    )
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    claim_ids = [c["claim_id"] for c in meth["claims"]]
    assert "bulk.contrast.treat_vs_ctrl" in claim_ids, (
        f"Non-dict design_intelligence caused build_run_ledger to wipe bulk claims: {claim_ids}"
    )


def test_d5_002_failed_run_ledger_pass_surfaces_error_banner_and_withheld_count(
    tmp_path, monkeypatch
):
    from aria.agents.narrative import run_ledger as rl_mod

    def _boom(_exp_ctx, _agent_results):
        raise RuntimeError("simulated run_ledger crash")

    monkeypatch.setattr(rl_mod, "build_run_ledger", _boom)
    agent = _make_agent(tmp_path)
    report = agent._render_html_report(
        experiment_id="d5_002_crash",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"bulk_rna_agent": _bulk_ok()},
        report_dir=tmp_path / "report",
    )
    html = report.read_text(encoding="utf-8")
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    assert meth["run_ledger"].get("status") == "error"
    assert "simulated run_ledger crash" in str(meth["run_ledger"].get("error", ""))
    assert meth["claim_compilation"]["n_withheld"] >= 1
    assert "withheld" in html.lower()
    assert "ledger" in html.lower() and "simulated run_ledger crash" in html


def test_d5_005_orchestrator_crashed_modality_agent_reaches_report_and_ledger(tmp_path):
    orch = OrchestratorAgent(ARIAMemory(":memory:"), llm=None)
    err_result = orch._run_agent(
        "bulk_rna_agent", "exp_d5_005", {"exp_context": ["malformed_ctx"]}
    )
    assert err_result["status"] == "error"
    assert err_result.get("error"), f"_run_agent crash dict missing 'error' key: {err_result}"
    agent = _make_agent(tmp_path)
    report = agent._render_html_report(
        experiment_id="d5_005",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"bulk_rna_agent": err_result},
        report_dir=tmp_path / "report",
    )
    html = report.read_text(encoding="utf-8")
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    bulk_blocks = [b for b in meth["narrative_blocks"] if b["id"].startswith("bulk.")]
    bulk_ledger = [e for e in meth["run_ledger"]["entries"] if e["modality"] == "bulk"]
    assert any(b.get("status") == "error" for b in bulk_blocks), (
        f"Crashed bulk_rna_agent did not produce an error block: {bulk_blocks}"
    )
    assert any(e.get("status") == "error" and e.get("reason") for e in bulk_ledger), (
        f"Crashed bulk_rna_agent not marked status='error' in RunLedger: {bulk_ledger}"
    )
    assert err_result["error_type"] in html


def test_d5_005_cross_modal_crash_reaches_report_and_ledger(tmp_path):
    agent = _make_agent(tmp_path)
    cm_err = {
        "status": "error",
        "ran": False,
        "analysis": "cross_modal_concordance",
        "error_type": "RuntimeError",
        "error": "RuntimeError: stack failed",
        "details": "stack failed",
    }
    report = agent._render_html_report(
        experiment_id="d5_005_cm",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"cross_modal_concordance": cm_err},
        report_dir=tmp_path / "report",
    )
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    cm_blocks = [b for b in meth["narrative_blocks"] if "cross_modal" in b["id"]]
    assert any(b.get("status") == "error" for b in cm_blocks), (
        f"Crashed cross_modal_concordance vanished from narrative_blocks: {meth['narrative_blocks']}"
    )


def test_d5_006_biological_synthesis_handles_nan_contrast_and_surfaces_slice_errors(
    monkeypatch,
):
    valid_contrasts = [
        {
            "name": "A_vs_ctrl",
            "status": "success",
            "n_significant": 5,
            "all_sig_gene_ids": ["G1", "G2", "G3", "G4", "G5"],
            "all_sig_genes": ["TP53", "MDM2", "CDKN1A", "BAX", "BBC3"],
            "up_gene_ids": ["G1", "G2", "G3"],
            "down_gene_ids": ["G4", "G5"],
            "pathways": {"GO_BP": [{"term": "Apoptosis", "padj": 1e-4}]},
        },
        {
            "name": "B_vs_ctrl",
            "status": "success",
            "n_significant": 4,
            "all_sig_gene_ids": ["G1", "G2", "G6", "G7"],
            "all_sig_genes": ["TP53", "MDM2", "MYC", "FOS"],
            "up_gene_ids": ["G1", "G2"],
            "down_gene_ids": ["G6", "G7"],
            "pathways": {"GO_BP": [{"term": "Apoptosis", "padj": 1e-3}]},
        },
    ]
    baseline = BiologicalSynthesisAgent().synthesize(
        {"bulk_rna_agent": {"status": "done", "findings": {"contrasts": valid_contrasts}}}
    )
    assert len(baseline) >= 2
    degraded = BiologicalSynthesisAgent().synthesize(
        {
            "bulk_rna_agent": {
                "status": "done",
                "findings": {
                    "contrasts": valid_contrasts
                    + [{"name": "C_vs_ctrl", "status": "success", "n_significant": float("nan")}]
                },
            }
        }
    )
    assert len(degraded) >= 2, (
        f"BiologicalSynthesisAgent dropped bulk synthesis blocks on NaN contrast: {degraded}"
    )

    # Even when a synthesis slice genuinely raises, an error NarrativeBlock must be emitted.
    from aria.agents import biological_synthesis_agent as bsa_mod

    monkeypatch.setattr(
        bsa_mod,
        "detect_bulk_patterns",
        lambda _c: (_ for _ in ()).throw(RuntimeError("slice1a boom")),
    )
    err_blocks = BiologicalSynthesisAgent().synthesize(
        {"bulk_rna_agent": {"status": "done", "findings": {"contrasts": valid_contrasts}}}
    )
    assert any(
        b.status == "error" and "slice1a boom" in str(b.error or "") for b in err_blocks
    ), f"Synthesis slice crash did not emit an error block: {err_blocks}"


def test_d5_008_chromatin_epigenomic_metadata_encoding_detection_and_fail_closed(tmp_path):
    meta_csv = tmp_path / "sample_metadata.csv"
    meta_csv.write_bytes(
        b"sample_id,role,is_control,spikein_reads,notes\n"
        b"ip1,target,false,50000,H3K27ac 1 \xb5g\n"
        b"ctrl1,input,true,48000,IgG control 1 \xb5g\n"
    )
    ip_bam = tmp_path / "ip1.bam"
    ct_bam = tmp_path / "ctrl1.bam"
    ip_bam.write_bytes(b"")
    ct_bam.write_bytes(b"")

    agent = ChromatinAgent(ARIAMemory(":memory:"), llm=None)
    out = agent._resolve_epigenomic_samples(
        {"sample_metadata_path": str(meta_csv)}, [str(ip_bam), str(ct_bam)]
    )
    assert out["all_target_sids"] == ["ip1"], (
        f"Control BAM was merged into target IPs: {out}"
    )
    assert out["control_files"] == [str(ct_bam)]
    assert out["spikein_counts"] == {"ip1": 50000.0, "ctrl1": 48000.0}
    assert out.get("sample_metadata_encoding") in {"cp1252", "latin-1"}

    # UTF-8 with BOM must strip BOM cleanly so sample_id matches.
    bom_csv = tmp_path / "sample_metadata_bom.csv"
    bom_csv.write_bytes(
        b"\xef\xbb\xbfsample_id,role,is_control,spikein_reads\n"
        b"ip1,target,false,51000\n"
        b"ctrl1,input,true,49000\n"
    )
    out_bom = agent._resolve_epigenomic_samples(
        {"sample_metadata_path": str(bom_csv)}, [str(ip_bam), str(ct_bam)]
    )
    assert out_bom["all_target_sids"] == ["ip1"]
    assert out_bom["control_files"] == [str(ct_bam)]
    assert out_bom.get("sample_metadata_encoding") == "utf-8-sig"

    # Truly corrupt / binary metadata file must fail closed (never silently discard).
    corrupt_csv = tmp_path / "corrupt_metadata.csv"
    corrupt_csv.write_bytes(b"\x00\x01\x02BINARY_GARBAGE\x00")
    with pytest.raises((ValueError, OSError)):
        agent._resolve_epigenomic_samples(
            {"sample_metadata_path": str(corrupt_csv)}, [str(ip_bam), str(ct_bam)]
        )
