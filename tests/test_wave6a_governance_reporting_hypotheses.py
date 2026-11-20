"""Wave 6A remediation tests: Governance, Reporting, Hypothesis Engine, Synthesis, Utils, and Top-Level Agents."""

from __future__ import annotations

import json
import math
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from aria.agents.data_audit_agent import DataAuditAgent
from aria.agents.design_agent import DesignAgent
from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.claim_compiler import classify_claim
from aria.agents.narrative.devils_advocate import build_devils_advocate
from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    SpeculativePromotionError,
    assert_no_speculative_promotion,
    build_evidence_index,
    build_proposer_prompt,
    build_signals_by_entity,
    bulk_rna_evidence,
    check_devils_advocate,
    inherited_caveat_codes,
    is_hypothesis_node,
    persist_speculative_manifest,
    rank_hypotheses,
    render_speculative_section_html,
    scrna_evidence,
    verify_hypothesis_grounding,
    visible_confounds,
)
from aria.agents.narrative.hypothesis.adapters._chromatin_common import (
    motif_signals,
    num as chromatin_num,
)
from aria.agents.narrative.hypothesis.ranking import _normalized_effect
from aria.agents.narrative.ledger_export import (
    write_reproducible_capsule,
    write_ro_crate,
)
from aria.agents.narrative.render_blocks import render_blocks
from aria.agents.narrative.report_builder import ReportBuilderMixin
from aria.agents.narrative.reporting import sections as report_sections_mod
from aria.agents.narrative import report_sections
from aria.agents.narrative.synthesis.cross_modal_triangulation import (
    build_cross_modal_evidence_graph,
)
from aria.agents.narrative.synthesis.footprint_rna import (
    _rna_favored_group,
    footprint_rna_concordance,
)
from aria.agents.narrative.synthesis.pattern_detector import (
    detect_bulk_patterns,
    detect_scrna_patterns,
)
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.utils import cross_modal_concordance as cm
from aria.utils.decision_manifest import (
    build_decision_manifest,
    commit_decision_manifest,
    verify_decision_manifest,
)
from aria.utils.power_estimation import bulk_power_estimate
from aria.utils.stats import contrast_family_significance


# ── 1. D5-H6: structured doublet report in devils_advocate ───────────────────


def test_d5h6_devils_advocate_reads_structured_doublet_flags():
    def _block():
        b = NarrativeBlock(
            id="scrna.pseudobulk.T.a_vs_b",
            modality="scRNA",
            analysis="pseudobulk_de",
            block_type="result",
            title="DE",
            status="success",
            confidence="medium",
            claim="Found 10 DE genes.",
            evidence=[EvidenceItem("DE genes", 10, "pseudobulk_de")],
        )
        b.metadata["claim"] = {"tier": "associative"}
        return b

    # Structured scrublet_report with ran=False MUST NOT flip doublets to True
    res_false = {
        "scrna_agent": {
            "findings": {
                "qc": {
                    "status": "success",
                    "scrublet_report": {"ran": False, "reason": "skipped"},
                }
            }
        }
    }
    m1 = build_devils_advocate([_block()], res_false, {})
    ch1 = {c["alternative"]: c["addressed"] for c in m1[0]["challenges"]}
    assert ch1["doublets"] is False

    # Structured doublet_detection with ran=False MUST also be False
    res_dbl_false = {
        "scrna_agent": {
            "findings": {
                "qc": {
                    "status": "success",
                    "doublet_detection": {"ran": False, "method": "scrublet"},
                }
            }
        }
    }
    m2 = build_devils_advocate([_block()], res_dbl_false, {})
    ch2 = {c["alternative"]: c["addressed"] for c in m2[0]["challenges"]}
    assert ch2["doublets"] is False

    # Structured scrublet_report with ran=True -> addressed=True
    res_true = {
        "scrna_agent": {
            "findings": {
                "qc": {
                    "status": "success",
                    "scrublet_report": {"ran": True, "n_doublets": 42},
                }
            }
        }
    }
    m3 = build_devils_advocate([_block()], res_true, {})
    ch3 = {c["alternative"]: c["addressed"] for c in m3[0]["challenges"]}
    assert ch3["doublets"] is True


# ── 2. D5-H8: footprint_rna finite, ties, case-insensitive, dimer support ────


def test_d5h8_footprint_rna_finite_ties_case_insensitive_and_dimers():
    # Non-finite, None, and bool means are skipped cleanly
    grp, means = _rna_favored_group(
        ["NAN_TF", "INF_TF", "NONE_TF", "BOOL_TF"],
        {
            "NAN_TF": {"A": float("nan"), "B": 1.0},
            "INF_TF": {"A": float("inf"), "B": 1.0},
            "NONE_TF": {"A": None, "B": 1.0},
            "BOOL_TF": {"A": True, "B": False},
        },
        "A",
        "B",
    )
    assert grp is None and means is None

    # Ties (e.g. 0.0 == 0.0 unexpressed TF) return neutral None, NOT group_a!
    summary = {
        "top_toward_A": [
            {"tf": "UNEXPR"},
            {"tf": "Cebpb"},  # case-insensitive lookup
            {"tf": "FOSL2::JUND"},  # dimer: FOSL2 strongly favors B, JUND supports A
        ],
        "top_toward_B": [
            {"tf": "TIED_B"},
        ],
    }
    rna = {
        "UNEXPR": {"A": 0.0, "B": 0.0},
        "TIED_B": {"A": 2.5, "B": 2.5},
        "cebpb": {"A": 4.0, "B": 1.0},
        "FosL2": {"A": 1.0, "B": 9.0},  # discordant with larger |diff| = 8
        "JunD": {"A": 2.0, "B": 1.0},   # concordant with A (|diff| = 1)
    }
    res = footprint_rna_concordance(summary, rna, "A", "B")
    toward_a = {row["tf"]: row for row in res["toward_A"]}
    toward_b = {row["tf"]: row for row in res["toward_B"]}

    # Ties must be neutral (rna_checked=False, rna_favored_group=None, rna_concordant=None)
    assert toward_a["UNEXPR"]["rna_favored_group"] is None
    assert toward_a["UNEXPR"]["rna_checked"] is False
    assert toward_a["UNEXPR"]["rna_concordant"] is None
    assert toward_b["TIED_B"]["rna_favored_group"] is None
    assert toward_b["TIED_B"]["rna_checked"] is False

    # Case-insensitive lookup matched "cebpb"
    assert toward_a["Cebpb"]["rna_checked"] is True
    assert toward_a["Cebpb"]["rna_concordant"] is True

    # Dimer FOSL2::JUND matched concordant subunit JunD even though FosL2 had larger |diff|
    assert toward_a["FOSL2::JUND"]["rna_checked"] is True
    assert toward_a["FOSL2::JUND"]["rna_favored_group"] == "A"
    assert toward_a["FOSL2::JUND"]["rna_concordant"] is True
    assert toward_a["FOSL2::JUND"]["rna_means"]["gene"] == "JUND"

    assert res["n_checked"] == 2
    assert res["n_concordant"] == 2
    assert res["concordance_rate"] == 1.0


# ── 3. D5-H11: canonical entity indexing + observed_claims caveat union ──────


def test_d5h11_canonical_variants_and_observed_claims_cannot_bypass_confounds():
    sig_nkx = EvidenceSignal(
        entity="NKX2-1",
        entity_kind="gene",
        modality="bulk_RNA",
        measure="log2fc",
        audited_node_ref="ledger://bulk/differential_expression",
        value=2.0,
        direction="up",
        context="nkx_ko_vs_wt",
        caveats_inherited=["batch"],
    )
    sig_nr1d1 = EvidenceSignal(
        entity="NR1D1",
        entity_kind="gene",
        modality="bulk_RNA",
        measure="log2fc",
        audited_node_ref="ledger://bulk/differential_expression",
        value=-2.5,
        direction="down",
        context="rev_erba_ko_vs_wt",
        caveats_inherited=["low_replication"],
    )
    sig_reverb = EvidenceSignal(
        entity="REV-ERBa",
        entity_kind="tf_motif",
        modality="bulk_ATAC",
        measure="motif_enrich",
        audited_node_ref="ledger://chromatin/motif_enrichment",
        value=1.8,
        direction="up",
        context="rev_erba_ko_vs_wt",
        caveats_inherited=["motif_not_binding"],
    )

    sigs = [sig_nkx, sig_nr1d1, sig_reverb]
    by_ent = build_signals_by_entity(sigs)
    by_id = {s.signal_id: s for s in sigs}

    # NKX2_1 <-> NKX2-1 canonical variant matches via both build_signals_by_entity AND raw dict
    hyp_nkx = Hypothesis(
        id="h_nkx",
        mechanism="NKX2_1 may regulate target expression",
        entities=["NKX2_1"],
        devils_advocate={"simpler_explanation": "shared regulator", "confounds": []},
    )
    assert "batch" in visible_confounds(hyp_nkx, by_ent)
    assert "batch" in visible_confounds(hyp_nkx, {"NKX2-1": sig_nkx})
    assert check_devils_advocate(hyp_nkx, by_ent).passed is False

    # REV-ERBα <-> REV-ERBa Greek variant matches motif_not_binding AND citing NR1D1 in observed_claims unions low_replication!
    hyp_reverb = Hypothesis(
        id="h_rev",
        mechanism="REV-ERBα loss may alter target transcription",
        entities=["REV-ERBα"],
        observation_refs=[
            "ledger://bulk/differential_expression",
            "ledger://chromatin/motif_enrichment",
        ],
        observed_claims=[
            {"signal_id": sig_nr1d1.signal_id, "stated_direction": "down"},
        ],
        experiment=DiscriminatingExperiment(
            "REV-ERBα rescue",
            "NR1D1 expression by RT-qPCR",
            "increase",
            "no change",
        ),
        devils_advocate={
            "simpler_explanation": "shared circadian shift",
            "confounds": ["motif_not_binding"],  # omitted low_replication!
        },
    )
    confs = visible_confounds(hyp_reverb, by_ent, signals_by_id=by_id)
    assert {"motif_not_binding", "low_replication"} <= confs
    assert check_devils_advocate(hyp_reverb, by_ent, signals_by_id=by_id).passed is False

    # When HypothesisAgent runs end-to-end, omitting low_replication rejects; owning both passes and stamps both!
    ledger = {
        "entries": [
            {"node_id": "ledger://bulk/differential_expression", "status": "ran"},
            {"node_id": "ledger://chromatin/motif_enrichment", "status": "ran"},
        ]
    }
    out_rej = HypothesisAgent(proposer=lambda s, c: [hyp_reverb]).generate(
        sigs, ledger, w_claim_passed=True, w_ledger_passed=True
    )
    assert out_rej["hypotheses"] == []
    assert any(
        f["gate"] == "devils_advocate" for f in out_rej["rejected"][0]["failures"]
    )

    hyp_reverb.devils_advocate["confounds"] = ["motif_not_binding", "low_replication"]
    out_ok = HypothesisAgent(proposer=lambda s, c: [hyp_reverb]).generate(
        sigs, ledger, w_claim_passed=True, w_ledger_passed=True
    )
    assert [h["id"] for h in out_ok["hypotheses"]] == ["h_rev"]
    assert set(out_ok["hypotheses"][0]["inherited_caveats"]) >= {
        "motif_not_binding",
        "low_replication",
    }
    assert build_evidence_index([sig_nkx])["nkx21"] is sig_nkx


# ── 4. D5-H12: Hi-C loop ranking & case-insensitive quarantine ───────────────


def test_d5h12_hic_loop_ranking_and_case_insensitive_quarantine():
    strong_loop = EvidenceSignal(
        entity="MYC",
        entity_kind="gene",
        modality="HiC",
        measure="hic_loop_anchor",
        audited_node_ref="ledger://hic/loops",
        value=1e-12,
        context="chr8:1-2<->chr8:3-4",
    )
    weak_loop = EvidenceSignal(
        entity="SOX2",
        entity_kind="gene",
        modality="HiC",
        measure="hic_loop_anchor",
        audited_node_ref="ledger://hic/loops",
        value=0.049,
        context="chr3:1-2<->chr3:3-4",
    )
    null_loop = EvidenceSignal(
        entity="ACTB",
        entity_kind="gene",
        modality="HiC",
        measure="hic_loop_anchor",
        audited_node_ref="ledger://hic/loops",
        value=1.0,
        context="chr7:1-2<->chr7:3-4",
    )

    e_strong = _normalized_effect(strong_loop)
    e_weak = _normalized_effect(weak_loop)
    e_null = _normalized_effect(null_loop)
    assert 0.0 <= e_null < e_weak < e_strong <= 1.0
    assert e_strong > 0.9 and e_weak < 0.35

    exp = DiscriminatingExperiment("KD", "qPCR", "decrease", "no change")
    h_weak = Hypothesis(
        id="weak",
        mechanism="SOX2 loop may act",
        entities=["SOX2"],
        observed_claims=[{"signal_id": weak_loop.signal_id, "stated_direction": "na"}],
        experiment=exp,
    )
    h_strong = Hypothesis(
        id="strong",
        mechanism="MYC loop may act",
        entities=["MYC"],
        observed_claims=[{"signal_id": strong_loop.signal_id, "stated_direction": "na"}],
        experiment=exp,
    )
    ranked = rank_hypotheses([h_weak, h_strong], [weak_loop, strong_loop])
    assert [h.id for h in ranked] == ["strong", "weak"]

    # Case-insensitive tier and anywhere-in-string hypothesis:// detection
    assert is_hypothesis_node("HYPOTHESIS://abc") is True
    assert is_hypothesis_node("derived from Hypothesis://h1 in note") is True
    assert is_hypothesis_node("ledger://bulk/differential_expression") is False

    for leaked_claim in (
        {"claim_id": "c1", "tier": "speculative"},
        {"claim_id": "c2", "tier": " Speculative "},
        {"claim_id": "c3", "tier": "associative", "note": "ref: HYPOTHESIS://h9"},
        {"claim_id": "c4", "tier": "associative", "evidence": [{"ref": "see hypothesis://x"}]},
    ):
        with pytest.raises(SpeculativePromotionError):
            assert_no_speculative_promotion([leaked_claim])


# ── 5. D5-H13, D5-M5, D4-M10b: QC fallbacks, XSS, staged table hrefs, render ─


def test_d5h13_qc_fallbacks_and_audit_xss_and_none_safety():
    class _Sec(report_sections_mod.ReportSectionsMixin):
        pass

    grouped = {"high": [], "medium": [], "low": [], "insufficient": []}
    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "bulk_ATAC": {
                    "findings": {
                        "peaks": {"n_peaks": 45000, "frip": 0.31},
                        "qc": {"tss_enrichment": 7.4},
                    }
                }
            },
        },
        "genome_arch_agent": {
            "status": "done",
            "findings": {"qc": {"status": "balanced"}},
        },
    }
    exp_ctx = {
        "audit_findings": [
            {
                "severity": "warning",
                "check": "<script>alert(1)</script>_check",
                "message": "<b>bad</b> & unsafe",
                "recommendation": None,
            },
            {
                "severity": "blocking",
                "check": None,
                "message": None,
                "recommendation": "<img src=x onerror=alert(2)>",
            },
        ]
    }
    html = _Sec()._build_qc_section(grouped, agent_results, exp_ctx)
    assert "TSS enrichment 7.4" in html
    assert "ICE balancing: balanced" in html
    assert "<script>" not in html
    assert "<img " not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html.lower()


def test_d5m5_staged_bulk_table_links_are_relative_and_reconciled(tmp_path: Path):
    work = tmp_path / "work" / "ko_vs_wt"
    (work / "tables").mkdir(parents=True)
    (work / "gsea").mkdir(parents=True)
    de_tsv = work / "tables" / "de_genes.tsv"
    gsea_csv = work / "gsea" / "gsea_results.csv"
    full_csv = work / "deseq2_full.csv"
    de_tsv.write_text("gene\tlog2FoldChange\nG1\t2.0\n", encoding="utf-8")
    gsea_csv.write_text("Term,NES,FDR\nP1,2.0,0.01\n", encoding="utf-8")
    full_csv.write_text("gene,log2FoldChange\nG1,2.0\n", encoding="utf-8")

    report_dir = tmp_path / "report"
    (report_dir / "figures").mkdir(parents=True)
    (report_dir / "tables").mkdir(parents=True)

    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "contrasts": [
                    {
                        "name": "ko_vs_wt",
                        "status": "success",
                        "contrast_dir": str(work),
                        "full_results_csv": str(full_csv),
                        "plots": {
                            "tables": {"de_genes": str(de_tsv)},
                            "gsea_table": str(gsea_csv),
                        },
                    }
                ]
            },
        }
    }

    class _Sec(report_sections_mod.ReportSectionsMixin):
        pass

    sec = _Sec()
    sec._stage_artifacts(agent_results, report_dir)

    assert (report_dir / "tables" / "ko_vs_wt_de_genes.tsv").is_file()
    assert (report_dir / "tables" / "ko_vs_wt_gsea_results.csv").is_file()
    assert (report_dir / "tables" / "ko_vs_wt_deseq2_full.csv").is_file()

    block = NarrativeBlock(
        id="bulk.contrast.ko_vs_wt",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="KO vs WT",
        status="success",
        confidence="medium",
        claim="Bulk contrast ko_vs_wt had 5 DE genes.",
        evidence=[EvidenceItem("DE genes", 5, "bulk_de")],
        tables=[
            {"id": "de_genes", "path": str(de_tsv), "label": "DE genes"},
            {"id": "gsea_table", "path": str(gsea_csv), "label": "GSEA"},
            {"id": "full_results", "path": str(full_csv), "label": "Full DE"},
        ],
    )
    html_out = render_blocks(
        [block],
        report_dir=report_dir,
        strict=True,
        staged_tables=sec._staged_table_map,
    )
    assert "href='tables/ko_vs_wt_de_genes.tsv'" in html_out
    assert "href='tables/ko_vs_wt_gsea_results.csv'" in html_out
    assert "href='tables/ko_vs_wt_deseq2_full.csv'" in html_out
    assert str(tmp_path) not in html_out

    # Also resolves cleanly without explicit staged_tables map via slug heuristics
    html_standalone = render_blocks([block], report_dir=report_dir, strict=True)
    assert "href='tables/ko_vs_wt_de_genes.tsv'" in html_standalone
    assert "href='tables/ko_vs_wt_gsea_results.csv'" in html_standalone
    assert "href='tables/ko_vs_wt_deseq2_full.csv'" in html_standalone


def test_d4m10b_render_filters_withheld_claims_from_devils_advocate(tmp_path: Path):
    from aria.agents.narrative_agent import NarrativeAgent

    def _make_agent(rdir: Path, blocks: list[NarrativeBlock]) -> NarrativeAgent:
        agent = NarrativeAgent.__new__(NarrativeAgent)
        agent.reports_dir = rdir
        agent.memory = type("M", (), {"db_path": ":memory:"})()
        agent._collect_narrative_blocks = lambda _ar, _ctx: list(blocks)
        agent._collect_tool_versions = lambda *a, **k: {}
        agent._fallback_executive_summary = lambda _g, _i: "ARIA completed the analysis."
        agent._summarize_agent_results_for_llm = lambda _ar: "10 DE genes"
        return agent

    good_block = NarrativeBlock(
        id="bulk.contrast.good",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="Good contrast",
        status="success",
        confidence="medium",
        claim="Good contrast found 10 DE genes.",
        evidence=[EvidenceItem("DE genes", 10, "bulk_de")],
        metrics={"n_significant": 10},
    )
    bad_block = NarrativeBlock(
        id="bulk.contrast.withheld_bad",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="Withheld contrast",
        status="success",
        confidence="medium",
        # Unsupported number 999 -> withheld by compile_public_claims
        claim="Withheld contrast found 999 DE genes.",
        evidence=[EvidenceItem("DE genes", 10, "bulk_de")],
        metrics={"n_significant": 10},
    )
    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {"contrasts": [{"name": "good", "status": "success"}]},
        }
    }
    harness = _make_agent(tmp_path, [good_block, bad_block])
    report_path = harness._render_html_report(
        experiment_id="d4m10b_test",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "test"},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="methods",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "rep",
    )
    html_text = report_path.read_text(encoding="utf-8")
    meth = json.loads((report_path.parent / "methodology.json").read_text(encoding="utf-8"))
    da_ids = {e["claim_id"] for e in meth["devils_advocate"]}
    claim_ids = {c["claim_id"] for c in meth["claims"]}
    assert "bulk.contrast.good" in claim_ids
    assert "bulk.contrast.good" in da_ids
    assert "bulk.contrast.withheld_bad" not in claim_ids
    assert "bulk.contrast.withheld_bad" not in da_ids
    assert "bulk.contrast.withheld_bad" not in html_text
    assert "Good contrast found 10 DE genes." in html_text

    # Also verify render-time withholding reconciles HTML Findings <-> methodology claims & devils_advocate
    good_block2 = NarrativeBlock(
        id="bulk.contrast.good",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="Good contrast",
        status="success",
        confidence="medium",
        claim="Good contrast found 10 DE genes.",
        evidence=[EvidenceItem("DE genes", 10, "bulk_de")],
        metrics={"n_significant": 10},
    )
    harness2 = _make_agent(tmp_path, [good_block2])
    harness2._build_findings_section = lambda *a, **k: (
        '<section class="narrative-block" data-block-id="bulk.contrast.good" '
        'data-withheld="true"><h4>withheld</h4></section>'
    )
    report_path2 = harness2._render_html_report(
        experiment_id="d4m10b_reconcile",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "test"},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="methods",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "rep2",
    )
    meth2 = json.loads((report_path2.parent / "methodology.json").read_text(encoding="utf-8"))
    assert "bulk.contrast.good" not in {c["claim_id"] for c in meth2["claims"]}
    assert "bulk.contrast.good" not in {e["claim_id"] for e in meth2["devils_advocate"]}
    assert any(w["claim_id"] == "bulk.contrast.good" for w in meth2["claim_compilation"]["withheld"])


# ── 6. D5-M1 & D5-M4: byte-identical replay, strict JSON allow_nan=False ────


def test_d5m1_d5m4_reproducible_manifest_byte_identity_and_strict_json(tmp_path: Path):
    ledger_path = tmp_path / "fork_ledger.jsonl"
    ws1, ws2 = tmp_path / "ws1", tmp_path / "ws2"
    ctx = {
        "organism": "Mus musculus",
        "genome": "mm39",
        "reproducible_mode": True,
        "modalities": {"bulk_RNA": ["/data/a.tsv"]},
        "input_files": [
            {"modality": "bulk_RNA", "path": "/data/a.tsv", "size_bytes": 10, "sha256": "a" * 64}
        ],
        "design": {"organism": "Mus musculus", "genome": "mm39", "main_factor": "cond"},
    }
    plan = {"steps": []}
    intent = {"analysis_type": "differential", "raw_user_question": "A vs B"}

    m1 = build_decision_manifest(
        experiment_id="exp_det",
        exp_context=ctx,
        plan=plan,
        biological_intent=intent,
        locked_parameter_decisions={"bad_nan": float("nan"), "bad_inf": float("inf")},
    )
    em1, fl1 = commit_decision_manifest(m1, workspace_dir=ws1, ledger_path=ledger_path)
    m2 = build_decision_manifest(
        experiment_id="exp_det",
        exp_context=ctx,
        plan=plan,
        biological_intent=intent,
        locked_parameter_decisions={"bad_nan": float("nan"), "bad_inf": float("inf")},
    )
    em2, fl2 = commit_decision_manifest(m2, workspace_dir=ws2, ledger_path=ledger_path)

    # Byte-identical decision_manifest.json across Run 1 and exact Replay Run 2!
    b1 = (ws1 / "decision_manifest.json").read_bytes()
    b2 = (ws2 / "decision_manifest.json").read_bytes()
    assert b1 == b2
    assert "ledger_runs_observed" not in em1["fork_ledger_summary"]
    assert "ledger_runs_observed" not in em2["fork_ledger_summary"]
    assert fl1["ledger_runs_observed"] == 1
    assert fl2["ledger_runs_observed"] == 2
    assert verify_decision_manifest(em1) == (True, "verified")
    # Non-finite values sanitized to null (valid RFC-8259 JSON)
    parsed_dm = json.loads(b1.decode("utf-8"))
    assert parsed_dm["locked_parameter_decisions"]["bad_nan"] is None
    assert parsed_dm["locked_parameter_decisions"]["bad_inf"] is None

    # External fork ledger still has both runs with ledger_runs_observed = 1, 2
    lines = [json.loads(ln) for ln in ledger_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert [ln["ledger_runs_observed"] for ln in lines] == [1, 2]

    # ledger_export sanitizes non-finite floats with allow_nan=False
    rd = tmp_path / "rep_nan"
    rd.mkdir()
    meth = {
        "provenance": {"version": "4.9", "git_commit": "abc123", "score": float("nan")},
        "inputs": [{"path": "/tmp/in.tsv", "sha256": "deadbeef", "bytes": 10}],
        "claims": [],
        "run_ledger": {"entries": [], "n_divergences": 0},
        "calibration": {"status": "pass", "metric": float("inf")},
    }
    (rd / "methodology.json").write_text(json.dumps(meth), encoding="utf-8")
    (rd / "report.html").write_text("<html></html>", encoding="utf-8")
    crate_path = write_ro_crate(rd)
    assert "NaN" not in crate_path.read_text(encoding="utf-8")
    cap_path = write_reproducible_capsule(rd, repo_root=tmp_path / "empty_repo")
    with zipfile.ZipFile(cap_path) as zf:
        for member in ("capsule_manifest.json", f"{rd.name}/methodology.json", f"{rd.name}/ro-crate-metadata.json"):
            text = zf.read(member).decode("utf-8")
            assert "NaN" not in text and "Infinity" not in text

    # hypothesis/report_section: sanitizes NaN in persist_speculative_manifest + renders blocked card for verification_evidence_absent
    section_nan = {
        "ran": True,
        "honest_null": False,
        "hypotheses": [{"id": "h1", "rank_evidence": {"mean_effect_norm": float("nan")}}],
        "quarantine": [],
        "rejected": [],
    }
    sp_path = persist_speculative_manifest(section_nan, tmp_path)
    sp_text = sp_path.read_text(encoding="utf-8")
    assert "NaN" not in sp_text
    assert json.loads(sp_text)["hypotheses"][0]["rank_evidence"]["mean_effect_norm"] is None

    blocked_html = render_speculative_section_html(
        {"ran": False, "reason": "verification_evidence_absent"}
    )
    assert "SPECULATIVE" in blocked_html
    assert "W-CLAIM/W-LEDGER" in blocked_html
    assert "withheld" in blocked_html


# ── 7. D5-M2 & D3-M3: stats, power_estimation, pattern_detector, triangulation ─


def test_d5m2_d3m3_stats_power_pattern_detector_and_triangulation():
    # stats.py: NaN/inf/bool LFC never corrupts n_sig/n_up/n_down
    group_stats = {
        "c1": {
            "G_OK_UP": {"pvalue": 1e-6, "log2fc": 2.5},
            "G_OK_DN": {"pvalue": 1e-6, "log2fc": -2.0},
            "G_NAN": {"pvalue": 1e-6, "log2fc": float("nan")},
            "G_INF": {"pvalue": 1e-6, "log2fc": float("inf")},
            "G_NINF": {"pvalue": 1e-6, "log2fc": float("-inf")},
            "G_BOOL": {"pvalue": 1e-6, "log2fc": True},
        }
    }
    res_lfc = contrast_family_significance(group_stats, padj_max=0.05, lfc_min=1.0)
    assert res_lfc["c1"]["sig_genes"] == ["G_OK_UP", "G_OK_DN"]
    assert (res_lfc["c1"]["n_sig"], res_lfc["c1"]["n_up"], res_lfc["c1"]["n_down"]) == (2, 1, 1)

    res_nol = contrast_family_significance(group_stats, padj_max=0.05, lfc_min=None)
    assert res_nol["c1"]["sig_genes"] == ["G_OK_UP", "G_OK_DN"]
    assert (res_nol["c1"]["n_sig"], res_nol["c1"]["n_up"], res_nol["c1"]["n_down"]) == (2, 1, 1)

    # power_estimation.py: negative/non-finite dispersion & <=0/non-finite mean_expression -> 0.0
    assert bulk_power_estimate(4, dispersion=-0.2) == 0.0
    assert bulk_power_estimate(4, dispersion=float("nan")) == 0.0
    assert bulk_power_estimate(4, dispersion=float("inf")) == 0.0
    assert bulk_power_estimate(4, mean_expression=0.0) == 0.0
    assert bulk_power_estimate(4, mean_expression=-10.0) == 0.0
    assert bulk_power_estimate(4, mean_expression=float("nan")) == 0.0
    assert bulk_power_estimate(4, mean_expression=float("inf")) == 0.0
    assert bulk_power_estimate(4, mean_expression=100.0, dispersion=0.0) > 0.0

    # pattern_detector.py: reference alias, fdr_strategy, finite non-bool powers
    bulk = detect_bulk_patterns(
        [
            {
                "name": "A_vs_WT",
                "reference": "WT",
                "status": "success",
                "n_significant": 2,
                "all_sig_gene_ids": ["g1", "g2"],
                "power_estimate_at_lfc_min": float("nan"),
            },
            {
                "name": "B_vs_WT",
                "reference": "WT",
                "status": "success",
                "n_significant": 2,
                "all_sig_gene_ids": ["g2", "g3"],
                "power_estimate_at_lfc_min": True,
            },
            {
                "name": "C_vs_WT",
                "denominator": "WT",
                "status": "success",
                "n_significant": 1,
                "all_sig_gene_ids": ["g2"],
                "power_estimate_at_lfc_min": 0.72,
            },
        ]
    )
    assert all(x["shared_reference"] == "WT" for x in bulk["cross_contrast"])
    assert bulk["reliability"]["min_power"] == 0.72
    assert bulk["reliability"]["max_power"] == 0.72

    scrna = detect_scrna_patterns(
        {
            "findings": {
                "qc": {"n_cells_after": 500},
                "pseudobulk_de": {
                    "multiple_testing": {"fdr_strategy": "per_cluster"},
                    "per_group": {
                        "Tcell": {
                            "per_comparison": {
                                "treat_vs_ctrl": {
                                    "status": "success",
                                    "n_significant_global": 0,
                                    "n_significant_local": 15,
                                    "n_up_local": 9,
                                    "n_down_local": 6,
                                    "power_estimate_at_lfc_min": float("inf"),
                                }
                            }
                        }
                    },
                },
            }
        }
    )
    assert scrna["n_pseudobulk_blocks"] == 1
    assert scrna["strongest_pseudobulk"]["n_de"] == 15
    assert scrna["strongest_pseudobulk"]["n_up"] == 9
    assert scrna["strongest_pseudobulk"]["n_down"] == 6
    assert scrna["reliability"]["min_power"] is None

    # cross_modal_triangulation.py: all "neg_" sign classes route to discordant_signal
    graph = build_cross_modal_evidence_graph(
        {
            "integration_agent": {
                "status": "done",
                "findings": {
                    "eregulons": {
                        "status": "success",
                        "top_eregulons": [
                            {"tf": "TF1", "sign_class": "neg_tf_pos_p2g", "top_target_genes": ["G1"]},
                            {"tf": "TF2", "sign_class": "pos_tf_neg_p2g", "top_target_genes": ["G2"]},
                            {"tf": "TF3", "sign_class": "neg_tf_neg_p2g", "top_target_genes": ["G3"]},
                            {"tf": "TF4", "sign_class": "pos_tf_pos_p2g", "top_target_genes": ["G4"]},
                        ],
                    }
                },
            },
            "chromatin_agent": {
                "status": "done",
                "findings": {
                    "motifs": {
                        "ran": True,
                        "per_group": {
                            "up": {
                                "top_motifs": [
                                    {"name": "TF1", "padj": 0.001},
                                    {"name": "TF2", "padj": 0.001},
                                    {"name": "TF3", "padj": 0.001},
                                    {"name": "TF4", "padj": 0.001},
                                ]
                            }
                        },
                    }
                },
            },
        }
    )
    by_id = {a["axis_id"]: a for a in graph["regulatory_axes"]}
    assert by_id["eregulon:TF1:neg_tf_pos_p2g"]["classification"] == "discordant_signal"
    assert by_id["eregulon:TF2:pos_tf_neg_p2g"]["classification"] == "discordant_signal"
    assert by_id["eregulon:TF3:neg_tf_neg_p2g"]["classification"] == "discordant_signal"
    assert by_id["eregulon:TF4:pos_tf_pos_p2g"]["classification"] == "convergent_cross_contrast"


# ── 8. D5-M3 & CLD-04: rna_not_tested, adapters, grounding, prompts ─────────


def test_d5m3_cld04_tf_rna_not_tested_adapters_grounding_and_questions():
    # cross_modal_concordance.tf_layer -> rna_not_tested when expressed-but-untested or padj=NaN
    rna_df = pd.DataFrame(
        [
            {"gene_key": "COOKS_TF", "rna_log2fc": 2.0, "rna_padj": float("nan"), "rna_direction": 1, "rna_known_sig": False},
            {"gene_key": "UNCHANGED_TF", "rna_log2fc": 0.1, "rna_padj": 0.8, "rna_direction": 1, "rna_known_sig": False},
        ]
    )
    atac = {
        "motifs": {
            "per_group": {
                "up": {
                    "top_motifs": [
                        {"name": "COOKS_TF", "padj": 0.001, "log2_enrichment": 1.5},
                        {"name": "EXPR_NO_DE", "padj": 0.001, "log2_enrichment": 1.5},
                        {"name": "ABSENT_TF", "padj": 0.001, "log2_enrichment": 1.5},
                        {"name": "UNCHANGED_TF", "padj": 0.001, "log2_enrichment": 1.5},
                    ]
                }
            }
        }
    }
    frame, summary = cm.tf_layer(
        atac,
        {"name": "A_vs_B"},
        rna_df,
        False,
        0.05,
        [],
        expression={"EXPR_NO_DE": 15.0, "ABSENT_TF": 0.0},
    )
    cls_map = dict(zip(frame["tf"], frame["classification"]))
    assert cls_map["COOKS_TF"] == "rna_not_tested"
    assert cls_map["EXPR_NO_DE"] == "rna_not_tested"
    assert cls_map["ABSENT_TF"] == "not_detected_in_rna"
    assert cls_map["UNCHANGED_TF"] == "chromatin_shift_rna_unchanged"
    assert summary["n_chromatin_shift_rna_unchanged"] == 1

    # Adapters: reject NaN/inf/bool and skip non-success bulk contrasts
    bulk_sigs = bulk_rna_evidence(
        {
            "bulk_rna_agent": {
                "findings": {
                    "contrasts": [
                        {
                            "name": "failed_c",
                            "status": "error",
                            "top_genes": [{"symbol": "ERRG", "log2fc": 3.0}],
                        },
                        {
                            "name": "ok_c",
                            "status": "success",
                            "top_genes": [
                                {"symbol": "GOODG", "log2fc": 2.0},
                                {"symbol": "NANG", "log2fc": float("nan")},
                                {"symbol": "INFG", "log2fc": float("inf")},
                                {"symbol": "BOOLG", "log2fc": True},
                            ],
                        },
                    ]
                }
            }
        }
    )
    assert [s.entity for s in bulk_sigs] == ["GOODG"]

    sc_sigs = scrna_evidence(
        {
            "scrna_agent": {
                "findings": {
                    "pseudobulk_de": {
                        "per_group": {
                            "T": {
                                "per_comparison": {
                                    "a_vs_b": {
                                        "status": "success",
                                        "top_genes": [
                                            {"gene": "SC_OK", "log2fc": -1.5},
                                            {"gene": "SC_NAN", "log2fc": float("nan")},
                                            {"gene": "SC_BOOL", "log2fc": False},
                                        ],
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    )
    assert [s.entity for s in sc_sigs] == ["SC_OK"]

    assert chromatin_num(float("nan")) is None
    assert chromatin_num(float("inf")) is None
    assert chromatin_num(True) is None
    assert chromatin_num(1.5) == 1.5
    msigs = motif_signals(
        {
            "padj_max": 0.05,
            "per_group": {
                "up": {
                    "top_motifs": [
                        {"name": "M1", "padj": 0.001, "log2_enrichment": float("nan")},
                    ]
                }
            },
        },
        modality="bulk_ATAC",
        node="ledger://chromatin/motif_enrichment",
        base_caveats=[],
        max_motifs=5,
    )
    assert msigs[0].value is None

    # grounding.py: short target MYC must NOT whitelist arbitrary unmeasured gene MYCN
    myc_sig = EvidenceSignal(
        entity="CDKN1A",
        entity_kind="gene",
        modality="bulk_RNA",
        measure="log2fc",
        audited_node_ref="ledger://bulk/differential_expression",
        value=2.0,
        direction="up",
        context="myc_ko_vs_wt",
    )
    hyp_smuggle = Hypothesis(
        id="smug",
        mechanism="MYC loss may act via unmeasured MYCN to elevate CDKN1A",
        entities=["MYC", "CDKN1A"],
        observation_refs=["ledger://bulk/differential_expression"],
        observed_claims=[{"signal_id": myc_sig.signal_id, "stated_direction": "up"}],
        experiment=DiscriminatingExperiment("MYC KD", "CDKN1A qPCR", "increase", "unchanged"),
    )
    gr = verify_hypothesis_grounding(hyp_smuggle, [myc_sig], None)
    assert gr.grounded is False
    assert "MYCN" in gr.ungrounded_prose_entities

    # proposer_llm.py: reads user_question fallback
    prompt = build_proposer_prompt([myc_sig], {"user_question": "Does MYC KO alter CDKN1A?"}, 2)
    assert "Does MYC KO alter CDKN1A?" in prompt

    # design_agent.py: reads user_question / raw_user_question / summary fallbacks
    da = DesignAgent.__new__(DesignAgent)
    captured: dict[str, str] = {}

    def _capture(prompt_txt, *args, **kwargs):
        captured["prompt"] = prompt_txt
        return {"groups": {"all": ["s1"]}, "confidence": "low", "reasoning": "stub"}

    da.think_structured = _capture
    for key in ("user_question", "raw_user_question", "summary"):
        captured.clear()
        da._propose_groups(
            [{"stem": "s1", "tokens": ["s1"], "paired_files": []}],
            {key: f"Compare treated vs control via {key}"},
        )
        assert f"Compare treated vs control via {key}" in captured["prompt"]


# ── 9. D5-M7 & D5-M8: claim_compiler fdr_strategy, executive negation, mouse h5ad ─


def test_d5m7_d5m8_claim_compiler_executive_negation_and_mouse_h5ad(tmp_path: Path):
    # claim_compiler: per_cluster strategy prefers n_significant_local over n_significant_global
    b_local = NarrativeBlock(
        id="scrna.pseudobulk.T.a_vs_b",
        modality="scRNA",
        analysis="pseudobulk_de",
        block_type="result",
        title="DE",
        status="success",
        confidence="medium",
        claim="T cells had 20 local-FDR DE genes.",
        evidence=[EvidenceItem("local-FDR DE genes", 20, "pseudobulk_de")],
        metrics={
            "fdr_strategy": "per_cluster",
            "n_significant_global": 0,
            "n_significant_local": 20,
        },
    )
    assert classify_claim(b_local).tier == "associative"

    b_global = NarrativeBlock(
        id="scrna.pseudobulk.B.a_vs_b",
        modality="scRNA",
        analysis="pseudobulk_de",
        block_type="result",
        title="DE",
        status="success",
        confidence="medium",
        claim="B cells had 0 global-FDR DE genes.",
        evidence=[EvidenceItem("global-FDR DE genes", 0, "pseudobulk_de")],
        metrics={
            "fdr_strategy": "global",
            "n_significant_global": 0,
            "n_significant_local": 20,
        },
    )
    assert classify_claim(b_global).tier == "descriptive"

    # executive.py: negated causal disclaimers pass; unnegated causal phrases fall back
    class _ExecHost(ReportBuilderMixin):
        def _fallback_executive_summary(self, grouped, intent):
            return "FALLBACK SUMMARY"

        def _summarize_agent_results_for_llm(self, agent_results):
            return "10 genes"

    host = _ExecHost()
    grouped = {"high": [], "medium": [], "low": [], "insufficient": []}

    ok_text = "Observed differences in 10 genes do not prove causation and do not establish causal regulation."
    out_ok, warn_ok = host._govern_executive_summary(ok_text, grouped, {}, {}, [])
    assert out_ok == ok_text
    assert warn_ok == ""

    bad_text = "Observed differences do not prove X; however GENEX drives disease."
    out_bad, warn_bad = host._govern_executive_summary(bad_text, grouped, {}, {}, [])
    assert out_bad == "FALLBACK SUMMARY"
    assert "unlicensed causal language" in warn_bad

    # data_audit_agent.py: mouse h5ad returns ("mm39", "Mus musculus")
    anndata = pytest.importorskip("anndata")
    mouse_h5ad = tmp_path / "mouse_only.h5ad"
    mouse_genes = ["Actb", "Gapdh", "Malat", "Sox", "Gm12345", "Gm23456", "Cebpb", "Gata"]
    ad_m = anndata.AnnData(
        X=np.ones((2, len(mouse_genes)), dtype=np.float32),
        var=pd.DataFrame(index=mouse_genes),
    )
    ad_m.write_h5ad(mouse_h5ad)
    assert DataAuditAgent._infer_h5ad_genome_organism([str(mouse_h5ad)]) == (
        "mm39",
        "Mus musculus",
    )

    human_h5ad = tmp_path / "human_only.h5ad"
    human_genes = ["SAMD11", "ISG15", "TMEM88B", "PRDM16", "MEGF6", "C1QA", "LINC00115"]
    ad_h = anndata.AnnData(
        X=np.ones((2, len(human_genes)), dtype=np.float32),
        var=pd.DataFrame(index=human_genes),
    )
    ad_h.write_h5ad(human_h5ad)
    assert DataAuditAgent._infer_h5ad_genome_organism([str(human_h5ad)]) == (
        "hg38",
        "Homo sapiens",
    )
