"""FAIL->PASS regression tests for D5 HIGH and MEDIUM findings (D5-003..D5-025)."""

from __future__ import annotations

import json
import sqlite3
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from aria.agents._narrative_chromatin import generate_figures as chrom_generate_figures
from aria.agents.audit_agent import AuditAgent
from aria.agents.bulk_rna_agent import BulkRNAAgent
from aria.agents.narrative.hypothesis.adapters.scrna import scrna_evidence
from aria.agents.narrative.hypothesis.caveats import AMBIENT, BATCH, DOUBLET
from aria.agents.narrative.ledger_export import (
    verify_reproducible_capsule,
    write_reproducible_capsule,
    write_ro_crate,
)
from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
from aria.agents.narrative.render_blocks import _claim_tier_badge
from aria.agents.narrative.run_ledger import build_run_ledger
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.agents.narrative_agent import NarrativeAgent
from aria.agents.scrna.agent import scRNAAgent
from aria.agents.scrna.de import DEMixin
from aria.memory.memory import ARIAMemory
from aria.runtime import experiment_view as ev


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


def test_d5_003_devils_advocate_handles_nan_and_surfaces_pass_failure(tmp_path, monkeypatch):
    agent = _make_agent(tmp_path)
    agent_results = {
        "bulk_rna_agent": _bulk_ok(),
        "integration_agent": {
            "status": "done",
            "findings": {"eregulons": {"status": "skipped", "n_depth_confounded": float("nan")}},
        },
    }
    report = agent._render_html_report(
        experiment_id="d5_003",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report_nan",
    )
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    da = meth["devils_advocate"]
    assert isinstance(da, list) and any(
        e.get("claim_id") == "bulk.contrast.treat_vs_ctrl" for e in da
    ), f"Associative claim missing from devils_advocate on NaN eregulons count: {da!r}"

    # When build_devils_advocate genuinely fails, surface error in HTML + methodology.json
    from aria.agents.narrative import devils_advocate as da_mod

    monkeypatch.setattr(
        da_mod,
        "build_devils_advocate",
        lambda *_a, **_kw: (_ for _ in ()).throw(RuntimeError("da pass boom")),
    )
    report2 = agent._render_html_report(
        experiment_id="d5_003_fail",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"bulk_rna_agent": _bulk_ok()},
        report_dir=tmp_path / "report_fail",
    )
    html2 = report2.read_text(encoding="utf-8")
    meth2 = json.loads((report2.parent / "methodology.json").read_text(encoding="utf-8"))
    assert isinstance(meth2["devils_advocate"], dict)
    assert meth2["devils_advocate"].get("status") == "error"
    assert "da pass boom" in html2


def test_d5_004_falsification_certificates_decoupled_from_robustness_and_surfaces_error(
    tmp_path, monkeypatch
):
    agent = _make_agent(tmp_path)
    agent_results = {
        "bulk_rna_agent": _bulk_ok(),
        "scrna_agent": {
            "status": "done",
            "findings": {
                "pseudobulk_de": {
                    "per_group": {
                        "T_cell": {
                            "per_comparison": {
                                "treat_vs_ctrl": {
                                    "status": "success",
                                    "n_significant_local": float("nan"),
                                    "n_significant": 5,
                                }
                            }
                        }
                    }
                }
            },
        },
    }
    report = agent._render_html_report(
        experiment_id="d5_004",
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
    html = report.read_text(encoding="utf-8")
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    assert meth["falsification_certificates"]["n_certificates"] > 0
    assert "Falsification Certificates" in html

    # When build_falsification_certificates itself fails, HTML must display error banner.
    from aria.agents.narrative import falsification_certificates as fc_mod

    monkeypatch.setattr(
        fc_mod,
        "build_falsification_certificates",
        lambda *_a, **_kw: (_ for _ in ()).throw(RuntimeError("fc pass boom")),
    )
    report2 = agent._render_html_report(
        experiment_id="d5_004_fail",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"bulk_rna_agent": _bulk_ok()},
        report_dir=tmp_path / "report_fail",
    )
    html2 = report2.read_text(encoding="utf-8")
    meth2 = json.loads((report2.parent / "methodology.json").read_text(encoding="utf-8"))
    assert meth2["falsification_certificates"].get("status") == "error"
    assert "Falsification Certificates" in html2 and "fc pass boom" in html2


def test_d5_007_scrna_qc_guardrails_record_error_findings_on_exception(tmp_path):
    h5_path = tmp_path / "qc_out.h5ad"
    ad.AnnData(
        X=np.ones((4, 3), dtype=np.float32),
        obs=pd.DataFrame(
            {"cell_type": ["T_cell", "T_cell", "B_cell", "B_cell"]},
            index=["c0", "c1", "c2", "c3"],
        ),
        var=pd.DataFrame(index=["CD3D", "MS4A1", "GAPDH"]),
    ).write_h5ad(h5_path)

    class _StubEnv:
        def run_in_stack(self, stack, script_path, params, **kw):
            if "rna_qc.py" in script_path:
                return {
                    "status": "success",
                    "output_path": str(h5_path),
                    "n_cells_before": 500,
                    "n_cells_after": 480,
                    "pct_removed": 4.0,
                    "batch_col": "batch",
                }
            if "rna_integration.py" in script_path:
                return {
                    "status": "success",
                    "output_path": str(h5_path),
                    "silhouette_before": 0.40,
                    "silhouette_after": -0.15,
                }
            if "rna_clustering.py" in script_path:
                return {
                    "status": "success",
                    "output_path": str(h5_path),
                    "n_clusters": 2,
                    "groupby": "cell_type",
                    "predef_clusters": False,
                    "silhouette": -0.10,
                    "cluster_sizes": [240, 240],
                    "top_markers": {"0": ["CD3D"], "1": ["MS4A1"]},
                }
            if "rna_celltypist.py" in script_path:
                return {
                    "status": "error",
                    "error_type": "InvalidExpressionMatrix",
                    "details": "scaled matrix",
                }
            return {"status": "skipped", "reason": "test_stub"}

    agent = scRNAAgent(ARIAMemory(":memory:"), llm=None)
    agent.env = _StubEnv()
    res = agent.run(
        "exp_d5_007",
        {
            "exp_context": {
                "modalities": {"scRNA": [str(h5_path)]},
                "design": {"batch_covariate": "batch"},
                "inferred_design": {
                    "groupby_col": "cell_type",
                    "obs_columns": {"sample.h5ad": 42},
                },
            },
            "biological_intent": {},
        },
    )
    findings = res["findings"]
    assert res["status"] == "done"
    for key in ("batch_qc", "annotation_qc"):
        assert key in findings, f"scRNAAgent omitted failed guardrail {key}"
        assert findings[key].get("status") == "error"
        assert findings[key].get("issues")


def test_d5_009_methodology_json_surfaces_compile_public_claims_error_and_redacts_blocks(
    tmp_path,
):
    agent = _make_agent(tmp_path)
    bad_block = NarrativeBlock(
        id="bulk.contrast.treat_vs_ctrl",
        modality="bulk RNA-seq",
        analysis="contrast",
        block_type="result",
        title="treat vs ctrl",
        status="success",
        confidence="high",
        claim="5 DE genes identified.",
        evidence=[EvidenceItem("DE genes", 5, "deseq2")],
        metadata=None,
    )
    meth = agent._build_methodology_json(
        provenance={"aria_version": "4.9.4"},
        exp_ctx={},
        agent_results={},
        decisions=[],
        narrative_blocks=[bad_block],
        run_ledger={"entries": [], "divergences": [], "n_divergences": 0},
    )
    cc = meth["claim_compilation"]
    assert cc.get("status") == "error"
    assert cc.get("error")
    assert len(cc.get("withheld", [])) == cc.get("n_withheld") == 1
    assert meth["narrative_blocks"] == []


def test_d5_010_withheld_claim_disclosed_in_html_and_compiler_internal_errors_tagged(
    tmp_path, monkeypatch
):
    agent = _make_agent(tmp_path)
    ar = _bulk_ok()
    ar["findings"]["contrasts"][0]["plots"] = {
        "volcano": str(tmp_path / "missing_volcano.png")
    }
    report = agent._render_html_report(
        experiment_id="d5_010",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={},
        executive_summary="ARIA completed the analysis.",
        findings_sections={},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="Methods.",
        decisions=[],
        agent_results={"bulk_rna_agent": ar},
        report_dir=tmp_path / "report",
    )
    html = report.read_text(encoding="utf-8")
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    assert meth["claim_compilation"]["n_withheld"] == 1
    assert "withheld" in html.lower()

    # Unexpected internal exception in compile_public_claims must be tagged distinctly.
    from aria.agents.narrative import claim_compiler as cc_mod

    monkeypatch.setattr(
        cc_mod,
        "compile_claim_manifest",
        lambda _b: (_ for _ in ()).throw(RuntimeError("internal compiler bug")),
    )
    blk = NarrativeBlock(
        id="bulk.contrast.c1",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="c1",
        status="success",
        confidence="medium",
        claim="5 DE genes identified.",
        evidence=[EvidenceItem("DE genes", 5, "deseq2")],
    )
    comp = cc_mod.compile_public_claims([blk], {}, run_ledger={"entries": []})
    assert len(comp.withheld) == 1
    assert "compiler_internal_error" in comp.withheld[0]["reason"]


def test_d5_011_synthesis_blocks_not_silently_dropped_before_claim_compiler():
    bad_syn = NarrativeBlock(
        id="integration.bulk_convergence",
        modality="integrated",
        analysis="convergent_evidence",
        block_type="result",
        title="Cross-contrast convergence",
        status="success",
        confidence="medium",
        claim="2 contrasts converge.",
        evidence=["malformed_evidence_item"],
    )
    kept = NarrativeAgent._drop_unverifiable_synthesis_blocks([bad_syn])
    assert len(kept) == 1


def test_d5_012_scrna_hypothesis_adapter_fails_conservative_on_signal_error(monkeypatch):
    agent_results = {
        "scrna_agent": {
            "status": "done",
            "findings": {
                "qc": {"status": "success"},
                "pseudobulk_de": {
                    "per_group": {
                        "T_cell": {
                            "per_comparison": {
                                "treat_vs_ctrl": {
                                    "status": "success",
                                    "top_genes": [
                                        {"gene": "TP53", "log2fc": 2.0, "padj": 1e-5}
                                    ],
                                }
                            }
                        }
                    }
                },
            },
        },
        "integration_agent": {
            "status": "done",
            "findings": {"eregulons": {"n_depth_confounded": float("nan")}},
        },
    }
    signals = scrna_evidence(agent_results, run_ledger=None, exp_ctx={})
    assert len(signals) == 1
    assert {AMBIENT, DOUBLET, BATCH}.issubset(set(signals[0].caveats_inherited))

    # Even if _run_signals raises directly, _scrna_caveats must fail conservative (all 3 caveats).
    from aria.agents.narrative.hypothesis.adapters import scrna as scrna_adapter

    monkeypatch.setattr(
        scrna_adapter,
        "_run_signals",
        lambda *_a, **_kw: (_ for _ in ()).throw(RuntimeError("signals crash")),
    )
    signals2 = scrna_evidence(agent_results, run_ledger=None, exp_ctx={})
    assert {AMBIENT, DOUBLET, BATCH}.issubset(set(signals2[0].caveats_inherited))


def test_d5_013_audit_agent_surfaces_check_exceptions_instead_of_clean_status(
    tmp_path, monkeypatch
):
    tsv = tmp_path / "counts.tsv"
    tsv.write_text("gene\tc1\tc2\tt1\tt2\nG1\t10\t12\t50\t52\n", encoding="utf-8")
    auditor = AuditAgent(ARIAMemory(":memory:"), llm=None)
    monkeypatch.setattr(
        auditor,
        "_check_pca_batch_dominance",
        lambda *_a, **_kw: (_ for _ in ()).throw(ValueError("simulated PCA crash")),
    )
    out = auditor.run_audit(
        {
            "files": [str(tsv)],
            "modalities": {"bulk_RNA": [str(tsv)]},
            "design": {
                "batch_factor": "batch",
                "groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]},
                "batch_map": {"c1": "B1", "c2": "B2", "t1": "B1", "t2": "B2"},
            },
        },
        experiment_id="d5_013",
    )
    checks = [f["check"] for f in out["findings"]]
    assert out["status"] != "clean"
    assert any("pca" in c and "failed" in c for c in checks)


def test_d5_014_scrna_expressed_gene_background_propagates_corrupt_h5ad_error(tmp_path):
    bad_h5ad = tmp_path / "corrupt.h5ad"
    bad_h5ad.write_bytes(b"NOT_AN_HDF5_FILE")
    with pytest.raises((OSError, ValueError, RuntimeError)):
        DEMixin._expressed_gene_background(str(bad_h5ad))


def test_d5_015_bulk_rna_metadata_read_failure_fails_closed_for_confirmed_covariates(
    tmp_path,
):
    bad_gz = tmp_path / "sample_sheet.tsv.gz"
    bad_gz.write_bytes(b"CORRUPT_GZIP_BYTES")
    with pytest.raises((OSError, ValueError, RuntimeError)):
        BulkRNAAgent._write_design_metadata(
            {"c1": "ctrl", "c2": "ctrl", "t1": "treat", "t2": "treat"},
            "condition",
            tmp_path / "out",
            covariates=["batch"],
            exp_ctx={"metadata_file": str(bad_gz)},
            design={"batch_covariate": "batch"},
        )


def test_d5_016_bulk_gsea_csv_parse_failure_emits_warning_caveat(tmp_path):
    gsea_csv = tmp_path / "gsea_report.csv"
    gsea_csv.write_bytes(b"Term,NES,FDR q-val\nIFN-\xb5 signaling,2.45,0.01\n")
    blocks = BulkRnaNarrator()._gsea_blocks_for_contrast(
        {"name": "treat_vs_ctrl", "plots": {"gsea_table": str(gsea_csv)}}
    )
    assert len(blocks) == 1
    blk = blocks[0]
    warn_cavs = [c for c in blk.caveats if c.severity in ("warning", "blocking")]
    assert blk.status != "success" or len(warn_cavs) > 0


def test_d5_017_executive_summary_governance_handles_none_log2fc():
    agent = NarrativeAgent.__new__(NarrativeAgent)
    ar = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "contrasts": [
                    {
                        "name": "treat_vs_ctrl",
                        "status": "success",
                        "n_significant": 42,
                        "n_upregulated": 25,
                        "n_downregulated": 17,
                        "top_genes": [
                            {"gene": "TP53", "symbol": "TP53", "log2fc": None, "padj": 1e-5}
                        ],
                    }
                ]
            },
        }
    }
    valid_summary = "Bulk RNA-seq identified 42 significant DE genes between treat and ctrl."
    text, warning = agent._govern_executive_summary(
        valid_summary,
        {"high": [], "medium": [], "low": [], "insufficient": []},
        {},
        ar,
        [],
    )
    assert text == valid_summary and "unsupported numeric claims" not in warning


def test_d5_018_run_ledger_records_triangulation_graph_exception_as_error():
    ar = {
        "integration_agent": {
            "status": "done",
            "findings": {
                "peak_to_gene": {
                    "status": "success",
                    "top_links": [
                        {"peak": "chr1:1-100", "gene": "TP53", "correlation": 0.8, "padj": 1e-4}
                    ],
                }
            },
        },
        "scrna_agent": {
            "status": "done",
            "findings": {
                "pseudobulk_de": {
                    "per_group": {"T_cell": {"per_comparison": ["malformed_list"]}}
                }
            },
        },
    }
    ledger = build_run_ledger({}, ar)
    node = next(e for e in ledger["entries"] if e["analysis"] == "cross_modal_triangulation")
    assert node["status"] == "error" and node["reason"]


def test_d5_019_corrupt_memory_db_fails_closed_in_get_decisions_log(tmp_path):
    corrupt_db = tmp_path / "corrupt.sqlite"
    mem = ARIAMemory(str(corrupt_db))
    mem._conn.close()
    corrupt_db.write_bytes(b"NOT_A_VALID_SQLITE_FILE")
    mem._conn = sqlite3.connect(str(corrupt_db), check_same_thread=False)
    agent = _make_agent(tmp_path, memory=mem)
    with pytest.raises((sqlite3.DatabaseError, RuntimeError)):
        agent._get_decisions_log("exp_d5_019")


def test_d5_020_speculative_section_surfaces_error_card_and_manifest(tmp_path):
    agent = _make_agent(tmp_path)
    out_dir = tmp_path / "spec_out"
    out_dir.mkdir()
    html = agent._build_speculative_section_html(
        {"bulk_rna_agent": _bulk_ok()},
        {"entries": [], "claim_ledger_verification": {"n_violations": "bad"}},
        {"enable_hypotheses": True},
        narrative_blocks=[],
        report_dir=out_dir,
    )
    assert html != "" and "error" in html.lower()
    spec_json = out_dir / "speculative_hypotheses.json"
    assert spec_json.exists()
    data = json.loads(spec_json.read_text(encoding="utf-8"))
    assert data.get("status") == "error"


def test_d5_021_verify_reproducible_capsule_does_not_pass_when_env_unverified(
    tmp_path, monkeypatch
):
    cap_dir = tmp_path / "report_cap"
    cap_dir.mkdir()
    (cap_dir / "report.html").write_text("<html>ok</html>", encoding="utf-8")
    mem = ARIAMemory(str(tmp_path / "lab.sqlite"))
    mem.create_wing("exp_cap_1", "Test", "Homo sapiens", "hg38")
    mem.export_experiment_snapshot("exp_cap_1", cap_dir / "memory_snapshot.sqlite")
    (cap_dir / "methodology.json").write_text(
        json.dumps(
            {
                "provenance": {
                    "aria_version": "4.9.4",
                    "git_commit": "abc123",
                    "environment": {},
                },
                "inputs": [],
                "claims": [],
                "run_ledger": {"entries": [], "divergences": [], "n_divergences": 0},
            }
        ),
        encoding="utf-8",
    )
    write_ro_crate(cap_dir)
    cap_zip = write_reproducible_capsule(cap_dir)

    from aria import version as ver_mod

    monkeypatch.setattr(
        ver_mod,
        "collect_environment_metadata",
        lambda *_a, **_kw: (_ for _ in ()).throw(RuntimeError("env probe failed")),
    )
    res = verify_reproducible_capsule(cap_zip)
    assert res["environment"]["status"] == "unverified"
    assert res["status"] == "fail"


def test_d5_022_experiment_view_surfaces_corrupt_methodology_and_errors(
    tmp_path, monkeypatch
):
    rep_dir = tmp_path / "rep"
    rep_dir.mkdir()
    (rep_dir / "report.html").write_text("<html></html>", encoding="utf-8")
    (rep_dir / "methodology.json").write_text("{CORRUPT_JSON", encoding="utf-8")

    artifacts = ev._artifact_views(str(rep_dir / "report.html"))
    meth_art = next(a for a in artifacts if a.category == "methodology")
    assert meth_art.status == "error"

    from aria.agents.narrative import run_ledger as rl_mod

    monkeypatch.setattr(
        rl_mod,
        "build_run_ledger",
        lambda *_a, **_kw: (_ for _ in ()).throw(RuntimeError("ledger view boom")),
    )
    sess = SimpleNamespace(
        exp_context={"organism": "Homo sapiens"},
        agent_results={"bulk_rna_agent": _bulk_ok()},
    )
    ledger_views = ev._ledger_views(sess)
    assert len(ledger_views) >= 1 and ledger_views[0].status == "error"

    # Privacy view must NOT fail open with enabled=False on air_gapped_enabled() error.
    from aria.utils import privacy as priv_mod

    monkeypatch.setattr(
        priv_mod,
        "air_gapped_enabled",
        lambda: (_ for _ in ()).throw(RuntimeError("privacy check boom")),
    )
    pview = ev._resource_privacy_view({})
    assert pview.status == "error"


def test_d5_023_claim_linkage_error_recorded_in_run_ledger(tmp_path, monkeypatch):
    from aria.agents.narrative import run_ledger as rl_mod

    orig_link = rl_mod.link_claims_to_ledger
    calls = {"n": 0}

    def _flaky_link(claims, run_ledger):
        calls["n"] += 1
        # compile_public_claims passes 1-element list; render/methodology passes full list
        if len(claims) >= 2:
            raise RuntimeError("linkage schema drift")
        return orig_link(claims, run_ledger)

    monkeypatch.setattr(rl_mod, "link_claims_to_ledger", _flaky_link)
    agent = _make_agent(tmp_path)
    report = agent._render_html_report(
        experiment_id="d5_023",
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
    meth = json.loads((report.parent / "methodology.json").read_text(encoding="utf-8"))
    cl = meth["run_ledger"].get("claim_linkage") or {}
    assert cl.get("status") == "error" and "linkage schema drift" in str(cl.get("error", ""))


def test_d5_024_chromatin_figure_error_populates_w_warn_errors_key(tmp_path):
    bad_csv = tmp_path / "corrupt_da.csv"
    bad_csv.write_bytes(b"\xff\xfeCORRUPT_BINARY")
    live = {
        "differential_accessibility": {
            "status": "success",
            "data_type": "bulk_ATAC",
            "full_results_csv": str(bad_csv),
            "comparisons": [
                {"test": "treat", "reference": "ctrl", "full_results_csv": str(bad_csv)}
            ],
        }
    }
    chrom_generate_figures(live, None, tmp_path / "figs", env_manager=None)
    assert (live.get("figures") or {}).get("_errors"), (
        f"Figure generation failed on corrupt DA CSV yet figures['_errors'] was not populated: "
        f"figures={live.get('figures')!r}"
    )


def test_d5_025_claim_tier_badge_renders_unclassified_badge_on_error(monkeypatch):
    from aria.agents.narrative import claim_compiler as cc_mod

    monkeypatch.setattr(
        cc_mod,
        "classify_claim",
        lambda _b: (_ for _ in ()).throw(RuntimeError("classifier boom")),
    )
    blk = NarrativeBlock(
        id="bulk.contrast.c1",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="c1",
        status="success",
        confidence="medium",
        claim="5 DE genes identified.",
        evidence=[EvidenceItem("DE genes", 5, "deseq2")],
    )
    badge = _claim_tier_badge(blk)
    assert "tier-unclassified" in badge or "unclassified" in badge.lower()
