"""End-to-end coverage for the W-LEDGER claim/ledger check on the REAL render
path (`NarrativeAgent._render_html_report`).

Two things are exercised: (1) the bulk GSEA-only case that used to make the
ledger read pathway_enrichment as not-run while the narrator still built an
associative GSEA block, and (2) C3 fail-closed publication — an unexpected
post-compilation linkage mismatch aborts before report publication. Both drive
the actual narrators through `_render_html_report` (not synthetic blocks).
"""

import json

import pytest

litellm = pytest.importorskip("litellm")  # NarrativeAgent import needs litellm

from aria.agents.narrative_agent import NarrativeAgent


def _agent(tmp_path):
    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    return agent


def _render(agent, tmp_path, agent_results, name="r"):
    return agent._render_html_report(
        experiment_id="exp_e2e",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "e2e"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="methods",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / name,
    )


def _gsea_plots(tmp_path, n=2):
    # The narrator turns contrast plots into figure refs that the renderer
    # validates as existing files, so write real (empty) plot files.
    paths = []
    for i in range(n):
        p = tmp_path / f"running_sum_{i}.png"
        p.write_bytes(b"\x89PNG\r\n")
        paths.append(str(p))
    return paths


def test_bulk_gsea_only_report_renders_and_pathway_node_is_ran(tmp_path):
    # A bulk contrast with GSEA running-sum plots only (no gsea_table, no ORA
    # pathways) builds a DE block AND an associative GSEA block. The report must
    # render, and the ledger must mark pathway_enrichment "ran" (the detection
    # fix), so there is NO spurious claim/ledger violation.
    agent_results = {"bulk_rna_agent": {"findings": {
        "sample_qc": {"n_samples": 6},
        "contrasts": [{
            "name": "treat_vs_ctrl", "status": "success",
            "n_significant": 5, "n_upregulated": 3, "n_downregulated": 2,
            "plots": {"gsea_running_sums": _gsea_plots(tmp_path)},
        }],
    }}}
    report = _render(_agent(tmp_path), tmp_path, agent_results)
    assert report.exists()                      # rendered, no crash

    methodology = json.loads((report.parent / "methodology.json").read_text())
    ledger = methodology["run_ledger"]
    by_node = {e["node_id"]: e for e in ledger["entries"]}
    assert by_node["ledger://bulk/pathway_enrichment"]["status"] == "ran"
    assert by_node["ledger://bulk/differential_expression"]["status"] == "ran"
    # The record-only verification ran and found no violation.
    verification = ledger.get("claim_ledger_verification", {})
    assert verification.get("n_violations", 0) == 0


def test_render_fails_closed_when_ledger_check_finds_a_violation(tmp_path, monkeypatch):
    # Force a post-compilation invariant failure: C3 must abort before writing a
    # public report, rather than publishing the invalid claim with a caveat.
    import aria.agents.narrative.run_ledger as rl

    def _fake_verify(blocks, run_ledger, strict=False):
        assert strict is True
        raise rl.LedgerLinkageError("forced C3 mismatch")

    monkeypatch.setattr(rl, "verify_blocks_against_ledger", _fake_verify)

    agent_results = {"bulk_rna_agent": {"findings": {
        "sample_qc": {"n_samples": 6},
        "contrasts": [{
            "name": "treat_vs_ctrl", "status": "success",
            "n_significant": 5, "n_upregulated": 3, "n_downregulated": 2,
            "plots": {"gsea_running_sums": _gsea_plots(tmp_path, 1)},
        }],
    }}}
    with pytest.raises(rl.LedgerLinkageError, match="forced C3 mismatch"):
        _render(_agent(tmp_path), tmp_path, agent_results, name="r2")
    assert not (tmp_path / "r2" / "report.html").exists()
