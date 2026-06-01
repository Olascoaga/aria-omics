"""P0-6 regression: no filename-fallback design in production.

When the confirmed DesignAgent design failed to apply, BulkRNAAgent silently fell
back to inferring experimental groups from file/column names ("for
compatibility") and ran DE on a guessed design. In production this must STOP; the
filename fallback is allowed only under an explicit ARIA_ALLOW_FILENAME_FALLBACK=1
opt-in, with a loud warning.
"""

import pytest


def _agent(monkeypatch, capture):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: capture["findings"].append((a, k))
    agent.memory = type("M", (), {"store_decision": staticmethod(lambda **k: None)})()
    # Spy: the production gate must NOT reach filename inference.
    agent._discover_groups = lambda *a, **k: capture.__setitem__("discovered", True) or ([], {})
    return BulkRNAAgent, agent


def _ctx_with_unmappable_design(tmp_path):
    counts = tmp_path / "counts.tsv"
    # Only 2 columns, but the design names 4 distinct samples -> neither name nor
    # positional mapping can succeed, so _apply_design raises.
    counts.write_text("gene_id\tX1\tX2\nGENE_1\t10\t80\n", encoding="utf-8")
    design = {"groups": {"A": ["S1", "S2"], "B": ["S3", "S4"]},
              "main_factor": "condition"}
    return {
        "exp_context": {"modalities": {"bulk_RNA": [str(counts)]}, "design": design},
        "biological_intent": {"summary": "compare A vs B"},
    }


def test_filename_fallback_allowed_reads_env(monkeypatch):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    monkeypatch.delenv("ARIA_ALLOW_FILENAME_FALLBACK", raising=False)
    assert BulkRNAAgent._filename_fallback_allowed() is False
    for truthy in ("1", "true", "YES", "on"):
        monkeypatch.setenv("ARIA_ALLOW_FILENAME_FALLBACK", truthy)
        assert BulkRNAAgent._filename_fallback_allowed() is True
    monkeypatch.setenv("ARIA_ALLOW_FILENAME_FALLBACK", "0")
    assert BulkRNAAgent._filename_fallback_allowed() is False


def test_failed_confirmed_design_stops_in_production(tmp_path, monkeypatch):
    monkeypatch.delenv("ARIA_ALLOW_FILENAME_FALLBACK", raising=False)
    capture = {"findings": [], "discovered": False}
    _, agent = _agent(monkeypatch, capture)

    result = agent.run("exp-p0-6", _ctx_with_unmappable_design(tmp_path))

    assert result["status"] == "failed"
    assert result["reason"] == "design_application_failed"
    # The analysis stopped — it never guessed groups from file/column names.
    assert capture["discovered"] is False


def test_flag_opt_in_allows_filename_fallback(tmp_path, monkeypatch):
    monkeypatch.setenv("ARIA_ALLOW_FILENAME_FALLBACK", "1")
    capture = {"findings": [], "discovered": False}
    _, agent = _agent(monkeypatch, capture)

    # With the opt-in set, the failed confirmed design falls through to filename
    # inference (which our spy stubs to "no groups", so the run still fails — but
    # crucially it WAS allowed to try).
    result = agent.run("exp-p0-6", _ctx_with_unmappable_design(tmp_path))

    assert capture["discovered"] is True
    assert result["status"] == "failed"  # spy returned no groups
