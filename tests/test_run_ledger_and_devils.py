"""Stage 4 P-LEDGER + R2/P-DEVIL (audit 2026-05-29):

- P-LEDGER: a deterministic planned-vs-run manifest flags any analysis the plan
  called for that did not run (the PBMC thin-report dispatch gap).
- P-DEVIL: a deterministic devil's advocate challenges every associative-or-
  stronger claim with standard technical confounders and records which were not
  addressed, reusing the Claim Compiler tiers (no LLM rounds).
"""

from aria.agents.narrative.types import NarrativeBlock, EvidenceItem
from aria.agents.narrative.run_ledger import build_run_ledger
from aria.agents.narrative.devils_advocate import build_devils_advocate


# ── P-LEDGER ─────────────────────────────────────────────────────────────────

def test_ledger_flags_planned_but_not_run():
    # Plan recommended pseudobulk DE + pathways, but only QC actually ran — the
    # PBMC thin-report shape.
    exp_ctx = {"design_intelligence": {
        "recommended": ["Donor-level pseudobulk DESeq2 between conditions."],
        "optional": ["Pathway/ORA enrichment on significant pseudobulk DE genes."],
    }}
    agent_results = {"scrna_agent": {"findings": {
        "qc": {"status": "success", "n_cells": 100},
    }}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    by_key = {e["analysis"]: e for e in ledger["entries"]}
    assert by_key["pseudobulk_de"]["planned"] is True
    assert by_key["pseudobulk_de"]["status"] == "not_run"
    assert by_key["pseudobulk_de"]["divergence"] is True
    assert by_key["pathway_enrichment"]["divergence"] is True
    assert by_key["qc"]["status"] == "ran"
    assert ledger["n_divergences"] >= 2


def test_ledger_no_divergence_when_plan_executed():
    exp_ctx = {"design_intelligence": {
        "recommended": ["Donor-level pseudobulk DESeq2 between conditions."],
        "optional": [],
    }}
    agent_results = {"scrna_agent": {"findings": {
        "qc": {"status": "success", "n_cells": 100},
        "pseudobulk_de": {"per_group": {"A": {}}},
    }}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    by_key = {e["analysis"]: e for e in ledger["entries"]}
    assert by_key["pseudobulk_de"]["status"] == "ran"
    assert by_key["pseudobulk_de"]["divergence"] is False
    assert ledger["n_divergences"] == 0


def test_ledger_records_skip_reason():
    exp_ctx = {"design_intelligence": {
        "recommended": ["LIANA cell-cell communication."], "optional": []}}
    agent_results = {"scrna_agent": {"findings": {
        "cell_communication": {"status": "skipped", "reason": "no_groupby_column"},
    }}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    cc = next(e for e in ledger["entries"] if e["analysis"] == "cell_communication")
    assert cc["status"] == "skipped"
    assert cc["reason"] == "no_groupby_column"
    assert cc["divergence"] is True   # planned but did not run


# ── P-DEVIL ──────────────────────────────────────────────────────────────────

def _assoc_block():
    b = NarrativeBlock(
        id="scrna.pseudobulk.Tcell.treat_vs_ctrl",
        modality="scRNA", analysis="pseudobulk_de", block_type="result",
        title="DE", status="success", confidence="medium",
        claim="Tcell shows differential expression between treat and ctrl.",
        evidence=[EvidenceItem(label="n_sig", value=42)],
        metrics={"corrected_for_composition": False, "low_power_warning": False},
    )
    b.metadata["claim"] = {"tier": "associative"}
    return b


def test_devils_advocate_flags_unaddressed_confounders():
    block = _assoc_block()
    # No integration, no ambient correction, no doublet signal, DA not run.
    agent_results = {"scrna_agent": {"findings": {"qc": {"status": "success"}}}}
    manifest = build_devils_advocate([block], agent_results, {})
    assert len(manifest) == 1
    entry = manifest[0]
    assert entry["tier"] == "associative"
    alts = {c["alternative"] for c in entry["challenges"]}
    assert {"batch_effect", "ambient_rna", "doublets", "composition_shift",
            "low_replication"} <= alts
    # ambient is never addressed (ARIA has no SoupX/decontX)
    assert "ambient_rna" in entry["unaddressed"]
    # an info caveat naming the unaddressed alternatives was attached
    assert any("not ruled out" in c.text for c in block.caveats)
    assert block.metadata["devils_advocate"]["unaddressed"]


def test_devils_advocate_addresses_confounders_with_evidence():
    block = _assoc_block()
    block.metrics["corrected_for_composition"] = True
    agent_results = {"scrna_agent": {"findings": {
        "qc": {"status": "success", "scrublet_doublet_scores": [0.1, 0.2]},
        "integration": {"status": "success", "method": "harmony"},
        "integration_qc": {"issues": []},
    }}}
    manifest = build_devils_advocate([block], agent_results, {})
    challenges = {c["alternative"]: c["addressed"] for c in manifest[0]["challenges"]}
    assert challenges["batch_effect"] is True       # integration ran, clean
    assert challenges["doublets"] is True           # scrublet signal present
    assert challenges["composition_shift"] is True  # block corrected
    assert challenges["ambient_rna"] is False       # still not in the pipeline


def test_devils_advocate_ambient_detector_is_not_a_correction():
    # P1-4: the ambient-RNA *detector* finding contains the word "ambient" but
    # is NOT an applied correction — it must not flip ambient_rna to addressed.
    block = _assoc_block()
    agent_results = {"scrna_agent": {"findings": {
        "qc": {"status": "success",
               "ambient_correction": {"ran": False, "reason": "not_requested"}},
        "ambient_qc": {"status": "warnings", "issues": [
            {"check": "possible_ambient_contamination",
             "message": "ambient (soup) RNA signature", "severity": "warning",
             "recommendation": "consider decontamination"}]},
    }}}
    manifest = build_devils_advocate([block], agent_results, {})
    challenges = {c["alternative"]: c["addressed"] for c in manifest[0]["challenges"]}
    assert challenges["ambient_rna"] is False


def test_devils_advocate_ambient_addressed_when_correction_ran():
    # P1-4: a real applied decontamination (ran=True) DOES address ambient_rna.
    block = _assoc_block()
    agent_results = {"scrna_agent": {"findings": {
        "qc": {"status": "success",
               "ambient_correction": {"ran": True, "method": "decontx"}},
    }}}
    manifest = build_devils_advocate([block], agent_results, {})
    challenges = {c["alternative"]: c["addressed"] for c in manifest[0]["challenges"]}
    assert challenges["ambient_rna"] is True


def test_devils_advocate_skips_descriptive_claims():
    qc = NarrativeBlock(
        id="scrna.qc", modality="scRNA", analysis="qc", block_type="qc",
        title="QC", status="success", confidence="high",
        claim="242k cells passed QC.",
        evidence=[EvidenceItem(label="n", value=242000)],
    )
    qc.metadata["claim"] = {"tier": "descriptive"}
    manifest = build_devils_advocate([qc], {"scrna_agent": {}}, {})
    assert manifest == []        # descriptive claims are not challenged


def test_devils_advocate_is_idempotent():
    block = _assoc_block()
    agent_results = {"scrna_agent": {"findings": {"qc": {"status": "success"}}}}
    build_devils_advocate([block], agent_results, {})
    n_caveats = len(block.caveats)
    build_devils_advocate([block], agent_results, {})   # second pass
    assert len(block.caveats) == n_caveats              # no duplicate caveat
