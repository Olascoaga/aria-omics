"""P0-7 regression: confirmed analysis thresholds propagate end-to-end.

The orchestrator records the user-confirmed CP3 thresholds in
exp_context["global_padj"]/["global_lfc"], and the per-cluster DE path honored
them, but scRNAAgent._run_pseudobulk dispatched HARDCODED padj_max=0.05 /
lfc_min=0.5 to rna_pseudobulk_de — silently overriding the confirmed thresholds.

The fix introduces an additive `AnalysisThresholds` value object resolved from the
experiment context and fed into the pseudobulk params, so no loose threshold
literal sits in the dispatch.
"""

import pytest

from aria.utils.thresholds import AnalysisThresholds


# ── The value object ─────────────────────────────────────────────────────────

def test_from_exp_context_reads_confirmed_thresholds():
    thr = AnalysisThresholds.from_exp_context(
        {"global_padj": 0.01, "global_lfc": 1.5}, modality="scRNA")
    assert thr.padj == 0.01
    assert thr.log2fc == 1.5
    assert thr.modality == "scRNA"


def test_from_exp_context_falls_back_to_defaults():
    thr = AnalysisThresholds.from_exp_context({}, modality="scRNA")
    assert thr.padj == 0.05
    assert thr.log2fc == 0.5          # scRNA pseudobulk default preserved
    assert thr.min_cells == 10
    assert thr.min_replicates == 3


def test_min_replicates_override_is_honored():
    thr = AnalysisThresholds.from_exp_context({}, min_replicates=2)
    assert thr.min_replicates == 2


def test_as_pseudobulk_params_shape():
    thr = AnalysisThresholds(padj=0.01, log2fc=1.5, min_cells=10,
                             min_replicates=2)
    assert thr.as_pseudobulk_params() == {
        "padj_max": 0.01,
        "lfc_min": 1.5,
        "min_cells_per_pseudosample": 10,
        "min_replicates_per_condition": 2,
    }


# ── End-to-end wiring: the agent dispatches the confirmed thresholds ─────────

def test_pseudobulk_dispatch_uses_confirmed_thresholds_not_hardcoded(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    agent = scRNAAgent.__new__(scRNAAgent)
    agent._workspace = lambda *a, **k: tmp_path
    agent._log_decision = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None

    captured = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            if script_path.endswith("rna_diff_abundance.py"):
                return {"status": "success", "any_significant": False}
            if script_path.endswith("rna_pseudobulk_de.py"):
                captured.update(params)
                # short-circuit after the dispatch we care about
                return {"status": "error", "error_type": "Stopped"}
            return {"status": "skipped"}

    agent.env = _Env()

    exp_ctx = {
        "global_padj": 0.01,
        "global_lfc": 1.5,
        "design": {
            "groups": {"A": ["r1", "r2", "r3"], "B": ["r4", "r5", "r6"]},
            "main_factor": "condition",
            "pseudobulk": {
                "from_obs": True,
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col": "cell_type",
                "comparisons": [["B", "A"]],
            },
        },
    }

    agent._run_pseudobulk("exp-p0-7", str(tmp_path / "x.h5ad"),
                          exp_ctx, {}, {})

    # The confirmed CP3 thresholds reached the pseudobulk script, not 0.05/0.5.
    assert captured["padj_max"] == 0.01
    assert captured["lfc_min"] == 1.5
