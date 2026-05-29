"""End-to-end coverage for the design checkpoint state machine.

This closes the F-ENG-E2E audit gap (2026-05-28): before this test the
orchestrator -> DesignAgent (CP2.1-2.6) -> CP2-plan control path was only
reachable through the interactive TUI and had no automated coverage.

It drives the DesignAgent state machine exactly as the headless runner /
interactive TUI would, using the deterministic answer policy, and asserts the
final confirmed design. Pure Python: no LLM, no subprocess.

Per ADR-011 the fixture uses neutral synthetic labels only (no real gene /
condition / dataset names).
"""

import pytest

from aria.headless import default_answer_policy


def _drive_design(agent, experiment_id, inferred_design):
    """Run start_design then resolve each published checkpoint via the policy.

    Returns (final_result, decisions) where decisions is [(cp, choice), ...].
    """
    from aria.bus.message_bus import bus

    exp_context = {
        "modalities": {"scRNA": ["/synthetic/path.h5ad"]},
        "inferred_design": inferred_design,
        "organism": inferred_design.get("organism"),
        "genome": inferred_design.get("genome"),
    }
    result = agent.start_design(
        experiment_id=experiment_id,
        exp_context=exp_context,
        biological_intent={"question": "compare A versus B in synthetic cells"},
    )

    decisions = []
    # The state machine advances one checkpoint per response. Resolve until
    # the agent reports it is done (or cancelled).
    for _ in range(20):
        pending = [
            m for m in bus.get_pending_checkpoints()
            if m.experiment_id == experiment_id and not m.payload.get("resolved")
        ]
        if not pending:
            break
        msg = pending[-1]  # newest unresolved design checkpoint
        cp_num = msg.payload.get("checkpoint")
        question = msg.payload.get("question", "")
        options = msg.payload.get("options", [])
        choice = default_answer_policy(cp_num, question, options)
        decisions.append((cp_num, choice))

        bus.resolve_checkpoint(msg.id, {"choice": choice})
        result = agent.handle_user_response(
            experiment_id=experiment_id,
            checkpoint_num=cp_num,
            choice=choice,
        )
        if result.get("status") in {"done", "cancelled"}:
            break
    return result, decisions


def test_design_state_machine_confirms_inferred_design(tmp_path):
    """Full CP2.1->2.6 walk on a high-confidence inferred design.

    Batch (2.4) and pseudorep (2.5) must auto-skip for clean donor-style
    replicates; the confirmed design must carry the inferred factor and
    formula, NOT options[0] defaults.
    """
    from aria.agents.design_agent import DesignAgent
    from aria.memory.memory import ARIAMemory

    inferred_design = {
        "source": "h5ad_obs",
        "confidence": "high",
        "organism": "Homo sapiens",
        "genome": "hg38",
        # donor-style ids that do NOT look like technical replicate suffixes
        # (r1/rep1/_1), so the pseudoreplication heuristic (CP2.5) auto-skips.
        "groups": {"COND_A": ["alpha", "beta", "gamma"],
                   "COND_B": ["alpha", "beta", "gamma"]},
        "main_factor": "grp",
        "condition_col": "grp",
        "replicate_col": "rep",
        "groupby_col": "ctype",
        "batch_covariate": None,
        "pseudobulk": {
            "from_obs": True,
            "condition_col": "grp",
            "replicate_col": "rep",
            "groupby_col": "ctype",
            "covariates": [],
            "comparisons": [["COND_A", "COND_B"]],
        },
    }

    agent = DesignAgent(memory=ARIAMemory(db_path=str(tmp_path / "mem.db")), llm=None)
    result, decisions = _drive_design(agent, "exp_e2e", inferred_design)

    assert result["status"] == "done", f"unexpected: {result}"
    design = result["design"]

    # Factor must be the inferred 'grp' (CP2.3 proposed 'grp'), not options[0]
    # which would be 'genotype'.
    assert design["main_factor"] == "grp"
    assert design["design_formula"] == "~ grp", design["design_formula"]
    assert design["organism"] == "Homo sapiens"
    assert design["genome"] == "hg38"
    assert set(design["groups"].keys()) == {"COND_A", "COND_B"}
    assert design["batch_covariate"] is None
    assert design["n_total_samples"] == 6
    assert design["pseudobulk"]["comparisons"] == [["COND_A", "COND_B"]]

    cps = [cp for cp, _ in decisions]
    # 2.1 groups, 2.2 organism, 2.3 factor, 2.6 confirm must all appear; the
    # batch/pseudorep checkpoints auto-skip for clean replicates.
    for required in (2.1, 2.2, 2.3, 2.6):
        assert required in cps, f"checkpoint {required} was not exercised: {cps}"
    assert 2.4 not in cps and 2.5 not in cps, f"batch/pseudorep should auto-skip: {cps}"


def test_design_state_machine_cancels_on_reject(tmp_path):
    """Rejecting the final design must cancel, not silently proceed."""
    from aria.agents.design_agent import DesignAgent
    from aria.memory.memory import ARIAMemory

    inferred_design = {
        "source": "h5ad_obs",
        "confidence": "high",
        "organism": "Mus musculus",
        "genome": "mm39",
        "groups": {"X": ["a", "b", "c"], "Y": ["a", "b", "c"]},
        "main_factor": "grp",
        "condition_col": "grp",
        "replicate_col": "rep",
        "groupby_col": "ctype",
        "batch_covariate": None,
        "pseudobulk": {"from_obs": True, "condition_col": "grp",
                       "replicate_col": "rep", "groupby_col": "ctype",
                       "covariates": [], "comparisons": [["X", "Y"]]},
    }

    def reject_at_confirm(cp_num, question, options):
        if cp_num == 2.6:
            return next((o for o in options if "cancel" in o.lower()), options[-1])
        return default_answer_policy(cp_num, question, options)

    from aria.bus.message_bus import bus

    agent = DesignAgent(memory=ARIAMemory(db_path=str(tmp_path / "mem.db")), llm=None)
    agent.start_design(
        experiment_id="exp_cancel",
        exp_context={"modalities": {"scRNA": ["/p.h5ad"]},
                     "inferred_design": inferred_design,
                     "organism": "Mus musculus", "genome": "mm39"},
        biological_intent={"question": "x vs y"},
    )
    result = {"status": "awaiting_user"}
    for _ in range(20):
        pending = [m for m in bus.get_pending_checkpoints()
                   if m.experiment_id == "exp_cancel" and not m.payload.get("resolved")]
        if not pending:
            break
        msg = pending[-1]
        cp_num = msg.payload.get("checkpoint")
        choice = reject_at_confirm(cp_num, msg.payload.get("question", ""),
                                   msg.payload.get("options", []))
        bus.resolve_checkpoint(msg.id, {"choice": choice})
        result = agent.handle_user_response(
            experiment_id="exp_cancel", checkpoint_num=cp_num, choice=choice)
        if result.get("status") in {"done", "cancelled"}:
            break

    assert result["status"] == "cancelled", result
