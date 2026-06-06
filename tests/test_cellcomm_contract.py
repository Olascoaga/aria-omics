"""P0-1 regression: LIANA must run through the agent path.

The scRNA agent dispatched `rna_cellcomm.py` with a `cell_type_col` param, but the
IPC contract for that script requires `groupby` (and the script itself read
`cell_type_col`). EnvironmentManager validates params against the contract BEFORE
the subprocess, so every agent-driven LIANA call failed at the contract gate with
"Missing required field 'groupby'" — cell-cell communication never ran from the
agent. These tests reproduce that failure at the agent and contract layers and
lock the canonical-`groupby` + `cell_type_col`-alias contract.
"""

from aria.agents.scrna_agent import scRNAAgent
from aria.utils.script_contracts import contract_for_script

_CELLCOMM = "aria/scripts/rna_cellcomm.py"


def _agent():
    agent = scRNAAgent.__new__(scRNAAgent)
    agent._log_decision = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    return agent


def test_agent_cellcomm_params_satisfy_the_contract(tmp_path):
    """Agent-level (not just script): the params the agent dispatches for LIANA
    must pass the rna_cellcomm IPC contract. Before P0-1 the agent sent
    `cell_type_col`, so the contract reported a missing `groupby` and LIANA was
    never executed."""
    agent = _agent()
    agent._workspace = lambda *a, **k: tmp_path

    clustered = tmp_path / "clustered.h5ad"
    clustered.write_bytes(b"")  # contract checks path existence, not contents

    captured = {}

    class _FakeEnv:
        def run_in_stack(self, *, stack, script_path, params):
            captured["script_path"] = script_path
            captured["params"] = params
            return {"status": "skipped", "reason": "fake env, no subprocess"}

    agent.env = _FakeEnv()

    agent._run_cell_communication(
        experiment_id="exp-p0-1",
        clustered_h5ad=str(clustered),
        exp_ctx={"organism": "Homo sapiens"},
        annotation={"label_col": "cell_type"},
    )

    params = captured["params"]
    # Canonical key is `groupby`, and it carries the chosen grouping column.
    assert "groupby" in params
    assert params["groupby"]

    # The real bug surface: the dispatched params must satisfy the contract.
    contract = contract_for_script(_CELLCOMM)
    issues = contract.validate_params(params)
    assert issues == [], f"agent params rejected by contract: {issues}"


def test_cellcomm_contract_accepts_cell_type_col_alias():
    """Backward compatibility: a legacy caller sending only `cell_type_col`
    (the pre-P0-1 alias) still satisfies the contract via the alias."""
    contract = contract_for_script(_CELLCOMM)
    issues = contract.validate_params({
        "data_path": __file__,           # an existing path
        "cell_type_col": "leiden",
    })
    assert issues == [], f"alias not accepted by contract: {issues}"


def test_cellcomm_contract_still_rejects_when_neither_key_present():
    """The grouping column is still mandatory; dropping both keys must fail."""
    contract = contract_for_script(_CELLCOMM)
    issues = contract.validate_params({"data_path": __file__})
    assert any(i.field == "groupby" for i in issues)


def test_script_resolves_groupby_then_cell_type_col_alias():
    """The script reads the canonical `groupby` first and falls back to the
    `cell_type_col` alias, then to the default."""
    from aria.scripts.rna_cellcomm import _resolve_groupby

    assert _resolve_groupby({"groupby": "cluster"}) == "cluster"
    assert _resolve_groupby({"cell_type_col": "leiden"}) == "leiden"
    assert _resolve_groupby({"groupby": "cluster",
                             "cell_type_col": "leiden"}) == "cluster"
    assert _resolve_groupby({}) == "cell_type"


def test_cellcomm_has_no_embedded_ligand_receptor_fallback():
    """When LIANA is absent, ARIA should skip honestly instead of emitting
    interactions from a baked-in ligand-receptor list."""
    import inspect
    import aria.scripts.rna_cellcomm as mod

    source = inspect.getsource(mod)
    assert "_LR_PAIRS" not in source
    assert "mean_expression_fallback" not in source
    assert "liana_not_installed" in source
