"""P1-12: every dispatchable script must declare an IPC contract.

EnvironmentManager.run_in_stack validates a script's params against its
ScriptContract before the subprocess. A dispatchable script with no contract
silently skips that gate (the FASTQ flow's rna_align / rna_quantify were the
original gap). This test AST-discovers every script dispatched via run_in_stack
across aria/agents and fails if one that exists on disk lacks a contract.
"""

import ast
import pathlib

from aria.utils.script_contracts import SCRIPT_CONTRACTS, contract_for_script

REPO = pathlib.Path(__file__).resolve().parents[1]
AGENTS = REPO / "aria" / "agents"


def _script_path_of(call: ast.Call) -> str | None:
    """Extract the script path from a run_in_stack(...) call.

    Supports both the keyword form (script_path="...") and the positional form
    run_in_stack(stack, script_path, params).
    """
    for kw in call.keywords:
        if kw.arg == "script_path" and isinstance(kw.value, ast.Constant):
            return kw.value.value
    if len(call.args) >= 2 and isinstance(call.args[1], ast.Constant):
        value = call.args[1].value
        if isinstance(value, str) and value.startswith("aria/scripts/"):
            return value
    return None


def _discover_dispatched_scripts() -> set[str]:
    found: set[str] = set()
    for path in AGENTS.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "run_in_stack"):
                sp = _script_path_of(node)
                if sp:
                    found.add(sp)
    return found


def test_every_dispatched_script_on_disk_has_a_contract():
    dispatched = _discover_dispatched_scripts()
    # Only scripts that actually exist must have a contract; planned-but-absent
    # scripts (e.g. v4.6 chromatin_differential/motifs) return
    # script_not_implemented and are intentionally not yet contracted.
    on_disk = {s for s in dispatched if (REPO / s).exists()}
    missing = sorted(s for s in on_disk if s not in SCRIPT_CONTRACTS)
    assert not missing, (
        "dispatchable scripts missing an IPC contract:\n  "
        + "\n  ".join(missing)
    )
    # Sanity: discovery is non-vacuous and the FASTQ-flow gap is covered.
    assert "aria/scripts/rna_align.py" in on_disk
    assert "aria/scripts/rna_quantify.py" in on_disk


def test_contract_lookup_resolves_for_each_contract_key():
    for key in SCRIPT_CONTRACTS:
        assert contract_for_script(key) is not None
        assert (REPO / key).exists(), f"contract for non-existent script: {key}"
