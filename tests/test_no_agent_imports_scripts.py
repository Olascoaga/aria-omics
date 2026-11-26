"""D1-A1 layering guard: modules under ``aria/agents/`` must not import ``aria.scripts``.

Compute scripts under ``aria/scripts/`` run in isolated Conda stacks and may
depend on heavy scientific packages; agents and narrators under
``aria/agents/`` must share pure contracts via ``aria/utils/`` instead.

Parses the tree with ``ast`` only — runs in any environment.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_AGENTS_DIR = _REPO_ROOT / "aria" / "agents"
_BULK_QC_CONTRACT = _REPO_ROOT / "aria" / "utils" / "bulk_qc_contract.py"

# Pre-existing function-local lazy script imports present on main (8a0cb68).
# Frozen as a strict ratchet: no top-level imports are allowed anywhere, and
# no new lazy imports from aria.scripts may be added under aria/agents/.
_LEGACY_LAZY_SCRIPT_IMPORTS: frozenset[tuple[str, str]] = frozenset({
    ("aria/agents/_narrative_chromatin.py", "aria.scripts.chromatin_figure_bulk"),
    ("aria/agents/chromatin_agent.py", "aria.scripts.chromatin_peak_counts"),
    ("aria/agents/chromatin_agent.py", "aria.scripts.chromatin_qc"),
    ("aria/agents/narrative/scrna/figures.py", "aria.scripts.rna_pathway_viz"),
})


def _resolve_import_targets(node: ast.AST, file_path: Path) -> list[str]:
    """Return fully-qualified imported module/symbol targets for an import node."""
    targets: list[str] = []
    if isinstance(node, ast.Import):
        for alias in node.names:
            targets.append(alias.name)
    elif isinstance(node, ast.ImportFrom):
        if node.level == 0:
            base = node.module or ""
        else:
            rel_parts = list(file_path.relative_to(_REPO_ROOT).with_suffix("").parts)
            if rel_parts and rel_parts[-1] == "__init__":
                pkg_parts = rel_parts[:-1]
            else:
                pkg_parts = rel_parts[:-1]
            up = node.level - 1
            prefix = pkg_parts[: len(pkg_parts) - up] if up < len(pkg_parts) else []
            base = ".".join(prefix + ([node.module] if node.module else []))
        if base == "aria.scripts":
            for alias in node.names:
                targets.append(f"aria.scripts.{alias.name}")
        elif base:
            targets.append(base)
    return targets


def _is_scripts_target(target: str) -> bool:
    return target == "aria.scripts" or target.startswith("aria.scripts.")


def _iter_module_level_nodes(stmts: list[ast.stmt]):
    """Yield AST statements at module import time (excluding def/class bodies)."""
    for stmt in stmts:
        yield stmt
        if isinstance(stmt, (ast.If, ast.For, ast.AsyncFor, ast.While, ast.With, ast.AsyncWith)):
            yield from _iter_module_level_nodes(stmt.body)
            yield from _iter_module_level_nodes(stmt.orelse)
        elif isinstance(stmt, ast.Try):
            yield from _iter_module_level_nodes(stmt.body)
            yield from _iter_module_level_nodes(stmt.orelse)
            yield from _iter_module_level_nodes(stmt.finalbody)
            for handler in stmt.handlers:
                yield from _iter_module_level_nodes(handler.body)


def _agent_py_files() -> list[Path]:
    return sorted(
        p for p in _AGENTS_DIR.rglob("*.py")
        if "__pycache__" not in p.parts
    )


def test_no_top_level_agent_imports_from_aria_scripts():
    violations: list[str] = []
    for path in _agent_py_files():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        rel = path.relative_to(_REPO_ROOT).as_posix()
        for node in _iter_module_level_nodes(tree.body):
            if not isinstance(node, (ast.Import, ast.ImportFrom)):
                continue
            for target in _resolve_import_targets(node, path):
                if _is_scripts_target(target):
                    violations.append(f"{rel}:{node.lineno} -> {target}")
    assert not violations, (
        "Module-level imports from aria.scripts are forbidden under aria/agents/ "
        "(move shared helpers to aria/utils/):\n  " + "\n  ".join(violations)
    )


def test_no_new_lazy_agent_imports_from_aria_scripts():
    violations: list[str] = []
    for path in _agent_py_files():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        rel = path.relative_to(_REPO_ROOT).as_posix()
        for node in ast.walk(tree):
            if not isinstance(node, (ast.Import, ast.ImportFrom)):
                continue
            for target in _resolve_import_targets(node, path):
                if _is_scripts_target(target) and (rel, target) not in _LEGACY_LAZY_SCRIPT_IMPORTS:
                    violations.append(f"{rel}:{node.lineno} -> {target}")
    assert not violations, (
        "Unexpected aria.scripts import(s) under aria/agents/:\n  "
        + "\n  ".join(violations)
    )


def test_bulk_qc_contract_is_stdlib_only():
    assert _BULK_QC_CONTRACT.is_file(), f"Missing {_BULK_QC_CONTRACT}"
    tree = ast.parse(_BULK_QC_CONTRACT.read_text(encoding="utf-8"), filename=str(_BULK_QC_CONTRACT))
    allowed = set(sys.stdlib_module_names) | {"__future__"}
    bad: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".", 1)[0]
                if root not in allowed:
                    bad.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level > 0:
                bad.append(f"relative:{node.module or ''}")
            else:
                root = (node.module or "").split(".", 1)[0]
                if root not in allowed:
                    bad.append(node.module or "")
    assert not bad, f"aria/utils/bulk_qc_contract.py must be stdlib-only; found: {bad}"
