"""D10 debt ratchet: structural debt in ``aria/`` may only go down.

Limits live in ``tests/debt_ratchet.json``. Each audit slice that closes lowers them.
Raising a limit requires an ADR reference in the JSON (``"adr"`` field per metric).

Metrics are computed with the standard-library ``ast`` module only, so the guard runs in
every CI job without extra dependencies.

- ``unannotated_broad_except``: ``except Exception`` / ``except BaseException`` / bare
  ``except`` handlers whose ``except`` line lacks ``# aria: allow-broad-except <reason>``.
- ``silent_broad_except``: the subset whose body only does
  ``pass``/``continue``/``return``/``return None``/``return <literal>``/``...``.
- ``zip_without_strict``: ``zip(...)`` calls without a ``strict=`` keyword in
  ``aria/scripts/`` and ``aria/utils/`` (scientific kernels).
- ``functions_over_300_lines``: functions/methods longer than 300 lines.
"""
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "aria"
LIMITS_PATH = Path(__file__).with_name("debt_ratchet.json")
ANNOTATION_RE = re.compile(r"#\s*aria:\s*allow-broad-except\s+\S")
KERNEL_DIRS = (PKG / "scripts", PKG / "utils")
_BROAD = {"Exception", "BaseException"}


def _is_broad(handler: ast.ExceptHandler) -> bool:
    t = handler.type
    if t is None:
        return True
    names = t.elts if isinstance(t, ast.Tuple) else [t]
    return any(isinstance(n, ast.Name) and n.id in _BROAD for n in names)


def _is_silent(handler: ast.ExceptHandler) -> bool:
    for stmt in handler.body:
        if isinstance(stmt, (ast.Pass, ast.Continue)):
            continue
        if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Constant):
            continue  # docstring / Ellipsis
        if isinstance(stmt, ast.Return) and (
            stmt.value is None or isinstance(stmt.value, ast.Constant)
            or (isinstance(stmt.value, (ast.List, ast.Dict, ast.Tuple, ast.Set))
                and not getattr(stmt.value, "elts", getattr(stmt.value, "keys", [1])))
        ):
            continue
        return False
    return True


def _iter_py(base: Path):
    yield from sorted(p for p in base.rglob("*.py") if "__pycache__" not in p.parts)


def collect_metrics() -> dict[str, dict]:
    """Return ``{metric: {"count": int, "sites": [...]}}`` for the package."""
    out = {k: {"count": 0, "sites": []} for k in (
        "unannotated_broad_except", "silent_broad_except",
        "zip_without_strict", "functions_over_300_lines")}
    for path in _iter_py(PKG):
        src = path.read_text(encoding="utf-8")
        lines = src.splitlines()
        tree = ast.parse(src, filename=str(path))
        rel = path.relative_to(ROOT).as_posix()
        in_kernel = any(path.is_relative_to(d) for d in KERNEL_DIRS)
        for node in ast.walk(tree):
            if isinstance(node, ast.ExceptHandler) and _is_broad(node):
                line = lines[node.lineno - 1] if node.lineno - 1 < len(lines) else ""
                if not ANNOTATION_RE.search(line):
                    out["unannotated_broad_except"]["sites"].append(f"{rel}:{node.lineno}")
                    if _is_silent(node):
                        out["silent_broad_except"]["sites"].append(f"{rel}:{node.lineno}")
            elif (in_kernel and isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Name) and node.func.id == "zip"
                  and not any(k.arg == "strict" for k in node.keywords)):
                out["zip_without_strict"]["sites"].append(f"{rel}:{node.lineno}")
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                n = (node.end_lineno or node.lineno) - node.lineno + 1
                if n > 300:
                    out["functions_over_300_lines"]["sites"].append(f"{rel}:{node.lineno}:{node.name}:{n}")
    for v in out.values():
        v["count"] = len(v["sites"])
    return out


def _limits() -> dict:
    return json.loads(LIMITS_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def metrics():
    return collect_metrics()


@pytest.mark.parametrize("metric", [
    "unannotated_broad_except", "silent_broad_except",
    "zip_without_strict", "functions_over_300_lines",
])
def test_debt_metric_does_not_exceed_ratchet(metrics, metric):
    limit = _limits()["limits"][metric]["max"]
    got = metrics[metric]["count"]
    assert got <= limit, (
        f"{metric} rose to {got} (ratchet max {limit}). New sites must be fixed, or "
        f"annotated with '# aria: allow-broad-except <reason>' where charter §1 applies. "
        f"First sites: {metrics[metric]['sites'][:10]}"
    )


def test_ratchet_limits_are_tight(metrics):
    """A limit far above reality lets debt creep back in silently; lower it when you fix debt."""
    limits = _limits()["limits"]
    slack = {m: limits[m]["max"] - metrics[m]["count"] for m in limits}
    loose = {m: s for m, s in slack.items() if s > limits[m].get("slack", 0)}
    assert not loose, f"Lower these ratchet limits in tests/debt_ratchet.json: {loose}"


def test_annotations_carry_a_reason():
    bad = []
    for path in _iter_py(PKG):
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if "aria: allow-broad-except" in line and not ANNOTATION_RE.search(line):
                bad.append(f"{path.relative_to(ROOT).as_posix()}:{i}")
    assert not bad, f"'# aria: allow-broad-except' needs a reason: {bad}"


if __name__ == "__main__":  # print current metrics to (re)seed the limits
    print(json.dumps({k: v["count"] for k, v in collect_metrics().items()}, indent=2))
