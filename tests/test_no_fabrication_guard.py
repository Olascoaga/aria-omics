"""P1-13 / P-FAKE-GUARD: mechanize ADR-002 (no silent fake science).

An AST scanner over aria/scripts and aria/agents that fails the build when it
finds the fabrication patterns prior audits had to remove by hand:

  A. returning an empty `AnnData()` as a placeholder data matrix (the P0-8
     integration_wnn pattern);
  B. a `_mock_*` / `mock_*` helper that returns a `status: "success"` result in a
     module that does NOT gate it through `mocks_allowed` (ungated fake success);
  C. returning a `hash(...)`-derived value as a metric (the B9 chromatin_qc TSS
     pattern).

The scanner is itself tested against synthetic good/bad modules so the guard is
proven to catch violations (failing-first), then applied to the real tree.
Documented, deliberate exceptions go in `_ALLOWLIST`.
"""

import ast
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCAN_DIRS = ["aria/scripts", "aria/agents"]

# (relpath, reason-substring) pairs that are deliberate, reviewed exceptions.
_ALLOWLIST: set[tuple[str, str]] = set()

_MOCK_NAME = re.compile(r"^_?mock_")


def _module_gates_mocks(tree: ast.AST) -> bool:
    for n in ast.walk(tree):
        if isinstance(n, ast.Name) and n.id == "mocks_allowed":
            return True
        if isinstance(n, ast.Attribute) and n.attr == "mocks_allowed":
            return True
    return False


def _call_name(call: ast.Call) -> str:
    f = call.func
    return f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")


def _returns_status_success(fn: ast.FunctionDef) -> bool:
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Dict):
            for k, v in zip(node.value.keys, node.value.values):
                if (isinstance(k, ast.Constant) and k.value == "status"
                        and isinstance(v, ast.Constant) and v.value == "success"):
                    return True
    return False


def scan_source(relpath: str, src: str) -> list[tuple[str, int, str]]:
    """Return a list of (relpath, lineno, reason) fabrication violations."""
    violations: list[tuple[str, int, str]] = []
    tree = ast.parse(src)
    gated = _module_gates_mocks(tree)

    for node in ast.walk(tree):
        # A. return empty AnnData()
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Call):
            if (_call_name(node.value) == "AnnData"
                    and not node.value.args and not node.value.keywords):
                violations.append((relpath, node.lineno,
                                   "returns empty AnnData() placeholder matrix"))
        # C. return a hash()-derived value (fabricated metric)
        if isinstance(node, ast.Return) and node.value is not None:
            for sub in ast.walk(node.value):
                if isinstance(sub, ast.Call) and _call_name(sub) == "hash":
                    violations.append((relpath, node.lineno,
                                       "returns a hash()-derived value (fabricated metric)"))
                    break
        # B. ungated _mock_* returning status:success
        if isinstance(node, ast.FunctionDef) and _MOCK_NAME.match(node.name):
            if _returns_status_success(node) and not gated:
                violations.append((relpath, node.lineno,
                                   f"{node.name}() returns status:success but the "
                                   f"module does not gate it via mocks_allowed"))

    return [v for v in violations
            if not any(v[0] == a and a_reason in v[2]
                       for a, a_reason in _ALLOWLIST)]


# ── The scanner catches the patterns (failing-first proof) ───────────────────

def test_scanner_flags_empty_anndata_return():
    bad = "import anndata as ad\ndef f():\n    return ad.AnnData()\n"
    v = scan_source("x.py", bad)
    assert any("AnnData" in r for _, _, r in v)


def test_scanner_flags_hash_derived_metric():
    bad = "def tss(f):\n    return {'tss': 1.0 + hash(f) % 30 / 10}\n"
    v = scan_source("x.py", bad)
    assert any("hash()" in r for _, _, r in v)


def test_scanner_flags_ungated_mock_success():
    bad = "def _mock_wnn(n):\n    return {'status': 'success', 'n': n}\n"
    v = scan_source("x.py", bad)
    assert any("mocks_allowed" in r for _, _, r in v)


def test_scanner_allows_gated_mock_success():
    ok = ("from aria.scripts._base import mocks_allowed\n"
          "def run(p):\n"
          "    if mocks_allowed(p):\n"
          "        return _mock(1)\n"
          "def _mock(n):\n"
          "    return {'status': 'success', 'n': n}\n")
    assert scan_source("x.py", ok) == []


# ── The real tree is clean ───────────────────────────────────────────────────

def test_repo_has_no_ungated_fabrication():
    violations: list[tuple[str, int, str]] = []
    for d in SCAN_DIRS:
        for path in (REPO / d).rglob("*.py"):
            rel = str(path.relative_to(REPO))
            violations.extend(scan_source(rel, path.read_text(encoding="utf-8")))
    assert not violations, "fabrication patterns found:\n" + "\n".join(
        f"  {r}:{ln} — {why}" for r, ln, why in violations)
