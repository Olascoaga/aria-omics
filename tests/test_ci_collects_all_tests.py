"""CI test-collection guard (D3-001, D3-015).

Parses `.github/workflows/ci.yml` and verifies that EVERY `tests/**/test_*.py`
file on disk is either:
  1. Collected by at least one CI job (via dynamic `pytest tests/` or an explicit
     target), OR
  2. Listed in `EXPLICIT_CI_EXCLUSIONS` with a concrete, non-empty reason.

Also enforces structural anti-drift invariants so CI can never silently regress
back to hand-maintained allowlists or 0-test ghost files:
  - Core jobs (`pr-checks`, `rna-env-tests`, `hic-env-tests`,
    `integration-env-tests`) collect `tests/` dynamically rather than hand-listing
    test filenames.
  - Every `--ignore=` target in `ci.yml` exists on disk.
  - No CI job references a non-existent test file, a `conftest.py:collect_ignore`
    entry, or a 0-test ghost file.
  - All statistical oracle suites (`tests/oracle/test_oracle_*.py`) are collected
    by both `pr-checks` and `rna-env-tests`, and modality-specific suites are
    matched by `hic-env-tests` and `integration-env-tests`.
"""
from __future__ import annotations

import ast
import fnmatch
import re
import shlex
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CI_YML_PATH = REPO_ROOT / ".github" / "workflows" / "ci.yml"
CONFTEST_PATH = REPO_ROOT / "tests" / "conftest.py"

# Explicit, reasoned exclusions for legacy top-level executable scripts that are
# kept in `tests/conftest.py:collect_ignore` because importing them executes
# top-level code and `sys.exit()`. Any new native pytest file MUST NOT be added
# here — it must be collected by CI automatically.
EXPLICIT_CI_EXCLUSIONS: dict[str, str] = {
    "tests/test_chromatin_agent.py": (
        "Legacy standalone executable script (top-level sys.exit, 0 pytest functions; "
        "excluded by tests/conftest.py:collect_ignore)."
    ),
    "tests/test_debate_council.py": (
        "Legacy standalone executable script (top-level sys.exit, 0 pytest functions; "
        "excluded by tests/conftest.py:collect_ignore)."
    ),
    "tests/test_genome_arch_agent.py": (
        "Legacy standalone executable script (top-level sys.exit, 0 pytest functions; "
        "excluded by tests/conftest.py:collect_ignore)."
    ),
    "tests/test_integration.py": (
        "Legacy standalone executable script (top-level sys.exit, 0 pytest functions; "
        "excluded by tests/conftest.py:collect_ignore)."
    ),
    "tests/test_narrative_agent.py": (
        "Legacy standalone executable script (top-level sys.exit, 0 pytest functions; "
        "excluded by tests/conftest.py:collect_ignore)."
    ),
    "tests/test_scrna_e2e.py": (
        "Legacy standalone executable script (top-level sys.exit, 0 pytest functions; "
        "excluded by tests/conftest.py:collect_ignore)."
    ),
}


def _load_conftest_collect_ignore() -> set[str]:
    """Parse `collect_ignore = [...]` from tests/conftest.py via AST."""
    tree = ast.parse(CONFTEST_PATH.read_text(encoding="utf-8"))
    ignored: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == "collect_ignore":
                    val = ast.literal_eval(node.value)
                    for item in val:
                        ignored.add(f"tests/{item.lstrip('./')}")
    return ignored


def _discover_test_files() -> set[str]:
    """Return all repo-relative `tests/**/test_*.py` paths."""
    return {
        p.relative_to(REPO_ROOT).as_posix()
        for p in (REPO_ROOT / "tests").rglob("test_*.py")
        if "__pycache__" not in p.parts
    }


def _count_ast_test_functions(rel_path: str) -> int:
    """Count top-level and Test* class `def test_*` functions in a test file."""
    tree = ast.parse((REPO_ROOT / rel_path).read_text(encoding="utf-8"))
    count = 0
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
            count += 1
        elif isinstance(node, ast.ClassDef) and node.name.startswith("Test"):
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name.startswith("test_"):
                    count += 1
    return count


def _file_searchable_tokens(rel_path: str) -> str:
    """Return lowercased file path + class/test function names for `-k` matching."""
    tree = ast.parse((REPO_ROOT / rel_path).read_text(encoding="utf-8"))
    tokens = [rel_path]
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
            tokens.append(node.name)
        elif isinstance(node, ast.ClassDef) and node.name.startswith("Test"):
            tokens.append(node.name)
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name.startswith("test_"):
                    tokens.append(item.name)
    return " ".join(tokens).lower()


def _matches_k_expr(k_expr: str | None, searchable_text: str) -> bool:
    """Evaluate simple pytest `-k` OR-of-terms expressions used in ci.yml."""
    if not k_expr:
        return True
    # Split on ` or ` (case-insensitive); all our CI -k expressions are unions of keywords.
    terms = [t.strip().lower() for t in re.split(r"\s+or\s+", k_expr, flags=re.IGNORECASE) if t.strip()]
    return any(term in searchable_text for term in terms)


def _extract_pytest_invocations(ci_data: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Parse every `pytest` command in every job of `.github/workflows/ci.yml`."""
    jobs_invocations: dict[str, list[dict[str, Any]]] = {}
    for job_name, job_cfg in (ci_data.get("jobs") or {}).items():
        invocations: list[dict[str, Any]] = []
        for step in job_cfg.get("steps") or []:
            run_block = step.get("run")
            if not run_block or "pytest" not in run_block:
                continue
            # Join shell line continuations (`\` + newline)
            joined = re.sub(r"\\\s*\n", " ", str(run_block))
            # Process line by line, unwrapping `bash -lc "..."` BEFORE splitting on && / ;
            for line in joined.splitlines():
                line = line.strip()
                if "pytest" not in line:
                    continue
                m_bash = re.search(r'bash\s+-lc\s+"([^"]+)"', line)
                body = m_bash.group(1) if m_bash else line
                for sub in re.split(r"(?:&&|;)", body):
                    if "pytest" not in sub:
                        continue
                    try:
                        tokens = shlex.split(sub)
                    except ValueError:
                        tokens = sub.split()
                    # Locate actual pytest invocation (either `pytest ...` or `python -m pytest ...`,
                    # NOT `pip install pytest ...`)
                    idx = None
                    for i, tok in enumerate(tokens):
                        if (tok == "pytest" or tok.endswith("/pytest")) and (
                            i == 0 or tokens[i - 1] == "-m"
                        ):
                            idx = i
                            break
                    if idx is None:
                        continue
                    args = tokens[idx + 1 :]
                    targets: list[str] = []
                    ignores: list[str] = []
                    ignore_globs: list[str] = []
                    k_expr: str | None = None
                    i = 0
                    while i < len(args):
                        a = args[i]
                        if a.startswith("--ignore="):
                            ignores.append(a.split("=", 1)[1].removeprefix("/opt/aria/").rstrip("/"))
                        elif a == "--ignore" and i + 1 < len(args):
                            ignores.append(args[i + 1].removeprefix("/opt/aria/").rstrip("/"))
                            i += 1
                        elif a.startswith("--ignore-glob="):
                            ignore_globs.append(a.split("=", 1)[1])
                        elif a == "--ignore-glob" and i + 1 < len(args):
                            ignore_globs.append(args[i + 1])
                            i += 1
                        elif a == "-k" and i + 1 < len(args):
                            k_expr = args[i + 1]
                            i += 1
                        elif a.startswith("-k"):
                            k_expr = a[2:]
                        elif a in ("--deselect", "-m", "-o", "-p", "-W", "-c", "--rootdir") and i + 1 < len(args):
                            i += 1
                        elif a.startswith("-"):
                            pass
                        else:
                            norm = a.removeprefix("/opt/aria/").split("::", 1)[0]
                            targets.append(norm)
                        i += 1
                    invocations.append({
                        "targets": targets,
                        "ignores": ignores,
                        "ignore_globs": ignore_globs,
                        "k_expr": k_expr,
                    })
        if invocations:
            jobs_invocations[job_name] = invocations
    return jobs_invocations


def _files_collected_by_invocation(
    inv: dict[str, Any],
    all_files: set[str],
    conftest_ignored: set[str],
    token_cache: dict[str, str],
) -> set[str]:
    """Compute the exact set of `tests/**/test_*.py` files collected by a pytest invocation."""
    selected: set[str] = set()
    ignores = set(inv["ignores"])
    ignore_globs = inv["ignore_globs"]
    k_expr = inv["k_expr"]

    for target in inv["targets"]:
        norm = target.rstrip("/")
        if norm in all_files:
            candidates = {norm}
        else:
            prefix = f"{norm}/"
            candidates = {f for f in all_files if f.startswith(prefix)}
        for f in candidates:
            if f in conftest_ignored:
                continue
            if f in ignores or any(f.startswith(f"{ig}/") for ig in ignores):
                continue
            if any(fnmatch.fnmatch(f, g) for g in ignore_globs):
                continue
            if k_expr and not _matches_k_expr(k_expr, token_cache[f]):
                continue
            selected.add(f)
    return selected


def _collect_ci_coverage(
    all_files: set[str],
    exclusions: dict[str, str],
) -> tuple[dict[str, set[str]], set[str]]:
    """Return (per_job_collected_files, uncovered_files)."""
    yaml = pytest.importorskip("yaml")
    ci_data = yaml.safe_load(CI_YML_PATH.read_text(encoding="utf-8"))
    conftest_ignored = _load_conftest_collect_ignore()
    token_cache = {f: _file_searchable_tokens(f) for f in all_files}
    job_invs = _extract_pytest_invocations(ci_data)

    per_job: dict[str, set[str]] = {}
    union_collected: set[str] = set()
    for job_name, invs in job_invs.items():
        job_files: set[str] = set()
        for inv in invs:
            job_files |= _files_collected_by_invocation(inv, all_files, conftest_ignored, token_cache)
        per_job[job_name] = job_files
        union_collected |= job_files

    uncovered = all_files - union_collected - set(exclusions)
    return per_job, uncovered


def test_every_test_file_is_collected_by_ci_or_explicitly_excluded():
    """No `tests/**/test_*.py` file may silently escape CI collection."""
    all_files = _discover_test_files()
    per_job, uncovered = _collect_ci_coverage(all_files, EXPLICIT_CI_EXCLUSIONS)

    assert not uncovered, (
        "Test file(s) neither collected by any CI job in .github/workflows/ci.yml "
        "nor listed in EXPLICIT_CI_EXCLUSIONS:\n  " + "\n  ".join(sorted(uncovered))
    )

    # Every exclusion must exist, carry a concrete reason, and match conftest collect_ignore.
    conftest_ignored = _load_conftest_collect_ignore()
    for excl_path, reason in EXPLICIT_CI_EXCLUSIONS.items():
        assert excl_path in all_files, f"Stale EXPLICIT_CI_EXCLUSIONS entry: {excl_path}"
        assert isinstance(reason, str) and len(reason.strip()) >= 15, (
            f"Exclusion reason for {excl_path} must be a descriptive explanation"
        )
        assert excl_path in conftest_ignored, (
            f"{excl_path} is excluded from CI but NOT in tests/conftest.py:collect_ignore"
        )
        assert all(excl_path not in files for files in per_job.values()), (
            f"{excl_path} is listed in EXPLICIT_CI_EXCLUSIONS yet collected by CI"
        )


def test_ci_uses_dynamic_directory_collection_and_valid_paths():
    """Core CI jobs must collect `tests/` dynamically and never reference ghost/missing files."""
    yaml = pytest.importorskip("yaml")
    ci_data = yaml.safe_load(CI_YML_PATH.read_text(encoding="utf-8"))
    conftest_ignored = _load_conftest_collect_ignore()
    job_invs = _extract_pytest_invocations(ci_data)

    for required_job in ("pr-checks", "rna-env-tests", "hic-env-tests", "integration-env-tests"):
        assert required_job in job_invs, f"Missing expected CI job: {required_job}"
        all_targets = [t.rstrip("/") for inv in job_invs[required_job] for t in inv["targets"]]
        assert "tests" in all_targets, (
            f"CI job '{required_job}' must collect `tests/` dynamically instead of "
            f"hand-maintaining a file allowlist; found targets={all_targets}"
        )

    for job_name, invs in job_invs.items():
        for inv in invs:
            for ig in inv["ignores"]:
                assert (REPO_ROOT / ig).exists(), (
                    f"CI job '{job_name}' has stale --ignore target: {ig}"
                )
            for tgt in inv["targets"]:
                norm = tgt.rstrip("/")
                assert (REPO_ROOT / norm).exists(), (
                    f"CI job '{job_name}' references non-existent test target: {tgt}"
                )
                assert norm not in conftest_ignored, (
                    f"CI job '{job_name}' explicitly references conftest-ignored ghost file: {norm}"
                )
                if norm.endswith(".py"):
                    assert _count_ast_test_functions(norm) > 0, (
                        f"CI job '{job_name}' explicitly references 0-test ghost file: {norm}"
                    )


def test_oracle_and_modality_suites_reach_their_target_ci_lanes():
    """Verify every oracle suite and modality stack reaches its matching CI job."""
    all_files = _discover_test_files()
    per_job, _ = _collect_ci_coverage(all_files, EXPLICIT_CI_EXCLUSIONS)

    oracle_files = {f for f in all_files if f.startswith("tests/oracle/test_oracle_")}
    assert len(oracle_files) == 10, f"Expected 10 oracle files, found {sorted(oracle_files)}"
    assert oracle_files <= per_job["pr-checks"]
    assert oracle_files <= per_job["rna-env-tests"]

    # Hi-C lane reaches both the new differential-3D oracle and Hi-C E2E suite
    assert "tests/oracle/test_oracle_hic_differential_3d.py" in per_job["hic-env-tests"]
    assert "tests/test_hic_production_e2e.py" in per_job["hic-env-tests"]

    # Multiome/integration lane reaches peak2gene, Layer D, SEACR/IDR, WNN, and MOFA+ suites
    for expected in (
        "tests/oracle/test_oracle_peak2gene_core.py",
        "tests/oracle/test_oracle_cross_modal_layer_d.py",
        "tests/oracle/test_oracle_seacr_idr.py",
        "tests/test_v48_i3_wnn_guardrail.py",
        "tests/test_v48_i4_mofa.py",
    ):
        assert expected in per_job["integration-env-tests"], (
            f"{expected} not selected by integration-env-tests"
        )

    # D3-015 regression guard: test_integration_agent.py is now a native pytest suite
    assert "tests/test_integration_agent.py" in per_job["pr-checks"]
    assert _count_ast_test_functions("tests/test_integration_agent.py") >= 7


def test_guard_parser_catches_simulated_uncollected_test_file():
    """Meta-test proving the guard actually fails if a test file escapes CI collection."""
    fake_files = _discover_test_files() | {"tests/test_cockpit_render.py"}
    # If pr-checks were stripped out, test_cockpit_render.py (ignored in conda lanes)
    # MUST be flagged as uncovered!
    yaml = pytest.importorskip("yaml")
    ci_data = yaml.safe_load(CI_YML_PATH.read_text(encoding="utf-8"))
    ci_data["jobs"].pop("pr-checks", None)
    conftest_ignored = _load_conftest_collect_ignore()
    token_cache = {f: _file_searchable_tokens(f) for f in fake_files}
    job_invs = _extract_pytest_invocations(ci_data)
    union = set()
    for invs in job_invs.values():
        for inv in invs:
            union |= _files_collected_by_invocation(inv, fake_files, conftest_ignored, token_cache)
    uncovered_without_pr = fake_files - union - set(EXPLICIT_CI_EXCLUSIONS)
    assert "tests/test_cockpit_render.py" in uncovered_without_pr
