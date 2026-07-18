"""P2-1: modern packaging + version-drift guard.

These checks are pure (stdlib `tomllib`, 3.11+) and import-light (`aria.version`
is stdlib-only), so they run in every CI lane. They pin the PEP 621 contract and
catch the README/env/setup `3.10` vs `3.11` drift from regressing.
"""

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"


def _proj():
    with open(PYPROJECT, "rb") as f:
        return tomllib.load(f)["project"]


def test_pyproject_exists_and_parses():
    assert PYPROJECT.exists(), "pyproject.toml must exist (PEP 621 packaging)"
    data = tomllib.loads(PYPROJECT.read_text())
    assert data["project"]["name"] == "aria-omics"


def test_requires_python_is_3_11():
    # The whole project (envs, install.sh, CI, README) standardized on 3.11;
    # packaging metadata must agree (was the lone ">=3.10").
    assert _proj()["requires-python"] == ">=3.11"


def test_license_uses_pep_639_spdx_expression():
    project = _proj()
    assert project["license"] == "MIT"
    assert not any(
        classifier.startswith("License ::")
        for classifier in project.get("classifiers", [])
    )


def test_no_python_310_classifier():
    classifiers = _proj().get("classifiers", [])
    assert not any("3.10" in c for c in classifiers)
    assert any("3.11" in c for c in classifiers)


def test_version_is_dynamic_from_single_source():
    data = tomllib.loads(PYPROJECT.read_text())
    assert "version" in data["project"].get("dynamic", [])
    attr = data["tool"]["setuptools"]["dynamic"]["version"]["attr"]
    assert attr == "aria.version.__version__"
    # And it actually resolves to the single source of truth.
    sys.path.insert(0, str(ROOT))
    from aria.version import __version__
    assert re.match(r"^\d+\.\d+", __version__)


def test_core_dependencies_have_version_ceilings():
    # P2-1: pin version ceilings so an unvetted major can't silently break a run.
    deps = _proj()["dependencies"]
    assert deps, "core dependencies must be declared"
    missing = [d for d in deps if "<" not in d]
    assert not missing, f"core deps without an upper bound: {missing}"


def test_console_scripts_present():
    scripts = _proj().get("scripts", {})
    assert scripts.get("aria") == "aria.tui:main"
    assert scripts.get("aria-doctor") == "aria.doctor:console_main"
