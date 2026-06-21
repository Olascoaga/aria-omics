"""T16 (tri-auditoría 2026-06-14): hermetic-reproducibility guard for the numpy
pin across the conda + pip lockfiles.

Gemini's "NumPy 2.0 ABI break / 91-of-117 fail" Blocker did NOT reproduce: aria-env
runs numpy 2.4.x green, and the three modality envs each import their full
scientific stack with no ABI error (verified live). This test is the durable fence
for what that verification established about the lockfiles themselves:

  * each env pins numpy in exactly one place (no conda-lock vs pip-lock conflict);
  * the numpy MAJOR per env matches intent — the RNA/report stack is numpy 2.x,
    while the chromatin (MACS2 needs numpy<2) and ingestion stacks are numpy 1.x.
  * the human-authored YAML constraints AGREE with the locked major (S1,
    pre-integration audit 2026-06-20). The RNA YAMLs used to cap numpy<2.0.0 while
    the validated lock + pydeseq2 0.5.4 require numpy 2.x — recreating the env from
    YAML produced a stack that broke at the pydeseq2 import. This fence now fails
    if a YAML constraint excludes the env's locked numpy major.

It parses files only (no scientific imports), so it runs in any env.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from packaging.specifiers import SpecifierSet
from packaging.version import Version

_ENVS_DIR = Path(__file__).resolve().parents[1] / "envs"

_CONDA_NUMPY_RE = re.compile(r"/numpy-(\d+)\.(\d+)\.(\d+)")
_PIP_NUMPY_RE = re.compile(r"(?mi)^numpy==(\d+)\.(\d+)\.(\d+)\b")
# A conda YAML numpy line, e.g. "  - numpy>=2.0,<3.0" (trailing comment ignored).
_YAML_NUMPY_RE = re.compile(r"(?m)^\s*-\s*numpy\s*([<>=!,\d.\s]*)")

# Expected numpy MAJOR per env (the T16 verified intent).
_EXPECTED_MAJOR = {"rna": 2, "chromatin": 1, "ingestion": 1}

# Every YAML that declares the env's numpy constraint (env + any CI/variant YAMLs
# that must stay in sync with it).
_ENV_YAMLS = {
    "rna": ["aria-rna-env.yml", "aria-rna-ci.yml"],
    "chromatin": ["aria-chromatin-env.yml"],
    "ingestion": ["aria-ingestion-env.yml"],
}


def _conda_numpy(env: str) -> tuple[int, int, int] | None:
    lock = _ENVS_DIR / f"aria-{env}-env.linux-64.lock"
    if not lock.exists():
        return None
    m = _CONDA_NUMPY_RE.search(lock.read_text())
    return tuple(int(x) for x in m.groups()) if m else None


def _pip_numpy(env: str) -> tuple[int, int, int] | None:
    lock = _ENVS_DIR / f"aria-{env}-env.pip.lock"
    if not lock.exists():
        return None
    m = _PIP_NUMPY_RE.search(lock.read_text())
    return tuple(int(x) for x in m.groups()) if m else None


@pytest.mark.parametrize("env", sorted(_EXPECTED_MAJOR))
def test_env_pins_numpy_exactly_once_and_consistently(env):
    conda_v = _conda_numpy(env)
    pip_v = _pip_numpy(env)

    # numpy must be pinned somewhere for a hermetic install.
    assert conda_v or pip_v, f"aria-{env}-env pins numpy in neither lockfile"

    # If both lockfiles pin numpy, they must agree (no conda/pip ABI split).
    if conda_v and pip_v:
        assert conda_v == pip_v, (
            f"aria-{env}-env numpy conflict: conda-lock {conda_v} vs pip-lock {pip_v}"
        )

    major = (conda_v or pip_v)[0]
    assert major == _EXPECTED_MAJOR[env], (
        f"aria-{env}-env numpy major {major} != expected {_EXPECTED_MAJOR[env]}"
    )


def test_rna_stack_is_numpy_2_and_chromatin_ingestion_are_numpy_1():
    # The green test stack (aria-env ≈ aria-rna-env) is numpy 2.x.
    rna = _conda_numpy("rna") or _pip_numpy("rna")
    assert rna is not None and rna[0] == 2
    # MACS2 (chromatin) and the kb ingestion stack require numpy 1.x — this is an
    # intentional pin, not a drift to chase.
    for env in ("chromatin", "ingestion"):
        v = _conda_numpy(env) or _pip_numpy(env)
        assert v is not None and v[0] == 1, f"aria-{env}-env should pin numpy 1.x"


def _yaml_numpy_spec(yaml_name: str) -> SpecifierSet | None:
    """Parse the numpy constraint from a conda env YAML, or None if absent."""
    path = _ENVS_DIR / yaml_name
    if not path.exists():
        return None
    m = _YAML_NUMPY_RE.search(path.read_text())
    if not m:
        return None
    spec = m.group(1).split("#")[0].strip()
    return SpecifierSet(spec) if spec else SpecifierSet("")


@pytest.mark.parametrize("env", sorted(_ENV_YAMLS))
def test_yaml_numpy_constraint_admits_locked_major(env):
    """S1: the YAML numpy constraint must ADMIT the env's locked numpy major and
    EXCLUDE the wrong major, so `conda env create -f <yaml>` reproduces the
    validated stack. This is the fence the RNA numpy<2 drift slipped past (the old
    test only checked the locks)."""
    locked = _conda_numpy(env) or _pip_numpy(env)
    assert locked is not None, f"aria-{env}-env has no locked numpy"
    major = locked[0]
    # Positive probe = the exact locked version (the YAML MUST admit what the lock
    # fixed). Negative probe = a realistic version of the opposite major the
    # correct constraint must exclude.
    good = Version(".".join(str(x) for x in locked))
    wrong = Version("2.0.0") if major == 1 else Version("1.26.4")

    for yaml_name in _ENV_YAMLS[env]:
        spec = _yaml_numpy_spec(yaml_name)
        assert spec is not None, f"{yaml_name} declares no numpy constraint"
        assert spec.contains(good), (
            f"{yaml_name} numpy constraint '{spec}' excludes the locked major "
            f"{major}.x — `conda env create` would diverge from the lock"
        )
        assert not spec.contains(wrong), (
            f"{yaml_name} numpy constraint '{spec}' admits the wrong major "
            f"(env is numpy {major}.x)"
        )
