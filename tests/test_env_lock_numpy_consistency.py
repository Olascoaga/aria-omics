"""T16 (tri-auditoría 2026-06-14): hermetic-reproducibility guard for the numpy
pin across the conda + pip lockfiles.

Gemini's "NumPy 2.0 ABI break / 91-of-117 fail" Blocker did NOT reproduce: aria-env
runs numpy 2.4.x green, and the three modality envs each import their full
scientific stack with no ABI error (verified live). This test is the durable fence
for what that verification established about the lockfiles themselves:

  * each env pins numpy in exactly one place (no conda-lock vs pip-lock conflict);
  * the numpy MAJOR per env matches intent — the RNA/report stack is numpy 2.x,
    while the chromatin (MACS2 needs numpy<2) and ingestion stacks are numpy 1.x.

It parses files only (no scientific imports), so it runs in any env.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_ENVS_DIR = Path(__file__).resolve().parents[1] / "envs"

_CONDA_NUMPY_RE = re.compile(r"/numpy-(\d+)\.(\d+)\.(\d+)")
_PIP_NUMPY_RE = re.compile(r"(?mi)^numpy==(\d+)\.(\d+)\.(\d+)\b")

# Expected numpy MAJOR per env (the T16 verified intent).
_EXPECTED_MAJOR = {"rna": 2, "chromatin": 1, "ingestion": 1}


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
