"""S3 (pre-integration audit): capture the active conda env's lock hash in provenance.

`image.env_lock_sha256` only fills inside a pinned Docker image. Locally (how every
benchmark artifact is produced) it was always None — the tool stack that produced an
artifact was not recorded. collect_environment_metadata now resolves the active conda
env -> its committed lockfile -> sha256, with honest nulls when no lock is committed
(never a fabricated hash). This is the base of the Reproducibility Capsule (S14).
"""

from __future__ import annotations

import hashlib
from pathlib import Path

from aria.version import (
    collect_environment_metadata,
    collect_version_metadata,
)

_ENVS = Path(__file__).resolve().parents[1] / "envs"


def test_active_env_with_lock_hashes_it(monkeypatch):
    lock = _ENVS / "aria-rna-env.linux-64.lock"
    assert lock.exists(), "expected the validated RNA lock to be committed"
    monkeypatch.setenv("CONDA_PREFIX", "/anyprefix/envs/aria-rna-env")
    meta = collect_environment_metadata()
    assert meta["env_name"] == "aria-rna-env"
    assert meta["source"] == "conda_lock"
    assert meta["env_lock_file"] == "envs/aria-rna-env.linux-64.lock"
    expected = hashlib.sha256(lock.read_bytes()).hexdigest()
    assert meta["env_lock_sha256"] == expected


def test_active_env_without_lock_is_honest_null(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/anyprefix/envs/aria-atacseq-env")
    meta = collect_environment_metadata()
    # atacseq env has no committed .lock yet -> honest null, not a fabricated hash.
    if (_ENVS / "aria-atacseq-env.linux-64.lock").exists():
        assert meta["env_lock_sha256"] is not None
    else:
        assert meta["env_lock_sha256"] is None
        assert meta["source"] == "no_lock"
        assert meta["env_name"] == "aria-atacseq-env"


def test_no_conda_env_is_honest(monkeypatch):
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    meta = collect_environment_metadata()
    assert meta["env_lock_sha256"] is None
    assert meta["source"] == "no_conda_env"
    assert meta["env_name"] is None


def test_deterministic(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/x/envs/aria-rna-env")
    a = collect_environment_metadata()
    b = collect_environment_metadata()
    assert a == b and a["env_lock_sha256"] is not None


def test_version_metadata_includes_environment_block(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/x/envs/aria-rna-env")
    vm = collect_version_metadata()
    assert "environment" in vm
    env = vm["environment"]
    for key in ("conda_prefix", "env_name", "env_lock_file",
                "env_lock_sha256", "source"):
        assert key in env
    assert env["env_lock_sha256"] is not None  # rna lock resolved


def test_provenance_propagates_environment(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/x/envs/aria-rna-env")
    from aria.utils.provenance import collect_provenance
    prov = collect_provenance()
    assert "environment" in prov
    assert prov["environment"]["env_lock_sha256"] is not None
