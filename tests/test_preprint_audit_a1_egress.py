"""Preprint-readiness audit — A1: air-gap must be a LIVE, per-call gate.

Two egress surfaces Codex flagged:
  - the LLM provider cached ``_air_gapped`` at construction (covered by the live
    ratchet in tests/test_preprint_audit_probes.py::test_probe_a1_...);
  - SetupAgent shelled out to curl/conda without consulting the gate.

This module locks the SetupAgent egress gate: under air-gap, the reference download
and the Miniforge install refuse WITHOUT shelling out (zero egress).

Tracker: memory/audit/ARIA_PLAN_AUDITORIA_preprint_journal_2026-07-09.md
"""
from __future__ import annotations

import os

import pytest

from aria.agents.setup_agent import SetupAgent
from aria.utils import privacy


class _Stub:
    """_download / _install_miniforge use no instance state, so a bare stub is a
    valid ``self`` and avoids constructing the full agent."""


@pytest.fixture
def _restore_airgap():
    prev = os.environ.get("ARIA_AIR_GAPPED")
    prev_reason = privacy._runtime_air_gapped_reason
    yield
    if prev is None:
        os.environ.pop("ARIA_AIR_GAPPED", None)
    else:
        os.environ["ARIA_AIR_GAPPED"] = prev
    privacy._runtime_air_gapped_reason = prev_reason


def test_setup_download_refuses_under_airgap_without_egress(tmp_path, monkeypatch,
                                                            _restore_airgap):
    import aria.agents.setup_agent as sa
    monkeypatch.setattr(sa.subprocess, "run",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("EGRESS: subprocess.run called")))
    privacy.enable_air_gapped_runtime(reason="a1_test")

    err = SetupAgent._download(_Stub(), "https://ftp.ensembl.org/x.fa.gz",
                               tmp_path / "x.fa.gz")
    assert err is not None and "AIR_GAPPED" in err
    assert not (tmp_path / "x.fa.gz").exists()


def test_setup_install_miniforge_refuses_under_airgap(monkeypatch, _restore_airgap):
    import aria.agents.setup_agent as sa
    monkeypatch.setattr(sa.subprocess, "run",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("EGRESS: subprocess.run called")))
    privacy.enable_air_gapped_runtime(reason="a1_test")

    res = SetupAgent._install_miniforge(_Stub())
    assert "AIR_GAPPED" in res


def test_setup_download_allowed_when_not_airgapped(monkeypatch, _restore_airgap):
    # Sanity: the gate does not block a normal (non-air-gapped) run — it proceeds
    # to the downloader loop (which we stub to a benign miss).
    os.environ.pop("ARIA_AIR_GAPPED", None)
    import aria.agents.setup_agent as sa
    monkeypatch.setattr(sa.shutil, "which", lambda *_a, **_k: None)  # no curl/wget
    err = SetupAgent._download(_Stub(), "https://example/x", "/tmp/nope")
    # Not the air-gap refusal — it fell through to the "no downloader" path.
    assert err is not None and "AIR_GAPPED" not in err
