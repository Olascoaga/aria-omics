"""P2-2: container image identity in report provenance.

`collect_image_metadata` reads the image stamp (ARIA_IMAGE_* env) + the run-time
registry digest (ARIA_IMAGE_DIGEST). It must report honest nulls when ARIA is
NOT running inside a pinned image (never fabricate a digest), and fold into
`collect_version_metadata` so a report cites the image it ran in.
"""

import os

from aria.version import collect_image_metadata, collect_version_metadata

_IMAGE_ENV = (
    "ARIA_IMAGE_KIND", "ARIA_IMAGE_DIGEST", "ARIA_IMAGE_REVISION",
    "ARIA_IMAGE_ENV_SHA", "ARIA_IMAGE_REF", "ARIA_IMAGE_VALIDATION",
)


def _clear(monkeypatch):
    for k in _IMAGE_ENV:
        monkeypatch.delenv(k, raising=False)


def test_image_metadata_null_when_not_containerized(monkeypatch):
    _clear(monkeypatch)
    img = collect_image_metadata()
    assert img["containerized"] is False
    assert img["digest"] is None
    assert img["kind"] is None


def test_image_metadata_populated_from_env(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("ARIA_IMAGE_KIND", "rna")
    monkeypatch.setenv("ARIA_IMAGE_DIGEST",
                       "sha256:abc123def4567890")
    monkeypatch.setenv("ARIA_IMAGE_REVISION", "deadbeef")
    monkeypatch.setenv("ARIA_IMAGE_ENV_SHA", "feed0001")
    img = collect_image_metadata()
    assert img["containerized"] is True
    assert img["kind"] == "rna"
    assert img["digest"] == "sha256:abc123def4567890"
    assert img["revision"] == "deadbeef"
    assert img["env_lock_sha256"] == "feed0001"


def test_blank_env_is_treated_as_null(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("ARIA_IMAGE_KIND", "   ")
    img = collect_image_metadata()
    assert img["kind"] is None
    assert img["containerized"] is False


def test_version_metadata_includes_image_block(monkeypatch):
    _clear(monkeypatch)
    meta = collect_version_metadata()
    assert "image" in meta
    assert meta["image"]["containerized"] is False


def test_report_provenance_cites_image_digest():
    # P2-2 acceptance: the report cites the image digest in-container, and says
    # so honestly when not containerized. Render directly (no LLM init).
    import pytest
    pytest.importorskip("litellm")
    from aria.agents.narrative_agent import NarrativeAgent
    na = NarrativeAgent.__new__(NarrativeAgent)
    html = na._build_provenance_section(
        {"aria_version": "4.5.4",
         "image": {"containerized": True, "kind": "rna",
                   "digest": "sha256:abcd1234", "reference": None,
                   "revision": "deadbeef", "env_lock_sha256": "feed01",
                   "validation": None}},
        [], {})
    assert "image_digest" in html and "sha256:abcd1234" in html

    html2 = na._build_provenance_section(
        {"aria_version": "4.5.4",
         "image": {"containerized": False, "kind": None, "digest": None}},
        [], {})
    assert "not containerized" in html2
