"""Preprint-readiness audit — Q1 (A8): public docs must not overclaim.

Codex A8 flagged four overclaims in the public docs that runtime does not back:
  1. "every LLM call ... fixed seed"      -> seed only on backends that accept one
  2. air-gap "governs **all** egress/Done" -> Partial; per-run enforcement pending (A1/A3)
  3. GEO "all four modalities"             -> classification/routing, not general fetch (E6)
  4. ATAC "full publishable workflow"      -> footprinting is descriptive-only (B7)

This fence keeps the corrected, honest phrasing from silently regressing. It parses
files only (runs in any env). The deeper A8 fix (capability matrix derived from the
runtime registry) is tracked separately; this guard locks in the interim honesty.

Tracker: memory/audit/ARIA_PLAN_AUDITORIA_preprint_journal_2026-07-09.md
"""
from __future__ import annotations

from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_README = (_ROOT / "README.md").read_text(encoding="utf-8")
_VALIDATION = (_ROOT / "docs" / "validation_status.md").read_text(encoding="utf-8")
_ROADMAP = (_ROOT / "docs" / "workflows" / "chromatin_roadmap.md").read_text(encoding="utf-8")


def test_seed_claim_is_qualified_per_backend():
    # The bare universal claim must be gone from both public surfaces.
    assert "with a fixed seed, and each report" not in _README
    assert "`temperature=0` + fixed seed;" not in _VALIDATION
    # ... and replaced with the honest per-backend qualifier.
    assert "backends that accept one" in _README
    assert "seed-deterministic" in _README
    assert "only to backends that accept one" in _VALIDATION


def test_air_gap_claim_is_not_overstated():
    # No "Done" row that claims the flag governs ALL egress.
    assert "governs **all** egress" not in _VALIDATION
    # The row is honestly marked Partial with the known enforcement gap.
    assert "| Air-gapped mode | Partial |" in _VALIDATION
    assert "resolved at provider construction" in _VALIDATION


def test_geo_connector_does_not_claim_general_four_modality_fetch():
    # GEO routing is classification, not general reproducible retrieval.
    assert "fastq_pending" in _VALIDATION
    assert "not" in _VALIDATION.lower()  # sanity: the caveat clause exists
    assert "not general reproducible retrieval" in _ROADMAP


def test_atac_is_not_claimed_publication_grade_end_to_end():
    # No unqualified "full publishable workflow end-to-end".
    assert "full publishable workflow end-to-end" not in _ROADMAP
    assert "full publishable ATAC workflow" not in _ROADMAP
    # Footprinting is explicitly called descriptive-only (B7).
    assert "FDR-controlled significance" in _ROADMAP
    assert "descriptive" in _ROADMAP.lower()
