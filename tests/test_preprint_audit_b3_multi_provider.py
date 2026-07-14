"""Preprint-readiness audit B3-multi (FASE 7 / Claim 6): reproducibility across
LLM providers and repetitions.

Claim 6 asserts that ARIA's public output is invariant under allowed LLM prose
variation: N providers x M repetitions over the SAME structured evidence must
yield identical deterministic statistics, identical public accepted claims, and
an equivalent reproducibility artifact. Only free narrative prose may differ.

This is the Fase A unit lane. It exercises the real claim compiler, run-ledger
linkage and methodology diff (``aria.benchmarks.multi_provider_repro``) with a
deterministic in-process prose-variation harness, so ordinary CI proves the
invariant without any billed live provider call. The heavy end-to-end lane over
``run_headless`` is a separate opt-in guard (Fase B).

The invariant rests on two structural facts, both checked below:

1. For a structured analysis, the RENDERED public prose is composed
   deterministically from the structured evidence and ignores the LLM-authored
   ``block.claim`` wording entirely, so provider prose variation is erased at the
   render boundary.
2. The public claim's evidence tier and ledger linkage are derived from that
   same structured evidence, never from prose, so ``diff_methodologies`` reports
   the cells identical regardless of which provider authored the prose.
"""
from __future__ import annotations

import copy
import json

import pytest

from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock, SemanticFact
from aria.benchmarks.multi_provider_repro import (
    ProseVariant,
    run_provider_matrix,
)


# ── Canonical structured payload (identical across every matrix cell) ─────────

_PROVENANCE = {
    "aria_version": "4.7.0",
    "git_commit": "b3multitestcommit",
    "git_sha": "b3multitestcommit",
    "git_dirty": False,
    "workflow_hash": "b3-multi-workflow-hash",
    "image_digest": None,
    "timestamp_utc": "2026-07-13T00:00:00+00:00",
}

_EXP_CTX = {
    "input_files": [{"path": "counts.tsv", "sha256": "0" * 64}],
    "design_intelligence": {
        "recommended": ["DESeq2 differential expression between conditions."],
        "optional": [],
    },
    "design": {},
}

# A bulk DE finding whose ``contrasts`` are truthy makes the run ledger record a
# typed, provenance-bearing ``ledger://bulk/differential_expression`` node, so a
# supported DE claim is eligible for public compilation.
_AGENT_RESULTS = {
    "bulk_rna_agent": {"findings": {"contrasts": [{"name": "cond_vs_ctrl"}]}}
}


def _canonical_blocks() -> list[NarrativeBlock]:
    """Fresh, structurally identical blocks for one matrix cell.

    A factory (not a shared template) because ``compile_public_claims`` mutates
    ``block.metadata`` in place; every cell must start from clean state.
    """
    return [
        NarrativeBlock(
            id="b3.bulk.de",
            modality="bulk RNA-seq",
            analysis="differential_expression",
            block_type="result",
            title="Differential expression (condition vs control)",
            status="success",
            confidence="medium",
            claim="GeneAlpha was differentially expressed between conditions.",
            evidence=[
                EvidenceItem("DE genes", 3, source="synthetic_de"),
                EvidenceItem(
                    "top gene GeneAlpha",
                    1.8,
                    source="synthetic_de",
                    facts=[
                        SemanticFact(
                            subject="GeneAlpha",
                            subject_type="gene",
                            predicate="differential_expression",
                            polarity="affirmed",
                            source="synthetic_de",
                        )
                    ],
                ),
            ],
            metrics={"n_significant": 3, "n_upregulated": 2, "n_downregulated": 1},
        )
    ]


# ── Prose transforms modelling how different providers word the SAME result ───

def _identity(text: str) -> str:
    return text


def _terse(text: str) -> str:
    return "GeneAlpha: DE between conditions."


def _verbose(text: str) -> str:
    return (
        "Our differential expression analysis clearly demonstrated that "
        "GeneAlpha was differentially expressed between the two conditions."
    )


def _causal_overreach(text: str) -> str:
    # Same supported subject (GeneAlpha, affirmed DE), but the wording asserts
    # causation and a speculative novel mechanism the evidence does not license.
    return (
        "GeneAlpha causally drives the phenotype and is the master regulator "
        "responsible for the condition; it was differentially expressed "
        "between conditions."
    )


_BENIGN_VARIANTS = [
    ProseVariant(provider="anthropic", claim_style=_identity,
                 free_narrative="Anthropic-style executive summary prose."),
    ProseVariant(provider="openai", claim_style=_terse,
                 free_narrative="OpenAI-style executive summary prose."),
    ProseVariant(provider="gemini", claim_style=_verbose,
                 free_narrative="Gemini-style executive summary prose."),
]


def _run(variants, repetitions=2):
    return run_provider_matrix(
        block_factory=_canonical_blocks,
        exp_ctx=copy.deepcopy(_EXP_CTX),
        agent_results=copy.deepcopy(_AGENT_RESULTS),
        provenance=copy.deepcopy(_PROVENANCE),
        variants=variants,
        repetitions=repetitions,
    )


# ── Guards ────────────────────────────────────────────────────────────────────

def test_baseline_publishes_the_supported_claim():
    """Sanity: the canonical payload actually yields a published public claim,
    otherwise an all-withheld matrix would be vacuously invariant."""
    result = _run(_BENIGN_VARIANTS)
    baseline = result.baseline
    assert baseline.published_claim_ids == ("b3.bulk.de",)
    assert baseline.withheld_claim_ids == ()
    tiers = {c["claim_id"]: c["tier"] for c in baseline.methodology["claims"]}
    assert tiers["b3.bulk.de"] == "associative"


def test_public_claims_invariant_across_providers_and_repetitions():
    """N providers x M repetitions -> every cell is identical to the baseline
    over tracked provenance, ledger, claims and calibration."""
    result = _run(_BENIGN_VARIANTS, repetitions=3)
    assert len(result.cells) == len(_BENIGN_VARIANTS) * 3
    assert result.invariant() is True
    for diff in result.pairwise_diffs():
        assert diff["identical"] is True, diff
    # Same published claim set and tiers in every cell.
    published = {c.published_claim_ids for c in result.cells}
    assert published == {("b3.bulk.de",)}
    tiers = {
        tuple(sorted((cl["claim_id"], cl["tier"])
                     for cl in c.methodology["claims"]))
        for c in result.cells
    }
    assert len(tiers) == 1


def test_rendered_public_prose_is_deterministic_regardless_of_claim_wording():
    """The rendered public sentence is composed from structured evidence and is
    byte-identical no matter how a provider worded ``block.claim``."""
    rendered = set()
    for style in (_identity, _terse, _verbose, _causal_overreach):
        block = _canonical_blocks()[0]
        block.claim = style(block.claim)
        rendered.add(compose_block_prose(block))
    assert len(rendered) == 1


def test_repetitions_are_byte_identical():
    """The same provider repeated M times yields byte-identical claim manifests
    (the compiler is a pure function of the structured evidence)."""
    result = _run(
        [ProseVariant(provider="anthropic", claim_style=_identity)],
        repetitions=4,
    )
    payloads = {
        json.dumps(c.methodology["claims"], sort_keys=True, default=str)
        for c in result.cells
    }
    assert len(payloads) == 1


def test_adversarial_prose_is_fail_closed_monotone_never_adds_or_elevates():
    """The public claim set is fail-closed monotone under prose variation: an
    adversarial provider that words the result as causal/speculative overreach
    can only ever REMOVE a claim (withhold it), never add a claim nor elevate a
    tier. ARIA never publishes more than the structured evidence licenses.

    This is the honest form of Claim 6: benign providers agree exactly (see
    ``test_public_claims_invariant_across_providers_and_repetitions``), while a
    misbehaving provider fails closed rather than smuggling a stronger claim.
    """
    adversarial = ProseVariant(provider="rogue", claim_style=_causal_overreach)
    result = _run(_BENIGN_VARIANTS + [adversarial], repetitions=1)

    baseline_published = set(result.baseline.published_claim_ids)
    assert baseline_published == {"b3.bulk.de"}

    rogue = next(c for c in result.cells if c.provider == "rogue")
    # Fail closed: the causal overreach is withheld, not published capped.
    assert set(rogue.published_claim_ids) <= baseline_published
    assert rogue.withheld_claim_ids == ("b3.bulk.de",)
    assert rogue.methodology["claims"] == []

    # No cell — benign or adversarial — ever publishes an elevated causal tier or
    # a claim id absent from the licensed baseline.
    for cell in result.cells:
        assert set(cell.published_claim_ids) <= baseline_published
        for claim in cell.methodology["claims"]:
            assert claim["tier"] != "causal_experimental"
            assert claim["licensed_language"] != "causal"


def test_matrix_detects_a_real_structural_divergence():
    """Negative control: the invariance check has teeth. If a cell's STRUCTURED
    evidence differs (no significant genes), the claim is withheld and the diff
    is not identical."""
    def _unsupported_blocks() -> list[NarrativeBlock]:
        # A positive DE count but NO semantic fact for the named gene: the
        # evidence does not support the entity, so the claim must be withheld.
        blocks = _canonical_blocks()
        blocks[0].claim = "GeneGhost was differentially expressed between conditions."
        blocks[0].evidence = [EvidenceItem("DE genes", 3, source="synthetic_de")]
        blocks[0].metrics = {"n_significant": 3}
        return blocks

    result = run_provider_matrix(
        block_factory=_canonical_blocks,
        exp_ctx=copy.deepcopy(_EXP_CTX),
        agent_results=copy.deepcopy(_AGENT_RESULTS),
        provenance=copy.deepcopy(_PROVENANCE),
        variants=_BENIGN_VARIANTS[:1],
        repetitions=1,
        extra_cells=[("unsupported_model", _unsupported_blocks)],
    )
    assert result.invariant() is False
    divergent = next(c for c in result.cells if c.provider == "unsupported_model")
    assert divergent.published_claim_ids == ()
    assert divergent.withheld_claim_ids == ("b3.bulk.de",)
