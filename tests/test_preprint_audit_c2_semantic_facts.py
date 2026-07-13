"""C2 guards: typed semantic facts and explicit polarity for W-CLAIM."""

import pytest

from aria.agents.narrative.evidence_verifier import (
    EvidenceVerificationError,
    verify_block_claim_support,
)
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock, SemanticFact


def _block(claim: str, evidence: list[EvidenceItem], metrics=None) -> NarrativeBlock:
    return NarrativeBlock(
        id="c2.synthetic",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="Synthetic C2 block",
        status="success",
        confidence="medium",
        claim=claim,
        evidence=evidence,
        metrics=metrics or {},
    )


def _fact(subject, polarity="affirmed", aliases=None):
    return SemanticFact(
        subject=subject,
        subject_type="gene",
        predicate="differential_expression",
        polarity=polarity,
        aliases=aliases or [],
        source="synthetic_de",
    )


def test_semantic_fact_round_trips_through_evidence_item():
    item = EvidenceItem(
        label="top gene GeneAlpha",
        value=1.5,
        source="synthetic_de",
        facts=[_fact("GeneAlpha", aliases=["AliasAlpha"])],
    )
    restored = EvidenceItem.from_dict(item.to_dict())
    assert restored.facts == item.facts


def test_semantic_fact_rejects_invalid_polarity():
    with pytest.raises(ValueError, match="polarity"):
        SemanticFact(
            subject="GeneAlpha",
            subject_type="gene",
            predicate="differential_expression",
            polarity="unknown",
        )


@pytest.mark.parametrize("symbol", ["Gfap", "GeneAlpha", "GENE_BETA"])
def test_typed_affirmed_fact_supports_mixed_symbol_styles(symbol):
    block = _block(
        f"{symbol} was differentially expressed between conditions.",
        [EvidenceItem("DE result", "supported", facts=[_fact(symbol)])],
    )
    manifest = verify_block_claim_support(block, strict=True)
    assert manifest["status"] == "supported"
    facts = manifest["evidence_card"]["semantic_facts"]
    assert facts[0]["predicate"] == "differential_expression"


def test_explicit_alias_supports_claim_without_token_overlap():
    block = _block(
        "AliasAlpha was differentially expressed between conditions.",
        [EvidenceItem(
            "canonical result",
            "supported",
            facts=[_fact("GeneAlpha", aliases=["AliasAlpha"])],
        )],
    )
    assert verify_block_claim_support(block, strict=True)["status"] == "supported"


def test_absent_entity_is_unsupported_even_with_positive_de_count():
    block = _block(
        "GeneMissing was differentially expressed between conditions.",
        [EvidenceItem("DE genes", 4, source="synthetic_de")],
        metrics={"n_significant": 4},
    )
    manifest = verify_block_claim_support(block, strict=False)
    assert manifest["status"] == "unsupported"
    assert "semantic fact" in " ".join(x["reason"] for x in manifest["unsupported"])


def test_global_null_evidence_rejects_affirmed_entity():
    block = _block(
        "Gfap was differentially expressed between conditions.",
        [EvidenceItem("Result", "No significant genes were detected")],
    )
    with pytest.raises(EvidenceVerificationError, match="explicit null evidence"):
        verify_block_claim_support(block, strict=True)


def test_public_claim_compiler_withholds_semantic_contradiction():
    from aria.agents.narrative.claim_compiler import compile_public_claims

    block = _block(
        "Gfap was differentially expressed between conditions.",
        [EvidenceItem("Result", "No significant genes were detected")],
    )
    compiled = compile_public_claims([block], exp_ctx={})
    assert compiled.blocks == []
    assert compiled.claims == []
    assert compiled.withheld == [{
        "claim_id": "c2.synthetic",
        "status": "withheld",
        "reason": "unsupported claim or invalid narrative block",
    }]


def test_global_null_claim_is_supported_by_zero_count():
    block = _block(
        "No significant genes were detected between conditions.",
        [EvidenceItem("DE genes", 0, source="synthetic_de")],
        metrics={"n_significant": 0},
    )
    assert verify_block_claim_support(block, strict=True)["status"] == "supported"


def test_global_null_claim_conflicts_with_positive_fact():
    block = _block(
        "No significant genes were detected between conditions.",
        [EvidenceItem("top gene GeneAlpha", 1.5, source="synthetic_de")],
    )
    manifest = verify_block_claim_support(block, strict=False)
    assert manifest["status"] == "unsupported"
    assert "contradicts affirmative" in " ".join(
        x["reason"] for x in manifest["unsupported"]
    )


def test_conflicting_null_and_positive_evidence_never_supports_a_null_claim():
    block = _block(
        "No significant genes were detected between conditions.",
        [
            EvidenceItem("DE genes", 0, source="synthetic_de"),
            EvidenceItem("top gene GeneAlpha", 1.5, source="synthetic_de"),
        ],
    )
    assert verify_block_claim_support(block, strict=False)["status"] == "unsupported"


def test_entity_negation_requires_matching_polarity():
    supported = _block(
        "GeneAlpha was not differentially expressed between conditions.",
        [EvidenceItem("tested gene", "GeneAlpha", facts=[_fact("GeneAlpha", "negated")])],
    )
    assert verify_block_claim_support(supported, strict=True)["status"] == "supported"

    contradicted = _block(
        "GeneAlpha was not differentially expressed between conditions.",
        [EvidenceItem("tested gene", "GeneAlpha", facts=[_fact("GeneAlpha")])],
    )
    assert verify_block_claim_support(contradicted, strict=False)["status"] == "unsupported"


def test_top_gene_evidence_is_upgraded_to_typed_affirmed_fact():
    block = _block(
        "GeneAlpha was differentially expressed between conditions.",
        [EvidenceItem("top gene GeneAlpha", 1.5, source="synthetic_de")],
    )
    manifest = verify_block_claim_support(block, strict=True)
    assert manifest["evidence_card"]["semantic_facts"][0]["subject"] == "GeneAlpha"


def test_bulk_narrator_emits_typed_top_gene_fact():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    blocks = BulkRnaNarrator().collect("bulk_rna_agent", {"findings": {
        "contrasts": [{
            "name": "condition_a_vs_condition_b",
            "status": "success",
            "n_significant": 1,
            "n_upregulated": 1,
            "n_downregulated": 0,
            "top_genes": [{"symbol": "GeneAlpha", "log2fc": 1.5}],
        }],
    }})
    top = next(
        item for block in blocks for item in block.evidence
        if item.label == "top gene GeneAlpha"
    )
    assert top.facts == [SemanticFact(
        subject="GeneAlpha",
        subject_type="gene",
        predicate="differential_expression",
        polarity="affirmed",
        value=1.5,
        source="bulk_de",
    )]


def test_scrna_narrator_emits_typed_top_gene_fact():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    blocks = ScrnaNarrator()._pseudobulk_blocks({
        "pseudobulk_de": {
            "per_group": {"GroupA": {"per_comparison": {
                "condition_a_vs_condition_b": {
                    "status": "success",
                    "n_significant": 1,
                    "n_significant_global": 1,
                    "n_significant_local": 1,
                    "n_up": 1,
                    "n_down": 0,
                    "top_genes": [{"gene": "GeneAlpha", "log2fc": 1.5}],
                },
            }}},
        },
    })
    top = next(
        item for block in blocks for item in block.evidence
        if item.label == "top gene GeneAlpha"
    )
    assert top.facts == [SemanticFact(
        subject="GeneAlpha",
        subject_type="gene",
        predicate="differential_expression",
        polarity="affirmed",
        value=1.5,
        source="pseudobulk_de",
    )]
