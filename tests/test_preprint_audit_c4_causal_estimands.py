"""C4: causal language is licensed per verified estimand, never per study."""

import pytest

from aria.agents.narrative.claim_compiler import (
    annotate_claim_tiers,
    compile_claim_manifest,
)
from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
from aria.agents.narrative.narrators.scrna import ScrnaNarrator
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock


def _block(
    block_id: str,
    analysis: str,
    *,
    estimand_id: str | None = None,
    contrast: str | None = None,
) -> NarrativeBlock:
    metadata = {}
    if estimand_id is not None or contrast is not None:
        metadata["estimand"] = {
            "id": estimand_id,
            "contrast": contrast,
        }
    return NarrativeBlock(
        id=block_id,
        modality="test",
        analysis=analysis,
        block_type="result",
        title=block_id,
        status="success",
        confidence="medium",
        claim="The contrast had 10 significant features.",
        evidence=[EvidenceItem(label="significant features", value=10, source="test")],
        metrics={"n_significant": 10},
        metadata=metadata,
    )


def _context(*estimands: dict, **design_overrides) -> dict:
    design = {
        # The legacy study-wide marker may describe the experiment, but C4
        # forbids using it as claim-level causal authority.
        "interventional": True,
        "causal_estimands": list(estimands),
        "unresolved_confounding": False,
        "confounded_covariates": [],
    }
    design.update(design_overrides)
    return {"design": design}


def _verified(estimand_id: str, contrast: str) -> dict:
    return {
        "id": estimand_id,
        "contrast": contrast,
        "interventional": True,
        "randomized": True,
        "confounding_verified": True,
    }


def test_study_contract_without_block_estimand_does_not_license_causality():
    blocks = [
        _block("bulk.contrast.drug", "differential_expression"),
        _block("scrna.composition.age", "differential_abundance"),
        _block("scrna.trajectory", "trajectory"),
    ]

    annotate_claim_tiers(
        blocks,
        _context(_verified("drug", "drug")),
    )

    assert [b.metadata["claim"]["licensed_language"] for b in blocks] == [
        "associative",
        "associative",
        "associative",
    ]


def test_mixed_experiment_licenses_only_the_exact_verified_estimand():
    blocks = [
        _block(
            "bulk.contrast.drug_vs_vehicle",
            "differential_expression",
            estimand_id="drug_effect_on_expression",
            contrast="drug vs vehicle",
        ),
        _block(
            "scrna.composition.old_vs_young",
            "differential_abundance",
            estimand_id="age_association_with_abundance",
            contrast="old vs young",
        ),
        _block(
            "scrna.trajectory",
            "trajectory",
            estimand_id="drug_effect_on_expression",
            contrast="drug vs vehicle",
        ),
    ]

    annotate_claim_tiers(
        blocks,
        _context(_verified("drug_effect_on_expression", "drug vs vehicle")),
    )

    assert blocks[0].metadata["claim"]["tier"] == "causal_experimental"
    assert blocks[1].metadata["claim"]["licensed_language"] == "associative"
    # A trajectory is exploratory context, not the effect estimand tested by DE.
    assert blocks[2].metadata["claim"]["licensed_language"] == "associative"


@pytest.mark.parametrize(
    "contract",
    [
        {
            "id": "drug_effect",
            "contrast": "drug vs vehicle",
            "interventional": False,
            "randomized": True,
            "confounding_verified": True,
        },
        {
            "id": "drug_effect",
            "contrast": "drug vs vehicle",
            "interventional": True,
            "randomized": False,
            "confounding_verified": True,
        },
        {
            "id": "drug_effect",
            "contrast": "drug vs vehicle",
            "interventional": True,
            "randomized": True,
            "confounding_verified": False,
        },
    ],
)
def test_each_causal_prerequisite_is_explicit_and_fail_closed(contract):
    block = _block(
        "bulk.contrast.drug_vs_vehicle",
        "differential_expression",
        estimand_id="drug_effect",
        contrast="drug vs vehicle",
    )

    annotate_claim_tiers([block], _context(contract))

    assert block.metadata["claim"]["licensed_language"] == "associative"


def test_ambiguous_duplicate_contracts_and_global_confounding_fail_closed():
    contract = _verified("drug_effect", "drug vs vehicle")
    ambiguous = _block(
        "bulk.contrast.drug_vs_vehicle",
        "differential_expression",
        estimand_id="drug_effect",
        contrast="drug vs vehicle",
    )
    confounded = _block(
        "bulk.contrast.drug_vs_vehicle.2",
        "differential_expression",
        estimand_id="drug_effect",
        contrast="drug vs vehicle",
    )

    annotate_claim_tiers([ambiguous], _context(contract, dict(contract)))
    annotate_claim_tiers(
        [confounded],
        _context(
            contract,
            unresolved_confounding=True,
            confounded_covariates=["batch"],
        ),
    )

    assert ambiguous.metadata["claim"]["licensed_language"] == "associative"
    assert confounded.metadata["claim"]["licensed_language"] == "associative"


def test_manifest_persists_estimand_license_decision():
    block = _block(
        "bulk.contrast.drug_vs_vehicle",
        "differential_expression",
        estimand_id="drug_effect",
        contrast="drug vs vehicle",
    )
    annotate_claim_tiers(
        [block],
        _context(_verified("drug_effect", "drug vs vehicle")),
    )

    manifest = compile_claim_manifest(block)

    assert manifest["causal_license"] == {
        "licensed": True,
        "estimand_id": "drug_effect",
        "contrast": "drug vs vehicle",
        "reason": "verified_estimand_contract",
    }


def test_bulk_narrator_carries_structured_contrast_to_the_license_resolver():
    block = BulkRnaNarrator()._contrast_blocks({
        "contrasts": [{
            "name": "drug vs vehicle",
            "status": "success",
            "n_significant": 10,
            "n_upregulated": 6,
            "n_downregulated": 4,
        }],
    })[0]

    annotate_claim_tiers(
        [block],
        _context(_verified("drug_effect", "drug vs vehicle")),
    )

    assert block.metadata["estimand"] == {
        "id": None,
        "contrast": "drug vs vehicle",
    }
    assert block.metadata["claim"]["causal_license"]["estimand_id"] == "drug_effect"
    assert block.metadata["claim"]["licensed_language"] == "causal"


def test_scrna_narrator_keeps_mixed_abundance_contrasts_separate():
    blocks = ScrnaNarrator()._composition_blocks({
        "differential_abundance": {
            "per_comparison": {
                "drug_vs_vehicle": {
                    "status": "success",
                    "n_significant": 1,
                    "per_cell_type": [],
                },
                "old_vs_young": {
                    "status": "success",
                    "n_significant": 1,
                    "per_cell_type": [],
                },
            },
        },
    })

    annotate_claim_tiers(
        blocks,
        _context(_verified("drug_abundance_effect", "drug_vs_vehicle")),
    )

    by_contrast = {
        block.metadata["estimand"]["contrast"]: block.metadata["claim"]
        for block in blocks
    }
    assert by_contrast["drug_vs_vehicle"]["licensed_language"] == "causal"
    assert by_contrast["old_vs_young"]["licensed_language"] == "associative"
