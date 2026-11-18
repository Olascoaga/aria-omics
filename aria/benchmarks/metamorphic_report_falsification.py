"""600-Case Deterministic Claim-Boundary & Metamorphic Report Falsification Suite (ADR-083).

Replaces subjective small-panel human grading in the automated Preprint `v1`
gate (`c4_metamorphic_600_falsification`) by systematically attacking ARIA's
live `compile_public_claims` + `W-CLAIM` `SemanticFact` verifier + `W-LEDGER`
`RunLedger` + `C4` causal licensing + `C7` `SPECULATIVE` quarantine wall across:

  - 5 Modality Templates (`bulk_rna`, `scrna`, `scatac`, `cut_and_run`, `hic`)
  - 6 Adversarial Error Families:
      1. `numeric_or_threshold_tampering`
      2. `entity_or_contrast_swap`
      3. `sign_or_negation_inversion`
      4. `unexecuted_or_errored_ledger_stage`
      5. `unlicensed_causal_verb_injection`
      6. `speculative_or_exploratory_promotion`
  - 10 Deterministic Seeds (`0..9`)

Total = `300` invalid metamorphic mutations + `300` matched valid controls
(`600` total end-to-end compilation & post-compilation verification trials).
"""

from __future__ import annotations

from typing import Any

from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.hypothesis import assert_no_speculative_promotion
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock, SemanticFact

METAMORPHIC_SUITE_SCHEMA = "aria.metamorphic_report_falsification.v1"

MODALITY_TEMPLATES: tuple[dict[str, Any], ...] = (
    {
        "key": "bulk_rna",
        "modality": "bulk_rna",
        "ledger_modality": "bulk",
        "analysis": "differential_expression",
        "predicate": "differential_expression",
        "ledger_node": "ledger://bulk/differential_expression",
        "entity_prefix": "GENE",
    },
    {
        "key": "scrna",
        "modality": "scrna",
        "ledger_modality": "scRNA",
        "analysis": "pseudobulk_de",
        "predicate": "differential_expression",
        "ledger_node": "ledger://scRNA/pseudobulk_de",
        "entity_prefix": "SCGENE",
    },
    {
        "key": "scatac",
        "modality": "chromatin",
        "ledger_modality": "chromatin",
        "analysis": "differential_accessibility",
        "predicate": "differential_accessibility",
        "ledger_node": "ledger://chromatin/differential_accessibility",
        "entity_prefix": "PEAK",
    },
    {
        "key": "cut_and_run",
        "modality": "chromatin",
        "ledger_modality": "chromatin",
        "analysis": "differential_accessibility",
        "predicate": "differential_accessibility",
        "ledger_node": "ledger://chromatin/differential_accessibility",
        "entity_prefix": "TARGET",
    },
    {
        "key": "hic",
        "modality": "hic",
        "ledger_modality": "hic",
        "analysis": "differential_3d",
        "predicate": "differential_3d",
        "ledger_node": "ledger://hic/differential_3d",
        "entity_prefix": "LOCUS",
    },
)

ERROR_FAMILIES: tuple[str, ...] = (
    "numeric_or_threshold_tampering",
    "entity_or_contrast_swap",
    "sign_or_negation_inversion",
    "unexecuted_or_errored_ledger_stage",
    "unlicensed_causal_verb_injection",
    "speculative_or_exploratory_promotion",
)


def _valid_case(template: dict[str, Any], family_idx: int, seed: int) -> tuple[NarrativeBlock, dict[str, Any], dict[str, Any]]:
    count_val = 12 + seed * 3 + family_idx
    entity = f"{template['entity_prefix']}{100 + seed * 7 + family_idx}"
    block_id = f"{template['key']}.valid.f{family_idx}.s{seed}"
    fact = SemanticFact(
        subject=entity,
        predicate=template["predicate"],
        polarity="affirmed",
        subject_type="gene",
        value=count_val,
        source="agent_results",
    )
    block = NarrativeBlock(
        id=block_id,
        modality=template["modality"],
        analysis=template["analysis"],
        block_type="result",
        title=f"{template['key']} contrast summary",
        status="success",
        confidence="high",
        claim=(
            f"Observational contrast identified {count_val} significant features "
            f"associated with {entity}."
        ),
        evidence=[
            EvidenceItem(
                label="significant features",
                value=count_val,
                source="agent_results",
                facts=[fact],
            ),
            EvidenceItem(
                label="lead feature",
                value=entity,
                source="agent_results",
            ),
        ],
        metrics={"power": 0.85, "n_significant": count_val},
    )
    ledger = {
        "entries": [
            {
                "node_id": template["ledger_node"],
                "node_type": "analysis_run",
                "provenance": {
                    "kind": "structured_agent_result",
                    "finding_keys": [template["analysis"]],
                },
                "modality": template["ledger_modality"],
                "analysis": template["analysis"],
                "status": "ran",
            }
        ]
    }
    exp_ctx = {"organism": "Homo sapiens", "genome": "hg38", "design": {"interventional": False}}
    return block, ledger, exp_ctx


def _mutated_case(
    template: dict[str, Any], family: str, family_idx: int, seed: int
) -> tuple[NarrativeBlock, dict[str, Any], dict[str, Any]]:
    block, ledger, exp_ctx = _valid_case(template, family_idx, seed)
    count_val = int(block.metrics["n_significant"])
    entity = str(block.evidence[1].value)
    block.id = f"{template['key']}.invalid.{family}.s{seed}"

    if family == "numeric_or_threshold_tampering":
        tampered_count = count_val + 77 + seed
        block.claim = (
            f"Observational contrast identified {tampered_count} significant features "
            f"associated with {entity}."
        )
    elif family == "entity_or_contrast_swap":
        decoy_entity = f"UNATTESTED_{template['entity_prefix']}_{900 + seed}"
        block.claim = (
            f"Observational contrast identified {count_val} significant features "
            f"associated with {decoy_entity}."
        )
    elif family == "sign_or_negation_inversion":
        block.claim = (
            f"Observational contrast did not identify {count_val} significant features "
            f"and {entity} was not differentially regulated."
        )
    elif family == "unexecuted_or_errored_ledger_stage":
        ledger["entries"][0]["status"] = "not_run" if seed % 2 == 0 else "failed"
    elif family == "unlicensed_causal_verb_injection":
        causal_verbs = (
            "directly drives",
            "causally induces",
            "causes",
            "mechanistically drives",
            "directly activates",
        )
        verb = causal_verbs[seed % len(causal_verbs)]
        block.claim = (
            f"Contrast {verb} {count_val} downstream targets through {entity}."
        )
    elif family == "speculative_or_exploratory_promotion":
        block.id = f"hypothesis://{template['key']}/s{seed}"
        block.metadata["tier"] = "SPECULATIVE"
        block.metadata["promotable"] = False
    else:
        raise ValueError(f"Unknown error family: {family}")

    return block, ledger, exp_ctx


def _evaluate_trial(block: NarrativeBlock, ledger: dict[str, Any], exp_ctx: dict[str, Any]) -> bool:
    """Return True iff the claim survived all C1-C4 & C7 public-publication gates."""
    compilation = compile_public_claims([block], exp_ctx=exp_ctx, run_ledger=ledger)
    if not compilation.claims or compilation.withheld:
        return False
    try:
        assert_no_speculative_promotion(compilation.claims)
    except Exception:
        return False
    # Also enforce C4 causal language violation absence on published associative claims
    for published in compilation.claims:
        if published.get("language_violation"):
            return False
        if published.get("tier") == "SPECULATIVE" or str(published.get("claim_id", "")).startswith("hypothesis://"):
            return False
    return True


def run_metamorphic_falsification_suite(n_seeds: int = 10) -> dict[str, Any]:
    """Execute all `6 * 5 * n_seeds` invalid + valid control pairs (`600` default)."""
    per_family: dict[str, dict[str, int]] = {
        fam: {"n_invalid_tested": 0, "n_invalid_caught": 0, "n_valid_tested": 0, "n_valid_published": 0}
        for fam in ERROR_FAMILIES
    }
    per_modality: dict[str, dict[str, int]] = {
        tpl["key"]: {"n_invalid_tested": 0, "n_invalid_caught": 0, "n_valid_tested": 0, "n_valid_published": 0}
        for tpl in MODALITY_TEMPLATES
    }

    failures: list[dict[str, Any]] = []
    total_invalid = 0
    caught_invalid = 0
    total_valid = 0
    published_valid = 0

    for fam_idx, fam in enumerate(ERROR_FAMILIES):
        for tpl in MODALITY_TEMPLATES:
            mkey = tpl["key"]
            for seed in range(int(n_seeds)):
                # 1. Valid matched control MUST publish cleanly
                v_block, v_ledger, v_ctx = _valid_case(tpl, fam_idx, seed)
                v_pub = _evaluate_trial(v_block, v_ledger, v_ctx)
                total_valid += 1
                per_family[fam]["n_valid_tested"] += 1
                per_modality[mkey]["n_valid_tested"] += 1
                if v_pub:
                    published_valid += 1
                    per_family[fam]["n_valid_published"] += 1
                    per_modality[mkey]["n_valid_published"] += 1
                else:
                    failures.append({"kind": "valid_control_rejected", "family": fam, "modality": mkey, "seed": seed})

                # 2. Invalid metamorphic mutation MUST be caught/blocked
                i_block, i_ledger, i_ctx = _mutated_case(tpl, fam, fam_idx, seed)
                i_pub = _evaluate_trial(i_block, i_ledger, i_ctx)
                total_invalid += 1
                per_family[fam]["n_invalid_tested"] += 1
                per_modality[mkey]["n_invalid_tested"] += 1
                if not i_pub:
                    caught_invalid += 1
                    per_family[fam]["n_invalid_caught"] += 1
                    per_modality[mkey]["n_invalid_caught"] += 1
                else:
                    failures.append({"kind": "invalid_mutation_leaked", "family": fam, "modality": mkey, "seed": seed})

    false_pub_rate = (total_invalid - caught_invalid) / max(1, total_invalid)
    valid_ret_rate = published_valid / max(1, total_valid)
    all_pass = (caught_invalid == total_invalid) and (published_valid == total_valid)

    return {
        "schema_version": METAMORPHIC_SUITE_SCHEMA,
        "status": "pass" if all_pass else "fail",
        "all_gates_pass": all_pass,
        "n_error_families": len(ERROR_FAMILIES),
        "n_modality_templates": len(MODALITY_TEMPLATES),
        "n_seeds": int(n_seeds),
        "total_cases": total_invalid + total_valid,
        "n_invalid_mutations": total_invalid,
        "n_invalid_caught": caught_invalid,
        "false_publication_rate": round(false_pub_rate, 6),
        "n_valid_controls": total_valid,
        "n_valid_published": published_valid,
        "valid_retention_rate": round(valid_ret_rate, 6),
        "per_family": per_family,
        "per_modality": per_modality,
        "failures": failures,
    }
