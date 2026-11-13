"""ADR-079 Track 2 — paired same-cell multiome integration promoted to production.

Pure-Python governance tests:
  1. `multiome_contracts`: multi-donor paired cohorts, mixed paired+unpaired,
     ATAC-only `.h5mu`, and raw Multiome FASTQ manifests (fail-closed paths).
  2. `IntegrationAuditAgent` / `build_capability_matrix`: production default is
     green/auto; measured risks downgrade to yellow/requires_ack; non-same-cell
     inputs are red/blocked for integration only; beta override is preserved.
  3. `registry_integrity`: production invariants are hard errors.
  4. `script_contracts`: all 8 integration contracts at production.
  5. Run ledger + `IntegrationNarrator`: artifact-gated preprocessing node,
     batch-harmonized WNN surfacing, zero causal verbs.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import pytest

from aria.agents.modality_audit import (
    IntegrationAuditAgent,
    _integration_production_readiness,
    build_capability_matrix,
)
from aria.agents.narrative.narrators.integration import IntegrationNarrator
from aria.agents.narrative.run_ledger import build_run_ledger, node_id_for
from aria.agents.narrative.validators import find_causal_language
from aria.utils.multiome_contracts import (
    enforce_same_cell_integration_gate,
    infer_multiome_contract,
    validate_multiome_contract,
)

ROOT = Path(__file__).resolve().parents[1]

PRODUCTION_VALIDATION: dict[str, Any] = {
    "level": "production",
    "dispatch_enabled": True,
    "requires_ack": False,
    "scope": "paired_same_cell_rna_atac",
    "production_promoted": True,
    "promotion_policy_adr": "ADR-056",
    "promotion_adr": "ADR-079",
    "reason": "Paired same-cell RNA+ATAC multiome integration (ADR-079).",
}

BETA_VALIDATION: dict[str, Any] = {
    "level": "beta",
    "dispatch_enabled": True,
    "requires_ack": True,
    "scope": "paired_same_cell_rna_atac",
    "production_promoted": False,
    "promotion_policy_adr": "ADR-056",
    "promotion_adr": "ADR-071",
}

SINGLE_ASSAY_VALIDATION: dict[str, dict[str, Any]] = {
    m: {"level": "production", "dispatch_enabled": True}
    for m in ("scRNA", "bulk_RNA", "bulk_RNA_raw", "scATAC", "bulk_ATAC")
}


def _paired_rec(path: str, fmt: str = "h5mu", paired: bool = True) -> dict[str, Any]:
    return {
        "path": path,
        "modality": "scATAC",
        "confidence": "high",
        "evidence": {
            "format": fmt,
            "modalities": ["rna", "atac"] if paired else ["atac"],
            "paired_rna": paired,
        },
    }


def _paired_ctx(*paths_fmt: tuple[str, str], **extra: Any) -> dict[str, Any]:
    paths = [p for p, _ in paths_fmt]
    ctx: dict[str, Any] = {
        "genome": "hg38",
        "organism": "Homo sapiens",
        "modalities": {"scATAC": list(paths)},
        "assay_detections": [_paired_rec(p, f) for p, f in paths_fmt],
    }
    ctx.update(extra)
    return ctx


SINGLE_H5MU = _paired_ctx(("/virtual/pbmc_paired.h5mu", "h5mu"))


def _readiness(ctx: dict[str, Any], validation: dict[str, Any] | None = PRODUCTION_VALIDATION):
    matrix = build_capability_matrix(
        ctx,
        modality_validation=SINGLE_ASSAY_VALIDATION,
        integration_validation=validation,
    )
    return matrix, matrix["integration_readiness"]


# ── 1. multiome_contracts ────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "fmts",
    [
        ("h5mu", "h5mu"),
        ("10x_multiome_h5", "10x_multiome_h5"),
        ("10x_multiome_mtx", "10x_multiome_mtx"),
        ("h5mu", "10x_multiome_h5", "10x_multiome_mtx"),
    ],
)
def test_multi_donor_paired_cohort_is_allowed(fmts):
    paths_fmt = tuple(
        (f"/virtual/donor{i}/{'outs' if f != 'h5mu' else 'mu'}.{f}", f)
        for i, f in enumerate(fmts, start=1)
    )
    ctx = _paired_ctx(*paths_fmt)
    contract = infer_multiome_contract(ctx)
    assert contract["object_type"] == "paired_mudata"
    assert contract["same_cell"] is True
    assert contract["cohort_merge_required"] is True
    assert contract["n_containers"] == len(fmts)
    assert contract["evidence_paths"] == [p for p, _ in paths_fmt]

    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is True
    assert gate["blockers"] == []


def test_single_paired_container_is_not_a_cohort():
    contract = infer_multiome_contract(SINGLE_H5MU)
    assert contract["cohort_merge_required"] is False
    assert enforce_same_cell_integration_gate(SINGLE_H5MU)["allowed"] is True


def test_cell_ranger_arc_raw_sibling_does_not_create_a_cohort():
    ctx = _paired_ctx(
        ("/virtual/run/outs/filtered_feature_bc_matrix.h5", "10x_multiome_h5"),
        ("/virtual/run/outs/raw_feature_bc_matrix.h5", "10x_multiome_h5"),
    )
    contract = infer_multiome_contract(ctx)
    assert contract["n_containers"] == 1
    assert contract["cohort_merge_required"] is False
    assert enforce_same_cell_integration_gate(ctx)["allowed"] is True


@pytest.mark.parametrize("unpaired_key", ["scRNA", "scATAC"])
def test_paired_container_mixed_with_unpaired_files_fails_closed(unpaired_key):
    ctx = _paired_ctx(("/virtual/donor1.h5mu", "h5mu"))
    ctx["modalities"].setdefault(unpaired_key, []).append(
        f"/virtual/unpaired_{unpaired_key}.h5ad"
    )
    contract = infer_multiome_contract(ctx)
    assert contract["same_cell"] is False
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert any(
        b["check"] == "multiome_mixed_paired_and_unpaired_containers"
        for b in gate["blockers"]
    )


def test_fragments_sidecar_does_not_count_as_unpaired_container():
    ctx = _paired_ctx(("/virtual/donor1.h5mu", "h5mu"))
    ctx["modalities"]["scATAC"].append("/virtual/atac_fragments.tsv.gz")
    assert enforce_same_cell_integration_gate(ctx)["allowed"] is True


def test_atac_only_h5mu_still_fails_closed():
    ctx = {
        "genome": "hg38",
        "modalities": {"scATAC": ["/virtual/atac_only.h5mu"]},
        "assay_detections": [_paired_rec("/virtual/atac_only.h5mu", paired=False)],
    }
    contract = infer_multiome_contract(ctx)
    assert contract["same_cell"] is False
    assert contract["object_type"] != "paired_mudata"
    assert enforce_same_cell_integration_gate(ctx)["allowed"] is False


def test_explicit_contract_listing_several_containers_without_cohort_flag_fails():
    ctx = {
        "genome": "hg38",
        "modalities": {"scATAC": ["/virtual/a.h5mu", "/virtual/b.h5mu"]},
        "multiome_contract": {
            "object_type": "paired_mudata",
            "same_cell": True,
            "cell_namespace": "shared_mudata_obs",
            "atac_feature_space": "peak_matrix",
            "evidence_paths": ["/virtual/a.h5mu", "/virtual/b.h5mu"],
        },
    }
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert any(
        b["check"] == "multiome_multiple_mudata_unsupported" for b in gate["blockers"]
    )
    ctx["multiome_contract"]["cohort_merge_required"] = True
    assert enforce_same_cell_integration_gate(ctx)["allowed"] is True


@pytest.mark.parametrize(
    "manifest,namespace",
    [
        (
            {
                "platform": "10x_multiome_arc_v1",
                "rna_fastqs": ["/raw/S1_GEX_R1.fastq.gz", "/raw/S1_GEX_R2.fastq.gz"],
                "atac_fastqs": ["/raw/S1_ATAC_R1.fastq.gz", "/raw/S1_ATAC_R2.fastq.gz"],
            },
            "10x_arc_whitelist_pair",
        ),
        (
            {
                "platform": "share_seq",
                "barcode_map_path": "/raw/barcode_map.tsv",
                "samples": [
                    {"sample_id": "d1", "rna_fastqs": ["/r/d1_rna.fq.gz"],
                     "atac_fastqs": ["/r/d1_atac.fq.gz"]},
                    {"sample_id": "d2", "rna_fastqs": ["/r/d2_rna.fq.gz"],
                     "atac_fastqs": ["/r/d2_atac.fq.gz"]},
                ],
            },
            "barcode_map_tsv",
        ),
        (
            {
                "platform": "dogma_seq",
                "rna_fastqs": ["/r/rna.fq.gz"],
                "atac_fastqs": ["/r/atac.fq.gz"],
            },
            "shared_barcodes",
        ),
        (
            {
                "gex_whitelist_path": "/wl/gex.txt",
                "atac_whitelist_path": "/wl/atac.txt",
                "rna_fastqs": ["/r/rna.fq.gz"],
                "atac_fastqs": ["/r/atac.fq.gz"],
            },
            "10x_arc_whitelist_pair",
        ),
    ],
)
def test_raw_multiome_fastq_manifest_passes_when_complete(manifest, namespace):
    ctx = {"genome": "hg38", "multiome_fastq_manifest": manifest}
    contract = infer_multiome_contract(ctx)
    assert contract["object_type"] == "raw_multiome_fastq"
    assert contract["same_cell"] is True
    assert contract["cell_namespace"] == namespace
    assert contract["atac_feature_space"] == "fragments"
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is True, gate["blockers"]


@pytest.mark.parametrize(
    "manifest,code",
    [
        (
            {"platform": "10x_multiome_arc_v1", "rna_fastqs": ["/r/rna.fq.gz"]},
            "multiome_raw_fastq_incomplete",
        ),
        (
            {
                "platform": "10x_multiome_arc_v1",
                "samples": [
                    {"sample_id": "d1", "rna_fastqs": ["/r/a"], "atac_fastqs": ["/r/b"]},
                    {"sample_id": "d2", "rna_fastqs": ["/r/c"], "atac_fastqs": []},
                ],
            },
            "multiome_raw_fastq_incomplete",
        ),
        (
            {
                "platform": "10x_multiome_arc_v1",
                "same_cell": False,
                "rna_fastqs": ["/r/rna.fq.gz"],
                "atac_fastqs": ["/r/atac.fq.gz"],
            },
            "multiome_not_same_cell",
        ),
        (
            # Unknown chemistry, no barcode translation: never assumed same-cell.
            {"rna_fastqs": ["/r/rna.fq.gz"], "atac_fastqs": ["/r/atac.fq.gz"]},
            "multiome_not_same_cell",
        ),
    ],
)
def test_raw_multiome_fastq_manifest_fails_closed(manifest, code):
    ctx = {"genome": "hg38", "multiome_fastq_manifest": manifest}
    contract = infer_multiome_contract(ctx)
    assert contract["same_cell"] is False
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert any(b["check"] == code for b in gate["blockers"])
    assert validate_multiome_contract(ctx)["status"] == "red"


def test_explicit_raw_contract_is_revalidated_against_its_fastqs():
    ctx = {
        "genome": "hg38",
        "multiome_contract": {
            "object_type": "raw_multiome_fastq",
            "same_cell": True,
            "cell_namespace": "barcode_map_tsv",
            "barcode_map": "/raw/map.tsv",
            "rna_fastqs": ["/r/rna.fq.gz"],
            "atac_fastqs": [],
        },
    }
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert any(
        b["check"] == "multiome_raw_fastq_incomplete" for b in gate["blockers"]
    )
    ctx["multiome_contract"]["atac_fastqs"] = ["/r/atac.fq.gz"]
    assert enforce_same_cell_integration_gate(ctx)["allowed"] is True


# ── 2. IntegrationAuditAgent / capability matrix ─────────────────────────────


def test_production_readiness_helper_mirrors_atac_helper():
    prod = _integration_production_readiness()
    assert prod["tier"] == "production"
    assert prod["requires_ack"] is False
    assert prod["production_promoted"] is True
    assert prod["governance_blockers"] == []
    assert prod["promotion_adr"] == "ADR-079"
    assert prod["caveated_sub_layers"] == {
        "peak_to_gene": "production_associative",
        "eregulons": "production_associative",
        "mofa_factors": "production_associative",
    }
    beta = _integration_production_readiness(tier_override="beta")
    assert beta["production_promoted"] is False
    assert beta["requires_ack"] is True
    assert beta["governance_blockers"]


def _assert_green_auto(ready: dict[str, Any]) -> None:
    assert ready["status"] == "green"
    assert ready["dispatch_policy"] == "auto"
    assert ready["requires_ack"] is False
    assert ready["requires_ack_for_dispatch"] is False
    assert ready["production_promoted"] is True
    assert ready["promotion_adr"] == "ADR-079"
    assert ready["promotion_policy_adr"] == "ADR-056"
    assert ready["governance_blockers"] == []
    assert ready["findings"] == []
    assert ready["caveated_sub_layers"]["peak_to_gene"] == "production_associative"


def test_production_single_donor_paired_multiome_is_green_auto():
    matrix, ready = _readiness(SINGLE_H5MU)
    _assert_green_auto(ready)
    assert ready["notes"] == []
    assert matrix["status"] == "green"


def test_production_is_the_default_when_no_validation_is_passed():
    _, ready = _readiness(SINGLE_H5MU, validation=None)
    _assert_green_auto(ready)


def test_production_multi_donor_cohort_is_green_with_advisory_note():
    ctx = _paired_ctx(
        ("/virtual/donor1.h5mu", "h5mu"),
        ("/virtual/donor2.h5mu", "h5mu"),
        n_shared_cells=8000,
    )
    _, ready = _readiness(ctx)
    _assert_green_auto(ready)
    assert any("sweep-line consensus ATAC peaks" in n for n in ready["notes"])
    assert ready["checks"]["donor_design"]["multi_sample_cohort"] is True


def test_production_raw_multiome_fastq_is_green_with_note():
    ctx = {
        "genome": "hg38",
        "multiome_fastq_manifest": {
            "platform": "10x_multiome_arc_v1",
            "rna_fastqs": ["/r/rna.fq.gz"],
            "atac_fastqs": ["/r/atac.fq.gz"],
        },
    }
    _, ready = _readiness(ctx)
    _assert_green_auto(ready)
    assert any("Raw Multiome FASTQs" in n for n in ready["notes"])


def _assert_yellow_ack(ready: dict[str, Any], check: str) -> None:
    assert ready["status"] == "yellow"
    assert ready["dispatch_policy"] == "requires_ack"
    assert ready["requires_ack_for_dispatch"] is True
    assert ready["production_promoted"] is True
    assert any(f["check"] == check for f in ready["findings"])


@pytest.mark.parametrize(
    "extra",
    [
        {"n_cells": 150},
        {"n_shared_cells": 199},
        {"design": {"n_cells": 120}},
    ],
)
def test_low_cell_count_requires_ack(extra):
    ctx = {**copy.deepcopy(SINGLE_H5MU), **extra}
    matrix, ready = _readiness(ctx)
    _assert_yellow_ack(ready, "integration_low_cell_count")
    assert matrix["status"] == "yellow"


def test_two_hundred_cells_is_enough():
    _, ready = _readiness({**copy.deepcopy(SINGLE_H5MU), "n_shared_cells": 200})
    _assert_green_auto(ready)


@pytest.mark.parametrize(
    "contract_extra",
    [{"allow_partial_overlap": True}, {"shared_barcode_fraction": 0.35}],
)
def test_partial_barcode_overlap_requires_ack(contract_extra):
    ctx = {
        "genome": "hg38",
        "modalities": {"scRNA": ["/virtual/rna.h5ad"], "scATAC": ["/virtual/atac.h5ad"]},
        "multiome_contract": {
            "object_type": "explicit_same_cell",
            "same_cell": True,
            "cell_namespace": "shared_10x_barcodes",
            "atac_feature_space": "peak_matrix",
            **contract_extra,
        },
    }
    _, ready = _readiness(ctx)
    _assert_yellow_ack(ready, "integration_partial_barcode_overlap")


def test_one_donor_per_condition_batch_confounded_cohort_requires_ack():
    ctx = _paired_ctx(
        ("/virtual/ctrl_donor1.h5mu", "h5mu"),
        ("/virtual/stim_donor2.h5mu", "h5mu"),
        design={
            "groups": {"ctrl": ["ctrl_donor1"], "stim": ["stim_donor2"]},
            "batch_map": {"ctrl_donor1": "run_A", "stim_donor2": "run_B"},
        },
    )
    _, ready = _readiness(ctx)
    _assert_yellow_ack(ready, "integration_donor_confounded_design")


def test_balanced_cohort_is_not_flagged_as_confounded():
    ctx = _paired_ctx(
        ("/virtual/c1.h5mu", "h5mu"),
        ("/virtual/c2.h5mu", "h5mu"),
        ("/virtual/s1.h5mu", "h5mu"),
        ("/virtual/s2.h5mu", "h5mu"),
        design={
            "groups": {"ctrl": ["c1", "c2"], "stim": ["s1", "s2"]},
            "batch_map": {"c1": "A", "s1": "A", "c2": "B", "s2": "B"},
        },
    )
    _, ready = _readiness(ctx)
    _assert_green_auto(ready)


def test_unpaired_scrna_scatac_is_blocked_while_single_lanes_are_unaffected():
    ctx = {
        "genome": "hg38",
        "modalities": {"scRNA": ["/virtual/rna.h5ad"], "scATAC": ["/virtual/atac.h5ad"]},
    }
    matrix, ready = _readiness(ctx)
    assert ready["status"] == "red"
    assert ready["dispatch_policy"] == "blocked"
    assert ready["requires_ack_for_dispatch"] is False
    assert ready["same_cell_gate_allowed"] is False
    assert matrix["dispatch"]["blocked"] == []
    assert matrix["cards"]["scATAC"]["status"] == "green"
    # Integration readiness never changes the per-assay cards.
    for modality in ("scRNA", "scATAC"):
        solo, _ = _readiness({"genome": "hg38", "modalities": {modality: ctx["modalities"][modality]}})
        assert matrix["cards"][modality]["status"] == solo["cards"][modality]["status"]


def test_bulk_rna_plus_bulk_atac_is_blocked_while_cards_stay_green():
    ctx = {
        "genome": "hg38",
        "modalities": {"bulk_RNA": ["/virtual/counts.tsv"], "bulk_ATAC": ["/virtual/s1.bam"]},
    }
    matrix, ready = _readiness(ctx)
    assert ready["status"] == "red"
    assert ready["dispatch_policy"] == "blocked"
    assert ready["requires_ack_for_dispatch"] is False
    assert matrix["cards"]["bulk_RNA"]["status"] == "green"
    assert matrix["cards"]["bulk_ATAC"]["status"] == "green"
    assert matrix["dispatch"]["blocked"] == []
    assert matrix["status"] == "green"


def test_single_assay_only_is_not_applicable():
    _, ready = _readiness({"genome": "hg38", "modalities": {"scATAC": ["/v/a.h5ad"]}})
    assert ready["status"] == "not_applicable"
    assert ready["dispatch_policy"] == "not_applicable"


def test_dispatch_disabled_or_scaffold_is_blocked():
    _, ready = _readiness(SINGLE_H5MU, {**PRODUCTION_VALIDATION, "dispatch_enabled": False})
    assert ready["status"] == "red" and ready["dispatch_policy"] == "blocked"
    _, ready = _readiness(SINGLE_H5MU, {"level": "scaffold", "dispatch_enabled": True})
    assert ready["status"] == "red" and ready["dispatch_policy"] == "blocked"
    assert ready["production_promoted"] is False


def test_beta_override_preserves_legacy_ack_behaviour():
    matrix, ready = _readiness(SINGLE_H5MU, BETA_VALIDATION)
    assert ready["tier"] == "beta"
    assert ready["status"] == "yellow"
    assert ready["dispatch_policy"] == "requires_ack"
    assert ready["requires_ack"] is True
    assert ready["requires_ack_for_dispatch"] is True
    assert ready["production_promoted"] is False
    assert ready["promotion_adr"] == "ADR-071"
    assert ready["governance_blockers"]
    assert any(f["check"] == "integration_beta_ack_required" for f in matrix["findings"])
    assert matrix["status"] == "yellow"


def test_production_tier_with_ack_flag_behaves_as_ack_gated():
    _, ready = _readiness(SINGLE_H5MU, {**PRODUCTION_VALIDATION, "requires_ack": True})
    assert ready["status"] == "yellow"
    assert ready["requires_ack_for_dispatch"] is True
    assert ready["production_promoted"] is False


def test_integration_audit_agent_direct_call_matches_matrix():
    multiome = validate_multiome_contract(SINGLE_H5MU)
    direct = IntegrationAuditAgent().audit(SINGLE_H5MU, multiome, PRODUCTION_VALIDATION)
    _, via_matrix = _readiness(SINGLE_H5MU)
    assert direct["status"] == via_matrix["status"] == "green"


# ── 3. registry_integrity ────────────────────────────────────────────────────


def _registry_errors(monkeypatch, validation: dict[str, Any]) -> list[str]:
    from aria.utils.registry_integrity import check_registry_integrity

    monkeypatch.setattr(
        "aria.agents.orchestrator_agent.INTEGRATION_VALIDATION", validation
    )
    _stub_new_integration_scripts(monkeypatch)
    return [i.code for i in check_registry_integrity() if i.severity == "error"]


_NEW_SCRIPTS = (
    "aria/scripts/integration_assemble_multiome.py",
    "aria/scripts/integration_merge_cohort.py",
)


def _stub_new_integration_scripts(monkeypatch) -> None:
    """Track 1 lands the two intake scripts; stub their existence until merge."""
    original_exists = Path.exists

    def _exists(self: Path, *args: Any, **kwargs: Any) -> bool:
        if any(str(self).endswith(s) for s in _NEW_SCRIPTS):
            return True
        return original_exists(self, *args, **kwargs)

    monkeypatch.setattr(Path, "exists", _exists)


def test_valid_production_integration_tier_passes_registry(monkeypatch):
    assert _registry_errors(monkeypatch, dict(PRODUCTION_VALIDATION)) == []


def test_current_registry_tier_passes(monkeypatch):
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION

    assert _registry_errors(monkeypatch, dict(INTEGRATION_VALIDATION)) == []


@pytest.mark.parametrize(
    "mutation,code",
    [
        ({"scope": "unpaired_any"}, "integration_production_scope_invalid"),
        ({"scope": "bulk_cross_modal"}, "integration_production_scope_invalid"),
        ({"dispatch_enabled": False}, "integration_production_dispatch_disabled"),
        ({"production_promoted": False}, "integration_production_flag_invalid"),
        ({"requires_ack": True}, "integration_production_ack_flag_invalid"),
        ({"promotion_adr": None}, "integration_production_missing_adr"),
        ({"promotion_adr": ""}, "integration_production_missing_adr"),
    ],
)
def test_malformed_production_tier_is_a_hard_error(monkeypatch, mutation, code):
    errors = _registry_errors(monkeypatch, {**PRODUCTION_VALIDATION, **mutation})
    assert code in errors, errors


@pytest.mark.parametrize(
    "validation,code",
    [
        ({"level": "scaffold", "dispatch_enabled": True}, "scaffold_dispatch_enabled"),
        ({**BETA_VALIDATION, "requires_ack": False}, "integration_beta_missing_ack_gate"),
        ({**BETA_VALIDATION, "production_promoted": True},
         "integration_beta_production_flag_invalid"),
        ({**BETA_VALIDATION, "scope": "unpaired_any"}, "integration_beta_scope_invalid"),
    ],
)
def test_scaffold_and_beta_registry_branches_are_intact(monkeypatch, validation, code):
    assert code in _registry_errors(monkeypatch, validation)


def test_integration_never_enters_modality_validation():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    assert "integration" not in MODALITY_VALIDATION
    assert len(MODALITY_VALIDATION) == 9


# ── 4. script_contracts ──────────────────────────────────────────────────────


_INTEGRATION_SCRIPTS = (
    "aria/scripts/integration_assemble_multiome.py",
    "aria/scripts/integration_merge_cohort.py",
    "aria/scripts/integration_paired_qc.py",
    "aria/scripts/integration_wnn.py",
    "aria/scripts/integration_wnn_diagnostics.py",
    "aria/scripts/integration_mofa.py",
    "aria/scripts/integration_peak2gene.py",
    "aria/scripts/integration_eregulons.py",
)


def test_all_integration_script_contracts_are_production():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    for key in _INTEGRATION_SCRIPTS:
        assert key in SCRIPT_CONTRACTS, key
        assert SCRIPT_CONTRACTS[key].validation_level == "production", key


def test_new_and_extended_contract_fields():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    def names(key: str, attr: str) -> set[str]:
        return {f.name for f in getattr(SCRIPT_CONTRACTS[key], attr)}

    assert {"batch_key", "harmonize_batches", "resume"} <= names(
        "aria/scripts/integration_wnn.py", "inputs")
    assert {"batch_key", "condition_col", "resume"} <= names(
        "aria/scripts/integration_mofa.py", "inputs")
    for key in ("integration_paired_qc.py", "integration_peak2gene.py",
                "integration_eregulons.py"):
        assert "resume" in names(f"aria/scripts/{key}", "inputs")
    assert {"output_path", "n_shared_cells", "n_genes", "n_peaks"} == names(
        "aria/scripts/integration_assemble_multiome.py", "success_outputs")
    assert {"output_path", "n_samples", "n_cells_total", "n_genes",
            "n_consensus_peaks"} == names(
        "aria/scripts/integration_merge_cohort.py", "success_outputs")
    merge_inputs = {f.name: f for f in SCRIPT_CONTRACTS[
        "aria/scripts/integration_merge_cohort.py"].inputs}
    assert merge_inputs["inputs"].type == "list[path]"
    assert merge_inputs["inputs"].required is True
    assert merge_inputs["inputs"].path_must_exist is True
    # All new optional inputs are optional.
    for key in ("aria/scripts/integration_assemble_multiome.py",):
        assert all(not f.required for f in SCRIPT_CONTRACTS[key].inputs)


def test_merge_cohort_contract_rejects_missing_inputs(tmp_path):
    from aria.utils.script_contracts import contract_for_script

    contract = contract_for_script("aria/scripts/integration_merge_cohort.py")
    assert contract is not None
    issues = contract.validate_params({"inputs": [str(tmp_path / "missing.h5mu")]})
    assert any(i.severity == "error" for i in issues)
    assert any(i.severity == "error" for i in contract.validate_params({"inputs": []}))
    real = tmp_path / "d1.h5mu"
    real.write_bytes(b"")
    issues_ok = contract.validate_params({"inputs": [str(real)], "resume": True})
    assert [i for i in issues_ok if i.severity == "error"] == []
    # Success output must point at a produced container.
    bad = contract.validate_result({
        "status": "success", "output_path": str(tmp_path / "nope.h5mu"),
        "n_samples": 2, "n_cells_total": 10, "n_genes": 5, "n_consensus_peaks": 7,
    })
    assert any(i.field == "output_path" for i in bad)


def test_agent_script_contracts_for_integration_agent(monkeypatch):
    from aria.utils.registry_integrity import AGENT_SCRIPT_CONTRACTS

    entry = AGENT_SCRIPT_CONTRACTS["integration_agent"]
    assert entry["required"] == (
        "aria/scripts/integration_paired_qc.py",
        "aria/scripts/integration_wnn.py",
    )
    assert set(entry["optional"]) == set(_INTEGRATION_SCRIPTS) - set(entry["required"])
    assert entry["planned"] == ()
    for script in entry["required"]:
        assert (ROOT / script).exists(), script
    # Optional scripts are existence-checked as well (stubbed for Track 1).
    _stub_new_integration_scripts(monkeypatch)
    for script in entry["optional"]:
        assert (ROOT / script).exists(), script


# ── 5. Run ledger + IntegrationNarrator ──────────────────────────────────────


def _integration_result(preprocessing: dict[str, Any] | None = None,
                        wnn: dict[str, Any] | None = None,
                        mofa: dict[str, Any] | None = None) -> dict[str, Any]:
    findings: dict[str, Any] = {}
    if preprocessing is not None:
        findings["preprocessing"] = preprocessing
    if wnn is not None:
        findings["wnn"] = wnn
    if mofa is not None:
        findings["mofa"] = mofa
    return {"status": "done", "findings": findings}


_COHORT_PRE = {
    "status": "success",
    "mode": "cohort_merge",
    "platform": "10x_multiome_arc_v1",
    "inputs": ["/virtual/donor1.h5mu", "/virtual/donor2.h5mu"],
    "output_path": "/virtual/out/cohort_merged.h5mu",
    "n_samples": 2,
    "n_cells_total": 9120,
    "n_genes": 36601,
    "n_consensus_peaks": 142233,
    "barcode_translation": "sample_prefixed",
}


def _ledger_entry(ledger: dict[str, Any], analysis: str) -> dict[str, Any]:
    return next(
        e for e in ledger["entries"]
        if e["modality"] == "integration" and e["analysis"] == analysis
    )


def test_ledger_preprocessing_executed_only_with_h5mu_artifact():
    ledger = build_run_ledger({}, {"integration_agent": _integration_result(_COHORT_PRE)})
    entry = _ledger_entry(ledger, "preprocessing")
    assert entry["node_id"] == "ledger://integration/preprocessing"
    assert entry["status"] == "executed"
    prov = entry["provenance"]
    assert prov["script"] == "aria/scripts/integration_merge_cohort.py"
    assert prov["environment_stack"] == "integration"
    assert prov["artifacts"] == ["/virtual/out/cohort_merged.h5mu"]
    assert prov["inputs"] == ["/virtual/donor1.h5mu", "/virtual/donor2.h5mu"]
    assert prov["validation_level"] == "production"


def test_ledger_preprocessing_success_without_artifact_is_skipped():
    pre = {k: v for k, v in _COHORT_PRE.items() if k != "output_path"}
    exp_ctx = {"design_intelligence": {"recommended": ["Cohort merge onto consensus peaks"]}}
    ledger = build_run_ledger(exp_ctx, {"integration_agent": _integration_result(pre)})
    entry = _ledger_entry(ledger, "preprocessing")
    assert entry["status"] == "skipped"
    assert entry["reason"] == "no assembled/merged .h5mu produced"
    assert entry["planned"] is True
    assert entry["divergence"] is True


def test_ledger_raw_assembly_uses_assembly_script():
    pre = {
        "status": "done",
        "mode": "raw_multiome_fastq",
        "rna_fastqs": ["/r/rna_R1.fq.gz"],
        "atac_fastqs": ["/r/atac_R1.fq.gz"],
        "output_path": "/virtual/out/assembled.h5mu",
    }
    ledger = build_run_ledger({}, {"integration_agent": _integration_result(pre)})
    entry = _ledger_entry(ledger, "preprocessing")
    assert entry["status"] == "executed"
    assert entry["provenance"]["script"] == "aria/scripts/integration_assemble_multiome.py"
    assert "/r/rna_R1.fq.gz" in entry["provenance"]["inputs"]


def test_ledger_absent_preprocessing_is_not_run():
    ledger = build_run_ledger({}, {"integration_agent": _integration_result(wnn={"status": "success"})})
    assert _ledger_entry(ledger, "preprocessing")["status"] == "not_run"


def test_node_id_aliases_map_intake_to_preprocessing():
    assert node_id_for("integration", "multiome_assembly") == "ledger://integration/preprocessing"
    assert node_id_for("integration", "cohort_merge") == "ledger://integration/preprocessing"
    assert node_id_for("multiome", "preprocessing") == "ledger://integration/preprocessing"


def _all_text(blocks) -> list[str]:
    out: list[str] = []
    for b in blocks:
        out.append(b.claim or "")
        out.extend(c.text for c in b.caveats)
        out.extend(str(e.label) for e in b.evidence)
    return out


def test_narrator_preprocessing_block_is_descriptive():
    narrator = IntegrationNarrator()
    result = _integration_result(_COHORT_PRE)
    assert narrator.accepts("integration_agent", {"findings": {"preprocessing": _COHORT_PRE}})
    blocks = narrator.collect("integration_agent", result, {"agent_results": {}})
    pre = next(b for b in blocks if b.id == "multiome.preprocessing")
    assert pre.analysis == "preprocessing"
    assert pre.modality == "integration"
    assert pre.status == "success"
    assert "2 sample(s)" in pre.claim
    assert "142233 consensus ATAC peaks" in pre.claim
    assert "sample_prefixed" in pre.claim
    assert pre.metrics["artifacts"] == ["/virtual/out/cohort_merged.h5mu"]
    for text in _all_text([pre]):
        assert find_causal_language(text) is None, text
    methods = " ".join(narrator.methods("integration_agent", result))
    assert "Multiome intake" in methods
    assert "142233 sweep-line consensus peaks" in methods
    assert find_causal_language(methods) is None


@pytest.mark.parametrize(
    "pre,expected_status",
    [
        ({"status": "error", "reason": "chromap index missing"}, "error"),
        ({"status": "skipped", "reason": "input already paired"}, "skipped"),
        ({k: v for k, v in _COHORT_PRE.items() if k != "output_path"}, "skipped"),
    ],
)
def test_narrator_preprocessing_non_executed_is_honest_diagnostic(pre, expected_status):
    narrator = IntegrationNarrator()
    result = _integration_result(pre)
    blocks = narrator.collect("integration_agent", result, {"agent_results": {}})
    block = next(b for b in blocks if b.id == "multiome.preprocessing")
    assert block.status == expected_status
    assert block.claim == ""
    assert not any(
        p.startswith("Multiome intake") for p in narrator.methods("integration_agent", result)
    )


_BATCH_WNN = {
    "status": "success",
    "n_cells": 9000,
    "n_joint_clusters": 14,
    "n_rna_only_clusters": 12,
    "n_atac_only_clusters": 11,
    "wnn_vs_rna_ari": 0.81,
    "wnn_vs_atac_ari": 0.74,
    "wnn_vs_rna_nmi": 0.85,
    "wnn_vs_atac_nmi": 0.79,
    "n_discordant_clusters": 1,
    "guardrail_status": "pass",
    "k": 20,
    "resolution": 0.5,
    "batch_integration": {
        "applied": True,
        "method": "Harmony",
        "batch_key": "donor",
        "n_batches": 2,
        "rna_batch_silhouette_before": 0.41,
        "rna_batch_silhouette_after": 0.07,
        "atac_batch_silhouette_before": 0.38,
        "atac_batch_silhouette_after": 0.09,
        "wnn_batch_ilisi_median": 1.82,
        "qc_assessment": {"status": "warn", "warnings": ["possible_overcorrection"]},
    },
}


def test_narrator_batch_harmonized_wnn_block_and_methods():
    narrator = IntegrationNarrator()
    result = _integration_result(wnn=copy.deepcopy(_BATCH_WNN))
    blocks = narrator.collect("integration_agent", result, {"agent_results": {}})
    wnn = next(b for b in blocks if b.id == "multiome.wnn")
    labels = {e.label: e.value for e in wnn.evidence}
    assert labels["batch key"] == "donor"
    assert labels["batches"] == 2
    assert labels["RNA batch silhouette before"] == 0.41
    assert labels["RNA batch silhouette after"] == 0.07
    assert labels["ATAC batch silhouette after"] == 0.09
    assert labels["WNN batch iLISI median"] == 1.82
    assert any("possible_overcorrection" in c.text for c in wnn.caveats)
    assert wnn.confidence == "medium"
    assert wnn.metrics["batch_integration"]["batch_key"] == "donor"
    for text in _all_text([wnn]):
        assert find_causal_language(text) is None, text

    methods = narrator.methods("integration_agent", result)
    batch_para = next(p for p in methods if p.startswith("Batch-aware WNN"))
    assert "Harmony" in batch_para
    assert "'donor' (2 batches)" in batch_para
    assert "0.41 before vs 0.07 after" in batch_para
    assert "1.82" in batch_para
    assert "<" not in batch_para and ">" not in batch_para
    assert find_causal_language(" ".join(methods)) is None


def test_narrator_wnn_without_batch_integration_is_unchanged():
    wnn = {k: v for k, v in _BATCH_WNN.items() if k != "batch_integration"}
    narrator = IntegrationNarrator()
    result = _integration_result(wnn=wnn)
    block = next(
        b for b in narrator.collect("integration_agent", result, {"agent_results": {}})
        if b.id == "multiome.wnn"
    )
    assert block.metrics["batch_integration"] is None
    assert not any(e.source == "wnn_batch" for e in block.evidence)
    assert not any(
        p.startswith("Batch-aware WNN") for p in narrator.methods("integration_agent", result)
    )
    # Not applied -> no batch surfacing either.
    wnn_na = {**wnn, "batch_integration": {"applied": False, "batch_key": "donor"}}
    block_na = next(
        b for b in narrator.collect(
            "integration_agent", _integration_result(wnn=wnn_na), {"agent_results": {}}
        )
        if b.id == "multiome.wnn"
    )
    assert block_na.metrics["batch_integration"] is None


def test_narrator_mofa_surfaces_batch_factors():
    mofa = {
        "status": "success",
        "n_cells": 9000,
        "n_factors": 15,
        "n_active_factors": 9,
        "excluded_factors": [{"factor_id": 2, "reason": "depth"}],
        "batch_factor_ids": [4, 7],
        "variance_explained": {"rna": 0.21, "atac": 0.12},
    }
    narrator = IntegrationNarrator()
    result = _integration_result(mofa=mofa)
    block = next(
        b for b in narrator.collect("integration_agent", result, {"agent_results": {}})
        if b.id == "multiome.mofa"
    )
    assert block.metrics["batch_factor_ids"] == [4, 7]
    assert any("Factor 4, 7" in c.text for c in block.caveats)
    labels = {e.label for e in block.evidence}
    assert {"batch-associated Factor 4", "batch-associated Factor 7"} <= labels
    methods = " ".join(narrator.methods("integration_agent", result))
    assert "Factor 4, 7" in methods
    for text in _all_text([block]) + [methods]:
        assert find_causal_language(text) is None, text
