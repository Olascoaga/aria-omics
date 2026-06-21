from __future__ import annotations

from aria.agents.modality_audit import build_capability_matrix
from aria.utils.multiome_contracts import (
    infer_multiome_contract,
    validate_multiome_contract,
)


def test_paired_h5mu_detection_infers_ready_multiome_contract():
    exp_context = {
        "genome": "hg38",
        "modalities": {"scATAC": ["/data/pbmc_multiome.h5mu"]},
        "assay_detections": [
            {
                "path": "/data/pbmc_multiome.h5mu",
                "modality": "scATAC",
                "confidence": "high",
                "evidence": {
                    "format": "h5mu",
                    "modalities": ["rna", "atac"],
                    "paired_rna": True,
                },
            }
        ],
    }

    contract = infer_multiome_contract(exp_context)
    assert contract["object_type"] == "paired_mudata"
    assert contract["same_cell"] is True
    assert contract["cell_namespace"] == "shared_mudata_obs"

    validation = validate_multiome_contract(exp_context)
    assert validation["status"] == "green"
    assert validation["findings"] == []
    assert validation["checks"]["cell_namespace"]["status"] == "pass"
    assert validation["checks"]["genome_feature_space"]["feature_space"] == "peak_matrix"


def test_split_scrna_scatac_without_pairing_degrades_not_blocks_modalities():
    exp_context = {
        "genome": "hg38",
        "modalities": {
            "scRNA": ["/data/rna.h5ad"],
            "scATAC": ["/data/fragments.tsv.gz"],
        },
        "design": {
            "groups": {"ctrl": ["c1", "c2"], "stim": ["s1", "s2"]},
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col": "cell_type",
            },
        },
    }

    validation = validate_multiome_contract(exp_context)
    assert validation["status"] == "yellow"
    assert validation["contract"]["degradation"] == "pairing_unconfirmed"
    assert any(
        finding["check"] == "multiome_contract_pairing_unconfirmed"
        for finding in validation["findings"]
    )

    matrix = build_capability_matrix(
        exp_context,
        modality_validation={
            "scRNA": {"level": "production", "dispatch_enabled": True},
            "scATAC": {"level": "beta", "dispatch_enabled": True},
        },
    )
    assert matrix["contracts"]["multiome"]["status"] == "yellow"
    assert "multiome" not in matrix["dispatch"]["blocked"]
    assert matrix["dispatch"]["blocked"] == []
    assert "scATAC" in matrix["dispatch"]["requires_ack"]
    assert {
        finding["check"] for finding in matrix["findings"]
    } >= {"multiome_contract_pairing_unconfirmed"}


def test_explicit_same_cell_contract_requires_barcode_namespace():
    validation = validate_multiome_contract(
        {
            "genome": "hg38",
            "modalities": {
                "scRNA": ["/data/rna.h5ad"],
                "scATAC": ["/data/atac.h5ad"],
            },
            "multiome_contract": {
                "object_type": "explicit_same_cell",
                "same_cell": True,
                "rna_modality": "scRNA",
                "atac_modality": "scATAC",
                "atac_feature_space": "peak_matrix",
            },
        }
    )

    assert validation["status"] == "red"
    assert validation["checks"]["cell_namespace"]["status"] == "blocked"
    assert any(
        finding["severity"] == "blocking"
        and finding["check"] == "multiome_contract_barcode_namespace_missing"
        for finding in validation["findings"]
    )


def test_explicit_same_cell_contract_accepts_shared_namespace_and_samples():
    validation = validate_multiome_contract(
        {
            "genome": "hg38",
            "modalities": {
                "scRNA": ["/data/rna.h5ad"],
                "scATAC": ["/data/atac.h5ad"],
            },
            "multiome_contract": {
                "object_type": "explicit_same_cell",
                "same_cell": True,
                "cell_namespace": "10x_multiome_barcodes",
                "rna_sample_ids": ["donor1", "donor2"],
                "atac_sample_ids": ["donor2", "donor1"],
                "atac_feature_space": "peak_matrix",
            },
        }
    )

    assert validation["status"] == "green"
    assert validation["checks"]["cell_namespace"]["status"] == "pass"
    assert validation["checks"]["sample_alignment"]["status"] == "pass"
    assert validation["findings"] == []


def test_explicit_same_cell_contract_blocks_mismatched_samples():
    validation = validate_multiome_contract(
        {
            "genome": "hg38",
            "modalities": {
                "scRNA": ["/data/rna.h5ad"],
                "scATAC": ["/data/atac.h5ad"],
            },
            "multiome_contract": {
                "object_type": "explicit_same_cell",
                "same_cell": True,
                "cell_namespace": "10x_multiome_barcodes",
                "rna_sample_ids": ["donor1", "donor2"],
                "atac_sample_ids": ["donor1", "donor3"],
                "atac_feature_space": "peak_matrix",
            },
        }
    )

    assert validation["status"] == "red"
    assert validation["checks"]["sample_alignment"]["missing_in_atac"] == ["donor2"]
    assert validation["checks"]["sample_alignment"]["missing_in_rna"] == ["donor3"]
    assert any(
        finding["check"] == "multiome_contract_sample_mismatch"
        for finding in validation["findings"]
    )


def test_data_audit_context_includes_inferred_multiome_contract(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent._last_assay_detections = [
        {
            "path": str(tmp_path / "paired.h5mu"),
            "modality": "scATAC",
            "confidence": "high",
            "evidence": {"format": "h5mu", "paired_rna": True},
        }
    ]
    p = tmp_path / "paired.h5mu"
    p.write_bytes(b"")

    exp_context = agent._build_context(
        "exp1",
        tmp_path,
        {"scATAC": [str(p)]},
        "hg38",
        "Homo sapiens",
        [],
        "multiome",
    )

    assert exp_context["multiome_contract"]["object_type"] == "paired_mudata"
    assert exp_context["multiome_contract"]["same_cell"] is True
