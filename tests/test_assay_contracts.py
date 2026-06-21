from __future__ import annotations

from aria.agents.modality_audit import build_capability_matrix
from aria.utils.assay_contracts import validate_assay_contract


def test_chromatin_contract_blocks_missing_genome_before_dispatch():
    matrix = build_capability_matrix(
        {"modalities": {"bulk_ATAC": ["/data/sample.bam"]}},
        modality_validation={"bulk_ATAC": {"level": "beta", "dispatch_enabled": True}},
    )

    card = matrix["cards"]["bulk_ATAC"]
    assert card["status"] == "red"
    assert card["dispatch_policy"] == "blocked"
    assert matrix["dispatch"]["blocked"] == ["bulk_ATAC"]
    assert "bulk_ATAC" not in matrix["dispatch"]["requires_ack"]
    assert card["checks"]["assay_contract"]["genome"]["required"] is True
    assert card["checks"]["assay_contract"]["genome"]["valid"] is False
    assert any(
        finding["check"] == "assay_contract_genome_missing"
        for finding in card["findings"]
    )


def test_scatac_fastq_contract_requires_barcode_namespace():
    matrix = build_capability_matrix(
        {"genome": "hg38", "modalities": {"scATAC": ["/data/R1.fastq.gz"]}},
        modality_validation={"scATAC": {"level": "beta", "dispatch_enabled": True}},
    )

    card = matrix["cards"]["scATAC"]
    assert card["status"] == "red"
    assert card["dispatch_policy"] == "blocked"
    assert card["checks"]["assay_contract"]["barcode_namespace"]["status"] == "blocked"
    assert any(
        finding["check"] == "assay_contract_scatac_fastq_barcode_missing"
        for finding in card["findings"]
    )


def test_scatac_fastq_contract_accepts_explicit_barcode_namespace():
    matrix = build_capability_matrix(
        {
            "genome": "hg38",
            "barcode_fastq": "/data/I1.fastq.gz",
            "barcode_whitelist": "/refs/737K-arc-v1.txt",
            "modalities": {"scATAC": ["/data/R1.fastq.gz"]},
        },
        modality_validation={"scATAC": {"level": "beta", "dispatch_enabled": True}},
    )

    card = matrix["cards"]["scATAC"]
    assert card["status"] == "yellow"
    assert card["dispatch_policy"] == "requires_ack"
    assert card["checks"]["assay_contract"]["barcode_namespace"]["status"] == "pass"
    assert not any(
        finding["check"] == "assay_contract_scatac_fastq_barcode_missing"
        for finding in card["findings"]
    )


def test_comparison_contract_warns_for_missing_sample_metadata():
    result = validate_assay_contract(
        {
            "genome": "hg38",
            "comparisons": [{"case": "stim", "control": "ctrl"}],
            "modalities": {"bulk_ATAC": ["/data/ctrl.bam", "/data/stim.bam"]},
        },
        "bulk_ATAC",
    )

    assert result["status"] == "yellow"
    contract = result["checks"]["assay_contract"]
    assert contract["sample_metadata"]["status"] == "warning"
    assert contract["replicate_structure"]["status"] == "warning"
    assert {
        "assay_contract_sample_metadata_incomplete",
        "assay_contract_replicates_low_or_missing",
    } <= {finding["check"] for finding in result["findings"]}


def test_coordinate_contract_blocks_invalid_bed_like_file(tmp_path):
    bad_bed = tmp_path / "bad.bed"
    bad_bed.write_text("chr1\t100\t50\n", encoding="utf-8")

    result = validate_assay_contract(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(bad_bed)]},
        },
        "bulk_ATAC",
    )

    assert result["status"] == "red"
    contract = result["checks"]["assay_contract"]
    assert contract["feature_coordinates"]["status"] == "blocked"
    assert contract["feature_coordinates"]["inspected"][0]["reason"] == "invalid_interval"
    assert any(
        finding["check"] == "assay_contract_coordinate_features_invalid"
        for finding in result["findings"]
    )


def test_coordinate_contract_records_valid_bed_convention(tmp_path):
    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text("chr1\t100\t250\tpeak1\n", encoding="utf-8")

    result = validate_assay_contract(
        {
            "genome": "hg38",
            "modalities": {"bulk_ATAC": [str(peaks)]},
        },
        "bulk_ATAC",
    )

    contract = result["checks"]["assay_contract"]
    assert result["status"] == "green"
    assert contract["feature_coordinates"]["status"] == "pass"
    assert contract["feature_coordinates"]["contig_style"] == "ucsc_chr"
    assert (
        contract["feature_coordinates"]["inspected"][0]["coordinate_convention"]
        == "bed_0_based_half_open"
    )
