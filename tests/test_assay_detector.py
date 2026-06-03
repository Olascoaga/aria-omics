from __future__ import annotations

import sys
import types

import pytest

litellm_stub = types.ModuleType("litellm")
litellm_stub.completion = lambda *args, **kwargs: None
sys.modules.setdefault("litellm", litellm_stub)

from aria.agents.data_audit_agent import DataAuditAgent


def _agent() -> DataAuditAgent:
    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.publish_status = lambda *args, **kwargs: None
    return agent


def test_data_audit_classifies_h5ad_by_content_with_nonstandard_suffix(tmp_path):
    h5py = pytest.importorskip("h5py")

    p = tmp_path / "processed_rna_payload.data"
    with h5py.File(p, "w") as h5:
        h5.create_group("obs")
        h5.create_group("var")
        h5.create_dataset("X", data=[[1, 0], [0, 1]])

    agent = _agent()
    assert p in agent._scan_directory(tmp_path)

    classified = agent._classify_files([p])

    assert classified == {"scRNA": [str(p)]}


def test_data_audit_classifies_paired_h5mu_by_internal_modalities(tmp_path):
    h5py = pytest.importorskip("h5py")

    p = tmp_path / "paired_multiome.payload"
    with h5py.File(p, "w") as h5:
        mod = h5.create_group("mod")
        mod.create_group("rna")
        mod.create_group("atac")

    agent = _agent()
    assert p in agent._scan_directory(tmp_path)

    classified = agent._classify_files([p])

    assert classified == {"scATAC": [str(p)]}


def test_data_audit_distinguishes_10x_atac_h5_by_feature_type(tmp_path):
    h5py = pytest.importorskip("h5py")

    p = tmp_path / "cellranger_payload.blob"
    with h5py.File(p, "w") as h5:
        matrix = h5.create_group("matrix")
        matrix.create_dataset("barcodes", data=[b"cell1", b"cell2"])
        matrix.create_dataset("data", data=[1, 2])
        features = matrix.create_group("features")
        features.create_dataset("feature_type", data=[b"Peaks", b"Gene Expression"])
        features.create_dataset("genome", data=[b"hg38", b"hg38"])

    agent = _agent()
    classified = agent._classify_files([p])

    assert classified == {"scATAC": [str(p)]}
    record = agent._last_assay_detections[0]
    assert record["confidence"] == "high"
    assert record["evidence"]["format"] == "10x_multiome_h5"
    assert record["possible_alternatives"] == ["scRNA"]


def test_data_audit_classifies_count_table_by_numeric_content(tmp_path):
    p = tmp_path / "matrix.payload"
    p.write_text(
        "gene_id,sample_a,sample_b,sample_c\n"
        "Gene1,10,12,9\n"
        "Gene2,0,3,1\n"
        "Gene3,15,14,16\n"
    )

    agent = _agent()
    assert p in agent._scan_directory(tmp_path)

    classified = agent._classify_files([p])

    assert classified == {"bulk_RNA": [str(p)]}


def test_data_audit_classifies_salmon_quant_with_opaque_name(tmp_path):
    p = tmp_path / "transcript_payload.out"
    p.write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "ENST1\t1000\t900\t12.5\t55.2\n"
        "ENST2\t2000\t1800\t0.0\t0.0\n"
    )

    agent = _agent()
    assert p in agent._scan_directory(tmp_path)

    classified = agent._classify_files([p])

    assert classified == {"bulk_RNA": [str(p)]}
    assert agent._last_assay_detections[0]["evidence"]["format"] == "salmon_quant"


def test_data_audit_classifies_misnamed_bam_from_header(tmp_path):
    p = tmp_path / "aligned_payload.bin"
    header = b"@HD\tVN:1.6\n@PG\tID:STAR\tPN:STAR\n"
    p.write_bytes(b"BAM\x01" + len(header).to_bytes(4, "little") + header)

    agent = _agent()
    assert p in agent._scan_directory(tmp_path)

    classified = agent._classify_files([p])

    assert classified == {"bulk_RNA": [str(p)]}
    record = agent._last_assay_detections[0]
    assert record["confidence"] == "medium"
    assert record["evidence"]["format"] == "bam_or_sam"


def test_data_audit_keeps_signature_fallback_when_content_is_uninformative(tmp_path):
    p = tmp_path / "sample_counts.tsv"
    p.write_text("placeholder\n")

    agent = _agent()
    classified = agent._classify_files([p])

    assert classified == {"bulk_RNA": [str(p)]}
