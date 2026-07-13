"""E5 guards for typed, multi-library 10x scATAC FASTQ manifests."""

from __future__ import annotations

import gzip
import json
from pathlib import Path


def _manifest(tmp_path: Path) -> dict:
    libraries = []
    for index, (library, sample, donor) in enumerate((
        ("lib_a", "sample_a", "donor_a"),
        ("lib_b", "sample_b", "donor_b"),
    ), start=1):
        roles = {}
        for role in ("R1", "R2", "R3"):
            path = tmp_path / f"{library}_S{index}_L001_{role}_001.fastq.gz"
            path.write_bytes(b"fastq")
            roles[role] = str(path)
        whitelist = tmp_path / f"{library}_whitelist.txt"
        whitelist.write_text("AAAA\n")
        libraries.append({
            "library_id": library,
            "sample_id": sample,
            "donor_id": donor,
            "fastqs": roles,
            "barcode_whitelist": str(whitelist),
            "metadata": {"condition": f"condition_{index}"},
        })
    return {"schema_version": "1", "libraries": libraries}


def test_e5_typed_manifest_loads_two_libraries_and_all_roles(tmp_path):
    from aria.utils.scatac_fastq_manifest import resolve_scatac_fastq_manifest

    raw = _manifest(tmp_path)
    manifest_path = tmp_path / "scatac_fastq_manifest.json"
    manifest_path.write_text(json.dumps(raw))

    result = resolve_scatac_fastq_manifest(
        {"data_dir": str(tmp_path),
         "scatac_fastq_manifest_path": str(manifest_path)},
        require_paths=True,
    )

    assert result["status"] == "valid"
    manifest = result["manifest"]
    assert manifest["schema_version"] == "1"
    assert [row["library_id"] for row in manifest["libraries"]] == [
        "lib_a", "lib_b"
    ]
    assert all(set(row["fastqs"]) == {"R1", "R2", "R3"}
               for row in manifest["libraries"])
    assert [row["metadata"]["donor_id"] for row in manifest["libraries"]] == [
        "donor_a", "donor_b"
    ]


def test_e5_manifest_is_authoritative_for_scatac_file_classification(tmp_path):
    from aria.utils.scatac_fastq_manifest import manifest_library_types

    manifest = _manifest(tmp_path)
    declarations = manifest_library_types(manifest)

    expected_fastqs = {
        path
        for library in manifest["libraries"]
        for path in library["fastqs"].values()
    }
    assert declarations == {path: "scATAC" for path in expected_fastqs}


def test_e5_data_audit_auto_discovers_manifest_for_tui_and_headless(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent

    manifest = _manifest(tmp_path)
    (tmp_path / "scatac_fastq_manifest.json").write_text(json.dumps(manifest))

    class _Memory:
        create_wing = lambda *args, **kwargs: None
        create_hall = lambda *args, **kwargs: None

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.memory = _Memory()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_escalation = lambda *args, **kwargs: None

    result = agent.run("exp", {
        "data_dir": str(tmp_path),
        "user_question": "compare conditions",
    })

    assert result["status"] == "awaiting_checkpoint"
    exp_context = result["exp_context"]
    assert len(exp_context["modalities"]["scATAC"]) == 6
    assert "unknown" not in exp_context["modalities"]
    assert len(exp_context["scatac_fastq_manifest"]["libraries"]) == 2


def test_e5_invalid_explicit_manifest_cannot_fall_through_to_bulk_rna(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent

    manifest = _manifest(tmp_path)
    del manifest["libraries"][1]["fastqs"]["R3"]
    (tmp_path / "scatac_fastq_manifest.json").write_text(json.dumps(manifest))

    class _Memory:
        create_wing = lambda *args, **kwargs: None
        create_hall = lambda *args, **kwargs: None

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.memory = _Memory()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_escalation = lambda *args, **kwargs: None
    result = agent.run("exp", {
        "data_dir": str(tmp_path), "user_question": "compare conditions"
    })

    exp_context = result["exp_context"]
    assert "scATAC" in exp_context["modalities"]
    assert "bulk_RNA_raw" not in exp_context["modalities"]
    assert exp_context["scatac_fastq_manifest_validation"]["status"] == "invalid"


def test_e5_assay_contract_blocks_incomplete_library_manifest(tmp_path):
    from aria.utils.assay_contracts import validate_assay_contract

    manifest = _manifest(tmp_path)
    del manifest["libraries"][1]["fastqs"]["R3"]
    files = [
        path
        for library in manifest["libraries"]
        for path in library["fastqs"].values()
    ]

    result = validate_assay_contract(
        {
            "genome": "hg38",
            "scatac_fastq_manifest": manifest,
            "modalities": {"scATAC": files},
        },
        "scATAC",
    )

    barcode = result["checks"]["assay_contract"]["barcode_namespace"]
    assert result["status"] == "red"
    assert barcode["status"] == "blocked"
    assert any("R3" in error for error in barcode["errors"])


def test_e5_manifest_rejects_unsafe_library_id(tmp_path):
    from aria.utils.scatac_fastq_manifest import resolve_scatac_fastq_manifest

    manifest = _manifest(tmp_path)
    manifest["libraries"][0]["library_id"] = "../escape"
    result = resolve_scatac_fastq_manifest(
        {"data_dir": str(tmp_path), "scatac_fastq_manifest": manifest},
        require_paths=True,
    )

    assert result["status"] == "invalid"
    assert any("library_id" in error for error in result["errors"])


def test_e5_agent_dispatches_canonical_two_library_manifest(tmp_path, monkeypatch):
    from aria.agents.chromatin_agent import ChromatinAgent

    manifest = _manifest(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")
    calls = []

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            calls.append((stack, script_path, params))
            return {
                "status": "success",
                "fragments_file": str(tmp_path / "combined.fragments.tsv.gz"),
                "library_manifest": params["libraries"],
                "n_libraries": 2,
            }

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = _Env()
    agent.publish_status = lambda *args, **kwargs: None

    result = agent._run_scatac_fastq_to_fragments(
        "exp",
        {
            "genome": "hg38",
            "genome_fasta": str(genome),
            "scatac_fastq_manifest": manifest,
        },
        [path for row in manifest["libraries"]
         for path in row["fastqs"].values()],
    )

    assert result["status"] == "success"
    assert len(calls) == 1
    stack, script, params = calls[0]
    assert stack == "atacseq"
    assert script.endswith("chromatin_scatac_align.py")
    assert [row["library_id"] for row in params["libraries"]] == [
        "lib_a", "lib_b"
    ]
    assert params["libraries"][0]["metadata"]["donor_id"] == "donor_a"
    assert "r1_fastq" not in params and "barcode_fastq" not in params


def test_e5_fragments_bridge_omits_manifest_for_legacy_input(tmp_path):
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.utils.script_contracts import contract_for_script

    fragments = tmp_path / "fragments.tsv.gz"
    fragments.write_bytes(b"fragments")
    calls = []

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            calls.append((script_path, params))
            return {"status": "success", "output_path": "matrix.h5ad"}

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = _Env()
    result = agent._run_fragments_to_matrix(
        "exp", {"genome": "hg38"}, [str(fragments)], {}
    )

    assert result["status"] == "success"
    script_path, params = calls[0]
    assert "library_manifest" not in params
    contract = contract_for_script(script_path)
    assert contract is not None
    assert contract.validate_params(params) == []


def test_e5_fragment_union_prefixes_barcodes_and_writes_manifest(tmp_path):
    from aria.scripts.chromatin_scatac_align import _merge_library_fragments

    manifest = _manifest(tmp_path)
    records = []
    for row in manifest["libraries"]:
        fragments = tmp_path / f"{row['library_id']}.fragments.tsv.gz"
        with gzip.open(fragments, "wt") as handle:
            handle.write("chr1\t10\t20\tAAAC-1\t1\n")
        records.append({**row, "fragments_file": str(fragments)})

    output = tmp_path / "combined.fragments.tsv"
    result = _merge_library_fragments(records, output)

    assert result["status"] == "success"
    lines = output.read_text().splitlines()
    assert lines == [
        "chr1\t10\t20\tlib_a#AAAC-1\t1",
        "chr1\t10\t20\tlib_b#AAAC-1\t1",
    ]
    assert result["barcode_prefixes"] == {
        "lib_a": "lib_a#", "lib_b": "lib_b#"
    }


def test_e5_multi_library_aligner_runs_each_library_then_unions(
    tmp_path, monkeypatch
):
    from aria.scripts import chromatin_scatac_align as align

    manifest = _manifest(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")
    calls = []

    def fake_single(params):
        library_id = Path(params["output_dir"]).name
        calls.append(library_id)
        fragments = tmp_path / f"{library_id}.fragments.tsv.gz"
        with gzip.open(fragments, "wt") as handle:
            handle.write("chr1\t10\t20\tAAAC-1\t1\n")
        return {"status": "success", "fragments_file": str(fragments)}

    monkeypatch.setattr(align, "chromatin_scatac_align", fake_single)
    monkeypatch.setattr(align, "_bgzip_tabix", lambda *args: "missing_tool")

    result = align._align_libraries({
        "libraries": manifest["libraries"],
        "genome_fasta": str(genome),
        "output_dir": str(tmp_path / "aligned"),
    })

    assert result["status"] == "success"
    assert calls == ["lib_a", "lib_b"]
    assert result["n_libraries"] == 2
    assert [row["sample_id"] for row in result["library_manifest"]] == [
        "sample_a", "sample_b"
    ]
    assert Path(result["library_manifest_path"]).is_file()


def test_e5_multi_library_aligner_is_all_or_fail(tmp_path, monkeypatch):
    from aria.scripts import chromatin_scatac_align as align

    manifest = _manifest(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")

    def fake_single(params):
        library_id = Path(params["output_dir"]).name
        return ({"status": "skipped", "reason": "alignment_failed"}
                if library_id == "lib_b" else
                {"status": "success", "fragments_file": "/tmp/lib_a.tsv"})

    monkeypatch.setattr(align, "chromatin_scatac_align", fake_single)
    result = align._align_libraries({
        "libraries": manifest["libraries"],
        "genome_fasta": str(genome),
        "output_dir": str(tmp_path / "aligned"),
    })

    assert result["status"] == "skipped"
    assert result["reason"] == "partial_library_alignment_failure"
    assert result["failed_libraries"] == ["lib_b"]


def test_e5_peak_matrix_annotation_recovers_library_sample_and_donor(tmp_path):
    import pandas as pd

    from aria.scripts.chromatin_fragments_to_matrix import (
        _annotate_library_metadata,
    )

    class _Adata:
        obs = pd.DataFrame(index=["lib_a#AAAC-1", "lib_b#AAAC-1"])

    manifest = _manifest(tmp_path)
    annotated = _annotate_library_metadata(_Adata(), manifest["libraries"])

    assert annotated == 2
    assert _Adata.obs["library_id"].tolist() == ["lib_a", "lib_b"]
    assert _Adata.obs["sample_id"].tolist() == ["sample_a", "sample_b"]
    assert _Adata.obs["donor_id"].tolist() == ["donor_a", "donor_b"]
    assert _Adata.obs["condition"].tolist() == ["condition_1", "condition_2"]


def test_e5_headless_context_and_cp1_corrections_accept_manifest(tmp_path):
    from aria.agents.data_audit_agent import apply_metadata_corrections
    from aria.headless import _build_headless_context

    manifest = _manifest(tmp_path)
    context = _build_headless_context(
        str(tmp_path), "question", False, False,
        context_overrides={"scatac_fastq_manifest": manifest},
    )
    assert context["scatac_fastq_manifest"]["libraries"][1]["donor_id"] == \
        "donor_b"

    exp_context = {"modalities": {"scATAC": []}}
    apply_metadata_corrections(
        exp_context, {"scatac_fastq_manifest": manifest}
    )
    corrected = exp_context["scatac_fastq_manifest"]
    assert corrected["libraries"][0]["metadata"]["donor_id"] == "donor_a"
    assert len(exp_context["modalities"]["scATAC"]) == 6


def test_e5_cockpit_editor_parses_and_validates_manifest(tmp_path):
    from aria.ui.scatac_manifest_editor import parse_scatac_manifest_edit

    manifest = _manifest(tmp_path)
    parsed, errors = parse_scatac_manifest_edit(
        json.dumps(manifest), str(tmp_path)
    )

    assert errors == []
    assert parsed is not None
    assert len(parsed["libraries"]) == 2
