"""W-LEDGER export (residual closure): a machine-readable reproducible capsule
(RO-Crate / W3C-PROV-style JSON-LD) per report, plus `aria diff` over two run
ledgers. Pure bookkeeping over methodology.json; no LLM, no biology.
"""

import json
import zipfile

import pytest

from aria.agents.narrative.ledger_export import (
    build_ro_crate,
    build_capsule_manifest,
    write_ro_crate,
    write_reproducible_capsule,
    verify_reproducible_capsule,
    diff_methodologies,
    format_diff,
)
from aria.agents.narrative.hypothesis import SpeculativePromotionError


def _methodology(version="4.5.4", commit="abc123", claims=None, ran=True,
                 calibration_status="pass"):
    return {
        "provenance": {
            "version": version, "git_commit": commit, "git_dirty": False,
            "workflow_hash": "wf" + commit,
            "timestamp_utc": "2026-06-04T00:00:00Z",
        },
        "inputs": [{"path": "/data/x.h5ad", "sha256": "deadbeef", "bytes": 100}],
        "seeds": {"global": 0},
        "claims": claims if claims is not None else [{
            "claim_id": "scrna.pb.A", "text": "A is DE", "modality": "scRNA-seq",
            "analysis": "pseudobulk_de", "tier": "associative",
            "evidence_card_id": "scrna.pb.A#evidence",
            "ledger_node_id": "ledger://scRNA/pseudobulk_de", "ledger_status": "ran",
        }],
        "run_ledger": {
            "entries": [{
                "node_id": "ledger://scRNA/pseudobulk_de",
                "analysis": "pseudobulk_de", "modality": "scRNA",
                "status": "ran" if ran else "skipped", "planned": True,
            }],
            "n_divergences": 0 if ran else 1,
        },
        "calibration": {"measured": True, "status": calibration_status,
                        "summary": {"bulk_recall": 0.9}},
    }


# ── RO-Crate / capsule ───────────────────────────────────────────────────────

def test_ro_crate_is_valid_jsonld_with_provenance_and_claims():
    crate = build_ro_crate(_methodology())
    assert crate["@context"]
    graph = crate["@graph"]
    ids = {e["@id"] for e in graph}
    assert "ro-crate-metadata.json" in ids
    assert "./" in ids                      # root dataset
    blob = json.dumps(crate)
    assert "4.5.4" in blob and "abc123" in blob          # provenance present
    assert "ledger://scRNA/pseudobulk_de" in blob        # claim -> ledger node
    assert "scrna.pb.A#evidence" in blob                 # claim -> evidence card


def test_write_ro_crate_and_capsule(tmp_path):
    rd = tmp_path / "report"
    rd.mkdir()
    (rd / "methodology.json").write_text(json.dumps(_methodology()))
    (rd / "report.html").write_text("<html></html>")
    repo = tmp_path / "repo"
    (repo / "envs").mkdir(parents=True)
    (repo / "envs" / "aria-rna-env.linux-64.lock").write_text("lock\n")
    graph_dir = repo / "docs" / "architecture" / "graphify"
    graph_dir.mkdir(parents=True)
    (repo / "docs" / "architecture").mkdir(parents=True, exist_ok=True)
    (repo / "docs" / "architecture" / "code_graph.md").write_text("graph\n")
    (graph_dir / "graph.json").write_text("{}\n")
    (graph_dir / "manifest.json").write_text("{}\n")
    (graph_dir / "README.md").write_text("readme\n")
    (graph_dir / "GRAPH_REPORT.md").write_text("report\n")

    crate_path = write_ro_crate(rd)
    assert crate_path.exists() and crate_path.name == "ro-crate-metadata.json"

    zip_path = write_reproducible_capsule(rd, repo_root=repo)
    assert zip_path.exists()
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        manifest = json.loads(z.read("capsule_manifest.json"))
    assert any(n.endswith("methodology.json") for n in names)
    assert any(n.endswith("ro-crate-metadata.json") for n in names)
    assert any(n.endswith("report.html") for n in names)
    assert "repository/envs/aria-rna-env.linux-64.lock" in names
    assert "repository/docs/architecture/graphify/graph.json" in names
    assert manifest["schema_version"].startswith("s14.")
    assert any(f["role"] == "lockfile" for f in manifest["files"])
    assert any(f["role"] == "graph_snapshot" for f in manifest["files"])


def test_build_capsule_manifest_records_reproduction_plan(tmp_path):
    rd = tmp_path / "report"
    rd.mkdir()
    (rd / "methodology.json").write_text(json.dumps(_methodology()))
    (rd / "report.html").write_text("<html></html>")

    manifest, sources = build_capsule_manifest(rd, repo_root=tmp_path / "empty_repo")

    assert manifest["reproduction"]["automatic_rerun"] is False
    assert manifest["reproduction"]["required_git_commit"] == "abc123"
    assert manifest["reproduction"]["input_hashes"][0]["sha256"] == "deadbeef"
    assert any(src.name == "methodology.json" for src, _arc, _role in sources)


def test_verify_reproducible_capsule_passes_and_warns_for_missing_inputs(tmp_path):
    rd = tmp_path / "report"
    rd.mkdir()
    (rd / "methodology.json").write_text(json.dumps(_methodology()))
    (rd / "report.html").write_text("<html></html>")
    zip_path = write_reproducible_capsule(rd, repo_root=tmp_path / "empty_repo")

    result = verify_reproducible_capsule(zip_path, repo_root=tmp_path / "empty_repo")

    assert result["status"] == "warning"  # portable input needs local relocation
    assert result["files"]["checked"] >= 3
    assert result["files"]["mismatched"] == []
    assert result["inputs"]["missing"] == []
    assert result["inputs"]["relocation_required"] == [
        "input://sha256/deadbeef"
    ]


def test_verify_reproducible_capsule_detects_file_hash_mismatch(tmp_path):
    rd = tmp_path / "report"
    rd.mkdir()
    (rd / "methodology.json").write_text(json.dumps(_methodology()))
    (rd / "report.html").write_text("<html></html>")
    zip_path = write_reproducible_capsule(rd, repo_root=tmp_path / "empty_repo")

    rewritten = tmp_path / "corrupted.zip"
    with zipfile.ZipFile(zip_path, "r") as src, zipfile.ZipFile(rewritten, "w") as dst:
        for name in src.namelist():
            payload = src.read(name)
            if name == "report/report.html":
                payload = b"<html>corrupted</html>"
            dst.writestr(name, payload)

    result = verify_reproducible_capsule(rewritten, repo_root=tmp_path / "empty_repo")

    assert result["status"] == "fail"
    assert any(m["path"] == "report/report.html" for m in result["files"]["mismatched"])


def test_verify_reproducible_capsule_can_diff_reproduced_report(tmp_path):
    rd = tmp_path / "report"
    rd.mkdir()
    (rd / "methodology.json").write_text(json.dumps(_methodology(commit="aaa")))
    (rd / "report.html").write_text("<html></html>")
    zip_path = write_reproducible_capsule(rd, repo_root=tmp_path / "empty_repo")

    reproduced = tmp_path / "reproduced"
    reproduced.mkdir()
    (reproduced / "methodology.json").write_text(json.dumps(_methodology(commit="bbb")))

    result = verify_reproducible_capsule(
        zip_path, compare_report=reproduced, repo_root=tmp_path / "empty_repo"
    )

    assert result["diff"]["identical"] is False
    assert result["diff"]["provenance"]["git_commit"] == {"a": "aaa", "b": "bbb"}


# ── aria diff ────────────────────────────────────────────────────────────────

def test_diff_detects_ledger_and_claim_and_provenance_changes():
    a = _methodology(version="4.5.4", commit="aaa", ran=True)
    b = _methodology(version="4.6.0", commit="bbb", ran=False, claims=[])
    d = diff_methodologies(a, b)
    assert d["identical"] is False
    assert d["provenance"]["version"] == {"a": "4.5.4", "b": "4.6.0"}
    assert d["provenance"]["git_commit"]["a"] == "aaa"
    assert any(c["node_id"] == "ledger://scRNA/pseudobulk_de"
               for c in d["ledger"]["status_changed"])
    assert "scrna.pb.A" in d["claims"]["removed"]
    text = format_diff(d)
    assert "4.5.4" in text and "4.6.0" in text


def test_diff_identical_is_empty():
    d = diff_methodologies(_methodology(), _methodology())
    assert d["identical"] is True
    assert format_diff(d).strip()           # still renders a human-readable line


# ── H17: the quarantine wall is enforced by the export/diff paths ────────────

def _contaminated_methodology():
    # A quarantined hypothesis node hidden inside a claim's nested evidence.
    return _methodology(claims=[{
        "claim_id": "leak", "tier": "associative",
        "ledger_node_id": "ledger://scRNA/pseudobulk_de", "ledger_status": "ran",
        "evidence": [{"source": "hypothesis://h1"}],
    }])


def test_build_ro_crate_rejects_quarantined_hypothesis_in_claims():
    with pytest.raises(SpeculativePromotionError):
        build_ro_crate(_contaminated_methodology())


def test_write_capsule_rejects_quarantined_hypothesis(tmp_path):
    rd = tmp_path / "report"
    rd.mkdir()
    (rd / "report.html").write_text("<html></html>")
    (rd / "methodology.json").write_text(json.dumps(_contaminated_methodology()))
    with pytest.raises(SpeculativePromotionError):
        write_reproducible_capsule(rd)


def test_diff_rejects_quarantined_hypothesis_in_either_report():
    clean, dirty = _methodology(), _contaminated_methodology()
    with pytest.raises(SpeculativePromotionError):
        diff_methodologies(clean, dirty)
    with pytest.raises(SpeculativePromotionError):
        diff_methodologies(dirty, clean)
