"""W-LEDGER export (residual closure): a machine-readable reproducible capsule
(RO-Crate / W3C-PROV-style JSON-LD) per report, plus `aria diff` over two run
ledgers. Pure bookkeeping over methodology.json; no LLM, no biology.
"""

import json
import zipfile

from aria.agents.narrative.ledger_export import (
    build_ro_crate,
    write_ro_crate,
    write_reproducible_capsule,
    diff_methodologies,
    format_diff,
)


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

    crate_path = write_ro_crate(rd)
    assert crate_path.exists() and crate_path.name == "ro-crate-metadata.json"

    zip_path = write_reproducible_capsule(rd)
    assert zip_path.exists()
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
    assert any(n.endswith("methodology.json") for n in names)
    assert any(n.endswith("ro-crate-metadata.json") for n in names)
    assert any(n.endswith("report.html") for n in names)


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
