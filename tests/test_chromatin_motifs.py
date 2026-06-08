"""v4.6 scATAC step 4 — chromatin motif enrichment guards.

Covers `aria/utils/motifs.py`, `scripts/fetch_motifs.py`, and
`aria/scripts/chromatin_motifs.py`:

- IPC contract registration;
- versioned local motif-collection resolution + manifest (W-PRIV mirror of ORA);
- honest skips when the motif collection or genome FASTA is absent (never a
  fabricated enrichment, ADR-002);
- the motif bootstrap refuses network egress under air-gapped mode (W-PRIV);
- the snapatac2 enrichment wiring (snapatac2-gated).

Synthetic labels only (ADR-011): no biological claims.
"""

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

_MEME = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0001.1 TF_ALPHA
letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0
 0.97 0.01 0.01 0.01
 0.01 0.97 0.01 0.01
 0.01 0.01 0.97 0.01
 0.01 0.01 0.01 0.97

MOTIF MA0002.1 TF_BETA
letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0
 0.25 0.25 0.25 0.25
 0.25 0.25 0.25 0.25
 0.25 0.25 0.25 0.25
 0.25 0.25 0.25 0.25
"""


def _stage_collection(motif_dir: Path, collection: str, manifest=True):
    base = motif_dir / collection
    base.mkdir(parents=True, exist_ok=True)
    (base / f"{collection}.meme").write_text(_MEME, encoding="utf-8")
    if manifest:
        import json
        (base / "manifest.json").write_text(json.dumps({
            "collection": collection, "source": "test",
            "release": "2024", "sha256": "deadbeef",
        }), encoding="utf-8")


# ── Contract ─────────────────────────────────────────────────────────────────

def test_motifs_has_ipc_contract():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS
    assert "aria/scripts/chromatin_motifs.py" in SCRIPT_CONTRACTS


# ── motifs.py resolver + manifest ─────────────────────────────────────────────

def test_motifs_dir_respects_env(tmp_path, monkeypatch):
    from aria.utils import motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path / "m"))
    assert motifs.motifs_dir() == tmp_path / "m"


def test_load_local_collection_absent_returns_none(tmp_path, monkeypatch):
    from aria.utils import motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    assert motifs.load_local_motif_collection("does_not_exist") is None


def test_load_local_collection_reads_meme_and_manifest(tmp_path, monkeypatch):
    from aria.utils import motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    _stage_collection(tmp_path, "JASPAR_TEST")
    loaded = motifs.load_local_motif_collection("JASPAR_TEST")
    assert loaded is not None
    meme_path, version = loaded
    assert meme_path.endswith("JASPAR_TEST.meme")
    assert version["n_motifs"] == 2           # counted from MOTIF records
    assert version["release"] == "2024"
    assert version["source"] == "test"


# ── chromatin_motifs honest skips (no snapatac2 needed) ───────────────────────

def test_motifs_skip_without_collection(tmp_path, monkeypatch):
    from aria.scripts.chromatin_motifs import chromatin_motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))   # empty
    res = chromatin_motifs({"genome_fasta": "x", "regions": {"0": ["chr1:1-9"]}})
    assert res["status"] == "success" and res["ran"] is False
    assert "fetch_motifs" in res["reason"]


def test_motifs_skip_without_genome(tmp_path, monkeypatch):
    from aria.scripts.chromatin_motifs import chromatin_motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    _stage_collection(tmp_path, "JASPAR_TEST")
    res = chromatin_motifs({
        "motif_collection": "JASPAR_TEST",
        "regions": {"0": ["chr1:1-9"]},
    })
    assert res["status"] == "success" and res["ran"] is False
    assert "genome_fasta" in res["reason"]
    # the motif provenance is still surfaced on the skip
    assert res["motif_source"]["collection"] == "JASPAR_TEST"


def test_motifs_skip_with_missing_genome_file(tmp_path, monkeypatch):
    from aria.scripts.chromatin_motifs import chromatin_motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    _stage_collection(tmp_path, "JASPAR_TEST")
    res = chromatin_motifs({
        "motif_collection": "JASPAR_TEST",
        "genome_fasta": str(tmp_path / "nope.fa"),
        "regions": {"0": ["chr1:1-9"]},
    })
    assert res["ran"] is False and "not found" in res["reason"]


# ── fetch_motifs governance (W-PRIV) ──────────────────────────────────────────

def _load_fetch_motifs():
    spec = importlib.util.spec_from_file_location(
        "aria_fetch_motifs", REPO_ROOT / "scripts" / "fetch_motifs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_fetch_motifs_refuses_egress_when_air_gapped(tmp_path, monkeypatch):
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    mod = _load_fetch_motifs()
    rc = mod.main(["--taxa", "vertebrates"])
    assert rc == 2                                    # refused, no download
    # nothing was written
    assert not any(tmp_path.rglob("*.meme"))


# ── snapatac2 enrichment wiring (gated) ───────────────────────────────────────

def test_motifs_no_regions_skips(tmp_path, monkeypatch):
    pytest.importorskip("snapatac2")
    from aria.scripts.chromatin_motifs import chromatin_motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    _stage_collection(tmp_path, "JASPAR_TEST")
    fa = tmp_path / "tiny.fa"
    fa.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    res = chromatin_motifs({
        "motif_collection": "JASPAR_TEST",
        "genome_fasta": str(fa),
        # no regions, no da_csv
    })
    assert res["status"] == "success" and res["ran"] is False
    assert "regions" in res["reason"] or "DA region" in res["reason"]
