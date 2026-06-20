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


# ── F7: motif threshold policy resolves by the REAL modality, not hardcoded ──

def test_motif_modality_thresholds_resolve_by_caller_not_hardcoded_scatac():
    from aria.scripts.chromatin_motifs import _resolve_modality_thresholds
    # scATAC default preserved when no modality is passed (existing behaviour)
    assert _resolve_modality_thresholds({}) == ("scATAC", 0.05)
    # bulk ATAC reuses the engine -> its own modality label, never scATAC
    modality, padj = _resolve_modality_thresholds({"modality": "bulk_ATAC"})
    assert modality == "bulk_ATAC"
    # data_type is accepted as the modality source too
    assert _resolve_modality_thresholds({"data_type": "bulk_ATAC"})[0] == "bulk_ATAC"
    # CP3-confirmed global_padj is honoured for the resolved modality
    modality, padj = _resolve_modality_thresholds(
        {"modality": "bulk_ATAC", "exp_context": {"global_padj": 0.1}})
    assert modality == "bulk_ATAC" and padj == 0.1
    # an explicit padj_max always wins over the resolved/default cutoff
    assert _resolve_modality_thresholds(
        {"modality": "bulk_ATAC", "padj_max": 0.2})[1] == 0.2


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
    monkeypatch.setenv("ARIA_GENOME_DIR", str(tmp_path / "genomes"))  # empty
    monkeypatch.delenv("ARIA_GENOME_FASTA", raising=False)  # deterministic
    _stage_collection(tmp_path, "JASPAR_TEST")
    res = chromatin_motifs({
        "motif_collection": "JASPAR_TEST",
        "regions": {"0": ["chr1:1-9"]},
    })
    # honest, user-facing skip (no env-var instruction); unknown assembly here
    assert res["status"] == "success" and res["ran"] is False
    assert "reference genome" in res["reason"]
    # the motif provenance is still surfaced on the skip
    assert res["motif_source"]["collection"] == "JASPAR_TEST"


def test_genome_fasta_from_env_resolves_existing_file(tmp_path, monkeypatch):
    from aria.utils import motifs
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\nACGT\n")
    monkeypatch.setenv("ARIA_GENOME_FASTA", str(fa))
    assert motifs.genome_fasta_from_env() == str(fa)
    # unset / missing file -> None (no fabrication)
    monkeypatch.setenv("ARIA_GENOME_FASTA", str(tmp_path / "nope.fa"))
    assert motifs.genome_fasta_from_env() is None
    monkeypatch.delenv("ARIA_GENOME_FASTA", raising=False)
    assert motifs.genome_fasta_from_env() is None


def test_motifs_uses_genome_fasta_env_fallback(tmp_path, monkeypatch):
    # With no genome_fasta param but ARIA_GENOME_FASTA set, the run gets PAST the
    # genome gate (it then needs snapatac2 to actually scan).
    pytest.importorskip("snapatac2")
    from aria.scripts.chromatin_motifs import chromatin_motifs
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(tmp_path))
    _stage_collection(tmp_path, "JASPAR_TEST")
    fa = tmp_path / "tiny.fa"
    fa.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    monkeypatch.setenv("ARIA_GENOME_FASTA", str(fa))
    res = chromatin_motifs({
        "motif_collection": "JASPAR_TEST",
        # no genome_fasta param -> env fallback
        "regions": {"0": ["chr1:1-9"]},
        "data_path": None, "background": ["chr1:1-9", "chr1:10-20"],
    })
    # past the genome gate: not the "no genome_fasta" skip
    assert not (res.get("ran") is False
                and "no genome_fasta" in str(res.get("reason", "")))


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


def test_regions_from_csv_caps_by_significance_not_csv_order(tmp_path):
    from aria.scripts.chromatin_motifs import _regions_from_csv
    import csv

    da_csv = tmp_path / "da.csv"
    with open(da_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["cluster", "peak", "log2fc", "padj", "significant"])
        writer.writerow(["0", "chr1:1-2", 9.0, 0.049, "True"])
        writer.writerow(["0", "chr1:3-4", 8.0, 0.048, "True"])
        for i in range(5000):
            writer.writerow([
                "0", f"chr2:{100 + i * 10}-{105 + i * 10}",
                1.0, 1e-8 + i * 1e-10, "True",
            ])
        writer.writerow(["0", "chr3:1-2", 99.0, 0.001, "False"])

    warnings = []
    groups = _regions_from_csv(str(da_csv), 5000, warnings)

    picked = groups["0"]
    assert len(picked) == 5000
    assert "chr1:1-2" not in picked
    assert "chr1:3-4" not in picked
    assert "chr3:1-2" not in picked
    assert picked[0] == "chr2:100-105"
    assert any("ranking by padj" in w for w in warnings)


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
