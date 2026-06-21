"""B2b: scATAC fragments -> cell x peak matrix bridge (snapatac2).

The bridge unlocks the full matrix pipeline from raw fragments. It is honest:
a missing fragments file / unknown assembly / air-gapped genome / missing
snapatac2 returns a structured not-run dict, never a fabricated matrix.
"""

from aria.scripts import chromatin_fragments_to_matrix as br


def test_missing_fragments_file(tmp_path):
    out = br.chromatin_fragments_to_matrix(
        {"fragments_file": str(tmp_path / "nope.tsv.gz"), "genome": "hg38"})
    assert out["status"] == "skipped"
    assert out["reason"] == "fragments_file_missing"
    assert out["ran"] is False


def test_unknown_assembly(tmp_path, monkeypatch):
    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"x")
    from aria.utils import genomes
    monkeypatch.setattr(genomes, "snapatac2_attr", lambda g: None)
    out = br.chromatin_fragments_to_matrix(
        {"fragments_file": str(frag), "genome": "zz99"})
    assert out["status"] == "skipped"
    assert out["reason"] == "unknown_assembly"


def test_air_gapped_genome(tmp_path, monkeypatch):
    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"x")
    from aria.utils import genomes, privacy
    monkeypatch.setattr(genomes, "snapatac2_attr", lambda g: "GRCh38")
    monkeypatch.setattr(privacy, "egress_allowed", lambda: False)
    out = br.chromatin_fragments_to_matrix(
        {"fragments_file": str(frag), "genome": "hg38"})
    assert out["status"] == "skipped"
    assert out["reason"] == "air_gapped_genome"


def test_snapatac2_unavailable(tmp_path, monkeypatch):
    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"x")
    from aria.utils import genomes, privacy
    monkeypatch.setattr(genomes, "snapatac2_attr", lambda g: "GRCh38")
    monkeypatch.setattr(privacy, "egress_allowed", lambda: True)
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "snapatac2":
            raise ImportError("no snapatac2 in this env")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    out = br.chromatin_fragments_to_matrix(
        {"fragments_file": str(frag), "genome": "hg38"})
    assert out["status"] == "skipped"
    assert out["reason"] == "snapatac2_unavailable"


def test_contract_registered():
    from aria.utils.script_contracts import contract_for_script
    c = contract_for_script("aria/scripts/chromatin_fragments_to_matrix.py")
    assert c is not None
    assert c.validation_level == "beta"
