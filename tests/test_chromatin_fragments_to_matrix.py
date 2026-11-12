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
    assert c.validation_level in ("beta", "production")


import gzip
from pathlib import Path
import anndata as ad
import pytest


def _local_fragment_inputs(tmp_path):
    fragments = tmp_path / "a.fragments.tsv"
    fragments.write_text("chrA\t10\t20\tCELL\t1\nchrA\t30\t40\tCELL\t1\n")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chrA\t0\t50\n")
    return {"fragments_file": str(fragments), "peak_file": str(peaks),
            "genome": "custom_assembly", "min_fragments": 1,
            "output_dir": str(tmp_path / "out")}


@pytest.mark.parametrize("source", ["genome_fai", "genome_fasta", "gzip_fasta", "chrom_sizes_dict", "chrom_sizes_file", "fasta_with_index"])
def test_local_assembly_succeeds_air_gapped(tmp_path, monkeypatch, source):
    pytest.importorskip("snapatac2")
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    from aria.utils import privacy
    assert privacy.egress_allowed() is False
    params = _local_fragment_inputs(tmp_path)
    if source in {"genome_fai", "chrom_sizes_file"}:
        sizes = tmp_path / "genome.fai"
        sizes.write_text("chrA\t1000\t6\t50\t51\n")
        params["genome_fai" if source == "genome_fai" else "chrom_sizes"] = str(sizes)
    elif source == "chrom_sizes_dict":
        params["chrom_sizes"] = {"chrA": 1000}
    elif source == "gzip_fasta":
        fasta = tmp_path / "genome.fa.gz"
        with gzip.open(fasta, "wt") as handle:
            handle.write(">chrA description\n" + "ACGT" * 250 + "\n")
        params["genome_fasta"] = str(fasta)
    else:
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chrA description\n" + "ACGT" * 250 + "\n")
        params["genome_fasta"] = str(fasta)
        if source == "fasta_with_index":
            Path(str(fasta) + ".fai").write_text("chrA\t1000\t6\t50\t51\n")
    result = br.chromatin_fragments_to_matrix(params)
    assert result["status"] == "success", result
    matrix = ad.read_h5ad(result["output_path"])
    assert matrix.shape == (1, 1)
    assert matrix.X.sum() == 2
    assert result["validation_level"] == "production"
    assert br.chromatin_fragments_to_matrix(params)["resumed"] is True


def test_no_implicit_genome_lookup(tmp_path, monkeypatch):
    from aria.utils import genomes, privacy
    params = _local_fragment_inputs(tmp_path)
    monkeypatch.setattr(genomes, "snapatac2_attr", lambda genome: "GRCh38")
    monkeypatch.setattr(privacy, "egress_allowed", lambda: False)
    assert br.chromatin_fragments_to_matrix(params)["reason"] == "air_gapped_genome"


def test_local_chrom_sizes_precedence_and_validation(tmp_path):
    fai = tmp_path / "genome.fai"
    fai.write_text("chrA\t120\t0\t80\t81\nchrB\t20\t0\t20\t21\n")
    sizes, source, warnings = br._resolve_local_chrom_sizes({"chrom_sizes": {"chrA": "40"}, "genome_fai": fai})
    assert sizes == {"chrA": 40} and source == "chrom_sizes" and not warnings
    sizes, source, warnings = br._resolve_local_chrom_sizes({"chrom_sizes": {"chrA": -1}, "genome_fai": fai})
    assert sizes == {"chrA": 120, "chrB": 20} and warnings
    assert br._resolve_local_chrom_sizes({}) == (None, None, [])
    assert br._resolve_local_chrom_sizes({"chrom_sizes": {"chrA": 1.5}})[0] is None


def test_multifile_fragment_union_and_metadata(tmp_path, monkeypatch):
    pytest.importorskip("snapatac2")
    from aria.scripts import chromatin_scatac_align as sa
    params = _local_fragment_inputs(tmp_path)
    first = Path(params.pop("fragments_file"))
    second = tmp_path / "b.fragments.tsv.gz"
    with gzip.open(second, "wt") as handle:
        handle.write("chrA\t25\t35\tsecond#CELL\t1\nchrA\t5\t15\tsecond#CELL\t1\n")
    params.update(fragments_files=[str(first), str(second)], sample_ids=["first", "second"],
                  sample_metadata={"first": {"donor_id": "donor_a", "condition": "group_a", "batch": "run1"},
                                   "second": {"replicate": "donor_b", "condition": "group_b", "batch": "run2"}},
                  chrom_sizes={"chrA": 1000})
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    monkeypatch.setattr(sa, "_bgzip_tabix", lambda *args: "missing_tool")
    result = br.chromatin_fragments_to_matrix(params)
    assert result["status"] == "success", result
    matrix = ad.read_h5ad(result["output_path"])
    assert set(matrix.obs_names) == {"first#CELL", "second#CELL"}
    assert matrix.obs.loc["first#CELL", "donor_id"] == "donor_a"
    assert matrix.obs.loc["second#CELL", "donor_id"] == "donor_b"
    assert matrix.obs.loc["first#CELL", "condition"] == "group_a"
    assert matrix.obs.loc["second#CELL", "batch"] == "run2"
    assert set(matrix.obs["library_id"]) == {"first", "second"}
    assert set(matrix.obs["sample_id"]) == {"first", "second"}
    rows = (tmp_path / "out" / "combined_cohort.fragments.tsv").read_text().splitlines()
    assert [int(row.split("\t")[1]) for row in rows] == [5, 10, 25, 30]
    assert "second#second#" not in "\n".join(rows)
    assert matrix.X.sum() == 4
    assert result["n_library_annotated"] == 2


def test_fragment_union_sort_python_fallback(tmp_path, monkeypatch):
    from aria.scripts import chromatin_scatac_align as sa
    first = tmp_path / "fragments.tsv"
    first.write_text("chrA\t30\t40\tCELL\t1\nchrA\t1\t5\tCELL\t1\n")
    def unavailable(*args, **kwargs):
        raise FileNotFoundError("sort")
    monkeypatch.setattr(sa.subprocess, "run", unavailable)
    destination = tmp_path / "combined.tsv"
    result = sa._merge_library_fragments([{"library_id": "sample", "fragments_file": str(first)}], destination)
    assert result["status"] == "success", result
    assert destination.read_text().splitlines()[0] == "chrA\t1\t5\tsample#CELL\t1"


def test_fragments_files_checks_every_file(tmp_path):
    params = _local_fragment_inputs(tmp_path)
    params["fragments_files"] = [params["fragments_file"], str(tmp_path / "missing.tsv")]
    assert br.chromatin_fragments_to_matrix(params)["reason"] == "fragments_file_missing"


def test_metadata_copied_when_peak_matrix_drops_obs(tmp_path, monkeypatch):
    snap = pytest.importorskip("snapatac2")
    original = snap.pp.make_peak_matrix
    def make_matrix(*args, **kwargs):
        matrix = original(*args, **kwargs)
        matrix.obs = matrix.obs.iloc[::-1, :0].copy()
        return matrix
    monkeypatch.setattr(snap.pp, "make_peak_matrix", make_matrix)
    params = _local_fragment_inputs(tmp_path)
    Path(params["fragments_file"]).write_text("chrA\t10\t20\tsample#CELL\t1\n")
    params.update(chrom_sizes={"chrA": 1000}, library_manifest=[{
        "library_id": "sample", "sample_id": "sample", "donor_id": "donor",
        "metadata": {"condition": "group_a"}}])
    result = br.chromatin_fragments_to_matrix(params)
    assert result["status"] == "success", result
    matrix = ad.read_h5ad(result["output_path"])
    assert matrix.obs.loc["sample#CELL", "condition"] == "group_a"
    assert matrix.obs.loc["sample#CELL", "donor_id"] == "donor"
    assert "n_fragment" in matrix.obs


@pytest.mark.parametrize("change", ["fragments", "peaks", "metadata", "min_fragments", "assembly"])
def test_matrix_stage_manifest_invalidation(tmp_path, change):
    pytest.importorskip("snapatac2")
    params = _local_fragment_inputs(tmp_path)
    params["chrom_sizes"] = {"chrA": 1000}
    assert br.chromatin_fragments_to_matrix(params)["status"] == "success"
    assert br.chromatin_fragments_to_matrix(params)["resumed"] is True
    if change == "fragments":
        with open(params["fragments_file"], "a") as handle:
            handle.write("chrA\t45\t48\tCELL\t1\n")
    elif change == "peaks":
        Path(params["peak_file"]).write_text("chrA\t0\t100\n")
    elif change == "metadata":
        params["sample_metadata"] = {"sample": {"condition": "group_b"}}
    elif change == "min_fragments":
        params["min_fragments"] = 2
    else:
        params["chrom_sizes"] = {"chrA": 2000}
    result = br.chromatin_fragments_to_matrix(params)
    assert result["status"] == "success", result
    assert result.get("resumed") is not True
