"""ADR-079 universal same-cell intake, donor consensus, and A4 regression tests."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture
def stack():
    np = pytest.importorskip("numpy")
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    sp = pytest.importorskip("scipy.sparse")
    pd = pytest.importorskip("pandas")
    return np, ad, md, sp, pd


def _adata(stack, barcodes=None, features=None, counts=None):
    np, ad, _, sp, pd = stack
    barcodes = barcodes if barcodes is not None else [f"BC{i}-1" for i in range(12)]
    features = features if features is not None else ["GENE1", "GENE2", "GENE3"]
    if counts is None:
        counts = [
            [1 + (i + j) % 4 for j in range(len(features))]
            for i in range(len(barcodes))
        ]
    obj = ad.AnnData(
        sp.csr_matrix(np.asarray(counts, dtype=np.int64)),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=features),
    )
    obj.layers["counts"] = obj.X.copy()
    return obj


def _paired(
    stack,
    path,
    *,
    barcodes=None,
    genes=None,
    peaks=None,
    rcounts=None,
    acounts=None,
    donor=None,
):
    rna = _adata(stack, barcodes, genes, rcounts)
    atac = _adata(
        stack,
        barcodes,
        peaks or ["chr1:100-200", "chr1:200-250", "chr2:100-200"],
        acounts,
    )
    rna.var["gene_ids"] = [f"id_{g}" for g in rna.var_names]
    if donor:
        rna.obs["donor_id"] = donor
        atac.obs["donor_id"] = donor
    mdata = stack[2].MuData({"rna": rna, "atac": atac})
    mdata.write_h5mu(path)
    return str(path)


def _arc_barcode(i):
    return "ACGTACGTACGT" + "".join("ACGT"[(i // (4**j)) % 4] for j in range(4)) + "-1"


def _tenx(stack, path, *, mex=False, compressed=False, types=None, duplicates=False):
    np, _, _, sp, _ = stack
    types = types or ["Gene Expression"] * 3 + ["Peaks"] * 3
    names = [
        "G1",
        "G1" if duplicates else "G2",
        "G3",
        "chr1_100_200",
        "chr1:200-250",
        "chr2:100-200",
    ][: len(types)]
    ids = [f"id{i}" for i in range(len(types))]
    barcodes = [_arc_barcode(i) for i in range(12)]
    counts = sp.csr_matrix(
        np.asarray(
            [[1 + (i + j) % 4 for j in range(len(types))] for i in range(12)],
            dtype=np.int64,
        )
    )
    matrix = counts.T.tocsc()
    if mex:
        from scipy.io import mmwrite

        path.mkdir()
        mmwrite(str(path / "matrix.mtx"), matrix)
        (path / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")
        (path / "features.tsv").write_text(
            "".join(f"{i}\t{n}\t{t}\n" for i, n, t in zip(ids, names, types))
        )
        if compressed:
            for f in path.iterdir():
                with gzip.open(str(f) + ".gz", "wb") as handle:
                    handle.write(f.read_bytes())
                f.unlink()
    else:
        h5py = pytest.importorskip("h5py")
        with h5py.File(path, "w") as handle:
            group = handle.create_group("matrix")
            for key, value in (
                ("data", matrix.data),
                ("indices", matrix.indices),
                ("indptr", matrix.indptr),
                ("shape", matrix.shape),
            ):
                group[key] = value
            group["barcodes"] = np.asarray(barcodes, dtype="S")
            features = group.create_group("features")
            for key, value in (("name", names), ("id", ids), ("feature_type", types)):
                features[key] = np.asarray(value, dtype="S")
    return str(path), counts


@pytest.mark.parametrize(
    "mex,compressed", [(False, False), (True, False), (True, True)]
)
def test_detector_and_universal_reader(stack, tmp_path, mex, compressed):
    from aria.utils.assay_detector import AssayDetector
    from aria.utils.mudata_io import (
        read_h5mu_paired,
        load_or_convert_paired_multiome_to_h5mu,
    )

    path, counts = _tenx(
        stack,
        tmp_path / ("mex" if mex else "arc.h5"),
        mex=mex,
        compressed=compressed,
        duplicates=True,
    )
    matrix = (
        Path(path) / ("matrix.mtx.gz" if compressed else "matrix.mtx")
        if mex
        else Path(path)
    )
    detection = AssayDetector().detect_file(matrix)
    assert detection.evidence["format"] == (
        "10x_multiome_mtx" if mex else "10x_multiome_h5"
    )
    assert detection.evidence["paired_rna"] is True
    for source in {path, str(matrix)}:
        result = read_h5mu_paired(source)
        assert result["status"] == "success", result
        assert result["n_shared_cells"] == 12
        assert result["rna_adata"].var_names.is_unique
        assert list(result["rna_adata"].var["gene_ids"]) == ["id0", "id1", "id2"]
        assert list(result["atac_adata"].var_names) == [
            "chr1:100-200",
            "chr1:200-250",
            "chr2:100-200",
        ]
        stack[0].testing.assert_array_equal(
            result["rna_counts_matrix"].toarray(), counts[:, :3].toarray()
        )
        assert result["atac_counts_matrix"].dtype.kind in "iu"
    converted = tmp_path / "converted.h5mu"
    out = load_or_convert_paired_multiome_to_h5mu(path, str(converted))
    assert out["status"] == "success", out
    reread = read_h5mu_paired(str(converted))
    assert reread["status"] == "success"
    assert reread["exact_same_order"] is True


@pytest.mark.parametrize("mex", [False, True])
@pytest.mark.parametrize("feature_type", ["Gene Expression", "Peaks"])
def test_single_modality_fails_closed_and_atac_detection(
    stack, tmp_path, mex, feature_type
):
    from aria.utils.assay_detector import AssayDetector
    from aria.utils.mudata_io import read_h5mu_paired

    path, _ = _tenx(
        stack,
        tmp_path / ("mex" if mex else "single.h5"),
        mex=mex,
        types=[feature_type] * 3,
    )
    # ATAC-only features must be genomic to isolate the paired-modality assertion.
    if feature_type == "Peaks":
        if mex:
            (Path(path) / "features.tsv").write_text(
                "".join(f"p{i}\tchr1:{i*100}-{i*100+50}\tPeaks\n" for i in range(3))
            )
        else:
            import h5py

            with h5py.File(path, "a") as h:
                del h["matrix/features/name"]
                h["matrix/features/name"] = stack[0].asarray(
                    [f"chr1:{i*100}-{i*100+50}" for i in range(3)], dtype="S"
                )
        detection = AssayDetector().detect_file(
            Path(path) / "matrix.mtx" if mex else path
        )
        assert detection.evidence["format"] == ("10x_mtx" if mex else "10x_atac_h5")
        assert detection.evidence["paired_rna"] is False
    result = read_h5mu_paired(path)
    assert result["status"] == "error"
    assert result["error_type"] == "MissingPairedModality"


def test_h5ad_rejected_even_if_renamed(stack, tmp_path):
    from aria.utils.mudata_io import read_h5mu_paired

    path = tmp_path / "single.h5ad"
    _adata(stack).write_h5ad(path)
    assert read_h5mu_paired(str(path))["error_type"] == "NotAnH5mu"
    path.rename(tmp_path / "renamed.h5")
    assert read_h5mu_paired(str(tmp_path / "renamed.h5"))["error_type"] == "NotAnH5mu"


def test_h5_interval_fallback_and_count_rejection(stack, tmp_path):
    import h5py
    from aria.utils.mudata_io import read_h5mu_paired

    path, _ = _tenx(stack, tmp_path / "arc.h5")
    with h5py.File(path, "a") as h:
        h["matrix/features/name"][3] = b"peak_1"
        h["matrix/features/interval"] = stack[0].asarray(
            ["", "", "", "chr1:100-200", "", ""], dtype="S"
        )
    assert read_h5mu_paired(path)["status"] == "success"
    with h5py.File(path, "a") as h:
        values = h["matrix/data"][:].astype(float)
        del h["matrix/data"]
        values[0] = 1.5
        h["matrix/data"] = values
    assert read_h5mu_paired(path)["error_type"] == "NonIntegerCounts"


@pytest.mark.parametrize("mex", [False, True])
def test_assemble_arc_fragments_and_resume(stack, tmp_path, mex):
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    path, _ = _tenx(stack, tmp_path / ("mex" if mex else "arc.h5"), mex=mex)
    fragments = tmp_path / "fragments.tsv.gz"
    with gzip.open(fragments, "wt") as handle:
        handle.write(
            f"chr1\t100\t150\t{_arc_barcode(0)}\t3\nchrM\t0\t10\t{_arc_barcode(0)}\t1\n"
        )
    params = {
        "arc_path": path,
        "fragments_path": str(fragments),
        "output_dir": str(tmp_path / "out"),
    }
    result = integration_assemble_multiome(params)
    assert result["status"] == "success", result
    assert result["mode"] == "cellranger_arc" and result["resumed"] is False
    out = stack[2].read_h5mu(result["output_path"])
    assert out.mod["atac"].obs.loc[_arc_barcode(0), "atac_fragments"] == 4
    assert out.mod["atac"].obs.loc[_arc_barcode(0), "atac_mito_frac"] == 0.25
    assert stack[0].isnan(out.mod["atac"].obs.loc[_arc_barcode(1), "atac_mito_frac"])
    for mod in out.mod.values():
        assert mod.X.dtype.kind in "iu"
        assert (mod.X != mod.layers["counts"]).nnz == 0
    assert integration_assemble_multiome(params)["resumed"] is True
    assert (
        integration_assemble_multiome(dict(params, min_shared_cells=11))["resumed"]
        is False
    )
    assert (
        integration_assemble_multiome(dict(params, min_shared_cells=13))["error_type"]
        == "InsufficientSharedCells"
    )
    assert integration_assemble_multiome(dict(params, resume=False))["resumed"] is False
    if mex:
        (Path(path) / "features.tsv").write_text(
            (Path(path) / "features.tsv").read_text().replace("G1", "GZ")
        )
        assert integration_assemble_multiome(params)["resumed"] is False


def _matrix_paths(stack, tmp_path, rbarcodes, abarcodes):
    rna, atac = tmp_path / "rna.h5ad", tmp_path / "atac.h5ad"
    _adata(stack, rbarcodes).write_h5ad(rna)
    _adata(stack, abarcodes, ["chr1:100-200", "chr2:100-200"]).write_h5ad(atac)
    return {
        "rna_path": str(rna),
        "atac_path": str(atac),
        "min_shared_cells": 2,
        "output_dir": str(tmp_path / "out"),
    }


@pytest.mark.parametrize("orientation", ["direct", "reverse_complement"])
def test_whitelist_translation_orientation(stack, tmp_path, orientation):
    from aria.scripts.integration_assemble_multiome import (
        integration_assemble_multiome,
        _revcomp,
    )

    gex, atac = ["AAAA", "CCCC", "ACAC"], ["AAAC", "CCCA", "ACCG"]
    observed = atac if orientation == "direct" else [_revcomp(b) for b in atac]
    params = _matrix_paths(
        stack, tmp_path, [b + "-1" for b in gex], [b + "-1" for b in observed]
    )
    gp, ap = tmp_path / "gex.txt.gz", tmp_path / "atac.txt"
    with gzip.open(gp, "wt") as handle:
        handle.write("\n".join(gex) + "\n")
    ap.write_text("\n".join(atac) + "\n")
    params.update(
        barcode_namespace="10x_arc_whitelist_pair",
        gex_whitelist_path=str(gp),
        atac_whitelist_path=str(ap),
    )
    result = integration_assemble_multiome(params)
    assert result["status"] == "success", result
    assert result["barcode_translation"] == {
        "mode": "10x_arc_whitelist_pair",
        "orientation": orientation,
        "n_translated": 3,
    }
    out = stack[2].read_h5mu(result["output_path"])
    assert (
        list(out.obs_names)
        == [b + "-1" for b in gex]
        == list(out.mod["atac"].obs_names)
    )
    result = integration_assemble_multiome(
        dict(params, atac_barcode_orientation=orientation)
    )
    assert result["status"] == "success"
    gp.unlink()
    with gzip.open(gp, "wt") as handle:
        handle.write("AAAA\nAAAA\nACAC\n")
    assert (
        integration_assemble_multiome(params)["error_type"] == "NonBijectiveBarcodeMap"
    )


@pytest.mark.parametrize("delimiter,header", [("\t", True), (",", False)])
def test_barcode_map_bijection_and_overlap_denominator(
    stack, tmp_path, delimiter, header
):
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    params = _matrix_paths(stack, tmp_path, ["R1", "R2"], ["A1", "A2", "unmapped"])
    mapping = tmp_path / "mapping.tsv"
    text = delimiter.join(["rna_barcode", "atac_barcode"]) + "\n" if header else ""
    text += f"R1{delimiter}A1\nR2{delimiter}A2\n"
    mapping.write_text(text)
    params.update(barcode_namespace="barcode_map_tsv", barcode_map_path=str(mapping))
    result = integration_assemble_multiome(params)
    assert result["status"] == "success", result
    assert result["n_atac_cells"] == 3
    assert result["shared_fraction_atac"] == 2 / 3
    assert (
        integration_assemble_multiome(dict(params, min_shared_barcode_fraction=0.9))[
            "error_type"
        ]
        == "LowBarcodeOverlap"
    )
    for extra in (f"R1{delimiter}A3\n", f"R3{delimiter}A1\n"):
        mapping.write_text(text + extra)
        result = integration_assemble_multiome(params)
        assert result["error_type"] == "NonBijectiveBarcodeMap"


def test_shared_suffix_alignment_unique_only(stack, tmp_path):
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    params = _matrix_paths(stack, tmp_path, ["A-1", "B-1"], ["B", "A"])
    result = integration_assemble_multiome(params)
    assert result["status"] == "success", result
    assert list(stack[2].read_h5mu(result["output_path"]).mod["atac"].obs_names) == [
        "A-1",
        "B-1",
    ]
    assert (
        integration_assemble_multiome(dict(params, normalize_10x_suffix=False))[
            "error_type"
        ]
        == "NoSharedBarcodes"
    )
    _adata(stack, ["A", "A-1"], ["chr1:1-2"]).write_h5ad(params["atac_path"])
    result = integration_assemble_multiome(
        dict(params, min_shared_cells=1, min_shared_barcode_fraction=0)
    )
    assert result["status"] == "success", result
    assert result["n_shared_cells"] == 1
    assert result["barcode_translation"]["n_translated"] == 0


def test_fragment_half_open_overlap_and_multiplicity(stack, tmp_path):
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    rna = tmp_path / "rna.h5ad"
    _adata(stack, ["A", "B"]).write_h5ad(rna)
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t100\t200\nchr1\t200\t300\nchr1\t120\t130\n")
    fragments = tmp_path / "fragments.tsv"
    fragments.write_text(
        "chr1\t90\t100\tA\t10\nchr1\t100\t120\tA\t2\nchr1\t199\t201\tA\t3\nchr1\t300\t310\tA\t5\nchr1\t125\t126\tB\t4\n"
    )
    result = integration_assemble_multiome(
        {
            "rna_path": str(rna),
            "fragments_path": str(fragments),
            "peaks_path": str(peaks),
            "min_shared_cells": 2,
            "output_dir": str(tmp_path / "out"),
        }
    )
    assert result["status"] == "success", result
    counts = stack[2].read_h5mu(result["output_path"]).mod["atac"].X
    stack[0].testing.assert_array_equal(counts.toarray(), [[5, 3, 0], [4, 0, 4]])


@pytest.mark.parametrize("kind", ["MissingDependency", "MissingReference"])
def test_raw_fastq_fail_closed(tmp_path, kind):
    pytest.importorskip("mudata")
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    params = {
        "mode": "raw_multiome_fastq",
        "platform": "10x_multiome_arc_v1",
        "rna_fastqs": ["R1", "R2"],
        "atac_fastqs": ["A1", "A2", "I1"],
        "kb_bin": "/nonexistent/kb",
        "chromap_bin": "/nonexistent/chromap",
        "output_dir": str(tmp_path),
    }
    if kind == "MissingReference":
        params.update(kb_bin=sys.executable, chromap_bin=sys.executable)
    result = integration_assemble_multiome(params)
    assert result["status"] == "error" and result["error_type"] == kind, result
    assert not list(tmp_path.glob("*.h5mu"))
    assert not list(tmp_path.glob("*.stage_manifest.json"))


def test_raw_fastq_delegated_outputs(stack, tmp_path):
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    rna = tmp_path / "rna.h5ad"
    _adata(stack, ["A", "B"]).write_h5ad(rna)
    fragments, peaks = tmp_path / "fragments.tsv", tmp_path / "peaks.bed"
    fragments.write_text("chr1\t100\t110\tA\t2\nchr1\t110\t120\tB\t3\n")
    peaks.write_text("chr1\t100\t200\n")
    result = integration_assemble_multiome(
        {
            "mode": "raw_multiome_fastq",
            "platform": "custom_paired",
            "rna_counts_out": str(rna),
            "atac_fragments_out": str(fragments),
            "peaks_path": str(peaks),
            "min_shared_cells": 2,
            "output_dir": str(tmp_path / "out"),
        }
    )
    assert result["status"] == "success", result
    assert result["mode"] == "raw_multiome_fastq"


def test_three_donor_merge_exact_consensus_counts_metadata_and_a4(stack, tmp_path):
    from aria.scripts.integration_merge_cohort import integration_merge_cohort

    peaks = [
        ["chr1:100-200", "chr1:200-250", "chr2:10-20"],
        ["chr1:150-210", "chr1:250-300", "chr2:20-30"],
        ["chr1:290-350", "chr2:15-25", "chr3:1-5"],
    ]
    paths = [
        _paired(
            stack,
            tmp_path / f"s{i}.h5mu",
            barcodes=["A", "B"],
            peaks=p,
            acounts=[[1, 2, 3], [4, 5, 6]],
            donor="cell_donor" if i == 0 else None,
        )
        for i, p in enumerate(peaks)
    ]
    params = {
        "inputs": paths,
        "sample_ids": ["s0", "s1", "s2"],
        "output_dir": str(tmp_path / "out"),
        "sample_metadata": {
            "s0": {"donor_id": "fallback", "condition": "control"},
            "s1": {"donor_id": "donor1", "condition": "treated"},
            "s2.h5mu": {"donor_id": "donor2", "condition": "treated"},
        },
    }
    result = integration_merge_cohort(params)
    assert result["status"] == "success", result
    assert result["n_samples"] == 3 and result["n_cells_total"] == 6
    assert result["n_consensus_peaks"] == 4
    out = stack[2].read_h5mu(result["output_path"])
    assert list(out.mod["atac"].var_names) == [
        "chr1:100-250",
        "chr1:250-350",
        "chr2:10-30",
        "chr3:1-5",
    ]
    # chr2 overlaps transitively via donor3; chr1 adjacency at 250 stays separate.
    expected = [
        [3, 0, 3, 0],
        [9, 0, 6, 0],
        [1, 2, 3, 0],
        [4, 5, 6, 0],
        [0, 1, 2, 3],
        [0, 4, 5, 6],
    ]
    stack[0].testing.assert_array_equal(out.mod["atac"].X.toarray(), expected)
    assert out.mod["atac"].X.dtype.kind in "iu"
    expected_bcs = [f"s{i}#{bc}" for i in range(3) for bc in ("A", "B")]
    for obs in (out.obs, out.mod["rna"].obs, out.mod["atac"].obs):
        assert list(obs.index) == expected_bcs
        assert (
            list(obs["donor_id"].astype(str))
            == ["cell_donor"] * 2 + ["donor1"] * 2 + ["donor2"] * 2
        )
        assert list(obs["condition"].astype(str)) == ["control"] * 2 + ["treated"] * 4
    manifest = json.loads(Path(result["stage_manifest_path"]).read_text())
    assert len(manifest["inputs"]) == 3
    assert all(len(entry["sha256"]) == 64 for entry in manifest["inputs"])
    assert integration_merge_cohort(params)["resumed"] is True
    assert (
        integration_merge_cohort(dict(params, rna_gene_join="outer"))["resumed"]
        is False
    )
    assert integration_merge_cohort(dict(params, resume=False))["resumed"] is False
    m = stack[2].read_h5mu(paths[0])
    m.mod["rna"].X[0, 0] += 1
    m.mod["rna"].layers["counts"] = m.mod["rna"].X.copy()
    m.write_h5mu(paths[0])
    assert integration_merge_cohort(params)["resumed"] is False
    Path(result["output_path"]).unlink()
    assert integration_merge_cohort(params)["resumed"] is False


@pytest.mark.parametrize("join", ["inner", "outer"])
def test_gene_alignment_and_empty_intersection(stack, tmp_path, join):
    from aria.scripts.integration_merge_cohort import integration_merge_cohort

    a = _paired(
        stack,
        tmp_path / "a.h5mu",
        barcodes=["a#A"],
        genes=["G1", "G2", "G3", "G4"],
        rcounts=[[1, 2, 3, 4]],
    )
    b = _paired(
        stack,
        tmp_path / "b.h5mu",
        barcodes=["A"],
        genes=["G3", "G2", "G1", "G5"],
        rcounts=[[5, 6, 7, 8]],
    )
    params = {
        "inputs": [a, b],
        "sample_ids": ["a", "b"],
        "rna_gene_join": join,
        "output_dir": str(tmp_path / "out"),
    }
    result = integration_merge_cohort(params)
    assert result["status"] == "success", result
    out = stack[2].read_h5mu(result["output_path"])
    assert list(out.obs_names) == ["a#A", "b#A"]
    rna = out.mod["rna"]
    assert list(rna.var_names) == (
        ["G1", "G2", "G3"] if join == "inner" else ["G1", "G2", "G3", "G4", "G5"]
    )
    stack[0].testing.assert_array_equal(
        rna.X.toarray(),
        (
            [[1, 2, 3], [7, 6, 5]]
            if join == "inner"
            else [[1, 2, 3, 4, 0], [7, 6, 5, 0, 8]]
        ),
    )
    assert list(rna.var["gene_ids"]) == [f"id_{g}" for g in rna.var_names]
    _paired(stack, Path(b), barcodes=["A"], genes=["OTHER1", "OTHER2", "G1"])
    result = integration_merge_cohort(dict(params, rna_gene_join="inner"))
    assert result["error_type"] == "EmptyGeneIntersection"


def _references(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "".join(
            f'chr1\ttest\tgene\t{200+i*200}\t{350+i*200}\t.\t+\t.\tgene_id "id{i}"; gene_name "GENE{i+1}";\n'
            for i in range(3)
        )
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 500 + "\n")
    Path(str(fasta) + ".fai").write_text("chr1\t2000\t6\t2000\t2001\n")
    motif = tmp_path / "motif.meme"
    motif.write_text(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF M1 GENE1\nletter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0\n0.97 0.01 0.01 0.01\n0.01 0.97 0.01 0.01\n0.01 0.01 0.97 0.01\n0.01 0.01 0.01 0.97\n"
    )
    return gtf, fasta, motif


@pytest.mark.parametrize(
    "stage", ["integration_paired_qc", "integration_peak2gene", "integration_eregulons"]
)
def test_existing_stage_a4_resume_input_parameter_reference_and_output_invalidation(
    stack, tmp_path, stage
):
    import importlib

    call = getattr(importlib.import_module("aria.scripts." + stage), stage)
    source = _paired(
        stack,
        tmp_path / "paired.h5mu",
        peaks=["chr1:200-250", "chr1:400-450", "chr1:600-650"],
    )
    gtf, fasta, motif = _references(tmp_path)
    params = {
        "h5mu_path": source,
        "output_dir": str(tmp_path / "out"),
        "min_shared_cells": 2,
        "min_cells_detected": 1,
        "min_genes": 0,
        "min_peaks": 0,
        "min_rna_counts": 0,
        "min_atac_counts": 0,
        "gtf_path": str(gtf),
        "genome_fasta": str(fasta),
        "motif_file": str(motif),
        "background_null": False,
        "use_fdr": False,
        "n_null_permutations": 2,
        "scanner_backend": "aria_numpy_dp",
    }
    fresh = call(params)
    assert fresh["status"] == "success", fresh
    assert fresh["resumed"] is False
    manifest = Path(params["output_dir"]) / (stage + ".stage_manifest.json")
    receipt = Path(params["output_dir"]) / ("." + stage + ".result.json")
    assert manifest.is_file() and receipt.is_file()
    assert call(params)["resumed"] is True
    assert call(dict(params, resume=False))["resumed"] is False
    changed = (
        dict(params, min_genes=1)
        if stage.endswith("qc")
        else dict(params, distance_kb=250)
    )
    assert call(changed)["resumed"] is False
    assert call(changed)["resumed"] is True
    m = stack[2].read_h5mu(source)
    m.mod["rna"].X[0, 0] += 1
    m.mod["rna"].layers["counts"] = m.mod["rna"].X.copy()
    m.write_h5mu(source)
    assert call(changed)["resumed"] is False
    if stage.endswith("qc"):
        reference = tmp_path / "reference.h5ad"
        ref = _adata(stack)
        ref.obs["subclass"] = "type_a"
        ref.write_h5ad(reference)
        changed = dict(changed, reference_h5ad=str(reference))
        assert call(changed)["status"] == "success"
        assert call(changed)["resumed"] is True
        ref.obs["subclass"] = "type_b"
        ref.write_h5ad(reference)
    else:
        gtf.write_text(gtf.read_text() + "# updated reference\n")
    assert call(changed)["resumed"] is False
    target = (
        fresh.get("triplets_path")
        if stage.endswith("eregulons")
        else fresh["output_path"]
    )
    Path(target).unlink()
    assert call(changed)["resumed"] is False
    manifest.write_text("corrupt")
    assert call(changed)["resumed"] is False


def test_qc_mex_component_hashing(stack, tmp_path):
    from aria.scripts.integration_paired_qc import integration_paired_qc

    source, _ = _tenx(stack, tmp_path / "mex", mex=True, compressed=True)
    params = {
        "h5mu_path": source,
        "output_dir": str(tmp_path / "out"),
        "min_genes": 0,
        "min_peaks": 0,
        "min_rna_counts": 0,
        "min_atac_counts": 0,
    }
    result = integration_paired_qc(params)
    assert result["status"] == "success", result
    assert integration_paired_qc(params)["resumed"] is True
    manifest = json.loads(
        (
            Path(params["output_dir"]) / "integration_paired_qc.stage_manifest.json"
        ).read_text()
    )
    assert len(manifest["inputs"]) == 3
    with gzip.open(Path(source) / "features.tsv.gz", "rt") as handle:
        text = handle.read()
    with gzip.open(Path(source) / "features.tsv.gz", "wt") as handle:
        handle.write(text.replace("G1", "GZ"))
    assert integration_paired_qc(params)["resumed"] is False


@pytest.mark.parametrize(
    "script", ["integration_assemble_multiome", "integration_merge_cohort"]
)
def test_ipc_scripts_bootstrap_from_other_directory(tmp_path, script):
    input_path, output_path = tmp_path / "input.json", tmp_path / "output.json"
    input_path.write_text("{}")
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "aria" / "scripts" / (script + ".py")),
            str(input_path),
            str(output_path),
        ],
        text=True,
        capture_output=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    receipt = json.loads(output_path.read_text())
    assert receipt["status"] == "error" and receipt["error_type"]


@pytest.mark.parametrize("mex", [False, True])
def test_detector_scans_feature_types_past_prefix(tmp_path, mex):
    np = pytest.importorskip("numpy")
    h5py = pytest.importorskip("h5py")
    from aria.utils.assay_detector import AssayDetector

    types = ["Gene Expression"] * 8200 + ["Peaks"]
    if mex:
        directory = tmp_path / "filtered_feature_bc_matrix"
        directory.mkdir()
        (directory / "barcodes.tsv").write_text(_arc_barcode(0) + "\n")
        (directory / "features.tsv").write_text(
            "".join(f"id{i}\tname{i}\t{t}\n" for i, t in enumerate(types))
        )
        path = directory / "matrix.mtx"
        path.write_text(
            "%%MatrixMarket matrix coordinate integer general\n8201 1 1\n8201 1 1\n"
        )
    else:
        path = tmp_path / "late_peaks.h5"
        with h5py.File(path, "w") as h:
            group = h.create_group("matrix")
            group["barcodes"] = [b"ACGTACGT-1"]
            group["data"] = [1]
            group.create_group("features")["feature_type"] = np.asarray(
                types, dtype="S"
            )
    result = AssayDetector().detect_file(path)
    assert result.evidence["format"] == (
        "10x_multiome_mtx" if mex else "10x_multiome_h5"
    )
    assert result.evidence["paired_rna"] is True


def test_merge_mixed_universal_containers(stack, tmp_path):
    from aria.scripts.integration_merge_cohort import integration_merge_cohort

    a, _ = _tenx(stack, tmp_path / "arc.h5")
    b, _ = _tenx(stack, tmp_path / "mex", mex=True)
    c = _paired(stack, tmp_path / "paired.h5mu", genes=["G1", "G2", "G3"])
    result = integration_merge_cohort(
        {"inputs": [a, b, c], "output_dir": str(tmp_path / "out")}
    )
    assert result["status"] == "success", result
    assert result["n_cells_total"] == 36 and result["n_genes"] == 3
    manifest = json.loads(Path(result["stage_manifest_path"]).read_text())
    assert len(manifest["inputs"]) == 5


def test_merge_verified_shared_subset_and_low_overlap(stack, tmp_path):
    from aria.scripts.integration_merge_cohort import integration_merge_cohort

    a = _paired(stack, tmp_path / "a.h5mu", barcodes=["A", "B"])
    rna = _adata(stack, ["A", "B", "C", "D", "E"])
    atac = _adata(stack, ["A", "X", "Y", "Z", "W"], ["chr1:100-200"])
    b = tmp_path / "b.h5mu"
    stack[2].MuData({"rna": rna, "atac": atac}).write_h5mu(b)
    params = {
        "inputs": [a, str(b)],
        "sample_ids": ["a", "b"],
        "output_dir": str(tmp_path / "out"),
    }
    assert integration_merge_cohort(params)["error_type"] == "LowBarcodeOverlap"
    result = integration_merge_cohort(dict(params, min_shared_barcode_fraction=0.2))
    assert result["status"] == "success", result
    out = stack[2].read_h5mu(result["output_path"])
    assert result["per_sample_cells"] == {"a": 2, "b": 1}
    for names in (out.obs_names, out.mod["rna"].obs_names, out.mod["atac"].obs_names):
        assert list(names) == ["a#A", "a#B", "b#A"]


def test_assembly_selects_real_counts_layer_and_rejects_normalized_only(
    stack, tmp_path
):
    from aria.scripts.integration_assemble_multiome import integration_assemble_multiome

    params = _matrix_paths(stack, tmp_path, ["A", "B"], ["A", "B"])
    rna = stack[1].read_h5ad(params["rna_path"])
    expected = rna.layers["counts"].copy()
    rna.X = rna.X.astype(float) / 3.5
    rna.write_h5ad(params["rna_path"])
    result = integration_assemble_multiome(params)
    assert result["status"] == "success", result
    stack[0].testing.assert_array_equal(
        stack[2].read_h5mu(result["output_path"]).mod["rna"].X.toarray(),
        expected.toarray(),
    )
    del rna.layers["counts"]
    rna.write_h5ad(params["rna_path"])
    assert integration_assemble_multiome(params)["error_type"] == "NonIntegerCounts"


@pytest.mark.parametrize("stage", ["integration_peak2gene", "integration_eregulons"])
def test_separate_matrix_stage_resume_and_hash_all_references(stack, tmp_path, stage):
    import importlib

    call = getattr(importlib.import_module("aria.scripts." + stage), stage)
    rna, atac = tmp_path / "rna.h5ad", tmp_path / "atac.h5ad"
    _adata(stack).write_h5ad(rna)
    _adata(stack, features=["chr1:200-250", "chr1:400-450", "chr1:600-650"]).write_h5ad(
        atac
    )
    gtf, fasta, motif = _references(tmp_path)
    params = {
        "rna_files": [str(rna)],
        "atac_files": [str(atac)],
        "gtf_path": str(gtf),
        "genome_fasta": str(fasta),
        "motif_file": str(motif),
        "min_shared_cells": 2,
        "min_cells_detected": 1,
        "n_null_permutations": 2,
        "scanner_backend": "aria_numpy_dp",
        "output_dir": str(tmp_path / "out"),
    }
    result = call(params)
    assert result["status"] == "success", result
    assert call(params)["resumed"] is True
    manifest = json.loads(Path(result["stage_manifest_path"]).read_text())
    roles = {entry["role"] for entry in manifest["inputs"]}
    assert {"rna:0", "atac:0", "gtf", "genome_fasta", "fasta_index"} <= roles
    if stage.endswith("eregulons"):
        assert "motifs" in roles
        motif.write_text(motif.read_text() + "\n# reference update\n")
        assert call(params)["resumed"] is False
    sequence = fasta.read_text()
    fasta.write_text(sequence.replace("ACGT", "AGGT", 1))
    assert call(params)["resumed"] is False
    r = stack[1].read_h5ad(rna)
    r.X[0, 0] += 1
    r.write_h5ad(rna)
    assert call(params)["resumed"] is False
