"""scATAC P4.3 — TOBIAS footprinting driver guards.

Pure helpers + the honest not-run gate are unit-tested in aria-env (no pysam, no
TOBIAS). The full ATACorrect->ScoreBigwig->BINDetect run is exercised separately in
the dedicated aria-tobias-env on real fragments. No fabrication: when TOBIAS or any
required asset is absent the driver returns ran=False with a concrete reason (ADR-002
/ W2.2: never an uncorrected footprint or an invented score).
"""

from __future__ import annotations

import pytest
import json
import sys
import types
from pathlib import Path


def test_chrom_sizes_from_fai(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _chrom_sizes_from_fai

    fai = tmp_path / "g.fa.fai"
    fai.write_text("chr1\t248956422\t112\t70\t71\nchr2\t242193529\t1\t70\t71\n",
                   encoding="utf-8")
    sizes = _chrom_sizes_from_fai(str(fai))
    assert sizes == {"chr1": 248956422, "chr2": 242193529}


def test_load_barcode_groups_skips_header(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _load_barcode_groups

    tsv = tmp_path / "groups.tsv"
    tsv.write_text("barcode\tgroup\nBC1-1\tMonocyte\nBC2-1\tT_cell\nBC3-1\tMonocyte\n",
                   encoding="utf-8")
    groups = _load_barcode_groups(str(tsv))
    assert groups == {"BC1-1": "Monocyte", "BC2-1": "T_cell", "BC3-1": "Monocyte"}


def test_fragment_cut_reads_two_reads_at_tn5_sites():
    from aria.scripts.chromatin_footprint_tobias import _fragment_cut_reads

    reads = list(_fragment_cut_reads("chr1", 100, 200, "BC:0", read_len=50))
    assert len(reads) == 2
    fwd, rev = reads
    assert fwd["is_reverse"] is False and fwd["pos"] == 100
    assert rev["is_reverse"] is True and rev["pos"] == 150  # end - read_len
    # degenerate fragment (end <= start) yields nothing
    assert list(_fragment_cut_reads("chr1", 200, 200, "BC:1")) == []


def test_fragment_cut_reads_clamps_read_len_to_fragment():
    from aria.scripts.chromatin_footprint_tobias import _fragment_cut_reads

    reads = list(_fragment_cut_reads("chr1", 100, 130, "BC:0", read_len=50))
    # read_len clamped to the 30bp fragment so reads never exceed the fragment
    assert all(r["length"] == 30 for r in reads)
    assert reads[1]["pos"] == 100  # max(start, end - 30)


def test_summarize_bindetect_top_differential_per_group(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import summarize_bindetect

    res = tmp_path / "bindetect_results.txt"
    res.write_text(
        "name\ttotal_tfbs\tMonocyte_T_cell_change\tMonocyte_T_cell_pvalue\n"
        "CEBPB\t97\t1.10\t1e-100\n"        # strongly toward Monocyte
        "TCF7\t500\t-0.40\t1e-50\n"          # toward T_cell
        "NOISE\t5\t0.90\t1e-80\n"            # below min_sites -> excluded
        "FLAT\t300\t0.01\t0.4\n",            # above p threshold -> excluded
        encoding="utf-8")
    s = summarize_bindetect(str(res), "Monocyte", "T_cell")
    assert s["parsed"] is True
    # B7 interim: descriptive candidate ranking, not an FDR-significance count.
    assert s["n_motifs_tested"] == 4 and s["n_ranked_candidates"] == 2
    assert "n_significant" not in s
    assert s["ranking_basis"]["fdr_controlled"] is False
    assert s["top_toward_Monocyte"][0]["tf"] == "CEBPB"
    assert s["top_toward_T_cell"][0]["tf"] == "TCF7"


def test_summarize_bindetect_honest_when_columns_absent(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import summarize_bindetect

    res = tmp_path / "bad.txt"
    res.write_text("name\tother\nX\t1\n", encoding="utf-8")
    s = summarize_bindetect(str(res), "Monocyte", "T_cell")
    assert s["parsed"] is False and "missing column" in s["reason"]


def test_top_tfs_both_directions_deduped():
    from aria.scripts.chromatin_footprint_tobias import _top_tfs

    summary = {
        "top_toward_Monocyte": [{"tf": "CEBPB"}, {"tf": "CEBPA"}, {"tf": "SPIC"}],
        "top_toward_T_cell": [{"tf": "EGR1"}, {"tf": "CEBPB"}],  # CEBPB dup
    }
    tfs = _top_tfs(summary, "Monocyte", "T_cell", n=2)
    # top-2 each direction, deduped, order preserved (no biology dropped one-sided)
    assert tfs == ["CEBPB", "CEBPA", "EGR1"]


def test_aggregate_plots_skips_tfs_without_beds(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _aggregate_plots

    # No BINDetect bed dirs present -> honest empty result, no TOBIAS call needed.
    bindetect = tmp_path / "bindetect"
    bindetect.mkdir()
    out = _aggregate_plots(bindetect, {"A": "a.bw", "B": "b.bw"},
                           ["CEBPB", "EGR1"], tmp_path)
    assert out == {}


def test_driver_honest_not_run_without_tobias(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    # TOBIAS absent (the aria-env reality) -> honest ran=False, never a fake footprint.
    monkeypatch.setattr(mod.shutil, "which", lambda name: None)
    res = mod.chromatin_footprint_tobias({
        "fragments_file": str(tmp_path / "f.tsv.gz"),
        "genome_fasta": str(tmp_path / "g.fa"),
        "peaks_bed": str(tmp_path / "p.bed"),
        "motif_meme": str(tmp_path / "m.meme"),
        "barcode_groups": str(tmp_path / "g.tsv"),
        "group_a": "Monocyte", "group_b": "T_cell",
        "output_dir": str(tmp_path / "out"),
    })
    assert res["ran"] is False
    assert "TOBIAS not installed" in res["reason"]


def test_driver_honest_not_run_missing_asset(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    # TOBIAS present but a required asset is missing -> concrete reason, still no run.
    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/TOBIAS")
    res = mod.chromatin_footprint_tobias({
        "fragments_file": str(tmp_path / "absent_fragments.tsv.gz"),
        "genome_fasta": str(tmp_path / "g.fa"),
        "peaks_bed": str(tmp_path / "p.bed"),
        "motif_meme": str(tmp_path / "m.meme"),
        "barcode_groups": str(tmp_path / "g.tsv"),
        "group_a": "Monocyte", "group_b": "T_cell",
        "output_dir": str(tmp_path / "out"),
    })
    assert res["ran"] is False
    assert "fragments_file missing" in res["reason"]


def test_dispatch_entrypoint_supports_bulk_mode_skip(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    monkeypatch.setattr(mod.shutil, "which", lambda name: None)
    res = mod.chromatin_footprint_tobias_dispatch({
        "mode": "bulk",
        "condition_bams": {"A": [str(tmp_path / "a.bam")],
                           "B": [str(tmp_path / "b.bam")]},
        "genome_fasta": str(tmp_path / "g.fa"),
        "peaks_bed": str(tmp_path / "p.bed"),
        "motif_meme": str(tmp_path / "m.meme"),
        "group_a": "A",
        "group_b": "B",
        "output_dir": str(tmp_path / "out"),
    })
    assert res["ran"] is False
    assert res["method"] == "tobias"
    assert "TOBIAS not installed" in res["reason"]


def test_footprint_tobias_contract_registered():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    key = "aria/scripts/chromatin_footprint_tobias.py"
    assert key in SCRIPT_CONTRACTS
    assert SCRIPT_CONTRACTS[key].validation_level in {"beta", "production"}


def _fake_fragment_pysam(monkeypatch):
    """Capture the real streaming writer without requiring the scientific env."""
    writes = {}

    class Writer:
        def __init__(self, path, *args, **kwargs):
            self.path = path
            writes[path] = []

        def write(self, read):
            writes[self.path].append(read)

        def close(self):
            pass

    def sort(*args):
        Path(args[1]).write_bytes(b"sorted")

    monkeypatch.setitem(sys.modules, "pysam", types.SimpleNamespace(
        AlignmentFile=Writer, AlignedSegment=types.SimpleNamespace,
        qualitystring_to_array=lambda value: value, sort=sort,
        index=lambda path: Path(path + ".bai").write_bytes(b"index")))
    return writes


@pytest.mark.parametrize("separator", ["#", "_", ":"])
@pytest.mark.parametrize("explicit_samples", [False, True])
def test_multisample_fragment_stream_matches_prefixed_barcodes(
    tmp_path, monkeypatch, separator, explicit_samples,
):
    from aria.scripts.chromatin_footprint_tobias import _write_group_bams

    writes = _fake_fragment_pysam(monkeypatch)
    import gzip
    paths = [tmp_path / "sample1.fragments.tsv.gz", tmp_path / "sample2.bed"]
    with gzip.open(paths[0], "wt") as fh:
        fh.write("chr1\t100\t200\tBC1\nchr1\t300\t400\tRAW\n")
    paths[1].write_text("chr1\t500\t600\tBC1\n")
    samples = ["custom1", "custom2"] if explicit_samples else ["sample1", "sample2"]
    barcodes = {f"{samples[0]}{separator}BC1": "A",
                f"{samples[1]}{separator}BC1": "B",
                "RAW": "A", f"{samples[0]}{separator}RAW": "B"}
    result = _write_group_bams(
        [str(p) for p in paths], barcodes, ["A", "B"], {"chr1": 1000}, tmp_path,
        reuse=False, sample_ids=samples if explicit_samples else None)
    assert result["A"]["n_fragments"] == 2  # unprefixed lookup has priority
    assert result["B"]["n_fragments"] == 1
    assert len(writes[str(tmp_path / "A.unsorted.bam")]) == 4
    assert len(writes[str(tmp_path / "B.unsorted.bam")]) == 2


def test_single_fragment_file_uses_explicit_sample_namespace(tmp_path, monkeypatch):
    from aria.scripts.chromatin_footprint_tobias import _write_group_bams

    _fake_fragment_pysam(monkeypatch)
    path = tmp_path / "fragments.tsv"
    path.write_text("chr1\t100\t200\tBC1\n")
    result = _write_group_bams(str(path), {"s#BC1": "A"}, ["A"], {"chr1": 1000},
                               tmp_path, sample_ids=["s"], reuse=False)
    assert result["A"]["n_fragments"] == 1


def _stage_params(tmp_path, mode):
    for name in ["genome.fa", "peaks.bed", "motifs.meme", "a.bam", "b.bam"]:
        (tmp_path / name).write_text("input\n")
    (tmp_path / "genome.fa.fai").write_text("chr1\t1000\n")
    fragments = tmp_path / "s.fragments.tsv"
    fragments.write_text("chr1\t100\t200\tBC1\n")
    groups = tmp_path / "groups.tsv"
    groups.write_text("barcode\tgroup\ns#BC1\tA\ns#BC2\tB\n")
    return {
        "mode": mode, "genome_fasta": str(tmp_path / "genome.fa"),
        "peaks_bed": str(tmp_path / "peaks.bed"),
        "motif_meme": str(tmp_path / "motifs.meme"),
        "group_a": "A", "group_b": "B", "output_dir": str(tmp_path / "out"),
        "fragments_files": [str(fragments)], "barcode_groups": str(groups),
        "sample_ids": ["s"], "make_aggregate_plots": False,
        "condition_bams": {"A": [str(tmp_path / "a.bam")],
                           "B": [str(tmp_path / "b.bam")]},
    }


@pytest.mark.parametrize("mode", ["scatac", "bulk"])
def test_stage_manifest_resume_and_peaks_invalidation(tmp_path, monkeypatch, mode):
    from aria.scripts import chromatin_footprint_tobias as mod

    params = _stage_params(tmp_path, mode)
    monkeypatch.setattr(mod.shutil, "which", lambda _: "/bin/TOBIAS")
    prep_calls, pipeline_calls = [], []

    def prepare(*args, **kwargs):
        prep_calls.append(kwargs)
        if mode == "scatac":
            assert args[0] == params["fragments_files"]
            assert kwargs["sample_ids"] == ["s"]
        else:
            assert kwargs["reference_fasta"] == params["genome_fasta"]
        return {g: {"bam": str(tmp_path / f"{g}.bam"), "n_fragments": 10}
                for g in ["A", "B"]}

    def pipeline(*args, **kwargs):
        pipeline_calls.append(True)
        assert not (Path(params["output_dir"]) /
                    "chromatin_footprint_tobias.stage_manifest.json").exists()
        result = Path(params["output_dir"]) / "bindetect_results.txt"
        result.write_text("name\ttotal_tfbs\tA_B_change\tA_B_pvalue\nTF1\t100\t2\t0.001\n")
        return {"bindetect_results": str(result)}

    monkeypatch.setattr(mod, "_write_group_bams" if mode == "scatac"
                        else "_merge_condition_bams", prepare)
    monkeypatch.setattr(mod, "_tobias_pipeline", pipeline)
    driver = mod.chromatin_footprint_tobias_dispatch
    first = driver(params)
    assert first["ran"] is True
    assert first["validation_level"] == "production"
    assert set(first["volcano_plot"]) == {"png", "svg"}
    manifest = Path(first["stage_manifest_path"])
    fingerprint = json.loads(manifest.read_text())["fingerprint"]
    assert len(fingerprint) == 64
    assert json.loads((manifest.parent / "chromatin_footprint_tobias.result.json").read_text()) == first
    cached = driver(params)
    assert cached == {**first, "reused_from_stage_manifest": True}
    assert len(prep_calls) == len(pipeline_calls) == 1
    Path(params["peaks_bed"]).write_text("chr1\t10\t80\nchanged peaks\n")
    rerun = driver(params)
    assert "reused_from_stage_manifest" not in rerun
    assert len(prep_calls) == len(pipeline_calls) == 2
    assert all(call["reuse"] is False for call in prep_calls)
    assert json.loads(manifest.read_text())["fingerprint"] != fingerprint
    Path(rerun["bindetect_results"]).unlink()
    driver(params)
    assert len(pipeline_calls) == 3
    driver({**params, "resume": False, "make_volcano_plot": False})
    assert len(pipeline_calls) == 4
    assert prep_calls[-1]["reuse"] is False


def test_multisample_driver_validates_every_fragment_path(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    params = _stage_params(tmp_path, "scatac")
    monkeypatch.setattr(mod.shutil, "which", lambda _: "/bin/TOBIAS")
    params["fragments_files"].append(str(tmp_path / "missing.tsv"))
    params["sample_ids"].append("missing")
    result = mod.chromatin_footprint_tobias(params)
    assert result["ran"] is False
    assert result["validation_level"] == "production"
    assert "missing.tsv" in result["reason"]


@pytest.mark.parametrize("sample_ids", [[], ["s", "extra"], [""]])
def test_fragment_sample_ids_require_one_nonempty_id_per_file(tmp_path, monkeypatch, sample_ids):
    from aria.scripts import chromatin_footprint_tobias as mod

    params = _stage_params(tmp_path, "scatac")
    monkeypatch.setattr(mod.shutil, "which", lambda _: "/bin/TOBIAS")
    result = mod.chromatin_footprint_tobias({**params, "sample_ids": sample_ids})
    assert result["ran"] is False
    assert "sample_ids" in result["reason"]


def test_empty_fragment_list_falls_back_to_legacy_file(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _fragment_inputs

    path = tmp_path / "legacy.tsv"
    path.touch()
    assert _fragment_inputs(str(path), []) == ([str(path)], ["legacy"])


@pytest.mark.parametrize("inferential", [True, False])
def test_render_footprint_volcano_exports_all_motifs_with_honest_labels(
    tmp_path, monkeypatch, inferential,
):
    from aria.scripts.chromatin_footprint_tobias import render_footprint_volcano
    from matplotlib.axes import Axes

    title, labels, annotated = [], [], []
    real_title, real_scatter, real_annotate = Axes.set_title, Axes.scatter, Axes.annotate

    def set_title(self, value, *args, **kwargs):
        title.append(value)
        return real_title(self, value, *args, **kwargs)

    def scatter(self, x, y, *args, **kwargs):
        labels.append((kwargs["label"], kwargs["c"], len(x)))
        return real_scatter(self, x, y, *args, **kwargs)

    def annotate(self, value, *args, **kwargs):
        annotated.append(value)
        return real_annotate(self, value, *args, **kwargs)

    monkeypatch.setattr(Axes, "set_title", set_title)
    monkeypatch.setattr(Axes, "scatter", scatter)
    monkeypatch.setattr(Axes, "annotate", annotate)
    table = tmp_path / "motifs.tsv"
    header = "tf\tchange\tpadj\n" if inferential else "name\tA_B_change\tA_B_pvalue\n"
    table.write_text(header + "UP\t2\t0\nDOWN\t-2\t0.001\nNULL\t0.1\t0.8\nBAD\tnan\t0\n")
    summary = {"inference": {"status": "success" if inferential else "descriptive_only"},
               "alpha": 0.05}
    result = render_footprint_volcano(summary, "A", "B", tmp_path,
                                      results_table=str(table), bindetect_results=str(table))
    assert Path(result["png"]).read_bytes().startswith(b"\x89PNG")
    assert "<svg" in Path(result["svg"]).read_text()
    assert all(Path(path).stat().st_size > 1000 for path in result.values())
    assert ("Higher in A", "#c0392b", 1) in labels
    assert ("Higher in B", "#2c7fb8", 1) in labels
    assert sum(item[2] for item in labels) == 3
    assert set(annotated) == {"UP", "DOWN"}
    assert ("replicate inference" if inferential else "descriptive — uncorrected site p") in title[-1]


def test_render_footprint_volcano_honestly_returns_none_without_finite_rows(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import render_footprint_volcano

    table = tmp_path / "empty.tsv"
    table.write_text("tf\tchange\tpadj\nBAD\tnan\t0.05\n")
    assert render_footprint_volcano({"inference": {"status": "success"}}, "A", "B",
                                    tmp_path, results_table=str(table)) is None
    assert render_footprint_volcano({}, "A", "B", tmp_path) is None
    assert not (tmp_path / "footprint_volcano.png").exists()


def test_volcano_falls_back_to_descriptive_when_inference_table_is_missing(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import render_footprint_volcano

    table = tmp_path / "bindetect_results.txt"
    table.write_text("name\tA_B_change\tA_B_pvalue\nTF\t1\t0.01\n")
    result = render_footprint_volcano({"inference": {"status": "success"}}, "A", "B",
                                     tmp_path, results_table=str(tmp_path / "missing.tsv"),
                                     bindetect_results=str(table))
    assert Path(result["svg"]).is_file()
    assert "descriptive" in Path(result["svg"]).read_text()


@pytest.mark.parametrize("key,value", [
    ("group_a", "other"), ("group_b", "other"), ("min_replicates_per_condition", 4),
    ("footprint_fdr", 0.01), ("make_aggregate_plots", True),
    ("replicate_batches", {"r1": "new_batch"}), ("sample_ids", ["other"]),
    ("max_label_permutations", 10), ("permutation_seed", 99),
])
def test_stage_fingerprint_binds_analysis_configuration(tmp_path, key, value):
    from aria.scripts.chromatin_footprint_tobias import _stage_fingerprint

    params = _stage_params(tmp_path, "scatac")
    paths = [params["peaks_bed"]]
    fingerprint = _stage_fingerprint(params, "scatac", paths)
    assert fingerprint == _stage_fingerprint(dict(reversed(list(params.items()))), "scatac", paths)
    assert fingerprint != _stage_fingerprint({**params, key: value}, "scatac", paths)
    assert fingerprint != _stage_fingerprint(params, "bulk", paths)
