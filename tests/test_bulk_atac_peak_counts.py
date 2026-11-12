"""V47 bulk ATAC peak-count matrix guards."""

from pathlib import Path


def test_peak_counts_contract_registered():
    from aria.utils.script_contracts import contract_for_script

    contract = contract_for_script("aria/scripts/chromatin_peak_counts.py")
    assert contract is not None
    assert contract.validation_level == "production"
    assert {f.name for f in contract.success_outputs} >= {
        "counts_matrix_path", "sample_metadata_path", "n_peaks", "n_samples",
    }


def test_peak_count_helpers_parse_bed_and_coverage(tmp_path):
    from aria.scripts.chromatin_peak_counts import (
        _parse_bedtools_counts,
        _read_peak_ids,
        _resolve_sample_ids,
    )

    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\tp1\nchr1\t20\t30\tp2\n")

    assert _read_peak_ids(peaks) == ["chr1:10-20", "chr1:20-30"]
    assert _parse_bedtools_counts(
        "chr1\t10\t20\tp1\t7\nchr1\t20\t30\tp2\t0\n",
        ["chr1:10-20", "chr1:20-30"],
    ) == [7, 0]
    assert _resolve_sample_ids(
        ["/data/A.bam", "/data/B.fragments.tsv.gz"]
    ) == ["A", "B"]


def test_peak_counts_missing_bedtools_is_honest_error(tmp_path, monkeypatch):
    from aria.scripts.chromatin_peak_counts import chromatin_peak_counts

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\n")
    monkeypatch.setattr("shutil.which", lambda name: None)

    res = chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam)],
        "peaks_path": str(peaks),
    })

    assert res["status"] == "error"
    assert res["error_type"] == "MissingDependency"


def test_peak_counts_writes_matrix_and_metadata(tmp_path, monkeypatch):
    import subprocess

    from aria.scripts.chromatin_peak_counts import chromatin_peak_counts

    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\nchr1\t20\t30\n")
    bam_a = tmp_path / "A.bam"
    bam_b = tmp_path / "B.bam"
    bam_a.write_bytes(b"")
    bam_b.write_bytes(b"")

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/bedtools")

    def fake_run(cmd, capture_output, text, check):
        assert cmd[:4] == ["/usr/bin/bedtools", "coverage", "-counts", "-a"]
        stdout = (
            "chr1\t10\t20\t5\nchr1\t20\t30\t1\n"
            if cmd[-1] == str(bam_a)
            else "chr1\t10\t20\t2\nchr1\t20\t30\t9\n"
        )
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    res = chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam_a), str(bam_b)],
        "sample_ids": ["ctrl_1", "treated_1"],
        "sample_metadata": {
            "ctrl_1": {"condition": "ctrl", "replicate": "r1"},
            "treated_1": {"condition": "treated", "replicate": "r1"},
        },
        "peaks_path": str(peaks),
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert res["validation_level"] == "production"
    assert res["n_peaks"] == 2
    assert res["n_samples"] == 2
    matrix = Path(res["counts_matrix_path"]).read_text().splitlines()
    assert matrix == [
        "peak_id\tctrl_1\ttreated_1",
        "chr1:10-20\t5\t2",
        "chr1:20-30\t1\t9",
    ]
    metadata = Path(res["sample_metadata_path"]).read_text()
    assert "condition" in metadata
    assert "replicate" in metadata
    # Empty stub BAMs have no readable header -> honest fallback warning.
    assert any("genome" in w for w in res["warnings"])


def test_peak_counts_uses_sorted_low_memory_coverage(tmp_path, monkeypatch):
    """Large bulk ATAC BAMs OOM `bedtools coverage` unless the low-memory
    sweeping algorithm (`-sorted -g`) is used. When a genome can be read from
    the alignment header, counting must engage `-sorted` and emit peaks in
    genome order. Counts are identical to the non-sorted call."""
    import subprocess

    from aria.scripts import chromatin_peak_counts as mod

    peaks = tmp_path / "peaks.bed"
    # Deliberately NOT in genome order, to prove we sort before -sorted counting.
    peaks.write_text("chr1\t20\t30\nchr1\t10\t20\n")
    bam = tmp_path / "A.bam"
    bam.write_bytes(b"")

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/bedtools")
    monkeypatch.setattr(mod, "_genome_sizes_from_bam", lambda f: [("chr1", 1000)])

    seen = {}

    def fake_run(cmd, capture_output, text, check):
        seen["cmd"] = cmd
        # bedtools emits one row per -a feature, in the (sorted) -a order.
        return subprocess.CompletedProcess(
            cmd, 0, stdout="chr1\t10\t20\t3\nchr1\t20\t30\t8\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    res = mod.chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam)],
        "sample_ids": ["s1"],
        "peaks_path": str(peaks),
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert "-sorted" in seen["cmd"] and "-g" in seen["cmd"]
    assert "sorted" in res["counting_method"]
    matrix = Path(res["counts_matrix_path"]).read_text().splitlines()
    assert matrix == ["peak_id\ts1", "chr1:10-20\t3", "chr1:20-30\t8"]


def test_peak_counts_falls_back_to_unsorted_without_genome(tmp_path, monkeypatch):
    """When no genome can be read from the header, counting falls back to the
    plain (non-sorted) call rather than failing."""
    import subprocess

    from aria.scripts import chromatin_peak_counts as mod

    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10\t20\n")
    bam = tmp_path / "A.bam"
    bam.write_bytes(b"")

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/bedtools")
    monkeypatch.setattr(mod, "_genome_sizes_from_bam", lambda f: None)

    def fake_run(cmd, capture_output, text, check):
        assert cmd[:4] == ["/usr/bin/bedtools", "coverage", "-counts", "-a"]
        return subprocess.CompletedProcess(cmd, 0, stdout="chr1\t10\t20\t4\n",
                                           stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    res = mod.chromatin_peak_counts({
        "data_type": "bulk_ATAC",
        "files": [str(bam)],
        "sample_ids": ["s1"],
        "peaks_path": str(peaks),
        "output_dir": str(tmp_path / "out"),
    })

    assert res["status"] == "success"
    assert res["counting_method"] == "bedtools coverage -counts"


def test_run_ledger_has_peak_count_matrix_node():
    from aria.agents.narrative.run_ledger import build_run_ledger

    ledger = build_run_ledger(
        {"design_intelligence": {"recommended": ["build peak count matrix"]}},
        {"chromatin_agent": {"status": "done", "findings": {
            "bulk_ATAC": {"status": "done", "findings": {
                "peak_counts": {"status": "success"},
            }}
        }}},
    )

    nodes = {entry["node_id"]: entry for entry in ledger["entries"]}
    assert nodes["ledger://chromatin/peak_count_matrix"]["planned"] is True
    assert nodes["ledger://chromatin/peak_count_matrix"]["status"] == "ran"


# Scientific fixtures execute in the chromatin stack even when pytest itself is
# hosted in aria-env, which deliberately does not include pysam.
_BAM_HELPERS = r'''
import json
from pathlib import Path
from unittest.mock import patch
import pysam

def make_bam(path, specs, marked=False):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chrM", "LN": 1000}]}
    if marked:
        header["PG"] = [{"ID": "markdup", "PN": "samtools", "CL": "samtools markdup"}]
    records = []
    for i, spec in enumerate(specs):
        read = pysam.AlignedSegment()
        read.query_name = spec.get("name", f"read{i}")
        read.query_sequence = "A" * 30
        read.query_qualities = pysam.qualitystring_to_array("I" * 30)
        read.flag = spec.get("flag", 0)
        read.reference_id = spec.get("tid", 0)
        read.reference_start = spec.get("start", 100)
        read.mapping_quality = 60
        read.cigarstring = "30M" if not read.is_unmapped else None
        if read.is_paired:
            read.next_reference_id = read.reference_id
            read.next_reference_start = read.reference_start + 50
        read.template_length = spec.get("tlen", 0)
        records.append(read)
    records.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in records:
            bam.write(read)
    return str(path)
'''


def _chromatin_python(tmp_path, code):
    import importlib.util
    import subprocess
    import sys

    source = _BAM_HELPERS + "\ntmp_path = Path(" + repr(str(tmp_path)) + ")\n" + code
    if importlib.util.find_spec("pysam") is not None:
        command = [sys.executable, "-c", source]
    else:
        command = ["conda", "run", "-n", "aria-chromatin-env", "python", "-c", source]
    proc = subprocess.run(command, capture_output=True, text=True, timeout=180)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_cpm_cache_hits_and_content_parameter_artifact_invalidation(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts.chromatin_peak_counts import chromatin_peak_counts
bam = make_bam(tmp_path / "a.bam", [{"start": 100}, {"start": 110}, {"start": 300}])
peaks = tmp_path / "p.bed"
peaks.write_text("chr1\t90\t150\nchr1\t290\t350\nchr1\t800\t900\n")
params = {"files": [bam], "peaks_path": str(peaks), "output_dir": str(tmp_path / "out"),
          "sample_ids": ["a"], "sample_metadata": {"a": {"condition": "ctrl"}}}
first = chromatin_peak_counts(params)
assert first["status"] == "success", first
assert first["validation_level"] == "production"
assert not first["reused_from_stage_manifest"]
assert Path(first["counts_cpm_path"]).name == "bulk_atac_peak_cpm.tsv"
rows = Path(first["counts_cpm_path"]).read_text().splitlines()
values = [float(row.split("\t")[1]) for row in rows[1:]]
assert abs(sum(values) - 1e6) < .01 and values[2] == 0
assert Path(first["counts_matrix_path"]).read_text().splitlines()[1:] == [
    "chr1:90-150\t2", "chr1:290-350\t1", "chr1:800-900\t0"]
with patch("aria.scripts.chromatin_peak_counts._count_one_sample", side_effect=AssertionError("cache miss")):
    assert chromatin_peak_counts(params)["reused_from_stage_manifest"]
for field in ("counts_cpm_path", "counts_matrix_path", "sample_metadata_path"):
    Path(first[field]).unlink()
    assert not chromatin_peak_counts(params)["reused_from_stage_manifest"]
params["sample_metadata"]["a"]["condition"] = "treated"
assert not chromatin_peak_counts(params)["reused_from_stage_manifest"]
assert "treated" in Path(first["sample_metadata_path"]).read_text()
peaks.write_text("chr1\t90\t150\n")
assert not chromatin_peak_counts(params)["reused_from_stage_manifest"]
assert chromatin_peak_counts(params)["reused_from_stage_manifest"]
make_bam(Path(bam), [{"start": 500}])
last = chromatin_peak_counts(params)
assert not last["reused_from_stage_manifest"]
assert float(Path(last["counts_cpm_path"]).read_text().splitlines()[1].split("\t")[1]) == 0
assert chromatin_peak_counts(params)["reused_from_stage_manifest"]
params["sample_ids"] = ["renamed"]
assert not chromatin_peak_counts(params)["reused_from_stage_manifest"]
assert Path(last["counts_matrix_path"]).read_text().splitlines()[0] == "peak_id\trenamed"
manifest = next((tmp_path / "out").glob("*.stage_manifest.json"))
manifest.write_text("not JSON")
assert not chromatin_peak_counts(params)["reused_from_stage_manifest"]
''')


def test_local_bam_qc_all_samples_duplicates_banding_frip_and_tss(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts.chromatin_qc import chromatin_qc, _estimate_dup_rate
specs = [{"start": 2996, "flag": 99, "tlen": 80 if i % 2 else 200} for i in range(20)]
specs += [{"start": 3200, "flag": 99, "tlen": 350},
          {"start": 1001}, {"start": 4991}, {"start": 50, "tid": 1},
          {"start": 2996, "flag": 1024}, {"start": 2996, "flag": 256},
          {"start": 2996, "flag": 2048}, {"start": 6000, "flag": 4}]
files = [make_bam(tmp_path / f"s{i}.bam", specs, marked=True) for i in range(9)]
peaks = tmp_path / "peaks.bed"
peaks.write_text("chr1\t2900\t3400\n")
gtf = tmp_path / "genes.gtf"
gtf.write_text('chr1\ttest\tgene\t3001\t4000\t.\t+\t.\tgene_id "g";\n'
               'chr1\ttest\ttranscript\t3001\t4000\t.\t+\t.\tgene_id "g";\n')
params = {"data_type": "bulk_ATAC", "files": files, "peaks_bed": str(peaks), "gtf": str(gtf)}
# No indexes exist: full counts must still include every primary alignment.
result = chromatin_qc(params)
assert result["status"] == "success", result
assert result["validation_level"] == "production"
assert result["n_samples"] == 9 and result["mapped_reads"] == 25 * 9
assert abs(result["dup_rate"] - 1 / 26) < 1e-10
assert len(result["dup_rate_per_sample"]) == 9
assert _estimate_dup_rate(files[0]) == 1 / 26
assert result["mito_fraction"] == .04
assert abs(result["frip"] - 21 / 24) < 1e-10
assert len(result["frip_per_sample"]) == 9 and len(result["frip_distribution"]) == 9
assert abs(result["tss_enrichment"] - (20 / 101) / (2 / 200)) < 1e-10
assert result["qc_complete"] and result["pass_qc"]
assert result["metrics_not_computed"] == []
f = result["fragment_sizes"]
assert f["n_fragments_sampled"] == 21 * 9
assert f["nucleosome_banding_detected"]
assert f["nfr_fraction"] == round(10 / 21, 4)
assert f["mononucleosomal_fraction"] == round(10 / 21, 4)
assert f["dinucleosomal_fraction"] == round(1 / 21, 4)
assert f["size_histogram"]["bin_edges"] == list(map(float, range(0, 1001, 10)))
assert len(f["size_histogram"]["counts"]) == 100
assert sum(f["size_histogram"]["counts"]) == 21 * 9
# The cap limits detail, never the exact read/mitochondrial/duplicate counters.
capped = chromatin_qc({**params, "max_reads_per_bam": 2})
assert capped["mapped_reads"] == result["mapped_reads"]
assert capped["dup_rate"] == result["dup_rate"]
assert capped["fragment_sizes"]["n_fragments_sampled"] <= 18
incomplete = chromatin_qc({"data_type": "bulk_ATAC", "files": [files[0]]})
assert incomplete["frip"] is None and incomplete["tss_enrichment"] is None
assert incomplete["metrics_not_computed"] == ["tss_enrichment", "frip"]
assert incomplete["pass_qc"] is None and not incomplete["qc_complete"]
failed = chromatin_qc({**params, "thresholds": {"min_frip": .99, "min_tss": 6}})
assert failed["qc_complete"] and failed["pass_qc"] is False
''')


def test_single_end_qc_missing_metrics_and_zero_dup_are_not_fabricated(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts.chromatin_qc import chromatin_qc, _estimate_dup_rate
bam = make_bam(tmp_path / "se.bam", [{"start": 100}])
assert _estimate_dup_rate(bam) is None
marked = make_bam(tmp_path / "marked.bam", [{"start": 100}], marked=True)
assert _estimate_dup_rate(marked) == 0.0
result = chromatin_qc({"files": [bam], "data_type": "bulk_ATAC"})
assert result["dup_rate"] is None and result.get("fragment_sizes") is None
assert result["frip"] is None and result["pass_qc"] is None
peaks = tmp_path / "p.bed"
peaks.write_text("chr1\t50\t200\n")
partial = chromatin_qc({"files": [bam], "data_type": "bulk_ATAC", "peaks_path": str(peaks)})
assert partial["frip"] == 1.0
assert partial["metrics_not_computed"] == ["tss_enrichment"] and partial["pass_qc"] is None
gtf = tmp_path / "g.gtf"
gtf.write_text('chr1\ttest\tgene\t8001\t9000\t.\t+\t.\tgene_id "g";\n')
empty = chromatin_qc({"files": [bam], "data_type": "bulk_ATAC", "gtf": str(gtf)})
assert empty["tss_enrichment"] is None
''')


def test_frip_without_samtools_bedtools_and_nested_interval_fallback(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts.chromatin_peaks import _compute_frip
bam = make_bam(tmp_path / "a.bam", [
    {"start": 100}, {"start": 400}, {"start": 900},
    {"start": 100, "flag": 1024}, {"start": 100, "flag": 256},
    {"start": 100, "flag": 2048}, {"start": 100, "flag": 4}])
peaks = tmp_path / "nested.bed"
peaks.write_text("chr1\t90\t800\nchr1\t95\t140\n")
assert abs(_compute_frip(bam, str(peaks)) - 2 / 3) < 1e-10
with patch("shutil.which", return_value=None), patch("subprocess.run", side_effect=AssertionError("external tool called")):
    assert abs(_compute_frip(bam, str(peaks)) - 2 / 3) < 1e-10
empty = tmp_path / "empty.bed"
empty.write_text("")
assert _compute_frip(bam, str(empty)) is None
''')


def test_peak_calling_se_detection_pe_preservation_blacklist_and_stage_cache(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts import chromatin_peaks as cp
se = make_bam(tmp_path / "se.bam", [{"start": 100}, {"start": 400}])
pe = make_bam(tmp_path / "pe.bam", [{"start": 100, "flag": 99, "tlen": 200}])
blacklist = tmp_path / "blacklist.bed"
blacklist.write_text("chr1\t390\t500\n")
calls = []
def fake_macs(cmd, timeout):
    calls.append(cmd)
    outdir = Path(cmd[cmd.index("--outdir") + 1])
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / (cmd[cmd.index("-n") + 1] + "_peaks.narrowPeak")).write_text(
        "chr1\t90\t150\tp1\nchr1\t390\t450\tp2\n")
    return {"status": "success"}
params = {"data_type": "bulk_ATAC", "files": [se], "output_dir": str(tmp_path / "out"),
          "blacklist_bed": str(blacklist), "macs3_params": {"format": "BAMPE", "qvalue": .01}}
with patch.object(cp, "_run_macs3", fake_macs):
    first = cp.chromatin_peaks(params)
    assert first["status"] == "success", first
    assert first["effective_format"] == "BAM_SE_TN5_SHIFT"
    assert first["n_peaks"] == 1 and first["n_blacklist_filtered"] == 1
    assert first["validation_level"] == "production"
    assert first["frip_per_sample"] == {"se": .5} and first["frip"] == .5
    assert any("Single-end" in w for w in first["warnings"])
    cmd = calls[-1]
    assert cmd[cmd.index("-f") + 1] == "BAM"
    assert "--nomodel" in cmd and cmd[cmd.index("--shift") + 1] == "-100"
    assert cmd[cmd.index("--extsize") + 1] == "200"
    assert cmd[cmd.index("-q") + 1] == "0.01"
    calls.clear()
    assert cp.chromatin_peaks(params)["reused_from_stage_manifest"] and not calls
    Path(first["peaks_path"]).unlink()
    assert not cp.chromatin_peaks(params)["reused_from_stage_manifest"]
    blacklist.write_text("chr1\t800\t900\n")
    assert not cp.chromatin_peaks(params)["reused_from_stage_manifest"]
    params["macs3_params"]["qvalue"] = .02
    assert not cp.chromatin_peaks(params)["reused_from_stage_manifest"]
    paired = cp.chromatin_peaks({**params, "files": [pe], "output_dir": str(tmp_path / "paired")})
    assert paired["effective_format"] == "BAMPE"
    assert "--shift" not in calls[-1]
    cohort = cp.chromatin_peaks({**params, "files": [se, pe], "output_dir": str(tmp_path / "cohort")})
    assert cohort["n_consensus_peaks"] == 2
    assert set(cohort["frip_per_sample"]) == {"se", "pe"}
    assert cp.chromatin_peaks({**params, "files": [se, pe], "output_dir": str(tmp_path / "cohort")})["reused_from_stage_manifest"]
    Path(cohort["consensus_peaks_path"]).unlink()
    assert not cp.chromatin_peaks({**params, "files": [se, pe], "output_dir": str(tmp_path / "cohort")})["reused_from_stage_manifest"]
''')


def test_cram_local_reference_qc_frip_and_peak_counts(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts.chromatin_peak_counts import chromatin_peak_counts
from aria.scripts.chromatin_qc import chromatin_qc
from aria.scripts.chromatin_peaks import _compute_frip
ref = tmp_path / "ref.fa"
ref.write_text(">chr1\n" + "A" * 10000 + "\n>chrM\n" + "A" * 1000 + "\n")
pysam.faidx(str(ref))
bam = make_bam(tmp_path / "a.bam", [{"start": 100}, {"start": 400}])
cram = tmp_path / "a.cram"
with pysam.AlignmentFile(bam, "rb") as source:
    with pysam.AlignmentFile(str(cram), "wc", template=source, reference_filename=str(ref)) as dest:
        for read in source:
            dest.write(read)
peaks = tmp_path / "p.bed"
peaks.write_text("chr1\t90\t150\nchr1\t390\t450\n")
assert _compute_frip(str(cram), str(peaks), reference_fasta=str(ref)) == 1.0
qc = chromatin_qc({"data_type": "bulk_ATAC", "files": [str(cram)],
                   "reference_fasta": str(ref), "peaks_bed": str(peaks)})
assert qc["mapped_reads"] == 2 and qc["frip"] == 1.0
params = {"files": [str(cram), bam], "sample_ids": ["cram", "bam"],
          "reference_fasta": str(ref), "peaks_path": str(peaks), "output_dir": str(tmp_path / "out")}
counts = chromatin_peak_counts(params)
assert counts["status"] == "success", counts
assert Path(counts["counts_matrix_path"]).read_text().splitlines() == [
    "peak_id\tcram\tbam", "chr1:90-150\t1\t1", "chr1:390-450\t1\t1"]
assert chromatin_peak_counts(params)["reused_from_stage_manifest"]
''')


def test_tss_reverse_cuts_strand_boundaries_and_unmeasured_background(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts.chromatin_qc import _read_local_tss, _scan_bam_metrics
# Minus-strand TSS is end - 1 = 4000, not gene start; annotation duplicates
# must not double-weight the same promoter.
gtf = tmp_path / "minus.gtf"
gtf.write_text('chr1\ttest\tgene\t3001\t4001\t.\t-\t.\tgene_id "g";\n'
               'chr1\ttest\ttranscript\t3001\t4001\t.\t-\t.\tgene_id "g";\n')
tss = _read_local_tss(gtf)
assert tss["chr1"][0] == [(4000, -1)]
# Reverse cut = reference_end - 5 = start + 25, with a 30M CIGAR.
# Two central cuts, two flank cuts, and cuts exactly at the excluded inner
# flank boundaries (+/-1900) ensure correct strand-aware window endpoints.
bam = make_bam(tmp_path / "minus.bam", [
    {"start": 3975, "flag": 16}, {"start": 3925, "flag": 16},
    {"start": 1975, "flag": 16}, {"start": 5975, "flag": 16},
    {"start": 2075, "flag": 16}, {"start": 5875, "flag": 16}])
metrics, _ = _scan_bam_metrics(bam, tss=tss)
assert abs(metrics["tss_enrichment"] - (2 / 101) / (2 / 200)) < 1e-10
center_only = make_bam(tmp_path / "center.bam", [{"start": 3975, "flag": 16}])
metrics, _ = _scan_bam_metrics(center_only, tss=tss)
assert metrics["tss_enrichment"] is None  # no fabricated denominator
''')


def test_failed_peak_rerun_cannot_reuse_partial_artifact(tmp_path):
    _chromatin_python(tmp_path, r'''
from aria.scripts import chromatin_peaks as cp
bam = make_bam(tmp_path / "a.bam", [{"start": 100}])
params = {"files": [bam], "data_type": "bulk_ATAC", "output_dir": str(tmp_path / "out")}
def successful(cmd, timeout):
    (tmp_path / "out" / "a_peaks.narrowPeak").write_text("chr1\t90\t150\n")
    return {"status": "success"}
with patch.object(cp, "_run_macs3", successful):
    result = cp.chromatin_peaks(params)
    assert result["status"] == "success"
Path(result["peaks_path"]).unlink()
def partial(cmd, timeout):
    Path(result["peaks_path"]).write_text("chr1\t0\t5\n")
    return {"status": "error", "error_type": "MACS3NonZero", "details": "failed"}
with patch.object(cp, "_run_macs3", partial):
    assert cp.chromatin_peaks(params)["status"] == "error"
with patch.object(cp, "_run_macs3", side_effect=successful) as call:
    rerun = cp.chromatin_peaks(params)
    assert rerun["status"] == "success" and not rerun["reused_from_stage_manifest"]
    assert call.call_count == 1
''')
