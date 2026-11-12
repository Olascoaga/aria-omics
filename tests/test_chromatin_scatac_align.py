"""B2b slice 2: scATAC FASTQ -> fragments front-end (chromap).

chromatin_scatac_align is the single-cell analogue of atac_align (bulk): chromap
FASTQ+barcode -> fragments.tsv.gz. It is binary/input-gated — missing inputs or a
missing aligner return a structured not-run, never a fabricated fragments file.
The agent only runs it when the barcode read + whitelist are explicit (the
barcode read cannot be inferred), then feeds the fragments into the bridge.
"""

from pathlib import Path

from aria.scripts import chromatin_scatac_align as sa
from aria.agents.chromatin_agent import _pick_read


# ── pure-input validation ─────────────────────────────────────────────────

def test_missing_inputs(tmp_path):
    out = sa.chromatin_scatac_align({"r1_fastq": str(tmp_path / "x.fq.gz")})
    assert out["status"] == "skipped"
    assert out["reason"] == "missing_input"


def test_missing_whitelist(tmp_path):
    for n in ("r1.fq.gz", "r3.fq.gz", "bc.fq.gz", "g.fa"):
        (tmp_path / n).write_text("x")
    out = sa.chromatin_scatac_align({
        "r1_fastq": str(tmp_path / "r1.fq.gz"),
        "r3_fastq": str(tmp_path / "r3.fq.gz"),
        "barcode_fastq": str(tmp_path / "bc.fq.gz"),
        "genome_fasta": str(tmp_path / "g.fa"),
    })
    assert out["status"] == "skipped"
    assert out["reason"] == "missing_barcode_whitelist"


def test_missing_genome(tmp_path):
    for n in ("r1.fq.gz", "r3.fq.gz", "bc.fq.gz", "wl.txt"):
        (tmp_path / n).write_text("x")
    out = sa.chromatin_scatac_align({
        "r1_fastq": str(tmp_path / "r1.fq.gz"),
        "r3_fastq": str(tmp_path / "r3.fq.gz"),
        "barcode_fastq": str(tmp_path / "bc.fq.gz"),
        "genome_fasta": str(tmp_path / "nope.fa"),
        "barcode_whitelist": str(tmp_path / "wl.txt"),
    })
    assert out["status"] == "skipped"
    assert out["reason"] == "missing_genome_fasta"


def test_chromap_unavailable(tmp_path, monkeypatch):
    for n in ("r1.fq.gz", "r3.fq.gz", "bc.fq.gz", "g.fa", "wl.txt"):
        (tmp_path / n).write_text("x")
    import subprocess

    def boom(*a, **k):
        raise FileNotFoundError("chromap")

    monkeypatch.setattr(subprocess, "run", boom)
    out = sa.chromatin_scatac_align({
        "r1_fastq": str(tmp_path / "r1.fq.gz"),
        "r3_fastq": str(tmp_path / "r3.fq.gz"),
        "barcode_fastq": str(tmp_path / "bc.fq.gz"),
        "genome_fasta": str(tmp_path / "g.fa"),
        "barcode_whitelist": str(tmp_path / "wl.txt"),
        "output_dir": str(tmp_path / "out"),
    })
    assert out["status"] == "skipped"
    assert out["reason"] == "chromap_unavailable"


def test_pick_read():
    fqs = ["/d/s_R1_001.fastq.gz", "/d/s_R2_001.fastq.gz", "/d/s_R3_001.fastq.gz"]
    assert _pick_read(fqs, ("_r1", "_R1")) == "/d/s_R1_001.fastq.gz"
    assert _pick_read(fqs, ("_r3", "_R3")) == "/d/s_R3_001.fastq.gz"
    assert _pick_read(fqs, ("_r9",)) is None


def test_contract_registered():
    from aria.utils.script_contracts import contract_for_script
    c = contract_for_script("aria/scripts/chromatin_scatac_align.py")
    assert c is not None and c.validation_level in ("beta", "production")


# ── agent wiring: FASTQ without barcode metadata skips honestly ────────────

class _Env:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, *, stack, script_path, params):
        self.calls.append(Path(script_path).name)
        return {"status": "skipped", "ran": False, "reason": "should_not_run"}


def _agent(env):
    from aria.agents.chromatin_agent import ChromatinAgent
    a = ChromatinAgent.__new__(ChromatinAgent)
    a.env = env
    a.publish_status = lambda *x, **k: None
    a.publish_finding = lambda *x, **k: None
    a.publish_escalation = lambda *x, **k: None
    a._publish_qc_finding = lambda *x, **k: None
    a._publish_peaks_finding = lambda *x, **k: None
    return a


def test_scatac_fastq_without_barcode_skips(tmp_path):
    fq = tmp_path / "s_R1_001.fastq.gz"
    fq.write_bytes(b"")
    env = _Env()
    agent = _agent(env)
    res = agent._run_scatac("exp", {"genome": "hg38"}, {}, [str(fq)])
    # No barcode read/whitelist -> honest skip, nothing dispatched.
    assert res["status"] == "done"
    assert res["findings"]["fastq_to_fragments"]["reason"] == "missing_barcode_inputs"
    assert env.calls == []


def test_scatac_fastq_aligns_then_continues(tmp_path, monkeypatch):
    r1 = tmp_path / "s_R1_001.fastq.gz"; r1.write_bytes(b"")
    r3 = tmp_path / "s_R3_001.fastq.gz"; r3.write_bytes(b"")
    bc = tmp_path / "s_R2_001.fastq.gz"; bc.write_bytes(b"")
    wl = tmp_path / "wl.txt"; wl.write_text("AAAA")

    class _OkEnv:
        def __init__(self):
            self.calls = []

        def run_in_stack(self, *, stack, script_path, params):
            name = Path(script_path).name
            self.calls.append((stack, name))
            if name == "chromatin_scatac_align.py":
                return {"status": "success", "ran": True, "aligner": "chromap",
                        "fragments_file": str(tmp_path / "fragments.tsv.gz")}
            if name == "chromatin_qc.py":
                return {"status": "success", "data_type": "scATAC", "warnings": []}
            if name == "chromatin_peaks.py":
                return {"status": "success", "data_type": "scATAC", "n_peaks": 5,
                        "peaks_path": "/tmp/p.narrowPeak"}
            if name == "chromatin_fragments_to_matrix.py":
                return {"status": "skipped", "ran": False, "reason": "snapatac2_unavailable"}
            if name == "chromatin_motifs.py":
                return {"status": "success"}
            return {"status": "success"}

    env = _OkEnv()
    agent = _agent(env)
    res = agent._run_scatac(
        "exp",
        {"genome": "hg38", "genome_fasta": str(tmp_path / "g.fa"),
         "barcode_fastq": str(bc), "barcode_whitelist": str(wl)},
        {},
        [str(r1), str(r3), str(bc)],
    )
    assert res["status"] == "done"
    assert res["findings"]["fastq_to_fragments"]["status"] == "success"
    names = [n for _s, n in env.calls]
    # chromap ran first (atacseq stack), then the fragments path (QC/peaks).
    assert names[0] == "chromatin_scatac_align.py"
    assert ("atacseq", "chromatin_scatac_align.py") in env.calls
    assert "chromatin_qc.py" in names and "chromatin_peaks.py" in names


import gzip
import pytest
from aria.utils.scatac_fastq_manifest import (
    resolve_scatac_fastq_manifest, manifest_fastq_files, manifest_library_types,
    first_manifest_fastq,
)


def _manifest_row(tmp_path, platform="10x_atac", lanes=1):
    row = {"library_id": "lib", "sample_id": "sample", "donor_id": "donor",
           "platform": platform, "fastqs": {}}
    for role in ("R1", "R2", "R3"):
        paths = []
        for lane in range(lanes):
            name = f"{role}_L{lane + 1:03d}.fq"
            (tmp_path / name).write_text("x")
            paths.append(name)
        row["fastqs"][role] = paths
    (tmp_path / "whitelist.txt").write_text("AACG\n")
    row["barcode_whitelist"] = "whitelist.txt"
    return row


def _validate_row(tmp_path, row, **defaults):
    return resolve_scatac_fastq_manifest(
        {"data_dir": str(tmp_path), "scatac_fastq_manifest": {
            "schema_version": "1", "libraries": [row], **defaults}}, require_paths=True)


def test_multilane_manifest_and_stable_flattening(tmp_path):
    row = _manifest_row(tmp_path, lanes=2)
    result = _validate_row(tmp_path, row)
    assert result["status"] == "valid", result
    canonical = result["manifest"]["libraries"][0]
    expected = [str(tmp_path / path) for role in ("R1", "R2", "R3") for path in row["fastqs"][role]]
    assert manifest_fastq_files(result["manifest"]) == expected
    assert list(manifest_library_types(result["manifest"])) == expected
    assert first_manifest_fastq(result["manifest"]) == expected[0]
    assert canonical["fastq_lanes"] == canonical["fastqs"]
    assert first_manifest_fastq({}) is None


def test_single_lane_manifest_preserves_v1_shape(tmp_path):
    row = _manifest_row(tmp_path)
    row["fastqs"]["R2"] = row["fastqs"]["R2"][0]
    result = _validate_row(tmp_path, row)
    canonical = result["manifest"]["libraries"][0]
    assert all(isinstance(value, str) for value in canonical["fastqs"].values())
    assert all(isinstance(value, list) for value in canonical["fastq_lanes"].values())


@pytest.mark.parametrize("bad", [[], [""], [123], [None]])
def test_manifest_rejects_invalid_lane_values(tmp_path, bad):
    row = _manifest_row(tmp_path)
    row["fastqs"]["R1"] = bad
    assert _validate_row(tmp_path, row)["status"] == "invalid"


def test_manifest_duplicate_lane_rejected(tmp_path):
    row = _manifest_row(tmp_path, lanes=2)
    row["fastqs"]["R2"][1] = row["fastqs"]["R1"][0]
    result = _validate_row(tmp_path, row)
    assert result["status"] == "invalid"
    assert "assigned more than once" in " ".join(result["errors"])


@pytest.mark.parametrize("platform", ["10x_atac", "10x_multiome_atac", "sci_atac", "share_seq", "surecell", "custom"])
def test_manifest_platform_presets(tmp_path, platform):
    row = _manifest_row(tmp_path, platform=platform)
    if platform in {"sci_atac", "surecell", "custom"}:
        row["fastqs"].pop("R3")
    if platform == "custom":
        row["fastqs"].pop("R2")
        row["read_format"] = "r1:0:-1,bc:0:15"
    if platform not in {"10x_atac", "10x_multiome_atac"}:
        row.pop("barcode_whitelist")
    result = _validate_row(tmp_path, row)
    assert result["status"] == "valid", result
    assert result["manifest"]["libraries"][0]["platform"] == platform


@pytest.mark.parametrize("change", ["platform", "orientation", "custom_format", "missing_R3"])
def test_manifest_rejects_platform_contract_violations(tmp_path, change):
    row = _manifest_row(tmp_path)
    if change == "platform":
        row["platform"] = "unknown"
    elif change == "orientation":
        row["barcode_orientation"] = "sideways"
    elif change == "custom_format":
        row["platform"] = "custom"
    else:
        row["fastqs"].pop("R3")
    assert _validate_row(tmp_path, row)["status"] == "invalid"


def test_manifest_defaults_and_orientation_alias(tmp_path):
    row = _manifest_row(tmp_path)
    row.pop("platform")
    row["fastqs"].pop("R3")
    result = _validate_row(tmp_path, row, platform="sci_atac", barcode_orientation="rc")
    assert result["status"] == "valid"
    canonical = result["manifest"]["libraries"][0]
    assert canonical["barcode_orientation"] == "reverse_complement"
    row["barcode_orientation"] = "fwd"
    assert _validate_row(tmp_path, row, platform="sci_atac")["manifest"]["libraries"][0]["barcode_orientation"] == "direct"


def _write_fastq(path, sequences):
    text = "".join(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(sequences))
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return str(path)


@pytest.mark.parametrize("sequence,orientation", [("AACG", "direct"), ("CGTT", "reverse_complement")])
def test_detect_whitelist_orientation(tmp_path, sequence, orientation):
    whitelist = tmp_path / "whitelist.txt.gz"
    with gzip.open(whitelist, "wt") as handle:
        handle.write("AACG\ninvalid\n")
    paths = [_write_fastq(tmp_path / "bc1.fq", [sequence] * 12),
             _write_fastq(tmp_path / "bc2.fq.gz", [sequence] * 12)]
    result = sa.detect_barcode_whitelist_orientation(paths, str(whitelist))
    assert result["orientation"] == orientation
    assert result["n_sampled"] == 24
    assert result["direct_rate"] == (1.0 if orientation == "direct" else 0.0)
    assert result["rc_rate"] == (1.0 if orientation == "reverse_complement" else 0.0)


def test_orientation_sampling_limit_slice_and_fallback(tmp_path):
    whitelist = tmp_path / "wl.txt"
    whitelist.write_text("AACG\n")
    path = _write_fastq(tmp_path / "bc.fq", ["GGAACGTT"] * 30)
    result = sa.detect_barcode_whitelist_orientation([path], str(whitelist), max_reads=7, bc_slice=(2, 6))
    assert result["n_sampled"] == 7 and result["direct_rate"] == 1
    Path(path).write_text("x")
    assert sa.detect_barcode_whitelist_orientation([path], str(whitelist)) == {
        "orientation": "direct", "n_sampled": 0, "direct_rate": None,
        "rc_rate": None, "status": "insufficient_reads_fallback"}
    whitelist.write_text("invalid\n")
    assert sa.detect_barcode_whitelist_orientation([], str(whitelist))["n_sampled"] == 0


def _alignment_params(tmp_path, barcode="AACG", lanes=1):
    params = {"output_dir": str(tmp_path / "out")}
    for key in ("r1_fastq", "r3_fastq", "barcode_fastq"):
        values = [_write_fastq(tmp_path / f"{key}, lane{i}.fq", [barcode] * 24) for i in range(lanes)]
        params[key] = values if lanes > 1 else values[0]
    (tmp_path / "g.fa").write_text(">contig\nACGT\n")
    (tmp_path / "wl.txt").write_text("AACG\n")
    params.update(genome_fasta=str(tmp_path / "g.fa"), barcode_whitelist=str(tmp_path / "wl.txt"))
    return params


def _capture_alignment(monkeypatch):
    calls = []

    def run(cmd, warnings, step):
        calls.append((step, cmd))
        if step == "align":
            Path(cmd[cmd.index("-o") + 1]).write_text("contig\t1\t3\tAACG\t1\n")
        return "ok"

    def compress(plain, gz, warnings):
        with gzip.open(gz, "wt") as handle:
            handle.write(plain.read_text())
        return "ok"

    monkeypatch.setattr(sa, "_run", run)
    monkeypatch.setattr(sa, "_bgzip_tabix", compress)
    monkeypatch.setattr(sa, "_chromap_version", lambda: "unit-test-tool")
    return calls


def test_auto_flip_writes_rc_whitelist(tmp_path, monkeypatch):
    params = _alignment_params(tmp_path, barcode="CGTT")
    calls = _capture_alignment(monkeypatch)
    result = sa.chromatin_scatac_align(params)
    assert result["status"] == "success", result
    assert result["barcode_orientation"] == "reverse_complement"
    rc_path = tmp_path / "out" / "barcode_whitelist.rc.txt"
    assert rc_path.read_text() == "CGTT\n"
    command = next(command for step, command in calls if step == "align")
    assert command[command.index("--barcode-whitelist") + 1] == str(rc_path)
    assert any("reverse-complement" in warning for warning in result["warnings"])


def test_orientation_mismatch_fails_before_tool(tmp_path, monkeypatch):
    params = _alignment_params(tmp_path, barcode="CCCC")
    monkeypatch.setattr(sa, "_run", lambda *args: pytest.fail("must reject before alignment"))
    result = sa.chromatin_scatac_align(params)
    assert result["reason"] == "barcode_whitelist_mismatch"
    assert result["n_sampled"] == 24
    assert result["direct_match_rate"] == result["rc_match_rate"] == 0


def test_repeated_lane_flags_and_stage_invalidation(tmp_path, monkeypatch):
    params = _alignment_params(tmp_path, lanes=2)
    calls = _capture_alignment(monkeypatch)
    assert sa.chromatin_scatac_align(params)["status"] == "success"
    command = next(command for step, command in calls if step == "align")
    for flag, key in (("-1", "r1_fastq"), ("-2", "r3_fastq"), ("-b", "barcode_fastq")):
        assert [command[i + 1] for i, value in enumerate(command) if value == flag] == params[key]
    assert (tmp_path / "out" / "chromatin_scatac_align.stage_manifest.json").is_file()
    assert sa.chromatin_scatac_align(params)["resumed"] is True
    assert sum(step == "align" for step, _ in calls) == 1
    _write_fastq(Path(params["r1_fastq"][1]), ["AACG"] * 25)
    assert sa.chromatin_scatac_align(params).get("resumed") is not True
    assert sum(step == "align" for step, _ in calls) == 2
    Path(params["barcode_whitelist"]).write_text("AACG\nTTTT\n")
    assert sa.chromatin_scatac_align(params).get("resumed") is not True
    assert sum(step == "align" for step, _ in calls) == 3
    params["barcode_orientation"] = "direct"
    assert sa.chromatin_scatac_align(params).get("resumed") is not True


@pytest.mark.parametrize("platform", ["sci_atac", "surecell", "custom"])
def test_inline_platform_omits_barcode_and_whitelist(tmp_path, monkeypatch, platform):
    params = _alignment_params(tmp_path)
    params.pop("barcode_fastq")
    params.pop("barcode_whitelist")
    params.update(platform=platform, read_format="r1:4:-1,r2:0:-1,bc:0:4")
    calls = _capture_alignment(monkeypatch)
    result = sa.chromatin_scatac_align(params)
    assert result["status"] == "success", result
    command = next(command for step, command in calls if step == "align")
    assert "-b" not in command and "--barcode-whitelist" not in command
    assert command[command.index("--read-format") + 1] == params["read_format"]


def test_normalize_fastqs_preserves_literal_comma_and_checks_all_lanes(tmp_path):
    literal = tmp_path / "lane, one.fq"
    literal.write_text("x")
    other = tmp_path / "two.fq"
    other.write_text("x")
    assert sa._normalize_fastq_list(str(literal)) == [str(literal)]
    assert sa._normalize_fastq_list(f"{other},{other}") == [str(other), str(other)]
    with pytest.raises(ValueError):
        sa._normalize_fastq_list([str(other), str(tmp_path / "missing")])
