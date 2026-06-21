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
    assert c is not None and c.validation_level == "beta"


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
