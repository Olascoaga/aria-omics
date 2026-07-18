"""E4 regression guards for the bulk-RNA FASTQ preprocessing contract."""

from pathlib import Path


def _touch_fastqs(tmp_path: Path, names: list[str]) -> list[str]:
    paths = []
    for name in names:
        path = tmp_path / name
        path.write_bytes(b"fastq")
        paths.append(str(path))
    return paths


def test_e4_multilane_paired_fastqs_form_one_sample_manifest(tmp_path):
    from aria.scripts.rna_fastq_qc import _detect_samples

    selected = _touch_fastqs(tmp_path, [
        "patient_a_S7_L001_R1_001.fastq.gz",
        "patient_a_S7_L001_R2_001.fastq.gz",
        "patient_a_S7_L002_R1_001.fastq.gz",
        "patient_a_S7_L002_R2_001.fastq.gz",
    ])
    warnings = []

    samples = _detect_samples(tmp_path, warnings, fastq_files=selected)

    assert warnings == []
    assert len(samples) == 1
    sample = samples[0]
    assert sample["name"] == "patient_a"
    assert sample["read_layout"] == "paired-end"
    assert sample["paired"] is True
    assert [lane["lane"] for lane in sample["lanes"]] == ["L001", "L002"]
    assert [Path(lane["r1"]).name for lane in sample["lanes"]] == [
        "patient_a_S7_L001_R1_001.fastq.gz",
        "patient_a_S7_L002_R1_001.fastq.gz",
    ]
    assert all(lane["r2"] for lane in sample["lanes"])


def test_e4_fastp_and_star_consume_every_lane(tmp_path, monkeypatch):
    from aria.scripts import rna_align as align
    from aria.scripts import rna_fastq_qc as qc

    raw = _touch_fastqs(tmp_path, [
        "sample_L001_R1_001.fastq.gz",
        "sample_L001_R2_001.fastq.gz",
        "sample_L002_R1_001.fastq.gz",
        "sample_L002_R2_001.fastq.gz",
    ])
    sample = qc._detect_samples(tmp_path, [], fastq_files=raw)[0]
    lane_calls = []

    def fake_fastp_lane(sample, trimmed_dir, **kwargs):
        lane_calls.append(sample)
        r1 = trimmed_dir / f"{sample['name']}_R1_trimmed.fq.gz"
        r2 = trimmed_dir / f"{sample['name']}_R2_trimmed.fq.gz"
        return {
            "name": sample["name"],
            "status": "success",
            "r1_trimmed": str(r1),
            "r2_trimmed": str(r2),
            "n_reads_raw": 100,
            "n_reads_trimmed": 90,
            "q30_rate": 95.0,
        }

    monkeypatch.setattr(qc, "_run_fastp_lane", fake_fastp_lane)
    processed = qc._run_fastp(
        sample, tmp_path / "trimmed", tmp_path / "reports", 4, 36, 20, []
    )

    assert [call["name"] for call in lane_calls] == [
        "sample_L001", "sample_L002"
    ]
    assert processed["n_reads_raw"] == 200
    assert processed["n_reads_trimmed"] == 180
    assert len(processed["r1_trimmed_files"]) == 2
    assert len(processed["r2_trimmed_files"]) == 2

    star_commands = []

    class _Proc:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(cmd, **kwargs):
        if cmd[0] == "STAR":
            star_commands.append(cmd)
            prefix = cmd[cmd.index("--outFileNamePrefix") + 1]
            Path(f"{prefix}Aligned.sortedByCoord.out.bam").write_bytes(b"bam")
            Path(f"{prefix}Log.final.out").write_text(
                "Number of input reads | 180\n"
                "Uniquely mapped reads % | 90.0%\n"
            )
        return _Proc()

    monkeypatch.setattr(align.subprocess, "run", fake_run)
    monkeypatch.setattr(align, "stage_is_current", lambda *a, **k: (False, ""))
    (tmp_path / "aligned").mkdir()
    aligned = align._align_sample(
        processed,
        genome_dir=tmp_path / "star",
        output_dir=tmp_path / "aligned",
        threads=4,
        two_pass=False,
        warnings=[],
    )

    assert aligned["status"] == "success"
    command = star_commands[0]
    reads_index = command.index("--readFilesIn")
    assert command[reads_index + 1].split(",") == \
        processed["r1_trimmed_files"]
    assert command[reads_index + 2].split(",") == \
        processed["r2_trimmed_files"]
    assert aligned["read_layout"] == "paired-end"


def test_e4_real_single_end_layout_reaches_featurecounts_without_pair_flags(
    tmp_path, monkeypatch
):
    from aria.scripts import rna_quantify as quantify

    bam = tmp_path / "single.bam"
    bam.write_bytes(b"bam")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("")
    commands = []

    class _Proc:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(cmd, **kwargs):
        commands.append(cmd)
        return _Proc()

    monkeypatch.setattr(quantify.subprocess, "run", fake_run)
    monkeypatch.setattr(
        quantify,
        "_clean_counts_matrix",
        lambda *args, **kwargs: (tmp_path / "counts_matrix.tsv", 10),
    )

    result = quantify.rna_quantify({
        "bam_files": [{
            "name": "single",
            "status": "success",
            "bam": str(bam),
            "read_layout": "single-end",
        }],
        "gtf_file": str(gtf),
        "output_dir": str(tmp_path / "counts"),
        "strand": 0,
    })

    assert result["status"] == "success"
    assert result["read_layout"] == "single-end"
    assert "-p" not in commands[0]
    assert "--countReadPairs" not in commands[0]


def test_e4_qc_is_all_or_fail_when_one_sample_fails(tmp_path, monkeypatch):
    from aria.scripts import rna_fastq_qc as qc

    selected = _touch_fastqs(tmp_path, [
        "ok_R1.fastq.gz",
        "bad_R1.fastq.gz",
    ])

    def fake_fastp(sample, **kwargs):
        return {
            "name": sample["name"],
            "status": "failed" if sample["name"] == "bad" else "success",
            "read_layout": sample["read_layout"],
        }

    monkeypatch.setattr(qc, "_run_fastp", fake_fastp)
    monkeypatch.setattr(qc, "_run_multiqc", lambda *args, **kwargs: None)

    result = qc.rna_fastq_qc({
        "fastq_dir": str(tmp_path),
        "fastq_files": selected,
        "output_dir": str(tmp_path / "qc"),
    })

    assert result["status"] == "error"
    assert result["error_type"] == "PartialFastqQCFailure"
    assert result["failed_samples"] == ["bad"]
    assert len(result["samples"]) == 2


def test_e4_alignment_is_all_or_fail_when_one_sample_fails(
    tmp_path, monkeypatch
):
    from aria.scripts import rna_align as align

    genome = tmp_path / "star"
    genome.mkdir()
    for name in ("SA", "SAindex", "Genome"):
        (genome / name).write_bytes(b"index")

    monkeypatch.setattr(align, "_star_version", lambda: "STAR test")
    monkeypatch.setattr(
        align,
        "_align_sample",
        lambda sample, **kwargs: {
            "name": sample["name"],
            "status": "failed" if sample["name"] == "bad" else "success",
            "bam": str(tmp_path / f"{sample['name']}.bam"),
            "read_layout": sample["read_layout"],
        },
    )

    result = align.rna_align({
        "samples": [
            {"name": "ok", "r1_trimmed": "ok.fq.gz",
             "read_layout": "single-end", "paired": False},
            {"name": "bad", "r1_trimmed": "bad.fq.gz",
             "read_layout": "single-end", "paired": False},
        ],
        "genome_dir": str(genome),
        "output_dir": str(tmp_path / "aligned"),
    })

    assert result["status"] == "error"
    assert result["error_type"] == "PartialAlignmentFailure"
    assert result["failed_samples"] == ["bad"]


def test_e4_bulk_agent_propagates_selection_and_does_not_force_paired(
    tmp_path, monkeypatch
):
    from aria.agents import bulk_rna_agent as bulk_module
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    monkeypatch.setattr(bulk_module.os, "cpu_count", lambda: 32)

    fastq = tmp_path / "single_R1.fastq.gz"
    fastq.write_bytes(b"fastq")
    calls = []

    class _Env:
        def run_in_stack(self, stack, script_path, params):
            calls.append((script_path, params))
            if script_path.endswith("rna_fastq_qc.py"):
                return {
                    "status": "success",
                    "n_samples": 1,
                    "samples": [{
                        "name": "single", "status": "success",
                        "read_layout": "single-end", "paired": False,
                        "r1_trimmed": "/tmp/single.fq.gz",
                        "pct_passed": 95.0,
                    }],
                }
            if script_path.endswith("rna_align.py"):
                return {
                    "status": "success",
                    "n_aligned": 1,
                    "bam_files": [{
                        "name": "single", "status": "success",
                        "bam": "/tmp/single.bam",
                        "read_layout": "single-end", "pct_unique": 90.0,
                    }],
                }
            return {
                "status": "success", "counts_matrix": "/tmp/counts.tsv",
                "n_genes": 10, "n_samples": 1,
                "read_layout": "single-end",
            }

    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    agent.env = _Env()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_finding = lambda *args, **kwargs: None

    counts, _ = agent._run_preprocessing(
        "exp", [str(fastq)], {"genome_config": {}, "n_cpus": 30}, {}
    )

    assert counts == ["/tmp/counts.tsv"]
    qc_params = calls[0][1]
    quant_params = calls[2][1]
    assert qc_params["fastq_files"] == [str(fastq)]
    assert [params["threads"] for _, params in calls] == [30, 30, 30]
    assert "paired" not in quant_params


def test_e4_cpu_requests_are_bounded_without_changing_ordinary_defaults(
    monkeypatch,
):
    from aria.agents import bulk_rna_agent as bulk_module
    from aria.scripts import rna_bulk_de

    monkeypatch.setattr(bulk_module.os, "cpu_count", lambda: 12)
    monkeypatch.setattr(rna_bulk_de.os, "cpu_count", lambda: 12)

    assert bulk_module._preprocessing_threads({"n_cpus": 30}) == 12
    assert bulk_module._preprocessing_threads({"n_cpus": 0}) == 1
    assert bulk_module._preprocessing_threads({}) == 8
    assert rna_bulk_de._bounded_n_cpus(30) == 12
    assert rna_bulk_de._bounded_n_cpus(0) == 1
    assert rna_bulk_de._bounded_n_cpus(None) is None
