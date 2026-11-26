"""Permanent regression tests for D1 silent-failure fixes in bulk RNA compute kernels."""
from __future__ import annotations

import gzip

import pytest


def test_d1_003_gtf_truncated_must_not_return_partial_biotype(tmp_path):
    """D1-003: Truncated GTF must fail closed instead of returning a partial biotype_map."""
    from aria.scripts.rna_bulk.gtf_io import _parse_gtf_biotype_and_length

    gtf_gz = tmp_path / "trunc.gtf.gz"
    lines = [
        f'chr1\tHAVANA\tgene\t{1000 + i * 100}\t{1050 + i * 100}\t.\t+\t.\t'
        f'gene_id "ENSG{i:011d}"; gene_biotype "protein_coding";\n'
        for i in range(600)
    ]
    raw_gz = gzip.compress("".join(lines).encode("utf-8"))
    gtf_gz.write_bytes(raw_gz[: len(raw_gz) // 2])
    with pytest.raises((EOFError, OSError, ValueError)):
        _parse_gtf_biotype_and_length(str(gtf_gz))


def test_d1_006_bulk_counts_multifile_must_fail_closed_on_corrupt_sample(tmp_path):
    """D1-006: Multi-file bulk counts must fail closed when a declared sample file is corrupt/empty."""
    pytest.importorskip("pydeseq2")
    import numpy as np
    from aria.scripts.rna_bulk_de import bulk_rna_de

    rng = np.random.default_rng(7)
    genes = [f"GENE{i}" for i in range(60)]
    specs = [
        ("ctrl_1", "ctrl"), ("ctrl_2", "ctrl"), ("ctrl_3", "ctrl"), ("ctrl_4", "ctrl"),
        ("trt_1", "treat"), ("trt_2", "treat"), ("trt_3", "treat"), ("trt_4", "treat"),
    ]
    files = []
    for sid, cond in specs:
        p = tmp_path / f"{sid}.tsv"
        if sid == "trt_4":
            p.write_text("")
        else:
            vals = rng.poisson(50 if cond == "ctrl" else 200, size=len(genes))
            p.write_text("\n".join(f"{g}\t{v}" for g, v in zip(genes, vals)) + "\n")
        files.append(str(p))
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tcondition\n" + "\n".join(f"{s}\t{c}" for s, c in specs) + "\n")
    res = bulk_rna_de({
        "files": files,
        "metadata_file": str(meta),
        "design_factor": "condition",
        "contrasts": [{"numerator": "treat", "denominator": "ctrl", "name": "t_vs_c"}],
        "run_pathways": False,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "error"


def test_d1_009_sample_qc_crash_must_not_emit_high_conf_zero_outliers(tmp_path, monkeypatch):
    """D1-009: Crashed _sample_qc must emit an error block, never high-conf '0 outlier(s) flagged'."""
    pytest.importorskip("litellm")
    import pandas as pd
    import aria.scripts.rna_bulk.qc as qc_mod
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    counts = pd.DataFrame(
        {"s1": [10, 20, 30], "s2": [15, 25, 35], "s3": [12, 22, 32], "s4": [14, 24, 34]},
        index=["G1", "G2", "G3"],
    )
    meta = pd.DataFrame({"condition": ["A", "A", "B", "B"]}, index=["s1", "s2", "s3", "s4"])
    monkeypatch.setattr(
        qc_mod, "_run_vst", lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("VST crashed"))
    )
    sqc = qc_mod._sample_qc(counts, meta, str(tmp_path), [], design_factor="condition")
    assert sqc.get("status") == "error" and sqc.get("outliers") is None
    blocks = BulkRnaNarrator().collect("bulk_rna_agent", {"findings": {"sample_qc": sqc}})
    qc_block = next(b for b in blocks if b.id == "bulk.qc")
    assert qc_block.status == "error" and "0 outlier(s)" not in qc_block.claim


def test_d1_009_failed_qc_preserves_real_lib_size_and_safe_across_all_consumers(
    tmp_path, monkeypatch
):
    """R1/D1-009: Failed sample_qc preserves real size_ratio when computable, omits it if not,
    and never crashes or emits '0 outlier' / '1.0x' / HIGH confidence in any consumer."""
    pytest.importorskip("litellm")
    import pandas as pd
    import aria.scripts.rna_bulk.qc as qc_mod
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.bus.message_bus import Confidence

    # 1. Producer: VST fails, but library-size stats from counts are real (4.0x: 60 vs 240).
    counts = pd.DataFrame(
        {"s1": [10, 20, 30], "s2": [40, 80, 120], "s3": [20, 40, 60], "s4": [20, 40, 60]},
        index=["G1", "G2", "G3"],
    )
    meta = pd.DataFrame({"condition": ["A", "A", "B", "B"]}, index=["s1", "s2", "s3", "s4"])
    monkeypatch.setattr(
        qc_mod, "_run_vst", lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("VST crashed"))
    )
    sqc_with_ratio = qc_mod._sample_qc(counts, meta, str(tmp_path), [], design_factor="condition")
    assert sqc_with_ratio["status"] == "error"
    assert sqc_with_ratio["outliers"] is None
    assert sqc_with_ratio["size_ratio"] == 4.0
    assert sqc_with_ratio["lib_size_range"] == [60, 240]

    # When even counts summation fails, size_ratio and lib_size_range must be omitted (not defaulted).
    bad_counts = pd.DataFrame({"s1": ["bad"], "s2": ["data"]}, index=["G1"])
    sqc_bad_counts = qc_mod._sample_qc(bad_counts, meta, str(tmp_path), [], design_factor="condition")
    assert sqc_bad_counts["status"] == "error"
    assert sqc_bad_counts["outliers"] is None
    assert "size_ratio" not in sqc_bad_counts and "lib_size_range" not in sqc_bad_counts

    # Simulate rna_bulk_de copying candidate_outliers = None and a variant without size_ratio.
    sqc_with_ratio["candidate_outliers"] = None
    sqc_no_ratio = {
        "status": "error",
        "n_samples": 4,
        "outliers": None,
        "candidate_outliers": None,
        "pca_variance": [],
        "error": "counts unreadable",
    }

    narr = NarrativeAgent.__new__(NarrativeAgent)
    for label, sqc in (("with_ratio", sqc_with_ratio), ("no_ratio", sqc_no_ratio)):
        bulk_findings = {
            "sample_qc": sqc,
            "contrasts": [{
                "name": "treat_vs_ctrl",
                "status": "success",
                "n_significant": 12,
                "n_upregulated": 7,
                "n_downregulated": 5,
                "top_genes": [{"gene": "G1", "symbol": "G1", "log2fc": 2.0, "padj": 1e-3}],
            }],
        }
        agent_results = {"bulk_rna_agent": {"status": "done", "findings": bulk_findings}}

        # Consumer 1: _deterministic_bulk_executive_summary
        exec_sum = narr._deterministic_bulk_executive_summary({}, {}, agent_results)
        assert "0 outlier" not in exec_sum.lower(), (label, exec_sum)
        assert "failed" in exec_sum.lower() or "not available" in exec_sum.lower(), (label, exec_sum)

        # Consumer 2: _summarize_agent_results_for_llm (must not raise TypeError on len(None))
        llm_sum = narr._summarize_agent_results_for_llm(agent_results)
        assert "0 outlier" not in llm_sum.lower(), (label, llm_sum)
        assert "failed" in llm_sum.lower() or "not available" in llm_sum.lower(), (label, llm_sum)
        if label == "no_ratio":
            assert "0×" not in llm_sum and "1.0×" not in llm_sum, llm_sum
            assert "not available" in llm_sum.lower(), llm_sum
        else:
            assert "4.0×" in llm_sum, llm_sum

        # Consumer 3: _summarize_bulk_rna
        bulk_sum = narr._summarize_bulk_rna(
            agent_results["bulk_rna_agent"],
            {"high": [], "medium": [], "low": [], "insufficient": []},
        )
        assert "0 outlier" not in bulk_sum.lower() and "no outliers" not in bulk_sum.lower(), (label, bulk_sum)
        assert "1.0×" not in bulk_sum and "1.0x" not in bulk_sum, (label, bulk_sum)
        assert "passed" not in bulk_sum.lower(), (label, bulk_sum)
        assert "failed" in bulk_sum.lower() or "not available" in bulk_sum.lower(), (label, bulk_sum)

        # Consumer 4: BulkRNAAgent._publish_findings
        published = []
        bulk_agent = BulkRNAAgent.__new__(BulkRNAAgent)
        bulk_agent.publish_finding = lambda exp_id, payload, conf: published.append((payload, conf))
        bulk_agent._publish_findings("exp1", bulk_findings)
        qc_pubs = [(p, c) for p, c in published if "qc" in str(p.get("summary", "")).lower()]
        assert len(qc_pubs) == 1, published
        qc_payload, qc_conf = qc_pubs[0]
        assert qc_conf != Confidence.HIGH, (label, qc_conf)
        assert "1.0x" not in qc_payload["summary"].lower(), (label, qc_payload)
        assert "failed" in qc_payload["summary"].lower() or "not available" in qc_payload["summary"].lower()

        # Consumers 5 & 6: BulkRnaNarrator._qc_blocks & ReportSectionsMixin._build_qc_section
        blocks = BulkRnaNarrator().collect("bulk_rna_agent", {"findings": bulk_findings})
        qc_block = next(b for b in blocks if b.id == "bulk.qc")
        assert qc_block.status == "error" and qc_block.confidence != "high"
        qc_html = narr._build_qc_section(
            {"high": [], "medium": [], "low": [], "insufficient": []},
            agent_results=agent_results,
        )
        assert "sample QC failed" in qc_html and "no outliers flagged" not in qc_html


def test_d1_014_lobo_error_must_not_be_masked_as_eligible_not_executed():
    """D1-014: LOBO jackknife status='error' must surface as status='error' in falsification certificates."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.falsification_certificates import build_falsification_certificates

    agent_results = {
        "bulk_rna_agent": {
            "findings": {
                "contrasts": [{
                    "name": "treat_vs_ctrl",
                    "status": "success",
                    "n_replicates": {"test": 4, "ref": 4},
                    "lobo_jackknife": {"status": "error", "details": "solver failed"},
                }]
            }
        }
    }
    for claim in (
        {"claim_id": "c1", "analysis": "differential_expression", "tier": "associative"},
        {"claim_id": "c2", "modality": "bulk RNA-seq", "analysis": "differential_expression", "tier": "associative"},
    ):
        certs = build_falsification_certificates([claim], agent_results, {})
        axis = certs["certificates"][0]["lobo_biological_unit_jackknife"]
        assert axis["status"] == "error", axis
        assert any("solver failed" in d for d in (axis.get("error_details") or []))


def test_d1_023_run_vst_inner_fallback_must_library_size_normalize(monkeypatch):
    """D1-023: Failing vst_fit + deseq2 must fall back to library-size-normalized log2, not raw log2."""
    pytest.importorskip("pydeseq2")
    import numpy as np
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from aria.scripts.rna_bulk.transforms import _run_vst

    counts = pd.DataFrame(
        {"s1": [10, 20, 30, 40], "s2": [100, 200, 300, 400], "s3": [10, 20, 30, 40], "s4": [100, 200, 300, 400]},
        index=["G1", "G2", "G3", "G4"],
    )
    meta = pd.DataFrame({"condition": ["A", "A", "B", "B"]}, index=counts.columns)
    monkeypatch.setattr(
        DeseqDataSet, "vst_fit", lambda self, **kw: (_ for _ in ()).throw(RuntimeError("vst fail"))
    )
    monkeypatch.setattr(
        DeseqDataSet, "deseq2", lambda self, **kw: (_ for _ in ()).throw(RuntimeError("deseq2 fail"))
    )
    vst_df = _run_vst(counts, meta, [])
    assert np.allclose(vst_df["s1"].values, vst_df["s2"].values, atol=0.1)


def test_d1_011_plot_failure_must_not_disable_pca_outlier_detection(tmp_path):
    """D1-011: Plot I/O failure must not zero PCA coords or suppress real replicate outliers."""
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")
    import numpy as np
    import pandas as pd
    from aria.scripts.rna_bulk.qc import _sample_qc

    rng = np.random.default_rng(1)
    base_a = rng.integers(50, 500, size=200)
    base_b = rng.integers(50, 500, size=200)
    s_bad = base_a.copy()
    s_bad[:25] = 50_000
    counts = pd.DataFrame(
        {
            "a1": base_a, "a2": base_a + rng.integers(0, 5, 200),
            "a3": base_a + rng.integers(0, 5, 200), "a4": base_a + rng.integers(0, 5, 200),
            "a5": base_a + rng.integers(0, 5, 200), "a6": base_a + rng.integers(0, 5, 200),
            "a_out": s_bad, "b1": base_b, "b2": base_b + rng.integers(0, 5, 200),
            "b3": base_b + rng.integers(0, 5, 200),
        },
        index=[f"G{i}" for i in range(200)],
    )
    meta = pd.DataFrame({"condition": ["A"] * 7 + ["B"] * 3}, index=list(counts.columns))
    ok_dir = tmp_path / "ok"
    ok_dir.mkdir()
    assert "a_out" in _sample_qc(counts, meta, str(ok_dir), [], design_factor="condition")["outliers"]
    bad_dir = tmp_path / "not_a_dir"
    bad_dir.write_text("x")
    assert "a_out" in _sample_qc(counts, meta, str(bad_dir), [], design_factor="condition")["outliers"]
