"""Regression tests for D7 layer-contract fixes in chromatin & scATAC."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def test_d7_006_cross_modal_triangulation_reads_scatac_da_peaks_by_cluster():
    """D7-006: chromatin_diffacc._per_cluster_accessibility writes per-cluster DA
    peaks under 'da_peaks_by_cluster'. cross_modal_triangulation must read that
    canonical key and emit atac_da_association edges."""
    from aria.agents.narrative.synthesis.cross_modal_triangulation import (
        build_cross_modal_evidence_graph,
    )
    from aria.scripts.chromatin_diffacc import _per_cluster_accessibility

    rng = np.random.default_rng(5)
    n_cells = 60
    X = rng.uniform(0.0, 0.2, size=(n_cells, 10)).astype(np.float32)
    X[:30, 0] += 5.0
    X[:30, 1] += 4.5
    X[30:, 2] += 5.0
    peaks = [f"chr1:{1000 + i * 1000}-{1500 + i * 1000}" for i in range(10)]
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"leiden": ["0"] * 30 + ["1"] * 30}, index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=peaks),
    )
    pc_res = _per_cluster_accessibility(adata, groupby="leiden", padj_max=0.05, lfc_min=0.5, top_n=10, warnings=[])
    assert pc_res["ran"] is True
    assert pc_res["n_da_total"] >= 2
    assert any(pc_res["da_peaks_by_cluster"].values())

    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "scATAC": {
                    "findings": {
                        "differential_accessibility": {
                            "status": "success",
                            "per_cluster": pc_res,
                        }
                    }
                }
            },
        },
        "integration_agent": {
            "status": "done",
            "findings": {
                "peak_to_gene": {
                    "status": "success",
                    "top_links": [
                        {"peak": peaks[0], "gene": "GENE_A", "correlation": 0.72, "padj": 1e-6}
                    ],
                }
            },
        },
    }
    graph = build_cross_modal_evidence_graph(agent_results, exp_ctx={})
    atac_edges = [e for e in graph.get("edges", []) if e.get("relation") == "atac_da_association"]
    assert len(atac_edges) > 0, (
        f"Expected atac_da_association edges for {peaks[0]}, got 0 (da_peaks_by_cluster ignored)"
    )


def test_d7_009_chromatin_tobias_forwards_paired_replicates_from_design_agent(tmp_path):
    """D7-009: DesignAgent._build_design encodes paired designs in
    replicate_handling['mode'] == 'paired_biological_block'.
    ChromatinAgent._run_bulk_atac_footprinting must forward paired_replicates=True
    to chromatin_footprint_tobias.py."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.agents.design_agent import DesignAgent
    from aria.memory.memory import ARIAMemory

    da = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da._organism = "Homo sapiens"
    da._genome = "hg38"
    da._confirmed_groups = {
        "ctrl": ["D1_ctrl", "D2_ctrl", "D3_ctrl"],
        "treat": ["D1_treat", "D2_treat", "D3_treat"],
    }
    da._main_factor = "condition"
    da._batch_covariate = "donor"
    da._paired_covariate_values = {
        "donor": {
            "D1_ctrl": "D1", "D1_treat": "D1",
            "D2_ctrl": "D2", "D2_treat": "D2",
            "D3_ctrl": "D3", "D3_treat": "D3",
        }
    }
    da._replicate_handling = da._build_replicate_handling(
        technical=False,
        paired_covariate="donor",
        sample_to_pair=da._paired_covariate_values["donor"],
    )
    da._inferred_design = {"source": "sample_sheet"}
    design = da._build_design()
    assert design["replicate_handling"]["mode"] == "paired_biological_block"

    bams = []
    for name in ["D1_ctrl", "D2_ctrl", "D3_ctrl", "D1_treat", "D2_treat", "D3_treat"]:
        p = tmp_path / f"{name}.bam"
        p.write_bytes(b"")
        bams.append(str(p))
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n")
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t0\t100\n")
    meme = tmp_path / "m.meme"
    meme.write_text("MEME version 4\n")

    captured = {}

    class _CapEnv:
        def run_in_stack(self, stack, script_path, params):
            captured["params"] = dict(params)
            return {"status": "success", "ran": False, "reason": "captured"}

    agent = ChromatinAgent(ARIAMemory(":memory:"), llm=None)
    agent.env = _CapEnv()
    exp_ctx = {
        "design": design,
        "genome": "hg38",
        "genome_fasta": str(fa),
        "peaks_bed": str(bed),
        "motif_meme": str(meme),
        "comparisons": [{"test": "treat", "reference": "ctrl"}],
    }
    agent._run_bulk_atac_footprinting(
        "exp1",
        exp_ctx,
        da_result={"comparisons": [{"test": "treat", "reference": "ctrl"}]},
        peaks_result={},
        files=bams,
    )
    assert captured["params"].get("paired_replicates") is True, (
        f"Paired design ({design['replicate_handling']['mode']}) was silently dropped: "
        f"paired_replicates={captured['params'].get('paired_replicates')!r}"
    )


def test_d7_016_warning_router_surfaces_chromatin_pass_qc_false_and_qc_failures():
    """D7-016: ChromatinAgent._merge_peak_frip_into_qc populates qc['pass_qc']=False
    and qc['qc_failures']. Through NarrativeRegistry.collect_blocks, warning_router
    must attach a structured warning caveat naming the failed ENCODE gates."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.registry import registry_with

    findings = {
        "qc": {
            "status": "success",
            "data_type": "bulk_ATAC",
            "n_samples": 4,
            "tss_enrichment": 2.1,
            "nrf": 0.42,
            "mito_fraction": 0.45,
            "metrics_not_computed": ["frip"],
            "qc_complete": False,
            "pass_qc": None,
            "qc_failures": [],
            "warnings": [],
        }
    }
    ChromatinAgent._merge_peak_frip_into_qc(
        findings,
        {"data_type": "bulk_ATAC", "frip": 0.04, "frip_per_sample": {"s1": 0.04}},
    )
    qc = findings["qc"]
    assert qc["pass_qc"] is False
    assert len(qc["qc_failures"]) >= 2
    # Clear free-text warnings to prove the STRUCTURED pass_qc + qc_failures
    # contract routes independently of free-text warnings.
    qc["warnings"] = []

    reg = registry_with((ChromatinNarrator(),))
    blocks = reg.collect_blocks({"chromatin_agent": {"status": "done", "findings": findings}})
    qc_block = next(b for b in blocks if b.id == "chromatin.qc.bulk_ATAC")
    warn_caveats = [c for c in qc_block.caveats if c.severity == "warning"]
    assert warn_caveats, f"Expected structured warning caveat for pass_qc=False, got {qc_block.caveats}"
    assert any("pass_qc=False" in c.text and "FRiP" in c.text for c in warn_caveats)
    sig_fields = {s["field"] for s in qc_block.metadata.get("degradation_signals", [])}
    assert "pass_qc" in sig_fields


def test_d7_017_warning_router_surfaces_chromatin_partial_pairing_and_dropped_donors():
    """D7-017 (chromatin): partial pairing with excluded_unpaired_samples (scATAC DA)
    and excluded_unpaired_replicates (TOBIAS footprinting) must route through
    warning_router as structured warning caveats naming dropped donors."""
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.registry import registry_with

    reg = registry_with((ChromatinNarrator(),))
    agent_results = {
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "differential_accessibility": {
                    "status": "success",
                    "data_type": "scATAC",
                    "pseudobulk": {
                        "ran": True,
                        "comparisons": [{
                            "test": "treat",
                            "reference": "ctrl",
                            "status": "success",
                            "n_sig": 4,
                            "pairing_status": "partial",
                            "excluded_unpaired_samples": ["D4_ctrl"],
                        }],
                    },
                },
                "footprinting": {
                    "status": "success",
                    "ran": True,
                    "group_a": "treat",
                    "group_b": "ctrl",
                    "differential_summary": {
                        "parsed": True,
                        "n_motifs_tested": 100,
                        "n_significant": 5,
                        "inference": {
                            "status": "success",
                            "test": "paired_t",
                            "replicates_per_condition": {"treat": 3, "ctrl": 3},
                            "pairing_status": "partial_complete_pair_subset",
                            "excluded_unpaired_replicates": ["D4"],
                            "null_label_controls": {"n_permutations": 10},
                        },
                    },
                },
            },
        }
    }
    blocks = reg.collect_blocks(agent_results)
    pb_block = next(b for b in blocks if b.id == "chromatin.differential_accessibility.pseudobulk")
    fp_block = next(b for b in blocks if b.id == "chromatin.differential_tf_footprinting")
    assert any(c.severity == "warning" and "D4_ctrl" in c.text for c in pb_block.caveats)
    assert any(s["field"] == "pairing_status" for s in pb_block.metadata.get("degradation_signals", []))
    assert any(c.severity == "warning" and "D4" in c.text for c in fp_block.caveats)
    assert any(s["field"] == "pairing_status" for s in fp_block.metadata.get("degradation_signals", []))
