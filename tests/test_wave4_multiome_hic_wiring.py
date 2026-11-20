"""Wave 4C regression tests: Multiome FASTQ/barcode wiring & Hi-C GenomeArchAgent wiring."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from aria.agents.genome_arch_agent import GenomeArchAgent
from aria.agents.integration_agent import IntegrationAgent
from aria.memory.memory import ARIAMemory
from aria.scripts.integration_assemble_multiome import (
    _translate_barcodes,
    integration_assemble_multiome,
)
from aria.utils.multiome_contracts import _unpaired_container_files


def _make_integration_agent(tmp_path: Path) -> IntegrationAgent:
    mem = ARIAMemory(db_path=str(tmp_path / "integ_mem.sqlite"))
    llm = MagicMock()
    agent = IntegrationAgent(mem, llm)
    agent.env = MagicMock()
    return agent


def test_integration_agent_raw_fastq_wires_barcode_namespace_and_barcode_map(tmp_path: Path):
    agent = _make_integration_agent(tmp_path)
    out_h5mu = tmp_path / "assembled.h5mu"
    out_h5mu.write_bytes(b"mock")
    bmap = tmp_path / "barcode_map.tsv"
    bmap.write_text("rna_barcode\tatac_barcode\nR1\tA1\n", encoding="utf-8")

    captured: list[dict] = []

    def _fake_run_in_stack(*, stack, script_path, params):
        captured.append({"stack": stack, "script_path": script_path, "params": dict(params)})
        return {"status": "success", "output_path": str(out_h5mu)}

    agent.env.run_in_stack.side_effect = _fake_run_in_stack

    exp_ctx = {
        "output_dir": str(tmp_path / "out"),
        "multiome_fastq_manifest": {
            "platform": "10x_multiome",
            "samples": [
                {
                    "sample_id": "S1",
                    "rna_fastqs": [str(tmp_path / "r1.fq.gz"), str(tmp_path / "r2.fq.gz")],
                    "atac_fastqs": [
                        str(tmp_path / "a1.fq.gz"),
                        str(tmp_path / "a2.fq.gz"),
                        str(tmp_path / "ab.fq.gz"),
                    ],
                }
            ],
        },
    }
    gate = {
        "allowed": True,
        "evidence_paths": [],
        "contract": {
            "object_type": "raw_multiome_fastq",
            "platform": "10x_multiome",
            "barcode_map": str(bmap),
        },
    }

    res = agent._run_multiome_preprocessing(
        "exp_w4c_raw", exp_ctx, {}, gate, exp_ctx.get("modalities") or {}
    )
    assert res is not None and res["status"] == "success"
    assert len(captured) == 1
    p = captured[0]["params"]
    assert p["barcode_map_path"] == str(bmap)
    assert p["barcode_namespace"] == "barcode_map_tsv"


def test_integration_agent_paired_matrices_sets_barcode_namespace_not_legacy_key(tmp_path: Path):
    agent = _make_integration_agent(tmp_path)
    rna_p = tmp_path / "rna.h5ad"
    atac_p = tmp_path / "atac.h5ad"
    gex_wl = tmp_path / "gex_wl.txt"
    atac_wl = tmp_path / "atac_wl.txt"
    for p in (rna_p, atac_p, gex_wl, atac_wl):
        p.write_text("x\n", encoding="utf-8")

    captured: list[dict] = []

    def _fake_run_in_stack(*, stack, script_path, params):
        captured.append(dict(params))
        return {"status": "success", "output_path": str(tmp_path / "paired.h5mu")}

    agent.env.run_in_stack.side_effect = _fake_run_in_stack

    exp_ctx = {
        "output_dir": str(tmp_path / "out"),
        "modalities": {"scRNA": [str(rna_p)], "scATAC": [str(atac_p)]},
    }
    gate = {
        "allowed": True,
        "evidence_paths": [],
        "contract": {
            "object_type": "explicit_same_cell",
            "cell_namespace": "10x_arc_whitelist_pair",
            "gex_whitelist_path": str(gex_wl),
            "atac_whitelist_path": str(atac_wl),
        },
    }

    res = agent._run_multiome_preprocessing(
        "exp_w4c_mat", exp_ctx, {}, gate, exp_ctx["modalities"]
    )
    assert res is not None and res["status"] == "success"
    assert len(captured) == 1
    assert captured[0]["barcode_namespace"] == "10x_arc_whitelist_pair"
    assert "barcode_translation" not in captured[0]


def test_store_eregulon_tunnels_keeps_exact_zero_padj_and_skips_nan(tmp_path: Path):
    agent = _make_integration_agent(tmp_path)
    agent.memory = MagicMock()

    ereg_result = {
        "top_eregulons": [
            {
                "tf": "SOX2",
                "regulon_id": "SOX2_pos",
                "sign_class": "pos_tf_pos_p2g",
                "n_enhancer_peaks": 12,
                "n_target_genes": 8,
                "regulon_padj": 0.0,  # Exact float64 underflow must NOT be dropped!
                "depth_confounded": False,
                "motif_family_ambiguity": False,
            },
            {
                "tf": "NAN_TF",
                "regulon_id": "NAN_pos",
                "regulon_padj": float("nan"),
                "depth_confounded": False,
            },
            {
                "tf": "INF_TF",
                "regulon_id": "INF_pos",
                "regulon_padj": float("inf"),
                "depth_confounded": False,
            },
            {
                "tf": "NEG_TF",
                "regulon_id": "NEG_pos",
                "regulon_padj": -0.01,
                "depth_confounded": False,
            },
            {
                "tf": "NS_TF",
                "regulon_id": "NS_pos",
                "regulon_padj": 0.25,
                "depth_confounded": False,
            },
        ]
    }
    agent._store_eregulon_tunnels("exp_w4c_reg", ereg_result)
    assert agent.memory.create_tunnel.call_count == 1
    kwargs = agent.memory.create_tunnel.call_args.kwargs
    assert kwargs["entity"] == "SOX2"
    assert kwargs["confidence"] == "high"
    assert "padj=0" in kwargs["description"]


def test_assemble_multiome_translates_fragments_and_accepts_platform_alias(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    mudata = pytest.importorskip("mudata")
    n_cells = 12
    rna_bcs = [f"RNA_CELL_{i:02d}" for i in range(n_cells)]
    atac_bcs = [f"ATAC_CELL_{i:02d}" for i in range(n_cells)]

    rna = ad.AnnData(
        X=sp.csr_matrix(np.ones((n_cells, 4), dtype=np.float32)),
        obs=pd.DataFrame(index=rna_bcs),
        var=pd.DataFrame(index=[f"G{i}" for i in range(4)]),
    )
    rna_path = tmp_path / "rna.h5ad"
    rna.write_h5ad(rna_path)

    peaks_path = tmp_path / "peaks.bed"
    peaks_path.write_text(
        "chr1\t100\t500\nchr1\t1000\t1500\n",
        encoding="utf-8",
    )

    frag_lines: list[str] = []
    for i, abc in enumerate(atac_bcs):
        # 4 nuclear fragments in peak1 + 1 mitochondrial fragment -> total=5, mito_frac=0.2
        frag_lines.append(f"chr1\t150\t350\t{abc}\t4")
        frag_lines.append(f"chrM\t100\t200\t{abc}\t1")
    frags_path = tmp_path / "fragments.tsv"
    frags_path.write_text("\n".join(frag_lines) + "\n", encoding="utf-8")

    bmap_path = tmp_path / "barcode_map.tsv"
    bmap_lines = ["rna_barcode\tatac_barcode"] + [
        f"{r}\t{a}" for r, a in zip(rna_bcs, atac_bcs)
    ]
    bmap_path.write_text("\n".join(bmap_lines) + "\n", encoding="utf-8")

    # 1) Auto-infer barcode_map_tsv even when barcode_namespace is omitted, AND translate fragment barcodes!
    res = integration_assemble_multiome(
        {
            "mode": "paired_matrices",
            "rna_path": str(rna_path),
            "fragments_path": str(frags_path),
            "peaks_path": str(peaks_path),
            "barcode_map_path": str(bmap_path),
            "output_dir": str(tmp_path / "assembled_paired"),
            "min_shared_cells": 10,
            "resume": False,
        }
    )
    assert res["status"] == "success", res
    assert res["barcode_translation"]["mode"] == "barcode_map_tsv"
    assert res["barcode_translation"]["n_translated"] == n_cells

    mdata = mudata.read_h5mu(res["output_path"])
    atac_obs = mdata["atac"].obs
    assert np.all(np.isfinite(atac_obs["atac_fragments"].to_numpy(dtype=float)))
    assert np.allclose(atac_obs["atac_fragments"].to_numpy(dtype=float), 5.0)
    assert np.allclose(atac_obs["atac_mito_frac"].to_numpy(dtype=float), 0.2)

    # 2) Raw FASTQ mode with delegated outputs accepts platform alias "10x_multiome"
    #    AND populates atac_fragments / atac_mito_frac from atac_fragments_out!
    res_raw = integration_assemble_multiome(
        {
            "mode": "raw_multiome_fastq",
            "platform": "10x_multiome",
            "rna_counts_out": str(rna_path),
            "atac_fragments_out": str(frags_path),
            "peaks_path": str(peaks_path),
            "barcode_map_path": str(bmap_path),
            "output_dir": str(tmp_path / "assembled_raw"),
            "min_shared_cells": 10,
            "resume": False,
        }
    )
    assert res_raw["status"] == "success", res_raw
    mdata_raw = mudata.read_h5mu(res_raw["output_path"])
    assert np.allclose(mdata_raw["atac"].obs["atac_fragments"].to_numpy(dtype=float), 5.0)
    assert np.allclose(mdata_raw["atac"].obs["atac_mito_frac"].to_numpy(dtype=float), 0.2)


def test_translate_barcodes_accepts_legacy_barcode_translation_alias(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    rna = ad.AnnData(X=np.ones((2, 2)), obs=pd.DataFrame(index=["R1", "R2"]))
    atac = ad.AnnData(X=np.ones((2, 2)), obs=pd.DataFrame(index=["A1", "A2"]))
    bmap = tmp_path / "map.tsv"
    bmap.write_text("rna_barcode\tatac_barcode\nR1\tA1\nR2\tA2\n", encoding="utf-8")
    tr, keep, names = _translate_barcodes(
        rna,
        atac,
        {"barcode_translation": "barcode_map_tsv", "barcode_map_path": str(bmap)},
    )
    assert tr["mode"] == "barcode_map_tsv"
    assert keep == [0, 1]
    assert names == ["R1", "R2"]


def test_unpaired_container_files_normalizes_relative_vs_absolute_paths(tmp_path: Path, monkeypatch):
    sub = tmp_path / "cohort"
    sub.mkdir()
    h5mu = sub / "sample1.h5mu"
    h5mu.write_bytes(b"mock")
    monkeypatch.chdir(tmp_path)

    rel_path = os.path.join("cohort", "sample1.h5mu")
    abs_path = str(h5mu.resolve())
    exp_ctx = {"modalities": {"scRNA": [rel_path], "scATAC": [abs_path]}}
    assert _unpaired_container_files(exp_ctx, [abs_path]) == []
    assert _unpaired_container_files(exp_ctx, [rel_path]) == []


def test_genome_arch_agent_microc_loop_res_svg_propagation_and_unphased_summary(tmp_path: Path):
    mem = ARIAMemory(db_path=str(tmp_path / "hic_mem.sqlite"))
    llm = MagicMock()
    agent = GenomeArchAgent(mem, llm)
    agent.env = MagicMock()

    mcool = tmp_path / "microc.mcool"
    mcool.write_bytes(b"HDF5")
    calls: list[tuple[str, dict]] = []

    def _fake_run(*, stack, script_path, params):
        calls.append((script_path, dict(params)))
        if script_path.endswith("hic_inspect.py"):
            return {
                "status": "success",
                "files": [str(mcool)],
                "formats": ["cooler"],
                "available_resolutions": [1_000, 2_000, 5_000, 10_000, 40_000, 100_000],
            }
        if script_path.endswith("hic_qc_and_balance.py"):
            return {
                "status": "success",
                "balanced_files": [str(mcool)],
                "ps_curve_tsv": str(tmp_path / "ps.tsv"),
                "per_file": [],
            }
        if script_path.endswith("hic_topology.py"):
            ana = params.get("analysis")
            if ana == "compartments":
                return {
                    "status": "success",
                    "ran": True,
                    "pct_A": 51.0,
                    "pct_B": 49.0,
                    "pc1_validated": False,
                    "compartments_tsv": str(tmp_path / "comp.tsv"),
                }
            if ana == "insulation_calibration":
                return {"status": "success", "boundary_strengths": {"200000": 0.5}}
            if ana == "tads":
                return {
                    "status": "success",
                    "ran": True,
                    "n_tads": 12,
                    "insulation_tsv": str(tmp_path / "ins.tsv"),
                    "tads_bed": str(tmp_path / "tads.bed"),
                }
            if ana == "loops":
                return {"status": "success", "ran": True, "n_loops": 5, "resolution": params["resolution"]}
        if script_path.endswith("hic_figure.py"):
            return {
                "status": "success",
                "figures": {"hic_contact_map": str(tmp_path / "map.png")},
                "figures_svg": {"hic_contact_map": str(tmp_path / "map.svg")},
                "figure_bundles": {
                    "hic_contact_map": {
                        "png": str(tmp_path / "map.png"),
                        "svg": str(tmp_path / "map.svg"),
                    },
                    "hic_contact_decay_ps": {
                        "png": str(tmp_path / "ps.png"),
                        "svg": str(tmp_path / "ps.svg"),
                    },
                },
            }
        return {"status": "success"}

    agent.env.run_in_stack.side_effect = _fake_run

    res = agent.run(
        "exp_microc",
        {
            "exp_context": {
                "genome": "hg38",
                "organism": "Homo sapiens",
                "modalities": {"HiC": [str(mcool)]},
                "output_dir": str(tmp_path / "work"),
            },
            "biological_intent": {"user_question": "Analyze TADs and loops"},
        },
    )
    assert res["status"] == "done"
    findings = res["findings"]

    # 1) Loop calling picked 5_000 bp (NOT 1_000 bp Micro-C resolution!)
    loop_calls = [p for s, p in calls if s.endswith("hic_topology.py") and p.get("analysis") == "loops"]
    assert len(loop_calls) == 1
    assert loop_calls[0]["resolution"] == 5_000

    # 2) hic_figure received pc1_validated=False and findings propagated figures_svg (including bundle fallback)
    fig_calls = [p for s, p in calls if s.endswith("hic_figure.py")]
    assert len(fig_calls) == 1
    assert fig_calls[0]["pc1_validated"] is False
    assert findings["figures_svg"]["hic_contact_map"].endswith("map.svg")
    assert findings["figures_svg"]["hic_contact_decay_ps"].endswith("ps.svg")

    # 3) RAM estimation is monotonic across Micro-C resolutions (never fallback 1.0 GB for 1kb/25kb!)
    dec = agent._advise_resolution(
        "exp_microc_ram",
        {"user_question": "Find fine chromatin loops"},
        {"available_resolutions": [1_000, 5_000, 25_000, 40_000]},
        {"organism": "Homo sapiens", "resolution": 1_000},
    )
    by_res = {c["resolution"]: c["ram_required_gb"] for c in dec["candidates"]}
    assert by_res[1_000] > by_res[5_000] > 10.0

    # 4) Unphased compartment finding never claims "A-compartment / B-compartment"
    published: list[dict] = []
    agent.publish_finding = lambda exp_id, payload, conf: published.append(payload)
    agent._publish_compartments_finding(
        "exp_microc",
        {"pct_A": 51.0, "pct_B": 49.0, "pc1_validated": False},
        {},
    )
    assert "unphased" in published[-1]["summary"].lower()
    assert "A-compartment" not in published[-1]["summary"]
