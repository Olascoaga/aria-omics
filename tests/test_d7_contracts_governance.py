"""Regression tests for D7 layer-contract fixes in governance / audit / manifest layers."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def test_d7_002_devils_advocate_reads_real_rna_qc_scrublet_ran_flag(tmp_path):
    """D7-002: rna_qc writes doublet detection status under qc['scrublet']={'ran': bool, ...}.
    devils_advocate must read qc['scrublet']['ran'] and never fall back to substring
    matching on JSON text (which falsely marked doublets addressed when ran=False)."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.narrative.devils_advocate import build_devils_advocate
    from aria.agents.narrative.hypothesis.adapters.scrna import scrna_evidence
    from aria.agents.narrative.hypothesis.caveats import DOUBLET
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.scripts.rna_qc import rna_qc

    rng = np.random.default_rng(42)
    n_cells, n_genes = 80, 400
    counts = rng.poisson(3, size=(n_cells, n_genes)).astype(np.float32)
    genes = [f"GENE{i}" for i in range(n_genes - 5)] + [f"MT-G{i}" for i in range(5)]
    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=genes),
    )
    in_h5ad = tmp_path / "raw.h5ad"
    adata.write_h5ad(str(in_h5ad))

    qc_res = rna_qc({
        "data_path": str(in_h5ad),
        "output_dir": str(tmp_path / "qc_out"),
        "run_scrublet": False,
        "min_genes": 10,
        "max_mt_pct": 50.0,
    })
    assert qc_res["status"] == "success", qc_res
    assert qc_res["scrublet"]["ran"] is False

    def _agent_results(qc_payload):
        return {
            "scrna_agent": {
                "status": "done",
                "findings": {
                    "qc": qc_payload,
                    "pseudobulk_de": {
                        "per_group": {
                            "T_cell": {
                                "per_comparison": {
                                    "treat_vs_ctrl": {
                                        "status": "success",
                                        "n_significant": 1,
                                        "top_genes": [{"gene": "GENE0", "log2fc": 2.1, "padj": 1e-4}],
                                    }
                                }
                            }
                        }
                    },
                },
            }
        }

    def _pb_block():
        return NarrativeBlock(
            id="scrna.pseudobulk.T_cell.treat_vs_ctrl",
            modality="scRNA",
            analysis="pseudobulk_de",
            block_type="result",
            title="T_cell treat vs ctrl",
            status="success",
            confidence="medium",
            claim="1 DE gene identified.",
            evidence=[EvidenceItem("DE genes", 1, source="rna_pseudobulk_de")],
            metadata={"claim": {"tier": "associative"}},
        )

    ar_false = _agent_results(qc_res)
    manifest_false = build_devils_advocate([_pb_block()], ar_false, exp_ctx={})
    dbl_false = next(c for c in manifest_false[0]["challenges"] if c["alternative"] == "doublets")
    sigs_false = scrna_evidence(ar_false, run_ledger=None, exp_ctx={})
    assert dbl_false["addressed"] is False, (
        f"Scrublet had ran=False ({qc_res['scrublet']}) yet devils_advocate reported {dbl_false}"
    )
    assert sigs_false and DOUBLET in sigs_false[0].caveats_inherited

    # Positive control on the same real rna_qc schema with scrublet['ran']=True:
    qc_ran = {**qc_res, "scrublet": {**qc_res["scrublet"], "ran": True, "n_doublets": 3}}
    ar_true = _agent_results(qc_ran)
    manifest_true = build_devils_advocate([_pb_block()], ar_true, exp_ctx={})
    dbl_true = next(c for c in manifest_true[0]["challenges"] if c["alternative"] == "doublets")
    sigs_true = scrna_evidence(ar_true, run_ledger=None, exp_ctx={})
    assert dbl_true["addressed"] is True
    assert sigs_true and DOUBLET not in sigs_true[0].caveats_inherited


def test_d7_008_falsification_certificates_reads_design_agent_technical_replicates():
    """D7-008: DesignAgent._build_replicate_handling(technical=True) writes
    mode='technical_aggregate' (never 'has_technical_replicates').
    falsification_certificates must resolve technical_replicates_aggregated=True."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.design_agent import DesignAgent
    from aria.agents.narrative.falsification_certificates import (
        build_falsification_certificates,
    )
    from aria.memory.memory import ARIAMemory

    da = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da._organism = "Homo sapiens"
    da._genome = "hg38"
    da._confirmed_groups = {
        "ctrl": ["cA_rep1", "cA_rep2", "cB_rep1", "cB_rep2", "cC_rep1", "cC_rep2"],
        "treat": ["tA_rep1", "tA_rep2", "tB_rep1", "tB_rep2", "tC_rep1", "tC_rep2"],
    }
    da._main_factor = "condition"
    da._replicate_handling = da._build_replicate_handling(technical=True)
    da._inferred_design = {"source": "sample_sheet"}
    design = da._build_design()
    assert design["replicate_handling"]["mode"] == "technical_aggregate"
    assert design["replicate_handling"]["aggregation"] == "sum_raw_counts"

    claims = [{
        "claim_id": "bulk.contrast.treat_vs_ctrl",
        "modality": "bulk_RNA",
        "analysis": "contrast",
        "tier": "associative",
    }]
    res = build_falsification_certificates(claims, agent_results={}, exp_ctx={"design": design})
    audit = res["certificates"][0]["design_aliasing_and_unit_audit"]
    assert audit["technical_replicates_aggregated"] is True

    # Independent biological design stays False:
    da._replicate_handling = da._build_replicate_handling(technical=False)
    design_ind = da._build_design()
    res_ind = build_falsification_certificates(claims, agent_results={}, exp_ctx={"design": design_ind})
    assert res_ind["certificates"][0]["design_aliasing_and_unit_audit"]["technical_replicates_aggregated"] is False


def test_d7_010_audit_agent_and_design_intelligence_detect_batch_on_design_agent_output(tmp_path):
    """D7-010: DesignAgent._build_design writes batch_covariate and
    paired_covariate_values (never batch_factor/batch_map). AuditAgent.run_audit
    must detect PC1 batch dominance and DesignIntelligence must detect batch
    confounding on real DesignAgent outputs."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.audit_agent import AuditAgent
    from aria.agents.design_agent import DesignAgent
    from aria.agents.design_intelligence import DesignIntelligence
    from aria.memory.memory import ARIAMemory

    rng = np.random.default_rng(23)
    samples = ["c1", "c2", "c3", "c4", "t1", "t2", "t3", "t4"]
    batch_for = {
        "c1": "B1", "c2": "B1", "c3": "B2", "c4": "B2",
        "t1": "B1", "t2": "B1", "t3": "B2", "t4": "B2",
    }
    base = rng.poisson(200, size=(80, 8)).astype(int)
    for j, s in enumerate(samples):
        if batch_for[s] == "B2":
            base[:40, j] *= 15
            base[40:, j] = np.maximum(1, base[40:, j] // 10)
    tsv = tmp_path / "counts.tsv"
    pd.DataFrame(base, index=[f"G{i}" for i in range(80)], columns=samples).to_csv(
        tsv, sep="\t", index_label="gene_id"
    )

    da = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da._organism = "Homo sapiens"
    da._genome = "hg38"
    da._confirmed_groups = {"ctrl": ["c1", "c2", "c3", "c4"], "treat": ["t1", "t2", "t3", "t4"]}
    da._main_factor = "condition"
    da._batch_covariate = "batch"
    da._paired_covariate_values = {"batch": batch_for}
    da._inferred_design = {"source": "sample_sheet"}
    design = da._build_design()

    auditor = AuditAgent(ARIAMemory(":memory:"), llm=None)
    out = auditor.run_audit(
        {"files": [str(tsv)], "modalities": {"bulk_RNA": [str(tsv)]}, "design": design},
        experiment_id="exp_audit",
    )
    checks = [f["check"] for f in out["findings"]]
    assert "pca_batch_dominance" in checks, (
        f"Strong batch PC1 dominance was missed on DesignAgent output; checks={checks}"
    )

    # DesignIntelligence also flags confounded batch from DesignAgent output:
    conf_batch = {
        "c1": "B1", "c2": "B1", "c3": "B1", "c4": "B1",
        "t1": "B2", "t2": "B2", "t3": "B2", "t4": "B2",
    }
    da._paired_covariate_values = {"batch": conf_batch}
    conf_design = da._build_design()
    di = DesignIntelligence()
    prof = di._bulk_rna_profile("bulk_RNA", {"design": conf_design}, intent={})
    assert any("confounded" in w.lower() for w in prof["warnings"]), prof["warnings"]


def test_d7_013_decision_manifest_preserves_design_formula_covariates_and_paired():
    """D7-013: decision_manifest._normalize_design must read DesignAgent's
    design_formula, batch_covariate, and replicate_handling.mode='paired_biological_block'
    instead of masking them to empty formula, empty covariates, and paired=False."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.design_agent import DesignAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.decision_manifest import build_decision_manifest, verify_decision_manifest

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
    assert design["design_formula"] == "~ donor + condition"
    assert design["replicate_handling"]["mode"] == "paired_biological_block"

    dm = build_decision_manifest(
        experiment_id="exp_dm",
        exp_context={"design": design},
        plan={"steps": []},
        biological_intent={},
    )
    assert verify_decision_manifest(dm) == (True, "verified")
    norm_d = dm["design"]
    assert norm_d["formula"] == "~ donor + condition"
    assert norm_d["paired"] is True
    assert "donor" in norm_d["covariates"]
    assert norm_d["donor_col"] == "donor"


def test_d7_004_devils_advocate_reads_design_agent_batch_covariate_and_formula():
    """D7-004: DesignAgent._build_design writes batch_covariate and design_formula
    (never top-level covariates). devils_advocate must resolve batch modeling via
    the real producer design contract and keep unmodeled designs flagged False."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.design_agent import DesignAgent
    from aria.agents.narrative.devils_advocate import build_devils_advocate
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.memory.memory import ARIAMemory

    def _bulk_block():
        return NarrativeBlock(
            id="bulk.contrast.treat_vs_ctrl",
            modality="bulk_RNA",
            analysis="contrast",
            block_type="result",
            title="treat vs ctrl",
            status="success",
            confidence="high",
            claim="42 DE genes identified.",
            evidence=[EvidenceItem("DE genes", 42, source="rna_bulk_de")],
            metadata={"claim": {"tier": "associative"}},
        )

    da = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da._organism = "Homo sapiens"
    da._genome = "hg38"
    da._confirmed_groups = {
        "ctrl": ["c1", "c2", "c3", "c4"],
        "treat": ["t1", "t2", "t3", "t4"],
    }
    da._main_factor = "condition"
    da._batch_covariate = "batch"
    da._inferred_design = {"source": "sample_sheet"}
    design_with_batch = da._build_design()
    assert design_with_batch["batch_covariate"] == "batch"
    assert design_with_batch["design_formula"] == "~ batch + condition"

    manifest = build_devils_advocate(
        [_bulk_block()], agent_results={}, exp_ctx={"design": design_with_batch}
    )
    batch_ch = next(c for c in manifest[0]["challenges"] if c["alternative"] == "batch_effect")
    assert batch_ch["addressed"] is True, (
        f"DesignAgent modeled {design_with_batch['design_formula']!r}, yet devils_advocate reported {batch_ch}"
    )

    # Non-'batch'-named blocking covariate confirmed by DesignAgent (e.g. 'donor'):
    da._batch_covariate = "donor"
    design_with_donor = da._build_design()
    assert design_with_donor["design_formula"] == "~ donor + condition"
    manifest_donor = build_devils_advocate(
        [_bulk_block()], agent_results={}, exp_ctx={"design": design_with_donor}
    )
    donor_ch = next(c for c in manifest_donor[0]["challenges"] if c["alternative"] == "batch_effect")
    assert donor_ch["addressed"] is True

    # Negative control: DesignAgent with no batch covariate stays addressed=False.
    da._batch_covariate = None
    design_no_batch = da._build_design()
    assert design_no_batch["design_formula"] == "~ condition"
    manifest_none = build_devils_advocate(
        [_bulk_block()], agent_results={}, exp_ctx={"design": design_no_batch}
    )
    none_ch = next(c for c in manifest_none[0]["challenges"] if c["alternative"] == "batch_effect")
    assert none_ch["addressed"] is False

