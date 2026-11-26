"""Regression tests for D7 layer-contract fixes in governance / audit / manifest layers."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def test_d7_002_devils_advocate_reads_real_rna_qc_scrublet_ran_flag(tmp_path):
    """D7-002: rna_qc writes doublet detection status under qc['scrublet']={'ran': bool, ...}.
    devils_advocate must read qc['scrublet']['ran'] and never fall back to substring
    matching on JSON text (which falsely marked doublets addressed when ran=False)."""
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
