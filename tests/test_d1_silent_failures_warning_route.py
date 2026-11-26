"""Permanent regression tests for D1 Track B findings: D1-005 (part 2) and D1-025."""
from __future__ import annotations

import sys
import types

import pandas as pd
import pytest


def test_d1_005_part2_scrna_methods_must_not_claim_apeglm_when_shrinkage_failed():
    """D1-005 (part 2): _lfc_shrinkage_clause and ScrnaNarrator.methods() must NOT claim
    apeGLM-shrunken LFCs when shrinkage failed on every comparison, and the block must be warned."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.scrna._common import _lfc_shrinkage_clause

    pb = {
        "status": "success",
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "lfc_shrinkage": {
            "requested": True,
            "applied": False,
            "n_applied": 0,
            "n_failed": 1,
            "failed_reasons": ["treat_vs_ctrl: lfc_shrink failed: solver error"],
        },
        "per_group": {
            "Tcell": {
                "per_comparison": {
                    "treat_vs_ctrl": {
                        "status": "success",
                        "n_significant": 5,
                        "n_up": 3,
                        "n_down": 2,
                        "n_replicates": {"test": 3, "ref": 3},
                        "lfc_shrinkage": {
                            "applied": False,
                            "method": None,
                            "reason": "lfc_shrink failed: solver error",
                        },
                        "top_genes": [{"gene": "G1", "log2fc": 2.1, "padj": 1e-4}],
                    }
                }
            }
        },
    }

    clause = _lfc_shrinkage_clause(pb)
    assert "are apeGLM-shrunken estimates" not in clause
    assert "unshrunken MLE" in clause

    agent_result = {"findings": {"pseudobulk_de": pb}}
    methods_text = " ".join(ScrnaNarrator().methods("scrna_agent", agent_result))
    assert "are apeGLM-shrunken estimates" not in methods_text
    assert "unshrunken MLE" in methods_text

    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    pb_block = next(b for b in blocks if b.id.startswith("scrna.pseudobulk."))
    assert pb_block.status == "warnings"
    assert any("shrink" in c.text.lower() and c.severity == "warning" for c in pb_block.caveats)


def test_d1_005_part2_pseudobulk_script_records_realized_shrink_counts(tmp_path, monkeypatch):
    """D1-005 (part 2 producer): rna_pseudobulk_de must record applied=False, n_applied=0,
    n_failed>0 when DeseqStats.lfc_shrink raises."""
    pytest.importorskip("pydeseq2")
    ad = pytest.importorskip("anndata")
    import numpy as np
    from pydeseq2.ds import DeseqStats
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    rng = np.random.default_rng(42)
    donors = ["d1", "d2", "d3", "d4", "d5", "d6"]
    conds = ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]
    X_rows, obs_rows = [], []
    for d, c in zip(donors, conds):
        base = rng.poisson(30, size=(12, 20))
        if c == "treat":
            base[:, :5] = rng.poisson(120, size=(12, 5))
        X_rows.append(base)
        obs_rows.extend({"cell_type": "Tcell", "condition": c, "donor": d} for _ in range(12))
    h5 = tmp_path / "pb.h5ad"
    ad.AnnData(
        X=np.vstack(X_rows).astype(np.int32),
        obs=pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(72)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(20)]),
    ).write_h5ad(h5)

    def _fail_shrink(self, coeff=None, **kw):
        raise RuntimeError("apeGLM solver failed")

    monkeypatch.setattr(DeseqStats, "lfc_shrink", _fail_shrink)
    res = rna_pseudobulk_de({
        "data_path": str(h5),
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["treat", "ctrl"]],
        "min_cells_per_pseudosample": 5,
        "min_replicates_per_condition": 2,
        "output_dir": str(tmp_path / "out"),
        "lfc_shrink": True,
    })
    sk = res["lfc_shrinkage"]
    assert sk["requested"] is True
    assert sk["applied"] is False
    assert sk["n_applied"] == 0
    assert sk["n_failed"] == 1
    assert any("apeGLM solver failed" in r for r in sk["failed_reasons"])


def test_d1_025_enrichr_partial_error_must_reach_narrator_caveats(monkeypatch):
    """D1-025 (partial outage): KEGG 503 while GO_BP succeeds must reach bulk.pathway caveats
    and downgrade status to 'warnings'."""
    pytest.importorskip("litellm")
    from aria.scripts.rna_bulk.ora import _enrichr_enrichment
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    fake_gp = types.ModuleType("gseapy")

    def _fake_enrichr(**kw):
        if "KEGG" in kw["gene_sets"]:
            raise RuntimeError("503 Service Unavailable")
        return types.SimpleNamespace(
            results=pd.DataFrame([{
                "Term": "Cell cycle",
                "Adjusted P-value": 0.001,
                "Overlap": "5/100",
                "Genes": "G1;G2",
            }])
        )

    fake_gp.enrichr = _fake_enrichr
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)
    monkeypatch.setattr("time.sleep", lambda *a, **kw: None)

    pathways, warns = _enrichr_enrichment(
        ["G1", "G2", "G3"],
        ["GO_BP", "KEGG"],
        {"GO_BP": "GO_Biological_Process_2023", "KEGG": "KEGG_2021_Human"},
        "Homo sapiens",
        ["G1", "G2", "G3", "G4"],
        False,
        ["G1", "G2", "G3"],
    )
    blocks = BulkRnaNarrator().collect("bulk_rna_agent", {
        "findings": {
            "warnings": warns,
            "contrasts": [{
                "name": "t_vs_c",
                "status": "success",
                "n_significant": 3,
                "pathways": pathways,
            }],
        }
    })
    pw_block = next(b for b in blocks if b.id.startswith("bulk.pathway."))
    assert pw_block.status == "warnings"
    assert any("KEGG" in c.text and c.severity == "warning" for c in pw_block.caveats)


def test_d1_025_enrichr_total_outage_synthesizes_degradation_block(monkeypatch):
    """D1-025 (total outage): when ALL Enrichr DBs fail and pathways=={}, router synthesizes
    a visible degradation block instead of silently emitting zero pathway blocks."""
    pytest.importorskip("litellm")
    from aria.scripts.rna_bulk.ora import _enrichr_enrichment
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    fake_gp = types.ModuleType("gseapy")
    fake_gp.enrichr = lambda **kw: (_ for _ in ()).throw(RuntimeError("503 Service Unavailable"))
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)
    monkeypatch.setattr("time.sleep", lambda *a, **kw: None)

    failed_dbs: dict[str, str] = {}
    pathways, warns = _enrichr_enrichment(
        ["G1", "G2", "G3"],
        ["GO_BP", "KEGG"],
        {"GO_BP": "GO_Biological_Process_2023", "KEGG": "KEGG_2021_Human"},
        "Homo sapiens",
        ["G1", "G2", "G3", "G4"],
        False,
        ["G1", "G2", "G3"],
        failed_databases=failed_dbs,
    )
    assert pathways == {}
    assert set(failed_dbs) == {"GO_BP", "KEGG"}

    blocks = BulkRnaNarrator().collect("bulk_rna_agent", {
        "findings": {
            "warnings": warns,
            "contrasts": [{
                "name": "t_vs_c",
                "status": "success",
                "n_significant": 3,
                "pathways": pathways,
                "pathway_ora": {"method": "enrichr", "failed_databases": failed_dbs},
            }],
        }
    })
    deg_blocks = [b for b in blocks if b.analysis == "pathway_enrichment"]
    assert deg_blocks, f"Expected synthesized pathway degradation block, got {[b.id for b in blocks]}"
    assert deg_blocks[0].status == "warnings"
    assert any("KEGG" in c.text for c in deg_blocks[0].caveats)
