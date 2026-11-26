"""Contract tests for centralized degradation/warning routing (Track B / W-WARN)."""
from __future__ import annotations

import pytest

pytest.importorskip("litellm")


def test_warning_router_contract_fields_and_status_transitions():
    """Every contract field (warnings, *_degraded, failed_databases, error, figures._errors)
    routes to caveats with exact status/confidence transitions and idempotency."""
    from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock
    from aria.agents.narrative.warning_router import route_degradation_signals

    # 1. Warning signal downgrades status success -> warnings and caps confidence high -> medium.
    blk_qc = NarrativeBlock(
        id="bulk.qc",
        modality="bulk RNA-seq",
        analysis="sample_qc",
        block_type="qc",
        title="Bulk RNA sample QC",
        status="success",
        confidence="high",
        claim="Bulk RNA sample QC evaluated 6 samples with 0 outlier(s) flagged.",
        evidence=[EvidenceItem("samples", 6, "sample_qc")],
    )
    blk_pw = NarrativeBlock(
        id="bulk.pathway.t_vs_c",
        modality="bulk RNA-seq",
        analysis="pathway_enrichment",
        block_type="result",
        title="Pathway enrichment t_vs_c",
        status="success",
        confidence="medium",
        claim="Pathway enrichment found 2 term(s) for t_vs_c.",
        evidence=[EvidenceItem("enriched terms", 2, "bulk_pathways")],
    )
    agent_result = {
        "findings": {
            "warnings": [
                "Gene annotation loaded from gencode.gtf: 19,000 biotypes, 19,000 lengths.",
                "pydeseq2 VST unavailable — using log2(size-factor normalized + 1).",
            ],
            "sample_qc": {"n_samples": 6, "transform_degraded": True},
            "contrasts": [{
                "name": "t_vs_c",
                "status": "success",
                "pathway_ora": {"failed_databases": {"KEGG": "503 Service Unavailable"}},
            }],
        }
    }

    routed = route_degradation_signals("bulk_rna_agent", agent_result, [blk_qc, blk_pw])
    qc_out = next(b for b in routed if b.id == "bulk.qc")
    pw_out = next(b for b in routed if b.id == "bulk.pathway.t_vs_c")

    assert qc_out.status == "warnings"
    assert qc_out.confidence == "medium"
    assert any("VST unavailable" in c.text and c.severity == "warning" for c in qc_out.caveats)
    assert any("transform_degraded" in c.text for c in qc_out.caveats)
    # Informational telemetry line stays info (does not masquerade as warning).
    assert any("Gene annotation loaded" in c.text and c.severity == "info" for c in qc_out.caveats)

    assert pw_out.status == "warnings"
    assert any("KEGG" in c.text and "503" in c.text and c.severity == "warning" for c in pw_out.caveats)

    # Idempotency: second pass adds zero duplicate caveats and does not downgrade confidence again.
    n_qc_cavs = len(qc_out.caveats)
    n_pw_cavs = len(pw_out.caveats)
    routed2 = route_degradation_signals("bulk_rna_agent", agent_result, routed)
    qc_out2 = next(b for b in routed2 if b.id == "bulk.qc")
    pw_out2 = next(b for b in routed2 if b.id == "bulk.pathway.t_vs_c")
    assert len(qc_out2.caveats) == n_qc_cavs
    assert len(pw_out2.caveats) == n_pw_cavs
    assert qc_out2.confidence == "medium"


def test_info_only_telemetry_does_not_downgrade_clean_block():
    """Routine informational messages in findings['warnings'] stay info and keep status='success'."""
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.agents.narrative.warning_router import route_degradation_signals

    blk = NarrativeBlock(
        id="bulk.qc",
        modality="bulk RNA-seq",
        analysis="sample_qc",
        block_type="qc",
        title="Bulk RNA sample QC",
        status="success",
        confidence="high",
        claim="Bulk RNA sample QC evaluated 6 samples with 0 outlier(s) flagged.",
        evidence=[EvidenceItem("samples", 6, "sample_qc")],
    )
    res = {
        "findings": {
            "warnings": [
                "Gene annotation loaded from gencode.gtf: 19,000 biotypes, 19,000 lengths.",
                "DR input: 2000 most-variable protein_coding genes (from 19000 total protein_coding).",
                "TPM table written (19000 genes) — supplementary only; not used for DE or PCA.",
                "120 low-count genes removed (< 10 counts in < 2 samples).",
            ]
        }
    }
    out = route_degradation_signals("bulk_rna_agent", res, [blk])
    assert out[0].status == "success"
    assert out[0].confidence == "high"
    assert all(c.severity == "info" for c in out[0].caveats)


def test_blocking_stage_error_fails_closed_even_if_narrator_emitted_success():
    """A primary stage error converts a mistakenly-successful block to status='error' with claim=''."""
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.agents.narrative.warning_router import route_degradation_signals

    blk = NarrativeBlock(
        id="multiome.mofa",
        modality="integration",
        analysis="mofa",
        block_type="result",
        title="MOFA+",
        status="success",
        confidence="high",
        claim="MOFA+ identified 3 factors.",
        evidence=[EvidenceItem("fitted factors", 3, "mofa")],
    )
    res = {"findings": {"mofa": {"error": "MOFA solver diverged"}}}
    out = route_degradation_signals("integration_agent", res, [blk])
    assert out[0].status == "error"
    assert out[0].block_type == "error"
    assert out[0].confidence == "insufficient"
    assert out[0].claim == ""
    assert "MOFA solver diverged" in (out[0].error or "")


def test_warnings_block_survives_full_publication_pipeline():
    """A block with status='warnings' keeps associative tier, prose, W-CLAIM, devil's advocate,
    W-LEDGER linkage, methodology limitations, and amber HTML badge."""
    from aria.agents.narrative.claim_compiler import compile_public_claims
    from aria.agents.narrative.compose_prose import compose_block_prose
    from aria.agents.narrative.devils_advocate import build_devils_advocate
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.run_ledger import build_run_ledger, ensure_report_ledger_nodes
    from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock, PUBLISHABLE_STATUSES
    from aria.agents.narrative.validators import validate_block

    assert "warnings" in PUBLISHABLE_STATUSES
    blk = NarrativeBlock(
        id="bulk.contrast.t_vs_c",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        block_type="result",
        title="Bulk contrast t_vs_c",
        status="warnings",
        confidence="medium",
        claim="Contrast t_vs_c identified 42 significant DE genes at padj < 0.05.",
        evidence=[
            EvidenceItem("significant genes", 42, "deseq2"),
            EvidenceItem("genes tested", 12000, "deseq2"),
            EvidenceItem("replicates", "4 vs 4", "deseq2"),
        ],
        caveats=[Caveat("Full DE table export failed: disk full.", severity="warning")],
        metrics={
            "n_significant": 42,
            "n_genes_tested": 12000,
            "n_replicates": {"test": 4, "ref": 4},
        },
    )
    validate_block(blk)
    prose = compose_block_prose(blk)
    assert "did not produce an interpretable result" not in prose
    assert "42" in prose

    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "contrasts": [{
                    "name": "t_vs_c",
                    "status": "success",
                    "n_significant": 42,
                    "n_genes_tested": 12000,
                }]
            },
        }
    }
    ledger = build_run_ledger({}, agent_results)
    ensure_report_ledger_nodes(ledger, [blk])
    comp = compile_public_claims([blk], {}, run_ledger=ledger)
    assert len(comp.claims) == 1, comp.withheld
    claim = comp.claims[0]
    assert claim["status"] == "warnings"
    assert claim["tier"] == "associative"
    assert any("Full DE table export failed" in lim for lim in claim["limitations"])

    da = build_devils_advocate(comp.blocks, agent_results, {})
    assert len(da) == 1 and da[0]["claim_id"] == "bulk.contrast.t_vs_c"

    rendered = render_blocks(comp.blocks, strict=True)
    assert "warnings / medium" in rendered
    assert "Full DE table export failed" in rendered


def test_cross_modal_parity_all_five_narrators_route_warnings():
    """All 5 modality narrators surface stage warnings via the single router in NarrativeRegistry."""
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.narrators.hic import HiCNarrator
    from aria.agents.narrative.narrators.integration import IntegrationNarrator
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.registry import registry_with

    reg = registry_with([
        BulkRnaNarrator(),
        ScrnaNarrator(),
        ChromatinNarrator(),
        HiCNarrator(),
        IntegrationNarrator(),
    ])
    agent_results = {
        "bulk_rna_agent": {
            "findings": {
                "warnings": ["de_all_genes.tsv export failed: permission denied"],
                "sample_qc": {"n_samples": 6, "outliers": []},
            }
        },
        "scrna_agent": {
            "findings": {
                "qc": {"status": "success", "n_cells_before": 1000, "n_cells_after": 900,
                       "warnings": ["ambient RNA estimation degraded"]},
            }
        },
        "chromatin_agent": {
            "findings": {
                "scATAC": {
                    "findings": {
                        "qc": {
                            "status": "success",
                            "data_type": "scATAC",
                            "n_cells": 500,
                            "n_fragments": 100000,
                            "warnings": ["Fragment scan capped at 100000 records"],
                        },
                        "figures": {"_errors": ["marker_heatmap failed: OOM"]},
                    }
                }
            }
        },
        "genome_arch_agent": {
            "findings": {
                "compartments": {
                    "status": "success",
                    "ran": True,
                    "pct_A": 51.0,
                    "pct_B": 49.0,
                    "resolution": 100000,
                    "pc1_validated": True,
                    "phasing_method": "gc_content",
                    "warnings": ["partial_unphased_chroms:chrY"],
                }
            }
        },
        "integration_agent": {
            "findings": {
                "mofa": {
                    "status": "success",
                    "ran": True,
                    "n_cells": 200,
                    "n_factors": 3,
                    "n_active_factors": 2,
                    "warnings": ["Seurat HVG selection fell back to variance ranking"],
                    "metrics_not_computed": ["cell_cycle_scores (parse failed)"],
                }
            }
        },
    }

    blocks = reg.collect_blocks(agent_results, {})
    by_id = {b.id: b for b in blocks}

    assert by_id["bulk.qc"].status == "warnings"
    assert any("de_all_genes.tsv" in c.text for c in by_id["bulk.qc"].caveats)

    assert by_id["scrna.qc"].status == "warnings"
    assert any("ambient RNA estimation degraded" in c.text for c in by_id["scrna.qc"].caveats)

    assert by_id["chromatin.qc.scATAC"].status == "warnings"
    assert any("Fragment scan capped" in c.text for c in by_id["chromatin.qc.scATAC"].caveats)
    assert any("marker_heatmap failed" in c.text for b in blocks for c in b.caveats)

    assert by_id["hic.compartments"].status == "warnings"
    assert any("partial_unphased_chroms:chrY" in c.text for c in by_id["hic.compartments"].caveats)

    assert by_id["multiome.mofa"].status == "warnings"
    assert any("variance ranking" in c.text for c in by_id["multiome.mofa"].caveats)
    assert any("cell_cycle_scores" in c.text for c in by_id["multiome.mofa"].caveats)
