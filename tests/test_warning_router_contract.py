"""Contract tests for centralized degradation & warning routing (Track B / W-WARN).

Covers D1_fix_design_warning_route.md rev. 2 + lead amendments A2 & A3:
  - purity (stdlib + narrative.types only, zero I/O, zero ``except Exception``,
    <=300 lines, single choke point in NarrativeRegistry.collect_blocks, no
    decorator);
  - structured signal routing (failed-stage blocking error, A3 success/absent-
    status error -> warning, *_error, figures._errors, failed_databases /
    failed_samples / failed_modalities, *_degraded, phasing_fallback_reason,
    lfc_shrinkage mismatch, metrics_not_computed, orphan limitation/error block
    synthesis, deduplication & idempotency);
  - free-text warnings structural-position-only routing (A2: severity="warning",
    zero keyword routing, never changes confidence/status, never synthesizes
    blocks, root warnings attach only to QC/summary block or fall back to
    first block metadata["degradation_signals"] with severity="info");
  - unchanged status contract (B2: no new status value, passes validate_blocks,
    compile_public_claims, and render_blocks).
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.narrators import (
    BulkRnaNarrator,
    IntegrationNarrator,
)
from aria.agents.narrative.registry import registry_with
from aria.agents.narrative.render_blocks import render_blocks
from aria.agents.narrative.run_ledger import build_run_ledger
from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock
from aria.agents.narrative.validators import validate_blocks

_REPO_ROOT = Path(__file__).resolve().parents[1]
_ROUTER_PATH = _REPO_ROOT / "aria" / "agents" / "narrative" / "warning_router.py"


def _sample_block(
    block_id: str,
    *,
    modality: str = "integration",
    analysis: str = "mofa",
    block_type: str = "result",
    status: str = "success",
    confidence: str = "high",
    claim: str = "MOFA+ retained 4 active biological factors across 600 cells.",
) -> NarrativeBlock:
    return NarrativeBlock(
        id=block_id,
        modality=modality,
        analysis=analysis,
        block_type=block_type,
        title=f"Title for {block_id}",
        status=status,
        confidence=confidence,
        claim=claim,
        evidence=[
            EvidenceItem("active biological factors", 4, source="mofa"),
            EvidenceItem("cells", 600, source="mofa"),
        ],
    )


# ── 1. Purity, size & single-choke-point architecture ────────────────────────

def test_warning_router_is_pure_small_and_has_no_broad_except():
    assert _ROUTER_PATH.is_file(), f"Missing {_ROUTER_PATH}"
    source = _ROUTER_PATH.read_text(encoding="utf-8")
    lines = source.splitlines()
    assert len(lines) <= 300, f"warning_router.py has {len(lines)} lines (target <=300)"

    tree = ast.parse(source, filename=str(_ROUTER_PATH))
    allowed_roots = set(sys.stdlib_module_names) | {"__future__"}
    disallowed_io = {"os", "io", "pathlib", "shutil", "subprocess", "socket", "urllib", "http"}

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".", 1)[0]
                assert root in allowed_roots and root not in disallowed_io, (
                    f"Disallowed import in warning_router.py: {alias.name}"
                )
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod != "aria.agents.narrative.types":
                root = mod.split(".", 1)[0]
                assert node.level == 0 and root in allowed_roots and root not in disallowed_io, (
                    f"Disallowed from-import in warning_router.py: {mod}"
                )
        elif isinstance(node, ast.ExceptHandler):
            if node.type is None:
                raise AssertionError("Bare except forbidden in warning_router.py")
            if isinstance(node.type, ast.Name) and node.type.id in {"Exception", "BaseException"}:
                raise AssertionError("except Exception forbidden in warning_router.py")
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            assert node.func.id != "open", "File I/O open() forbidden in warning_router.py"

    # Single choke point: no @routed_collect decorator anywhere under narrative/.
    narrative_dir = _REPO_ROOT / "aria" / "agents" / "narrative"
    for p in narrative_dir.rglob("*.py"):
        if "__pycache__" in p.parts:
            continue
        assert "routed_collect" not in p.read_text(encoding="utf-8"), (
            f"Forbidden decorator routed_collect found in {p}"
        )


# ── 2. Structured signals & Amendment A3 ─────────────────────────────────────

def test_structured_blocking_requires_failed_stage_a3():
    from aria.agents.narrative.warning_router import route_degradation_signals

    # (a) Stage with status="error" + error key -> blocking error block.
    blk_fail = _sample_block("multiome.mofa", analysis="mofa", confidence="high")
    out_fail = route_degradation_signals(
        "integration_agent",
        {"findings": {"mofa": {"status": "error", "error": "MOFA convergence failed"}}},
        [blk_fail],
    )
    assert len(out_fail) == 1
    assert out_fail[0].status == "error"
    assert out_fail[0].block_type == "error"
    assert out_fail[0].claim == ""
    assert out_fail[0].confidence == "insufficient"
    assert out_fail[0].error == "MOFA convergence failed"
    assert any(c.severity == "blocking" and "MOFA convergence failed" in c.text for c in out_fail[0].caveats)

    # (b) A3 negative test: error key on a stage with status="success" -> warning, NOT blocking!
    blk_ok = _sample_block("multiome.mofa", analysis="mofa", confidence="high")
    out_ok = route_degradation_signals(
        "integration_agent",
        {"findings": {"mofa": {"status": "success", "error": "partial view 2 read error"}}},
        [blk_ok],
    )
    assert out_ok[0].status == "success"
    assert out_ok[0].block_type == "result"
    assert out_ok[0].claim != ""
    assert out_ok[0].confidence == "medium"  # capped high -> medium
    assert any(c.severity == "warning" and "partial view 2 read error" in c.text for c in out_ok[0].caveats)
    assert not any(c.severity == "blocking" for c in out_ok[0].caveats)

    # (c) A3 negative test: error key with status absent -> warning, NOT blocking!
    blk_absent = _sample_block("hic.qc", modality="hic", analysis="qc", confidence="high")
    out_absent = route_degradation_signals(
        "genome_arch_agent",
        {"findings": {"inspection": {"error": "sample2.mcool header warning"}}},
        [blk_absent],
    )
    assert out_absent[0].status == "success"
    assert out_absent[0].confidence == "medium"
    assert any(c.severity == "warning" and "sample2.mcool" in c.text for c in out_absent[0].caveats)
    assert not any(c.severity == "blocking" for c in out_absent[0].caveats)


def test_structured_auxiliary_errors_degraded_flags_shrinkage_and_mnc():
    from aria.agents.narrative.warning_router import route_degradation_signals

    # *_error + figures["_errors"] -> structured warning + confidence cap.
    blk_qc = _sample_block("chromatin.qc.bulk_ATAC", modality="chromatin", analysis="qc", block_type="qc", confidence="high")
    out_chrom = route_degradation_signals(
        "chromatin_agent",
        {
            "findings": {
                "qc": {"status": "success", "pca_error": "SVD did not converge"},
                "figures": {"_errors": ["heatmap: empty matrix"]},
            }
        },
        [blk_qc],
    )
    assert out_chrom[0].status == "success"
    assert out_chrom[0].confidence == "medium"
    texts = [(c.severity, c.text) for c in out_chrom[0].caveats]
    assert any(sev == "warning" and "SVD did not converge" in t for sev, t in texts)
    assert any(sev == "warning" and "heatmap: empty matrix" in t for sev, t in texts)

    # *_degraded + phasing_fallback_reason -> structured warning + confidence cap.
    blk_comp = _sample_block("hic.compartments", modality="hic", analysis="compartments", confidence="high")
    out_hic = route_degradation_signals(
        "genome_arch_agent",
        {
            "findings": {
                "topology": {
                    "compartments": {
                        "status": "success",
                        "phasing_fallback_reason": "FASTA open failed; fell back to GTF phasing",
                        "warnings": ["FASTA open failed; fell back to GTF phasing"],
                    }
                }
            }
        },
        [blk_comp],
    )
    assert out_hic[0].status == "success"
    assert out_hic[0].confidence == "medium"
    # Identical text in phasing_fallback_reason + warnings is deduplicated to 1 caveat!
    phasing_cavs = [c for c in out_hic[0].caveats if "FASTA open failed" in c.text]
    assert len(phasing_cavs) == 1
    assert phasing_cavs[0].severity == "warning"

    # background_degraded + lfc_shrinkage mismatch -> structured warnings + confidence cap.
    blk_pb = _sample_block("scrna.pseudobulk.T.A_vs_B", modality="scRNA-seq", analysis="pseudobulk_de", confidence="high")
    out_sc = route_degradation_signals(
        "scrna_agent",
        {
            "findings": {
                "pseudobulk_de": {
                    "background_degraded": True,
                    "background_source": "expressed_fallback",
                    "lfc_shrinkage": {
                        "requested": True,
                        "applied": False,
                        "n_applied": 0,
                        "n_failed": 2,
                        "failed_reasons": ["singular Hessian"],
                    },
                }
            }
        },
        [blk_pb],
    )
    assert out_sc[0].status == "success"
    assert out_sc[0].confidence == "medium"
    assert any(c.severity == "warning" and "background degraded" in c.text.lower() for c in out_sc[0].caveats)
    assert any(c.severity == "warning" and "shrinkage" in c.text.lower() for c in out_sc[0].caveats)

    # metrics_not_computed -> severity="info", NEVER caps confidence, NEVER synthesizes orphan blocks.
    blk_mofa = _sample_block("multiome.mofa", analysis="mofa", confidence="high")
    out_mnc = route_degradation_signals(
        "integration_agent",
        {
            "findings": {
                "mofa": {"status": "success", "metrics_not_computed": ["batch_factor_association"]},
                "wnn": {"status": "success", "metrics_not_computed": ["atac_frip"]},
            }
        },
        [blk_mofa],
    )
    assert len(out_mnc) == 1  # wnn had 0 blocks and metrics_not_computed must NOT synthesize an orphan
    assert out_mnc[0].confidence == "high"
    assert any(c.severity == "info" and "batch_factor_association" in c.text for c in out_mnc[0].caveats)


def test_orphan_block_synthesis_only_for_structured_warning_or_blocking():
    from aria.agents.narrative.warning_router import route_degradation_signals

    # Structured warning on orphan stage -> limitation block with claim="", confidence="insufficient".
    out_warn = route_degradation_signals(
        "bulk_rna_agent",
        {
            "findings": {
                "contrasts": [
                    {
                        "name": "A_vs_B",
                        "pathway_ora": {"failed_databases": {"KEGG": "HTTP 503"}},
                    }
                ]
            }
        },
        [],
    )
    assert len(out_warn) == 1
    assert out_warn[0].block_type == "limitation"
    assert out_warn[0].status == "limitation"
    assert out_warn[0].claim == ""
    assert out_warn[0].confidence == "insufficient"
    assert any(c.severity == "warning" and "KEGG" in c.text for c in out_warn[0].caveats)
    validate_blocks(out_warn)

    # Structured blocking error on orphan stage -> error block with claim="", confidence="insufficient".
    out_err = route_degradation_signals(
        "integration_agent",
        {"findings": {"mofa": {"status": "error", "error": "rank-deficient matrix"}}},
        [],
    )
    assert len(out_err) == 1
    assert out_err[0].block_type == "error"
    assert out_err[0].status == "error"
    assert out_err[0].claim == ""
    assert out_err[0].confidence == "insufficient"
    assert out_err[0].error == "rank-deficient matrix"
    validate_blocks(out_err)

    # Auxiliary stage error (plots) on bulk contrast must NOT turn primary contrast into error.
    blk_c = _sample_block(
        "bulk.contrast.A_vs_B",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        confidence="high",
    )
    out_aux = route_degradation_signals(
        "bulk_rna_agent",
        {
            "findings": {
                "contrasts": [
                    {
                        "name": "A_vs_B",
                        "status": "success",
                        "plots": {"status": "error", "error": "matplotlib backend error"},
                    }
                ]
            }
        },
        [blk_c],
    )
    assert out_aux[0].status == "success"
    assert out_aux[0].confidence == "medium"
    assert any(c.severity == "warning" and "matplotlib backend error" in c.text for c in out_aux[0].caveats)


# ── 3. Free-text warnings & Amendments A2 / B3 negative tests ────────────────

def test_freetext_stage_warnings_keep_warning_severity_never_change_confidence_or_synthesize():
    from aria.agents.narrative.warning_router import route_degradation_signals

    blk_mofa = _sample_block("multiome.mofa", analysis="mofa", confidence="high")
    blk_wnn = _sample_block("multiome.wnn", analysis="wnn", confidence="high")

    out = route_degradation_signals(
        "integration_agent",
        {
            "findings": {
                "mofa": {
                    "status": "success",
                    "warnings": ["<b>Subsampled</b> to 1000 cells\nfor training."],
                },
                "eregulons": {
                    "status": "success",
                    "warnings": ["Orphan warning that must NOT synthesize a block"],
                },
            }
        },
        [blk_mofa, blk_wnn],
    )
    # No orphan block synthesized for eregulons!
    assert [b.id for b in out] == ["multiome.mofa", "multiome.wnn"]
    # mofa gets sanitized severity="warning" caveat (A2), confidence stays "high", status stays "success".
    mofa_blk = out[0]
    assert mofa_blk.status == "success"
    assert mofa_blk.confidence == "high"
    assert any(
        c.severity == "warning" and c.text == "Subsampled to 1000 cells for training."
        for c in mofa_blk.caveats
    )
    # wnn untouched (structural position only).
    assert out[1].caveats == []
    assert out[1].confidence == "high"


def test_root_warnings_no_keyword_routing_and_qc_or_metadata_fallback_a2():
    from aria.agents.narrative.warning_router import route_degradation_signals

    # (a) With QC block present: root warnings mentioning Enrichr/GTF/VST/PCA/TPM/outlier
    # must NOT be keyword-routed to pathway/contrast blocks; attach ONLY to bulk.qc with
    # severity="warning" and WITHOUT changing bulk.qc.confidence ("high" stays "high").
    qc_blk = _sample_block(
        "bulk.qc",
        modality="bulk RNA-seq",
        analysis="sample_qc",
        block_type="qc",
        confidence="high",
    )
    contrast_blk = _sample_block(
        "bulk.contrast.A_vs_B",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        confidence="high",
    )
    pw_blk = _sample_block(
        "bulk.pathway.A_vs_B",
        modality="bulk RNA-seq",
        analysis="pathway_enrichment",
        confidence="high",
    )
    root_warn = "Enrichr GTF VST PCA TPM outlier diagnostic note"
    out = route_degradation_signals(
        "bulk_rna_agent",
        {"findings": {"warnings": [root_warn]}},
        [qc_blk, contrast_blk, pw_blk],
    )
    assert out[0].confidence == "high"
    assert out[0].status == "success"
    assert any(c.severity == "warning" and c.text == root_warn for c in out[0].caveats)
    assert out[1].caveats == []
    assert out[2].caveats == []

    # (b) A2 fallback when NO QC/summary block exists: attach NOTHING to caveats;
    # record in metadata["degradation_signals"] of the first block with severity="info".
    c_only = _sample_block(
        "bulk.contrast.A_vs_B",
        modality="bulk RNA-seq",
        analysis="differential_expression",
        confidence="high",
    )
    p_only = _sample_block(
        "bulk.pathway.A_vs_B",
        modality="bulk RNA-seq",
        analysis="pathway_enrichment",
        confidence="high",
    )
    out_no_qc = route_degradation_signals(
        "bulk_rna_agent",
        {"findings": {"warnings": [root_warn]}},
        [c_only, p_only],
    )
    assert len(out_no_qc) == 2
    assert out_no_qc[0].caveats == []
    assert out_no_qc[1].caveats == []
    assert out_no_qc[0].confidence == "high"
    sigs = out_no_qc[0].metadata.get("degradation_signals") or []
    assert any(
        s["scope"] == "__root__" and s["field"] == "warnings" and s["severity"] == "info" and s["message"] == root_warn
        for s in sigs
    )
    assert not out_no_qc[1].metadata.get("degradation_signals")


# ── 4. Registry choke point, deduplication (D1-009/D1-025) & unchanged status ─

def test_registry_collect_blocks_routes_deduplicates_and_preserves_status_contract():
    from aria.agents.narrative.warning_router import route_degradation_signals

    reg = registry_with((BulkRnaNarrator(), IntegrationNarrator()))
    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "sample_qc": {
                    "status": "error",
                    "n_samples": 6,
                    "outliers": None,
                    "error": "PCA SVD crash",
                },
                "contrasts": [
                    {
                        "name": "treat_vs_ctrl",
                        "status": "success",
                        "n_significant": 12,
                        "n_upregulated": 7,
                        "n_downregulated": 5,
                        "pathways": {
                            "GO_BP": [{"term": "Immune response", "padj": 0.01}],
                        },
                        "pathway_ora": {
                            "method": "enrichr",
                            "failed_databases": {"KEGG": "ConnectionError"},
                        },
                        "plots": {"pca_error": "figure render skipped"},
                    }
                ],
            },
        }
    }
    blocks = reg.collect_blocks(agent_results)
    by_id = {b.id: b for b in blocks}

    # D1-009 sample_qc error is NOT duplicated.
    qc_blk = by_id["bulk.qc"]
    assert qc_blk.status == "error"
    assert len(qc_blk.caveats) == 1

    # D1-025 failed_databases is NOT duplicated.
    pw_blk = by_id["bulk.pathway.treat_vs_ctrl"]
    assert pw_blk.status == "success"
    kegg_cavs = [c for c in pw_blk.caveats if "KEGG" in c.text]
    assert len(kegg_cavs) == 1

    # plots.pca_error routed to bulk.contrast.treat_vs_ctrl with standard status="success".
    c_blk = by_id["bulk.contrast.treat_vs_ctrl"]
    assert c_blk.status == "success"
    assert any(c.severity == "warning" and "figure render skipped" in c.text for c in c_blk.caveats)

    # Idempotency check: second pass changes nothing.
    snap_before = [(b.id, b.status, b.confidence, [(c.severity, c.text) for c in b.caveats], b.metadata.get("degradation_signals")) for b in blocks]
    blocks_again = route_degradation_signals("bulk_rna_agent", agent_results["bulk_rna_agent"], blocks)
    snap_after = [(b.id, b.status, b.confidence, [(c.severity, c.text) for c in b.caveats], b.metadata.get("degradation_signals")) for b in blocks_again]
    assert snap_before == snap_after

    # Unchanged status contract: compile_public_claims + render_blocks work cleanly.
    ledger = build_run_ledger({}, agent_results)
    compiled = compile_public_claims(blocks, {}, run_ledger=ledger)
    assert compiled.withheld == []
    assert any(
        (c.get("claim_id") or c.get("block_id") if isinstance(c, dict) else getattr(c, "claim_id", None))
        == "bulk.contrast.treat_vs_ctrl"
        for c in compiled.claims
    ), f"Expected bulk.contrast.treat_vs_ctrl in claims; got={compiled.claims}"
    html_out = render_blocks(blocks, strict=False)
    assert "figure render skipped" in html_out
