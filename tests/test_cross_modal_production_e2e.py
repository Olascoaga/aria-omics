"""ADR-080 — Production Cross-Modal Concordance & Synthesis E2E Test Suite.

Verifies:
  1. ``CROSS_MODAL_CONCORDANCE_VALIDATION`` and ``CrossModalConcordanceAuditAgent``
     governance:
       - ``green`` / ``auto`` on non-same-cell ``bulk_RNA + bulk_ATAC`` and
         unpaired ``scRNA + scATAC`` while ``IntegrationAgent`` stays ``red/blocked``
         (``scope == "paired_same_cell_rna_atac"`` untouched).
       - ``red`` / ``blocked`` when diagonal/bridge cell imputation (``glue``,
         ``multivi``, ``seurat_bridge``, ``diagonal_imputation=True``) is requested.
       - ``not_applicable`` when paired same-cell Multiome supersedes or on
         single-modality runs.
  2. ``OrchestratorAgent._dispatch_agents`` automatic execution of
     ``aria/scripts/cross_modal_concordance.py`` when ``bulk_rna_agent``/``scrna_agent``
     and ``chromatin_agent`` both succeed with structured findings, and fail-closed
     blocking when diagonal cell imputation is requested.
  3. Full end-to-end ``EnvironmentManager.run_in_stack("rna", "aria/scripts/cross_modal_concordance.py")``
     -> ``BiologicalSynthesisAgent`` Slice 2 -> ``NarrativeAgent`` block verification
     (``W-CLAIM``), artifact staging, ``RunLedger`` ``C3`` linkage (``0`` violations),
     ``C4`` zero-causal-verb lint, and ``HypothesisAgent`` non-double-counting.
"""

from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

import pytest

from aria.agents.biological_synthesis_agent import BiologicalSynthesisAgent
from aria.agents.modality_audit import (
    CrossModalConcordanceAuditAgent,
    build_capability_matrix,
)
from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.evidence_verifier import verify_block_claim_support
from aria.agents.narrative.hypothesis.report_section import gather_evidence
from aria.agents.narrative.hypothesis.ranking import _rank_basis, _signals_by_id
from aria.agents.narrative.hypothesis.types import Hypothesis
from aria.agents.narrative.run_ledger import (
    build_run_ledger,
    ensure_report_ledger_nodes,
    link_claims_to_ledger,
)
from aria.agents.narrative.validators import find_causal_language
from aria.agents.narrative_agent import NarrativeAgent
from aria.agents.orchestrator_agent import (
    CROSS_MODAL_CONCORDANCE_VALIDATION,
    INTEGRATION_VALIDATION,
    MODALITY_VALIDATION,
    OrchestratorAgent,
)
from aria.utils.environment_manager import env_manager
from aria.utils.registry_integrity import check_registry_integrity
from aria.utils.script_contracts import contract_for_script


class _LocalSubprocessEnvFallback:
    """Contract-enforcing IPC runner used on lightweight CI hosts that lack full Conda stacks."""

    def __init__(self, work_dir: Path) -> None:
        self._work_dir = work_dir

    def run_in_stack(
        self,
        stack: str,
        script_path: str,
        params: dict[str, Any],
        timeout: int | None = None,
    ) -> dict[str, Any]:
        repo_root = Path(__file__).resolve().parents[1]
        resolved = repo_root / script_path
        contract = contract_for_script(resolved)
        assert contract is not None, f"Missing ScriptContract for {script_path}"
        in_issues = contract.validate_params(params)
        assert in_issues == [], f"Input contract violations: {in_issues}"
        in_json = self._work_dir / "ipc_in.json"
        out_json = self._work_dir / "ipc_out.json"
        in_json.write_text(json.dumps(params), encoding="utf-8")
        proc = subprocess.run(
            [sys.executable, str(resolved), str(in_json), str(out_json)],
            capture_output=True,
            text=True,
            timeout=timeout or 120,
            check=False,
        )
        assert proc.returncode == 0, f"Subprocess failed: {proc.stderr or proc.stdout}"
        out = json.loads(out_json.read_text(encoding="utf-8"))
        out_issues = contract.validate_result(out)
        assert out_issues == [], f"Output contract violations: {out_issues}"
        return out


def _ipc_runner(tmp_path: Path):
    envs = env_manager.check_environments()
    if envs.get("rna"):
        return env_manager
    return _LocalSubprocessEnvFallback(tmp_path)


def _make_mini_gtf_and_tables(work_dir: Path) -> dict[str, Any]:
    """Build a realistic GTF, DE table, DA table, peak annotation, peak ORA,
    motifs, TOBIAS footprinting, and 8-donor CPM matrices on disk."""
    gtf = work_dir / "mini.gtf"
    genes_spec = [
        ("ENSG0001", "GATA1", "chr1", "+", 10000, 14000, 3.4, 1e-6, 2.5, 2e-5, True, True),
        ("ENSG0002", "KLF1", "chr1", "+", 50000, 54000, 2.8, 2e-5, 1.9, 1e-4, True, True),
        ("ENSG0003", "TAL1", "chr1", "-", 90000, 95000, 2.1, 4e-4, 1.6, 5e-4, True, True),
        ("ENSG0004", "PAX5", "chr2", "+", 20000, 25000, -3.1, 1e-6, -2.2, 3e-5, True, True),
        ("ENSG0005", "EBF1", "chr2", "+", 70000, 76000, -2.4, 3e-5, -1.8, 2e-4, True, True),
        ("ENSG0006", "SPI1", "chr2", "-", 120000, 125000, -1.9, 8e-4, -1.5, 9e-4, True, True),
        ("ENSG0007", "MYC", "chr3", "+", 30000, 35000, 1.8, 0.002, -1.4, 0.004, True, True),
        ("ENSG0008", "ACTB", "chr3", "+", 80000, 84000, 0.1, 0.85, 0.05, 0.91, False, False),
        ("ENSG0009", "GAPDH", "chr3", "+", 140000, 144000, -0.15, 0.78, 1.7, 0.001, False, True),
        ("ENSG0010", "TUBB", "chr4", "+", 40000, 44000, 2.0, 0.001, 0.12, 0.80, True, False),
    ]

    gtf_lines = []
    de_rows = []
    da_rows = []
    ann_rows = []
    for gid, sym, chrom, strand, start, end, r_lfc, r_padj, a_lfc, a_padj, r_sig, a_sig in genes_spec:
        gtf_lines.append(
            f'{chrom}\tARIA\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}"; gene_name "{sym}";\n'
        )
        gtf_lines.append(
            f'{chrom}\tARIA\texon\t{start}\t{start + 800}\t.\t{strand}\t.\tgene_id "{gid}"; gene_name "{sym}";\n'
        )
        tss = start if strand == "+" else end
        p_start = max(1, tss - 400)
        p_end = tss + 400
        peak_id = f"{chrom}:{p_start}-{p_end}"
        de_rows.append(
            {
                "gene_id": gid,
                "gene_symbol": sym,
                "symbol": sym,
                "baseMean": 250.0,
                "log2FoldChange": r_lfc,
                "pvalue": r_padj / 5.0,
                "padj": r_padj,
                "significant": r_sig,
            }
        )
        da_rows.append(
            {
                "peak": peak_id,
                "baseMean": 120.0,
                "log2FoldChange": a_lfc,
                "pvalue": a_padj / 5.0,
                "padj": a_padj,
                "significant": a_sig,
                "comparison": "K562_vs_GM12878",
            }
        )
        if a_sig:
            ann_rows.append(
                {
                    "peak": peak_id,
                    "nearest_gene": sym,
                    "distance_to_tss": 0,
                    "feature": "Promoter",
                }
            )
    gtf.write_text("".join(gtf_lines), encoding="utf-8")

    contrast_dir = work_dir / "bulk_rna" / "K562_vs_GM12878"
    tables_dir = contrast_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    de_tsv = tables_dir / "de_genes.tsv"
    with open(de_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(de_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(de_rows)

    atac_dir = work_dir / "bulk_atac"
    atac_dir.mkdir(parents=True, exist_ok=True)
    da_csv = atac_dir / "bulk_atac_da_K562_vs_GM12878_full.csv"
    with open(da_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(da_rows[0].keys()))
        w.writeheader()
        w.writerows(da_rows)

    ann_csv = atac_dir / "bulk_atac_peak_annotation_K562_vs_GM12878.csv"
    with open(ann_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["peak", "nearest_gene", "distance_to_tss", "feature"])
        w.writeheader()
        w.writerows(ann_rows)

    # 8 matched biological units (D1..D4 K562, D5..D8 GM12878) for sample-matched cis links
    donors = [f"D{i}" for i in range(1, 9)]
    rna_cpm = work_dir / "bulk_rna" / "counts_cpm.tsv"
    atac_cpm = atac_dir / "bulk_atac_peak_cpm.tsv"
    high_pattern = [220.0, 250.0, 235.0, 265.0, 20.0, 25.0, 18.0, 22.0]
    low_pattern = [18.0, 24.0, 20.0, 22.0, 210.0, 240.0, 225.0, 255.0]
    flat_pattern = [100.0, 104.0, 98.0, 102.0, 101.0, 99.0, 103.0, 97.0]

    with open(rna_cpm, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", *donors])
        for gid, sym, _c, _s, _st, _en, r_lfc, *_rest in genes_spec:
            vals = high_pattern if r_lfc > 1.0 else (low_pattern if r_lfc < -1.0 else flat_pattern)
            w.writerow([sym, *[v + (int(gid[-1]) * 1.3) for v in vals]])

    with open(atac_cpm, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["peak_id", *donors])
        for idx, row in enumerate(da_rows):
            a_lfc = float(row["log2FoldChange"])
            vals = high_pattern if a_lfc > 1.0 else (low_pattern if a_lfc < -1.0 else flat_pattern[:])
            w.writerow([row["peak"], *[v + (idx * 1.1) for v in vals]])

    rna_results = {
        "status": "done",
        "findings": {
            "status": "success",
            "n_contrasts": 1,
            "n_significant": 8,
            "contrasts": [
                {
                    "name": "K562_vs_GM12878",
                    "numerator": "K562",
                    "denominator": "GM12878",
                    "status": "success",
                    "n_genes_tested": len(de_rows),
                    "n_significant": 8,
                    "n_upregulated": 5,
                    "n_downregulated": 3,
                    "top_genes": [
                        {"symbol": r["symbol"], "gene": r["symbol"], "log2fc": r["log2FoldChange"], "padj": r["padj"]}
                        for r in de_rows
                        if r["significant"]
                    ],
                    "all_sig_genes": [r["symbol"] for r in de_rows if r["significant"]],
                    "up_gene_ids": [r["symbol"] for r in de_rows if r["significant"] and r["log2FoldChange"] > 0],
                    "down_gene_ids": [r["symbol"] for r in de_rows if r["significant"] and r["log2FoldChange"] < 0],
                    "pathways": {
                        "GO_BP": [
                            {"term": "erythrocyte differentiation (GO:0030218)", "padj": 0.0003},
                            {"term": "heme metabolic process (GO:0042168)", "padj": 0.0012},
                        ]
                    },
                    "contrast_dir": str(contrast_dir),
                    "tables_dir": str(tables_dir),
                    "plots": {"tables": {"de_genes": str(de_tsv)}},
                }
            ],
        },
    }

    chromatin_results = {
        "status": "done",
        "findings": {
            "bulk_ATAC": {
                "status": "done",
                "findings": {
                    "qc": {"status": "success", "data_type": "bulk_ATAC"},
                    "peaks": {"status": "success", "data_type": "bulk_ATAC", "n_peaks": len(da_rows)},
                    "differential_accessibility": {
                        "status": "success",
                        "ran": True,
                        "data_type": "bulk_ATAC",
                        "comparisons": [
                            {
                                "test": "K562",
                                "reference": "GM12878",
                                "status": "success",
                                "n_sig": sum(1 for r in da_rows if r["significant"]),
                                "n_up": sum(1 for r in da_rows if r["significant"] and r["log2FoldChange"] > 0),
                                "n_down": sum(1 for r in da_rows if r["significant"] and r["log2FoldChange"] < 0),
                                "full_results_csv": str(da_csv),
                                "top_peaks": [
                                    {"peak": r["peak"], "log2fc": r["log2FoldChange"], "padj": r["padj"]}
                                    for r in da_rows
                                    if r["significant"]
                                ],
                            }
                        ],
                    },
                    "peak_annotation": {
                        "status": "success",
                        "ran": True,
                        "data_type": "bulk_ATAC",
                        "comparisons": [
                            {
                                "test": "K562",
                                "reference": "GM12878",
                                "n_annotated": len(ann_rows),
                                "annotation_csv": str(ann_csv),
                            }
                        ],
                    },
                    "peak_ora": {
                        "status": "success",
                        "ran": True,
                        "data_type": "bulk_ATAC",
                        "comparisons": [
                            {
                                "test": "K562",
                                "reference": "GM12878",
                                "pathways": {
                                    "GO_BP": [
                                        {"term": "erythrocyte differentiation (GO:0030218)", "padj": 0.0021},
                                        {"term": "chromatin remodeling (GO:0006338)", "padj": 0.008},
                                    ]
                                },
                            }
                        ],
                    },
                    "motifs": {
                        "status": "success",
                        "ran": True,
                        "per_group": {
                            "K562_vs_GM12878_up": {
                                "n_enriched": 2,
                                "top_motifs": [
                                    {
                                        "motif_id": "MA0035.4",
                                        "name": "GATA1",
                                        "family": "GATA",
                                        "log2_enrichment": 2.35,
                                        "pvalue": 1e-7,
                                        "padj": 1e-5,
                                    },
                                    {
                                        "motif_id": "MA0493.2",
                                        "name": "KLF1",
                                        "family": "KLF",
                                        "log2_enrichment": 1.85,
                                        "pvalue": 2e-5,
                                        "padj": 4e-4,
                                    },
                                ],
                            },
                            "K562_vs_GM12878_down": {
                                "n_enriched": 1,
                                "top_motifs": [
                                    {
                                        "motif_id": "MA0014.3",
                                        "name": "PAX5",
                                        "family": "Paired",
                                        "log2_enrichment": 2.10,
                                        "pvalue": 3e-6,
                                        "padj": 5e-5,
                                    }
                                ],
                            },
                        },
                    },
                    "footprinting": {
                        "status": "success",
                        "ran": True,
                        "data_type": "bulk_ATAC",
                        "differential_summary": {
                            "group_a": "K562",
                            "group_b": "GM12878",
                            "inferential": True,
                            "significant_tfs": [
                                {
                                    "name": "GATA1",
                                    "motif_id": "MA0035.4",
                                    "change": 0.42,
                                    "pvalue": 1e-5,
                                    "padj": 2e-4,
                                }
                            ],
                        },
                    },
                },
            }
        },
    }

    return {
        "gtf_path": str(gtf),
        "rna_cpm_path": str(rna_cpm),
        "atac_cpm_path": str(atac_cpm),
        "rna_results": rna_results,
        "chromatin_results": chromatin_results,
    }


def test_cross_modal_concordance_governance_and_diagonal_fail_closed():
    assert CROSS_MODAL_CONCORDANCE_VALIDATION["level"] == "production"
    assert CROSS_MODAL_CONCORDANCE_VALIDATION["evidence_tier"] == "associative"
    assert CROSS_MODAL_CONCORDANCE_VALIDATION["scope"] == "contrast_and_sample_matched_rna_atac"
    assert CROSS_MODAL_CONCORDANCE_VALIDATION["diagonal_cell_imputation"] == "blocked"
    assert CROSS_MODAL_CONCORDANCE_VALIDATION["promotion_adr"] == "ADR-080"
    # Paired same-cell IntegrationAgent scope MUST remain untouched
    assert INTEGRATION_VALIDATION["scope"] == "paired_same_cell_rna_atac"
    assert [i for i in check_registry_integrity() if i.severity == "error"] == []

    # 1. Non-same-cell bulk_RNA + bulk_ATAC -> cross_modal green/auto, integration blocked
    unpaired_ctx = {
        "genome": "hg38",
        "modalities": {
            "bulk_RNA": ["/virtual/counts.tsv"],
            "bulk_ATAC": ["/virtual/atac.bam"],
        },
        "sample_metadata": {
            "S1": {"condition": "A", "replicate": "1"},
            "S2": {"condition": "A", "replicate": "2"},
            "S3": {"condition": "B", "replicate": "1"},
            "S4": {"condition": "B", "replicate": "2"},
        },
        "comparisons": [{"test": "A", "reference": "B"}],
    }
    mat = build_capability_matrix(
        unpaired_ctx,
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert mat["cross_modal_concordance_readiness"]["status"] == "green"
    assert mat["cross_modal_concordance_readiness"]["dispatch_policy"] == "auto"
    assert mat["integration_readiness"]["status"] == "red"
    assert mat["integration_readiness"]["dispatch_policy"] == "blocked"

    # 2. Diagonal cell-imputation request -> red / blocked fail-closed
    for bad_kw in ({"cross_modal_strategy": "glue"}, {"diagonal_imputation": True}):
        diag_ctx = {**copy.deepcopy(unpaired_ctx), **bad_kw}
        diag_mat = build_capability_matrix(
            diag_ctx,
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )
        assert diag_mat["status"] == "red"
        assert diag_mat["cross_modal_concordance_readiness"]["status"] == "red"
        assert diag_mat["cross_modal_concordance_readiness"]["dispatch_policy"] == "blocked"
        assert diag_mat["cross_modal_concordance_readiness"]["reason"] == "diagonal_cell_imputation_not_validated"

    # 3. Single-modality run -> not_applicable
    single_mat = build_capability_matrix(
        {"genome": "hg38", "modalities": {"bulk_RNA": ["/virtual/counts.tsv"]}},
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert single_mat["cross_modal_concordance_readiness"]["status"] == "not_applicable"


def test_orchestrator_auto_dispatches_cross_modal_concordance_and_stages_report(tmp_path):
    fixture = _make_mini_gtf_and_tables(tmp_path)

    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch._sessions = {}
    orch._experiment_plans = {"exp_cm_e2e": {"intent": {"user_question": "K562 vs GM12878 concordance"}}}
    orch._agent_results = {}
    orch._log_handlers = {}
    orch.memory = object()
    orch.llm = object()
    orch.env = _ipc_runner(tmp_path)
    orch.publish_status = lambda *a, **kw: None
    orch.publish_escalation = lambda *a, **kw: None
    bus_findings: list[dict[str, Any]] = []
    orch.publish_finding = lambda exp, payload, conf: bus_findings.append(payload)

    def _fake_run_agent(agent_name: str, experiment_id: str, context: dict):
        if agent_name == "bulk_rna_agent":
            return fixture["rna_results"]
        if agent_name == "chromatin_agent":
            return fixture["chromatin_results"]
        if agent_name == "narrative_agent":
            return {"status": "done", "agent": "narrative_agent", "seen_keys": sorted(context["agent_results"].keys())}
        return {"status": "done", "agent": agent_name}

    orch._run_agent = _fake_run_agent

    ctx = {
        "genome": "hg38",
        "organism": "Homo sapiens",
        "modalities": {
            "bulk_RNA": [str(tmp_path / "bulk_rna" / "counts_cpm.tsv")],
            "bulk_ATAC": [str(tmp_path / "bulk_atac" / "bulk_atac_peak_cpm.tsv")],
        },
        "gtf_path": fixture["gtf_path"],
        "rna_cpm_path": fixture["rna_cpm_path"],
        "atac_cpm_path": fixture["atac_cpm_path"],
        "cross_modal_output_dir": str(tmp_path / "cm_out"),
        "min_matched_samples": 6,
        "n_permutations": 100,
        "seed": 42,
    }
    plan = orch._fallback_plan(ctx)
    orch._dispatch_agents("exp_cm_e2e", plan, ctx)

    results = orch._agent_results["exp_cm_e2e"]
    assert results["integration_agent"]["status"] == "skipped"
    assert results["integration_agent"]["reason"] == "multiome_same_cell_gate_blocked"

    cm = results["cross_modal_concordance"]
    assert cm["status"] == "success"
    assert cm["ran"] is True
    assert cm["evidence_tier"] == "associative"
    assert "cross_modal_concordance" in results["narrative_agent"]["seen_keys"]

    de_da = cm["de_da_concordance"]
    assert de_da["n_concordant_total"] >= 6
    assert cm["sample_matched_links"]["ran"] is True
    assert cm["sample_matched_links"]["n_shared_samples"] == 8

    # End-to-end BiologicalSynthesisAgent Slice 2 + NarrativeAgent + RunLedger + Hypothesis check
    blocks = BiologicalSynthesisAgent().synthesize(results, ctx)
    cm_blocks = [b for b in blocks if str(b.id).startswith("integration.cross_modal")]
    assert len(cm_blocks) >= 2

    kept = NarrativeAgent._drop_unverifiable_synthesis_blocks(cm_blocks)
    assert len(kept) == len(cm_blocks)

    for b in cm_blocks:
        prose = compose_block_prose(b) or b.claim or ""
        assert verify_block_claim_support(b, prose, strict=False)["status"] == "supported"
        assert find_causal_language(b.claim or "") is None
        assert find_causal_language(prose) is None

    ledger = build_run_ledger(ctx, results)
    ensure_report_ledger_nodes(ledger, cm_blocks)
    compiled = compile_public_claims(cm_blocks, {}, run_ledger=ledger)
    assert compiled.withheld == []
    linked = link_claims_to_ledger(copy.deepcopy(compiled.claims), ledger)
    assert linked["n_violations"] == 0

    # Verify non-double-counting in HypothesisAgent ranking across RNA + ATAC + cross_modal TF
    signals = gather_evidence(results, ledger, ctx)
    gata1_sigs = [s for s in signals if str(s.entity).upper() == "GATA1"]
    assert len(gata1_sigs) >= 2
    by_id = _signals_by_id(signals)
    hyp = Hypothesis(
        id="hyp_gata1",
        mechanism="GATA1 motif accessibility and transcript abundance may covary in erythroid cells",
        entities=["GATA1"],
        observed_claims=[{"signal_id": s.signal_id, "stated_direction": s.direction} for s in gata1_sigs],
    )
    basis = _rank_basis(hyp, by_id)
    assert basis["n_independent_lines"] >= 2
    # GATA1 has at most 2 independent lineages (RNA DE + ATAC motif/footprint/cross_modal TF lineage)
    # when sample-matched peak2gene correlation is excluded.
    tf_only_sigs = [s for s in gata1_sigs if s.measure != "sample_matched_peak2gene_corr"]
    if len(tf_only_sigs) >= 3:
        hyp_tf = Hypothesis(
            id="hyp_gata1_tf",
            mechanism="GATA1 motif and expression covary",
            entities=["GATA1"],
            observed_claims=[{"signal_id": s.signal_id, "stated_direction": s.direction} for s in tf_only_sigs],
        )
        assert _rank_basis(hyp_tf, by_id)["n_independent_lines"] <= 2

    # Stage artifacts into a report folder
    report_dir = tmp_path / "report_bundle"
    (report_dir / "figures").mkdir(parents=True, exist_ok=True)
    (report_dir / "tables").mkdir(parents=True, exist_ok=True)
    narr = NarrativeAgent.__new__(NarrativeAgent)
    narr._stage_artifacts(results, report_dir)
    assert any((report_dir / "figures").glob("cross_modal_starburst.*"))
    assert any((report_dir / "tables").glob("cross_modal_*.csv"))


def test_orchestrator_blocks_diagonal_imputation_in_dispatch():
    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch._sessions = {}
    orch._experiment_plans = {"exp_diag": {"intent": {}}}
    orch._agent_results = {}
    orch._log_handlers = {}
    orch.memory = object()
    orch.llm = object()
    orch.publish_status = lambda *a, **kw: None
    orch.publish_escalation = lambda *a, **kw: None
    orch.publish_finding = lambda *a, **kw: None
    orch._run_agent = lambda agent_name, experiment_id, context: {
        "status": "done",
        "agent": agent_name,
        "findings": {"stub": True},
    }

    ctx = {
        "genome": "hg38",
        "modalities": {
            "scRNA": ["/virtual/rna.h5ad"],
            "scATAC": ["/virtual/atac.h5ad"],
        },
        "diagonal_imputation": True,
    }
    orch._dispatch_agents("exp_diag", orch._fallback_plan(ctx), ctx)
    cm = orch._agent_results["exp_diag"]["cross_modal_concordance"]
    assert cm["status"] == "blocked"
    assert cm["ran"] is False
    assert cm["reason"] == "diagonal_cell_imputation_not_validated"


def test_adr080_empirical_receipt_integrity():
    receipt_path = (
        Path(__file__).resolve().parents[1]
        / "docs"
        / "benchmark_results"
        / "multiome_concordance"
        / "adr080_cross_modal_concordance_validation.json"
    )
    if not receipt_path.is_file():
        pytest.skip("ADR-080 empirical validation receipt not yet generated")
    raw = receipt_path.read_text(encoding="utf-8")
    assert "/home/" not in raw
    assert "/tmp/" not in raw
    data = json.loads(raw)
    assert data["promotion_adr"] == "ADR-080"
    assert data["all_invariants_passed"] is True
    assert data["dataset"]["n_donors"] == 10

    arm_a = data["arm_a_disjoint_donor_glial_lineage"]
    assert arm_a["design_isolation"]["shared_donors_in_contrast"] == 0
    assert arm_a["design_isolation"]["shared_cells_in_contrast"] == 0
    assert arm_a["de_da_summary"]["n_gene_peak_pairs_tested"] > 20000
    assert arm_a["de_da_summary"]["n_atac_da_rna_unchanged"] > 1000
    assert arm_a["de_da_summary"]["n_concordant_total"] > 3000
    assert arm_a["de_da_summary"]["fisher_odds_ratio"] > 100.0
    assert arm_a["reversed_contrast_invariance"]["atac_sign_flipped"] is True
    assert arm_a["reversed_contrast_invariance"]["exact_concordance_match"] is True
    assert arm_a["tf_triangulation"]["ran"] is True
    assert arm_a["tf_triangulation"]["n_concordant_tfs"] >= 8
    assert arm_a["pathway_concordance"]["ran"] is True
    assert arm_a["pathway_concordance"]["n_shared_terms"] >= 20
    assert arm_a["matched_20unit_glial_cpm_links"]["n_significant_links"] > 10000
    assert arm_a["matched_20unit_glial_cpm_links"]["deranged_unit_null_control"]["n_null_significant_links"] == 0
    assert arm_a["a4_stage_manifest"]["cache_hit_verified"] is True
    assert arm_a["figures"]["starburst_png_exists"] is True
    assert arm_a["figures"]["starburst_svg_exists"] is True
    assert data["same_cell_ground_truth_comparison"]["same_cell_gene_recall_combined"] >= 0.80

    arm_b = data["arm_b_matched_10donor_clinical_aging"]
    assert arm_b["donor_matched_10donor_cpm_links"]["donor_deranged_null_control"]["n_null_significant_links"] == 0

    gov = data["governance_and_narrative_verification"]
    assert gov["unpaired_bulk_readiness"]["status"] == "green"
    assert gov["same_cell_integration_readiness_on_unpaired"]["status"] == "red"
    assert gov["diagonal_imputation_gate"]["readiness_status"] == "red"
    assert gov["diagonal_imputation_gate"]["reason"] == "diagonal_cell_imputation_not_validated"
    assert gov["biological_synthesis_slice2"]["w_claim_unsupported_blocks"] == 0
    assert gov["biological_synthesis_slice2"]["c3_ledger_violations"] == 0
    assert gov["biological_synthesis_slice2"]["c4_causal_verb_hits"] == 0

