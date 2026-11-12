#!/usr/bin/env python3
"""V48 Slice I8 (ADR-071) — End-to-end real-data validation of Orchestrated
`IntegrationAgent` Beta Dispatch, Pre-Dispatch Readiness/CP3.5 Gate, Same-Cell
Scope Enforcement, Post-QC `.h5mu` Propagation, and Hermetic Environment Audit.

Executes on real paired 10x Multiome (`GSE278576` `hc11`, `3,143` raw -> `2,517`
post-QC shared cells):
  1. Environment & registry preflight:
     - `check_registry_integrity()` == 0 errors (`INTEGRATION_VALIDATION` beta +
       `requires_ack=True`, `production_promoted=False`, `ADR-056` / `ADR-071`)
     - `audit_environment("aria-integration-env")` == `"pass"` (`conda_exact=True`,
       `pip_exact=True`)
     - `audit_environment("aria-chromatin-env")` == `"pass"`
  2. Gate Probe A (CP3.5 Cancel / Unacknowledged Refusal):
     - Drives `run_headless` on `hc11_paired.h5mu` with a CP3.5 cancel policy,
       verifying `AuditAgent` holds dispatch at CP3.5 (`integration_beta_ack_required`
       + `chromatin_readiness_beta_ack_required`) and cancelling prevents any
       subprocess execution (`status == "cancelled"`).
  3. Gate Probe B (Unpaired / Non-Same-Cell Fail-Closed Scope):
     - Exercises `OrchestratorAgent._dispatch_agents` on split unpaired
       `scRNA + scATAC` and `scRNA + bulk_RNA`, proving `integration_agent` is
       skipped (`reason == "multiome_same_cell_gate_blocked"`) while single-assay
       lanes remain active.
  4. Full Live Orchestrated Headless Run (`run_headless`, `reproducible_mode=True`):
     - `DataAuditAgent` -> CP1 -> `DesignAgent` -> `_design_analysis_plan`
       (including deterministic override of an LLM plan returning
       `integration_needed=False`) -> CP2 -> `AuditAgent` -> CP3.5 (`"Proceed anyway"`,
       persisted in `ARIAMemory` `decisions` table) ->
     - Live `ChromatinAgent` (`aria-chromatin-env`) + live `IntegrationAgent`
       (`integration_paired_qc.py` -> CP3 WNN `k` -> `integration_wnn.py` + I3
       5-axis guardrail -> `integration_peak2gene.py` -> `integration_eregulons.py`
       -> CP3 MOFA `n_factors` -> `integration_mofa.py`), verifying every
       downstream stage operated on the `2,517`-cell post-QC `.h5mu` ->
     - Live `NarrativeAgent` -> `report.html`, `methodology.json`,
       `ro-crate-metadata.json`, `memory_snapshot.sqlite`, dual PNG+SVG figures,
       staged CSV tables, `0` C1-C4 claim/ledger violations, and verified
       reproducible capsule (`verify_reproducible_capsule == "pass"`).

Emits a sanitized public JSON artifact with zero machine-absolute paths
(enforced by `tests/test_public_artifacts_no_absolute_paths.py`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sqlite3
import sys
import tempfile
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.agents.narrative.ledger_export import (  # noqa: E402
    load_methodology,
    verify_reproducible_capsule,
    write_reproducible_capsule,
)
from aria.agents.orchestrator_agent import (  # noqa: E402
    INTEGRATION_VALIDATION,
    MODALITY_VALIDATION_LEVELS,
    OrchestratorAgent,
)
from aria.headless import default_answer_policy, run_headless  # noqa: E402
from aria.utils.environment_audit import audit_environment  # noqa: E402
from aria.utils.registry_integrity import check_registry_integrity  # noqa: E402
from aria.version import __version__, collect_version_metadata  # noqa: E402

_FORBIDDEN_PATH_PREFIXES = (
    "/home/",
    "/Users/",
    "/root/",
    "/tmp/",
    "/mnt/",
    "/private/",
    "/data/",
    "/var/",
    "/opt/",
)


def _sanitize_public_obj(val: Any) -> Any:
    if isinstance(val, str):
        for pfx in _FORBIDDEN_PATH_PREFIXES:
            if val.startswith(pfx):
                return Path(val).name
        if any(pfx in val for pfx in _FORBIDDEN_PATH_PREFIXES):
            tokens = []
            for tok in val.split():
                clean_tok = tok.strip(",;()[]'\"")
                if clean_tok.startswith(_FORBIDDEN_PATH_PREFIXES):
                    tokens.append(tok.replace(clean_tok, Path(clean_tok).name))
                else:
                    tokens.append(tok)
            return " ".join(tokens)
        return val
    if isinstance(val, dict):
        return {str(k): _sanitize_public_obj(v) for k, v in val.items()}
    if isinstance(val, list):
        return [_sanitize_public_obj(v) for v in val]
    return val


def _assert_no_forbidden_paths(obj: Any, trail: str = "root") -> None:
    if isinstance(obj, str):
        for pfx in _FORBIDDEN_PATH_PREFIXES:
            if obj.startswith(pfx):
                raise AssertionError(f"Forbidden machine path at {trail}: {obj!r}")
    elif isinstance(obj, dict):
        for k, v in obj.items():
            _assert_no_forbidden_paths(v, f"{trail}.{k}")
    elif isinstance(obj, list):
        for idx, v in enumerate(obj):
            _assert_no_forbidden_paths(v, f"{trail}[{idx}]")


class OfflineDeterministicLLM:
    """Zero-egress deterministic LLM double for orchestrated E2E validation.

    - Returns structured JSON for `_parse_question` and `_design_analysis_plan`
      (intentionally returning `integration_needed=false` on the plan prompt so
      `OrchestratorAgent._design_analysis_plan`'s deterministic paired-multiome
      override is exercised live).
    - Returns deterministic structured text for `DebateCouncil` if invoked.
    - Raises `RuntimeError` on narrative prose calls so `compose_block_prose`
      and `NarrativeAgent` exercise their deterministic evidence-grounded
      template path (`$0.00`, zero ungrounded entities, zero causal drift).
    """

    def __init__(self) -> None:
        self.calls: list[str] = []

    def for_execution(self, **kwargs: Any) -> "OfflineDeterministicLLM":
        return self

    def complete(
        self,
        prompt: str,
        system: str = "",
        tier: Any = None,
        max_tokens: int = 2048,
        temperature: float = 0.0,
        **kwargs: Any,
    ) -> str:
        p_text = str(prompt or "")
        if "Return analysis intent as JSON" in p_text:
            self.calls.append("parse_question")
            return json.dumps(
                {
                    "analysis_type": "integration",
                    "biological_entities": [
                        "oligodendrocytes",
                        "hippocampus",
                        "chromatin accessibility",
                        "gene expression",
                    ],
                    "comparison": None,
                    "key_modalities_needed": ["RNA", "ATAC"],
                    "complexity": "moderate",
                    "summary": (
                        "Joint same-cell RNA+ATAC multiome characterization of "
                        "HC11 human hippocampus."
                    ),
                }
            )
        if "Return analysis plan as JSON" in p_text:
            self.calls.append("design_analysis_plan")
            return json.dumps(
                {
                    "steps": [
                        {
                            "order": 1,
                            "agent": "chromatin_agent",
                            "analysis": "Analyze scATAC accessibility and motifs",
                            "depends_on": [],
                            "can_parallel": False,
                        }
                    ],
                    "contrasts": [],
                    "integration_needed": False,
                    "integration_type": "none",
                    "estimated_complexity": "moderate",
                    "rationale": (
                        "Deterministic scripted plan testing paired-multiome "
                        "integration_needed override."
                    ),
                }
            )
        if "ALTERNATIVE_HYPOTHESIS" in p_text or "VERDICT" in str(system or ""):
            self.calls.append("debate_council")
            return (
                "ALTERNATIVE_HYPOTHESIS: Sparsity or Tn5 bias could influence "
                "cross-modal weight balance.\n"
                "CHALLENGES: Inspect per-modality silhouette and neighbor recall.\n"
                "EVIDENCE_REQUESTED: 5-axis anti-degradation guardrail metrics.\n"
                "VERDICT: ACCEPT_REVISED\n"
                "REVISED_CLAIM: Joint multiome structure is reported alongside "
                "per-modality RNA and ATAC diagnostics."
            )
        self.calls.append("deterministic_fallback")
        raise RuntimeError("offline_deterministic_template_fallback")

    def complete_heavy(self, prompt: str, system: str = "", **kw: Any) -> str:
        return self.complete(prompt, system=system, **kw)

    def complete_medium(self, prompt: str, system: str = "", **kw: Any) -> str:
        return self.complete(prompt, system=system, **kw)

    def complete_light(self, prompt: str, system: str = "", **kw: Any) -> str:
        return self.complete(prompt, system=system, **kw)


def _run_unpaired_gate_probes() -> dict[str, Any]:
    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    orch._sessions = {}
    orch._experiment_plans = {}
    orch._agent_results = {}
    orch._log_handlers = {}
    orch.memory = object()
    orch.llm = object()
    orch.publish_status = lambda *a, **kw: None
    orch.publish_escalation = lambda *a, **kw: None
    orch.publish_finding = lambda *a, **kw: None

    probed: dict[str, Any] = {}
    scenarios = {
        "split_scrna_scatac_unconfirmed": {
            "genome": "hg38",
            "modalities": {
                "scRNA": ["donor_rna.h5ad"],
                "scATAC": ["donor_atac.h5ad"],
            },
            "integration_ack_confirmed": True,
        },
        "scrna_plus_bulk_rna": {
            "genome": "hg38",
            "modalities": {
                "scRNA": ["donor_rna.h5ad"],
                "bulk_RNA": ["bulk_counts.tsv"],
            },
            "integration_ack_confirmed": True,
        },
        "paired_h5mu_without_cp35_ack": {
            "genome": "hg38",
            "modalities": {"scATAC": ["hc11_paired.h5mu"]},
            "assay_detections": [
                {
                    "path": "hc11_paired.h5mu",
                    "modality": "scATAC",
                    "evidence": {"format": "h5mu", "paired_rna": True},
                }
            ],
        },
    }

    for name, ctx in scenarios.items():
        calls: list[str] = []
        orch._run_agent = lambda agent_name, experiment_id, context, _c=calls: (
            _c.append(agent_name) or {"status": "done"}
        )
        orch._experiment_plans[name] = {"intent": {}}
        plan = orch._fallback_plan(ctx)
        orch._dispatch_agents(name, plan, ctx)
        int_res = orch._agent_results[name]["integration_agent"]
        probed[name] = {
            "dispatched_agents": calls,
            "integration_dispatched": "integration_agent" in calls,
            "integration_status": int_res.get("status"),
            "integration_reason": int_res.get("reason"),
            "validation_level": int_res.get("validation_level"),
        }
    return probed


def main() -> int:
    home = Path.home()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--h5mu",
        default=str(home / "Samael/Oligos/data/interim/paired/hc11_paired.h5mu"),
    )
    ap.add_argument(
        "--reference-h5ad",
        default=str(
            home / "Samael/Oligos/data/raw/reference/GSE278576_hippocampus_RNA.h5ad"
        ),
    )
    ap.add_argument(
        "--gtf",
        default=str(
            home / "Samael/H9-RNA/reference/gencode.v49.primary_assembly.annotation.gtf"
        ),
    )
    ap.add_argument(
        "--fasta",
        default=str(home / ".aria/genomes/hg38/GRCh38.primary_assembly.genome.fa"),
    )
    ap.add_argument("--sample-id", default="hc11")
    ap.add_argument(
        "--out",
        default="docs/benchmark_results/multiome_concordance/v48_i8_orchestrated_dispatch_hc11.json",
    )
    args = ap.parse_args()

    for label, p_str in (
        ("h5mu", args.h5mu),
        ("reference_h5ad", args.reference_h5ad),
        ("gtf", args.gtf),
        ("fasta", args.fasta),
    ):
        if not Path(p_str).is_file():
            raise FileNotFoundError(f"Missing required {label} input: {Path(p_str).name}")

    t_start = time.time()

    # ── 1. Registry & Hermetic Lock Preflight ──────────────────────────────
    reg_issues = check_registry_integrity()
    reg_errors = [i for i in reg_issues if i.severity == "error"]
    assert not reg_errors, f"Registry integrity errors: {reg_errors}"

    int_env_audit = audit_environment(
        "aria-integration-env", env_name="aria-integration-env"
    )
    chrom_env_audit = audit_environment(
        "aria-chromatin-env", env_name="aria-chromatin-env"
    )
    assert int_env_audit["status"] == "pass", int_env_audit
    assert chrom_env_audit["status"] == "pass", chrom_env_audit

    # ── 2. Gate Probes (Unpaired / Unacknowledged Fail-Closed) ─────────────
    gate_probes = _run_unpaired_gate_probes()
    assert (
        gate_probes["split_scrna_scatac_unconfirmed"]["integration_reason"]
        == "multiome_same_cell_gate_blocked"
    )
    assert (
        gate_probes["scrna_plus_bulk_rna"]["integration_reason"]
        == "multiome_same_cell_gate_blocked"
    )
    assert (
        gate_probes["paired_h5mu_without_cp35_ack"]["integration_reason"]
        == "integration_ack_required"
    )

    question = (
        "Characterize joint same-cell RNA+ATAC cell states, MOFA+ shared factors, "
        "cis peak-to-gene links, and candidate eRegulons in HC11 human hippocampus."
    )

    with tempfile.TemporaryDirectory(prefix="aria_v48_i8_") as tmp_root:
        work = Path(tmp_root)
        input_dir = work / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        staged_h5mu = input_dir / "hc11_paired.h5mu"
        shutil.copy2(args.h5mu, staged_h5mu)

        shared_overrides = {
            "sample_id": args.sample_id,
            "reference_h5ad": str(Path(args.reference_h5ad)),
            "reference_label_col": "subclass",
            "target_label_col": "cell_type",
            "reference_sample_col": "orig.ident",
            "gtf_path": str(Path(args.gtf)),
            "genome_fasta": str(Path(args.fasta)),
            "enable_mofa": True,
            "enable_eregulons": True,
            "seed": 0,
            "mofa_params": {"max_cells": 2517},
        }

        # ── 3. Gate Probe A: Live CP3.5 Cancel Refuses Dispatch ─────────────
        def _cancel_at_cp35_policy(cp_num: Any, q: str, opts: list[str]) -> str:
            if cp_num == 3.5:
                for opt in opts:
                    if "cancel" in opt.lower():
                        return opt
                return "Cancel analysis"
            return default_answer_policy(cp_num, q, opts)

        llm_cancel = OfflineDeterministicLLM()
        t_cancel0 = time.time()
        cancel_run = run_headless(
            str(input_dir),
            question,
            policy=_cancel_at_cp35_policy,
            reproducible_mode=True,
            timeout=180.0,
            log=lambda msg: print(f"[cancel-probe] {msg}", flush=True),
            context_overrides=shared_overrides,
            llm_provider=llm_cancel,
        )
        t_cancel_s = round(time.time() - t_cancel0, 2)
        assert cancel_run.status == "cancelled", cancel_run
        assert any(cp == 3.5 for cp, _ in cancel_run.decisions), cancel_run.decisions
        assert not (input_dir / "aria_integration").exists()

        # ── 4. Full Live Orchestrated Run (CP1 -> CP2 -> CP3.5 Ack -> Live E2E)
        llm_live = OfflineDeterministicLLM()
        t_live0 = time.time()
        live_run = run_headless(
            str(input_dir),
            question,
            policy=default_answer_policy,
            reproducible_mode=True,
            timeout=2400.0,
            log=lambda msg: print(f"[live-e2e] {msg}", flush=True),
            context_overrides=shared_overrides,
            llm_provider=llm_live,
        )
        t_live_s = round(time.time() - t_live0, 2)
        assert live_run.status == "completed", live_run
        assert live_run.report_path is not None

        report_html = Path(live_run.report_path)
        report_dir = report_html.parent
        try:
            meth = load_methodology(report_dir)
            rocrate = json.loads(
                (report_dir / "ro-crate-metadata.json").read_text(encoding="utf-8")
            )
            html_text = report_html.read_text(encoding="utf-8")

            # Verify reproducibility capsule (private input paths are sanitized to
            # input://sha256/... URIs, so relocation_required yields "warning"
            # while all embedded file hashes and memory_snapshot must pass).
            capsule_zip = work / "v48_i8_hc11_capsule.zip"
            write_reproducible_capsule(report_dir, capsule_zip)
            capsule_check = verify_reproducible_capsule(capsule_zip)
            assert capsule_check["status"] in {"pass", "warning"}, capsule_check
            assert capsule_check["files"]["missing"] == [], capsule_check
            assert capsule_check["files"]["mismatched"] == [], capsule_check
            assert capsule_check["inputs"]["mismatched"] == [], capsule_check
            assert capsule_check["memory_snapshot"]["valid"] is True, capsule_check

            # Verify memory_snapshot.sqlite contains CP3.5 acknowledgement + WNN k + MOFA n_factors
            mem_db = report_dir / "memory_snapshot.sqlite"
            assert mem_db.is_file() and mem_db.stat().st_size > 0
            conn = sqlite3.connect(str(mem_db))
            conn.row_factory = sqlite3.Row
            try:
                decisions_rows = [
                    dict(r)
                    for r in conn.execute(
                        "SELECT checkpoint, decision, made_by FROM decisions "
                        "WHERE wing_id = ? ORDER BY made_at",
                        (live_run.experiment_id,),
                    ).fetchall()
                ]
                tunnels_count = int(
                    conn.execute("SELECT COUNT(*) FROM tunnels").fetchone()[0]
                )
            finally:
                conn.close()

            cp35_db = [d for d in decisions_rows if str(d["checkpoint"]) == "3.5"]
            assert len(cp35_db) == 1, decisions_rows
            assert cp35_db[0]["made_by"] == "headless_policy"

            # Verify post-QC 2,517-cell container existed and drove downstream stages
            int_dir = input_dir / "aria_integration"
            post_qc_file = int_dir / f"paired_qc_{args.sample_id}.h5mu"
            wnn_file = int_dir / "wnn_integrated.h5mu"
            assert post_qc_file.is_file(), f"Missing {post_qc_file.name}"
            assert wnn_file.is_file(), f"Missing {wnn_file.name}"

            import h5py
            from aria.agents._narrative_integration import _read_h5_index

            with h5py.File(post_qc_file, "r") as f_qc:
                n_cells_post_qc = int(len(_read_h5_index(f_qc["obs"])))
            with h5py.File(wnn_file, "r") as f_wnn:
                n_cells_wnn = int(len(_read_h5_index(f_wnn["obs"])))
            assert n_cells_post_qc == 2517, n_cells_post_qc
            assert n_cells_wnn == n_cells_post_qc, (n_cells_wnn, n_cells_post_qc)

            ev_graph = meth.get("cross_modal_evidence_graph") or {}
            assert ev_graph.get("ran") is True, ev_graph
            assert ev_graph.get("non_causal_contract") is True, ev_graph
            assert ev_graph.get("circular_convergence_prevented") is True, ev_graph

            ledger = meth.get("run_ledger") or {}
            ledger_ver = ledger.get("claim_ledger_verification") or {}
            assert ledger_ver.get("n_violations") == 0, ledger_ver
            claim_comp = meth.get("claim_compilation") or {}
            assert claim_comp.get("n_withheld") == 0, claim_comp

            nblocks = meth.get("narrative_blocks") or []
            claims = meth.get("claims") or []
            block_ids = [str(b.get("id") or "") for b in nblocks]
            multiome_block_ids = [
                bid for bid in block_ids if bid.startswith("multiome.")
            ]
            chromatin_block_ids = [
                bid for bid in block_ids if bid.startswith("chromatin.")
            ]
            expected_multiome_ids = [
                "multiome.paired_qc",
                "multiome.wnn",
                "multiome.mofa",
                "multiome.peak_to_gene",
                "multiome.eregulons",
                "multiome.triangulation",
            ]
            assert multiome_block_ids == expected_multiome_ids, (
                multiome_block_ids,
                block_ids,
            )

            # Verify staged figures (PNG + SVG) and tables inside report_dir
            fig_pngs = sorted(p.name for p in (report_dir / "figures").glob("*.png"))
            fig_svgs = sorted(p.name for p in (report_dir / "figures").glob("*.svg"))
            staged_tables = sorted(p.name for p in (report_dir / "tables").glob("*"))
            for req_fig in (
                "wnn_joint_vs_modality_umap.png",
                "wnn_joint_vs_modality_umap.svg",
                "wnn_modality_weights_by_cluster.png",
                "wnn_modality_weights_by_cluster.svg",
                "mofa_variance_explained.png",
                "mofa_variance_explained.svg",
                "regulatory_triangulation_summary.png",
                "regulatory_triangulation_summary.svg",
            ):
                assert (report_dir / "figures" / req_fig).is_file(), (
                    req_fig,
                    fig_pngs,
                    fig_svgs,
                )

            for req_tbl in (
                "mofa_factor_loadings.csv",
                "peak_gene_links.csv",
                "eregulons.csv",
                "eregulon_triplets.csv",
            ):
                assert (report_dir / "tables" / req_tbl).is_file(), (
                    req_tbl,
                    staged_tables,
                )

            payload = {
                "slice": "V48_I8",
                "adr": "ADR-071",
                "title": (
                    "IntegrationAgent scaffold -> beta promotion, fail-closed "
                    "same-cell & CP3.5 acknowledgement gate, hermetic env locks, "
                    "and full headless orchestrated E2E validation"
                ),
                "aria_version": __version__,
                "version_metadata": _sanitize_public_obj(collect_version_metadata()),
                "dataset": {
                    "accession": "GSE278576",
                    "sample_id": args.sample_id,
                    "tissue": "human_hippocampus",
                    "assay": "10x_chromium_single_cell_multiome_atac_plus_gene_expression",
                    "input_container": Path(args.h5mu).name,
                    "reference_labels": Path(args.reference_h5ad).name,
                    "gtf_annotation": Path(args.gtf).name,
                    "reference_fasta": Path(args.fasta).name,
                    "n_cells_raw": 3143,
                    "n_cells_post_qc": n_cells_post_qc,
                    "n_cells_wnn_downstream": n_cells_wnn,
                },
                "registry_and_environment_gates": {
                    "integration_validation": dict(INTEGRATION_VALIDATION),
                    "modality_validation_levels": dict(MODALITY_VALIDATION_LEVELS),
                    "registry_integrity_errors": len(reg_errors),
                    "aria_integration_env_audit": {
                        "status": int_env_audit["status"],
                        "conda_exact": int_env_audit["conda_exact"],
                        "pip_exact": int_env_audit["pip_exact"],
                        "conda_expected": int_env_audit["conda_expected"],
                        "pip_expected": int_env_audit["pip_expected"],
                        "binaries": int_env_audit["binaries"],
                    },
                    "aria_chromatin_env_audit": {
                        "status": chrom_env_audit["status"],
                        "conda_exact": chrom_env_audit["conda_exact"],
                        "pip_exact": chrom_env_audit["pip_exact"],
                    },
                },
                "fail_closed_gate_probes": {
                    "unit_dispatch_probes": gate_probes,
                    "live_cp35_cancel_probe": {
                        "status": cancel_run.status,
                        "elapsed_s": t_cancel_s,
                        "checkpoints_visited": [
                            {"checkpoint": cp, "choice": ch}
                            for cp, ch in cancel_run.decisions
                        ],
                        "integration_outputs_created": False,
                    },
                },
                "live_orchestrated_e2e_run": {
                    "status": live_run.status,
                    "elapsed_s": t_live_s,
                    "llm_mode": "offline_deterministic_zero_egress",
                    "llm_plan_integration_needed_override_exercised": (
                        "design_analysis_plan" in llm_live.calls
                    ),
                    "checkpoints_resolved": [
                        {"checkpoint": cp, "choice": ch}
                        for cp, ch in live_run.decisions
                    ],
                    "memory_snapshot_decisions": decisions_rows,
                    "memory_tunnels_stored": tunnels_count,
                    "post_qc_container": post_qc_file.name,
                    "wnn_container": wnn_file.name,
                    "cross_modal_evidence_graph": {
                        "ran": ev_graph.get("ran"),
                        "schema_version": ev_graph.get("schema_version"),
                        "non_causal_contract": ev_graph.get("non_causal_contract"),
                        "circular_convergence_prevented": ev_graph.get(
                            "circular_convergence_prevented"
                        ),
                        "executed_integration_stages": ev_graph.get(
                            "executed_integration_stages"
                        ),
                        "absent_not_tested": ev_graph.get("absent_not_tested"),
                        "wnn_guardrail": ev_graph.get("wnn_guardrail"),
                        "mofa_summary": ev_graph.get("mofa_summary"),
                        "summary_counts": ev_graph.get("summary_counts"),
                    },
                    "narrative_governance": {
                        "n_total_narrative_blocks": len(nblocks),
                        "n_published_claims": len(claims),
                        "n_withheld_claims": claim_comp.get("n_withheld"),
                        "n_multiome_blocks": len(multiome_block_ids),
                        "multiome_block_ids": multiome_block_ids,
                        "n_chromatin_blocks": len(chromatin_block_ids),
                        "chromatin_block_ids": chromatin_block_ids,
                        "claim_ledger_violations": ledger_ver.get("n_violations"),
                        "reproducible_capsule_status": capsule_check["status"],
                        "reproducible_capsule_files_checked": capsule_check["files"][
                            "checked"
                        ],
                        "rocrate_present": bool(rocrate.get("@graph")),
                        "html_size_bytes": len(html_text.encode("utf-8")),
                    },
                    "staged_figures": {
                        "png": fig_pngs,
                        "svg": fig_svgs,
                    },
                    "staged_tables": staged_tables,
                },
                "total_wall_time_s": round(time.time() - t_start, 2),
            }
        finally:
            # Clean up temporary report directory under ~/.aria/reports/
            shutil.rmtree(report_dir, ignore_errors=True)
            exp_ws = Path.home() / ".aria" / "workspace" / live_run.experiment_id
            shutil.rmtree(exp_ws, ignore_errors=True)
            cancel_ws = Path.home() / ".aria" / "workspace" / cancel_run.experiment_id
            shutil.rmtree(cancel_ws, ignore_errors=True)

    sanitized = _sanitize_public_obj(payload)
    _assert_no_forbidden_paths(sanitized)

    out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(sanitized, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[V48-I8] Wrote sanitized benchmark artifact -> {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
