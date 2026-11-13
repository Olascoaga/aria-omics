"""ADR-079 — Paired Same-Cell Multiome Production Promotion End-to-End Suite.

Verifies end-to-end production readiness for ``scope="paired_same_cell_rna_atac"``:
  1. Governance & ``IntegrationAuditAgent`` readiness (``green``/``auto`` at
     ``production`` vs ``yellow``/``requires_ack`` on ``<200`` cells, partial
     barcode overlap, or 1-donor-per-condition confounding vs ``red``/``blocked``
     on unpaired/diagonal/bulk cross-modal inputs while single-assay lanes stay
     unaffected).
  2. Live ``DataAuditAgent`` -> ``SetupAgent`` -> ``OrchestratorAgent``
     (autonomous ``green`` dispatch without CP3.5 pause) -> ``IntegrationAgent``
     -> ``NarrativeAgent`` (``C3`` ``ledger://integration/preprocessing ==
     "executed"``, ``claim_ledger_verification.n_violations == 0``) on a Cell
     Ranger ARC dual-modality MEX directory (``10x_multiome_mtx``).
  3. Live multi-donor paired ``.h5mu`` cohort through ``IntegrationAgent``
     (``integration_merge_cohort.py`` -> ``integration_paired_qc.py`` ->
     donor-aware Harmony ``integration_wnn.py`` -> donor/condition ANOVA
     ``eta^2`` ``integration_mofa.py``) plus ``A4`` SHA-256 stage-manifest
     instant resume (``resumed=True``) and parameter-change invalidation.
  4. Quantitative MOFA+ batch/donor factor exclusion & ``IntegrationNarrator`` /
     ``IntegrationAgent._resolve_conflicts`` disclosures with zero causal verbs.
  5. Empirical receipt invariants on ``adr079_multidonor_multiome_validation.json``
     and AST anti-fabrication / causal-verb audit across all ADR-079 modules.
"""

from __future__ import annotations

import ast
import copy
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

_CAUSAL_VERBS = {"drives", "regulates", "activates", "represses", "causes", "controls"}


def _env_python(env_name: str) -> Path | None:
    exe = Path(sys.executable).resolve()
    for parent in exe.parents:
        cand = parent / "envs" / env_name / "bin" / "python"
        if cand.is_file():
            return cand
    home_cand = Path.home() / "anaconda3" / "envs" / env_name / "bin" / "python"
    return home_cand if home_cand.is_file() else None


def _require_multiome_envs() -> tuple[Path, Path]:
    chrom_py = _env_python("aria-chromatin-env")
    integ_py = _env_python("aria-integration-env")
    if chrom_py is None or integ_py is None:
        pytest.skip("aria-chromatin-env and aria-integration-env are required for live multiome E2E")
    return chrom_py, integ_py


def _write_arc_dual_mex(target_dir: Path, n_cells: int = 240, seed: int = 7) -> Path:
    """Write a realistic 10x Cell Ranger ARC filtered_feature_bc_matrix directory
    containing both ``Gene Expression`` and ``Peaks`` features over shared barcodes."""
    target_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    bases = ("A", "C", "G", "T")
    barcodes: list[str] = []
    idx = 0
    while len(barcodes) < n_cells:
        val = idx
        chars = []
        for _ in range(16):
            chars.append(bases[val & 3])
            val >>= 2
        bc = "".join(chars) + "-1"
        if bc not in barcodes:
            barcodes.append(bc)
        idx += 17

    n_genes = 120
    n_peaks = 160
    genes = [f"MT-CO{i + 1}" if i < 3 else f"GENE_{i:03d}" for i in range(n_genes)]
    peaks = [f"chr1:{100_000 + i * 2_000}-{100_600 + i * 2_000}" for i in range(n_peaks)]

    half = n_cells // 2
    rna_mat = rng.poisson(1.4, size=(n_cells, n_genes)).astype(np.int64)
    rna_mat[:half, 3:35] += rng.poisson(6.0, size=(half, 32))
    rna_mat[half:, 35:70] += rng.poisson(6.0, size=(n_cells - half, 35))
    rna_mat[:, 0:3] = rng.poisson(0.4, size=(n_cells, 3))
    # Guarantee every cell passes paired QC floors (>=200 RNA genes? paired_qc adaptive/default)
    rna_mat[:, 70:] += 1

    atac_mat = rng.poisson(1.1, size=(n_cells, n_peaks)).astype(np.int64)
    atac_mat[:half, :50] += rng.poisson(5.0, size=(half, 50))
    atac_mat[half:, 50:105] += rng.poisson(5.0, size=(n_cells - half, 55))
    atac_mat[:, 105:] += 1

    full = np.hstack([rna_mat, atac_mat])  # shape (n_cells, n_features)

    with gzip.open(target_dir / "barcodes.tsv.gz", "wt", encoding="utf-8") as fh:
        for bc in barcodes:
            fh.write(f"{bc}\n")

    with gzip.open(target_dir / "features.tsv.gz", "wt", encoding="utf-8") as fh:
        for j, g in enumerate(genes):
            fh.write(f"ENSG{j:08d}\t{g}\tGene Expression\tchr1\t{1000 + j * 100}\t{1500 + j * 100}\n")
        for j, pk in enumerate(peaks):
            chrom, span = pk.split(":")
            st, en = span.split("-")
            fh.write(f"{pk}\t{pk}\tPeaks\t{chrom}\t{st}\t{en}\n")

    nz_cells, nz_feats = np.nonzero(full)
    nnz = int(len(nz_cells))
    with gzip.open(target_dir / "matrix.mtx.gz", "wt", encoding="utf-8") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        fh.write(f"{n_genes + n_peaks} {n_cells} {nnz}\n")
        for c_i, f_i in zip(nz_cells, nz_feats):
            fh.write(f"{int(f_i) + 1} {int(c_i) + 1} {int(full[c_i, f_i])}\n")

    return target_dir


def _build_two_donor_h5mu_cohort(chrom_py: Path, out_dir: Path) -> tuple[Path, Path]:
    """Create 2 synthetic paired .h5mu donor files in aria-chromatin-env with
    distinct overlapping ATAC peak sets, shared 2-cluster biology, and donor shift."""
    out_dir.mkdir(parents=True, exist_ok=True)
    d1 = out_dir / "donorA_paired.h5mu"
    d2 = out_dir / "donorB_paired.h5mu"
    code = r"""
import sys
from pathlib import Path
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp

d1_path, d2_path = Path(sys.argv[1]), Path(sys.argv[2])
rng = np.random.default_rng(42)

def make_donor(donor_id, condition, peak_shift, donor_axis_sign, seed_off):
    r = np.random.default_rng(42 + seed_off)
    n_cells = 130
    half = n_cells // 2
    bcs = [f"BC{i:04d}AA-1" for i in range(n_cells)]
    genes = [f"MT-ND{i+1}" if i < 2 else f"G_{i:03d}" for i in range(140)]
    peaks = [
        f"chr1:{200_000 + i * 1_500 + (peak_shift if i % 2 == 0 else 0)}-{200_700 + i * 1_500 + (peak_shift if i % 2 == 0 else 0)}"
        for i in range(180)
    ]
    if seed_off == 2:
        peaks.extend([f"chr2:{50_000 + j * 2_000}-{50_600 + j * 2_000}" for j in range(5)])

    rna = r.poisson(1.2, size=(n_cells, len(genes))).astype(np.int64)
    rna[:half, 2:40] += r.poisson(8.5, size=(half, 38))
    rna[half:, 40:80] += r.poisson(8.5, size=(n_cells - half, 40))
    if donor_axis_sign > 0:
        rna[:, 80:92] += r.poisson(1.6, size=(n_cells, 12))
    else:
        rna[:, 92:104] += r.poisson(1.6, size=(n_cells, 12))
    rna[:, 130:] += 1

    atac = r.poisson(1.1, size=(n_cells, len(peaks))).astype(np.int64)
    atac[:half, :55] += r.poisson(7.5, size=(half, 55))
    atac[half:, 55:110] += r.poisson(7.5, size=(n_cells - half, 55))
    if donor_axis_sign > 0:
        atac[:, 110:122] += r.poisson(1.5, size=(n_cells, 12))
    else:
        atac[:, 122:134] += r.poisson(1.5, size=(n_cells, 12))
    atac[:, 170:] += 1

    obs = pd.DataFrame(
        {
            "donor_id": donor_id,
            "condition": condition,
            "cell_type": ["TypeA"] * half + ["TypeB"] * (n_cells - half),
        },
        index=bcs,
    )
    rna_ad = ad.AnnData(sp.csr_matrix(rna), obs=obs.copy(), var=pd.DataFrame(index=genes))
    atac_ad = ad.AnnData(sp.csr_matrix(atac), obs=obs.copy(), var=pd.DataFrame(index=peaks))
    rna_ad.layers["counts"] = rna_ad.X.copy()
    atac_ad.layers["counts"] = atac_ad.X.copy()
    m = md.MuData({"rna": rna_ad, "atac": atac_ad})
    m.obs = obs.copy()
    return m

make_donor("donorA", "ctrl", 0, +1, 1).write_h5mu(str(d1_path))
make_donor("donorB", "ctrl", 120, -1, 2).write_h5mu(str(d2_path))
"""
    env = dict(os.environ)
    env["NUMBA_CACHE_DIR"] = str(out_dir / "numba")
    proc = subprocess.run([str(chrom_py), "-c", code, str(d1), str(d2)], capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    return d1, d2


def _make_offline_integration_agent(tmp_path: Path):
    from aria.agents.integration_agent import IntegrationAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.environment_manager import EnvironmentManager

    mem = ARIAMemory(db_path=str(tmp_path / "integ.sqlite"))
    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.memory = mem
    agent.env = EnvironmentManager()
    agent.llm = SimpleNamespace(
        complete=lambda *a, **kw: "Descriptive cross-modal summary.",
        complete_heavy=lambda *a, **kw: "Descriptive cross-modal summary.",
        complete_medium=lambda *a, **kw: "Descriptive cross-modal summary.",
        complete_light=lambda *a, **kw: "Descriptive cross-modal summary.",
    )
    agent.advisor = SimpleNamespace(
        advise_wnn_k=lambda **kw: SimpleNamespace(
            decision_id="wnn_k",
            chosen_value=15,
            justification="Deterministic k=15 for paired cohort.",
            candidates=[],
            chosen_by="advisor",
            approved_by_user=False,
        ),
        advise_mofa_factors=lambda **kw: SimpleNamespace(
            decision_id="mofa_factors",
            chosen_value=5,
            justification="Deterministic 5 factors.",
            candidates=[],
            chosen_by="advisor",
            approved_by_user=False,
        ),
    )
    agent.publish_status = lambda *a, **kw: None
    agent.publish_finding = lambda *a, **kw: None
    agent.publish_escalation = lambda *a, **kw: None
    agent.publish_blocking_escalation = lambda **kw: ("ack", {"choice": kw["options"][0]})
    return agent, mem


# ══════════════════════════════════════════════════════════════════════════════
# 1. Production Governance, Docker Validation & CP3.5 Risk Escalation
# ══════════════════════════════════════════════════════════════════════════════


def test_adr079_production_governance_and_risk_escalation_matrix(tmp_path, monkeypatch):
    from aria.agents.audit_agent import AuditAgent
    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.orchestrator_agent import (
        INTEGRATION_VALIDATION,
        MODALITY_VALIDATION,
        OrchestratorAgent,
    )
    from aria.memory.memory import ARIAMemory
    from aria.utils.registry_integrity import AGENT_SCRIPT_CONTRACTS, check_registry_integrity
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    assert INTEGRATION_VALIDATION == {
        "level": "production",
        "dispatch_enabled": True,
        "requires_ack": False,
        "scope": "paired_same_cell_rna_atac",
        "production_promoted": True,
        "promotion_policy_adr": "ADR-056",
        "promotion_adr": "ADR-079",
        "reason": INTEGRATION_VALIDATION["reason"],
    }
    assert "integration" not in MODALITY_VALIDATION
    assert check_registry_integrity() == []

    all_integration_scripts = (
        *AGENT_SCRIPT_CONTRACTS["integration_agent"]["required"],
        *AGENT_SCRIPT_CONTRACTS["integration_agent"]["optional"],
    )
    assert len(all_integration_scripts) == 8
    for rel in all_integration_scripts:
        assert SCRIPT_CONTRACTS[rel].validation_level == "production", rel

    dockerfile_txt = (ROOT / "docker" / "Dockerfile.integration").read_text(encoding="utf-8")
    assert "ARIA_IMAGE_VALIDATION=validated" in dockerfile_txt

    # 1. Balanced 3x3 multi-donor paired cohort -> green / auto across both scATAC and integration
    balanced_ctx = {
        "genome": "hg38",
        "modalities": {"scATAC": [f"/v/d{i}.h5mu" for i in range(1, 7)]},
        "assay_detections": [
            {"path": f"/v/d{i}.h5mu", "modality": "scATAC", "evidence": {"format": "h5mu", "paired_rna": True}}
            for i in range(1, 7)
        ],
        "design": {
            "groups": {"ctrl": ["d1", "d2", "d3"], "treat": ["d4", "d5", "d6"]},
            "sample_metadata": {
                "d1": {"donor_id": "d1", "condition": "ctrl"},
                "d2": {"donor_id": "d2", "condition": "ctrl"},
                "d3": {"donor_id": "d3", "condition": "ctrl"},
                "d4": {"donor_id": "d4", "condition": "treat"},
                "d5": {"donor_id": "d5", "condition": "treat"},
                "d6": {"donor_id": "d6", "condition": "treat"},
            },
        },
    }
    mat_green = build_capability_matrix(
        balanced_ctx,
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert mat_green["integration_readiness"]["status"] == "green"
    assert mat_green["integration_readiness"]["dispatch_policy"] == "auto"
    assert mat_green["integration_readiness"]["requires_ack_for_dispatch"] is False

    # 2. 1-donor-per-condition confounded multiome cohort -> yellow / requires_ack at CP3.5
    confounded_ctx = {
        "genome": "hg38",
        "modalities": {"scATAC": ["/v/d1.h5mu", "/v/d2.h5mu"]},
        "assay_detections": balanced_ctx["assay_detections"][:2],
        "design": {
            "groups": {"ctrl": ["d1"], "treat": ["d2"]},
            "sample_metadata": {
                "d1": {"donor_id": "d1", "condition": "ctrl"},
                "d2": {"donor_id": "d2", "condition": "treat"},
            },
        },
    }
    mat_conf = build_capability_matrix(
        confounded_ctx,
        modality_validation=MODALITY_VALIDATION,
        integration_validation=INTEGRATION_VALIDATION,
    )
    assert mat_conf["integration_readiness"]["status"] == "yellow"
    assert mat_conf["integration_readiness"]["requires_ack_for_dispatch"] is True
    assert any(
        f["check"] == "integration_donor_confounded_design"
        for f in mat_conf["integration_readiness"]["findings"]
    )

    # 3. Orchestrator _gate_and_dispatch dispatches green cohort WITHOUT CP3.5 pause,
    #    while pausing confounded yellow cohort at CP3.5!
    mem = ARIAMemory(db_path=str(tmp_path / "gov.sqlite"))
    try:
        mem.create_wing("exp_green", name="green", organism="Homo sapiens", genome="hg38")
        mem.create_wing("exp_conf", name="conf", organism="Homo sapiens", genome="hg38")
        orch = OrchestratorAgent.__new__(OrchestratorAgent)
        orch.memory = mem
        orch.llm = object()
        orch._sessions = {}
        orch._pending_checkpoints = {}
        orch._pending_dispatch = {}
        orch._experiment_plans = {"exp_green": {"intent": {}}, "exp_conf": {"intent": {}}}
        orch._agent_results = {}
        orch._log_handlers = {}
        escalations: list[dict[str, Any]] = []
        dispatched: list[Any] = []
        orch.publish_status = lambda *a, **kw: None
        orch.publish_escalation = lambda **kw: escalations.append(kw)

        class _NoopThread:
            def __init__(self, target, args, daemon):
                self.args = args

            def start(self):
                dispatched.append(self.args)

        monkeypatch.setattr("aria.agents.orchestrator_agent.threading.Thread", _NoopThread)

        res_green = orch._gate_and_dispatch(
            "exp_green",
            {"steps": [], "integration_needed": True},
            copy.deepcopy(balanced_ctx),
        )
        assert res_green["status"] == "analysis_running"
        assert escalations == []
        assert len(dispatched) == 1

        res_conf = orch._gate_and_dispatch(
            "exp_conf",
            {"steps": [], "integration_needed": True},
            copy.deepcopy(confounded_ctx),
        )
        assert res_conf["status"] == "audit_blocking"
        assert len(escalations) == 1
        assert escalations[0]["checkpoint"] == 3.5
    finally:
        mem.close()


# ══════════════════════════════════════════════════════════════════════════════
# 2. Live Cell Ranger ARC Dual MEX -> DataAudit -> Integration -> Narrative E2E
# ══════════════════════════════════════════════════════════════════════════════


def test_live_cellranger_arc_mex_autonomous_integration_and_c3_narrative_e2e(tmp_path):
    _require_multiome_envs()
    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION, MODALITY_VALIDATION
    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.setup_agent import SetupAgent
    from aria.memory.memory import ARIAMemory

    arc_dir = _write_arc_dual_mex(tmp_path / "sample_arc" / "filtered_feature_bc_matrix", n_cells=240)

    mem = ARIAMemory(db_path=str(tmp_path / "arc_e2e.sqlite"))
    try:
        da = DataAuditAgent(mem, llm=SimpleNamespace(complete_light=lambda *a, **kw: "ok"))
        da.publish_status = lambda *a, **kw: None
        da.publish_escalation = lambda *a, **kw: None
        audit_out = da.run(
            "exp_arc_mex",
            {
                "data_dir": str(arc_dir.parent),
                "user_question": "Characterize joint RNA+ATAC cell states",
            },
        )
        exp_ctx = audit_out["exp_context"]
        assert "scRNA" not in exp_ctx["modalities"], (
            "Cell Ranger ARC dual MEX barcodes/features must not leak into scRNA"
        )
        assert "scATAC" in exp_ctx["modalities"]
        assert exp_ctx["multiome_contract"]["same_cell"] is True
        assert exp_ctx["multiome_contract"]["object_type"] == "paired_mudata"
        assert any(
            (d.get("evidence") or {}).get("format") == "10x_multiome_mtx"
            for d in exp_ctx.get("assay_detections", [])
        )

        setup = SetupAgent.__new__(SetupAgent)
        envs = setup._needed_envs(
            exp_ctx["modalities"],
            has_fastq=False,
            exp_ctx=exp_ctx,
        )
        assert "aria-chromatin-env" in envs
        assert "aria-integration-env" in envs

        cap = build_capability_matrix(
            exp_ctx,
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )
        assert cap["integration_readiness"]["status"] == "green"
        assert cap["integration_readiness"]["requires_ack_for_dispatch"] is False

        integ_agent, integ_mem = _make_offline_integration_agent(tmp_path)
        try:
            integ_res = integ_agent.run(
                "exp_arc_mex",
                {
                    "exp_context": {
                        **exp_ctx,
                        "enable_mofa": False,
                        "paired_qc_params": {
                            "min_genes": 10,
                            "min_peaks": 10,
                            "min_rna_counts": 20,
                            "min_atac_counts": 20,
                        },
                        "wnn_params": {"k": 15, "n_pcs": 20, "n_lsi_comps": 20},
                    },
                    "biological_intent": {"user_question": "Joint RNA+ATAC states"},
                },
            )
        finally:
            integ_mem.close()

        assert integ_res["status"] == "done", integ_res
        findings = integ_res["findings"]
        assert findings["preprocessing"]["status"] == "success"
        assert Path(findings["preprocessing"]["output_path"]).is_file()
        assert findings["paired_qc"]["status"] == "success"
        assert findings["wnn"]["status"] == "success"
        assert findings["wnn"]["n_joint_clusters"] >= 2

        narr_agent = NarrativeAgent(memory=mem, llm=None)
        narr_agent.reports_dir = tmp_path / "reports"
        narr_agent.reports_dir.mkdir(parents=True, exist_ok=True)
        narr_res = narr_agent.run(
            experiment_id="exp_arc_mex",
            context={
                "exp_context": {
                    **exp_ctx,
                    "organism": "Homo sapiens",
                    "genome": "hg38",
                    "reproducible_mode": True,
                },
                "biological_intent": {"summary": "Cell Ranger ARC paired multiome"},
                "agent_results": {"integration_agent": integ_res},
                "findings": [],
            },
        )
        assert narr_res["status"] == "done", narr_res
        report_dir = Path(narr_res["report_path"]).parent
        meth = json.loads((report_dir / "methodology.json").read_text(encoding="utf-8"))
        assert meth["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0
        pre_entries = [
            e
            for e in meth["run_ledger"]["entries"]
            if e["modality"] == "integration" and e["analysis"] == "preprocessing"
        ]
        assert len(pre_entries) == 1
        assert pre_entries[0]["status"] == "executed"
    finally:
        mem.close()


# ══════════════════════════════════════════════════════════════════════════════
# 3. Live Multi-Donor Cohort Merger + Harmony WNN + MOFA+ ANOVA eta^2 + A4 Resume
# ══════════════════════════════════════════════════════════════════════════════


def test_live_multidonor_cohort_harmony_wnn_mofa_eta2_and_a4_resume(tmp_path):
    chrom_py, _ = _require_multiome_envs()
    d1, d2 = _build_two_donor_h5mu_cohort(chrom_py, tmp_path / "donors")

    exp_ctx = {
        "genome": "hg38",
        "organism": "Homo sapiens",
        "modalities": {"scATAC": [str(d1), str(d2)]},
        "assay_detections": [
            {"path": str(d1), "modality": "scATAC", "evidence": {"format": "h5mu", "paired_rna": True}},
            {"path": str(d2), "modality": "scATAC", "evidence": {"format": "h5mu", "paired_rna": True}},
        ],
        "enable_mofa": True,
        "batch_key": "donor_id",
        "condition_col": "condition",
        "paired_qc_params": {
            "min_genes": 10,
            "min_peaks": 10,
            "min_rna_counts": 20,
            "min_atac_counts": 20,
        },
        "wnn_params": {
            "k": 15,
            "n_pcs": 20,
            "n_lsi_comps": 20,
            "harmonize_batches": True,
            "audit_group_col": "cell_type",
        },
        "mofa_params": {
            "n_factors": 5,
            "n_top_genes": 120,
            "n_top_peaks": 150,
            "max_iter": 150,
            "write_h5mu": False,
        },
    }

    integ_agent, integ_mem = _make_offline_integration_agent(tmp_path)
    try:
        res1 = integ_agent.run(
            "exp_cohort_1",
            {
                "exp_context": copy.deepcopy(exp_ctx),
                "biological_intent": {"user_question": "Multi-donor joint states and factors"},
            },
        )
        assert res1["status"] == "done", res1
        f1 = res1["findings"]
        assert f1["preprocessing"]["mode"] == "paired_cohort_merge"
        assert f1["preprocessing"]["n_samples"] == 2
        assert f1["preprocessing"]["n_cells_total"] == 260
        assert f1["preprocessing"]["n_consensus_peaks"] > 180  # union across shifted donor peaks
        assert f1["wnn"]["status"] == "success", f1["wnn"]
        assert f1["mofa"]["status"] == "success", f1["mofa"]

        bi = f1["wnn"]["batch_integration"]
        assert bi["applied"] is True
        assert bi["method"] == "harmony"
        assert bi["batch_key"] == "donor_id"
        assert bi["wnn_batch_ilisi_median"] > 1.0
        assert isinstance(bi["quality_assessment"], dict)

        assert f1["mofa"]["batch_evaluated"] is True
        assert f1["mofa"]["batch_key"] == "donor_id"
        assert any(f.get("batch_eta2") is not None for f in (f1["mofa"].get("top_factors") or []))

        # Second identical call hits A4 SHA-256 stage manifests across cohort_merge, paired_qc, wnn, and mofa!
        res2 = integ_agent.run(
            "exp_cohort_2",
            {
                "exp_context": copy.deepcopy(exp_ctx),
                "biological_intent": {"user_question": "Multi-donor joint states and factors"},
            },
        )
        assert res2["status"] == "done", res2
        f2 = res2["findings"]
        assert f2["preprocessing"].get("resumed") is True
        assert f2["paired_qc"].get("resumed") is True
        assert f2["wnn"].get("resumed") is True
        assert f2["mofa"].get("resumed") is True

        # Changing WNN resolution invalidates WNN A4 stage manifest while keeping upstream cache hits!
        exp_ctx_mut = copy.deepcopy(exp_ctx)
        exp_ctx_mut["wnn_params"]["resolution"] = 0.8
        res3 = integ_agent.run(
            "exp_cohort_3",
            {
                "exp_context": exp_ctx_mut,
                "biological_intent": {"user_question": "Multi-donor joint states and factors"},
            },
        )
        assert res3["status"] == "done", res3
        assert res3["findings"]["paired_qc"].get("resumed") is True
        assert not res3["findings"]["wnn"].get("resumed")
    finally:
        integ_mem.close()


# ══════════════════════════════════════════════════════════════════════════════
# 4. Conflict Resolution & Narrator Disclosures (Zero Causal Language)
# ══════════════════════════════════════════════════════════════════════════════


def test_wnn_overcorrection_and_mofa_batch_confounding_disclosures_are_non_causal():
    from aria.agents.integration_agent import IntegrationAgent
    from aria.agents.narrative.narrators.integration import IntegrationNarrator

    agent = IntegrationAgent.__new__(IntegrationAgent)
    synthetic_findings = {
        "preprocessing": {
            "status": "success",
            "mode": "paired_matrices",
            "output_path": "/v/assembled.h5mu",
            "n_samples": 1,
            "n_cells_total": 500,
            "barcode_translation": {"mode": "10x_arc_whitelist_pair", "orientation": "reverse_complement"},
        },
        "wnn": {
            "status": "success",
            "output_path": "/v/wnn.h5mu",
            "n_cells": 500,
            "n_joint_clusters": 4,
            "n_rna_only_clusters": 4,
            "n_discordant_clusters": 0,
            "mean_rna_weight": 0.52,
            "mean_atac_weight": 0.48,
            "batch_integration": {
                "batch_key": "donor_id",
                "n_batches": 2,
                "harmonized": True,
                "method": "harmony",
                "wnn_batch_ilisi_median": 1.85,
                "wnn_batch_ilisi_norm_median": 0.85,
                "quality_assessment": {
                    "possible_overcorrection": True,
                    "warning": "Harmony batch integration reduced reference group purity below 0.70.",
                },
            },
        },
        "mofa": {
            "status": "success",
            "output_path": "/v/mofa.hdf5",
            "n_cells": 500,
            "n_factors": 4,
            "variance_explained": {"rna": 0.22, "atac": 0.15},
            "batch_factor_flag": True,
            "batch_factor_ids": [2],
            "batch_screening": {
                "ran": True,
                "batch_key": "donor_id",
                "condition_col": "condition",
                "batch_condition_confounded": True,
                "n_batches": 2,
            },
            "excluded_factors": [
                {"factor_id": 2, "exclusion_reasons": ["batch_confounded_eta2=0.81"]}
            ],
        },
    }

    conflicts_res = agent._resolve_conflicts("exp_disc", {}, {}, synthetic_findings)
    conflict_types = {c["type"] for c in conflicts_res["conflicts"]}
    assert "wnn_harmony_possible_overcorrection" in conflict_types
    assert "mofa_confounder_factor_excluded" in conflict_types
    assert "mofa_batch_confounded_factors" in conflict_types

    blocks = IntegrationNarrator().collect(
        "integration_agent",
        {"status": "done", "findings": synthetic_findings},
        {"agent_results": {}},
    )
    by_id = {b.id: b for b in blocks}
    assert "multiome.preprocessing" in by_id
    assert "10x_arc_whitelist_pair (reverse_complement)" in by_id["multiome.preprocessing"].claim

    all_words: set[str] = set()
    for b in blocks:
        for txt in [b.claim or "", *(c.text for c in b.caveats)]:
            for token in txt.lower().replace(".", " ").replace(",", " ").split():
                all_words.add(token)
    assert all_words.isdisjoint(_CAUSAL_VERBS), f"Causal verbs leaked: {all_words & _CAUSAL_VERBS}"


# ══════════════════════════════════════════════════════════════════════════════
# 5. Empirical Receipt & AST Anti-Fabrication Invariants
# ══════════════════════════════════════════════════════════════════════════════


def test_adr079_multidonor_empirical_receipt_and_no_fabrication_scan():
    receipt_path = (
        ROOT
        / "docs"
        / "benchmark_results"
        / "multiome_concordance"
        / "adr079_multidonor_multiome_validation.json"
    )
    if receipt_path.is_file():
        raw_txt = receipt_path.read_text(encoding="utf-8")
        for forbidden in ("/home/", "/tmp/", "/Users/", "/private/"):
            assert forbidden not in raw_txt, f"Absolute path {forbidden!r} leaked into receipt"
        receipt = json.loads(raw_txt)
        assert receipt["adr"] == "ADR-079"
        assert receipt["dataset"]["accession"] == "GSE278576"
        assert receipt["dataset"]["n_donors"] == 4
        assert receipt["stage_1_cohort_merge"]["a4_resume_verified"] is True
        assert receipt["stage_2_paired_qc"]["a4_resume_verified"] is True
        assert receipt["stage_3_donor_aware_wnn"]["harmony_integrated"]["harmonized"] is True
        assert receipt["stage_3_donor_aware_wnn"]["harmony_integrated"]["a4_resume_verified"] is True
        assert receipt["stage_4_donor_screened_mofa"]["a4_resume_verified"] is True
        p2g = receipt["stage_5_cis_peak2gene_vs_barcode_permuted_null"]
        assert p2g["real_paired"]["n_links"] >= 10
        assert p2g["barcode_shuffled_null_control"]["empirical_false_discovery_rate_at_q005"] <= 0.05
        assert p2g["barcode_shuffled_null_control"]["per_test_false_positive_rate"] < 1e-4
        assert p2g["barcode_shuffled_null_control"]["derangement_fixed_points"] == 0
        gov = receipt["stage_6_governance_and_c3_ledger"]
        assert gov["c3_preprocessing_stage_status"] == "executed"
        assert gov["claim_ledger_violations"] == 0

    spec = importlib.util.spec_from_file_location(
        "_fab_guard_adr079",
        ROOT / "tests" / "test_no_fabrication_guard.py",
    )
    assert spec and spec.loader
    fab_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fab_mod)

    adr079_targets = (
        "aria/scripts/integration_assemble_multiome.py",
        "aria/scripts/integration_merge_cohort.py",
        "aria/scripts/integration_paired_qc.py",
        "aria/scripts/integration_wnn.py",
        "aria/scripts/integration_mofa.py",
        "aria/scripts/integration_peak2gene.py",
        "aria/scripts/integration_eregulons.py",
        "aria/utils/mudata_io.py",
        "aria/utils/multiome_contracts.py",
        "aria/agents/integration_agent.py",
        "aria/agents/modality_audit.py",
        "aria/agents/orchestrator_agent.py",
    )
    for rel in adr079_targets:
        src = (ROOT / rel).read_text(encoding="utf-8")
        ast.parse(src)
        hits = fab_mod.scan_source(rel, src)
        assert hits == [], f"No-fabrication guard flagged {rel}: {hits}"
