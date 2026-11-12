"""Integrity checks for ARIA agent and script registries."""

from __future__ import annotations

from dataclasses import dataclass
import importlib
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class IntegrityIssue:
    code: str
    message: str
    severity: str = "error"


AGENT_SCRIPT_CONTRACTS = {
    "bulk_rna_agent": {
        "required": (
            "aria/scripts/rna_bulk_de.py",
            "aria/scripts/rna_fastq_qc.py",
            "aria/scripts/rna_align.py",
            "aria/scripts/rna_umi_dedup.py",
            "aria/scripts/rna_quantify.py",
            "aria/scripts/rna_starsolo_bulk.py",
        ),
        "planned": (),
    },
    "scrna_agent": {
        "required": (
            "aria/scripts/rna_qc.py",
            "aria/scripts/rna_concat.py",
            "aria/scripts/rna_integration.py",
            "aria/scripts/rna_clustering.py",
            "aria/scripts/rna_celltypist.py",
            "aria/scripts/rna_apply_cluster_labels.py",
            "aria/scripts/rna_diff_abundance.py",
            "aria/scripts/rna_pseudobulk_de.py",
            "aria/scripts/rna_pathway_per_cluster.py",
            "aria/scripts/rna_cellcomm.py",
            "aria/scripts/rna_trajectory.py",
        ),
        "planned": (),
    },
    "chromatin_agent": {
        "required": (
            "aria/scripts/chromatin_qc.py",
            "aria/scripts/chromatin_peaks.py",
            "aria/scripts/chromatin_scatac_align.py",
            "aria/scripts/chromatin_fragments_to_matrix.py",
            "aria/scripts/chromatin_merge_matrices.py",
            "aria/scripts/chromatin_lsi_clustering.py",
            "aria/scripts/chromatin_diffacc.py",
            "aria/scripts/chromatin_motifs.py",
            "aria/scripts/chromatin_peak_annotation.py",
            "aria/scripts/chromatin_peak_ora.py",
        ),
        "planned": (
            "aria/scripts/chromatin_differential.py",
        ),
    },
    "genome_arch_agent": {
        "required": (
            "aria/scripts/hic_inspect.py",
            "aria/scripts/hic_qc_and_balance.py",
            "aria/scripts/hic_topology.py",
        ),
        "planned": (),
    },
    "integration_agent": {
        "required": (
            "aria/scripts/integration_paired_qc.py",
            "aria/scripts/integration_wnn.py",
            "aria/scripts/integration_wnn_diagnostics.py",
            "aria/scripts/integration_peak2gene.py",
            "aria/scripts/integration_eregulons.py",
            "aria/scripts/integration_mofa.py",
        ),
        "planned": (),
    },
}


def check_registry_integrity(repo_root: Path | None = None) -> list[IntegrityIssue]:
    root = repo_root or Path(__file__).resolve().parents[2]
    issues: list[IntegrityIssue] = []

    from aria.agents.base_agent import BaseAgent
    from aria.agents.orchestrator_agent import (
        AGENT_REGISTRY,
        INTEGRATION_VALIDATION,
        MODALITY_TO_AGENT,
        MODALITY_VALIDATION,
    )
    from aria.utils.script_contracts import (
        CONTRACT_VERSION,
        SCRIPT_CONTRACTS as IPC_SCRIPT_CONTRACTS,
    )

    for agent_name, entry in AGENT_REGISTRY.items():
        try:
            module = importlib.import_module(entry["module"])
            cls = getattr(module, entry["class"])
        except Exception as exc:
            issues.append(IntegrityIssue(
                "agent_import_failed",
                f"{agent_name}: {entry} could not be imported ({exc})",
            ))
            continue
        if not issubclass(cls, BaseAgent):
            issues.append(IntegrityIssue(
                "agent_not_baseagent",
                f"{agent_name}: {entry['class']} is not a BaseAgent subclass",
            ))

    for modality, agent_name in MODALITY_TO_AGENT.items():
        if agent_name not in AGENT_REGISTRY:
            issues.append(IntegrityIssue(
                "modality_unknown_agent",
                f"{modality}: maps to missing agent {agent_name}",
            ))
        meta = MODALITY_VALIDATION.get(modality)
        if not meta or not meta.get("level"):
            issues.append(IntegrityIssue(
                "modality_missing_validation_level",
                f"{modality}: missing validation level",
            ))
        elif meta["level"] == "scaffold" and meta.get("dispatch_enabled", True):
            issues.append(IntegrityIssue(
                "scaffold_dispatch_enabled",
                f"{modality}: scaffold modality is dispatch-enabled",
                severity="warning",
            ))

    int_level = INTEGRATION_VALIDATION.get("level")
    int_enabled = INTEGRATION_VALIDATION.get("dispatch_enabled", True)
    if not int_level:
        issues.append(IntegrityIssue(
            "integration_missing_validation_level",
            "integration_agent: missing validation level",
        ))
    elif int_level == "scaffold" and int_enabled:
        issues.append(IntegrityIssue(
            "scaffold_dispatch_enabled",
            "integration_agent: scaffold agent is dispatch-enabled",
        ))
    elif int_level == "beta":
        if INTEGRATION_VALIDATION.get("requires_ack") is not True:
            issues.append(IntegrityIssue(
                "integration_beta_missing_ack_gate",
                "integration_agent: beta tier must declare requires_ack=True (ADR-056/ADR-071)",
            ))
        if INTEGRATION_VALIDATION.get("production_promoted") is not False:
            issues.append(IntegrityIssue(
                "integration_beta_production_flag_invalid",
                "integration_agent: beta tier must declare production_promoted=False (ADR-056)",
            ))
        if INTEGRATION_VALIDATION.get("scope") != "paired_same_cell_rna_atac":
            issues.append(IntegrityIssue(
                "integration_beta_scope_invalid",
                "integration_agent: beta tier scope must be 'paired_same_cell_rna_atac' (ADR-071)",
            ))
    elif int_level == "production":
        issues.append(IntegrityIssue(
            "integration_premature_production_tier",
            "integration_agent: production promotion is blocked by ADR-056 governance criteria",
        ))

    for agent_name, contract in AGENT_SCRIPT_CONTRACTS.items():
        for script in _scripts(contract.get("required", ())):
            if not (root / script).exists():
                issues.append(IntegrityIssue(
                    "required_script_missing",
                    f"{agent_name}: required script missing: {script}",
                ))
        for script in _scripts(contract.get("planned", ())):
            if (root / script).exists():
                continue
            if agent_name != "chromatin_agent":
                issues.append(IntegrityIssue(
                    "planned_script_missing_unblocked",
                    f"{agent_name}: planned script missing without scaffold gate: {script}",
                ))

    for script, contract in IPC_SCRIPT_CONTRACTS.items():
        if not (root / script).exists():
            issues.append(IntegrityIssue(
                "contract_script_missing",
                f"IPC contract references missing script: {script}",
            ))
        if contract.contract_version != CONTRACT_VERSION:
            issues.append(IntegrityIssue(
                "contract_version_mismatch",
                (
                    f"{script}: contract version {contract.contract_version} "
                    f"does not match runtime {CONTRACT_VERSION}"
                ),
            ))

    return issues


def _scripts(paths: Iterable[str]) -> tuple[str, ...]:
    return tuple(str(path) for path in paths)
