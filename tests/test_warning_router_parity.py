"""Cross-modal warning-routing parity and manual-loop retirement tests (W-WARN §9.4).

Verifies:
  1. The manual ``for ... in <stage>.get("warnings")`` loops in ``HiCNarrator``
     and ``ChromatinNarrator`` are retired (no duplicated per-stage warning
     plumbing inside narrators).
  2. Every caveat that those retired loops previously produced is still
     produced through ``NarrativeRegistry.collect_blocks`` with the exact same
     ``severity="warning"`` and sanitized text, without altering confidence or
     status.
  3. All 5 modality narrators (``bulk_rna_agent``, ``scrna_agent``,
     ``chromatin_agent``, ``genome_arch_agent``, ``integration_agent``) achieve
     cross-modal parity for both free-text stage warnings and structured
     degradation signals.
"""

from __future__ import annotations

import ast
from pathlib import Path
import pytest

from aria.agents.narrative.narrators import (
    BulkRnaNarrator,
    ChromatinNarrator,
    HiCNarrator,
    IntegrationNarrator,
    ScrnaNarrator,
)
from aria.agents.narrative.registry import registry_with

_REPO_ROOT = Path(__file__).resolve().parents[1]
_HIC_PATH = _REPO_ROOT / "aria" / "agents" / "narrative" / "narrators" / "hic.py"
_CHROM_PATH = _REPO_ROOT / "aria" / "agents" / "narrative" / "narrators" / "chromatin.py"


def _full_registry():
    return registry_with((
        ScrnaNarrator(),
        BulkRnaNarrator(),
        ChromatinNarrator(),
        HiCNarrator(),
        IntegrationNarrator(),
    ))


def _manual_warning_for_loops(path: Path) -> list[int]:
    """Return line numbers of any ``for ... in *.get("warnings")`` loops in ``path``."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    lines: list[int] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.For):
            continue
        for sub in ast.walk(node.iter):
            if (
                isinstance(sub, ast.Call)
                and isinstance(sub.func, ast.Attribute)
                and sub.func.attr == "get"
                and sub.args
                and isinstance(sub.args[0], ast.Constant)
                and sub.args[0].value == "warnings"
            ):
                lines.append(node.lineno)
                break
    return lines


def test_manual_warning_loops_are_retired_in_hic_and_chromatin_narrators():
    hic_loops = _manual_warning_for_loops(_HIC_PATH)
    chrom_loops = _manual_warning_for_loops(_CHROM_PATH)
    assert hic_loops == [], f"Unretired manual warning loop(s) in hic.py at line(s): {hic_loops}"
    assert chrom_loops == [], f"Unretired manual warning loop(s) in chromatin.py at line(s): {chrom_loops}"


def test_hic_all_five_stages_warning_parity_through_registry():
    agent_result = {
        "status": "done",
        "findings": {
            "balancing": {
                "status": "success",
                "n_samples": 4,
                "n_valid_pairs": 120000,
                "cis_trans_ratio": 3.2,
                "pass_qc": True,
                "warnings": ["ICE weight convergence note on chrM"],
            },
            "topology": {
                "compartments": {
                    "status": "success",
                    "pct_A": 51.0,
                    "pct_B": 49.0,
                    "resolution_bp": 50000,
                    "pc1_validated": True,
                    "phasing_method": "gc_content",
                    "warnings": ["Compartment phasing skipped unplaced contig"],
                },
                "tads": {
                    "status": "success",
                    "n_tads": 42,
                    "median_size_kb": 450.0,
                    "warnings": ["TAD insulation boundary prominence relaxed on short arm"],
                },
                "loops": {
                    "status": "success",
                    "n_loops": 15,
                    "resolution_bp": 10000,
                    "warnings": ["Loop caller skipped sparse chromosome chrY"],
                },
                "differential_3d": {
                    "status": "success",
                    "ran": True,
                    "descriptive_only": False,
                    "pc1_validated": True,
                    "n_comparisons": 1,
                    "n_significant_compartment_bins": 8,
                    "n_ab_switches": 2,
                    "n_differential_insulation_bins": 5,
                    "warnings": ["Differential 3D omitted 1 zero-variance bin"],
                },
            },
        },
    }

    # Direct narrator call must NOT emit stage warning caveats (loops retired).
    raw_by_id = {b.id: b for b in HiCNarrator().collect("genome_arch_agent", agent_result)}
    assert not any("chrM" in c.text for c in raw_by_id["hic.qc"].caveats)
    assert not any("unplaced contig" in c.text for c in raw_by_id["hic.compartments"].caveats)

    # Through NarrativeRegistry.collect_blocks: all 5 stage warnings arrive with
    # severity="warning" and identical confidence/status.
    routed_by_id = {
        b.id: b for b in _full_registry().collect_blocks({"genome_arch_agent": agent_result})
    }
    expected = {
        "hic.qc": ("high", "ICE weight convergence note on chrM"),
        "hic.compartments": ("medium", "Compartment phasing skipped unplaced contig"),
        "hic.tads": ("medium", "TAD insulation boundary prominence relaxed on short arm"),
        "hic.loops": ("medium", "Loop caller skipped sparse chromosome chrY"),
        "hic.differential_3d": ("medium", "Differential 3D omitted 1 zero-variance bin"),
    }
    for bid, (exp_conf, exp_text) in expected.items():
        blk = routed_by_id[bid]
        assert blk.status == "success"
        assert blk.confidence == exp_conf
        assert any(c.severity == "warning" and c.text == exp_text for c in blk.caveats), (
            f"Missing parity caveat on {bid}: got {blk.caveats}"
        )


def test_chromatin_retired_and_newly_covered_stages_parity_through_registry():
    agent_result = {
        "status": "done",
        "findings": {
            "bulk_ATAC": {
                "status": "done",
                "findings": {
                    "qc": {
                        "status": "success",
                        "data_type": "bulk_ATAC",
                        "n_samples": 4,
                        "mito_fraction": 0.05,
                        "warnings": ["QC sample s4 had borderline NRF"],
                    },
                    "peaks": {
                        "status": "success",
                        "data_type": "bulk_ATAC",
                        "n_peaks": 2500,
                        "genome": "hg38",
                        "warnings": ["MACS3 broad-peak cutoff defaulted"],
                    },
                    "peak_counts": {
                        "status": "success",
                        "data_type": "bulk_ATAC",
                        "n_peaks": 2500,
                        "n_samples": 4,
                        "warnings": ["Single-end fragment extension assumed 200 bp"],
                    },
                    "differential_accessibility": {
                        "status": "success",
                        "data_type": "bulk_ATAC",
                        "ran": True,
                        "n_comparisons_success": 1,
                        "n_sig_total": 35,
                        "n_replicate_samples": 4,
                        "n_peaks_tested": 2500,
                        "warnings": ["Low-power contrast with n=2 replicates"],
                    },
                    "motifs": {
                        "status": "success",
                        "ran": True,
                        "motif_source": {"collection": "JASPAR2024", "n_motifs": 100},
                        "per_group": {"g1": {"n_enriched": 3}},
                        "warnings": ["Motif scan omitted 2 masked peaks"],
                    },
                    "figures": {"_errors": ["bulk_atac_pca: legend clipped"]},
                },
            }
        },
    }

    # Direct narrator call must NOT append raw stage warnings (loops retired).
    raw_by_id = {b.id: b for b in ChromatinNarrator().collect("chromatin_agent", agent_result)}
    assert not any("MACS3 broad-peak" in c.text for c in raw_by_id["chromatin.peak_calling.bulk_ATAC"].caveats)
    assert not any("Single-end fragment" in c.text for c in raw_by_id["chromatin.peak_count_matrix.bulk_ATAC"].caveats)
    assert not any("Low-power contrast" in c.text for c in raw_by_id["chromatin.differential_accessibility.bulk_atac"].caveats)

    # Through NarrativeRegistry.collect_blocks: retired loops AND newly covered stages
    # (qc, motifs, figures._errors) all surface with severity="warning".
    routed_by_id = {
        b.id: b for b in _full_registry().collect_blocks({"chromatin_agent": agent_result})
    }
    expected = {
        "chromatin.qc.bulk_ATAC": "QC sample s4 had borderline NRF",
        "chromatin.peak_calling.bulk_ATAC": "MACS3 broad-peak cutoff defaulted",
        "chromatin.peak_count_matrix.bulk_ATAC": "Single-end fragment extension assumed 200 bp",
        "chromatin.differential_accessibility.bulk_atac": "Low-power contrast with n=2 replicates",
        "chromatin.motifs": "Motif scan omitted 2 masked peaks",
    }
    for bid, exp_text in expected.items():
        blk = routed_by_id[bid]
        assert blk.status == "success"
        assert blk.confidence == "medium"
        assert any(c.severity == "warning" and c.text == exp_text for c in blk.caveats), (
            f"Missing parity caveat on {bid}: got {blk.caveats}"
        )
    assert any(
        c.severity == "warning" and "bulk_atac_pca: legend clipped" in c.text
        for c in routed_by_id["chromatin.qc.bulk_ATAC"].caveats
    )

    # Skipped bulk DA limitation block also preserves its stage warnings.
    skipped_res = {
        "status": "done",
        "findings": {
            "bulk_ATAC": {
                "status": "done",
                "findings": {
                    "differential_accessibility": {
                        "status": "skipped",
                        "data_type": "bulk_ATAC",
                        "ran": False,
                        "reason": "fewer than 2 replicates",
                        "warnings": ["Skipped DA replicate floor warning"],
                    }
                },
            }
        },
    }
    skipped_blocks = {
        b.id: b for b in _full_registry().collect_blocks({"chromatin_agent": skipped_res})
    }
    sk_blk = skipped_blocks["chromatin.differential_accessibility.bulk_atac_skipped"]
    assert sk_blk.status == "limitation"
    assert any(c.severity == "warning" and c.text == "Skipped DA replicate floor warning" for c in sk_blk.caveats)


@pytest.mark.parametrize(
    "agent_name,agent_result,target_block_id,warn_text,structured_sub",
    [
        (
            "bulk_rna_agent",
            {
                "status": "done",
                "findings": {
                    "sample_qc": {
                        "n_samples": 6,
                        "outliers": [],
                        "size_ratio": 1.4,
                        "warnings": ["Bulk sample QC free-text warning"],
                        "qc_error": "secondary MDS axis skipped",
                    }
                },
            },
            "bulk.qc",
            "Bulk sample QC free-text warning",
            "secondary MDS axis skipped",
        ),
        (
            "scrna_agent",
            {
                "status": "done",
                "findings": {
                    "qc": {
                        "status": "success",
                        "n_cells_before": 1000,
                        "n_cells_after": 900,
                        "warnings": ["scRNA QC free-text warning"],
                        "doublet_error": "Scrublet skipped on 1 lane",
                    }
                },
            },
            "scrna.qc",
            "scRNA QC free-text warning",
            "Scrublet skipped on 1 lane",
        ),
        (
            "chromatin_agent",
            {
                "status": "done",
                "findings": {
                    "qc": {
                        "status": "success",
                        "data_type": "scATAC",
                        "n_cells": 800,
                        "n_peaks": 12000,
                        "warnings": ["scATAC QC free-text warning"],
                        "tss_error": "GTF chrom mismatch",
                    }
                },
            },
            "chromatin.qc.scATAC",
            "scATAC QC free-text warning",
            "GTF chrom mismatch",
        ),
        (
            "genome_arch_agent",
            {
                "status": "done",
                "findings": {
                    "balancing": {
                        "status": "success",
                        "n_samples": 3,
                        "n_valid_pairs": 50000,
                        "cis_trans_ratio": 2.5,
                        "pass_qc": True,
                        "warnings": ["Hi-C balancing free-text warning"],
                        "scc_error": "HiCRep skipped chrY",
                    }
                },
            },
            "hic.qc",
            "Hi-C balancing free-text warning",
            "HiCRep skipped chrY",
        ),
        (
            "integration_agent",
            {
                "status": "done",
                "findings": {
                    "mofa": {
                        "status": "success",
                        "n_cells": 700,
                        "n_factors": 5,
                        "n_active_factors": 3,
                        "shared_factors": [1, 2],
                        "modality_specific_factors": [3],
                        "excluded_factors": [],
                        "variance_explained": {"rna": 0.25, "atac": 0.18},
                        "warnings": ["MOFA free-text warning"],
                        "projection_error": "posterior residual high",
                    }
                },
            },
            "multiome.mofa",
            "MOFA free-text warning",
            "posterior residual high",
        ),
    ],
)
def test_cross_modal_parity_across_all_five_narrators(
    agent_name, agent_result, target_block_id, warn_text, structured_sub
):
    blocks = _full_registry().collect_blocks({agent_name: agent_result})
    by_id = {b.id: b for b in blocks}
    blk = by_id[target_block_id]
    assert blk.status == "success"
    # Structured *_error caps high -> medium while preserving status="success".
    assert blk.confidence == "medium"
    assert any(c.severity == "warning" and c.text == warn_text for c in blk.caveats)
    assert any(c.severity == "warning" and structured_sub in c.text for c in blk.caveats)
    sigs = blk.metadata.get("degradation_signals") or []
    assert any(s["field"] == "warnings" and s["message"] == warn_text for s in sigs)
    assert any(s["field"].endswith("_error") and structured_sub in s["message"] for s in sigs)
