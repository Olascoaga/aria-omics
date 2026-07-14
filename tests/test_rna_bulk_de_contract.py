"""A7 characterization: pin the rna_bulk_de public contract through extraction.

`aria/scripts/rna_bulk_de.py` (2.4k lines) is the third A7 giant — a bulk DE
script that already has a partial helper subpackage (`aria.scripts.rna_bulk`:
gtf_io/ora/plots, P2-8). This slice continues that split, moving the remaining
counts/contrasts/transforms/QC/DESeq2 helper groups into the subpackage and
re-exporting them from `rna_bulk_de.py`, so callers and tests are unaffected.

This file locks two things the extraction must preserve:

1. Re-export surface — every public function the script exposes (incl. the ones
   other modules import: `_run_deseq2`, `_load_counts`, `bulk_rna_de`, …) stays
   importable from `aria.scripts.rna_bulk_de`.
2. Behavior of the pure helpers most at risk when moved (slug/formula/group
   inference).

Heavy DESeq2/QC paths (pydeseq2, sklearn) stay covered by the bulk-DE suites; the
module imports cleanly without them (they are lazy).
"""
from __future__ import annotations

import importlib

import pytest

_mod = importlib.import_module("aria.scripts.rna_bulk_de")


# ── 1. Re-export surface: every public name survives the split ────────────────

_PUBLIC_NAMES = [
    # entry
    "bulk_rna_de",
    # counts / metadata
    "_load_counts", "_enforce_metadata_correspondence",
    "_metadata_inference_allowed", "_load_or_infer_metadata",
    "_aggregate_technical_replicates", "_infer_groups", "_resolve_comparison",
    # contrasts
    "_slugify", "_format_top_genes", "_suggest_contrasts", "_auto_contrasts",
    "_contrast_overlap",
    # transforms
    "_run_vst", "_select_variable_genes", "_compute_tpm",
    # qc / outliers
    "_sample_qc", "_prune_outliers_for_design", "_run_outlier_sensitivity",
    # deseq2
    "_build_design_formula", "_resolve_covariates", "_shrink_coeff",
    "_serialize_fit_warnings", "_run_deseq2", "_mock_de_result",
    # already-extracted subpackage re-exports (P2-8) must remain
    "_run_pathway_enrichment", "_load_symbol_map", "_generate_plots",
]


@pytest.mark.parametrize("name", _PUBLIC_NAMES)
def test_public_name_is_importable_from_facade(name):
    assert hasattr(_mod, name), f"{name} no longer importable from rna_bulk_de"


def test_external_consumers_still_import_run_deseq2():
    # benchmarks/chromatin import this by name; keep the exact path working.
    from aria.scripts.rna_bulk_de import _run_deseq2, _load_counts  # noqa: F401


# ── 2. Pure helper behavior ───────────────────────────────────────────────────

def test_slugify_sanitises_and_defaults():
    assert _mod._slugify("Treated vs Control!") == "treated_vs_control"
    assert _mod._slugify("") == "contrast"
    assert _mod._slugify("!!!") == "contrast"


def test_build_design_formula_puts_factor_last_and_dedupes():
    assert _mod._build_design_formula("condition", []) == "~ condition"
    assert _mod._build_design_formula("condition", ["batch"]) == "~ batch + condition"
    # factor never repeated as a covariate; covariates deduped, order kept
    assert _mod._build_design_formula(
        "condition", ["batch", "condition", "batch", "sex"]
    ) == "~ batch + sex + condition"


def test_infer_groups_detects_condition_replicate_naming():
    groups = _mod._infer_groups(["ctrl_1", "treat_1", "ctrl_2", "treat_2"])
    assert groups == {
        "ctrl_1": "ctrl", "treat_1": "treat",
        "ctrl_2": "ctrl", "treat_2": "treat",
    }
    # a single group / unmatchable names -> None (no fabricated design)
    assert _mod._infer_groups(["onlyname"]) is None
