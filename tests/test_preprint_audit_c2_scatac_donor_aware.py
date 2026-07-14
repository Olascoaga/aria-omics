"""Plumbing guards for the c2_scatac_donor_aware preprint-freeze lane.

These validate the lane registration, the frozen donor cohort split and the
artifact contract cheaply, without the real multi-donor DA concordance run. The
real receipt is regenerated once, against the final clean source snapshot,
during the freeze regeneration step (see memory/NEXT_SESSION.md).

The lane wires the already real-validated
``scripts/run_scatac_multisample_da_concordance.py`` runner (ARIA pseudobulk
CONDITION-DA vs R DESeq2/edgeR/limma on the real GSE278576 per-donor scATAC
matrices) to the freeze inventory; it is fully deterministic and LLM-free.
"""

from __future__ import annotations

import shlex
from pathlib import Path

from aria.benchmarks.preprint_freeze import LANES, _resource_defaults


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = REPO_ROOT / "scripts" / "run_scatac_multisample_da_concordance.py"

# The frozen young/old donor split, pinned from the real-validated prior run
# (docs/benchmark_results/scatac_concordance/p3_multisample_aging_da_concordance.json).
FROZEN_YOUNG = ("hc937", "hc5579", "hc13344", "hc935", "hc5614")
FROZEN_OLD = ("hc26", "hc98", "hc35", "hc40", "hc73")


def _lane():
    return next(item for item in LANES if item["lane_id"] == "c2_scatac_donor_aware")


def _command_opts():
    """Parse the lane command into a {flag: value} map for its long options."""
    tokens = shlex.split(_lane()["command"])
    opts: dict[str, str] = {}
    for i, tok in enumerate(tokens):
        if tok.startswith("--") and i + 1 < len(tokens):
            opts[tok] = tokens[i + 1]
    return opts


def test_lane_is_registered_and_executable():
    lane = _lane()
    assert lane["claims"] == ["claim_2"]
    assert lane["implementation"] == "available"
    assert lane["evidence_kind"] == "external_concordance"
    assert lane["required_for_freeze"] is True
    assert "run_scatac_multisample_da_concordance.py" in lane["command"]


def test_runner_and_referenced_scripts_exist():
    assert RUNNER.is_file()
    # The runner dispatches into these two helpers; a missing one would only
    # surface at run time, so pin them here.
    assert (REPO_ROOT / "scripts" / "aria_pseudobulk_da_from_tsv.py").is_file()
    assert (REPO_ROOT / "aria" / "scripts"
            / "benchmark_a1_external_comparators.R").is_file()


def test_lane_declares_single_manifest_artifact():
    lane = _lane()
    assert tuple(lane["expected_artifacts"]) == (
        "claim_2/scatac_donor/scatac_donor_da_concordance.json",
    )


def test_lane_binds_both_scientific_environments():
    lane = _lane()
    assert lane["resources"] == [
        "env:aria-rna-env", "env:aria-bench-env", "data:gse278576"
    ]
    opts = _command_opts()
    # The nested conda-run dispatch targets must match the declared envs.
    assert opts["--rna-env"] == "aria-rna-env"
    assert opts["--bench-env"] == "aria-bench-env"


def test_command_pins_the_frozen_donor_split():
    opts = _command_opts()
    young = tuple(opts["--young"].split(","))
    old = tuple(opts["--old"].split(","))
    assert young == FROZEN_YOUNG
    assert old == FROZEN_OLD
    # Balanced 5v5 design; no donor may appear in both groups.
    assert len(young) == len(old) == 5
    assert not (set(young) & set(old))


def test_command_output_matches_declared_artifact():
    opts = _command_opts()
    out_dir = opts["--output-dir"].rstrip("/")
    manifest = opts["--manifest-name"]
    relative = f"{out_dir}/{manifest}".removeprefix(
        "docs/benchmark_results/preprint_v1/"
    )
    assert relative == _lane()["expected_artifacts"][0]


def test_data_probe_points_at_the_real_scatac_donor_matrix():
    defaults = _resource_defaults(
        REPO_ROOT / "docs/benchmark_results/preprint_v1"
    )
    path_ref, path = defaults["data:gse278576"]
    # The probe must reflect scATAC donor availability, not the RNA-only h5ad.
    assert "scatac" in path_ref
    assert path is not None and path.name.endswith("_paired.h5mu")
