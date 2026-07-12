"""Preprint-readiness audit B4: no post-confirmation column-name inference.

P0-6/F5 already ban inferring the experimental design from file/column names in the
production path. The residual back door was inside ``BulkRNAAgent._apply_design``:
when a CONFIRMED design's sample count did not match the count-matrix columns, the
mapping silently fell back to ``_infer_col_groups`` and rebuilt the design from
``ctrl_*`` / ``treat_*`` column-name patterns — bypassing the gate and the explicit
``ARIA_ALLOW_FILENAME_FALLBACK`` opt-in.

After B4 an unmappable confirmed design must fail closed (raise), which routes the
run through the existing P0-6 gate: block, or infer only under the explicit opt-in
via the sample-name path — never silently from column names.
"""
from __future__ import annotations

import pytest


def _agent():
    pytest.importorskip("litellm")  # importing the agent pulls aria.llm.provider
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    return object.__new__(BulkRNAAgent)


def _write_counts(path, header_cols):
    path.write_text(
        "gene\t" + "\t".join(header_cols) + "\n"
        + "g1\t" + "\t".join("1" for _ in header_cols) + "\n"
    )
    return str(path)


def test_infer_col_groups_backdoor_removed():
    pytest.importorskip("litellm")
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    # The column-name inference helper (and its only call site) is gone.
    assert not hasattr(BulkRNAAgent, "_infer_col_groups")


def test_apply_design_blocks_column_name_inference_on_mismatch(tmp_path):
    agent = _agent()
    # Confirmed design carries two GSM samples; the matrix has four ctrl_*/treat_*
    # columns (incomplete metadata / sample-count mismatch). Previously this was
    # silently rescued by inferring {ctrl, treat} from the column names.
    design = {
        "groups": {"KO": ["GSM_KO"], "WT": ["GSM_WT"]},
        "main_factor": "genotype",
    }
    counts = _write_counts(
        tmp_path / "counts.tsv", ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]
    )
    with pytest.raises(ValueError):
        agent._apply_design(design, [counts], "exp-b4")


def test_apply_design_blocks_even_when_counts_match_but_names_do_not(tmp_path):
    """A same-count mismatch where the confirmed stems/aliases match NOTHING in the
    columns must also block, not be rescued by ctrl_*/treat_* name inference."""
    agent = _agent()
    design = {
        "groups": {"KO": ["GSM_KO_a", "GSM_KO_b"], "WT": ["GSM_WT_a", "GSM_WT_b"]},
        "main_factor": "genotype",
    }
    # Four columns, four confirmed samples, but the NAMES share no signal with the
    # confirmed stems. The positional fallback uses the confirmed design (order),
    # so this must NOT silently switch to column-name groups; if the confirmed
    # design maps positionally it keeps KO/WT, and it must never yield ctrl/treat.
    counts = _write_counts(
        tmp_path / "counts.tsv", ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]
    )
    sample_names, group_labels, factor, contrasts = agent._apply_design(
        design, [counts], "exp-b4"
    )
    assert set(group_labels.values()) == {"KO", "WT"}
    assert "ctrl" not in set(group_labels.values())
    assert "treat" not in set(group_labels.values())


def test_apply_design_still_maps_confirmed_design_via_aliases(tmp_path):
    """Positive control: a confirmed design whose aliases match the columns still
    maps cleanly — B4 removes the name-inference rescue, not confirmed mapping."""
    agent = _agent()
    design = {
        "groups": {"KO": ["GSM1"], "WT": ["GSM2"]},
        "main_factor": "genotype",
        "sample_aliases": {"GSM1": ["GSM1", "mut_a"], "GSM2": ["GSM2", "wt_a"]},
        "contrasts": [["KO", "WT"]],
    }
    counts = _write_counts(tmp_path / "counts.tsv", ["mut_a", "wt_a"])
    sample_names, group_labels, factor, contrasts = agent._apply_design(
        design, [counts], "exp-b4"
    )
    assert set(group_labels.values()) == {"KO", "WT"}
    assert contrasts and contrasts[0]["numerator"] == "KO"
