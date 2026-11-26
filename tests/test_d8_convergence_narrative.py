"""D8 convergence regression tests for narrative helpers and gene symbol validation (D8-003, D8-017, D8-019)."""

from __future__ import annotations

import pytest


def test_d8_003_narrator_stage_state_fail_closed():
    """D8-003: HiC and Integration narrators must share fail-closed stage_state."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.hic import HiCNarrator, _stage_state as hic_state
    from aria.agents.narrative.narrators.integration import (
        IntegrationNarrator,
        _stage_state as int_state,
    )

    stage = {"status": "cancelled", "n_cells_before": 500, "n_cells_after": 500}
    assert hic_state(stage) == "skipped"
    assert int_state(stage) == "skipped"
    hic_blocks = HiCNarrator().collect(
        "GenomeArchAgent",
        {"status": "success", "findings": {"qc": {"status": "cancelled", "n_samples": 2}}},
    )
    int_blocks = IntegrationNarrator().collect(
        "IntegrationAgent",
        {"status": "success", "findings": {"paired_qc": stage}},
    )
    assert hic_blocks and hic_blocks[0].status == "limitation" and hic_blocks[0].claim == ""
    assert int_blocks and int_blocks[0].status != "success" and int_blocks[0].claim == ""

    from aria.agents.narrative.protocols import stage_state

    assert hic_state is stage_state
    assert int_state is stage_state


def test_d8_017_narrator_as_int_inf_and_diagnostic_block_converge():
    """D8-017: narrator _as_int must not crash on inf and _diagnostic_block must unpack dict details with error=None on skipped."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.hic import (
        _as_int as hic_int,
        _diagnostic_block as hic_diag,
    )
    from aria.agents.narrative.narrators.integration import (
        _as_int as int_int,
        _diagnostic_block as int_diag,
    )

    assert hic_int(float("inf")) is None
    assert int_int(float("inf")) is None
    assert int_int(float("nan")) is None

    stage = {"status": "skipped", "details": {"reason": "low_barcode_overlap"}}
    b_hic = hic_diag("b1", "qc", "QC", stage, "skipped")
    b_int = int_diag("b1", "qc", "QC", stage, "skipped")
    assert b_hic.error is None and "low_barcode_overlap" in b_hic.caveats[0].text
    assert b_int.error is None and "low_barcode_overlap" in b_int.caveats[0].text
    assert "{" not in b_int.caveats[0].text

    from aria.agents.narrative.protocols import as_int

    assert hic_int is as_int
    assert int_int is as_int


def test_d8_019_gene_symbol_validators_converge():
    """D8-019: all 5 narrative gene-symbol validators must reject raw Ensembl/Entrez IDs, peak coords, and sentinels."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.hypothesis.adapters.cross_modal import _gene_symbol
    from aria.agents.narrative.hypothesis.adapters.hic import _split_gene_field
    from aria.agents.narrative.hypothesis.adapters.integration import _valid_gene_symbol
    from aria.agents.narrative.synthesis.cross_modal_synthesis import _clean_symbol
    from aria.agents.narrative.synthesis.pattern_detector import _looks_like_raw_id

    for bad in (
        "ENSG00000141510",
        "ENSMUSG00000059552",
        "7157",
        "nan",
        "None",
        "chr1:100-500",
        "chr1-100-500",
    ):
        assert _gene_symbol(bad) is None, bad
        assert _clean_symbol(bad) is None, bad
        assert _valid_gene_symbol(bad) is None, bad
        assert _split_gene_field(bad) == [], bad
        assert _looks_like_raw_id(bad) is True, bad

    for good in ("TP53", "SOX2", "ENSA", "HLA-DRA"):
        assert _gene_symbol(good) == good
        assert _clean_symbol(good) == good
        assert _valid_gene_symbol(good) == good
        assert _split_gene_field(good) == [good]
        assert _looks_like_raw_id(good) is False

    from aria.utils.text_matching import clean_gene_symbol

    assert _gene_symbol is clean_gene_symbol
    assert _valid_gene_symbol is clean_gene_symbol
