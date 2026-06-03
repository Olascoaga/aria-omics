"""P3-1: ChromatinNarrator skeleton (pre-4.6 polish).

The skeleton must (1) be registered as a built-in narrator from day one, (2)
surface only the structured QC metrics the scaffold chromatin_qc.py actually
measures, (3) mark not-computed metrics and the scaffold limitation honestly,
and (4) fabricate no chromatin science (no peaks/DA/motif claims). It only
fires on chromatin_agent results, so it is a no-op for the validated RNA paths.
"""

from aria.agents.narrative.narrators import ChromatinNarrator
from aria.agents.narrative.narrators.chromatin import ChromatinNarrator as Direct


def test_chromatin_narrator_is_a_registered_builtin():
    assert ChromatinNarrator is Direct
    n = ChromatinNarrator()
    assert n.name == "chromatin"


def test_chromatin_narrator_only_accepts_chromatin_agent():
    n = ChromatinNarrator()
    assert n.accepts("chromatin_agent", {"findings": {"qc": {}}}) is True
    assert n.accepts("scrna_agent", {"findings": {"qc": {}}}) is False
    assert n.accepts("chromatin_agent", {}) is False


def test_chromatin_narrator_surfaces_measured_qc_and_marks_not_computed():
    n = ChromatinNarrator()
    result = {"findings": {"qc": {
        "status": "success",
        "data_type": "scATAC",
        "n_cells": 500,
        "n_fragments": 1_000_000,
        "mito_fraction": 0.04,
        "frip": None,            # not computed until peaks called
        "tss_enrichment": None,  # not computed until a reference TSS QC runs
        "qc_complete": False,
        "pass_qc": None,
        "metrics_not_computed": ["frip", "tss_enrichment"],
    }}}
    blocks = n.collect("chromatin_agent", result, {})
    assert len(blocks) == 1
    qc = blocks[0]
    assert qc.modality == "chromatin"
    assert qc.metadata.get("validation_level") == "scaffold"
    # measured metrics appear as evidence; uncomputed metrics never fabricated
    labels = {e.label for e in qc.evidence}
    assert any("barcode" in lbl.lower() for lbl in labels)
    ev_values = {str(e.value) for e in qc.evidence}
    assert "None" not in ev_values
    # the scaffold limitation + not-computed list are disclosed
    caveat_text = " ".join(c.text.lower() for c in qc.caveats)
    assert "scaffold" in caveat_text
    assert "frip" in caveat_text and "tss_enrichment" in caveat_text


def test_chromatin_narrator_fabricates_no_science_without_qc():
    n = ChromatinNarrator()
    # No QC finding -> no blocks (never invents peaks/DA/motif narration).
    blocks = n.collect("chromatin_agent", {"findings": {}}, {})
    assert blocks == []
