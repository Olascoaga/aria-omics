"""Causal-language guard hardening (audit 2026-05-28, F-SCI-CAUSAL).

Covers the broadened lexicon and the new render-level prose scan, which the
claim/evidence-only validator did not previously reach.
"""

from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.agents.narrative.validators import find_causal_language, validate_block
from aria.agents.narrative.render_blocks import render_blocks


def test_find_causal_language_detects_broadened_terms():
    for term in ("induces", "triggers", "master regulator", "causes",
                 "leads to", "controls", "upstream of", "mediates"):
        sentence = f"This module {term} the observed program."
        assert find_causal_language(sentence) is not None, term


def test_find_causal_language_passes_associative_text():
    assert find_causal_language(
        "Expression of this set is associated with the stimulated state and "
        "correlates with the response signature."
    ) is None


def _causal_block():
    return NarrativeBlock(
        id="x.causal",
        modality="scRNA",
        analysis="freeform",   # not a recognized analysis -> prose echoes claim
        block_type="result",
        title="Synthetic block",
        status="success",
        confidence="high",
        claim="Factor M induces the differentiation program in cluster 3.",
        evidence=[EvidenceItem(label="n", value=10, source="test")],
    )


def test_validate_block_downgrades_and_caveats_on_causal_claim():
    block = validate_block(_causal_block())
    # high -> medium downgrade and an associative caveat must be attached.
    assert block.confidence == "medium"
    assert any("associative" in c.text.lower() for c in block.caveats)


def test_render_warns_when_composed_prose_is_causal():
    html = render_blocks([_causal_block()])
    assert "correlation, not causation" in html


def test_render_suppresses_warning_when_causal_evidence_declared():
    block = _causal_block()
    block.metadata["causal_evidence"] = True
    html = render_blocks([block])
    assert "correlation, not causation" not in html
