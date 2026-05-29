"""Causal-language guard hardening (audit 2026-05-28, F-SCI-CAUSAL).

Covers the broadened lexicon and the new render-level prose scan, which the
claim/evidence-only validator did not previously reach.
"""

from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock
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


# --- B5 (audit 2026-05-29): a database term name that contains a causal verb
# is data ARIA reports, not a claim ARIA asserts. It must not downgrade or warn.


def _ora_block_with_regulatory_term_name():
    """Honest ORA block whose enriched term name embeds a causal verb."""
    return NarrativeBlock(
        id="bulk.pathway.x",
        modality="bulk RNA-seq",
        analysis="pathway_enrichment",
        block_type="result",
        title="Pathway enrichment cluster 3",
        status="success",
        confidence="medium",
        claim="Pathway enrichment found 2 term(s) for cluster 3.",
        evidence=[
            EvidenceItem(label="enriched terms", value=2, source="bulk_pathways"),
            EvidenceItem(
                label="Reactome term",
                value="TP53 Regulates Transcription of Cell Cycle Genes",
                source="bulk_pathways",
            ),
            EvidenceItem(label="Reactome term FDR", value=0.01,
                         source="bulk_pathways"),
        ],
        caveats=[Caveat("ORA summarizes DE gene sets, not causal proof.",
                        "info")],
        metrics={"n_terms": 2},
    )


def test_term_name_with_causal_verb_does_not_trip_validator():
    block = _ora_block_with_regulatory_term_name()
    validate_block(block)
    # confidence is unchanged and no causal caveat was attached
    assert block.confidence == "medium"
    assert not any("Causal language pattern" in c.text for c in block.caveats)


def test_term_name_with_causal_verb_does_not_warn_at_render():
    # the composed prose interpolates the term name "...Regulates..."
    html = render_blocks([_ora_block_with_regulatory_term_name()])
    assert "correlation, not causation" not in html


def test_causal_claim_still_caught_despite_term_name_evidence():
    """A real causal CLAIM is still downgraded even when evidence has term
    names — only the named entities are exempt, not ARIA's assertion."""
    block = _ora_block_with_regulatory_term_name()
    block.claim = "This pathway drives the differentiation program."
    validate_block(block)
    assert block.confidence == "low"  # medium -> low
    assert any("Causal language pattern" in c.text and "associative" in c.text
               for c in block.caveats)


def test_find_causal_language_exclude_redacts_named_entities():
    text = ("Leading terms include TP53 Regulates Transcription of Cell Cycle "
            "Genes and the immune response signature.")
    assert find_causal_language(text) is not None  # without exclude, fires
    assert find_causal_language(
        text, exclude=["TP53 Regulates Transcription of Cell Cycle Genes"]
    ) is None
