"""Structured claim-to-evidence verification for narrative blocks.

This is the deterministic W-CLAIM layer: it checks ARIA-authored claim text and
rendered prose against the block's structured evidence card. It is deliberately
conservative and dependency-light; it verifies numeric support, named evidence
tokens, analysis-family support, and causal licensing without asking an LLM to
judge prose.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Any

from aria.agents.narrative.types import NarrativeBlock
from aria.agents.narrative.validators import find_causal_language


class EvidenceVerificationError(ValueError):
    """Raised when report prose is not supported by structured evidence."""


_EXPRESSION_TERMS = (
    "de gene", "de genes", "differentially expressed",
    "differential expression", "global-fdr de", "local-fdr de",
)
_PATHWAY_TERMS = (
    "pathway", "pathways", "gene-set", "gene set", "enriched term",
    "enriched terms", "gsea", "ora",
)
_ABUNDANCE_TERMS = (
    "abundance", "composition", "cell type shift", "cell-type shift",
)
_COMMUNICATION_TERMS = (
    "ligand-receptor", "ligand receptor", "interaction", "interactions",
    "cell communication",
)
_TRAJECTORY_TERMS = ("paga", "dpt", "trajectory", "pseudotime")

_ANALYSIS_FAMILY = {
    "differential_expression": "expression",
    "contrast": "expression",
    "pseudobulk_de": "expression",
    "pathway_enrichment": "pathway",
    "pseudobulk_pathways": "pathway",
    "gsea_preranked": "pathway",
    "differential_abundance": "abundance",
    "cell_communication": "communication",
    "trajectory": "trajectory",
    "qc": "descriptive",
    "sample_qc": "descriptive",
    "power": "descriptive",
}

_GENE_STOPWORDS = {
    "ARIA", "ATAC", "BAM", "CTRL", "DE", "DNA", "DPT", "FDR", "FRIP",
    "GEO", "GO", "GSEA", "HTML", "IF", "KEGG", "LLM", "NES", "ORA",
    "PAGA", "PBMC", "QC", "RNA", "SRA", "TAD", "UMAP", "VST",
}

_NUMBER_RE = re.compile(r"(?<![A-Za-z0-9_])-?\d+(?:\.\d+)?%?")
_GENE_LIKE_RE = re.compile(r"\b[A-Z][A-Z0-9_-]{2,}\b")


@dataclass
class VerificationIssue:
    location: str
    sentence: str
    reason: str

    def as_dict(self) -> dict[str, str]:
        return {
            "location": self.location,
            "sentence": self.sentence,
            "reason": self.reason,
        }


@dataclass
class EvidenceCard:
    evidence_card_id: str
    refs: list[dict[str, Any]] = field(default_factory=list)
    support_text: str = ""
    numbers: set[str] = field(default_factory=set)

    def as_dict(self) -> dict[str, Any]:
        return {
            "evidence_card_id": self.evidence_card_id,
            "refs": self.refs,
            "n_refs": len(self.refs),
        }


def verify_block_claim_support(
    block: NarrativeBlock,
    rendered_prose: str | None = None,
    *,
    strict: bool = True,
) -> dict[str, Any]:
    """Verify that a block's claim/prose is supported by structured evidence.

    The returned manifest is serializable and can be embedded in
    ``methodology.json``. When ``strict`` is true, any unsupported sentence raises
    ``EvidenceVerificationError`` so unsupported prose cannot be rendered.
    """
    card = build_evidence_card(block)
    issues: list[VerificationIssue] = []
    checked = 0

    seen_text: set[str] = set()
    for location, text in (
        ("claim", block.claim),
        ("rendered_prose", rendered_prose or ""),
    ):
        normalized_text = str(text or "").strip()
        if not normalized_text or normalized_text in seen_text:
            continue
        seen_text.add(normalized_text)
        for sentence in _sentences(normalized_text):
            checked += 1
            issues.extend(_issues_for_sentence(block, card, location, sentence))

    manifest = {
        "status": "supported" if not issues else "unsupported",
        "checked_sentences": checked,
        "unsupported": [issue.as_dict() for issue in issues],
        "evidence_card": card.as_dict(),
        "policy": (
            "Every ARIA-authored claim sentence must be backed by structured "
            "evidence, metrics, tables, figures, or an explicitly licensed "
            "analysis family."
        ),
    }
    if strict and issues:
        first = issues[0]
        raise EvidenceVerificationError(
            f"{block.id}: unsupported {first.location} sentence: "
            f"{first.reason}: {first.sentence}"
        )
    return manifest


def build_evidence_card(block: NarrativeBlock) -> EvidenceCard:
    """Return the structured evidence card used to verify one block."""
    refs: list[dict[str, Any]] = []
    support_parts = [
        block.id,
        block.title,
        block.analysis,
        block.modality,
    ]
    numbers: set[str] = set()

    for idx, ev in enumerate(block.evidence or []):
        ref = {
            "ref_id": f"{block.id}:evidence:{idx}",
            "kind": "evidence",
            "label": ev.label,
            "value": ev.value,
            "source": ev.source,
            "path": ev.path,
        }
        refs.append(ref)
        support_parts.extend(
            str(x) for x in (ev.label, ev.value, ev.source, ev.path) if x is not None
        )
        numbers.update(_numbers_in(" ".join(str(x) for x in (ev.label, ev.value))))

    for key, value in _flatten_mapping(block.metrics or {}):
        refs.append({
            "ref_id": f"{block.id}:metric:{key}",
            "kind": "metric",
            "label": key,
            "value": value,
        })
        support_parts.extend((str(key), str(value)))
        numbers.update(_numbers_in(f"{key} {value}"))

    for idx, table in enumerate(block.tables or []):
        if not isinstance(table, dict):
            continue
        refs.append({
            "ref_id": f"{block.id}:table:{idx}",
            "kind": "table",
            "label": table.get("label") or table.get("id"),
            "path": table.get("path"),
        })
        support_parts.extend(str(v) for v in table.values() if v is not None)

    for idx, figure in enumerate(block.figures or []):
        if not isinstance(figure, dict):
            continue
        refs.append({
            "ref_id": f"{block.id}:figure:{idx}",
            "kind": "figure",
            "label": figure.get("caption") or figure.get("id"),
            "path": figure.get("path"),
        })
        support_parts.extend(str(v) for v in figure.values() if v is not None)

    return EvidenceCard(
        evidence_card_id=f"{block.id}#evidence",
        refs=refs,
        support_text=" ".join(support_parts).lower(),
        numbers=numbers,
    )


def _issues_for_sentence(
    block: NarrativeBlock,
    card: EvidenceCard,
    location: str,
    sentence: str,
) -> list[VerificationIssue]:
    issues: list[VerificationIssue] = []
    low = sentence.lower()
    family = _ANALYSIS_FAMILY.get(block.analysis)

    for number in _claim_numbers(sentence):
        if number not in card.numbers:
            issues.append(VerificationIssue(
                location,
                sentence,
                f"numeric value {number!r} is absent from the evidence card",
            ))

    for token in _claim_entities(sentence):
        if token.lower() not in card.support_text:
            issues.append(VerificationIssue(
                location,
                sentence,
                f"named entity {token!r} is absent from the evidence card",
            ))

    required = _required_family(low)
    if required and required != family:
        # Pathway prose often explains that terms summarize a DE list. That is
        # still pathway evidence, not a second unsupported DE claim.
        pathway_explains_de = (
            required == "expression"
            and family == "pathway"
            and any(term in low for term in _PATHWAY_TERMS)
        )
        if not pathway_explains_de:
            issues.append(VerificationIssue(
                location,
                sentence,
                f"sentence asserts {required!r} evidence but block analysis is "
                f"{block.analysis!r}",
            ))

    claim_meta = block.metadata.get("claim") or {}
    licensed = claim_meta.get("licensed_language")
    if (
        licensed
        and licensed != "causal"
        and not block.metadata.get("causal_evidence")
    ):
        hit = find_causal_language(sentence)
        if hit and "not " not in low[max(0, low.find(hit.lower()) - 8):low.find(hit.lower())]:
            issues.append(VerificationIssue(
                location,
                sentence,
                f"causal language {hit!r} exceeds the licensed "
                f"{licensed!r} claim tier",
            ))

    return issues


def _required_family(sentence_lower: str) -> str | None:
    if "composition covariate" in sentence_lower:
        return None
    checks = (
        ("expression", _EXPRESSION_TERMS),
        ("pathway", _PATHWAY_TERMS),
        ("abundance", _ABUNDANCE_TERMS),
        ("communication", _COMMUNICATION_TERMS),
        ("trajectory", _TRAJECTORY_TERMS),
    )
    for family, terms in checks:
        if any(term in sentence_lower for term in terms):
            return family
    return None


def _sentences(text: str) -> list[str]:
    return [
        sentence.strip()
        for sentence in re.split(r"(?<=[.!?])\s+", str(text or "").strip())
        if sentence.strip()
    ]


def _claim_numbers(sentence: str) -> set[str]:
    numbers = set()
    for raw in _numbers_in(sentence):
        # Very small values usually represent thresholds (FDR < 0.05), signs,
        # or ordinal prose. Large integers/counts are the support-critical facts.
        try:
            value = float(raw.rstrip("%"))
        except ValueError:
            continue
        if abs(value) >= 2:
            numbers.add(raw)
    return numbers


def _numbers_in(text: str) -> set[str]:
    return {_normalize_number(m.group(0)) for m in _NUMBER_RE.finditer(str(text or ""))}


def _normalize_number(raw: str) -> str:
    raw = str(raw).strip()
    pct = raw.endswith("%")
    if pct:
        raw = raw[:-1]
    try:
        value = float(raw)
    except ValueError:
        return raw
    if value.is_integer():
        out = str(int(value))
    else:
        out = f"{value:.6g}"
    return out + ("%" if pct else "")


def _claim_entities(sentence: str) -> set[str]:
    entities = set()
    for match in _GENE_LIKE_RE.finditer(sentence):
        token = match.group(0)
        if token in _GENE_STOPWORDS:
            continue
        if token.startswith("X_"):
            continue
        entities.add(token)
    return entities


def _flatten_mapping(data: dict[str, Any], prefix: str = ""):
    for key, value in data.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            yield from _flatten_mapping(value, name)
        elif isinstance(value, list):
            for idx, item in enumerate(value[:10]):
                item_name = f"{name}.{idx}"
                if isinstance(item, dict):
                    yield from _flatten_mapping(item, item_name)
                else:
                    yield item_name, item
        else:
            yield name, value
