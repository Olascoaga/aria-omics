"""Claim Compiler (audit item X14).

Every biological claim ARIA reports is classified into an evidence tier from the
STRUCTURED evidence that actually exists for it — never from the LLM's prose.
The tier caps the language the report is allowed to use, so an associative
result can never be dressed up as a causal one.

Tiers (ascending strength):

  descriptive          counts / QC / composition; no between-condition claim.
  associative          a correlative contrast (DE, abundance, communication);
                       observational — no mechanism, no intervention.
  weak_mechanistic     associative + one independent mechanistic line
                       (regulatory motif / accessibility / regulon).
  strong_mechanistic   associative + multiple converging mechanistic lines,
                       still observational.
  causal_experimental  the design is interventional (perturbation / KO / KD /
                       controlled treatment) and the contrast tests it.

Licensed language per tier:

  descriptive / associative / weak_mechanistic / strong_mechanistic  -> at most
      ASSOCIATIVE language (a mechanism may be proposed as a hypothesis, never
      asserted as causation) unless the design is interventional.
  causal_experimental -> CAUSAL language is licensed.

This is deterministic and dependency-light. It builds on the existing
`NarrativeBlock` evidence cards and the causal-language guard in `validators`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from aria.agents.narrative.types import NarrativeBlock
from aria.agents.narrative.evidence_verifier import verify_block_claim_support
from aria.agents.narrative.validators import find_causal_language


TIER_ORDER = (
    "descriptive",
    "associative",
    "weak_mechanistic",
    "strong_mechanistic",
    "causal_experimental",
)

# Map a block's analysis type to the evidence category it contributes.
_ANALYSIS_EVIDENCE = {
    "qc": "descriptive",
    "sample_qc": "descriptive",
    "power": "descriptive",
    "executive_summary": "descriptive",
    # B-DD1 (scRNA-lane audit 2026-06-11): per-cluster marker discovery is a
    # cluster-vs-rest test run on the SAME cells/embedding that defined the
    # clusters (Leiden on the same PCA/HVG). Its p-values are anti-conservative
    # by selection-then-test (double-dipping), so they rank candidate markers
    # but are NOT valid inferential significance. Cap it to descriptive — the
    # between-condition inferential weight is the pseudobulk DE, never markers.
    "marker_discovery": "descriptive",
    "differential_expression": "expression_change",
    "contrast": "expression_change",
    "pseudobulk_de": "expression_change",
    "differential_abundance": "abundance_change",
    "pathway_enrichment": "functional_annotation",
    "pseudobulk_pathways": "functional_annotation",
    "gsea_preranked": "functional_annotation",
    "cell_communication": "communication",
    "trajectory": "trajectory",
    # forward-looking (chromatin / multiome, v4.6+):
    "motif_enrichment": "regulatory",
    "differential_accessibility": "regulatory",
    "peak2gene": "regulatory",
    "regulon": "regulatory",
}

# Evidence categories that count as observational/associative support.
_ASSOCIATIVE = {
    "expression_change", "abundance_change", "functional_annotation",
    "communication", "trajectory",
}
# Evidence categories that count as independent mechanistic lines.
_MECHANISTIC = {"regulatory"}
_STATS_GATE_ANALYSES = {"pseudobulk_de", "differential_expression"}
_POWER_MIN = 0.5


@dataclass
class ClaimClassification:
    tier: str
    licensed_language: str          # "descriptive" | "associative" | "causal"
    evidence_categories: list[str]
    rationale: str
    limitations: list[str] = field(default_factory=list)
    language_violation: str | None = None   # causal phrase found above the tier

    def as_dict(self) -> dict[str, Any]:
        return {
            "tier": self.tier,
            "licensed_language": self.licensed_language,
            "evidence_categories": self.evidence_categories,
            "rationale": self.rationale,
            "limitations": self.limitations,
            "language_violation": self.language_violation,
        }


def design_is_interventional(exp_ctx: dict | None) -> bool:
    """Conservative: observational unless the design is EXPLICITLY interventional.

    ARIA cannot tell a controlled perturbation from an observational label by
    name alone, so causal licensing requires an explicit marker. This keeps the
    default honest (associative) and lets interventional designs opt in.
    """
    design = (exp_ctx or {}).get("design", {}) or {}
    if design.get("interventional") is True:
        return True
    di = (exp_ctx or {}).get("design_intelligence", {}) or {}
    if di.get("interventional") is True:
        return True
    return False


def _block_evidence_category(block: NarrativeBlock) -> str | None:
    if block.block_type == "qc":
        return "descriptive"
    cat = _ANALYSIS_EVIDENCE.get(block.analysis)
    if cat:
        return cat
    # Unknown successful result with evidence is at least associative.
    if block.status == "success" and block.evidence:
        return "expression_change"
    return None


def _block_subject(block: NarrativeBlock) -> str:
    """Group key for converging evidence (e.g. a cell type x comparison)."""
    parts = str(block.id).split(".")
    # scrna.pseudobulk.<group>.<comparison> -> group.comparison
    if len(parts) >= 4:
        return f"{parts[2]}.{parts[3]}"
    if len(parts) >= 2:
        return parts[-1]
    return str(block.id)


def _numeric_metric(block: NarrativeBlock, *names: str) -> float | None:
    for name in names:
        value = (block.metrics or {}).get(name)
        if value is None:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def _flag_metric(block: NarrativeBlock, name: str) -> bool:
    if name in (block.metrics or {}):
        return bool(block.metrics.get(name))
    text = " ".join(c.text for c in (block.caveats or [])).lower()
    if name == "low_power_warning":
        return "low replicate" in text or "low power" in text
    if name == "lognorm_recovered":
        return "log-normalized" in text or "reverse-engineered" in text
    return False


def _quantitative_stats_gate(block: NarrativeBlock) -> tuple[bool, list[str]]:
    """P-CLAIM2: DE claims need numerical support, not only analysis type."""
    if block.analysis not in _STATS_GATE_ANALYSES:
        return True, []
    failures: list[str] = []
    n_sig = _numeric_metric(
        block, "n_significant", "n_significant_global", "n_significant_local"
    )
    if n_sig is not None and n_sig <= 0:
        failures.append("no significant features passed the selected FDR rule")
    power = _numeric_metric(
        block, "power_estimate_at_effective_alpha",
        "power_estimate_at_lfc_min",
    )
    if power is not None and power < _POWER_MIN:
        failures.append(
            f"estimated power ({power:.2f}) is below {_POWER_MIN:.2f}"
        )
    if _flag_metric(block, "low_power_warning"):
        failures.append("low replicate/power warning is active")
    if _flag_metric(block, "lognorm_recovered"):
        failures.append("DE used recovered log-normalized counts")
    return not failures, failures


def classify_claim(block: NarrativeBlock,
                   *,
                   interventional: bool = False,
                   converging_categories: set[str] | None = None) -> ClaimClassification:
    """Classify one block's claim into an evidence tier.

    ``converging_categories`` are the evidence categories available from OTHER
    blocks about the same subject, so independent mechanistic lines can raise
    the tier without the LLM's help.
    """
    own = _block_evidence_category(block)
    categories: set[str] = set()
    if own:
        categories.add(own)
    categories |= (converging_categories or set())

    limitations: list[str] = []

    # Pure descriptive (QC / counts / power): no contrast claim.
    contentful = categories - {"descriptive"}
    if not contentful:
        tier = "descriptive"
        licensed = "descriptive"
        rationale = "Descriptive summary (QC/counts); no between-condition claim."
        return _finalize(block, tier, licensed, sorted(categories), rationale,
                         limitations, interventional)

    has_assoc = bool(categories & _ASSOCIATIVE)
    n_mech = len(categories & _MECHANISTIC)

    if interventional and has_assoc:
        tier = "causal_experimental"
        licensed = "causal"
        rationale = ("Interventional design: the contrast tests a perturbation, "
                     "so the effect of the intervention is causally licensed.")
    elif n_mech >= 2 and has_assoc:
        tier = "strong_mechanistic"
        licensed = "associative"
        rationale = ("Associative contrast plus multiple converging mechanistic "
                     "lines (e.g. regulatory evidence); still observational.")
        limitations.append(
            "Observational design: converging mechanism is a strong hypothesis, "
            "not experimental proof of causation.")
    elif n_mech >= 1 and has_assoc:
        tier = "weak_mechanistic"
        licensed = "associative"
        rationale = ("Associative contrast plus one independent mechanistic line; "
                     "mechanism is a hypothesis, not established.")
        limitations.append(
            "Observational design with limited mechanistic support; "
            "causation is not established.")
    else:
        tier = "associative"
        licensed = "associative"
        rationale = ("Correlative contrast (expression/abundance/communication) "
                     "without independent mechanistic or perturbation evidence.")
        limitations.append(
            "Observational association only; no mechanistic or interventional "
            "evidence supports a causal interpretation.")

    stats_ok, stats_failures = _quantitative_stats_gate(block)
    if not stats_ok and tier != "causal_experimental":
        tier = "descriptive"
        licensed = "descriptive"
        rationale = (
            "Quantitative stats-evidence gate did not license an associative "
            "claim for this DE result."
        )
        limitations.append(
            "Stats-evidence gate: " + "; ".join(stats_failures) + "."
        )

    return _finalize(block, tier, licensed, sorted(categories), rationale,
                     limitations, interventional)


def _finalize(block, tier, licensed, categories, rationale, limitations,
              interventional) -> ClaimClassification:
    # Inherit block-level limitations already known to the kernel.
    for cav in block.caveats or []:
        if cav.text and cav.text not in limitations:
            limitations.append(cav.text)
    if "trajectory" in categories:
        limitations.append(
            "Trajectory ordering (PAGA/DPT) is exploratory without velocity or "
            "time-course data.")

    # Language check: does the claim/prose exceed the licensed tier?
    violation = None
    if licensed != "causal" and not block.metadata.get("causal_evidence"):
        text = " ".join(filter(None, [
            block.claim,
            *[str(c.text) for c in (block.caveats or [])],
        ]))
        hit = find_causal_language(text)
        if hit:
            violation = hit
            limitations.append(
                f"Claim language ('{hit}') implies causation beyond the "
                f"'{tier}' evidence tier; treat as associative.")

    # De-duplicate while preserving order.
    seen = set()
    deduped = [x for x in limitations if not (x in seen or seen.add(x))]
    return ClaimClassification(
        tier=tier,
        licensed_language=licensed,
        evidence_categories=categories,
        rationale=rationale,
        limitations=deduped,
        language_violation=violation,
    )


def annotate_claim_tiers(blocks: list[NarrativeBlock],
                         exp_ctx: dict | None = None) -> list[NarrativeBlock]:
    """Classify every block in place, storing the result under
    ``block.metadata['claim']``. Returns the same list for chaining."""
    interventional = design_is_interventional(exp_ctx)

    # Converging evidence: union of categories per subject.
    by_subject: dict[str, set[str]] = {}
    for b in blocks:
        cat = _block_evidence_category(b)
        if cat:
            by_subject.setdefault(_block_subject(b), set()).add(cat)

    for b in blocks:
        subject = _block_subject(b)
        own = _block_evidence_category(b)
        converging = set(by_subject.get(subject, set()))
        if own:
            converging.discard(own)  # the rest are the converging lines
        classification = classify_claim(
            b, interventional=interventional, converging_categories=converging,
        )
        b.metadata["claim"] = classification.as_dict()
    return blocks


def compile_claim_manifest(block: NarrativeBlock) -> dict[str, Any]:
    """Build a per-claim evidence manifest for methodology.json."""
    claim_meta = block.metadata.get("claim") or {}
    verification = block.metadata.get("claim_verification")
    if not verification:
        verification = verify_block_claim_support(block, strict=False)
        block.metadata["claim_verification"] = verification
    evidence = []
    for ev in block.evidence or []:
        evidence.append({
            "label": ev.label,
            "value": ev.value,
            "source": ev.source,
            "path": ev.path,
        })
    tables = [t.get("path") for t in (block.tables or []) if t.get("path")]
    figures = [f.get("path") for f in (block.figures or []) if f.get("path")]
    return {
        "claim_id": block.id,
        "text": block.claim,
        "modality": block.modality,
        "analysis": block.analysis,
        "status": block.status,
        "confidence": block.confidence,
        "tier": claim_meta.get("tier"),
        "licensed_language": claim_meta.get("licensed_language"),
        "evidence_categories": claim_meta.get("evidence_categories", []),
        "rationale": claim_meta.get("rationale"),
        "limitations": claim_meta.get("limitations", []),
        "language_violation": claim_meta.get("language_violation"),
        "verification": verification,
        "evidence_card_id": (
            verification.get("evidence_card", {}).get("evidence_card_id")
            if isinstance(verification, dict) else None
        ),
        "evidence": evidence,
        "tables": tables,
        "figures": figures,
    }


def compile_claims(blocks: list[NarrativeBlock],
                   exp_ctx: dict | None = None) -> list[dict[str, Any]]:
    """Annotate tiers (if not already) and return the claim manifests for all
    contentful (non-descriptive) blocks plus successful descriptive ones."""
    annotate_claim_tiers(blocks, exp_ctx)
    return [compile_claim_manifest(b) for b in blocks]


@dataclass
class PublicClaimCompilation:
    """The only block/claim set eligible for public report rendering."""

    blocks: list[NarrativeBlock] = field(default_factory=list)
    claims: list[dict[str, Any]] = field(default_factory=list)
    withheld: list[dict[str, str]] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        return {
            "compiler": "compile_public_claims",
            "n_published": len(self.claims),
            "n_withheld": len(self.withheld),
            "withheld": list(self.withheld),
        }


def compile_public_claims(
    blocks: list[NarrativeBlock],
    exp_ctx: dict | None = None,
) -> PublicClaimCompilation:
    """Compile the exclusive public claim set, withholding unsupported blocks.

    MessageBus payloads and legacy prose are deliberately not accepted by this
    API: callers must provide typed ``NarrativeBlock`` instances with structured
    evidence. Each block is validated, tiered and checked against the exact prose
    the renderer will show. A failed block contributes only an opaque withheld
    record; its unverified text never reaches HTML or the claim manifest.
    """
    from aria.agents.narrative.compose_prose import compose_block_prose
    from aria.agents.narrative.validators import validate_block

    candidates = list(blocks or [])
    annotate_claim_tiers(candidates, exp_ctx or {})
    compiled = PublicClaimCompilation()
    for block in candidates:
        try:
            validate_block(block)
            prose = compose_block_prose(block) or block.claim or ""
            verification = verify_block_claim_support(
                block, prose, strict=True
            )
            block.metadata["claim_verification"] = verification
            manifest = compile_claim_manifest(block)
            if manifest.get("verification", {}).get("status") != "supported":
                raise ValueError("claim verification did not return supported")
            compiled.blocks.append(block)
            compiled.claims.append(manifest)
        except Exception:
            compiled.withheld.append({
                "claim_id": str(getattr(block, "id", "unknown")),
                "status": "withheld",
                "reason": "unsupported claim or invalid narrative block",
            })
    return compiled
