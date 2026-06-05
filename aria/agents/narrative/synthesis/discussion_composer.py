"""Compose the Integrated Biological Discussion as governed NarrativeBlocks.

The composer turns the deterministic pattern manifest into ``NarrativeBlock``s
(``integration.*``). Because they are blocks, they flow through the existing
governance: claim tiering (Claim Compiler), STRICT evidence verification (every
number/entity must be on the block's evidence card), devil's advocate, and
run-ledger linkage. Observational omics caps at the associative tier — no causal
language. Every number a claim states is attached as an evidence item, so the
strict verifier can confirm it.

Slice 1: bulk single-modality (within-contrast convergence, cross-contrast
convergence/divergence over a shared reference, reliability limits). Cross-modal
(RNA+ATAC) is Slice 2.
"""

from __future__ import annotations

from typing import Any

from aria.agents.narrative.types import NarrativeBlock, EvidenceItem, Caveat

_ASSOCIATIVE_CAVEAT = Caveat(
    "This is an observational association across analyses, not a demonstration of "
    "causation; shared regulators are inferred from co-differential expression, "
    "not from a perturbational rescue.",
    "info",
)


def _ev(label: str, value: Any) -> EvidenceItem:
    return EvidenceItem(label=label, value=value, source="biological_synthesis")


def compose_discussion_blocks(patterns: dict) -> list[NarrativeBlock]:
    """Return the integration NarrativeBlocks the data supports (may be empty)."""
    blocks: list[NarrativeBlock] = []
    within = patterns.get("within_contrast", []) or []
    cross = patterns.get("cross_contrast", []) or []
    reliability = patterns.get("reliability", {}) or {}
    n_contrasts = int(patterns.get("n_contrasts", 0) or 0)
    if n_contrasts == 0:
        return blocks

    # 1) Integrated signal: how many contrasts produced DE + concordant pathways.
    converging = [w for w in within if w.get("converges")]
    if converging:
        names = "; ".join(w["name"] for w in converging)
        blocks.append(NarrativeBlock(
            id="integration.signal",
            modality="integrated synthesis",
            analysis="integrated_signal",
            block_type="integration",
            title="Integrated biological signal",
            status="success",
            confidence="medium",
            claim=(
                f"Across {n_contrasts} analyzable contrast(s), {len(converging)} "
                f"produced significant differential expression accompanied by "
                f"functional pathway enrichment, indicating a coordinated "
                f"transcriptional response."
            ),
            evidence=[
                _ev("analyzable contrasts", n_contrasts),
                _ev("contrasts with DE + enrichment", len(converging)),
                _ev("contrasts", names),
            ],
            caveats=[_ASSOCIATIVE_CAVEAT],
            metrics={"n_contrasts": n_contrasts, "n_converging": len(converging)},
        ))

    # 2) Convergent evidence: only for pairs that share a reference level, where
    #    "same direction" is biologically meaningful (the detector flags this).
    for x in cross:
        if not x.get("shared_reference") or x.get("n_shared_genes", 0) <= 0:
            continue
        a, b = x["contrast_a"], x["contrast_b"]
        shared = x["n_shared_genes"]
        ev = [
            _ev("contrasts compared", f"{a}; {b}"),
            _ev("shared reference", x["shared_reference"]),
            _ev("shared DE genes", shared),
            _ev("shared enriched terms", x["n_shared_terms"]),
        ]
        direction_clause = ""
        if x.get("direction_known") and x.get("n_direction_concordant", -1) >= 0:
            conc = x["n_direction_concordant"]
            ev.append(_ev("shared genes, same direction", conc))
            direction_clause = (
                f", of which {conc} change in the same direction in both"
            )
        terms_clause = (
            f"; {x['n_shared_terms']} enriched term(s) are common to both"
            if x.get("n_shared_terms") else ""
        )
        blocks.append(NarrativeBlock(
            id=f"integration.convergent.{_safe(a)}__{_safe(b)}",
            modality="integrated synthesis",
            analysis="convergent_evidence",
            block_type="integration",
            title="Convergent evidence across analyses",
            status="success",
            confidence="medium",
            claim=(
                f"{a} and {b} share {shared} differentially expressed gene(s)"
                f"{direction_clause}{terms_clause}, consistent with a shared "
                f"transcriptional program relative to {x['shared_reference']}."
            ),
            evidence=ev,
            caveats=[_ASSOCIATIVE_CAVEAT],
            metrics={"n_shared_genes": shared,
                     "n_shared_terms": x["n_shared_terms"]},
        ))

    # 3) Divergent / contrast-specific signal (also reference-anchored).
    for x in cross:
        if not x.get("shared_reference"):
            continue
        a, b = x["contrast_a"], x["contrast_b"]
        spec_a, spec_b = x["n_specific_a"], x["n_specific_b"]
        if spec_a <= 0 and spec_b <= 0:
            continue
        discord_clause = ""
        if x.get("direction_known") and x.get("n_direction_discordant", 0) > 0:
            discord_clause = (
                f" {x['n_direction_discordant']} shared gene(s) move in opposite "
                f"directions (discordant evidence)."
            )
        blocks.append(NarrativeBlock(
            id=f"integration.divergent.{_safe(a)}__{_safe(b)}",
            modality="integrated synthesis",
            analysis="divergent_evidence",
            block_type="integration",
            title="Divergent or condition-specific signal",
            status="success",
            confidence="medium",
            claim=(
                f"{a} additionally shows {spec_a} contrast-specific differentially "
                f"expressed gene(s) and {b} shows {spec_b}, indicating partially "
                f"distinct programs.{discord_clause}"
            ),
            evidence=[
                _ev("contrasts compared", f"{a}; {b}"),
                _ev(f"{a} specific genes", spec_a),
                _ev(f"{b} specific genes", spec_b),
                _ev("discordant shared genes", x.get("n_direction_discordant", 0)),
            ],
            caveats=[_ASSOCIATIVE_CAVEAT],
            metrics={"n_specific_a": spec_a, "n_specific_b": spec_b},
        ))

    # 4) Limitations of interpretation — MANDATORY whenever any integrated claim
    #    was made.
    if len(blocks) > 0:
        ev = []
        bound = ""
        if reliability.get("min_power") is not None:
            ev.append(_ev("minimum power", reliability["min_power"]))
            ev.append(_ev("maximum power", reliability["max_power"]))
            bound = (
                f" Approximate statistical power ranged from "
                f"{reliability['min_power']} to {reliability['max_power']} "
                f"across contrasts."
            )
        low = reliability.get("low_power_contrasts") or []
        if low:
            ev.append(_ev("low-power contrasts", "; ".join(low)))
        blocks.append(NarrativeBlock(
            id="integration.limitations",
            modality="integrated synthesis",
            analysis="integration_limitations",
            block_type="limitation",
            title="Limits of interpretation",
            status="success",
            confidence="high",
            claim=(
                "These integrated results are associative and do not establish "
                "causality or direct regulation; cross-analysis agreement "
                "prioritizes candidates for independent experimental validation "
                "rather than proving a mechanism." + bound
            ),
            evidence=ev or [_ev("scope", "association_only")],
            caveats=[_ASSOCIATIVE_CAVEAT],
        ))

    return blocks


def _safe(name: str) -> str:
    out = "".join(ch if ch.isalnum() else "_" for ch in str(name))
    return out.strip("_") or "x"
