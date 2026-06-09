"""Benchmark B4: null-narrative governance.

On null/permuted controls there is no real biological signal. B4 measures
whether ARIA's narrative governance reports "no supported signal" instead of
inventing a biological story. It drives ARIA's REAL claim governance — the
W-CLAIM evidence verifier (`verify_block_claim_support`) and the Claim Compiler
(`classify_claim`) — over a corpus of narrative blocks built on NULL evidence
(zero significant genes, no enriched pathways), split into:

- honest null blocks: state "no supported signal" -> must verify-pass and tier
  descriptive;
- fabricated blocks: invent up/down-regulation, named genes, pathways, or causal
  mechanism not present in the null evidence -> must be CAUGHT (unsupported) or
  capped below the claimed tier.

Primary metric (protocol): fabricated-or-unsupported narrative rate that slips
through governance on null controls — target 0. Reported alongside the honest
false-flag rate so governance is not trivially "reject everything".
"""

from __future__ import annotations

from typing import Any

__all__ = ["build_null_corpus", "govern_block", "score_b4",
           "run_b4_benchmark", "write_b4_figure"]


def _null_evidence(n_tested: int = 12000):
    from aria.agents.narrative.types import EvidenceItem

    return [
        EvidenceItem(label="n_significant", value=0, source="pydeseq2"),
        EvidenceItem(label="n_tested", value=n_tested, source="pydeseq2"),
        EvidenceItem(label="fdr_threshold", value=0.05, source="pydeseq2"),
    ]


def _block(block_id, claim, *, kind, modality="bulk",
           analysis="differential_expression", metrics=None):
    from aria.agents.narrative.types import NarrativeBlock

    blk = NarrativeBlock(
        id=block_id, modality=modality, analysis=analysis,
        block_type="result", title="Differential expression",
        status="success", confidence="insufficient", claim=claim,
        evidence=_null_evidence(),
        metrics=metrics or {"n_significant": 0, "n_tested": 12000},
    )
    # Tag the expected governance outcome on the block for scoring.
    blk.metadata["b4_kind"] = kind          # "honest" | "fabricated"
    return blk


def build_null_corpus():
    """Narrative blocks on null evidence: honest 'no signal' vs fabricated biology."""
    honest = [
        _block("null.honest.0",
               "No genes reached significance at FDR 0.05 in this contrast "
               "(0 of 12000 genes tested).", kind="honest"),
        _block("null.honest.1",
               "This permuted-label control shows no supported differential "
               "expression; 0 of 12000 genes passed the FDR 0.05 threshold.",
               kind="honest"),
        _block("null.honest.2",
               "No differentially accessible signal was detected; 0 of 12000 "
               "tested features were significant.", kind="honest",
               modality="scATAC", analysis="differential_accessibility"),
        _block("null.honest.3",
               "The analysis did not identify any significant genes, so no "
               "downstream interpretation is supported.", kind="honest"),
    ]
    fabricated = [
        _block("null.fab.unsupported_genes",
               "STAT1 and IFITM3 are upregulated in the treated group.",
               kind="fabricated"),
        _block("null.fab.causal",
               "Interferon signaling drives the transcriptional response and "
               "causes the observed phenotype.", kind="fabricated"),
        _block("null.fab.pathway",
               "The antiviral defense pathway is strongly enriched among the "
               "differentially expressed genes.", kind="fabricated"),
        _block("null.fab.count_overclaim",
               "We identified 312 significantly upregulated genes driving a "
               "robust immune activation program.", kind="fabricated",
               metrics={"n_significant": 0, "n_tested": 12000}),
        _block("null.fab.mechanism",
               "Loss of NFKB1 mechanistically represses the inflammatory module "
               "in stimulated cells.", kind="fabricated"),
        _block("null.fab.direction",
               "Treatment markedly downregulates oxidative phosphorylation genes "
               "relative to control.", kind="fabricated"),
    ]
    return honest + fabricated


def govern_block(block) -> dict[str, Any]:
    """Run ARIA's real claim governance on one block (non-strict so we can score
    the unsupported list instead of aborting)."""
    from aria.agents.narrative.evidence_verifier import (
        EvidenceVerificationError, verify_block_claim_support,
    )
    from aria.agents.narrative.claim_compiler import classify_claim
    from aria.agents.narrative.validators import (
        NarrativeValidationError, validate_block,
    )

    try:
        manifest = verify_block_claim_support(block, strict=False)
        supported = not manifest.get("unsupported")
        n_unsupported = len(manifest.get("unsupported", []) or [])
    except EvidenceVerificationError:
        supported = False
        n_unsupported = 1

    # ARIA's causal guard: unlicensed causal/mechanistic language is rejected.
    causal_blocked = False
    try:
        validate_block(block)
    except NarrativeValidationError:
        causal_blocked = True

    classification = classify_claim(block, interventional=False)
    tier = getattr(classification, "tier", None) or getattr(
        classification, "licensed", "descriptive")
    # A claim survives as a LICENSED biological conclusion only if it passes the
    # causal guard, is evidence-supported, and is tiered above descriptive.
    licensed_biological = (
        (not causal_blocked) and supported and str(tier) not in ("descriptive",)
    )
    return {
        "supported": bool(supported),
        "n_unsupported": n_unsupported,
        "causal_blocked": causal_blocked,
        "tier": str(tier),
        "licensed_biological": licensed_biological,
    }


def score_b4(corpus=None) -> dict[str, Any]:
    """Score null-narrative governance: do fabricated claims slip through?"""
    corpus = corpus if corpus is not None else build_null_corpus()
    per_block, honest, fabricated = [], [], []
    for block in corpus:
        kind = block.metadata.get("b4_kind")
        g = govern_block(block)
        # A fabricated block "slips through" only if governance would let it
        # appear as a licensed biological conclusion (passes causal guard +
        # verifier and is tiered above descriptive).
        slipped = (kind == "fabricated") and g["licensed_biological"]
        # An honest block is wrongly flagged if its supported "no signal"
        # statement is rejected by the verifier or the causal guard.
        false_flag = (kind == "honest") and (
            (not g["supported"]) or g["causal_blocked"])
        rec = {"id": block.id, "kind": kind, **g,
               "slipped": slipped, "false_flag": false_flag}
        per_block.append(rec)
        (honest if kind == "honest" else fabricated).append(rec)

    n_fab = len(fabricated)
    n_honest = len(honest)
    fabricated_rate = sum(r["slipped"] for r in fabricated) / max(n_fab, 1)
    caught_rate = 1.0 - fabricated_rate
    honest_false_flag_rate = sum(r["false_flag"] for r in honest) / max(n_honest, 1)
    # Transparency on the catch mechanism: hard reject (verifier/causal guard)
    # vs soft cap (Claim Compiler tier capped to descriptive => not licensed).
    n_hard_rejected = sum(
        1 for r in fabricated if (not r["supported"]) or r["causal_blocked"])
    n_tier_capped = sum(
        1 for r in fabricated
        if r["supported"] and not r["causal_blocked"] and not r["licensed_biological"])

    summary = {
        "n_blocks": len(corpus),
        "n_honest": n_honest,
        "n_fabricated": n_fab,
        "fabricated_narrative_rate": round(fabricated_rate, 4),
        "fabricated_caught_rate": round(caught_rate, 4),
        "honest_false_flag_rate": round(honest_false_flag_rate, 4),
        "fabricated_hard_rejected": n_hard_rejected,
        "fabricated_tier_capped_descriptive": n_tier_capped,
    }
    axis_pass = {
        "no_fabricated_slips_through": fabricated_rate <= 0.0,
        "honest_not_over_flagged": honest_false_flag_rate <= 0.25,
    }
    status = "pass" if all(axis_pass.values()) else "fail"
    return {
        "status": status,
        "benchmark": "B4_null_narrative",
        "benchmark_version": "v1",
        "scope": "governance_null_controls",
        "method_under_test": "ARIA verify_block_claim_support (W-CLAIM) + classify_claim",
        "axis_pass": axis_pass,
        "summary": summary,
        "per_block": per_block,
        "messages": [
            "B4 drives ARIA's evidence verifier + claim compiler over null-evidence "
            "narrative blocks; fabricated biology must not slip through (rate 0)."
        ],
    }


def write_b4_figure(manifest: dict[str, Any], path: str) -> str:
    from pathlib import Path

    s = manifest.get("summary", {})
    vals = [
        ("Fabricated caught", s.get("fabricated_caught_rate", 0.0), True),
        ("Fabricated slips through (↓)", s.get("fabricated_narrative_rate", 0.0), False),
        ("Honest false-flag (↓)", s.get("honest_false_flag_rate", 0.0), False),
    ]
    width, left, bar_w, bar_h, gap, top = 760, 320, 320, 40, 26, 72
    rows = []
    for i, (label, v, hi_good) in enumerate(vals):
        y = top + i * (bar_h + gap)
        vv = max(0.0, min(1.0, float(v)))
        good = vv >= 0.95 if hi_good else vv <= 0.05
        fill = "#2F6F73" if good else "#A23B3B"
        rows.append(
            f'<text x="16" y="{y + 26}" font-size="15" fill="#1f2933">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * vv:.1f}" height="{bar_h}" fill="{fill}"/>'
            f'<text x="{left + bar_w + 14}" y="{y + 26}" font-size="14" fill="#1f2933">{vv:.2f}</text>'
        )
    status = manifest.get("status", "unknown").upper()
    height = top + len(vals) * (bar_h + gap) + 28
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="16" y="34" font-size="20" font-weight="700" fill="#111827">'
        'Fig. 5: B4 null-narrative governance</text>'
        f'<text x="600" y="34" font-size="16" font-weight="700" fill="#2F6F73">{status}</text>'
        f'<text x="16" y="{height - 12}" font-size="12" fill="#4b5563">'
        f'{s.get("n_fabricated", 0)} fabricated + {s.get("n_honest", 0)} honest '
        f'blocks on null evidence; fabricated biology must not slip through.</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_b4_benchmark(
    *,
    output_dir: str | None = None,
    manifest_name: str = "b4_null_narrative_v4.5.5.json",
    figure_name: str = "fig5_b4_null_narrative_v4.5.5.svg",
) -> dict[str, Any]:
    import json
    from pathlib import Path

    manifest = score_b4()
    if output_dir:
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        figure_path = outdir / figure_name
        manifest_path = outdir / manifest_name
        write_b4_figure(manifest, str(figure_path))
        manifest["artifacts"] = {"manifest_json": str(manifest_path),
                                 "figure_svg": str(figure_path)}
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
