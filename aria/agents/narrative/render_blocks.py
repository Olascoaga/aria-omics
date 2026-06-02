"""HTML composer for structured narrative blocks."""

from __future__ import annotations

import base64
import html
from pathlib import Path

from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.types import NarrativeBlock
from aria.agents.narrative.validators import (
    validate_blocks,
    NarrativeValidationError,
    find_causal_language,
    collect_named_entities,
)
from aria.agents.narrative.evidence_verifier import (
    EvidenceVerificationError,
    verify_block_claim_support,
)


BLOCK_ORDER = {
    "qc": 0,
    "result": 1,
    "exploratory": 2,
    "limitation": 3,
    "error": 4,
    "method": 5,
}


def render_blocks(blocks: list[NarrativeBlock],
                  report_dir: str | Path | None = None) -> str:
    blocks = validate_blocks(
        list(blocks or []),
        base_dir=report_dir,
        check_files=report_dir is not None,
    )
    ordered = sorted(
        blocks,
        key=lambda b: (
            _group_key(b),
            BLOCK_ORDER.get(b.block_type, 99),
            b.id,
        ),
    )
    return "\n".join(_render_block(block, report_dir=report_dir)
                     for block in ordered)


def group_blocks_by_prefix(blocks: list[NarrativeBlock]) -> dict[str, list[NarrativeBlock]]:
    grouped: dict[str, list[NarrativeBlock]] = {}
    for block in blocks or []:
        grouped.setdefault(_group_key(block), []).append(block)
    return grouped


def _group_key(block: NarrativeBlock) -> str:
    return str(block.id).split(".", 1)[0]


def _render_block(block: NarrativeBlock,
                  report_dir: str | Path | None = None) -> str:
    status_class = "low" if block.status != "success" else block.confidence
    if status_class == "insufficient":
        status_class = "insuff"
    evidence_html = _render_evidence(block)
    caveats_html = _render_caveats(block)
    figures_html = _render_figures(block, report_dir=report_dir)
    tables_html = _render_tables(block)
    error_html = (
        f"<p class='warning'><strong>Error:</strong> {html.escape(str(block.error))}</p>"
        if block.error else ""
    )
    raw_prose = compose_block_prose(block) or block.claim or ""

    # F-SCI-CAUSAL (audit 2026-05-28): the causal guard in validate_blocks only
    # scans claim+evidence. Scan the FINAL composed prose too — that is what the
    # reader sees — and surface a visible warning if unlicensed causal phrasing
    # slipped in (unless the block declares explicit causal evidence).
    block_warnings = list(block.warnings or [])
    if not block.metadata.get("causal_evidence"):
        # B5 (audit 2026-05-29): redact external named entities (DB term names,
        # gene symbols) that compose_prose interpolates into the prose so a GO/
        # Reactome term like "TP53 Regulates Transcription…" does not trip the
        # guard. Only ARIA's authored phrasing should be evaluated.
        causal_hit = find_causal_language(
            raw_prose, exclude=collect_named_entities(block)
        )
        if causal_hit:
            warn = (
                f"Associative result described with causal language "
                f"('{causal_hit}'); interpret as correlation, not causation."
            )
            if warn not in block_warnings:
                block_warnings.append(warn)

    warnings_html = ""
    if block_warnings:
        warnings_html = (
            "<ul style='margin:0.4rem 0 0.4rem 1.1rem;color:var(--amber)'>"
            + "".join(f"<li>{html.escape(str(w))}</li>" for w in block_warnings)
            + "</ul>"
        )
    prose = html.escape(raw_prose)
    try:
        block.metadata["claim_verification"] = verify_block_claim_support(
            block, raw_prose, strict=True
        )
    except EvidenceVerificationError as exc:
        raise NarrativeValidationError(str(exc)) from exc
    tier_badge = _claim_tier_badge(block)
    return f"""
<section class="narrative-block" data-block-id="{html.escape(block.id)}"
         style="border-top:1px solid var(--border);padding-top:0.9rem;margin-top:0.9rem">
  <h4>{html.escape(block.title)}
    <span class="badge {status_class}" style="margin-left:0.35rem">
      {html.escape(block.status)} / {html.escape(block.confidence)}
    </span>{tier_badge}
  </h4>
  {f"<p>{prose}</p>" if prose else ""}
  {error_html}
  {evidence_html}
  {caveats_html}
  {warnings_html}
  {figures_html}
  {tables_html}
</section>"""


def _claim_tier_badge(block: NarrativeBlock) -> str:
    """X14: render the evidence-tier badge for a block's claim.

    Uses the precomputed `block.metadata['claim']` when present (set by
    annotate_claim_tiers); falls back to classifying on the fly so render is
    correct even when called standalone (e.g. tests, offline harness).
    """
    claim = block.metadata.get("claim")
    if not claim:
        try:
            from aria.agents.narrative.claim_compiler import classify_claim
            claim = classify_claim(block).as_dict()
        except Exception:
            return ""
    tier = claim.get("tier")
    if not tier or tier == "descriptive":
        return ""
    licensed = claim.get("licensed_language", "associative")
    label = tier.replace("_", " ")
    title = html.escape(str(claim.get("rationale") or ""))
    return (
        f"<span class='badge' title='{title}' "
        f"style='margin-left:0.35rem;background:var(--border);color:var(--muted)'>"
        f"evidence: {html.escape(label)} · {html.escape(licensed)}</span>"
    )


def _render_evidence(block: NarrativeBlock) -> str:
    if not block.evidence:
        return ""
    rows = []
    for ev in block.evidence:
        value = "" if ev.value is None else ev.value
        source = ev.source or ""
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(ev.label))}</td>"
            f"<td><code>{html.escape(str(value))}</code></td>"
            f"<td>{html.escape(str(source))}</td>"
            "</tr>"
        )
    return (
        "<details style='margin-top:0.45rem'>"
        "<summary style='font-size:0.85rem;color:var(--muted);cursor:pointer'>"
        "Structured evidence</summary>"
        "<table style='font-size:0.84rem;margin-top:0.35rem'>"
        "<tr><th>Evidence</th><th>Value</th><th>Source</th></tr>"
        + "".join(rows)
        + "</table>"
        "</details>"
    )


def _render_caveats(block: NarrativeBlock) -> str:
    if not block.caveats:
        return ""
    items = "".join(
        f"<li><strong>{html.escape(c.severity)}</strong>: "
        f"{html.escape(c.text)}</li>"
        for c in block.caveats
    )
    return (
        "<div class='warning'><strong>Caveats</strong>"
        f"<ul style='margin:0.4rem 0 0 1.1rem'>{items}</ul></div>"
    )


def _render_figures(block: NarrativeBlock,
                    report_dir: str | Path | None = None) -> str:
    html_parts = []
    for fig in block.figures or []:
        raw = fig.get("path")
        if not raw:
            continue
        uri = _image_uri(raw, report_dir=report_dir)
        caption = html.escape(str(fig.get("caption") or fig.get("id") or raw))
        if uri:
            html_parts.append(
                f"<figure><img src='{uri}' alt='{caption}'><figcaption>"
                f"{caption}</figcaption></figure>"
            )
    return "\n".join(html_parts)


def _render_tables(block: NarrativeBlock) -> str:
    if not block.tables:
        return ""
    links = []
    for table in block.tables:
        raw = table.get("path")
        if not raw:
            continue
        label = html.escape(str(table.get("label") or table.get("id") or raw))
        links.append(f"<a href='{html.escape(str(raw))}'>{label}</a>")
    if not links:
        return ""
    return (
        "<p style='font-size:0.85rem;color:var(--muted)'>Tables: "
        + " &middot; ".join(links)
        + "</p>"
    )


def _image_uri(raw: str, report_dir: str | Path | None = None) -> str:
    path = Path(raw)
    if not path.is_absolute() and report_dir is not None:
        path = Path(report_dir) / path
    if not path.exists():
        return ""
    suffix = path.suffix.lower()
    mime = "image/png"
    if suffix == ".svg":
        mime = "image/svg+xml"
    elif suffix in {".jpg", ".jpeg"}:
        mime = "image/jpeg"
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"
