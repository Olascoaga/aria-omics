"""SPECULATIVE report section (ADR-057 S9): the visible face of the tier.

Gathers audited evidence from whichever single-modality adapters apply to the
run, runs the HypothesisAgent (opt-in, downstream of W-CLAIM + W-LEDGER), and
renders a clearly-walled report section: an explicit machine-generated header, a
model-provenance line, and each ranked hypothesis with its audited observation,
hedged mechanism, discriminating experiment, devils_advocate, and quarantine
node. This section is NEVER part of the audited claim manifest.
"""

from __future__ import annotations

import html as _html
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from .adapters import (
    bulk_atac_evidence,
    bulk_rna_evidence,
    scatac_evidence,
    scrna_evidence,
)

log = logging.getLogger(__name__)

SPECULATIVE_MANIFEST_SCHEMA = "aria.speculative_hypotheses.v1"
SPECULATIVE_MANIFEST_FILENAME = "speculative_hypotheses.json"

SPECULATIVE_HEADER = (
    "Machine-generated hypotheses. Not findings. Not verified against evidence. "
    "Each is grounded in audited results, but the proposed mechanism is "
    "speculative and must be tested by the experiment it names."
)


def _chromatin_signals(agent_results, run_ledger, exp_ctx):
    """Route the chromatin lane to the bulk ATAC or scATAC adapter (not both)."""
    try:
        from aria.agents.narrative.narrators.chromatin import (
            _is_bulk_atac,
            unwrap_chromatin_findings,
        )

        findings = unwrap_chromatin_findings(
            (agent_results or {}).get("chromatin_agent", {})
        )
        if not findings:
            return []
        if _is_bulk_atac(findings):
            return bulk_atac_evidence(agent_results, run_ledger, exp_ctx)
        return scatac_evidence(agent_results, run_ledger, exp_ctx)
    except Exception as exc:
        # A chromatin-adapter failure must NOT silently masquerade as "no
        # evidence"; log it so the missing speculative input is visible.
        log.warning(
            "Speculative chromatin evidence adapter failed; treating as no "
            "chromatin evidence: %s",
            exc,
            exc_info=True,
        )
        return []


def gather_evidence(agent_results, run_ledger=None, exp_ctx=None) -> list:
    """Collect EvidenceSignal items from every applicable audited modality."""
    agent_results = agent_results or {}
    signals: list = []
    if "bulk_rna_agent" in agent_results:
        signals += bulk_rna_evidence(agent_results, run_ledger, exp_ctx)
    if "scrna_agent" in agent_results:
        signals += scrna_evidence(agent_results, run_ledger, exp_ctx)
    if "chromatin_agent" in agent_results:
        signals += _chromatin_signals(agent_results, run_ledger, exp_ctx)
    return signals


def build_speculative_section(
    agent_results,
    run_ledger=None,
    exp_ctx=None,
    *,
    proposer=None,
    w_claim_passed: bool = False,
    w_ledger_passed: bool = False,
):
    """Build the SPECULATIVE section result, or None when there is nothing to show.

    Opt-in is the caller's responsibility (e.g. ``exp_ctx['enable_hypotheses']``).
    Returns the HypothesisAgent output augmented with a ``header``; returns None
    when no audited evidence is available so no empty section is emitted.

    H13/F5: ``w_claim_passed`` / ``w_ledger_passed`` are fail-closed (default
    ``False``) — the real caller (``report_builder._speculative_verification_state``)
    passes the run's actual verification state; a caller that omits them gets no
    speculation rather than a silent fail-open.
    """
    from aria.agents.hypothesis_agent import HypothesisAgent

    signals = gather_evidence(agent_results, run_ledger, exp_ctx)
    if not signals:
        return None
    out = HypothesisAgent(proposer=proposer).generate(
        signals,
        run_ledger,
        exp_ctx,
        w_claim_passed=w_claim_passed,
        w_ledger_passed=w_ledger_passed,
    )
    out["header"] = SPECULATIVE_HEADER
    return out


def build_speculative_manifest(
    section: dict | None, *, generated_utc: str | None = None
) -> dict | None:
    """Build the auditable, NON-PROMOTABLE speculative manifest from a section.

    A structured record of the speculative layer for a run: the quarantine nodes,
    the model provenance, the ranked hypotheses, and — when nothing was emitted —
    the honest reason (gate withheld / per-gate rejections). It is deliberately a
    SEPARATE artifact from the audited claim manifest and is stamped
    ``promotable: False`` + the ``hypothesis://`` tier so it can never be mistaken
    for an audited claim. Returns None when there was no speculative section.
    """
    if not section:
        return None
    hypotheses = section.get("hypotheses") or []
    provenance = (hypotheses[0].get("provenance") if hypotheses else None) or {}
    return {
        "schema": SPECULATIVE_MANIFEST_SCHEMA,
        "tier": "SPECULATIVE",
        "promotable": False,
        "note": (
            "Machine-generated speculative hypotheses (ADR-057). NOT part of the "
            "audited claim manifest and mechanically non-promotable "
            "(hypothesis:// namespace)."
        ),
        "generated_utc": (
            generated_utc
            if generated_utc is not None
            else datetime.now(timezone.utc).isoformat()
        ),
        "ran": bool(section.get("ran")),
        "honest_null": bool(section.get("honest_null")),
        "reason": section.get("reason"),
        "null_reason": section.get("null_reason"),
        "null_summary": section.get("null_summary") or {},
        "proposer_diagnostics": section.get("proposer_diagnostics"),
        "n_evidence": section.get("n_evidence"),
        "n_candidates": section.get("n_candidates"),
        "provenance": provenance,
        "quarantine": section.get("quarantine") or [],
        "hypotheses": hypotheses,
        "rejected": section.get("rejected") or [],
    }


def persist_speculative_manifest(
    section: dict | None, report_dir, *, reproducible: bool = False
) -> Path | None:
    """Write ``speculative_hypotheses.json`` next to the report. Returns its path.

    Persists whenever a speculative section was engaged (including honest-null and
    gate-withheld runs — the attempt + reason is itself auditable provenance). No
    file is written when there was no section (no evidence / opt-out). Under
    ``reproducible`` the timestamp is redacted so the artifact stays byte-identical.
    """
    if not section or report_dir is None:
        return None
    generated = (
        "<timestamp redacted for byte-identity>" if reproducible else None
    )
    manifest = build_speculative_manifest(section, generated_utc=generated)
    path = Path(report_dir) / SPECULATIVE_MANIFEST_FILENAME
    path.write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8"
    )
    return path


def _esc(value) -> str:
    return _html.escape(str(value if value is not None else ""))


def _render_hypothesis(hyp: dict, rank: int) -> str:
    exp = hyp.get("experiment") or {}
    da = hyp.get("devils_advocate") or {}
    confounds = ", ".join(str(c) for c in (da.get("confounds") or [])) or "none"
    rank_ev = hyp.get("rank_evidence") or {}
    competing = ", ".join(str(c) for c in (hyp.get("competing_with") or []))
    competing_html = (
        f"<p style='color:var(--muted);font-size:0.85rem'>competes with: "
        f"{_esc(competing)}</p>"
        if competing
        else ""
    )
    return (
        "<div class='card' style='margin:0.5rem 0'>"
        f"<h3>{rank}. {_esc(hyp.get('mechanism'))}</h3>"
        f"<p><strong>Arises from:</strong> {_esc(', '.join(hyp.get('observation_refs') or []))} "
        f"(entities: {_esc(', '.join(hyp.get('entities') or []))})</p>"
        f"<p><strong>Discriminating experiment:</strong> "
        f"{_esc(exp.get('perturbation'))} &rarr; {_esc(exp.get('readout'))}; "
        f"predicts <em>{_esc(exp.get('predicted_direction'))}</em>; "
        f"refuted if {_esc(exp.get('refuting_outcome'))}</p>"
        f"<p><strong>Devil's advocate:</strong> simpler explanation &mdash; "
        f"{_esc(da.get('simpler_explanation'))}; confounds &mdash; {_esc(confounds)}</p>"
        f"{competing_html}"
        f"<p style='color:var(--muted);font-size:0.85rem'>ranking basis: "
        f"{_esc(rank_ev.get('n_independent_lines'))} independent line(s), "
        f"mean normalized effect {_esc(rank_ev.get('mean_effect_norm'))}, "
        f"confound load {_esc(rank_ev.get('evidence_caveat_load'))} &middot; "
        f"quarantine node {_esc(hyp.get('ledger_node'))}</p>"
        "</div>"
    )


def _render_gate_blocked(section: dict) -> str:
    """Visible note when the causal gate withheld the section (rail #1).

    The agent did not run because the run's audited claims failed W-CLAIM/
    W-LEDGER verification. That is a governance signal, not nothing — render it
    instead of silently dropping the section.
    """
    return "\n".join(
        [
            "<h2>Speculative Hypotheses "
            "<span style='font-size:0.7rem;background:var(--amber);color:#000;"
            "padding:0.1rem 0.4rem;border-radius:0.25rem'>SPECULATIVE</span></h2>",
            "<div class='card' style='border-left:4px solid var(--amber)'>",
            f"<p><strong>{_esc(SPECULATIVE_HEADER)}</strong></p>",
            "<p>No hypotheses were generated: the run's audited claims did not "
            "pass W-CLAIM/W-LEDGER verification, so the speculative layer is "
            "withheld (a hypothesis must not build on an unverified claim).</p>",
            "</div>",
        ]
    )


def _render_null_summary(section: dict) -> str:
    """Per-gate breakdown of why every candidate was rejected (honest-null)."""
    summary = section.get("null_summary") or {}
    if not summary:
        return ""
    parts = ", ".join(
        f"{_esc(gate)}: {_esc(count)}"
        for gate, count in sorted(summary.items())
    )
    n = section.get("n_candidates")
    lead = (
        f"{_esc(n)} candidate(s) proposed; none survived the publication gates"
        if n is not None
        else "No candidate survived the publication gates"
    )
    return (
        f"<p style='color:var(--muted);font-size:0.85rem'>{lead} "
        f"(rejections by gate &mdash; {parts}).</p>"
    )


def render_speculative_section_html(section: dict | None) -> str:
    """Render the SPECULATIVE section to HTML, or '' when there is nothing to show."""
    if not section:
        return ""
    if not section.get("ran"):
        # The causal gate withheld generation (verification did not pass). Make
        # it visible rather than silently emitting nothing.
        if section.get("reason") == "verification_gate_not_passed":
            return _render_gate_blocked(section)
        return ""
    hypotheses = section.get("hypotheses") or []
    prov = (hypotheses[0].get("provenance") if hypotheses else None) or {}
    model = prov.get("model_label") or "machine"
    body = [
        "<h2>Speculative Hypotheses "
        "<span style='font-size:0.7rem;background:var(--amber);color:#000;"
        "padding:0.1rem 0.4rem;border-radius:0.25rem'>SPECULATIVE</span></h2>",
        "<div class='card' style='border-left:4px solid var(--amber)'>",
        f"<p><strong>{_esc(section.get('header'))}</strong></p>",
        f"<p style='color:var(--muted);font-size:0.85rem'>generator: {_esc(model)}; "
        "this section is not part of the audited claim manifest and cannot be "
        "promoted to a finding.</p>",
    ]
    if not hypotheses:
        # Distinguish a true honest-null (model declined / nothing survived the
        # gates) from a generation failure (the model answered but the response
        # could not be parsed, almost always a token-budget truncation). Never
        # present a generation failure as "no defensible hypothesis".
        reason = section.get("null_reason")
        if reason in ("parse_error", "malformed_items"):
            body.append(
                "<p>The generator returned a response that could not be parsed "
                "(likely truncated); no hypothesis is shown. This is a generation "
                "issue, not a conclusion that the evidence supports none. "
                "Re-running may resolve it.</p>"
            )
        elif reason in ("empty_response", "no_evidence"):
            body.append(
                "<p>The generator returned no usable response; no hypothesis is "
                "shown (generation issue, not a conclusion about the evidence).</p>"
            )
        else:
            body.append(
                "<p>No defensible hypothesis from the audited evidence "
                "(honest-null).</p>"
            )
            # Surface WHY: the per-gate rejection breakdown, so an honest-null is
            # explainable (e.g. grounding rejected an un-measured entity) and not
            # an opaque "nothing here".
            summary_html = _render_null_summary(section)
            if summary_html:
                body.append(summary_html)
    else:
        for i, hyp in enumerate(hypotheses, start=1):
            body.append(_render_hypothesis(hyp, i))
    body.append("</div>")
    return "\n".join(body)
