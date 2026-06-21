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

from .adapters import (
    bulk_atac_evidence,
    bulk_rna_evidence,
    scatac_evidence,
    scrna_evidence,
)

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
    except Exception:
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
    w_claim_passed: bool = True,
    w_ledger_passed: bool = True,
):
    """Build the SPECULATIVE section result, or None when there is nothing to show.

    Opt-in is the caller's responsibility (e.g. ``exp_ctx['enable_hypotheses']``).
    Returns the HypothesisAgent output augmented with a ``header``; returns None
    when no audited evidence is available so no empty section is emitted.
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


def _esc(value) -> str:
    return _html.escape(str(value if value is not None else ""))


def _render_hypothesis(hyp: dict, rank: int) -> str:
    exp = hyp.get("experiment") or {}
    da = hyp.get("devils_advocate") or {}
    confounds = ", ".join(str(c) for c in (da.get("confounds") or [])) or "none"
    rank_ev = hyp.get("rank_evidence") or {}
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
        f"<p style='color:var(--muted);font-size:0.85rem'>ranking basis: "
        f"{_esc(rank_ev.get('n_independent_lines'))} independent line(s), "
        f"mean |effect| {_esc(rank_ev.get('mean_abs_effect'))}, "
        f"confound load {_esc(rank_ev.get('evidence_caveat_load'))} &middot; "
        f"quarantine node {_esc(hyp.get('ledger_node'))}</p>"
        "</div>"
    )


def render_speculative_section_html(section: dict | None) -> str:
    """Render the SPECULATIVE section to HTML, or '' when there is nothing to show."""
    if not section or not section.get("ran"):
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
        body.append(
            "<p>No defensible hypothesis from the audited evidence "
            "(honest-null).</p>"
        )
    else:
        for i, hyp in enumerate(hypotheses, start=1):
            body.append(_render_hypothesis(hyp, i))
    body.append("</div>")
    return "\n".join(body)
