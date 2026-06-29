"""Grounding verifier (ADR-057 rail #7): zero invention over the facts.

Every entity a hypothesis names must resolve to a real audited
``EvidenceSignal``, and every observation it arises from must cite a run-ledger
node that actually ran. A hypothesis that fails is REJECTED, never caveated.
This is the mechanical wall that lets the LLM be free over the *connection*
while the *facts* stay real — "LLM proposes, code guarantees" on the most
dangerous layer.

Reuses the existing W-LEDGER machinery (``run_ledger._node_index`` +
``_NOT_RUN_STATUSES``) so a hypothesis cannot arise from an analysis the run
marked not-run/skipped/error — the same contradiction W-LEDGER catches for
claims. S1 does not modify the ledger; it only reads it.

The wall also guards the rendered ``mechanism`` prose, not only the structured
``entities`` field: an LLM can list grounded entities yet still smuggle an
un-measured entity into the free-text mechanism the reader actually sees. This
reuses W-CLAIM's own named-entity check (``evidence_verifier._claim_entities``)
so a gene-like token named in the mechanism but absent from the audited evidence
is rejected, the same way W-CLAIM rejects it for an audited claim.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field

from aria.agents.narrative.evidence_verifier import _claim_entities
from aria.agents.narrative.run_ledger import _NOT_RUN_STATUSES, _node_index

from .types import EvidenceSignal, Hypothesis


def _norm(entity: str) -> str:
    return str(entity or "").strip().lower()


def build_evidence_index(
    signals: list[EvidenceSignal],
) -> dict[str, EvidenceSignal]:
    """Index audited evidence by normalized entity — the grounding universe.

    Only real ``EvidenceSignal`` items with a non-empty entity enter the index;
    anything else is ignored (no fabrication of a grounding target).
    """
    index: dict[str, EvidenceSignal] = {}
    for sig in signals or []:
        if isinstance(sig, EvidenceSignal) and _norm(sig.entity):
            # First-wins: deterministic representative. An entity can now carry
            # several context-distinct signals (H4); use build_signals_by_entity
            # when ALL of them matter (caveat union, independent-line counting).
            index.setdefault(_norm(sig.entity), sig)
    return index


def build_signals_by_entity(
    signals: list[EvidenceSignal],
) -> dict[str, list[EvidenceSignal]]:
    """Index audited evidence by entity, KEEPING every context-distinct signal.

    Unlike :func:`build_evidence_index` (one representative per entity), this
    preserves the full list so a gene measured in two contrasts contributes both
    its confounds and both its converging lines (H4: no destructive dedup).
    """
    out: dict[str, list[EvidenceSignal]] = {}
    for sig in signals or []:
        if isinstance(sig, EvidenceSignal) and _norm(sig.entity):
            out.setdefault(_norm(sig.entity), []).append(sig)
    return out


@dataclass
class GroundingResult:
    """Outcome of grounding one hypothesis against audited evidence."""

    grounded: bool
    missing_entities: list[str] = field(default_factory=list)
    not_run_refs: list[dict] = field(default_factory=list)
    ungrounded_prose_entities: list[str] = field(default_factory=list)
    misattributed_refs: list[str] = field(default_factory=list)
    # H15: signal-level grounding of the hypothesis's observed_claims.
    unknown_signals: list[str] = field(default_factory=list)
    misattributed_signals: list[str] = field(default_factory=list)
    contradicting_claims: list[dict] = field(default_factory=list)
    missing_observed_claims: bool = False
    vacuous: bool = False
    reason: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


# H15 (Codex blocker 2): direction vocabularies to reconstruct, from code, what a
# stated reading of a signal means and whether it contradicts the audited
# direction. Best-effort and deterministic: a stated_direction carrying BOTH an
# up- and a down-token (or neither) normalises to ``na`` (ambiguous → not a
# contradiction). Only a DEFINITE opposite (up vs down) is rejected, so a
# speculative downstream effect is never caught here — only a false restatement of
# the audited measurement.
_UP_DIRECTION_TOKENS = {
    "up", "increase", "increased", "increases", "increasing", "higher", "high",
    "elevated", "elevation", "gain", "gained", "enriched", "enrichment",
    "upregulated", "upregulation", "induced", "induction", "activated", "more",
    "accumulate", "accumulation", "accumulated", "positive",
}
_DOWN_DIRECTION_TOKENS = {
    "down", "decrease", "decreased", "decreases", "decreasing", "lower", "low",
    "reduced", "reduction", "depleted", "depletion", "loss", "lost",
    "downregulated", "downregulation", "repressed", "suppressed", "silenced",
    "diminished", "less", "fewer", "negative",
}


def _normalize_direction(text: str) -> str:
    """Map a free-text stated direction to ``up`` / ``down`` / ``na``.

    Ambiguous (both polarities, or neither) → ``na``, which never contradicts an
    audited direction. Hyphens are split so ``down-regulated`` matches.
    """
    toks = set(re.split(r"[^a-z0-9]+", str(text or "").lower()))
    up = bool(toks & _UP_DIRECTION_TOKENS)
    down = bool(toks & _DOWN_DIRECTION_TOKENS)
    if up and not down:
        return "up"
    if down and not up:
        return "down"
    return "na"


def _verify_observed_claims(
    hypothesis: Hypothesis,
    signals: list[EvidenceSignal],
    declared: set[str],
) -> tuple[list[str], list[str], list[dict], bool]:
    """Validate the signal-level observed_claims (H15).

    Returns ``(unknown_signals, misattributed_signals, contradicting, missing)``:

    - ``unknown_signals``: a cited ``signal_id`` not in the audited universe — the
      hypothesis claims to read a measurement that does not exist.
    - ``misattributed_signals``: a cited signal whose entity the hypothesis never
      named — a reading attributed to an entity it is not about.
    - ``contradicting``: a stated_direction that is the DEFINITE opposite of the
      cited signal's audited direction (e.g. "increased" over a ``down`` signal).
    - ``missing``: the hypothesis cites NO audited signal at all — it must declare
      at least one measurement it reads (faithful), so a directional claim can
      never live only in un-anchored prose.
    """
    by_id = {
        sig.signal_id: sig
        for sig in (signals or [])
        if isinstance(sig, EvidenceSignal) and sig.signal_id
    }
    unknown: list[str] = []
    misattributed: list[str] = []
    contradicting: list[dict] = []
    claims = hypothesis.observed_claims or []
    for claim in claims:
        sid = str((claim or {}).get("signal_id", ""))
        sig = by_id.get(sid)
        if sig is None:
            unknown.append(sid)
            continue
        if _norm(sig.entity) not in declared:
            misattributed.append(sid)
            continue
        stated = _normalize_direction((claim or {}).get("stated_direction", ""))
        audited = sig.direction if sig.direction in ("up", "down") else "na"
        if stated in ("up", "down") and audited in ("up", "down") and stated != audited:
            contradicting.append(
                {
                    "signal_id": sid,
                    "entity": sig.entity,
                    "stated": stated,
                    "audited": audited,
                }
            )
    return unknown, misattributed, contradicting, not claims


# Methodological/assay tokens that look gene-like but are NOT biological
# entities (perturbation systems, readout assays, reporters, stains, reagents).
# Same rationale as W-CLAIM's ``_GENE_STOPWORDS``: fixed methodological
# vocabulary, not biological content (ADR-011 does not apply). Compared in an
# alphanumeric-normalized, upper-cased form, so ``RT-PCR`` / ``RT-qPCR`` /
# ``RNA-seq`` all resolve to a listed token. Two-letter tokens (KO/KD/WB/IF/IP)
# are already dropped by the length>=3 filter. This list is intentionally
# extensible — readout grounding (H11) is best-effort; an exotic assay acronym
# not listed here will surface as a visible (tunable) rejection, never a silent
# pass of a smuggled gene.
_METHOD_STOPWORDS = {
    # perturbation systems
    "CRISPR", "CRISPRI", "CRISPRA", "CAS9", "CAS12", "CAS13",
    "SHRNA", "SIRNA", "SGRNA", "GRNA",
    # PCR / sequencing assays
    "PCR", "QPCR", "RTPCR", "RTQPCR", "DDPCR", "RNASEQ", "SCRNASEQ",
    "ATACSEQ", "CHIPSEQ", "CUTRUN", "CUTTAG", "NGS", "WES", "WGS", "UMI", "CHIP",
    # protein / blot / interaction assays
    "ELISA", "WESTERN", "NORTHERN", "SOUTHERN", "IHC", "ICC", "COIP", "FRET",
    # microscopy / imaging / flow cytometry
    "FACS", "MFI", "FSC", "SSC", "FMO", "TEM", "SEM", "DAPI", "FISH",
    "SMFISH", "TUNEL",
    # reporters / tags
    "GFP", "EGFP", "RFP", "YFP", "CFP", "BFP", "LUC", "LACZ",
    # common reagents / buffers
    "BSA", "PBS", "FBS", "DMEM", "EDTA", "DMSO",
}


def _stop_key(token: str) -> str:
    """Alphanumeric-normalized, upper-cased key for method-stopword matching."""
    return re.sub(r"[^A-Za-z0-9]", "", token).upper()


def _prose_entity_tokens(text: str) -> set[str]:
    """Gene-like tokens from prose, cleaned of edge punctuation and method noise.

    Reuses W-CLAIM's ``_claim_entities`` (which already drops ``_GENE_STOPWORDS``),
    then strips edge punctuation (so ``RT-qPCR`` does not leak a ``RT-`` fragment)
    and drops short fragments and methodological acronyms (compared in an
    alphanumeric-normalized form so hyphenated assays like ``RT-PCR`` resolve).
    """
    tokens: set[str] = set()
    for raw in _claim_entities(text):
        cleaned = raw.strip("-_.")
        if len(cleaned) < 3 or _stop_key(cleaned) in _METHOD_STOPWORDS:
            continue
        tokens.add(cleaned)
    return tokens


# EVERY generated free-text field the reader sees and the model is free to
# author — including the experiment ``readout`` (H11, reversing the H1 carve-out).
# An entity smuggled into the readout ("TP53 protein abundance") is just as much
# an invented fact as one in the mechanism; the assay vocabulary that made the
# readout noisy (RT-qPCR / FACS / ELISA / GFP / DAPI ...) is now handled by the
# alphanumeric-normalized ``_METHOD_STOPWORDS`` + length filter, so honest assay
# descriptions pass while a smuggled gene is caught.
def _generated_prose(hypothesis: Hypothesis) -> str:
    exp = getattr(hypothesis, "experiment", None)
    da = getattr(hypothesis, "devils_advocate", None) or {}
    parts = [str(getattr(hypothesis, "mechanism", "") or "")]
    if exp is not None:
        parts += [
            str(getattr(exp, "perturbation", "") or ""),
            str(getattr(exp, "readout", "") or ""),
            str(getattr(exp, "predicted_direction", "") or ""),
            str(getattr(exp, "refuting_outcome", "") or ""),
        ]
    parts.append(str(da.get("simpler_explanation", "") or ""))
    return " . ".join(p for p in parts if p)


def verify_hypothesis_grounding(
    hypothesis: Hypothesis,
    signals: list[EvidenceSignal],
    run_ledger: dict | None = None,
) -> GroundingResult:
    """Reject a hypothesis that invents facts or is anchored to nothing.

    Four mechanical checks:

    1. Every entity in ``hypothesis.entities`` must exist in the audited
       evidence index (built from ``signals``). An entity not measured by any
       audited signal is invented — the hypothesis is rejected.
    2. Every gene-like entity named in ANY generated free-text field — the
       mechanism, the discriminating experiment, AND the devils-advocate
       "simpler explanation" — must also resolve to the audited evidence. The
       structured ``entities`` field is not the only surface the reader sees; an
       entity smuggled only into the prose would otherwise evade check (1).
       Reuses W-CLAIM's ``_claim_entities`` extractor (which already drops
       non-gene acronyms via ``_GENE_STOPWORDS``); an undeclared prose entity
       absent from evidence is rejected.
    3. Non-vacuity: the hypothesis must name at least one GROUNDED entity AND
       cite at least one observation. A hypothesis anchored to nothing (no
       entities, no refs) is the ultimate ungrounded case — it is rejected, not
       vacuously accepted (H1, bug 1).
    4. Evidence↔citation lineage (H10): every cited ``observation_ref`` must be
       the ``audited_node_ref`` of at least one of the hypothesis's named
       entities. Existing + run is not enough — a DE gene cited to the pathway
       node is a misattribution (the cited analysis did not produce that entity),
       and is rejected.
    5. Every ``hypothesis.observation_refs`` ledger node must exist AND have run
       (status not in not-run/skipped/error). Reusing W-LEDGER, a hypothesis
       cannot arise from an analysis that did not actually produce results.

    The check on (5) only runs when a ``run_ledger`` is supplied; entity
    grounding (1), prose grounding (2), non-vacuity (3) and lineage (4) are
    always enforced.
    """
    evidence = build_evidence_index(signals)
    declared = {_norm(ent) for ent in (hypothesis.entities or [])}
    missing = [
        ent
        for ent in (hypothesis.entities or [])
        if _norm(ent) not in evidence
    ]

    # (2) Prose grounding over EVERY generated field. Gene-like tokens that are
    # neither in the evidence universe nor among the declared entities are the
    # evasion path the structured-only check (1) misses. Declared-but-missing
    # entities are already reported by (1), so exclude them here to keep the
    # signals distinct.
    ungrounded_prose = sorted(
        {
            token
            for token in _prose_entity_tokens(_generated_prose(hypothesis))
            if _norm(token) not in evidence and _norm(token) not in declared
        }
    )

    # (3) Non-vacuity: at least one grounded entity AND at least one cited
    # observation. Closes the "grounded by naming nothing" acceptance.
    grounded_entities = [
        ent for ent in (hypothesis.entities or []) if _norm(ent) in evidence
    ]
    has_observation = bool(hypothesis.observation_refs)
    vacuous = (not grounded_entities) or (not has_observation)

    # (4) Evidence↔citation lineage: the nodes that actually produced the
    # hypothesis's named entities. A cited observation_ref that is NOT one of
    # them is a misattribution — the analysis it points to did not produce any
    # entity the hypothesis names.
    entity_nodes = {
        sig.audited_node_ref
        for sig in (signals or [])
        if isinstance(sig, EvidenceSignal)
        and _norm(sig.entity) in declared
        and sig.audited_node_ref
    }
    misattributed = [
        ref
        for ref in (hypothesis.observation_refs or [])
        if ref not in entity_nodes
    ]

    # (6) Signal-level grounding (H15): the hypothesis must cite the exact audited
    # signals it reads (observed_claims), and a stated reading must not contradict
    # the audited direction/entity of the cited signal.
    (
        unknown_signals,
        misattributed_signals,
        contradicting_claims,
        missing_observed,
    ) = _verify_observed_claims(hypothesis, signals, declared)

    not_run: list[dict] = []
    if run_ledger is not None:
        index = _node_index(run_ledger)
        for ref in hypothesis.observation_refs or []:
            node = index.get(ref)
            if node is None:
                not_run.append({"node_id": ref, "status": "no_ledger_node"})
            elif node.get("status") in _NOT_RUN_STATUSES:
                not_run.append(
                    {
                        "node_id": ref,
                        "status": node.get("status"),
                        "reason": node.get("reason"),
                    }
                )

    grounded = (
        not missing
        and not not_run
        and not ungrounded_prose
        and not vacuous
        and not misattributed
        and not unknown_signals
        and not misattributed_signals
        and not contradicting_claims
        and not missing_observed
    )
    reason = None
    if not grounded:
        parts: list[str] = []
        if vacuous:
            parts.append(
                "vacuous: requires at least one grounded entity and one "
                "cited observation"
            )
        if missing_observed:
            parts.append(
                "no observed_claims: must cite at least one audited signal it "
                "reads (signal_id + stated_direction)"
            )
        if contradicting_claims:
            parts.append(
                "stated direction contradicts the audited signal: "
                f"{contradicting_claims}"
            )
        if unknown_signals:
            parts.append(
                f"unknown cited signal_ids: {sorted(set(unknown_signals))}"
            )
        if misattributed_signals:
            parts.append(
                "observed_claims cite a signal whose entity is not named: "
                f"{sorted(set(misattributed_signals))}"
            )
        if misattributed:
            parts.append(
                "misattributed observation refs (did not produce a named "
                f"entity): {misattributed}"
            )
        if missing:
            parts.append(f"ungrounded entities: {sorted(set(missing))}")
        if ungrounded_prose:
            parts.append(f"ungrounded prose entities: {ungrounded_prose}")
        if not_run:
            parts.append(
                "observations not run: "
                f"{[n['node_id'] for n in not_run]}"
            )
        reason = "; ".join(parts)
    return GroundingResult(
        grounded=grounded,
        missing_entities=missing,
        not_run_refs=not_run,
        ungrounded_prose_entities=ungrounded_prose,
        misattributed_refs=misattributed,
        unknown_signals=unknown_signals,
        misattributed_signals=misattributed_signals,
        contradicting_claims=contradicting_claims,
        missing_observed_claims=missing_observed,
        vacuous=vacuous,
        reason=reason,
    )
