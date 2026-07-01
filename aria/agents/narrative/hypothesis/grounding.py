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


# ── entity resolution (round-4 grounding calibration) ────────────────────────
# The audited evidence uses OFFICIAL gene symbols (NR1D1, ARNTL); a model reasons
# in common names and formatting variants (REV-ERBα, REV-ERBa, BMAL1). Literal
# string grounding then rejects a common name for a gene that IS measured. These
# helpers canonicalise names for MATCHING only (never for display), and resolve
# the declared perturbation targets from the run's own contrast labels so a
# hypothesis about the knocked-out gene is grounded in the DESIGN, not "invented".
# No gene names are hardcoded: the design targets come from the data.

_GREEK_TO_LATIN = str.maketrans({
    "α": "a", "β": "b", "γ": "g", "δ": "d", "ε": "e", "κ": "k", "λ": "l",
    "μ": "u", "σ": "s", "ω": "w",
})

# Condition markers that identify an explicitly PERTURBED gene in a contrast label
# (``bmal1_ko`` -> target ``bmal1``). Baselines (wt/ctrl/control) are NOT targets.
_PERTURBATION_MARKERS = (
    "ko", "kd", "oe", "ki", "mut", "knockout", "knockdown", "overexpression",
)
# All condition suffixes stripped from a prose token before grounding
# (``BMAL1_KO`` -> ``BMAL1``); includes the baselines so a labelled prose token
# resolves to its gene.
_CONDITION_SUFFIXES = _PERTURBATION_MARKERS + ("wt", "ctrl", "control")

# 3-letter amino-acid codes: a residue+position token (Ser15, Thr183) is a
# phospho-site, not a gene — it must not read as an ungrounded entity.
_AMINO_ACIDS = {
    "SER", "THR", "TYR", "CYS", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN",
    "GLN", "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TRP", "MET",
}


def _canon(name: str) -> str:
    """Canonical match key: lowercase, Greek→latin, alphanumeric only.

    ``REV-ERBα`` / ``REV-ERBa`` / ``rev_erba`` all collapse to ``reverba`` so a
    model's common name matches a design target or a measured symbol. For MATCHING
    only — the original string is always what is displayed.
    """
    lowered = str(name or "").strip().lower().translate(_GREEK_TO_LATIN)
    return re.sub(r"[^a-z0-9]", "", lowered)


def _strip_condition_suffix(token: str) -> str:
    """``BMAL1_KO`` / ``BMAL1-KO`` -> ``BMAL1`` (drop a trailing condition marker)."""
    m = re.match(
        r"^(.*?)[_\- ]?(" + "|".join(_CONDITION_SUFFIXES) + r")$",
        str(token or "").strip(),
        re.IGNORECASE,
    )
    return m.group(1) if (m and m.group(1)) else str(token or "").strip()


def _design_targets(signals: list[EvidenceSignal]) -> set[str]:
    """Canonical perturbation targets declared by the run's own contrast labels.

    Each signal's ``context`` is its contrast label (``bmal1_ko_vs_wt``). The gene
    a group explicitly perturbs (``bmal1_ko`` -> ``bmal1``) is part of the study's
    DESIGN — the most grounded fact of the experiment — so a hypothesis naming it
    is grounded even when the gene is not itself a differentially expressed row.
    Only explicitly perturbed groups (ko/kd/oe/...) yield a target; baselines do
    not. Data-driven; no gene name is hardcoded.
    """
    targets: set[str] = set()
    for sig in signals or []:
        ctx = str(getattr(sig, "context", "") or "")
        for group in re.split(r"_vs_|\bvs\b|_versus_|\bversus\b", ctx, flags=re.IGNORECASE):
            group = group.strip(" _-")
            m = re.match(
                r"^(.+?)[_\- ]?(" + "|".join(_PERTURBATION_MARKERS) + r")$",
                group, re.IGNORECASE,
            )
            if m and m.group(1):
                canon = _canon(m.group(1))
                if len(canon) >= 2:
                    targets.add(canon)
    return targets


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
    - ``misattributed_signals``: retained for shape compatibility but no longer
      populated (round-4). The old rule rejected a cited signal whose gene was not
      in ``entities``; because the evidence uses official symbols while a model
      cites in common names (NR1D1 vs REV-ERBα), that literal check rejected a
      FAITHFUL citation of a real measurement. Invention is still blocked: the
      signal must EXIST (``unknown``), its direction must not be contradicted
      (``contradicting``), and the hypothesis's own entities are grounded
      separately. Citing a real audited signal is evidence, not fabrication.
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
    misattributed: list[str] = []  # round-4: intentionally never populated
    contradicting: list[dict] = []
    claims = hypothesis.observed_claims or []
    for claim in claims:
        sid = str((claim or {}).get("signal_id", ""))
        sig = by_id.get(sid)
        if sig is None:
            unknown.append(sid)
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
    # bare assay/modality acronyms + mass-spec / genomics + biological concepts
    # that read gene-like but are not entities (round-4 calibration)
    "ATAC", "RNA", "DNA", "CHROMATIN", "ICPMS", "ICP", "MS", "LCMS", "GCMS",
    "MALDI", "SPR", "ITC", "NMR", "EMSA", "SELEX", "GSEA", "ORA", "GO", "KEGG",
    "ECM", "SASP", "ROS", "TCA", "OXPHOS", "EMT",
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
        # Round-4: drop a trailing condition marker (BMAL1_KO -> BMAL1) so a
        # labelled perturbation resolves to its gene, and drop residue+position
        # phospho-site tokens (Ser15) which are not entities.
        cleaned = _strip_condition_suffix(cleaned)
        if len(cleaned) < 3 or _stop_key(cleaned) in _METHOD_STOPWORDS:
            continue
        if _is_residue_token(cleaned):
            continue
        tokens.add(cleaned)
    return tokens


def _is_residue_token(token: str) -> bool:
    """True for a residue+position phospho-site token (Ser15, Thr183, S15)."""
    m = re.match(r"^([A-Za-z]{1,3})(\d+)$", str(token or "").strip())
    if not m:
        return False
    aa = m.group(1).upper()
    return aa in _AMINO_ACIDS or len(aa) == 1


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
    # Round-4 grounding universe: audited entities PLUS the run's declared
    # perturbation targets, matched by canonical key so a common name / formatting
    # variant (REV-ERBα) resolves to a measured symbol (NR1D1) or a design target
    # (rev_erba_ko). ``_grounded`` is the single membership test used everywhere.
    evidence_canon = {
        _canon(sig.entity)
        for sig in (signals or [])
        if isinstance(sig, EvidenceSignal) and _norm(sig.entity)
    }
    design_targets = _design_targets(signals)
    universe = evidence_canon | design_targets

    def _grounded(name: str) -> bool:
        return _canon(name) in universe

    declared_canon = {_canon(ent) for ent in (hypothesis.entities or [])}
    missing = [
        ent for ent in (hypothesis.entities or []) if not _grounded(ent)
    ]

    # (2) Prose grounding over EVERY generated field. A gene-like token neither in
    # the universe nor among the declared entities is the evasion path check (1)
    # misses. A token that is a FRAGMENT of a design target (``REV`` split from
    # ``REV-ERBα``) is not invention — it is tolerated. Declared-but-missing
    # entities are already reported by (1), so exclude them here.
    def _fragment_of_target(canon_tok: str) -> bool:
        return len(canon_tok) >= 3 and any(
            canon_tok in tgt or tgt in canon_tok for tgt in design_targets
        )

    ungrounded_prose = sorted(
        {
            token
            for token in _prose_entity_tokens(_generated_prose(hypothesis))
            if not _grounded(token)
            and _canon(token) not in declared_canon
            and not _fragment_of_target(_canon(token))
        }
    )

    # (3) Non-vacuity: at least one grounded entity AND at least one cited
    # observation. Closes the "grounded by naming nothing" acceptance.
    grounded_entities = [
        ent for ent in (hypothesis.entities or []) if _grounded(ent)
    ]
    has_observation = bool(hypothesis.observation_refs)
    vacuous = (not grounded_entities) or (not has_observation)

    # (4) Evidence↔citation lineage (round-4 softening of H10): a cited
    # observation_ref is valid when it is the node of ANY audited signal in this
    # run — citing the DE node for genes AND the enriched-pathway node for the
    # pathway context is legitimate provenance, not misattribution. What stays
    # rejected is a ref that points at NO audited node (a fabricated citation).
    # Node existence-and-ran is still enforced by check (5) below.
    audited_nodes = {
        sig.audited_node_ref
        for sig in (signals or [])
        if isinstance(sig, EvidenceSignal) and sig.audited_node_ref
    }
    misattributed = [
        ref
        for ref in (hypothesis.observation_refs or [])
        if ref not in audited_nodes
    ]

    # (6) Signal-level grounding (H15): the hypothesis must cite the exact audited
    # signals it reads (observed_claims), and a stated reading must not contradict
    # the audited direction of the cited signal.
    (
        unknown_signals,
        misattributed_signals,
        contradicting_claims,
        missing_observed,
    ) = _verify_observed_claims(hypothesis, signals, declared_canon)

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
