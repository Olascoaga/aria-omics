"""Schemas for the SPECULATIVE hypothesis layer (ADR-057).

Pure data: these dataclasses assert no biology and run no LLM. ``EvidenceSignal``
is the modality-agnostic Evidence Interface — per-modality adapters (S5-S8)
flatten audited artifacts into a uniform, already-grounded vocabulary the agent
reasons over. ``Hypothesis`` is the agent's structured output; the grounding
verifier (grounding.py) enforces that every entity it names is real.

Note: ``EvidenceSignal`` is deliberately NOT the narrative ``EvidenceItem`` in
``aria/agents/narrative/types.py`` (that one is a report-block evidence row). The
distinct name avoids the collision.
"""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass, field
from typing import Literal

EntityKind = Literal[
    "gene", "peak", "tf_motif", "pathway", "cluster", "celltype", "other"
]
Modality = Literal["bulk_RNA", "scRNA", "bulk_ATAC", "scATAC", "integration"]
Direction = Literal["up", "down", "na"]

_ENTITY_KINDS = {
    "gene", "peak", "tf_motif", "pathway", "cluster", "celltype", "other",
}
_DIRECTIONS = {"up", "down", "na"}


@dataclass
class EvidenceSignal:
    """One grounded, audited measurement the agent may reason over.

    ``entity`` must already be a real entity from the audited results and
    ``audited_node_ref`` a run-ledger node id (``ledger://<modality>/<analysis>``)
    that actually ran. The grounding verifier rejects anything else, so the LLM
    is free over the *connection* while these *facts* stay real.
    """

    entity: str
    entity_kind: EntityKind
    modality: Modality
    measure: str  # LFC | padj | motif_enrich | footprint_score | composition_delta | ...
    audited_node_ref: str  # ledger://<modality>/<analysis> that produced it
    value: float | int | str | None = None
    direction: Direction = "na"
    caveats_inherited: list[str] = field(default_factory=list)
    # The analysis context the signal was measured in (a contrast / comparison /
    # cluster / group label). The SAME entity measured in two contexts is two
    # distinct, independent signals — the adapters must NOT collapse them, or a
    # contrast is lost and the ranking under-counts the converging lines (H4).
    context: str = ""
    # Stable identity. Two signals are the same iff they share modality, node,
    # entity, measure, direction AND context. Derived deterministically when not
    # supplied so the index keys by signal, not by bare entity.
    signal_id: str = ""

    def __post_init__(self) -> None:
        if not self.signal_id:
            self.signal_id = self._derive_signal_id()

    def _derive_signal_id(self) -> str:
        basis = "|".join(
            str(part)
            for part in (
                self.modality,
                self.audited_node_ref,
                str(self.entity).strip().lower(),
                self.measure,
                self.direction,
                self.context,
            )
        )
        return "sig_" + hashlib.sha1(basis.encode("utf-8")).hexdigest()[:12]

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "EvidenceSignal":
        kind = str(data.get("entity_kind", "other"))
        if kind not in _ENTITY_KINDS:
            kind = "other"
        direction = str(data.get("direction", "na"))
        if direction not in _DIRECTIONS:
            direction = "na"
        return cls(
            entity=str(data.get("entity", "")),
            entity_kind=kind,  # type: ignore[arg-type]
            modality=str(data.get("modality", "")),  # type: ignore[arg-type]
            measure=str(data.get("measure", "")),
            audited_node_ref=str(data.get("audited_node_ref", "")),
            value=data.get("value"),
            direction=direction,  # type: ignore[arg-type]
            caveats_inherited=list(data.get("caveats_inherited", []) or []),
            context=str(data.get("context", "") or ""),
            signal_id=str(data.get("signal_id", "") or ""),
        )


@dataclass
class DiscriminatingExperiment:
    """Falsifiability schema (ADR-057 rail #8).

    A concrete experiment/analysis that would confirm or refute the hypothesis.
    The S2 gate enforces completeness; in S1 the schema simply exists so the
    output shape is stable.
    """

    perturbation: str = ""
    readout: str = ""
    predicted_direction: str = ""
    refuting_outcome: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "DiscriminatingExperiment":
        data = data or {}
        return cls(
            perturbation=str(data.get("perturbation", "")),
            readout=str(data.get("readout", "")),
            predicted_direction=str(data.get("predicted_direction", "")),
            refuting_outcome=str(data.get("refuting_outcome", "")),
        )

    def is_complete(self) -> bool:
        """True iff every falsifiability field is populated (used by the S2 gate)."""
        return all(
            str(v).strip()
            for v in (
                self.perturbation,
                self.readout,
                self.predicted_direction,
                self.refuting_outcome,
            )
        )


@dataclass
class Hypothesis:
    """A machine-generated, SPECULATIVE hypothesis.

    ``entities`` are the audited entities the hypothesis names (the grounding
    target); ``observation_refs`` are the run-ledger node ids the hypothesis
    arises from (a). ``mechanism`` is the proposed connection (b), and
    ``experiment`` is the discriminating experiment (c). Fields populated by
    later slices (competing set S3-ranking, devils_advocate S3, provenance S9,
    ``ledger_node`` quarantine S4) default empty.
    """

    id: str
    mechanism: str
    entities: list[str] = field(default_factory=list)
    observation_refs: list[str] = field(default_factory=list)
    experiment: DiscriminatingExperiment = field(
        default_factory=DiscriminatingExperiment
    )
    competing_with: list[str] = field(default_factory=list)
    devils_advocate: dict = field(default_factory=dict)
    rank_evidence: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    tier: str = "SPECULATIVE"
    ledger_node: str | None = None

    def to_dict(self) -> dict:
        data = asdict(self)
        data["experiment"] = self.experiment.to_dict()
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "Hypothesis":
        return cls(
            id=str(data.get("id", "")),
            mechanism=str(data.get("mechanism", "")),
            entities=list(data.get("entities", []) or []),
            observation_refs=list(data.get("observation_refs", []) or []),
            experiment=DiscriminatingExperiment.from_dict(
                data.get("experiment", {}) or {}
            ),
            competing_with=list(data.get("competing_with", []) or []),
            devils_advocate=dict(data.get("devils_advocate", {}) or {}),
            rank_evidence=dict(data.get("rank_evidence", {}) or {}),
            provenance=dict(data.get("provenance", {}) or {}),
            tier=str(data.get("tier", "SPECULATIVE")),
            ledger_node=data.get("ledger_node"),
        )
