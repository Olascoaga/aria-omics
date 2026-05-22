"""Serializable narrative block schema used by modality narrators."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal


Confidence = Literal["high", "medium", "low", "insufficient"]
CaveatSeverity = Literal["info", "warning", "blocking"]


@dataclass
class EvidenceItem:
    label: str
    value: str | int | float | None = None
    source: str | None = None
    path: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "EvidenceItem":
        return cls(
            label=str(data.get("label", "")),
            value=data.get("value"),
            source=data.get("source"),
            path=data.get("path"),
        )


@dataclass
class Caveat:
    text: str
    severity: CaveatSeverity = "warning"

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Caveat":
        severity = data.get("severity", "warning")
        if severity not in {"info", "warning", "blocking"}:
            severity = "warning"
        return cls(text=str(data.get("text", "")), severity=severity)


@dataclass
class NarrativeBlock:
    id: str
    modality: str
    analysis: str
    block_type: str
    title: str
    status: str
    confidence: Confidence
    claim: str
    evidence: list[EvidenceItem] = field(default_factory=list)
    caveats: list[Caveat] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    figures: list[dict] = field(default_factory=list)
    tables: list[dict] = field(default_factory=list)
    methods: list[str] = field(default_factory=list)
    error: str | None = None
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.status == "success":
            if not str(self.claim or "").strip():
                raise ValueError("successful NarrativeBlock requires claim")
            if not self.evidence:
                raise ValueError("successful NarrativeBlock requires evidence")

    def to_dict(self) -> dict:
        data = asdict(self)
        data["evidence"] = [item.to_dict() for item in self.evidence]
        data["caveats"] = [item.to_dict() for item in self.caveats]
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "NarrativeBlock":
        return cls(
            id=str(data.get("id", "")),
            modality=str(data.get("modality", "")),
            analysis=str(data.get("analysis", "")),
            block_type=str(data.get("block_type", "")),
            title=str(data.get("title", "")),
            status=str(data.get("status", "")),
            confidence=data.get("confidence", "insufficient"),
            claim=str(data.get("claim", "")),
            evidence=[
                EvidenceItem.from_dict(item)
                for item in data.get("evidence", []) or []
            ],
            caveats=[
                Caveat.from_dict(item)
                for item in data.get("caveats", []) or []
            ],
            metrics=dict(data.get("metrics", {}) or {}),
            figures=list(data.get("figures", []) or []),
            tables=list(data.get("tables", []) or []),
            methods=list(data.get("methods", []) or []),
            error=data.get("error"),
            warnings=list(data.get("warnings", []) or []),
            metadata=dict(data.get("metadata", {}) or {}),
        )
