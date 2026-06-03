"""Chromatin (scATAC / bulk-ATAC / ChIP) narrative plugin — SKELETON (P3-1).

A from-day-one ``ChromatinNarrator`` so that, the moment the v4.6 chromatin stack
lands, its report sections compose from validated ``NarrativeBlock`` objects in
this package (exactly like scRNA / bulk RNA) instead of growing
``narrative_agent.py`` — and so a thin chromatin report is reconciled by the run
ledger as a divergence rather than vanishing silently.

This is deliberately a SKELETON: it asserts no chromatin biology and fabricates
nothing. It surfaces only the structured QC metrics the scaffold
``chromatin_qc.py`` actually measures, marks the metrics that need the v4.6 stack
(FRiP, TSS enrichment) as not computed, and states the scaffold /
not-publication-grade limitation. Richer chromatin narration (peaks, differential
accessibility, motifs) is intentionally deferred to v4.6, where it is added HERE,
never back in ``narrative_agent.py``. The narrator only fires on
``chromatin_agent`` results, so it is a no-op on the validated RNA paths.
"""

from __future__ import annotations

import re

from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock


def _safe_id(value) -> str:
    text = str(value or "chromatin")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "chromatin"


class ChromatinNarrator:
    name = "chromatin"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        return (
            agent_name == "chromatin_agent"
            and isinstance(agent_result, dict)
            and bool(agent_result.get("findings", agent_result))
        )

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        findings = agent_result.get("findings", {}) or {}
        blocks: list[NarrativeBlock] = []
        qc = findings.get("qc")
        if isinstance(qc, dict) and str(qc.get("status")) in ("success", "done"):
            blocks.append(self._qc_block(qc))
        return blocks

    def _qc_block(self, qc: dict) -> NarrativeBlock:
        # Only measured metrics become evidence; None-valued (uncomputed) metrics
        # are never turned into evidence — that is the no-fabrication contract.
        evidence: list[EvidenceItem] = []
        for label, key in (
            ("Unique barcodes", "n_cells"),
            ("Fragments scanned", "n_fragments"),
            ("Mitochondrial fraction", "mito_fraction"),
        ):
            val = qc.get(key)
            if val is not None:
                evidence.append(
                    EvidenceItem(label=label, value=val, source="chromatin_qc")
                )

        caveats = [Caveat(
            "scATAC/chromatin analysis is a v4.6 scaffold: ARIA reports only the "
            "QC metrics actually measured and does not yet produce "
            "publication-grade chromatin results.",
            severity="warning",
        )]
        not_computed = qc.get("metrics_not_computed") or []
        if not_computed:
            caveats.append(Caveat(
                "Not computed (needs the v4.6 stack / called peaks): "
                + ", ".join(str(m) for m in not_computed),
                severity="info",
            ))

        status = "success" if evidence else "limitation"
        claim = (
            "Chromatin QC reported the measured fragment and barcode metrics below."
            if evidence else ""
        )
        return NarrativeBlock(
            id=f"chromatin.qc.{_safe_id(qc.get('data_type'))}",
            modality="chromatin",
            analysis="qc",
            block_type="data_quality",
            title="Chromatin QC (scaffold)",
            status=status,
            confidence="low",
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            metrics={
                "qc_complete": qc.get("qc_complete"),
                "pass_qc": qc.get("pass_qc"),
            },
            metadata={"validation_level": "scaffold"},
        )

    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        findings = agent_result.get("findings", {}) or {}
        if not findings.get("qc"):
            return []
        return [
            "Chromatin QC (scaffold): only measured fragment/barcode metrics are "
            "reported; TSS enrichment and FRiP remain uncomputed until the v4.6 "
            "chromatin stack and called peaks are available."
        ]
