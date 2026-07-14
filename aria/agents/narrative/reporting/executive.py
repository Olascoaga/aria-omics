"""Executive-summary governance and evidence-block construction."""

from __future__ import annotations

from aria.agents.narrative.reporting._base import *


class ExecutiveSummaryMixin:
    def _govern_executive_summary(self, executive_summary: str,
                                  grouped_findings: dict,
                                  intent: dict,
                                  agent_results: dict,
                                  narrative_blocks: list) -> tuple[str, str]:
        """Guard free-text executive summaries before they reach HTML."""
        text = str(executive_summary or "")
        fallback = ""
        try:
            fallback = self._fallback_executive_summary(grouped_findings, intent)
        except Exception:
            fallback = ""
        if fallback and text.strip() == str(fallback).strip():
            return text, ""

        named_entities: list[str] = []
        for block in narrative_blocks or []:
            try:
                named_entities.extend(collect_named_entities(block))
            except Exception:
                continue

        violations: list[str] = []
        causal_hit = find_causal_language(text, exclude=named_entities)
        if causal_hit:
            violations.append("unlicensed causal language")

        unsupported_numbers = self._unsupported_executive_summary_numbers(
            text, agent_results
        )
        if unsupported_numbers:
            violations.append("unsupported numeric claims")

        if not violations:
            return text, ""

        if not fallback:
            fallback = (
                "ARIA completed the analysis. See the findings table below for "
                "the governed results and limitations."
            )
        warning = (
            "Executive summary governance: free-text summary failed "
            f"{' and '.join(violations)} checks; deterministic fallback shown."
        )
        return fallback, warning

    def _unsupported_executive_summary_numbers(self, text: str,
                                               agent_results: dict) -> set[str]:
        observed = _executive_summary_numbers(text)
        if not observed:
            return set()
        concrete = ""
        try:
            concrete = self._summarize_agent_results_for_llm(agent_results)
        except Exception:
            concrete = ""
        allowed = _executive_summary_numbers(concrete)
        return observed - allowed

    def _build_executive_summary_block(self, executive_summary: str,
                                       executive_summary_warning: str,
                                       grouped_findings: dict,
                                       intent: dict,
                                       exp_ctx: dict,
                                       agent_results: dict,
                                       narrative_blocks: list
                                       ) -> tuple[str, str, NarrativeBlock]:
        """Represent the executive summary as a governed narrative block."""
        text = str(executive_summary or "")
        warning = str(executive_summary_warning or "")
        block = self._make_executive_summary_block(
            text, warning, grouped_findings, intent, exp_ctx, agent_results,
            narrative_blocks,
        )
        try:
            from aria.agents.narrative.evidence_verifier import (
                verify_block_claim_support,
            )
            block.metadata["claim_verification"] = verify_block_claim_support(
                block, strict=True
            )
            return text, warning, block
        except Exception as exc:
            log.warning(
                "Executive summary W-CLAIM verification failed; using "
                "deterministic fallback: %s",
                exc,
            )
            try:
                fallback = self._fallback_executive_summary(
                    grouped_findings, intent
                )
            except Exception:
                fallback = (
                    "ARIA completed the analysis. See the governed findings "
                    "and limitations below."
                )
            if warning:
                warning = (
                    f"{warning} W-CLAIM verification also failed; "
                    "deterministic fallback shown."
                )
            else:
                warning = (
                    "Executive summary governance: W-CLAIM verification failed; "
                    "deterministic fallback shown."
                )
            block = self._make_executive_summary_block(
                fallback, warning, grouped_findings, intent, exp_ctx,
                agent_results, narrative_blocks,
            )
            try:
                from aria.agents.narrative.evidence_verifier import (
                    verify_block_claim_support,
                )
                block.metadata["claim_verification"] = verify_block_claim_support(
                    block, strict=True
                )
            except Exception as fallback_exc:
                log.warning(
                    "Fallback executive summary verification failed: %s",
                    fallback_exc,
                )
            return fallback, warning, block

    def _make_executive_summary_block(self, text: str, warning: str,
                                      grouped_findings: dict,
                                      intent: dict,
                                      exp_ctx: dict,
                                      agent_results: dict,
                                      narrative_blocks: list) -> NarrativeBlock:
        concrete = ""
        try:
            concrete = self._summarize_agent_results_for_llm(agent_results)
        except Exception:
            concrete = "(concrete result summary unavailable)"
        total_findings = sum(len(v) for v in (grouped_findings or {}).values())
        evidence = [
            EvidenceItem(
                label="concrete pipeline results",
                value=concrete,
                source="agent_results",
            ),
            EvidenceItem(
                label="findings recorded",
                value=total_findings,
                source="message_bus",
            ),
            EvidenceItem(
                label="high-confidence findings",
                value=len((grouped_findings or {}).get("high", [])),
                source="message_bus",
            ),
            EvidenceItem(
                label="narrative blocks summarized",
                value=len(narrative_blocks or []),
                source="narrative_registry",
            ),
            EvidenceItem(
                label="allowed confidence labels",
                value="HIGH MEDIUM LOW INSUFFICIENT",
                source="report_policy",
            ),
        ]
        caveats = [Caveat(warning, "warning")] if warning else []
        return NarrativeBlock(
            id="executive_summary",
            modality="report",
            analysis="executive_summary",
            block_type="summary",
            title="Executive Summary",
            status="success",
            confidence="medium",
            claim=text,
            evidence=evidence,
            caveats=caveats,
            metrics={
                "n_findings": total_findings,
                "n_high_confidence": len(
                    (grouped_findings or {}).get("high", [])
                ),
                "n_narrative_blocks_summarized": len(narrative_blocks or []),
            },
            metadata={"render_surface": "executive_summary"},
        )

