"""Paired multiome (RNA + ATAC) narrator (V48 Slice I7).

Translates ``IntegrationAgent`` structured outputs (I1-I6) into governed
``NarrativeBlock`` instances under the ``multiome.*`` block namespace with
``modality="integration"``:

  - ``multiome.gate`` (diagnostic when same-cell gate blocked execution)
  - ``multiome.paired_qc`` (``analysis="paired_qc"``)
  - ``multiome.wnn`` (``analysis="wnn"``; surfaces 5-axis guardrail)
  - ``multiome.mofa`` (``analysis="mofa"``; surfaces excluded confounders)
  - ``multiome.peak_to_gene`` (``analysis="peak_to_gene"``; associative cis)
  - ``multiome.eregulons`` (``analysis="eregulons"``; candidate eRegulons)
  - ``multiome.triangulation`` (``analysis="cross_modal_triangulation"``)

Strict invariants (ADR-064, ADR-070, C1-C4):
  - Zero causal verbs in claims, evidence items, caveats, or methods.
  - Non-executed / skipped / failed stages emit honest diagnostic blocks with
    ``claim=""`` so failures remain visible without fabricating public claims.
  - ``methods()`` emits copy-pasteable manuscript paragraphs strictly for stages
    that executed, with markup-stripped dynamic tokens safe for plain-text export
    and single-pass HTML escaping.
"""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

from aria.agents.narrative.synthesis.cross_modal_triangulation import (
    build_cross_modal_evidence_graph,
    build_triangulation_block,
)
from aria.agents.narrative.types import (
    Caveat,
    EvidenceItem,
    NarrativeBlock,
    SemanticFact,
)

_OK_STATUSES = {"success", "ok", "done"}
_NON_RUN_STATUSES = {"skipped", "not_run", "blocked", "error", "failed"}
_TAG_STRIP_RE = re.compile(r"<[^>]+>")


def unwrap_integration_findings(agent_result: dict | None) -> dict[str, Any]:
    """Return the flat integration findings dict from an ``IntegrationAgent`` envelope."""
    if not isinstance(agent_result, dict):
        return {}
    findings = agent_result.get("findings")
    if isinstance(findings, dict):
        return findings
    return agent_result


def _clean_token(val: Any, default: str = "unspecified") -> str:
    """Sanitize dynamic identifiers for copy-pasteable Methods (no HTML tags/control chars)."""
    if val is None:
        return default
    txt = _TAG_STRIP_RE.sub("", str(val))
    txt = re.sub(r"[\r\n\t]+", " ", txt).strip()
    return txt or default


def _format_motif_source(val: Any, default: str = "JASPAR2024_CORE_vertebrates") -> str:
    if isinstance(val, dict):
        cand = val.get("collection") or val.get("source") or val.get("release")
        return _clean_token(cand, default)
    return _clean_token(val, default)


def _count_factor_collection(val: Any, *, exclude_confounders: bool = True) -> int:
    items: list[Any] = []
    if isinstance(val, dict):
        for sub in val.values():
            if isinstance(sub, list):
                items.extend(sub)
    elif isinstance(val, list):
        items.extend(val)
    if not exclude_confounders:
        return len(items)
    return sum(
        1
        for item in items
        if not (
            isinstance(item, dict)
            and bool(item.get("excluded_from_biological_interpretation"))
        )
    )


def _finite_round(val: Any, digits: int = 4) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _stage_state(stage: Any) -> str:
    """Classify stage dict into ``'ok'``, ``'skipped'``, ``'error'``, or ``'absent'``."""
    if not isinstance(stage, dict) or not stage:
        return "absent"
    st = str(stage.get("status") or "").lower().strip()
    if st in {"error", "failed"}:
        return "error"
    if st in {"skipped", "not_run", "blocked"} or stage.get("ran") is False:
        return "skipped"
    if st in _OK_STATUSES or stage.get("ran") is True:
        return "ok"
    return "ok"


def _diagnostic_block(
    block_id: str,
    analysis: str,
    title: str,
    stage: dict[str, Any],
    state: str,
    block_type: str = "result",
) -> NarrativeBlock:
    reason = _clean_token(
        stage.get("reason")
        or stage.get("details")
        or stage.get("error_type")
        or f"stage reported status={stage.get('status')}",
        default="stage did not execute",
    )
    status = "error" if state == "error" else "skipped"
    return NarrativeBlock(
        id=block_id,
        modality="integration",
        analysis=analysis,
        block_type=block_type,
        title=title,
        status=status,
        confidence="insufficient",
        claim="",
        error=reason,
        caveats=[
            Caveat(
                text=f"{title} status is {status} ({reason}); no inferred result was published.",
                severity="warning" if status == "error" else "info",
            )
        ],
    )


def _figure_entry(findings: dict[str, Any], key: str, caption: str) -> list[dict[str, str]]:
    figs = findings.get("figures") or {}
    figs_svg = findings.get("figures_svg") or {}
    raw = figs.get(key)
    if not raw or not isinstance(raw, str):
        return []
    rec: dict[str, str] = {"id": f"multiome.{key}", "path": raw, "caption": caption}
    svg_p = figs_svg.get(key)
    if isinstance(svg_p, str) and svg_p:
        rec["svg_path"] = svg_p
    return [rec]


def _table_entry(findings: dict[str, Any], key: str, fallback_path: Any, label: str) -> list[dict[str, str]]:
    staged = (findings.get("staged_tables") or {}).get(key)
    chosen = staged if isinstance(staged, str) and staged else None
    if not chosen and fallback_path:
        p = Path(str(fallback_path))
        if p.exists():
            chosen = str(p)
    if not chosen:
        return []
    return [{"id": f"multiome.{key}", "path": chosen, "label": label}]


class IntegrationNarrator:
    """Modality narrator for paired same-cell multiome integration (I1-I7)."""

    name = "integration"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        if agent_name != "integration_agent" or not isinstance(agent_result, dict):
            return False
        st = str(agent_result.get("status") or "").lower()
        if st in {"done", "success", "ok", "partial", "error", "blocked"}:
            return True
        findings = unwrap_integration_findings(agent_result)
        return any(
            k in findings
            for k in ("multiome_gate", "paired_qc", "wnn", "mofa", "peak_to_gene", "eregulons")
        )

    def collect(
        self,
        agent_name: str,
        agent_result: dict,
        context: dict | None = None,
    ) -> list[NarrativeBlock]:
        findings = unwrap_integration_findings(agent_result)
        blocks: list[NarrativeBlock] = []

        # 0. Same-cell multiome hard gate diagnostic when blocked/failed.
        gate = findings.get("multiome_gate")
        if isinstance(gate, dict) and gate.get("passed") is False:
            reason = _clean_token(
                gate.get("reason") or gate.get("details"),
                default="same-cell barcode pairing gate did not pass",
            )
            blocks.append(
                NarrativeBlock(
                    id="multiome.gate",
                    modality="integration",
                    analysis="paired_qc",
                    block_type="qc",
                    title="Same-cell multiome pairing gate",
                    status="skipped",
                    confidence="insufficient",
                    claim="",
                    error=reason,
                    caveats=[
                        Caveat(
                            text=(
                                "Same-cell multiome integration requires verified "
                                "shared cell barcodes across RNA and ATAC; cross-sample "
                                "pseudo-pairing is disallowed."
                            ),
                            severity="blocking",
                        )
                    ],
                )
            )

        # 1. Paired QC (I1).
        pqc = findings.get("paired_qc")
        pqc_state = _stage_state(pqc)
        if pqc_state in {"skipped", "error"} and isinstance(pqc, dict):
            blocks.append(
                _diagnostic_block(
                    "multiome.paired_qc", "paired_qc",
                    "Paired multiome quality control", pqc, pqc_state,
                    block_type="qc",
                )
            )
        elif pqc_state == "ok" and isinstance(pqc, dict):
            blocks.append(self._paired_qc_block(pqc))

        # 2. Joint WNN + 5-axis anti-degradation guardrail (I2 + I3).
        wnn = findings.get("wnn")
        wnn_state = _stage_state(wnn)
        if wnn_state in {"skipped", "error"} and isinstance(wnn, dict):
            blocks.append(
                _diagnostic_block(
                    "multiome.wnn", "wnn",
                    "Weighted Nearest Neighbors (WNN) joint manifold",
                    wnn, wnn_state,
                )
            )
        elif wnn_state == "ok" and isinstance(wnn, dict):
            blocks.append(self._wnn_block(wnn, findings))

        # 3. MOFA+ variance decomposition & confounder quarantine (I4).
        mofa = findings.get("mofa")
        mofa_state = _stage_state(mofa)
        if mofa_state in {"skipped", "error"} and isinstance(mofa, dict):
            blocks.append(
                _diagnostic_block(
                    "multiome.mofa", "mofa",
                    "MOFA+ multi-omics factor decomposition",
                    mofa, mofa_state,
                )
            )
        elif mofa_state == "ok" and isinstance(mofa, dict):
            blocks.append(self._mofa_block(mofa, findings))

        # 4. Signac-faithful cis peak_to_gene associations (I5).
        p2g = findings.get("peak_to_gene")
        p2g_state = _stage_state(p2g)
        ereg = findings.get("eregulons")
        ereg_state = _stage_state(ereg)
        if p2g_state in {"skipped", "error"} and isinstance(p2g, dict):
            blocks.append(
                _diagnostic_block(
                    "multiome.peak_to_gene", "peak_to_gene",
                    "Candidate cis peak-to-gene associations",
                    p2g, p2g_state,
                )
            )
        elif p2g_state == "ok" and isinstance(p2g, dict):
            blocks.append(
                self._peak_to_gene_block(
                    p2g, findings, attach_triangulation_fig=(ereg_state != "ok")
                )
            )

        # 5. SCENIC+-style candidate eRegulons (I6).
        if ereg_state in {"skipped", "error"} and isinstance(ereg, dict):
            blocks.append(
                _diagnostic_block(
                    "multiome.eregulons", "eregulons",
                    "Candidate eRegulon regulatory modules",
                    ereg, ereg_state,
                )
            )
        elif ereg_state == "ok" and isinstance(ereg, dict):
            blocks.append(self._eregulons_block(ereg, findings))

        # 6. Deterministic Cross-Modal Evidence Triangulation (I7).
        exp_ctx = (context or {}).get("exp_context") or {}
        all_agent_results = (context or {}).get("agent_results")
        if not isinstance(all_agent_results, dict):
            all_agent_results = {"integration_agent": agent_result}
        elif "integration_agent" not in all_agent_results:
            all_agent_results = {**all_agent_results, "integration_agent": agent_result}

        graph = build_cross_modal_evidence_graph(all_agent_results, exp_ctx)
        tri_block = build_triangulation_block(graph)
        if tri_block is not None:
            blocks.append(tri_block)

        return blocks

    @staticmethod
    def _paired_qc_block(pqc: dict[str, Any]) -> NarrativeBlock:
        n_before = int(pqc.get("n_cells_before") or 0)
        n_after = int(pqc.get("n_cells_after") or 0)
        pct_rem = _finite_round(pqc.get("pct_removed"), 2) or 0.0
        n_genes = int(pqc.get("n_genes") or 0)
        n_peaks = int(pqc.get("n_peaks") or 0)

        evidence: list[EvidenceItem] = [
            EvidenceItem("cells before QC", n_before, source="paired_qc"),
            EvidenceItem("cells after QC", n_after, source="paired_qc"),
            EvidenceItem("removed fraction pct", f"{pct_rem}%", source="paired_qc"),
            EvidenceItem("RNA genes", n_genes, source="paired_qc"),
            EvidenceItem("ATAC peaks", n_peaks, source="paired_qc"),
        ]
        lj = pqc.get("label_join")
        if isinstance(lj, dict) and lj.get("n_matched") is not None:
            evidence.append(
                EvidenceItem(
                    "label join matched cells",
                    int(lj["n_matched"]),
                    source="paired_qc",
                )
            )

        caveats: list[Caveat] = []
        not_computed = [str(m) for m in (pqc.get("metrics_not_computed") or []) if m]
        if not_computed:
            caveats.append(
                Caveat(
                    text=(
                        "Metrics not computed due to absent fragment or annotation "
                        "inputs: " + ", ".join(not_computed) + "."
                    ),
                    severity="info",
                )
            )

        confidence = "high" if n_after >= 500 else ("medium" if n_after >= 100 else "low")
        claim = (
            f"Paired multiome QC retained {n_after} of {n_before} shared barcodes "
            f"({pct_rem}% removed) across {n_genes} RNA genes and {n_peaks} ATAC peaks."
        )
        return NarrativeBlock(
            id="multiome.paired_qc",
            modality="integration",
            analysis="paired_qc",
            block_type="qc",
            title="Paired multiome quality control",
            status="success",
            confidence=confidence,
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            metrics={
                "n_cells_before": n_before,
                "n_cells_after": n_after,
                "pct_removed": pct_rem,
                "n_genes": n_genes,
                "n_peaks": n_peaks,
                "metrics_not_computed": not_computed,
            },
        )

    @staticmethod
    def _wnn_block(wnn: dict[str, Any], findings: dict[str, Any]) -> NarrativeBlock:
        n_cells = int(wnn.get("n_cells") or 0)
        n_joint = int(wnn.get("n_joint_clusters") or 0)
        n_rna_cl = int(wnn.get("n_rna_only_clusters") or 0)
        n_atac_cl = int(wnn.get("n_atac_only_clusters") or 0)
        mean_rw = _finite_round(wnn.get("mean_rna_weight"), 4)
        mean_aw = _finite_round(wnn.get("mean_atac_weight"), 4)
        ari_r = _finite_round(wnn.get("wnn_vs_rna_ari"), 4) or 0.0
        ari_a = _finite_round(wnn.get("wnn_vs_atac_ari"), 4) or 0.0
        nmi_r = _finite_round(wnn.get("wnn_vs_rna_nmi"), 4) or 0.0
        nmi_a = _finite_round(wnn.get("wnn_vs_atac_nmi"), 4) or 0.0
        n_disc = int(wnn.get("n_discordant_clusters") or 0)
        guard = str(wnn.get("guardrail_status") or "skipped")

        evidence: list[EvidenceItem] = [
            EvidenceItem("shared cells", n_cells, source="wnn"),
            EvidenceItem("WNN joint clusters", n_joint, source="wnn"),
            EvidenceItem("RNA-only clusters", n_rna_cl, source="wnn"),
            EvidenceItem("ATAC-only clusters", n_atac_cl, source="wnn"),
            EvidenceItem("WNN vs RNA ARI", ari_r, source="wnn"),
            EvidenceItem("WNN vs ATAC ARI", ari_a, source="wnn"),
            EvidenceItem("WNN vs RNA NMI", nmi_r, source="wnn"),
            EvidenceItem("WNN vs ATAC NMI", nmi_a, source="wnn"),
            EvidenceItem("discordant clusters", n_disc, source="wnn"),
            EvidenceItem("WNN guardrail status", guard, source="wnn_diagnostics"),
        ]
        if mean_rw is not None and mean_aw is not None:
            evidence.append(EvidenceItem("mean RNA weight", mean_rw, source="wnn"))
            evidence.append(EvidenceItem("mean ATAC weight", mean_aw, source="wnn"))

        caveats: list[Caveat] = [
            Caveat(
                text=(
                    "Joint WNN manifold integration is an exploratory multi-modal "
                    "lens rather than a replacement for single-modality RNA and "
                    "ATAC clusterings."
                ),
                severity="info",
            )
        ]
        anti = wnn.get("anti_degradation") or {}
        if isinstance(anti, dict):
            for issue in (anti.get("issues") or [])[:4]:
                if not isinstance(issue, dict):
                    continue
                code = _clean_token(issue.get("code"), "guardrail_notice")
                sev = "warning" if issue.get("severity") == "degraded" else "info"
                caveats.append(
                    Caveat(
                        text=(
                            f"WNN anti-degradation check [{code}] flagged status "
                            f"{issue.get('severity', 'warn')}; inspect per-modality "
                            "embeddings alongside the joint graph."
                        ),
                        severity=sev,
                    )
                )

        if guard == "degraded":
            confidence = "low"
        elif guard in {"warn", "skipped"} or n_cells < 500:
            confidence = "medium"
        else:
            confidence = "high"

        weight_clause = (
            f", mean RNA weight {mean_rw} and mean ATAC weight {mean_aw},"
            if (mean_rw is not None and mean_aw is not None)
            else ""
        )
        claim = (
            f"WNN integration across {n_cells} shared cells identified {n_joint} "
            f"joint clusters (vs {n_rna_cl} RNA-only and {n_atac_cl} ATAC-only){weight_clause} "
            f"with WNN vs RNA ARI {ari_r} (NMI {nmi_r}), WNN vs ATAC ARI {ari_a} "
            f"(NMI {nmi_a}), {n_disc} discordant cluster(s), and guardrail status {guard}."
        )

        figs = (
            _figure_entry(
                findings,
                "wnn_joint_vs_modality_umap",
                "Joint WNN UMAP compared against single-modality RNA and ATAC UMAPs.",
            )
            + _figure_entry(
                findings,
                "wnn_modality_weights_by_cluster",
                "Per-cluster WNN RNA vs ATAC modality weight distribution and dominance thresholds.",
            )
        )
        return NarrativeBlock(
            id="multiome.wnn",
            modality="integration",
            analysis="wnn",
            block_type="result",
            title="Weighted Nearest Neighbors (WNN) joint manifold",
            status="success",
            confidence=confidence,
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            figures=figs,
            metrics={
                "n_cells": n_cells,
                "n_joint_clusters": n_joint,
                "n_rna_only_clusters": n_rna_cl,
                "n_atac_only_clusters": n_atac_cl,
                "mean_rna_weight": mean_rw,
                "mean_atac_weight": mean_aw,
                "wnn_vs_rna_ari": ari_r,
                "wnn_vs_atac_ari": ari_a,
                "wnn_vs_rna_nmi": nmi_r,
                "wnn_vs_atac_nmi": nmi_a,
                "n_discordant_clusters": n_disc,
                "guardrail_status": guard,
            },
        )

    @staticmethod
    def _mofa_block(mofa: dict[str, Any], findings: dict[str, Any]) -> NarrativeBlock:
        n_cells = int(mofa.get("n_cells") or 0)
        n_trained = int(mofa.get("n_cells_trained") or n_cells)
        n_factors = int(mofa.get("n_factors") or 0)
        n_active = int(mofa.get("n_active_factors") or 0)
        n_shared = _count_factor_collection(mofa.get("shared_factors"))
        n_spec = _count_factor_collection(mofa.get("modality_specific_factors"))
        excluded = list(mofa.get("excluded_factors") or [])
        n_excl = len(excluded)

        ve = mofa.get("variance_explained") or {}
        ve_add = mofa.get("variance_explained_additive") or {}
        rna_r2_pct = _finite_round(float(ve.get("rna") or 0.0) * 100.0, 2) or 0.0
        atac_r2_pct = _finite_round(float(ve.get("atac") or 0.0) * 100.0, 2) or 0.0
        rna_add_pct = _finite_round(float(ve_add.get("rna") or 0.0) * 100.0, 2) or 0.0
        atac_add_pct = _finite_round(float(ve_add.get("atac") or 0.0) * 100.0, 2) or 0.0

        evidence: list[EvidenceItem] = [
            EvidenceItem("MOFA cells evaluated", n_cells, source="mofa"),
            EvidenceItem("MOFA cells trained", n_trained, source="mofa"),
            EvidenceItem("fitted factors", n_factors, source="mofa"),
            EvidenceItem("active factors", n_active, source="mofa"),
            EvidenceItem("shared factors", n_shared, source="mofa"),
            EvidenceItem("modality-specific factors", n_spec, source="mofa"),
            EvidenceItem("excluded confounder factors", n_excl, source="mofa"),
            EvidenceItem("RNA non-additive R2 pct", f"{rna_r2_pct}%", source="mofa"),
            EvidenceItem("ATAC non-additive R2 pct", f"{atac_r2_pct}%", source="mofa"),
            EvidenceItem("RNA additive R2 sum pct", f"{rna_add_pct}%", source="mofa"),
            EvidenceItem("ATAC additive R2 sum pct", f"{atac_add_pct}%", source="mofa"),
        ]

        excluded_ids: list[int] = []
        for ex in excluded:
            if not isinstance(ex, dict):
                continue
            raw_fid = (
                ex.get("factor_id")
                if ex.get("factor_id") is not None
                else ex.get("factor")
            )
            if raw_fid is None:
                continue
            fid = int(raw_fid)
            excluded_ids.append(fid)
            evidence.append(
                EvidenceItem(
                    f"excluded confounder Factor {fid}",
                    fid,
                    source="mofa",
                    facts=[
                        SemanticFact(
                            subject=f"Factor_{fid}",
                            subject_type="latent_factor",
                            predicate="confounder_excluded",
                            polarity="affirmed",
                            value=fid,
                            source="mofa",
                        )
                    ],
                )
            )

        caveats: list[Caveat] = [
            Caveat(
                text=(
                    "MOFA+ latent factors summarize linear co-variation across "
                    "log-normalized RNA and TF-IDF ATAC views; factor loadings "
                    "are associative signatures rather than causal programs."
                ),
                severity="info",
            )
        ]
        if n_excl > 0:
            caveats.append(
                Caveat(
                    text=(
                        f"{n_excl} MOFA+ factor(s) correlated with sequencing depth "
                        "or cell-cycle scores were quarantined from biological "
                        "factor rankings."
                    ),
                    severity="warning",
                )
            )
        if mofa.get("subsampled"):
            caveats.append(
                Caveat(
                    text=(
                        f"MOFA+ model fitting used Horvitz-Thompson stratified "
                        f"subsampling ({n_trained} of {n_cells} cells) followed by "
                        "closed-form ridge posterior projection to all cells."
                    ),
                    severity="info",
                )
            )

        confidence = "high" if (n_active > 0 and n_cells >= 500) else "medium"
        claim = (
            f"MOFA+ decomposition across {n_cells} cells ({n_factors} fitted "
            f"factors) retained {n_active} active factors ({n_shared} shared and "
            f"{n_spec} modality-specific), quarantined {n_excl} technical or "
            f"cell-cycle confounder factor(s), and explained non-additive R2 of "
            f"{rna_r2_pct}% in RNA and {atac_r2_pct}% in ATAC."
        )
        return NarrativeBlock(
            id="multiome.mofa",
            modality="integration",
            analysis="mofa",
            block_type="result",
            title="MOFA+ multi-omics factor decomposition",
            status="success",
            confidence=confidence,
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            figures=_figure_entry(
                findings,
                "mofa_variance_explained",
                "MOFA+ per-factor variance explained across RNA and ATAC views with excluded confounders highlighted.",
            ),
            tables=_table_entry(
                findings,
                "mofa_factor_loadings",
                mofa.get("factor_loadings_path"),
                "MOFA+ top feature loadings per active factor",
            ),
            metrics={
                "n_cells": n_cells,
                "n_cells_trained": n_trained,
                "n_factors": n_factors,
                "n_active_factors": n_active,
                "n_shared_factors": n_shared,
                "n_modality_specific_factors": n_spec,
                "n_excluded_factors": n_excl,
                "excluded_factor_ids": excluded_ids,
                "rna_r2_pct": rna_r2_pct,
                "atac_r2_pct": atac_r2_pct,
                "rna_additive_r2_pct": rna_add_pct,
                "atac_additive_r2_pct": atac_add_pct,
            },
        )

    @staticmethod
    def _peak_to_gene_block(
        p2g: dict[str, Any],
        findings: dict[str, Any],
        *,
        attach_triangulation_fig: bool = False,
    ) -> NarrativeBlock:
        n_cells = int(p2g.get("n_cells") or 0)
        dist_kb = int(p2g.get("distance_kb") or 150)
        n_tested = int(p2g.get("n_tested") or 0)
        n_links = int(p2g.get("n_links") or 0)
        n_pos = int(p2g.get("n_positive_correlations") or 0)
        n_neg = int(p2g.get("n_negative_correlations") or 0)
        n_prom = int(p2g.get("n_promoter_proximal") or 0)
        n_dist = int(p2g.get("n_distal_cis") or 0)

        evidence: list[EvidenceItem] = [
            EvidenceItem("shared cells", n_cells, source="peak_to_gene"),
            EvidenceItem("cis window kb", dist_kb, source="peak_to_gene"),
            EvidenceItem("tested peak-gene pairs", n_tested, source="peak_to_gene"),
            EvidenceItem("BH-FDR cis links", n_links, source="peak_to_gene"),
            EvidenceItem("positive cis links", n_pos, source="peak_to_gene"),
            EvidenceItem("negative cis links", n_neg, source="peak_to_gene"),
            EvidenceItem("promoter-proximal links", n_prom, source="peak_to_gene"),
            EvidenceItem("distal cis links", n_dist, source="peak_to_gene"),
        ]

        top_genes_listed: list[str] = []
        for row in (p2g.get("top_links") or [])[:5]:
            if not isinstance(row, dict):
                continue
            gene = str(row.get("gene") or "").strip()
            peak = str(row.get("peak") or "").strip()
            raw_corr = (
                row.get("correlation")
                if row.get("correlation") is not None
                else row.get("corr")
            )
            corr = _finite_round(raw_corr, 4)
            if not gene or corr is None:
                continue
            if gene not in top_genes_listed:
                top_genes_listed.append(gene)
            evidence.append(
                EvidenceItem(
                    f"candidate cis link {gene}",
                    corr,
                    source="peak_to_gene",
                    facts=[
                        SemanticFact(
                            subject=gene,
                            subject_type="gene",
                            predicate="cis_regulatory_association",
                            polarity="affirmed",
                            value=corr,
                            source="peak_to_gene",
                        )
                    ],
                )
            )
            if peak:
                evidence.append(
                    EvidenceItem(
                        f"candidate cis peak for {gene}",
                        peak,
                        source="peak_to_gene",
                    )
                )

        caveats: list[Caveat] = [
            Caveat(
                text=(
                    "Peak-to-gene cis links measure observational accessibility-expression "
                    "co-variation across shared cells against GC- and accessibility-matched "
                    "trans-chromosome nulls; physical 3D looping and enhancer activity "
                    "remain unverified without chromatin conformation or perturbation assays."
                ),
                severity="info",
            )
        ]

        confidence = (
            "high" if (n_links > 0 and n_cells >= 500)
            else ("medium" if n_cells >= 100 else "low")
        )
        figs = (
            _figure_entry(
                findings,
                "regulatory_triangulation_summary",
                "Candidate cis peak-to-gene associations and regulatory triangulation summary.",
            )
            if attach_triangulation_fig
            else []
        )
        claim = (
            f"Cis peak-to-gene linkage within {dist_kb} kb tested {n_tested} "
            f"candidate pairs across {n_cells} shared cells and identified "
            f"{n_links} BH-FDR cis associations ({n_pos} positive and {n_neg} "
            f"negative; {n_prom} promoter-proximal and {n_dist} distal cis)."
        )
        return NarrativeBlock(
            id="multiome.peak_to_gene",
            modality="integration",
            analysis="peak_to_gene",
            block_type="result",
            title="Candidate cis peak-to-gene associations",
            status="success",
            confidence=confidence,
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            figures=figs,
            tables=_table_entry(
                findings,
                "peak_gene_links",
                p2g.get("output_path"),
                "Signac-faithful BH-FDR cis peak-to-gene links",
            ),
            metrics={
                "n_cells": n_cells,
                "distance_kb": dist_kb,
                "n_tested": n_tested,
                "n_links": n_links,
                "n_positive_correlations": n_pos,
                "n_negative_correlations": n_neg,
                "n_promoter_proximal": n_prom,
                "n_distal_cis": n_dist,
                "top_genes": top_genes_listed[:4],
            },
        )

    @staticmethod
    def _eregulons_block(
        ereg: dict[str, Any],
        findings: dict[str, Any],
    ) -> NarrativeBlock:
        n_cells = int(ereg.get("n_shared_cells") or ereg.get("n_cells") or 0)
        n_tested = int(ereg.get("n_candidate_modules_tested") or 0)
        n_ereg = int(ereg.get("n_eregulons") or 0)
        n_trip = int(ereg.get("n_triplets") or 0)
        n_tfs = int(ereg.get("n_unique_tfs") or 0)
        n_pos = int(ereg.get("n_pos_tf_pos_p2g") or 0)
        n_neg = int(ereg.get("n_neg_tf_pos_p2g") or 0)
        n_paralog = int(ereg.get("n_paralog_ambiguous") or 0)
        n_depth = int(ereg.get("n_depth_confounded") or 0)
        motif_src = _format_motif_source(ereg.get("motif_source"))

        evidence: list[EvidenceItem] = [
            EvidenceItem("shared cells", n_cells, source="eregulons"),
            EvidenceItem("candidate modules tested", n_tested, source="eregulons"),
            EvidenceItem("BH-FDR candidate eRegulons", n_ereg, source="eregulons"),
            EvidenceItem("TF-enhancer-gene triplets", n_trip, source="eregulons"),
            EvidenceItem("unique TFs", n_tfs, source="eregulons"),
            EvidenceItem("concordant pos_tf_pos_p2g modules", n_pos, source="eregulons"),
            EvidenceItem("inverse neg_tf_pos_p2g modules", n_neg, source="eregulons"),
            EvidenceItem("paralog-ambiguous modules", n_paralog, source="eregulons"),
            EvidenceItem("depth-confounded modules", n_depth, source="eregulons"),
            EvidenceItem("motif collection", motif_src, source="eregulons"),
        ]

        top_tfs_listed: list[str] = []
        for mod in (ereg.get("top_eregulons") or [])[:5]:
            if not isinstance(mod, dict):
                continue
            tf = str(mod.get("tf") or "").strip()
            raw_tf_r = (
                mod.get("mean_r_tf_gene")
                if mod.get("mean_r_tf_gene") is not None
                else mod.get("tf_gene_corr")
            )
            tf_r = _finite_round(raw_tf_r, 4)
            n_tg = int(mod.get("n_target_genes") or 0)
            if not tf or tf_r is None:
                continue
            if tf not in top_tfs_listed:
                top_tfs_listed.append(tf)
            evidence.append(
                EvidenceItem(
                    f"candidate eRegulon {tf}",
                    tf_r,
                    source="eregulons",
                    facts=[
                        SemanticFact(
                            subject=tf,
                            subject_type="tf_regulator",
                            predicate="candidate_eregulon",
                            polarity="affirmed",
                            value=tf_r,
                            source="eregulons",
                        )
                    ],
                )
            )
            evidence.append(
                EvidenceItem(
                    f"eRegulon {tf} target genes",
                    n_tg,
                    source="eregulons",
                )
            )

        caveats: list[Caveat] = [
            Caveat(
                text=(
                    "Candidate eRegulons combine motif sequence matches, GC-stratified "
                    "cistrome enrichment, and shared-cell TF-target co-expression "
                    "(pos_tf_pos_p2g and neg_tf_pos_p2g); motif presence and "
                    "co-variation do not establish physical TF occupancy or "
                    "transcriptional causation."
                ),
                severity="info",
            )
        ]
        if n_paralog > 0:
            caveats.append(
                Caveat(
                    text=(
                        f"{n_paralog} candidate eRegulon module(s) share position "
                        "weight matrices across paralogous TF family members "
                        "(motif_family_ambiguity) and require ChIP/CUT&RUN or "
                        "perturbation disambiguation."
                    ),
                    severity="warning",
                )
            )
        if n_depth > 0:
            caveats.append(
                Caveat(
                    text=(
                        f"{n_depth} candidate eRegulon module(s) showed elevated "
                        "correlation with library sequencing depth (depth_confounded)."
                    ),
                    severity="warning",
                )
            )

        confidence = (
            "medium"
            if (n_ereg == 0 or n_depth > 0 or n_cells < 500)
            else "high"
        )
        tables = (
            _table_entry(
                findings,
                "eregulons",
                ereg.get("output_path"),
                "Candidate eRegulon modules summary",
            )
            + _table_entry(
                findings,
                "eregulon_triplets",
                ereg.get("triplets_path"),
                "Pruned TF-enhancer-gene regulatory triplets",
            )
        )
        claim = (
            f"Candidate eRegulon inference across {n_cells} shared cells tested "
            f"{n_tested} motif-linked modules and retained {n_ereg} BH-FDR "
            f"candidate eRegulons ({n_pos} pos_tf_pos_p2g and {n_neg} "
            f"neg_tf_pos_p2g) spanning {n_tfs} unique TFs and {n_trip} "
            f"TF-enhancer-gene triplets ({n_paralog} paralog-ambiguous, "
            f"{n_depth} depth-confounded)."
        )
        return NarrativeBlock(
            id="multiome.eregulons",
            modality="integration",
            analysis="eregulons",
            block_type="result",
            title="Candidate eRegulon regulatory modules",
            status="success",
            confidence=confidence,
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            figures=_figure_entry(
                findings,
                "regulatory_triangulation_summary",
                "Cross-modal regulatory triangulation of cis peak-to-gene links and candidate eRegulons.",
            ),
            tables=tables,
            metrics={
                "n_shared_cells": n_cells,
                "n_candidate_modules_tested": n_tested,
                "n_eregulons": n_ereg,
                "n_triplets": n_trip,
                "n_unique_tfs": n_tfs,
                "n_pos_tf_pos_p2g": n_pos,
                "n_neg_tf_pos_p2g": n_neg,
                "n_paralog_ambiguous": n_paralog,
                "n_depth_confounded": n_depth,
                "top_tfs": top_tfs_listed[:4],
            },
        )

    def methods(
        self,
        agent_name: str,
        agent_result: dict,
        context: dict | None = None,
    ) -> list[str]:
        """Return copy-pasteable manuscript Methods paragraphs strictly for executed stages."""
        findings = unwrap_integration_findings(agent_result)
        paragraphs: list[str] = []

        pqc = findings.get("paired_qc")
        if _stage_state(pqc) == "ok" and isinstance(pqc, dict):
            rna_s = pqc.get("rna_qc_summary") or {}
            atac_s = pqc.get("atac_qc_summary") or {}
            min_genes = rna_s.get("min_genes", 200)
            max_mt = rna_s.get("max_pct_mt", 20.0)
            min_peaks = atac_s.get("min_peaks", 200)
            uncomp = [
                _clean_token(x) for x in (pqc.get("metrics_not_computed") or []) if x
            ]
            uncomp_txt = (
                f" Uncomputed fragment-dependent metrics ({', '.join(uncomp)}) "
                "were recorded explicitly as null rather than imputed."
                if uncomp
                else ""
            )
            paragraphs.append(
                f"Paired multiome quality control: Same-cell barcode pairing was "
                f"verified prior to integration, retaining {pqc.get('n_cells_after', '?')} "
                f"of {pqc.get('n_cells_before', '?')} shared barcodes "
                f"(min_genes >= {min_genes}, mitochondrial fraction <= {max_mt}%, "
                f"min_peaks >= {min_peaks}) across {pqc.get('n_genes', '?')} RNA genes "
                f"and {pqc.get('n_peaks', '?')} ATAC peaks.{uncomp_txt}"
            )

        wnn = findings.get("wnn")
        if _stage_state(wnn) == "ok" and isinstance(wnn, dict):
            k = wnn.get("k", 20)
            res = wnn.get("resolution", 0.5)
            guard = _clean_token(wnn.get("guardrail_status"), "skipped")
            paragraphs.append(
                f"Weighted Nearest Neighbors (WNN) integration: RNA counts were "
                f"library-size log-normalized, subset to highly variable genes, "
                f"projected via PCA, and L2-normalized per cell; ATAC peak counts "
                f"were transformed via log-TF-IDF and truncated SVD/LSI, dropping "
                f"LSI components strongly correlated with log10 sequencing depth "
                f"before row-wise L2 normalization. Joint WNN graph construction "
                f"(muon, k={k}) and Leiden clustering (resolution={res}) were "
                f"audited by ARIA's 5-axis anti-degradation guardrail "
                f"(guardrail_status={guard}), treating WNN as a complementary "
                f"multi-modal lens alongside single-modality clusterings."
            )

        mofa = findings.get("mofa")
        if _stage_state(mofa) == "ok" and isinstance(mofa, dict):
            n_fac = mofa.get("n_factors", "?")
            n_act = mofa.get("n_active_factors", "?")
            n_excl = len(mofa.get("excluded_factors") or [])
            sub_clause = ""
            if mofa.get("subsampled"):
                sub_clause = (
                    f" Stratified Horvitz-Thompson subsampling trained on "
                    f"{mofa.get('n_cells_trained', '?')} of {mofa.get('n_cells', '?')} "
                    f"cells followed by closed-form ridge posterior projection to "
                    f"all cells."
                )
            paragraphs.append(
                f"MOFA+ latent factor decomposition: Joint automatic relevance "
                f"determination factor analysis (mofapy2, {n_fac} initial factors) "
                f"was fit over feature-standardized log1p-CP10k highly variable "
                f"RNA genes and TF-IDF hypervariable ATAC peaks.{sub_clause} "
                f"Factors exhibiting sequencing-depth association (|r| >= 0.70) "
                f"or cell-cycle enrichment were quarantined ({n_excl} excluded), "
                f"yielding {n_act} active biological factors with non-additive "
                f"total variance explained reported per view."
            )

        p2g = findings.get("peak_to_gene")
        if _stage_state(p2g) == "ok" and isinstance(p2g, dict):
            dkb = p2g.get("distance_kb", 150)
            min_r = p2g.get("min_corr", 0.10)
            max_q = p2g.get("max_padj", 0.05)
            bg_m = _clean_token(p2g.get("background_method"), "trans_chromosome_gc_accessibility")
            paragraphs.append(
                f"Cis peak-to-gene linkage: Pearson correlations between log1p-CP10k "
                f"gene expression and log-TF-IDF peak accessibility were evaluated "
                f"across shared cells for peaks within +-{dkb} kb of canonical "
                f"transcription start sites. Empirical Z-scores and one-sided "
                f"normal tail p-values were calibrated against reference-FASTA "
                f"GC-fraction and log10-accessibility matched trans-chromosome "
                f"background peaks ({bg_m}), with Benjamini-Hochberg FDR <= {max_q} "
                f"and |r| >= {min_r} defining associative candidate cis links."
            )

        ereg = findings.get("eregulons")
        if _stage_state(ereg) == "ok" and isinstance(ereg, dict):
            win_bp = ereg.get("scan_window_bp", 500)
            mpval = ereg.get("motif_pvalue", 1e-4)
            msrc = _format_motif_source(ereg.get("motif_source"))
            backend = _clean_token(ereg.get("scanner_backend"), "pwm_scanner")
            paragraphs.append(
                f"Candidate eRegulon inference: Summit-centered {win_bp} bp cis "
                f"peak windows were scanned on both DNA strands ({backend}, "
                f"{msrc}, motif p <= {mpval}) and tested for 5-stratum "
                f"GC-matched Cochran-Mantel-Haenszel cistrome enrichment. "
                f"TF-to-target expression associations were calibrated against a "
                f"5x5 (mean, detection-rate) stratified null, controlled by "
                f"Benjamini-Hochberg FDR across all tested (TF, motif, sign_class) "
                f"hypotheses prior to paralog collapse, pruned via weighted "
                f"leading-edge recovery, and scored per cell by RNA and ATAC AUCell "
                f"rankings as non-causal associative candidate modules."
            )

        return paragraphs

    def figures(
        self,
        agent_name: str,
        agent_result: dict,
        report_dir: Path | None = None,
    ) -> list[dict]:
        findings = unwrap_integration_findings(agent_result)
        figs = findings.get("figures") or {}
        figs_svg = findings.get("figures_svg") or {}
        out: list[dict] = []
        for key, path in figs.items():
            if not isinstance(path, str) or not path:
                continue
            rec = {
                "id": f"multiome.{key}",
                "path": path,
                "caption": key.replace("_", " "),
            }
            if isinstance(figs_svg.get(key), str):
                rec["svg_path"] = figs_svg[key]
            out.append(rec)
        return out

    def tables(
        self,
        agent_name: str,
        agent_result: dict,
        report_dir: Path | None = None,
    ) -> list[dict]:
        findings = unwrap_integration_findings(agent_result)
        out: list[dict] = []
        out.extend(
            _table_entry(
                findings,
                "mofa_factor_loadings",
                (findings.get("mofa") or {}).get("factor_loadings_path"),
                "MOFA+ top feature loadings",
            )
        )
        out.extend(
            _table_entry(
                findings,
                "peak_gene_links",
                (findings.get("peak_to_gene") or {}).get("output_path"),
                "Cis peak-to-gene links",
            )
        )
        out.extend(
            _table_entry(
                findings,
                "eregulons",
                (findings.get("eregulons") or {}).get("output_path"),
                "Candidate eRegulon modules",
            )
        )
        out.extend(
            _table_entry(
                findings,
                "eregulon_triplets",
                (findings.get("eregulons") or {}).get("triplets_path"),
                "Pruned TF-enhancer-gene triplets",
            )
        )
        return out
