"""Deterministic Cross-Modal Evidence Graph & Regulatory Triangulation (V48 I7).

Builds an auditable, zero-LLM evidence graph linking structured outputs from:
  - ``integration_agent`` (paired QC, WNN + 5-axis anti-degradation guardrail,
    MOFA+ active & excluded confounder factors, Signac-faithful ``peak_to_gene``
    cis links, SCENIC+-style candidate ``eregulons``)
  - optional ``scrna_agent`` / ``bulk_rna_agent`` differential expression
  - optional ``chromatin_agent`` differential accessibility & motif enrichment

Core governance invariants (ADR-064, ADR-070, C1-C4):
  1. Associative, never causal: every edge carries ``non_causal=True``; claims,
     evidence items, and caveats use strictly observational candidate phrasing.
  2. Non-circular convergence: ``eregulons`` consumes ``peak_to_gene`` links on
     the same shared cell barcodes. That lineage is recorded explicitly via
     ``shared_upstream=["peak_to_gene"]`` and ``shared_cell_universe=True`` so
     pipeline refinement is never double-counted as two independent cohorts.
  3. Absence != discordance: unexecuted upstream contrasts are recorded under
     ``absent_not_tested`` and never misclassified as discordant biology.
  4. Lens, not replacement: WNN 5-axis guardrail status (``pass``, ``warn``,
     ``degraded``, ``skipped``), paralog motif-sharing ambiguity, sequencing-depth
     confounding, and excluded MOFA+ technical/cell-cycle factors remain visible.
  5. Portable serialization: zero machine-absolute filesystem paths enter the
     exported graph manifest.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from aria.agents.narrative.types import (
    Caveat,
    EvidenceItem,
    NarrativeBlock,
    SemanticFact,
)

SCHEMA_VERSION = "1.0.0"

_EXECUTED_STATUSES = {"success", "ok", "done"}
_NON_RUN_STATUSES = {"not_run", "skipped", "blocked", "error", "failed"}


def _is_executed_stage(stage: Any) -> bool:
    if not isinstance(stage, dict) or not stage:
        return False
    st = str(stage.get("status") or "").lower().strip()
    if st in _NON_RUN_STATUSES:
        return False
    if stage.get("ran") is False:
        return False
    if st in _EXECUTED_STATUSES or stage.get("ran") is True:
        return True
    # Structured payload without explicit failure status.
    return bool(stage)


def _finite_float(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        out = float(val)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _round_num(val: Any, digits: int = 4) -> float | None:
    f = _finite_float(val)
    return round(f, digits) if f is not None else None


def _link_corr(link: dict[str, Any]) -> float | None:
    raw = link.get("correlation") if link.get("correlation") is not None else link.get("corr")
    return _finite_float(raw)


def _portable_name(path_val: Any) -> str | None:
    if not path_val:
        return None
    text = str(path_val).strip()
    if not text:
        return None
    return Path(text).name


def _unwrap_integration(agent_results: dict | None) -> dict[str, Any]:
    integ = (agent_results or {}).get("integration_agent")
    if not isinstance(integ, dict):
        return {}
    findings = integ.get("findings")
    if isinstance(findings, dict):
        return findings
    return integ


def _extract_rna_de_genes(agent_results: dict | None) -> tuple[bool, dict[str, dict[str, Any]]]:
    """Return ``(ran, {UPPER_GENE: {'symbol': str, 'log2fc': float, 'source': str}})``."""
    out: dict[str, dict[str, Any]] = {}
    ran = False

    bulk = (agent_results or {}).get("bulk_rna_agent")
    if isinstance(bulk, dict) and str(bulk.get("status") or "").lower() == "done":
        bf = bulk.get("findings", bulk)
        if isinstance(bf, dict):
            for contrast in bf.get("contrasts") or []:
                if not isinstance(contrast, dict):
                    continue
                if str(contrast.get("status") or "").lower() != "success":
                    continue
                ran = True
                for row in contrast.get("top_genes") or []:
                    if not isinstance(row, dict):
                        continue
                    sym = str(row.get("symbol") or row.get("gene") or "").strip()
                    lfc = _finite_float(row.get("log2fc"))
                    if sym and lfc is not None:
                        key = sym.upper()
                        prev = out.get(key)
                        if prev is None or abs(lfc) > abs(prev["log2fc"]):
                            out[key] = {
                                "symbol": sym,
                                "log2fc": round(lfc, 4),
                                "source": "bulk_rna_de",
                            }

    scrna = (agent_results or {}).get("scrna_agent")
    if isinstance(scrna, dict) and str(scrna.get("status") or "").lower() == "done":
        try:
            from aria.agents import _narrative_scrna
            sf = _narrative_scrna.unwrap_scrna_findings(scrna) or {}
        except Exception:
            sf = scrna.get("findings", scrna) if isinstance(scrna, dict) else {}
        pb = (sf.get("pseudobulk_de") or {}) if isinstance(sf, dict) else {}
        per_group = pb.get("per_group") if isinstance(pb, dict) else None
        if isinstance(per_group, dict):
            for grp_data in per_group.values():
                if not isinstance(grp_data, dict):
                    continue
                per_comp = grp_data.get("per_comparison") or {}
                for comp in per_comp.values():
                    if not isinstance(comp, dict):
                        continue
                    if str(comp.get("status") or "").lower() != "success":
                        continue
                    ran = True
                    for row in comp.get("top_genes") or []:
                        if not isinstance(row, dict):
                            continue
                        sym = str(row.get("gene") or row.get("symbol") or "").strip()
                        lfc = _finite_float(row.get("log2fc"))
                        if sym and lfc is not None:
                            key = sym.upper()
                            prev = out.get(key)
                            if prev is None or abs(lfc) > abs(prev["log2fc"]):
                                out[key] = {
                                    "symbol": sym,
                                    "log2fc": round(lfc, 4),
                                    "source": "scrna_pseudobulk_de",
                                }
    return ran, out


def _extract_chromatin_signals(
    agent_results: dict | None,
) -> tuple[bool, dict[str, dict[str, Any]], bool, set[str]]:
    """Return ``(da_ran, da_peaks_by_id, motifs_ran, enriched_tf_upper)``."""
    da_peaks: dict[str, dict[str, Any]] = {}
    enriched_tfs: set[str] = set()
    da_ran = False
    motifs_ran = False

    chrom = (agent_results or {}).get("chromatin_agent")
    if not isinstance(chrom, dict) or str(chrom.get("status") or "").lower() != "done":
        return False, da_peaks, False, enriched_tfs

    try:
        from aria.agents.narrative.narrators.chromatin import (
            unwrap_chromatin_findings,
        )
        cf = unwrap_chromatin_findings(chrom) or {}
    except Exception:
        cf = chrom.get("findings", chrom) if isinstance(chrom, dict) else {}

    da = cf.get("differential_accessibility") if isinstance(cf, dict) else None
    if isinstance(da, dict):
        pb = da.get("pseudobulk") or {}
        if isinstance(pb, dict) and str(pb.get("status") or "").lower() == "success":
            da_ran = True
            for row in pb.get("top_peaks") or []:
                if not isinstance(row, dict):
                    continue
                pid = str(row.get("peak") or row.get("feature") or "").strip()
                lfc = _finite_float(row.get("log2fc"))
                if pid and lfc is not None:
                    da_peaks[pid] = {"peak": pid, "log2fc": round(lfc, 4)}
        pc = da.get("per_cluster") or {}
        if isinstance(pc, dict) and pc.get("ran"):
            da_ran = True
            top_by_cl = pc.get("top_peaks_by_cluster") or {}
            if isinstance(top_by_cl, dict):
                for rows in top_by_cl.values():
                    for row in rows or []:
                        if not isinstance(row, dict):
                            continue
                        pid = str(row.get("peak") or row.get("feature") or "").strip()
                        lfc = _finite_float(row.get("log2fc") or row.get("score"))
                        if pid and pid not in da_peaks:
                            da_peaks[pid] = {
                                "peak": pid,
                                "log2fc": round(lfc, 4) if lfc is not None else 1.0,
                            }

    motifs = cf.get("motifs") if isinstance(cf, dict) else None
    if isinstance(motifs, dict) and motifs.get("ran"):
        motifs_ran = True
        per_group = motifs.get("per_group") or {}
        if isinstance(per_group, dict):
            for grp in per_group.values():
                if not isinstance(grp, dict):
                    continue
                for row in grp.get("top_motifs") or []:
                    if not isinstance(row, dict):
                        continue
                    tf_name = str(
                        row.get("tf") or row.get("motif_name") or ""
                    ).strip()
                    if tf_name:
                        for part in tf_name.replace("::", ",").split(","):
                            clean = part.split("(")[0].strip().upper()
                            if clean:
                                enriched_tfs.add(clean)
    return da_ran, da_peaks, motifs_ran, enriched_tfs


def build_cross_modal_evidence_graph(
    agent_results: dict | None,
    exp_ctx: dict | None = None,
) -> dict[str, Any]:
    """Construct a deterministic, non-circular cross-modal evidence graph.

    Parameters
    ----------
    agent_results
        Mapping of agent names to structured result envelopes.
    exp_ctx
        Optional experiment context (used only for assembly/organism metadata,
        never for biological priors).
    """
    findings = _unwrap_integration(agent_results)
    pqc = findings.get("paired_qc") if isinstance(findings, dict) else None
    wnn = findings.get("wnn") if isinstance(findings, dict) else None
    mofa = findings.get("mofa") if isinstance(findings, dict) else None
    p2g = findings.get("peak_to_gene") if isinstance(findings, dict) else None
    ereg = findings.get("eregulons") if isinstance(findings, dict) else None

    stage_flags = {
        "paired_qc": _is_executed_stage(pqc),
        "wnn": _is_executed_stage(wnn),
        "mofa": _is_executed_stage(mofa),
        "peak_to_gene": _is_executed_stage(p2g),
        "eregulons": _is_executed_stage(ereg),
    }
    executed_stages = [k for k, v in stage_flags.items() if v]

    rna_de_ran, rna_de_genes = _extract_rna_de_genes(agent_results)
    atac_da_ran, atac_da_peaks, atac_motifs_ran, atac_enriched_tfs = (
        _extract_chromatin_signals(agent_results)
    )

    absent_not_tested: list[str] = []
    for stage_name, did_run in stage_flags.items():
        if not did_run:
            absent_not_tested.append(f"integration.{stage_name}")
    if not rna_de_ran:
        absent_not_tested.append("modality.rna_differential_expression")
    if not atac_da_ran:
        absent_not_tested.append("modality.atac_differential_accessibility")
    if not atac_motifs_ran:
        absent_not_tested.append("modality.atac_motif_enrichment")

    # Graph runs iff integration executed at least two stages or one regulatory
    # integration stage alongside an upstream modality layer.
    ran = (
        len(executed_stages) >= 2
        or (
            bool(stage_flags["peak_to_gene"] or stage_flags["eregulons"])
            and (rna_de_ran or atac_da_ran or atac_motifs_ran)
        )
    )

    nodes: dict[str, dict[str, Any]] = {}
    edges: list[dict[str, Any]] = []
    regulatory_axes: list[dict[str, Any]] = []

    def _put_node(node_id: str, node_type: str, label: str, **attrs: Any) -> None:
        existing = nodes.get(node_id)
        if existing is None:
            nodes[node_id] = {
                "id": node_id,
                "node_type": node_type,
                "label": label,
                **attrs,
            }
        else:
            for k, v in attrs.items():
                if v is not None and k not in existing:
                    existing[k] = v

    for stage_name, did_run in stage_flags.items():
        raw_stage = findings.get(stage_name) if isinstance(findings, dict) else None
        raw_status = (
            str(raw_stage.get("status") or "executed")
            if isinstance(raw_stage, dict) and raw_stage
            else "absent_not_tested"
        )
        _put_node(
            f"state:{stage_name}",
            "modality_state",
            stage_name,
            status="executed" if did_run else raw_status,
            artifact=_portable_name(
                raw_stage.get("output_path") if isinstance(raw_stage, dict) else None
            ),
        )

    # WNN 5-axis anti-degradation guardrail node ("Lens, not replacement").
    guardrail_status = "not_evaluated"
    wnn_discordant_clusters = 0
    if stage_flags["wnn"] and isinstance(wnn, dict):
        guardrail_status = str(wnn.get("guardrail_status") or "skipped")
        wnn_discordant_clusters = int(wnn.get("n_discordant_clusters") or 0)
        anti_deg = wnn.get("anti_degradation") or {}
        axis_statuses = {}
        if isinstance(anti_deg, dict):
            for ax_name, ax_payload in (anti_deg.get("axes") or {}).items():
                if isinstance(ax_payload, dict):
                    axis_statuses[str(ax_name)] = str(
                        ax_payload.get("status") or "unknown"
                    )
        _put_node(
            "guardrail:wnn_anti_degradation",
            "guardrail",
            "WNN 5-axis anti-degradation guardrail",
            guardrail_status=guardrail_status,
            n_discordant_clusters=wnn_discordant_clusters,
            wnn_vs_rna_ari=_round_num(wnn.get("wnn_vs_rna_ari")),
            wnn_vs_atac_ari=_round_num(wnn.get("wnn_vs_atac_ari")),
            axis_statuses=axis_statuses,
            lens_not_replacement=True,
        )

    # MOFA+ active & excluded confounder factor nodes.
    n_mofa_shared = 0
    n_mofa_excluded = 0
    if stage_flags["mofa"] and isinstance(mofa, dict):
        for fac in mofa.get("top_factors") or []:
            if not isinstance(fac, dict):
                continue
            if fac.get("active") is False or bool(
                fac.get("excluded_from_biological_interpretation")
            ):
                continue
            fid = fac.get("factor_id") if fac.get("factor_id") is not None else fac.get("factor")
            if fid is None:
                continue
            cls_name = str(
                fac.get("category") or fac.get("classification") or "active"
            )
            if cls_name == "shared":
                n_mofa_shared += 1
            r2_view = fac.get("r2_by_view") if isinstance(fac.get("r2_by_view"), dict) else {}
            rna_v = fac.get("rna_r2") if fac.get("rna_r2") is not None else r2_view.get("rna")
            atac_v = fac.get("atac_r2") if fac.get("atac_r2") is not None else r2_view.get("atac")
            _put_node(
                f"mofa:Factor_{fid}",
                "mofa_factor",
                f"Factor {fid}",
                status="active",
                classification=cls_name,
                rna_r2=_round_num(rna_v, 6),
                atac_r2=_round_num(atac_v, 6),
            )
        for ex in mofa.get("excluded_factors") or []:
            if not isinstance(ex, dict):
                continue
            fid = ex.get("factor_id") if ex.get("factor_id") is not None else ex.get("factor")
            if fid is None:
                continue
            reasons = [str(r) for r in (ex.get("exclusion_reasons") or [])]
            r2_view = ex.get("r2_by_view") if isinstance(ex.get("r2_by_view"), dict) else {}
            rna_v = ex.get("rna_r2") if ex.get("rna_r2") is not None else r2_view.get("rna")
            atac_v = ex.get("atac_r2") if ex.get("atac_r2") is not None else r2_view.get("atac")
            n_mofa_excluded += 1
            _put_node(
                f"mofa:Factor_{fid}",
                "mofa_factor",
                f"Factor {fid}",
                status="excluded_confounder",
                classification="excluded_confounder",
                exclusion_reasons=reasons,
                rna_r2=_round_num(rna_v, 6),
                atac_r2=_round_num(atac_v, 6),
            )
            regulatory_axes.append({
                "axis_id": f"mofa_confounder:Factor_{fid}",
                "axis_kind": "latent_factor_confounder",
                "classification": "caveated_ambiguous",
                "subject": f"Factor {fid}",
                "target_gene": None,
                "cis_peak": None,
                "tf_regulator": None,
                "sign": "confounded",
                "shared_cell_universe": True,
                "shared_upstream": [],
                "independent_contrast_support": False,
                "reasons": reasons or ["technical_or_cell_cycle_confounder"],
                "non_causal": True,
            })

    # Index top peak_to_gene links by gene (upper) and (gene, peak).
    p2g_by_gene: dict[str, list[dict[str, Any]]] = {}
    if stage_flags["peak_to_gene"] and isinstance(p2g, dict):
        for link in p2g.get("top_links") or []:
            if not isinstance(link, dict):
                continue
            gene = str(link.get("gene") or "").strip()
            peak = str(link.get("peak") or "").strip()
            corr = _link_corr(link)
            padj = _finite_float(link.get("padj"))
            if not gene or not peak or corr is None:
                continue
            region_class = str(link.get("region_class") or "cis")
            dist_bp = link.get("distance_bp")
            _put_node(f"gene:{gene}", "gene", gene)
            _put_node(
                f"peak:{peak}",
                "cis_peak",
                peak,
                region_class=region_class,
                distance_bp=int(dist_bp) if isinstance(dist_bp, (int, float)) else None,
            )
            sign = "positive" if corr >= 0 else "negative"
            edges.append({
                "source": f"peak:{peak}",
                "target": f"gene:{gene}",
                "relation": "peak_to_gene_association",
                "sign": sign,
                "corr": round(corr, 4),
                "padj": _round_num(padj, 6),
                "region_class": region_class,
                "shared_cell_universe": True,
                "shared_upstream": [],
                "non_causal": True,
            })
            p2g_by_gene.setdefault(gene.upper(), []).append(link)

            if rna_de_ran and gene.upper() in rna_de_genes:
                de_info = rna_de_genes[gene.upper()]
                edges.append({
                    "source": "state:rna_de",
                    "target": f"gene:{gene}",
                    "relation": "rna_de_association",
                    "sign": "positive" if de_info["log2fc"] >= 0 else "negative",
                    "log2fc": de_info["log2fc"],
                    "shared_cell_universe": False,
                    "shared_upstream": [],
                    "non_causal": True,
                })
            if atac_da_ran and peak in atac_da_peaks:
                da_info = atac_da_peaks[peak]
                edges.append({
                    "source": "state:atac_da",
                    "target": f"peak:{peak}",
                    "relation": "atac_da_association",
                    "sign": "positive" if da_info["log2fc"] >= 0 else "negative",
                    "log2fc": da_info["log2fc"],
                    "shared_cell_universe": False,
                    "shared_upstream": [],
                    "non_causal": True,
                })

    # Candidate eRegulons (explicitly recording shared_upstream=["peak_to_gene"]).
    ereg_genes_covered: set[str] = set()
    if stage_flags["eregulons"] and isinstance(ereg, dict):
        for mod in ereg.get("top_eregulons") or []:
            if not isinstance(mod, dict):
                continue
            tf = str(mod.get("tf") or "").strip()
            if not tf:
                continue
            sign_class = str(mod.get("sign_class") or "pos_tf_pos_p2g")
            motif_id = str(
                mod.get("best_motif_id")
                or mod.get("motif_id")
                or ((mod.get("supporting_motif_ids") or [""])[0] if mod.get("supporting_motif_ids") else "")
            )
            paralog_ambig = bool(mod.get("motif_family_ambiguity"))
            depth_conf = bool(mod.get("depth_confounded"))
            raw_tf_r = (
                mod.get("mean_r_tf_gene")
                if mod.get("mean_r_tf_gene") is not None
                else mod.get("tf_gene_corr")
            )
            raw_padj = (
                mod.get("regulon_padj")
                if mod.get("regulon_padj") is not None
                else mod.get("padj")
            )
            raw_aucell = (
                mod.get("rna_atac_activity_corr")
                if mod.get("rna_atac_activity_corr") is not None
                else mod.get("rna_atac_aucell_corr")
            )
            tf_corr = _finite_float(raw_tf_r)
            padj = _finite_float(raw_padj)
            aucell_r = _finite_float(raw_aucell)
            target_genes = [
                str(g).strip()
                for g in (mod.get("top_target_genes") or mod.get("target_genes") or [])
                if str(g).strip()
            ]
            enhancer_peaks = [
                str(p).strip()
                for p in (mod.get("enhancer_peaks") or [])
                if str(p).strip()
            ]
            _put_node(
                f"tf:{tf}",
                "tf_regulator",
                tf,
                motif_id=motif_id or None,
                motif_family_ambiguity=paralog_ambig,
                depth_confounded=depth_conf,
            )

            for ep in enhancer_peaks[:4]:
                _put_node(f"peak:{ep}", "cis_peak", ep, region_class="cis")
                edges.append({
                    "source": f"tf:{tf}",
                    "target": f"peak:{ep}",
                    "relation": "eregulon_tf_to_peak",
                    "sign": "positive",
                    "motif_id": motif_id or None,
                    "shared_cell_universe": True,
                    "shared_upstream": ["peak_to_gene"],
                    "non_causal": True,
                })

            if atac_motifs_ran and tf.upper() in atac_enriched_tfs:
                edges.append({
                    "source": "state:atac_motifs",
                    "target": f"tf:{tf}",
                    "relation": "atac_motif_enrichment",
                    "sign": "positive",
                    "shared_cell_universe": False,
                    "shared_upstream": [],
                    "non_causal": True,
                })

            for tg in target_genes[:5]:
                ereg_genes_covered.add(tg.upper())
                _put_node(f"gene:{tg}", "gene", tg)
                edges.append({
                    "source": f"tf:{tf}",
                    "target": f"gene:{tg}",
                    "relation": "eregulon_tf_to_gene",
                    "sign": (
                        "negative"
                        if sign_class.startswith("neg_")
                        else "positive"
                    ),
                    "sign_class": sign_class,
                    "tf_gene_corr": _round_num(tf_corr, 4),
                    "padj": _round_num(padj, 6),
                    "rna_atac_aucell_corr": _round_num(aucell_r, 4),
                    "shared_cell_universe": True,
                    "shared_upstream": ["peak_to_gene"],
                    "non_causal": True,
                })

            # Classify candidate eRegulon axis without circular inflation.
            reasons: list[str] = []
            cross_contrast = False
            discordant_contrast = False
            if rna_de_ran:
                tf_de = rna_de_genes.get(tf.upper())
                tg_de_list = [
                    rna_de_genes[g.upper()]
                    for g in target_genes
                    if g.upper() in rna_de_genes
                ]
                if tf_de is not None and tg_de_list:
                    same_dir = any(
                        (tf_de["log2fc"] * g_de["log2fc"]) > 0
                        for g_de in tg_de_list
                    )
                    if sign_class == "pos_tf_pos_p2g" and same_dir:
                        cross_contrast = True
                        reasons.append("concordant_rna_de_tf_and_target")
                    elif sign_class == "pos_tf_pos_p2g" and not same_dir:
                        discordant_contrast = True
                        reasons.append("opposite_rna_de_tf_vs_target")
            if atac_motifs_ran and tf.upper() in atac_enriched_tfs:
                cross_contrast = True
                reasons.append("concordant_standalone_atac_motif_enrichment")

            if paralog_ambig:
                reasons.append("motif_family_ambiguity")
            if depth_conf:
                reasons.append("depth_confounded")
            if guardrail_status == "degraded":
                reasons.append("wnn_guardrail_degraded")
            elif guardrail_status == "warn":
                reasons.append("wnn_guardrail_warn")

            if sign_class == "neg_tf_pos_p2g" or discordant_contrast:
                axis_cls = "discordant_signal"
                if sign_class == "neg_tf_pos_p2g":
                    reasons.insert(0, "inverse_tf_expression_vs_accessible_peak2gene")
            elif paralog_ambig or depth_conf or guardrail_status == "degraded":
                axis_cls = "caveated_ambiguous"
            elif cross_contrast:
                axis_cls = "convergent_cross_contrast"
            elif (
                sign_class == "pos_tf_pos_p2g"
                and aucell_r is not None
                and aucell_r >= 0.20
            ):
                axis_cls = "convergent_paired_multiome"
                reasons.append("concordant_rna_atac_aucell_on_shared_cells")
            else:
                axis_cls = "associative_cis_candidate"
                reasons.append("candidate_eregulon_shared_cell_association")

            regulatory_axes.append({
                "axis_id": f"eregulon:{tf}:{sign_class}",
                "axis_kind": "candidate_eregulon",
                "classification": axis_cls,
                "subject": f"{tf} ({sign_class})",
                "tf_regulator": tf,
                "motif_id": motif_id or None,
                "target_gene": target_genes[0] if target_genes else None,
                "n_target_genes": int(mod.get("n_target_genes") or len(target_genes)),
                "n_enhancer_peaks": int(mod.get("n_enhancer_peaks") or 0),
                "cis_peak": None,
                "sign": (
                    "negative" if sign_class.startswith("neg_") else "positive"
                ),
                "tf_gene_corr": _round_num(tf_corr, 4),
                "padj": _round_num(padj, 6),
                "rna_atac_aucell_corr": _round_num(aucell_r, 4),
                "motif_family_ambiguity": paralog_ambig,
                "depth_confounded": depth_conf,
                "shared_cell_universe": True,
                "shared_upstream": ["peak_to_gene"],
                "independent_contrast_support": cross_contrast,
                "reasons": reasons,
                "non_causal": True,
            })

    # Classify standalone top peak_to_gene links not already covered by an eRegulon.
    if stage_flags["peak_to_gene"] and isinstance(p2g, dict):
        for link in (p2g.get("top_links") or [])[:15]:
            if not isinstance(link, dict):
                continue
            gene = str(link.get("gene") or "").strip()
            peak = str(link.get("peak") or "").strip()
            corr = _link_corr(link)
            padj = _finite_float(link.get("padj"))
            if not gene or not peak or corr is None:
                continue
            if gene.upper() in ereg_genes_covered:
                continue
            reasons = []
            de_info = rna_de_genes.get(gene.upper()) if rna_de_ran else None
            da_info = atac_da_peaks.get(peak) if atac_da_ran else None
            cross_contrast = False
            discordant_contrast = False
            if de_info is not None and da_info is not None:
                if (de_info["log2fc"] * da_info["log2fc"] > 0) and corr > 0:
                    cross_contrast = True
                    reasons.append("concordant_rna_de_and_atac_da")
                elif (de_info["log2fc"] * da_info["log2fc"] < 0):
                    discordant_contrast = True
                    reasons.append("opposite_rna_de_vs_atac_da")

            if corr < 0 or discordant_contrast:
                axis_cls = "discordant_signal"
                if corr < 0:
                    reasons.insert(0, "negative_peak_to_gene_correlation")
            elif guardrail_status == "degraded":
                axis_cls = "caveated_ambiguous"
                reasons.append("wnn_guardrail_degraded")
            elif cross_contrast:
                axis_cls = "convergent_cross_contrast"
            else:
                axis_cls = "associative_cis_candidate"
                reasons.append("cis_accessibility_expression_correlation")

            regulatory_axes.append({
                "axis_id": f"peak_to_gene:{peak}->{gene}",
                "axis_kind": "cis_peak_to_gene",
                "classification": axis_cls,
                "subject": f"{gene} ({peak})",
                "tf_regulator": None,
                "target_gene": gene,
                "cis_peak": peak,
                "region_class": str(link.get("region_class") or "cis"),
                "sign": "positive" if corr >= 0 else "negative",
                "corr": round(corr, 4),
                "padj": _round_num(padj, 6),
                "shared_cell_universe": True,
                "shared_upstream": [],
                "independent_contrast_support": cross_contrast,
                "reasons": reasons,
                "non_causal": True,
            })

    # Surface WNN-vs-RNA cluster discordance or degraded guardrail explicitly.
    if stage_flags["wnn"] and (
        wnn_discordant_clusters > 0 or guardrail_status == "degraded"
    ):
        reasons = []
        if wnn_discordant_clusters > 0:
            reasons.append("wnn_cluster_low_rna_purity")
        if guardrail_status == "degraded":
            reasons.append("wnn_anti_degradation_guardrail_degraded")
        regulatory_axes.append({
            "axis_id": "manifold:wnn_vs_modality_discordance",
            "axis_kind": "manifold_concordance",
            "classification": "discordant_signal",
            "subject": "WNN vs single-modality clustering",
            "tf_regulator": None,
            "target_gene": None,
            "cis_peak": None,
            "sign": "discordant",
            "n_discordant_clusters": wnn_discordant_clusters,
            "guardrail_status": guardrail_status,
            "shared_cell_universe": True,
            "shared_upstream": [],
            "independent_contrast_support": False,
            "reasons": reasons,
            "non_causal": True,
        })

    counts_by_class = {
        "convergent_cross_contrast": 0,
        "convergent_paired_multiome": 0,
        "associative_cis_candidate": 0,
        "discordant_signal": 0,
        "caveated_ambiguous": 0,
    }
    for ax in regulatory_axes:
        cls_key = str(ax.get("classification") or "")
        if cls_key in counts_by_class:
            counts_by_class[cls_key] += 1

    n_convergent_total = (
        counts_by_class["convergent_cross_contrast"]
        + counts_by_class["convergent_paired_multiome"]
    )

    sorted_nodes = sorted(
        nodes.values(), key=lambda n: (str(n.get("node_type")), str(n.get("id")))
    )
    sorted_edges = sorted(
        edges,
        key=lambda e: (
            str(e.get("relation")),
            str(e.get("source")),
            str(e.get("target")),
        ),
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "ran": ran,
        "non_causal_contract": True,
        "circular_convergence_prevented": True,
        "organism": str((exp_ctx or {}).get("organism") or "unspecified"),
        "genome": str((exp_ctx or {}).get("genome") or "unspecified"),
        "executed_integration_stages": executed_stages,
        "external_modality_contrasts": {
            "rna_differential_expression": rna_de_ran,
            "atac_differential_accessibility": atac_da_ran,
            "atac_motif_enrichment": atac_motifs_ran,
        },
        "absent_not_tested": sorted(set(absent_not_tested)),
        "wnn_guardrail": {
            "status": guardrail_status,
            "n_discordant_clusters": wnn_discordant_clusters,
            "lens_not_replacement": True,
        },
        "mofa_summary": {
            "n_shared_factors": n_mofa_shared,
            "n_excluded_confounder_factors": n_mofa_excluded,
        },
        "summary_counts": {
            "n_nodes": len(sorted_nodes),
            "n_edges": len(sorted_edges),
            "n_axes_total": len(regulatory_axes),
            "n_convergent_total": n_convergent_total,
            **counts_by_class,
        },
        "nodes": sorted_nodes,
        "edges": sorted_edges,
        "regulatory_axes": regulatory_axes,
    }


def build_triangulation_block(graph: dict[str, Any]) -> NarrativeBlock | None:
    """Emit a governed ``NarrativeBlock`` summarizing cross-modal triangulation."""
    if not isinstance(graph, dict) or not graph.get("ran"):
        return None

    counts = graph.get("summary_counts") or {}
    stages = list(graph.get("executed_integration_stages") or [])
    n_stages = len(stages)
    n_nodes = int(counts.get("n_nodes") or 0)
    n_edges = int(counts.get("n_edges") or 0)
    n_axes = int(counts.get("n_axes_total") or 0)
    n_conv = int(counts.get("n_convergent_total") or 0)
    n_cross = int(counts.get("convergent_cross_contrast") or 0)
    n_paired_conv = int(counts.get("convergent_paired_multiome") or 0)
    n_assoc = int(counts.get("associative_cis_candidate") or 0)
    n_disc = int(counts.get("discordant_signal") or 0)
    n_cav = int(counts.get("caveated_ambiguous") or 0)

    wnn_guard = graph.get("wnn_guardrail") or {}
    guard_status = str(wnn_guard.get("status") or "not_evaluated")

    evidence: list[EvidenceItem] = [
        EvidenceItem("integration stages executed", n_stages, source="cross_modal_graph"),
        EvidenceItem("evidence graph nodes", n_nodes, source="cross_modal_graph"),
        EvidenceItem("evidence graph edges", n_edges, source="cross_modal_graph"),
        EvidenceItem("triangulated axes", n_axes, source="cross_modal_graph"),
        EvidenceItem("convergent axes", n_conv, source="cross_modal_graph"),
        EvidenceItem("cross-contrast convergent axes", n_cross, source="cross_modal_graph"),
        EvidenceItem("paired multiome concordant axes", n_paired_conv, source="cross_modal_graph"),
        EvidenceItem("associative cis candidates", n_assoc, source="cross_modal_graph"),
        EvidenceItem("discordant axes", n_disc, source="cross_modal_graph"),
        EvidenceItem("caveated ambiguous axes", n_cav, source="cross_modal_graph"),
        EvidenceItem("WNN guardrail status", guard_status, source="wnn_diagnostics"),
    ]

    # Attach top non-confounded convergent or candidate axes as typed SemanticFacts.
    axes = graph.get("regulatory_axes") or []
    highlighted: list[str] = []
    for ax in axes:
        if not isinstance(ax, dict):
            continue
        tf = ax.get("tf_regulator")
        gene = ax.get("target_gene")
        cls_name = str(ax.get("classification") or "")
        if tf and cls_name in {
            "convergent_cross_contrast",
            "convergent_paired_multiome",
            "associative_cis_candidate",
        }:
            if tf not in highlighted:
                highlighted.append(str(tf))
                evidence.append(
                    EvidenceItem(
                        f"triangulated regulator {tf}",
                        cls_name,
                        source="cross_modal_graph",
                        facts=[
                            SemanticFact(
                                subject=str(tf),
                                subject_type="tf_regulator",
                                predicate="cross_modal_triangulation",
                                polarity="affirmed",
                                source="cross_modal_graph",
                            )
                        ],
                    )
                )
        elif gene and cls_name in {
            "convergent_cross_contrast",
            "associative_cis_candidate",
        }:
            if gene not in highlighted:
                highlighted.append(str(gene))
                evidence.append(
                    EvidenceItem(
                        f"triangulated target {gene}",
                        cls_name,
                        source="cross_modal_graph",
                        facts=[
                            SemanticFact(
                                subject=str(gene),
                                subject_type="gene",
                                predicate="cross_modal_triangulation",
                                polarity="affirmed",
                                source="cross_modal_graph",
                            )
                        ],
                    )
                )
        if len(highlighted) >= 4:
            break

    caveats: list[Caveat] = [
        Caveat(
            text=(
                "Cross-modal triangulation summarizes observational co-variation "
                "across paired modalities; candidate eRegulons share upstream "
                "peak-to-gene cis links on the same cells and are not counted as "
                "independent cohort replication or experimental proof of "
                "transcriptional causation."
            ),
            severity="info",
        )
    ]
    if n_cav > 0:
        caveats.append(
            Caveat(
                text=(
                    f"{n_cav} axis/factor entry(ies) carry motif-paralog "
                    "sharing, sequencing-depth association, or latent-factor "
                    "confounder flags and remain quarantined as ambiguous."
                ),
                severity="warning",
            )
        )
    if guard_status == "degraded":
        caveats.append(
            Caveat(
                text=(
                    "WNN anti-degradation guardrail status is degraded; "
                    "per-modality RNA and ATAC analyses remain the primary "
                    "reference."
                ),
                severity="warning",
            )
        )

    absent = graph.get("absent_not_tested") or []
    if absent:
        caveats.append(
            Caveat(
                text=(
                    "Unrun upstream contrasts are tracked as absent_not_tested "
                    "and are not interpreted as cross-modal discordance."
                ),
                severity="info",
            )
        )

    confidence = "medium"
    if guard_status == "degraded":
        confidence = "low"
    elif n_conv > 0 and n_stages >= 4 and guard_status == "pass":
        confidence = "high"

    claim = (
        f"Cross-modal evidence triangulation across {n_stages} integration "
        f"stage(s) linked {n_nodes} nodes and {n_edges} associative edges into "
        f"{n_axes} regulatory or latent axes ({n_conv} convergent, "
        f"{n_assoc} associative cis candidates, {n_disc} discordant, and "
        f"{n_cav} caveated ambiguous; WNN guardrail {guard_status})."
    )

    return NarrativeBlock(
        id="multiome.triangulation",
        modality="integration",
        analysis="cross_modal_triangulation",
        block_type="integration",
        title="Cross-modal evidence triangulation",
        status="success",
        confidence=confidence,
        claim=claim,
        evidence=evidence,
        caveats=caveats,
        metrics={
            "n_stages": n_stages,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "n_axes_total": n_axes,
            "n_convergent_total": n_conv,
            "convergent_cross_contrast": n_cross,
            "convergent_paired_multiome": n_paired_conv,
            "associative_cis_candidate": n_assoc,
            "discordant_signal": n_disc,
            "caveated_ambiguous": n_cav,
            "wnn_guardrail_status": guard_status,
            "highlighted_entities": highlighted,
        },
    )
