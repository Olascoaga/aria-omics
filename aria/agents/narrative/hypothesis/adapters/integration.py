"""Paired cross-modal integration evidence adapter (ADR-057 S10; ADR-072 V48 I9).

Pure and deterministic — zero LLM, zero biological priors. Flattens audited
paired 10x Multiome (same-cell RNA + ATAC) ``IntegrationAgent`` findings into
``EvidenceSignal`` items across three associative cross-modal layers plus WNN
anti-degradation guardrail propagation:

  1. ``peak_to_gene`` (``ledger://integration/peak_to_gene``) -> cis-linked
     ``gene`` signals (``peak2gene_correlation``), caveated with
     ``PEAK2GENE_ASSOCIATIVE``. Raw ``chr:start-end`` peak coordinates are never
     emitted as entities.
  2. ``eregulons`` (``ledger://integration/eregulons``) -> candidate regulator
     ``tf_motif`` signals (``eregulon_mean_tf_gene_correlation``, signed by
     observed TF–target correlation) + top pruned ``gene`` membership signals
     (``eregulon_target_membership``, ``value=None``, ``direction="na"`` so
     module-mean correlations are never misattributed to individual targets).
     Inherits ``MOTIF_NOT_BINDING`` + ``PEAK2GENE_ASSOCIATIVE``, plus
     ``PARALOG_MOTIF_AMBIGUITY`` and/or ``DEPTH_CONFOUNDED`` strictly when
     flagged on that module.
  3. ``mofa`` (``ledger://integration/mofa``) -> signed RNA ``gene`` loadings
     (``mofa_shared_rna_loading``, ``direction="na"`` respecting latent factor
     sign indeterminacy ``Z W^T = (-Z)(-W)^T``) ONLY from active ``shared``
     RNA+ATAC factors that are NOT flagged in ``excluded_factors`` or technical/
     cell-cycle confounder fields. Bare unsigned ``factor1_top_features`` lists
     are never used.
  4. ``wnn`` (``ledger://integration/wnn``) -> numeric cluster IDs are never
     emitted as entities; if WNN ran and ``guardrail_status`` is ``warn``,
     ``degraded``, or ``fail`` (ADR-066 "lens, not replacement"), every emitted
     integration signal inherits ``WNN_GUARDRAIL_DEGRADED``.

Non-circularity (ADR-070 + ADR-072):
  ``eregulons`` directly consumes ``peak_to_gene`` cis links on the same shared
  cell barcodes. Both layers carry ``independence_key="paired_cis_regulatory"``
  so competitive ranking counts shared-upstream cis representations as ONE
  lineage line rather than inflating convergence.
"""

from __future__ import annotations

import json
import logging
import math
import re
from pathlib import Path
from typing import Any

from ..caveats import (
    DEPTH_CONFOUNDED,
    MOFA_LATENT_ASSOCIATIVE,
    MOTIF_NOT_BINDING,
    PARALOG_MOTIF_AMBIGUITY,
    PEAK2GENE_ASSOCIATIVE,
    WNN_GUARDRAIL_DEGRADED,
)
from ..types import EvidenceSignal
from ._chromatin_common import node_ran

log = logging.getLogger(__name__)

_WNN_NODE = "ledger://integration/wnn"
_MOFA_NODE = "ledger://integration/mofa"
_P2G_NODE = "ledger://integration/peak_to_gene"
_EREG_NODE = "ledger://integration/eregulons"

PAIRED_CIS_REGULATORY_KEY = "paired_cis_regulatory"
PAIRED_MOFA_LATENT_KEY = "paired_mofa_latent"

_MAX_P2G_GENES = 20
_MAX_EREGULONS = 10
_MAX_TARGETS_PER_EREGULON = 3
_MAX_MOFA_FACTORS = 3
_MAX_GENES_PER_FACTOR = 5

_NON_RUN_STATUSES = frozenset(
    {"not_run", "skipped", "blocked", "error", "failed"}
)
_DEGRADED_GUARDRAIL_STATUSES = frozenset({"warn", "degraded", "fail"})
_PEAK_COORD_RE = re.compile(
    r"^(?:chr[\w.\-]+|[\dXYM]+|MT):\d+[-:]\d+$", re.IGNORECASE
)
_COMPOSITE_TF_MARKERS = ("::", ",", "/", ";", "(", ")", " ")


def _finite_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _finite_corr(value: Any) -> float | None:
    out = _finite_float(value)
    if out is None:
        return None
    if abs(out) > 1.0 + 1e-6:
        return None
    return max(-1.0, min(1.0, out))


def _is_peak_coordinate(token: str) -> bool:
    s = str(token or "").strip()
    if not s:
        return False
    if _PEAK_COORD_RE.match(s):
        return True
    # Also guard Signac-style chr1-10000-10500 coordinate triples.
    parts = s.split("-")
    if (
        len(parts) == 3
        and parts[0].lower().startswith("chr")
        and parts[1].isdigit()
        and parts[2].isdigit()
    ):
        return True
    return False


def _valid_gene_symbol(raw: Any) -> str | None:
    sym = str(raw or "").strip()
    if not sym or sym == "?":
        return None
    if _is_peak_coordinate(sym):
        return None
    if any(ch in sym for ch in (" ", "\t", "\n", ";", ",")):
        return None
    return sym


def _valid_single_tf_symbol(raw: Any) -> str | None:
    """Return a clean single TF gene symbol, rejecting composite motif labels.

    ``eregulon_core.py`` already resolves candidate eRegulons 1-to-1 to expressed
    TF gene symbols. If an unsplit heterodimer/family string (e.g. ``JUN::FOS``)
    ever reaches the adapter, discard it rather than inventing component genes.
    """
    sym = _valid_gene_symbol(raw)
    if sym is None:
        return None
    if any(marker in sym for marker in _COMPOSITE_TF_MARKERS):
        return None
    return sym


def _is_executed_stage(stage: Any) -> bool:
    if not isinstance(stage, dict) or not stage:
        return False
    status = str(stage.get("status") or "").strip().lower()
    if status in _NON_RUN_STATUSES:
        return False
    if stage.get("ran") is False:
        return False
    return True


def _unwrap_integration_findings(agent_results: dict | None) -> dict[str, Any]:
    if not isinstance(agent_results, dict):
        return {}
    integ = agent_results.get("integration_agent")
    if not isinstance(integ, dict) or not integ:
        return {}
    top_status = str(integ.get("status") or "").strip().lower()
    if top_status in _NON_RUN_STATUSES:
        return {}
    findings = integ.get("findings")
    if isinstance(findings, dict):
        return findings
    return integ


def _dedup_codes(codes: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for code in codes:
        if code and code not in seen:
            seen.add(code)
            out.append(code)
    return out


def _wnn_guardrail_caveats(
    findings: dict[str, Any], run_ledger: dict | None
) -> list[str]:
    wnn = findings.get("wnn")
    if not _is_executed_stage(wnn) or not node_ran(run_ledger, _WNN_NODE):
        return []
    guardrail = str(wnn.get("guardrail_status") or "").strip().lower()
    if guardrail in _DEGRADED_GUARDRAIL_STATUSES:
        return [WNN_GUARDRAIL_DEGRADED]
    return []


def _peak_to_gene_signals(
    p2g: dict[str, Any],
    *,
    wnn_caveats: list[str],
    max_genes: int,
) -> list[EvidenceSignal]:
    padj_cutoff = _finite_float(p2g.get("max_padj"))
    if padj_cutoff is None:
        padj_cutoff = 0.05

    candidates: list[tuple[float, float, str, str, str, float, str]] = []
    for row in p2g.get("top_links") or []:
        if not isinstance(row, dict):
            continue
        gene = _valid_gene_symbol(row.get("gene"))
        if gene is None:
            continue
        peak_raw = str(row.get("peak") or "").strip()
        corr_raw = (
            row.get("correlation")
            if row.get("correlation") is not None
            else row.get("corr")
        )
        corr = _finite_corr(corr_raw)
        if corr is None or corr == 0.0:
            continue
        if "padj" in row and row.get("padj") is not None:
            padj = _finite_float(row.get("padj"))
            if padj is None or padj < 0.0 or padj > padj_cutoff + 1e-12:
                continue
        else:
            padj = None

        direction = "up" if corr > 0.0 else "down"
        declared_dir = str(row.get("direction") or "").strip().lower()
        if declared_dir:
            if ("pos" in declared_dir or declared_dir == "up") and direction != "up":
                continue
            if ("neg" in declared_dir or declared_dir == "down") and direction != "down":
                continue

        region_class = str(row.get("region_class") or "cis").strip() or "cis"
        ctx = (
            f"cis:{region_class}:{peak_raw}"
            if peak_raw
            else f"cis:{region_class}"
        )
        candidates.append(
            (
                padj if padj is not None else 0.0,
                -abs(corr),
                gene.lower(),
                peak_raw,
                gene,
                corr,
                ctx,
            )
        )

    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

    signals: list[EvidenceSignal] = []
    seen_gene_dir: set[tuple[str, str]] = set()
    caveats = _dedup_codes([PEAK2GENE_ASSOCIATIVE, *wnn_caveats])
    for _, _, gene_low, _, gene, corr, ctx in candidates:
        direction = "up" if corr > 0.0 else "down"
        key = (gene_low, direction)
        if key in seen_gene_dir:
            continue
        seen_gene_dir.add(key)
        signals.append(
            EvidenceSignal(
                entity=gene,
                entity_kind="gene",
                modality="integration",
                measure="peak2gene_correlation",
                audited_node_ref=_P2G_NODE,
                value=round(corr, 4),
                direction=direction,
                caveats_inherited=list(caveats),
                context=ctx,
                independence_key=PAIRED_CIS_REGULATORY_KEY,
            )
        )
        if len(signals) >= max_genes:
            break
    return signals


def _eregulon_signals(
    ereg: dict[str, Any],
    *,
    wnn_caveats: list[str],
    max_eregulons: int,
    max_targets: int,
) -> list[EvidenceSignal]:
    padj_cutoff = _finite_float(ereg.get("max_padj"))
    if padj_cutoff is None:
        padj_cutoff = 0.05

    ranked_mods: list[tuple[float, float, str, str, dict[str, Any], str, float]] = []
    for mod in ereg.get("top_eregulons") or []:
        if not isinstance(mod, dict):
            continue
        tf = _valid_single_tf_symbol(mod.get("tf"))
        if tf is None:
            continue
        sign_class = str(mod.get("sign_class") or "pos_tf_pos_p2g").strip()
        raw_r = (
            mod.get("mean_r_tf_gene")
            if mod.get("mean_r_tf_gene") is not None
            else mod.get("tf_gene_corr")
        )
        mean_r = _finite_corr(raw_r)
        if mean_r is None or mean_r == 0.0:
            continue
        if sign_class.startswith("pos_") and mean_r <= 0.0:
            continue
        if sign_class.startswith("neg_") and mean_r >= 0.0:
            continue

        raw_padj = (
            mod.get("regulon_padj")
            if mod.get("regulon_padj") is not None
            else mod.get("padj")
        )
        if raw_padj is not None:
            padj = _finite_float(raw_padj)
            if padj is None or padj < 0.0 or padj > padj_cutoff + 1e-12:
                continue
        else:
            padj = None

        ranked_mods.append(
            (
                padj if padj is not None else 0.0,
                -abs(mean_r),
                tf.lower(),
                sign_class,
                mod,
                tf,
                mean_r,
            )
        )

    ranked_mods.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

    signals: list[EvidenceSignal] = []
    seen_tf_ctx: set[tuple[str, str]] = set()
    mods_emitted = 0

    for _, _, tf_low, sign_class, mod, tf, mean_r in ranked_mods:
        if mods_emitted >= max_eregulons:
            break
        regulon_id = str(mod.get("regulon_id") or f"{tf}_{sign_class}").strip()
        ctx = (
            f"eregulon:{regulon_id}"
            if sign_class and sign_class in regulon_id
            else f"eregulon:{regulon_id}:{sign_class}"
        )
        if (tf_low, ctx) in seen_tf_ctx:
            continue
        seen_tf_ctx.add((tf_low, ctx))

        caveats = [MOTIF_NOT_BINDING, PEAK2GENE_ASSOCIATIVE]
        if bool(mod.get("motif_family_ambiguity")) or bool(
            mod.get("ambiguous_paralog_tfs")
        ):
            caveats.append(PARALOG_MOTIF_AMBIGUITY)
        if bool(mod.get("depth_confounded")):
            caveats.append(DEPTH_CONFOUNDED)
        caveats.extend(wnn_caveats)
        mod_caveats = _dedup_codes(caveats)

        tf_dir = "up" if mean_r > 0.0 else "down"
        signals.append(
            EvidenceSignal(
                entity=tf,
                entity_kind="tf_motif",
                modality="integration",
                measure="eregulon_mean_tf_gene_correlation",
                audited_node_ref=_EREG_NODE,
                value=round(mean_r, 4),
                direction=tf_dir,
                caveats_inherited=list(mod_caveats),
                context=ctx,
                independence_key=PAIRED_CIS_REGULATORY_KEY,
            )
        )

        raw_targets = (
            mod.get("top_target_genes")
            if mod.get("top_target_genes") is not None
            else mod.get("target_genes")
        ) or []
        seen_tg: set[str] = set()
        n_tg = 0
        for raw_tg in raw_targets:
            if n_tg >= max_targets:
                break
            tg = _valid_gene_symbol(raw_tg)
            if tg is None or tg.lower() == tf_low or tg.lower() in seen_tg:
                continue
            seen_tg.add(tg.lower())
            signals.append(
                EvidenceSignal(
                    entity=tg,
                    entity_kind="gene",
                    modality="integration",
                    measure="eregulon_target_membership",
                    audited_node_ref=_EREG_NODE,
                    value=None,
                    direction="na",
                    caveats_inherited=list(mod_caveats),
                    context=ctx,
                    independence_key=PAIRED_CIS_REGULATORY_KEY,
                )
            )
            n_tg += 1

        mods_emitted += 1

    return signals


def _resolve_mofa_loadings_map(
    mofa: dict[str, Any],
) -> dict[str, dict[str, list[dict[str, Any]]]]:
    """Resolve signed per-factor, per-view loadings from inline dict or JSON path."""
    for key in ("top_loadings", "factor_loadings", "loadings"):
        inline = mofa.get(key)
        if isinstance(inline, dict) and inline:
            factors = inline.get("factors")
            if isinstance(factors, dict):
                return factors
            return inline

    loadings_path = mofa.get("factor_loadings_path")
    if isinstance(loadings_path, str) and loadings_path.strip():
        p = Path(loadings_path.strip())
        try:
            if p.is_file():
                payload = json.loads(p.read_text(encoding="utf-8"))
                if isinstance(payload, dict):
                    factors = payload.get("factors")
                    if isinstance(factors, dict):
                        return factors
        except Exception as exc:
            log.warning(
                "Could not read MOFA+ factor_loadings_path for hypothesis "
                "adapter (%s); omitting MOFA+ loadings.",
                exc,
            )
    return {}


def _rna_view_keys(
    mofa: dict[str, Any], factor_views: dict[str, Any]
) -> list[str]:
    views_meta = mofa.get("views")
    preferred: list[str] = []
    if isinstance(views_meta, dict):
        for v_name, v_info in views_meta.items():
            if (
                isinstance(v_info, dict)
                and str(v_info.get("kind") or "").strip().lower() == "rna"
                and v_name in factor_views
            ):
                preferred.append(str(v_name))
    if preferred:
        return sorted(preferred)

    matched: list[str] = []
    for v_name in factor_views:
        low = str(v_name).strip().lower()
        if any(bad in low for bad in ("atac", "peak", "chrom")):
            continue
        if "rna" in low or "expr" in low or "gex" in low:
            matched.append(str(v_name))
    return sorted(matched)


def _mofa_excluded_ids(mofa: dict[str, Any]) -> set[str]:
    excluded: set[str] = set()
    for ex in mofa.get("excluded_factors") or []:
        if not isinstance(ex, dict):
            continue
        fid = (
            ex.get("factor_id")
            if ex.get("factor_id") is not None
            else ex.get("factor")
        )
        if fid is not None:
            excluded.add(str(fid).strip())
    if bool(mofa.get("technical_factor_flag")) and (
        mofa.get("technical_factor_id") is not None
    ):
        excluded.add(str(mofa.get("technical_factor_id")).strip())
    for cc_id in mofa.get("cell_cycle_factor_ids") or []:
        if cc_id is not None:
            excluded.add(str(cc_id).strip())
    return excluded


def _factor_sort_id(fid_str: str) -> tuple[int, int | str]:
    try:
        return (0, int(fid_str))
    except ValueError:
        return (1, fid_str)


def _mofa_signals(
    mofa: dict[str, Any],
    *,
    wnn_caveats: list[str],
    max_factors: int,
    max_genes_per_factor: int,
) -> list[EvidenceSignal]:
    factors_map = _resolve_mofa_loadings_map(mofa)
    if not factors_map:
        return []

    excluded_ids = _mofa_excluded_ids(mofa)
    eligible: list[tuple[float, tuple[int, int | str], str]] = []
    seen_fids: set[str] = set()

    for fac in mofa.get("top_factors") or []:
        if not isinstance(fac, dict):
            continue
        fid_raw = (
            fac.get("factor_id")
            if fac.get("factor_id") is not None
            else fac.get("factor")
        )
        if fid_raw is None:
            continue
        fid_str = str(fid_raw).strip()
        if not fid_str or fid_str in seen_fids:
            continue
        seen_fids.add(fid_str)

        # Exclusion always wins on contradiction.
        if (
            fid_str in excluded_ids
            or fac.get("active") is not True
            or bool(fac.get("excluded_from_biological_interpretation"))
            or bool(fac.get("technical_depth_flag"))
            or bool(fac.get("cell_cycle_flag"))
            or bool(fac.get("exclusion_reasons"))
        ):
            continue

        cls_name = str(
            fac.get("category") or fac.get("classification") or ""
        ).strip().lower()
        if cls_name != "shared":
            continue

        r2_view = (
            fac.get("r2_by_view")
            if isinstance(fac.get("r2_by_view"), dict)
            else {}
        )
        total_r2 = 0.0
        for v_val in r2_view.values():
            fv = _finite_float(v_val)
            if fv is not None and fv > 0.0:
                total_r2 += fv
        if total_r2 == 0.0:
            for rk in ("rna_r2", "atac_r2", "r2_rna", "r2_atac", "r2_scRNA", "r2_scATAC"):
                fv = _finite_float(fac.get(rk))
                if fv is not None and fv > 0.0:
                    total_r2 += fv

        eligible.append((-total_r2, _factor_sort_id(fid_str), fid_str))

    eligible.sort(key=lambda item: (item[0], item[1]))

    signals: list[EvidenceSignal] = []
    caveats = _dedup_codes([MOFA_LATENT_ASSOCIATIVE, *wnn_caveats])
    factors_emitted = 0

    for _, _, fid_str in eligible:
        if factors_emitted >= max_factors:
            break
        factor_views = (
            factors_map.get(f"Factor{fid_str}")
            or factors_map.get(f"Factor_{fid_str}")
            or factors_map.get(fid_str)
        )
        if not isinstance(factor_views, dict):
            continue
        rna_keys = _rna_view_keys(mofa, factor_views)
        if not rna_keys:
            continue

        gene_rows: list[tuple[float, str, str, float]] = []
        for rna_key in rna_keys:
            rows = factor_views.get(rna_key)
            if not isinstance(rows, list):
                continue
            for row in rows:
                if not isinstance(row, dict):
                    continue
                gene = _valid_gene_symbol(
                    row.get("feature") or row.get("gene") or row.get("symbol")
                )
                if gene is None:
                    continue
                weight = _finite_float(row.get("weight"))
                if weight is None or weight == 0.0:
                    continue
                gene_rows.append((-abs(weight), gene.lower(), gene, weight))

        if not gene_rows:
            continue

        gene_rows.sort(key=lambda item: (item[0], item[1]))
        seen_genes: set[str] = set()
        n_genes = 0
        for _, gene_low, gene, weight in gene_rows:
            if n_genes >= max_genes_per_factor:
                break
            if gene_low in seen_genes:
                continue
            seen_genes.add(gene_low)
            pole = "positive_loading" if weight > 0.0 else "negative_loading"
            signals.append(
                EvidenceSignal(
                    entity=gene,
                    entity_kind="gene",
                    modality="integration",
                    measure="mofa_shared_rna_loading",
                    audited_node_ref=_MOFA_NODE,
                    value=round(weight, 4),
                    direction="na",
                    caveats_inherited=list(caveats),
                    context=f"Factor{fid_str}:rna:{pole}",
                    independence_key=PAIRED_MOFA_LATENT_KEY,
                )
            )
            n_genes += 1

        if n_genes > 0:
            factors_emitted += 1

    return signals


def integration_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
    *,
    max_p2g_genes: int = _MAX_P2G_GENES,
    max_eregulons: int = _MAX_EREGULONS,
    max_targets_per_eregulon: int = _MAX_TARGETS_PER_EREGULON,
    max_mofa_factors: int = _MAX_MOFA_FACTORS,
    max_genes_per_factor: int = _MAX_GENES_PER_FACTOR,
) -> list[EvidenceSignal]:
    """Flatten audited paired multiome ``IntegrationAgent`` findings into ``EvidenceSignal`` items."""
    del exp_ctx  # Reserved for uniform adapter signature; never injects priors.
    findings = _unwrap_integration_findings(agent_results)
    if not findings:
        return []

    wnn_caveats = _wnn_guardrail_caveats(findings, run_ledger)
    signals: list[EvidenceSignal] = []

    p2g = findings.get("peak_to_gene")
    if _is_executed_stage(p2g) and node_ran(run_ledger, _P2G_NODE):
        signals += _peak_to_gene_signals(
            p2g,
            wnn_caveats=wnn_caveats,
            max_genes=max_p2g_genes,
        )

    ereg = findings.get("eregulons")
    if _is_executed_stage(ereg) and node_ran(run_ledger, _EREG_NODE):
        signals += _eregulon_signals(
            ereg,
            wnn_caveats=wnn_caveats,
            max_eregulons=max_eregulons,
            max_targets=max_targets_per_eregulon,
        )

    mofa = findings.get("mofa")
    if _is_executed_stage(mofa) and node_ran(run_ledger, _MOFA_NODE):
        signals += _mofa_signals(
            mofa,
            wnn_caveats=wnn_caveats,
            max_factors=max_mofa_factors,
            max_genes_per_factor=max_genes_per_factor,
        )

    return signals
