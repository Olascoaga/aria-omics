"""Centralized degradation and warning router for modality narrators (W-WARN).

Ensures that warnings, stage errors, partial database/sample failures,
``*_degraded`` flags, ``metrics_not_computed`` entries, ``figures._errors``,
and ``lfc_shrinkage`` requested-vs-applied discrepancies emitted by compute
scripts always reach the governed ``NarrativeBlock`` stream before claim
compilation (C1-C4) and HTML/JSON rendering.

Pure and dependency-light: zero LLM calls, zero I/O, strictly idempotent.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import wraps
import re
from typing import Any

from aria.agents.narrative.types import (
    PUBLISHABLE_STATUSES,
    Caveat,
    CaveatSeverity,
    NarrativeBlock,
)

_TAG_STRIP_RE = re.compile(r"<[^>]+>")
_CONTRAST_PREFIX_RE = re.compile(r"^\[([^\]]+)\]\s*(.+)$")
_LOW_COUNT_INFO_RE = re.compile(r"^\d[\d,]*\s+low-count genes removed\b", re.IGNORECASE)

_INFO_WARNING_PREFIXES: tuple[str, ...] = (
    "gene annotation loaded from ",
    "dr input: ",
    "tpm table written ",
    "tpm export suppressed for protocol ",
    "technical libraries were summed within biological units",
    "scatac peak-matrix qc on .h5ad:",
    "qc metrics not computed:",
    "no significant ",
)

_DEGRADED_NOT_COMPUTED_MARKERS: tuple[str, ...] = (
    "error", "fail", "corrupt", "unparseable", "no finite numeric", "invalid", "crash",
)

_STAGE_TO_ANALYSES: dict[str, tuple[str, ...]] = {
    "preprocessing": ("preprocessing", "scatac_preprocessing", "bulk_atac_preprocessing"),
    "inspection": ("qc",),
    "balancing": ("qc",),
    "qc": ("qc",),
    "sample_qc": ("sample_qc",),
    "outlier_sensitivity": ("sample_qc", "differential_expression"),
    "peaks": ("peak_calling",),
    "peak_counts": ("peak_count_matrix",),
    "lsi": ("dimensionality_reduction",),
    "lsi_clustering": ("dimensionality_reduction",),
    "differential_accessibility": ("differential_accessibility",),
    "diffacc": ("differential_accessibility",),
    "bulk_diffacc": ("differential_accessibility",),
    "motifs": ("motif_enrichment",),
    "peak_annotation": ("peak_annotation",),
    "peak_ora": ("peak_ora",),
    "regulatory": ("regulatory_layers",),
    "footprinting": ("differential_tf_footprinting", "footprinting"),
    "compartments": ("compartments",),
    "tads": ("tads",),
    "loops": ("loops",),
    "differential_3d": ("differential_3d",),
    "paired_qc": ("paired_qc",),
    "wnn": ("wnn",),
    "mofa": ("mofa",),
    "peak_to_gene": ("peak_to_gene",),
    "eregulons": ("eregulons",),
    "pseudobulk_de": ("pseudobulk_de",),
    "pseudobulk_pathways": ("pseudobulk_pathways", "pathway_enrichment"),
    "pathways": ("pathway_enrichment", "pseudobulk_pathways"),
    "pathway_ora": ("pathway_enrichment",),
    "contrast": ("differential_expression", "contrast"),
    "plots": ("differential_expression", "sample_qc"),
    "figures": ("qc", "sample_qc", "dimensionality_reduction", "differential_accessibility", "peak_calling"),
    "integration": ("qc", "data_quality"),
    "clustering": ("qc", "marker_discovery"),
    "cell_types": ("qc", "data_quality"),
    "trajectory": ("trajectory",),
    "cell_communication": ("cell_communication",),
}

_FLAT_STAGE_KEYS: tuple[str, ...] = (
    "preprocessing", "inspection", "balancing", "qc", "sample_qc",
    "outlier_sensitivity", "peaks", "peak_counts", "lsi", "lsi_clustering",
    "differential_accessibility", "diffacc", "bulk_diffacc", "motifs",
    "peak_annotation", "peak_ora", "regulatory", "footprinting",
    "compartments", "tads", "loops", "differential_3d",
    "paired_qc", "wnn", "mofa", "peak_to_gene", "eregulons",
    "pseudobulk_de", "pseudobulk_pathways", "pathways",
    "integration", "clustering", "cell_types", "trajectory",
    "cell_communication", "figures",
)

_PRIMARY_ERROR_SCOPES: frozenset[str] = frozenset({
    "qc", "sample_qc", "paired_qc", "wnn", "mofa", "peak_to_gene", "eregulons",
    "compartments", "tads", "loops", "differential_3d", "peaks", "peak_counts",
    "differential_accessibility", "diffacc", "bulk_diffacc", "pseudobulk_de",
    "contrast", "motifs",
})

_AGENT_META: dict[str, tuple[str, str]] = {
    "bulk_rna_agent": ("bulk", "bulk RNA-seq"),
    "scrna_agent": ("scrna", "scRNA-seq"),
    "rna_agent": ("scrna", "scRNA-seq"),
    "chromatin_agent": ("chromatin", "chromatin"),
    "genome_arch_agent": ("hic", "hic"),
    "integration_agent": ("multiome", "integration"),
}

_CHROMATIN_ASSAYS: tuple[str, ...] = (
    "scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG",
)


@dataclass(frozen=True)
class _DegradationSignal:
    scope: str
    target_analyses: tuple[str, ...]
    contrast_hint: str | None
    modality_hint: str | None
    field: str
    severity: CaveatSeverity
    message: str
    synthesize_if_missing: bool = False
    broadcast_all_targets: bool = False


def _safe_id(value: Any) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "unknown"))
    return text.strip("_") or "unknown"


def _clean_text(value: Any, max_len: int = 300) -> str:
    txt = _TAG_STRIP_RE.sub("", str(value or ""))
    txt = re.sub(r"[\r\n\t]+", " ", txt).strip()
    return re.sub(r" {2,}", " ", txt)[:max_len]


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip().lower())


def _classify_warning_severity(msg: str) -> CaveatSeverity:
    low = msg.lower().strip()
    if _LOW_COUNT_INFO_RE.match(low):
        return "info"
    if low.startswith("pca flagged ") and "likely real biological variation" in low:
        return "info"
    if any(low.startswith(p) for p in _INFO_WARNING_PREFIXES):
        return "info"
    return "warning"


def _infer_root_targets(msg: str, contrast_hint: str | None) -> tuple[tuple[str, ...], bool]:
    low = msg.lower()
    if any(k in low for k in ("enrichr", "pathway", "ora ", "ora:", "gene set", "gmt")):
        return ("pathway_enrichment", "pseudobulk_pathways"), True
    if any(k in low for k in (
        "vst", "rlog", "pca", "mds", "sample qc", "gtf", "gene annotation",
        "tpm", "cpm", "outlier", "low-count", "technical libraries",
    )):
        return ("sample_qc", "qc"), False
    if contrast_hint:
        return ("differential_expression", "contrast"), False
    return ("sample_qc", "qc", "differential_expression", "pseudobulk_de"), False


def _warnings_signals(
    scope: str, payload: dict[str, Any], default_targets: tuple[str, ...],
    c_hint: str | None, m_hint: str | None,
) -> list[_DegradationSignal]:
    raw = payload.get("warnings")
    if not isinstance(raw, list):
        return []
    out: list[_DegradationSignal] = []
    for raw_w in raw:
        msg = _clean_text(raw_w) if raw_w else ""
        if not msg:
            continue
        cur_c, body = c_hint, msg
        m = _CONTRAST_PREFIX_RE.match(msg)
        if m and not cur_c:
            cur_c, body = m.group(1).strip(), m.group(2).strip()
        sev = _classify_warning_severity(body)
        if scope == "__root__":
            targets, synth = _infer_root_targets(body, cur_c)
            synth = synth and (sev != "info")
            broadcast = targets == ("pathway_enrichment", "pseudobulk_pathways")
        else:
            targets, synth = default_targets, False
            broadcast = scope in {"contrast", "pathway_ora", "pseudobulk_de"}
        out.append(_DegradationSignal(
            scope, targets, cur_c, m_hint, "warnings", sev, msg,
            synthesize_if_missing=synth, broadcast_all_targets=broadcast,
        ))
    return out


def _error_and_degraded_signals(
    scope: str, payload: dict[str, Any], default_targets: tuple[str, ...],
    c_hint: str | None, m_hint: str | None,
) -> list[_DegradationSignal]:
    out: list[_DegradationSignal] = []
    err_val = payload.get("error")
    if isinstance(err_val, str) and err_val.strip():
        sev: CaveatSeverity = "blocking" if scope in _PRIMARY_ERROR_SCOPES else "warning"
        out.append(_DegradationSignal(
            scope, default_targets, c_hint, m_hint, "error", sev, _clean_text(err_val),
            synthesize_if_missing=(scope != "__root__"),
        ))
    for k, v in payload.items():
        if k == "error" or not isinstance(k, str):
            continue
        if k.endswith("_error") and isinstance(v, str) and v.strip():
            out.append(_DegradationSignal(
                scope, default_targets, c_hint, m_hint, k, "warning", _clean_text(f"{k}: {v}"),
            ))
        elif k == "_errors" and isinstance(v, list):
            out.extend(
                _DegradationSignal(
                    scope, default_targets, c_hint, m_hint, "_errors", "warning",
                    _clean_text(f"Figure rendering error: {item}"),
                )
                for item in v if item
            )
        elif k.endswith("_degraded") and v is True:
            base = k[:-len("_degraded")]
            ctx_val = payload.get(f"{base}_source") or payload.get(f"{base}_reason")
            desc = f"{k}=True ({ctx_val})" if ctx_val else f"{k}=True"
            out.append(_DegradationSignal(
                scope, default_targets, c_hint, m_hint, k, "warning",
                _clean_text(f"Degraded computation: {desc}."),
                broadcast_all_targets=scope in {"pseudobulk_de", "contrast"},
            ))
    return out


def _partial_failure_signals(
    scope: str, payload: dict[str, Any], default_targets: tuple[str, ...],
    c_hint: str | None, m_hint: str | None,
) -> list[_DegradationSignal]:
    out: list[_DegradationSignal] = []
    targets = default_targets or ("pathway_enrichment",)
    for fail_key, label in (
        ("failed_databases", "Pathway ORA database"),
        ("failed_samples", "Cohort sample"),
        ("failed_modalities", "Modality"),
    ):
        obj = payload.get(fail_key)
        if isinstance(obj, dict):
            items = [f"{label} {name} failed: {reason}." for name, reason in obj.items()]
        elif isinstance(obj, list):
            items = [f"{label} failed: {item}." for item in obj if item]
        else:
            continue
        out.extend(
            _DegradationSignal(
                scope, targets, c_hint, m_hint, fail_key, "warning", _clean_text(txt),
                synthesize_if_missing=True, broadcast_all_targets=True,
            )
            for txt in items
        )
    return out


def _not_computed_signals(
    scope: str, payload: dict[str, Any], default_targets: tuple[str, ...],
    c_hint: str | None, m_hint: str | None,
) -> list[_DegradationSignal]:
    mnc = payload.get("metrics_not_computed")
    if not isinstance(mnc, list):
        return []
    out: list[_DegradationSignal] = []
    for item in mnc:
        txt = _clean_text(item) if item else ""
        if not txt:
            continue
        low = txt.lower()
        sev: CaveatSeverity = (
            "warning" if any(tok in low for tok in _DEGRADED_NOT_COMPUTED_MARKERS) else "info"
        )
        out.append(_DegradationSignal(
            scope, default_targets, c_hint, m_hint, "metrics_not_computed", sev,
            f"Metric not computed: {txt}",
        ))
    return out


def _shrinkage_signals(
    scope: str, payload: dict[str, Any], default_targets: tuple[str, ...],
    c_hint: str | None, m_hint: str | None,
) -> list[_DegradationSignal]:
    shrink = payload.get("lfc_shrinkage")
    if not isinstance(shrink, dict):
        return []
    requested = shrink.get("requested")
    applied = shrink.get("applied")
    n_failed = shrink.get("n_failed")
    n_applied = shrink.get("n_applied")
    reason = shrink.get("reason")
    targets = default_targets or ("pseudobulk_de", "differential_expression")

    if requested is True and (
        applied is False
        or (isinstance(n_failed, int) and n_failed > 0)
        or (isinstance(n_applied, int) and n_applied == 0 and n_failed is not None)
    ):
        failed_reasons = shrink.get("failed_reasons") or []
        detail = (
            "; ".join(str(r) for r in failed_reasons[:3])
            if failed_reasons
            else str(reason or "solver failed or reverted to raw MLE")
        )
        return [_DegradationSignal(
            scope, targets, c_hint, m_hint, "lfc_shrinkage", "warning",
            _clean_text(
                f"apeGLM LFC shrinkage did not apply to all comparisons ({detail}); "
                "affected log2FC values are unshrunken MLE estimates."
            ),
            broadcast_all_targets=True,
        )]
    if requested is None and applied is False and reason and str(reason).lower() != "disabled":
        return [_DegradationSignal(
            scope, targets, c_hint, m_hint, "lfc_shrinkage", "warning",
            _clean_text(
                f"apeGLM LFC shrinkage was not applied ({reason}); "
                "reported log2FC is the unshrunken MLE."
            ),
            broadcast_all_targets=True,
        )]
    return []


def _extract_payload_signals(
    scope: str,
    payload: dict[str, Any],
    *,
    contrast_hint: str | None = None,
    modality_hint: str | None = None,
) -> list[_DegradationSignal]:
    targets = _STAGE_TO_ANALYSES.get(scope, ())
    return (
        _warnings_signals(scope, payload, targets, contrast_hint, modality_hint)
        + _error_and_degraded_signals(scope, payload, targets, contrast_hint, modality_hint)
        + _partial_failure_signals(scope, payload, targets, contrast_hint, modality_hint)
        + _not_computed_signals(scope, payload, targets, contrast_hint, modality_hint)
        + _shrinkage_signals(scope, payload, targets, contrast_hint, modality_hint)
    )


def _walk_flat_scopes(
    container: dict[str, Any], mod_hint: str | None = None
) -> list[tuple[str, dict[str, Any], str | None, str | None]]:
    out: list[tuple[str, dict[str, Any], str | None, str | None]] = [
        ("__root__", container, None, mod_hint)
    ]
    for key in _FLAT_STAGE_KEYS:
        sub = container.get(key)
        if isinstance(sub, dict):
            out.append((key, sub, None, mod_hint))
    topo = container.get("topology")
    if isinstance(topo, dict):
        for tkey in ("compartments", "tads", "loops", "differential_3d"):
            sub = topo.get(tkey)
            if isinstance(sub, dict):
                out.append((tkey, sub, None, mod_hint))
    for c in container.get("contrasts") or []:
        if not isinstance(c, dict):
            continue
        cname = str(c.get("name") or "") or None
        out.append(("contrast", c, cname, mod_hint))
        for sub_key in ("pathway_ora", "plots"):
            if isinstance(c.get(sub_key), dict):
                out.append((sub_key, c[sub_key], cname, mod_hint))
    return out


def _iter_agent_scopes(
    agent_name: str, agent_result: dict[str, Any]
) -> list[tuple[str, dict[str, Any], str | None, str | None]]:
    if not isinstance(agent_result, dict):
        return []
    root = agent_result.get("findings")
    if not isinstance(root, dict):
        root = agent_result
    if agent_name in {"scrna_agent", "rna_agent"}:
        sc_wrap = root.get("scRNA")
        inner = sc_wrap.get("findings") if isinstance(sc_wrap, dict) else None
        return _walk_flat_scopes(inner if isinstance(inner, dict) and inner else root)
    if agent_name == "chromatin_agent":
        nested: list[tuple[str, dict[str, Any], str | None, str | None]] = []
        for assay in _CHROMATIN_ASSAYS:
            entry = root.get(assay)
            if isinstance(entry, dict):
                inner = entry.get("findings") if isinstance(entry.get("findings"), dict) else entry
                nested.extend(_walk_flat_scopes(inner, mod_hint=assay))
        return nested or _walk_flat_scopes(root)
    return _walk_flat_scopes(root)


def _select_target_blocks(
    sig: _DegradationSignal, blocks: list[NarrativeBlock]
) -> list[NarrativeBlock]:
    candidates = list(blocks)
    if sig.modality_hint:
        mod_tok = _safe_id(sig.modality_hint).lower()
        mod_matches = [
            b for b in candidates
            if mod_tok in b.id.lower()
            or (sig.modality_hint == "bulk_ATAC" and "bulk_atac" in b.id.lower())
        ]
        if mod_matches:
            candidates = mod_matches

    by_analysis = (
        [b for b in candidates if b.analysis in sig.target_analyses]
        if sig.target_analyses else []
    )
    if sig.contrast_hint and by_analysis:
        ctok = _safe_id(sig.contrast_hint).lower()
        by_contrast = [b for b in by_analysis if ctok in b.id.lower()]
        if by_contrast:
            return by_contrast
        if sig.scope == "pathway_ora":
            return []

    if by_analysis:
        return by_analysis if sig.broadcast_all_targets else [by_analysis[0]]
    if sig.scope == "figures" and candidates:
        return [candidates[0]]
    return []


def _already_present(block: NarrativeBlock, msg: str) -> bool:
    norm_msg = _norm(msg)
    if not norm_msg:
        return True
    for cav in block.caveats or []:
        norm_c = _norm(cav.text)
        if norm_msg == norm_c or norm_msg in norm_c or norm_c in norm_msg:
            return True
    return any(
        isinstance(rec, dict) and _norm(str(rec.get("message") or "")) == norm_msg
        for rec in (block.metadata.get("degradation_signals") or [])
    )


def _apply_signal_to_block(block: NarrativeBlock, sig: _DegradationSignal) -> None:
    if not _already_present(block, sig.message):
        block.caveats.append(Caveat(text=sig.message, severity=sig.severity))
    records = block.metadata.setdefault("degradation_signals", [])
    if not any(
        isinstance(r, dict) and _norm(str(r.get("message") or "")) == _norm(sig.message)
        for r in records
    ):
        records.append({
            "scope": sig.scope,
            "field": sig.field,
            "severity": sig.severity,
            "message": sig.message,
        })

    if sig.severity == "blocking" and block.status in PUBLISHABLE_STATUSES:
        block.status = "error"
        block.block_type = "error"
        block.confidence = "insufficient"
        block.claim = ""
        block.error = block.error or sig.message
    elif sig.severity == "warning":
        if block.status == "success":
            block.status = "warnings"
        if block.confidence == "high":
            block.confidence = "medium"


def _synthesize_orphan_block(
    synth_id: str, modality_label: str, analysis: str, sigs: list[_DegradationSignal]
) -> NarrativeBlock:
    has_blocking = any(s.severity == "blocking" for s in sigs)
    caveats: list[Caveat] = []
    seen_norm: set[str] = set()
    meta_records: list[dict[str, str]] = []
    for s in sigs:
        nm = _norm(s.message)
        if nm in seen_norm:
            continue
        seen_norm.add(nm)
        caveats.append(Caveat(text=s.message, severity=s.severity))
        meta_records.append({
            "scope": s.scope, "field": s.field, "severity": s.severity, "message": s.message,
        })
    return NarrativeBlock(
        id=synth_id,
        modality=modality_label,
        analysis=analysis,
        block_type="error" if has_blocking else "limitation",
        title=f"{modality_label} {analysis.replace('_', ' ')} degradation notice",
        status="error" if has_blocking else "warnings",
        confidence="insufficient",
        claim="",
        error=next((s.message for s in sigs if s.severity == "blocking"), None) if has_blocking else None,
        caveats=caveats,
        metadata={"degradation_signals": meta_records, "degradation_routed": True},
    )


def route_degradation_signals(
    agent_name: str,
    agent_result: dict[str, Any],
    blocks: list[NarrativeBlock],
) -> list[NarrativeBlock]:
    """Attach degradation/warning signals from ``agent_result`` onto ``blocks``.

    Idempotent: calling this function repeatedly on the same blocks produces
    identical caveats, status, confidence, and synthesized degradation blocks.
    """
    out = list(blocks or [])
    signals = [
        sig
        for scope, payload, c_hint, m_hint in _iter_agent_scopes(agent_name, agent_result)
        for sig in _extract_payload_signals(
            scope, payload, contrast_hint=c_hint, modality_hint=m_hint
        )
    ]
    prefix, modality_label = _AGENT_META.get(
        agent_name, (_safe_id(agent_name), agent_name)
    )
    orphans: dict[tuple[str, str], list[_DegradationSignal]] = {}

    for sig in signals:
        targets = _select_target_blocks(sig, out)
        if targets:
            for blk in targets:
                _apply_signal_to_block(blk, sig)
        elif sig.severity in {"warning", "blocking"} and (sig.synthesize_if_missing or not out):
            analysis = sig.target_analyses[0] if sig.target_analyses else sig.scope
            group_tok = _safe_id(f"{analysis}_{sig.contrast_hint}" if sig.contrast_hint else analysis)
            orphans.setdefault((analysis, group_tok), []).append(sig)

    for (analysis, group_tok), sigs in orphans.items():
        synth_id = f"{prefix}.degradation.{group_tok}"
        existing = next((b for b in out if b.id == synth_id), None)
        if existing is not None:
            for sig in sigs:
                _apply_signal_to_block(existing, sig)
        else:
            out.append(_synthesize_orphan_block(synth_id, modality_label, analysis, sigs))

    for blk in out:
        blk.metadata["degradation_routed"] = True
    return out


def routed_collect(fn):
    """Decorator ensuring direct ``narrator.collect(...)`` calls run the shared router."""

    @wraps(fn)
    def _wrapped(self, agent_name: str, agent_result: dict, context: dict | None = None):
        blocks = fn(self, agent_name, agent_result, context)
        return route_degradation_signals(agent_name, agent_result or {}, blocks)

    return _wrapped


def reconcile_methods_claims(
    agent_name: str,
    agent_result: dict[str, Any],
    method_lines: list[str],
) -> list[str]:
    """Ensure Methods text never omits realized degradation flags."""
    lines = list(method_lines or [])
    if not isinstance(agent_result, dict):
        return lines
    root = agent_result.get("findings")
    if not isinstance(root, dict):
        root = agent_result
    if agent_name in {"scrna_agent", "rna_agent"}:
        sc_wrap = root.get("scRNA")
        inner = sc_wrap.get("findings") if isinstance(sc_wrap, dict) else None
        sc_f = inner if isinstance(inner, dict) and inner else root
        pb = sc_f.get("pseudobulk_de") or {}
        if isinstance(pb, dict) and pb.get("background_degraded") is True:
            note = (
                "Note: pseudobulk ORA background fell back to all dataset genes "
                "(background_degraded=True) rather than strictly expressed genes."
            )
            if not any("background_degraded" in ln for ln in lines):
                lines.append(note)
    return lines
