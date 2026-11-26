"""Centralized degradation and warning router for modality narrators (W-WARN).

Pure post-collection pass called solely from ``NarrativeRegistry.collect_blocks``:
structured fields cap confidence / mark failed stages (A3) / synthesize orphans;
free-text ``warnings`` route by structural position only with ``severity="warning"``
(A2) and never change confidence, status, or block count. Pure: stdlib + types only.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any

from aria.agents.narrative.types import Caveat, CaveatSeverity, NarrativeBlock

_TAG_RE = re.compile(r"<[^>]+>")
_WS_RE = re.compile(r"\s+")
_SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")
_FAILED_STATUSES = frozenset({"error", "failed"})
_AUXILIARY_SCOPES = frozenset({"plots", "figures", "lobo_jackknife", "outlier_sensitivity"})
_QC_SUMMARY_ANALYSES = frozenset({"qc", "sample_qc", "paired_qc", "summary"})
_CHROMATIN_ASSAYS = ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG")
_SEV_RANK = {"info": 0, "warning": 1, "blocking": 2}
_AGENT_META = {
    "bulk_rna_agent": ("bulk", "bulk RNA-seq"), "scrna_agent": ("scrna", "scRNA-seq"),
    "rna_agent": ("scrna", "scRNA-seq"), "chromatin_agent": ("chromatin", "chromatin"),
    "genome_arch_agent": ("hic", "hic"), "integration_agent": ("multiome", "integration"),
}
_STAGE_TO_ANALYSES: dict[str, tuple[str, ...]] = {
    "fastq_alignment": ("fastq_alignment",),
    "preprocessing": ("preprocessing", "scatac_preprocessing", "bulk_atac_preprocessing"),
    "fastq_to_fragments": ("scatac_preprocessing",), "fragments_to_matrix": ("scatac_preprocessing",),
    "cohort_matrix_union": ("scatac_preprocessing",), "inspection": ("qc",), "balancing": ("qc",),
    "qc": ("qc",), "sample_qc": ("sample_qc",), "outlier_sensitivity": ("sample_qc", "differential_expression"),
    "peaks": ("peak_calling",), "peak_counts": ("peak_count_matrix",),
    "lsi": ("dimensionality_reduction",), "lsi_clustering": ("dimensionality_reduction",),
    "differential_accessibility": ("differential_accessibility",), "diffacc": ("differential_accessibility",),
    "bulk_diffacc": ("differential_accessibility",), "motifs": ("motif_enrichment",),
    "peak_annotation": ("peak_annotation",), "peak_ora": ("peak_ora",),
    "regulatory": ("regulatory_layers",), "footprint": ("differential_tf_footprinting", "footprinting"),
    "footprinting": ("differential_tf_footprinting", "footprinting"), "compartments": ("compartments",),
    "tads": ("tads",), "loops": ("loops",), "differential_3d": ("differential_3d",),
    "paired_qc": ("paired_qc",), "wnn": ("wnn",), "mofa": ("mofa",),
    "peak_to_gene": ("peak_to_gene",), "eregulons": ("eregulons",),
    "differential_expression": ("marker_discovery", "differential_expression"),
    "differential_abundance": ("differential_abundance",), "pseudobulk_de": ("pseudobulk_de",),
    "pseudobulk_pathways": ("pseudobulk_pathways", "pathway_enrichment"),
    "pathways": ("pathway_enrichment", "pseudobulk_pathways"), "pathway_ora": ("pathway_enrichment",),
    "contrast": ("differential_expression",), "plots": ("differential_expression", "sample_qc"),
    "lobo_jackknife": ("differential_expression", "pseudobulk_de", "differential_accessibility"),
    "figures": ("qc", "sample_qc", "paired_qc", "dimensionality_reduction", "differential_accessibility", "peak_calling"),
    "integration": ("qc", "data_quality"), "clustering": ("qc", "marker_discovery"),
    "cell_types": ("qc", "data_quality"), "trajectory": ("trajectory",), "cell_communication": ("cell_communication",),
}
_CHROM_KEYS = (
    "fastq_to_fragments", "fragments_to_matrix", "cohort_matrix_union", "preprocessing", "qc",
    "peaks", "peak_counts", "lsi", "lsi_clustering", "differential_accessibility", "diffacc",
    "bulk_diffacc", "motifs", "peak_annotation", "peak_ora", "footprint", "footprinting", "regulatory", "figures",
)
_SCRNA_KEYS = (
    "qc", "integration", "clustering", "cell_types", "differential_expression", "differential_abundance",
    "pseudobulk_de", "pseudobulk_pathways", "pathways", "trajectory", "cell_communication", "figures",
)


@dataclass(frozen=True)
class _Signal:
    scope: str
    field: str
    severity: CaveatSeverity
    message: str
    structured: bool
    id_hint: str | None = None
    dedup_tokens: tuple[str, ...] = ()


def _safe_id(val: Any) -> str:
    return _SAFE_ID_RE.sub("_", str(val or "unknown")).strip("_") or "unknown"

def _sanitize(text: Any, max_len: int = 300) -> str:
    cleaned = _WS_RE.sub(" ", _TAG_RE.sub("", str(text or ""))).strip()
    return cleaned[:max_len] if len(cleaned) > max_len else cleaned

def _norm(text: str) -> str:
    return _WS_RE.sub(" ", text).strip().lower()


def _iter_stage_dicts(agent_name: str, agent_result: dict) -> list[tuple[str, dict[str, Any], str | None]]:
    if not isinstance(agent_result, dict):
        return []
    root = agent_result.get("findings") if isinstance(agent_result.get("findings"), dict) else agent_result
    out: list[tuple[str, dict[str, Any], str | None]] = []

    def _add(container: dict[str, Any], keys: tuple[str, ...], hint: str | None = None) -> None:
        out.extend((k, container[k], hint) for k in keys if isinstance(container.get(k), dict) and container[k])

    if root is not agent_result and isinstance(agent_result.get("warnings"), (list, tuple)):
        out.append(("__root__", {"warnings": agent_result["warnings"]}, None))

    if agent_name == "bulk_rna_agent":
        out.append(("__root__", root, None))
        _add(root, ("preprocessing", "sample_qc", "outlier_sensitivity", "figures"))
        for c in (x for x in (root.get("contrasts") or []) if isinstance(x, dict)):
            cname = _safe_id(c.get("name", "contrast"))
            cid = f"bulk.contrast.{cname}"
            out.append(("contrast", c, cid))
            for sub_key, tid in (("pathway_ora", f"bulk.pathway.{cname}"), ("plots", cid), ("lobo_jackknife", cid)):
                if isinstance(c.get(sub_key), dict) and c[sub_key]:
                    out.append((sub_key, c[sub_key], tid))
    elif agent_name in {"scrna_agent", "rna_agent"}:
        sc = root.get("scRNA") if isinstance(root.get("scRNA"), dict) else None
        inner = (sc.get("findings") if isinstance(sc.get("findings"), dict) else sc) if sc else root
        if inner is not root and isinstance(root.get("warnings"), (list, tuple)):
            out.append(("__root__", {"warnings": root["warnings"]}, None))
        out.append(("__root__", inner, None))
        _add(inner, _SCRNA_KEYS)
        pb = inner.get("pseudobulk_de") if isinstance(inner.get("pseudobulk_de"), dict) else {}
        for grp, ginfo in (pb.get("per_group") or {}).items():
            for ckey, comp in ((ginfo.get("per_comparison") or {}) if isinstance(ginfo, dict) else {}).items():
                if isinstance(comp, dict) and comp:
                    out.append(("pseudobulk_de", comp, f"scrna.pseudobulk.{_safe_id(grp)}.{_safe_id(ckey)}"))
    elif agent_name == "chromatin_agent":
        assays = [(k, root[k]) for k in _CHROMATIN_ASSAYS if isinstance(root.get(k), dict)]
        out.append(("__root__", root, None))
        if assays:
            for mod_key, sub in assays:
                nested = sub.get("findings") if isinstance(sub.get("findings"), dict) else sub
                hint = _safe_id(mod_key) if len(assays) > 1 else None
                out.append(("__root__", nested, hint))
                _add(nested, _CHROM_KEYS, hint=hint)
        else:
            _add(root, _CHROM_KEYS)
    elif agent_name == "genome_arch_agent":
        out.append(("__root__", root, None))
        _add(root, ("fastq_alignment", "preprocessing", "inspection", "balancing", "qc", "figures"))
        topo = root.get("topology") if isinstance(root.get("topology"), dict) else root
        _add(topo, ("compartments", "tads", "loops", "differential_3d"))
    elif agent_name == "integration_agent":
        out.append(("__root__", root, None))
        _add(root, ("preprocessing", "paired_qc", "wnn", "mofa", "peak_to_gene", "eregulons", "figures"))
    else:
        out.append(("__root__", root, None))
        _add(root, tuple(_STAGE_TO_ANALYSES))
    return out


def _extract_signals(scope: str, stage: dict[str, Any], id_hint: str | None) -> list[_Signal]:
    sigs: list[_Signal] = []
    stage_failed = str(stage.get("status") or "").strip().lower() in _FAILED_STATUSES

    if scope == "figures":
        errs = stage.get("_errors")
        raw_items = errs if isinstance(errs, (list, tuple)) else (
            [f"{k}: {v}" if v else str(k) for k, v in errs.items()] if isinstance(errs, dict) else ()
        )
        sigs.extend(_Signal(scope, "figures._errors", "warning", msg, True, id_hint) for item in raw_items if (msg := _sanitize(item)))

    if (err_val := stage.get("error")) and not isinstance(err_val, bool) and (msg := _sanitize(err_val)):
        sev: CaveatSeverity = "blocking" if (stage_failed and scope not in _AUXILIARY_SCOPES and scope != "__root__") else "warning"
        sigs.append(_Signal(scope, "error", sev, msg, True, id_hint))

    for k, v in stage.items():
        if isinstance(k, str) and k.endswith("_error") and v and not isinstance(v, bool) and (msg := _sanitize(f"{k}: {v}")):
            sigs.append(_Signal(scope, k, "warning", msg, True, id_hint))

    for key, noun in (("failed_databases", "database(s)"), ("failed_samples", "sample(s)"), ("failed_modalities", "modality(ies)")):
        raw = stage.get(key)
        if isinstance(raw, dict) and (names := tuple(str(k).strip() for k in raw if str(k).strip())):
            details = "; ".join(f"{k} ({_sanitize(v, 100)})" if v else str(k) for k, v in raw.items() if str(k).strip())
            sigs.append(_Signal(scope, key, "warning", _sanitize(f"Partial failure in {scope}: failed {noun} {details}."), True, id_hint, names))
        elif isinstance(raw, (list, tuple)) and (names := tuple(_sanitize(x, 80) for x in raw if str(x or "").strip())):
            sigs.append(_Signal(scope, key, "warning", _sanitize(f"Partial failure in {scope}: failed {noun} {', '.join(names)}."), True, id_hint, names))

    for k, v in stage.items():
        if isinstance(k, str) and k.endswith("_degraded") and v is True:
            prefix = k[: -len("_degraded")]
            src, reason = stage.get(f"{prefix}_source"), stage.get(f"{prefix}_reason") or stage.get("degradation_reason")
            parts = [f"{prefix.replace('_', ' ')} degraded"] + ([f"(source: {src})"] if src else []) + ([f": {reason}"] if reason else [])
            sigs.append(_Signal(scope, k, "warning", _sanitize(" ".join(parts) + "."), True, id_hint))

    if isinstance(pfr := stage.get("phasing_fallback_reason"), str) and pfr.strip():
        sigs.append(_Signal(scope, "phasing_fallback_reason", "warning", _sanitize(pfr), True, id_hint))

    for item in (stage.get("degradations") if isinstance(stage.get("degradations"), list) else ()):
        if isinstance(item, dict) and (msg := _sanitize(item.get("message") or "")):
            d_sev = item.get("severity") if item.get("severity") in _SEV_RANK else "warning"
            sigs.append(_Signal(str(item.get("scope") or scope), str(item.get("field") or "degradation"), d_sev, msg, True, id_hint))

    ls = stage.get("lfc_shrinkage")
    if isinstance(ls, dict) and ls.get("requested") is True:
        n_failed = ls.get("n_failed")
        has_failed = isinstance(n_failed, int) and not isinstance(n_failed, bool) and n_failed > 0
        if ls.get("applied") is False or has_failed:
            reason = ls.get("reason") or "; ".join(str(r) for r in (ls.get("failed_reasons") or []) if r)
            msg = (
                f"LFC shrinkage was requested but not applied{f' ({reason})' if reason else ''}; reported log2FC values are unshrunken MLE estimates."
                if ls.get("applied") is False
                else f"LFC shrinkage was partially applied ({ls.get('n_applied', 0)} succeeded, {n_failed} failed); failed comparisons retain unshrunken MLE estimates."
            )
            sigs.append(_Signal(scope, "lfc_shrinkage", "warning", _sanitize(msg), True, id_hint))

    if isinstance(mnc := stage.get("metrics_not_computed"), (list, tuple)) and (items := tuple(_sanitize(x, 120) for x in mnc if str(x or "").strip())):
        sigs.append(_Signal(scope, "metrics_not_computed", "info", _sanitize("Metrics not computed: " + ", ".join(items) + "."), True, id_hint, items))

    if isinstance(raw_w := stage.get("warnings"), (list, tuple)):
        sigs.extend(_Signal(scope, "warnings", "warning", txt, False, id_hint) for w in raw_w if (txt := _sanitize(w)))
    return sigs


def _match_blocks(sig: _Signal, blocks: list[NarrativeBlock]) -> list[NarrativeBlock]:
    if sig.scope == "__root__":
        qc_blocks = [b for b in blocks if b.analysis in _QC_SUMMARY_ANALYSES]
        if sig.id_hint:
            qc_blocks = [b for b in qc_blocks if sig.id_hint.lower() in b.id.lower()] or qc_blocks
        return [qc_blocks[0]] if qc_blocks else ([] if not sig.structured else blocks[:1])
    if sig.id_hint and "." in sig.id_hint:
        return [b for b in blocks if b.id == sig.id_hint]
    matched = [b for b in blocks if b.analysis in _STAGE_TO_ANALYSES.get(sig.scope, (sig.scope,))]
    if sig.id_hint and matched:
        matched = [b for b in matched if sig.id_hint.lower() in b.id.lower()] or matched
    return [matched[0]] if (sig.scope == "figures" and matched) else matched


def _already_surfaced(block: NarrativeBlock, sig: _Signal) -> bool:
    norm_msg = _norm(sig.message)
    if any(_norm(c.text) == norm_msg and _SEV_RANK.get(c.severity, 1) >= _SEV_RANK.get(sig.severity, 1) for c in block.caveats):
        return True
    if sig.field == "error" and block.status == "error":
        return norm_msg in _norm(block.error or "") or any(norm_msg in _norm(c.text) for c in block.caveats)
    if sig.field == "failed_databases" and sig.dedup_tokens:
        return any(c.severity in {"warning", "blocking"} and all(t.lower() in c.text.lower() for t in sig.dedup_tokens) for c in block.caveats)
    if sig.field == "metrics_not_computed" and sig.dedup_tokens:
        return any("not computed" in c.text.lower() and all(t.lower() in c.text.lower() for t in sig.dedup_tokens) for c in block.caveats)
    return False


def _record_audit(block: NarrativeBlock, scope: str, field: str, severity: str, message: str) -> None:
    trail = block.metadata.setdefault("degradation_signals", [])
    key = (scope, field, severity, _norm(message))
    if not any((e.get("scope"), e.get("field"), e.get("severity"), _norm(str(e.get("message", "")))) == key for e in trail):
        trail.append({"scope": scope, "field": field, "severity": severity, "message": message})


def _apply_signal(block: NarrativeBlock, sig: _Signal) -> None:
    _record_audit(block, sig.scope, sig.field, sig.severity, sig.message)
    if not _already_surfaced(block, sig):
        block.caveats.append(Caveat(text=sig.message, severity=sig.severity))
    if not sig.structured:
        return
    if sig.severity == "blocking":
        block.status, block.block_type, block.claim, block.confidence = "error", "error", "", "insufficient"
        if not block.error:
            block.error = sig.message
    elif sig.severity == "warning" and block.confidence == "high":
        block.confidence = "medium"


def route_degradation_signals(agent_name: str, agent_result: dict, blocks: list[NarrativeBlock]) -> list[NarrativeBlock]:
    """Route structured degradation signals and free-text warnings onto ``blocks``."""
    out = list(blocks or [])
    orphans: dict[tuple[str, str | None], list[_Signal]] = {}
    for scope, stage, id_hint in _iter_stage_dicts(agent_name, agent_result):
        for sig in _extract_signals(scope, stage, id_hint):
            targets = _match_blocks(sig, out)
            if targets:
                for block in targets:
                    _apply_signal(block, sig)
            elif sig.scope == "__root__" and not sig.structured:
                if out:
                    _record_audit(out[0], sig.scope, sig.field, "info", sig.message)
            elif sig.structured and sig.severity in {"warning", "blocking"}:
                orphans.setdefault((sig.scope, sig.id_hint), []).append(sig)

    prefix, modality = _AGENT_META.get(agent_name, (_safe_id(agent_name), str(agent_name)))
    existing_ids = {b.id for b in out}
    for (scope, id_hint), group in orphans.items():
        is_blocking = any(s.severity == "blocking" for s in group)
        base_id = id_hint if (id_hint and "." in id_hint) else f"{prefix}.{_safe_id(scope)}"
        bid = base_id if base_id not in existing_ids else f"{base_id}.degradation"
        existing_ids.add(bid)
        caveats, seen_cav, trail = [], set(), []
        for s in group:
            if (k := (s.severity, _norm(s.message))) not in seen_cav:
                seen_cav.add(k)
                caveats.append(Caveat(text=s.message, severity=s.severity))
            trail.append({"scope": s.scope, "field": s.field, "severity": s.severity, "message": s.message})
        out.append(NarrativeBlock(
            id=bid, modality=modality, analysis=_STAGE_TO_ANALYSES.get(scope, (scope,))[0],
            block_type="error" if is_blocking else "limitation",
            title=f"{modality} {scope.replace('_', ' ')}",
            status="error" if is_blocking else "limitation", confidence="insufficient", claim="",
            caveats=caveats, error=next((s.message for s in group if s.severity == "blocking"), None),
            metadata={"degradation_signals": trail},
        ))
    return out
