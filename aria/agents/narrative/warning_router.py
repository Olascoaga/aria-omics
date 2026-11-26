"""Centralized degradation and warning router for modality narrators (W-WARN).

Pure post-collection pass called from ``NarrativeRegistry.collect_blocks`` (and
cross-modal synthesis): structured fields cap confidence / mark failed stages
(A3) / synthesize orphans; free-text ``warnings`` route by structural position
only with ``severity="warning"`` (A2). Pure: stdlib + narrative.types only.
"""
from __future__ import annotations

import re
from typing import Any, NamedTuple

from aria.agents.narrative.types import Caveat, CaveatSeverity, NarrativeBlock

_TAG_RE, _WS_RE, _SAFE_ID_RE = re.compile(r"<[^>]+>"), re.compile(r"\s+"), re.compile(r"[^A-Za-z0-9_.-]+")
_FAILED_STATUSES = frozenset({"error", "failed"})
_AUXILIARY_SCOPES = frozenset({"plots", "figures", "lobo_jackknife", "outlier_sensitivity"})
_QC_SUMMARY_ANALYSES = frozenset({"qc", "sample_qc", "paired_qc", "summary"})
_CHROMATIN_ASSAYS = ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG")
_SEV_RANK = {"info": 0, "warning": 1, "blocking": 2}
SEX_CHR_SIG_FRAC_THRESHOLD = 0.20  # >=20% (>4x ~4.5% chrX+chrY gene baseline)
_AGENT_META = {
    "bulk_rna_agent": ("bulk", "bulk RNA-seq"), "scrna_agent": ("scrna", "scRNA-seq"),
    "rna_agent": ("scrna", "scRNA-seq"), "chromatin_agent": ("chromatin", "chromatin"),
    "genome_arch_agent": ("hic", "hic"), "integration_agent": ("multiome", "integration"),
    "cross_modal_concordance": ("integration.cross_modal", "cross_modal"),
}
_STAGE_TO_ANALYSES: dict[str, tuple[str, ...]] = {
    "preprocessing": ("preprocessing", "scatac_preprocessing", "bulk_atac_preprocessing"),
    "fastq_to_fragments": ("scatac_preprocessing",), "fragments_to_matrix": ("scatac_preprocessing",),
    "cohort_matrix_union": ("scatac_preprocessing",), "inspection": ("qc",), "balancing": ("qc",),
    "outlier_sensitivity": ("sample_qc", "differential_expression"), "peaks": ("peak_calling",),
    "peak_counts": ("peak_count_matrix",), "lsi": ("dimensionality_reduction",), "lsi_clustering": ("dimensionality_reduction",),
    "diffacc": ("differential_accessibility",), "bulk_diffacc": ("differential_accessibility",), "motifs": ("motif_enrichment",),
    "regulatory": ("regulatory_layers",), "footprint": ("differential_tf_footprinting", "footprinting"),
    "footprinting": ("differential_tf_footprinting", "footprinting"), "differential_expression": ("marker_discovery", "differential_expression"),
    "pseudobulk_pathways": ("pseudobulk_pathways", "pathway_enrichment"), "pathways": ("pathway_enrichment", "pseudobulk_pathways"),
    "pathway_ora": ("pathway_enrichment",), "contrast": ("differential_expression",), "plots": ("differential_expression", "sample_qc"),
    "lobo_jackknife": ("differential_expression", "pseudobulk_de", "differential_accessibility"),
    "figures": ("qc", "sample_qc", "paired_qc", "dimensionality_reduction", "differential_accessibility", "peak_calling"),
    "integration": ("qc", "data_quality"), "clustering": ("qc", "marker_discovery"),
    "cell_types": ("qc", "data_quality"), "sample_matched_links": ("concordance", "cross_modal_synthesis"),
}
_AGENT_STAGES: dict[str, tuple[str, ...]] = {
    "bulk_rna_agent": ("preprocessing", "sample_qc", "outlier_sensitivity", "figures"),
    "scrna_agent": ("qc", "integration", "clustering", "cell_types", "differential_expression", "differential_abundance", "pseudobulk_de", "pseudobulk_pathways", "pathways", "trajectory", "cell_communication", "figures"),
    "chromatin_agent": ("fastq_to_fragments", "fragments_to_matrix", "cohort_matrix_union", "preprocessing", "qc", "peaks", "peak_counts", "lsi", "lsi_clustering", "differential_accessibility", "diffacc", "bulk_diffacc", "motifs", "peak_annotation", "peak_ora", "footprint", "footprinting", "regulatory", "figures"),
    "genome_arch_agent": ("fastq_alignment", "preprocessing", "inspection", "balancing", "qc", "figures"),
    "integration_agent": ("preprocessing", "paired_qc", "wnn", "mofa", "peak_to_gene", "eregulons", "figures"),
}
_FALLBACK_STAGES = tuple(dict.fromkeys((*_STAGE_TO_ANALYSES, *(s for v in _AGENT_STAGES.values() for s in v), "compartments", "tads", "loops", "differential_3d")))


class _Signal(NamedTuple):
    scope: str; field: str; severity: CaveatSeverity; message: str
    structured: bool; id_hint: str | None = None; dedup_tokens: tuple[str, ...] = ()

def _safe_id(val: Any) -> str:
    return _SAFE_ID_RE.sub("_", str(val or "unknown")).strip("_") or "unknown"

def _sanitize(text: Any, max_len: int = 300) -> str:
    c = _WS_RE.sub(" ", _TAG_RE.sub("", str(text or ""))).strip()
    return c[:max_len] if len(c) > max_len else c

def _norm(text: str) -> str:
    return _WS_RE.sub(" ", text).strip().lower()

def _pos_int(val: Any) -> bool:
    return isinstance(val, int) and not isinstance(val, bool) and val > 0

def _str_items(raw: Any, max_len: int = 120) -> tuple[str, ...]:
    return tuple(_sanitize(x, max_len) for x in raw if str(x or "").strip()) if isinstance(raw, (list, tuple)) else ()


def _iter_stage_dicts(agent_name: str, agent_result: dict) -> list[tuple[str, dict[str, Any], str | None]]:
    if not isinstance(agent_result, dict):
        return []
    root = agent_result.get("findings") if isinstance(agent_result.get("findings"), dict) else agent_result
    out: list[tuple[str, dict[str, Any], str | None]] = []

    def _push(sc: str, d: Any, h: str | None = None) -> None:
        if isinstance(d, dict) and d:
            out.append((sc, d, h))

    def _add(container: dict[str, Any], keys: tuple[str, ...], hint: str | None = None) -> None:
        for k in keys:
            _push(k, container.get(k), hint)

    def _add_chrom(container: dict[str, Any], hint: str | None = None) -> None:
        _add(container, _AGENT_STAGES["chromatin_agent"], hint=hint)
        pb_h = f"chromatin.differential_accessibility.pseudobulk.{hint}" if hint else "chromatin.differential_accessibility.pseudobulk"
        for da_k in ("differential_accessibility", "diffacc", "bulk_diffacc"):
            if isinstance(da := container.get(da_k), dict):
                for sub, sub_h in ((da.get("pseudobulk"), pb_h), (da, hint)):
                    if isinstance(sub, dict):
                        if sub is not da:
                            _push("differential_accessibility", sub, sub_h)
                        for comp in (sub.get("comparisons") or []):
                            _push("differential_accessibility", comp, sub_h)
                            _push("lobo_jackknife", comp.get("lobo_jackknife") if isinstance(comp, dict) else None, sub_h)
        for fp_k in ("footprinting", "footprint"):
            if isinstance(fp := container.get(fp_k), dict):
                _push("footprinting", (fp.get("differential_summary") or {}).get("inference"), hint)

    if root is not agent_result and isinstance(agent_result.get("warnings"), (list, tuple)):
        out.append(("__root__", {"warnings": agent_result["warnings"]}, None))
    if agent_name in {"scrna_agent", "rna_agent"}:
        sc = root.get("scRNA") if isinstance(root.get("scRNA"), dict) else None
        inner = (sc.get("findings") if isinstance(sc.get("findings"), dict) else sc) if sc else root
        if inner is not root and isinstance(root.get("warnings"), (list, tuple)):
            out.append(("__root__", {"warnings": root["warnings"]}, None))
        out.append(("__root__", inner, None))
        _add(inner, _AGENT_STAGES["scrna_agent"])
        pb = inner.get("pseudobulk_de") if isinstance(inner.get("pseudobulk_de"), dict) else {}
        for grp, ginfo in (pb.get("per_group") or {}).items():
            for ckey, comp in ((ginfo.get("per_comparison") or {}) if isinstance(ginfo, dict) else {}).items():
                if isinstance(comp, dict) and comp:
                    tid = f"scrna.pseudobulk.{_safe_id(grp)}.{_safe_id(ckey)}"
                    out.append(("pseudobulk_de", comp, tid))
                    _push("lobo_jackknife", comp.get("lobo_jackknife"), tid)
        return out

    out.append(("__root__", root, None))
    if agent_name == "bulk_rna_agent":
        _add(root, _AGENT_STAGES["bulk_rna_agent"])
        for c in (x for x in (root.get("contrasts") or []) if isinstance(x, dict)):
            cn = _safe_id(c.get("name", "contrast"))
            cid = f"bulk.contrast.{cn}"
            out.append(("contrast", c, cid))
            for sk, tid in (("pathway_ora", f"bulk.pathway.{cn}"), ("plots", cid), ("lobo_jackknife", cid)):
                _push(sk, c.get(sk), tid)
    elif agent_name == "chromatin_agent":
        assays = [(k, root[k]) for k in _CHROMATIN_ASSAYS if isinstance(root.get(k), dict)]
        for mod_key, sub in assays:
            nested = sub.get("findings") if isinstance(sub.get("findings"), dict) else sub
            hint = _safe_id(mod_key) if len(assays) > 1 else None
            out.append(("__root__", nested, hint))
            _add_chrom(nested, hint=hint)
        if not assays:
            _add_chrom(root)
    elif agent_name == "cross_modal_concordance":
        _push("sample_matched_links", root.get("sample_matched_links"), "integration.cross_modal_sample_links")
    elif agent_name in _AGENT_STAGES:
        _add(root, _AGENT_STAGES[agent_name])
        if agent_name == "genome_arch_agent":
            _add(root.get("topology") if isinstance(root.get("topology"), dict) else root, ("compartments", "tads", "loops", "differential_3d"))
    else:
        _add(root, _FALLBACK_STAGES)
    return out


def _extract_signals(scope: str, stage: dict[str, Any], id_hint: str | None) -> list[_Signal]:
    sigs: list[_Signal] = []
    stage_failed = str(stage.get("status") or "").strip().lower() in _FAILED_STATUSES
    qc_fail_norms: set[str] = set()

    def _emit(field: str, sev: CaveatSeverity, msg: Any, structured: bool = True, tok: tuple[str, ...] = ()) -> None:
        if (txt := _sanitize(msg)):
            sigs.append(_Signal(scope, field, sev, txt, structured, id_hint, tok))

    if scope == "figures":
        errs = stage.get("_errors")
        for item in (errs if isinstance(errs, (list, tuple)) else ([f"{k}: {v}" if v else str(k) for k, v in errs.items()] if isinstance(errs, dict) else ())):
            _emit("figures._errors", "warning", item)
    if (err_val := stage.get("error")) and not isinstance(err_val, bool):
        _emit("error", "blocking" if (stage_failed and scope not in _AUXILIARY_SCOPES and scope != "__root__") else "warning", err_val)
    for k, v in stage.items():
        if isinstance(k, str) and k.endswith("_error") and v and not isinstance(v, bool):
            _emit(k, "warning", f"{k}: {v}")
        elif isinstance(k, str) and k.endswith("_degraded") and v is True:
            p = k[: -len("_degraded")]
            src, rsn = stage.get(f"{p}_source"), stage.get(f"{p}_reason") or stage.get("degradation_reason")
            _emit(k, "warning", " ".join([f"{p.replace('_', ' ')} degraded"] + ([f"(source: {src})"] if src else []) + ([f": {rsn}"] if rsn else [])) + ".")
    for key, noun in (("failed_databases", "database(s)"), ("failed_samples", "sample(s)"), ("failed_modalities", "modality(ies)")):
        raw = stage.get(key)
        if isinstance(raw, dict) and (names := tuple(str(k).strip() for k in raw if str(k).strip())):
            det = "; ".join(f"{k} ({_sanitize(v, 100)})" if v else str(k) for k, v in raw.items() if str(k).strip())
            _emit(key, "warning", f"Partial failure in {scope}: failed {noun} {det}.", True, names)
        elif (names := _str_items(raw, 80)):
            _emit(key, "warning", f"Partial failure in {scope}: failed {noun} {', '.join(names)}.", True, names)
    if isinstance(pfr := stage.get("phasing_fallback_reason"), str) and pfr.strip():
        _emit("phasing_fallback_reason", "warning", pfr)
    for item in (stage.get("degradations") if isinstance(stage.get("degradations"), list) else ()):
        if isinstance(item, dict) and (msg := _sanitize(item.get("message") or "")):
            sigs.append(_Signal(str(item.get("scope") or scope), str(item.get("field") or "degradation"), item.get("severity") if item.get("severity") in _SEV_RANK else "warning", msg, True, id_hint))
    if isinstance(ls := stage.get("lfc_shrinkage"), dict) and ls.get("requested") is True and (ls.get("applied") is False or _pos_int(ls.get("n_failed"))):
        rsn = ls.get("reason") or "; ".join(str(r) for r in (ls.get("failed_reasons") or []) if r)
        _emit("lfc_shrinkage", "warning", f"LFC shrinkage was requested but not applied{f' ({rsn})' if rsn else ''}; reported log2FC values are unshrunken MLE estimates." if ls.get("applied") is False else f"LFC shrinkage was partially applied ({ls.get('n_applied', 0)} succeeded, {ls.get('n_failed')} failed); failed comparisons retain unshrunken MLE estimates.")
    if isinstance(cov := stage.get("covariate_adjustment"), dict) and cov.get("unadjusted_confounding_risk") is True:
        c_rsn = _sanitize(cov.get("reason") or "", 100)
        _emit("covariate_adjustment.unadjusted_confounding_risk", "warning", f"Sample-matched correlations are unadjusted for sample covariates (unadjusted confounding risk{f': {c_rsn}' if c_rsn else ''}).")
    sfrac, n_sex = stage.get("sex_chromosome_significant_fraction"), stage.get("n_sex_chromosome_significant_links")
    if isinstance(sfrac, (int, float)) and not isinstance(sfrac, bool) and sfrac >= SEX_CHR_SIG_FRAC_THRESHOLD and (n_sex is None or n_sex > 0):
        cnt = f" ({n_sex}/{stage.get('n_significant_links')})" if (n_sex is not None and stage.get("n_significant_links")) else ""
        _emit("sex_chromosome_significant_fraction", "warning", f"Sex-chromosome dominance in significant links: {sfrac:.0%}{cnt} lie on chrX/chrY (threshold >= {SEX_CHR_SIG_FRAC_THRESHOLD:.0%}), indicating possible sex confounding.")
    if stage.get("pass_qc") is False:
        fails = _str_items(stage.get("qc_failures"))
        qc_fail_norms = {_norm(f) for f in fails}
        _emit("pass_qc", "warning", f"QC gate failed (pass_qc=False): {'; '.join(fails)}." if fails else "QC gate failed (pass_qc=False).", True, fails)
    p_st = str(stage.get("pairing_status") or "").strip().lower()
    ex_names = _str_items(stage.get("excluded_unpaired_samples") or stage.get("excluded_unpaired_replicates"), 80)
    n_unp = stage.get("n_unpaired_units_excluded")
    if p_st == "partial" or p_st.startswith("partial_") or ex_names or _pos_int(n_unp):
        det = f"excluded {len(ex_names)} unpaired donor/unit(s) ({', '.join(ex_names[:8])})" if ex_names else (f"excluded {n_unp} unpaired biological unit(s)" if _pos_int(n_unp) else "unpaired donors were excluded from the paired model")
        _emit("pairing_status", "warning", f"Partial paired design: {det}.", True, ex_names)
    if stage.get("batch_condition_confounded") is True:
        _emit("batch_condition_confounded", "warning", "Batch/donor covariate is 1:1 confounded with experimental condition (batch_condition_confounded=True); donor-specific technical effects cannot be disentangled from condition.")
    if _pos_int(n_nc := stage.get("n_nonconverged_excluded")):
        cmp_lbl = f" in {stage['test']} vs {stage['reference']}" if (stage.get("test") and stage.get("reference")) else ""
        _emit("n_nonconverged_excluded", "warning", f"Excluded {n_nc} non-converged apeGLM fit peak(s){cmp_lbl} (n_nonconverged_excluded={n_nc}).")
    if (items := _str_items(stage.get("metrics_not_computed"))):
        _emit("metrics_not_computed", "info", "Metrics not computed: " + ", ".join(items) + ".", True, items)
    for w in (stage.get("warnings") if isinstance(stage.get("warnings"), (list, tuple)) else ()):
        if (txt := _sanitize(w)) and _norm(txt) not in qc_fail_norms:
            sigs.append(_Signal(scope, "warnings", "warning", txt, False, id_hint))
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
    if sig.structured and sig.severity == "blocking":
        block.status, block.block_type, block.claim, block.confidence = "error", "error", "", "insufficient"
        block.error = block.error or sig.message
    elif sig.structured and sig.severity == "warning" and block.confidence == "high":
        block.confidence = "medium"


def route_degradation_signals(agent_name: str, agent_result: dict, blocks: list[NarrativeBlock]) -> list[NarrativeBlock]:
    """Route structured degradation signals and free-text warnings onto ``blocks``."""
    out, orphans = list(blocks or []), {}
    for scope, stage, id_hint in _iter_stage_dicts(agent_name, agent_result):
        for sig in _extract_signals(scope, stage, id_hint):
            if (targets := _match_blocks(sig, out)):
                for block in targets:
                    _apply_signal(block, sig)
            elif sig.scope == "__root__" and not sig.structured and out:
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
        seen_cav: set[tuple[str, str]] = set()
        caveats = [Caveat(text=s.message, severity=s.severity) for s in group if not ((k := (s.severity, _norm(s.message))) in seen_cav or seen_cav.add(k))]
        out.append(NarrativeBlock(
            id=bid, modality=modality, analysis=_STAGE_TO_ANALYSES.get(scope, (scope,))[0],
            block_type="error" if is_blocking else "limitation", title=f"{modality} {scope.replace('_', ' ')}",
            status="error" if is_blocking else "limitation", confidence="insufficient", claim="", caveats=caveats,
            error=next((s.message for s in group if s.severity == "blocking"), None),
            metadata={"degradation_signals": [{"scope": s.scope, "field": s.field, "severity": s.severity, "message": s.message} for s in group]},
        ))
    return out
