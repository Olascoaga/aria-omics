"""ADR-080 production cross-modal associative concordance entrypoint."""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import pandas as pd

from aria.scripts._base import run_script
from aria.scripts.chromatin_peak_annotation import _resolve_gtf
from aria.utils.csv_artifacts import persist_csv
from aria.utils.stage_manifest import build_stage_manifest, stage_is_current, write_stage_manifest
from aria.utils.cross_modal_concordance import (
    aggregate_pathways, align_contrasts, atac_contrasts, contrast_arms, contrast_name,
    de_da_layer, gtf_symbols, load_atac_table, load_rna_table, normalize_inputs,
    parse_gtf_gene_models, pathway_layer, read_cpm, render_starburst,
    resolve_cpm_paths, rna_contrasts, rna_table_paths, sample_matched_layer,
    safe_label, summarize_de_da, summarize_tfs, tf_layer,
)

VERSION = "ADR-080.1"
DEFAULTS = {"padj_threshold": .05, "cis_distance_bp": 100_000, "promoter_bp": 2_000,
            "min_matched_samples": 6, "n_permutations": 200, "seed": 0}


def _json_value(value):
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def _input_files(params, rna, atac, rpath, apath, gtf):
    found = set()
    def visit(value):
        if isinstance(value, dict):
            for item in value.values():
                visit(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                visit(item)
        elif isinstance(value, (str, Path)):
            try:
                if Path(value).is_file():
                    found.add(str(value))
            except (OSError, ValueError):
                pass
    visit([params.get("rna_results"), params.get("chromatin_results"), rpath, apath, gtf])
    for comp in rna_contrasts(rna):
        visit(rna_table_paths(comp))
    for comp in atac_contrasts(atac):
        visit(comp.get("full_results_csv"))
    return [(f"input_{i}", path) for i, path in enumerate(sorted(found))]


def _empty_result(rm, am, sample_summary):
    de_da = summarize_de_da(pd.DataFrame())
    de_da.update(best_fisher_odds_ratio=None, best_fisher_pvalue=None, best_spearman_r=None, best_spearman_pvalue=None)
    return {"status": "success", "ran": False, "reason": None, "analysis": "cross_modal_concordance",
            "validation_level": "production", "evidence_tier": "associative", "associative_only": True,
            "rna_modality": rm, "atac_modality": am, "n_matched_contrasts": 0, "matched_contrasts": [],
            "de_da_concordance": de_da, "tf_triangulation": summarize_tfs(pd.DataFrame()),
            "pathway_concordance": aggregate_pathways({}, {}, pd.DataFrame(), False), "sample_matched_links": sample_summary,
            "tables": {"de_da_concordance_csv": None, "tf_triangulation_csv": None, "shared_pathways_csv": None, "sample_matched_links_csv": None},
            "figures": {}, "figures_svg": {}, "summary_json_path": None, "stage_manifest_path": None,
            "resumed_from_manifest": False, "warnings": []}


def _valid_params(params):
    try:
        alpha = float(params["padj_threshold"])
        if not 0 < alpha <= 1:
            return False
        for key in ("cis_distance_bp", "promoter_bp", "min_matched_samples", "n_permutations", "seed"):
            value = params[key]
            if isinstance(value, bool) or int(value) != float(value):
                return False
            params[key] = int(value)
        params["padj_threshold"] = alpha
        return (params["cis_distance_bp"] >= 0 and params["promoter_bp"] >= 0 and params["min_matched_samples"] >= 3
                and params["n_permutations"] >= 1 and params["seed"] >= 0
                and isinstance(params.get("sample_id_map", {}), dict)
                and all(isinstance(k, str) and isinstance(v, str) for k, v in params.get("sample_id_map", {}).items()))
    except (ValueError, TypeError, OverflowError):
        return False


def cross_modal_concordance(params: dict) -> dict:
    """Run available evidence layers and persist A4-addressed artifacts.

    RNA-to-ATAC sample maps are explicit, one-to-one, exact string mappings.
    Condition-level comparisons of scRNA cell groups retain the group for every
    gene rather than pooling donors or presenting cells as biological replicates.
    """
    params = {**DEFAULTS, **params}
    warnings = []
    rna, atac, rm, am = normalize_inputs(params.get("rna_results"), params.get("chromatin_results"))
    if not _valid_params(params):
        result = _empty_result(rm, am, {"ran": False, "reason": "invalid_parameters", "n_shared_samples": 0, "shared_sample_ids": []})
        result["reason"] = "invalid_parameters"
        return result
    out_dir = Path(params.get("output_dir") or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    gtf = _resolve_gtf(params, out_dir)
    gtf = str(gtf) if gtf and Path(gtf).is_file() else None
    rpath, apath = resolve_cpm_paths(params, rna, atac)
    inputs = _input_files(params, rna, atac, rpath, apath, gtf)
    # Include every in-memory upstream value, including missing-file hints. File
    # contents have independent signatures; resume/output locations are excluded.
    signature_params = {k: _json_value(v) for k, v in params.items() if k not in ("resume", "output_dir")}
    digest = hashlib.sha256(json.dumps(signature_params, sort_keys=True, default=str).encode()).hexdigest()
    manifest_params = {**{k: params[k] for k in DEFAULTS}, "input_payload_sha256": digest}
    summary_path = out_dir / "cross_modal_concordance_summary.json"
    manifest_path = out_dir / "cross_modal_concordance.stage_manifest.json"
    if params.get("resume") and summary_path.is_file():
        current, _ = stage_is_current(manifest_path, inputs=inputs, params=manifest_params, tool_version=VERSION)
        if current:
            try:
                cached = json.loads(summary_path.read_text())
                outputs = [p for key in ("tables", "figures", "figures_svg") for p in cached.get(key, {}).values() if p]
                if all(Path(p).is_file() and Path(p).stat().st_size > 0 for p in outputs):
                    cached["resumed_from_manifest"] = True
                    return cached
            except (OSError, ValueError, TypeError):
                warnings.append("cached_summary_unreadable")
    models = parse_gtf_gene_models(gtf) if gtf else {}
    aliases = gtf_symbols(gtf)
    rna_cpm, atac_cpm = read_cpm(rpath, aliases, warnings), read_cpm(apath, {}, warnings)
    expression = {str(k).upper(): float(v) for k, v in rna_cpm.mean(axis=1).items() if np.isfinite(v)}
    sample_frame, sample_summary = sample_matched_layer(rna_cpm, atac_cpm, models, params, warnings)
    result = _empty_result(rm, am, sample_summary)
    result["warnings"] = warnings
    pairs, had_contrasts = align_contrasts(rna, atac, warnings)
    da_frames, tf_frames, pw_frames = [], [], []
    rterms, aterms, pw_ran = {}, {}, False
    for rcomp, acomp, flipped, granularity in pairs:
        name = contrast_name(rcomp)
        canonical_name = "_vs_".join(contrast_arms(rcomp))
        contrast_id = f"{rcomp['group']}:{canonical_name}" if rcomp.get("group") is not None else canonical_name
        rframe, rscope = load_rna_table(rcomp, aliases, warnings, params["padj_threshold"])
        aframe, ascope = load_atac_table(acomp, atac, models, params, flipped, warnings)
        df, da = de_da_layer(rframe, aframe, params, contrast_id)
        da.update(annotation_scope=ascope, rna_scope=rscope, cis_distance_bp=params["cis_distance_bp"], promoter_bp=params["promoter_bp"])
        tf, tf_summary = tf_layer(atac, acomp, rframe, flipped, params["padj_threshold"], warnings, expression)
        pw, pw_summary, rt, at = pathway_layer(rna, atac, rcomp, acomp, params["padj_threshold"])
        for frame, target in ((df, da_frames), (tf, tf_frames), (pw, pw_frames)):
            if not frame.empty:
                frame["contrast_id"] = contrast_id
                target.append(frame)
        rterms.update(rt)
        aterms.update(at)
        pw_ran = pw_ran or pw_summary["ran"]
        result["matched_contrasts"].append({"contrast_id": contrast_id, "rna_contrast_name": name,
            "atac_contrast_name": contrast_name(acomp), "atac_sign_flipped": flipped,
            "matching_granularity": granularity, "cell_group": rcomp.get("group"),
            "de_da_summary": da, "tf_triangulation": tf_summary, "pathway_concordance": pw_summary})
    combine = lambda frames: pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    da_frame, tf_frame, pw_frame = combine(da_frames), combine(tf_frames), combine(pw_frames)
    aggregate = summarize_de_da(da_frame, params["promoter_bp"])
    # Pooled effects can mix contrasts; best per-contrast statistics are explicit.
    per_da = [x["de_da_summary"] for x in result["matched_contrasts"]]
    for kind in ("fisher", "spearman"):
        candidates = [x for x in per_da if x.get(f"{kind}_pvalue") is not None]
        best = min(candidates, key=lambda x: x[f"{kind}_pvalue"]) if candidates else {}
        aggregate[f"best_{kind}_pvalue"] = best.get(f"{kind}_pvalue")
        aggregate[f"best_{kind}_odds_ratio" if kind == "fisher" else "best_spearman_r"] = best.get("fisher_odds_ratio" if kind == "fisher" else "spearman_r")
    result.update(n_matched_contrasts=len(pairs), de_da_concordance=aggregate, tf_triangulation=summarize_tfs(tf_frame),
                  pathway_concordance=aggregate_pathways(rterms, aterms, pw_frame, pw_ran))
    result["ran"] = any(result[key]["ran"] for key in ("de_da_concordance", "tf_triangulation", "pathway_concordance", "sample_matched_links"))
    if not result["ran"]:
        result["reason"] = "incompatible_contrast_denominators" if had_contrasts and not pairs else "no_cross_modal_evidence_available"
    specs = (("de_da_concordance_csv", "cross_modal_de_da_concordance.csv", da_frame, aggregate["ran"]),
             ("tf_triangulation_csv", "cross_modal_tf_triangulation.csv", tf_frame, result["tf_triangulation"]["ran"]),
             ("shared_pathways_csv", "cross_modal_shared_pathways.csv", pw_frame, pw_ran),
             ("sample_matched_links_csv", "cross_modal_matched_sample_links.csv", sample_frame, sample_summary["ran"]))
    for key, filename, frame, ran in specs:
        if ran:
            path = str(out_dir / filename)
            if frame.empty and key == "shared_pathways_csv":
                frame = pd.DataFrame(columns=["term", "atac_term", "safe_label", "rna_padj", "atac_padj", "geometric_mean_padj", "contrast_id"])
            result["tables"][key] = persist_csv(path, filename, lambda f=frame, p=path: f.to_csv(p, index=False), warnings)
    if result["ran"]:
        png, svg = render_starburst(da_frame, tf_frame, pw_frame, aggregate, out_dir)
        result["figures"], result["figures_svg"] = {"cross_modal_starburst": png}, {"cross_modal_starburst": svg}
    portable_warnings = []
    for warning in warnings:
        warning = str(warning).replace(str(out_dir) + "/", "output_dir/")
        for _, path in inputs:
            warning = warning.replace(path, Path(path).name)
        portable_warnings.append(safe_label(warning))
    result["warnings"] = list(dict.fromkeys(portable_warnings))
    result.update(summary_json_path=str(summary_path), stage_manifest_path=str(manifest_path))
    result = _json_value(result)
    # Strict JSON disallows NaN/Inf, including undefined Fisher odds ratios.
    tmp = summary_path.with_name(summary_path.name + ".tmp")
    tmp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    tmp.replace(summary_path)
    write_stage_manifest(manifest_path, build_stage_manifest(stage="cross_modal_concordance", inputs=inputs,
                         params=manifest_params, tool_version=VERSION, extra={"evidence_tier": "associative"}))
    return result


if __name__ == "__main__":
    run_script(cross_modal_concordance)
