"""Cross-modal (RNA x ATAC) biological synthesis — ADR-028 Slice 2 (ADR-080).

Non-same-cell multi-omics: bulk RNA x bulk ATAC on separate libraries, or
unpaired scRNA x scATAC cohorts. Paired same-cell Multiome is owned by
``IntegrationNarrator`` (``multiome.triangulation``) and is not duplicated here.

Two deterministic sources, never an LLM:

1. **Concordance engine** (preferred): when ``cross_modal_concordance`` ran
   (``aria/scripts/cross_modal_concordance.py``), patterns are read directly from
   its structured result and every block routes to the IPC ledger node
   ``ledger://cross_modal/concordance`` (``modality="cross_modal"``,
   ``analysis="concordance"``), carrying the starburst figure.
2. **In-memory fallback**: when the engine result is absent but BOTH an RNA lane
   and the chromatin lane finished, matched contrasts are inspected from the
   agents' own structured outputs (DE genes x nearest-TSS-annotated DA peaks,
   motif / footprint direction x TF gene RNA, RNA ORA x peak ORA). Blocks route
   to the derived report node ``ledger://report/biological_synthesis``
   (``modality="integrated"``, ``analysis="cross_modal_synthesis"``).

Contrast matching is orientation-aware: an ATAC comparison run in the reverse
direction (``B_vs_A`` for RNA ``A_vs_B``) is sign-aligned (multiplied by -1); any
comparison whose arms differ is rejected, never force-matched.

Every number and named symbol a claim states is attached verbatim as an
``EvidenceItem`` so W-CLAIM strict verification passes; no claim, caveat or
evidence string carries causal vocabulary (C4).
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any

from aria.agents.narrative.run_ledger import (
    cross_modal_findings,
    cross_modal_status,
)
from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock
from aria.agents.narrative.validators import find_causal_language

_DONE = {"done", "success", "ok"}
_GO_ID_RE = re.compile(r"\s*\((?:GO|R-[A-Z]{3}|hsa|mmu|WP)[:\-]?[\w.\-]*\)\s*$")
_RAW_ID_PREFIXES = ("ENSG", "ENSMUSG", "ENST", "ENSMUST", "ENSR")
_PEAK_COORD_RE = re.compile(
    r"^(?:chr[\w.\-]+|[\dXYM]+|MT)[:\-_]\d+[-:_]\d+$", re.IGNORECASE
)
_MAX_NAMED = 3
_ENGINE_SOURCE = "cross_modal_concordance"
_MEMORY_SOURCE = "biological_synthesis"
_STARBURST_KEY = "cross_modal_starburst"

_CROSS_MODAL_CAVEAT = Caveat(
    "Cross-modal agreement between separately profiled RNA and chromatin "
    "libraries is an observational association at the condition or donor "
    "level; it is not same-cell evidence and does not demonstrate a mechanism.",
    "info",
)
_TF_CAVEAT = Caveat(
    "An enriched or footprinted motif is a sequence match with protected "
    "accessibility, not proof of TF binding; paralog TFs share motifs, so the "
    "specific family member remains ambiguous without direct binding data.",
    "info",
)
_SAMPLE_CAVEAT = Caveat(
    "Across-sample peak-gene correlation reflects co-variation between "
    "separately profiled libraries of the same biological units, not "
    "single-cell cis coupling.",
    "info",
)


# ── small pure helpers ────────────────────────────────────────────────────────

def _finite(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _int(value: Any) -> int | None:
    f = _finite(value)
    if f is None:
        return None
    return int(round(f))


def _sign(value: Any) -> int:
    f = _finite(value)
    if f is None or f == 0:
        return 0
    return 1 if f > 0 else -1


def _fmt_r(value: float) -> str:
    return f"{value:.2f}"


def _fmt_p(value: float) -> str:
    return f"{value:.2g}"


def _fmt_ratio(value: float) -> str:
    return f"{value:.3g}"


def _join(items: list[str]) -> str:
    if len(items) <= 1:
        return items[0] if items else ""
    sep = ", and " if len(items) > 2 else " and "
    return ", ".join(items[:-1]) + sep + items[-1]


def _safe(name: str) -> str:
    out = "".join(ch if ch.isalnum() else "_" for ch in str(name))
    return out.strip("_") or "x"


def _status_done(result: Any) -> bool:
    if not isinstance(result, dict) or not result:
        return False
    status = str(result.get("status") or "").lower().strip()
    if status:
        return status in _DONE
    findings = result.get("findings")
    return isinstance(findings, dict) and bool(findings)


def _clean_symbol(raw: Any) -> str | None:
    """A reportable gene/TF symbol: no raw Ensembl/numeric ids or coordinates."""
    sym = str(raw or "").strip()
    if not sym or sym in {"?", "nan", "None", "NA"} or sym.isdigit():
        return None
    if sym.upper().startswith(_RAW_ID_PREFIXES) or _PEAK_COORD_RE.match(sym):
        return None
    if any(ch in sym for ch in (" ", "\t", "\n", ";", ",", "<", ">")):
        return None
    if find_causal_language(sym):
        return None
    return sym


def _clean_term(term: Any) -> str:
    """Strip a trailing database accession so a term reads as a process name."""
    return _GO_ID_RE.sub("", str(term or "")).strip()


def _name_terms(rows: list, k: int = _MAX_NAMED) -> list[str]:
    """Up to k accession-free, causal-vocabulary-free shared term names."""
    out: list[str] = []
    for row in rows or []:
        if isinstance(row, dict):
            raw = row.get("clean_term") or row.get("term") or row.get("name")
        else:
            raw = row
        name = _clean_term(raw)
        if (
            name
            and name not in out
            and find_causal_language(name) is None
            and not any(ch in name for ch in "<>")
        ):
            out.append(name)
        if len(out) >= k:
            break
    return out


def _name_symbols(rows: list, key: str, k: int = _MAX_NAMED) -> list[str]:
    """Up to k clean single symbols; composite motif labels (``A::B``) and
    parenthesised variants are skipped, never split into invented genes."""
    out: list[str] = []
    for row in rows or []:
        raw = row.get(key) if isinstance(row, dict) else row
        if any(m in str(raw or "") for m in ("::", "(", ")", "/")):
            continue
        sym = _clean_symbol(raw)
        if sym and sym not in out:
            out.append(sym)
        if len(out) >= k:
            break
    return out


def _safe_reason(reason: Any) -> str | None:
    """A code-authored skip reason fit for a public sentence (or None)."""
    text = re.sub(r"\s+", " ", str(reason or "")).strip().rstrip(".")
    if not text or find_causal_language(text):
        return None
    text = text.replace("<", "fewer than ").replace(">", "more than ")
    return text[:200]


def _existing(path: Any) -> str | None:
    if isinstance(path, str) and path and Path(path).is_file():
        return path
    return None


# ── activation + engine-backed patterns ───────────────────────────────────────

def _paired_integration_done(agent_results: dict) -> bool:
    return _status_done((agent_results or {}).get("integration_agent"))


def _inactive(reason: str) -> dict[str, Any]:
    return {
        "active": False,
        "backed_by_concordance_engine": False,
        "skip_reason": reason,
        "contrasts": [],
        "rejected_contrasts": [],
        "sample_links": None,
        "sample_links_skip_reason": None,
        "figures": [],
        "tables": [],
    }


def _engine_contrasts(cm: dict) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    raw = [mc for mc in (cm.get("matched_contrasts") or []) if isinstance(mc, dict)]
    if not raw and any(
        isinstance(cm.get(k), dict)
        for k in ("de_da_concordance", "tf_triangulation", "pathway_concordance")
    ):
        raw = [{
            "contrast_id": "matched_contrast",
            "de_da_summary": cm.get("de_da_concordance"),
            "tf_triangulation": cm.get("tf_triangulation"),
            "pathway_concordance": cm.get("pathway_concordance"),
        }]
    for idx, mc in enumerate(raw):
        cid = str(mc.get("contrast_id") or mc.get("rna_contrast_name")
                  or f"contrast_{idx + 1}")
        de_da = mc.get("de_da_summary")
        tf = mc.get("tf_triangulation")
        pw = mc.get("pathway_concordance")
        out.append({
            "contrast_id": cid,
            "rna_contrast_name": str(mc.get("rna_contrast_name") or cid),
            "atac_contrast_name": str(mc.get("atac_contrast_name") or cid),
            "atac_sign_flipped": bool(mc.get("atac_sign_flipped")),
            "matching_granularity": mc.get("matching_granularity"),
            "de_da": de_da if isinstance(de_da, dict) else None,
            "tf": tf if isinstance(tf, dict) and tf.get("ran") is not False else None,
            "pathways": pw if isinstance(pw, dict) and pw.get("ran") is not False else None,
        })
    return out


def _engine_figures(cm: dict) -> list[dict[str, str]]:
    figs = cm.get("figures") or {}
    figs_svg = cm.get("figures_svg") or {}
    png = _existing(figs.get(_STARBURST_KEY)) if isinstance(figs, dict) else None
    svg = _existing(figs_svg.get(_STARBURST_KEY)) if isinstance(figs_svg, dict) else None
    if not png and not svg:
        return []
    rec = {
        "id": f"integration.{_STARBURST_KEY}",
        "path": png or svg,
        "caption": (
            "Cross-modal starburst: RNA log2 fold change versus ATAC log2 fold "
            "change for co-tested gene-peak pairs, coloured by concordance "
            "quadrant (associative, non-same-cell)."
        ),
    }
    if png and svg:
        rec["svg_path"] = svg
    return [rec]


def _engine_tables(cm: dict) -> list[dict[str, str]]:
    labels = {
        "de_da_concordance_csv": "DE×DA cis concordance pairs",
        "tf_triangulation_csv": "TF motif/footprint × TF RNA triangulation",
        "shared_pathways_csv": "Shared RNA and peak ORA terms",
        "sample_matched_links_csv": "Sample-matched cis peak-gene correlations",
    }
    out = []
    tables = cm.get("tables") or {}
    if not isinstance(tables, dict):
        return out
    for key, label in labels.items():
        path = _existing(tables.get(key))
        if path:
            out.append({"id": f"integration.{key}", "path": path, "label": label})
    return out


def _engine_patterns(cm: dict) -> dict[str, Any]:
    links = cm.get("sample_matched_links")
    links_ok = isinstance(links, dict) and links.get("ran") is True
    links_reason = None
    if isinstance(links, dict) and links.get("ran") is False:
        links_reason = _safe_reason(links.get("reason")) or "not_reported"
    return {
        "active": True,
        "backed_by_concordance_engine": True,
        "skip_reason": None,
        "rna_modality": str(cm.get("rna_modality") or "RNA"),
        "atac_modality": str(cm.get("atac_modality") or "ATAC"),
        "contrasts": _engine_contrasts(cm),
        "rejected_contrasts": [],
        "sample_links": links if links_ok else None,
        "sample_links_skip_reason": links_reason,
        "figures": _engine_figures(cm),
        "tables": _engine_tables(cm),
    }


# ── in-memory fallback extraction ─────────────────────────────────────────────

def _levels(contrast: dict, name: str) -> tuple[str, str] | None:
    for a_key, b_key in (("numerator", "denominator"), ("test", "reference"),
                         ("group_a", "group_b")):
        a, b = contrast.get(a_key), contrast.get(b_key)
        if a not in (None, "") and b not in (None, ""):
            return str(a).strip(), str(b).strip()
    parts = re.split(r"_vs_", str(name or ""))
    if len(parts) == 2 and all(p.strip() for p in parts):
        return parts[0].strip(), parts[1].strip()
    return None


def _norm_level(level: str) -> str:
    return str(level or "").strip().lower()


def _rna_terms(contrast: dict) -> list[str]:
    out: list[str] = []
    pw = contrast.get("pathways")
    groups = pw.values() if isinstance(pw, dict) else [pw or []]
    for rows in groups:
        for row in rows or []:
            raw = (row.get("term") or row.get("Term") or row.get("name")
                   if isinstance(row, dict) else row)
            name = _clean_term(raw)
            if name and name not in out:
                out.append(name)
    return out


def _bulk_gene_signs(contrast: dict) -> dict[str, tuple[str, float]]:
    """``{UPPER: (symbol, signed effect)}`` from the contrast's DE gene lists."""
    out: dict[str, tuple[str, float]] = {}
    for row in contrast.get("top_genes") or []:
        if not isinstance(row, dict):
            continue
        sym = _clean_symbol(row.get("symbol") or row.get("gene"))
        lfc = _finite(row.get("log2fc"))
        if sym and lfc:
            out.setdefault(sym.upper(), (sym, lfc))
    ids = [str(x) for x in (contrast.get("all_sig_gene_ids") or [])]
    syms = [str(x) for x in (contrast.get("all_sig_genes") or [])]
    id_to_sym = dict(zip(ids, syms)) if ids and len(ids) == len(syms) else {}
    for key, sign in (("up_gene_ids", 1.0), ("down_gene_ids", -1.0)):
        for gid in contrast.get(key) or []:
            sym = _clean_symbol(id_to_sym.get(str(gid), str(gid)))
            if sym:
                out.setdefault(sym.upper(), (sym, sign))
    for row in contrast.get("all_sig") or []:
        if isinstance(row, dict):
            sym = _clean_symbol(row.get("gene") or row.get("symbol"))
            lfc = _finite(row.get("log2fc"))
            if sym and lfc:
                out.setdefault(sym.upper(), (sym, lfc))
    return out


def _rna_contrasts(agent_results: dict) -> tuple[str | None, list[dict]]:
    """Structured RNA contrasts (bulk first, else scRNA pseudobulk)."""
    bulk = agent_results.get("bulk_rna_agent")
    if _status_done(bulk) or (
        isinstance(bulk, dict) and isinstance(bulk.get("contrasts"), list)
    ):
        findings = bulk.get("findings", bulk) if isinstance(bulk, dict) else {}
        out = []
        for idx, c in enumerate((findings or {}).get("contrasts") or []):
            if not isinstance(c, dict) or c.get("status", "success") != "success":
                continue
            name = str(c.get("name") or f"contrast_{idx + 1}")
            out.append({
                "name": name,
                "levels": _levels(c, name),
                "genes": _bulk_gene_signs(c),
                "terms": _rna_terms(c),
            })
        if out:
            return "bulk_RNA", out

    sc = agent_results.get("scrna_agent") or agent_results.get("rna_agent")
    if _status_done(sc):
        try:
            from aria.agents import _narrative_scrna
            sf = _narrative_scrna.unwrap_scrna_findings(sc) or {}
        except Exception:
            sf = sc.get("findings", sc) if isinstance(sc, dict) else {}
        pb = sf.get("pseudobulk_de") if isinstance(sf, dict) else None
        per_group = pb.get("per_group") if isinstance(pb, dict) else None
        out = []
        for group, gdata in (per_group or {}).items():
            if not isinstance(gdata, dict):
                continue
            for comp_key, comp in (gdata.get("per_comparison") or {}).items():
                if not isinstance(comp, dict):
                    continue
                if str(comp.get("status") or "").lower() != "success":
                    continue
                out.append({
                    "name": f"{group}:{comp_key}",
                    "levels": _levels(comp, str(comp_key)),
                    "genes": _bulk_gene_signs(comp),
                    "terms": _rna_terms(comp),
                })
        if out:
            return "scRNA", out
    return None, []


def _chromatin(agent_results: dict) -> tuple[str | None, dict]:
    chrom = agent_results.get("chromatin_agent")
    if not isinstance(chrom, dict):
        return None, {}
    try:
        from aria.agents.narrative.narrators.chromatin import (
            _is_bulk_atac,
            unwrap_chromatin_findings,
        )
        cf = unwrap_chromatin_findings(chrom) or {}
        is_bulk = _is_bulk_atac(cf)
    except Exception:
        cf = chrom.get("findings", chrom) or {}
        is_bulk = False
    if not cf or not (_status_done(chrom) or chrom.get("status") is None):
        return None, {}
    return ("bulk_ATAC" if is_bulk else "scATAC"), cf


def _atac_comparisons(cf: dict) -> list[dict]:
    da = cf.get("differential_accessibility") or {}
    if not isinstance(da, dict):
        return []
    comps = da.get("comparisons")
    if not comps and isinstance(da.get("pseudobulk"), dict):
        comps = da["pseudobulk"].get("comparisons")
    return [
        c for c in (comps or [])
        if isinstance(c, dict) and c.get("status", "success") == "success"
        and c.get("test") not in (None, "") and c.get("reference") not in (None, "")
    ]


def _read_csv_rows(path: Any) -> list[dict]:
    p = _existing(path)
    if not p:
        return []
    try:
        with open(p, newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))
    except (OSError, csv.Error, UnicodeDecodeError):
        return []


def _da_peaks(comp: dict) -> dict[str, float]:
    """Significant DA peaks -> log2FC (test vs reference orientation)."""
    out: dict[str, float] = {}
    for row in _read_csv_rows(comp.get("full_results_csv")):
        if str(row.get("significant", "")).strip().lower() != "true":
            continue
        peak = str(row.get("peak") or "").strip()
        lfc = _finite(row.get("log2FoldChange") or row.get("log2fc"))
        if peak and lfc:
            out[peak] = lfc
    if not out:
        for row in comp.get("top_peaks") or []:
            if isinstance(row, dict):
                peak = str(row.get("peak") or row.get("feature") or "").strip()
                lfc = _finite(row.get("log2fc"))
                if peak and lfc:
                    out[peak] = lfc
    return out


def _same_comparison(rec: dict, test: str, ref: str) -> bool:
    return (_norm_level(rec.get("test")) == _norm_level(test)
            and _norm_level(rec.get("reference")) == _norm_level(ref))


def _peak_genes(cf: dict, test: str, ref: str) -> dict[str, str]:
    ann = cf.get("peak_annotation") or {}
    if not isinstance(ann, dict) or not ann.get("ran"):
        return {}
    out: dict[str, str] = {}
    for rec in ann.get("comparisons") or []:
        if not isinstance(rec, dict) or not _same_comparison(rec, test, ref):
            continue
        rows = list(rec.get("rows") or []) + _read_csv_rows(rec.get("annotation_csv"))
        for row in rows:
            if not isinstance(row, dict):
                continue
            peak = str(row.get("peak") or "").strip()
            gene = _clean_symbol(row.get("nearest_gene") or row.get("gene"))
            if peak and gene:
                out.setdefault(peak, gene)
    return out


def _peak_ora_terms(cf: dict, test: str, ref: str) -> list[str] | None:
    ora = cf.get("peak_ora") or {}
    if not isinstance(ora, dict) or not ora.get("ran"):
        return None
    for rec in ora.get("comparisons") or []:
        if isinstance(rec, dict) and _same_comparison(rec, test, ref):
            return _rna_terms(rec)
    return None


def _tf_symbols(name: Any) -> tuple[list[str], bool]:
    raw = str(name or "").split("(")[0].strip()
    if not raw:
        return [], False
    dimer = "::" in raw
    parts = raw.split("::") if dimer else [raw]
    syms = [s for s in (_clean_symbol(p.strip()) for p in parts) if s]
    return syms, dimer


def _chromatin_tf_directions(cf: dict, test: str, ref: str) -> dict[str, dict]:
    """``{UPPER_TF: {tf, direction(+1/-1 test-vs-ref), layers, paralog}}``."""
    out: dict[str, dict] = {}
    comp_key = f"{test}_vs_{ref}"

    def _add(name: Any, direction: int, layer: str, family: Any = None) -> None:
        syms, dimer = _tf_symbols(name)
        for sym in syms:
            rec = out.setdefault(sym.upper(), {
                "tf": sym, "directions": set(), "layers": set(),
                "paralog": False, "families": set(),
            })
            rec["directions"].add(direction)
            rec["layers"].add(layer)
            rec["paralog"] = rec["paralog"] or dimer
            if family:
                rec["families"].add(str(family).strip().lower())

    motifs = cf.get("motifs") or {}
    if isinstance(motifs, dict) and motifs.get("ran"):
        padj_max = _finite(motifs.get("padj_max"))
        for group, gdata in (motifs.get("per_group") or {}).items():
            g = str(group)
            if not isinstance(gdata, dict):
                continue
            if g.lower().startswith(comp_key.lower() + "::up_in_"):
                level = g.split("::up_in_", 1)[1]
                if _norm_level(level) == _norm_level(test):
                    direction = 1
                elif _norm_level(level) == _norm_level(ref):
                    direction = -1
                else:
                    continue
            elif g.lower() == f"{comp_key.lower()}_up":
                direction = 1
            elif g.lower() == f"{comp_key.lower()}_down":
                direction = -1
            else:
                continue
            for motif in gdata.get("top_motifs") or []:
                if not isinstance(motif, dict):
                    continue
                padj = _finite(motif.get("padj"))
                if padj is not None and padj_max is not None and padj >= padj_max:
                    continue
                _add(motif.get("name") or motif.get("tf"), direction, "motif",
                     motif.get("family"))

    fp = cf.get("footprinting") or {}
    if isinstance(fp, dict) and fp.get("ran") is True:
        summary = fp.get("differential_summary") or {}
        inference = summary.get("inference") if isinstance(summary, dict) else None
        if isinstance(inference, dict) and inference.get("status") == "success":
            ga, gb = str(fp.get("group_a") or ""), str(fp.get("group_b") or "")
            alpha = _finite(inference.get("alpha")) or 0.05
            if {_norm_level(ga), _norm_level(gb)} == {_norm_level(test), _norm_level(ref)}:
                for favored in (ga, gb):
                    direction = 1 if _norm_level(favored) == _norm_level(test) else -1
                    for row in summary.get(f"top_toward_{favored}") or []:
                        if not isinstance(row, dict):
                            continue
                        padj = _finite(row.get("padj"))
                        if padj is None or padj > alpha:
                            continue
                        _add(row.get("tf"), direction, "footprint")

    # A motif family carried by more than one evaluated TF is paralog-ambiguous.
    fam_count: dict[str, int] = {}
    for rec in out.values():
        for fam in rec["families"]:
            fam_count[fam] = fam_count.get(fam, 0) + 1
    for rec in out.values():
        if any(fam_count.get(f, 0) > 1 for f in rec["families"]):
            rec["paralog"] = True
    return out


def _memory_de_da(da_peaks: dict[str, float], peak_genes: dict[str, str],
                  rna_genes: dict[str, tuple[str, float]], flip: int) -> dict | None:
    annotated = {p: g for p, g in peak_genes.items() if p in da_peaks}
    if not annotated or not rna_genes:
        return None
    pairs = []
    for peak, gene in annotated.items():
        rna = rna_genes.get(gene.upper())
        if rna is None:
            continue
        atac_sign = _sign(da_peaks[peak]) * flip
        rna_sign = _sign(rna[1])
        if not atac_sign or not rna_sign:
            continue
        pairs.append((abs(da_peaks[peak]), rna[0], atac_sign, rna_sign))
    pairs.sort(key=lambda t: (-t[0], t[1]))
    up = sum(1 for _, _, a, r in pairs if a > 0 and r > 0)
    down = sum(1 for _, _, a, r in pairs if a < 0 and r < 0)
    top = []
    for _, gene, a, r in pairs:
        if a == r and gene not in top:
            top.append(gene)
    return {
        "n_da_peaks_annotated": len(annotated),
        "n_both_significant": len(pairs),
        "n_concordant_up": up,
        "n_concordant_down": down,
        "n_concordant_total": up + down,
        "n_discordant_total": len(pairs) - up - down,
        "top_concordant_genes": [{"gene": g} for g in top[:10]],
    }


def _memory_tf(tf_dirs: dict[str, dict], rna_genes: dict[str, tuple[str, float]],
               flip: int) -> dict | None:
    if not tf_dirs:
        return None
    conc, disc, unchanged, ambiguous_dir = [], 0, 0, 0
    n_par = 0
    for key in sorted(tf_dirs):
        rec = tf_dirs[key]
        if len(rec["directions"]) != 1:
            ambiguous_dir += 1
            continue
        chrom_dir = next(iter(rec["directions"])) * flip
        n_par += int(bool(rec["paralog"]))
        rna = rna_genes.get(key)
        if rna is None:
            unchanged += 1
        elif _sign(rna[1]) == chrom_dir:
            conc.append({"tf": rec["tf"], "status": "concordant",
                         "paralog_motif_ambiguity": rec["paralog"]})
        else:
            disc += 1
    n_eval = len(conc) + disc + unchanged
    if n_eval == 0:
        return None
    return {
        "ran": True,
        "n_tfs_evaluated": n_eval,
        "n_concordant_tfs": len(conc),
        "n_discordant_tfs": disc,
        "n_chromatin_shift_rna_unchanged": unchanged,
        "n_paralog_ambiguous": n_par,
        "top_concordant_tfs": conc,
    }


def _memory_pathways(rna_terms: list[str], atac_terms: list[str] | None) -> dict | None:
    if atac_terms is None or not rna_terms or not atac_terms:
        return None
    atac_idx = {t.lower(): t for t in atac_terms}
    shared = [t for t in rna_terms if t.lower() in atac_idx]
    union = {t.lower() for t in rna_terms} | set(atac_idx)
    return {
        "ran": True,
        "n_rna_terms": len(rna_terms),
        "n_atac_terms": len(atac_terms),
        "n_shared_terms": len(shared),
        "jaccard_index": round(len(shared) / len(union), 4) if union else None,
        "shared_terms": [{"term": t, "clean_term": t} for t in shared],
    }


def _memory_patterns(agent_results: dict) -> dict[str, Any]:
    rna_mod, rna_contrasts = _rna_contrasts(agent_results)
    atac_mod, cf = _chromatin(agent_results)
    if not rna_contrasts or not cf:
        return _inactive("single_modality_or_missing_lane")
    comps = _atac_comparisons(cf)
    if not comps:
        return _inactive("no_chromatin_comparisons")

    contrasts: list[dict] = []
    rejected: list[dict] = []
    for rc in rna_contrasts:
        levels = rc["levels"]
        if levels is None:
            rejected.append({"rna_contrast_name": rc["name"],
                             "reason": "contrast_levels_unparseable"})
            continue
        num, den = (_norm_level(x) for x in levels)
        match = None
        for comp in comps:
            t, r = _norm_level(comp["test"]), _norm_level(comp["reference"])
            if (t, r) == (num, den):
                match = (comp, False)
                break
            if (t, r) == (den, num):
                match = (comp, True)
        if match is None:
            rejected.append({"rna_contrast_name": rc["name"],
                             "reason": "no_orientation_compatible_chromatin_comparison"})
            continue
        comp, flipped = match
        test, ref = str(comp["test"]), str(comp["reference"])
        flip = -1 if flipped else 1
        contrasts.append({
            "contrast_id": rc["name"],
            "rna_contrast_name": rc["name"],
            "atac_contrast_name": f"{test}_vs_{ref}",
            "atac_sign_flipped": flipped,
            "matching_granularity": "contrast_levels",
            "de_da": _memory_de_da(_da_peaks(comp), _peak_genes(cf, test, ref),
                                   rc["genes"], flip),
            "tf": _memory_tf(_chromatin_tf_directions(cf, test, ref),
                             rc["genes"], flip),
            "pathways": _memory_pathways(rc["terms"],
                                         _peak_ora_terms(cf, test, ref)),
        })
    if not contrasts:
        out = _inactive("no_orientation_compatible_contrasts")
        out["rejected_contrasts"] = rejected
        return out
    return {
        "active": True,
        "backed_by_concordance_engine": False,
        "skip_reason": None,
        "rna_modality": rna_mod,
        "atac_modality": atac_mod,
        "contrasts": contrasts,
        "rejected_contrasts": rejected,
        "sample_links": None,
        "sample_links_skip_reason": None,
        "figures": [],
        "tables": [],
    }


def detect_cross_modal_patterns(agent_results: dict,
                                exp_ctx: dict | None = None) -> dict[str, Any]:
    """Deterministic non-same-cell RNA x ATAC pattern manifest (Slice 2).

    Returns ``{"active": False, "skip_reason": ...}`` for single-modality runs,
    for paired same-cell Multiome without an engine run (IntegrationNarrator
    already owns that triangulation), and when the engine was present but did
    not run (skipped/error/no matched contrasts) — never a fallback that would
    contradict the engine's own honest skip.
    """
    ar = agent_results or {}
    if "cross_modal_concordance" in ar:
        cm = cross_modal_findings(ar)
        status, _reason = cross_modal_status(cm)
        if status != "ran":
            return _inactive(f"cross_modal_concordance_{status}")
        return _engine_patterns(cm)
    if _paired_integration_done(ar):
        return _inactive("paired_same_cell_integration_owns_triangulation")
    return _memory_patterns(ar)


# ── block composition ─────────────────────────────────────────────────────────

def _route(patterns: dict) -> dict[str, str]:
    if patterns.get("backed_by_concordance_engine"):
        return {"modality": "cross_modal", "analysis": "concordance",
                "limit_analysis": "concordance", "source": _ENGINE_SOURCE}
    return {"modality": "integrated", "analysis": "cross_modal_synthesis",
            "limit_analysis": "cross_modal_limitations", "source": _MEMORY_SOURCE}


def _ev(route: dict, label: str, value: Any) -> EvidenceItem:
    return EvidenceItem(label=label, value=value, source=route["source"])


def _block_id(base: str, contrast_id: str, multi: bool) -> str:
    return f"{base}.{_safe(contrast_id)}" if multi else base


def _de_da_block(route, mc, multi, engine) -> NarrativeBlock | None:
    d = mc.get("de_da") or {}
    cid = mc["contrast_id"]
    n_conc = _int(d.get("n_concordant_total"))
    n_disc = _int(d.get("n_discordant_total"))
    if n_conc is None or n_disc is None:
        return None
    n_up = _int(d.get("n_concordant_up"))
    n_down = _int(d.get("n_concordant_down"))
    ev = [_ev(route, "matched contrast", cid),
          _ev(route, "sign-concordant co-significant pairs", n_conc),
          _ev(route, "sign-discordant co-significant pairs", n_disc)]
    split = ""
    if n_up is not None and n_down is not None:
        ev += [_ev(route, "concordant pairs, both up", n_up),
               _ev(route, "concordant pairs, both down", n_down)]
        split = f" ({n_up} up, {n_down} down)"
    metrics: dict[str, Any] = {"n_concordant_total": n_conc,
                               "n_discordant_total": n_disc}
    if engine:
        n_pairs = _int(
            d.get("n_gene_peak_pairs_tested")
            if d.get("n_gene_peak_pairs_tested") is not None
            else d.get("n_pairs_tested")
        )
        if n_pairs is None or n_pairs <= 0:
            return None
        ev.append(_ev(route, "co-tested gene-peak pairs", n_pairs))
        metrics["n_gene_peak_pairs_tested"] = n_pairs
        claim = (
            f"For {cid}, DE×DA cis concordance over {n_pairs} co-tested "
            f"gene-peak pair(s) found {n_conc} sign-concordant{split} and "
            f"{n_disc} sign-discordant co-significant pair(s)"
        )
        n_atac_only = _int(d.get("n_atac_da_rna_unchanged"))
        n_rna_only = _int(d.get("n_rna_de_atac_unchanged"))
        if n_atac_only is not None:
            ev.append(_ev(route, "accessibility-only shifts without RNA change",
                          n_atac_only))
            tail = ""
            if n_rna_only is not None:
                ev.append(_ev(route, "RNA-only changes without accessibility shift",
                              n_rna_only))
                tail = (f", and {n_rna_only} changed in RNA without an "
                        f"accessibility shift")
            claim += (f", while {n_atac_only} pair(s) shifted in accessibility "
                      f"without an RNA change{tail}")
        claim += "."
        r = _finite(d.get("spearman_r"))
        r_p = _finite(d.get("spearman_pvalue"))
        r_name = "Spearman"
        if r is None:
            r, r_p, r_name = (_finite(d.get("pearson_r")),
                              _finite(d.get("pearson_pvalue")), "Pearson")
        stats = []
        if r is not None:
            ev.append(_ev(route, f"{r_name} r (RNA vs ATAC log2FC)", _fmt_r(r)))
            txt = f"{r_name} r = {_fmt_r(r)}"
            if r_p is not None:
                ev.append(_ev(route, f"{r_name} p-value", _fmt_p(r_p)))
                txt += f", p = {_fmt_p(r_p)}"
            stats.append(f"RNA and ATAC log2 fold changes across co-tested "
                         f"pairs gave {txt}")
            metrics[f"{r_name.lower()}_r"] = r
        odds = _finite(d.get("fisher_odds_ratio"))
        fp = _finite(d.get("fisher_pvalue"))
        if odds is not None:
            ev.append(_ev(route, "directional Fisher odds ratio", _fmt_ratio(odds)))
            txt = (f"the directional Fisher test gave an odds ratio of "
                   f"{_fmt_ratio(odds)}")
            if fp is not None:
                ev.append(_ev(route, "directional Fisher p-value", _fmt_p(fp)))
                txt += f" (p = {_fmt_p(fp)})"
            stats.append(txt)
            metrics["fisher_odds_ratio"] = odds
        if stats:
            joined = "; ".join(stats)
            claim += " " + joined[0].upper() + joined[1:] + "."
    else:
        n_ann = _int(d.get("n_da_peaks_annotated"))
        n_both = _int(d.get("n_both_significant"))
        if not n_ann or n_both is None:
            return None
        ev += [_ev(route, "differentially accessible peaks with nearest-TSS gene",
                   n_ann),
               _ev(route, "DA peak-DE gene pairs", n_both)]
        metrics.update({"n_da_peaks_annotated": n_ann, "n_both_significant": n_both})
        claim = (
            f"For {cid}, {n_both} of {n_ann} differentially accessible peak(s) "
            f"with a nearest-TSS gene mapped to a gene that was also "
            f"differentially expressed in the matched RNA contrast: {n_conc} "
            f"sign-concordant{split} and {n_disc} sign-discordant."
        )
    genes = _name_symbols(d.get("top_concordant_genes") or [], "gene")
    if genes and n_conc:
        ev.append(_ev(route, "top concordant genes", ", ".join(genes)))
        claim += f" Top concordant genes include {_join(genes)}."
    return NarrativeBlock(
        id=_block_id("integration.cross_modal_de_da", cid, multi),
        modality=route["modality"], analysis=route["analysis"],
        block_type="integration",
        title="RNA × chromatin accessibility concordance (DE×DA)",
        status="success", confidence="medium", claim=claim, evidence=ev,
        caveats=[_CROSS_MODAL_CAVEAT], metrics=metrics,
        metadata={"contrast_id": cid, "cross_modal_layer": "de_da",
                  "backed_by_concordance_engine": engine},
    )


def _tf_block(route, mc, multi, engine) -> NarrativeBlock | None:
    t = mc.get("tf") or {}
    cid = mc["contrast_id"]
    n_eval = _int(t.get("n_tfs_evaluated"))
    n_conc = _int(t.get("n_concordant_tfs"))
    if not n_eval or n_conc is None:
        return None
    n_disc = _int(t.get("n_discordant_tfs")) or 0
    n_shift = _int(t.get("n_chromatin_shift_rna_unchanged")) or 0
    n_par = _int(t.get("n_paralog_ambiguous")) or 0
    ev = [_ev(route, "matched contrast", cid),
          _ev(route, "TFs evaluated", n_eval),
          _ev(route, "sign-concordant TFs", n_conc),
          _ev(route, "sign-discordant TFs", n_disc),
          _ev(route, "chromatin shift without RNA change (TFs)", n_shift),
          _ev(route, "paralog-ambiguous TF motifs", n_par)]
    no_rna = ("without an RNA change" if engine
              else "without a significant change of the TF gene in the RNA DE set")
    claim = (
        f"For {cid}, TF motif or footprint shifts in chromatin were compared "
        f"with each TF's own RNA change across {n_eval} evaluated TF(s): "
        f"{n_conc} were sign-concordant, {n_disc} sign-discordant, and "
        f"{n_shift} showed a chromatin shift {no_rna}."
    )
    tfs = _name_symbols(
        [r for r in (t.get("top_concordant_tfs") or [])
         if isinstance(r, dict)
         and str(r.get("status") or r.get("classification") or "concordant").lower() == "concordant"],
        "tf",
    )
    if tfs and n_conc:
        ev.append(_ev(route, "concordant TFs", ", ".join(tfs)))
        claim += f" Concordant TFs include {_join(tfs)}."
    if n_par > 0:
        claim += (
            f" {n_par} evaluated TF motif(s) are shared with paralog family "
            f"members, so the specific TF behind each motif signal is ambiguous "
            f"and the concordance is associative."
        )
    return NarrativeBlock(
        id=_block_id("integration.cross_modal_tf", cid, multi),
        modality=route["modality"], analysis=route["analysis"],
        block_type="integration",
        title="TF motif/footprint × TF RNA triangulation",
        status="success", confidence="low" if n_par else "medium",
        claim=claim, evidence=ev,
        caveats=[_CROSS_MODAL_CAVEAT, _TF_CAVEAT],
        metrics={"n_tfs_evaluated": n_eval, "n_concordant_tfs": n_conc,
                 "n_paralog_ambiguous": n_par},
        metadata={"contrast_id": cid, "cross_modal_layer": "tf_triangulation",
                  "backed_by_concordance_engine": engine},
    )


def _pathway_block(route, mc, multi, engine) -> NarrativeBlock | None:
    p = mc.get("pathways") or {}
    cid = mc["contrast_id"]
    n_rna = _int(p.get("n_rna_terms"))
    n_atac = _int(p.get("n_atac_terms"))
    n_shared = _int(p.get("n_shared_terms"))
    if not n_rna or not n_atac or n_shared is None:
        return None
    ev = [_ev(route, "matched contrast", cid),
          _ev(route, "RNA enriched terms", n_rna),
          _ev(route, "chromatin peak-gene enriched terms", n_atac),
          _ev(route, "shared enriched terms", n_shared)]
    if n_shared == 0:
        claim = (
            f"For {cid}, none of the {n_rna} RNA enriched term(s) was shared "
            f"with the {n_atac} enriched term(s) of genes near differentially "
            f"accessible peaks."
        )
    else:
        claim = (
            f"For {cid}, {n_shared} enriched biological term(s) were shared "
            f"between the {n_rna} RNA enriched term(s) and the {n_atac} "
            f"enriched term(s) of genes near differentially accessible peaks"
        )
        jac = _finite(p.get("jaccard_index"))
        if jac is not None:
            ev.append(_ev(route, "Jaccard index of term sets", _fmt_r(jac)))
            claim += f" (Jaccard index {_fmt_r(jac)})"
        claim += "."
        names = _name_terms(p.get("shared_terms") or [])
        if names:
            ev.append(_ev(route, "shared enriched term", "; ".join(names)))
            claim += f" Shared terms include {_join(names)}."
    return NarrativeBlock(
        id=_block_id("integration.cross_modal_pathways", cid, multi),
        modality=route["modality"], analysis=route["analysis"],
        block_type="integration",
        title="Convergent enriched terms across RNA and chromatin",
        status="success", confidence="medium", claim=claim, evidence=ev,
        caveats=[_CROSS_MODAL_CAVEAT], metrics={"n_shared_terms": n_shared},
        metadata={"contrast_id": cid, "cross_modal_layer": "pathways",
                  "backed_by_concordance_engine": engine},
    )


def _sample_links_block(route, links: dict) -> NarrativeBlock | None:
    n_samples = _int(links.get("n_shared_samples"))
    n_tested = _int(links.get("n_pairs_tested"))
    n_sig = _int(links.get("n_significant_links"))
    if not n_samples or not n_tested or n_sig is None:
        return None
    n_pos = _int(links.get("n_positive_links")) or 0
    n_neg = _int(links.get("n_negative_links")) or 0
    ev = [_ev(route, "sample-matched biological units", n_samples),
          _ev(route, "cis peak-gene pairs tested", n_tested),
          _ev(route, "significant sample-matched links", n_sig),
          _ev(route, "positive links", n_pos),
          _ev(route, "negative links", n_neg)]
    window = _int(links.get("cis_distance_bp"))
    win_txt = ""
    if window:
        ev.append(_ev(route, "cis window (bp from TSS)", window))
        win_txt = f" within {window} bp of the TSS"
    claim = (
        f"Across {n_samples} sample- or donor-matched biological units profiled "
        f"in both assays, {n_tested} cis peak-gene pair(s){win_txt} were tested "
        f"and {n_sig} showed a significant across-sample correlation between "
        f"accessibility and expression ({n_pos} positive, {n_neg} negative)."
    )
    perm = links.get("permutation_null") or {}
    efdr = _finite(perm.get("empirical_fdr_at_padj")) if isinstance(perm, dict) else None
    n_perm = _int(perm.get("n_permutations")) if isinstance(perm, dict) else None
    if efdr is not None and n_perm:
        ev += [_ev(route, "sample permutations", n_perm),
               _ev(route, "empirical permutation FDR", _fmt_p(efdr))]
        claim += (f" A sample-label permutation null ({n_perm} permutations) "
                  f"estimated an empirical FDR of {_fmt_p(efdr)}.")
    genes = _name_symbols(links.get("top_links") or [], "gene")
    if genes and n_sig:
        ev.append(_ev(route, "top linked genes", ", ".join(genes)))
        claim += f" Top linked genes include {_join(genes)}."
    return NarrativeBlock(
        id="integration.cross_modal_sample_links",
        modality=route["modality"], analysis=route["analysis"],
        block_type="integration",
        title="Sample-matched cis peak-gene correlation",
        status="success", confidence="medium", claim=claim, evidence=ev,
        caveats=[_CROSS_MODAL_CAVEAT, _SAMPLE_CAVEAT],
        metrics={"n_shared_samples": n_samples, "n_pairs_tested": n_tested,
                 "n_significant_links": n_sig},
        metadata={"cross_modal_layer": "sample_matched_links",
                  "backed_by_concordance_engine": True},
    )


def _limitations_block(route, patterns: dict, engine: bool) -> NarrativeBlock:
    rna_mod = str(patterns.get("rna_modality") or "RNA")
    atac_mod = str(patterns.get("atac_modality") or "ATAC")
    ev = [_ev(route, "scope", "non_same_cell_association_only"),
          _ev(route, "RNA assay", rna_mod),
          _ev(route, "chromatin assay", atac_mod)]
    claim = (
        f"This cross-modal synthesis compares {rna_mod} and {atac_mod} measured "
        f"in separate libraries or unpaired cohorts rather than in the same "
        f"cells, so agreement reflects condition- or donor-level co-variation "
        f"and is associative; it does not establish single-cell cis coupling or "
        f"a mechanism."
    )
    flipped = [mc for mc in patterns.get("contrasts") or []
               if mc.get("atac_sign_flipped")]
    for mc in flipped:
        ev.append(_ev(route, f"reverse-oriented chromatin comparison for "
                             f"{mc['contrast_id']}", mc["atac_contrast_name"]))
        claim += (
            f" Chromatin effects for {mc['contrast_id']} were sign-aligned to "
            f"the RNA orientation because accessibility was contrasted in the "
            f"reverse direction ({mc['atac_contrast_name']})."
        )
    rejected = patterns.get("rejected_contrasts") or []
    if rejected:
        ev.append(_ev(route, "RNA contrasts without a compatible chromatin "
                             "comparison", len(rejected)))
        claim += (
            f" {len(rejected)} RNA contrast(s) had no orientation-compatible "
            f"chromatin comparison and were not compared."
        )
    reason = patterns.get("sample_links_skip_reason")
    if engine and reason:
        ev.append(_ev(route, "sample-matched correlation skip reason", reason))
        claim += f" Sample-matched peak-gene correlation was not reported ({reason})."
    elif not engine:
        claim += (
            " No sample-matched peak-gene correlation was computed in this "
            "in-memory synthesis."
        )
    return NarrativeBlock(
        id="integration.cross_modal_limitations",
        modality=route["modality"], analysis=route["limit_analysis"],
        block_type="limitation",
        title="Limits of cross-modal interpretation",
        status="success", confidence="high", claim=claim, evidence=ev,
        caveats=[_CROSS_MODAL_CAVEAT],
        metadata={"cross_modal_layer": "limitations",
                  "backed_by_concordance_engine": engine},
    )


def compose_cross_modal_discussion_blocks(
    patterns: dict[str, Any],
    agent_results: dict | None = None,
) -> list[NarrativeBlock]:
    """Governed ``integration.cross_modal_*`` blocks for the active patterns.

    Only layers that actually ran and carry data produce a block; the mandatory
    limitations block is appended whenever any cross-modal claim was made.
    ``agent_results`` is accepted for API symmetry with the other composers.
    """
    if not isinstance(patterns, dict) or not patterns.get("active"):
        return []
    engine = bool(patterns.get("backed_by_concordance_engine"))
    route = _route(patterns)
    contrasts = patterns.get("contrasts") or []
    multi = len(contrasts) > 1
    blocks: list[NarrativeBlock] = []
    for mc in contrasts:
        for builder in (_de_da_block, _tf_block, _pathway_block):
            block = builder(route, mc, multi, engine)
            if block is not None:
                blocks.append(block)
    links = patterns.get("sample_links")
    if engine and isinstance(links, dict):
        block = _sample_links_block(route, links)
        if block is not None:
            blocks.append(block)
    if not blocks:
        return []
    primary = next(
        (b for b in blocks if b.id.startswith("integration.cross_modal_de_da")),
        blocks[0],
    )
    if patterns.get("figures"):
        primary.figures = [dict(f) for f in patterns["figures"]]
    if patterns.get("tables"):
        primary.tables = [dict(t) for t in patterns["tables"]]
    blocks.append(_limitations_block(route, patterns, engine))
    return blocks
