"""Deterministic associative concordance for independently assayed RNA and ATAC.

Contrast effects, genomic proximity, and sample correlations are distinct evidence
layers. A cis assignment is a positional association, never a molecular mechanism.
Missing effects or significance values are retained as missing, never invented.
"""
from __future__ import annotations

import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, pearsonr, rankdata, spearmanr, t

from aria.benchmarks.bulk_atac_marker_mapping import parse_peak_interval
from aria.scripts.chromatin_peak_annotation import annotate_peaks, parse_gtf_gene_models

# Pure-stdlib mirror of aria.agents.narrative.validators.CAUSAL_PATTERNS so
# subprocess execution inside isolated scientific stacks (aria-rna-env,
# aria-chromatin-env) never imports aria.agents -> litellm.
CAUSAL_PATTERNS = (
    "drives", "enforces", "acts as", "governs", "position as", "positions as",
    "directly regulates", "binds to", "hierarchical gatekeeper",
    "hierarchical gatekeepers", "cause", "causes", "caused by", "causally",
    "induces", "induced by",
    "triggers", "leads to", "results in", "gives rise to", "controls",
    "regulates", "drives the", "promotes", "inhibits", "suppresses",
    "activates", "represses", "responsible for", "master regulator",
    "master regulators", "upstream of", "downstream effector", "orchestrates",
    "dictates", "determines", "mediates",
)
_CAUSAL_RE = re.compile(
    r"\b(" + "|".join(re.escape(p) for p in CAUSAL_PATTERNS) + r")\b",
    flags=re.IGNORECASE,
)


def find_causal_language(text: str) -> str | None:
    """Return the first causal pattern found in ``text``, or None."""
    if not text:
        return None
    match = _CAUSAL_RE.search(text)
    return match.group(0) if match else None

CATEGORIES = (
    "concordant_up", "concordant_down", "discordant_rna_up_atac_down",
    "discordant_rna_down_atac_up", "atac_da_rna_unchanged",
    "rna_de_atac_unchanged", "both_unchanged",
)


def finite(value):
    """A JSON-compatible finite number, or None."""
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError):
        return None


def safe_label(value):
    """Remove database identifiers and replace assertion-like label fragments."""
    label = re.sub(r"\s*\(?GO:\d+\)?|\s*\(?R-HSA-\d+\)?|\s*\(?hsa\d+\)?", "", str(value), flags=re.I).strip()
    while (word := find_causal_language(label)) is not None:
        label = re.sub(re.escape(word), "association", label, flags=re.I)
    return label


def unwrap(value):
    """Peel agent findings envelopes without guessing absent biological labels."""
    value = value if isinstance(value, dict) else {}
    while isinstance(value.get("findings"), dict):
        value = value["findings"]
    return value


def normalize_inputs(rna_results, chromatin_results):
    rna, atac = unwrap(rna_results), unwrap(chromatin_results)
    rm = None
    for key, modality in (("bulk_RNA", "bulk_RNA"), ("scRNA", "scRNA")):
        if isinstance(rna.get(key), dict):
            rna, rm = unwrap(rna[key]), modality
            break
    rm = rm or ("scRNA" if "pseudobulk_de" in rna or "per_group" in rna else "bulk_RNA" if "contrasts" in rna else None)
    am = None
    for key in ("bulk_ATAC", "scATAC"):
        if isinstance(atac.get(key), dict):
            atac, am = unwrap(atac[key]), key
            break
    da = unwrap(atac.get("differential_accessibility", atac))
    am = am or (atac.get("data_type") if atac.get("data_type") in ("bulk_ATAC", "scATAC") else None)
    am = am or ("scATAC" if "pseudobulk" in da or "per_group" in da else "bulk_ATAC" if "comparisons" in da else None)
    return rna, atac, rm, am


def contrast_arms(record):
    numerator = record.get("test") or record.get("numerator")
    denominator = record.get("reference") or record.get("denominator") or record.get("ref")
    name = record.get("name") or record.get("comparison") or record.get("contrast") or record.get("contrast_id")
    if (not numerator or not denominator) and isinstance(name, str) and "_vs_" in name:
        left, right = name.split("_vs_", 1)
        numerator, denominator = numerator or left, denominator or right
    return (str(numerator), str(denominator)) if numerator and denominator and numerator != denominator else None


def contrast_name(record):
    arms = contrast_arms(record)
    return str(record.get("name") or record.get("comparison") or ("_vs_".join(arms) if arms else ""))


def _eligible(record):
    return isinstance(record, dict) and record.get("status") not in ("skipped", "error") and record.get("ran") is not False


def rna_contrasts(rna):
    if "contrasts" in rna:
        return [dict(x) for x in rna.get("contrasts", []) if _eligible(x)]
    pb = unwrap(rna.get("pseudobulk_de", rna))
    rows = []
    for group, data in (pb.get("per_group") or {}).items():
        for name, comp in (data.get("per_comparison") or {}).items():
            if _eligible(comp):
                rows.append({
                    **comp,
                    "name": name,
                    "group": str(group),
                    "output_csv": comp.get("output_csv") or pb.get("output_csv"),
                    "full_results_csv": comp.get("full_results_csv") or pb.get("full_results_csv"),
                })
    return rows


def atac_contrasts(atac):
    da = unwrap(atac.get("differential_accessibility", atac))
    da = unwrap(da.get("pseudobulk", da))
    rows = []
    for comp in da.get("comparisons") or []:
        if _eligible(comp):
            rows.append({**comp, "full_results_csv": comp.get("full_results_csv") or da.get("full_results_csv")})
    for group, data in (da.get("per_group") or {}).items():
        for comp in data.get("comparisons") or []:
            if _eligible(comp):
                rows.append({**comp, "group": str(group), "full_results_csv": comp.get("full_results_csv") or data.get("full_results_csv")})
        for name, comp in (data.get("per_comparison") or {}).items():
            if _eligible(comp):
                rows.append({**comp, "name": name, "group": str(group), "full_results_csv": comp.get("full_results_csv") or data.get("full_results_csv")})
    return rows


def align_contrasts(rna, atac, warnings):
    """Exact arm and optional cell-group matching, including reversed effects."""
    pairs = []
    rc, ac = rna_contrasts(rna), atac_contrasts(atac)
    stratified = any(x.get("group") is not None for x in ac)
    for r in rc:
        ra = contrast_arms(r)
        matched = False
        for a in ac:
            aa = contrast_arms(a)
            if not ra or not aa or (stratified and r.get("group") != a.get("group")):
                continue
            if aa not in (ra, ra[::-1]):
                continue
            granularity = ("cell_group_matched" if stratified else "condition_pseudobulk_vs_cell_groups" if r.get("group") is not None else "condition_matched")
            pairs.append((r, a, aa == ra[::-1], granularity))
            matched = True
        if not matched:
            warnings.append("contrast_skipped: incompatible_contrast_denominators_or_cell_groups")
    return pairs, bool(rc and ac)


def read_table(path, warnings):
    if not path or not Path(path).is_file():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t" if str(path).endswith((".tsv", ".txt", ".tsv.gz")) else ",")
    except (OSError, ValueError, pd.errors.ParserError):
        warnings.append("input_table_unreadable")
        return pd.DataFrame()


def _column(df, *names):
    lookup = {str(c).lower(): c for c in df.columns}
    return next((lookup[n.lower()] for n in names if n.lower() in lookup), None)


def _value(record, *names):
    lookup = {str(key).lower(): value for key, value in record.items()}
    for name in names:
        value = lookup.get(name.lower())
        if value is None or (isinstance(value, str) and not value.strip()):
            continue
        if isinstance(value, (bool, np.bool_)):
            continue
        try:
            if not math.isfinite(float(value)):
                continue
        except (TypeError, ValueError):
            if not isinstance(value, str) or pd.isna(value):
                continue
        return value
    return None


def filter_comparison(df, comp):
    for col in ("comparison", "contrast"):
        if col in df:
            df = df[df[col].astype(str) == contrast_name(comp)]
    if "group" in df and comp.get("group") is not None:
        df = df[df["group"].astype(str) == str(comp["group"])]
    return df


def rna_table_paths(comp):
    paths = [comp.get("full_results_csv")]
    if comp.get("tables_dir"):
        paths.append(str(Path(comp["tables_dir"]) / "de_all_genes.tsv"))
    paths.append(((comp.get("plots") or {}).get("tables") or {}).get("de_genes"))
    if comp.get("tables_dir"):
        paths.append(str(Path(comp["tables_dir"]) / "de_genes.tsv"))
    paths.append(comp.get("output_csv"))
    return paths


def gtf_symbols(path):
    """ID-to-symbol aliases from the same GTF as the genomic assignments."""
    aliases = {}
    if not path:
        return aliases
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                gid = re.search(r'gene_id "([^"]+)"', fields[8])
                symbol = re.search(r'gene_name "([^"]+)"', fields[8])
                if gid and symbol:
                    aliases[gid[1].split(".")[0].upper()] = symbol[1]
    except (OSError, EOFError):
        return {}
    return aliases


def load_rna_table(comp, aliases, warnings, alpha=.05):
    aliases = dict(aliases)
    for row in comp.get("top_genes") or []:
        if isinstance(row, dict) and row.get("gene") and row.get("symbol"):
            aliases[str(row["gene"]).split(".")[0].upper()] = str(row["symbol"])
    df = pd.DataFrame()
    full = False
    for path in rna_table_paths(comp):
        if path and Path(path).is_file():
            df = filter_comparison(read_table(path, warnings), comp)
            if df.empty:
                continue
            # ARIA pseudobulk output_csv contains significant rows only.
            full = not (path == comp.get("output_csv") and comp.get("group") is not None)
            break
    if df.empty:
        records = []
        for key in ("top_genes", "all_sig", "all_sig_genes"):
            for row in comp.get(key) or []:
                records.append({**row, "_known_sig": key != "top_genes"} if isinstance(row, dict) else {"gene": str(row), "_known_sig": key != "top_genes"})
        for key, sign in (("up_gene_ids", 1), ("down_gene_ids", -1)):
            records.extend({"gene": str(g), "_known_sig": True, "_direction_sign": sign} for g in comp.get(key) or [])
        df = pd.DataFrame(records)
        full = False
    # export_de_table writes significant rows, despite the generic TSV name.
    # Some callers supply an extended all-tested TSV under the same name; a
    # non-significant adjusted p-value establishes that extended scope.
    if full and path != comp.get("full_results_csv"):
        pcol = _column(df, "padj", "FDR", "adjusted_pvalue", "Adjusted P-value")
        if pcol is not None and not (pd.to_numeric(df[pcol], errors="coerce") >= alpha).any():
            full = False
    if not full:
        warnings.append("rna_scope: reported_genes_only; unreported_genes_are_not_tested_background")
    gene_col = _column(df, "gene_symbol", "symbol", "hgnc_symbol", "HGNC Symbol", "gene", "gene_id", "Gene ID", "Unnamed: 0")
    if gene_col is None:
        return pd.DataFrame(), "all_tested_genes" if full else "reported_genes_only"
    rows = []
    for raw in df.to_dict("records"):
        gene = next((raw[c] for c in (gene_col, "symbol", "gene", "gene_id", "Gene ID", "Unnamed: 0")
                     if c in raw and pd.notna(raw[c]) and str(raw[c]).strip()), None)
        if gene is None:
            continue
        gene = str(gene)
        gene = aliases.get(gene.split(".")[0].upper(), gene)
        lfc = finite(_value(raw, "log2fc", "log2FoldChange", "log2_fold_change", "Log2 Fold Change"))
        padj = finite(_value(raw, "padj", "FDR", "adjusted_pvalue", "Adjusted P-value"))
        direction = finite(raw.get("_direction_sign"))
        if direction is None and raw.get("direction") in ("up", "down"):
            direction = 1 if raw["direction"] == "up" else -1
        raw_gid = next(
            (str(raw[c]).split(".")[0].strip() for c in ("gene_id", "Gene ID")
             if c in raw and pd.notna(raw[c]) and str(raw[c]).strip()),
            None,
        )
        rows.append({"gene": gene, "gene_key": gene.upper(), "rna_log2fc": lfc, "rna_padj": padj,
                     "rna_direction": int(np.sign(lfc)) if lfc is not None else direction,
                     "rna_known_sig": raw.get("_known_sig") is True,
                     "best_cell_group": comp.get("group"),
                     "_locus_id": raw_gid})
    result = pd.DataFrame(rows)
    if not result.empty:
        # Row-preserving deduplication per gene_key (D4-011): never splice
        # columns across distinct loci via groupby().first() and never OR
        # rna_known_sig across opposing directions or finite non-sig padj.
        deduped_rows = []
        ambiguous_keys: list[str] = []
        dup_locus_keys: list[str] = []
        for gkey, grp in result.groupby("gene_key", sort=False):
            if len(grp) == 1:
                deduped_rows.append(grp.iloc[0].to_dict())
                continue
            ordered = grp.sort_values("rna_padj", na_position="last", kind="mergesort")
            dirs = {
                int(d) for d in ordered["rna_direction"].dropna().tolist()
                if int(d) != 0
            }
            loci = {
                str(x) for x in ordered["_locus_id"].dropna().tolist()
                if str(x).strip()
            }
            finite_padjs = ordered["rna_padj"].dropna().tolist()
            has_sig = bool(
                any(p < alpha for p in finite_padjs)
                or ordered["rna_known_sig"].any()
            )
            has_nonsig = bool(any(p >= alpha for p in finite_padjs))
            if len(dirs) > 1 or (has_sig and has_nonsig):
                ambiguous_keys.append(str(gkey))
                continue
            if len(loci) > 1:
                dup_locus_keys.append(str(gkey))
            winner = ordered.iloc[0].to_dict()
            # Only propagate rna_known_sig from compatible in-memory lists when
            # the winning row itself has no conflicting non-significant padj.
            w_padj = winner.get("rna_padj")
            if (
                ordered["rna_known_sig"].any()
                and (w_padj is None or pd.isna(w_padj) or float(w_padj) < alpha)
            ):
                winner["rna_known_sig"] = True
            deduped_rows.append(winner)
        if ambiguous_keys:
            warnings.append(
                f"rna_ambiguous_symbol_collisions_excluded: n={len(ambiguous_keys)} "
                f"({', '.join(sorted(ambiguous_keys)[:10])})"
            )
        if dup_locus_keys:
            warnings.append(
                f"rna_duplicate_symbol_loci_deduplicated: n={len(dup_locus_keys)} "
                f"({', '.join(sorted(dup_locus_keys)[:10])})"
            )
        cols = ["gene", "gene_key", "rna_log2fc", "rna_padj", "rna_direction",
                "rna_known_sig", "best_cell_group"]
        result = (
            pd.DataFrame(deduped_rows)[cols]
            if deduped_rows
            else pd.DataFrame(columns=cols)
        )
    return result, "all_tested_genes" if full else "reported_genes_only"


def _annotation_records(atac, comp, warnings):
    ann = unwrap(atac.get("peak_annotation", {}))
    paths = [ann.get("annotation_csv")]
    for entry in ann.get("comparisons") or []:
        if contrast_arms(entry) == contrast_arms(comp) and (entry.get("group") is None or entry.get("group") == comp.get("group")):
            paths.append(entry.get("annotation_csv"))
    frames = [read_table(p, warnings) for p in paths if p]
    return pd.concat(frames, ignore_index=True).drop_duplicates("peak") if frames and "peak" in pd.concat(frames).columns else pd.DataFrame()


def load_atac_table(comp, atac, models, params, flipped, warnings):
    path = comp.get("full_results_csv")
    full = bool(path and Path(path).is_file())
    df = filter_comparison(read_table(path, warnings), comp) if full else pd.DataFrame(comp.get("top_peaks") or [])
    if "peak" not in df or df.empty:
        return pd.DataFrame(), "unavailable"
    df = df.drop_duplicates("peak").copy()
    if models:
        annotations = pd.DataFrame(annotate_peaks(df["peak"].astype(str), models,
                                    params["promoter_bp"], params["promoter_bp"])["rows"])
        scope = "all_tested_peaks" if full else "reported_peaks_only"
        df = df.drop(columns=[c for c in ("nearest_gene", "distance_to_tss", "feature") if c in df])
        df = df.merge(annotations, on="peak", how="inner") if not annotations.empty else pd.DataFrame()
    elif {"nearest_gene", "distance_to_tss"}.issubset(df.columns):
        scope = "all_tested_peaks" if full else "reported_peaks_only"
    else:
        annotations = _annotation_records(atac, comp, warnings)
        df = df.merge(annotations, on="peak", how="inner", suffixes=("_input", "")) if not annotations.empty else pd.DataFrame()
        scope = "significant_peaks_only"
        warnings.append("annotation_scope: significant_peaks_only; continuous_associations_have_selected_peak_background")
    rows = []
    for raw in df.to_dict("records"):
        gene = raw.get("nearest_gene")
        distance = finite(raw.get("distance_to_tss"))
        lfc = finite(_value(raw, "log2fc", "log2FoldChange"))
        padj = finite(raw.get("padj"))
        if not isinstance(gene, str) or not gene or distance is None or abs(distance) > params["cis_distance_bp"]:
            continue
        conv = raw.get("lfc_converged")
        unconverged = (isinstance(conv, (int, float, np.number)) and conv == 0) or str(conv).strip().lower() in ("false", "0", "0.0")
        tested = padj is not None and lfc is not None and not unconverged
        if not tested:
            continue
        rows.append({"peak": str(raw["peak"]), "gene_key": gene.upper(), "atac_log2fc": -lfc if flipped else lfc,
                     "atac_padj": padj, "distance_to_tss": distance, "feature": raw.get("feature"),
                     "atac_direction": ("up" if (-lfc if flipped else lfc) > 0 else "down" if lfc != 0 else "unchanged")})
    return pd.DataFrame(rows), scope


def _category(row, alpha):
    rs = (finite(row.get("rna_padj")) is not None and row["rna_padj"] < alpha) or row.get("rna_known_sig") is True
    ats = row["atac_padj"] < alpha
    rd, ad = row.get("rna_direction"), np.sign(row["atac_log2fc"])
    if rs and ats:
        if rd == 1 and ad == 1:
            return "concordant_up"
        if rd == -1 and ad == -1:
            return "concordant_down"
        if rd == 1 and ad == -1:
            return "discordant_rna_up_atac_down"
        if rd == -1 and ad == 1:
            return "discordant_rna_down_atac_up"
        return "both_significant_direction_unavailable"
    return "rna_de_atac_unchanged" if rs else "atac_da_rna_unchanged" if ats else "both_unchanged"


def _correlations(frame):
    result = {k: None for k in ("spearman_r", "spearman_pvalue", "pearson_r", "pearson_pvalue")}
    if frame.empty:
        return result
    valid = frame[["rna_log2fc", "atac_log2fc"]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(valid) < 3 or (valid.nunique() < 2).any():
        return result
    for name, func in (("spearman", spearmanr), ("pearson", pearsonr)):
        stat, p = func(valid.iloc[:, 0], valid.iloc[:, 1])
        result[f"{name}_r"], result[f"{name}_pvalue"] = finite(stat), finite(p)
    return result


def _counts(frame):
    counts = Counter(frame["classification"]) if not frame.empty else Counter()
    out = {f"n_{c}": int(counts[c]) for c in CATEGORIES}
    out["n_concordant_total"] = counts["concordant_up"] + counts["concordant_down"]
    out["n_discordant_total"] = counts["discordant_rna_up_atac_down"] + counts["discordant_rna_down_atac_up"]
    out["n_both_significant"] = out["n_concordant_total"] + out["n_discordant_total"] + counts["both_significant_direction_unavailable"]
    return out


def summarize_de_da(frame, promoter_bp=2000):
    out = {"ran": not frame.empty, "n_pairs_tested": len(frame), "n_gene_peak_pairs_tested": len(frame),
           "n_promoter_pairs_tested": int((frame.distance_to_tss.abs() <= promoter_bp).sum()) if not frame.empty else 0,
           "n_genes_tested": 0,
           "gene_counting_unit": "gene_by_matched_contrast",
           **_counts(frame), "concordance_fraction": None, "fisher_odds_ratio": None, "fisher_pvalue": None,
           "enrichment_odds_ratio": None, "enrichment_pvalue": None, **_correlations(pd.DataFrame()),
           "promoter_spearman_r": None, "promoter_spearman_pvalue": None,
           "top_concordant_genes": [], "top_discordant_genes": []}
    if frame.empty:
        out["gene_counts"] = _counts(frame)
        return out
    keys = [c for c in ("contrast_id", "gene_key") if c in frame]
    sorted_frame = frame.assign(_distance=frame.distance_to_tss.abs(), _effect=frame.atac_log2fc.abs()).sort_values(
        ["atac_padj", "_distance", "_effect", "peak"], ascending=[True, True, False, True])
    primary = sorted_frame.drop_duplicates(keys) if keys else sorted_frame
    gene_counts = _counts(primary)
    out["n_genes_tested"] = len(primary)
    out["gene_counts"] = gene_counts
    out.update(_correlations(primary))
    # Enrichment is based on any DA peak per gene, not the representative peak.
    gene_sig = frame.groupby(keys)[["rna_significant", "atac_significant"]].any()
    both = int((gene_sig.rna_significant & gene_sig.atac_significant).sum())
    rna_only = int((gene_sig.rna_significant & ~gene_sig.atac_significant).sum())
    atac_only = int((~gene_sig.rna_significant & gene_sig.atac_significant).sum())
    neither = int((~gene_sig.rna_significant & ~gene_sig.atac_significant).sum())
    out["gene_enrichment_contingency"] = [[both, rna_only], [atac_only, neither]]
    if neither > 0:
        odds, p = fisher_exact(out["gene_enrichment_contingency"], alternative="two-sided")
        out.update(enrichment_odds_ratio=finite(odds), enrichment_pvalue=finite(p))
    if gene_counts["n_both_significant"]:
        out["concordance_fraction"] = gene_counts["n_concordant_total"] / gene_counts["n_both_significant"]
    if gene_counts["n_concordant_total"] + gene_counts["n_discordant_total"] >= 2:
        odds, p = fisher_exact([[gene_counts["n_concordant_up"], gene_counts["n_discordant_rna_up_atac_down"]],
                               [gene_counts["n_discordant_rna_down_atac_up"], gene_counts["n_concordant_down"]]], alternative="two-sided")
        out.update(fisher_odds_ratio=finite(odds), fisher_pvalue=finite(p))
    promoter_sorted = sorted_frame[sorted_frame.distance_to_tss.abs() <= promoter_bp]
    promoter_primary = promoter_sorted.drop_duplicates(keys) if (keys and not promoter_sorted.empty) else promoter_sorted
    promoter = _correlations(promoter_primary)
    out.update(promoter_spearman_r=promoter["spearman_r"], promoter_spearman_pvalue=promoter["spearman_pvalue"])
    for kind in ("concordant", "discordant"):
        subset = primary[primary.classification.str.startswith(kind)].sort_values(["rna_padj", "atac_padj"], na_position="last")
        out[f"top_{kind}_genes"] = subset[["gene", "rna_log2fc", "atac_log2fc", "classification", "best_cell_group"]].head(20).to_dict("records")
    return out


def de_da_layer(rna_df, atac_df, params, contrast_id):
    if rna_df.empty or atac_df.empty:
        frame = pd.DataFrame()
    else:
        # NaN adjusted p-values from independently filtered genes are not an
        # unchanged comparison background.
        eligible = rna_df.rna_padj.notna() | rna_df.rna_known_sig.fillna(False)
        frame = rna_df[eligible].merge(atac_df, on="gene_key", how="inner")
    if not frame.empty:
        alpha = params["padj_threshold"]
        frame["classification"] = [_category(x, alpha) for x in frame.to_dict("records")]
        frame["rna_significant"] = (frame.rna_padj < alpha) | frame.rna_known_sig.fillna(False)
        frame["atac_significant"] = frame.atac_padj < alpha
        frame["contrast_id"] = contrast_id
        rank = frame.assign(_dist=frame.distance_to_tss.abs(), _effect=frame.atac_log2fc.abs()).sort_values(
            ["atac_padj", "_dist", "_effect", "peak"], ascending=[True, True, False, True])
        frame["is_primary_peak"] = frame.index.isin(rank.drop_duplicates("gene_key").index)
    return frame, summarize_de_da(frame, params["promoter_bp"])


def tf_subunits(name):
    name = re.sub(r"\([^)]*\)", "", str(name)).strip()
    name = re.sub(r"^MA\d+\.\d+\s+", "", name)
    return [x.strip() for x in name.split("::") if x.strip()]


def _stem(tf):
    if re.match(r"^(FOS|JUN)", tf, flags=re.I):
        return "AP1"
    return re.sub(r"[\d]+[A-Z]?$", "", tf.upper())


def _group_direction(group, comp):
    """Only explicit directional groups with the exact contrast are eligible."""
    group = str(group)
    if group in ("up", "down"):
        return 1 if group == "up" else -1
    arms = contrast_arms(comp)
    names = {contrast_name(comp)} | ({"_vs_".join(arms)} if arms else set())
    if "::up_in_" in group and arms:
        prefix, level = group.split("::up_in_", 1)
        valid_prefixes = names | ({f"{comp['group']}_{n}" for n in names} if comp.get("group") is not None else set())
        if prefix in valid_prefixes:
            if level == arms[0]:
                return 1
            if level == arms[1]:
                return -1
    for sign, suffix in ((1, "_up"), (-1, "_down")):
        prefix = group[:-len(suffix)] if group.endswith(suffix) else None
        if prefix in names or (comp.get("group") is not None and prefix in {f"{comp['group']}_{name}" for name in names}):
            return sign
    return None


def _footprint_records(foot, comp, warnings):
    if not foot.get("ran"):
        return [], 1
    ga, gb = foot.get("group_a"), foot.get("group_b")
    arms = contrast_arms(comp)
    foot_arms = (str(ga), str(gb)) if ga and gb else contrast_arms(foot)
    if arms is None or foot_arms is None:
        warnings.append("footprint_skipped: contrast_orientation_unavailable")
        return [], 1
    if foot_arms not in (arms, arms[::-1]):
        warnings.append("footprint_skipped: incompatible_contrast_denominators")
        return [], 1
    orientation = -1 if foot_arms and foot_arms == arms[::-1] else 1
    summary = foot.get("differential_summary") or {}
    if summary.get("inference", {}).get("status") == "descriptive_only" or summary.get("ranking_basis", {}).get("fdr_controlled") is False:
        warnings.append("footprint_skipped: descriptive_ranking_without_replicate_inference")
        return [], orientation
    rows = []
    for container in (foot, summary):
        for key in ("significant_tfs", "top_differential_motifs", "all_results"):
            rows.extend(x for x in container.get(key) or [] if isinstance(x, dict))
        for key, values in container.items():
            if key.startswith("top_toward_") and isinstance(values, list):
                rows.extend(x for x in values if isinstance(x, dict))
    if not rows:
        path = foot.get("bindetect_results_tsv") or (foot.get("tables") or {}).get("bindetect_results_tsv")
        rows = read_table(path, warnings).to_dict("records")
    normalized = []
    for row in rows:
        row = dict(row)
        # Raw BINDetect per-site p-values are not biological-replicate FDR.
        change_key = f"{ga}_{gb}_change"
        padj = finite(_value(row, "padj", "qvalue", "fdr", f"{ga}_{gb}_padj"))
        change = finite(_value(row, "change", "log2fc", "log2FoldChange", "score", change_key))
        if padj is not None and change is not None:
            normalized.append({**row, "padj": padj, "change": change})
    return normalized, orientation


def tf_layer(atac, comp, rna_df, flipped, alpha, warnings, expression=None):
    motifs = unwrap(atac.get("motifs", {}))
    motif_rows = read_table(motifs.get("output_csv"), warnings).to_dict("records")
    if not motif_rows:
        for group, data in (motifs.get("per_group") or {}).items():
            motif_rows.extend({**row, "group": group} for row in data.get("top_motifs") or [])
    evidence = []
    for motif in motif_rows:
        padj = finite(motif.get("padj"))
        direction = _group_direction(motif.get("group", motif.get("direction", "")), comp)
        if padj is None or not 0.0 <= padj < alpha or direction is None:
            continue
        score = finite(motif.get("log2_enrichment"))
        if score is not None and score <= 0:
            continue
        evidence.append({**motif, "source": "motif", "direction_sign": -direction if flipped else direction,
                         # Zero FDR alone provides no finite enrichment magnitude.
                         "chromatin_score": score if score is not None else -math.log10(padj) if padj > 0 else None})
    foot_rows, orientation = _footprint_records(unwrap(atac.get("footprinting", {})), comp, warnings)
    for row in foot_rows:
        if 0.0 <= row["padj"] < alpha and row["change"] != 0:
            change = row["change"] * orientation * (-1 if flipped else 1)
            evidence.append({**row, "name": row.get("tf") or row.get("name"), "source": "footprint",
                             "direction_sign": int(np.sign(change)), "chromatin_score": abs(change)})
    lookup = rna_df.set_index("gene_key").to_dict("index") if not rna_df.empty else {}
    families = defaultdict(set)
    stems = defaultdict(set)
    expanded = []
    for item in evidence:
        name = item.get("name") or item.get("tf")
        if not name:
            continue
        subunits = tf_subunits(name)
        signs = {lookup[x.upper()].get("rna_direction") for x in subunits if x.upper() in lookup}
        divergent = 1 in signs and -1 in signs
        for tf in subunits:
            families[str(item.get("family") or "").lower()].add(tf.upper())
            stems[_stem(tf)].add(tf.upper())
            expanded.append({**item, "tf": tf, "motif_name": str(name), "dimer_direction_divergence": divergent})
    rows = []
    for item in expanded:
        tf = item["tf"]
        rna = lookup.get(tf.upper())
        expr = (expression or {}).get(tf.upper())
        rna_padj = finite(rna.get("rna_padj")) if rna is not None else None
        rna_known_sig = bool(rna.get("rna_known_sig") is True) if rna is not None else False
        significant = rna is not None and ((rna_padj is not None and rna_padj < alpha) or rna_known_sig)
        if rna is None and not (expr is not None and expr > 0):
            classification = "not_detected_in_rna"
        elif rna is None or (rna_padj is None and not rna_known_sig):
            classification = "rna_not_tested"
        elif significant and rna.get("rna_direction") in (-1, 1):
            classification = "concordant" if rna["rna_direction"] == item["direction_sign"] else "discordant"
        else:
            classification = "chromatin_shift_rna_unchanged"
        family = str(item.get("family") or "").lower()
        ambiguous = item["dimer_direction_divergence"] or len(stems[_stem(tf)]) > 1 or bool(family and len(families[family]) > 1)
        rows.append({"tf": tf, "motif_name": item["motif_name"], "motif_id": item.get("motif_id"),
                     "family": item.get("family"), "source": item["source"], "chromatin_direction": "up" if item["direction_sign"] > 0 else "down",
                     "chromatin_score": item["chromatin_score"] * item["direction_sign"] if item["chromatin_score"] is not None else None, "chromatin_padj": item["padj"],
                     "rna_log2fc": rna.get("rna_log2fc") if rna else None, "rna_padj": rna.get("rna_padj") if rna else None,
                     "rna_mean_cpm": expr, "classification": classification, "paralog_motif_ambiguity": bool(ambiguous)})
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["chromatin_score"] = pd.to_numeric(frame["chromatin_score"], errors="coerce")
        frame = frame.drop_duplicates(["tf", "source", "motif_name", "chromatin_direction"])
    return frame, summarize_tfs(frame)


def summarize_tfs(frame):
    out = {"ran": not frame.empty, "n_tfs_evaluated": 0, "n_concordant_tfs": 0, "n_discordant_tfs": 0,
           "n_chromatin_shift_rna_unchanged": 0, "n_paralog_ambiguous": 0, "top_concordant_tfs": []}
    if frame.empty:
        return out
    work = frame.copy()
    work["_tf_key"] = work["tf"].astype(str).str.upper()
    work["_abs_score"] = pd.to_numeric(work.get("chromatin_score"), errors="coerce").abs().fillna(-1.0)
    keys = [x for x in ("contrast_id", "_tf_key") if x in work]
    unique = (
        work.sort_values(
            ["chromatin_padj", "_abs_score", "_tf_key"],
            ascending=[True, False, True],
            kind="stable",
        )
        .drop_duplicates(keys)
    )
    counts = Counter(unique.classification)
    drop_cols = ["_tf_key", "_abs_score"]
    top_df = unique[unique.classification == "concordant"].drop(columns=drop_cols).head(20)
    out.update(n_tfs_evaluated=len(unique), n_concordant_tfs=counts["concordant"], n_discordant_tfs=counts["discordant"],
               n_chromatin_shift_rna_unchanged=counts["chromatin_shift_rna_unchanged"],
               n_paralog_ambiguous=int(work.groupby(keys).paralog_motif_ambiguity.any().sum()),
               top_concordant_tfs=[
                   {key: None if isinstance(value, float) and not math.isfinite(value) else value
                    for key, value in row.items()}
                   for row in top_df.to_dict("records")])
    return out


def pathway_terms(value, alpha):
    """Flatten ARIA ORA dictionaries, retaining only finite significant terms."""
    terms = {}
    def visit(node):
        if isinstance(node, list):
            for child in node:
                visit(child)
        elif isinstance(node, dict):
            raw = _value(node, "term", "Term", "name")
            padj = finite(_value(node, "padj", "Adjusted P-value", "adjusted_pvalue"))
            if raw is not None and padj is not None:
                label = safe_label(raw)
                # Match the database term after identifier stripping, before
                # display sanitization, so unrelated labels cannot collapse.
                clean = re.sub(r"\s*\(?GO:\d+\)?|\s*\(?R-HSA-\d+\)?|\s*\(?hsa\d+\)?", "", str(raw), flags=re.I).strip()
                key = re.sub(r"\s+", " ", clean).casefold()
                if key and 0.0 <= padj < alpha and (key not in terms or padj < terms[key]["padj"]):
                    terms[key] = {"term": str(raw), "safe_label": label, "padj": padj}
            else:
                for key, child in node.items():
                    if key not in ("warnings", "figures", "gene_set_versions", "background_genes"):
                        visit(child)
    visit(value)
    return terms


def pathway_layer(rna, atac, rcomp, acomp, alpha):
    rp = rcomp.get("pathways")
    if rp is None and rcomp.get("group") is not None:
        rp = (unwrap(rna.get("pathways", {})).get("per_cluster") or {}).get(str(rcomp["group"]))
    ap = None
    ora = unwrap(atac.get("peak_ora", {}))
    arms = contrast_arms(acomp)
    for comp in ora.get("comparisons") or []:
        comp_arms = contrast_arms(comp)
        if arms is not None and comp_arms is not None and comp_arms in (arms, arms[::-1]):
            if comp.get("group") is None or comp.get("group") == acomp.get("group"):
                ap = comp.get("pathways")
                break
    if ap is None and "comparisons" not in ora:
        ap = ora.get("pathways")
    rt, at = pathway_terms(rp, alpha), pathway_terms(ap, alpha)
    shared = []
    for key in rt.keys() & at.keys():
        rr, aa = rt[key], at[key]
        shared.append({"term": rr["term"], "atac_term": aa["term"], "safe_label": rr["safe_label"],
                       "rna_padj": rr["padj"], "atac_padj": aa["padj"],
                       "geometric_mean_padj": math.sqrt(rr["padj"]) * math.sqrt(aa["padj"]), "term_key": key})
    shared.sort(key=_pathway_sort_key)
    union = rt.keys() | at.keys()
    summary = {"ran": bool(rp) and bool(ap), "n_rna_terms": len(rt), "n_atac_terms": len(at),
               "n_shared_terms": len(shared), "jaccard_index": len(shared) / len(union) if union else 0.0,
               "shared_terms": [{**{k: v for k, v in row.items() if k not in ("term", "atac_term", "term_key")},
                                  "term": row["safe_label"]} for row in shared],
               "rna_only_top_terms": sorted([rt[k]["safe_label"] for k in rt.keys() - at.keys()], key=str.casefold)[:10],
               "atac_only_top_terms": sorted([at[k]["safe_label"] for k in at.keys() - rt.keys()], key=str.casefold)[:10]}
    return pd.DataFrame(shared), summary, rt, at


def _pathway_sort_key(row):
    larger, smaller = max(row["rna_padj"], row["atac_padj"]), min(row["rna_padj"], row["atac_padj"])
    term_key = row.get("term_key")
    if term_key is None:
        term_key = str(row.get("safe_label") or row.get("term") or "").casefold()
    # A reported zero has unknown precision; rank it conservatively without
    # changing the measured p-values or their geometric mean.
    return (row["geometric_mean_padj"] if smaller > 0 else larger, larger, smaller, term_key)


def aggregate_pathways(rna_terms, atac_terms, frame, ran):
    if frame.empty:
        rows = []
    else:
        rows = sorted(frame.to_dict("records"), key=_pathway_sort_key)
        rows = pd.DataFrame(rows).drop_duplicates("term_key").to_dict("records")
    shared = [{**{k: v for k, v in row.items() if k not in ("term", "atac_term", "term_key")},
               "term": row["safe_label"]} for row in rows]
    union = set(rna_terms) | set(atac_terms)
    return {"ran": ran, "n_rna_terms": len(rna_terms), "n_atac_terms": len(atac_terms),
            "n_shared_terms": len(shared), "jaccard_index": len(shared) / len(union) if union else 0.0, "shared_terms": shared}


def resolve_cpm_paths(params, rna, atac):
    rpath, apath = params.get("rna_cpm_path"), params.get("atac_cpm_path")
    def search(node, keys):
        if isinstance(node, dict):
            for key in keys:
                value = node.get(key)
                if isinstance(value, (str, Path)) and Path(value).is_file():
                    return str(value)
            for value in node.values():
                found = search(value, keys)
                if found:
                    return found
        elif isinstance(node, list):
            for value in node:
                found = search(value, keys)
                if found:
                    return found
        return None
    keys = ("counts_cpm_path", "cpm_path", "pseudobulk_cpm_path", "donor_cpm_path", "counts_tpm_path")
    rpath = rpath or search(rna, keys)
    apath = apath or search(atac.get("peak_counts", atac), keys)
    if not rpath:
        for comp in rna_contrasts(rna):
            bases = []
            for path in rna_table_paths(comp):
                if path:
                    bases.extend(list(Path(path).parents)[:3])
            if comp.get("contrast_dir"):
                bases.extend([Path(comp["contrast_dir"]), Path(comp["contrast_dir"]).parent])
            for base in bases:
                for name in ("counts_cpm.tsv", "counts_tpm.tsv", "pseudobulk_cpm.tsv", "pseudobulk_cpm.csv"):
                    if (base / name).is_file():
                        rpath = str(base / name)
                        break
                if rpath:
                    break
            if rpath:
                break
    if not apath:
        for comp in atac_contrasts(atac):
            if comp.get("full_results_csv"):
                candidate = Path(comp["full_results_csv"]).parent / "bulk_atac_peak_cpm.tsv"
                if candidate.is_file():
                    apath = str(candidate)
                    break
    return rpath, apath


def read_cpm(path, aliases, warnings):
    """Feature-by-biological-unit matrices or explicit long sample_id matrices."""
    if not path or not Path(path).is_file():
        return pd.DataFrame()
    sep = "\t" if str(path).endswith((".tsv", ".txt", ".tsv.gz")) else ","
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            header = next(csv.reader(handle, delimiter=sep), [])
        if len(set(header)) != len(header):
            warnings.append("cpm_matrix_skipped: duplicate_sample_headers")
            return pd.DataFrame()
        # Explicit unit names such as 001 must retain their exact string form.
        df = pd.read_csv(path, sep=sep, dtype=str)
    except (OSError, ValueError, pd.errors.ParserError):
        warnings.append("cpm_matrix_unreadable")
        return pd.DataFrame()
    if df.empty:
        return pd.DataFrame()
    # Prefer unique gene_id over non-unique gene_symbol (D4-012), while
    # preserving any inline gene_symbol column as a local alias map so distinct
    # Ensembl loci sharing an HGNC symbol sum cleanly instead of aborting.
    feature = _column(df, "gene_id", "gene_symbol", "gene", "peak", "feature", "Unnamed: 0")
    feature = feature if feature is not None else df.columns[0]
    sym_col = _column(df, "gene_symbol", "symbol")
    table_aliases: dict[str, str] = {}
    if sym_col is not None and sym_col != feature:
        for fid, sym in zip(df[feature], df[sym_col], strict=True):
            if pd.notna(fid) and pd.notna(sym):
                fid_s = str(fid).split(".")[0].strip().upper()
                sym_s = str(sym).strip()
                if fid_s and sym_s and sym_s.lower() != "nan":
                    table_aliases.setdefault(fid_s, sym_s)
    unit = _column(df, "sample_id", "donor_id")
    amount = _column(df, "cpm", "tpm", "expression")
    if unit is not None and amount is not None:
        if df.duplicated([feature, unit]).any():
            warnings.append("cpm_matrix_skipped: duplicate_feature_sample_rows")
            return pd.DataFrame()
        df = df.pivot(index=feature, columns=unit, values=amount)
    else:
        df = df.set_index(feature)
        df = df.drop(columns=[c for c in ("gene_id", "gene_symbol", "symbol") if c in df], errors="ignore")
    df.columns = df.columns.astype(str)
    if df.columns.duplicated().any() or df.index.duplicated().any():
        warnings.append("cpm_matrix_skipped: duplicate_feature_or_sample_ids")
        return pd.DataFrame()
    df = df.apply(pd.to_numeric, errors="coerce")
    eff_aliases = {**table_aliases, **aliases} if table_aliases else aliases
    df.index = [eff_aliases.get(str(x).split(".")[0].upper(), str(x)) for x in df.index]
    if df.index.duplicated().any():
        warnings.append("cpm_feature_alias_collisions: summed_mapped_loci")
        df = df.groupby(level=0, sort=False).sum(min_count=1)
    return df


from aria.utils.stats import bh_correct as bh_adjust


def _norm_chrom(chrom):
    chrom = str(chrom)
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return "MT" if chrom in ("M", "MT") else chrom


def _standardize(matrix):
    centered = matrix - matrix.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    return np.divide(centered, norms, out=np.zeros_like(centered), where=norms > 0)


def _residualize_covariates(matrix, shared_samples, sample_covariates):
    """OLS-residualize rows of `matrix` (n_features x n_samples) against covariates.

    Returns (residual_matrix, adjustment_meta). Fails closed (returns original
    matrix and adjusted=False) if covariates are missing, incomplete, constant,
    or leave fewer than 3 residual degrees of freedom.
    """
    n_samples = len(shared_samples)
    if not isinstance(sample_covariates, dict) or not sample_covariates:
        return matrix, {
            "adjusted": False,
            "covariates_adjusted": [],
            "reason": "no_covariates_supplied",
            "unadjusted_confounding_risk": True,
        }

    cols = [np.ones(n_samples, dtype=float)]
    used_covs: list[str] = []
    skipped_covs: dict[str, str] = {}
    for cov_name in sorted(sample_covariates.keys()):
        raw_map = sample_covariates.get(cov_name)
        if not isinstance(raw_map, dict):
            skipped_covs[str(cov_name)] = "invalid_mapping"
            continue
        vals = [raw_map.get(s) for s in shared_samples]
        if any(v is None or (isinstance(v, float) and not math.isfinite(v)) or str(v).strip() == "" for v in vals):
            skipped_covs[str(cov_name)] = "missing_values_for_shared_samples"
            continue
        # Try numeric first; fall back to categorical dummy encoding (drop first level).
        is_num = True
        num_arr = []
        for v in vals:
            if isinstance(v, bool):
                is_num = False
                break
            try:
                fv = float(v)
                if not math.isfinite(fv):
                    is_num = False
                    break
                num_arr.append(fv)
            except (TypeError, ValueError):
                is_num = False
                break
        if is_num:
            arr = np.asarray(num_arr, dtype=float)
            std = float(np.std(arr))
            if std <= 1e-12:
                skipped_covs[str(cov_name)] = "constant_covariate"
                continue
            cols.append((arr - float(np.mean(arr))) / std)
            used_covs.append(str(cov_name))
        else:
            cats = [str(v).strip() for v in vals]
            levels = sorted(set(cats))
            if len(levels) < 2:
                skipped_covs[str(cov_name)] = "constant_covariate"
                continue
            for lvl in levels[1:]:
                cols.append(np.asarray([1.0 if c == lvl else 0.0 for c in cats], dtype=float))
            used_covs.append(str(cov_name))

    if not used_covs:
        return matrix, {
            "adjusted": False,
            "covariates_adjusted": [],
            "skipped_covariates": skipped_covs,
            "reason": "no_estimable_covariates",
            "unadjusted_confounding_risk": True,
        }

    design = np.column_stack(cols)
    u, s, _ = np.linalg.svd(design, full_matrices=False)
    tol = max(design.shape) * np.finfo(float).eps * (float(s[0]) if s.size else 1.0)
    rank = int(np.sum(s > tol))
    residual_df = n_samples - rank
    if residual_df < 3:
        return matrix, {
            "adjusted": False,
            "covariates_adjusted": [],
            "candidate_covariates": used_covs,
            "design_rank": rank,
            "residual_df": residual_df,
            "reason": "insufficient_residual_df_for_covariates",
            "unadjusted_confounding_risk": True,
        }

    # Orthogonal projector onto Col(design)^perp (rank-safe via thin SVD).
    # matrix: (n_features, n_samples), proj_perp: (n_samples, n_samples).
    u_orth = u[:, :rank]
    proj_perp = np.eye(n_samples, dtype=float) - (u_orth @ u_orth.T)
    resid = matrix @ proj_perp
    return resid, {
        "adjusted": True,
        "covariates_adjusted": used_covs,
        "skipped_covariates": skipped_covs,
        "design_rank": rank,
        "residual_df": residual_df,
        "method": "ols_residualization_on_log2_cpm",
        "permutation_method": "freedman_lane_residual_reprojection",
        "unadjusted_confounding_risk": False,
        "_proj_perp": proj_perp,
    }


# Minimum sample-matched biological units below which across-sample cis
# correlation power is low (df < 8 -> power < 80% even for strong effects after
# BH adjustment across tested peak-gene pairs; D7-014).
LOW_POWER_MIN_MATCHED_SAMPLES = 10


def sample_matched_layer(rna_cpm, atac_cpm, models, params, warnings):
    """Cis correlations with a pooled empirical sample-label permutation null.

    Each null draw applies one common permutation to RNA biological units,
    preserving feature covariance. When ``params['sample_covariates']`` is
    provided with sufficient residual df, both RNA and ATAC log2(CPM+1) matrices
    are OLS-residualized before correlation; otherwise an explicit unadjusted
    covariate warning and sex-chromosome link audit are emitted.
    """
    result = {"ran": False, "reason": "cpm_matrices_unavailable", "n_shared_samples": 0, "shared_sample_ids": [],
              "low_power_warning": None,
              "cis_distance_bp": params["cis_distance_bp"], "min_required": params["min_matched_samples"],
              "n_pairs_tested": 0, "n_significant_links": 0, "n_positive_links": 0, "n_negative_links": 0,
              "n_sex_chromosome_significant_links": 0, "sex_chromosome_significant_fraction": 0.0,
              "covariate_adjustment": {"adjusted": False, "covariates_adjusted": [], "reason": "layer_not_run",
                                       "unadjusted_confounding_risk": True},
              "permutation_null": {}, "top_links": []}
    if rna_cpm.empty or atac_cpm.empty:
        return pd.DataFrame(), result
    mapping = params.get("sample_id_map") or {}
    aligned = [(x, str(mapping.get(x, x))) for x in rna_cpm.columns if str(mapping.get(x, x)) in atac_cpm.columns]
    if len({y for _, y in aligned}) != len(aligned):
        result["reason"] = "non_unique_sample_id_map"
        return pd.DataFrame(), result
    shared = [x for x, _ in aligned]
    result.update(n_shared_samples=len(shared), shared_sample_ids=shared, sample_id_pairs=[{"rna": x, "atac": y} for x, y in aligned])
    if not aligned:
        result["reason"] = "unmatched_sample_namespaces"
        return pd.DataFrame(), result
    if len(aligned) < params["min_matched_samples"]:
        result["reason"] = "insufficient_matched_biological_units"
        return pd.DataFrame(), result
    if not models:
        result["reason"] = "gtf_gene_models_unavailable"
        return pd.DataFrame(), result
    rna = rna_cpm[shared]
    atac = atac_cpm[[y for _, y in aligned]]
    rna = rna[np.isfinite(rna).all(axis=1) & (rna >= 0).all(axis=1)]
    atac = atac[np.isfinite(atac).all(axis=1) & (atac >= 0).all(axis=1)]
    rna, atac = np.log2(rna + 1.0), np.log2(atac + 1.0)
    rna, atac = rna[rna.var(axis=1) > 0], atac[atac.var(axis=1) > 0]
    if len(atac) > 50_000:
        atac = atac.loc[atac.var(axis=1).sort_values(ascending=False, kind="stable").head(50_000).index]
        warnings.append("sample_link_candidate_cap: top_50000_variable_peaks")

    sample_covariates = params.get("sample_covariates")
    rna_arr, cov_adj = _residualize_covariates(rna.to_numpy(dtype=float), shared, sample_covariates)
    proj_perp = cov_adj.pop("_proj_perp", None)
    if cov_adj.get("adjusted") and proj_perp is not None:
        atac_arr = atac.to_numpy(dtype=float) @ proj_perp
        rna_keep = np.var(rna_arr, axis=1) > 1e-12
        atac_keep = np.var(atac_arr, axis=1) > 1e-12
        rna = rna.iloc[rna_keep]
        atac = atac.iloc[atac_keep]
        rna_arr = rna_arr[rna_keep]
        atac_arr = atac_arr[atac_keep]
    else:
        atac_arr = atac.to_numpy(dtype=float)
        warnings.append("cpm_sample_matched_layer_unadjusted_for_covariates")

    result["covariate_adjustment"] = cov_adj
    rlookup = {str(g).upper(): i for i, g in enumerate(rna.index)}
    # Search sorted peak centers for every measured gene TSS, retaining every
    # cis candidate rather than only the nearest gene.
    chrom_peaks = defaultdict(list)
    for i, peak in enumerate(atac.index):
        parsed = parse_peak_interval(str(peak))
        if parsed:
            chrom, start, end = parsed
            chrom_peaks[_norm_chrom(chrom)].append(((start + end) // 2, i))
    rows, ri, ai, seen = [], [], [], set()
    window = params["cis_distance_bp"]
    for chrom, genes in models.items():
        norm_c = _norm_chrom(chrom)
        is_sex_chrom = norm_c in {"X", "Y"}
        peaks = sorted(chrom_peaks.get(norm_c, []))
        centers = np.array([x[0] for x in peaks])
        for gene in genes:
            gene_index = rlookup.get(gene["gene"].upper())
            if gene_index is None:
                continue
            lo = np.searchsorted(centers, gene["tss"] - window, side="left")
            hi = np.searchsorted(centers, gene["tss"] + window, side="right")
            for center, peak_index in peaks[lo:hi]:
                if (gene_index, peak_index) in seen:
                    continue
                seen.add((gene_index, peak_index))
                ri.append(gene_index)
                ai.append(peak_index)
                distance = center - gene["tss"]
                rows.append({
                    "gene": str(rna.index[gene_index]),
                    "peak": str(atac.index[peak_index]),
                    "chrom": norm_c,
                    "is_sex_chromosome": is_sex_chrom,
                    "distance_to_tss": distance if gene["strand"] != "-" else -distance,
                })
    if not rows:
        result["reason"] = "no_variable_cis_pairs"
        return pd.DataFrame(), result
    ri, ai = np.asarray(ri), np.asarray(ai)
    rz, az = _standardize(rna_arr), _standardize(atac_arr)
    rrank, arank = _standardize(rankdata(rna_arr, axis=1)), _standardize(rankdata(atac_arr, axis=1))
    observed = np.clip(np.einsum("ij,ij->i", rz[ri], az[ai]), -1, 1)
    rho = np.clip(np.einsum("ij,ij->i", rrank[ri], arank[ai]), -1, 1)
    rng = np.random.default_rng(params["seed"])
    draws = params["n_permutations"]
    def permutation():
        perm = rng.permutation(len(shared))
        while 2 <= len(shared) <= 8 and np.array_equal(perm, np.arange(len(shared))):
            perm = rng.permutation(len(shared))
        return perm
    # Chunk pairs so temporary sample arrays do not scale with the full family.
    # Under covariate adjustment, re-project permuted RNA residuals back onto
    # Col(design)^perp and re-normalize (Freedman-Lane) so null draws inhabit
    # the exact same (n - rank) subspace as observed partial correlations.
    null = np.empty((draws, len(rows)), dtype=float)
    for draw in range(draws):
        perm = permutation()
        rz_draw = (
            _standardize(rz[:, perm] @ proj_perp)
            if proj_perp is not None
            else rz[:, perm]
        )
        for start in range(0, len(rows), 4096):
            end = start + 4096
            null[draw, start:end] = np.abs(np.einsum("ij,ij->i", rz_draw[ri[start:end]], az[ai[start:end]]))
    pool = np.sort(null.ravel())
    empirical = lambda values: (len(pool) - np.searchsorted(pool, np.maximum(0, np.abs(values) - 1e-12), side="left") + 1) / (len(pool) + 1)
    p_emp = empirical(observed)
    adjusted = bh_adjust(p_emp)
    held_perm = permutation()
    rz_held = (
        _standardize(rz[:, held_perm] @ proj_perp)
        if proj_perp is not None
        else rz[:, held_perm]
    )
    held_r = np.empty(len(rows), dtype=float)
    for start in range(0, len(rows), 4096):
        end = start + 4096
        held_r[start:end] = np.einsum("ij,ij->i", rz_held[ri[start:end]], az[ai[start:end]])
    held_padj = bh_adjust(empirical(held_r))
    df_corr = max(1, (cov_adj["residual_df"] - 1) if cov_adj.get("adjusted") else (len(shared) - 2))
    def analytical(r):
        stat = np.abs(r) * np.sqrt(df_corr / np.maximum(1 - r * r, np.finfo(float).eps))
        return 2 * t.sf(stat, df=df_corr)
    frame = pd.DataFrame(rows)
    frame["pearson_r"], frame["spearman_r"] = observed, rho
    frame["pearson_pvalue"], frame["spearman_pvalue"] = analytical(observed), analytical(rho)
    frame["empirical_pvalue"], frame["padj"] = p_emp, adjusted
    frame["n_samples"] = len(shared)
    frame = frame.sort_values(["padj", "empirical_pvalue", "gene", "peak"], kind="stable")
    significant = frame[frame.padj < params["padj_threshold"]]
    n_sex_sig = int(significant["is_sex_chromosome"].sum()) if len(significant) > 0 else 0
    sex_sig_frac = round(float(n_sex_sig) / float(len(significant)), 4) if len(significant) > 0 else 0.0
    if not cov_adj.get("adjusted") and n_sex_sig > 0:
        warnings.append(f"unadjusted_sex_chromosome_links_detected:{n_sex_sig}/{len(significant)}")
    held_n = int((held_padj < params["padj_threshold"]).sum())
    empirical_fdr_at_padj = (
        min(1.0, float(held_n) / float(len(significant)))
        if len(significant) > 0
        else None
    )
    result.update(ran=True, reason=None, n_pairs_tested=len(frame), n_significant_links=len(significant),
                  low_power_warning=bool(len(shared) < LOW_POWER_MIN_MATCHED_SAMPLES),
                  n_positive_links=int((significant.pearson_r > 0).sum()), n_negative_links=int((significant.pearson_r < 0).sum()),
                  n_sex_chromosome_significant_links=n_sex_sig,
                  sex_chromosome_significant_fraction=sex_sig_frac,
                  covariate_adjustment=cov_adj,
                  permutation_null={"n_permutations": draws, "seed": params["seed"], "method": "pooled_absolute_pearson_sample_permutation",
                                    "n_null_statistics": len(pool), "held_out_n_significant": held_n,
                                    "held_out_discovery_fraction": float(np.mean(held_padj < params["padj_threshold"])),
                                    "empirical_fdr_at_padj": empirical_fdr_at_padj,
                                    "held_out_draw_independent": True}, top_links=frame.head(20).to_dict("records"))
    return frame, result


def _adaptive_neglog10(pvals):
    """Return (-log10(p) with data-adaptive cap for p==0, cap_value, has_zero_cap)."""
    arr = np.asarray(pvals, dtype=float)
    valid = np.isfinite(arr) & (arr >= 0.0) & (arr <= 1.0)
    pos = arr[valid & (arr > 0.0)]
    if pos.size > 0:
        max_fin = float(-np.log10(np.min(pos)))
        cap = max(max_fin * 1.15, max_fin + 1.0, 5.0)
    else:
        cap = 10.0
    zero_mask = valid & (arr == 0.0)
    out = np.full(arr.shape, np.nan, dtype=float)
    out[valid & (arr > 0.0)] = -np.log10(arr[valid & (arr > 0.0)])
    out[zero_mask] = cap
    return out, cap, bool(np.any(zero_mask))


def render_starburst(de_da, tfs, pathways, summary, output_dir):
    """Render editable, headless publication panels with associative labels."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib import rc_context

    png, svg = Path(output_dir) / "cross_modal_starburst.png", Path(output_dir) / "cross_modal_starburst.svg"
    with rc_context({"svg.fonttype": "none", "svg.hashsalt": "aria_cross_modal_concordance",
                     "font.size": 9, "axes.spines.top": False, "axes.spines.right": False}):
        fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), layout="constrained")
        a, b, c = axes
        a.set_title("A  Cross-modal effect concordance", loc="left", fontweight="bold")
        a.set(xlabel="RNA log2 fold change", ylabel="ATAC cis peak log2 fold change")
        a.axhline(0, color="#aaaaaa", lw=.7)
        a.axvline(0, color="#aaaaaa", lw=.7)
        colors = {"concordant_up": "#bc4934", "concordant_down": "#346da2", "discordant": "#9467a6", "single_modality": "#b9b9b9"}
        if not de_da.empty:
            plot = de_da.dropna(subset=["rna_log2fc", "atac_log2fc"]).copy()
            r_lfc = pd.to_numeric(plot["rna_log2fc"], errors="coerce").to_numpy(dtype=float)
            a_lfc = pd.to_numeric(plot["atac_log2fc"], errors="coerce").to_numpy(dtype=float)
            plot = plot[np.isfinite(r_lfc) & np.isfinite(a_lfc)].copy()
            if not plot.empty:
                plot["plot_group"] = plot.classification.map(lambda x: x if x in ("concordant_up", "concordant_down") else "discordant" if str(x).startswith("discordant") else "single_modality")
                for key, color in colors.items():
                    part = plot[plot.plot_group == key]
                    if not part.empty:
                        a.scatter(part.rna_log2fc, part.atac_log2fc, s=21, alpha=.65, color=color, label=key.replace("_", " "), edgecolors="none")
                labels = plot[plot.classification.astype(str).str.startswith(("concordant", "discordant"))].sort_values(["rna_padj", "atac_padj", "gene_key"], kind="stable").drop_duplicates("gene_key").head(8)
                for i, row in enumerate(labels.itertuples()):
                    a.annotate(safe_label(row.gene), (row.rna_log2fc, row.atac_log2fc), xytext=(4, 5 + (i % 2) * 7), textcoords="offset points", fontsize=7)
                a.legend(frameon=False, fontsize=7, loc="lower right")
            else:
                a.text(.5, .5, "No finite RNA/ATAC log2FC pairs available", ha="center", transform=a.transAxes)
        else:
            a.text(.5, .5, "No DE × DA layer available", ha="center", transform=a.transAxes)
        rho = summary.get("spearman_r")
        fp = summary.get("fisher_pvalue")
        fraction = summary.get("concordance_fraction")
        box = f"Spearman ρ = {rho:.2f}" if (rho is not None and math.isfinite(float(rho))) else "Spearman ρ unavailable"
        if fp is not None and math.isfinite(float(fp)):
            fp_val = float(fp)
            box += "\nFisher p < 1e-300" if fp_val == 0.0 else f"\nFisher p = {fp_val:.2g}"
        if fraction is not None and math.isfinite(float(fraction)):
            box += f"\nConcordant fraction = {float(fraction):.2f}"
        a.text(.03, .97, box, va="top", transform=a.transAxes, fontsize=8, bbox={"facecolor": "white", "alpha": .85, "edgecolor": "#dddddd"})
        b.set_title("B  TF motif / footprint × TF RNA", loc="left", fontweight="bold")
        if tfs.empty:
            b.text(.5, .5, "No TF motif/footprint layer available", ha="center", va="center", transform=b.transAxes)
            b.set_axis_off()
        else:
            tf_work = tfs.copy()
            tf_work["_tf_key"] = tf_work["tf"].astype(str).str.upper()
            tf_work["_abs_score"] = pd.to_numeric(tf_work["chromatin_score"], errors="coerce").abs().fillna(-1.0)
            rows = (
                tf_work.sort_values(
                    ["chromatin_padj", "_abs_score", "_tf_key"],
                    ascending=[True, False, True],
                    kind="stable",
                )
                .drop_duplicates("_tf_key")
                .head(10)
                .iloc[::-1]
                .copy()
            )
            c_scores = pd.to_numeric(rows["chromatin_score"], errors="coerce").to_numpy(dtype=float)
            r_scores = pd.to_numeric(rows["rna_log2fc"], errors="coerce").to_numpy(dtype=float)
            # Keep chromatin scores visually comparable to RNA log2FC without 300.0 spikes.
            fin_c = c_scores[np.isfinite(c_scores)]
            fin_r = r_scores[np.isfinite(r_scores)]
            ref_scale = max(
                float(np.percentile(np.abs(fin_r), 95)) if fin_r.size else 2.0,
                float(np.median(np.abs(fin_c))) if fin_c.size else 2.0,
                2.0,
            )
            score_cap = max(ref_scale * 4.0, 10.0)
            c_capped = False
            if np.any(np.isfinite(c_scores) & (np.abs(c_scores) > score_cap)):
                c_scores = np.clip(c_scores, -score_cap, score_cap)
                c_capped = True
            # Fill missing chromatin_score (e.g., motif padj==0 without log2_enrichment) with signed cap.
            missing_c = ~np.isfinite(c_scores)
            if np.any(missing_c):
                fill_mag = min(score_cap, max(float(np.max(np.abs(fin_c))) * 1.15 if fin_c.size else 3.0, 2.0))
                dirs = np.where(rows["chromatin_direction"].astype(str).to_numpy() == "down", -1.0, 1.0)
                c_scores = np.where(missing_c, dirs * fill_mag, c_scores)
                c_capped = True

            y = np.arange(len(rows))
            b.hlines(y, 0, c_scores, color="#759c94", lw=2)
            b.scatter(c_scores, y, color="#326e62", label="Chromatin score", s=36)
            fin_r_mask = np.isfinite(r_scores)
            if np.any(fin_r_mask):
                b.scatter(r_scores[fin_r_mask], y[fin_r_mask], color="#cf7546", label="RNA log2 fold change", marker="D", s=26)
            b.set_yticks(y, [safe_label(x.tf) + (" *" if x.paralog_motif_ambiguity else "") for x in rows.itertuples()])
            b.axvline(0, color="#aaaaaa", lw=.7)
            xlabel_b = "Signed chromatin score / RNA log2 fold change"
            if c_capped:
                xlabel_b += f" (capped ±{score_cap:.0f})"
            if bool(rows["paralog_motif_ambiguity"].any()):
                xlabel_b += "\n* Motif-family or dimer ambiguity"
            b.set_xlabel(xlabel_b)
            b.legend(frameon=False, fontsize=7)
        c.set_title("C  Shared functional pathways", loc="left", fontweight="bold")
        if not pathways.empty:
            ordered = sorted(pathways.to_dict("records"), key=_pathway_sort_key)
            rows = (
                pd.DataFrame(ordered)
                .drop_duplicates("safe_label")
                .head(7)
                .iloc[::-1]
            )
            y = np.arange(len(rows))
            both_p = np.concatenate([
                pd.to_numeric(rows.rna_padj, errors="coerce").to_numpy(dtype=float),
                pd.to_numeric(rows.atac_padj, errors="coerce").to_numpy(dtype=float),
            ])
            both_neg, cap_c, capped_c = _adaptive_neglog10(both_p)
            n_r = len(rows)
            rna_y = np.nan_to_num(both_neg[:n_r], nan=0.0)
            atac_y = np.nan_to_num(both_neg[n_r:], nan=0.0)
            c.barh(y - .18, rna_y, .34, color="#cf7546", label="RNA ORA")
            c.barh(y + .18, atac_y, .34, color="#326e62", label="ATAC peak ORA")
            import textwrap
            c.set_yticks(y, [textwrap.fill(safe_label(x), 27) for x in rows.safe_label])
            c_xlab = "−log10 adjusted p-value"
            if capped_c:
                c_xlab += f" (padj = 0 capped at {cap_c:.1f})"
            c.set_xlabel(c_xlab)
            c.legend(frameon=False, fontsize=7)
        else:
            keys = ("concordant_up", "concordant_down", "discordant_rna_up_atac_down", "discordant_rna_down_atac_up", "atac_da_rna_unchanged", "rna_de_atac_unchanged")
            vals = [int(summary.get(f"n_{k}", 0) or 0) for k in keys]
            if any(v > 0 for v in vals):
                c.barh(np.arange(len(keys)), vals, color=[colors["concordant_up"], colors["concordant_down"], colors["discordant"], colors["discordant"], "#b9b9b9", "#b9b9b9"])
                c.set_yticks(np.arange(len(keys)), [k.replace("_", " ") for k in keys], fontsize=7)
                c.set_xlabel("Co-tested gene–peak pairs")
                c.set_title("C  Cross-modal quadrant counts", loc="left", fontweight="bold")
            else:
                c.text(.5, .5, "No shared pathways or quadrant counts available", ha="center", va="center", transform=c.transAxes)
                c.set_axis_off()
        fig.suptitle("Cross-modal concordance • associative evidence", fontsize=13, fontweight="bold")
        try:
            fig.savefig(png, dpi=300, facecolor="white")
            fig.savefig(svg, format="svg", facecolor="white", metadata={"Date": None})
        finally:
            plt.close(fig)
    return str(png), str(svg)
