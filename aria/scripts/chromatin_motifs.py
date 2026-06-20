"""
ARIA Chromatin Motif Enrichment Script (v4.6 scATAC, step 4)
------------------------------------------------------------
Transcription-factor motif over-representation in differentially accessible
(DA) peak sets, using snapatac2's native enrichment against a versioned,
LOCAL, offline MEME motif collection.

Pipeline:
  1. Read DA peak groups (from step 3's per-cluster CSV, or explicit `regions`);
  2. Read a versioned MEME motif collection (default JASPAR2024 CORE
     vertebrates) resolved from `ARIA_MOTIF_DIR` (W-PRIV: offline, versioned);
  3. Scan each DA peak group for motif hits against the peak universe as
     background and test enrichment (binomial / hypergeometric) with
     `snapatac2.tl.motif_enrichment`;
  4. Report the enriched motifs per group with the exact motif release used.

Executed inside aria-chromatin-env by EnvironmentManager.

Honesty contract (ADR-002 / ADR-011 / W-PRIV):
  - Motif enrichment runs ONLY with a local versioned motif collection AND a
    local genome FASTA. Either missing => honest skip with a concrete reason
    (e.g. "run scripts/fetch_motifs.py"), never a fabricated enrichment.
  - This script performs NO network egress; the only governed download lives in
    `scripts/fetch_motifs.py`.
  - The motif release (collection, release, sha256) is recorded for
    reproducibility. A motif id is a database identifier, not a biological
    claim about the cluster — the report stays descriptive/associative.
  - Per-cell chromVAR-style motif ACTIVITY is out of scope here (a heavier later
    increment); this is over-representation in DA peak sets only, and that limit
    is reported, not hidden.

Input params:
    data_path:        str   — clustered scATAC .h5ad (peak universe = var_names,
                              used as the enrichment background).
    da_csv:           str   (optional) — step-3 `chromatin_da_per_cluster.csv`
                              (columns include `cluster`, `peak`); grouped into
                              per-cluster DA peak sets.
    regions:          dict  (optional) — explicit {group: [peak, ...]} override.
    background:       list  (optional) — explicit background peaks (else all
                              var_names from data_path).
    genome_fasta:     str   — local genome FASTA path (REQUIRED for a real run).
    motif_collection: str   (optional) — collection name under ARIA_MOTIF_DIR
                              (default JASPAR2024_CORE_vertebrates).
    motif_file:       str   (optional) — explicit MEME path (overrides the
                              collection resolution).
    method:           str   (optional) — "hypergeometric" (default; DA peaks are
                              a subset of background) or "binomial".
    padj_max:         float (optional) — enrichment significance cutoff.
    top_n:            int   (optional) — top enriched motifs per group (def 25).
    max_regions_per_group: int (optional) — cap on peaks scanned per group
                              (default 5000; excess is dropped with a warning).
    exp_context:      dict  (optional) — confirmed thresholds (global_padj).
    output_dir:       str   (optional).

Output:
    {
      "status":          "success" | "error",
      "ran":             bool,
      "reason":          str | None,         — why it skipped, if it did
      "method":          str,
      "genome_fasta":    str | None,
      "motif_source":    {collection, release, sha256, n_motifs, source} | None,
      "n_groups":        int,
      "per_group": {
        group: {n_regions, n_regions_dropped, n_motifs_tested,
                n_enriched, top_motifs: [{motif_id, name, family,
                                          log2_enrichment, pvalue, padj}, ...]}
      },
      "output_csv":      str | None,
      "warnings":        [str],
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script, mocks_allowed


def _skip(reason, **extra):
    out = {"status": "success", "ran": False, "reason": reason,
           "per_group": {}, "n_groups": 0, "warnings": []}
    out.update(extra)
    return out


def _regions_from_csv(da_csv, max_per_group, warnings):
    """Build ranked {cluster: [peak, ...]} from a DA CSV.

    When adjusted p-values/effect sizes are available, oversized groups are
    capped after ranking by padj ascending and then absolute LFC descending.
    """
    import csv
    groups: dict = {}
    with open(da_csv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or "peak" not in reader.fieldnames:
            return {}
        cols = reader.fieldnames
        gcol = "cluster" if "cluster" in cols else None
        lfc_col = _detect_col(cols, "log2fold", "log2fc", "logfold")
        padj_col = _detect_col(cols, "padj", "pvals_adj", "adjusted", "fdr",
                               "q-value", "qval")
        sig_col = "significant" if "significant" in cols else None
        for row in reader:
            g = str(row.get(gcol, "all")) if gcol else "all"
            peak = str(row.get("peak", "")).strip()
            if not peak:
                continue
            if sig_col and str(row.get(sig_col, "")).strip().lower() != "true":
                continue
            padj = _as_float(row.get(padj_col)) if padj_col else None
            lfc = _as_float(row.get(lfc_col)) if lfc_col else None
            groups.setdefault(g, []).append(
                (peak, padj if padj is not None else float("inf"),
                 abs(lfc) if lfc is not None else 0.0))
    return {
        g: _rank_da_region_group(rows, max_per_group, g, warnings)
        for g, rows in groups.items()
    }


def _rank_da_region_group(rows, max_per_group, group, warnings):
    ranked = sorted(rows, key=lambda x: (x[1], -x[2], x[0]))
    if len(ranked) > max_per_group:
        warnings.append(
            f"group '{group}': capped {len(ranked)} -> {max_per_group} peaks "
            "for motif scanning after ranking by padj then |log2FC|.")
    return [peak for peak, _padj, _abs_lfc in ranked[:max_per_group]]


def _detect_col(cols, *needles):
    low = {c.lower(): c for c in cols}
    for n in needles:
        for lc, orig in low.items():
            if n in lc:
                return orig
    return None


def _genome_chroms(genome_fasta):
    """Set of contig names in the reference FASTA, or None if undiscoverable.

    Bulk ATAC peak universes (MACS3 on the raw BAM) can include scaffold/alt
    contigs absent from a primary-assembly FASTA; snapatac2 sequence fetching
    raises a KeyError on those. Returning the valid contig set lets the caller
    drop such peaks (disclosed) instead of crashing. A snapatac2 Genome object
    (auto-fetch path) or any non-path input yields None (no filtering). Prefers
    the cheap `.fai` index; falls back to pysam, which builds it.
    """
    from pathlib import Path
    if not isinstance(genome_fasta, str):
        return None
    fai = Path(genome_fasta + ".fai")
    if fai.is_file():
        chroms = set()
        try:
            with open(fai, encoding="utf-8") as fh:
                for line in fh:
                    name = line.split("\t", 1)[0].strip()
                    if name:
                        chroms.add(name)
        except OSError:
            return None
        return chroms or None
    try:
        import pysam
        with pysam.FastaFile(genome_fasta) as fa:
            return set(fa.references) or None
    except Exception:
        return None


def _resolve_modality_thresholds(params: dict):
    """Resolve the analysis modality + padj cutoff for motif enrichment.

    The modality is taken from the caller (``modality`` or ``data_type``), NOT
    hardcoded: bulk ATAC reuses this engine via ``_run_bulk_atac_motifs`` and must
    not inherit the scATAC threshold-policy label. scATAC stays the default for the
    scATAC lane (which does not pass an explicit modality), so existing behaviour is
    unchanged. An explicit ``padj_max`` always wins over the resolved CP3/default
    cutoff. Returns ``(modality, padj_max)``.
    """
    from aria.utils.thresholds import AnalysisThresholds
    modality = str(params.get("modality") or params.get("data_type") or "scATAC")
    thr = AnalysisThresholds.from_exp_context(
        params.get("exp_context"), modality=modality)
    padj_max = float(params.get("padj_max", thr.padj))
    return modality, padj_max


def chromatin_motifs(params: dict) -> dict:
    import json
    from pathlib import Path

    data_path = params.get("data_path")
    da_csv = params.get("da_csv")
    regions_param = params.get("regions")
    genome_fasta = params.get("genome_fasta")
    assembly = str(params.get("assembly") or params.get("genome") or "").strip()
    allow_genome_fetch = bool(params.get("allow_genome_fetch", False))
    motif_collection = params.get("motif_collection")
    motif_file = params.get("motif_file")
    method = params.get("method", "hypergeometric")
    top_n = int(params.get("top_n", 25))
    max_per_group = int(params.get("max_regions_per_group", 5000))
    output_dir = params.get("output_dir")
    allow_mock = mocks_allowed(params)
    truncation_strategy = str(
        params.get("foreground_truncation_strategy")
        or "rank_by_padj_then_abs_log2fc_before_cap_when_scores_available"
    )

    modality, padj_max = _resolve_modality_thresholds(params)

    warnings: list = []

    # ── Resolve the local versioned motif collection (offline, W-PRIV) ───────
    from aria.utils.motifs import (
        load_local_motif_collection, DEFAULT_COLLECTION,
    )
    if motif_file:
        if not Path(motif_file).is_file():
            return _skip(f"motif_file not found: {motif_file}")
        meme_path = str(motif_file)
        motif_version = {"collection": Path(motif_file).stem, "source": "explicit",
                         "release": "unknown"}
    else:
        coll = motif_collection or DEFAULT_COLLECTION
        loaded = load_local_motif_collection(coll)
        if loaded is None:
            return _skip(
                f"motif collection '{coll}' not staged under ARIA_MOTIF_DIR; "
                f"run scripts/fetch_motifs.py to download it (one-time, "
                f"governed). Motif enrichment is offline-only and never "
                f"fabricated.")
        meme_path, motif_version = loaded

    # ── Reference genome (resolved automatically; no env var required) ───────
    # Resolution policy (aria.utils.genomes): explicit path → ARIA_GENOME_FASTA
    # → managed store ARIA_GENOME_DIR/<assembly> → (governed, opt-in) snapatac2
    # auto-fetch by assembly. A missing genome is an honest skip — never an
    # env-var instruction to the user.
    from aria.utils import genomes
    genome_source = "explicit_param" if genome_fasta else None
    if genome_fasta and not Path(str(genome_fasta)).is_file():
        return _skip(f"genome_fasta not found: {genome_fasta}",
                     motif_source=motif_version)
    if not genome_fasta:
        local, src = genomes.resolve_local_genome_fasta(assembly)
        if local:
            genome_fasta, genome_source = local, src

    # Decide whether a governed snapatac2 auto-fetch is warranted. It downloads +
    # caches the assembly FASTA (~hundreds of MB), so it requires an explicit
    # opt-in AND network egress (ARIA_AIR_GAPPED off) — ARIA never fetches
    # silently. The decision is light; the actual fetch happens after the deps.
    attr = genomes.snapatac2_attr(assembly)
    will_fetch = bool((not genome_fasta) and attr and allow_genome_fetch)
    if will_fetch:
        from aria.utils import privacy
        if not privacy.egress_allowed():
            will_fetch = False

    # No genome and not going to fetch → honest skip NOW (no heavy deps needed),
    # with a user-facing reason (never an env-var instruction).
    if not genome_fasta and not will_fetch:
        if attr and not allow_genome_fetch:
            reason = (f"reference genome for assembly '{assembly}' is not staged "
                      f"locally; ARIA can auto-download it (~hundreds of MB) once "
                      f"the genome checkpoint is approved. Skipping motif "
                      f"enrichment for now.")
        elif attr:
            reason = (f"reference genome for assembly '{assembly}' is unavailable "
                      f"(air-gapped, no local copy). Stage one under "
                      f"ARIA_GENOME_DIR/<assembly> to enable motif enrichment.")
        else:
            reason = (f"could not resolve a reference genome for assembly "
                      f"'{assembly or 'unknown'}'. Provide one via the run or "
                      f"stage it under ARIA_GENOME_DIR/<assembly>.")
        return _skip(reason, motif_source=motif_version)

    # ── Dependencies ─────────────────────────────────────────────────────────
    try:
        import snapatac2 as snap
        from aria.utils.safe_h5ad import read_h5ad
    except ImportError as e:
        if not allow_mock:
            return {"status": "error", "error_type": "MissingDependency",
                    "details": (f"motif enrichment requires snapatac2 "
                                f"(aria-chromatin-env). ImportError: {e}")}
        return _skip(f"mock (deps missing: {e})", motif_source=motif_version)

    # Governed auto-fetch (after deps): snap.genome.<attr> downloads + caches.
    if not genome_fasta and will_fetch:
        try:
            gobj = getattr(snap.genome, attr, None)
            if gobj is not None:
                genome_fasta = gobj       # snapatac2 Genome object (cached fetch)
                genome_source = f"snapatac2_auto_fetch:{attr}"
        except Exception as e:
            warnings.append(f"snapatac2 genome auto-fetch failed: {e}")
        if not genome_fasta:
            return _skip(f"snapatac2 auto-fetch did not yield a genome for "
                         f"assembly '{assembly}'.", motif_source=motif_version)

    # ── Assemble region groups + background ──────────────────────────────────
    if regions_param:
        groups = {str(k): [str(p) for p in v] for k, v in regions_param.items()}
    elif da_csv and Path(da_csv).is_file():
        groups = _regions_from_csv(da_csv, max_per_group, warnings)
    else:
        return _skip("no DA regions supplied (need `regions` or `da_csv` from "
                     "chromatin_diffacc step 3).", motif_source=motif_version)
    groups = {g: peaks for g, peaks in groups.items() if peaks}
    if not groups:
        return _skip("DA region groups are empty.", motif_source=motif_version)

    background = params.get("background")
    if not background:
        if not data_path or not Path(data_path).exists():
            return _skip("no background peaks: supply `background` or a "
                         "`data_path` clustered .h5ad.",
                         motif_source=motif_version)
        adata = read_h5ad(str(data_path))
        background = [str(p) for p in adata.var_names]
    background = [str(p) for p in background]

    # Bulk ATAC peaks (MACS3 on the raw BAM) may sit on scaffold/alt contigs
    # absent from the reference FASTA; snapatac2 sequence fetching would crash
    # on them. Drop out-of-reference peaks up front and disclose the count.
    # Regions are intersected with this background below, so off-contig DA
    # peaks are removed (and disclosed) transitively.
    valid_chroms = _genome_chroms(genome_fasta)
    if valid_chroms is not None:
        before_bg = len(background)
        background = [p for p in background
                      if p.split(":", 1)[0] in valid_chroms]
        dropped_bg = before_bg - len(background)
        if dropped_bg:
            warnings.append(
                f"dropped {dropped_bg} background peak(s) on contigs absent "
                f"from the reference genome before motif enrichment.")
        if not background:
            return _skip("no background peaks remain after restricting to "
                         "reference-genome contigs.",
                         motif_source=motif_version)
    bg_set = set(background)

    # For hypergeometric, regions MUST be a subset of background; intersect and
    # cap per group. Dropped peaks are disclosed, never silently ignored.
    clean_groups: dict = {}
    dropped_counts: dict = {}
    for g, peaks in groups.items():
        in_bg = [p for p in peaks if p in bg_set]
        dropped = len(peaks) - len(in_bg)
        if len(in_bg) > max_per_group:
            warnings.append(
                f"group '{g}': capped {len(in_bg)} -> {max_per_group} peaks "
                f"for motif scanning using {truncation_strategy}.")
            in_bg = in_bg[:max_per_group]
        if dropped:
            warnings.append(
                f"group '{g}': dropped {dropped} peak(s) not in the background "
                f"peak universe before enrichment.")
        if in_bg:
            clean_groups[g] = in_bg
            dropped_counts[g] = dropped
    if not clean_groups:
        return _skip("no DA peaks overlap the background peak universe.",
                     motif_source=motif_version)

    # ── Read motifs + run enrichment ─────────────────────────────────────────
    motifs = snap.datasets.read_motifs(meme_path)
    n_motifs_loaded = len(motifs)
    try:
        results = snap.tl.motif_enrichment(
            motifs=motifs,
            regions=clean_groups,
            genome_fasta=genome_fasta,
            background=background,
            method=method,
        )
    except Exception as e:
        return {"status": "error", "error_type": "MotifEnrichmentFailed",
                "details": f"{type(e).__name__}: {str(e)[:400]}",
                "motif_source": motif_version}

    # ── Normalize polars results -> plain dicts ──────────────────────────────
    per_group: dict = {}
    rows_csv: list = []
    for g, df in results.items():
        try:
            recs = df.to_dicts()
        except Exception:
            recs = []
        cols = list(recs[0].keys()) if recs else []
        id_c = _detect_col(cols, "id")
        name_c = _detect_col(cols, "name")
        fam_c = _detect_col(cols, "family")
        p_c = _detect_col(cols, "p-value", "pvalue", "p_value", "pval")
        padj_c = _detect_col(cols, "adjusted", "padj", "fdr", "q-value", "qval")
        lfc_c = _detect_col(cols, "log2", "fold", "enrichment", "odds")

        norm = []
        for r in recs:
            padj_v = _as_float(r.get(padj_c)) if padj_c else None
            mid = str(r.get(id_c)) if id_c and r.get(id_c) is not None else None
            nm = (str(r.get(name_c))
                  if name_c and r.get(name_c) is not None else None)
            # JASPAR MEME records carry id + TF name in one token
            # (e.g. "MA0752.2 ZNF410"); split so the matrix id and the TF name
            # are reported separately. The TF name is a database fact, not an
            # inferred biological claim about the cluster.
            if nm is None and mid and " " in mid:
                matrix_id, nm = mid.split(None, 1)
            else:
                matrix_id = mid
            norm.append({
                "motif_id": matrix_id,
                "name": nm,
                "family": str(r.get(fam_c)) if fam_c else None,
                "log2_enrichment": _as_float(r.get(lfc_c)) if lfc_c else None,
                "pvalue": _as_float(r.get(p_c)) if p_c else None,
                "padj": padj_v,
            })
        # sort by padj (then pvalue) ascending; count enriched at threshold
        norm.sort(key=lambda x: (x["padj"] if x["padj"] is not None else 1.0,
                                 x["pvalue"] if x["pvalue"] is not None else 1.0))
        n_enriched = sum(1 for x in norm
                         if x["padj"] is not None and x["padj"] < padj_max)
        top = norm[:top_n]
        for x in top:
            row = {"group": g}
            row.update(x)
            rows_csv.append(row)
        per_group[g] = {
            "n_regions": len(clean_groups[g]),
            "n_regions_dropped": dropped_counts.get(g, 0),
            "n_motifs_tested": n_motifs_loaded,
            "n_enriched": n_enriched,
            "top_motifs": top,
        }

    out_dir = Path(output_dir) if output_dir else (
        Path(data_path).parent if data_path else Path("."))
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = None
    if rows_csv:
        import csv as _csv
        csv_path = str(out_dir / "chromatin_motif_enrichment.csv")
        try:
            with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(rows_csv[0].keys()))
                w.writeheader()
                w.writerows(rows_csv)
        except Exception:
            csv_path = None

    return {
        "status": "success",
        "ran": True,
        "reason": None,
        "modality": modality,
        "padj_max": padj_max,
        "method": method,
        "genome_fasta": str(genome_fasta),
        "genome_source": genome_source,
        "assembly": assembly or None,
        "motif_source": motif_version,
        "n_groups": len(per_group),
        "per_group": per_group,
        "output_csv": csv_path,
        "truncation_strategy": truncation_strategy,
        "warnings": warnings,
    }


def _as_float(v):
    try:
        f = float(v)
        return f
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    run_script(chromatin_motifs)
