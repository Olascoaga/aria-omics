"""Local, offline over-representation analysis (ORA) for ARIA (P1-7, W-PRIV).

Default pathway enrichment runs a **hypergeometric test locally** against
versioned GMT gene-set libraries and the dataset's expressed-gene background, so
the DE gene list never leaves the machine. Enrichr (which ships the gene list to
an external server) becomes **opt-in** via ``ARIA_ALLOW_ENRICHR=1``.

Gene-set libraries are read from ``ARIA_GMT_DIR`` (default ``~/.aria/genesets``),
one sub-directory per library:

    <genesets_dir>/<library>/<library>.gmt
    <genesets_dir>/<library>/manifest.json   # {library, source, release, date, sha256}

The manifest is surfaced as ``gene_set_version`` so ``methodology.json`` records
exactly which release of each library produced the enrichment (reproducibility).

This module is dependency-light on purpose (stdlib ``math`` + ``aria.utils.stats``
for BH) so the enrichment engine is testable in the light CI lane without gseapy
or scipy. ``scripts/fetch_genesets.py`` is the explicit, one-time online
bootstrap that materializes the versioned ``.gmt`` files from Enrichr libraries.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

from aria.utils.reference_integrity import (
    public_integrity_result, reference_is_usable, verify_reference_file,
)
from aria.utils.stats import bh_correct

ENRICHR_OPT_IN_ENV = "ARIA_ALLOW_ENRICHR"
GMT_DIR_ENV = "ARIA_GMT_DIR"


def enrichr_opt_in() -> bool:
    """True only when the user explicitly opts into Enrichr network egress.

    Local hypergeometric ORA is the default; Enrichr is never used unless this
    flag is set (and, separately, air-gapped mode is off — see
    ``aria.utils.privacy.egress_allowed``)."""
    return os.environ.get(ENRICHR_OPT_IN_ENV, "0").strip().lower() in {
        "1", "true", "yes", "on",
    }


def genesets_dir() -> Path:
    """Resolve the versioned GMT library directory.

    ``ARIA_GMT_DIR`` overrides; otherwise ``$ARIA_HOME/genesets`` and finally
    ``~/.aria/genesets``."""
    override = os.environ.get(GMT_DIR_ENV)
    if override:
        return Path(override).expanduser()
    home = os.environ.get("ARIA_HOME")
    base = Path(home).expanduser() if home else (Path.home() / ".aria")
    return base / "genesets"


def parse_gmt(path) -> dict:
    """Parse a GMT file into ``{term: [gene, ...]}``.

    GMT format: ``term<TAB>description<TAB>gene1<TAB>gene2<TAB>...``. Gene
    symbols are upper-cased so matching is case-insensitive (Enrichr stores
    human symbols upper-cased; mouse symbols are mixed-case)."""
    sets: dict = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            genes = [g.strip().upper() for g in parts[2:] if g and g.strip()]
            if term and genes:
                sets[term] = genes
    return sets


def load_local_library(library_name: str):
    """Return ``(gene_sets, version)`` for a local versioned GMT, or ``None``.

    ``version`` is the parsed ``manifest.json`` (plus a fallback ``n_terms``),
    so the report can state the exact gene-set release used."""
    base = genesets_dir() / library_name
    gmt = base / f"{library_name}.gmt"
    if not gmt.is_file():
        return None
    gene_sets = parse_gmt(gmt)
    if not gene_sets:
        return None
    version: dict = {"library": library_name, "n_terms": len(gene_sets)}
    manifest = base / "manifest.json"
    integrity = verify_reference_file(gmt, manifest_path=manifest)
    if not reference_is_usable(integrity):
        return None
    if manifest.is_file():
        manifest_data = integrity.get("manifest")
        if isinstance(manifest_data, dict):
            version.update(manifest_data)
        else:
            try:
                version.update(json.loads(manifest.read_text(encoding="utf-8")))
            except Exception:
                pass
    version.setdefault("source", "unknown")
    version.setdefault("release", "unknown")
    version["n_terms"] = len(gene_sets)
    version["integrity"] = public_integrity_result(integrity)
    return gene_sets, version


def _logcomb(n: int, k: int) -> float:
    """log(n choose k); -inf when out of range so exp() -> 0."""
    if k < 0 or k > n:
        return float("-inf")
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    )


def hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """Upper-tail hypergeometric survival ``P(X >= k)``.

    ``M`` population size (background), ``n`` successes in population (genes in
    the term ∩ background), ``N`` draws (query genes ∩ background), ``k``
    observed overlap. This is the standard over-representation p-value."""
    if k <= 0:
        return 1.0
    hi = min(n, N)
    if k > hi:
        return 0.0
    log_denom = _logcomb(M, N)
    if log_denom == float("-inf"):
        return 1.0
    total = 0.0
    for i in range(k, hi + 1):
        total += math.exp(_logcomb(n, i) + _logcomb(M - n, N - i) - log_denom)
    return min(1.0, max(0.0, total))


def _odds_ratio(k: int, M: int, n: int, N: int) -> float:
    """2x2 odds ratio for the overlap, Haldane-corrected to avoid div-by-zero."""
    a = k                       # in query & in term
    b = N - k                   # in query, not in term
    c = n - k                   # not in query, in term
    d = M - n - b               # not in query, not in term
    a += 0.5
    b += 0.5
    c += 0.5
    d += 0.5
    return (a * d) / (b * c)


def run_ora(query_genes, gene_sets: dict, background_genes,
            *, padj_max: float = 0.05, min_overlap: int = 1,
            top: int = 20) -> list:
    """Local hypergeometric ORA of ``query_genes`` against ``gene_sets``.

    The universe is ``background_genes`` (the dataset's expressed genes); gene
    sets and the query are restricted to it. When no background is given the
    universe falls back to the union of all gene-set genes ∪ query (standard
    ORA fallback). BH-corrects across all tested terms, returns the significant
    terms (``padj < padj_max``) sorted by padj, capped at ``top``.

    Output dicts match the Enrichr-path schema so downstream viz/narrators are
    unaffected: ``{term, padj, pvalue, overlap, odds_ratio, combined_score,
    genes}``."""
    query = {str(g).upper() for g in (query_genes or []) if g and str(g).lower() != "nan"}
    bg = {str(g).upper() for g in (background_genes or []) if g and str(g).lower() != "nan"}

    if bg:
        query = query & bg
    else:
        bg = set(query)
        for genes in gene_sets.values():
            bg.update(g.upper() for g in genes)

    M = len(bg)
    N = len(query)
    if N == 0 or M == 0:
        return []

    tested = []
    for term, genes in gene_sets.items():
        gs = {g.upper() for g in genes}
        if background_genes:
            gs = gs & bg
        n = len(gs)
        if n == 0:
            continue
        overlap = query & gs
        k = len(overlap)
        if k < min_overlap:
            continue
        p = hypergeom_sf(k, M, n, N)
        odds = _odds_ratio(k, M, n, N)
        tested.append({
            "term": term,
            "pvalue": p,
            "_k": k,
            "_n": n,
            "_genes": sorted(overlap),
            "odds_ratio": round(odds, 2),
        })

    if not tested:
        return []

    padj = bh_correct([t["pvalue"] for t in tested])
    out = []
    for t, q in zip(tested, padj):
        p = t["pvalue"]
        combined = round(-math.log10(max(p, 1e-300)) * max(t["odds_ratio"], 0.0), 1)
        out.append({
            "term": t["term"],
            "padj": round(float(q), 5),
            "pvalue": round(float(p), 8),
            "overlap": f"{t['_k']}/{t['_n']}",
            "odds_ratio": t["odds_ratio"],
            "combined_score": combined,
            "genes": t["_genes"][:10],
        })

    out = [r for r in out if r["padj"] < padj_max]
    out.sort(key=lambda r: r["padj"])
    return out[:top]


def local_ora_for_databases(query_genes, db_library_map: dict, background_genes,
                            *, padj_max: float = 0.05, min_overlap: int = 1,
                            top: int = 20):
    """Run local ORA for several databases.

    ``db_library_map``: ``{display_label: gmt_library_name}`` (e.g.
    ``{"GO_BP": "GO_Biological_Process_2021"}``). Returns
    ``(results_by_db, versions_by_db, missing_labels)`` where ``missing_labels``
    lists databases with no local versioned GMT available."""
    results: dict = {}
    versions: dict = {}
    missing: list = []
    for label, library_name in db_library_map.items():
        lib = load_local_library(library_name)
        if lib is None:
            missing.append(label)
            continue
        gene_sets, version = lib
        results[label] = run_ora(
            query_genes, gene_sets, background_genes,
            padj_max=padj_max, min_overlap=min_overlap, top=top,
        )
        versions[label] = version
    return results, versions, missing
