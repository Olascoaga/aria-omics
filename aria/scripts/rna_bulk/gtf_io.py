"""Bulk RNA-seq GTF / symbol-map / IO helpers (P2-8 split from rna_bulk_de.py).

Behavior-preserving extraction: these functions are re-exported from
`aria.scripts.rna_bulk_de` so the public surface is unchanged."""

from __future__ import annotations

from pathlib import Path


def _load_symbol_map(files: list, warnings: list,
                       gtf_hint: str = None) -> dict:
    """
    Load Ensembl-ID → HGNC-symbol mapping.

    Strategy:
      1. Look for counts_with_symbols.tsv next to the counts file.
      2. If missing, auto-regenerate from a GTF (search standard locations
         + hint provided by the agent).
      3. If still no GTF found, fall back to {} and warn loudly.

    The auto-regeneration matters because resume logic (v3.1+) skips
    featureCounts when counts_matrix.tsv exists — but if the matrix was
    written by an older version of ARIA (before v3.3), the symbols file
    won't exist. This function patches that gap without re-running
    featureCounts.
    """
    if not files:
        return {}
    try:
        import pandas as pd
        counts_path = Path(files[0])
        sym_path    = counts_path.parent / "counts_with_symbols.tsv"

        # ── Path A: symbols file already exists ──────────────────────
        if sym_path.exists():
            df = pd.read_csv(sym_path, sep="\t", index_col=0)
            if "gene_symbol" not in df.columns:
                warnings.append(
                    "counts_with_symbols.tsv missing 'gene_symbol' "
                    "column — regenerating."
                )
            else:
                clean = {}
                for k, v in df["gene_symbol"].items():
                    if isinstance(k, str) and isinstance(v, str):
                        clean[k.split(".")[0]] = v
                return clean

        # ── Path B: regenerate from GTF ──────────────────────────────
        gtf_path = _locate_gtf(counts_path, gtf_hint)
        if not gtf_path:
            warnings.append(
                "counts_with_symbols.tsv missing AND no GTF found in "
                "standard locations (~/.aria/genomes/*/annotation.gtf*). "
                "Pathway enrichment will use Ensembl IDs (likely 0 matches)."
            )
            return {}

        warnings.append(
            f"Auto-generating gene-symbol map from GTF: {gtf_path.name}"
        )
        mapping = _gtf_to_symbol_map(str(gtf_path))
        if not mapping:
            warnings.append(
                "GTF parse returned 0 mappings — check gene_name "
                "annotations are present."
            )
            return {}

        # Persist the regenerated symbols file for next time
        try:
            counts_df = pd.read_csv(counts_path, sep="\t", index_col=0)
            counts_with_sym = counts_df.copy()
            counts_with_sym.insert(
                0, "gene_symbol",
                [mapping.get(str(g).split(".")[0], str(g))
                 for g in counts_with_sym.index]
            )
            counts_with_sym.to_csv(str(sym_path), sep="\t")
            warnings.append(
                f"Wrote {sym_path.name} ({len(mapping):,} symbol mappings) "
                f"for future runs."
            )
        except Exception as e:
            warnings.append(f"Could not persist symbols file: {e}")

        return mapping

    except Exception as e:
        warnings.append(f"Symbol map load failed: {e}")
        return {}

def _locate_gtf(counts_path: Path, hint: str = None) -> Path | None:
    """
    Find a GTF file in standard locations.

    Search order:
      1. Explicit hint from the caller
      2. Sibling of the counts file (~/.../counts/annotation.gtf)
      3. ~/.aria/genomes/*/annotation.gtf{,.gz}
      4. ~/.aria/genomes/hg38/annotation.gtf etc.
    """
    candidates = []
    if hint:
        candidates.append(Path(hint))

    # Sibling of counts file
    candidates.extend([
        counts_path.parent / "annotation.gtf",
        counts_path.parent / "annotation.gtf.gz",
        counts_path.parent.parent / "annotation.gtf",
    ])

    # ~/.aria/genomes/*/annotation.gtf
    aria_genomes = Path.home() / ".aria" / "genomes"
    if aria_genomes.exists():
        for genome_dir in aria_genomes.iterdir():
            if genome_dir.is_dir():
                candidates.extend([
                    genome_dir / "annotation.gtf",
                    genome_dir / "annotation.gtf.gz",
                ])

    for c in candidates:
        if c.exists() and c.stat().st_size > 100_000:
            return c
    return None

def _gtf_to_symbol_map(gtf_file: str) -> dict:
    """
    Parse a GTF for {ensembl_id (no version): gene_symbol}.
    Streaming, low memory. Same logic as rna_quantify._build_ensembl_to_symbol_map.
    """
    import gzip as _gzip
    import re as _re

    mapping = {}
    opener = _gzip.open if str(gtf_file).endswith(".gz") else open
    try:
        with opener(gtf_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attrs = fields[8]
                gid_match = _re.search(r'gene_id "([^"]+)"', attrs)
                sym_match = _re.search(r'gene_name "([^"]+)"', attrs)
                if gid_match and sym_match:
                    gid = gid_match.group(1).split(".")[0]
                    mapping[gid] = sym_match.group(1)
    except Exception:
        return {}

    return mapping

def _load_gene_annotation(files: list, warnings: list) -> dict:
    """
    Load {biotype_map, length_map} from the GTF — used for DR filtering
    (protein_coding only) and TPM computation.

    biotype_map: {ensembl_id_no_version: biotype_string}
    length_map:  {ensembl_id_no_version: union_exon_length_bp}

    Exon-union length is the gene-model-appropriate length for TPM
    (sums non-overlapping exon regions). Falls back to gene coordinate
    length if exon parsing fails.

    Returns empty dicts gracefully if GTF not locatable.
    """
    if not files:
        return {"biotype": {}, "length": {}}
    try:
        counts_path = Path(files[0])
        gtf_path    = _locate_gtf(counts_path)
        if not gtf_path:
            warnings.append(
                "GTF not found for biotype/length annotation — "
                "DR will use all biotypes; TPM will not be computed."
            )
            return {"biotype": {}, "length": {}}

        biotype_map, length_map = _parse_gtf_biotype_and_length(str(gtf_path))
        warnings.append(
            f"Gene annotation loaded from {gtf_path.name}: "
            f"{len(biotype_map):,} biotypes, {len(length_map):,} lengths."
        )
        return {"biotype": biotype_map, "length": length_map}
    except Exception as e:
        warnings.append(f"Gene annotation load failed (non-fatal): {e}")
        return {"biotype": {}, "length": {}}

def _parse_gtf_biotype_and_length(gtf_file: str) -> tuple:
    """
    Single pass through GTF to extract both biotype (from gene lines)
    and exon-union length (sum of non-overlapping exon spans per gene).

    Returns (biotype_map, length_map).
    """
    import gzip as _gzip
    import re as _re
    from collections import defaultdict

    biotype_map = {}
    exons_by_gene = defaultdict(list)  # gid -> list of (start, end)

    opener = _gzip.open if str(gtf_file).endswith(".gz") else open
    try:
        with opener(gtf_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                feature = fields[2]
                if feature not in ("gene", "exon"):
                    continue

                attrs = fields[8]
                gid_match = _re.search(r'gene_id "([^"]+)"', attrs)
                if not gid_match:
                    continue
                gid = gid_match.group(1).split(".")[0]

                if feature == "gene":
                    bio_match = (_re.search(r'gene_biotype "([^"]+)"', attrs)
                                  or _re.search(r'gene_type "([^"]+)"', attrs))
                    if bio_match:
                        biotype_map[gid] = bio_match.group(1)

                elif feature == "exon":
                    try:
                        start = int(fields[3])
                        end   = int(fields[4])
                        exons_by_gene[gid].append((start, end))
                    except ValueError:
                        continue
    except Exception:
        return biotype_map, {}

    # Compute union-exon length per gene (merge overlapping intervals)
    length_map = {}
    for gid, intervals in exons_by_gene.items():
        if not intervals:
            continue
        # Merge overlapping intervals
        intervals.sort()
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        length_map[gid] = sum(e - s + 1 for s, e in merged)

    return biotype_map, length_map

def _to_symbols(gene_ids: list, symbol_map: dict) -> list:
    """
    Convert Ensembl IDs to HGNC symbols. Genes without a symbol are
    dropped (Enrichr won't match them anyway). Returns deduplicated list.
    """
    if not symbol_map:
        return list(gene_ids)
    seen = set()
    out  = []
    for gid in gene_ids:
        if not isinstance(gid, str):
            continue
        clean = gid.split(".")[0]   # strip version
        sym = symbol_map.get(clean)
        # If not in map, check if input might already be a symbol
        if not sym:
            if not clean.startswith(("ENSG", "ENSMUSG", "ENS")):
                sym = clean   # already a symbol
            else:
                continue       # unmapped Ensembl ID — skip
        if sym not in seen:
            seen.add(sym)
            out.append(sym)
    return out
