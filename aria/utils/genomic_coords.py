"""Canonical genomic contig normalization and BED/GTF coordinate conversions.

Conventions
-----------
- **BED / internal interval convention:** 0-based, half-open ``[start, end)``
  (``0 <= start < end``, width ``end - start``). Used by BED, narrowPeak /
  broadPeak, BAM/CRAM ``reference_start``/``reference_end``, FASTA ``.fai``
  byte offsets, and AnnData peak coordinates (``chr1:start-end``).
- **GTF / GFF convention:** 1-based, closed ``[start_1based, end_1based]``
  (``1 <= start_1based <= end_1based``, width ``end_1based - start_1based + 1``).
- **Single-conversion rule:** Convert GTF 1-based closed coordinates to 0-based
  half-open ``[start_1based - 1, end_1based)`` *once* at parse time via
  :func:`gtf_to_bed_interval` and :func:`gtf_tss_0based`. Never subtract 1 from
  coordinates that are already 0-based BED intervals (use
  :func:`bed_interval_clamped` when clamping to a contig length).
- **Contig canonicalization:** :func:`canonical_contig` reconciles UCSC
  (``chr1``..``chr22``, ``chrX``, ``chrY``, ``chrM``/``chrMT``) and Ensembl
  (``1``..``22``, ``X``, ``Y``, ``M``/``MT``) aliases into canonical ``chr*``
  names while preserving unplaced scaffolds verbatim. :func:`ensembl_contig`
  returns the stripped Ensembl/``MT`` key for callers whose public schema uses
  bare chromosome names.

Boundary invariant: pure Python standard library only; safe to import from any
ARIA script, utility, or agent environment.
"""

from __future__ import annotations

from collections import defaultdict
import gzip
from pathlib import Path
import re
from typing import Any

_STANDARD_AUTOSOMES = frozenset(str(i) for i in range(1, 23))
_PEAK_HYPHEN_RE = re.compile(r"^(.+)-(\d+)-(\d+)$")
_PEAK_UNDERSCORE_RE = re.compile(r"^(.+)_(\d+)_(\d+)$")
_ENS_VERSION_RE = re.compile(r"^(ENS[^.]+)\.\d+$")
_GTF_ATTR_RE = re.compile(r'(?:^|;)\s*([\w]+)\s*(?:=|\s)\s*(?:"([^"]*)"|([^;]+))')


def canonical_contig(chrom: str | None) -> str | None:
    """Normalize UCSC (``chr1``, ``chrM``) and Ensembl (``1``, ``MT``) aliases.

    Standard nuclear chromosomes (``1..22``, ``X``, ``Y``) and mitochondrial
    aliases (``M``, ``MT``, ``chrM``, ``chrMT``) map to canonical ``chr*`` keys
    (e.g. ``"chr1"``, ``"chrX"``, ``"chrM"``), while unplaced scaffolds and
    non-standard contigs are preserved stripped.
    """
    if chrom is None:
        return None
    raw = str(chrom).strip()
    if not raw:
        return None
    low = raw.lower()
    if low in ("m", "mt", "chrm", "chrmt"):
        return "chrM"
    stripped = raw[3:] if low.startswith("chr") else raw
    up = stripped.upper()
    if stripped in _STANDARD_AUTOSOMES or up in ("X", "Y"):
        return f"chr{up}"
    if stripped.isdigit() and 1 <= int(stripped) <= 99:
        return f"chr{int(stripped)}"
    return raw


def is_standard_contig(chrom: str | None) -> bool:
    """Return True if ``chrom`` maps to ``chr1..chr22``, ``chrX``, or ``chrY``."""
    canon = canonical_contig(chrom)
    if canon is None or not canon.startswith("chr"):
        return False
    body = canon[3:]
    return body in _STANDARD_AUTOSOMES or body in ("X", "Y")


def ensembl_contig(chrom: str | None) -> str:
    """Normalize standard UCSC/Ensembl contig aliases to bare Ensembl form (``1``, ``X``, ``MT``).

    Unplaced scaffolds and alternate contigs (e.g. ``chrUn_KI270442v1``) are
    preserved verbatim rather than having their ``chr`` prefix stripped.
    """
    canon = canonical_contig(chrom)
    if canon is None:
        return ""
    if canon == "chrM":
        return "MT"
    if canon.startswith("chr"):
        body = canon[3:]
        if (
            body in _STANDARD_AUTOSOMES
            or body in ("X", "Y")
            or (body.isdigit() and 1 <= int(body) <= 99)
        ):
            return body
    return canon


def contig_aliases(chrom: str | None) -> tuple[str, ...]:
    """Return ordered unique UCSC/Ensembl aliases for ``chrom``."""
    raw = str(chrom or "").strip()
    if not raw:
        return ()
    canon = canonical_contig(raw)
    ens = ensembl_contig(raw)
    cands: list[str] = [raw]
    if canon == "chrM":
        cands.extend(("chrM", "chrMT", "MT", "M"))
    else:
        if canon:
            cands.append(canon)
        if ens:
            cands.append(ens)
    return tuple(dict.fromkeys(c for c in cands if c))


def resolve_contig_in_valid(chrom: str | None, valid_chroms) -> str | None:
    """Resolve ``chrom`` against ``valid_chroms`` using UCSC/Ensembl canonical aliasing."""
    if not chrom or not valid_chroms:
        return None
    for cand in contig_aliases(chrom):
        if cand in valid_chroms:
            return cand
    target = canonical_contig(chrom)
    if target is None:
        return None
    for existing in valid_chroms:
        if canonical_contig(existing) == target:
            return existing
    return None


def gtf_to_bed_interval(start_1based: int, end_1based: int) -> tuple[int, int]:
    """Convert a 1-based closed GTF interval ``[start, end]`` to 0-based half-open BED ``[start-1, end)``."""
    s1 = int(start_1based)
    e1 = int(end_1based)
    if s1 < 1 or e1 < s1:
        raise ValueError(
            f"Invalid 1-based closed GTF interval: [{start_1based}, {end_1based}]"
        )
    return s1 - 1, e1


def gtf_tss_0based(start_1based: int, end_1based: int, strand: str) -> int:
    """Return the 0-based TSS coordinate for a 1-based closed GTF gene interval."""
    s0, e0 = gtf_to_bed_interval(start_1based, end_1based)
    return e0 - 1 if str(strand).strip() == "-" else s0


def bed_interval_clamped(
    start_0based: int, end_0based: int, seq_len: int
) -> tuple[int, int]:
    """Clamp an already 0-based half-open BED interval ``[start, end)`` to ``[0, seq_len]``."""
    length = max(0, int(seq_len))
    s0 = max(0, min(length, int(start_0based)))
    e0 = max(s0, min(length, int(end_0based)))
    return s0, e0


def parse_peak_coord(
    name: Any,
    *,
    allow_reversed: bool = True,
    require_str: bool = False,
) -> tuple[str, int, int] | None:
    """Parse a genomic peak coordinate string into ``(chrom, start, end)``.

    Supports UCSC/Ensembl ``chr1:100-200`` / ``1:100-200`` / ``chr1:100_200``,
    GRCh38 unplaced scaffolds (``chr1_KI270706v1_random:1000-1500``), Signac/10x
    hyphen-delimited ``chr1-1100-1500``, underscore-delimited ``chr1_1000_1500``,
    and BED tab lines (``chr1\\t100\\t200``). Rejects negative coordinates,
    whitespace contigs, and malformed extra-colon strings.
    """
    if name is None:
        return None
    if require_str and not isinstance(name, str):
        return None
    text = str(name).strip()
    if not text:
        return None
    try:
        if "\t" in text:
            parts = text.split("\t")
            if len(parts) < 3:
                return None
            chrom, s, e = parts[0].strip(), int(parts[1]), int(parts[2])
        elif ":" in text:
            if text.count(":") != 1:
                return None
            chrom, pos = text.split(":", 1)
            chrom = chrom.strip()
            pos = pos.strip()
            if "-" in pos:
                s_txt, e_txt = pos.split("-", 1)
            elif "_" in pos:
                s_txt, e_txt = pos.split("_", 1)
            else:
                return None
            s, e = int(s_txt), int(e_txt)
        else:
            m = _PEAK_HYPHEN_RE.match(text) or _PEAK_UNDERSCORE_RE.match(text)
            if not m:
                return None
            chrom = m.group(1).strip()
            s, e = int(m.group(2)), int(m.group(3))
    except ValueError:
        return None
    if not chrom or any(ch.isspace() for ch in chrom):
        return None
    if s < 0 or e < 0:
        return None
    if e < s:
        if not allow_reversed:
            return None
        s, e = e, s
    return chrom, s, e


def strip_gene_version(gene: str) -> str:
    """Strip numeric Ensembl versions (e.g. ``ENSG00000141510.17``) without truncating dotted IDs."""
    return _ENS_VERSION_RE.sub(r"\1", str(gene).strip())


def parse_gtf_attributes(text: str) -> dict[str, str]:
    """Parse quoted GTF2.2 (``key "val"``) or GFF3 (``key=val``) attribute strings."""
    return {
        m.group(1): (m.group(2) if m.group(2) is not None else m.group(3)).strip()
        for m in _GTF_ATTR_RE.finditer(text)
    }


def gtf_has_tag(attrs: str, tag_value: str) -> bool:
    """True when GTF/GFF3 attributes contain ``tag "<tag_value>"`` or ``tag=<tag_value>``."""
    target = tag_value.lower()
    for chunk in attrs.split(";"):
        item = chunk.strip()
        low = item.lower()
        if low.startswith("tag ") or low.startswith("tag="):
            val = item.split(None, 1)[-1] if " " in item else item.split("=", 1)[-1]
            if val.strip().strip('"').lower() == target:
                return True
    return False


def is_par_y_record(chrom: str, gene_id: str, attrs: str) -> bool:
    """True for pseudoautosomal chrY duplicates (``_PAR_Y`` suffix or ``PAR`` tag on ``chrY``)."""
    return gene_id.upper().endswith("_PAR_Y") or (
        gtf_has_tag(attrs, "PAR") and canonical_contig(chrom) == "chrY"
    )


def parse_gtf_symbol_map(gtf_path: str | Path | None) -> dict[str, str]:
    """Stream a GTF/GFF3 into ``{gene_id_without_ensembl_version: gene_symbol}``.

    Prefers ``gene`` rows over ``transcript``/``exon`` rows while supporting
    transcript-only RefSeq/UCSC GTFs and GFF3 ``=`` attributes. Preserves dotted
    non-Ensembl IDs (e.g. ``LOC100.1``) verbatim.
    """
    if not gtf_path:
        return {}
    p = Path(str(gtf_path))
    if not p.is_file():
        return {}
    opener = gzip.open if str(p).lower().endswith(".gz") else open
    mapping: dict[str, str] = {}
    try:
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feat = fields[2].lower()
                if feat not in ("gene", "transcript", "mrna", "exon"):
                    continue
                attrs = parse_gtf_attributes(fields[8])
                gid = attrs.get("gene_id") or attrs.get("ID")
                symbol = attrs.get("gene_name") or attrs.get("gene_symbol") or attrs.get("Name")
                if gid and symbol:
                    gid_clean = strip_gene_version(gid)
                    if feat == "gene" or gid_clean not in mapping:
                        mapping[gid_clean] = symbol
    except (OSError, EOFError):
        return {}
    return mapping


def parse_gtf_tss_records(
    gtf_path: str | Path | None,
    wanted: set[str] | None = None,
    *,
    allow_unnamed: bool = False,
    require_strand: bool = False,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Parse canonical 0-based per-gene TSS records from a local GTF/GFF3 (``.gz`` ok).

    Single-conversion invariant (GTF2.2/GFF3 -> 0-based half-open BED):
      - ``start``, ``end`` are 0-based half-open ``[start_1based - 1, end_1based)``.
      - ``tss`` is 0-based ``start_1based - 1`` (``+`` strand) or ``end_1based - 1`` (``-``).
      - Falls back to ``transcript`` / ``mRNA`` rows when ``gene`` rows are absent.
      - Drops pseudoautosomal Y duplicates (``_PAR_Y`` / ``PAR`` tag on ``chrY``).
      - Deterministically resolves duplicate symbols by
        ``(MANE_Select/Ensembl_canonical, protein_coding, basic, span, smallest gene_id)``.
    """
    stats: dict[str, Any] = {
        "n_gene_rows_scanned": 0,
        "n_par_y_dropped": 0,
        "n_duplicate_symbols_resolved": 0,
        "n_unique_genes": 0,
    }
    if not gtf_path or not Path(str(gtf_path)).is_file():
        return {}, stats

    opener = gzip.open if str(gtf_path).lower().endswith(".gz") else open
    best_by_gene: dict[str, tuple[tuple[int, int, int, int], str, dict[str, Any]]] = {}

    wanted_set = set(wanted) if wanted is not None else None
    wanted_no_ver: dict[str, str] = {}
    wanted_upper_unique: dict[str, str] = {}
    if wanted_set is not None:
        no_ver_buckets: dict[str, list[str]] = {}
        upper_buckets: dict[str, list[str]] = {}
        for w in wanted_set:
            ws = str(w).strip()
            if not ws:
                continue
            base = ws.split(".", 1)[0]
            no_ver_buckets.setdefault(base, []).append(ws)
            upper_buckets.setdefault(ws.upper(), []).append(ws)
        wanted_no_ver = {k: v[0] for k, v in no_ver_buckets.items() if len(v) == 1}
        wanted_upper_unique = {k: v[0] for k, v in upper_buckets.items() if len(v) == 1}

    def _scan_features(target_features: frozenset[str]) -> None:
        with opener(str(gtf_path), "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2].lower() not in target_features:
                    continue
                stats["n_gene_rows_scanned"] += 1
                attrs_raw = parts[8]
                attrs = parse_gtf_attributes(attrs_raw)
                gene_id = attrs.get("gene_id") or attrs.get("ID") or ""
                gene_name = attrs.get("gene_name") or attrs.get("gene_symbol") or attrs.get("Name") or gene_id
                try:
                    start_1, end_1 = int(parts[3]), int(parts[4])
                    start_0, end_0 = gtf_to_bed_interval(start_1, end_1)
                except ValueError:
                    continue
                if require_strand and parts[6] not in {"+", "-"}:
                    continue
                strand = "-" if parts[6] == "-" else "+"
                tss_0 = gtf_tss_0based(start_1, end_1, strand)
                if not gene_name:
                    if not allow_unnamed:
                        continue
                    gene_name = f"{parts[0].strip()}:{tss_0}:{strand}"

                if wanted_set is not None:
                    if gene_name in wanted_set:
                        key_symbol = gene_name
                    elif gene_id and gene_id in wanted_set:
                        key_symbol = gene_id
                    elif gene_id and gene_id.split(".", 1)[0] in wanted_no_ver:
                        key_symbol = wanted_no_ver[gene_id.split(".", 1)[0]]
                    elif gene_name.upper() in wanted_upper_unique:
                        key_symbol = wanted_upper_unique[gene_name.upper()]
                    elif gene_id and gene_id.upper() in wanted_upper_unique:
                        key_symbol = wanted_upper_unique[gene_id.upper()]
                    else:
                        continue
                else:
                    key_symbol = gene_name

                if is_par_y_record(parts[0], gene_id, attrs_raw):
                    stats["n_par_y_dropped"] += 1
                    continue

                span = max(0, end_0 - start_0)
                biotype = (attrs.get("gene_type") or attrs.get("gene_biotype") or "").lower()
                is_mane = int(
                    gtf_has_tag(attrs_raw, "MANE_Select")
                    or gtf_has_tag(attrs_raw, "Ensembl_canonical")
                )
                is_pc = int(biotype == "protein_coding")
                is_basic = int(gtf_has_tag(attrs_raw, "basic"))
                prio = (is_mane, is_pc, is_basic, span)

                rec = {
                    "gene": key_symbol,
                    "gene_id": gene_id or key_symbol,
                    "chrom": parts[0].strip(),
                    "canon_chrom": canonical_contig(parts[0]),
                    "start": start_0,
                    "end": end_0,
                    "tss": tss_0,
                    "strand": strand,
                    "gene_biotype": biotype or None,
                }
                prev = best_by_gene.get(key_symbol)
                if prev is None:
                    best_by_gene[key_symbol] = (prio, rec["gene_id"], rec)
                else:
                    stats["n_duplicate_symbols_resolved"] += 1
                    prev_prio, prev_gid, _ = prev
                    if prio > prev_prio or (prio == prev_prio and rec["gene_id"] < prev_gid):
                        best_by_gene[key_symbol] = (prio, rec["gene_id"], rec)

    _scan_features(frozenset({"gene"}))
    if not best_by_gene:
        _scan_features(frozenset({"transcript", "mrna"}))

    records = {k: v[2] for k, v in best_by_gene.items()}
    stats["n_unique_genes"] = len(records)
    return records, stats


def parse_gtf_gene_models(gtf_path: str | Path) -> dict[str, list[dict[str, Any]]]:
    """Stream a GTF/GFF3 into per-chromosome 0-based BED gene models with exons.

    Raises on read/decompression errors (fail-closed, D1-021). Falls back to
    ``transcript`` / ``mRNA`` rows when explicit ``gene`` rows are absent, and
    drops pseudoautosomal ``_PAR_Y`` duplicates on ``chrY``.
    """
    genes: dict[str, dict[str, Any]] = {}
    tx_genes: dict[str, dict[str, Any]] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    opener = gzip.open if str(gtf_path).lower().endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feat = fields[2].lower()
            if feat not in ("gene", "transcript", "mrna", "exon"):
                continue
            attrs_raw = fields[8]
            attrs = parse_gtf_attributes(attrs_raw)
            gid = attrs.get("gene_id") or attrs.get("ID")
            if not gid:
                continue
            if is_par_y_record(fields[0], gid, attrs_raw):
                continue
            try:
                start_1, end_1 = int(fields[3]), int(fields[4])
                start_0, end_0 = gtf_to_bed_interval(start_1, end_1)
            except ValueError:
                continue
            if feat == "exon":
                exons[gid].append((start_0, end_0))
                continue
            strand = fields[6]
            gname = attrs.get("gene_name") or attrs.get("gene_symbol") or attrs.get("Name") or strip_gene_version(gid)
            rec = {
                "gene": gname,
                "chrom": ensembl_contig(fields[0]),
                "strand": strand,
                "start": start_0,
                "end": end_0,
                "tss": gtf_tss_0based(start_1, end_1, strand),
            }
            if feat == "gene":
                genes[gid] = rec
            else:
                prev = tx_genes.get(gid)
                if prev is None:
                    tx_genes[gid] = rec
                else:
                    s0 = min(prev["start"], start_0)
                    e0 = max(prev["end"], end_0)
                    prev["start"], prev["end"] = s0, e0
                    prev["tss"] = (e0 - 1) if prev["strand"] == "-" else s0

    merged_genes = dict(tx_genes)
    merged_genes.update(genes)

    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for gid, gene in merged_genes.items():
        gene["exons"] = sorted(exons.get(gid, []))
        by_chrom[gene["chrom"]].append(gene)
    return {
        chrom: sorted(glist, key=lambda g: g["tss"])
        for chrom, glist in by_chrom.items()
    }
