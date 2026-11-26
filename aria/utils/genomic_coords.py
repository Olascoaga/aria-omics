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

import re
from typing import Any

_STANDARD_AUTOSOMES = frozenset(str(i) for i in range(1, 23))
_PEAK_HYPHEN_RE = re.compile(r"^(.+)-(\d+)-(\d+)$")
_PEAK_UNDERSCORE_RE = re.compile(r"^(.+)_(\d+)_(\d+)$")


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
