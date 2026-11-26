"""Deterministic word- and phrase-boundary keyword matching helpers.

Replaces naive ``any(k in question for k in keywords)`` substring checks that
caused false positives across biological prompts (e.g. ``"fine"`` matching
``"define"``, ``"type"`` matching ``"genotype"`` / ``"wild-type"``, ``"dot"``
matching ``"dotplot"``, ``"factor"`` matching ``"transcription factor"``, and
``"tf "`` missing ``"TF,"`` or matching ``"wtf "``).
"""

from __future__ import annotations

from functools import lru_cache
import re
from typing import Iterable

# Characters treated as part of a single token (including hyphen so compound
# tokens like "wild-type" do not match standalone "type", while explicit
# hyphenated keywords like "sub-tad" or "hi-c loop" still match verbatim).
_TOKEN_CHAR_CLASS = r"0-9A-Za-z_À-ÖØ-öø-ÿ-"


@lru_cache(maxsize=256)
def _compile_boundary_pattern(keywords: tuple[str, ...]) -> re.Pattern[str] | None:
    cleaned = []
    for raw in keywords:
        kw = str(raw or "").strip().lower()
        if not kw:
            continue
        # Collapse internal whitespace to flexible \s+ while escaping literals.
        parts = [re.escape(part) for part in kw.split()]
        cleaned.append(r"\s+".join(parts))
    if not cleaned:
        return None
    # Sort longest-first so multi-word / hyphenated phrases match cleanly.
    cleaned.sort(key=len, reverse=True)
    alternation = "|".join(cleaned)
    pattern = (
        rf"(?:^|(?<![{_TOKEN_CHAR_CLASS}]))"
        rf"(?:{alternation})"
        rf"(?:$|(?![{_TOKEN_CHAR_CLASS}]))"
    )
    return re.compile(pattern, flags=re.IGNORECASE)


def has_keyword_boundary(text: str | None, keywords: Iterable[str]) -> bool:
    """Return True iff *text* contains any keyword/phrase at token boundaries."""
    if not text:
        return False
    key_tuple = tuple(str(k) for k in keywords if k)
    if not key_tuple:
        return False
    regex = _compile_boundary_pattern(key_tuple)
    if regex is None:
        return False
    return bool(regex.search(str(text)))


_GENE_SKIP_TOKENS = frozenset({
    "",
    "?",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "<na>",
    "intergenic",
    "distal",
    "unannotated",
    "unspecified",
    ".",
    "-",
})
_RAW_ENSEMBL_ID_RE = re.compile(
    r"^(?:ENSG|ENSMUSG|ENSMUST|ENSMUSP|ENST|ENSP|ENSRNOG|ENSRNOT|ENS[A-Z]{1,4}[GTPE]|ENS)\d+(?:\.\d+)?$",
    re.IGNORECASE,
)
_PEAK_COORD_SYMBOL_RE = re.compile(
    r"^(?:(?:chr[\w.\-]+|[\dXYM]+|MT):\d+[-:]\d+|chr[\w.\-]+-\d+-\d+)$",
    re.IGNORECASE,
)


def clean_gene_symbol(raw) -> str | None:
    """Return a clean HGNC/MGI gene symbol, or ``None`` for raw IDs, peak coords, or sentinels.

    Rejects:
      - missing / whitespace-only / sentinel strings (``?``, ``nan``, ``None``, ``NA``, ``null``, etc.);
      - bare numeric Entrez/integer IDs (e.g. ``"7157"``);
      - raw Ensembl identifiers (e.g. ``"ENSG00000141510"``, ``"ENSMUSG00000059552"``);
      - genomic peak coordinates (``chr1:100-200``, ``1:100-200``, ``chr1-100-200``);
      - tokens containing whitespace, delimiters, or HTML angle brackets.
    """
    if raw is None or isinstance(raw, bool):
        return None
    sym = str(raw).strip()
    if not sym or sym.lower() in _GENE_SKIP_TOKENS or sym.isdigit():
        return None
    if _RAW_ENSEMBL_ID_RE.match(sym) or _PEAK_COORD_SYMBOL_RE.match(sym):
        return None
    if any(ch in sym for ch in (" ", "\t", "\n", "\r", ";", ",", "<", ">")):
        return None
    return sym
