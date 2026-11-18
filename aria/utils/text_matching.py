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
