"""Ledger quarantine for the SPECULATIVE tier (ADR-057 rail #6).

A machine hypothesis gets its own ledger namespace, ``hypothesis://<id>``, linked
to the audited evidence it cites — but mechanically marked so NOTHING downstream
(the audited claim manifest, the RO-Crate export, ``aria diff``, a future run)
can ever promote it to a claim. Zero feedback loop: a speculation never becomes a
premise.

The wall is by construction, and this module makes it explicit:
  - ``hypothesis://`` is a distinct namespace from the audited ``ledger://``
    nodes and from any claim id, so claim<->ledger linkage
    (``run_ledger.link_claims_to_ledger``) cannot resolve a claim to a hypothesis
    node.
  - ``SPECULATIVE`` is not a Claim Compiler tier (``claim_compiler.TIER_ORDER``)
    and not in ``run_ledger._NON_DESCRIPTIVE_TIERS``, so it can never be rendered
    as an associative-or-stronger claim.
  - ``assert_no_speculative_promotion`` is the active enforcer the report /
    RO-Crate assembly (S9) calls to prove no audited claim carried a hypothesis
    node or the SPECULATIVE tier.
"""

from __future__ import annotations

import re
from typing import Any, Iterable

from .types import Hypothesis

SPECULATIVE_TIER = "SPECULATIVE"
_HYPOTHESIS_SCHEME = "hypothesis://"


class SpeculativePromotionError(ValueError):
    """A speculative hypothesis leaked into the audited claim path."""


# A quarantine node id is part of a URL-like namespace; keep it to a safe,
# stable slug so an LLM-authored id cannot inject scheme/path characters.
_UNSAFE_ID_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def _slug_id(raw: Any) -> str:
    """Normalize a hypothesis id to a safe, non-empty slug."""
    slug = _UNSAFE_ID_CHARS.sub("-", str(raw or "").strip()).strip("-")
    return slug or "unnamed"


def quarantine_node_id(hyp: Hypothesis) -> str:
    """Deterministic ``hypothesis://<slug>`` quarantine node id for a hypothesis.

    Slugged but NOT deduplicated here — uniqueness across a batch is enforced by
    :func:`quarantine_hypotheses`, which is the only place that can see the other
    hypotheses' ids.
    """
    return f"{_HYPOTHESIS_SCHEME}{_slug_id(getattr(hyp, 'id', ''))}"


def is_hypothesis_node(node_id: Any) -> bool:
    """True iff ``node_id`` is a quarantined hypothesis node."""
    return isinstance(node_id, str) and node_id.startswith(_HYPOTHESIS_SCHEME)


def quarantine_hypotheses(hypotheses: list[Hypothesis]) -> list[dict]:
    """Assign each hypothesis its ``hypothesis://`` node and build the manifest.

    Stamps ``Hypothesis.ledger_node`` and returns a non-promotable quarantine
    manifest linking each hypothesis node to the audited evidence it cites
    (``observation_refs`` / ``entities``). Every entry is explicitly
    ``promotable: False``.
    """
    manifest: list[dict] = []
    seen: set[str] = set()
    for hyp in hypotheses or []:
        node = quarantine_node_id(hyp)
        # Enforce uniqueness across the batch: two hypotheses with the same id
        # (or ids that slug to the same value) must NOT collapse onto one
        # quarantine node, or a hypothesis would silently inherit another's
        # provenance/cites (H1, bug 3). Suffix deterministically on collision.
        if node in seen:
            base, suffix = node, 2
            while node in seen:
                node = f"{base}-{suffix}"
                suffix += 1
        seen.add(node)
        hyp.ledger_node = node
        manifest.append(
            {
                "hypothesis_node": node,
                "tier": SPECULATIVE_TIER,
                "cites": list(hyp.observation_refs or []),
                "entities": list(hyp.entities or []),
                "promotable": False,
            }
        )
    return manifest


def assert_no_speculative_promotion(claims: Iterable[dict]) -> None:
    """Raise if any audited claim carries a hypothesis node or the SPECULATIVE tier.

    The active wall the report / RO-Crate assembly calls so a speculation can
    never be emitted as an audited claim. Raises
    :class:`SpeculativePromotionError` on the first contamination.
    """
    for claim in claims or []:
        if not isinstance(claim, dict):
            continue
        if str(claim.get("tier") or "") == SPECULATIVE_TIER:
            raise SpeculativePromotionError(
                f"claim {claim.get('claim_id')!r} carries the SPECULATIVE tier; "
                "machine hypotheses may not enter the audited claim manifest"
            )
        for key in ("node_id", "ledger_node_id", "ledger_node"):
            if is_hypothesis_node(claim.get(key)):
                raise SpeculativePromotionError(
                    f"claim {claim.get('claim_id')!r} cites quarantined node "
                    f"{claim.get(key)!r}; hypotheses are non-promotable"
                )
