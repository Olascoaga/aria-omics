"""Non-same-cell cross-modal concordance evidence adapter (ADR-080; ADR-028 S2).

Pure and deterministic — zero LLM, zero biological priors. Flattens the audited
``cross_modal_concordance`` IPC result (bulk RNA x bulk ATAC, or unpaired
scRNA x scATAC cohorts) into ``EvidenceSignal`` items, gated strictly on the
``ledger://cross_modal/concordance`` node having actually run.

Two layers are emitted:

  1. ``sample_matched_links`` -> ``gene`` signals
     (``measure="sample_matched_peak2gene_corr"``): sample/donor-matched cis
     peak-gene correlations across separately profiled libraries of the SAME
     biological units. This is a genuinely new observation (it is not a
     re-expression of either single-assay contrast), so it carries no lineage
     key; it inherits ``PEAK2GENE_ASSOCIATIVE`` + ``CROSS_MODAL_UNPAIRED_COHORT``.
  2. ``tf_triangulation`` -> ``tf_motif`` signals
     (``measure="cross_modal_tf_concordance"``) for TFs whose chromatin motif /
     footprint shift agrees in sign with the TF gene's own RNA change.

Non-double-counting (ADR-072 ``independence_key`` discipline): a TF
concordance signal is ALGEBRAICALLY DERIVED from the upstream single-assay RNA
DE and chromatin motif/footprint signals of the same contrast — it is not a
third wet-lab assay. It therefore carries
``independence_key = chromatin_tf_lineage_key(tf)``, and ``gather_evidence``
stamps the same key on the chromatin motif/footprint signals of that TF, so a
hypothesis citing ``[RNA DE, ATAC motif, cross-modal TF concordance]`` counts
exactly two independent lines (RNA + ATAC), never three.
"""

from __future__ import annotations

import math
import re
from typing import Any

from aria.agents.narrative.run_ledger import (
    cross_modal_findings,
    cross_modal_status,
)

from .. import devils_advocate as _devils_advocate
from ..caveats import (
    CROSS_MODAL_UNPAIRED_COHORT,
    EXTRA_CONFOUND_TOKENS,
    FOOTPRINT_ASSOCIATIVE,
    LOW_REPLICATION,
    MOTIF_NOT_BINDING,
    PARALOG_MOTIF_AMBIGUITY,
    PEAK2GENE_ASSOCIATIVE,
)
from ..types import EvidenceSignal
from ._chromatin_common import node_ran

# Additive, idempotent registration of the new caveat's recognition tokens in
# the adversarial gate, so a hypothesis can acknowledge it in free text.
for _code, _tokens in EXTRA_CONFOUND_TOKENS.items():
    _devils_advocate._CONFOUND_TOKENS.setdefault(_code, _tokens)

CROSS_MODAL_NODE = "ledger://cross_modal/concordance"
SAMPLE_LINK_MEASURE = "sample_matched_peak2gene_corr"
TF_CONCORDANCE_MEASURE = "cross_modal_tf_concordance"
CHROMATIN_TF_LINEAGE_PREFIX = "chromatin_tf_lineage::"

_MAX_SAMPLE_LINK_GENES = 20
_MAX_TFS_PER_CONTRAST = 15
_DEFAULT_LINK_ALPHA = 0.05

_PEAK_COORD_RE = re.compile(
    r"^(?:chr[\w.\-]+|[\dXYM]+|MT):\d+[-:]\d+$", re.IGNORECASE
)
_RAW_ID_PREFIXES = ("ENSG", "ENSMUSG", "ENST", "ENSMUST")
_COMPOSITE_MARKERS = ("::", ",", "/", ";", "(", ")", " ")


def chromatin_tf_lineage_key(tf: str) -> str:
    """Shared lineage key for one TF's chromatin-derived evidence (ADR-072)."""
    return f"{CHROMATIN_TF_LINEAGE_PREFIX}{str(tf or '').strip().lower()}"


def _finite(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _corr(value: Any) -> float | None:
    out = _finite(value)
    if out is None or abs(out) > 1.0 + 1e-6:
        return None
    return max(-1.0, min(1.0, out))


def _gene_symbol(raw: Any) -> str | None:
    sym = str(raw or "").strip()
    if not sym or sym == "?" or sym.isdigit():
        return None
    if _PEAK_COORD_RE.match(sym) or sym.upper().startswith(_RAW_ID_PREFIXES):
        return None
    parts = sym.split("-")
    if len(parts) == 3 and parts[0].lower().startswith("chr") and all(
        p.isdigit() for p in parts[1:]
    ):
        return None
    if any(ch in sym for ch in (" ", "\t", "\n", ";", ",")):
        return None
    return sym


def _single_tf(raw: Any) -> str | None:
    """A clean single TF symbol; composite motif labels are discarded."""
    sym = str(raw or "").strip()
    if not sym or any(m in sym for m in _COMPOSITE_MARKERS):
        return None
    return _gene_symbol(sym)


def _engine_ran(findings: dict, run_ledger: dict | None) -> bool:
    if cross_modal_status(findings)[0] != "ran":
        return False
    if run_ledger is None:
        return True
    return node_ran(run_ledger, CROSS_MODAL_NODE)


def _sample_link_signals(findings: dict) -> list[EvidenceSignal]:
    links = findings.get("sample_matched_links") or {}
    if not isinstance(links, dict) or links.get("ran") is not True:
        return []
    n_shared = links.get("n_shared_samples")
    try:
        n_shared_int = int(n_shared)
    except (TypeError, ValueError):
        return []
    if n_shared_int <= 0:
        return []
    alpha = _finite(links.get("padj_threshold")) or _DEFAULT_LINK_ALPHA
    caveats = [PEAK2GENE_ASSOCIATIVE, CROSS_MODAL_UNPAIRED_COHORT]
    if links.get("low_power_warning") and LOW_REPLICATION not in caveats:
        caveats.append(LOW_REPLICATION)
    context = f"matched_samples_n={n_shared_int}"

    best: dict[str, tuple[float, str]] = {}
    order: list[str] = []
    for row in links.get("top_links") or []:
        if not isinstance(row, dict):
            continue
        gene = _gene_symbol(row.get("gene"))
        r = _corr(row.get("spearman_r"))
        if gene is None or r is None or r == 0.0:
            continue
        padj = _finite(row.get("padj"))
        if padj is not None and padj > alpha:
            continue
        key = gene.lower()
        if key not in best:
            order.append(key)
            best[key] = (r, gene)
        elif abs(r) > abs(best[key][0]):
            best[key] = (r, gene)

    out: list[EvidenceSignal] = []
    for key in order[:_MAX_SAMPLE_LINK_GENES]:
        r, gene = best[key]
        out.append(
            EvidenceSignal(
                entity=gene,
                entity_kind="gene",
                modality="cross_modal",  # type: ignore[arg-type]
                measure=SAMPLE_LINK_MEASURE,
                audited_node_ref=CROSS_MODAL_NODE,
                value=round(r, 4),
                direction="up" if r > 0 else "down",
                caveats_inherited=list(caveats),
                context=context,
            )
        )
    return out


def _tf_triangulation_blocks(findings: dict) -> list[tuple[str, dict]]:
    """(contrast_id, tf_triangulation) pairs; per-contrast first, else top-level."""
    out: list[tuple[str, dict]] = []
    for idx, mc in enumerate(findings.get("matched_contrasts") or []):
        if not isinstance(mc, dict):
            continue
        tri = mc.get("tf_triangulation")
        if isinstance(tri, dict) and tri.get("ran") is True:
            cid = str(mc.get("contrast_id") or f"contrast_{idx}")
            out.append((cid, tri))
    if not out:
        tri = findings.get("tf_triangulation")
        if isinstance(tri, dict) and tri.get("ran") is True:
            out.append((str(tri.get("contrast_id") or "matched_contrast"), tri))
    return out


def _tf_signals(findings: dict) -> list[EvidenceSignal]:
    out: list[EvidenceSignal] = []
    seen: set[tuple[str, str]] = set()
    for contrast_id, tri in _tf_triangulation_blocks(findings):
        n = 0
        for row in tri.get("top_concordant_tfs") or []:
            if not isinstance(row, dict):
                continue
            if str(row.get("status") or row.get("classification") or "").lower() != "concordant":
                continue
            tf = _single_tf(row.get("tf"))
            effect = _finite(
                row.get("chromatin_effect")
                if row.get("chromatin_effect") is not None
                else row.get("chromatin_score")
            )
            if tf is None or effect is None or effect == 0.0:
                continue
            key = (tf.lower(), contrast_id)
            if key in seen:
                continue
            seen.add(key)
            caveats = [MOTIF_NOT_BINDING, CROSS_MODAL_UNPAIRED_COHORT]
            if row.get("paralog_motif_ambiguity"):
                caveats.append(PARALOG_MOTIF_AMBIGUITY)
            if "footprint" in str(row.get("source_layer") or row.get("source") or "").lower():
                caveats.append(FOOTPRINT_ASSOCIATIVE)
            out.append(
                EvidenceSignal(
                    entity=tf,
                    entity_kind="tf_motif",
                    modality="cross_modal",  # type: ignore[arg-type]
                    measure=TF_CONCORDANCE_MEASURE,
                    audited_node_ref=CROSS_MODAL_NODE,
                    value=round(effect, 4),
                    direction="up" if effect > 0 else "down",
                    caveats_inherited=caveats,
                    context=contrast_id,
                    independence_key=chromatin_tf_lineage_key(tf),
                )
            )
            n += 1
            if n >= _MAX_TFS_PER_CONTRAST:
                break
    return out


def cross_modal_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited cross-modal concordance into EvidenceSignal items.

    Emits nothing unless the concordance engine stamped ``ran is True`` AND (when
    a ledger is supplied) ``ledger://cross_modal/concordance`` has status
    ``ran``. ``exp_ctx`` is accepted for adapter-signature parity only.
    """
    findings = cross_modal_findings(agent_results)
    if not findings or not _engine_ran(findings, run_ledger):
        return []
    return _sample_link_signals(findings) + _tf_signals(findings)


def tag_chromatin_tf_lineage(signals: list, cross_modal_signals: list) -> int:
    """Stamp the shared chromatin TF lineage key on single-assay ATAC TF signals.

    For every TF that has a ``cross_modal_tf_concordance`` signal, the bulk
    ATAC / scATAC motif-enrichment and footprint ``tf_motif`` signals naming
    that TF (a ``::`` dimer matches by subunit) receive the SAME
    ``independence_key`` — the derived concordance never counts as a third
    line. Signals that already carry a lineage key are left untouched (no
    destructive re-keying). Returns the number of signals tagged.
    """
    lineage_tfs = {
        str(s.entity).strip().lower()
        for s in cross_modal_signals or []
        if getattr(s, "measure", None) == TF_CONCORDANCE_MEASURE
    }
    if not lineage_tfs:
        return 0
    n = 0
    for sig in signals or []:
        if getattr(sig, "modality", None) not in (
            "bulk_ATAC",
            "scATAC",
            "ChIP",
            "CUT_AND_RUN",
            "CUT_AND_TAG",
        ):
            continue
        if getattr(sig, "entity_kind", None) != "tf_motif":
            continue
        if getattr(sig, "independence_key", ""):
            continue
        raw = str(getattr(sig, "entity", "") or "").split("(")[0]
        parts = sorted(
            p.strip().lower() for p in raw.split("::") if p.strip()
        )
        match = next((p for p in parts if p in lineage_tfs), None)
        if match is None:
            continue
        sig.independence_key = chromatin_tf_lineage_key(match)
        n += 1
    return n
