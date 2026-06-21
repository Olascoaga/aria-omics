"""scRNA evidence adapter (ADR-057 S6).

Pure and deterministic — no LLM. Flattens the AUDITED scRNA findings into
``EvidenceSignal`` items across the modality's three associative/inferential
layers:

  - pseudobulk donor-level DE  -> gene signals (the inferential primary).
  - differential abundance      -> cell-type signals (composition shifts).
  - LIANA cell-cell comms       -> ligand/receptor gene + source/target cell-type
                                    signals (transcript-supported candidates).

Every signal inherits the scRNA technical confounds the run did NOT address
(unaddressed ambient RNA, doublets, residual batch — read from the same
``devils_advocate._run_signals`` the report uses) so the devils_advocate gate
forces the LLM to own them. A layer is emitted only when its run-ledger node
actually ran; honest-empty otherwise. Raw per-cluster markers are intentionally
excluded (double-dipping / descriptive, not inferential).
"""

from __future__ import annotations

from typing import Any

from aria.agents.narrative.devils_advocate import _run_signals
from aria.agents.narrative.run_ledger import _node_index, _scrna_findings

from ..types import EvidenceSignal

_DE_NODE = "ledger://scRNA/pseudobulk_de"
_DA_NODE = "ledger://scRNA/differential_abundance"
_CCC_NODE = "ledger://scRNA/cell_communication"
_MAX_GENES_PER_COMPARISON = 15
_MAX_INTERACTIONS = 10

_UP_TOKENS = ("up", "incre", "higher", "gain", "enrich", "elevat", "expand")
_DOWN_TOKENS = ("down", "decre", "lower", "loss", "deplet", "reduc", "contract")


def _node_ran(run_ledger: dict | None, node_id: str) -> bool:
    if run_ledger is None:
        return True
    node = _node_index(run_ledger).get(node_id)
    return bool(node) and node.get("status") == "ran"


def _num(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _num_direction(value: Any) -> str:
    v = _num(value)
    if v is None:
        return "na"
    return "up" if v > 0 else "down" if v < 0 else "na"


def _str_direction(value: Any) -> str:
    v = _num(value)
    if v is not None:
        return "up" if v > 0 else "down" if v < 0 else "na"
    s = str(value or "").lower()
    if any(t in s for t in _UP_TOKENS):
        return "up"
    if any(t in s for t in _DOWN_TOKENS):
        return "down"
    return "na"


def _scrna_caveats(agent_results: dict, exp_ctx: dict | None) -> list[str]:
    """Global scRNA technical confounds the run did not address."""
    try:
        sig = _run_signals(agent_results, exp_ctx or {})
    except Exception:
        return []
    caveats: list[str] = []
    if not sig.get("ambient_corrected"):
        caveats.append("ambient RNA contamination not corrected")
    if not sig.get("doublets_handled"):
        caveats.append("doublets not ruled out")
    if not sig.get("integration_clean"):
        caveats.append("residual batch effect (integration not clean)")
    return caveats


def scrna_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited scRNA findings into EvidenceSignal items."""
    findings = _scrna_findings(agent_results)
    if not findings:
        return []

    base = _scrna_caveats(agent_results, exp_ctx)
    signals: list[EvidenceSignal] = []
    seen: set[tuple[str, str]] = set()

    def _add(entity, kind, measure, node, *, value=None, direction="na", extra=None):
        ent = str(entity or "").strip()
        if not ent or ent == "?":
            return
        key = (ent.lower(), node)
        if key in seen:
            return
        seen.add(key)
        caveats = list(base)
        if extra:
            caveats.append(extra)
        signals.append(
            EvidenceSignal(
                entity=ent,
                entity_kind=kind,
                modality="scRNA",
                measure=measure,
                audited_node_ref=node,
                value=value,
                direction=direction,
                caveats_inherited=caveats,
            )
        )

    # pseudobulk donor-level DE -> gene signals
    if _node_ran(run_ledger, _DE_NODE):
        pb = findings.get("pseudobulk_de") or {}
        for info in (pb.get("per_group") or {}).values():
            if not isinstance(info, dict):
                continue
            for comp in (info.get("per_comparison") or {}).values():
                if not isinstance(comp, dict) or comp.get("status") != "success":
                    continue
                n = 0
                for gene in comp.get("top_genes") or []:
                    if not isinstance(gene, dict):
                        continue
                    symbol = gene.get("gene") or gene.get("symbol")
                    lfc = gene.get("log2fc")
                    if not symbol or lfc is None:
                        continue
                    _add(symbol, "gene", "log2fc", _DE_NODE,
                         value=_num(lfc), direction=_num_direction(lfc))
                    n += 1
                    if n >= _MAX_GENES_PER_COMPARISON:
                        break

    # differential abundance -> cell-type signals
    if _node_ran(run_ledger, _DA_NODE):
        da = findings.get("differential_abundance") or {}
        for comp in (da.get("per_comparison") or {}).values():
            if not isinstance(comp, dict) or comp.get("status", "success") != "success":
                continue
            for row in comp.get("per_cell_type") or []:
                if not isinstance(row, dict) or not row.get("significant"):
                    continue
                _add(row.get("name"), "celltype", "abundance", _DA_NODE,
                     direction=_str_direction(row.get("direction")),
                     extra="cell-type composition shift")

    # LIANA cell-cell communication -> ligand/receptor genes + cell-type signals
    if _node_ran(run_ledger, _CCC_NODE):
        ccc = findings.get("cell_communication") or {}
        if ccc.get("status") in {"done", "success"}:
            extra = "ligand-receptor candidate (transcript-supported, needs review)"
            for row in (ccc.get("top_interactions") or [])[:_MAX_INTERACTIONS]:
                if not isinstance(row, dict):
                    continue
                pair = (
                    f"{row.get('source', '?')}->{row.get('target', '?')} "
                    f"({row.get('ligand', '?')}-{row.get('receptor', '?')})"
                )
                _add(row.get("ligand"), "gene", "cell_communication", _CCC_NODE,
                     value=pair, extra=extra)
                _add(row.get("receptor"), "gene", "cell_communication", _CCC_NODE,
                     value=pair, extra=extra)
                _add(row.get("source"), "celltype", "cell_communication", _CCC_NODE,
                     extra=extra)
                _add(row.get("target"), "celltype", "cell_communication", _CCC_NODE,
                     extra=extra)

    return signals
