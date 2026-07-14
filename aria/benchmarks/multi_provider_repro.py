"""B3-multi harness: reproducibility across LLM providers and repetitions.

Claim 6 (preprint audit, FASE 7): ARIA's public output is invariant under
allowed LLM prose variation. Given ONE structured payload (typed narrative
blocks + agent findings + provenance), running it through the real claim
compiler and run-ledger linkage under N provider prose styles x M repetitions
must yield an identical public claim set and an identical methodology diff. Only
free narrative prose may differ between cells.

This is the Fase A, in-process lane. It does not call any provider: a
:class:`ProseVariant` models *how* a given provider would word the same result,
and the harness proves that this wording cannot move the public boundary. The
compiler (`compile_public_claims`) tiers every claim from structured evidence,
never from prose, and `diff_methodologies` compares provenance, ledger nodes,
claim tiers/linkage and calibration — not prose text. The heavy end-to-end lane
over ``run_headless`` is a separate opt-in guard (Fase B).

The harness deliberately reuses production components (`compile_public_claims`,
`build_run_ledger`, `ensure_report_ledger_nodes`, `link_claims_to_ledger`,
`diff_methodologies`) so a regression that makes any public field prose-derived
is caught here rather than hidden behind a bespoke reimplementation.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Callable

from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.ledger_export import diff_methodologies
from aria.agents.narrative.run_ledger import (
    build_run_ledger,
    ensure_report_ledger_nodes,
    link_claims_to_ledger,
)
from aria.agents.narrative.types import NarrativeBlock


def _identity(text: str) -> str:
    return text


@dataclass(frozen=True)
class ProseVariant:
    """One provider's stylistic rendering of the SAME structured evidence.

    ``claim_style`` rewrites a block's free ``claim`` prose (a paraphrase, or an
    adversarial causal/speculative overreach). It must never touch structured
    evidence, metrics, facts or estimands — those are the analysis output, not
    prose. ``free_narrative`` is optional operational prose a provider might emit
    alongside the report; it is never a public-claim source and is carried only
    so a caller can record what varied.
    """

    provider: str
    claim_style: Callable[[str], str] = _identity
    free_narrative: str = ""


@dataclass
class CellResult:
    """The compiled outcome for one (provider, repetition) matrix cell."""

    provider: str
    repetition: int
    methodology: dict[str, Any]
    published_claim_ids: tuple[str, ...]
    withheld_claim_ids: tuple[str, ...]
    free_narrative: str = ""


@dataclass
class MatrixResult:
    """Every matrix cell plus the invariance verdict against the baseline."""

    cells: list[CellResult] = field(default_factory=list)

    @property
    def baseline(self) -> CellResult:
        if not self.cells:
            raise ValueError("matrix produced no cells")
        return self.cells[0]

    def pairwise_diffs(self) -> list[dict[str, Any]]:
        """Diff every non-baseline cell's methodology against the baseline."""
        base = self.baseline.methodology
        return [
            diff_methodologies(base, cell.methodology)
            for cell in self.cells[1:]
        ]

    def invariant(self) -> bool:
        """True iff every cell is identical to the baseline over tracked fields."""
        return all(diff.get("identical") for diff in self.pairwise_diffs())

    def summary(self) -> dict[str, Any]:
        diffs = self.pairwise_diffs()
        return {
            "harness": "multi_provider_repro",
            "n_cells": len(self.cells),
            "providers": sorted({c.provider for c in self.cells}),
            "repetitions": sorted({c.repetition for c in self.cells}),
            "invariant": all(d.get("identical") for d in diffs),
            "n_divergent_cells": sum(
                0 if d.get("identical") else 1 for d in diffs
            ),
            "baseline_published": list(self.baseline.published_claim_ids),
        }


def build_cell_methodology(
    *,
    blocks: list[NarrativeBlock],
    exp_ctx: dict[str, Any],
    agent_results: dict[str, Any],
    provenance: dict[str, Any],
) -> tuple[dict[str, Any], Any]:
    """Compile one cell's public claims and assemble a diff-able methodology.

    Mirrors the diff-relevant subset of
    ``ReportBuilderMixin._build_methodology_json``: a typed run ledger, the
    public claim compilation, fixed seeds and the input hashes. Returns
    ``(methodology, compilation)``.
    """
    run_ledger = build_run_ledger(exp_ctx, agent_results)
    ensure_report_ledger_nodes(run_ledger, blocks)
    compilation = compile_public_claims(blocks, exp_ctx, run_ledger=run_ledger)
    run_ledger["claim_linkage"] = link_claims_to_ledger(
        compilation.claims, run_ledger
    )
    methodology = {
        "provenance": provenance,
        "inputs": exp_ctx.get("input_files", []),
        "claims": compilation.claims,
        "run_ledger": run_ledger,
        "seeds": {"global": 0, "scanpy": 0, "harmony": 0},
    }
    return methodology, compilation


def _apply_variant(
    blocks: list[NarrativeBlock], variant: ProseVariant
) -> list[NarrativeBlock]:
    """Reword each block's free ``claim`` prose; never touch structured fields."""
    styled = []
    for block in blocks:
        clone = copy.deepcopy(block)
        clone.claim = variant.claim_style(clone.claim)
        styled.append(clone)
    return styled


def _cell(
    *,
    provider: str,
    repetition: int,
    blocks: list[NarrativeBlock],
    exp_ctx: dict[str, Any],
    agent_results: dict[str, Any],
    provenance: dict[str, Any],
    free_narrative: str = "",
) -> CellResult:
    methodology, compilation = build_cell_methodology(
        blocks=blocks,
        exp_ctx=copy.deepcopy(exp_ctx),
        agent_results=copy.deepcopy(agent_results),
        provenance=copy.deepcopy(provenance),
    )
    return CellResult(
        provider=provider,
        repetition=repetition,
        methodology=methodology,
        published_claim_ids=tuple(
            str(c.get("claim_id")) for c in compilation.claims
        ),
        withheld_claim_ids=tuple(
            str(w.get("claim_id")) for w in compilation.withheld
        ),
        free_narrative=free_narrative,
    )


def run_provider_matrix(
    *,
    block_factory: Callable[[], list[NarrativeBlock]],
    exp_ctx: dict[str, Any],
    agent_results: dict[str, Any],
    provenance: dict[str, Any],
    variants: list[ProseVariant],
    repetitions: int = 1,
    extra_cells: list[tuple[str, Callable[[], list[NarrativeBlock]]]] | None = None,
) -> MatrixResult:
    """Run the provider x repetition matrix over one structured payload.

    ``block_factory`` returns FRESH, structurally identical blocks per cell
    (the compiler mutates ``block.metadata`` in place, so cells must not share
    state). Each :class:`ProseVariant` rewords the free claim prose; the
    structured evidence is untouched. ``extra_cells`` injects labelled cells with
    a different structured payload — used by tests as a negative control to prove
    the invariance check is sensitive to real structural divergence.
    """
    if repetitions < 1:
        raise ValueError("repetitions must be >= 1")
    if not variants:
        raise ValueError("at least one ProseVariant is required")

    cells: list[CellResult] = []
    for variant in variants:
        for rep in range(repetitions):
            cells.append(
                _cell(
                    provider=variant.provider,
                    repetition=rep,
                    blocks=_apply_variant(block_factory(), variant),
                    exp_ctx=exp_ctx,
                    agent_results=agent_results,
                    provenance=provenance,
                    free_narrative=variant.free_narrative,
                )
            )

    for label, factory in extra_cells or []:
        cells.append(
            _cell(
                provider=label,
                repetition=0,
                blocks=factory(),
                exp_ctx=exp_ctx,
                agent_results=agent_results,
                provenance=provenance,
            )
        )

    return MatrixResult(cells=cells)
