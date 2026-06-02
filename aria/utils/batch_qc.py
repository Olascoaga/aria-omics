"""Hidden (unmodeled) batch red-flags for scRNA (audit item P1-4).

ARIA only runs Harmony when a batch column is declared in the confirmed design
(`scrna_agent._resolve_batch_column`). If the data carries a technical/batch
column that the user did NOT declare — a sequencing lane, library, 10x chip,
processing run, flowcell — nothing corrects it and nothing models it in DE, so
its variance can masquerade as biology. `integration_qc` only fires AFTER
integration ran, so it cannot see this case.

This module WARNS, it never corrects (P1-4 acceptance: "detector que advierte,
no corrige a ciegas"). It is name + design based — it inspects obs column NAMES
and the confirmed design (condition / replicate / declared batch), never
cell-level values. The batch-token list is generic technical vocabulary (an
ADR-011 technical-detection exception, like `aria/utils/sensitivity.py`), not
biological content.

Limitation (stated, not hidden): a technical factor whose column name does not
match a known batch token, or technical structure with no recorded column at
all, cannot be detected here. Declaring the column in the design (so Harmony
runs / DE adjusts for it) is the recommended remedy — ARIA does not guess and
correct on its own.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# Generic technical / batch column-name tokens (substring match, lowercased).
# Donor/subject/patient/individual are intentionally EXCLUDED: those are the
# biological replicate, which the pseudobulk path already models — flagging them
# as "hidden batch" would be wrong.
_BATCH_TOKENS = (
    "batch", "lane", "library", "flowcell", "flow_cell", "run_id", "run",
    "chip", "pool", "10x", "chemistry", "processing", "process_date",
    "seq_date", "sequencing", "capture", "channel",
)


@dataclass(frozen=True)
class QCIssue:
    severity: str           # "warning" | "blocking"
    check: str
    message: str
    recommendation: str

    def as_dict(self) -> dict[str, str]:
        return {
            "severity": self.severity,
            "check": self.check,
            "message": self.message,
            "recommendation": self.recommendation,
        }


def _looks_like_batch(col: str) -> bool:
    name = str(col).lower()
    return any(tok in name for tok in _BATCH_TOKENS)


def assess_hidden_batch(
    obs_columns: list[str] | None,
    *,
    condition_col: str | None = None,
    replicate_col: str | None = None,
    declared_batch: str | None = None,
    integration_ran: bool = False,
) -> dict[str, Any]:
    """Flag candidate technical/batch columns that are present but unmodeled.

    Args:
        obs_columns: obs column NAMES present in the data (no values).
        condition_col: the confirmed condition/main factor (modeled in DE).
        replicate_col: the confirmed replicate/donor column (modeled in
            pseudobulk); never treated as a hidden batch.
        declared_batch: the batch column declared in the design (corrected by
            integration and/or available as a DE covariate).
        integration_ran: whether batch integration (e.g. Harmony) actually ran.

    Returns: {"status": "clean"|"warnings", "issues": [...],
              "candidate_batch_columns": [...]}.
    """
    cols = [str(c) for c in (obs_columns or [])]
    issues: list[QCIssue] = []

    # Columns that are already accounted for and must not be re-flagged.
    modeled = {c for c in (condition_col, replicate_col, declared_batch) if c}
    # If integration ran on the declared batch, that column is corrected.
    declared_lower = str(declared_batch).lower() if declared_batch else None

    candidates: list[str] = []
    for col in cols:
        if col in modeled:
            continue
        if not _looks_like_batch(col):
            continue
        # A declared batch under a different casing is still modeled.
        if declared_lower and str(col).lower() == declared_lower:
            continue
        candidates.append(col)

    if candidates:
        joined = ", ".join(candidates)
        plural = "columns" if len(candidates) > 1 else "column"
        issues.append(QCIssue(
            "warning",
            "unmodeled_batch",
            (f"Potential technical/batch {plural} present but not modeled: "
             f"{joined}. No batch correction or DE covariate uses {'them' if len(candidates) > 1 else 'it'}, "
             f"so the associated technical variance is not accounted for and "
             f"could be mistaken for biological signal. ARIA does not correct "
             f"this automatically."),
            (f"Declare the relevant column in the experimental design (batch "
             f"factor) so integration corrects it and/or DE adjusts for it; or "
             f"confirm it is biologically irrelevant. Do not interpret "
             f"cluster/DE structure that aligns with {joined} as biology."),
        ))

    return {
        "status": "warnings" if issues else "clean",
        "issues": [i.as_dict() for i in issues],
        "candidate_batch_columns": candidates,
    }
