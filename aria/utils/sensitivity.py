"""Input sensitivity classification (P1-8a, W-PRIV).

Before any analysis, ARIA assesses whether the input looks human and/or carries
clinical / PHI-like signals, and surfaces that assessment at the first
checkpoint so the user can decide whether to enable air-gapped mode (which blocks
ALL network egress: cloud LLM, Enrichr ORA, GEO/SRA connectors).

ARIA never auto-disables egress on its own (per the adopted policy): it
classifies, recommends, and lets the user choose at the checkpoint. The detection
uses generic *technical* tokens (column/path naming conventions), not dataset- or
biology-specific content, consistent with ADR-011's technical-detection
exceptions (e.g. organism inference).

This module is dependency-free and pure so it is trivially testable.
"""

from __future__ import annotations

import re

# Clinical / PHI-like field or path naming conventions. Generic identifiers a
# clinical dataset tends to carry — NOT biological content.
_PHI_TOKENS = (
    "patient", "subject", "mrn", "medicalrecord", "diagnosis", "icd",
    "dateofbirth", "birthdate", "dob", "ssn", "clinical", "phi", "hospital",
    "biopsy", "specimen", "consent", "accession", "ethnicity", "zipcode",
)
# Quasi-identifiers: meaningful only in combination (age alone is common and not
# sensitive; age + sex/gender starts to re-identify).
_QUASI_AGE = ("age", "ageyears", "ageatdiagnosis")
_QUASI_SEX = ("sex", "gender")

_HUMAN_TOKENS = ("human", "sapiens", "homosapiens", "hsapiens", "hg38", "hg19", "grch")


def _norm(text) -> str:
    """Lowercase and strip non-alphanumerics so `patient_id`/`patientID` match."""
    return re.sub(r"[^a-z0-9]", "", str(text).lower())


def classify_sensitivity(organism: str = "",
                         field_names=None,
                         path_hints=None) -> dict:
    """Classify input sensitivity from organism + field/path naming.

    Returns a JSON-serializable assessment:
      {level, is_human, phi_signals, quasi_identifier_signals,
       recommend_air_gapped, summary}

    `level` is one of `low` / `elevated` / `high`:
      * high     — clinical/PHI-like signals present (recommend air-gapped);
      * elevated — human data without PHI signals (privacy worth a decision);
      * low      — neither.
    """
    fields = [str(f) for f in (field_names or [])]
    hints = [str(h) for h in (path_hints or [])]
    normed = [_norm(x) for x in (fields + hints) if x]

    is_human = any(tok in _norm(organism) for tok in _HUMAN_TOKENS) or \
        any(any(tok in n for tok in _HUMAN_TOKENS) for n in normed)

    phi_signals = sorted({
        tok for tok in _PHI_TOKENS
        for n in normed if tok in n
    })
    has_age = any(any(tok == n or tok in n for tok in _QUASI_AGE) for n in normed)
    has_sex = any(any(tok == n or tok in n for tok in _QUASI_SEX) for n in normed)
    quasi = []
    if has_age and has_sex:
        quasi = ["age+sex"]

    if phi_signals:
        level = "high"
    elif is_human and quasi:
        level = "high"
    elif is_human:
        level = "elevated"
    else:
        level = "low"

    recommend = level == "high"

    if level == "high":
        reasons = []
        if phi_signals:
            reasons.append("clinical/PHI-like fields (" + ", ".join(phi_signals) + ")")
        if quasi:
            reasons.append("quasi-identifiers (" + ", ".join(quasi) + ")")
        summary = (
            "Input looks sensitive: " + "; ".join(reasons) + ". "
            "Air-gapped mode is recommended to block ALL network egress "
            "(cloud LLM, Enrichr ORA, GEO/SRA)."
        )
    elif level == "elevated":
        summary = (
            "Input appears to be human data with no obvious PHI fields. "
            "Consider air-gapped mode if the data are not de-identified."
        )
    else:
        summary = "No human or clinical/PHI signals detected in the input metadata."

    return {
        "level": level,
        "is_human": is_human,
        "phi_signals": phi_signals,
        "quasi_identifier_signals": quasi,
        "recommend_air_gapped": recommend,
        "summary": summary,
    }


# ── Checkpoint contract (kept here so the wiring is unit-testable) ────────────

AIR_GAPPED_OPTION = "Confirm and enable air-gapped mode"


def checkpoint_options(sensitivity: dict) -> list[str]:
    """CP1 options. 'Confirm and continue' stays FIRST so the default (incl.
    headless) is never silently changed; the air-gapped option is always offered
    and recommended in the question text when the input looks sensitive."""
    return [
        "Confirm and continue",
        AIR_GAPPED_OPTION,
        "Correct metadata",
        "Cancel",
    ]


def annotate_checkpoint_question(question: str, sensitivity: dict) -> str:
    """Append the sensitivity assessment + air-gapped guidance to the CP1 text."""
    rec = " (RECOMMENDED)" if sensitivity.get("recommend_air_gapped") else ""
    return question + (
        f"\n\n— Data sensitivity ({sensitivity.get('level', 'low').upper()}): "
        f"{sensitivity.get('summary', '')}"
        f"\n  To block ALL network egress (cloud LLM, Enrichr ORA, GEO/SRA), "
        f"choose '{AIR_GAPPED_OPTION}'{rec}."
    )


def decision_enables_air_gapped(decision: str) -> bool:
    """True when a CP1 decision string selects the air-gapped option."""
    return "air-gap" in str(decision).lower()
