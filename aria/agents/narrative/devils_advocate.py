"""Deterministic devil's advocate on the validated path (R2 / P-DEVIL,
audit 2026-05-29).

DebateCouncil's "alternative hypothesis first" reasoning only ran on the
scaffold modalities, so the science ARIA actually ships (scRNA / bulk RNA) had
no adversarial check. Rather than pay for 3 LLM rounds, this is a single-shot,
deterministic challenge: for every associative-or-stronger claim, enumerate the
standard TECHNICAL alternative explanations (batch, ambient RNA, doublets,
composition shift, low replication) and state, from the run's own structured
evidence, whether each was addressed.

The confounder list is methodological/technical, not biological content — the
same category as the causal-language lexicon (ADR-011 does not apply to fixed
QC/technical vocabulary). It reuses the Claim Compiler tiers and adds an
``info`` caveat naming the unaddressed alternatives, plus a manifest for
``methodology.json``.
"""

from __future__ import annotations

import json
from typing import Any

from aria.agents.narrative.types import Caveat, NarrativeBlock


_CHALLENGEABLE_TIERS = {"associative", "weak_mechanistic", "strong_mechanistic"}


def _scrna_findings(agent_results: dict) -> dict:
    sc = (agent_results or {}).get("scrna_agent", {})
    try:
        from aria.agents import _narrative_scrna
        return _narrative_scrna.unwrap_scrna_findings(sc) or {}
    except Exception:
        if isinstance(sc, dict):
            return sc.get("findings", sc) or {}
        return {}


def _ran(val: Any) -> bool:
    if isinstance(val, dict):
        return str(val.get("status") or "").lower() not in ("skipped", "error") \
            and bool(val)
    return bool(val)


def _contains(blob: Any, *needles: str) -> bool:
    try:
        text = json.dumps(blob, default=str).lower()
    except Exception:
        text = str(blob).lower()
    return any(n in text for n in needles)


def _run_signals(agent_results: dict, exp_ctx: dict) -> dict:
    f = _scrna_findings(agent_results)
    integration_ran = _ran(f.get("integration"))
    iqc = f.get("integration_qc") or {}
    integration_issue = bool(iqc.get("issues") or iqc.get("flags")
                             or iqc.get("residual_batch"))
    da = f.get("differential_abundance") or {}
    covariates = [str(c).lower()
                  for c in ((exp_ctx or {}).get("design", {}) or {}).get(
                      "covariates", []) or []]
    return {
        "integration_clean": integration_ran and not integration_issue,
        "integration_ran": integration_ran,
        "doublets_handled": _contains(f.get("qc") or {}, "doublet", "scrublet"),
        "ambient_corrected": _contains(f, "soupx", "decontx", "ambient"),
        "da_ran": _ran(da),
        "da_significant": bool(da.get("any_significant")),
        "batch_covariate": any("batch" in c for c in covariates),
    }


def _block_low_power(block: NarrativeBlock) -> bool:
    if block.metrics.get("low_power_warning") or block.metrics.get("low_power"):
        return True
    return _contains([c.text for c in (block.caveats or [])],
                     "low replicate", "low power")


def _block_composition_ok(block: NarrativeBlock, signals: dict) -> bool:
    # Either the cell-type proportion did not shift, or the block modeled it.
    if block.metrics.get("corrected_for_composition"):
        return True
    if signals["da_ran"] and not signals["da_significant"]:
        return True
    return False


def _challenges_for(block: NarrativeBlock, signals: dict) -> list[dict]:
    modality = str(block.modality or "").lower()
    is_scrna = "scrna" in modality or modality in ("rna", "")
    out: list[dict] = []

    # Batch effect (both modalities).
    if is_scrna:
        addressed = signals["integration_clean"]
        basis = ("batch correction (integration) ran and residual batch was "
                 "within tolerance" if addressed
                 else "no batch correction ran, or residual batch was flagged")
    else:
        addressed = signals["batch_covariate"]
        basis = ("a batch covariate was included in the design" if addressed
                 else "no batch covariate was modeled in the design")
    out.append({"alternative": "batch_effect", "addressed": addressed,
                "basis": basis})

    if is_scrna:
        out.append({
            "alternative": "ambient_rna",
            "addressed": signals["ambient_corrected"],
            "basis": ("ambient-RNA correction was applied"
                      if signals["ambient_corrected"] else
                      "no ambient-RNA correction (SoupX/decontX) is in the "
                      "pipeline — contamination is not ruled out"),
        })
        out.append({
            "alternative": "doublets",
            "addressed": signals["doublets_handled"],
            "basis": ("doublet detection ran in QC"
                      if signals["doublets_handled"] else
                      "no doublet detection signal found in QC"),
        })
        comp_ok = _block_composition_ok(block, signals)
        out.append({
            "alternative": "composition_shift",
            "addressed": comp_ok,
            "basis": ("composition was modeled or the cell-type proportion did "
                      "not shift" if comp_ok else
                      "cell-type composition shift was not ruled out for this "
                      "block"),
        })

    low_power = _block_low_power(block)
    out.append({
        "alternative": "low_replication",
        "addressed": not low_power,
        "basis": ("replicate support met the floor" if not low_power else
                  "low replicate support — effect-size/FDR estimates are noisy"),
    })
    return out


def build_devils_advocate(blocks: list[NarrativeBlock],
                          agent_results: dict,
                          exp_ctx: dict | None = None) -> list[dict]:
    """Annotate each challengeable block in place with a devil's-advocate caveat
    and ``metadata['devils_advocate']``, and return a manifest for
    ``methodology.json``. Requires claim tiers to be annotated first
    (``claim_compiler.annotate_claim_tiers``)."""
    signals = _run_signals(agent_results, exp_ctx or {})
    manifest: list[dict] = []

    for block in blocks:
        if block.status != "success":
            continue
        tier = ((block.metadata.get("claim") or {}).get("tier"))
        if tier not in _CHALLENGEABLE_TIERS:
            continue
        # Idempotent: this pass may run once before HTML render and once during
        # methodology build. Don't re-append the caveat on the second call.
        existing = block.metadata.get("devils_advocate")
        if existing is not None:
            manifest.append({
                "claim_id": block.id, "tier": tier,
                "challenges": existing.get("challenges", []),
                "unaddressed": existing.get("unaddressed", []),
            })
            continue
        challenges = _challenges_for(block, signals)
        unaddressed = [c["alternative"] for c in challenges if not c["addressed"]]
        block.metadata["devils_advocate"] = {
            "challenges": challenges,
            "unaddressed": unaddressed,
        }
        if unaddressed:
            block.caveats.append(Caveat(
                text=("Alternative technical explanations not ruled out: "
                      + ", ".join(unaddressed)
                      + ". Interpret the association with these confounders in "
                        "mind."),
                severity="info",
            ))
        manifest.append({
            "claim_id": block.id,
            "tier": tier,
            "challenges": challenges,
            "unaddressed": unaddressed,
        })
    return manifest
