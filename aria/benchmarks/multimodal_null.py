"""C5 multimodal label-permutation null — RNA + ATAC through the public compiler.

Claim 5 is that ARIA reports an honest null instead of inventing biology. The B4
lane proves this on a synthetic null corpus; this lane proves it END-TO-END and
MULTIMODALLY on real analyses: a controlled RNA count matrix and a controlled
scATAC pseudobulk peak matrix, each carrying a real condition signal, are run
through the shared real DESeq2 core under the TRUE labels (a positive control
that must detect the signal) and under many SEEDED label permutations (a genuine
null in which the signal is destroyed). Each run's real findings become a real
NarrativeBlock and pass through ARIA's real ``compile_public_claims``.

The headline metric per modality is the false-positive narrative rate: the
fraction of null permutations whose emitted public claim asserts a significant
result — target ~0, bounded by the FDR level. No fabrication: the block is built
from the real DESeq2 output; a permutation that finds nothing yields an honest
"no significant features" claim, never an invented one.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

RNA_CONTRAST = "COND_B_vs_COND_A"


# ── controlled RNA matrix + label permutation (pure, seedable) ──────────────

def synthesize_rna_counts(n_genes: int = 60, n_per_group: int = 6, *,
                          seed: int = 13):
    """Controlled gene x sample counts with a real DE signal, analysed for real."""
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    samples = ([f"A{i + 1}" for i in range(n_per_group)]
               + [f"B{i + 1}" for i in range(n_per_group)])
    conditions = (["COND_A"] * n_per_group) + (["COND_B"] * n_per_group)
    n_de = max(6, n_genes // 5)
    base = rng.integers(80, 400, size=n_genes).astype(float)
    mat = np.zeros((n_genes, len(samples)), dtype=int)
    for j, cond in enumerate(conditions):
        means = base.copy()
        if cond == "COND_B":
            means[:n_de // 2] *= 4.0
            means[n_de // 2:n_de] /= 4.0
        mat[:, j] = rng.poisson(np.maximum(means, 1.0))
    counts = pd.DataFrame(mat, index=[f"GENE{i:03d}" for i in range(n_genes)],
                          columns=samples)
    counts.insert(0, "gene", counts.index)
    metadata = pd.DataFrame({"sample": samples, "condition": conditions})
    return counts, metadata


def _grouping(samples, labels) -> frozenset:
    """The COND_B sample set — the identity of a two-level grouping."""
    return frozenset(s for s, c in zip(samples, labels) if c == "COND_B")


def permute_conditions(metadata, *, seed: int, max_tries: int = 64):
    """Return metadata with the condition column shuffled across samples.

    A genuine label-permutation null must BREAK the true grouping. The true
    labeling and its exact complement (swapping COND_A<->COND_B yields the same
    grouping with a sign-symmetric |LFC|) both recover the real signal, so they
    are not null draws; they are rejected and re-drawn. Every returned grouping
    therefore differs from the observed one and from its complement.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    samples = list(metadata["sample"])
    true_labels = list(metadata["condition"])
    true_group = _grouping(samples, true_labels)
    complement = frozenset(samples) - true_group
    for _ in range(max_tries):
        shuffled = list(rng.permutation(true_labels))
        group = _grouping(samples, shuffled)
        if group != true_group and group != complement:
            permuted = metadata.copy()
            permuted["condition"] = shuffled
            return permuted
    raise RuntimeError(
        "could not draw a null permutation that breaks the true grouping; "
        "increase the replicate count")


# ── real findings -> real NarrativeBlock, per modality ──────────────────────

def rna_block_from_de(de: Mapping[str, Any]):
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    lfc = dict(de.get("lfc_by_peak") or {})
    sig = list(de.get("sig_peaks") or [])
    top = sorted(sig, key=lambda g: abs(lfc.get(g, 0.0)), reverse=True)[:6]
    findings = {"design_used": "~condition", "contrasts": [{
        "name": RNA_CONTRAST, "status": "success",
        "n_significant": int(de.get("n_sig", 0)),
        "n_upregulated": int(de.get("n_up", 0)),
        "n_downregulated": int(de.get("n_down", 0)),
        "estimand_id": "c5_rna_est",
        "top_genes": [{"symbol": g, "log2fc": float(lfc.get(g, 0.0))} for g in top],
    }]}
    agent_result = {"findings": findings}
    exp_ctx = {"design": {"groups": {"COND_A": ["A1"], "COND_B": ["B1"]}},
               "modalities": {"bulk_RNA": ["c5"]}}
    blocks = BulkRnaNarrator().collect("bulk_rna_agent", agent_result, exp_ctx)
    de_blocks = [b for b in blocks if b.analysis == "differential_expression"]
    return de_blocks[0], _agent_results_bulk(findings), exp_ctx


def _agent_results_bulk(findings):
    return {"bulk_rna_agent": {"findings": findings}}


def atac_block_from_de(de: Mapping[str, Any]):
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    n_sig = int(de.get("n_sig", 0))
    da = {
        "data_type": "scATAC",
        "padj_max": 0.05, "lfc_min": 0.5,
        "pseudobulk": {"ran": True, "comparisons": [{
            "test": "COND_B", "reference": "COND_A", "status": "success",
            "n_sig": n_sig, "n_up": int(de.get("n_up", 0)),
            "n_down": int(de.get("n_down", 0)),
        }]},
    }
    blocks = ChromatinNarrator()._diffacc_blocks(da)
    pb_blocks = [b for b in blocks
                 if b.id.endswith("pseudobulk") or "pseudobulk" in b.title.lower()]
    block = (pb_blocks or blocks)[0]
    agent_results = {"chromatin_agent": {"findings": {"differential_accessibility": da}}}
    exp_ctx = {"design": {"groups": {"COND_A": ["r1"], "COND_B": ["r2"]}},
               "modalities": {"scATAC": ["c5"]}}
    return block, agent_results, exp_ctx


# ── real public compiler: did an emitted claim assert significance? ─────────

def emits_significant_claim(block, agent_results: Mapping[str, Any],
                            exp_ctx: Mapping[str, Any], n_sig: int) -> dict[str, Any]:
    """Compile through the real boundary and classify the emitted narrative."""
    from aria.agents.narrative.claim_compiler import compile_public_claims
    from aria.agents.narrative.run_ledger import (
        build_run_ledger, ensure_report_ledger_nodes,
    )

    ledger = build_run_ledger(dict(exp_ctx), dict(agent_results))
    ensure_report_ledger_nodes(ledger, [block])
    comp = compile_public_claims([block], dict(exp_ctx), run_ledger=dict(ledger))
    emitted = any(str(b.id) == str(block.id) for b in comp.blocks)
    # An emitted claim asserts significance only when the real findings had a
    # non-empty significant set; an honest "0 features" claim is not a false
    # positive.
    asserts_significant = emitted and n_sig > 0
    return {"emitted": emitted, "n_sig": n_sig,
            "asserts_significant": asserts_significant,
            "withheld": [w["claim_id"] for w in comp.withheld]}


def classify_run(modality: str, de: Mapping[str, Any],
                 *, is_permuted: bool) -> dict[str, Any]:
    n_sig = int(de.get("n_sig", 0))
    if modality == "rna":
        block, agent_results, exp_ctx = rna_block_from_de(de)
    elif modality == "atac":
        block, agent_results, exp_ctx = atac_block_from_de(de)
    else:
        raise ValueError(f"unknown modality {modality!r}")
    verdict = emits_significant_claim(block, agent_results, exp_ctx, n_sig)
    return {
        "modality": modality, "is_permuted": is_permuted, "n_sig": n_sig,
        "emitted": verdict["emitted"],
        "asserts_significant": verdict["asserts_significant"],
        "false_positive": bool(is_permuted and verdict["asserts_significant"]),
    }


# ── scoring ─────────────────────────────────────────────────────────────────

def score_modality(modality: str, true_run: Mapping[str, Any],
                   perm_runs: Sequence[Mapping[str, Any]],
                   *, fdr_tolerance: float = 0.1) -> dict[str, Any]:
    n_perms = len(perm_runs)
    n_fp = sum(1 for r in perm_runs if r["false_positive"])
    fp_rate = n_fp / max(n_perms, 1)
    return {
        "modality": modality,
        "positive_control_detects_signal": bool(true_run["asserts_significant"]),
        "n_permutations": n_perms,
        "n_false_positive_narratives": n_fp,
        "false_positive_narrative_rate": round(fp_rate, 4),
        "max_null_n_sig": max((r["n_sig"] for r in perm_runs), default=0),
        "axis_pass": {
            "positive_control": bool(true_run["asserts_significant"]),
            "null_false_positive_bounded": fp_rate <= fdr_tolerance,
        },
    }


def score_multimodal_null(modality_scores: Mapping[str, Mapping[str, Any]]
                          ) -> dict[str, Any]:
    axes = [m["axis_pass"] for m in modality_scores.values()]
    status = "pass" if all(
        all(a.values()) for a in axes) and modality_scores else "fail"
    return {
        "benchmark": "C5_multimodal_null_permutations",
        "benchmark_version": "v1",
        "scope": "multimodal_label_permutation_null_through_public_compiler",
        "modalities": dict(modality_scores),
        "status": status,
        "caveats": [
            "Controlled matrices carry a real condition signal (positive "
            "control); label permutation destroys it (the null). DESeq2 and the "
            "public-claim compiler are ARIA's real code — the null verdict is "
            "computed, not asserted.",
            "The false-positive narrative rate is bounded by the FDR level; the "
            "lane tests governance + FDR control, not effect-size accuracy.",
        ],
    }
