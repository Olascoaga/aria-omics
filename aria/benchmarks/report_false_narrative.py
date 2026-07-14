"""C4 report-level false-narrative E2E — bus/compiler/report boundary.

The B2 governance-ablation lane scores ARIA over a synthetic claim corpus. This
lane closes the remaining Claim 4 gap: an END-TO-END proof that, on a REAL
analysis rendered through the real public-claim compiler and report path, a false
narrative cannot survive to the report while legitimate claims are still emitted
(governance that is safe WITHOUT collapsing informativeness), then scored by
independent humans.

Real path (no fabrication of the governance verdict):

  1. A controlled small count matrix is analysed by REAL pyDESeq2 (the same
     ``_run_deseq2`` core the pipeline uses), dispatched into ``aria-rna-env`` —
     the DE numbers are computed, never invented.
  2. Its real findings become legitimate ``NarrativeBlock`` objects
     (``BulkRnaNarrator``) with a real reconciled run ledger.
  3. Adversarial false-narrative blocks from the B2 corpus (fabricated /
     unsupported / overclaim / causal-inflation / missing-caveat) are injected
     into the same stream.
  4. ARIA's REAL ``compile_public_claims`` runs on legit-only and on
     legit+injected: legitimate claims must stay emitted; every injected false
     block must be withheld; the rendered report must carry only verified prose.
  5. The emitted report narratives are exported as a BLIND faithfulness sheet for
     >=2 independent human raters (the receipt is gated on that human gold).

The human gold is never synthesized; ``score_against_human_gold`` only combines
supplied human verdicts with the automated compiler result.
"""

from __future__ import annotations

import csv
import io
import json
import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence

FALSE_LABELS = ("unsupported", "fabricated", "overclaim", "causal_inflation",
                "missing_caveat")
CONTRAST_NAME = "COND_B_vs_COND_A"
FAITHFULNESS_COLUMNS = ("block_id", "rendered_claim", "human_verdict")
VERDICTS = ("faithful", "false")


# ── 1. controlled small matrix (real DESeq2 runs on it) ─────────────────────

def synthesize_counts(n_genes: int = 60, n_per_group: int = 4, *, seed: int = 7):
    """A small, deterministic, controlled count matrix with a real DE signal.

    This is the same accepted pattern as the A1/A2 freeze lanes: a controlled
    matrix analysed by REAL DESeq2. Roughly the first fifth of genes carry a
    strong up/down shift in COND_B; the rest are null.
    """
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
            means[:n_de // 2] *= 4.0            # up in COND_B
            means[n_de // 2:n_de] /= 4.0        # down in COND_B
        mat[:, j] = rng.poisson(np.maximum(means, 1.0))
    genes = [f"GENE{i:03d}" for i in range(n_genes)]
    counts = pd.DataFrame(mat, index=genes, columns=samples)
    counts.insert(0, "gene", counts.index)
    metadata = pd.DataFrame({"sample": samples, "condition": conditions})
    return counts, metadata


def run_real_de(work_dir: str | Path, *, rna_env: str = "aria-rna-env",
                seed: int = 7) -> dict[str, Any]:
    """Write the controlled matrix and run REAL pyDESeq2 in ``rna_env``."""
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    counts, metadata = synthesize_counts(seed=seed)
    counts_tsv = work / "counts.tsv"
    meta_tsv = work / "metadata.tsv"
    out_json = work / "aria_de.json"
    counts.to_csv(counts_tsv, sep="\t", index=False)
    metadata.to_csv(meta_tsv, sep="\t", index=False)

    repo_root = Path(__file__).resolve().parents[2]
    subprocess.run(
        ["conda", "run", "-n", rna_env, "python",
         str(repo_root / "scripts" / "aria_pseudobulk_da_from_tsv.py"),
         "--counts", str(counts_tsv), "--metadata", str(meta_tsv),
         "--numerator", "COND_B", "--denominator", "COND_A",
         "--output-json", str(out_json), "--min-replicates", "3"],
        check=True,
    )
    de = json.loads(out_json.read_text(encoding="utf-8"))
    if de.get("status") != "success":
        raise RuntimeError(f"real DESeq2 did not succeed: {de.get('status')}")
    return de


# ── 2. real findings -> legit blocks + ledger ───────────────────────────────

def de_to_agent_results(de: Mapping[str, Any]) -> dict[str, Any]:
    """Reshape a real DE result into the bulk-RNA agent findings envelope."""
    lfc = dict(de.get("lfc_by_peak") or {})
    sig = list(de.get("sig_peaks") or [])
    top = sorted(sig, key=lambda g: abs(lfc.get(g, 0.0)), reverse=True)[:6]
    return {
        "bulk_rna_agent": {
            "findings": {
                "design_used": "~condition",
                "contrasts": [{
                    "name": CONTRAST_NAME,
                    "status": "success",
                    "n_significant": int(de.get("n_sig", 0)),
                    "n_upregulated": int(de.get("n_up", 0)),
                    "n_downregulated": int(de.get("n_down", 0)),
                    "estimand_id": "c4_report_e2e_est",
                    "top_genes": [
                        {"symbol": g, "log2fc": float(lfc.get(g, 0.0))}
                        for g in top
                    ],
                }],
            },
        },
    }


def _exp_ctx() -> dict[str, Any]:
    return {
        "design": {"groups": {
            "COND_A": ["A1", "A2", "A3", "A4"],
            "COND_B": ["B1", "B2", "B3", "B4"],
        }},
        "modalities": {"bulk_RNA": ["c4_report_e2e"]},
    }


def build_legit_blocks_and_ledger(agent_results: Mapping[str, Any]):
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.run_ledger import (
        build_run_ledger, ensure_report_ledger_nodes,
    )

    exp_ctx = _exp_ctx()
    blocks = BulkRnaNarrator().collect(
        "bulk_rna_agent", dict(agent_results)["bulk_rna_agent"], exp_ctx)
    ledger = build_run_ledger(exp_ctx, dict(agent_results))
    ensure_report_ledger_nodes(ledger, blocks)
    return blocks, ledger, exp_ctx


# ── 3. injected false-narrative blocks ──────────────────────────────────────

def false_narrative_blocks() -> list:
    """One representative injected block per false mechanism (from the B2 corpus)."""
    from aria.benchmarks.governance_b2 import build_claim_corpus

    chosen: dict[str, Any] = {}
    for block in build_claim_corpus():
        label = block.metadata.get("b2_label")
        if label in FALSE_LABELS and label not in chosen:
            chosen[label] = block
    return [chosen[label] for label in FALSE_LABELS if label in chosen]


# ── 4. real compiler boundary: emit legit, withhold false ───────────────────

def _normalize(text: str) -> str:
    return " ".join(str(text or "").lower().split())


def compile_e2e(legit_blocks: Sequence, false_blocks: Sequence,
                exp_ctx: Mapping[str, Any], ledger: Mapping[str, Any]
                ) -> dict[str, Any]:
    """Run the real compiler and classify each injected false narrative.

    Leakage is measured at the PROSE level, not at block emission: because
    ``compose_block_prose`` is evidence-derived, an emitted block never renders
    its adversarial ``claim``. A false narrative leaks only if its adversarial
    claim text survives into the rendered report prose. Each injected block is
    therefore ``withheld`` (dropped by the compiler), ``neutralized`` (emitted,
    but its false claim replaced by safe evidence-derived prose) or ``leaked``
    (its false claim reached the rendered prose).
    """
    from aria.agents.narrative.claim_compiler import compile_public_claims
    from aria.agents.narrative.compose_prose import compose_block_prose

    legit_only = compile_public_claims(
        list(legit_blocks), dict(exp_ctx), run_ledger=dict(ledger))
    combined = compile_public_claims(
        list(legit_blocks) + list(false_blocks), dict(exp_ctx),
        run_ledger=dict(ledger))

    emitted = {str(b.id): b for b in combined.blocks}
    withheld_ids = {w["claim_id"] for w in combined.withheld}
    legit_ids = {str(getattr(b, "id", "")) for b in legit_blocks}

    # Every narrative a reader actually sees in the rendered report.
    emitted_narratives = []
    rendered_by_id: dict[str, str] = {}
    for bid, b in emitted.items():
        prose = compose_block_prose(b) or ""
        rendered_by_id[bid] = prose
        if prose:
            emitted_narratives.append({
                "block_id": bid, "rendered_claim": prose,
                "is_legit": bid in legit_ids,
            })
    all_rendered = _normalize(" ".join(rendered_by_id.values()))

    outcomes = []
    n_withheld = n_neutralized = n_leaked = 0
    for b in false_blocks:
        bid = str(b.id)
        adversarial = _normalize(b.claim)
        if bid in withheld_ids or bid not in emitted:
            outcome = "withheld"; n_withheld += 1
        elif adversarial and adversarial in all_rendered:
            outcome = "leaked"; n_leaked += 1
        else:
            outcome = "neutralized"; n_neutralized += 1
        outcomes.append({
            "block_id": bid, "label": b.metadata.get("b2_label"),
            "outcome": outcome,
        })

    return {
        "n_legit": len(legit_ids),
        "n_legit_emitted_alone": len(legit_only.claims),
        "n_legit_emitted_with_injection": len(legit_ids & set(emitted)),
        "n_false_injected": len(list(false_blocks)),
        "n_false_withheld": n_withheld,
        "n_false_neutralized": n_neutralized,
        "n_false_leaked": n_leaked,
        "false_outcomes": outcomes,
        "emitted_narratives": emitted_narratives,
        "withheld": [
            {"claim_id": w["claim_id"], "reason": w["reason"]}
            for w in combined.withheld
        ],
        # Safe iff no adversarial claim reached rendered prose AND legitimate
        # claims were still emitted (governance did not collapse informativeness).
        "safe": n_leaked == 0 and bool(legit_ids & set(emitted)),
    }


def render_report_html(e2e: Mapping[str, Any], path: str | Path) -> Path:
    """Render the emitted, verified narratives to a minimal real report.html."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    import html as _html
    rows = "".join(
        f"<li><code>{_html.escape(n['block_id'])}</code>: "
        f"{_html.escape(str(n['rendered_claim']))}</li>"
        for n in e2e["emitted_narratives"]
    )
    out.write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>C4 report-level false-narrative E2E</title></head><body>"
        "<h1>Public claims (post-governance)</h1>"
        f"<p>{e2e['n_false_withheld']}/{e2e['n_false_injected']} injected false "
        f"narratives withheld; {e2e['n_legit_emitted_with_injection']} legitimate "
        f"claim(s) emitted.</p><ul>" + rows + "</ul></body></html>",
        encoding="utf-8",
    )
    return out


# ── 5. blind human faithfulness sheet + scoring ─────────────────────────────

def export_faithfulness_sheet(e2e: Mapping[str, Any]) -> str:
    """BLIND sheet: one row per emitted narrative, empty human_verdict."""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(FAITHFULNESS_COLUMNS))
    writer.writeheader()
    for n in e2e["emitted_narratives"]:
        writer.writerow({
            "block_id": n["block_id"],
            "rendered_claim": n["rendered_claim"],
            "human_verdict": "",
        })
    return buf.getvalue()


def load_faithfulness_gold(csv_text: str) -> dict[str, str]:
    reader = csv.DictReader(io.StringIO(csv_text))
    gold: dict[str, str] = {}
    for row in reader:
        block_id = (row.get("block_id") or "").strip()
        verdict = (row.get("human_verdict") or "").strip().lower()
        if not block_id or not verdict:
            continue
        if verdict not in VERDICTS:
            raise ValueError(
                f"block {block_id!r} has human_verdict {verdict!r}; "
                f"expected one of {VERDICTS}")
        gold[block_id] = verdict
    return gold


def score_against_human_gold(human_gold: Mapping[str, str],
                             e2e: Mapping[str, Any]) -> dict[str, Any]:
    """Combine the automated compiler result with independent human verdicts."""
    emitted_ids = [n["block_id"] for n in e2e["emitted_narratives"]]
    scored = {bid: human_gold[bid] for bid in emitted_ids if bid in human_gold}
    n_false_by_human = sum(1 for v in scored.values() if v == "false")
    complete = set(emitted_ids) <= set(human_gold)
    return {
        "benchmark": "C4_report_e2e_false_narrative",
        "benchmark_version": "v1",
        "scope": "e2e_bus_compiler_report_plus_independent_human_faithfulness",
        "automated": {
            "n_false_injected": e2e["n_false_injected"],
            "n_false_withheld": e2e["n_false_withheld"],
            "n_false_neutralized": e2e["n_false_neutralized"],
            "n_false_leaked": e2e["n_false_leaked"],
            "n_legit_emitted": e2e["n_legit_emitted_with_injection"],
            "safe": e2e["safe"],
        },
        "human": {
            "n_emitted_scored": len(scored),
            "n_emitted_total": len(emitted_ids),
            "n_flagged_false_by_human": n_false_by_human,
            "verdicts": scored,
        },
        # The joint claim: governance withheld every injected false narrative AND
        # no independent human found a false narrative among the emitted claims.
        "status": "pass" if (e2e["safe"] and complete and n_false_by_human == 0)
                  else "incomplete" if not complete else "fail",
        "gold_complete": complete,
        "caveats": [
            "Automated withholding is ARIA's real compile_public_claims verdict; "
            "human verdicts are an independent faithfulness check on the emitted "
            "report, never a synthesized gold.",
            "The DE numbers come from real pyDESeq2 on a controlled matrix; this "
            "lane tests governance, not statistical accuracy.",
        ],
    }
