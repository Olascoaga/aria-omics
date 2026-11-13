"""3D Genome Architecture (`HiC` / Micro-C) evidence adapter (ADR-082; ADR-057).

Pure and deterministic — zero LLM, zero biological priors. Flattens audited
``GenomeArchAgent`` chromatin-loop anchor annotations into ``EvidenceSignal``
items strictly downstream of ``W-CLAIM`` and ``W-LEDGER``.

Only named gene symbols anchored to statistically significant focal contact dots
(``ledger://hic/loops``) are emitted here; raw genomic bin coordinates from
compartments or TAD boundaries are not groundable biological entity names and
therefore never enter the ``SPECULATIVE`` entity pool. Every emitted signal
inherits ``HIC_LOOP_ASSOCIATIVE`` (physical 3D proximity does not establish
direct regulatory targeting or transcriptional activation without perturbation
validation) plus ``LOW_REPLICATION`` when fewer than two replicates were
profiled.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
from typing import Any

from aria.agents.narrative.narrators.hic import unwrap_hic_findings

from ..caveats import HIC_LOOP_ASSOCIATIVE, LOW_REPLICATION
from ..types import EvidenceSignal
from ._chromatin_common import node_ran

HIC_LOOPS_NODE = "ledger://hic/loops"
HIC_LOOP_MEASURE = "hic_loop_anchor"
_MAX_LOOP_GENES = 25

_COORD_RE = re.compile(r"^(?:chr[\w.\-]+|[\dXYM]+|MT):\d+[-:]\d+$", re.IGNORECASE)
_RAW_ID_PREFIXES = ("ENSG0", "ENSMUSG0", "ENST0", "ENSMUST0")
_SKIP_TOKENS = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "intergenic",
    "distal",
    "unannotated",
    "unspecified",
    ".",
    "-",
}


def _finite_float(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    return round(f, 6) if math.isfinite(f) else None


def _split_gene_field(raw: Any) -> list[str]:
    if raw is None:
        return []
    if isinstance(raw, (list, tuple, set)):
        tokens: list[str] = []
        for item in raw:
            tokens.extend(_split_gene_field(item))
        return tokens
    out: list[str] = []
    for part in re.split(r"[,;|/]+", str(raw)):
        gene = part.strip()
        if not gene or gene.lower() in _SKIP_TOKENS:
            continue
        if _COORD_RE.match(gene):
            continue
        if gene.upper().startswith(_RAW_ID_PREFIXES):
            continue
        if any(ch.isspace() for ch in gene):
            continue
        out.append(gene)
    return out


def _low_replication(findings: dict[str, Any]) -> bool:
    qc = findings.get("qc") or findings.get("balancing") or {}
    if isinstance(qc, dict):
        n_samp = qc.get("n_samples")
        try:
            if n_samp is not None and int(n_samp) >= 2:
                return False
            if n_samp is not None and int(n_samp) < 2:
                return True
        except (TypeError, ValueError):
            pass
    topo = findings.get("topology") if isinstance(findings.get("topology"), dict) else findings
    diff3d = (topo.get("differential_3d") if isinstance(topo, dict) else None) or findings.get("differential_3d") or {}
    if isinstance(diff3d, dict) and diff3d.get("descriptive_only") is True:
        return True
    samples = findings.get("samples")
    if isinstance(samples, list):
        return len(samples) < 2
    return True


def _loop_context(row: dict[str, Any]) -> str:
    c1 = row.get("chrom1")
    s1 = row.get("start1")
    e1 = row.get("end1")
    c2 = row.get("chrom2")
    s2 = row.get("start2")
    e2 = row.get("end2")
    if c1 is not None and s1 is not None and e1 is not None and c2 is not None and s2 is not None and e2 is not None:
        return f"{c1}:{s1}-{e1}<->{c2}:{s2}-{e2}"
    loop_id = row.get("loop_id") or row.get("id")
    return str(loop_id) if loop_id else "hic_loop"


def _candidate_loop_rows(loops: dict[str, Any]) -> list[dict[str, Any]]:
    for key in ("top_loops", "top_annotated_loops", "annotated_loops"):
        rows = loops.get(key)
        if isinstance(rows, list) and rows:
            valid = [r for r in rows if isinstance(r, dict)]
            if valid:
                return valid

    ann_path = loops.get("annotated_loops_path")
    if isinstance(ann_path, str) and ann_path:
        p = Path(ann_path)
        if p.is_file():
            out: list[dict[str, Any]] = []
            try:
                with p.open("r", encoding="utf-8", errors="replace") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for idx, r in enumerate(reader):
                        if idx >= 100:
                            break
                        out.append(dict(r))
            except OSError:
                return []
            return out
    return []


def hic_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited Hi-C loop-anchor genes into ``EvidenceSignal`` items."""
    del exp_ctx
    raw = (agent_results or {}).get("genome_arch_agent")
    findings = unwrap_hic_findings(raw)
    if not findings:
        return []

    topo = findings.get("topology") if isinstance(findings.get("topology"), dict) else findings
    loops = (topo.get("loops") if isinstance(topo, dict) else None) or findings.get("loops")
    if not isinstance(loops, dict):
        return []
    status = str(loops.get("status") or "").lower().strip()
    if status in {"error", "failed", "skipped", "not_run", "blocked"} or loops.get("ran") is False:
        return []
    if not node_ran(run_ledger, HIC_LOOPS_NODE):
        return []

    fdr_cut = _finite_float(loops.get("fdr_threshold"))
    caveats = [HIC_LOOP_ASSOCIATIVE]
    if _low_replication(findings):
        caveats.append(LOW_REPLICATION)

    signals: list[EvidenceSignal] = []
    seen_genes: set[str] = set()

    for row in _candidate_loop_rows(loops):
        qval = _finite_float(
            row.get("qval")
            if row.get("qval") is not None
            else row.get("la_exp.donut.qval", row.get("fdr"))
        )
        if qval is not None and fdr_cut is not None and qval > fdr_cut:
            continue
        val = qval if qval is not None else _finite_float(row.get("enrichment", row.get("count")))
        ctx = _loop_context(row)
        genes = (
            _split_gene_field(row.get("anchor1_genes"))
            + _split_gene_field(row.get("anchor1_gene"))
            + _split_gene_field(row.get("anchor2_genes"))
            + _split_gene_field(row.get("anchor2_gene"))
        )
        for gene in genes:
            gkey = gene.lower()
            if gkey in seen_genes:
                continue
            seen_genes.add(gkey)
            signals.append(
                EvidenceSignal(
                    entity=gene,
                    entity_kind="gene",
                    modality="HiC",
                    measure=HIC_LOOP_MEASURE,
                    audited_node_ref=HIC_LOOPS_NODE,
                    value=val,
                    direction="na",
                    caveats_inherited=list(caveats),
                    context=ctx,
                )
            )
            if len(signals) >= _MAX_LOOP_GENES:
                return signals

    return signals
