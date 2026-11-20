"""3D Genome Architecture (`HiC` / Micro-C) narrative plugin (ADR-082).

Translates ``GenomeArchAgent`` structured outputs into governed
``NarrativeBlock`` instances under the ``hic.*`` block namespace with
``modality="hic"``:

  - ``hic.preprocessing`` (``analysis="preprocessing"``; raw FASTQ / BAM /
    ``.pairs`` / ``.hic`` intake into ``.cool`` / ``.mcool``)
  - ``hic.qc`` (``analysis="qc"``; measured valid pairs, cis/trans ratio, ICE
    balancing, P(s) decay slope, and multi-replicate HiCRep SCC)
  - ``hic.compartments`` (``analysis="compartments"``; cis eigenvector E1 A/B
    decomposition with explicit GC / gene-density sign phasing audit)
  - ``hic.tads`` (``analysis="tads"``; Crane et al. 2015 insulation score and
    prominence-filtered TAD boundaries)
  - ``hic.loops`` (``analysis="loops"``; expected-calibrated donut focal contact
    enrichment and positional anchor gene annotation)
  - ``hic.differential_3d`` (``analysis="differential_3d"``; replicate-level
    empirical-Bayes moderated compartment/insulation testing, or honest
    limitation block when under-replicated / skipped)

Strict invariants (ADR-002, ADR-011, ADR-082, C1-C4):
  - Every number and proper noun in a public ``claim`` is grounded in that
    block's ``evidence`` list (``W-CLAIM``).
  - Zero causal verbs in claims, evidence items, caveats, or methods.
  - Non-executed, skipped, descriptive-only, or failed stages emit honest
    diagnostic/limitation blocks with ``claim=""``.
"""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock

_OK_STATUSES = {"success", "ok", "done"}
_TAG_STRIP_RE = re.compile(r"<[^>]+>")


def unwrap_hic_findings(agent_result: dict | None) -> dict[str, Any]:
    """Return the flat Hi-C findings dict from a ``GenomeArchAgent`` envelope."""
    if not isinstance(agent_result, dict):
        return {}
    findings = agent_result.get("findings")
    if isinstance(findings, dict):
        return findings
    return agent_result


def _clean_token(val: Any, default: str = "unspecified") -> str:
    if val is None:
        return default
    txt = _TAG_STRIP_RE.sub("", str(val))
    txt = re.sub(r"[\r\n\t]+", " ", txt).strip()
    return txt or default


def _finite_round(val: Any, digits: int = 4) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _as_int(val: Any) -> int | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        if isinstance(val, float) and not math.isfinite(val):
            return None
        return int(val)
    except (TypeError, ValueError, OverflowError):
        # Wave 4A: int(inf) raises OverflowError; int("nan") ValueError.
        return None


def _first_key(d: dict[str, Any], *keys: str) -> Any:
    """First non-None value among canonical script keys and legacy aliases."""
    if not isinstance(d, dict):
        return None
    for key in keys:
        if d.get(key) is not None:
            return d.get(key)
    return None


def _svg_map(findings: dict[str, Any]) -> dict[str, str]:
    """SVG paths from ``figures_svg`` and ``figure_bundles[k]["svg"]``."""
    out: dict[str, str] = {}
    figs_svg = findings.get("figures_svg")
    if isinstance(figs_svg, dict):
        for k, v in figs_svg.items():
            if isinstance(v, str) and v:
                out[str(k)] = v
    bundles = findings.get("figure_bundles")
    if isinstance(bundles, dict):
        for k, b in bundles.items():
            if isinstance(b, dict) and isinstance(b.get("svg"), str) and b["svg"]:
                out.setdefault(str(k), b["svg"])
    return out


def _png_map(findings: dict[str, Any]) -> dict[str, str]:
    """PNG paths from ``figures`` and ``figure_bundles[k]["png"]``."""
    out: dict[str, str] = {}
    figs = findings.get("figures")
    if isinstance(figs, dict):
        for k, v in figs.items():
            if isinstance(v, str) and v:
                out[str(k)] = v
    bundles = findings.get("figure_bundles")
    if isinstance(bundles, dict):
        for k, b in bundles.items():
            if isinstance(b, dict) and isinstance(b.get("png"), str) and b["png"]:
                out.setdefault(str(k), b["png"])
    return out


def _base_resolution(pre: dict[str, Any]) -> int | None:
    base = _as_int(_first_key(pre, "base_resolution", "base_res", "resolution"))
    if base is not None:
        return base
    res_vals: list[int] = []
    for rec in pre.get("samples") or []:
        if isinstance(rec, dict):
            for r in rec.get("resolutions") or []:
                ri = _as_int(r)
                if ri is not None:
                    res_vals.append(ri)
    return min(res_vals) if res_vals else None


def _loop_pe_count(loops: dict[str, Any]) -> int | None:
    n = _as_int(loops.get("n_promoter_enhancer_candidates"))
    if n is not None:
        return n
    counts = loops.get("loop_class_counts")
    if isinstance(counts, dict) and counts:
        return _as_int(counts.get("promoter_enhancer"))
    return None


def _stage_state(stage: Any) -> str:
    if not isinstance(stage, dict) or not stage:
        return "absent"
    st = str(stage.get("status") or "").lower().strip()
    if st in {"error", "failed"}:
        return "error"
    if st in {"skipped", "not_run", "blocked"} or stage.get("ran") is False:
        return "skipped"
    if st in _OK_STATUSES or stage.get("ran") is True:
        return "ok"
    return "skipped"


def _ev(label: str, value: Any, source: str = "genome_arch_agent", path: str | None = None) -> EvidenceItem:
    return EvidenceItem(label=label, value=value, source=source, path=path)


def _figure_entries(findings: dict[str, Any], key: str, caption: str) -> list[dict[str, str]]:
    raw = _png_map(findings).get(key)
    if not raw or not isinstance(raw, str):
        return []
    rec: dict[str, str] = {"id": f"hic.{key}", "path": raw, "caption": caption}
    svg_p = _svg_map(findings).get(key)
    if isinstance(svg_p, str) and svg_p:
        rec["svg_path"] = svg_p
    return [rec]


def _diagnostic_block(
    block_id: str,
    analysis: str,
    title: str,
    stage: dict[str, Any],
    state: str,
) -> NarrativeBlock:
    raw_details = stage.get("details")
    if isinstance(raw_details, dict):
        raw_details = (
            raw_details.get("details")
            or raw_details.get("reason")
            or raw_details.get("error_type")
            or raw_details.get("message")
            or stage.get("stage")
        )
    reason = _clean_token(
        stage.get("reason")
        or stage.get("message")
        or raw_details
        or stage.get("error_type")
        or f"stage reported status={stage.get('status')}",
        default="stage did not execute",
    )
    status = "error" if state == "error" else "limitation"
    return NarrativeBlock(
        id=block_id,
        modality="hic",
        analysis=analysis,
        block_type="limitation",
        title=title,
        status=status,
        confidence="insufficient",
        claim="",
        error=reason if status == "error" else None,
        caveats=[
            Caveat(
                text=f"{title} did not produce inferential output ({reason}); no unmeasured values were substituted.",
                severity="warning" if status == "error" else "info",
            )
        ],
        metadata={"validation_level": stage.get("validation_level") or "production"},
    )


class HiCNarrator:
    """Governed C1-C4 narrator for ``GenomeArchAgent`` (`HiC` / Micro-C)."""

    name = "hic"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        return (
            agent_name == "genome_arch_agent"
            and isinstance(agent_result, dict)
            and bool(unwrap_hic_findings(agent_result))
        )

    def collect(
        self,
        agent_name: str,
        agent_result: dict,
        context: dict | None = None,
    ) -> list[NarrativeBlock]:
        findings = unwrap_hic_findings(agent_result)
        blocks: list[NarrativeBlock] = []

        fq_align = findings.get("fastq_alignment")
        fq_state = _stage_state(fq_align)
        if fq_state in {"skipped", "error"} and isinstance(fq_align, dict):
            blocks.append(
                _diagnostic_block(
                    "hic.fastq_alignment",
                    "fastq_alignment",
                    "Hi-C FASTQ chimeric alignment",
                    fq_align,
                    fq_state,
                )
            )

        pre = findings.get("preprocessing")
        pre_state = _stage_state(pre)
        if pre_state == "ok" and isinstance(pre, dict) and (pre.get("cooler_files") or pre.get("n_samples")):
            blocks.append(self._preprocessing_block(pre))
        elif pre_state in {"skipped", "error"} and isinstance(pre, dict):
            blocks.append(
                _diagnostic_block(
                    "hic.preprocessing",
                    "preprocessing",
                    "3D contact preprocessing",
                    pre,
                    pre_state,
                )
            )

        bal = findings.get("balancing") or findings.get("qc")
        bal_state = _stage_state(bal)
        if bal_state == "ok" and isinstance(bal, dict):
            blocks.append(self._qc_block(bal, findings))
        elif bal_state in {"skipped", "error"} and isinstance(bal, dict):
            blocks.append(
                _diagnostic_block(
                    "hic.qc", "qc", "Hi-C contact QC & ICE balancing", bal, bal_state
                )
            )

        topo = findings.get("topology") if isinstance(findings.get("topology"), dict) else findings

        comp = topo.get("compartments") if isinstance(topo, dict) else None
        comp_state = _stage_state(comp)
        if comp_state == "ok" and isinstance(comp, dict):
            blocks.append(self._compartments_block(comp, findings))
        elif comp_state in {"skipped", "error"} and isinstance(comp, dict):
            blocks.append(
                _diagnostic_block(
                    "hic.compartments",
                    "compartments",
                    "A/B compartment decomposition",
                    comp,
                    comp_state,
                )
            )

        tads = topo.get("tads") if isinstance(topo, dict) else None
        tads_state = _stage_state(tads)
        if tads_state == "ok" and isinstance(tads, dict):
            blocks.append(
                self._tads_block(tads, findings, attach_track_fig=(comp_state != "ok"))
            )
        elif tads_state in {"skipped", "error"} and isinstance(tads, dict):
            blocks.append(
                _diagnostic_block(
                    "hic.tads",
                    "tads",
                    "Topologically associating domains (insulation)",
                    tads,
                    tads_state,
                )
            )

        loops = topo.get("loops") if isinstance(topo, dict) else None
        loops_state = _stage_state(loops)
        if loops_state == "ok" and isinstance(loops, dict):
            blocks.append(self._loops_block(loops, findings))
        elif loops_state in {"skipped", "error"} and isinstance(loops, dict):
            blocks.append(
                _diagnostic_block(
                    "hic.loops",
                    "loops",
                    "Chromatin loop calling",
                    loops,
                    loops_state,
                )
            )

        diff3d = topo.get("differential_3d") if isinstance(topo, dict) else None
        diff_state = _stage_state(diff3d)
        if isinstance(diff3d, dict):
            if diff_state == "ok" and diff3d.get("ran") is True and not diff3d.get("descriptive_only"):
                blocks.append(self._differential_3d_block(diff3d, findings))
            elif diff3d.get("descriptive_only") or diff_state in {"skipped", "error"}:
                blocks.append(
                    _diagnostic_block(
                        "hic.differential_3d",
                        "differential_3d",
                        "Differential 3D compartment & insulation testing",
                        diff3d,
                        "error" if diff_state == "error" else "skipped",
                    )
                )

        return blocks

    def _preprocessing_block(self, pre: dict[str, Any]) -> NarrativeBlock:
        coolers = [str(c) for c in (pre.get("cooler_files") or []) if c]
        n_samples = _as_int(pre.get("n_samples"))
        if n_samples is None:
            n_samples = len(coolers)
        base_res = _base_resolution(pre)
        evidence = [
            _ev("Assay", "Hi-C", "hic_preprocess"),
            _ev("Samples preprocessed", n_samples, "hic_preprocess"),
            _ev("Cooler artifacts produced", len(coolers), "hic_preprocess"),
        ]
        if base_res is not None:
            evidence.append(_ev("Base bin resolution (bp)", base_res, "hic_preprocess"))
        for idx, cpath in enumerate(coolers, start=1):
            evidence.append(_ev(f"Cooler artifact {idx}", "cooler matrix", "hic_preprocess", path=cpath))

        claim = (
            f"Hi-C contact preprocessing assembled {len(coolers)} cooler matrix "
            f"artifact(s) across {n_samples} sample(s)."
        )
        return NarrativeBlock(
            id="hic.preprocessing",
            modality="hic",
            analysis="preprocessing",
            block_type="result",
            title="3D contact preprocessing & matrix build",
            status="success",
            confidence="high",
            claim=claim,
            evidence=evidence,
            caveats=[
                Caveat(
                    "Contact binning and deduplication counts describe technical "
                    "data intake and matrix construction rather than 3D biological "
                    "differences.",
                    severity="info",
                )
            ],
            metrics={"n_samples": n_samples, "n_coolers": len(coolers), "base_resolution": base_res},
            metadata={"validation_level": pre.get("validation_level") or "production"},
        )

    def _qc_block(self, bal: dict[str, Any], findings: dict[str, Any]) -> NarrativeBlock:
        # Wave 4A: an explicit 0 must stay 0 (no falsy ``or``); the script
        # emits ``per_file`` (``per_sample`` is a legacy alias).
        n_samples = _as_int(bal.get("n_samples"))
        if n_samples is None:
            per = bal.get("per_file") or bal.get("per_sample") or []
            n_ok = sum(
                1 for r in per
                if not isinstance(r, dict) or r.get("status") in (None, "success")
            )
            n_samples = max(1, n_ok)
        n_valid = _as_int(bal.get("n_valid_pairs"))
        cis_trans = _finite_round(bal.get("cis_trans_ratio"), 3)
        long_cis = _finite_round(bal.get("long_range_cis_fraction", bal.get("long_range_cis")), 3)
        ps_slope = _finite_round(bal.get("ps_decay_slope"), 3)
        mean_scc = _finite_round(bal.get("mean_replicate_scc"), 3)

        evidence: list[EvidenceItem] = [
            _ev("Assay", "Hi-C", "hic_qc_and_balance"),
            _ev("Normalization", "ICE", "hic_qc_and_balance"),
            _ev("Samples evaluated", n_samples, "hic_qc_and_balance"),
        ]
        if n_valid is not None:
            evidence.append(_ev("Valid contact pairs", n_valid, "hic_qc_and_balance"))
        if cis_trans is not None:
            evidence.append(_ev("Cis/trans contact ratio", cis_trans, "hic_qc_and_balance"))
        if long_cis is not None:
            evidence.append(_ev("Long-range cis fraction (>=20 kb)", long_cis, "hic_qc_and_balance"))
        if ps_slope is not None:
            evidence.append(_ev("Contact decay P(s) log-log slope", ps_slope, "hic_qc_and_balance"))
        if mean_scc is not None:
            evidence.append(_ev("Mean replicate stratum-adjusted correlation (SCC)", mean_scc, "hic_qc_and_balance"))
        if bal.get("ps_curve_tsv"):
            evidence.append(_ev("Contact probability P(s) table", "tsv", "hic_qc_and_balance", path=str(bal["ps_curve_tsv"])))
        if bal.get("scc_tsv"):
            evidence.append(_ev("HiCRep SCC matrix", "tsv", "hic_qc_and_balance", path=str(bal["scc_tsv"])))

        if n_valid is not None and cis_trans is not None:
            claim = (
                f"Hi-C contact QC and ICE balancing evaluated {n_samples} sample(s) "
                f"with {n_valid} valid contact pairs and cis/trans ratio {cis_trans}."
            )
        else:
            claim = (
                f"Hi-C contact QC and ICE balancing completed across {n_samples} sample(s)."
            )

        contact_map_caption = (
            "ICE-balanced Hi-C cis contact frequency matrix."
            if bool(bal.get("balanced_files"))
            else "Raw Hi-C cis contact frequency matrix (unbalanced fallback)."
        )
        figs = _figure_entries(
            findings,
            "hic_contact_decay_ps",
            "Distance-dependent contact probability P(s) decay curve across genomic separation.",
        ) + _figure_entries(
            findings,
            "hic_contact_map",
            contact_map_caption,
        )
        caveats = [
            Caveat(
                "ICE matrix balancing equalizes marginal bin visibility across "
                "mappable genomic bins; low-coverage centromeric or repetitive "
                "bins are masked as NaN rather than imputed.",
                severity="info",
            )
        ]
        for w in bal.get("warnings") or []:
            caveats.append(Caveat(_clean_token(w), severity="warning"))

        return NarrativeBlock(
            id="hic.qc",
            modality="hic",
            analysis="qc",
            block_type="qc",
            title="Hi-C contact QC, ICE balancing & P(s) decay",
            status="success",
            confidence="high" if (bal.get("pass_qc") is not False) else "low",
            claim=claim,
            evidence=evidence,
            figures=figs,
            caveats=caveats,
            metrics={
                "n_samples": n_samples,
                "n_valid_pairs": n_valid,
                "cis_trans_ratio": cis_trans,
                "ps_decay_slope": ps_slope,
                "mean_replicate_scc": mean_scc,
            },
            metadata={"validation_level": bal.get("validation_level") or "production"},
        )

    def _compartments_block(self, comp: dict[str, Any], findings: dict[str, Any]) -> NarrativeBlock:
        pct_a = _finite_round(comp.get("pct_A"), 2)
        pct_b = _finite_round(comp.get("pct_B"), 2)
        res_bp = _as_int(comp.get("resolution_bp") or comp.get("resolution"))
        phased = comp.get("pc1_validated") is True
        phasing_track = _clean_token(
            comp.get("phasing_method") or comp.get("phasing_track"), "unphased"
        )
        n_bins = _as_int(comp.get("n_valid_bins") or comp.get("n_bins"))

        evidence: list[EvidenceItem] = [
            _ev("Assay", "Hi-C", "hic_topology"),
            _ev("Analysis", "cis eigenvector compartment decomposition", "hic_topology"),
            _ev("Sign orientation validated", "yes" if phased else "no", "hic_topology"),
            _ev("Phasing reference", phasing_track, "hic_topology"),
        ]
        if res_bp is not None:
            evidence.append(_ev("Compartment bin resolution (bp)", res_bp, "hic_topology"))
        if n_bins is not None:
            evidence.append(_ev("Valid genomic bins", n_bins, "hic_topology"))
        # Wave 4A: "active"/"inactive" only when the E1 sign was phased
        # (pc1_validated is True); otherwise report raw E1 polarity.
        if phased:
            a_label, b_label = "Active compartment fraction", "Inactive compartment fraction"
        else:
            a_label = "Positive-E1 bin fraction (unphased sign)"
            b_label = "Negative-E1 bin fraction (unphased sign)"
        if pct_a is not None:
            evidence.append(_ev(a_label, f"{pct_a}%", "hic_topology"))
        if pct_b is not None:
            evidence.append(_ev(b_label, f"{pct_b}%", "hic_topology"))
        if comp.get("compartments_tsv"):
            evidence.append(_ev("Compartment eigenvector table", "tsv", "hic_topology", path=str(comp["compartments_tsv"])))

        if pct_a is not None and pct_b is not None and res_bp is not None and phased:
            claim = (
                f"Hi-C cis eigenvector compartment decomposition at {res_bp} bp "
                f"resolution classified {pct_a}% active and {pct_b}% "
                f"inactive compartment bins."
            )
        elif pct_a is not None and pct_b is not None and res_bp is not None:
            claim = (
                f"Hi-C cis eigenvector compartment decomposition at {res_bp} bp "
                f"resolution found {pct_a}% positive-E1 and {pct_b}% negative-E1 "
                f"bins; the eigenvector sign is unphased, so these are not "
                f"active/inactive compartment assignments."
            )
        else:
            claim = "Hi-C cis eigenvector compartment decomposition completed."

        caveats = [
            Caveat(
                "A/B compartment calls reflect long-range plaid contact "
                "co-segregation only; B-compartment assignment indicates relative "
                "chromatin compaction and is not standalone proof of transcriptional "
                "inactivity without matched RNA measurements.",
                severity="info",
            )
        ]
        if not phased:
            caveats.append(
                Caveat(
                    "Eigenvector sign was not phased against GC content or TSS gene "
                    "density (`pc1_validated=False`), so per-chromosome A vs B polarity "
                    "is arbitrary until phased against a reference FASTA or GTF.",
                    severity="warning",
                )
            )
        for w in comp.get("warnings") or []:
            caveats.append(Caveat(_clean_token(w), severity="warning"))

        return NarrativeBlock(
            id="hic.compartments",
            modality="hic",
            analysis="compartments",
            block_type="result",
            title="A/B chromatin compartment decomposition",
            status="success",
            confidence="medium" if phased else "low",
            claim=claim,
            evidence=evidence,
            figures=_figure_entries(
                findings,
                "hic_compartment_insulation_track",
                "Cis compartment E1 eigenvector and log2 insulation score profile.",
            ),
            caveats=caveats,
            metrics={
                "pct_A": pct_a,
                "pct_B": pct_b,
                "resolution_bp": res_bp,
                "pc1_validated": phased,
                "phasing_track": phasing_track,
            },
            metadata={"validation_level": comp.get("validation_level") or "production"},
        )

    def _tads_block(
        self,
        tads: dict[str, Any],
        findings: dict[str, Any],
        *,
        attach_track_fig: bool = False,
    ) -> NarrativeBlock:
        n_tads = _as_int(tads.get("n_tads"))
        n_bounds = _as_int(tads.get("n_boundaries"))
        med_kb = _finite_round(tads.get("median_size_kb"), 1)
        win_bp = _as_int(_first_key(tads, "window_size", "window_bp"))

        evidence: list[EvidenceItem] = [
            _ev("Assay", "Hi-C", "hic_topology"),
            _ev("Method", "diamond insulation score", "hic_topology"),
        ]
        if n_tads is not None:
            evidence.append(_ev("TAD intervals identified", n_tads, "hic_topology"))
        if n_bounds is not None:
            evidence.append(_ev("Prominent insulation boundaries", n_bounds, "hic_topology"))
        if med_kb is not None:
            evidence.append(_ev("Median TAD size (kb)", med_kb, "hic_topology"))
        if win_bp is not None:
            evidence.append(_ev("Insulation diamond window (bp)", win_bp, "hic_topology"))
        if tads.get("tads_bed"):
            evidence.append(_ev("TAD domain BED", "bed", "hic_topology", path=str(tads["tads_bed"])))
        if tads.get("insulation_tsv"):
            evidence.append(_ev("Insulation score track", "tsv", "hic_topology", path=str(tads["insulation_tsv"])))

        if n_tads is not None and med_kb is not None:
            claim = (
                f"Hi-C diamond insulation score analysis identified {n_tads} TAD "
                f"intervals with median size {med_kb} kb."
            )
        elif n_tads is not None:
            claim = f"Hi-C diamond insulation score analysis identified {n_tads} TAD intervals."
        else:
            claim = "Hi-C diamond insulation score analysis completed."

        caveats = [
            Caveat(
                "Insulation-defined TAD intervals depend on matrix resolution and "
                "sliding-diamond window scale; boundaries denote local contact "
                "depletion rather than impermeable enhancer-blocking insulators.",
                severity="info",
            )
        ]
        for w in tads.get("warnings") or []:
            caveats.append(Caveat(_clean_token(w), severity="warning"))

        figs = (
            _figure_entries(
                findings,
                "hic_compartment_insulation_track",
                "Crane log2 insulation score profile and prominent TAD boundaries.",
            )
            if attach_track_fig
            else []
        )
        return NarrativeBlock(
            id="hic.tads",
            modality="hic",
            analysis="tads",
            block_type="result",
            title="Topologically associating domains (insulation score)",
            status="success",
            confidence="medium",
            claim=claim,
            evidence=evidence,
            figures=figs,
            caveats=caveats,
            metrics={
                "n_tads": n_tads,
                "n_boundaries": n_bounds,
                "median_size_kb": med_kb,
                "window_bp": win_bp,
            },
            metadata={"validation_level": tads.get("validation_level") or "production"},
        )

    def _loops_block(self, loops: dict[str, Any], findings: dict[str, Any]) -> NarrativeBlock:
        n_loops = _as_int(loops.get("n_loops"))
        res_bp = _as_int(loops.get("resolution_bp") or loops.get("resolution"))
        n_pe = _loop_pe_count(loops)
        loops_path = _first_key(loops, "loops_path", "loops_bedpe")
        ann_path = _first_key(loops, "annotated_loops_path", "annotated_loops_tsv")

        evidence: list[EvidenceItem] = [
            _ev("Assay", "Hi-C", "hic_topology"),
            _ev("Method", "donut local background enrichment", "hic_topology"),
        ]
        if n_loops is not None:
            evidence.append(_ev("Chromatin loops called", n_loops, "hic_topology"))
        if res_bp is not None:
            evidence.append(_ev("Loop calling resolution (bp)", res_bp, "hic_topology"))
        if n_pe is not None:
            evidence.append(_ev("Promoter-enhancer loop candidates", n_pe, "hic_topology"))
        if loops_path:
            evidence.append(_ev("Called chromatin loops", "bedpe", "hic_topology", path=str(loops_path)))
        if ann_path:
            evidence.append(_ev("Anchor-annotated loops", "tsv", "hic_topology", path=str(ann_path)))

        if n_loops is not None and res_bp is not None:
            claim = (
                f"Hi-C donut local background enrichment at {res_bp} bp resolution "
                f"identified {n_loops} chromatin loop(s)."
            )
        elif n_loops is not None:
            claim = f"Hi-C donut local background enrichment identified {n_loops} chromatin loop(s)."
        else:
            claim = "Hi-C chromatin loop analysis completed."

        caveats = [
            Caveat(
                "Focal 3D chromatin loops represent elevated physical contact "
                "frequency above local distance-stratified background; physical "
                "proximity with a TSS is associative and does not establish "
                "transcriptional regulation without perturbation evidence.",
                severity="info",
            )
        ]
        for w in loops.get("warnings") or []:
            caveats.append(Caveat(_clean_token(w), severity="warning"))

        return NarrativeBlock(
            id="hic.loops",
            modality="hic",
            analysis="loops",
            block_type="result",
            title="Focal chromatin loop calling",
            status="success",
            confidence="medium",
            claim=claim,
            evidence=evidence,
            caveats=caveats,
            metrics={
                "n_loops": n_loops,
                "resolution_bp": res_bp,
                "n_promoter_enhancer_candidates": n_pe,
            },
            metadata={"validation_level": loops.get("validation_level") or "production"},
        )

    def _differential_3d_block(self, diff3d: dict[str, Any], findings: dict[str, Any]) -> NarrativeBlock:
        comps = [c for c in (diff3d.get("comparisons") or []) if isinstance(c, dict)]

        def _contrast_tested(c: dict[str, Any]) -> bool:
            if c.get("ran") is False or c.get("descriptive_only") is True:
                return False
            st = str(c.get("status") or "").strip().lower()
            if st in {"skipped", "error", "failed"}:
                return False
            return True

        tested_comps = [c for c in comps if _contrast_tested(c)]
        skipped_comps = [c for c in comps if not _contrast_tested(c)]
        if comps:
            n_comps = len(tested_comps)
        else:
            n_comps = _as_int(diff3d.get("n_comparisons")) or 0

        topo = findings.get("topology") if isinstance(findings.get("topology"), dict) else findings
        comp_stage = (topo or {}).get("compartments") if isinstance(topo, dict) else None
        raw_pc1_val = diff3d.get("pc1_validated")
        if raw_pc1_val is True:
            is_phased = True
        elif raw_pc1_val is False:
            is_phased = False
        else:
            is_phased = bool(
                isinstance(comp_stage, dict) and comp_stage.get("pc1_validated") is True
            )

        n_sig_comp = _as_int(
            diff3d.get("n_significant_compartment_bins")
            if diff3d.get("n_significant_compartment_bins") is not None
            else diff3d.get("n_differential_compartment_bins")
        )
        if n_sig_comp is None:
            n_sig_comp = sum(
                _as_int(
                    _first_key(c, "n_differential_compartment_bins",
                               "n_sig_compartments", "n_significant_e1")
                ) or 0
                for c in tested_comps
            )
        n_switches = _as_int(diff3d.get("n_ab_switches"))
        if n_switches is None:
            n_switches = sum(
                (
                    _as_int(c.get("n_ab_switches"))
                    if _as_int(c.get("n_ab_switches")) is not None
                    else (_as_int(c.get("n_a_to_b")) or 0) + (_as_int(c.get("n_b_to_a")) or 0)
                )
                for c in tested_comps
            )
        n_unphased_flips = _as_int(diff3d.get("n_unphased_sign_flips"))
        if n_unphased_flips is None:
            n_unphased_flips = sum(
                _as_int(c.get("n_unphased_sign_flips")) or 0
                for c in tested_comps
            )
        n_sig_ins = _as_int(
            diff3d.get("n_differential_insulation_bins")
            if diff3d.get("n_differential_insulation_bins") is not None
            else diff3d.get("n_significant_boundaries", diff3d.get("n_differential_boundaries"))
        )
        if n_sig_ins is None:
            n_sig_ins = sum(
                _as_int(
                    _first_key(c, "n_differential_insulation_bins",
                               "n_differential_boundaries", "n_sig_insulation")
                ) or 0
                for c in tested_comps
            )

        evidence: list[EvidenceItem] = [
            _ev("Assay", "Hi-C", "hic_topology"),
            _ev("Method", "empirical-Bayes moderated linear model with BH-FDR", "hic_topology"),
            _ev("Contrasts tested", n_comps, "hic_topology"),
            _ev("Differential compartment bins", n_sig_comp, "hic_topology"),
        ]
        if is_phased:
            evidence.append(_ev("Compartment switches", n_switches, "hic_topology"))
        else:
            if n_switches > 0:
                evidence.append(_ev("Compartment switches (phased chroms)", n_switches, "hic_topology"))
            evidence.append(_ev("Unphased E1 sign flips", n_unphased_flips, "hic_topology"))
        evidence.append(_ev("Differential insulation bins", n_sig_ins, "hic_topology"))
        for c in tested_comps:
            label = _clean_token(c.get("contrast") or c.get("comparison") or f"{c.get('test')}_vs_{c.get('reference')}")
            comp_csv = c.get("differential_compartments_csv") or c.get("compartments_csv")
            ins_csv = c.get("differential_insulation_csv") or c.get("insulation_csv")
            if comp_csv:
                evidence.append(_ev(f"Differential compartments ({label})", "csv", "hic_topology", path=str(comp_csv)))
            if ins_csv:
                evidence.append(_ev(f"Differential insulation ({label})", "csv", "hic_topology", path=str(ins_csv)))

        if is_phased:
            claim = (
                f"Hi-C replicate-level differential topology testing across {n_comps} "
                f"contrast(s) identified {n_sig_comp} differential compartment bin(s), "
                f"including {n_switches} compartment switch(es), and {n_sig_ins} "
                f"differential insulation bin(s)."
            )
        elif n_switches > 0:
            claim = (
                f"Hi-C replicate-level differential topology testing across {n_comps} "
                f"contrast(s) identified {n_sig_comp} differential compartment bin(s), "
                f"including {n_switches} phased compartment switch(es) and "
                f"{n_unphased_flips} unphased E1 sign flip(s), and {n_sig_ins} "
                f"differential insulation bin(s)."
            )
        else:
            claim = (
                f"Hi-C replicate-level differential topology testing across {n_comps} "
                f"contrast(s) identified {n_sig_comp} differential compartment bin(s), "
                f"including {n_unphased_flips} unphased E1 sign flip(s) (A/B polarity unvalidated), "
                f"and {n_sig_ins} differential insulation bin(s)."
            )

        caveats = [
            Caveat(
                "Differential compartment and insulation shifts are observational "
                "3D reorganization measurements across biological replicates; "
                "co-occurring transcriptional changes remain associative.",
                severity="info",
            )
        ]
        if not is_phased:
            caveats.append(
                Caveat(
                    "Eigenvector E1 was not phased against a verified GC or gene-density track "
                    "across all analyzed chromosomes; directional A-to-B / B-to-A compartment "
                    "switch labels are suppressed and sign reversals are reported as unphased sign flips.",
                    severity="warning",
                )
            )
        for sc in skipped_comps:
            sc_label = _clean_token(sc.get("contrast") or sc.get("comparison") or f"{sc.get('test')}_vs_{sc.get('reference')}")
            sc_reason = _clean_token(sc.get("reason") or sc.get("status") or "skipped")
            caveats.append(
                Caveat(
                    f"Contrast {sc_label} was skipped from statistical testing ({sc_reason}).",
                    severity="warning",
                )
            )
        for w in diff3d.get("warnings") or []:
            caveats.append(Caveat(_clean_token(w), severity="warning"))

        return NarrativeBlock(
            id="hic.differential_3d",
            modality="hic",
            analysis="differential_3d",
            block_type="result",
            title="Differential 3D compartment & insulation reorganization",
            status="success",
            confidence="medium" if is_phased else "low",
            claim=claim,
            evidence=evidence,
            figures=_figure_entries(
                findings,
                "hic_differential_compartments",
                "Replicate-moderated differential compartment E1 shifts (delta E1 vs -log10 BH-FDR).",
            ),
            caveats=caveats,
            metrics={
                "n_comparisons": n_comps,
                "n_significant_compartment_bins": n_sig_comp,
                "n_ab_switches": n_switches,
                "n_unphased_sign_flips": n_unphased_flips,
                "n_significant_boundaries": n_sig_ins,
                "pc1_validated": is_phased,
            },
            metadata={"validation_level": diff3d.get("validation_level") or "production"},
        )

    def methods(
        self,
        agent_name: str,
        agent_result: dict,
        context: dict | None = None,
    ) -> list[str]:
        findings = unwrap_hic_findings(agent_result)
        out: list[str] = []
        pre = findings.get("preprocessing") or {}
        if _stage_state(pre) == "ok" and pre.get("cooler_files"):
            out.append(
                f"3D contact preprocessing: contact inputs across "
                f"{_clean_token(pre.get('n_samples'), '?')} sample(s) were "
                f"standardized into multi-resolution cooler matrices"
                + (f" (base bin {_base_resolution(pre)} bp)"
                   if _base_resolution(pre) is not None else "")
                + " via pairtools/cooler/hicstraw in aria-hic-env."
            )
        bal = findings.get("balancing") or findings.get("qc") or {}
        if _stage_state(bal) == "ok":
            mean_scc = _finite_round(bal.get("mean_replicate_scc"), 3)
            scc_clause = (
                f"; multi-sample reproducibility was quantified by distance-"
                f"stratified HiCRep stratum-adjusted correlation (mean SCC "
                f"{mean_scc})"
                if mean_scc is not None
                else ""
            )
            out.append(
                "Hi-C contact QC and normalization: contact matrices were "
                "iteratively balanced with cooler ICE normalization on staged "
                "copies, exact valid-pair sums and cis/trans contact fractions "
                "were tabulated from matrix pixels, and distance-dependent "
                f"contact decay P(s) slopes were estimated with cooltools{scc_clause}."
            )
        topo = findings.get("topology") if isinstance(findings.get("topology"), dict) else findings
        comp = (topo or {}).get("compartments") or {}
        if _stage_state(comp) == "ok":
            phase_ref = _first_key(comp, "phasing_method", "phasing_track")
            if comp.get("pc1_validated") is True:
                phase_txt = (
                    f"phased against {_clean_token(phase_ref)}"
                    if phase_ref else "sign-phased against a reference track"
                )
            else:
                phase_txt = "reported unphased with explicit sign-arbitrariness disclosure"
            out.append(
                f"A/B compartments: cis eigenvector decomposition (cooltools "
                f"eigs_cis, E1) was computed at {_clean_token(comp.get('resolution_bp') or comp.get('resolution'), '?')} "
                f"bp resolution and {phase_txt}."
            )
        tads = (topo or {}).get("tads") or {}
        if _stage_state(tads) == "ok":
            out.append(
                f"Topologically associating domains: diamond log2 insulation "
                f"scores (Crane et al., 2015"
                + (f"; window {_as_int(_first_key(tads, 'window_size', 'window_bp'))} bp"
                   if _as_int(_first_key(tads, 'window_size', 'window_bp')) is not None else "")
                + ") "
                "and prominence-filtered boundary minima were called with cooltools insulation."
            )
        loops = (topo or {}).get("loops") or {}
        if _stage_state(loops) == "ok":
            out.append(
                f"Chromatin loops: focal contact enrichment over local donut/vertical/"
                f"horizontal backgrounds was tested with cooltools dots at "
                f"{_clean_token(loops.get('resolution_bp') or loops.get('resolution'), '?')} bp "
                "resolution after cis expected-contact calibration."
            )
        diff3d = (topo or {}).get("differential_3d") or {}
        if _stage_state(diff3d) == "ok" and diff3d.get("ran") is True and not diff3d.get("descriptive_only"):
            diff_e1_desc = (
                "per-replicate phased E1 compartment scores"
                if diff3d.get("pc1_validated") is True
                else "per-replicate unphased E1 compartment scores (sign-aligned across replicates without A/B polarity validation)"
            )
            out.append(
                f"Differential 3D organization: {diff_e1_desc} "
                "and log2 insulation profiles were tested across biological "
                "conditions using empirical-Bayes variance-moderated linear models "
                "with Benjamini-Hochberg FDR control."
            )
        return out

    def figures(
        self,
        agent_name: str,
        agent_result: dict,
        report_dir: Path | None = None,
    ) -> list[dict]:
        findings = unwrap_hic_findings(agent_result)
        figs_svg = _svg_map(findings)
        out: list[dict] = []
        for key, path in _png_map(findings).items():
            rec = {
                "id": f"hic.{key}",
                "path": path,
                "caption": str(key).replace("_", " "),
            }
            if figs_svg.get(key):
                rec["svg_path"] = figs_svg[key]
            out.append(rec)
        return out

    def tables(
        self,
        agent_name: str,
        agent_result: dict,
        report_dir: Path | None = None,
    ) -> list[dict]:
        findings = unwrap_hic_findings(agent_result)
        tables: list[dict] = []
        seen_paths: set[str] = set()

        def _add_table(tid: str, path: Any, label: str) -> None:
            if not path:
                return
            p_str = str(path)
            if p_str in seen_paths:
                return
            seen_paths.add(p_str)
            tables.append({"id": tid, "path": p_str, "label": label})

        bal = findings.get("balancing") or findings.get("qc") or {}
        if isinstance(bal, dict):
            _add_table("hic.ps_curve", bal.get("ps_curve_tsv"), "Hi-C P(s) contact decay curve")
            per_file = bal.get("per_file") or bal.get("per_sample") or []
            for idx, rec in enumerate(per_file, start=1):
                if not isinstance(rec, dict) or not rec.get("ps_curve_tsv"):
                    continue
                sid = _clean_token(rec.get("sample_id"), f"sample_{idx}")
                _add_table(
                    f"hic.ps_curve.{sid}",
                    rec.get("ps_curve_tsv"),
                    f"Hi-C P(s) contact decay curve ({sid})",
                )
            _add_table("hic.scc_matrix", bal.get("scc_tsv"), "HiCRep replicate SCC matrix")
        topo = findings.get("topology") if isinstance(findings.get("topology"), dict) else findings
        if isinstance(topo, dict):
            comp = topo.get("compartments") or {}
            if isinstance(comp, dict):
                _add_table("hic.compartments", comp.get("compartments_tsv"), "Hi-C E1 compartment track")
                _add_table(
                    "hic.compartments_bedgraph",
                    comp.get("compartments_bedgraph"),
                    "Hi-C E1 compartment track (bedGraph)",
                )
            tads = topo.get("tads") or {}
            if isinstance(tads, dict):
                _add_table("hic.tads", tads.get("tads_bed"), "Hi-C TAD intervals (BED)")
                _add_table("hic.insulation", tads.get("insulation_tsv"), "Hi-C diamond insulation track")
            loops = topo.get("loops") or {}
            if isinstance(loops, dict):
                loops_path = _first_key(loops, "loops_path", "loops_bedpe")
                ann_path = _first_key(loops, "annotated_loops_path", "annotated_loops_tsv")
                _add_table("hic.loops", loops_path, "Hi-C chromatin loops (BEDPE)")
                _add_table("hic.annotated_loops", ann_path, "Hi-C anchor-annotated loops")
            diff3d = topo.get("differential_3d") or {}
            if isinstance(diff3d, dict):
                for c in diff3d.get("comparisons") or []:
                    if not isinstance(c, dict):
                        continue
                    label = _clean_token(c.get("contrast") or c.get("comparison") or f"{c.get('test')}_vs_{c.get('reference')}")
                    comp_csv = _first_key(c, "differential_compartments_csv", "compartments_csv")
                    ins_csv = _first_key(c, "differential_insulation_csv", "insulation_csv")
                    desc_csv = c.get("descriptive_csv")
                    _add_table(
                        f"hic.diff_compartments.{label}",
                        comp_csv,
                        f"Differential compartments ({label})",
                    )
                    _add_table(
                        f"hic.diff_insulation.{label}",
                        ins_csv,
                        f"Differential insulation ({label})",
                    )
                    _add_table(
                        f"hic.descriptive_compartments.{label}",
                        desc_csv,
                        f"Descriptive compartment deltas ({label})",
                    )
                _add_table(
                    "hic.diff_compartments",
                    diff3d.get("differential_compartments_csv"),
                    "Differential compartments",
                )
                _add_table(
                    "hic.diff_insulation",
                    diff3d.get("differential_insulation_csv"),
                    "Differential insulation",
                )
                _add_table(
                    "hic.descriptive_compartments",
                    diff3d.get("descriptive_csv"),
                    "Descriptive compartment deltas",
                )
        return tables
