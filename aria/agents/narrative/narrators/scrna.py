"""scRNA narrative plugin wrapping the legacy scRNA narrative helpers."""

from __future__ import annotations

import re
from pathlib import Path

from aria.agents import _narrative_scrna
from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock


def _safe_id(value) -> str:
    text = str(value or "unknown")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "unknown"


def _evidence(label: str, value, source: str,
              path: str | None = None) -> EvidenceItem:
    return EvidenceItem(label=label, value=value, source=source, path=path)


def _design_issues(comp: dict, severities: set[str] | None = None) -> list[dict]:
    issues = ((comp.get("design_check") or {}).get("issues") or [])
    if severities is None:
        return [i for i in issues if isinstance(i, dict)]
    return [
        i for i in issues
        if isinstance(i, dict) and i.get("severity") in severities
    ]


def _first_design_issue(comp: dict) -> str:
    issues = _design_issues(comp)
    if not issues:
        return ""
    issue = issues[0]
    return issue.get("message") or issue.get("check") or ""


class ScrnaNarrator:
    name = "scrna"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        if agent_name not in {"scrna_agent", "rna_agent"}:
            return False
        if not isinstance(agent_result, dict):
            return False
        findings = _narrative_scrna.unwrap_scrna_findings(agent_result)
        return bool(findings)

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        findings = _narrative_scrna.unwrap_scrna_findings(agent_result)
        blocks: list[NarrativeBlock] = []
        blocks.extend(self._qc_blocks(findings))
        blocks.extend(self._data_quality_blocks(findings))
        blocks.extend(self._error_blocks(findings))
        blocks.extend(self._composition_blocks(findings))
        blocks.extend(self._pseudobulk_blocks(findings))
        blocks.extend(self._pathway_blocks(findings))
        blocks.extend(self._cellcomm_blocks(findings))
        blocks.extend(self._trajectory_blocks(findings))
        self._attach_artifacts(blocks, findings)
        return blocks

    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        findings = _narrative_scrna.unwrap_scrna_findings(agent_result)
        text = _narrative_scrna.build_scrna_methods(findings)
        return [text] if text else []

    def figures(self, agent_name: str, agent_result: dict,
                report_dir: Path | None = None) -> list[dict]:
        findings = _narrative_scrna.unwrap_scrna_findings(agent_result)
        figs = findings.get("figures") or {}
        out = []
        for key, value in figs.items():
            if isinstance(value, str):
                out.append({"id": key, "path": value, "caption": key})
            elif isinstance(value, dict):
                for sub_key, paths in value.items():
                    for path in paths or []:
                        out.append({
                            "id": f"{key}.{sub_key}",
                            "path": path,
                            "caption": f"{key} {sub_key}",
                        })
        return out

    def tables(self, agent_name: str, agent_result: dict,
               report_dir: Path | None = None) -> list[dict]:
        findings = _narrative_scrna.unwrap_scrna_findings(agent_result)
        return [
            {"id": key, "path": path, "label": key.replace("_", " ")}
            for key, path in (findings.get("tables") or {}).items()
        ]

    def _qc_blocks(self, findings: dict) -> list[NarrativeBlock]:
        qc = findings.get("qc") or {}
        if not qc:
            return []
        n_after = qc.get("n_cells_after")
        n_before = qc.get("n_cells_before")
        claim = (
            f"scRNA QC retained {n_after:,} of {n_before:,} cells."
            if n_after and n_before else
            f"scRNA QC retained {n_after:,} cells."
            if n_after else
            "scRNA QC completed."
        )
        evidence = []
        if n_before is not None:
            evidence.append(_evidence("cells before QC", n_before, "qc"))
        if n_after is not None:
            evidence.append(_evidence("cells after QC", n_after, "qc"))
        if qc.get("pct_removed") is not None:
            evidence.append(_evidence("percent removed", qc.get("pct_removed"), "qc"))
        return [NarrativeBlock(
            id="scrna.qc",
            modality="scRNA-seq",
            analysis="qc",
            block_type="qc",
            title="scRNA quality control",
            status="success",
            confidence="high",
            claim=claim,
            evidence=evidence or [_evidence("QC status", "completed", "qc")],
            metrics={k: v for k, v in qc.items() if isinstance(v, (int, float, str))},
        )]

    def _data_quality_blocks(self, findings: dict) -> list[NarrativeBlock]:
        """X8/X9: surface integration-overcorrection and annotation-coherence
        red-flags as a visible limitation block instead of buried numbers."""
        caveats = []
        for key, label in (("integration_qc", "Integration"),
                            ("annotation_qc", "Annotation")):
            qc = findings.get(key) or {}
            for issue in qc.get("issues", []) or []:
                sev = "warning" if issue.get("severity") != "blocking" else "blocking"
                caveats.append(Caveat(
                    f"{label}: {issue.get('message', '')} "
                    f"{issue.get('recommendation', '')}".strip(),
                    sev,
                ))
        if not caveats:
            return []
        return [NarrativeBlock(
            id="scrna.data_quality",
            modality="scRNA-seq",
            analysis="data_quality",
            block_type="limitation",
            title="Data quality & integration checks",
            status="warnings",
            confidence="medium",
            claim="Automated quality checks flagged issues to review before "
                  "interpreting downstream cell-type-resolved results.",
            caveats=caveats,
        )]

    def _error_blocks(self, findings: dict) -> list[NarrativeBlock]:
        blocks = []
        de = findings.get("differential_expression") or {}
        if de.get("status") and de.get("status") != "success":
            err = de.get("details") or de.get("reason") or de.get("error_type")
            blocks.append(NarrativeBlock(
                id="scrna.marker_discovery",
                modality="scRNA-seq",
                analysis="marker_discovery",
                block_type="error",
                title="Per-cluster marker discovery",
                status=de.get("status", "error"),
                confidence="insufficient",
                claim="",
                error=str(err or "marker discovery did not complete"),
                caveats=[Caveat(
                    "Per-cluster marker output is unavailable; completed "
                    "pseudobulk results remain separately interpretable.",
                    "warning",
                )],
            ))
        return blocks

    def _composition_blocks(self, findings: dict) -> list[NarrativeBlock]:
        da = findings.get("differential_abundance") or {}
        blocks = []
        for comp_key, comp in (da.get("per_comparison") or {}).items():
            status = comp.get("status", "success")
            if status != "success":
                blocks.append(NarrativeBlock(
                    id=f"scrna.composition.{_safe_id(comp_key)}",
                    modality="scRNA-seq",
                    analysis="differential_abundance",
                    block_type="error",
                    title=f"Cell-type abundance {comp_key}",
                    status=status,
                    confidence="insufficient",
                    claim="",
                    error=comp.get("reason", "differential abundance did not complete"),
                    caveats=[Caveat("Cell-type abundance could not be evaluated.")],
                ))
                continue
            n_sig = comp.get("n_significant", 0)
            rows = comp.get("per_cell_type", []) or []
            top = [r for r in rows if r.get("significant")][:5]
            evidence = [
                _evidence("cell types tested", len(rows), "differential_abundance"),
                _evidence("significant shifts", n_sig, "differential_abundance"),
            ]
            for row in top:
                evidence.append(_evidence(
                    f"{row.get('name')} abundance direction",
                    row.get("direction"),
                    "differential_abundance",
                ))
            blocks.append(NarrativeBlock(
                id=f"scrna.composition.{_safe_id(comp_key)}",
                modality="scRNA-seq",
                analysis="differential_abundance",
                block_type="result",
                title=f"Cell-type abundance {comp_key}",
                status="success",
                confidence="medium",
                claim=(
                    f"{n_sig} cell type(s) shifted in abundance for {comp_key}."
                ),
                evidence=evidence,
                metrics={"n_significant": n_sig},
            ))
        return blocks

    def _pseudobulk_blocks(self, findings: dict) -> list[NarrativeBlock]:
        pb = findings.get("pseudobulk_de") or {}
        blocks = []
        for group, info in (pb.get("per_group", {}) or {}).items():
            for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
                block_id = (
                    f"scrna.pseudobulk.{_safe_id(group)}.{_safe_id(comp_key)}"
                )
                if comp.get("status") != "success":
                    issue_text = _first_design_issue(comp)
                    blocks.append(NarrativeBlock(
                        id=block_id,
                        modality="scRNA-seq",
                        analysis="pseudobulk_de",
                        block_type="error",
                        title=f"{group} {comp_key}",
                        status=comp.get("status", "skipped"),
                        confidence="insufficient",
                        claim="",
                        error=issue_text or comp.get(
                            "reason", "pseudobulk block did not complete"
                        ),
                        caveats=[Caveat("No DE conclusion is drawn for this block.")],
                    ))
                    continue
                # Primary significance count follows the FDR strategy the script
                # applied. Legacy results without fdr_strategy keep global wording.
                strategy = comp.get("fdr_strategy", "global")
                fdr_label = ("per-cluster FDR" if strategy == "per_cluster"
                             else "global-FDR")
                n_sig = comp.get("n_significant", comp.get("n_significant_global", 0))
                if not n_sig:
                    continue
                evidence = [
                    _evidence("global-FDR DE genes",
                              comp.get("n_significant_global", n_sig), "pseudobulk_de"),
                    _evidence(
                        "local-FDR DE genes",
                        comp.get("n_significant_local", comp.get("n_significant", 0)),
                        "pseudobulk_de",
                    ),
                    _evidence("up genes", comp.get("n_up", comp.get("n_up_global", 0)),
                              "pseudobulk_de"),
                    _evidence(
                        "down genes",
                        comp.get("n_down", comp.get("n_down_global", 0)),
                        "pseudobulk_de",
                    ),
                ]
                for gene in (comp.get("top_genes") or [])[:6]:
                    evidence.append(_evidence(
                        f"top gene {gene.get('gene', '?')}",
                        gene.get("log2fc"),
                        "pseudobulk_de",
                    ))
                caveats = []
                # F-SCI-LOGNORM (audit 2026-05-28): counts recovered from
                # log-normalized values are quantitatively unreliable for the NB
                # dispersion estimate DESeq2 depends on. Surface the caveat AND
                # cap this block's confidence so the report never presents
                # recovered-count DE at the same trust level as raw-count DE.
                lognorm_recovered = bool(pb.get("lognorm_recovered"))
                if lognorm_recovered:
                    caveats.append(Caveat(
                        "Counts were reverse-engineered from log-normalized "
                        "values (raw counts were not available); DESeq2 inputs "
                        "are approximate integer reconstructions, so effect "
                        "sizes and dispersion are only directionally reliable.",
                        "warning",
                    ))
                if comp.get("low_power_warning"):
                    caveats.append(Caveat("Low replicate support; interpret cautiously."))
                for issue in _design_issues(comp, severities={"warning"}):
                    caveats.append(Caveat(
                        f"Design-matrix warning: {issue.get('message', '')}",
                        "warning",
                    ))
                if comp.get("corrected_for_composition"):
                    caveats.append(Caveat(
                        "DESeq2 design included a log-proportion composition covariate.",
                        "info",
                    ))
                else:
                    caveats.append(Caveat(
                        "No composition covariate was used for this DE block.",
                        "info",
                    ))
                blocks.append(NarrativeBlock(
                    id=block_id,
                    modality="scRNA-seq",
                    analysis="pseudobulk_de",
                    block_type="result",
                    title=f"{group} {comp_key}",
                    status="success",
                    confidence="low" if lognorm_recovered else "medium",
                    claim=f"{group} {comp_key} had {n_sig} {fdr_label} DE genes.",
                    evidence=evidence,
                    caveats=caveats,
                    metrics={
                        "n_significant_global": comp.get(
                            "n_significant_global", n_sig
                        ),
                        "n_significant_local": comp.get(
                            "n_significant_local", comp.get("n_significant", 0)
                        ),
                        "power_estimate_at_lfc_min": comp.get(
                            "power_estimate_at_lfc_min"
                        ),
                        "power_estimate_at_effective_alpha": comp.get(
                            "power_estimate_at_effective_alpha"
                        ),
                        "n_significant": n_sig,
                        "low_power_warning": bool(comp.get("low_power_warning")),
                        "lognorm_recovered": lognorm_recovered,
                        "corrected_for_composition": bool(
                            comp.get("corrected_for_composition")
                        ),
                    },
                ))
        return blocks

    def _pathway_blocks(self, findings: dict) -> list[NarrativeBlock]:
        pwp = findings.get("pseudobulk_pathways") or {}
        blocks = []
        for block_key, block in (pwp.get("per_cluster", {}) or {}).items():
            if block.get("status") and block.get("status") != "success":
                continue
            terms = []
            for db_name, rows in (block.get("results", {}) or {}).items():
                for row in rows or []:
                    term = row.get("term") or row.get("Term")
                    if term:
                        terms.append((db_name, term, row))
            if not terms:
                continue
            evidence = [
                _evidence("enriched terms", block.get("n_significant", len(terms)),
                          "pseudobulk_pathways")
            ]
            for db_name, term, row in terms[:5]:
                evidence.append(_evidence(
                    f"{db_name} term",
                    term,
                    "pseudobulk_pathways",
                ))
                if row.get("adjusted_p") is not None:
                    evidence.append(_evidence(
                        f"{term} FDR",
                        row.get("adjusted_p"),
                        "pseudobulk_pathways",
                    ))
            blocks.append(NarrativeBlock(
                id=f"scrna.pathway.{_safe_id(block_key)}",
                modality="scRNA-seq",
                analysis="pseudobulk_pathways",
                block_type="result",
                title=f"Pathway support {block_key}",
                status="success",
                confidence="medium",
                claim=(
                    f"ORA found {block.get('n_significant', len(terms))} "
                    f"enriched term(s) for {block_key}."
                ),
                evidence=evidence,
                caveats=[Caveat(
                    "ORA is an over-representation summary of DE genes, not "
                    "proof of pathway activity.",
                    "info",
                )],
                metrics={"n_significant": block.get("n_significant", len(terms))},
            ))
        return blocks

    def _cellcomm_blocks(self, findings: dict) -> list[NarrativeBlock]:
        ccc = findings.get("cell_communication") or {}
        if ccc.get("status") not in {"done", "success"}:
            return []
        evidence = [
            _evidence("interactions", ccc.get("n_interactions", 0),
                      "cell_communication"),
            _evidence("cell types", ccc.get("n_cell_types", 0),
                      "cell_communication"),
            _evidence("autocrine pairs excluded",
                      ccc.get("n_autocrine_dropped", 0),
                      "cell_communication"),
        ]
        for row in (ccc.get("top_interactions") or [])[:5]:
            evidence.append(_evidence(
                f"{row.get('source', '?')}->{row.get('target', '?')}",
                f"{row.get('ligand', '?')}-{row.get('receptor', '?')}",
                "cell_communication",
            ))
        return [NarrativeBlock(
            id="scrna.cellcomm",
            modality="scRNA-seq",
            analysis="cell_communication",
            block_type="exploratory",
            title="Cell-cell communication",
            status="success",
            confidence="low",
            claim=(
                f"LIANA reported {ccc.get('n_interactions', 0)} "
                "non-autocrine ligand-receptor interaction candidates."
            ),
            evidence=evidence,
            caveats=[Caveat(
                "Ligand-receptor scores are transcript-supported candidates "
                "requiring manual review.",
                "warning",
            )],
            metrics={"n_interactions": ccc.get("n_interactions", 0)},
        )]

    def _trajectory_blocks(self, findings: dict) -> list[NarrativeBlock]:
        traj = findings.get("trajectory") or {}
        if traj.get("status") not in {"done", "success"}:
            return []
        paga = traj.get("paga", {}) or {}
        pt = traj.get("pseudotime", {}) or {}
        vel = traj.get("velocity", {}) or {}
        evidence = [
            _evidence("PAGA connections", paga.get("n_connections", 0),
                      "trajectory"),
            _evidence("strong PAGA edges", paga.get("n_strong", 0),
                      "trajectory"),
            _evidence("DPT computed", bool(pt.get("computed")), "trajectory"),
        ]
        if pt.get("pseudotime_by_group"):
            ordered = sorted(pt["pseudotime_by_group"].items(), key=lambda kv: kv[1])
            evidence.append(_evidence(
                "DPT ordered groups",
                " -> ".join(str(g) for g, _ in ordered),
                "trajectory",
            ))
        return [NarrativeBlock(
            id="scrna.trajectory",
            modality="scRNA-seq",
            analysis="trajectory",
            block_type="exploratory",
            title="Trajectory context",
            status="success",
            confidence="low",
            claim=(
                f"PAGA/DPT evaluated {paga.get('n_connections', 0)} "
                "cluster-pair connections as exploratory manifold context."
            ),
            evidence=evidence,
            metadata={"velocity_computed": bool(vel.get("computed"))},
        )]

    def _attach_artifacts(self, blocks: list[NarrativeBlock],
                          findings: dict) -> None:
        figs = findings.get("figures") or {}
        tables = [
            {"id": key, "path": path, "label": key.replace("_", " ")}
            for key, path in (findings.get("tables") or {}).items()
        ]
        for block in blocks:
            if block.id == "scrna.qc":
                for key, path in figs.items():
                    if isinstance(path, str) and key.startswith("umap_"):
                        block.figures.append({
                            "id": key,
                            "path": path,
                            "caption": key.replace("_", " "),
                        })
            elif block.id.startswith("scrna.pseudobulk"):
                path = figs.get("per_celltype_de_bar")
                if path:
                    block.figures.append({
                        "id": "per_celltype_de_bar",
                        "path": path,
                        "caption": "Pseudobulk DE counts per group",
                    })
            elif block.id.startswith("scrna.pathway"):
                for block_key, paths in (figs.get("pathway_dotplots") or {}).items():
                    if _safe_id(block_key) in block.id:
                        for path in paths or []:
                            block.figures.append({
                                "id": f"pathway.{block_key}",
                                "path": path,
                                "caption": f"Pathway enrichment {block_key}",
                            })
            elif block.id == "scrna.cellcomm":
                for key in ("cellcomm_heatmap", "cellcomm_top_pairs"):
                    if figs.get(key):
                        block.figures.append({
                            "id": key,
                            "path": figs[key],
                            "caption": key.replace("_", " "),
                        })
            elif block.id == "scrna.trajectory":
                for key in ("paga_graph", "paga_log10_graph", "umap_dpt_pseudotime"):
                    if figs.get(key):
                        block.figures.append({
                            "id": key,
                            "path": figs[key],
                            "caption": key.replace("_", " "),
                        })
            block.tables.extend(tables)
