"""Bulk RNA narrative plugin."""

from __future__ import annotations

import re
from pathlib import Path

from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock


def _safe_id(value) -> str:
    text = str(value or "unknown")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "unknown"


def _evidence(label: str, value, source: str,
              path: str | None = None) -> EvidenceItem:
    return EvidenceItem(label=label, value=value, source=source, path=path)


class BulkRnaNarrator:
    name = "bulk_rna"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        return (
            agent_name == "bulk_rna_agent"
            and isinstance(agent_result, dict)
            and bool(agent_result.get("findings", agent_result))
        )

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        findings = agent_result.get("findings", {}) or {}
        blocks: list[NarrativeBlock] = []
        blocks.extend(self._qc_blocks(findings))
        blocks.extend(self._contrast_blocks(findings))
        blocks.extend(self._power_blocks(findings))
        self._attach_artifacts(blocks, findings)
        return blocks

    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        findings = agent_result.get("findings", {}) or {}
        contrasts = findings.get("contrasts", []) or []
        if not contrasts and not findings:
            return []
        design = findings.get("design_used", "~condition")
        padj = findings.get("padj_threshold", 0.05)
        lfc = findings.get("lfc_threshold", 1.0)
        return [
            "Bulk RNA-seq differential expression used DESeq2/pyDESeq2 "
            f"with design {design}; significance used adjusted p-value < "
            f"{padj} and |log2FC| > {lfc}."
        ]

    def figures(self, agent_name: str, agent_result: dict,
                report_dir: Path | None = None) -> list[dict]:
        findings = agent_result.get("findings", {}) or {}
        out = []
        sqc = findings.get("sample_qc", {}) or {}
        for key in ("pca_plot", "mds_plot"):
            if sqc.get(key):
                out.append({"id": key, "path": sqc[key], "caption": key})
        for contrast in findings.get("contrasts", []) or []:
            for key, value in (contrast.get("plots", {}) or {}).items():
                if isinstance(value, str):
                    out.append({
                        "id": f"{contrast.get('name', 'contrast')}.{key}",
                        "path": value,
                        "caption": f"{contrast.get('name')} {key}",
                    })
                elif isinstance(value, dict):
                    for sub_key, path in value.items():
                        out.append({
                            "id": f"{contrast.get('name', 'contrast')}.{key}.{sub_key}",
                            "path": path,
                            "caption": f"{contrast.get('name')} {sub_key}",
                        })
        return out

    def tables(self, agent_name: str, agent_result: dict,
               report_dir: Path | None = None) -> list[dict]:
        findings = agent_result.get("findings", {}) or {}
        tables = []
        for contrast in findings.get("contrasts", []) or []:
            for key, path in ((contrast.get("plots", {}) or {}).get("tables", {}) or {}).items():
                tables.append({
                    "id": f"{contrast.get('name', 'contrast')}.{key}",
                    "path": path,
                    "label": f"{contrast.get('name')} {key}",
                })
        return tables

    def _qc_blocks(self, findings: dict) -> list[NarrativeBlock]:
        sqc = findings.get("sample_qc", {}) or {}
        preprocessing = findings.get("preprocessing", {}) or {}
        if not sqc and not preprocessing:
            return []
        n_samples = sqc.get("n_samples") or (preprocessing.get("qc", {}) or {}).get(
            "n_samples"
        )
        outliers = sqc.get("outliers", []) or []
        evidence = []
        if n_samples is not None:
            evidence.append(_evidence("samples", n_samples, "sample_qc"))
        if sqc.get("size_ratio") is not None:
            evidence.append(_evidence("library-size range", sqc.get("size_ratio"),
                                      "sample_qc"))
        evidence.append(_evidence("outliers removed", len(outliers), "sample_qc"))
        claim = (
            f"Bulk RNA sample QC evaluated {n_samples} samples with "
            f"{len(outliers)} outlier(s) removed."
        )
        return [NarrativeBlock(
            id="bulk.qc",
            modality="bulk RNA-seq",
            analysis="sample_qc",
            block_type="qc",
            title="Bulk RNA sample QC",
            status="success",
            confidence="high",
            claim=claim,
            evidence=evidence or [_evidence("QC status", "completed", "sample_qc")],
            metrics={"n_samples": n_samples, "n_outliers": len(outliers)},
        )]

    def _contrast_blocks(self, findings: dict) -> list[NarrativeBlock]:
        blocks = []
        contrasts = findings.get("contrasts", []) or []
        for contrast in contrasts:
            name = contrast.get("name", "contrast")
            block_id = f"bulk.contrast.{_safe_id(name)}"
            status = contrast.get("status", "success")
            if status != "success":
                blocks.append(NarrativeBlock(
                    id=block_id,
                    modality="bulk RNA-seq",
                    analysis="differential_expression",
                    block_type="error",
                    title=f"Bulk contrast {name}",
                    status=status,
                    confidence="insufficient",
                    claim="",
                    error=contrast.get("reason", "contrast did not complete"),
                    caveats=[Caveat("No DE conclusion is drawn for this contrast.")],
                ))
                continue
            n_sig = contrast.get("n_significant", 0)
            evidence = [
                _evidence("DE genes", n_sig, "bulk_de"),
                _evidence("up genes", contrast.get("n_upregulated", 0), "bulk_de"),
                _evidence(
                    "down genes", contrast.get("n_downregulated", 0), "bulk_de"
                ),
            ]
            for gene in (contrast.get("top_genes") or [])[:6]:
                evidence.append(_evidence(
                    f"top gene {gene.get('symbol') or gene.get('gene', '?')}",
                    gene.get("log2fc"),
                    "bulk_de",
                ))
            caveats = []
            if contrast.get("low_power_warning"):
                caveats.append(Caveat(
                    contrast.get("low_power_reason")
                    or "Low replicate support; interpret cautiously."
                ))
            blocks.append(NarrativeBlock(
                id=block_id,
                modality="bulk RNA-seq",
                analysis="differential_expression",
                block_type="result",
                title=f"Bulk contrast {name}",
                status="success",
                confidence="medium",
                claim=f"Bulk contrast {name} had {n_sig} DE genes.",
                evidence=evidence,
                caveats=caveats,
                metrics={
                    "n_significant": n_sig,
                    "n_upregulated": contrast.get("n_upregulated", 0),
                    "n_downregulated": contrast.get("n_downregulated", 0),
                },
            ))
            blocks.extend(self._pathway_blocks_for_contrast(contrast))
        return blocks

    def _pathway_blocks_for_contrast(self, contrast: dict) -> list[NarrativeBlock]:
        blocks = []
        name = contrast.get("name", "contrast")
        pathways = contrast.get("pathways", {}) or {}
        terms = []
        for db_name, rows in pathways.items():
            for row in rows or []:
                term = row.get("term") or row.get("Term")
                if term:
                    terms.append((db_name, term, row))
        if not terms:
            return blocks
        evidence = [_evidence("enriched terms", len(terms), "bulk_pathways")]
        for db_name, term, row in terms[:5]:
            evidence.append(_evidence(f"{db_name} term", term, "bulk_pathways"))
            if row.get("adjusted_p") is not None:
                evidence.append(_evidence(f"{term} FDR", row.get("adjusted_p"),
                                          "bulk_pathways"))
        blocks.append(NarrativeBlock(
            id=f"bulk.pathway.{_safe_id(name)}",
            modality="bulk RNA-seq",
            analysis="pathway_enrichment",
            block_type="result",
            title=f"Pathway enrichment {name}",
            status="success",
            confidence="medium",
            claim=f"Pathway enrichment found {len(terms)} term(s) for {name}.",
            evidence=evidence,
            caveats=[Caveat(
                "ORA/GSEA terms summarize DE gene sets and are not causal "
                "pathway proof.",
                "info",
            )],
            metrics={"n_terms": len(terms)},
        ))
        return blocks

    def _power_blocks(self, findings: dict) -> list[NarrativeBlock]:
        powers = [
            c.get("power_estimate_at_lfc_min")
            for c in findings.get("contrasts", []) or []
            if c.get("status") == "success"
            and isinstance(c.get("power_estimate_at_lfc_min"), (int, float))
        ]
        if not powers:
            return []
        return [NarrativeBlock(
            id="bulk.power",
            modality="bulk RNA-seq",
            analysis="power",
            block_type="limitation",
            title="Bulk RNA power",
            status="success",
            confidence="medium",
            claim=(
                f"Approximate bulk RNA power ranged from {min(powers):.0%} "
                f"to {max(powers):.0%} across analyzable contrasts."
            ),
            evidence=[
                _evidence("minimum power", min(powers), "power_estimation"),
                _evidence("maximum power", max(powers), "power_estimation"),
            ],
            caveats=[Caveat(
                "Power is an approximation from replicate count, expression, "
                "and dispersion.",
                "info",
            )],
            metrics={"min_power": min(powers), "max_power": max(powers)},
        )]

    def _attach_artifacts(self, blocks: list[NarrativeBlock],
                          findings: dict) -> None:
        by_id = {block.id: block for block in blocks}
        sqc = findings.get("sample_qc", {}) or {}
        qc = by_id.get("bulk.qc")
        if qc:
            for key in ("pca_plot", "mds_plot"):
                if sqc.get(key):
                    qc.figures.append({
                        "id": key,
                        "path": sqc[key],
                        "caption": key.replace("_", " "),
                    })
        for contrast in findings.get("contrasts", []) or []:
            cid = f"bulk.contrast.{_safe_id(contrast.get('name', 'contrast'))}"
            block = by_id.get(cid)
            pathway = by_id.get(
                f"bulk.pathway.{_safe_id(contrast.get('name', 'contrast'))}"
            )
            plots = contrast.get("plots", {}) or {}
            if block:
                for key in ("volcano", "heatmap_padj", "heatmap_lfc", "heatmap"):
                    if plots.get(key):
                        block.figures.append({
                            "id": key,
                            "path": plots[key],
                            "caption": f"{contrast.get('name')} {key}",
                        })
                for key, path in (plots.get("tables", {}) or {}).items():
                    block.tables.append({
                        "id": key,
                        "path": path,
                        "label": f"{contrast.get('name')} {key}",
                    })
            if pathway:
                for db, path in (plots.get("ora_dotplots", {}) or {}).items():
                    pathway.figures.append({
                        "id": f"ora.{db}",
                        "path": path,
                        "caption": f"{contrast.get('name')} ORA {db}",
                    })
