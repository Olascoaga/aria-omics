"""Chromatin (scATAC) narrative plugin (v4.6 steps 5-6).

Composes the chromatin report section from validated ``NarrativeBlock`` objects
in this package — exactly like the scRNA / bulk RNA narrators — instead of
growing ``narrative_agent.py``. The narrator only fires on ``chromatin_agent``
results, so it is a no-op on the validated RNA paths.

It narrates the v4.6 scATAC peak-matrix pipeline:
  - QC (``chromatin_qc`` — measured metrics only; FRiP/TSS stay not-computed);
  - LSI dimensionality reduction + Leiden clustering (``chromatin_lsi_clustering``);
  - differential accessibility — per-cluster marker peaks + replicate-gated
    pseudobulk DA (``chromatin_diffacc``);
  - TF motif enrichment in DA peak sets (``chromatin_motifs``).

Honesty contract (ADR-002 / ADR-011): only measured quantities become evidence;
not-run / skipped lanes become honest limitation blocks with the concrete
reason, never fabricated results. Claims stay descriptive/associative — clusters
and peaks carry NO biological identity, and a JASPAR motif id/TF name is a
database fact, not a mechanistic claim about the cluster.
"""

from __future__ import annotations

import re
from pathlib import Path

from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock


def _safe_id(value) -> str:
    text = str(value or "chromatin")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "chromatin"


def _ev(label, value, source="chromatin"):
    return EvidenceItem(label=label, value=value, source=source)


def _ok(result) -> bool:
    return isinstance(result, dict) and str(result.get("status")) in (
        "success", "done")


class ChromatinNarrator:
    name = "chromatin"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        return (
            agent_name == "chromatin_agent"
            and isinstance(agent_result, dict)
            and bool(agent_result.get("findings", agent_result))
        )

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        findings = agent_result.get("findings", {}) or {}
        blocks: list[NarrativeBlock] = []

        qc = findings.get("qc")
        if isinstance(qc, dict) and _ok(qc):
            blocks.append(self._qc_block(qc))

        lsi = findings.get("lsi") or findings.get("lsi_clustering")
        if isinstance(lsi, dict) and _ok(lsi):
            blocks.append(self._clustering_block(lsi))

        da = findings.get("differential_accessibility")
        if isinstance(da, dict) and _ok(da):
            blocks.extend(self._diffacc_blocks(da))

        motifs = findings.get("motifs")
        if isinstance(motifs, dict) and _ok(motifs):
            blocks.append(self._motif_block(motifs))

        return blocks

    # ── QC ────────────────────────────────────────────────────────────────
    def _qc_block(self, qc: dict) -> NarrativeBlock:
        evidence: list[EvidenceItem] = []
        for label, key in (
            ("Cells (barcodes)", "n_cells"),
            ("Peaks", "n_peaks"),
            ("Fragments scanned", "n_fragments"),
            ("Mitochondrial fraction", "mito_fraction"),
        ):
            val = qc.get(key)
            if val is not None:
                evidence.append(_ev(label, val, "chromatin_qc"))

        caveats = [Caveat(
            "scATAC QC is a v4.6 scaffold: it reports only the metrics actually "
            "measured from the input; FRiP and TSS enrichment stay not-computed "
            "until called peaks and a reference TSS annotation are available.",
            severity="info",
        )]
        not_computed = qc.get("metrics_not_computed") or []
        if not_computed:
            caveats.append(Caveat(
                "Not computed: " + ", ".join(str(m) for m in not_computed),
                severity="info"))

        status = "success" if evidence else "limitation"
        claim = (
            f"scATAC QC measured {qc.get('n_cells', '?')} cells and "
            f"{qc.get('n_peaks', '?')} peaks."
            if evidence else "")
        return NarrativeBlock(
            id=f"chromatin.qc.{_safe_id(qc.get('data_type'))}",
            modality="chromatin", analysis="qc", block_type="qc",
            title="Chromatin QC",
            status=status, confidence="medium" if evidence else "low",
            claim=claim, evidence=evidence, caveats=caveats,
            metrics={"qc_complete": qc.get("qc_complete"),
                     "pass_qc": qc.get("pass_qc")},
            metadata={"validation_level": "scaffold"},
        )

    # ── LSI + clustering ────────────────────────────────────────────────────
    def _clustering_block(self, lsi: dict) -> NarrativeBlock:
        n_clusters = lsi.get("n_clusters")
        n_cells = lsi.get("n_cells_used")
        n_peaks = lsi.get("n_peaks")
        dropped = lsi.get("dropped_components") or []
        n_used = lsi.get("n_components_used")

        evidence = [
            _ev("Method", "TF-IDF/LSI", "chromatin_lsi_clustering"),
            _ev("Clusters", n_clusters, "chromatin_lsi_clustering"),
            _ev("Cells", n_cells, "chromatin_lsi_clustering"),
            _ev("Peaks", n_peaks, "chromatin_lsi_clustering"),
            _ev("LSI components used", n_used, "chromatin_lsi_clustering"),
            _ev("Depth-associated components removed", len(dropped),
                "chromatin_lsi_clustering"),
        ]
        evidence = [e for e in evidence if e.value is not None]

        caveats = [Caveat(
            "Clusters are accessibility-defined groupings only; ARIA assigns no "
            "cell-type identity to them.", severity="info")]
        if lsi.get("sketch_used"):
            caveats.append(Caveat(
                "A random cell sketch was used before LSI; cluster assignments "
                "cover the sketch.", severity="warning"))

        claim = (
            f"TF-IDF/LSI over {n_peaks} peaks in {n_cells} cells produced "
            f"{n_clusters} clusters after removing {len(dropped)} "
            f"depth-associated component(s).")
        return NarrativeBlock(
            id="chromatin.clustering", modality="chromatin",
            analysis="dimensionality_reduction", block_type="result",
            title="Chromatin LSI clustering",
            status="success", confidence="medium",
            claim=claim, evidence=evidence, caveats=caveats,
            metrics={"resolution": lsi.get("resolution"),
                     "rep_used": lsi.get("rep_used"),
                     "cluster_sizes": lsi.get("cluster_sizes")},
        )

    # ── Differential accessibility ──────────────────────────────────────────
    def _diffacc_blocks(self, da: dict) -> list[NarrativeBlock]:
        blocks: list[NarrativeBlock] = []
        padj_max = da.get("padj_max")
        lfc_min = da.get("lfc_min")

        pc = da.get("per_cluster") or {}
        if pc.get("ran"):
            n_total = pc.get("n_da_total", 0)
            n_with = sum(1 for v in (pc.get("n_da_by_cluster") or {}).values()
                         if v)
            evidence = [
                _ev("Differentially accessible peaks", n_total,
                    "chromatin_diffacc"),
                _ev("Clusters with DA peaks", n_with, "chromatin_diffacc"),
                _ev("padj threshold", padj_max, "chromatin_diffacc"),
                _ev("|log2FC| threshold", lfc_min, "chromatin_diffacc"),
            ]
            evidence = [e for e in evidence if e.value is not None]
            if pc.get("output_csv"):
                evidence.append(EvidenceItem(
                    label="Per-cluster DA table", value="csv",
                    source="chromatin_diffacc", path=pc["output_csv"]))
            claim = (
                f"Per-cluster accessibility testing identified {n_total} "
                f"differentially accessible peaks across {n_with} cluster(s) "
                f"at padj<{padj_max}, |log2FC|>{lfc_min}.")
            blocks.append(NarrativeBlock(
                id="chromatin.differential_accessibility.per_cluster",
                modality="chromatin", analysis="differential_accessibility",
                block_type="result", title="Per-cluster differential accessibility",
                status="success", confidence="medium", claim=claim,
                evidence=evidence,
                caveats=[Caveat(
                    "Marker peaks distinguish clusters by accessibility; both "
                    "the effect size and the in/out detection fraction are "
                    "reported. No transcription-factor or gene identity is "
                    "inferred from peak coordinates.", severity="info")],
                metrics={"n_da_by_cluster": pc.get("n_da_by_cluster")},
            ))

        pb = da.get("pseudobulk") or {}
        if pb.get("ran"):
            comps = [c for c in (pb.get("comparisons") or [])
                     if c.get("status") == "success"]
            total_sig = sum(int(c.get("n_sig", 0)) for c in comps)
            evidence = [
                _ev("Method", "pseudobulk DESeq2", "chromatin_diffacc"),
                _ev("Pseudobulk comparisons run", len(comps),
                    "chromatin_diffacc"),
                _ev("Differentially accessible peaks (pseudobulk)", total_sig,
                    "chromatin_diffacc"),
            ]
            claim = (
                f"Replicate-level pseudobulk DA ran {len(comps)} comparison(s) "
                f"and found {total_sig} differentially accessible peaks via the "
                f"shared DESeq2 core.")
            blocks.append(NarrativeBlock(
                id="chromatin.differential_accessibility.pseudobulk",
                modality="chromatin", analysis="differential_accessibility",
                block_type="result",
                title="Pseudobulk differential accessibility",
                status="success", confidence="medium", claim=claim,
                evidence=evidence,
                caveats=[Caveat(
                    "Pseudobulk DA aggregates peak counts per biological "
                    "replicate; single cells are not treated as replicates.",
                    severity="info")],
                metrics={"comparisons": pb.get("comparisons")},
            ))
        elif pb.get("reason"):
            blocks.append(NarrativeBlock(
                id="chromatin.differential_accessibility.pseudobulk_skipped",
                modality="chromatin", analysis="differential_accessibility",
                block_type="limitation",
                title="Pseudobulk differential accessibility not run",
                status="limitation", confidence="insufficient", claim="",
                evidence=[],
                caveats=[Caveat(
                    f"Cross-condition pseudobulk DA was not run: "
                    f"{pb.get('reason')}.", severity="info")],
            ))

        return blocks

    # ── Motif enrichment ────────────────────────────────────────────────────
    def _motif_block(self, motifs: dict) -> NarrativeBlock:
        if not motifs.get("ran"):
            return NarrativeBlock(
                id="chromatin.motifs.skipped", modality="chromatin",
                analysis="motif_enrichment", block_type="limitation",
                title="TF motif enrichment not run",
                status="limitation", confidence="insufficient", claim="",
                evidence=[],
                caveats=[Caveat(
                    f"Motif enrichment was not run: "
                    f"{motifs.get('reason', 'prerequisites missing')}.",
                    severity="info")],
            )

        src = motifs.get("motif_source") or {}
        per_group = motifs.get("per_group") or {}
        n_groups = len(per_group)
        n_with = sum(1 for g in per_group.values()
                     if (g or {}).get("n_enriched"))
        n_motifs = src.get("n_motifs")
        collection = src.get("collection", "motif collection")
        release = src.get("release")
        method = motifs.get("method", "enrichment")

        evidence = [
            _ev("Analysis", "TF motif enrichment", "chromatin_motifs"),
            _ev("Method", str(method), "chromatin_motifs"),
            _ev("Motif collection", str(collection), "chromatin_motifs"),
            _ev("Motifs tested", n_motifs, "chromatin_motifs"),
            _ev("Peak groups tested", n_groups, "chromatin_motifs"),
            _ev("Groups with enriched motifs", n_with, "chromatin_motifs"),
        ]
        if release is not None:
            evidence.append(_ev("Motif release", str(release), "chromatin_motifs"))
        evidence = [e for e in evidence if e.value is not None]
        if motifs.get("output_csv"):
            evidence.append(EvidenceItem(
                label="Motif enrichment table", value="csv",
                source="chromatin_motifs", path=motifs["output_csv"]))

        claim = (
            f"TF motif enrichment ({collection}, {method}) tested {n_motifs} "
            f"motifs and found enriched motifs in {n_with} of {n_groups} "
            f"peak group(s).")
        return NarrativeBlock(
            id="chromatin.motifs", modality="chromatin",
            analysis="motif_enrichment", block_type="result",
            title="TF motif enrichment", status="success",
            confidence="medium", claim=claim, evidence=evidence,
            caveats=[Caveat(
                "Motif enrichment is association-only: an enriched JASPAR motif "
                "(matrix id + TF name) is a database match in the peak set, not "
                "evidence that the factor is active or regulates a gene. "
                "Per-cell motif activity (chromVAR-style) is out of scope.",
                severity="info")],
            metrics={"genome_fasta": motifs.get("genome_fasta"),
                     "motif_source": src},
        )

    # ── Methods + tables ────────────────────────────────────────────────────
    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        findings = agent_result.get("findings", {}) or {}
        out: list[str] = []
        if _ok(findings.get("qc")):
            out.append(
                "Chromatin QC reports only measured fragment/barcode/peak "
                "metrics; TSS enrichment and FRiP remain uncomputed until "
                "called peaks and a reference TSS annotation are available.")
        lsi = findings.get("lsi") or findings.get("lsi_clustering")
        if _ok(lsi):
            out.append(
                f"scATAC dimensionality reduction: TF-IDF (log-TF-log-IDF) on "
                f"the peak matrix, truncated SVD/LSI "
                f"({lsi.get('n_components_computed', '?')} components), removal "
                f"of LSI components correlated with sequencing depth, then "
                f"neighbors/UMAP/Leiden (resolution {lsi.get('resolution','?')}).")
        if _ok(findings.get("differential_accessibility")):
            out.append(
                "Differential accessibility: per-cluster Wilcoxon marker peaks "
                "(cluster vs rest) with in/out accessibility fractions, plus "
                "replicate-level pseudobulk testing through the shared DESeq2 "
                "core when condition + replicate metadata and an explicit "
                "comparison are present.")
        motifs = findings.get("motifs")
        if _ok(motifs) and motifs.get("ran"):
            src = motifs.get("motif_source") or {}
            out.append(
                f"TF motif enrichment: {motifs.get('method','enrichment')} test "
                f"of DA peak sets against background using the versioned local "
                f"{src.get('collection','motif')} collection "
                f"(release {src.get('release','?')}); offline, no network egress.")
        return out

    def figures(self, agent_name: str, agent_result: dict,
                report_dir: Path | None = None) -> list[dict]:
        return []

    def tables(self, agent_name: str, agent_result: dict,
               report_dir: Path | None = None) -> list[dict]:
        findings = agent_result.get("findings", {}) or {}
        tables = []
        da = findings.get("differential_accessibility") or {}
        pc = (da.get("per_cluster") or {})
        if pc.get("output_csv"):
            tables.append({"id": "chromatin.da_per_cluster",
                           "path": pc["output_csv"],
                           "label": "Per-cluster differential accessibility"})
        pb = (da.get("pseudobulk") or {})
        if pb.get("output_csv"):
            tables.append({"id": "chromatin.da_pseudobulk",
                           "path": pb["output_csv"],
                           "label": "Pseudobulk differential accessibility"})
        motifs = findings.get("motifs") or {}
        if motifs.get("output_csv"):
            tables.append({"id": "chromatin.motif_enrichment",
                           "path": motifs["output_csv"],
                           "label": "TF motif enrichment"})
        return tables
