"""Chromatin narrative plugin (scATAC + beta bulk ATAC slices).

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

It also narrates the V47 bulk ATAC beta slice:
  - measured bulk chromatin QC;
  - MACS3 peak calling.
  - scaffolded peak-count matrix output when a comparison is requested;
  - replicate-gated DESeq2 DA when explicit condition/replicate/comparison
    metadata are present.

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


def _modality_label(data_type: str | None) -> str:
    labels = {
        "scATAC": "scATAC",
        "bulk_ATAC": "bulk ATAC",
        "ChIP": "ChIP-seq",
        "CUT_AND_RUN": "CUT&RUN",
        "CUT_AND_TAG": "CUT&TAG",
    }
    return labels.get(str(data_type or ""), "chromatin")


def _is_bulk_atac(findings: dict) -> bool:
    """True when the (unwrapped) findings come from the bulk ATAC lane.

    Bulk ATAC and scATAC are separate dispatch lanes, so any one
    ``data_type == "bulk_ATAC"`` stamp on QC/peaks/counts/DA is decisive.
    """
    for key in ("qc", "peaks", "peak_counts", "differential_accessibility"):
        sub = findings.get(key)
        if isinstance(sub, dict) and sub.get("data_type") == "bulk_ATAC":
            return True
    return False


# ChromatinAgent.run() returns findings keyed by modality, with the analysis
# results nested one level deeper, e.g.
#   {"findings": {"scATAC": {"status": "done",
#                            "findings": {"qc": ..., "lsi": ..., ...}}}}.
# Both the narrator and the run-ledger reconcile against the FLAT analysis keys
# (qc/lsi/differential_accessibility/motifs), so unwrap the per-modality wrapper
# first. A flat findings dict (used by older tests / a direct script feed) is
# returned unchanged.
_CHROMATIN_MODALITY_KEYS = (
    "scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG",
)


def unwrap_chromatin_findings(agent_result: dict) -> dict:
    """Flatten ChromatinAgent's per-modality finding wrapper to analysis keys."""
    if not isinstance(agent_result, dict):
        return {}
    findings = agent_result.get("findings", agent_result)
    if not isinstance(findings, dict):
        return {}
    modality_subs = {
        k: v for k, v in findings.items()
        if k in _CHROMATIN_MODALITY_KEYS and isinstance(v, dict)
    }
    if not modality_subs:
        return findings
    out: dict = {}
    for sub in modality_subs.values():
        nested = sub.get("findings", sub)
        if isinstance(nested, dict):
            out.update(nested)
    return out


class ChromatinNarrator:
    name = "chromatin"

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        return (
            agent_name == "chromatin_agent"
            and isinstance(agent_result, dict)
            and bool(unwrap_chromatin_findings(agent_result))
        )

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        findings = unwrap_chromatin_findings(agent_result)
        blocks: list[NarrativeBlock] = []

        qc = findings.get("qc")
        if isinstance(qc, dict) and _ok(qc):
            blocks.append(self._qc_block(qc))

        peaks = findings.get("peaks")
        if isinstance(peaks, dict) and _ok(peaks):
            blocks.append(self._peaks_block(peaks))

        peak_counts = findings.get("peak_counts")
        if isinstance(peak_counts, dict) and _ok(peak_counts):
            blocks.append(self._peak_counts_block(peak_counts))

        lsi = findings.get("lsi") or findings.get("lsi_clustering")
        if isinstance(lsi, dict) and _ok(lsi):
            blocks.append(self._clustering_block(lsi))

        da = findings.get("differential_accessibility")
        if isinstance(da, dict) and _ok(da):
            blocks.extend(self._diffacc_blocks(da))

        motifs = findings.get("motifs")
        if isinstance(motifs, dict) and _ok(motifs):
            blocks.append(self._motif_block(motifs))

        regulatory = findings.get("regulatory")
        if isinstance(regulatory, dict) and _ok(regulatory):
            blocks.append(self._regulatory_block(regulatory))

        footprinting = findings.get("footprinting")
        if (isinstance(footprinting, dict) and footprinting.get("ran")
                and isinstance(footprinting.get("differential_summary"), dict)
                and footprinting["differential_summary"].get("parsed")):
            blocks.append(self._footprint_block(footprinting))

        return blocks

    # ── QC ────────────────────────────────────────────────────────────────
    def _qc_block(self, qc: dict) -> NarrativeBlock:
        data_type = qc.get("data_type")
        modality_label = _modality_label(data_type)
        evidence: list[EvidenceItem] = []
        for ev_label, key in (
            ("Modality", "data_type"),
            ("Cells (barcodes)", "n_cells"),
            ("Samples", "n_samples"),
            ("Peaks", "n_peaks"),
            ("Fragments scanned", "n_fragments"),
            ("Mitochondrial fraction", "mito_fraction"),
            ("Duplicate rate", "dup_rate"),
        ):
            val = qc.get(key)
            if val is not None:
                evidence.append(_ev(ev_label, val, "chromatin_qc"))

        caveats = [Caveat(
            "Chromatin QC reports only the metrics actually measured from the "
            "input; metrics that require additional assets or post-peak-calling "
            "steps stay not-computed, and scaffolded layers remain explicitly "
            "caveated.",
            severity="info",
        )]
        not_computed = qc.get("metrics_not_computed") or []
        if not_computed:
            caveats.append(Caveat(
                "Not computed: " + ", ".join(str(m) for m in not_computed),
                severity="info"))

        status = "success" if evidence else "limitation"
        if data_type == "bulk_ATAC" and evidence:
            claim = (
                f"Bulk ATAC QC measured {qc.get('n_samples', '?')} sample(s) "
                f"with mean mitochondrial fraction {qc.get('mito_fraction', '?')}.")
        elif evidence:
            claim = (
                f"{modality_label} QC measured {qc.get('n_cells', '?')} cells and "
                f"{qc.get('n_peaks', '?')} peaks.")
        else:
            claim = ""
        return NarrativeBlock(
            id=f"chromatin.qc.{_safe_id(qc.get('data_type'))}",
            modality="chromatin", analysis="qc", block_type="qc",
            title="Chromatin QC",
            status=status, confidence="medium" if evidence else "low",
            claim=claim, evidence=evidence, caveats=caveats,
            metrics={"qc_complete": qc.get("qc_complete"),
                     "pass_qc": qc.get("pass_qc")},
            metadata={"validation_level": "beta"
                      if data_type == "bulk_ATAC" else "scaffold"},
        )

    # ── Peak calling ──────────────────────────────────────────────────────
    def _peaks_block(self, peaks: dict) -> NarrativeBlock:
        data_type = peaks.get("data_type")
        label = _modality_label(data_type)
        n_peaks = peaks.get("n_peaks")
        frip = peaks.get("frip")
        evidence = [
            _ev("Analysis", "MACS3 peak calling", "chromatin_peaks"),
            _ev("Modality", data_type, "chromatin_peaks"),
            _ev("Peaks called", n_peaks, "chromatin_peaks"),
            _ev("Genome", peaks.get("genome"), "chromatin_peaks"),
        ]
        if frip is not None:
            evidence.append(_ev("FRiP", frip, "chromatin_peaks"))
        if peaks.get("peaks_path"):
            evidence.append(EvidenceItem(
                label="Peak table", value="peak file",
                source="chromatin_peaks", path=peaks["peaks_path"]))
        if peaks.get("consensus_peaks_path"):
            evidence.append(EvidenceItem(
                label="Consensus peaks", value="peak file",
                source="chromatin_peaks", path=peaks["consensus_peaks_path"]))
        evidence = [e for e in evidence if e.value is not None]

        caveats = [Caveat(
            "Peak calling identifies accessible genomic intervals. These peaks "
            "are not differential accessibility results and do not establish "
            "regulatory mechanism on their own.",
            severity="info",
        )]
        if frip is None:
            caveats.append(Caveat(
                "FRiP was not computed after peak calling; no default FRiP value "
                "was substituted.",
                severity="info"))
        for warning in peaks.get("warnings") or []:
            caveats.append(Caveat(str(warning), severity="warning"))

        claim = f"MACS3 peak calling for {label} identified {n_peaks} peaks."
        return NarrativeBlock(
            id=f"chromatin.peak_calling.{_safe_id(data_type)}",
            modality="chromatin", analysis="peak_calling",
            block_type="result", title="Chromatin peak calling",
            status="success", confidence="medium",
            claim=claim, evidence=evidence, caveats=caveats,
            metrics={
                "n_peaks": n_peaks,
                "frip": frip,
                "consensus_peaks_path": peaks.get("consensus_peaks_path"),
            },
            metadata={"validation_level": "beta"
                      if data_type == "bulk_ATAC" else "scaffold"},
        )

    def _peak_counts_block(self, counts: dict) -> NarrativeBlock:
        data_type = counts.get("data_type")
        n_peaks = counts.get("n_peaks")
        n_samples = counts.get("n_samples")
        evidence = [
            _ev("Analysis", "Peak-count matrix", "chromatin_peak_counts"),
            _ev("Method", counts.get("counting_method"),
                "chromatin_peak_counts"),
            _ev("Modality", data_type, "chromatin_peak_counts"),
            _ev("Peaks counted", n_peaks, "chromatin_peak_counts"),
            _ev("Samples", n_samples, "chromatin_peak_counts"),
        ]
        if counts.get("counts_matrix_path"):
            evidence.append(EvidenceItem(
                label="Count matrix", value="peak x sample TSV",
                source="chromatin_peak_counts",
                path=counts["counts_matrix_path"]))
        if counts.get("sample_metadata_path"):
            evidence.append(EvidenceItem(
                label="Sample metadata", value="sample TSV",
                source="chromatin_peak_counts",
                path=counts["sample_metadata_path"]))
        evidence = [e for e in evidence if e.value is not None]

        caveats = [Caveat(
            "The peak-count matrix is a technical artifact for bulk ATAC. It is "
            "not a differential accessibility result; replicate-aware DA "
            "remains scaffolded until separately validated.",
            severity="warning",
        )]
        for warning in counts.get("warnings") or []:
            caveats.append(Caveat(str(warning), severity="warning"))

        claim = (
            f"Bulk ATAC peak counting produced a matrix with {n_peaks} peaks "
            f"across {n_samples} samples."
        )
        return NarrativeBlock(
            id=f"chromatin.peak_count_matrix.{_safe_id(data_type)}",
            modality="chromatin", analysis="peak_count_matrix",
            block_type="result", title="Bulk ATAC peak-count matrix",
            status="success", confidence="low",
            claim=claim, evidence=evidence, caveats=caveats,
            metrics={
                "n_peaks": n_peaks,
                "n_samples": n_samples,
                "counts_matrix_path": counts.get("counts_matrix_path"),
            },
            metadata={"validation_level": "scaffold"},
        )

    # ── LSI + clustering ────────────────────────────────────────────────────
    def _clustering_block(self, lsi: dict) -> NarrativeBlock:
        n_clusters = lsi.get("n_clusters")
        n_cells = lsi.get("n_cells_used")
        n_peaks = lsi.get("n_peaks")
        dropped = lsi.get("dropped_components") or []
        n_used = lsi.get("n_components_used")
        doublets = lsi.get("doublets") or {}
        batch_qc = lsi.get("batch_qc") or {}
        consensus = lsi.get("consensus_peaks") or {}

        evidence = [
            _ev("Method", "TF-IDF/LSI", "chromatin_lsi_clustering"),
            _ev("Clusters", n_clusters, "chromatin_lsi_clustering"),
            _ev("Cells", n_cells, "chromatin_lsi_clustering"),
            _ev("Peaks", n_peaks, "chromatin_lsi_clustering"),
            _ev("LSI components used", n_used, "chromatin_lsi_clustering"),
            _ev("Depth-associated components removed", len(dropped),
                "chromatin_lsi_clustering"),
        ]
        if doublets:
            evidence.extend([
                _ev("Doublet detector", doublets.get("method"),
                    "chromatin_lsi_clustering"),
                _ev("Predicted doublets removed", doublets.get("removed", 0),
                    "chromatin_lsi_clustering"),
            ])
        if consensus:
            evidence.append(_ev(
                "Consensus peak provenance", consensus.get("status"),
                "chromatin_lsi_clustering"))
        evidence = [e for e in evidence if e.value is not None]

        caveats = [Caveat(
            "Clusters are accessibility-defined groupings only; ARIA assigns no "
            "cell-type identity to them.", severity="info")]
        if lsi.get("sketch_used"):
            caveats.append(Caveat(
                "A random cell sketch was used before LSI; cluster assignments "
                "cover the sketch.", severity="warning"))
        if doublets and not doublets.get("ran"):
            caveats.append(Caveat(
                "scATAC doublet detection was not run: "
                f"{doublets.get('reason', 'prerequisites missing')}.",
                severity="warning"))
        if batch_qc.get("issues"):
            checks = ", ".join(
                str(i.get("check")) for i in batch_qc.get("issues", [])
            )
            caveats.append(Caveat(
                f"scATAC batch QC raised warning(s): {checks}.",
                severity="warning"))
        if consensus and consensus.get("status") != "verified":
            caveats.append(Caveat(
                "Consensus peak provenance is "
                f"{consensus.get('status')}; peak reproducibility and rare-peak "
                "preservation are not fully verified from this matrix alone.",
                severity="warning"))

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
                     "cluster_sizes": lsi.get("cluster_sizes"),
                     "doublets": doublets,
                     "batch_qc": batch_qc,
                     "consensus_peaks": consensus},
        )

    # ── Differential accessibility ──────────────────────────────────────────
    def _diffacc_blocks(self, da: dict) -> list[NarrativeBlock]:
        if da.get("data_type") == "bulk_ATAC" and da.get("ran"):
            return [self._bulk_atac_diffacc_block(da)]

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

    def _bulk_atac_diffacc_block(self, da: dict) -> NarrativeBlock:
        n_comps = da.get("n_comparisons_success")
        n_sig = da.get("n_sig_total")
        n_reps = da.get("n_replicate_samples")
        n_peaks = da.get("n_peaks_tested")
        evidence = [
            _ev("Analysis", "Bulk ATAC differential accessibility",
                "chromatin_bulk_diffacc"),
            _ev("Method", da.get("method"), "chromatin_bulk_diffacc"),
            _ev("Modality", da.get("data_type"), "chromatin_bulk_diffacc"),
            _ev("Comparisons run", n_comps, "chromatin_bulk_diffacc"),
            _ev("Replicate samples", n_reps, "chromatin_bulk_diffacc"),
            _ev("Peaks tested", n_peaks, "chromatin_bulk_diffacc"),
            _ev("Differentially accessible peaks", n_sig,
                "chromatin_bulk_diffacc"),
            _ev("padj threshold", da.get("padj_max"),
                "chromatin_bulk_diffacc"),
            _ev("|log2FC| threshold", da.get("lfc_min"),
                "chromatin_bulk_diffacc"),
        ]
        if da.get("output_csv"):
            evidence.append(EvidenceItem(
                label="Bulk ATAC DA summary", value="csv",
                source="chromatin_bulk_diffacc", path=da["output_csv"]))
        if da.get("replicate_counts_path"):
            evidence.append(EvidenceItem(
                label="Replicate peak counts", value="peak x replicate TSV",
                source="chromatin_bulk_diffacc",
                path=da["replicate_counts_path"]))
        evidence = [e for e in evidence if e.value is not None]

        caveats = [Caveat(
            "Bulk ATAC DA uses biological replicate-level peak counts and the "
            "shared DESeq2 core. It requires explicit condition, replicate, and "
            "comparison metadata; ARIA does not infer contrasts from filenames.",
            severity="info",
        )]
        for warning in da.get("warnings") or []:
            caveats.append(Caveat(str(warning), severity="warning"))

        claim = (
            f"Bulk ATAC DESeq2 differential accessibility ran {n_comps} "
            f"comparison(s) over {n_reps} replicate samples and found {n_sig} "
            f"differentially accessible peaks.")
        return NarrativeBlock(
            id="chromatin.differential_accessibility.bulk_atac",
            modality="chromatin", analysis="differential_accessibility",
            block_type="result",
            title="Bulk ATAC differential accessibility",
            status="success", confidence="medium",
            claim=claim, evidence=evidence, caveats=caveats,
            metrics={
                "comparisons": da.get("comparisons"),
                "n_peaks_tested": n_peaks,
                "n_replicate_samples": n_reps,
                "n_sig_total": n_sig,
            },
            metadata={"validation_level": da.get("validation_level", "beta")},
        )

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

    # ── P2 regulatory layers ────────────────────────────────────────────────
    def _regulatory_block(self, regulatory: dict) -> NarrativeBlock:
        layers = {
            "motif_activity": "Motif activity",
            "gene_scores": "Gene activity scores",
            "footprinting": "Tn5 footprinting",
            "peak_to_gene": "Peak-to-gene links",
            "label_transfer": "scRNA label transfer",
        }
        ran = [key for key in layers if (regulatory.get(key) or {}).get("ran")]
        skipped = [key for key in layers if key not in ran]

        evidence = [
            _ev("Analysis", "scATAC regulatory layers", "chromatin_regulatory"),
            _ev("Optional layers", len(layers), "chromatin_regulatory"),
            _ev("Layers run", len(ran), "chromatin_regulatory"),
            _ev("Layers skipped", len(skipped), "chromatin_regulatory"),
        ]
        for key in ran:
            sub = regulatory.get(key) or {}
            label = layers[key]
            if key == "motif_activity":
                evidence.append(_ev(label, sub.get("n_motifs"),
                                    "chromatin_regulatory"))
            elif key == "gene_scores":
                evidence.append(_ev(label, sub.get("n_genes_scored"),
                                    "chromatin_regulatory"))
            elif key == "peak_to_gene":
                evidence.append(_ev(label, sub.get("n_links"),
                                    "chromatin_regulatory"))
            elif key == "label_transfer":
                evidence.append(_ev(label, sub.get("n_labeled_cells"),
                                    "chromatin_regulatory"))
            else:
                evidence.append(_ev(label, sub.get("method"),
                                    "chromatin_regulatory"))
        evidence = [e for e in evidence if e.value is not None]

        caveats = [
            Caveat(
                "Regulatory layers are associative. Motif activity, gene scores, "
                "footprints, peak-to-gene links, and transferred labels do not "
                "establish causal regulation.",
                severity="info",
            ),
            Caveat(
                "Transferred scRNA labels are report hypotheses only and are "
                "not trusted as inferential groupby columns.",
                severity="info",
            ),
        ]
        for key in skipped:
            sub = regulatory.get(key) or {}
            caveats.append(Caveat(
                f"{layers[key]} not run: "
                f"{sub.get('reason', 'prerequisites missing')}.",
                severity="info"))

        if (regulatory.get("peak_to_gene") or {}).get("validation_level") == "beta":
            caveats.append(Caveat(
                "Peak-to-gene link recovery is beta-grade: externally validated "
                "for concordance against a canonical single-cell peak-gene linker "
                "(see benchmark artifact). The other regulatory layers remain "
                "scaffold/exploratory. Beta refers to link recovery, not causal "
                "regulation — the links stay associative.",
                severity="info"))

        claim = (
            f"scATAC regulatory-layer analysis ran {len(ran)} of "
            f"{len(layers)} optional layer(s); skipped layers are reported with "
            f"their prerequisites.")
        status = "success" if ran else "limitation"
        block_type = "result" if ran else "limitation"
        return NarrativeBlock(
            id="chromatin.regulatory_layers", modality="chromatin",
            analysis="regulatory_layers", block_type=block_type,
            title="scATAC regulatory layers", status=status,
            confidence="low" if ran else "insufficient",
            claim=claim if ran else "", evidence=evidence, caveats=caveats,
            metrics={key: regulatory.get(key) for key in layers},
        )

    # ── P4.3 differential TF footprinting (governed, cross-modal) ────────────
    def _footprint_block(self, footprinting: dict) -> NarrativeBlock:
        """Governed associative block for TOBIAS differential TF binding, with optional
        footprint<->RNA cross-evidence. Claim carries COUNTS only (no TF names, so the
        W-CLAIM named-entity guard is satisfied); specific TFs live in evidence/metrics.
        Capped associative: a differential footprint is an occupancy difference, not
        regulation."""
        summary = footprinting.get("differential_summary") or {}
        cross = footprinting.get("rna_cross_evidence") or {}
        ga = footprinting.get("group_a", "group A")
        gb = footprinting.get("group_b", "group B")
        n_sig = summary.get("n_significant")
        n_tested = summary.get("n_motifs_tested")

        evidence = [
            _ev("Analysis", "differential TF footprinting", "chromatin_footprint_tobias"),
            _ev("Method", "TOBIAS ATACorrect/ScoreBigwig/BINDetect (Tn5-corrected)",
                "chromatin_footprint_tobias"),
            _ev("Cell-type groups", f"{ga} vs {gb}", "chromatin_footprint_tobias"),
            _ev("Motifs tested", n_tested, "chromatin_footprint_tobias"),
            _ev("Significant differential motifs", n_sig, "chromatin_footprint_tobias"),
        ]
        n_checked = cross.get("n_checked")
        n_conc = cross.get("n_concordant")
        if n_checked:
            evidence.append(_ev("Top TFs cross-checked vs RNA", n_checked,
                                "footprint_rna"))
            evidence.append(_ev("RNA-concordant (TF gene up in same group)", n_conc,
                                "footprint_rna"))
        evidence = [e for e in evidence if e.value is not None]

        cross_txt = ""
        if n_checked:
            cross_txt = (f" Of {n_checked} top differential factors cross-checked, "
                         f"{n_conc} have the TF gene's RNA concordantly higher in the "
                         f"same group (associative cross-modal support).")
        claim = (
            f"Tn5-bias-corrected footprinting found {n_sig} of {n_tested} TF motifs "
            f"with differential binding occupancy between {ga} and {gb}.{cross_txt}")

        figures = []
        for tf, paths in (footprinting.get("aggregate_plots") or {}).items():
            png = (paths or {}).get("png")
            if png:
                figures.append({"id": f"footprint_{tf}", "path": png,
                                "caption": f"{tf} aggregate footprint ({ga} vs {gb}); "
                                           "central dip = protected TF binding site."})
        return NarrativeBlock(
            id="chromatin.differential_tf_footprinting", modality="chromatin",
            analysis="differential_tf_footprinting", block_type="result",
            title="Differential TF footprinting", status="success",
            confidence="medium", claim=claim, evidence=evidence, figures=figures,
            caveats=[Caveat(
                "Differential TF footprinting is associative: a footprint-occupancy "
                "difference between cell-type groups is not evidence that the factor "
                "regulates a gene or drives the cell state. Footprints are pseudobulk "
                "and Tn5-bias-corrected; concordant RNA is supporting association, not "
                "causation.", severity="info")],
            metrics={"differential_summary": summary,
                     "rna_cross_evidence": cross or None},
        )

    # ── Methods + tables ────────────────────────────────────────────────────
    def _bulk_atac_methods(self, findings: dict) -> list[str]:
        """Copy-pasteable bulk ATAC methods (aligner-input/caller/counting/
        model parameters), sourced ONLY from measured findings so the section
        never fabricates a parameter that was not actually used."""
        out: list[str] = []
        if _ok(findings.get("qc")):
            out.append(
                "Bulk ATAC QC reports only measured library metrics "
                "(per-sample fragment/read counts, mitochondrial fraction, "
                "duplicate rate); FRiP is reported only when post-peak-calling "
                "read counting succeeds and is never default-filled.")
        peaks = findings.get("peaks") or {}
        if _ok(peaks):
            line = (
                f"Bulk ATAC peak calling: MACS3 on the provided aligned "
                f"BAM/CRAM input(s) (genome {peaks.get('genome', '?')}); ARIA "
                f"consumes existing alignments and does not realign in this "
                f"lane. {peaks.get('n_peaks', '?')} peaks were called")
            if peaks.get("consensus_peaks_path"):
                line += (", and per-replicate peaks were merged into a "
                         "consensus peak universe")
            cmd = peaks.get("macs3_cmd")
            line += f". Exact MACS3 command: {cmd}." if cmd else "."
            out.append(line)
        counts = findings.get("peak_counts") or {}
        if _ok(counts):
            out.append(
                f"Peak quantification: per-peak counts for each sample using "
                f"{counts.get('counting_method', 'bedtools coverage')} over "
                f"the called peak universe ({counts.get('n_peaks', '?')} peaks "
                f"x {counts.get('n_samples', '?')} samples), producing a "
                f"peak-by-sample count matrix.")
        da = findings.get("differential_accessibility") or {}
        if da.get("data_type") == "bulk_ATAC" and da.get("ran"):
            covariates = da.get("covariates_adjusted") or []
            cov_txt = (f", adjusting for {', '.join(map(str, covariates))}"
                       if covariates else "")
            dropped_covariates = [
                d.get("covariate")
                for d in (da.get("covariates_dropped") or [])
                if isinstance(d, dict) and d.get("covariate")
            ]
            comparisons = da.get("comparisons") or []
            design = next((c.get("fitted_design_formula") for c in comparisons
                           if isinstance(c, dict)
                           and c.get("fitted_design_formula")), None)
            line = (
                f"Bulk ATAC differential accessibility: technical samples were "
                f"aggregated to biological replicate-level counts "
                f"(condition::replicate), all-zero peaks dropped, then tested "
                f"with {da.get('method', 'replicate-level DESeq2')} on "
                f"{da.get('n_replicate_samples', '?')} replicate samples "
                f"comparing {da.get('condition_col', 'condition')} levels "
                f"(replicate factor {da.get('replicate_col', 'replicate')}"
                f"{cov_txt})")
            if design:
                line += f" with design {design}"
            if dropped_covariates:
                line += (
                    "; requested covariate(s) not adjusted for after "
                    "replicate aggregation: "
                    f"{', '.join(sorted(set(map(str, dropped_covariates))))}"
                )
            line += (
                f". Peaks with |log2FC| >= {da.get('lfc_min', '?')} and "
                f"adjusted p <= {da.get('padj_max', '?')} were called "
                f"differentially accessible, and non-converged log2 "
                f"fold-changes were gated out. At least "
                f"{da.get('min_replicates_per_condition', '?')} replicates per "
                f"condition are required")
            if any(isinstance(c, dict) and c.get("low_power_warning")
                   for c in comparisons):
                line += ("; this run carried a low-power warning because a "
                         "condition had fewer replicates than the default "
                         "floor")
            line += "."
            out.append(line)
        motifs = findings.get("motifs") or {}
        if _ok(motifs) and motifs.get("ran"):
            src = motifs.get("motif_source") or {}
            out.append(
                f"TF motif enrichment: DA peaks were split by accessibility "
                f"direction (both conditions, no one-sided pruning) and tested "
                f"for {motifs.get('method', 'hypergeometric')} motif "
                f"over-representation against the tested-peak background using "
                f"the versioned local {src.get('collection', 'motif')} "
                f"collection (release {src.get('release', '?')}); offline, no "
                f"network egress. An enriched motif is an associative database "
                f"match in the peak set, not evidence of TF activity or "
                f"regulation.")
        return out

    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        findings = unwrap_chromatin_findings(agent_result)
        if _is_bulk_atac(findings):
            return self._bulk_atac_methods(findings)
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
        regulatory = findings.get("regulatory")
        if _ok(regulatory):
            out.append(
                "scATAC regulatory layers are optional and input-gated: motif "
                "activity requires an explicit motif-to-peak map, gene scores "
                "and peak-to-gene links require local gene coordinates, Tn5 "
                "footprinting requires fragments plus a Tn5 bias model, and "
                "scRNA label transfer is report-only rather than an inferential "
                "grouping.")
        return out

    def figures(self, agent_name: str, agent_result: dict,
                report_dir: Path | None = None) -> list[dict]:
        # W0.1 (scATAC P0): surface figures rendered by
        # `_narrative_chromatin.generate_figures` into findings["figures"]
        # ({key: png_path}). Empty until the figure pipeline runs (honest).
        findings = unwrap_chromatin_findings(agent_result)
        figs = findings.get("figures") or {}
        out = []
        for key, value in figs.items():
            if isinstance(value, str):
                out.append({
                    "id": key,
                    "path": value,
                    "caption": str(key).replace("_", " "),
                })
        return out

    def tables(self, agent_name: str, agent_result: dict,
               report_dir: Path | None = None) -> list[dict]:
        findings = unwrap_chromatin_findings(agent_result)
        tables = []
        peaks = findings.get("peaks") or {}
        if peaks.get("peaks_path"):
            tables.append({"id": "chromatin.peaks",
                           "path": peaks["peaks_path"],
                           "label": "Called peaks"})
        if peaks.get("consensus_peaks_path"):
            tables.append({"id": "chromatin.consensus_peaks",
                           "path": peaks["consensus_peaks_path"],
                           "label": "Consensus peaks"})
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
        regulatory = findings.get("regulatory") or {}
        for key, label in (
            ("motif_activity", "Motif activity"),
            ("gene_scores", "Gene activity scores"),
            ("peak_to_gene", "Peak-to-gene links"),
            ("label_transfer", "scRNA label-transfer hypotheses"),
        ):
            sub = regulatory.get(key) or {}
            if sub.get("output_csv"):
                tables.append({"id": f"chromatin.{key}",
                               "path": sub["output_csv"],
                               "label": label})
        return tables
