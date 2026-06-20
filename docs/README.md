# ARIA Documentation

ARIA is documented by validation level. The project intentionally separates
workflows that are validated on controlled + small real datasets from beta,
alpha, and scaffolded roadmap modules. "Validated" never means
publication-grade for a specific study — a domain expert must still review the
design, the fitted model, and the conclusions.

## Start Here

- [Architecture Overview](architecture/overview.md)
- [Code Dependency Graph](architecture/code_graph.md)
- [Generated Graphify Map](architecture/graphify/README.md)
- [Design Principles](architecture/design_principles.md)
- [Reporting and Outputs](architecture/reporting_and_outputs.md)
- [Biological Synthesis (Integrated Discussion)](architecture/biological_synthesis.md)
- [Installation Guide](INSTALLATION.md)
- [Validation Status](validation_status.md)
- [v4.7.0 Release Notes](release_notes_v4.7.0.md)
- [v4.6 Release Notes](release_notes_v4.6.md)
- [v4.6.0-alpha Release Notes](release_notes_v4.6.0-alpha.md)
- [v4.5.5 Release Notes](release_notes_v4.5.5.md)
- [v4.5.4 Release Notes](release_notes_v4.5.4.md)
- [v4.3.19 Release Notes](release_notes_v4.3.19.md)
- [v4.3.18 Release Notes](release_notes_v4.3.18.md)
- [v4.3.17 Release Notes](release_notes_v4.3.17.md)
- [v4.3.16 Release Notes](release_notes_v4.3.16.md)
- [v4.3.15 Release Notes](release_notes_v4.3.15.md)
- [v4.3.14 Release Notes](release_notes_v4.3.14.md)
- [v4.3.13 Release Notes](release_notes_v4.3.13.md)
- [v4.3.12 Release Notes](release_notes_v4.3.12.md)

## Validated Workflows (controlled + small real datasets)

- [Bulk RNA-seq](workflows/bulk_rna.md)
- [Single-cell RNA-seq](workflows/scrna.md)
- [Pseudobulk scRNA-seq from h5ad obs metadata](workflows/pseudobulk_scrna.md)

## Validated / Beta Workflows

- [Trajectory analysis: PAGA + root-gated DPT](workflows/trajectory.md)
- [Cell-cell communication: LIANA](workflows/cell_communication.md)

## Chromatin Workflows (beta, dispatch behind explicit acknowledgement)

These chromatin paths are beta-grade and remain review-required; they dispatch
only after an explicit acknowledgement.

- [Chromatin roadmap: scATAC matrix workflow](workflows/chromatin_roadmap.md) —
  complete beta (de-alpha v4.7.0 / ADR-048): QC/clustering, motif enrichment,
  replicate-gated pseudobulk DA, peak-to-gene links (ADR-050), TOBIAS footprinting,
  publication figures, gene-activity (caveated).
- [Chromatin roadmap: bulk ATAC-seq](workflows/chromatin_roadmap.md) — complete V47
  beta lane: QC + MACS3 peaks + peak×sample counts + replicate-gated DESeq2 DA +
  TF-motif interpretation (validated on ENCODE K562 vs GM12878).

## Alpha Workflows

No modalities are currently at alpha — scATAC was promoted to scoped **beta** in
v4.7.0 (ADR-048). Any future alpha modality dispatches only behind explicit
acknowledgement and stays review-required.

## Scaffolded Roadmap Workflows

These modules exist in code, but should not yet be described as stable
production workflows.

- [Chromatin roadmap: ChIP, CUT&RUN, CUT&TAG](workflows/chromatin_roadmap.md)
- [Genome architecture roadmap: Hi-C / Micro-C](workflows/hic_roadmap.md)
- [Multimodal integration roadmap: WNN, MOFA+ single-cell integration](workflows/integration_roadmap.md)
  (standalone peak-to-gene link recovery is beta; cross-modal WNN/MOFA+ stays scaffold)

## Diagrams

- [ARIA overview](diagrams/aria_overview.mmd)
- [Bulk RNA flow](diagrams/bulk_rna_flow.mmd)
- [scRNA flow](diagrams/scrna_flow.mmd)
- [Pseudobulk h5ad obs flow](diagrams/pseudobulk_scrna_flow.mmd)
- [Validation boundary](diagrams/validation_boundary.mmd)
