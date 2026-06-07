# ARIA Documentation

ARIA is documented by validation level. The project intentionally separates
workflows that are validated on controlled + small real datasets from beta
paths and scaffolded roadmap modules. "Validated" never means publication-grade
for a specific study — a domain expert must still review the design, the fitted
model, and the conclusions.

## Start Here

- [Architecture Overview](architecture/overview.md)
- [Code Dependency Graph](architecture/code_graph.md)
- [Generated Graphify Map](architecture/graphify/README.md)
- [Design Principles](architecture/design_principles.md)
- [Reporting and Outputs](architecture/reporting_and_outputs.md)
- [Biological Synthesis (Integrated Discussion)](architecture/biological_synthesis.md)
- [Installation Guide](INSTALLATION.md)
- [Validation Status](validation_status.md)
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

- [Trajectory analysis: PAGA + DPT](workflows/trajectory.md)
- [Cell-cell communication: LIANA](workflows/cell_communication.md)

## Scaffolded Roadmap Workflows

These modules exist in code, but should not yet be described as stable
production workflows.

- [Chromatin roadmap: scATAC, ATAC, ChIP, CUT&RUN, CUT&TAG](workflows/chromatin_roadmap.md)
- [Genome architecture roadmap: Hi-C / Micro-C](workflows/hic_roadmap.md)
- [Multimodal integration roadmap: WNN, MOFA+, peak-to-gene](workflows/integration_roadmap.md)

## Diagrams

- [ARIA overview](diagrams/aria_overview.mmd)
- [Bulk RNA flow](diagrams/bulk_rna_flow.mmd)
- [scRNA flow](diagrams/scrna_flow.mmd)
- [Pseudobulk h5ad obs flow](diagrams/pseudobulk_scrna_flow.mmd)
- [Validation boundary](diagrams/validation_boundary.mmd)
