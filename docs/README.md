# ARIA Documentation

ARIA is documented by validation level. The project intentionally separates
workflows that are already production-like from beta paths and scaffolded
roadmap modules.

## Start Here

- [Architecture Overview](architecture/overview.md)
- [Design Principles](architecture/design_principles.md)
- [Reporting and Outputs](architecture/reporting_and_outputs.md)
- [Installation Guide](INSTALLATION.md)
- [Validation Status](validation_status.md)
- [v4.3.19 Release Notes](release_notes_v4.3.19.md)
- [v4.3.18 Release Notes](release_notes_v4.3.18.md)
- [v4.3.17 Release Notes](release_notes_v4.3.17.md)
- [v4.3.16 Release Notes](release_notes_v4.3.16.md)
- [v4.3.15 Release Notes](release_notes_v4.3.15.md)
- [v4.3.14 Release Notes](release_notes_v4.3.14.md)
- [v4.3.13 Release Notes](release_notes_v4.3.13.md)
- [v4.3.12 Release Notes](release_notes_v4.3.12.md)

## Production-Like Validated Workflows

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
