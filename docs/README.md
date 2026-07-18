# ARIA Documentation

This documentation is organized around current user tasks and explicit
validation boundaries. Historical release notes describe the code at a past
tag; they are not the source of truth for current `main`.

“Production”, “beta”, and “scaffold” are ARIA readiness tiers. None of them
replace expert review of a study's experimental design, fitted model, power,
quality control, or biological interpretation.

## Start here

| Goal | Document |
|---|---|
| Install or reproduce an exact checkout | [Installation](INSTALLATION.md) |
| Use the current terminal interface | [Control Center](CONTROL_CENTER.md) |
| See what is enabled and how strongly it is validated | [Validation Status](validation_status.md) |
| Understand the system flow | [Architecture Overview](architecture/overview.md) |
| Review reporting and evidence rules | [Reporting and Outputs](architecture/reporting_and_outputs.md) |
| Assess the impact of a code change | [Code Dependency Graph](architecture/code_graph.md) |
| Navigate the repository broadly | [Generated Graphify Map](architecture/graphify/README.md) |

## Workflow guides

### Production-tier runtime modalities

- [Bulk RNA-seq](workflows/bulk_rna.md)
- [Single-cell RNA-seq](workflows/scrna.md)
- [Donor-level pseudobulk scRNA-seq](workflows/pseudobulk_scrna.md)

Bulk RNA FASTQ preprocessing is a beta entry path into the production bulk RNA
analysis lane.

### Optional beta RNA layers

- [Trajectory: PAGA + root-gated DPT](workflows/trajectory.md)
- [Cell-cell communication: LIANA](workflows/cell_communication.md)

These are exploratory or associative layers and must not be presented as causal
evidence.

### Acknowledgement-gated beta chromatin

- [scATAC and bulk ATAC workflows](workflows/chromatin_roadmap.md)

Both runtime modalities are beta and require explicit acknowledgement. The same
guide distinguishes beta steps from descriptive or caveated sub-analyses.

### Dispatch-disabled scaffolds

- [ChIP-seq, CUT&RUN, and CUT&TAG](workflows/chromatin_roadmap.md)
- [Hi-C / Micro-C](workflows/hic_roadmap.md)
- [WNN / MOFA+ multimodal integration](workflows/integration_roadmap.md)

Scaffold code is not a supported production workflow. Standalone scATAC
peak-to-gene links use the beta chromatin path; that does not validate the
cross-modal integration scaffold.

## Architecture and governance

- [Design Principles](architecture/design_principles.md)
- [Biological Synthesis](architecture/biological_synthesis.md)
- [Preprint-v1 Evidence Freeze](architecture/preprint_freeze_v1.md)
- [Frozen v4.5 Benchmarking Protocol](architecture/benchmarking_v45.md)

The preprint freeze is fail-closed. Receipts are tied to a clean indexed source
snapshot, including tracked documentation, and become stale after source changes.

## Diagrams

- [ARIA overview](diagrams/aria_overview.mmd)
- [Bulk RNA flow](diagrams/bulk_rna_flow.mmd)
- [scRNA flow](diagrams/scrna_flow.mmd)
- [Pseudobulk scRNA flow](diagrams/pseudobulk_scrna_flow.mmd)
- [Validation boundary](diagrams/validation_boundary.mmd)

## Releases

- [v4.7.0](release_notes_v4.7.0.md) — current release tag
- [v4.6.1](release_notes_v4.6.1.md)
- [v4.6](release_notes_v4.6.md)
- [v4.6.0-alpha](release_notes_v4.6.0-alpha.md)
- [v4.5.5](release_notes_v4.5.5.md)
- [v4.5.4](release_notes_v4.5.4.md)
- [v4.5.2](release_notes_v4.5.2.md)

Older v4.3 release notes remain in this directory as historical records.
