# ARIA Documentation

This documentation is organized around concrete analytical workflows, explicit
validation boundaries, and auditable governance contracts. Historical release
notes under [`releases/`](releases/) describe past tags; [`validation_status.md`](validation_status.md)
and the runtime registry are the source of truth for current `main` (`v4.8.0`).

“Production”, “beta”, and “scaffold” are ARIA validation tiers — none of them
replace domain-expert review of experimental design, fitted models, statistical
power, quality control, or biological conclusions.

---

## Start here

| Goal | Document |
| :--- | :--- |
| Install ARIA or reconstruct an exact lock-pinned checkout | [Installation Guide](INSTALLATION.md) |
| Run the CLI, flags, capsules, diagnostics, or `run_headless` | [CLI & Headless Guide](CLI_AND_HEADLESS.md) |
| Navigate the interactive Textual terminal cockpit | [Control Center Guide](CONTROL_CENTER.md) |
| Inspect enabled modalities, dispatch gates, and readiness tiers | [Validation Status](validation_status.md) |
| Understand the agent orchestration and Conda stack isolation | [Architecture Overview](architecture/overview.md) |
| Audit the `C1–C4` claim compiler, `RunLedger`, and report capsules | [Reporting and Outputs](architecture/reporting_and_outputs.md) |
| Understand the opt-in `HypothesisAgent` (`S0–S10`) & `SPECULATIVE` wall | [HypothesisAgent Architecture](architecture/hypothesis_agent.md) |
| Assess downstream impact before modifying shared code | [Code Dependency & Contract Graph](architecture/code_graph.md) |

---

## Workflow guides

### Production-tier runtime modalities
- [Bulk RNA-seq](workflows/bulk_rna.md) *(count-matrix DE and universal multi-platform raw FASTQ preprocessing are Production)*
- [Single-cell RNA-seq](workflows/scrna.md)
- [Donor-level pseudobulk scRNA-seq](workflows/pseudobulk_scrna.md)

### Optional beta RNA layers
- [Trajectory: PAGA + root-gated DPT](workflows/trajectory.md)
- [Cell-cell communication: LIANA](workflows/cell_communication.md)

### Acknowledgement-gated beta chromatin & paired 10x Multiome
- [scATAC-seq and Bulk ATAC-seq workflows](workflows/chromatin_roadmap.md)
- [Paired same-cell 10x Multiome integration (WNN / MOFA+ / peak2gene / eRegulons)](workflows/integration_roadmap.md)

`scATAC`, `bulk_ATAC`, and paired same-cell RNA+ATAC multiome integration are
scoped `beta` lanes requiring explicit `CP3.5` acknowledgement persisted in
`ARIAMemory`.

### Dispatch-disabled scaffolds & blocked cross-modal scopes
- [ChIP-seq, CUT&RUN, and CUT&TAG](workflows/chromatin_roadmap.md)
- [Hi-C / Micro-C](workflows/hic_roadmap.md)

Unpaired, diagonal, or bulk cross-modal integration combinations fail closed at
`enforce_same_cell_integration_gate`.

---

## Architecture, benchmarks & governance

- [Design Principles](architecture/design_principles.md)
- [Biological Synthesis & Cross-Modal Triangulation](architecture/biological_synthesis.md)
- [HypothesisAgent & `SPECULATIVE` Quarantine (`S0–S10`)](architecture/hypothesis_agent.md)
- [Preprint-v1 Evidence Freeze Protocol](architecture/preprint_freeze_v1.md)
- [Frozen v4.5 Benchmarking Protocol](architecture/benchmarking_v45.md)
- [Benchmark & Concordance Artifacts Index](benchmark_results/README.md)
- [Track B3 Frozen Open-Weight Reagent Ablation Supplement](benchmark_results/b3_frozen_reagent/B3_FROZEN_REAGENT_SUPPLEMENT.md)
- [Maintenance & Benchmark Runner Scripts (`scripts/`)](../scripts/README.md)

---

## Diagrams

- [ARIA system overview](diagrams/aria_overview.mmd)
- [Bulk RNA flow](diagrams/bulk_rna_flow.mmd)
- [scRNA flow](diagrams/scrna_flow.mmd)
- [Pseudobulk scRNA flow](diagrams/pseudobulk_scrna_flow.mmd)
- [Validation boundary](diagrams/validation_boundary.mmd)

---

## Releases & changelog

- [Master Changelog (`CHANGELOG.md`)](../CHANGELOG.md)
- [v4.8.0 Release Notes](release_notes_v4.8.0.md) — **current release tag**
- [Historical Release Archive (`v4.3.12` – `v4.7.0`)](releases/)
