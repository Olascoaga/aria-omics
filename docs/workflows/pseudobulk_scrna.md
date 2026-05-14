# Pseudobulk scRNA-seq From h5ad obs Metadata

Validation level: production-like validated analytically for `v4.3.12`.

## Goal

Run between-condition differential expression at the biological-replicate level
inside each cell type or cluster.

This is the preferred path for condition contrasts in scRNA-seq because it
avoids treating individual cells as independent biological replicates.

## Required Metadata

The `.h5ad.obs` table should contain:

- condition column, such as `age_group`, `condition`, `treatment`, or genotype;
- biological replicate column, such as `donor_id`, `sample_id`, or `orig.ident`;
- group/cell-type column, such as `subclass`, `cell_type`, or CellTypist labels;
- optional covariates, such as sex, batch, donor attributes, or chemistry.

## Flow

```mermaid
flowchart TD
    H[h5ad with obs metadata] --> AUDIT[DataAuditAgent inspects obs]
    AUDIT --> DESIGN[DesignAgent confirms factor, groups, replicates, covariates]
    DESIGN --> QC[rna_qc.py processed-h5ad aware QC]
    QC --> ANNO[Annotation: CellTypist or existing/fallback labels]
    ANNO --> OBS{Use native obs design?}
    OBS -->|yes| PB[rna_pseudobulk_de.py]
    OBS -->|no| INJ[rna_inject_condition.py]
    INJ --> PB
    PB --> ORA[rna_pathway_per_cluster.py per group x comparison]
    PB --> N[NarrativeAgent]
    ORA --> N
    N --> R[Report + DE/pathway TSV supplements]
```

The same diagram is stored as
[pseudobulk_scrna_flow.mmd](../diagrams/pseudobulk_scrna_flow.mmd).

## Design Rules

- Conditions must map to biological replicates.
- Each condition should have enough biological replicates for DE.
- Covariates must be recorded in the design formula.
- If ARIA cannot identify a condition or replicate column, it should ask or
  skip pseudobulk rather than guessing silently.
- If condition metadata already exists in `obs`, ARIA should use it directly
  instead of injecting filename-derived labels.
- If metadata is not available in `obs`, ARIA may inject condition labels from
  the user-confirmed sample-to-group design before pseudobulk DE.

## Outputs

- pseudobulk DE summaries per group and comparison;
- significant genes with log2FC and adjusted p values;
- ORA pathway hits per group and comparison;
- summary bar plots;
- supplementary TSV tables;
- methods section with groupby, condition, replicate, covariate, and threshold
  details.

## Current Evidence

The GSE278576 consolidated hippocampus `.h5ad` path completed analytically on a
large processed object:

- processed-h5ad QC retained cells using existing obs metrics;
- pseudobulk DE ran across 18 subclasses;
- pathway ORA, LIANA, and PAGA/DPT outputs were present;
- narrative/report fixes were added after detecting a contradiction between
  executive summary and body sections.

Release policy: donor-level pseudobulk is the primary evidence for
between-condition scRNA claims. `rna_de_per_cluster.py` is optional on
atlas-scale inputs and may time out without invalidating a completed
pseudobulk/ORA/LIANA/trajectory run.
