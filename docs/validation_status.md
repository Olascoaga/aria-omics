# Validation Status

Last updated: May 2026.

ARIA uses explicit validation boundaries so users can distinguish mature
workflows from beta analysis paths and implementation scaffolds.

## Production-Like Validated

| Area | Status | Evidence |
|---|---|---|
| Bulk RNA-seq count matrix | Validated | Synthetic regression suite and H9 three-condition workflow |
| scRNA single-sample | Validated | PBMC 3k and GSE278576 single-sample runs |
| scRNA multi-sample | Validated | GSE278576 3-donor concat + Harmony workflow |
| Processed h5ad pseudobulk | Validated pending final rerun review | 40-donor hippocampus run completed analytically; narrative fixes require final report review |

## Validated / Beta

| Area | Status | Notes |
|---|---|---|
| Bulk RNA FASTQ preprocessing | Beta | Scripted path exists; dependency and real-data coverage should expand |
| Trajectory: PAGA + DPT | Beta | Validated on hippocampus subset; exploratory, not causal |
| Cell-cell communication: LIANA | Beta | Validated on GSE278576 annotated h5ad |
| GEO/SRA connector | Beta | GSE183948 path validated; public metadata remains heterogeneous |

## Scaffolded / Roadmap

| Area | Status | Required before stable |
|---|---|---|
| scATAC-seq | Scaffolded | LSI, clustering, DA, motifs, report section, fixtures |
| Bulk ATAC-seq | Scaffolded | Peak count matrix, DA, QC summaries, report section |
| ChIP-seq / CUT&RUN / CUT&TAG | Scaffolded | Clear assay-specific QC and peak interpretation |
| Hi-C / Micro-C | Scaffolded | End-to-end validation, memory-safe fixtures, report integration |
| WNN / MOFA+ / peak-to-gene | Scaffolded | Stable standalone RNA + ATAC paths first |

## Current Release Gate

Before tagging `v4.3.12`, rerun or regenerate the hippocampus report and verify:

- executive summary does not contradict body sections;
- pseudobulk DE, pathway, LIANA, and trajectory outputs are represented when present;
- missing outputs stay explicitly missing;
- supplementary TSV tables are non-empty when structured outputs exist;
- methods record design formula, thresholds, grouping columns, covariates, and warnings.
