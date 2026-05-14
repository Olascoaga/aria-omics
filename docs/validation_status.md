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
| Processed h5ad pseudobulk | Validated | 40-donor hippocampus rerun completed analytically; report review confirmed pseudobulk, ORA, LIANA, trajectory, figures, and TSV exports |

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

## v4.3.12 Release Gate

The `v4.3.12` closeout checks are:

- executive summary does not contradict body sections;
- pseudobulk DE, pathway, LIANA, and trajectory outputs are represented when present;
- missing outputs stay explicitly missing;
- supplementary TSV tables are non-empty when structured outputs exist;
- methods record design formula, thresholds, grouping columns, covariates, and warnings;
- existing `obs` annotations are reported as reused groupings, not as newly
  inferred Leiden clusters;
- `rna_de_per_cluster.py` is treated as optional on atlas-scale inputs. If it
  times out, the pipeline may still be scientifically valid when donor-level
  pseudobulk, pathway, communication, and trajectory outputs complete.

Latest reviewed hippocampus rerun:
`/home/medusa/.aria/reports/aria_20260514_143352_oligodendrocytes_opcs_microglia_-009/report.html`.
The reviewed run produced real output tables and figures. A narrative wording
bug around reused `obs['subclass']` labels was fixed after this report; new
reports generated from the same code path will state that Leiden was skipped
when input annotations are reused.
