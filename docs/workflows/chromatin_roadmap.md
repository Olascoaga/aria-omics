# Chromatin Workflow Roadmap

Validation level: scaffolded.

This module should not yet be described as production-ready.

## Target Assays

- scATAC-seq;
- bulk ATAC-seq;
- ChIP-seq;
- CUT&RUN;
- CUT&TAG.

## Existing Pieces

- `ChromatinAgent`;
- `chromatin_qc.py`;
- `chromatin_peaks.py`;
- MACS3 parameter profiles for assay types;
- basic report hooks.

## Current v4.6 Entry Path

The planned same-cell RNA+ATAC `.h5mu` input path is now validated at the
reader/QC-scaffold level. On the local HC11 validation file, `chromatin_qc.py`
reads ATAC modality `atac` and reports real dimensions of 3,143 cells x 60,990
peaks. It still does not compute FRiP, TSS enrichment, complete QC, or pass/fail;
those require the v4.6 chromatin stack, peak calling, and reference-backed TSS
enrichment. This module remains scaffolded.

## Required Before Stable

For scATAC:

- fragments / peak matrix input detection beyond the validated `.h5mu` entry
  path;
- TSS enrichment and FRiP QC;
- TF-IDF + LSI;
- clustering;
- differential accessibility;
- motif enrichment;
- peak-to-gene handoff contract;
- small deterministic fixtures;
- NarrativeAgent chromatin section.

For bulk ATAC / ChIP / CUT&RUN / CUT&TAG:

- robust BAM/fragment validation;
- assay-specific QC;
- peak calling and consensus peak strategy;
- count matrix generation;
- differential accessibility or binding;
- motif/pathway summary;
- report methods.

## Proposed Flow

```mermaid
flowchart TD
    I[Fragments / BAM / peak matrix] --> QC[Assay-specific QC]
    QC --> PEAKS[Peak calling or peak matrix]
    PEAKS --> FEAT[Feature matrix]
    FEAT --> DA[Differential accessibility / binding]
    DA --> MOTIF[Motif enrichment]
    MOTIF --> N[NarrativeAgent chromatin section]
```

## Release Recommendation

Do not start integration work before standalone chromatin workflows have
validated fixtures and report sections.
