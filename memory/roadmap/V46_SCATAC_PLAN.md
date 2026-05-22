---
status: active
source_of_truth_for: v46_scatac_plan
last_updated: 2026-05-22
---

# v4.6 scATAC Plan

Moved from `v4.4` on 2026-05-15. `v4.4` is now the Publication Readiness
sprint (see `V44_PUBLICATION_READINESS.md`). scATAC must inherit the new
guarantees the v4.4 sprint introduces (composition correction patterns,
global FDR, provenance block, conda lockfiles, methodology.json,
reproducible mode). Moved again from `v4.5` to `v4.6` on 2026-05-20 so
`v4.5` can add Raw Ingestion first (see `V45_RAW_INGESTION_PLAN.md`).
Do not start scATAC implementation until Raw Ingestion is tagged and
validated.

## Input

Primary validation input:

`/home/medusa/Samael/Erosion/data_inputs/muon_processed/hc11_paired.h5mu`

Known context:

- same-cell paired RNA + ATAC;
- approximately 3,143 cells x 60,990 peaks;
- single-donor validation scale;
- related to the GSE278576 hippocampus project.

## Current Code

Existing scaffold:

- `aria/agents/chromatin_agent.py`
- `aria/scripts/chromatin_qc.py`
- `aria/scripts/chromatin_peaks.py`

## Required First Step

Create `aria-chromatin-env` and a lockfile before adding new analysis scripts.
Do not mix chromatin dependencies into RNA unless there is a deliberate
documented decision.

## Implementation Order

1. Validate `chromatin_qc.py` standalone on the `.h5mu`.
2. Add `chromatin_lsi_clustering.py`:
   - TF-IDF;
   - SVD/LSI;
   - drop depth-associated component if needed;
   - neighbors/UMAP/Leiden.
3. Add `chromatin_diffacc.py`:
   - per-cluster accessibility;
   - per-condition pseudobulk DA when replicates exist.
4. Add `chromatin_motifs.py`.
5. Add `_narrative_chromatin.py`.
6. Wire chromatin findings into `NarrativeAgent`.
7. Add focused smoke tests before broad refactors.

Narrative expectation: follow the post-v4.5.1 scRNA narrative-depth pattern
from `ADR-009`. The chromatin narrative should not stop at aggregate counts;
it should describe each major QC, accessibility, and motif result with local
evidence, limitations, and how it supports or fails to support the biological
question.

## Done Criteria

- QC summary in report: FRiP, TSS, fragment distribution when available.
- LSI/UMAP by cluster/group.
- Peak or accessibility summary.
- Motif enrichment summary.
- Methods section with exact parameters.
- Clear confidence levels and limitations.
