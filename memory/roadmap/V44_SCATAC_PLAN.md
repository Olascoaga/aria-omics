---
status: active
source_of_truth_for: v44_scatac_plan
last_updated: 2026-05-14
---

# v4.4 scATAC Plan

This is the next possible milestone after the 4.3 closeout, but it should not
start unless Samael explicitly asks.

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

## Done Criteria

- QC summary in report: FRiP, TSS, fragment distribution when available.
- LSI/UMAP by cluster/group.
- Peak or accessibility summary.
- Motif enrichment summary.
- Methods section with exact parameters.
- Clear confidence levels and limitations.
