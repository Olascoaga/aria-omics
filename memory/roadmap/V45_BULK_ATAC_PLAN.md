---
status: active
source_of_truth_for: v45_bulk_atac_plan
last_updated: 2026-05-14
---

# v4.5 Bulk ATAC Plan

Deferred until v4.4 scATAC is stable unless Samael explicitly reprioritizes it.

## Existing Code

- `chromatin_qc.py` has a bulk chromatin QC branch.
- `chromatin_peaks.py` wraps MACS3 narrow/broad peak calling.

## Missing Pieces

- Peak count matrix handling.
- Differential accessibility over peak counts.
- Replicate-aware DESeq2-style model for bulk ATAC peak counts.
- Narrative/report section in `_narrative_chromatin.py`.
- Supplementary tables for peaks and DA results.

## Done Criteria

- QC summary suitable for ATAC libraries.
- Peak-calling summary.
- Differential accessibility table.
- Motif/pathway-style interpretation only when data supports it.
- Methods section with aligner/caller/counting/model parameters.
