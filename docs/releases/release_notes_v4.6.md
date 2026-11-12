# ARIA v4.6 Release Notes

`v4.6` promotes the scATAC release line out of pre-release (`4.6.0-alpha` → `4.6`).
This is a **version-line** promotion: the pre-4.6 polish gate is closed and the
scATAC dispatch path is live- and multi-sample-validated end to end. It is **not**
a modality promotion — scATAC stays `alpha` + `requires_ack` (see Honesty Boundary).

## What this release line delivers

- **scATAC dispatch lane (alpha):** `ChromatinAgent` dispatches the implemented
  scATAC path after explicit CP3.5 acknowledgement — QC (honest `None` for
  FRiP/TSS on a pre-called `.h5mu`), LSI clustering, differential accessibility
  (per-cluster descriptive + pseudobulk inferential via the shared validated
  DESeq2 core), and governed motif enrichment.
- **Live end-to-end validation (ADR-034):** the scATAC path ran through the real
  headless/orchestrator control path on HC11, closing four bugs no fake-env unit
  test reached.
- **Pseudobulk DA correctness (ADR-041):** numerically validated against
  synthetic truth (recall scales with replication; empirical FDR ≤ 0.01).
- **Multi-sample DA, honestly diagnosed and fixed (ADR-042 → ADR-043):** the
  real multi-donor lane runs on the Erosion consensus peak universe. ADR-043
  corrected the ADR-042 root cause: the per-donor consensus peaks are not
  biologically disjoint (~66% genomic overlap); the apparent block structure was
  exact-string peak matching. Genomic-overlap unification (fragment-free) merges
  the boundary-shifted per-donor duplicates of each region — 142,228 peaks →
  16,057 intervals, n_sig 10,922 → 113 on quantitatively-comparable peaks.

## Honesty Boundary

- **scATAC remains `alpha` + `requires_ack`** in `MODALITY_VALIDATION` and the
  readiness matrix. It is not declared fully autonomous or publication-grade;
  de-alpha is reserved for expert review (ADR-033/042/043 unchanged).
- bulk ATAC, ChIP, CUT&RUN, CUT&TAG remain scaffolded and blocked from dispatch.
- The multi-sample DA result (n_sig=113) is execution validation on comparable
  peaks, not a vetted age-DA biological finding.
- The version badge and this release note are synchronized from
  `aria.version.__version__`.

## Validation

- `test_version_stamps_are_derived_from_single_source` passes after the bump
  (README badge + this release note synchronized from `__version__`).
- scATAC multi-sample harness guards pass (6 passed / 1 dataset-gated skip).
- Pre-4.6 polish gate (§4 master plan) was verified PASS before this promotion;
  the RNA core, governance spine (W-PRIV/W-CLAIM/W-LEDGER/W-CALIB), and benchmark
  spine (A1/A2/B1/B2/B4) remain validated and unchanged.
