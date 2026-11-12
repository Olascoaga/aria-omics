# ARIA v4.6.1 Release Notes

`v4.6.1` is a **patch release** that consolidates the post-`v4.6` tri-audit
hardening phase (ADR-046). It is **disclosure / correctness / test-debt** work on
the validated RNA + reporting baseline — **no DE/DA math changed**, and it is **not**
a modality promotion: scATAC stays `alpha` + `requires_ack`.

## Tri-audit remediation (T1–T14 + T16 closed; T15 parked)

Three independent adversarial audits (internal Claude + Codex + Gemini) were
consolidated into one tracker and remediated one finding at a time, failing-test
first. Highlights:

- **T1 (blocker) — annotation↔inference boundary.** LLM-only / marker-only
  cell-type labels are report-only hypotheses and can no longer define an
  inferential pseudobulk groupby; automatic grouping comes only from a
  CellTypist-backed or input-obs trusted label.
- **T2 — integration overcorrection is operative.** A `blocking` integration-QC
  finding now stops differential abundance / pseudobulk dispatch and caps legacy
  downstream blocks, instead of being cosmetic.
- **T3 / T4 — CellTypist governance.** Model-hub downloads obey W-PRIV egress
  (structured `EgressBlocked` under `ARIA_AIR_GAPPED=1`), and the immune-default
  fallback requires an explicit acknowledgement before any import/download/read.
- **T5 / T6 — ingestion & QC.** Only complete 10x MEX triplets collapse to a
  sample (an incomplete MEX is reported at CP1 with the missing component named),
  and a degenerate `MAD=0` no longer collapses the QC bounds to zero width.
- **T7 — annotation-confidence honesty.** A failed extraction of CellTypist
  per-cell confidence is surfaced as a visible degradation, distinct from
  "no probabilities available", without fabricating label uncertainty.
- **T8 — role-aware executive-summary verification.** The W-CLAIM evidence card
  for the executive summary no longer lets an identifier number ("cluster 7")
  back a measured-quantity claim ("7 DE genes").
- **T9 — compositional-dependence disclosure.** Differential-abundance reports
  disclose that per-cell-type CLR tests are not independent (sum-to-zero), naming
  scCODA/propeller; the runtime CLR-OLS method is unchanged.
- **T10 — test debt cleared.** Gemini's "NumPy 2.0 ABI break / 91-of-117 fail"
  Blocker did **not** reproduce; the real state was 6 stale/inaccurate tests, all
  fixed (incl. the scATAC DA simulator now honoring `n_cells_per_condition`).
- **T11–T14 — reviewer-defensibility polish.** Mixed-case (mouse/rat) gene
  verification, the differential-abundance FDR threshold named in the claim,
  whole-word progenitor root matching for DPT, and a single-source PAGA
  strong-connectivity constant.
- **T16 — reproducibility verified.** Each conda/pip lockfile pins numpy once and
  consistently with the installed env; all three modality stacks import with no
  ABI error (RNA on numpy 2.x; chromatin/ingestion on numpy 1.26.4 for MACS2/kb).

**Parked:** T15 (scATAC apeGLM LFC≈0 with a significant Wald p on non-converged
peaks) — fixed at the scATAC de-alpha via pydeseq2 convergence filtering. Alpha
lane; not a preprint blocker.

## Honesty boundary

- The DE/DA inferential core was never touched (`_run_deseq2`, pre-registered FDR,
  lfcThreshold-in-Wald, apeGLM, donor-level CLR-OLS).
- scATAC remains `alpha` + `requires_ack`; bulk ATAC / ChIP / CUT&RUN / CUT&TAG
  stay scaffolded; Hi-C stays opt-in/experimental.
- This is a version-line patch only; no scientific claim was upgraded.

## Validation

- Full test suite: **686 passed / 0 failed / 40 skipped** in `aria-env`
  (numpy 2.4.4).
- `v4.5.4` remains the last stable RNA/reporting tag and is unchanged; `v4.6`
  and `v4.6.0-alpha` are retained as history.
