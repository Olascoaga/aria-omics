# ARIA v4.5.5 Release Notes

`v4.5.5` closes the v4.5 RNA governance freeze with the first executable
Benchmark A1 artifact for the preprint methods lane.

## Benchmark A1

- Added a repeatable A1 runner for preliminary synthetic bulk-DE validation:
  `scripts/run_a1_bulk_de_benchmark.py`.
- The A1 manifest reports the four frozen axes:
  FDR calibration, LFC concordance, ranking concordance, and significant-call
  concordance at FDR 0.05.
- The benchmark uses ARIA's real bulk DESeq2 path (`_run_deseq2`) with apeGLM
  enabled against a neutral synthetic truth set.
- The A1 output includes a simple Fig. 1 SVG summary and a structured JSON
  manifest under `docs/benchmark_results/`.

## Scope Boundary

- This is the ARIA-path A1 preliminary lane from the frozen v4.5 benchmarking
  plan.
- External DESeq2/edgeR/limma comparator execution remains assigned to the
  separate `aria-bench-env` lane, so the A1 artifact does not overclaim method
  superiority or identity.

## Validation

- `scripts/run_a1_bulk_de_benchmark.py` executed in `aria-rna-env`.
- Focused benchmark/unit tests were run for the A1 scorer and synthetic-DE
  benchmark support.
- Version metadata is derived from `aria.version.__version__ == "4.5.5"`.
