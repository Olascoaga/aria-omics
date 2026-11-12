# ARIA v4.5.4 Release Notes

`v4.5.4` is a scientific-honesty hardening release on top of the `v4.5.3`
pre-ATAC integrity freeze.

## Scientific Behavior Change

- scRNA pseudobulk DE now defaults to `fdr_strategy="per_cluster"`.
- `rna_pseudobulk_de.py` still computes and reports both `padj_local` and
  `padj_global`.
- The selected strategy controls the primary significant set (`n_significant`,
  `top_genes`, `all_sig`, up/down counts) and downstream ORA inputs.
- Legacy results without `fdr_strategy` keep global-FDR wording when rendered.

This is intentional: significant-gene counts can differ from pre-`v4.5.4`
reports because the primary FDR family changed.

## Report Honesty

- Pseudobulk power now includes the effective global-BH cutoff
  (`effective_alpha_global`) and `power_estimate_at_effective_alpha`.
- Nominal-alpha power remains visible but is labeled as an upper bound.
- DE based on counts reverse-engineered from log-normalized values is capped at
  low confidence and carries an explicit caveat.

## Validation

- `python -m compileall -q aria` passed.
- aria-env targeted suite: 128 passed / 4 skipped.
- aria-rna-env pydeseq2-gated synthetic DE benchmark: 3 passed.
- `python -c "import aria; print(aria.__version__)"` returned `4.5.4`.
- `git diff --check` passed.
