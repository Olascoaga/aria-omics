# ARIA v4.3.13 Release Notes

`v4.3.13` is superseded by `v4.3.14`.

`v4.3.13` was a pre-v4.4 maintenance patch for the terminal interface and
scRNA design checkpoint visibility.

## Fixes

- Long pasted TUI prompts are captured as one biological question instead of
  leaking remaining lines into `Launch ARIA?` and checkpoint prompts.
- Checkpoint 1 now displays h5ad-derived covariates, such as `Gender`, when
  DataAudit infers them from `obs`.
- `Homo sapiens` in the user question now maps to the default human genome
  hint (`hg38`) during DataAudit.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
