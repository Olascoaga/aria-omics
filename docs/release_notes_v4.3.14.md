# ARIA v4.3.14 Release Notes

`v4.3.14` supersedes `v4.3.13` as the pre-v4.4 TUI/DataAudit maintenance
patch.

## Fixes

- The TUI biological-question field now uses explicit multi-line entry. Paste
  the full prompt and finish with a line containing only `END`.
- Stale queued pasted input is discarded before `Launch ARIA?` and checkpoint
  choice prompts.
- DataAudit can infer `Homo sapiens` / `hg38` from human-style h5ad gene
  symbols when explicit organism metadata is absent.
- Checkpoint 1 continues to show h5ad-derived covariates such as `Gender`.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
