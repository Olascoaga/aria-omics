# ARIA v4.5.2 Release Notes

`v4.5.2` introduces the Narrative Kernel, a structured reporting layer that
turns modality outputs into validated `NarrativeBlock` objects before HTML
composition.

## Highlights

- Added `aria.agents.narrative` with block schema, narrator protocol, registry,
  validators, and HTML block renderer.
- Added scRNA and bulk RNA narrators.
- scRNA reports now compose QC, composition, pseudobulk DE, ORA, LIANA, and
  PAGA/DPT sections from validated blocks.
- Bulk RNA reports now compose QC, contrast, pathway, and power blocks.
- `methodology.json` includes serialized narrative blocks.
- The offline scRNA E2E adapter now records input file SHA-256 hashes so
  harness-rendered reports retain reviewer-grade input provenance.
- Validators enforce evidence-backed claims, explicit failed-analysis blocks,
  low/insufficient visibility, deterministic causal-language caveats, PAGA/DPT
  caveats, and figure/table existence checks at render time.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_narrative_types.py tests/test_narrative_validators.py tests/test_narrator_scrna.py tests/test_narrator_bulk.py tests/test_narrative_render_blocks.py`
  -> 16 passed
- `python -m pytest -q tests/test_pytest_smoke.py`
  -> 86 passed, 4 skipped
