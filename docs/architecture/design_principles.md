# Design Principles

## LLM Proposes, Code Guarantees

LLMs can interpret intent, summarize, and help draft biological prose. They
must not be trusted to enforce scientific invariants.

Code must guarantee:

- sample-to-group mappings;
- contrast generation;
- threshold application;
- pathway input gene lists;
- report input schemas;
- missing-result handling;
- production-vs-development mock behavior.

## No Silent Fake Science

Production runs must fail loudly when a required tool or dependency is absent.
Mock results are allowed only with explicit development opt-in:

- `allow_mock=true`
- `allow_mocks=true`
- `ARIA_ALLOW_MOCKS=1`
- `ARIA_DEV_MODE=true`

Reports must never silently substitute fake analyses for missing dependencies.

## Missing Results Stay Missing

If enrichment, annotation, figures, integration, trajectory, or communication
analysis fails, the exact warning should propagate forward. Narrative code must
not infer plausible biology from absent outputs.

## Resume Is File- and Parameter-Validated

Heavy steps can resume only when outputs are present and valid for the current
parameters. ARIA should validate:

- file existence;
- parseable JSON summaries;
- expected columns;
- non-empty outputs;
- input modification times;
- matching parameter or manifest signatures.

## Design Before Compute

The most dangerous failures usually happen before statistics:

- wrong condition labels;
- technical replicates treated as biological replicates;
- hidden batch factors;
- wrong organism or genome assembly;
- unclear control/reference group.

ARIA therefore treats experimental design confirmation as a first-class
checkpoint.

## Honest Uncertainty

Every major finding should have a confidence level:

- `HIGH`: strong structured evidence;
- `MEDIUM`: plausible, but with caveats;
- `LOW`: weak or partial support;
- `INSUFFICIENT`: data or outputs do not support a claim.

Low-confidence and insufficient-evidence findings must remain visible in the
report.
