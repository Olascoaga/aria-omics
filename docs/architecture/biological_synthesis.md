# Biological Synthesis

`BiologicalSynthesisAgent` produces the report's **Integrated Biological
Discussion** — the section that answers "what does it all mean together, what
converges, what conflicts, what can and cannot be claimed". It is an **evidence
integrator, not a writer**: it organizes the structured results ARIA already
produced, detects cross-analysis patterns deterministically, and writes only what
the counts support. It never reads raw data, searches external literature,
proposes unmeasured mechanisms, or adds general knowledge.

Durable decision: **ADR-028**. Code: `aria/agents/biological_synthesis_agent.py`
and `aria/agents/narrative/synthesis/`.

## Where it sits

```
structured agent results
        ↓
Claim Compiler (claim tiers)
        ↓
BiologicalSynthesisAgent
   ├─ synthesis/pattern_detector.py   (deterministic — code guarantees)
   └─ synthesis/discussion_composer.py (emits integration.* NarrativeBlocks)
        ↓
[ render_blocks · evidence_verifier · devil's advocate · run-ledger linkage ]
        ↓
NarrativeAgent → "Integrated Biological Discussion" section
```

It runs inside `NarrativeAgent._collect_narrative_blocks`, after the per-modality
narrators, and appends `integration.*` blocks to the same list. Because the
output is `NarrativeBlock`s, the synthesis **inherits** the existing governance
instead of adding a parallel one.

## Design rules

1. **The LLM does not decide patterns (ADR-002: LLM proposes, code guarantees).**
   `pattern_detector.py` is pure set/sign math over the structured DE + pathway
   results. No LLM, no hardcoded biology. The numbers are reproducible.

2. **Emit blocks, inherit governance.** Every integration block passes through:
   - **Claim Compiler** — observational omics caps at the *associative* tier; the
     spec's E1–E5 inference levels map onto the existing tiers (E5
     `causal_experimental` requires an interventional design).
   - **Strict evidence verification** (`evidence_verifier`) — every number a claim
     states must be on the block's evidence card, or render rejects it.
   - **Devil's advocate** and **run-ledger linkage**.

3. **Claim boundaries.** Causal vocabulary (`CAUSAL_PATTERNS`) is kept out of
   `block.claim`; the "this is not causal" disclaimers live in caveats (which are
   not claim-verified). A limitations block is mandatory whenever any integrated
   claim is made. Conflicts are surfaced as *discordant evidence*, never smoothed.

4. **Reference-anchored comparisons only.** Convergence / direction-concordance
   claims are made only for contrast pairs that share a reference level — because
   "same direction" is not comparable across pairs that swap numerator and
   denominator.

5. **`data_only=True` is the only mode.** A `literature_augmented` mode is
   intentionally not built; a data-supported discussion is the point of ARIA.

6. **Fail-safe.** A synthesis block that cannot be verified is dropped
   (`_drop_unverifiable_synthesis_blocks`), and report rendering is resilient
   (`render_blocks(..., strict=False)`): one unverifiable block is withheld, never
   aborts the whole report.

## What it detects

`pattern_detector.detect_bulk_patterns(contrasts)` reads
`agent_results['bulk_rna_agent']['findings']['contrasts']` — each carries
`all_sig_gene_ids`, `up_gene_ids`/`down_gene_ids`, and `pathways` — and computes:

- **within-contrast convergence** — a contrast produced both DE and pathway
  enrichment;
- **cross-contrast convergence** — shared DE genes between two contrasts, how many
  move in the **same direction**, and the shared enriched processes;
- **cross-contrast divergence** — contrast-specific genes and processes, plus
  discordant (opposite-direction) shared genes;
- **top shared genes** — the shared DE genes ranked by significance, named by
  symbol (raw Ensembl/numeric ids are skipped);
- **reliability** — power and low-power flags that bound interpretation.

## What it writes

The Integrated Biological Discussion section (Slice 1, bulk single-modality):

1. **Main integrated pattern** — the headline: the shared-vs-distinct program, the
   processes it points to, and each condition's specific program.
2. **Convergent evidence** — shared genes (count + same-direction count), the
   shared enriched processes the convergence points to, and the top shared genes
   by symbol.
3. **Divergent / condition-specific signal** — contrast-specific gene counts and
   the processes each condition uniquely engages; discordant shared genes.
4. **Limits of interpretation** (mandatory) — associative scope, no causal/direct
   regulation, candidates prioritized for independent validation; power range.

Sections appear only when the data supports them — e.g. the cross-modal section is
omitted (not faked) when there is no ATAC.

## Scope

- **Slice 1 (current):** bulk RNA-seq, single modality (within-contrast,
  cross-contrast, reliability).
- **Slice 2 (deferred):** cross-modal patterns (RNA expression + ATAC
  accessibility + motifs + peak-to-gene), to be designed against a real
  multi-omic run, not blind.

## Tests

`tests/test_biological_synthesis.py` enforces the spec guards: no fabricated
pathways, no causal language, cross-modal only when the modality exists, mandatory
limitations, conflicts surfaced, every claim mapped to evidence, the main-pattern
headline, named processes and genes — plus a real-data golden that runs against a
generated report on the machine and skips in CI.
