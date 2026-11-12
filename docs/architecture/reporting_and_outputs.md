# Reporting, Claim Governance & Reproducible Outputs

ARIA reports are compiled by `NarrativeAgent` (`aria/agents/narrative_agent.py` + `aria/agents/narrative/reporting/`), governed mechanically so that **public scientific claims are derived strictly from verified structured outputs**, never from unconstrained LLM prose.

---

## 1. The `C1–C4` Public-Claim Boundary

Every public sentence in `report.html` and every claim entry in `methodology.json["claims"]` passes through a single pre-render compiler (`claim_compiler.compile_public_claims`):

1. **`C1` Untrusted-Data Envelope & Exclusive Claim Compiler (`prompt_boundary.py`, `claim_compiler.py`):**
   - User questions, filenames, sample labels, gene symbols, and pathway terms enter LLMs exclusively inside a single escaped `<untrusted_data>` JSON envelope.
   - Raw `MessageBus` status text and legacy free-text summaries cannot enter public claim sections; unsupported `NarrativeBlock`s are withheld opaquely before rendering.
2. **`C2` Typed Semantic Fact & Polarity Verification (`SemanticFact`, `evidence_verifier.py`, `compose_prose.py`):**
   - Narrators attach structured `EvidenceCard` metrics and `SemanticFact` records (typed entity/aliases, predicate, polarity, and numeric values).
   - `W-CLAIM` verifies numeric fidelity, entity presence, and directional/null polarity (catching negation inversions).
   - For structured analyses, `compose_block_prose` renders the canonical public sentence deterministically from the `EvidenceCard`, guaranteeing **Claim 6 fail-closed monotone invariance** across LLM providers (provider phrasing can cause an overreaching block to be withheld, never added or elevated).
3. **`C3` Fail-Closed Run-Ledger Linkage (`run_ledger.py`, `W-LEDGER`):**
   - Every non-empty result claim must resolve to an executed `analysis_run` (or allowlisted `derived_report_analysis`) node with provenance (`ledger://<modality>/<analysis>`) plus a typed `evidence_card` node. Claims referencing skipped, failed, or unexecuted stages are withheld as violations.
4. **`C4` Per-Estimand Causal Licensing (`resolve_causal_estimand_license`, `validators.py`):**
   - Observational omics caps at the `associative` tier; causal vocabulary (`drives`, `causes`, `regulates directly`) is rejected or downgraded unless the block's structured contrast uniquely matches a verified `design.causal_estimands` entry with `interventional=True`, `randomized=True`, and `confounding_verified=True`.

---

## 2. Modality Narrators & Cross-Modal Triangulation

`ModalityNarrator` plugins under `aria/agents/narrative/narrators/` translate structured agent findings into governed `NarrativeBlock`s:

- **`BulkRnaNarrator` (`bulk_rna.py`):** Sample QC + outlier-sensitivity retention (`P1-5`), contrast-level `pyDESeq2` DE, local GMT ORA, and GSEA blocks (`F3`: free-text LLM DE prose disabled).
- **`ScRnaNarrator` (`scrna.py` + `narrative/scrna/`):** Per-sample QC, Harmony integration + overcorrection diagnostics, clustering, CellTypist/fallback annotation caveats, donor-aware pseudobulk DE (disclosing pre-registered per-cluster vs global BH-FDR), ORA, exploratory LIANA communication, and root-gated PAGA/DPT trajectory blocks.
- **`ChromatinNarrator` (`chromatin.py`):** Governed `scATAC` and `bulk_ATAC` blocks covering measured QC, consensus peaks, TF-IDF/LSI clustering, replicate-gated `pyDESeq2` differential accessibility, peak annotation, peak ORA, offline motif enrichment, cis `peak2gene`, and replicate-aware/descriptive `TOBIAS` footprinting.
- **`IntegrationNarrator` (`integration.py` + `synthesis/cross_modal_triangulation.py`):** Governed `multiome.{paired_qc,wnn,mofa,peak_to_gene,eregulons,triangulation}` blocks, HTML-safe Methods, dual `300 DPI` PNG + editable SVG figures (`_narrative_integration.py`), CSV table staging (`mofa_factor_loadings.csv`, `peak_gene_links.csv`, `eregulons.csv`, `eregulon_triplets.csv`), and the deterministic **Cross-Modal Evidence Graph** (`methodology.json["cross_modal_evidence_graph"]`, stamping `shared_upstream=["peak_to_gene"]` so shared cis lineage never self-inflates into independent contrast replication).

---

## 3. Report Bundle & Reproducibility Capsule Anatomy

Each execution creates a unique, non-overwritable bundle directory (`exist_ok=False`):

```text
aria_<timestamp>_<slug>_<id>/
├── report.html                    # Self-contained governed HTML report + copy-pasteable Methods
├── methodology.json               # Public claims, evidence cards, RunLedger, caveats, Cross-Modal Graph & provenance
├── figures/                       # Publication figures (PNG + editable SVG pairs where supported)
├── tables/                        # Supplementary CSV/TSV tables (DE, DA, pathways, MOFA+ loadings, peak2gene, eRegulons)
├── memory_snapshot.sqlite         # WAL-consistent single-experiment ARIAMemory subtree (reproducible runs)
├── ro-crate-metadata.json         # W3C/RO-Crate 1.1 provenance graph linking claims -> ledger -> evidence cards
├── capsule.zip                    # Atomic SHA-256 reproducibility capsule (via --reproducible or `aria export`)
└── speculative_hypotheses.json    # Present ONLY when `--hypotheses` is enabled (promotable: false)
```

Completed bundles are inspected and verified with `aria diff`, `aria export`, and `aria reproduce` (see [CLI & Headless Guide](../CLI_AND_HEADLESS.md) and [HypothesisAgent Architecture](hypothesis_agent.md)).

---

## 4. Confidence Visibility & Caveat Honesty

ARIA uses four explicit confidence tiers — `HIGH`, `MEDIUM`, `LOW`, and `INSUFFICIENT` — and **never hides null, skipped, low-power, or degraded findings**:

- Under-replicated or confounded contrasts remain visible as structured skips or `LOW`/`INSUFFICIENT` blocks.
- Marker-panel fallback annotations warn downstream pseudobulk, trajectory, and communication sections that labels are curation targets.
- Degraded WNN anti-degradation guardrails (`warn`/`fail`), paralog TF motif sharing (`motif_family_ambiguity`), library-depth attenuation (`depth_confounded`), and excluded technical MOFA+ factors are surfaced explicitly in both report prose and `methodology.json`.
