# ARIA v4.7.0 Release Notes

`v4.7.0` is a **minor release** over `v4.6.1`. It ships the interactive
"ARIA Control Center" (opt-in Textual TUI over the canonical headless path), the
scATAC external-concordance benchmarks (P3a/P3b), and the **scATAC de-alpha to a
scoped beta** (ADR-048). The **DE/DA estimation core is unchanged**; scATAC keeps
`requires_ack`. This is a normal release, **not** the frozen preprint-artifact tag
(B2 multi-annotator + a clean-checkout benchmark regeneration remain preprint work).

## Headline: scATAC de-alpha alpha -> beta (scoped) — ADR-048

scATAC is promoted from `alpha` to a **scoped `beta`** in the single-source
registry `orchestrator_agent.MODALITY_VALIDATION` (+ `chromatin_agent` + the
`ChromatinAuditAgent` readiness card). `requires_ack` is **retained** (beta ->
status yellow -> CP3.5 acknowledgement still gates dispatch). The promotion is
evidence-backed and scoped:

- **Beta-grade:** QC/clustering, motif enrichment, and the **pseudobulk
  condition-DA** path (donor-level, with biological replication; ARIA enforces
  `min_replicates_per_condition>=3`).
- **Caveated:** the per-cluster **Wilcoxon marker** path stays flagged as
  single-sample-fragile (exploratory, not the primary DA evidence).

### Evidence

- **Cross-tool clustering/motif (SnapATAC2, HC11).** P3b wires a real SnapATAC2
  2.9.0 reference pipeline (dedicated `aria-bench-atac-env`) on the SAME ARIA peak
  matrix; the P3a scorer reports cluster **ARI 0.532 / NMI 0.669** (85-99%
  bidirectional cluster purity — the ARI penalty is granularity, not
  disagreement). Motif null on HC11 (CIS-BP vs JASPAR2024 namespace + local genome
  naming, disclosed honestly).
- **Multi-replicate pseudobulk DA (edgeR/limma/DESeq2).** ARIA's pseudobulk DA
  (`rna_bulk_de._run_deseq2`) was validated against edgeR-QLF / R-DESeq2 /
  limma-voom on **5 young (20-39) vs 5 old (80-100) GSE278576** hippocampus donors
  over a genomic-overlap consensus peak universe. Over 134,276 shared tested peaks,
  ARIA's 2,052 DA peaks are a clean **subset of R-DESeq2's 6,646 (recall 1.000)** and
  a superset of the conservative methods (all 9 edgeR + 248/249 limma calls), with
  **LFC Spearman 0.61-0.76** and sign agreement 0.80-0.95. The earlier weak HC11 DA
  concordance was a single-sample artifact of the Wilcoxon marker path, not the
  pseudobulk math.

This validation is reproducible in-repo: `aria/benchmarks/consensus_pseudobulk.py`,
`concordance_atac.score_da_lfc_concordance`, `scripts/aria_pseudobulk_da_from_tsv.py`,
and the CLI `scripts/run_scatac_multisample_da_concordance.py`. Artifacts:
`docs/benchmark_results/scatac_concordance/p3a_concordance_snapatac2.json` and
`p3_multisample_aging_da_concordance.json`.

## ARIA Control Center (opt-in Textual TUI) — ADR-047

A UI-agnostic read-model (`aria/runtime/experiment_view.py`) drives a Textual
"control center" (`aria/ui/cockpit.py`) plus a front-door intake. It surfaces run
state, per-agent progress, checkpoint decisions, findings by confidence, the run
ledger, modality readiness cards, a resource center, an artifact browser, and
cross-experiment resume/history. Textual is an **opt-in `tui` extra**; the cockpit
launches only on a real TTY and never becomes load-bearing for a claim —
`aria/headless.py` stays the canonical reproducible path. Pure Rich renderers in
`aria/ui/render.py` keep the read-model UI-toolkit-free.

## scATAC external-concordance benchmarks (P3a/P3b)

- **P3a** — pure-numpy concordance scorer (`aria/benchmarks/concordance_atac.py`:
  cluster ARI/NMI, exact + genomic DA-peak overlap, motif top-k/RBO, seed
  stability) + dispatchable runner + CLI. Missing comparator -> honest `not_run`.
- **P3b** — wired SnapATAC2 reference driver in the dedicated `aria-bench-atac-env`
  (numpy>=2, kept apart from the R `aria-bench-env`). Honest `motifs=None` when
  genome/DB assets are unavailable (no fabrication, ADR-002).

## Other changes since v4.6.1

- **Startup egress fix.** Importing LiteLLM no longer fetches the remote model cost
  map; ARIA sets `LITELLM_LOCAL_MODEL_COST_MAP=True` before import.
- **Raw ingestion hardening.** Bulk MatrixMarket sidecars no longer misclassify as
  scRNA without 10X cell-barcode evidence; `RawIngestionAgent` honors CP1-confirmed
  modalities; scans publish visible progress counters.
- **Documentation freeze sync.** README, `docs/validation_status.md`, and the
  capability matrix now report scATAC as beta consistently with ADR-048.

## Scope guard

DE/DA estimation core untouched (`_run_deseq2`, pre-registered FDR,
lfcThreshold-in-Wald, apeGLM). Hi-C / bulk ATAC / ChIP / CUT&RUN / CUT&TAG remain
scaffolded. Existing tags are not moved. scATAC stays `requires_ack`.

## Not in this release (preprint freeze, tracked separately)

Per the external readiness review, the frozen preprint depends on: (1) a complete
multi-annotator **B2** evidence-governance study, and (2) regenerating all cited
benchmark artifacts from a single clean checkout into a dedicated
`preprint_v1/` directory. Those are deliberately **not** bundled here; `v4.7.0` is
a normal release, not the immutable preprint-artifact tag.
