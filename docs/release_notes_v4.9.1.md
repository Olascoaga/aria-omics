# ARIA v4.9.1 Release Notes

`v4.9.1` is a **patch release** over `v4.9.0` delivering a full-modality
scientific, statistical, visualization, and contract hardening sweep across all
**9 / 9 single-assay modalities** and **3 cross-cutting multi-modal engines**,
alongside **ORA / GSEA precision and publication-renderer upgrades**, **CI /
Docker / environment hygiene fixes**, and a freshly re-verified **`19 / 19`
`automated_v1_gate` preprint-v1 machine receipt freeze** bound to `v4.9.1`.
The **`pyDESeq2` statistical core (`aria/scripts/rna_bulk/deseq2.py`),
`aria/utils/stage_manifest.py`, and human-gold reader-study CSVs
(`docs/benchmark_results/preprint_v1/human/*.csv`) remain untouched (`0` diff
relative to `v4.9.0`)**.

## Headline 1: ORA Precision Fix & Publication-Grade GSEA Upgrade

- **Full-Precision `padj`/`pvalue`/`qval` Preservation:** Eliminated
  fixed-decimal truncation (`round(q, 5)`) that collapsed `padj < 5e-6` to
  `0.0` and produced artificial `-log10(FDR) = 300` dotplot artifacts across
  `aria/utils/ora.py`, `aria/scripts/rna_bulk/ora.py`,
  `aria/scripts/rna_pathway_per_cluster.py`, `aria/scripts/chromatin_bulk_diffacc.py`,
  `aria/scripts/chromatin_diffacc.py`, `aria/scripts/rna_de_per_cluster.py`,
  `aria/scripts/hic_topology.py`, and `aria/utils/lobo_jackknife.py`.
- **Native 300-DPI PNG + Editable SVG GSEA Renderers (`aria/scripts/rna_pathway_viz.py`):**
  Fixed `NES=inf` on large GO sets by adopting Broad standard bounds
  (`min_size=15, max_size=500, seed=0`) with strict non-finite filtering, added
  offline/air-gap local GMT loading and mouse/human case-insensitive gene-symbol
  matching, fixed ORA dotplot overlap gene counts (`k/n`) and colorbar
  `hue_norm` alignment, and replaced fragile `blitz.plot.*` wrappers with native
  running-sum and multi-pathway top-table renderers displaying full red/blue hit
  barcodes across every contrast.

## Headline 2: Full-Modality 4-Wave Audit Remediation

Following a systematic five-track audit of the entire codebase (`31 CRITICAL`,
`53 HIGH`, `33 MEDIUM` verified defects), four atomic remediation waves hardened
every analytical layer:

- **Wave 1 — Core Statistics, Precision & C1–C4 Governance (`edb079a`, 23 files):**
  NaN/Inf-safe Benjamini–Hochberg FDR, scientific-notation parsing in
  `evidence_verifier.py`, deterministic non-leaking fallback prose in
  `compose_prose.py`, hardened `claim_compiler.py`, `run_ledger.py`,
  `falsification_certificates.py`, `lobo_jackknife.py`, `power_estimation.py`,
  `sensitivity.py`, `reporting/*`, `cross_modal_triangulation.py`, and
  full-precision loop q-values in `hypothesis/adapters/hic.py`.
- **Wave 2 — Visualization, Finite Scales & Dual PNG/SVG Parity (`55df872`, 15 files):**
  Data-adaptive caps for underflowed `p == 0` (eliminating fake `y = 300`
  spikes), case-insensitive MOFA+ view rendering, mappable contact-map window
  selection avoiding the mouse `mm10`/`mm39` `0–3 Mb` telocentric N-gap,
  unphased E1 track labeling, finite-guarded volcano/PCA/PAGA/UMAP/LIANA/TOBIAS
  plots, and dual PNG + SVG parity.
- **Wave 3 — Chromatin, RNA, Multiome & Hi-C Scientific Engines (`2407348`, 41 files):**
  - *Chromatin & Epigenomics (13 files):* finite `log2_enrichment > 0` motif
    filtering (separating depleted motifs), sign-pure TOBIAS TF rankings,
    post-doublet LSI sketching cell counts, full-file scATAC mito-fraction
    streaming, UCSC/Ensembl contig normalization, fail-closed heterologous
    spike-in normalization, `Promoter > Exonic > Intronic > Distal` annotation
    priority without cross-species GTF fallback, preserved overlap peak counts
    post-blacklist, and `O(N)` barcode set lookups.
  - *Bulk RNA & scRNA Engines (12 files):* hardened differential abundance,
    CellTypist species/model handling, CellPhoneDB validation, disconnected
    PAGA/DPT trajectory handling, non-finite QC MAD guards (bumping QC cache
    schema `6 -> 7`), pseudobulk/per-cluster DE guards, and bulk count/transform/QC
    finiteness.
  - *Multiome Integration Engines (8 files):* FASTA contig aliasing and
    CMH skipped-stratum exclusion in `eregulon_core.py`, unversioned Ensembl &
    case-folded GTF matching in `peak2gene_core.py`, singleton/NaN-weight/falsy-zero
    guards in `wnn_diagnostics.py`, NaN-safe MOFA+ ANOVA $\eta^2$, self-loop/Inf
    filtering in graph iLISI, and paired QC safeguards.
  - *Hi-C / Micro-C 3D Genome Engines (5 files):* out-of-core sparse banded
    2D-smoothed **Yang et al. (2017) HiCRep SCC** with `vstran` rank weights
    (eliminating dense `N × N` OOM), integer-multiple insulation window snapping,
    Micro-C `< 4 kb` loop-calling skip guard, per-chromosome phasing validation
    with contig aliasing, sorted viewframes, `n_valid`-weighted $P(s)$ decay,
    bounded `.hic` conversion, and ICE convergence gating on `pass_qc`.
- **Wave 4 — Narrators, Mouse/Species Parity & Agent Wiring (`bbb7b6b`, 23 files):**
  - Aligned all five narrators (`bulk_rna`, `scrna` + helpers, `chromatin`,
    `integration`, `hic`) with canonical script keys (`padj`, `figures_svg`,
    `window_size`, `phasing_method`, `base_resolution`, `jsd_vs_control`,
    `thresholds_used`, primary `fdr_strategy`, exact `cellphone_pval == 0.0`,
    and `pc1_validated`-gated active/inactive wording), and added
    `batch_factor_ids` / `batch_flag` quarantine to `hypothesis/adapters/integration.py`.
  - Unified mouse/human species & assembly normalization (`mouse`, `mm10`,
    `mm39`, `GRCm39`, `hg38`, etc.) across `rna_bulk/ora.py`, `gtf_io.py`,
    `rna_pathway_per_cluster.py`, and `chromatin_peak_ora.py`; wired
    `design_factor` into `_sample_qc`, gated GSEA behind `run_pathways`, and
    implemented two-pass contrast-family global BH-FDR in `rna_bulk_de.py`.
  - Fixed `barcode_namespace`/`barcode_map` forwarding, `"10x_multiome"`
    aliasing, translated ATAC fragment barcode accounting, and exact
    `regulon_padj == 0.0` tunnel retention across `integration_agent.py`,
    `integration_assemble_multiome.py`, and `multiome_contracts.py`; restricted
    `GenomeArchAgent` Micro-C loop resolution selection to `[4 kb, 25 kb]` and
    propagated `figures_svg` + `pc1_validated`.

## Headline 3: CI, Dockerfiles & Preprint `v1` Receipt Refresh

- **CI & Supply-Chain Hygiene:** Synced `CAUSAL_PATTERNS` across concordance
  and report validators, returned structured `missing_conda` from
  `environment_audit.py` when `conda` is absent, updated `envs/aria-rna-env.yml`
  (`decoupler>=1.4`, `pytest>=7.4`) and Python 3.11-compatible Dockerfiles
  (`atacseq`, `tobias`) so both `ARIA CI` and the tag-triggered `release-repro`
  Docker + CycloneDX SBOM lane run green.
- **Fresh `automated_v1_gate` Preprint Freeze (`19 / 19` Verified on `v4.9.1`):**
  All `19 / 19` machine-verifiable lanes across `Claims 1–7` are refreshed and
  cryptographically bound to `v4.9.1` under
  `docs/benchmark_results/preprint_v1/inventory.json`, while the `3` independent
  `human_gold` reader-study lanes remain untouched (`0` diff) as
  `deferred_to_v2_reader_study`.

## Compatibility & Invariants

- **Zero diff against `v4.9.0`** on `aria/scripts/rna_bulk/deseq2.py`,
  `aria/utils/stage_manifest.py`, and `docs/benchmark_results/preprint_v1/human/*.csv`.
- **11 new regression test suites** added (`tests/test_wave{1,2,3,4}_*.py`);
  verified with `1,038 passed, 67 skipped` (`aria-env`), `260 passed`
  (`aria-rna-env`), and live checks in `aria-hic-env` and `aria-chromatin-env`.
- External Tier 1 blind benchmarks (`docs/benchmark_results/external_benchmarks_v49/`)
  remain explicitly attributed to immutable tag `v4.9.0`.
