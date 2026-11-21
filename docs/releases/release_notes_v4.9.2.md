# ARIA v4.9.2 Release Notes

`v4.9.2` is a **patch release** over `v4.9.1` closing every verified defect
surfaced by an independent **tripartite blind post-Wave-4 re-audit** (cross-checked
at `file:line` across Gemini, Claude, and Codex across all five analytical
domains), delivered in two atomic remediation waves (**Wave 5** and **Wave 6**)
alongside refreshed **`ADR-080`** and **`ADR-082`** empirical certifications, a
fixed **`docker/Dockerfile.rna`** supply-chain build, and a freshly re-verified
**`19 / 19` `automated_v1_gate` preprint-v1 machine receipt freeze** bound to
`v4.9.2`. **`aria/utils/stage_manifest.py` and human-gold reader-study CSVs
(`docs/benchmark_results/preprint_v1/human/*.csv`) remain untouched (`0` diff
relative to `v4.9.0`)**, and all external Tier 1 blind benchmarks
(`docs/benchmark_results/external_benchmarks_v49/`) remain explicitly attributed
to immutable tag **`v4.9.0`**.

## Headline 1: Wave 5 — Consensus, Critical, Regression & Statistical Core Fixes (`f215667`, 76 files)

- **`pydeseq2 >= 0.5.4` Explicit Reference Anchoring & Sign-Flip Guard (`BI-01`):**
  Updated `aria/scripts/rna_bulk/deseq2.py` and `aria/scripts/rna_pseudobulk_de.py`
  to pass explicit `C(factor, Treatment(reference=...))` formulas to `pydeseq2`
  so apeGLM `lfc_shrink` always targets the requested reference level even when
  `test < reference` lexicographically, added exact `_shrink_coeff` resolution
  and a post-shrinkage sign-flip guard, and synchronized `docs/architecture/code_graph.md`.
- **Full-Universe Benjamini–Hochberg FDR in Local ORA (`CDX-01`):**
  Fixed `aria/utils/ora.py` so BH adjustment is computed across the full tested
  gene-set universe (including zero-overlap sets) rather than only `k > 0` hits,
  eliminating anti-conservative q-value bias.
- **Gene-Level Primary Inference in Cross-Modal DE $\times$ DA Concordance (`CLD-01`):**
  Updated `aria/utils/cross_modal_concordance.py` (`summarize_de_da`) so
  `concordance_fraction`, directional/co-significance Fisher exact tests,
  Spearman $\rho$, Pearson $r$, and promoter $\rho$ are computed on deduplicated
  gene-level `primary` rows (`gene_counts`), eliminating multi-peak-per-gene
  pseudoreplication while preserving full pair counts alongside gene counts.
- **Honest Per-Pair FPR vs. Empirical FDR Labeling & Statistical Hardening (`TRI-01`):**
  Separated `empirical_fdr_at_padj` from per-pair permutation FPR, fixed
  `run_ledger.py` & `evidence_verifier.py` edge cases, hardened `lobo_jackknife.py`,
  `falsification_certificates.py`, `sensitivity.py`, and `power_estimation.py`,
  and restored full precision in TOBIAS `change` / `raw_change` outputs.
- **Symmetric Sparse Yang et al. (2017) HiCRep SCC (`k_op @ m_sym @ k_op.T`):**
  Fixed upper-triangle boundary asymmetry in `aria/scripts/hic_qc_and_balance.py`
  by symmetrizing banded cis contact matrices prior to 2D box-kernel smoothing,
  plus hardened Multiome (`integration_wnn.py`, `integration_mofa.py`,
  `integration_peak2gene.py`, `eregulon_core.py`), Chromatin/TOBIAS, and scRNA
  pipelines (`rna_qc` cache schema `7 -> 8`).
- **Docker `release-repro` Build Fix (`CI-01`):**
  Installed prebuilt `conda-forge::python-annoy` in `docker/Dockerfile.rna`
  before solving `envs/aria-rna-env.yml`, resolving the source-compile `g++`
  failure encountered on tag `v4.9.1`'s `release-repro` build.

## Headline 2: Wave 6 — Domain-Wide HIGH & MEDIUM Hardening (`4899f8c`, 96 files)

- **Governance, Reporting & Hypotheses:** Added HTML XSS escaping and staged
  relative bulk table links in report renderers, universal `allow_nan=False`
  JSON serialization guards, canonical `_canon` / q-value ranking and
  case-insensitive quarantine matching in `HypothesisAgent`, and strict
  non-empty / finite validation across reporting and ledger builders.
- **Bulk RNA & scRNA Polish:** Achieved complete dual 300-DPI PNG + editable SVG
  figure parity across scRNA renderers, hardened numeric contrast level coercion,
  fixed `per_read` dual-UMI accounting in `rna_umi_dedup.py`, and fixed per-chrom
  strand exon capping in `rna_quantify.py`.
- **Chromatin, TOBIAS, Multiome & Hi-C Polish:** Hardened cohort-aware LSI
  doublet estimation, long-gene peak-to-TSS annotation, depth-normalized motif
  deviations, SHA-256 + PID/UUID-safe temp files in Multiome and Hi-C stage
  caches, and Hi-C $P(s)$ scaling, sorted loop anchors, `chromosomes` filtering,
  and `switch_regions` exports.

## Headline 3: Refreshed Empirical Certifications & `19 / 19` Preprint `v1` Freeze

- **Refreshed `ADR-080` Cross-Modal Concordance Receipt (`docs/benchmark_results/multiome_concordance/adr080_cross_modal_concordance_validation.json`):**
  Re-executed on the 10-donor balanced `GSE278576` human hippocampal cohort
  (`all_invariants_passed=True`): the disjoint-donor split (`5` brains RNA-only
  $\times$ `5` disjoint brains ATAC-only, `Oligo_vs_Astro`) achieves **`93.19%`
  gene-level directional agreement** across `11,613` tested genes (`3,215`
  concordant cis pairs spanning `2,007` unique concordant genes, gene-level
  **Fisher `OR = 218.71`**, **Spearman $\rho = 0.5686$**, **Pearson $r = 0.6056$**),
  **`15` concordant `JASPAR2024` TFs** (`SOX8`, `SOX10`, `NFIA`, `MEIS1`,
  `OLIG2`, `TCF12`, `RFX2`, `RFX3`, `NFIB`, `ZIC2`, `BNC2`, `JUN`, `NFATC1`,
  `FEZF2`, `RORA`), **`30` shared full-precision ORA pathways (`Jaccard = 0.3333`)**,
  **`100%` (`10/10`)** same-cell `Signac` `peak2gene` linked-gene recall, exact
  reversed-contrast (`Astro_vs_Oligo`) sign-flip invariance, `18,371` matched
  20-unit CPM cis links, and **`0` false discoveries** across both 10-donor
  (`143,266` pairs) and 20-unit (`140,506` pairs) deranged-label null controls
  (`0` `W-CLAIM` / `C3` / `C4` violations).
- **Refreshed `ADR-082` 3D Genome Architecture Receipt (`docs/benchmark_results/hic/adr082_hic_3d_genome_validation.json`):**
  Re-certified in `aria-hic-env` (`all_gates_pass=True`, `aria_version: 4.9.2`):
  symmetric sparse Yang et al. (2017) HiCRep SCC elevates within-condition
  replicate concordance to **`0.9845`** vs **`0.5937`** between-condition
  (separation gap `+0.3908`, up from `0.8213` vs `0.7273`), while preserving
  `r = 1.0000` vs standalone `cooltools 0.7.1` on GC-phased `E1` and Crane
  insulation, `100%` TAD boundary, loop, and `A_to_B`/`B_to_A` switch recall,
  and **`0` false switches / boundaries** across balanced `4v4` null permutations.
- **Fresh `automated_v1_gate` Preprint Freeze (`19 / 19` Verified on `v4.9.2`):**
  All `19 / 19` machine-verifiable lanes across `Claims 1–7` are regenerated
  and cryptographically bound to `v4.9.2` under
  `docs/benchmark_results/preprint_v1/inventory.json`, while the `3` independent
  `human_gold` reader-study lanes remain untouched (`0` diff) as
  `deferred_to_v2_reader_study`.

## Compatibility & Verification

- **8 new regression test suites** added (`tests/test_wave5{a,b,c,d}_*.py`,
  `tests/test_wave6{a,b,c,d}_*.py`); verified with **`2,931+` passed** in
  `aria-env`, **`254` passed (`0` failed)** in `aria-rna-env`, live checks in
  `aria-chromatin-env` and `aria-hic-env`, green `ARIA CI`, and green
  tag-triggered `release-repro` Docker + CycloneDX SBOM builds.
- **Zero diff against `v4.9.0`** on `aria/utils/stage_manifest.py` and
  `docs/benchmark_results/preprint_v1/human/*.csv`.
- External Tier 1 blind benchmarks (`docs/benchmark_results/external_benchmarks_v49/`)
  remain explicitly attributed to immutable tag `v4.9.0`.
