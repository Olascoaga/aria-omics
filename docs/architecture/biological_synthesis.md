# Biological Synthesis & Cross-Modal Triangulation

ARIA separates **observed multi-analysis integration** from both unconstrained LLM prose and quarantined speculation through two deterministic, zero-LLM synthesis layers under `aria/agents/narrative/synthesis/`:

1. **`BiologicalSynthesisAgent` (`ADR-028`):** within-modality and cross-contrast RNA synthesis (`pattern_detector.py` + `discussion_composer.py`), producing the **Integrated Biological Discussion**.
2. **`Cross-Modal Evidence Graph` (`ADR-070`, V48 Slice `I7`):** paired same-cell RNA + ATAC regulatory triangulation (`cross_modal_triangulation.py`), serialized into `methodology.json["cross_modal_evidence_graph"]` and rendered via `IntegrationNarrator` (`multiome.triangulation`).

Neither component reads raw matrices directly, queries external literature, invents unmeasured biology, or bypasses `C1–C4` claim governance.

---

## 1. Governance Pipeline

```text
Structured Agent Results (Bulk RNA / scRNA / Chromatin / Paired Multiome)
        │
        ├─► synthesis/pattern_detector.py + discussion_composer.py  (integration.* NarrativeBlocks)
        └─► synthesis/cross_modal_triangulation.py + IntegrationNarrator (multiome.triangulation block)
                │
                ▼
[ C1–C4 Claim Compiler · W-CLAIM SemanticFact Verifier · Devil's Advocate · W-LEDGER Linkage ]
                │
                ▼
      report.html + methodology.json
```

---

## 2. Single-Modality RNA Synthesis (`pattern_detector.py` + `discussion_composer.py`)

- **Pure Set & Sign Math (`ADR-002`):** Detects within-contrast DE + pathway convergence, reference-anchored cross-contrast shared DE genes (same-direction vs discordant opposite-direction counts), contrast-specific programs, dominant donor-aware scRNA pseudobulk DE blocks, and observational context layers (composition, LIANA, PAGA/DPT).
- **Reference-Anchored Comparisons Only:** Cross-contrast direction concordance is computed exclusively when two contrasts share the exact same denominator/reference level.
- **Mandatory Limitations & Conflict Visibility:** Discordant genes are surfaced explicitly (never smoothed away), `data_only=True` is the sole mode, causal vocabulary is excluded on observational designs, and every integrated section emits an explicit limitations block.

---

## 3. Non-Circular Cross-Modal Triangulation (`cross_modal_triangulation.py`)

In paired 10x Multiome (`RNA + ATAC`) analyses, candidate `eRegulons` (`integration_eregulons.py`) consume `integration_peak2gene.py` BH-FDR cis links on the identical cell universe. Treating their co-occurrence as "independent cross-modal replication" would be circular double-counting.

`cross_modal_triangulation.py` enforces a strict, deterministic 4-tier non-causal taxonomy (`schema_version="1.0.0"`):

| Tier | Definition & Lineage Stamp |
| :--- | :--- |
| **`convergent_paired_multiome`** | `pos_tf_pos_p2g` candidate eRegulons exhibiting concordant per-cell RNA and ATAC AUCell activity (`rna_atac_aucell_corr >= 0.15`), explicitly stamped `shared_upstream=["peak_to_gene"]`, `shared_cell_universe=True`, `independent_contrast_support=False`. |
| **`convergent_cross_contrast`** | Reserved strictly for regulatory axes corroborated by **independent differential contrasts** (`RNA DE` + concordant `ATAC DA` and/or `ATAC motif enrichment`); unexecuted lanes are recorded transparently in `absent_not_tested`. |
| **`associative_cis_candidate`** | BH-FDR `peak2gene` cis links or positive candidate eRegulons without independent condition-contrast corroboration. |
| **`discordant_signal` / `caveated_ambiguous`** | Inverse `neg_tf_pos_p2g` candidates, opposite `RNA LFC` vs `ATAC DA` directions, paralog motif-family ambiguity (`motif_family_ambiguity`), library-depth attenuation (`depth_confounded`), excluded technical MOFA+ factors, or degraded WNN guardrails. |

---

## 4. Verification Anchors

- `tests/test_biological_synthesis.py` (within/cross-contrast bulk & scRNA synthesis, reference-anchoring, conflict visibility, zero causal overreach).
- `tests/test_v48_i7_integration_narrative.py` (bit-exact `Cross-Modal Evidence Graph` determinism, non-circular `shared_upstream` lineage, `C1–C4` + Claim 6 monotone withholding).
