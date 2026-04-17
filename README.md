# ARIA v3 Updates — Bulk RNA-seq Pipeline

## What changed

### Biology-critical fixes
1. **Label-aware intent parsing** — BMAL1 → label B, REV-ERBα → label R, wildtype → WT
2. **Multiple contrasts** — B vs WT AND R vs WT AND B vs R in one run
3. **TF-aware LFC threshold** — 0.58 for TF knockouts (BMAL1, etc.), 1.0 for others
4. **Cross-contrast overlap** — reports shared vs unique DE genes between contrasts

### Infrastructure fixes (caught during your test run)
5. **STAR FASTA decompression** — SetupAgent now runs gunzip on `.fa.gz` since STAR cannot read gzipped FASTAs. Also decompresses GTF for universal compatibility.
6. **Dispatcher dedup** — pipeline used to run 3× because LLM plans list multiple steps per agent. Now each agent runs once.
7. **featureCounts `-t exon`** — was `-t gene` (counts introns as noise).
8. **Strandedness auto-detection** — real pysam-based detection. Protects against ~50% read loss on stranded libraries.
9. **Replicate concordance QC** — Spearman correlation within groups. Flags samples <0.85 mean r with their replicates.
10. **PCA threshold** — was 2 SD (too aggressive). Now 2.5 SD + combined with replicate concordance for outlier confirmation.

### NarrativeAgent
- `_summarize_bulk_rna` — proper multi-contrast prose
- Embeds SVG volcanos + PCA inline in HTML (no external files)
- Bulk RNA methods section with exact params
- Handles both `bulk_rna_agent` and `scrna_agent` (was hardcoded to `rna_agent`)

## How to apply

```bash
cd ~/Samael/ARIA
bash APPLY.sh
```

The script backs up your current files to `/tmp/aria_backup_*` before overwriting.
Then runs integration tests to verify.

## Before relaunching the pipeline

```bash
# Force STAR re-index with the new decompressed FASTA
rm -rf ~/.aria/genomes/hg38/star_index/

# SetupAgent will decompress genome.fa.gz → genome.fa on next run
# (adds ~30 sec but STAR can actually use it now)
```

## Tests
- `test_integration.py`: 21/21 ✓
- `test_bulk_rna.py`: 22/25 (3 pending — tests for features still rolling in)
