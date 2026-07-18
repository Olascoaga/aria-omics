# ARIA

### Agentic Research Intelligence for supervised omics analysis

> Ask the biological question. Confirm the design. Keep every conclusion tied
> to evidence.

![Version](https://img.shields.io/badge/version-4.7.0-blue)
![Status](https://img.shields.io/badge/status-active-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.11-blue)

ARIA is an open-source analysis system that places a supervised, auditable
control layer around omics workflows. It inspects the inputs, asks the user to
confirm experimental design, runs deterministic modality-specific tools in
isolated environments, and builds an evidence-linked report instead of filling
gaps with plausible prose.

The package version is **4.7.0**. The `main` branch contains post-`v4.7.0`
hardening, the current Control Center, the bulk ATAC V47 beta lane, and
preprint-readiness work. A report records the exact commit and dirty state, so
the package version alone is never treated as a complete reproducibility pin.

## What is ready today?

The runtime registry is the source of truth for modality readiness. This public
summary reflects that registry:

| Modality | Readiness | Dispatch policy |
|---|---|---|
| scRNA-seq | Production | Enabled |
| Bulk RNA-seq count matrices | Production | Enabled |
| Bulk RNA-seq FASTQ | Beta | Enabled with visible limitations |
| scATAC-seq | Beta | Explicit acknowledgement required |
| Bulk ATAC-seq | Beta | Explicit acknowledgement required |
| ChIP-seq, CUT&RUN, CUT&TAG | Scaffold | Disabled |
| Hi-C / Micro-C | Scaffold | Disabled by default; experimental opt-in only |
| WNN / MOFA+ multimodal integration | Scaffold | Disabled |

“Production” and “beta” describe ARIA's validation boundary, not whether a
particular study is ready to publish. A domain expert must still review the
design, fitted model, power, QC, and biological conclusions. The detailed and
test-guarded boundary is in [Validation Status](docs/validation_status.md).

### Validated RNA paths

- Bulk count-matrix DE with explicit metadata, covariates, and contrasts.
- scRNA single-sample and multi-sample QC, integration, clustering, annotation,
  and governed reporting.
- Donor-level pseudobulk DE from `.h5ad` metadata, with raw-count handoff when
  available.
- Optional exploratory trajectory and cell-cell communication layers.
- Atomic, manifest-verified GEO/SRA retrieval for supported payloads.

### Acknowledgement-gated chromatin paths

- scATAC QC, TF-IDF/LSI clustering, accessibility analysis, local motif
  enrichment, beta peak-to-gene links, TOBIAS footprinting, gene-activity
  scoring, and report figures. Single-sample footprint comparisons and
  gene-activity remain explicitly caveated.
- Bulk ATAC alignment or BAM intake, measured QC, MACS3 peaks, reproducibility
  consensus, peak-by-sample counts, replicate-gated DESeq2 differential
  accessibility, annotation, ORA, motif interpretation, and descriptive
  footprinting.

Chromatin beta does not mean autonomous or universally publication-grade.
Missing resources, insufficient replication, and unsupported inference return
structured skips or limitations.

## Why ARIA is different

- **Design before compute.** Groups, organism, factor, batch, biological
  replicates, reference levels, and contrasts are reviewed at checkpoints.
- **No silent fake science.** Missing tools and missing evidence stay missing;
  production paths cannot silently substitute mock outputs.
- **Evidence-linked reporting.** Public claims are compiled from structured
  evidence cards and linked to the planned-vs-run ledger.
- **Honest uncertainty.** Findings remain visible as `HIGH`, `MEDIUM`, `LOW`,
  or `INSUFFICIENT`, including null, skipped, and low-power results.
- **Reproducible provenance.** Reports stamp the version, commit, dirty state,
  workflow hash, input hashes, parameter hashes, environment locks, and LLM
  usage.
- **Deterministic narrative boundary.** LLM calls use `temperature=0`; a fixed
  seed is sent only to backends that accept one, and reports disclose whether
  the answering backend is seed-deterministic.
- **Isolated science stacks.** The orchestrator never substitutes its own
  environment for a missing RNA, ingestion, alignment, chromatin, or
  footprinting environment.

## Quick start

ARIA currently targets Ubuntu and WSL2 with Python 3.11. Exact scientific locks
are published for Linux `linux-64`.

```bash
git clone https://github.com/Olascoaga/aria-omics.git
cd aria-omics
bash install.sh
conda activate aria-env
aria doctor --smoke
aria
```

The installer sets up the core environment, installs the current Textual
Control Center, offers LLM-provider configuration, and downloads the optional
PBMC 3k smoke dataset. Scientific workflows use additional modality-specific
Conda environments; see the [Installation Guide](docs/INSTALLATION.md) before a
real run.

For a lightweight source install when Conda is already available:

```bash
conda create -n aria-env python=3.11 -y
conda activate aria-env
python -m pip install -e '.[tui]'
aria doctor --smoke
```

This convenience install resolves the declared version ranges. For an archival
or shared execution, install from the committed locks and pin the Git commit as
described in the installation guide.

## Control Center and launch modes

Running `aria` on an interactive terminal opens the Textual **ARIA Control
Center** when the `tui` extra is installed. It provides a single intake screen,
checkpoint decisions, live agent progress, findings, readiness, resources,
run-ledger state, and generated artifacts.

```bash
aria                 # Control Center when available
aria --classic-tui   # force the Rich terminal flow
aria --reproducible  # interactive reproducible mode; disables the cockpit
aria --hypotheses    # explicitly enable the speculative hypotheses section
```

The Control Center is a presentation layer over the same orchestrator and
checkpoint resolver. It does not own scientific state and cannot turn a
scaffolded modality into a validated one. See the
[Control Center Guide](docs/CONTROL_CENTER.md) for its views and shortcuts.

The non-interactive runner is available as `aria.headless.run_headless` for
automation and test harnesses. It records every policy-selected checkpoint
answer; callers that need study-specific decisions should provide an explicit
answer policy instead of accepting the default policy blindly.

## Provider and privacy boundary

ARIA supports Anthropic, OpenAI, Gemini, and local Ollama-compatible endpoints.
Provider credentials are usage-billed independently of consumer chat
subscriptions; consult the provider's current terms and pricing rather than
assuming a free tier.

Raw analysis files are processed locally, but a configured cloud model can
receive the biological question and structured context needed by an agent. For
a run that must prohibit egress from process start, configure a local model and
set the flag before launching ARIA:

```bash
export ARIA_AIR_GAPPED=1
aria
```

The in-run sensitivity checkpoint is useful, but it occurs after the initial
question parse and cannot retroactively protect that first call. Run
`aria doctor --llm` and `aria doctor --secrets` before sensitive work.

## Scientific environments

`aria.utils.environment_specs` owns the routing and lock policy.

| Environment | Purpose |
|---|---|
| `aria-rna-env` | scRNA, bulk DE, pseudobulk, and shared DESeq2-backed DA |
| `aria-ingestion-env` | Raw scRNA FASTQ ingestion with kallisto/bustools |
| `aria-rnaseq-env` | Raw bulk RNA FASTQ QC, alignment, and quantification |
| `aria-atacseq-env` | Bulk/scATAC raw-read alignment |
| `aria-chromatin-env` | Chromatin QC, peaks, clustering, and regulatory layers |
| `aria-tobias-env` | Tn5-bias-corrected footprinting |
| `aria-bench-env` | External benchmark comparators |
| `aria-hic-env`, `aria-integration-env` | Dispatch-disabled scaffold stacks |

Do not regenerate locks just to install ARIA. Use the committed lock artifacts;
lock generation is a maintainer/release operation.

## Outputs you can audit

A completed report bundle can contain:

```text
report.html                 evidence-linked analysis report
methodology.json            structured methods, claims, caveats, and provenance
figures/                    rendered analysis figures
tables/                     supplementary TSV outputs
ro-crate-metadata.json      optional reproducibility metadata
```

ARIA also exposes ledger-aware utilities:

```bash
aria diff REPORT_A REPORT_B
aria export REPORT_DIR
aria reproduce CAPSULE.zip
```

`aria reproduce` verifies a capsule and can compare provenance and public claim
sets; it is not a command that reruns an entire study automatically.

## Documentation map

- [Documentation home](docs/README.md)
- [Installation](docs/INSTALLATION.md)
- [Control Center](docs/CONTROL_CENTER.md)
- [Validation status](docs/validation_status.md)
- [Architecture overview](docs/architecture/overview.md)
- [Code dependency graph](docs/architecture/code_graph.md)
- [Generated Graphify map](docs/architecture/graphify/README.md)
- [Reporting and outputs](docs/architecture/reporting_and_outputs.md)
- [Workflow guides](docs/workflows/)
- [v4.7.0 release notes](docs/release_notes_v4.7.0.md)

Historical release notes describe the state at those releases. They are not the
source of truth for the current `main` branch.

## Development

```bash
python -m pip install -e '.[dev,tui]'
python -m pytest -q tests/test_docs_drift_guard.py \
  tests/test_preprint_audit_q1_doc_claims.py \
  tests/test_packaging.py
```

Before changing shared runtime paths, consult the curated
[Code Dependency Graph](docs/architecture/code_graph.md) and the generated
[Graphify map](docs/architecture/graphify/README.md). Keep modality readiness in
the runtime registry first; documentation mirrors that source and is guarded in
tests.

## Versioning and reproducibility

- `aria/version.py` is the single package-version source.
- `v4.7.0` is the current release tag; `main` is newer.
- Clone `main` for current development, a tag for a release baseline, or an
  explicit commit for an exact collaboration snapshot.
- The preprint evidence freeze is fail-closed. Evidence receipts are bound to a
  clean indexed source snapshot, so any tracked source or documentation change
  makes existing receipts stale until they are regenerated.

## License

MIT — free for academic and commercial use.
