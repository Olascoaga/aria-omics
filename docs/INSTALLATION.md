# Installing ARIA

ARIA currently targets **Ubuntu or WSL2**, **Python 3.11**, and Linux `linux-64`
Conda environments. Other Linux distributions may work, but the committed exact
locks and automated validation are Linux-specific. Native Windows and macOS are
not release-validated installation targets.

## Choose the snapshot first

ARIA's package version and source revision are separate reproducibility fields.
The current package version is `4.7.0`, while `main` contains post-tag work.

```bash
git clone https://github.com/Olascoaga/aria-omics.git
cd aria-omics

# Current development snapshot
git switch main

# Or a release baseline
git checkout v4.7.0

# Or the exact revision supplied by a collaborator
git checkout <commit-sha>
```

Record the resolved revision:

```bash
git rev-parse HEAD
git status --short
```

An archival run should start from an empty `git status --short`. Generated
reports record the commit and dirty state, but a clean checkout makes that
provenance easier to review.

## Requirements

| Requirement | Minimum guidance |
|---|---|
| OS | Ubuntu 22.04+ or WSL2 Ubuntu |
| CPU architecture | `linux-64` / x86-64 for committed exact locks |
| Python | 3.11 |
| RAM | 8 GB for small examples; 16 GB+ for typical analysis |
| Disk | 10 GB for core setup; substantially more for raw reads and all science stacks |
| Environment manager | Conda, Miniforge, or compatible Mamba installation |
| Network | Needed to clone/install and for optional public-data/provider access |

LLM API access is optional when a supported local endpoint is configured.
Consumer chat subscriptions do not normally include API usage. ARIA does not
promise provider pricing or free-tier availability; check the provider's current
terms before configuring a cloud model.

## Option A: guided Ubuntu/WSL2 installer

For a workstation setup:

```bash
bash install.sh
```

The installer:

1. checks Ubuntu/WSL2 resources and network access;
2. installs required system packages with `apt`;
3. installs Miniforge if Conda is absent;
4. creates `aria-env` with Python 3.11;
5. installs the committed pinned core requirements, ARIA, and the Textual
   Control Center;
6. offers Anthropic/Gemini provider configuration;
7. downloads an optional PBMC 3k smoke dataset;
8. runs core doctor and mock integration checks.

It is interactive and needs `sudo` for system packages. Its Python core uses the
committed lock, but the full workstation setup still depends on the host OS and
the separately installed science environments. It creates the core orchestrator
environment, not every modality-specific science environment.

After it finishes:

```bash
conda activate aria-env
aria doctor --smoke
aria doctor --llm
aria
```

## Option B: manual source install

Use this route when system dependencies and Conda already exist.

### Convenience development install

```bash
conda create -n aria-env python=3.11 -y
conda activate aria-env
python -m pip install --upgrade pip
python -m pip install -e '.[tui]'
aria doctor --smoke
```

This resolves the dependency ranges declared in `pyproject.toml`. Add the
development extra only on a development machine:

```bash
python -m pip install -e '.[dev,tui]'
```

### Exact core reconstruction

`requirements.lock` is the pinned Linux core/orchestrator fallback, including
the Control Center dependency:

```bash
conda create -n aria-env python=3.11 -y
conda activate aria-env
python -m pip install --requirement requirements.lock
python -m pip install --no-deps -e .
aria doctor --smoke
```

The last command installs the selected source revision without asking pip to
re-resolve dependencies.

## Install the scientific environments

ARIA dispatches heavy scripts to isolated environments and fails if the exact
registered stack is missing. It never falls back to `aria-env` for scientific
compute.

Install only the stacks required by the intended workflow. Each exact Conda
lock is an `@EXPLICIT` Linux artifact list; a matching `.pip.lock` adds only the
PyPI packages that were part of the validated environment.

```bash
conda create --name aria-rna-env \
  --file envs/aria-rna-env.linux-64.lock
conda run --name aria-rna-env python -m pip install --no-deps \
  --requirement envs/aria-rna-env.pip.lock
```

Use the same two-command pattern for any registered environment that has both
lock files:

| Workflow need | Environment |
|---|---|
| Bulk count matrices, scRNA, pseudobulk, shared DESeq2 DA | `aria-rna-env` |
| Raw scRNA FASTQ ingestion | `aria-ingestion-env` |
| Raw bulk RNA FASTQ preprocessing | `aria-rnaseq-env` |
| Raw bulk/scATAC alignment | `aria-atacseq-env` |
| Chromatin QC, peaks, clustering, regulatory layers | `aria-chromatin-env` |
| TOBIAS footprinting | `aria-tobias-env` |
| External RNA benchmark comparators | `aria-bench-env` |

Hi-C and integration environments are explicit scaffolds and are not active
release-lock targets. Installing those YAMLs does not promote their runtime
readiness or enable normal dispatch.

Audit installed active stacks against the registry and locks:

```bash
python -m aria.utils.environment_audit --env aria-rna-env
```

Maintainers can validate a clean reconstruction below an empty temporary root:

```bash
python -m aria.utils.environment_audit \
  --clean-root /tmp/aria-clean-envs \
  --env aria-rna-env
```

Do **not** run `scripts/generate_locks.sh` as an installation step. That script
regenerates release artifacts from already validated maintainer environments.

## Configure an LLM provider

ARIA can use Anthropic, OpenAI, Gemini, or a local Ollama-compatible endpoint.
Credentials are read from the process environment or `~/.aria/.env`; model
routing can be overridden in `~/.aria/config.yaml`.

Example private environment file:

```bash
mkdir -p ~/.aria
chmod 700 ~/.aria
printf '%s\n' 'ANTHROPIC_API_KEY=replace-me' > ~/.aria/.env
chmod 600 ~/.aria/.env
```

Equivalent variables are `OPENAI_API_KEY`, `GEMINI_API_KEY`, and
`GOOGLE_API_KEY`. Never commit credentials to the repository. Load the file
when needed:

```bash
set -a
source ~/.aria/.env
set +a
aria doctor --secrets
aria doctor --llm
```

A minimal local-provider configuration looks like:

```yaml
llm:
  heavy:
    provider: ollama
    model: ollama/llama3:70b
    api_base: http://localhost:11434
  medium:
    provider: ollama
    model: ollama/llama3:8b
    api_base: http://localhost:11434
  light:
    provider: ollama
    model: ollama/llama3:8b
    api_base: http://localhost:11434
```

Model identifiers evolve. Prefer a model confirmed by `aria doctor --llm` and
your provider account over copying an old release note.

## Privacy and air-gapped launch

Raw omics files are processed locally, but cloud LLM providers can receive the
question and structured agent context. To prohibit egress from the beginning of
a run, configure a local model and set:

```bash
export ARIA_AIR_GAPPED=1
aria
```

Set the flag **before** launch. The interactive sensitivity checkpoint occurs
after the initial question parse; enabling air gap there cannot retroactively
protect that earlier call. Air-gapped mode also blocks governed connector,
resource-download, enrichment, and child-process egress paths.

## Launch and verify the interface

```bash
conda activate aria-env
python -c 'import aria; print(aria.__version__)'
python -c 'from aria.ui.cockpit import cockpit_available; print(cockpit_available())'
aria
```

`cockpit_available()` should print `True` for the current Control Center. If it
prints `False`, install the UI extra in the same environment:

```bash
python -m pip install -e '.[tui]'
```

ARIA automatically falls back to the classic Rich interface when Textual is
missing, the process is not attached to a TTY, `ARIA_NO_TUI` is set,
`--classic-tui` is passed, or `--reproducible` is active.

See [Control Center](CONTROL_CENTER.md) for views and shortcuts.

## Data intake boundaries

Stable entry points include bulk count matrices, 10x `.h5`, a complete 10x MEX
directory, and `.h5ad`. Beta raw-read and chromatin entry points require typed
metadata or manifests where the assay contract needs them.

Do not rely on filenames to establish the biological design. ARIA intentionally
refuses silent production inference for groups, references, replicates, and
contrasts. Prepare explicit sample metadata and review every design checkpoint.
For scATAC FASTQ, provide the required read roles, whitelist, and genome assets;
for bulk ATAC DA, provide condition, biological replicate, and comparison
metadata.

## Updating a checkout

An update changes the analysis source snapshot. Do not update in the middle of
an archival run.

```bash
git switch main
git pull --ff-only
conda activate aria-env
python -m pip install --no-deps -e .
aria doctor --smoke
git rev-parse HEAD
```

If dependency or environment locks changed, reconstruct or audit the affected
environment instead of adding packages ad hoc.

## Troubleshooting

### The old interface appears

```bash
conda activate aria-env
python -c 'import aria, aria.tui; print(aria.__version__, aria.__file__)'
python -c 'from aria.ui.cockpit import cockpit_available; print(cockpit_available())'
python -m pip install -e '.[tui]'
```

Also verify that `which aria` points into the expected environment and that the
checkout is on the intended commit.

### `conda: command not found`

Initialize the Conda installation, then open a new shell. For Miniforge:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate aria-env
```

### A scientific environment is missing

Do not install its packages into `aria-env`. Create the registered environment
from its exact lock, then run `python -m aria.utils.environment_audit --env ...`.

### A provider cannot be reached

```bash
aria doctor --llm
aria doctor --secrets
```

Check that the configured provider has a matching credential or that the local
endpoint is running. ARIA will not silently replace an unavailable configured
model with fabricated analysis output.

### Installation was interrupted

The guided installer can reuse an existing `aria-env`, but it is not a fully
transactional resume system. Run `aria doctor --smoke` first; if the core
environment is inconsistent, remove only that named environment and recreate it
after confirming no other project depends on it.
