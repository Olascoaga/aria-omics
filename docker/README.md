# ARIA container images (P2-2)

One image per modality on a shared base layer. The base carries the OS,
micromamba, and the ARIA source; each modality image adds only its conda
science env. This keeps the heavy per-modality conda solve isolated and the
common layer cached, mirroring ARIA's conda-per-env design.

| Image | Dockerfile | Status |
|-------|-----------|--------|
| `aria-base` | `Dockerfile.base` | common layer (no science env) |
| `aria-rna` | `Dockerfile.rna` | **validated** — built + benchmarked in release CI |
| `aria-chromatin` | `Dockerfile.chromatin` | scaffold (v4.6; not built in CI) |
| `aria-integration` | `Dockerfile.integration` | scaffold (v4.7; not built in CI) |

## Hermeticity

`.dockerignore` (repo root) keeps `.git`, the **private** `memory/`, build
caches, and reports out of the build context, so they never enter an image
layer. Do not `COPY` those back in.

## Build (RNA)

```bash
REV=$(git rev-parse HEAD)
ENV_SHA=$(sha256sum envs/aria-rna-env.yml | cut -d' ' -f1)

docker build -f docker/Dockerfile.base -t aria-base:local \
  --build-arg ARIA_IMAGE_REVISION="$REV" .
docker build -f docker/Dockerfile.rna -t aria-rna:local \
  --build-arg ARIA_IMAGE_ENV_SHA="$ENV_SHA" .

# Run the TUI:
docker run --rm -it -v "$PWD":/work -v "$HOME/.aria":/root/.aria aria-rna:local aria
```

## Provenance / digest

Each image stamps `ARIA_IMAGE_KIND`, `ARIA_IMAGE_REVISION` (git SHA), and
`ARIA_IMAGE_ENV_SHA` (sha256 of the solved env yml) as env vars. When you run a
**published** image, pass its registry digest so reports cite it:

```bash
docker run --rm -e ARIA_IMAGE_DIGEST="$(docker inspect --format '{{index .RepoDigests 0}}' aria-rna:published)" \
  -e ARIA_IMAGE_REF=aria-rna:published aria-rna:published aria
```

`aria/version.py:collect_image_metadata()` reads these and folds them into every
report's provenance and `methodology.json`. Outside a container all fields are
honestly `null` (`containerized: false`) — ARIA never fabricates a digest.

## Supply chain (CI)

- **Secret scanning:** `gitleaks` runs on every PR/push (`secret-scan` job).
- **SBOM:** the release lane generates a CycloneDX SBOM of the built `aria-rna`
  image (syft) and uploads it as a build artifact.

## Multi-platform / non-RNA env locks

Multi-platform (osx) locks and the chromatin/integration env locks are produced
by building these images on the corresponding CI runners — they are NOT
snapshot-fabricated locally (conda-lock's solver hangs on the bioconda+pip mix;
those envs are not installed on the maintainer's machine). See P2-1.
