# ARIA v4.9.3 Release Notes

`v4.9.3` is a focused **container supply-chain patch release** over `v4.9.2`
making the tag-triggered **`release-repro` Docker & CycloneDX SBOM pipeline
100% green** across all published release images (`Dockerfile.base`,
`Dockerfile.rna`, `Dockerfile.atacseq`, `Dockerfile.chromatin`, and
`Dockerfile.tobias`), paired with a freshly re-verified **`19 / 19`
`automated_v1_gate` preprint-v1 machine receipt freeze** bound to `v4.9.3`.

## What Changed

- **Prebuilt `bioconda::macs3=3.0.4` in `docker/Dockerfile.chromatin`:**
  On tag `v4.9.2`, `Dockerfile.rna` passed for the first time after the
  `conda-forge::python-annoy` fix, unmasking a pre-existing source-compile
  failure in `docker/Dockerfile.chromatin` where `pip` attempted to compile
  `MACS3==3.0.4` and its `cykhash` build dependency from source without `gcc`
  in the thin `mambaorg/micromamba:1.5.8` base image. `Dockerfile.chromatin`
  now pre-installs prebuilt `bioconda::macs3=3.0.4` (and `conda-forge::cykhash`)
  before solving `/opt/aria/envs/aria-chromatin-env.yml`, eliminating C
  compilation during image builds while keeping `envs/aria-chromatin-env.yml`
  and its lockfiles completely untouched.
- **Full `release-repro` + CycloneDX SBOM Verification:**
  Validated across `Dockerfile.base`, `Dockerfile.rna` (including the live
  in-container `tests/test_benchmark_synthetic_de.py` recovery gate and
  CycloneDX SBOM artifact), `Dockerfile.atacseq`, `Dockerfile.chromatin`, and
  `Dockerfile.tobias`.
- **Fresh `19 / 19` `automated_v1_gate` Preprint Freeze Bound to `v4.9.3`:**
  Regenerated all `19 / 19` machine-verifiable freeze lanes across `Claims 1–7`
  under `docs/benchmark_results/preprint_v1/inventory.json` (`git_dirty=false`,
  `0` machine blockers).

## Compatibility & Invariants

- **Zero analytical or statistical code changes** relative to `v4.9.2`
  (only `docker/Dockerfile.chromatin`, version/release metadata, and fresh
  `preprint_v1` receipts).
- **Zero diff against `v4.9.0`** on `aria/utils/stage_manifest.py` and
  `docs/benchmark_results/preprint_v1/human/*.csv`.
- External Tier 1 blind benchmarks (`docs/benchmark_results/external_benchmarks_v49/`)
  remain explicitly attributed to immutable tag `v4.9.0`.
