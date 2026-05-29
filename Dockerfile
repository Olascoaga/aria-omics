# ARIA reproducible RNA stack
# ---------------------------
# Closes part of F-ENG-REPRO (senior audit 2026-05-28): pins the scRNA/bulk-RNA
# analytical environment in an image so a report's provenance (git SHA + input
# SHA-256 + this image digest) is bit-reproducible off Samael's machine.
#
# Build:   docker build -t aria-rna:local .
# Run TUI: docker run --rm -it -v "$PWD":/work -v "$HOME/.aria":/root/.aria aria-rna:local aria
#
# NOTE: only the RNA stack is containerized here. chromatin / hic / integration
# stacks get their own images as those modalities (v4.6+) are productionized.

FROM mambaorg/micromamba:1.5.8

# Pinned numpy<2 and the strict conda-forge>bioconda priority live in the yml.
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/aria-rna-env.yml /tmp/aria-rna-env.yml
RUN micromamba install -y -n base -f /tmp/aria-rna-env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install ARIA itself into the solved environment.
WORKDIR /opt/aria
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/aria
RUN python -m pip install --no-deps -e .

# Deterministic numba cache location (matches the local validation gates).
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["aria"]
