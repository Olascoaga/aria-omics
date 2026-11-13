# ARIA 3D genome architecture (Hi-C / Micro-C) image — PRODUCTION / VALIDATED (ADR-082).
# Builds FROM the common base and adds ONLY the hermetic aria-hic-env stack
# pinned in envs/aria-hic-env.yml (pairtools, cooler, cooltools, bioframe,
# hic-straw, pysam, bedtools) used by aria/scripts/hic_*.py. Raw Hi-C FASTQ
# alignment (-SP5M) runs in docker/Dockerfile.atacseq (bwa-mem2), keeping this
# 3D contact-map analysis image isolated and reproducible.
#
# Build:
#   docker build -f docker/Dockerfile.base -t aria-base:local \
#     --build-arg ARIA_IMAGE_REVISION=$(git rev-parse HEAD) .
#   docker build -f docker/Dockerfile.hic -t aria-hic:local \
#     --build-arg ARIA_IMAGE_ENV_SHA=$(sha256sum envs/aria-hic-env.yml | cut -d' ' -f1) .
ARG ARIA_BASE_IMAGE=aria-base:local
FROM ${ARIA_BASE_IMAGE}

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -n base -f /opt/aria/envs/aria-hic-env.yml && \
    micromamba clean --all --yes
RUN python -m pip install --no-deps -e /opt/aria

ARG ARIA_IMAGE_ENV_SHA=unknown
ENV ARIA_IMAGE_KIND=hic \
    ARIA_IMAGE_VALIDATION=validated \
    ARIA_IMAGE_ENV_SHA=${ARIA_IMAGE_ENV_SHA} \
    NUMBA_CACHE_DIR=/tmp/numba_cache

LABEL org.opencontainers.image.title="aria-hic" \
      org.opencontainers.image.description="ARIA 3D genome architecture Hi-C / Micro-C image (ADR-082)"

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["aria"]
