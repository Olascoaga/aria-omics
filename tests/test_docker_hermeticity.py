"""P2-2: Docker hermeticity + per-modality image structure.

Pure repo-layout checks (no docker daemon needed). They guard that the build
context never leaks the private `memory/`, `.git`, or build caches into image
layers, and that each modality image builds FROM the common base layer.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKER = ROOT / "docker"
DOCKERIGNORE = ROOT / ".dockerignore"


def _ignore_lines():
    return [l.strip() for l in DOCKERIGNORE.read_text().splitlines()
            if l.strip() and not l.strip().startswith("#")]


def test_dockerignore_exists():
    assert DOCKERIGNORE.exists(), ".dockerignore must exist for a hermetic build"


def test_dockerignore_excludes_private_and_vcs_and_caches():
    # Normalize trailing slashes (both `memory/` and `memory` are valid).
    lines = {l.rstrip("/") for l in _ignore_lines()}
    # The private operational memory must NEVER enter an image layer.
    assert "memory" in lines, "must ignore memory/"
    assert ".git" in lines, "must ignore .git"
    for pat in ("__pycache__", "*.egg-info"):
        assert pat in lines, f"must ignore {pat}"


def test_per_modality_dockerfiles_exist():
    for name in ("Dockerfile.base", "Dockerfile.rna",
                 "Dockerfile.chromatin", "Dockerfile.integration"):
        assert (DOCKER / name).exists(), f"missing docker/{name}"


def test_base_image_has_no_science_env_solve():
    # The common base layer must stay thin (no micromamba env solve) so the
    # heavy conda layer is per-modality and cacheable.
    base = (DOCKER / "Dockerfile.base").read_text()
    assert "micromamba install" not in base


def test_modality_images_build_from_base():
    for name in ("Dockerfile.rna", "Dockerfile.chromatin",
                 "Dockerfile.integration"):
        text = (DOCKER / name).read_text()
        assert "ARIA_BASE_IMAGE" in text and "FROM ${ARIA_BASE_IMAGE}" in text, \
            f"docker/{name} must build FROM the common base"
        assert "micromamba install" in text, \
            f"docker/{name} must add its own conda env layer"


def test_scaffold_modalities_are_marked_unvalidated():
    for name in ("Dockerfile.chromatin", "Dockerfile.integration"):
        text = (DOCKER / name).read_text()
        assert "ARIA_IMAGE_VALIDATION=scaffold" in text, \
            f"docker/{name} must stamp itself as a scaffold (not a validated image)"
