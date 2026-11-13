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


_MODALITY_IMAGES = ("Dockerfile.rna", "Dockerfile.rnaseq", "Dockerfile.chromatin",
                    "Dockerfile.atacseq", "Dockerfile.tobias",
                    "Dockerfile.integration", "Dockerfile.hic")


def test_per_modality_dockerfiles_exist():
    for name in ("Dockerfile.base",) + _MODALITY_IMAGES:
        assert (DOCKER / name).exists(), f"missing docker/{name}"


def test_base_image_has_no_science_env_solve():
    # The common base layer must stay thin (no micromamba env solve) so the
    # heavy conda layer is per-modality and cacheable.
    base = (DOCKER / "Dockerfile.base").read_text()
    assert "micromamba install" not in base


def test_modality_images_build_from_base():
    for name in _MODALITY_IMAGES:
        text = (DOCKER / name).read_text()
        assert "ARIA_BASE_IMAGE" in text and "FROM ${ARIA_BASE_IMAGE}" in text, \
            f"docker/{name} must build FROM the common base"
        assert "micromamba install" in text, \
            f"docker/{name} must add its own conda env layer"


def test_stack_images_install_their_own_env():
    expected = {
        "Dockerfile.chromatin": ("aria-chromatin-env.yml", "chromatin"),
        "Dockerfile.atacseq": ("aria-atacseq-env.yml", "atacseq"),
        "Dockerfile.tobias": ("aria-tobias-env.yml", "tobias"),
        "Dockerfile.integration": ("aria-integration-env.yml", "integration"),
        "Dockerfile.hic": ("aria-hic-env.yml", "hic"),
    }
    for name, (env_yml, kind) in expected.items():
        text = (DOCKER / name).read_text()
        assert f"/opt/aria/envs/{env_yml}" in text, \
            f"docker/{name} must install envs/{env_yml}"
        assert (ROOT / "envs" / env_yml).exists(), f"missing envs/{env_yml}"
        assert f"ARIA_IMAGE_KIND={kind}" in text, \
            f"docker/{name} must stamp ARIA_IMAGE_KIND={kind}"


def test_no_modality_image_is_stamped_scaffold_or_beta():
    # ADR-079 promoted the last scaffold image (integration); every modality
    # image is now a validated release image. A future scaffold image must be
    # added here deliberately, never silently.
    for name in _MODALITY_IMAGES:
        text = (DOCKER / name).read_text()
        assert "ARIA_IMAGE_VALIDATION=scaffold" not in text, name
        assert "ARIA_IMAGE_VALIDATION=beta" not in text, name


def test_validated_integration_image_is_stamped_validated():
    text = (DOCKER / "Dockerfile.integration").read_text()
    assert "ARIA_IMAGE_VALIDATION=validated" in text
    assert "ADR-079" in text
    assert "SCAFFOLD" not in text


def test_validated_chromatin_images_are_stamped_validated():
    for name in ("Dockerfile.chromatin", "Dockerfile.atacseq", "Dockerfile.tobias",
                 "Dockerfile.integration"):
        text = (DOCKER / name).read_text()
        assert "ARIA_IMAGE_VALIDATION=validated" in text, \
            f"docker/{name} must stamp ARIA_IMAGE_VALIDATION=validated"
        assert "ARIA_IMAGE_VALIDATION=scaffold" not in text
        assert "ARIA_IMAGE_VALIDATION=beta" not in text
