"""S7 (pre-integration audit): chromatin/ATAC env specs are present + declare their
key tools.

Codex finding #3: CI validated RNA/core but not the chromatin/ATAC stacks, so an
env spec could drift (a missing aligner, a renamed package) and only surface at run
time on a reviewer's machine. A full conda solve of these envs in CI is heavy and
flaky; this fence is the cheap half — it parses each YAML and asserts it exists and
declares the binaries the dispatch lanes call. The honest-skip + dispatch behaviour
itself is exercised by the chromatin test battery the CI now runs (S7 ci.yml step).

Parses files only — runs in any env.
"""

from __future__ import annotations

from pathlib import Path

import pytest

try:
    import yaml
except Exception:  # pragma: no cover - PyYAML ships with the toolchain
    yaml = None

_ENVS = Path(__file__).resolve().parents[1] / "envs"

# Each chromatin/ATAC env -> the tokens that MUST appear in its dependency specs
# (the tools the dispatch lanes actually invoke).
_REQUIRED_TOOLS = {
    "aria-atacseq-env.yml": ("bwa-mem2", "samtools", "chromap"),
    "aria-chromatin-env.yml": ("snapatac2",),
    "aria-tobias-env.yml": ("tobias", "samtools"),
}


def _dep_text(yaml_name: str) -> str:
    path = _ENVS / yaml_name
    assert path.exists(), f"missing env spec: envs/{yaml_name}"
    text = path.read_text(encoding="utf-8")
    if yaml is not None:
        # Parse to confirm the YAML is structurally valid (not just grep-able).
        doc = yaml.safe_load(text)
        assert isinstance(doc, dict) and doc.get("dependencies"), (
            f"envs/{yaml_name} has no dependencies block"
        )
    return text


@pytest.mark.parametrize("yaml_name", sorted(_REQUIRED_TOOLS))
def test_chromatin_env_declares_its_tools(yaml_name):
    text = _dep_text(yaml_name).lower()
    for tool in _REQUIRED_TOOLS[yaml_name]:
        assert tool.lower() in text, (
            f"envs/{yaml_name} must declare '{tool}' (a tool its dispatch lane calls)"
        )


def test_atacseq_env_name_matches_environment_manager():
    # The atacseq stack must resolve to the env the YAML defines (S2/B2a wiring).
    from aria.utils.environment_manager import EnvironmentManager
    assert EnvironmentManager.STACKS.get("atacseq") == "aria-atacseq-env"
    assert EnvironmentManager.STACKS.get("tobias") == "aria-tobias-env"
    assert EnvironmentManager.STACKS.get("chromatin") == "aria-chromatin-env"
