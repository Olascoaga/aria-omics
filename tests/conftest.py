"""
Pytest collection guard for legacy executable test scripts.

Most ARIA tests were originally written as standalone validation scripts with
top-level execution and sys.exit(). Importing them during pytest collection runs
the whole script and can abort collection. Keep them runnable directly, and use
pytest wrappers for the scripts we want in CI.
"""

import os

import pytest

collect_ignore = [
    # Converted to native pytest (P1-11 + follow-up + D3) — collected directly:
    #   test_bulk_rna.py, test_scrna.py, test_environment_manager.py,
    #   test_pbmc_e2e.py, test_integration_agent.py.
    "test_chromatin_agent.py",
    "test_debate_council.py",
    "test_genome_arch_agent.py",
    "test_integration.py",
    "test_narrative_agent.py",
    "test_scrna_e2e.py",
]

_STACK_STATUS_CACHE: dict[str, dict[str, bool]] = {}


def require_aria_stack(stack: str) -> None:
    """Skip the current test if the named ARIA conda stack is not ready."""
    from aria.utils.environment_manager import EnvironmentManager

    cache_key = os.environ.get("PATH", "")
    envs = _STACK_STATUS_CACHE.get(cache_key)
    if envs is None:
        envs = EnvironmentManager().check_environments()
        _STACK_STATUS_CACHE[cache_key] = envs
    if not envs.get(stack, False):
        env_name = EnvironmentManager.STACKS.get(stack, stack)
        pytest.skip(
            f"ARIA stack '{stack}' (conda env '{env_name}') is not available"
        )


@pytest.fixture(name="require_aria_stack")
def _require_aria_stack_fixture():
    return require_aria_stack
