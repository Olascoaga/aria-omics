"""
Pytest collection guard for legacy executable test scripts.

Most ARIA tests were originally written as standalone validation scripts with
top-level execution and sys.exit(). Importing them during pytest collection runs
the whole script and can abort collection. Keep them runnable directly, and use
pytest wrappers for the scripts we want in CI.
"""

collect_ignore = [
    # Converted to native pytest (P1-11 + follow-up) — collected directly:
    #   test_bulk_rna.py, test_scrna.py, test_environment_manager.py,
    #   test_pbmc_e2e.py.
    "test_chromatin_agent.py",
    "test_debate_council.py",
    "test_genome_arch_agent.py",
    "test_integration.py",
    "test_integration_agent.py",
    "test_narrative_agent.py",
    "test_scrna_e2e.py",
]
