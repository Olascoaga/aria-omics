# ARIA packaging moved to PEP 621 metadata in pyproject.toml (P2-1).
# This shim remains only so legacy `python setup.py` invocations and editable
# installs on older tooling keep working; all metadata (name, dynamic version
# from aria.version.__version__, dependencies + ceilings, extras, scripts,
# classifiers) lives in pyproject.toml — do NOT duplicate it here.
from setuptools import setup

setup()
