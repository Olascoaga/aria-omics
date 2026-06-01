#!/usr/bin/env python3
"""Synchronize docs stamped from aria.version."""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.version import __version__, version_badge_url  # noqa: E402


def update_readme_badge() -> None:
    readme = ROOT / "README.md"
    text = readme.read_text(encoding="utf-8")
    pattern = re.compile(
        r"!\[Version\]\(https://img\.shields\.io/badge/version-[^)]+-blue\)"
    )
    replacement = f"![Version]({version_badge_url()})"
    new_text, n_replaced = pattern.subn(replacement, text, count=1)
    if n_replaced != 1:
        raise SystemExit("README version badge not found or ambiguous")
    readme.write_text(new_text, encoding="utf-8")


def validate_current_release_notes() -> None:
    release_notes = ROOT / "docs" / f"release_notes_v{__version__}.md"
    if not release_notes.exists():
        raise SystemExit(f"Missing release notes for aria.version {__version__}")
    title = release_notes.read_text(encoding="utf-8", errors="replace").splitlines()[0]
    if f"v{__version__}" not in title:
        raise SystemExit(
            f"Release notes title does not match aria.version {__version__}"
        )


def main() -> int:
    update_readme_badge()
    validate_current_release_notes()
    print(f"Synchronized version metadata from aria.version {__version__}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
