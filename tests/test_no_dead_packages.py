"""S6 (pre-integration audit): no dead package directories under aria/.

aria/core/ was an empty leftover (only a __pycache__, no source, untracked,
imported nowhere) from an old refactor. This fence keeps the package tree clean:
every directory under aria/ must carry at least one .py source file, so a dead
package (or a dir left holding only stale bytecode) cannot accumulate again.

Parses the tree only — runs in any env, clean checkout or local.
"""

from __future__ import annotations

from pathlib import Path

_ARIA = Path(__file__).resolve().parents[1] / "aria"


def _source_dirs():
    for d in _ARIA.rglob("*"):
        if not d.is_dir():
            continue
        if "__pycache__" in d.parts:
            continue
        yield d


def test_no_dir_under_aria_lacks_python_source():
    dead = []
    for d in _source_dirs():
        has_py = any(
            p.suffix == ".py"
            for p in d.iterdir()
            if p.is_file()
        )
        if not has_py:
            dead.append(str(d.relative_to(_ARIA.parent)))
    assert not dead, (
        "dead package directories (no .py source) under aria/: "
        + ", ".join(sorted(dead))
        + " — remove them (S6) or add the missing module"
    )


def test_aria_core_is_gone():
    # The specific leftover this slice removed must not return.
    assert not (_ARIA / "core").exists(), "aria/core/ reappeared (dead package)"
