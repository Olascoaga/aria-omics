"""F9 (preprint audit 2026-06-19): persist a public artifact CSV, disclosing
write failures instead of swallowing them.

Public DA / motif result CSVs were written under a bare ``except Exception:
path = None``, so a failed write became an invisible downstream skip — a consumer
(figure builder, narrator, benchmark) just saw no table with no reason. ``persist_csv``
keeps the honest-skip semantics (it still returns ``None`` on failure, never a
fabricated path) but appends a structured ``CsvWriteFailed`` warning so the failure is
visible in the result's ``warnings``.
"""
from __future__ import annotations

from typing import Callable, Optional


def persist_csv(path: str, label: str, write_fn: Callable[[], object],
                warnings: list) -> Optional[str]:
    """Run ``write_fn`` to write ``path``; on failure record a warning, return None.

    ``label`` is the human-facing artifact name used in the warning. ``write_fn`` is
    the actual write (e.g. ``lambda: df.to_csv(path, index=False)``) so the CSV format
    of each call site is unchanged. Returns ``path`` on success, ``None`` on failure.
    """
    try:
        write_fn()
        return path
    except Exception as e:  # noqa: BLE001 — disclosed via a warning, not swallowed
        warnings.append(
            f"CsvWriteFailed: could not write {label} ({type(e).__name__}: "
            f"{str(e)[:200]}); the table is unavailable to downstream consumers.")
        return None
