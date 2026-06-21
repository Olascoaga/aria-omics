"""
ARIA motif-collection resolver (v4.6 scATAC motif enrichment, W-PRIV).
---------------------------------------------------------------------
Transcription-factor motif enrichment runs LOCALLY and OFFLINE against a
versioned MEME-format motif collection, exactly like the pathway-ORA GMT story
(P1-7/W-PRIV). Motif files are read from ``ARIA_MOTIF_DIR`` (default
``~/.aria/motifs``):

    <motifs_dir>/<collection>/<collection>.meme    # MEME-format PWMs
    <motifs_dir>/<collection>/manifest.json        # {collection, source,
                                                    #  release, date, sha256, url}

The manifest is surfaced so ``methodology.json`` records exactly which motif
release produced the enrichment (reproducibility). This module deliberately has
NO heavy dependencies (no snapatac2) so it is importable anywhere for path
resolution and tests; the actual MEME parsing / scanning happens in
``aria/scripts/chromatin_motifs.py`` inside the chromatin env.

The ONLY place that talks to the network is ``scripts/fetch_motifs.py`` (a
one-time, explicit, governed bootstrap). At analysis time a missing collection
is an honest skip, never a fabricated enrichment (ADR-002).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from aria.utils.reference_integrity import (
    public_integrity_result, reference_is_usable, verify_reference_file,
)

MOTIF_DIR_ENV = "ARIA_MOTIF_DIR"

# Curated, CC0, semantically versioned, field-standard for ATAC reporting.
DEFAULT_COLLECTION = "JASPAR2024_CORE_vertebrates"

# Genome resolution lives in aria.utils.genomes (the full policy: explicit →
# env → managed store → governed auto-fetch). Re-exported here for backward
# compatibility — `ARIA_GENOME_FASTA` is now only the power-user env fallback.
from aria.utils.genomes import (  # noqa: E402,F401
    GENOME_FASTA_ENV, genome_fasta_from_env,
)


def motifs_dir() -> Path:
    """Resolve the versioned motif directory.

    ``ARIA_MOTIF_DIR`` overrides; otherwise ``$ARIA_HOME/motifs`` and finally
    ``~/.aria/motifs``."""
    override = os.environ.get(MOTIF_DIR_ENV)
    if override:
        return Path(override).expanduser()
    home = os.environ.get("ARIA_HOME")
    base = Path(home).expanduser() if home else (Path.home() / ".aria")
    return base / "motifs"


def collection_paths(collection: str) -> tuple[Path, Path]:
    """Return ``(meme_path, manifest_path)`` for a collection (may not exist)."""
    base = motifs_dir() / collection
    return base / f"{collection}.meme", base / "manifest.json"


def _count_meme_motifs(meme_path: Path) -> int:
    """Count ``MOTIF`` records in a MEME file (cheap line scan)."""
    n = 0
    try:
        with open(meme_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("MOTIF"):
                    n += 1
    except Exception:
        return 0
    return n


def load_local_motif_collection(collection: str = DEFAULT_COLLECTION):
    """Return ``(meme_path: str, version: dict)`` for a local versioned motif
    collection, or ``None`` when it is not staged.

    ``version`` is the parsed ``manifest.json`` (plus a counted ``n_motifs``), so
    the report can state the exact motif release used. No network, no snapatac2.
    """
    meme_path, manifest_path = collection_paths(collection)
    if not meme_path.is_file():
        return None
    n_motifs = _count_meme_motifs(meme_path)
    if n_motifs == 0:
        return None
    version: dict = {"collection": collection, "n_motifs": n_motifs}
    integrity = verify_reference_file(meme_path, manifest_path=manifest_path)
    if not reference_is_usable(integrity):
        return None
    if manifest_path.is_file():
        manifest_data = integrity.get("manifest")
        if isinstance(manifest_data, dict):
            version.update(manifest_data)
        else:
            try:
                version.update(json.loads(manifest_path.read_text(encoding="utf-8")))
            except Exception:
                pass
    version.setdefault("source", "unknown")
    version.setdefault("release", "unknown")
    version["n_motifs"] = n_motifs  # counted value is authoritative
    version["integrity"] = public_integrity_result(integrity)
    return str(meme_path), version
