#!/usr/bin/env python
"""One-time, explicit, ONLINE bootstrap of versioned TF motif collections.

ARIA's scATAC motif enrichment runs locally and offline (W-PRIV) against
versioned MEME-format motif files. This helper is the *only* place that talks to
the network for motifs: it downloads a JASPAR CORE MEME collection once and
writes it as ``<collection>.meme`` plus a ``manifest.json`` recording the source,
release, fetch date, SHA-256, and URL — so every analysis afterwards is
reproducible and egress-free.

    python scripts/fetch_motifs.py                       # JASPAR2024 CORE vertebrates
    python scripts/fetch_motifs.py --taxa plants
    ARIA_MOTIF_DIR=/data/motifs python scripts/fetch_motifs.py

The download is governed by ``aria.utils.privacy``: if ``ARIA_AIR_GAPPED`` is on,
it refuses egress and exits without writing anything (no fabrication).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from aria.utils.motifs import motifs_dir  # noqa: E402
from aria.utils import privacy  # noqa: E402

JASPAR_RELEASE = "2024"
# JASPAR CORE non-redundant PWMs, MEME format, by taxonomic group.
_URL_TMPL = (
    "https://jaspar.elixir.no/download/data/{release}/CORE/"
    "JASPAR{release}_CORE_{taxa}non-redundant_pfms_meme.txt"
)
# taxa token -> (collection name suffix, URL infix)
_TAXA = {
    "vertebrates": ("vertebrates", "vertebrates_"),
    "all": ("all", ""),
    "plants": ("plants", "plants_"),
    "insects": ("insects", "insects_"),
    "fungi": ("fungi", "fungi_"),
    "nematodes": ("nematodes", "nematodes_"),
}


def _collection_name(taxa: str) -> str:
    suffix, _ = _TAXA[taxa]
    return f"JASPAR{JASPAR_RELEASE}_CORE_{suffix}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxa", choices=sorted(_TAXA), default="vertebrates",
                    help="JASPAR CORE taxonomic group (default: vertebrates — "
                         "the correct set for human/mouse).")
    ap.add_argument("--timeout", type=int, default=120)
    args = ap.parse_args(argv)

    # W-PRIV: this is the single governed egress point for motifs.
    try:
        privacy.assert_egress_allowed("jaspar")
    except privacy.EgressBlocked as exc:
        print(f"[fetch_motifs] refusing download: {exc}", file=sys.stderr)
        return 2

    _, infix = _TAXA[args.taxa]
    url = _URL_TMPL.format(release=JASPAR_RELEASE, taxa=infix)
    collection = _collection_name(args.taxa)
    out_base = motifs_dir() / collection
    out_base.mkdir(parents=True, exist_ok=True)
    meme_path = out_base / f"{collection}.meme"

    print(f"[fetch_motifs] downloading {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ARIA-fetch-motifs"})
        with urllib.request.urlopen(req, timeout=args.timeout) as resp:
            data = resp.read()
    except Exception as exc:
        print(f"[fetch_motifs] download failed: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 1

    text = data.decode("utf-8", errors="ignore")
    n_motifs = sum(1 for ln in text.splitlines() if ln.startswith("MOTIF"))
    if n_motifs == 0:
        print("[fetch_motifs] downloaded file has no MOTIF records; not writing.",
              file=sys.stderr)
        return 1

    meme_path.write_text(text, encoding="utf-8")
    sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    manifest = {
        "collection": collection,
        "source": "JASPAR CORE (MEME format)",
        "release": JASPAR_RELEASE,
        "taxa": args.taxa,
        "date": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d"),
        "sha256": sha,
        "n_motifs": n_motifs,
        "url": url,
        "license": "CC0 1.0 (JASPAR)",
    }
    (out_base / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[fetch_motifs] wrote {meme_path} ({n_motifs} motifs)")
    print(f"[fetch_motifs] manifest: release={JASPAR_RELEASE} sha256={sha[:16]}…")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
