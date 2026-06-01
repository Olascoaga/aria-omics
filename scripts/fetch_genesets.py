#!/usr/bin/env python
"""One-time, explicit, ONLINE bootstrap of versioned GMT gene-set libraries.

ARIA's pathway ORA runs locally and offline (P1-7/W-PRIV) against versioned GMT
files. This helper is the *only* place that talks to the network for gene sets:
it downloads named Enrichr libraries once via ``gseapy.get_library`` and writes
each one as a ``.gmt`` plus a ``manifest.json`` recording the source, release,
fetch date, and SHA-256 — so every analysis afterwards is reproducible and
egress-free.

    python scripts/fetch_genesets.py                 # default human libraries
    python scripts/fetch_genesets.py --organism mouse
    python scripts/fetch_genesets.py --libraries GO_Biological_Process_2021 Reactome_2022
    ARIA_GMT_DIR=/data/genesets python scripts/fetch_genesets.py

Run inside an env with gseapy installed (e.g. aria-rna-env).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from aria.utils.ora import genesets_dir  # noqa: E402

# Defaults mirror aria/scripts/rna_bulk_de.py:_get_gene_sets so the local GMT
# library names match the labels the ORA engine asks for.
DEFAULT_LIBRARIES = {
    "human": [
        "GO_Biological_Process_2021",
        "KEGG_2021_Human",
        "Reactome_2022",
    ],
    "mouse": [
        "GO_Biological_Process_2021",
        "KEGG_2019_Mouse",
        "Reactome_2022",
    ],
}


def _release_from_name(name: str) -> str:
    """Best-effort release tag: the year embedded in the Enrichr library name."""
    m = re.search(r"(\d{4})", name)
    return m.group(1) if m else "unknown"


def _write_library(name: str, gene_sets: dict, out_dir: Path) -> dict:
    base = out_dir / name
    base.mkdir(parents=True, exist_ok=True)
    gmt_path = base / f"{name}.gmt"
    lines = []
    for term, genes in gene_sets.items():
        clean = [str(g).strip() for g in genes if g and str(g).strip()]
        if clean:
            lines.append("\t".join([str(term), "", *clean]))
    text = "\n".join(lines) + "\n"
    gmt_path.write_text(text, encoding="utf-8")

    sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    manifest = {
        "library": name,
        "source": "Enrichr (gseapy.get_library)",
        "release": _release_from_name(name),
        "date": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d"),
        "sha256": sha,
        "n_terms": len(lines),
        "url": f"https://maayanlab.cloud/Enrichr/geneSetLibrary?libraryName={name}",
    }
    (base / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", choices=sorted(DEFAULT_LIBRARIES),
                    default="human",
                    help="Which default library set to fetch (default: human).")
    ap.add_argument("--libraries", nargs="+", default=None,
                    help="Explicit Enrichr library names (overrides --organism).")
    ap.add_argument("--dir", default=None,
                    help="Destination dir (default: ARIA_GMT_DIR or ~/.aria/genesets).")
    args = ap.parse_args(argv)

    try:
        import gseapy as gp
    except ImportError:
        print("ERROR: gseapy is not installed. Run inside aria-rna-env.",
              file=sys.stderr)
        return 2

    out_dir = Path(args.dir).expanduser() if args.dir else genesets_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    libraries = args.libraries or DEFAULT_LIBRARIES[args.organism]

    print(f"Writing versioned GMT libraries to: {out_dir}")
    failures = 0
    for name in libraries:
        try:
            gene_sets = gp.get_library(name=name)
        except Exception as exc:  # network / unknown-library
            print(f"  [FAIL] {name}: {type(exc).__name__}: {exc}", file=sys.stderr)
            failures += 1
            continue
        if not gene_sets:
            print(f"  [FAIL] {name}: empty library returned", file=sys.stderr)
            failures += 1
            continue
        manifest = _write_library(name, gene_sets, out_dir)
        print(f"  [OK]   {name}: {manifest['n_terms']} terms "
              f"(release {manifest['release']}, sha256 {manifest['sha256'][:12]}…)")

    if failures:
        print(f"Done with {failures} failure(s).", file=sys.stderr)
        return 1
    print("Done. ORA will now run locally and offline against these libraries.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
