#!/usr/bin/env python3
"""One-time bootstrap of the SEQC/MAQC A1 reference bundle.

Runs the R extractor in ``aria-bench-env`` (needs the ``seqc`` Bioconductor data
package installed there) and writes a versioned bundle with provenance hashes:

    counts.tsv  samples.tsv  taqman.tsv  manifest.json

Then point the benchmark at it:

    python scripts/fetch_seqc_maqc_reference.py --out ~/.aria/benchmarks/seqc_maqc
    ARIA_SEQC_MAQC_BUNDLE=~/.aria/benchmarks/seqc_maqc \\
        conda run -n aria-rna-env python scripts/run_a1_seqc_maqc_benchmark.py

This is a governed online bootstrap (like scripts/fetch_genesets.py); the bundle
is reference benchmark data, kept out of the repo. ARIA fabricates nothing.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
R_SCRIPT = ROOT / "aria" / "scripts" / "fetch_seqc_maqc_reference.R"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default=os.path.expanduser("~/.aria/benchmarks/seqc_maqc"),
        help="Output bundle directory (kept out of the repo).",
    )
    parser.add_argument("--env", default="aria-bench-env")
    parser.add_argument("--count-table", default="ILM_refseq_gene_BGI")
    args = parser.parse_args(argv)

    out = Path(os.path.expanduser(args.out))
    out.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "--name", args.env, "--no-capture-output",
        "Rscript", str(R_SCRIPT), str(out), args.count_table,
    ]
    print(f"[fetch] {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0 or "SEQC_BUNDLE_DONE" not in proc.stdout:
        sys.stderr.write(proc.stderr[-4000:])
        print("[fetch] FAILED — is the `seqc` package installed in the env?")
        return 1

    files = ["counts.tsv", "samples.tsv", "taqman.tsv"]
    missing = [f for f in files if not (out / f).exists()]
    if missing:
        print(f"[fetch] FAILED — extractor did not write: {missing}")
        return 1

    manifest = {
        "source": "seqc Bioconductor data package",
        "references": [
            "SEQC/MAQC-III Consortium, Nature Biotechnology 2014",
            "MAQC Consortium (TaqMan), Nature Biotechnology 2006",
        ],
        "count_table": args.count_table,
        "truth": "TaqMan qPCR log2(A/B) by gene Symbol (mean of 4 reps each)",
        "fetched_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "env": args.env,
        "sha256": {f: _sha256(out / f) for f in files},
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(f"[fetch] wrote bundle to {out}")
    print(f"[fetch] sha256 counts.tsv = {manifest['sha256']['counts.tsv'][:16]}…")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
