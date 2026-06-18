#!/usr/bin/env python3
"""scATAC P4.2 — export stage for the peak-to-gene concordance benchmark.

Runs in ``aria-chromatin-env`` (muon/anndata). Reads a paired RNA+ATAC ``.h5mu``,
then:

  1. runs ARIA's own peak-to-gene linker (``chromatin_regulatory._peak_to_gene``)
     on the SAME paired cells + local GTF and writes the ARIA links CSV;
  2. exports the ATAC peak-count matrix and the RNA count matrix (+ peaks / genes /
     barcodes) so the Signac ``LinkPeaks`` reference driver consumes the IDENTICAL
     cells, peaks, and genes.

No fabrication: peaks that do not parse to a coordinate are dropped from the export
(and were never linkable by ARIA either). Nothing is downloaded.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

_PEAK_RE = re.compile(r"^([^:\-_]+)[:\-_](\d+)[\-_](\d+)$")


def _parse_peak(name: str):
    m = _PEAK_RE.match(str(name).strip())
    if not m:
        return None
    s, e = int(m.group(2)), int(m.group(3))
    return (m.group(1), min(s, e), max(s, e))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5mu", required=True)
    p.add_argument("--gtf", required=True, help="local gene-annotation GTF (matches the peak naming, e.g. UCSC chr*)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--atac-mod", default="atac")
    p.add_argument("--rna-mod", default="rna")
    p.add_argument("--distance-bp", type=int, default=500_000)
    p.add_argument("--min-corr", type=float, default=0.3)
    p.add_argument("--min-shared-cells", type=int, default=100)
    p.add_argument("--standard-chroms-only", action="store_true", default=True)
    args = p.parse_args(argv)

    import muon as mu
    import numpy as np
    import scipy.io as sio
    import scipy.sparse as sp

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    mdata = mu.read(args.h5mu)
    atac = mdata.mod[args.atac_mod]
    rna = mdata.mod[args.rna_mod]

    # Align both modalities to the shared barcodes (paired by construction).
    shared = [c for c in atac.obs_names if c in set(rna.obs_names)]
    atac = atac[shared, :].copy()
    rna = rna[shared, :].copy()

    # Drop peaks that do not parse to a coordinate (and optionally non-standard
    # chromosomes) so ARIA and Signac see exactly the same peak universe.
    std = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    keep_idx, peak_rows = [], []
    for i, name in enumerate(map(str, atac.var_names)):
        parsed = _parse_peak(name)
        if parsed is None:
            continue
        chrom, start, end = parsed
        if args.standard_chroms_only and chrom not in std:
            continue
        keep_idx.append(i)
        peak_rows.append((name, chrom, start, end))
    atac = atac[:, keep_idx].copy()
    print(f"shared cells={len(shared)} peaks={len(peak_rows)} genes={rna.n_vars}", flush=True)

    # 1) ARIA peak-to-gene on these paired cells.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from aria.scripts.chromatin_regulatory import _peak_to_gene

    rna_h5ad = out / "rna.h5ad"
    rna.write_h5ad(rna_h5ad)
    aria_res = _peak_to_gene(
        atac,
        {
            "rna_data_path": str(rna_h5ad),
            "gtf_path": args.gtf,
            "peak_gene_distance_bp": args.distance_bp,
            "min_peak_gene_corr": args.min_corr,
            "min_shared_cells": args.min_shared_cells,
        },
        out,
    )
    print(f"ARIA peak2gene: ran={aria_res.get('ran')} "
          f"n_links={aria_res.get('n_links')} n_tested={aria_res.get('n_tested')} "
          f"reason={aria_res.get('reason')}", flush=True)
    # _peak_to_gene writes chromatin_peak_to_gene.csv into out/.

    # 2) Export matrices for Signac (features x cells, MatrixMarket).
    def _csc(mat):
        return sp.csc_matrix(mat) if not sp.issparse(mat) else mat.tocsc()

    sio.mmwrite(str(out / "atac_counts.mtx"), _csc(atac.X).T)   # peaks x cells
    sio.mmwrite(str(out / "rna_counts.mtx"), _csc(rna.X).T)     # genes x cells
    (out / "peaks.tsv").write_text(
        "\n".join(f"{c}\t{s}\t{e}\t{n}" for (n, c, s, e) in peak_rows) + "\n")
    (out / "genes.tsv").write_text("\n".join(map(str, rna.var_names)) + "\n")
    (out / "barcodes.tsv").write_text("\n".join(map(str, atac.obs_names)) + "\n")
    print(f"wrote export to {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
