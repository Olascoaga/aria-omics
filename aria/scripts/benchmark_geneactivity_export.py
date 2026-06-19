#!/usr/bin/env python3
"""scATAC P4.3 — export stage for the gene-activity concordance benchmark.

Runs in ``aria-chromatin-env``. Reads a 10x ARC peak-count matrix (or any ATAC
``.h5ad``/10x ``.h5`` with ``chr:start-end`` peaks), then:

  1. runs ARIA's own gene-activity scorer
     (``chromatin_regulatory._gene_scores``: peak accessibility summed in a window
     around each gene's TSS) + a local GTF -> ARIA per-gene activity CSV;
  2. writes the barcodes so the Signac ``GeneActivity`` reference driver scores the
     SAME cells from the fragments file.

No fabrication: peaks that do not parse to a coordinate are dropped (and were never
scorable by ARIA either). Nothing is downloaded.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True, help="10x ARC .h5 or ATAC .h5ad (peaks)")
    p.add_argument("--gtf", required=True, help="local gene-annotation GTF (chr* naming)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--window-bp", type=int, default=100_000)
    p.add_argument("--standard-chroms-only", action="store_true", default=True)
    args = p.parse_args(argv)

    import numpy as np
    import scanpy as sc

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if args.matrix.endswith(".h5"):
        ad = sc.read_10x_h5(args.matrix, gex_only=False)
        ad.var_names_make_unique()
        atac = ad[:, ad.var["feature_types"] == "Peaks"].copy()
        # ARC peak var_names are the gene-id column (e.g. chr1:9797-10686); fall back
        # to the interval column when present.
        if "interval" in atac.var:
            atac.var_names = [str(x) for x in atac.var["interval"]]
    else:
        atac = sc.read_h5ad(args.matrix)

    std = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    keep = []
    for name in map(str, atac.var_names):
        head = name.split(":")[0]
        if ":" in name and "-" in name and (not args.standard_chroms_only or head in std):
            keep.append(name)
    atac = atac[:, keep].copy()
    print(f"cells={atac.n_obs} peaks={atac.n_vars}", flush=True)

    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from aria.scripts.chromatin_regulatory import _gene_scores

    res = _gene_scores(atac, {"gtf_path": args.gtf,
                              "gene_score_window_bp": args.window_bp}, out)
    print(f"ARIA gene scores: ran={res.get('ran')} "
          f"n_genes={res.get('n_genes_scored')} reason={res.get('reason')}", flush=True)
    # _gene_scores writes chromatin_gene_scores.csv into out/.

    (out / "barcodes.tsv").write_text(
        "\n".join(map(str, atac.obs_names)) + "\n", encoding="utf-8")
    print(f"wrote export to {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
