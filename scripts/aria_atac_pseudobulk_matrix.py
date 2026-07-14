#!/usr/bin/env python3
"""Build a controlled scATAC AnnData and aggregate it to a pseudobulk peak matrix.

Runs inside aria-chromatin-env (anndata/snapatac2 stack). Emits a pseudobulk
peak x replicate counts TSV (``gene`` column of peak ids + one column per
replicate) plus a metadata TSV (``sample``/``condition``), the exact shape the
shared DESeq2 core consumes in aria-rna-env. This is the ATAC arm of the C5
multimodal label-permutation null: the matrix carries a real condition signal so
the TRUE-label run is a positive control, and label permutation destroys it.

No fabrication of a governance verdict: this only prepares the controlled matrix;
the differential test and its null permutations run downstream through the real
DESeq2 core and the real public-claim compiler.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-counts", required=True)
    parser.add_argument("--out-metadata", required=True)
    parser.add_argument("--n-peaks", type=int, default=120)
    parser.add_argument("--n-reps-per-cond", type=int, default=6)
    parser.add_argument("--cells-per-rep", type=int, default=200)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args(argv)

    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(args.seed)
    conditions = ["COND_A", "COND_B"]
    n_de = max(8, args.n_peaks // 6)
    peaks = [f"chr1:{1000 + i * 500}-{1000 + i * 500 + 300}"
             for i in range(args.n_peaks)]

    rows = []
    obs_cond, obs_rep = [], []
    base = rng.integers(3, 20, size=args.n_peaks).astype(float)
    for cond in conditions:
        for r in range(args.n_reps_per_cond):
            rep = f"{cond}_rep{r + 1}"
            for _ in range(args.cells_per_rep):
                means = base.copy()
                if cond == "COND_B":
                    means[:n_de] *= 3.0            # real accessibility signal
                rows.append(rng.poisson(np.maximum(means, 0.2)))
                obs_cond.append(cond)
                obs_rep.append(rep)
    X = np.asarray(rows, dtype="float32")
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"condition": obs_cond, "replicate": obs_rep,
                          "cell_type": "ctype0"}),
        var=pd.DataFrame(index=peaks),
    )

    # Pseudobulk: sum peak counts per biological replicate (the only defensible
    # cross-condition unit; single cells are not replicates).
    reps = list(dict.fromkeys(adata.obs["replicate"]))
    pb = np.zeros((adata.n_vars, len(reps)), dtype=int)
    rep_to_cond = {}
    for j, rep in enumerate(reps):
        mask = (adata.obs["replicate"] == rep).to_numpy()
        pb[:, j] = np.rint(np.asarray(adata.X[mask].sum(axis=0)).ravel())
        rep_to_cond[rep] = adata.obs.loc[mask, "condition"].iloc[0]

    counts = pd.DataFrame(pb, index=peaks, columns=reps)
    counts.insert(0, "gene", counts.index)
    Path(args.out_counts).parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(args.out_counts, sep="\t", index=False)
    pd.DataFrame({"sample": reps,
                  "condition": [rep_to_cond[r] for r in reps]}).to_csv(
        args.out_metadata, sep="\t", index=False)
    print(f"ATAC pseudobulk: {adata.n_obs} cells -> {len(reps)} replicate "
          f"pseudobulk columns x {adata.n_vars} peaks -> {args.out_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
