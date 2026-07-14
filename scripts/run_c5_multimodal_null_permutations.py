#!/usr/bin/env python3
"""C5 multimodal label-permutation null lane.

Runs, for RNA and ATAC, a TRUE-label positive control plus many seeded label
permutations of a controlled matrix through the SHARED real DESeq2 core and
ARIA's real public-claim compiler, and writes the false-positive-narrative-rate
manifest (target ~0 per modality).

Environments (dispatched like the other multi-env lanes):
  * aria-chromatin-env — builds the controlled scATAC AnnData + pseudobulk peaks;
  * aria-rna-env       — the shared DESeq2 core for every RNA and ATAC contrast;
  * aria-env (this)    — narrative block building, the public compiler, scoring.

No fabrication: every verdict comes from real DESeq2 output through the real
compiler; a permutation that finds nothing yields an honest null claim.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aria.benchmarks import multimodal_null as mn  # noqa: E402


def _run_deseq2(counts_tsv: Path, meta_tsv: Path, out_json: Path,
                rna_env: str) -> dict:
    subprocess.run(
        ["conda", "run", "-n", rna_env, "python",
         str(ROOT / "scripts" / "aria_pseudobulk_da_from_tsv.py"),
         "--counts", str(counts_tsv), "--metadata", str(meta_tsv),
         "--numerator", "COND_B", "--denominator", "COND_A",
         "--output-json", str(out_json), "--min-replicates", "3"],
        check=True,
    )
    de = json.loads(out_json.read_text(encoding="utf-8"))
    if de.get("status") != "success":
        raise RuntimeError(f"DESeq2 did not succeed: {de.get('status')}")
    return de


def _null_runs(modality: str, counts_tsv: Path, metadata, work: Path,
               n_perms: int, seed: int, rna_env: str) -> tuple[dict, list[dict]]:
    import pandas as pd  # noqa: F401

    true_meta = work / f"{modality}_meta_true.tsv"
    metadata.to_csv(true_meta, sep="\t", index=False)
    true_de = _run_deseq2(counts_tsv, true_meta, work / f"{modality}_de_true.json",
                          rna_env)
    true_run = mn.classify_run(modality, true_de, is_permuted=False)

    perm_runs: list[dict] = []
    for i in range(n_perms):
        pmeta = mn.permute_conditions(metadata, seed=seed + i + 1)
        pmeta_tsv = work / f"{modality}_meta_perm{i}.tsv"
        pmeta.to_csv(pmeta_tsv, sep="\t", index=False)
        de = _run_deseq2(counts_tsv, pmeta_tsv,
                         work / f"{modality}_de_perm{i}.json", rna_env)
        perm_runs.append(mn.classify_run(modality, de, is_permuted=True))
    return true_run, perm_runs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="docs/benchmark_results/preprint_v1/claim_5/multimodal_null")
    parser.add_argument("--manifest-name", default="multimodal_null.json")
    parser.add_argument("--n-perms", type=int, default=20)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--work-dir", default="/tmp/c5_multimodal_null")
    parser.add_argument("--rna-env", default="aria-rna-env")
    parser.add_argument("--chromatin-env", default="aria-chromatin-env")
    args = parser.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)

    # RNA controlled matrix (this env) + its DESeq2 null permutations (rna-env).
    rna_counts, rna_meta = mn.synthesize_rna_counts(seed=args.seed)
    rna_counts_tsv = work / "rna_counts.tsv"
    rna_counts.to_csv(rna_counts_tsv, sep="\t", index=False)
    rna_true, rna_perms = _null_runs(
        "rna", rna_counts_tsv, rna_meta, work, args.n_perms, args.seed,
        args.rna_env)

    # ATAC controlled pseudobulk matrix (chromatin-env) + DESeq2 null (rna-env).
    atac_counts_tsv = work / "atac_counts.tsv"
    atac_meta_tsv = work / "atac_metadata.tsv"
    subprocess.run(
        ["conda", "run", "-n", args.chromatin_env, "python",
         str(ROOT / "scripts" / "aria_atac_pseudobulk_matrix.py"),
         "--out-counts", str(atac_counts_tsv),
         "--out-metadata", str(atac_meta_tsv), "--seed", str(args.seed)],
        check=True,
    )
    import pandas as pd
    atac_meta = pd.read_csv(atac_meta_tsv, sep="\t")
    atac_true, atac_perms = _null_runs(
        "atac", atac_counts_tsv, atac_meta, work, args.n_perms, args.seed,
        args.rna_env)

    modality_scores = {
        "rna": mn.score_modality("rna", rna_true, rna_perms),
        "atac": mn.score_modality("atac", atac_true, atac_perms),
    }
    manifest = mn.score_multimodal_null(modality_scores)

    from aria.version import __version__, collect_version_metadata
    manifest["aria_version"] = __version__
    manifest["provenance"] = collect_version_metadata()
    manifest["config"] = {"n_permutations": args.n_perms, "seed": args.seed}

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / args.manifest_name
    out_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    for mod, s in modality_scores.items():
        print(f"{mod}: positive_control="
              f"{s['positive_control_detects_signal']} fp_rate="
              f"{s['false_positive_narrative_rate']} "
              f"(max null n_sig={s['max_null_n_sig']})")
    print(f"status={manifest['status']} -> wrote {out_path}")
    return 0 if manifest["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
