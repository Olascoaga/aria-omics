#!/usr/bin/env python3
"""Run ARIA's real pseudobulk DESeq2 (rna_bulk_de._run_deseq2 — the same core the
scATAC pseudobulk DA lane uses) on a counts/metadata TSV pair, and write a compact
DA JSON. This is the ARIA-side mirror of benchmark_a1_external_comparators.R: same
``gene``/``sample``/``condition`` TSV schema, so the orchestrator can run ARIA and
the R comparators on the IDENTICAL matrix and score concordance.

Runs in the RNA stack (e.g. ``conda run -n aria-rna-env``) which has pydeseq2.
Output JSON: {status, n_tested, n_sig, n_up, n_down, sig_peaks, lfc_by_peak,
padj_by_peak}. No fabrication: a non-success DESeq2 result is written as-is.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True, help="TSV with a 'gene' column + one column per sample")
    p.add_argument("--metadata", required=True, help="TSV with 'sample' + 'condition' columns")
    p.add_argument("--output-json", required=True)
    p.add_argument("--numerator", default="COND_B", help="condition level tested vs --denominator")
    p.add_argument("--denominator", default="COND_A", help="reference condition level")
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=0.5)
    p.add_argument("--min-replicates", type=int, default=3)
    args = p.parse_args(argv)

    import pandas as pd
    from aria.scripts.rna_bulk_de import _run_deseq2

    counts = pd.read_csv(args.counts, sep="\t")
    if "gene" not in counts.columns:
        raise SystemExit("counts TSV must contain a 'gene' column")
    counts = counts.set_index("gene")
    counts.index = counts.index.astype(str)
    meta = pd.read_csv(args.metadata, sep="\t").set_index("sample")
    meta = meta.loc[list(counts.columns)]

    res, warns = _run_deseq2(
        counts, meta, design_factor="condition",
        numerator=args.numerator, denominator=args.denominator,
        padj_thr=args.padj, lfc_thr=args.lfc,
        min_replicates_per_condition=args.min_replicates, lfc_shrink=True)

    out: dict = {"status": res.get("status"), "warnings": list(warns)}
    if res.get("status") == "success":
        df = res["results"]
        df.index = df.index.astype(str)
        out.update({
            "n_tested": int(len(df)),
            "n_sig": int(res.get("n_sig", 0)),
            "n_up": int(res.get("n_up", 0)),
            "n_down": int(res.get("n_down", 0)),
            "sig_peaks": [str(g) for g in res.get("sig_genes", [])],
            "lfc_by_peak": {str(g): float(v) for g, v in df["log2FoldChange"].items()},
            "padj_by_peak": {str(g): float(v) for g, v in df["padj"].items()},
        })
    else:
        out.update({"error_type": res.get("error_type"), "details": res.get("details")})

    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_json).write_text(json.dumps(out), encoding="utf-8")
    print(f"ARIA pseudobulk DA: status={out['status']} "
          f"n_sig={out.get('n_sig')} (of {out.get('n_tested')} tested)")
    return 0 if out["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
