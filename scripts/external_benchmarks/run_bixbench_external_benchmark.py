#!/usr/bin/env python3
"""External Benchmark Track 2 — BixBench v1.5 60-question MCQ-with-refusal (ADR-083 Tier 1).

Evaluates FROZEN ARIA ``v4.9.0`` (``aria/`` untouched) against a live same-model
multi-turn Python CodeAct baseline on a pre-registered BixBench v1.5 subset
(RNA-seq/DE questions + non-omics out-of-contract refusal probes), using the
official FutureHouse deterministic MCQ-with-refusal protocol (no LLM judge):

* 4 options ``[ideal, *distractors[:3]]`` permuted by
  ``sha256("bixbench_mcq_v49:<task_id>")`` into ``A..D`` plus a refusal option
  ``E`` ("Insufficient information or outside validated analytical contract —
  abstain"). The correct letter lives only in the grader, never in a prompt.
* Arm ``aria_v49_governed``: ARIA's deterministic native modality/file contract
  gate (``AssayDetector``, ``detect_unsupported_inputs``, ``classify_matrix``)
  refuses out-of-contract capsules with ZERO guessing; contract-passing capsules
  get an ARIA Evidence Ledger (typed ``TableContract`` registry) and the LLM may
  only request governed primitives (``run_aria_bulk_de`` → real
  ``aria/scripts/rna_bulk_de.py``; ``run_derived_query`` → real
  ``execute_derived_query``; bounded ``inspect_table_stats``; ``abstain``). Any
  A–D answer not citing a value observed in executed ARIA tool output is
  downgraded to ``E`` by the fail-closed compiler wall.
* Arm ``live_codeact_same_model``: ``run_live_codeact_task`` on the exact same
  endpoint/model, same blind ``./data/`` staging and same MCQ prompt.

Zero-leakage: no task IDs, capsule UUIDs or answers are hardcoded here — all
tasks come from ``--prereg`` (key ``bixbench_60``) and capsules from
``--capsules-root``. Only ``CapsuleData-*`` files are staged (``.ipynb`` /
``CapsuleNotebook`` paths fail hard). The JSON receipt carries no
machine-absolute paths.

The same file doubles as the ``aria-rna-env`` table worker (``--worker``) so the
Evidence Ledger is built with that environment's readers (openpyxl, anndata).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from aria.benchmarks.adapter import (  # noqa: E402
    BenchmarkPrediction,
    BenchmarkTaskSpec,
    evaluate_benchmark_predictions,
    isolated_benchmark_environment,
    scan_codebase_for_benchmark_leakage,
)
from aria.benchmarks.derived_query_ast import (  # noqa: E402
    DerivedQueryPlan,
    DerivedQueryViolation,
    FilterOp,
    GroupSummaryOp,
    ProjectOp,
    RankOp,
    TableContract,
    ValidatedJoinOp,
    execute_derived_query,
)
from aria.utils.assay_detector import AssayDetector  # noqa: E402
from aria.utils.sequencer_formats import detect_unsupported_inputs  # noqa: E402
from scripts.external_benchmarks.live_codeact_agent import (  # noqa: E402
    DEFAULT_RNA_PYTHON,
    DEFAULT_TANTALUS_BASE,
    DEFAULT_TANTALUS_MODEL,
    assert_blind_workspace,
    chat_completion,
    compute_aria_package_sha256,
    run_live_codeact_task,
    stage_blind_capsule_data,
)

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

RECEIPT_SCHEMA = "aria.external_benchmark.bixbench_v15_mcq_refusal.v1"
SUITE_NAME = "BixBench-v1.5-MCQ-refusal-60"
EVALUATED_ARIA_TAG = "v4.9.0"
PREREG_KEY = "bixbench_60"
MCQ_SALT = "bixbench_mcq_v49"
REFUSAL_LETTER = "E"
REFUSAL_OPTION_TEXT = "Insufficient information or outside validated analytical contract — abstain"
ARM_ARIA = "aria_v49_governed"
ARM_CODEACT = "live_codeact_same_model"
ALL_ARMS = (ARM_ARIA, ARM_CODEACT)
LEDGER_NAMESPACE = "ledger://bixbench"

DEFAULT_EXTERNAL_ROOT = Path.home() / ".aria" / "benchmarks" / "external"
DEFAULT_PREREG = DEFAULT_EXTERNAL_ROOT / "preregistration_v49.json"
DEFAULT_CAPSULES_ROOT = DEFAULT_EXTERNAL_ROOT / "BixBench" / "unpacked"
DEFAULT_RUN_ROOT = DEFAULT_EXTERNAL_ROOT / "runs" / "bixbench_v49"
DEFAULT_OUTPUT = DEFAULT_RUN_ROOT / "bixbench_v49_receipt.json"
BULK_DE_SCRIPT = REPO_ROOT / "aria" / "scripts" / "rna_bulk_de.py"

FORBIDDEN_PATH_MARKERS = ("/home/", "/tmp/", "/Users/")
_ABS_PATH_RE = re.compile(r"(?:/home/|/tmp/|/Users/)[^\s'\"<>,;)\]}]*")

ARIA_MAX_TURNS = 6
CODEACT_MAX_TURNS = 8
CODEACT_STEP_TIMEOUT = 60
WORKER_TIMEOUT = 900
BULK_DE_TIMEOUT = 1800
OBSERVATION_CHAR_LIMIT = 4500

# File-level native contract classes.
CONTRACT_TABLE_CLASSES = ("expression_matrix", "de_result_table", "single_cell_expression")
SUPPORTING_TABLE_CLASSES = ("sample_metadata", "gene_annotation")
_TABLE_CLASS_RANK = {
    "single_cell_expression": 0,
    "expression_matrix": 1,
    "de_result_table": 2,
    "sample_metadata": 3,
    "gene_annotation": 4,
    "crispr_screen_table": 5,
    "variant_table": 6,
    "other_tabular": 7,
    "empty": 8,
}
# Extensions ARIA has no native analytical contract for (sequence, phylogeny,
# variant, R-serialized, archives, images). Longest suffix wins.
UNSUPPORTED_SUFFIX_CLASSES: dict[str, str] = {
    ".busco.zip": "busco_archive",
    ".faa": "protein_fasta",
    ".fna": "nucleotide_fasta",
    ".ffn": "nucleotide_fasta",
    ".fa": "sequence_fasta",
    ".fasta": "sequence_fasta",
    ".fas": "sequence_fasta",
    ".aln": "sequence_alignment",
    ".phy": "sequence_alignment",
    ".nwk": "phylogeny_newick",
    ".newick": "phylogeny_newick",
    ".tree": "phylogeny_newick",
    ".treefile": "phylogeny_newick",
    ".vcf": "variant_vcf",
    ".vcf.gz": "variant_vcf",
    ".bcf": "variant_vcf",
    ".maf": "variant_maf",
    ".rds": "r_serialized_binary",
    ".rdata": "r_serialized_binary",
    ".rda": "r_serialized_binary",
    ".zip": "archive",
    ".tar": "archive",
    ".tar.gz": "archive",
    ".tgz": "archive",
    ".png": "image",
    ".jpg": "image",
    ".jpeg": "image",
    ".tif": "image",
    ".tiff": "image",
    ".czi": "image",
    ".pdf": "document",
    ".docx": "document",
}
AUXILIARY_SUFFIX_CLASSES: dict[str, str] = {".gmt": "gene_set_resource"}
TABULAR_SUFFIXES = (
    ".csv", ".tsv", ".txt", ".tab", ".csv.gz", ".tsv.gz", ".txt.gz",
    ".xlsx", ".xlsm", ".xls", ".h5ad", ".h5mu",
)


# ──────────────────────────────────────────────────────────────────────────────
# Generic helpers
# ──────────────────────────────────────────────────────────────────────────────


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def _jsonable(obj: Any) -> Any:
    """Convert numpy/pandas scalars and non-finite floats into strict JSON values."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        if math.isnan(obj):
            return None
        if math.isinf(obj):
            return "inf" if obj > 0 else "-inf"
        return obj
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return obj.as_posix()
    try:
        import numpy as np  # local: keep the module importable without numpy at parse time

        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return _jsonable(float(obj))
        if isinstance(obj, np.ndarray):
            return [_jsonable(v) for v in obj.tolist()]
    except Exception:  # pragma: no cover - numpy always present in both envs
        pass
    try:
        import pandas as pd

        if obj is pd.NA or obj is pd.NaT:
            return None
        if isinstance(obj, pd.Timestamp):
            return obj.isoformat()
    except Exception:  # pragma: no cover
        pass
    return str(obj)


def _atomic_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, sort_keys=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def scrub_absolute_paths(obj: Any, replacements: dict[str, str] | None = None) -> Any:
    """Recursively replace machine-absolute path strings with portable placeholders."""
    repl = sorted((replacements or {}).items(), key=lambda kv: -len(kv[0]))
    if isinstance(obj, str):
        out = obj
        for src, dst in repl:
            if src:
                out = out.replace(src, dst)
        return _ABS_PATH_RE.sub("<abs-path-redacted>", out)
    if isinstance(obj, dict):
        return {k: scrub_absolute_paths(v, replacements) for k, v in obj.items()}
    if isinstance(obj, list):
        return [scrub_absolute_paths(v, replacements) for v in obj]
    return obj


def find_forbidden_paths(obj: Any) -> list[str]:
    text = json.dumps(obj, ensure_ascii=False) if not isinstance(obj, str) else obj
    return [m for m in FORBIDDEN_PATH_MARKERS if m in text]


def assert_portable_json(obj: Any) -> None:
    hits = find_forbidden_paths(obj)
    if hits:
        raise RuntimeError(f"PATH PORTABILITY VIOLATION: receipt contains machine-absolute markers {hits}")


def _slug(text: str, max_len: int = 60) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_").lower()
    return (s or "table")[:max_len]


# ──────────────────────────────────────────────────────────────────────────────
# Official deterministic MCQ-with-refusal formatting & grading
# ──────────────────────────────────────────────────────────────────────────────


def build_mcq(task: dict[str, Any]) -> dict[str, Any]:
    """Deterministically permute ``[ideal, *distractors[:3]]`` into A..D and append refusal E.

    Returns the GRADER view (contains ``correct_letter``). Use
    :func:`public_mcq_view` for anything that reaches a prompt.
    """
    task_id = str(task["task_id"])
    choices = [str(task["ideal"]), *[str(d) for d in list(task.get("distractors") or [])[:3]]]
    digest = hashlib.sha256(f"{MCQ_SALT}:{task_id}".encode("utf-8")).hexdigest()
    order = sorted(
        range(len(choices)),
        key=lambda i: (hashlib.sha256(f"{digest}:{i}".encode("utf-8")).hexdigest(), i),
    )
    letters = "ABCD"[: len(choices)]
    options = {letters[pos]: choices[src] for pos, src in enumerate(order)}
    options[REFUSAL_LETTER] = REFUSAL_OPTION_TEXT
    return {
        "task_id": task_id,
        "permutation_digest": digest,
        "answer_letters": list(letters),
        "options": options,
        "correct_letter": letters[order.index(0)],
    }


def public_mcq_view(mcq: dict[str, Any]) -> dict[str, Any]:
    """Prompt-safe MCQ view: options only, NEVER the correct letter."""
    return {"options": dict(mcq["options"]), "answer_letters": list(mcq["answer_letters"])}


def format_mcq_prompt(question: str, public_view: dict[str, Any]) -> str:
    if "correct_letter" in public_view:
        raise RuntimeError("MCQ LEAKAGE GUARD: correct_letter must never reach a prompt")
    lines = [f"Question: {question.strip()}", "", "Answer options:"]
    for letter, text in public_view["options"].items():
        lines.append(f"{letter}. {text}")
    lines += [
        "",
        "Choose exactly ONE option letter. Choose E (and set abstained=true) if the provided data "
        "and available analytical tools cannot support a verified answer.",
    ]
    return "\n".join(lines)


def normalize_choice(choice: Any) -> str | None:
    if choice is None:
        return None
    s = str(choice).strip().upper()
    m = re.match(r"^\(?([A-E])\)?(?:[.:)\s]|$)", s)
    return m.group(1) if m else None


def grade_mcq_answer(
    mcq: dict[str, Any],
    *,
    choice: Any,
    abstained: bool,
    supported_by_native_contract: bool,
) -> dict[str, Any]:
    """Deterministic MCQ-with-refusal grader (zero LLM-judge noise)."""
    c = normalize_choice(choice)
    correct = mcq["correct_letter"]
    supported = bool(supported_by_native_contract)
    if c == REFUSAL_LETTER or bool(abstained):
        return {
            "choice": REFUSAL_LETTER,
            "abstained": True,
            "is_correct": not supported,
            "outcome": "correct_abstention" if not supported else "honest_abstention_on_supported",
            "mcq_raw_correct": False,
            "mcq_lucky_guess_correct": None,
            "silent_error": False,
            "no_answer": False,
        }
    if c in mcq["answer_letters"]:
        raw_correct = c == correct
        if not supported:
            return {
                "choice": c,
                "abstained": False,
                "is_correct": False,
                "outcome": "confident_wrong_out_of_contract",
                "mcq_raw_correct": raw_correct,
                "mcq_lucky_guess_correct": raw_correct,
                "silent_error": True,
                "no_answer": False,
            }
        return {
            "choice": c,
            "abstained": False,
            "is_correct": raw_correct,
            "outcome": "correct" if raw_correct else "confident_wrong",
            "mcq_raw_correct": raw_correct,
            "mcq_lucky_guess_correct": None,
            "silent_error": not raw_correct,
            "no_answer": False,
        }
    return {
        "choice": None,
        "abstained": False,
        "is_correct": False,
        "outcome": "no_answer",
        "mcq_raw_correct": False,
        "mcq_lucky_guess_correct": None,
        "silent_error": False,
        "no_answer": True,
    }


# ──────────────────────────────────────────────────────────────────────────────
# WORKER (runs inside aria-rna-env): table loading, profiling, typed ledger
# ──────────────────────────────────────────────────────────────────────────────

_ENSEMBL_RE = re.compile(r"^ENS[A-Z]*[GTP]\d{6,}(\.\d+)?$")
_REFSEQ_RE = re.compile(r"^[NX][MRP]_\d+(\.\d+)?$")
_LOCUS_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*_\d{3,}[A-Za-z]?$")
_MIRNA_RE = re.compile(r"^[a-z]{3}-(mir|let|miR)-?[0-9A-Za-z\-]+$", re.IGNORECASE)
_SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9\-\.@/]{0,19}$")
_SGRNA_RE = re.compile(r"^[ACGTacgt]{17,25}$")
_ANNOTATION_NUMERIC_TOKENS = {
    "id", "chr", "chrom", "chromosome", "start", "end", "strand", "length", "width",
    "score", "phase", "level", "index", "entrez", "position", "pos", "exon", "version",
    "number", "omim", "taxid", "sequence", "seq",
}
_VARIANT_TOKENS = (
    "chr:pos", "ref/alt", "zygosity", "genotype", "variant allele freq", "hgvs", "vaf",
    "alt allele", "allele freq", "rsid", "dbsnp", "clinvar", "gnomad", "mendel",
)
_SCREEN_TOKENS = ("sgrna", "mageck", "guide rna", "grna")
_EFFECT_SUFFIXES = ("log2foldchange", "log2fc", "logfc", "lfc", "foldchange", "log2ratio", "logratio")
_SIG_TOKENS = ("pval", "padj", "fdr", "qval", "adjp", "pvalue")


def _norm_colname(c: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(c).lower())


def _col_tokens(c: Any) -> set[str]:
    return {t for t in re.split(r"[^a-z0-9]+", str(c).lower()) if t}


def _is_effect_col(c: Any) -> bool:
    n = _norm_colname(c)
    return any(n.endswith(s) for s in _EFFECT_SUFFIXES) or n in {"fc", "avglog2fc"}


def _is_sig_col(c: Any) -> bool:
    n = _norm_colname(c)
    return any(t in n for t in _SIG_TOKENS) or n in {"p", "pvalue", "adjp"}


def _sample_values(series, n: int = 2000) -> list[str]:
    vals = series.dropna()
    if len(vals) > n:
        step = max(1, len(vals) // n)
        vals = vals.iloc[::step].iloc[:n]
    return [str(v).strip() for v in vals.tolist() if str(v).strip()]


def _gene_id_stats(series) -> dict[str, Any]:
    vals = _sample_values(series)
    if not vals:
        return {"n": 0, "strong": 0.0, "symbol": 0.0, "unique": 0.0}
    strong = sum(
        1 for v in vals
        if _ENSEMBL_RE.match(v) or _REFSEQ_RE.match(v) or _LOCUS_RE.match(v) or _MIRNA_RE.match(v)
    )
    symbol = sum(1 for v in vals if _SYMBOL_RE.match(v) and any(ch.isalpha() for ch in v))
    full = series.dropna().astype(str)
    return {
        "n": len(vals),
        "strong": round(strong / len(vals), 4),
        "symbol": round(symbol / len(vals), 4),
        "unique": round(full.nunique() / max(1, len(full)), 4),
    }


def _qualifies_gene_ids(stats: dict[str, Any], n_rows: int, *, relaxed: bool) -> bool:
    if stats["n"] == 0:
        return False
    if stats["strong"] >= 0.6:
        return True
    if relaxed:
        return stats["symbol"] >= 0.85 and n_rows >= 5
    return stats["symbol"] >= 0.85 and stats["unique"] >= 0.9 and n_rows >= 200


def _coerce_numeric_columns(df):
    import pandas as pd

    for c in list(df.columns):
        col = df[c]
        if col.dtype == object:
            non_null = int(col.notna().sum())
            if non_null == 0:
                continue
            conv = pd.to_numeric(col, errors="coerce")
            if int(conv.notna().sum()) == non_null:
                df[c] = conv
    return df


def _dedupe_columns(cols: list[Any]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for j, c in enumerate(cols):
        valid = c is not None and str(c).strip() and str(c) != "nan"
        name = str(c).strip() if valid else ("feature_id" if j == 0 else f"col_{j}")
        if j == 0 and name.lower().startswith("unnamed:"):
            name = "feature_id"
        if name in seen:
            seen[name] += 1
            name = f"{name}.{seen[name]}"
        else:
            seen[name] = 0
        out.append(name)
    return out


def _pick_header_row(raw, max_rows: int = 15) -> int:
    best_idx, best_score = 0, -1
    for i in range(min(max_rows, len(raw))):
        row = raw.iloc[i].tolist()
        score = 0
        for v in row:
            if v is None or (isinstance(v, float) and math.isnan(v)):
                continue
            s = str(v).strip()
            if not s:
                continue
            try:
                float(s)
            except ValueError:
                score += 1
        if score > best_score:
            best_idx, best_score = i, score
    return best_idx


def _load_delimited(path: Path):
    import pandas as pd

    opener_gz = path.name.lower().endswith(".gz")
    import gzip

    with (gzip.open(path, "rt", errors="replace") if opener_gz else open(path, "r", errors="replace")) as fh:
        head_lines = []
        for _ in range(50):
            line = fh.readline()
            if not line:
                break
            head_lines.append(line)
    skip = 0
    for line in head_lines:
        if line.startswith("#"):
            skip += 1
        else:
            break
    header = head_lines[skip] if len(head_lines) > skip else ""
    if "\t" in header:
        sep: str = "\t"
    elif "," in header:
        sep = ","
    elif ";" in header:
        sep = ";"
    else:
        sep = r"\s+"
    try:
        df = pd.read_csv(path, sep=sep, skiprows=skip, low_memory=False)
    except UnicodeDecodeError:
        df = pd.read_csv(path, sep=sep, skiprows=skip, low_memory=False, encoding="latin-1")
    df.columns = _dedupe_columns(list(df.columns))
    return [(None, df)]


def _load_excel(path: Path, max_sheets: int = 40):
    import pandas as pd

    suffix = path.suffix.lower()
    engine = "openpyxl" if suffix in {".xlsx", ".xlsm"} else None
    sheets = pd.read_excel(path, sheet_name=None, header=None, engine=engine)
    out = []
    for sheet_name, raw in list(sheets.items())[:max_sheets]:
        raw = raw.dropna(axis=0, how="all").dropna(axis=1, how="all")
        if raw.empty:
            out.append((str(sheet_name), raw))
            continue
        raw = raw.reset_index(drop=True)
        h = _pick_header_row(raw)
        df = raw.iloc[h + 1:].reset_index(drop=True)
        df.columns = _dedupe_columns(raw.iloc[h].tolist())
        df = df.dropna(axis=0, how="all")
        out.append((str(sheet_name), df))
    return out


def _profile_h5ad(path: Path) -> dict[str, Any]:
    import anndata as ad
    import numpy as np

    from aria.utils.count_classifier import classify_matrix

    adata = ad.read_h5ad(path, backed="r")
    try:
        n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
        x_slice = adata.X[: min(200, n_obs)]
        if hasattr(x_slice, "toarray"):
            x_slice = x_slice.toarray()
        x_slice = np.asarray(x_slice, dtype=float)
        x_cls = classify_matrix(x_slice.T, gene_ids=[str(v) for v in list(adata.var_names[:200])],
                                source_hint=path.name)
        raw_cls = None
        if adata.raw is not None:
            try:
                r = adata.raw.X[: min(200, n_obs)]
                r = r.toarray() if hasattr(r, "toarray") else r
                raw_cls = classify_matrix(np.asarray(r, dtype=float).T, source_hint=path.name)
            except Exception as exc:  # backed raw access can be unavailable
                raw_cls = {"error": f"{type(exc).__name__}: {str(exc)[:120]}"}
        obs = adata.obs.copy()
        var = adata.var.copy()
        return {
            "n_obs": n_obs,
            "n_vars": n_vars,
            "x_classification": x_cls,
            "raw_classification": raw_cls,
            "layers": sorted(list(adata.layers.keys())) if hasattr(adata, "layers") else [],
            "obs": obs,
            "var": var,
        }
    finally:
        try:
            adata.file.close()
        except Exception:
            pass


def _classify_frame(df, *, source_name: str, workbook_is_screen: bool) -> dict[str, Any]:
    """Classify one table into an ARIA native contract class (no question text used)."""
    import numpy as np

    from aria.utils.count_classifier import classify_matrix

    n_rows, n_cols = int(df.shape[0]), int(df.shape[1])
    cols = [str(c) for c in df.columns]
    info: dict[str, Any] = {"n_rows": n_rows, "n_cols": n_cols}
    if n_rows == 0 or n_cols == 0:
        info["contract_class"] = "empty"
        return info
    header_l = " | ".join(c.lower() for c in cols)
    seq_col = None
    for c in cols[:4]:
        vals = _sample_values(df[c], 300)
        if vals and sum(1 for v in vals if _SGRNA_RE.match(v)) / len(vals) >= 0.8:
            seq_col = c
            break
    if workbook_is_screen or seq_col is not None or any(t in header_l for t in _SCREEN_TOKENS):
        info["contract_class"] = "crispr_screen_table"
        info["reason"] = "sgRNA library / screen signature (no native ARIA contract)"
        return info
    variant_hits = [t for t in _VARIANT_TOKENS if t in header_l]
    norm_cols = {_norm_colname(c) for c in cols}
    if {"ref", "alt"} <= norm_cols and (norm_cols & {"pos", "position", "start", "chrom", "chr"}):
        variant_hits.append("ref+alt+pos")
    if len(variant_hits) >= 2:
        info["contract_class"] = "variant_table"
        info["reason"] = f"variant-call signature {variant_hits[:4]}"
        return info

    id_col, id_stats = None, None
    for c in cols[:3]:
        st = _gene_id_stats(df[c])
        if _qualifies_gene_ids(st, n_rows, relaxed=True) and (id_stats is None or st["strong"] > id_stats["strong"]):
            id_col, id_stats = c, st
    effect_cols = [c for c in cols if _is_effect_col(c)]
    sig_cols = [c for c in cols if _is_sig_col(c)]
    info["de_columns"] = {"effect": effect_cols[:8], "significance": sig_cols[:8]}
    if id_col is not None:
        info["id_column"] = id_col
        info["id_stats"] = id_stats
    if effect_cols and sig_cols and id_col is not None:
        info["contract_class"] = "de_result_table"
        info["reason"] = "gene identifiers + effect-size + significance columns"
        return info
    numeric_cols = [
        c for c in cols
        if c != id_col
        and str(df[c].dtype).startswith(("int", "float", "uint"))
        and not (_col_tokens(c) & _ANNOTATION_NUMERIC_TOKENS)
    ]
    info["n_numeric_columns"] = len(numeric_cols)
    info["numeric_columns_preview"] = numeric_cols[:10]
    if (
        id_col is not None
        and _qualifies_gene_ids(id_stats, n_rows, relaxed=False)
        and len(numeric_cols) >= 2
        and n_rows >= 100
    ):
        block = df[numeric_cols].to_numpy(dtype=float)
        cls = classify_matrix(block, gene_ids=[str(v) for v in df[id_col].head(200).tolist()],
                              source_hint=source_name)
        info["count_classification"] = {
            k: cls.get(k) for k in ("kind", "is_raw_counts", "raw_count_score", "confidence", "max", "min")
        }
        if not np.isfinite(block).any():
            info["contract_class"] = "other_tabular"
            return info
        info["contract_class"] = "expression_matrix"
        info["sample_columns"] = numeric_cols
        info["reason"] = f"gene x sample numeric matrix (classify_matrix kind={cls.get('kind')})"
        return info
    if id_col is not None and _qualifies_gene_ids(id_stats, n_rows, relaxed=False):
        info["contract_class"] = "gene_annotation"
        info["reason"] = "gene identifiers without an expression block"
        return info
    info["contract_class"] = "other_tabular"
    info["reason"] = "no gene-level expression/count/DE structure"
    return info


def _canon_sample(s: Any) -> str:
    return re.sub(r"[\s\-_.]+", "_", str(s).strip().lower())


class _Ledger:
    """Typed TableContract registry persisted under ``<ledger_dir>`` (worker side)."""

    def __init__(self, ledger_dir: Path):
        self.dir = Path(ledger_dir)
        self.tables_dir = self.dir / "tables"
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.registry_path = self.dir / "registry.json"
        self.registry: dict[str, Any] = (
            json.loads(self.registry_path.read_text()) if self.registry_path.exists() else {}
        )

    def save(self) -> None:
        _atomic_write_json(self.registry_path, self.registry)

    def unique_id(self, base: str) -> str:
        tid = _slug(base)
        k = 2
        out = tid
        while out in self.registry:
            out = f"{tid}_{k}"
            k += 1
        return out

    def register(self, table_id: str, df, *, biological_unit: str, source_sha256: str | None,
                 profile: dict[str, Any], key_candidates: list[str] | None = None) -> dict[str, Any]:
        pk = None
        for c in (key_candidates or []):
            if c in df.columns and df[c].notna().all() and not df[c].duplicated().any():
                pk = (c,)
                break
        if pk is None:
            df = df.reset_index(drop=True)
            df.insert(0, "__row_id", range(len(df)))
            pk = ("__row_id",)
        contract = TableContract(
            table_id=table_id,
            biological_unit=biological_unit,
            primary_key=pk,
            ledger_node_id=f"{LEDGER_NAMESPACE}/{table_id}",
            source_path=None,
            source_sha256=source_sha256,
        )
        rel_pickle = f"tables/{table_id}.pkl"
        df.to_pickle(self.dir / rel_pickle)
        self.registry[table_id] = {
            "table_id": table_id,
            "biological_unit": contract.biological_unit,
            "primary_key": list(contract.primary_key),
            "ledger_node_id": contract.ledger_node_id,
            "source_sha256": source_sha256,
            "pickle": rel_pickle,
            "columns": [str(c) for c in df.columns],
            "profile": profile,
        }
        return self.registry[table_id]

    def load(self, table_id: str):
        import pandas as pd

        if table_id not in self.registry:
            raise DerivedQueryViolation(
                f"Unknown table_id {table_id!r}; registered: {sorted(self.registry)[:40]}"
            )
        ent = self.registry[table_id]
        df = pd.read_pickle(self.dir / ent["pickle"])
        contract = TableContract(
            table_id=ent["table_id"],
            biological_unit=ent["biological_unit"],
            primary_key=tuple(ent["primary_key"]),
            ledger_node_id=ent["ledger_node_id"],
            source_path=None,
            source_sha256=ent.get("source_sha256"),
        )
        return contract, df


def _unit_for_class(cls: str) -> str:
    return {
        "expression_matrix": "gene",
        "de_result_table": "gene",
        "gene_annotation": "gene",
        "sample_metadata": "sample",
    }.get(cls, "sample")


def _worker_profile(req: dict[str, Any]) -> dict[str, Any]:
    from aria.utils.design_matrix import validate_design_matrix

    data_dir = Path(req["data_dir"])
    ledger = _Ledger(Path(req["ledger_dir"]))
    max_files = int(req.get("max_files", 150))
    rels = [str(r) for r in req.get("files") or []]
    file_reports: list[dict[str, Any]] = []
    frames: list[dict[str, Any]] = []
    for rel in rels[:max_files]:
        p = data_dir / rel
        low = rel.lower()
        rep: dict[str, Any] = {"file": rel, "tables": []}
        try:
            sha = _sha256_file(p)
            if low.endswith((".h5ad", ".h5mu")):
                prof = _profile_h5ad(p)
                x_cls = prof["x_classification"]
                ok_x = prof["n_vars"] >= 200 and x_cls.get("kind") != "empty"
                cls = "single_cell_expression" if ok_x else "other_tabular"
                stem = Path(rel).name.rsplit(".", 1)[0]
                summary = {
                    "contract_class": cls,
                    "n_obs": prof["n_obs"],
                    "n_vars": prof["n_vars"],
                    "x_classification": {k: x_cls.get(k) for k in ("kind", "is_raw_counts", "raw_count_score", "max")},
                    "raw_classification": (
                        {k: prof["raw_classification"].get(k) for k in ("kind", "is_raw_counts", "error")}
                        if prof["raw_classification"] else None
                    ),
                    "layers": prof["layers"],
                    "reason": "AnnData cell x gene matrix",
                }
                obs = prof["obs"].copy()
                obs.insert(0, "cell_id", obs.index.astype(str))
                obs = obs.reset_index(drop=True)
                for c in obs.columns:
                    if str(obs[c].dtype) == "category":
                        obs[c] = obs[c].astype(str)
                obs_prof = dict(summary)
                obs_prof["table_role"] = "obs (cell metadata)"
                obs_prof["factor_levels"] = _factor_levels(obs, max_levels=40)
                tid = ledger.unique_id(f"{stem}__obs")
                ledger.register(tid, obs, biological_unit="cell", source_sha256=sha,
                                profile={"file": rel, "sheet": None, **obs_prof}, key_candidates=["cell_id"])
                rep["tables"].append(tid)
                var = prof["var"].copy()
                var.insert(0, "gene_id", var.index.astype(str))
                var = var.reset_index(drop=True)
                for c in var.columns:
                    if str(var[c].dtype) == "category":
                        var[c] = var[c].astype(str)
                tid2 = ledger.unique_id(f"{stem}__var")
                ledger.register(tid2, var, biological_unit="gene", source_sha256=sha,
                                profile={"file": rel, "sheet": None, "table_role": "var (gene metadata)",
                                         "contract_class": "gene_annotation", "n_rows": int(len(var)),
                                         "n_cols": int(var.shape[1])},
                                key_candidates=["gene_id"])
                rep["tables"].append(tid2)
                rep["contract_class"] = cls
                rep.update({k: summary[k] for k in ("n_obs", "n_vars", "x_classification")})
                file_reports.append(rep)
                continue
            if low.endswith((".xlsx", ".xlsm", ".xls")):
                loaded = _load_excel(p)
            else:
                loaded = _load_delimited(p)
            sheet_names = " ".join(str(s or "").lower() for s, _ in loaded)
            workbook_is_screen = any(t in sheet_names for t in _SCREEN_TOKENS)
            for sheet, df in loaded:
                df = _coerce_numeric_columns(df)
                info = _classify_frame(df, source_name=Path(rel).name, workbook_is_screen=workbook_is_screen)
                frames.append({"file": rel, "sheet": sheet, "sha": sha, "df": df, "info": info})
        except Exception as exc:
            rep["contract_class"] = "unreadable"
            rep["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
        file_reports.append(rep)

    # Cross-link sample metadata to expression-matrix sample columns.
    sample_names: dict[str, set[str]] = {}
    for fr in frames:
        if fr["info"].get("contract_class") == "expression_matrix":
            key = f"{fr['file']}::{fr['sheet']}"
            sample_names[key] = {str(s) for s in fr["info"].get("sample_columns", [])}
    for fr in frames:
        info = fr["info"]
        if info.get("contract_class") not in {"other_tabular", "gene_annotation"} or not sample_names:
            continue
        df = fr["df"]
        if len(df) > 5000:
            continue
        best = None
        for c in list(df.columns)[:25]:
            vals = {str(v).strip() for v in df[c].dropna().tolist()}
            if len(vals) < 2:
                continue
            for mkey, names in sample_names.items():
                exact = len(vals & names)
                canon_names = {_canon_sample(n) for n in names}
                norm = len({_canon_sample(v) for v in vals} & canon_names)
                score = max(exact, norm)
                frac = score / max(1, min(len(vals), len(names)))
                if score >= 2 and frac >= 0.5 and (best is None or score > best[0]):
                    best = (score, c, mkey, "exact" if exact >= norm else "normalized", round(frac, 3))
        if best is not None:
            info["contract_class"] = "sample_metadata"
            info["sample_id_column"] = best[1]
            info["matched_matrix"] = best[2]
            info["sample_id_match_mode"] = best[3]
            info["sample_id_match_fraction"] = best[4]
            info["reason"] = "sample-level metadata whose IDs match expression-matrix columns"

    for fr in frames:
        info = fr["info"]
        df = fr["df"]
        cls = info.get("contract_class", "other_tabular")
        if cls in {"sample_metadata", "other_tabular"} and len(df) <= 5000:
            info["factor_levels"] = _factor_levels(df, exclude=[info.get("sample_id_column")])
        if cls == "sample_metadata":
            checks: dict[str, Any] = {}
            for col in list(info.get("factor_levels", {}).keys())[:6]:
                try:
                    res = validate_design_matrix(df, col, min_replicates_per_condition=2)
                    checks[col] = {
                        "status": res.get("status"),
                        "condition_counts": res.get("condition_counts"),
                        "issues": [i.get("check") for i in res.get("issues", [])][:4],
                    }
                except Exception as exc:
                    checks[col] = {"status": "error", "error": str(exc)[:120]}
            info["design_checks"] = checks
        base = Path(fr["file"]).name
        for suf in (".gz",):
            if base.lower().endswith(suf):
                base = base[: -len(suf)]
        base = base.rsplit(".", 1)[0]
        tid = ledger.unique_id(f"{base}__{fr['sheet']}" if fr["sheet"] else base)
        prof = {"file": fr["file"], "sheet": fr["sheet"], **{k: v for k, v in info.items() if k != "sample_columns"}}
        if "sample_columns" in info:
            prof["n_sample_columns"] = len(info["sample_columns"])
            prof["sample_columns_preview"] = info["sample_columns"][:12]
        key_candidates = [c for c in [info.get("id_column"), info.get("sample_id_column")] if c]
        ledger.register(tid, df, biological_unit=_unit_for_class(cls), source_sha256=fr["sha"],
                        profile=prof, key_candidates=key_candidates)
        for rep in file_reports:
            if rep["file"] == fr["file"]:
                rep["tables"].append(tid)
                prev = rep.get("contract_class")
                if prev is None or _TABLE_CLASS_RANK.get(cls, 99) < _TABLE_CLASS_RANK.get(prev, 99):
                    rep["contract_class"] = cls
    for rep in file_reports:
        rep.setdefault("contract_class", "empty")
    ledger.save()
    return {"status": "ok", "files": file_reports, "n_tables": len(ledger.registry)}


def _factor_levels(df, *, max_levels: int = 12, exclude: list[Any] | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for c in df.columns:
        if c in (exclude or []) or c == "__row_id":
            continue
        s = df[c].dropna()
        if s.empty:
            continue
        nun = int(s.astype(str).nunique())
        if 2 <= nun <= max_levels:
            vc = s.astype(str).value_counts().sort_index()
            out[str(c)] = {str(k): int(v) for k, v in vc.items()}
        if len(out) >= 15:
            break
    return out


def _parse_step(d: dict[str, Any]):
    kind = str(d.get("step") or d.get("type") or "").lower()
    if not kind and str(d.get("op", "")).lower() in {"filter", "project", "rank", "group_summary", "validated_join"}:
        kind = str(d["op"]).lower()
    if kind == "filter":
        op = d.get("operator") if "operator" in d else d.get("op")
        if str(op).lower() == "filter":
            op = d.get("operator")
        return FilterOp(column=str(d["column"]), op=str(op), value=d.get("value"))
    if kind == "project":
        return ProjectOp(columns=tuple(str(c) for c in d["columns"]))
    if kind == "rank":
        by = d["by"] if isinstance(d["by"], list) else [d["by"]]
        asc = d.get("ascending", [True] * len(by))
        asc = asc if isinstance(asc, list) else [bool(asc)] * len(by)
        lim = d.get("limit")
        return RankOp(by=tuple(str(b) for b in by), ascending=tuple(bool(a) for a in asc),
                      limit=None if lim is None else int(lim))
    if kind == "group_summary":
        by = d["by"] if isinstance(d["by"], list) else [d["by"]]
        return GroupSummaryOp(by=tuple(str(b) for b in by),
                              aggregations={str(k): str(v) for k, v in dict(d["aggregations"]).items()})
    if kind == "validated_join":
        on = d["on"] if isinstance(d["on"], list) else [d["on"]]
        return ValidatedJoinOp(
            right_table_id=str(d["right_table_id"]),
            on=tuple(str(c) for c in on),
            how=str(d.get("how", "inner")),
            expected_cardinality=str(d.get("expected_cardinality", "1:1")),
            allow_unit_projection=bool(d.get("allow_unit_projection", False)),
        )
    raise DerivedQueryViolation(f"Unsupported DerivedQuery step: {d!r}")


def _coerce_filter_values(steps, df_columns_dtypes: dict[str, str]):
    out = []
    for s in steps:
        if isinstance(s, FilterOp) and isinstance(s.value, str):
            dt = df_columns_dtypes.get(s.column, "")
            if dt.startswith(("int", "float", "uint")):
                try:
                    s = FilterOp(column=s.column, op=s.op, value=float(s.value))
                except ValueError:
                    pass
        out.append(s)
    return out


def _run_plan(ledger: _Ledger, plan_d: dict[str, Any]) -> dict[str, Any]:
    base_id = str(plan_d["base_table_id"])
    steps = [_parse_step(s) for s in plan_d.get("steps") or []]
    needed = {base_id} | {s.right_table_id for s in steps if isinstance(s, ValidatedJoinOp)}
    tables = {tid: ledger.load(tid) for tid in needed}
    dtypes = {str(c): str(t) for c, t in tables[base_id][1].dtypes.items()}
    steps = _coerce_filter_values(steps, dtypes)
    plan = DerivedQueryPlan(base_table_id=base_id, steps=tuple(steps))
    return execute_derived_query(plan, tables)


def _register_result_frame(ledger: _Ledger, name: str, df, unit: str, parent_ids: list[str],
                           key_candidates: list[str]) -> str:
    tid = ledger.unique_id(f"derived_{name}")
    ledger.register(tid, df, biological_unit=unit, source_sha256=None,
                    profile={"file": None, "sheet": None, "derived_from": parent_ids,
                             "contract_class": "derived_table", "n_rows": int(len(df)),
                             "n_cols": int(df.shape[1])},
                    key_candidates=key_candidates)
    ledger.save()
    return tid


def _worker_derived_query(req: dict[str, Any]) -> dict[str, Any]:
    import pandas as pd

    ledger = _Ledger(Path(req["ledger_dir"]))
    res = _run_plan(ledger, req["plan"])
    max_rows = max(1, min(int(req.get("max_rows") or 25), 50))
    out = {k: res[k] for k in ("schema_version", "ast_sha256", "ledger_node_ids", "biological_unit",
                                "input_denominator", "output_row_count", "columns")}
    out["rows_shown"] = res["rows"][:max_rows]
    out["rows_truncated"] = len(res["rows"]) > max_rows
    if req.get("register_as"):
        df = pd.DataFrame(res["rows"], columns=res["columns"])
        out["registered_table_id"] = _register_result_frame(
            ledger, str(req["register_as"]), df, res["biological_unit"],
            list(res["ledger_node_ids"]), [c for c in res["columns"][:3]],
        )
    return {"status": "ok", "result": out}


def _numeric(series):
    import pandas as pd

    return pd.to_numeric(series, errors="coerce")


def _worker_inspect_stats(req: dict[str, Any]) -> dict[str, Any]:
    import numpy as np
    import pandas as pd
    from scipy import stats as sps

    from aria.utils.stats import bh_correct

    ledger = _Ledger(Path(req["ledger_dir"]))
    a = dict(req.get("args") or {})
    op = str(a.get("op", "")).lower()
    out: dict[str, Any] = {"op": op}

    if op == "fisher_exact":
        tbl = np.asarray(a["table"], dtype=float)
        if tbl.shape != (2, 2) or (tbl < 0).any():
            raise ValueError("fisher_exact requires a non-negative 2x2 table")
        odds, p = sps.fisher_exact(tbl.astype(int), alternative=str(a.get("alternative", "two-sided")))
        out.update({"table": tbl.astype(int).tolist(), "odds_ratio": float(odds), "p_value": float(p)})
        return {"status": "ok", "result": out}
    if op == "hypergeometric":
        k, M, n, N = (int(a[x]) for x in ("k", "M", "n", "N"))
        p = float(sps.hypergeom.sf(k - 1, M, n, N))
        out.update({"k": k, "M": M, "n": n, "N": N, "p_value_upper_tail": p,
                    "expected_overlap": float(n * N / M) if M else None,
                    "fold_enrichment": float(k / (n * N / M)) if M and n and N else None})
        return {"status": "ok", "result": out}
    if op == "count_where":
        ids = a.get("table_ids")
        if not ids and a.get("table_prefix"):
            ids = [t for t in sorted(ledger.registry) if t.startswith(str(a["table_prefix"]))]
        ids = [str(t) for t in (ids or [])][:80]
        filters = [{"step": "filter", **f} if "step" not in f else f for f in (a.get("filters") or [])]
        counts: dict[str, Any] = {}
        for tid in ids:
            try:
                r = _run_plan(ledger, {"base_table_id": tid, "steps": filters})
                counts[tid] = {"n_rows_matching": r["output_row_count"], "n_rows_total": r["input_denominator"]}
            except Exception as exc:
                counts[tid] = {"error": f"{type(exc).__name__}: {str(exc)[:160]}"}
        out["counts"] = counts
        return {"status": "ok", "result": out}

    tid = str(a["table_id"])
    if a.get("pre_steps"):
        r = _run_plan(ledger, {"base_table_id": tid, "steps": a["pre_steps"]})
        df = pd.DataFrame(r["rows"], columns=r["columns"])
        out["pre_steps_ast_sha256"] = r["ast_sha256"]
    else:
        _, df = ledger.load(tid)
    out["table_id"] = tid
    out["n_rows_input"] = int(len(df))

    def _need(col: str) -> str:
        if col not in df.columns:
            raise DerivedQueryViolation(f"Column {col!r} not in table {tid!r} (columns: {list(df.columns)[:40]})")
        return col

    if op == "head":
        n = max(1, min(int(a.get("n", 10)), 25))
        cols = [_need(c) for c in (a.get("columns") or list(df.columns)[:30])]
        out["rows"] = df[cols].head(n).to_dict(orient="records")
    elif op == "describe":
        cols = a.get("columns") or [c for c in df.columns if str(df[c].dtype).startswith(("int", "float"))][:20]
        desc: dict[str, Any] = {}
        for c in cols[:40]:
            _need(c)
            s = _numeric(df[c]).dropna()
            if len(s) and len(s) >= 0.5 * df[c].notna().sum():
                desc[c] = {
                    "n": int(len(s)), "mean": float(s.mean()), "std": float(s.std(ddof=1)) if len(s) > 1 else None,
                    "min": float(s.min()), "q25": float(s.quantile(0.25)), "median": float(s.median()),
                    "q75": float(s.quantile(0.75)), "max": float(s.max()),
                    "skewness": float(sps.skew(s)) if len(s) > 2 else None,
                    "excess_kurtosis": float(sps.kurtosis(s)) if len(s) > 3 else None,
                    "n_positive": int((s > 0).sum()), "n_negative": int((s < 0).sum()), "n_zero": int((s == 0).sum()),
                }
            else:
                vc = df[c].astype(str).value_counts().head(10)
                desc[c] = {"n": int(df[c].notna().sum()), "n_unique": int(df[c].nunique()),
                           "top_values": {str(k): int(v) for k, v in vc.items()}}
        out["describe"] = desc
    elif op == "value_counts":
        c = _need(str(a["column"]))
        top = max(1, min(int(a.get("top", 30)), 100))
        vc = df[c].astype(str).value_counts()
        out.update({"column": c, "n_unique": int(vc.size), "n_total": int(df[c].notna().sum()),
                    "counts": {str(k): int(v) for k, v in vc.head(top).items()}})
    elif op == "corr":
        x, y = _need(str(a["x"])), _need(str(a["y"]))
        method = str(a.get("method", "pearson")).lower()
        pair = pd.DataFrame({"x": _numeric(df[x]), "y": _numeric(df[y])}).dropna()
        if len(pair) < 3:
            raise ValueError("corr needs >= 3 complete pairs")
        res = sps.spearmanr(pair["x"], pair["y"]) if method == "spearman" else sps.pearsonr(pair["x"], pair["y"])
        out.update({"x": x, "y": y, "method": method, "n": int(len(pair)),
                    "r": float(res[0]), "p_value": float(res[1])})
    elif op in {"ttest", "wilcoxon"}:
        v, g = _need(str(a["value_column"])), _need(str(a["group_column"]))
        ga, gb = str(a["group_a"]), str(a["group_b"])
        xa = _numeric(df.loc[df[g].astype(str) == ga, v]).dropna()
        xb = _numeric(df.loc[df[g].astype(str) == gb, v]).dropna()
        if len(xa) < 2 or len(xb) < 2:
            raise ValueError(f"need >=2 values per group (got {len(xa)}, {len(xb)})")
        if op == "ttest":
            equal_var = bool(a.get("equal_var", False))
            r = sps.ttest_ind(xa, xb, equal_var=equal_var)
            out.update({"test": "student_t" if equal_var else "welch_t", "statistic": float(r.statistic),
                        "p_value": float(r.pvalue)})
        else:
            r = sps.mannwhitneyu(xa, xb, alternative=str(a.get("alternative", "two-sided")))
            out.update({"test": "wilcoxon_rank_sum_mann_whitney_u", "statistic_U": float(r.statistic),
                        "p_value": float(r.pvalue)})
        out.update({"group_a": ga, "group_b": gb, "n_a": int(len(xa)), "n_b": int(len(xb)),
                    "mean_a": float(xa.mean()), "mean_b": float(xb.mean()),
                    "median_a": float(xa.median()), "median_b": float(xb.median())})
    elif op == "anova":
        v, g = _need(str(a["value_column"])), _need(str(a["group_column"]))
        levels = [str(x) for x in (a.get("groups") or sorted(df[g].dropna().astype(str).unique()))][:50]
        samples = [_numeric(df.loc[df[g].astype(str) == lv, v]).dropna() for lv in levels]
        kept = [(lv, s) for lv, s in zip(levels, samples) if len(s) >= 2]
        if len(kept) < 2:
            raise ValueError("anova needs >= 2 groups with >= 2 values")
        r = sps.f_oneway(*[s for _, s in kept])
        out.update({"test": "one_way_anova", "F": float(r.statistic), "p_value": float(r.pvalue),
                    "groups": {lv: int(len(s)) for lv, s in kept}})
    elif op == "bh_adjust":
        c = _need(str(a["column"]))
        alpha = float(a.get("alpha", 0.05))
        p = _numeric(df[c])
        mask = p.notna()
        adj = np.full(len(df), np.nan)
        adj[mask.to_numpy()] = bh_correct(p[mask].to_numpy())
        new_col = f"{c}_bh"
        df = df.copy()
        df[new_col] = adj
        out.update({"column": c, "new_column": new_col, "alpha": alpha, "n_tested": int(mask.sum()),
                    "n_below_alpha": int((adj < alpha).sum()), "n_raw_below_alpha": int((p < alpha).sum())})
        if a.get("register_as"):
            out["registered_table_id"] = _register_result_frame(
                ledger, str(a["register_as"]), df, ledger.registry[tid]["biological_unit"], [tid],
                list(ledger.registry[tid]["primary_key"]))
    elif op == "row_aggregate":
        cols = [_need(str(c)) for c in (a.get("columns") or [])]
        if not cols:
            raise ValueError("row_aggregate requires an explicit non-empty 'columns' list")
        func = str(a.get("func", "mean")).lower()
        if func not in {"mean", "median", "sum", "var", "std", "min", "max"}:
            raise ValueError(f"unsupported row_aggregate func {func!r}")
        key = _need(str(a.get("key_column") or ledger.registry[tid]["primary_key"][0]))
        vals = df[cols].apply(pd.to_numeric, errors="coerce")
        agg = getattr(vals, func)(axis=1)
        name = str(a.get("value_name") or f"row_{func}")
        res_df = pd.DataFrame({key: df[key].to_numpy(), name: agg.to_numpy()})
        s = res_df[name].dropna()
        out.update({"func": func, "n_columns_aggregated": len(cols), "value_name": name,
                    "summary": {"n": int(len(s)), "mean": float(s.mean()) if len(s) else None,
                                "median": float(s.median()) if len(s) else None}})
        if a.get("register_as"):
            out["registered_table_id"] = _register_result_frame(
                ledger, str(a["register_as"]), res_df, ledger.registry[tid]["biological_unit"], [tid], [key])
        else:
            out["rows_preview"] = res_df.head(10).to_dict(orient="records")
    else:
        raise ValueError(
            f"Unsupported inspect_table_stats op {op!r}; allowed: head, describe, value_counts, corr, ttest, "
            "wilcoxon, anova, fisher_exact, hypergeometric, bh_adjust, row_aggregate, count_where"
        )
    return {"status": "ok", "result": out}


def _worker_prepare_de(req: dict[str, Any]) -> dict[str, Any]:
    """Assemble a validated counts CSV + metadata TSV for ARIA's bulk DE script."""
    import pandas as pd

    from aria.utils.count_classifier import validate_raw_count_matrix
    from aria.utils.design_matrix import validate_design_matrix

    ledger = _Ledger(Path(req["ledger_dir"]))
    a = dict(req["args"])
    out_dir = Path(req["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    _, counts_df = ledger.load(str(a["counts_table_id"]))
    _, meta_df = ledger.load(str(a["metadata_table_id"]))
    sample_col = str(a["sample_column"])
    factor = str(a["design_factor"])
    covariates = [str(c) for c in (a.get("covariates") or [])]
    for c in [sample_col, factor, *covariates]:
        if c not in meta_df.columns:
            return {"status": "error", "error_type": "MissingMetadataColumn",
                    "details": f"{c!r} not in metadata columns {list(meta_df.columns)[:40]}"}
    meta = meta_df.copy()
    for col, allowed in dict(a.get("subset") or {}).items():
        if col not in meta.columns:
            return {"status": "error", "error_type": "MissingSubsetColumn", "details": f"{col!r} not in metadata"}
        allowed_l = [str(v) for v in (allowed if isinstance(allowed, list) else [allowed])]
        meta = meta[meta[col].astype(str).isin(allowed_l)]
    meta = meta.dropna(subset=[sample_col])
    meta[sample_col] = meta[sample_col].astype(str).str.strip()
    if meta[sample_col].duplicated().any():
        return {"status": "error", "error_type": "DuplicateSampleIds",
                "details": f"duplicated sample IDs in {sample_col!r}"}
    centry = ledger.registry[str(a["counts_table_id"])]
    id_col = centry["profile"].get("id_column") or centry["primary_key"][0]
    sample_cols = [c for c in counts_df.columns if c not in {id_col, "__row_id"}]
    normalize = bool(a.get("normalize_sample_ids", False))
    mapping: dict[str, str] = {}
    meta_ids = list(meta[sample_col])
    for sid in meta_ids:
        if sid in sample_cols:
            mapping[sid] = sid
    if normalize:
        canon_cols: dict[str, list[str]] = {}
        for c in sample_cols:
            canon_cols.setdefault(_canon_sample(c), []).append(str(c))
        for sid in meta_ids:
            if sid in mapping:
                continue
            hits = canon_cols.get(_canon_sample(sid), [])
            if len(hits) == 1:
                mapping[sid] = hits[0]
            elif len(hits) > 1:
                return {"status": "error", "error_type": "AmbiguousSampleIdNormalization",
                        "details": f"{sid!r} maps to multiple count columns {hits[:4]}"}
    unmatched = [s for s in meta_ids if s not in mapping]
    if unmatched:
        return {"status": "error", "error_type": "MetadataCountColumnMismatch",
                "details": (f"{len(unmatched)} metadata sample(s) have no count column "
                            f"(examples {unmatched[:5]}; count columns e.g. {sample_cols[:5]}); "
                            "set normalize_sample_ids=true only if IDs differ by separators/case")}
    ids = counts_df[id_col].astype(str)
    if ids.duplicated().any():
        return {"status": "error", "error_type": "DuplicateFeatureIds",
                "details": f"{int(ids.duplicated().sum())} duplicated feature IDs in {id_col!r}"}
    counts = counts_df[[mapping[s] for s in meta_ids]].apply(pd.to_numeric, errors="coerce")
    counts.columns = meta_ids
    counts.index = ids.to_numpy()
    counts.index.name = "gene_id"
    raw_check = validate_raw_count_matrix(counts.to_numpy(dtype=float))
    meta_out = meta.set_index(sample_col)[[factor, *covariates]].astype(str)
    meta_out.index.name = "sample"
    prefix = str(a.get("categorical_label_prefix") or "")
    relabel = None
    if prefix:
        # Audited application of ARIA's own design recommendation: encode a
        # numeric-coded condition as categorical labels (never implicit).
        meta_out[factor] = [v if v.startswith(prefix) else f"{prefix}{v}" for v in meta_out[factor]]
        relabel = {"design_factor": factor, "categorical_label_prefix": prefix,
                   "levels": sorted(set(meta_out[factor]))}
    design = validate_design_matrix(meta_out, factor, covariates,
                                    min_replicates_per_condition=int(a.get("min_replicates_per_condition", 3)))
    counts_path = out_dir / "counts.csv"
    meta_path = out_dir / "metadata.tsv"
    counts.to_csv(counts_path)
    meta_out.to_csv(meta_path, sep="\t")
    return {
        "status": "ok",
        "counts_path": str(counts_path),
        "metadata_path": str(meta_path),
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "sample_id_mapping_mode": "normalized" if normalize and any(k != v for k, v in mapping.items()) else "exact",
        "condition_relabel": relabel,
        "raw_count_check": {k: raw_check.get(k) for k in ("valid", "reason", "n_offending")},
        "design_check": {k: design.get(k) for k in ("status", "issues", "condition_counts", "rank", "n_columns",
                                                     "residual_degrees_of_freedom")},
    }


def _worker_register_table(req: dict[str, Any]) -> dict[str, Any]:
    import pandas as pd

    ledger = _Ledger(Path(req["ledger_dir"]))
    df = pd.read_csv(req["path"], sep="\t")
    tid = ledger.unique_id(str(req["name"]))
    ledger.register(tid, df, biological_unit=str(req.get("biological_unit", "gene")),
                    source_sha256=_sha256_file(Path(req["path"])),
                    profile={"file": None, "sheet": None, "contract_class": "aria_de_result",
                             "derived_from": req.get("derived_from") or [], "n_rows": int(len(df)),
                             "n_cols": int(df.shape[1])},
                    key_candidates=["gene_id"])
    ledger.save()
    return {"status": "ok", "table_id": tid, "n_rows": int(len(df)), "columns": [str(c) for c in df.columns]}


_WORKER_HANDLERS: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    "profile": _worker_profile,
    "derived_query": _worker_derived_query,
    "inspect_stats": _worker_inspect_stats,
    "prepare_de": _worker_prepare_de,
    "register_table": _worker_register_table,
}


def _worker_main(req_path: str, out_path: str) -> int:
    import warnings

    warnings.filterwarnings("ignore")
    req = json.loads(Path(req_path).read_text())
    try:
        res = _WORKER_HANDLERS[req["kind"]](req)
    except Exception as exc:
        res = {"status": "error", "error_type": type(exc).__name__, "details": str(exc)[:600]}
    Path(out_path).write_text(json.dumps(_jsonable(res), ensure_ascii=False))
    return 0


# ──────────────────────────────────────────────────────────────────────────────
# HARNESS side: worker invocation, contract gate, governed loop
# ──────────────────────────────────────────────────────────────────────────────


def _resolve_rna_python(rna_python: str | None) -> str:
    cand = rna_python or DEFAULT_RNA_PYTHON
    return cand if Path(cand).exists() else sys.executable


def _subprocess_env(scratch: Path) -> dict[str, str]:
    base = {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "HOME": os.environ.get("HOME", str(scratch)),
        "TMPDIR": str(scratch),
        "PYTHONPATH": str(REPO_ROOT),
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(scratch / ".mpl"),
        "NUMBA_CACHE_DIR": str(scratch / ".numba"),
        "OMP_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "4",
        "MKL_NUM_THREADS": "4",
        "NO_PROXY": "127.0.0.1,localhost",
        "HTTP_PROXY": "http://127.0.0.1:9",
        "HTTPS_PROXY": "http://127.0.0.1:9",
    }
    return isolated_benchmark_environment(scratch, base_env=base)


class WorkerClient:
    """Invoke the aria-rna-env table worker (this same file in ``--worker`` mode)."""

    def __init__(self, scratch: Path, rna_python: str | None = None, timeout: int = WORKER_TIMEOUT):
        self.scratch = Path(scratch).resolve()
        self.python = _resolve_rna_python(rna_python)
        self.timeout = timeout
        self.env = _subprocess_env(self.scratch)
        self._n = 0

    def call(self, kind: str, payload: dict[str, Any], timeout: int | None = None) -> dict[str, Any]:
        self._n += 1
        req_path = self.scratch / f"_worker_req_{self._n:03d}.json"
        out_path = self.scratch / f"_worker_out_{self._n:03d}.json"
        req_path.write_text(json.dumps(_jsonable({"kind": kind, **payload}), ensure_ascii=False))
        try:
            proc = subprocess.run(
                [self.python, str(Path(__file__).resolve()), "--worker", str(req_path), str(out_path)],
                cwd=str(self.scratch), env=self.env, capture_output=True, text=True,
                timeout=timeout or self.timeout, check=False,
            )
        except subprocess.TimeoutExpired:
            return {"status": "error", "error_type": "WorkerTimeout", "details": f"{kind} exceeded timeout"}
        if not out_path.exists():
            return {"status": "error", "error_type": "WorkerCrashed",
                    "details": (proc.stderr or "")[-600:]}
        return json.loads(out_path.read_text())


def _classify_suffix(rel: str) -> tuple[str | None, str | None]:
    low = rel.lower()
    for suf in sorted(UNSUPPORTED_SUFFIX_CLASSES, key=len, reverse=True):
        if low.endswith(suf):
            return "unsupported", UNSUPPORTED_SUFFIX_CLASSES[suf]
    for suf, cls in AUXILIARY_SUFFIX_CLASSES.items():
        if low.endswith(suf):
            return "auxiliary", cls
    if low.endswith(TABULAR_SUFFIXES):
        return "tabular", None
    return "other", "unrecognized_format"


def run_contract_gate(data_dir: Path, worker: WorkerClient, ledger_dir: Path) -> dict[str, Any]:
    """ARIA Stage 1 — deterministic native modality & file contract gate (before any LLM call)."""
    data_dir = Path(data_dir)
    assert_blind_workspace(data_dir)
    rels = sorted(p.relative_to(data_dir).as_posix() for p in data_dir.rglob("*") if p.is_file())
    seq_hits = detect_unsupported_inputs(rels)
    detector = AssayDetector()
    file_rows: dict[str, dict[str, Any]] = {}
    tabular: list[str] = []
    for rel in rels:
        group, cls = _classify_suffix(rel)
        row: dict[str, Any] = {"file": rel, "suffix_group": group}
        try:
            det = detector.detect_file(data_dir / rel)
        except Exception as exc:  # detector must never crash the gate
            det = None
            row["assay_detector_error"] = type(exc).__name__
        if det is not None:
            row["assay_detector"] = {"modality": det.modality, "confidence": det.confidence,
                                     "reason": det.reason[:160]}
        if group == "tabular":
            tabular.append(rel)
        else:
            row["contract_class"] = cls
        file_rows[rel] = row
    profile: dict[str, Any] = {"status": "skipped", "files": []}
    if tabular:
        profile = worker.call("profile", {"data_dir": str(data_dir), "ledger_dir": str(ledger_dir),
                                          "files": tabular})
        for rep in profile.get("files", []) or []:
            if rep["file"] in file_rows:
                file_rows[rep["file"]]["contract_class"] = rep.get("contract_class", "unreadable")
                file_rows[rep["file"]]["tables"] = rep.get("tables", [])
                if rep.get("error"):
                    file_rows[rep["file"]]["error"] = rep["error"]
        if profile.get("status") != "ok":
            for rel in tabular:
                file_rows[rel].setdefault("contract_class", "unreadable")
    for rel in tabular:
        file_rows[rel].setdefault("contract_class", "unreadable")

    omics_ok_detector = {"bulk_RNA", "scRNA", "single_cell_RNA", "scATAC", "bulk_ATAC", "chromatin"}
    omics_files: list[str] = []
    non_omics: list[str] = []
    for rel, row in file_rows.items():
        cls = row.get("contract_class")
        det_mod = (row.get("assay_detector") or {}).get("modality")
        if cls in CONTRACT_TABLE_CLASSES or cls in SUPPORTING_TABLE_CLASSES:
            row["group"] = "omics_tabular_or_h5ad"
            omics_files.append(rel)
        elif cls in {"other_tabular", "unreadable", "empty"} and det_mod in omics_ok_detector:
            row["group"] = "omics_tabular_or_h5ad"
            row["contract_class"] = "assay_detector_verified_omics"
            omics_files.append(rel)
        elif cls in AUXILIARY_SUFFIX_CLASSES.values():
            row["group"] = "auxiliary_resource"
        else:
            row["group"] = "unsupported_native_or_non_omics"
            non_omics.append(rel)
    contract_files = [
        r for r in omics_files
        if file_rows[r]["contract_class"] in CONTRACT_TABLE_CLASSES
        or file_rows[r]["contract_class"] == "assay_detector_verified_omics"
    ]
    passed = (not seq_hits) and bool(contract_files)
    detected = [f"{Path(r).name}[{file_rows[r].get('contract_class')}]" for r in (non_omics or rels)]
    reason = None
    if not passed:
        listing = ",".join(detected[:10]) + (f",+{len(detected) - 10}_more" if len(detected) > 10 else "")
        if seq_hits:
            listing = ",".join(f"{h['path']}[{h['format']}]" for h in seq_hits[:10]) + ("," + listing if listing else "")
        reason = f"outside_native_omics_contract:{listing or 'no_files'}"
    return {
        "passed": passed,
        "abstain_reason": reason,
        "n_files": len(rels),
        "sequencer_native_hits": seq_hits,
        "contract_files": contract_files,
        "omics_tabular_or_h5ad": omics_files,
        "unsupported_native_or_non_omics": non_omics,
        "files": [file_rows[r] for r in rels][:200],
        "profile_status": profile.get("status"),
        "profile_error": profile.get("details"),
    }


def _load_registry(ledger_dir: Path) -> dict[str, Any]:
    p = Path(ledger_dir) / "registry.json"
    return json.loads(p.read_text()) if p.exists() else {}


def render_ledger_for_prompt(registry: dict[str, Any], max_chars: int = 14000) -> str:
    lines: list[str] = []
    order = sorted(registry.values(), key=lambda e: (
        _TABLE_CLASS_RANK.get(str(e.get("profile", {}).get("contract_class")), 50), e["table_id"]))
    for ent in order:
        prof = ent.get("profile", {})
        cols = ent.get("columns", [])
        d: dict[str, Any] = {
            "table_id": ent["table_id"],
            "ledger_node": ent["ledger_node_id"],
            "source": {"file": prof.get("file"), "sheet": prof.get("sheet")},
            "contract_class": prof.get("contract_class"),
            "biological_unit": ent["biological_unit"],
            "primary_key": ent["primary_key"],
            "shape": [prof.get("n_rows", prof.get("n_obs")), len(cols)],
            "columns": cols[:30] + ([f"...(+{len(cols) - 30} more)"] if len(cols) > 30 else []),
        }
        for k in ("id_column", "count_classification", "de_columns", "n_sample_columns", "sample_id_column",
                  "sample_id_match_mode", "factor_levels", "design_checks", "x_classification", "n_vars",
                  "table_role", "derived_from"):
            if prof.get(k) not in (None, {}, []):
                d[k] = prof[k]
        lines.append(json.dumps(_jsonable(d), ensure_ascii=False))
    text = "\n".join(lines)
    if len(text) > max_chars:
        text = text[:max_chars] + "\n...(ledger listing truncated; use inspect_table_stats head/describe)"
    return text


ARIA_GOVERNED_SYSTEM_PROMPT = """\
You are the planner inside ARIA v4.9.0, a governed bioinformatics system. You CANNOT write or run code.
You may only request ARIA's verified, deterministic primitives over the registered ARIA Evidence Ledger tables.
Each turn, emit exactly ONE JSON object inside <ACTION>...</ACTION> tags and nothing else of substance.

Allowed actions:
1. inspect_table_stats — bounded deterministic statistics on a registered table:
   {"action":"inspect_table_stats","table_id":"...","op":"<op>", ...op args..., "pre_steps":[optional DerivedQuery steps]}
   ops: head{n,columns} | describe{columns} | value_counts{column,top} | corr{x,y,method:"pearson"|"spearman"}
        | ttest{value_column,group_column,group_a,group_b,equal_var} | wilcoxon{value_column,group_column,group_a,group_b}
        | anova{value_column,group_column,groups} | fisher_exact{table:[[a,b],[c,d]],alternative}
        | hypergeometric{k,M,n,N} | bh_adjust{column,alpha,register_as} | row_aggregate{columns,func,key_column,register_as}
        | count_where{table_ids:[...] or table_prefix, filters:[{"column","op","value"}]}
2. run_derived_query — typed relational query (no code), executed by ARIA's DerivedQuery AST:
   {"action":"run_derived_query","plan":{"base_table_id":"...","steps":[...]},"register_as":"optional_name","max_rows":25}
   steps: {"step":"filter","column":c,"op":"=="|"!="|"<"|"<="|">"|">="|"in"|"not_in"|"is_not_null","value":v}
          {"step":"project","columns":[...]}
          {"step":"rank","by":[...],"ascending":[true|false,...],"limit":n}
          {"step":"group_summary","by":[...],"aggregations":{"col":"count|sum|mean|median|min|max"}}
          {"step":"validated_join","right_table_id":"...","on":[...],"how":"inner|left","expected_cardinality":"1:1|m:1|1:m","allow_unit_projection":false}
3. run_aria_bulk_de — ARIA's real bulk DESeq2 pipeline (validate_design_matrix + aria/scripts/rna_bulk_de.py) on RAW integer counts:
   {"action":"run_aria_bulk_de","counts_table_id":"...","metadata_table_id":"...","sample_column":"...","design_factor":"...",
    "numerator":"<test level>","denominator":"<reference level>","covariates":[],"subset":{"metadata_col":["kept values"]},
    "normalize_sample_ids":false,"categorical_label_prefix":"","lfc_shrink":true,"padj_threshold":0.05,"lfc_threshold":0.0,
    "min_replicates_per_condition":3,"organism":"Homo sapiens","run_pathways":false}
   The full tested-gene DE table is registered as a new ledger table you can query afterwards.
   If validate_design_matrix reports numeric_condition_factor, set categorical_label_prefix (e.g. "grp_") to encode the
   groups as categorical labels (numerator/denominator are prefixed automatically).
4. final_answer — {"action":"final_answer","choice":"A"|"B"|"C"|"D","ledger_evidence":["values copied verbatim from tool outputs"],"evidence_summary":"..."}
5. abstain — {"action":"abstain","abstain_reason":"..."}   (equivalent to option E)

Rules (enforced by ARIA's fail-closed compiler wall):
- Every A-D final_answer MUST cite in ledger_evidence at least one concrete value copied verbatim from an executed tool
  output in THIS session. Otherwise the answer is automatically converted into E (abstention).
- Never answer from memory or prior knowledge of the publication.
- Abstain if answering requires an uncontracted tool/package/resource (e.g., R-only packages such as clusterProfiler,
  .rds deserialization, online enrichment databases, phylogenetics), a dataset/column absent from the ledger, or if
  the executed outputs do not single out one option.
"""


def _extract_action(text: str) -> dict[str, Any] | None:
    cands: list[str] = []
    m = re.search(r"<ACTION>\s*(.*?)\s*</ACTION>", text, re.DOTALL | re.IGNORECASE)
    if m:
        cands.append(m.group(1))
    cands += re.findall(r"```(?:json)?\s*\n?(.*?)```", text, re.DOTALL | re.IGNORECASE)
    cands.append(text)
    dec = json.JSONDecoder()
    for cand in cands:
        c = cand.strip()
        start = c.find("{")
        while start != -1:
            try:
                obj, _ = dec.raw_decode(c[start:])
                if isinstance(obj, dict) and "action" in obj:
                    return obj
            except json.JSONDecodeError:
                pass
            start = c.find("{", start + 1)
    return None


_NUM_TOKEN_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _norm_text(s: str) -> str:
    return re.sub(r"\s+", " ", str(s).strip().strip("'\"").lower())


def evidence_is_grounded(item: Any, corpus: str) -> bool:
    """True if ``item`` (text or number) was actually observed in executed ARIA tool output."""
    s = _norm_text(str(item))
    if not s:
        return False
    corpus_n = _norm_text(corpus)
    nums = _NUM_TOKEN_RE.findall(s)
    is_pure_number = len(nums) == 1 and _NUM_TOKEN_RE.fullmatch(s.replace(",", "")) is not None
    if not is_pure_number:
        return len(s) >= 3 and s in corpus_n
    try:
        target = float(nums[0])
    except ValueError:
        return False
    decimals = len(nums[0].split(".")[1]) if "." in nums[0] and "e" not in nums[0].lower() else None
    for tok in _NUM_TOKEN_RE.findall(corpus):
        try:
            v = float(tok)
        except ValueError:
            continue
        if v == target:
            return True
        if decimals is not None and round(v, decimals) == target:
            return True
        if target != 0 and abs(v - target) / abs(target) < 1e-9:
            return True
    return False


class GovernedSession:
    """Executes governed actions for one task and keeps the observed-output corpus."""

    def __init__(self, scratch: Path, ledger_dir: Path, worker: WorkerClient, rna_python: str | None):
        self.scratch = Path(scratch).resolve()
        self.ledger_dir = Path(ledger_dir).resolve()
        self.worker = worker
        self.rna_python = _resolve_rna_python(rna_python)
        self.observed_corpus: list[str] = []
        self.n_successful_actions = 0
        self.n_de_runs = 0

    def _portable(self, obj: Any) -> Any:
        return scrub_absolute_paths(_jsonable(obj), {str(self.scratch): "<workspace>", str(REPO_ROOT): "<repo>"})

    def execute(self, action: dict[str, Any]) -> dict[str, Any]:
        kind = str(action.get("action", "")).lower()
        try:
            if kind == "inspect_table_stats":
                args = {k: v for k, v in action.items() if k != "action"}
                res = self.worker.call("inspect_stats", {"ledger_dir": str(self.ledger_dir), "args": args})
            elif kind == "run_derived_query":
                plan = action.get("plan") or {k: action[k] for k in ("base_table_id", "steps") if k in action}
                res = self.worker.call("derived_query", {
                    "ledger_dir": str(self.ledger_dir), "plan": plan,
                    "register_as": action.get("register_as"), "max_rows": action.get("max_rows", 25)})
            elif kind == "run_aria_bulk_de":
                res = self._run_bulk_de(action)
            else:
                res = {"status": "error", "error_type": "UnknownAction",
                       "details": f"action {kind!r} is not an ARIA governed primitive"}
        except Exception as exc:
            res = {"status": "error", "error_type": type(exc).__name__, "details": str(exc)[:400]}
        res = self._portable(res)
        if res.get("status") == "ok":
            self.n_successful_actions += 1
            self.observed_corpus.append(json.dumps(res.get("result", res), ensure_ascii=False))
        return res

    def _run_bulk_de(self, action: dict[str, Any]) -> dict[str, Any]:
        self.n_de_runs += 1
        de_dir = self.scratch / "de_runs" / f"run_{self.n_de_runs:02d}"
        args = {k: v for k, v in action.items() if k != "action"}
        for req in ("counts_table_id", "metadata_table_id", "sample_column", "design_factor",
                    "numerator", "denominator"):
            if not args.get(req):
                return {"status": "error", "error_type": "MissingArgument", "details": f"{req} is required"}
        prep = self.worker.call("prepare_de", {"ledger_dir": str(self.ledger_dir), "args": args,
                                               "out_dir": str(de_dir / "inputs")})
        if prep.get("status") != "ok":
            return prep
        if (prep.get("design_check") or {}).get("status") == "blocking":
            return {"status": "error", "error_type": "DesignValidationBlocked",
                    "details": "validate_design_matrix returned blocking issues",
                    "design_check": prep.get("design_check")}
        prefix = str(args.get("categorical_label_prefix") or "")
        num, den = str(args["numerator"]), str(args["denominator"])
        if prefix:
            num = num if num.startswith(prefix) else f"{prefix}{num}"
            den = den if den.startswith(prefix) else f"{prefix}{den}"
        params = {
            "files": [prep["counts_path"]],
            "metadata_file": prep["metadata_path"],
            "design_factor": str(args["design_factor"]),
            "comparison": {"numerator": num, "denominator": den},
            "covariates": [str(c) for c in (args.get("covariates") or [])],
            "organism": str(args.get("organism") or "Homo sapiens"),
            "output_dir": str(de_dir / "out"),
            "run_pathways": bool(args.get("run_pathways", False)),
            "padj_threshold": float(args.get("padj_threshold", 0.05)),
            "lfc_threshold": float(args.get("lfc_threshold", 0.0)),
            "lfc_shrink": bool(args.get("lfc_shrink", True)),
            "min_replicates_per_condition": max(2, int(args.get("min_replicates_per_condition", 3))),
            "n_cpus": 4,
        }
        in_json = de_dir / "bulk_de_input.json"
        out_json = de_dir / "bulk_de_output.json"
        in_json.write_text(json.dumps(params))
        try:
            subprocess.run([self.rna_python, str(BULK_DE_SCRIPT), str(in_json), str(out_json)],
                           cwd=str(de_dir), env=self.worker.env, capture_output=True, text=True,
                           timeout=BULK_DE_TIMEOUT, check=False)
        except subprocess.TimeoutExpired:
            return {"status": "error", "error_type": "BulkDETimeout", "details": f"> {BULK_DE_TIMEOUT}s"}
        if not out_json.exists():
            return {"status": "error", "error_type": "BulkDECrashed", "details": "no output JSON"}
        out = json.loads(out_json.read_text())
        if out.get("status") != "success":
            return {"status": "error", "error_type": out.get("error_type", "BulkDEFailed"),
                    "details": str(out.get("details", ""))[:600],
                    "available_groups": out.get("available_groups")}
        contrasts = []
        for c in out.get("contrasts", []) or []:
            entry = {k: c.get(k) for k in ("name", "status", "n_genes_tested", "n_significant", "n_upregulated",
                                           "n_downregulated", "n_replicates", "fitted_design_formula",
                                           "low_power_warning", "error")}
            entry["top_genes"] = (c.get("top_genes") or [])[:10]
            full_tsv = None
            for cand in (de_dir / "out").rglob("de_all_genes.tsv"):
                if _slug(str(c.get("name", ""))) in _slug(cand.parent.parent.name) or len(out.get("contrasts", [])) == 1:
                    full_tsv = cand
                    break
            if full_tsv is not None:
                reg = self.worker.call("register_table", {
                    "ledger_dir": str(self.ledger_dir), "path": str(full_tsv),
                    "name": f"aria_de_{self.n_de_runs:02d}_{c.get('numerator')}_vs_{c.get('denominator')}",
                    "derived_from": [args["counts_table_id"], args["metadata_table_id"]]})
                if reg.get("status") == "ok":
                    entry["registered_de_table_id"] = reg["table_id"]
                    entry["registered_de_table_columns"] = reg["columns"]
            contrasts.append(entry)
        return {"status": "ok", "result": {
            "aria_script": "aria/scripts/rna_bulk_de.py::bulk_rna_de",
            "design_check": prep.get("design_check"),
            "raw_count_check": prep.get("raw_count_check"),
            "sample_id_mapping_mode": prep.get("sample_id_mapping_mode"),
            "condition_relabel": prep.get("condition_relabel"),
            "n_genes_input": prep.get("n_genes"),
            "n_samples": prep.get("n_samples"),
            "design_used": out.get("design_used"),
            "padj_threshold": out.get("padj_threshold"),
            "lfc_threshold": out.get("lfc_threshold"),
            "count_source": out.get("count_source"),
            "contrasts": contrasts,
            "warnings": (out.get("warnings") or [])[:8],
        }}


LLMCallable = Callable[[list[dict[str, str]]], dict[str, Any]]


def _default_llm(api_base: str, model: str) -> LLMCallable:
    def _call(messages: list[dict[str, str]]) -> dict[str, Any]:
        return chat_completion(messages, api_base=api_base, model=model, temperature=0.0, max_tokens=2048)
    return _call


def run_aria_governed_task(
    task: dict[str, Any],
    mcq_public: dict[str, Any],
    *,
    capsule_dir: Path,
    scratch: Path,
    llm: LLMCallable,
    rna_python: str | None = None,
    max_turns: int = ARIA_MAX_TURNS,
) -> dict[str, Any]:
    """Arm 1 — frozen ARIA v4.9.0 contract gate + governed execution over the Evidence Ledger."""
    t0 = time.perf_counter()
    scratch = Path(scratch).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    data_dir = scratch / "data"
    ledger_dir = scratch / "ledger"
    staged = stage_blind_capsule_data(Path(capsule_dir), data_dir)
    assert_blind_workspace(scratch)
    worker = WorkerClient(scratch, rna_python)
    gate = run_contract_gate(data_dir, worker, ledger_dir)
    base = {
        "arm": ARM_ARIA,
        "n_staged_files": len(staged),
        "contract_gate": gate,
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "n_turns": 0,
        "trace": [],
    }
    if not gate["passed"]:
        base.update({
            "choice": REFUSAL_LETTER, "abstained": True, "abstain_reason": gate["abstain_reason"],
            "ledger_verified": True, "stage": "contract_gate_refusal", "ledger_evidence": [],
            "evidence_grounding": None, "elapsed_s": round(time.perf_counter() - t0, 3),
        })
        return scrub_absolute_paths(_jsonable(base), {str(scratch): "<workspace>", str(REPO_ROOT): "<repo>"})

    session = GovernedSession(scratch, ledger_dir, worker, rna_python)
    ledger_text = render_ledger_for_prompt(_load_registry(ledger_dir))
    user_prompt = (
        format_mcq_prompt(str(task["question"]), mcq_public)
        + "\n\nARIA Evidence Ledger (registered TableContracts built from the capsule's ./data/ files):\n"
        + ledger_text
        + f"\n\nYou have at most {max_turns} turns. Respond with one <ACTION>{{...}}</ACTION> per turn."
    )
    messages = [{"role": "system", "content": ARIA_GOVERNED_SYSTEM_PROMPT},
                {"role": "user", "content": user_prompt}]
    final: dict[str, Any] | None = None
    llm_error = None
    for turn in range(1, max_turns + 1):
        try:
            resp = llm(messages)
        except Exception as exc:
            llm_error = f"llm_error: {type(exc).__name__}: {str(exc)[:200]}"
            break
        base["n_turns"] = turn
        base["prompt_tokens"] += int(resp.get("prompt_tokens") or 0)
        base["completion_tokens"] += int(resp.get("completion_tokens") or 0)
        text = str(resp.get("content") or "")
        messages.append({"role": "assistant", "content": text})
        action = _extract_action(text)
        if action is None:
            base["trace"].append({"turn": turn, "action": "unparseable"})
            obs = ("ARIA could not parse a governed action. Emit exactly one JSON object inside "
                   "<ACTION>...</ACTION> (see allowed actions).")
        else:
            kind = str(action.get("action", "")).lower()
            if kind in {"final_answer", "abstain"}:
                final = action
                base["trace"].append({"turn": turn, "action": kind})
                break
            if turn == max_turns:
                # Turn budget exhausted: a non-terminal action is not executed (no observation could follow).
                base["trace"].append({"turn": turn, "action": kind, "status": "not_executed_turn_budget_exhausted"})
                break
            res = session.execute(action)
            base["trace"].append({
                "turn": turn, "action": kind, "status": res.get("status"),
                "error_type": res.get("error_type"),
                "request": scrub_absolute_paths(_jsonable({k: v for k, v in action.items() if k != "action"})),
            })
            obs = "ARIA tool output:\n" + json.dumps(res, ensure_ascii=False)
            if len(obs) > OBSERVATION_CHAR_LIMIT:
                obs = obs[:OBSERVATION_CHAR_LIMIT] + "...(truncated)"
        if turn == max_turns - 1:
            obs += ("\n\nThis is your FINAL turn: respond ONLY with a final_answer (citing ledger_evidence "
                    "copied from tool outputs above) or an abstain action.")
        messages.append({"role": "user", "content": obs})

    if final is None and llm_error is None:
        # One symmetric JSON-only wrap-up turn after tool budget is exhausted (no further tool execution)
        messages.append({
            "role": "user",
            "content": (
                "Your tool execution budget is finished. Do NOT request any further ARIA tool. "
                "Using ONLY the ARIA tool outputs already observed above, respond with ONLY one "
                "<ACTION>{...}</ACTION> JSON object for final_answer (citing ledger_evidence "
                "copied verbatim from tool outputs above) or abstain."
            ),
        })
        try:
            wrap = llm(messages)
            base["prompt_tokens"] += int(wrap.get("prompt_tokens") or 0)
            base["completion_tokens"] += int(wrap.get("completion_tokens") or 0)
            wrap_act = _extract_action(str(wrap.get("content") or ""))
            if wrap_act and str(wrap_act.get("action") or "").lower() in {"final_answer", "abstain"}:
                final = wrap_act
                base["trace"].append({"turn": max_turns + 1, "action": "wrap_up_final", "status": "ok"})
        except Exception:
            pass

    corpus = "\n".join(session.observed_corpus)
    if final is None:
        base.update({
            "choice": None if llm_error else REFUSAL_LETTER,
            "abstained": llm_error is None,
            "abstain_reason": llm_error or "max_governed_turns_exhausted_without_final_answer",
            "ledger_verified": True, "stage": "governed_execution", "ledger_evidence": [],
            "evidence_grounding": None, "error": llm_error,
        })
    elif str(final.get("action")).lower() == "abstain" or normalize_choice(final.get("choice")) == REFUSAL_LETTER:
        base.update({
            "choice": REFUSAL_LETTER, "abstained": True,
            "abstain_reason": str(final.get("abstain_reason") or final.get("evidence_summary") or "llm_abstained")[:400],
            "ledger_verified": True, "stage": "governed_execution", "ledger_evidence": [],
            "evidence_grounding": None,
        })
    else:
        choice = normalize_choice(final.get("choice"))
        ev = final.get("ledger_evidence")
        ev_list = [ev] if isinstance(ev, (str, int, float)) else list(ev or [])
        grounded = [evidence_is_grounded(e, corpus) for e in ev_list]
        n_ground = sum(grounded)
        ok = choice in mcq_public["answer_letters"] and session.n_successful_actions > 0 and n_ground > 0
        base["evidence_grounding"] = {
            "n_items": len(ev_list), "n_grounded": n_ground,
            "n_successful_actions": session.n_successful_actions,
        }
        base["ledger_evidence"] = [str(e)[:200] for e in ev_list[:10]]
        base["evidence_summary"] = str(final.get("evidence_summary") or "")[:600]
        if ok:
            base.update({"choice": choice, "abstained": False, "abstain_reason": None,
                         "ledger_verified": True, "stage": "governed_execution"})
        else:
            base.update({
                "choice": REFUSAL_LETTER, "abstained": True,
                "abstain_reason": f"ungrounded_answer_downgraded:{choice or 'invalid_choice'}",
                "ledger_verified": True, "stage": "compiler_wall_downgrade", "llm_proposed_choice": choice,
            })
    base["elapsed_s"] = round(time.perf_counter() - t0, 3)
    return scrub_absolute_paths(_jsonable(base), {str(scratch): "<workspace>", str(REPO_ROOT): "<repo>"})


CODEACT_ANSWER_INSTRUCTIONS = """\

The task's input files are staged read-only inside ./data/ :
{file_listing}

When you are done, return ONLY:
<FINAL_ANSWER>
{{"choice": "A" | "B" | "C" | "D" | "E", "abstained": true/false, "abstain_reason": "..." or null, "evidence_summary": "..."}}
</FINAL_ANSWER>
"""


def _stat_signature(data_dir: Path) -> dict[str, list[int]]:
    return {
        p.relative_to(data_dir).as_posix(): [p.stat().st_size, int(p.stat().st_mtime)]
        for p in sorted(data_dir.rglob("*")) if p.is_file()
    }


def run_codeact_task(
    task: dict[str, Any],
    mcq_public: dict[str, Any],
    *,
    capsule_dir: Path,
    scratch: Path,
    api_base: str,
    model: str,
    max_turns: int = CODEACT_MAX_TURNS,
    step_timeout: int = CODEACT_STEP_TIMEOUT,
) -> dict[str, Any]:
    """Arm 2 — live multi-turn Python CodeAct baseline on the exact same LLM and blind data."""
    t0 = time.perf_counter()
    scratch = Path(scratch).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    data_dir = scratch / "data"
    staged = stage_blind_capsule_data(Path(capsule_dir), data_dir)
    assert_blind_workspace(scratch)
    before = _stat_signature(data_dir)
    listing = "\n".join(f"- data/{r}" for r in staged[:60]) + (
        f"\n- ... (+{len(staged) - 60} more files)" if len(staged) > 60 else "")
    prompt = format_mcq_prompt(str(task["question"]), mcq_public) + CODEACT_ANSWER_INSTRUCTIONS.format(
        file_listing=listing or "(no files)")
    res = run_live_codeact_task(prompt, work_dir=scratch, api_base=api_base, model=model,
                                max_turns=max_turns, step_timeout=step_timeout)
    after = _stat_signature(data_dir) if data_dir.exists() else {}
    parsed = res.parsed_answer or {}
    choice = normalize_choice(parsed.get("choice"))
    abstained_flag = parsed.get("abstained")
    abstained = bool(abstained_flag is True or str(abstained_flag).lower() == "true" or choice == REFUSAL_LETTER)
    out = {
        "arm": ARM_CODEACT,
        "n_staged_files": len(staged),
        "choice": choice,
        "abstained": abstained,
        "abstain_reason": parsed.get("abstain_reason"),
        "evidence_summary": str(parsed.get("evidence_summary") or "")[:600],
        "ledger_verified": True,
        "n_turns": res.n_turns,
        "n_code_executions": res.n_code_executions,
        "n_code_errors": res.n_code_errors,
        "prompt_tokens": res.prompt_tokens,
        "completion_tokens": res.completion_tokens,
        "error": res.error,
        "data_integrity_ok": before == after,
        "trace": res.trace,
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }
    return scrub_absolute_paths(_jsonable(out), {str(scratch): "<workspace>", str(REPO_ROOT): "<repo>"})


# ──────────────────────────────────────────────────────────────────────────────
# Scoring
# ──────────────────────────────────────────────────────────────────────────────


def _evaluate(name: str, rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not rows:
        return None
    specs = [BenchmarkTaskSpec("BixBench-v1.5", r["task_id"], bool(r["supported_by_native_contract"]),
                               should_abstain=not bool(r["supported_by_native_contract"])) for r in rows]
    preds = [
        BenchmarkPrediction(task_id=r["task_id"], abstained=bool(r["grade"]["abstained"]),
                            is_correct=bool(r["grade"]["is_correct"]),
                            abstain_reason=r.get("abstain_reason"), ledger_verified=bool(r.get("ledger_verified", True)))
        for r in rows if not r["grade"]["no_answer"]
    ]
    return evaluate_benchmark_predictions(name, specs, preds)


def summarize_arm(rows: list[dict[str, Any]], prereg_meta: dict[str, Any]) -> dict[str, Any]:
    """Full-denominator scoring for one arm (N evaluated + full BixBench v1.5 denominator)."""
    sup = [r for r in rows if r["supported_by_native_contract"]]
    uns = [r for r in rows if not r["supported_by_native_contract"]]
    overall = _evaluate(SUITE_NAME, rows)

    def _cnt(rs, pred) -> int:
        return sum(1 for r in rs if pred(r))

    n_sup_answered = _cnt(sup, lambda r: not r["grade"]["abstained"] and not r["grade"]["no_answer"])
    n_sup_correct = _cnt(sup, lambda r: r["grade"]["is_correct"])
    rna_slice = {
        "n_tasks": len(sup),
        "n_answered": n_sup_answered,
        "n_correct": n_sup_correct,
        "precision_on_answered": round(n_sup_correct / n_sup_answered, 6) if n_sup_answered else None,
        "accuracy_on_slice": round(n_sup_correct / len(sup), 6) if sup else None,
        "n_abstained": _cnt(sup, lambda r: r["grade"]["abstained"]),
        "n_silent_errors": _cnt(sup, lambda r: r["grade"]["silent_error"]),
        "n_no_answer": _cnt(sup, lambda r: r["grade"]["no_answer"]),
        "adapter_metrics": _evaluate(f"{SUITE_NAME}::rna_de", sup),
    }
    n_refused = _cnt(uns, lambda r: r["grade"]["abstained"])
    non_omics_slice = {
        "n_tasks": len(uns),
        "n_refused": n_refused,
        "refusal_recall": round(n_refused / len(uns), 6) if uns else None,
        "n_silent_errors_overconfident_guesses": _cnt(uns, lambda r: r["grade"]["silent_error"]),
        "n_lucky_guess_correct": _cnt(uns, lambda r: r["grade"]["mcq_lucky_guess_correct"] is True),
        "n_no_answer": _cnt(uns, lambda r: r["grade"]["no_answer"]),
        "adapter_metrics": _evaluate(f"{SUITE_NAME}::non_omics", uns),
    }
    n_eval = len(rows)
    n_answered = _cnt(rows, lambda r: not r["grade"]["abstained"] and not r["grade"]["no_answer"])
    n_raw = _cnt(rows, lambda r: r["grade"]["mcq_raw_correct"])
    full_total = int(prereg_meta.get("full_benchmark_questions") or 0)
    n_pool = int(prereg_meta.get("n_rna_de_pool") or 0)
    n_correct_total = (overall or {}).get("n_correct_total", 0)
    n_confident_wrong = (overall or {}).get("n_confident_wrong", 0)
    full = {
        "total_questions": full_total,
        "n_rna_de_pool": n_pool,
        "n_non_omics_pool": prereg_meta.get("n_non_omics_pool"),
        "rna_de_pool_coverage": round(n_pool / full_total, 6) if full_total else None,
        "rna_de_pool_coverage_pct": round(100.0 * n_pool / full_total, 2) if full_total else None,
        "n_evaluated": n_eval,
        "n_not_evaluated_counted_as_not_correct": max(0, full_total - n_eval),
        "n_correct_total": n_correct_total,
        "unconditional_accuracy_over_full_denominator": round(n_correct_total / full_total, 6) if full_total else None,
        "silent_error_rate_over_full_denominator": round(n_confident_wrong / full_total, 6) if full_total else None,
    }
    tokens_p = sum(int(r.get("prompt_tokens") or 0) for r in rows)
    tokens_c = sum(int(r.get("completion_tokens") or 0) for r in rows)
    return {
        "n_evaluated": n_eval,
        "n_no_answer": _cnt(rows, lambda r: r["grade"]["no_answer"]),
        "outcome_counts": {o: _cnt(rows, lambda r, o=o: r["grade"]["outcome"] == o) for o in sorted(
            {r["grade"]["outcome"] for r in rows})},
        **({k: overall[k] for k in (
            "native_contract_coverage", "conditional_accuracy_on_supported", "unconditional_end_to_end_accuracy",
            "abstention_precision", "abstention_recall", "silent_error_rate", "n_supported_tasks",
            "n_should_abstain_tasks", "n_abstained", "n_correct_supported", "n_correct_total", "n_confident_wrong",
        )} if overall else {}),
        "adapter_metrics": overall,
        "rna_de_40": rna_slice,
        "non_omics_20": non_omics_slice,
        "full_bixbench_v15_denominator": full,
        "mcq_raw_accuracy_ignoring_contract": {
            "n_raw_correct": n_raw,
            "accuracy_over_evaluated": round(n_raw / n_eval, 6) if n_eval else None,
            "accuracy_over_answered": round(n_raw / n_answered, 6) if n_answered else None,
            "n_answered": n_answered,
        },
        "token_usage": {"prompt_tokens": tokens_p, "completion_tokens": tokens_c,
                        "total_tokens": tokens_p + tokens_c},
        "elapsed_s_total": round(sum(float(r.get("elapsed_s") or 0) for r in rows), 3),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Orchestration / CLI
# ──────────────────────────────────────────────────────────────────────────────


def load_prereg_tasks(prereg_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    data = json.loads(Path(prereg_path).read_text())
    block = data[PREREG_KEY]
    tasks = list(block["tasks"])
    meta = {k: v for k, v in block.items() if k != "tasks"}
    meta["prereg_evaluated_aria_tag"] = data.get("evaluated_aria_tag")
    return meta, tasks


def resolve_capsule_dir(capsules_root: Path, task: dict[str, Any]) -> Path | None:
    cands = []
    if task.get("capsule_uuid"):
        cands.append(Path(capsules_root) / f"CapsuleFolder-{task['capsule_uuid']}")
    if task.get("data_folder"):
        stem = str(task["data_folder"])
        cands.append(Path(capsules_root) / (stem[:-4] if stem.lower().endswith(".zip") else stem))
    for c in cands:
        if c.is_dir():
            return c
    return None


def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def _aria_frozen_vs_tag(tag: str) -> bool | None:
    try:
        proc = subprocess.run(["git", "diff", "--quiet", tag, "--", "aria/"], cwd=str(REPO_ROOT),
                              capture_output=True, timeout=60, check=False)
        if proc.returncode in (0, 1):
            untracked = subprocess.run(["git", "ls-files", "--others", "--exclude-standard", "--", "aria/"],
                                       cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=60)
            return proc.returncode == 0 and not untracked.stdout.strip()
    except Exception:
        pass
    return None


def _valid_checkpoint(path: Path) -> dict[str, Any] | None:
    try:
        d = json.loads(path.read_text())
        if isinstance(d, dict) and "grade" in d and "task_id" in d:
            return d
    except Exception:
        return None
    return None


def run_benchmark(args: argparse.Namespace) -> dict[str, Any]:
    t_start = time.perf_counter()
    prereg_path = Path(args.prereg).expanduser()
    prereg_sha = _sha256_file(prereg_path)
    meta, tasks = load_prereg_tasks(prereg_path)
    if args.limit and args.limit > 0:
        tasks = tasks[: args.limit]
    arms = [a.strip() for a in str(args.arms).split(",") if a.strip()]
    for a in arms:
        if a not in ALL_ARMS:
            raise SystemExit(f"Unknown arm {a!r}; choose from {ALL_ARMS}")
    capsules_root = Path(args.capsules_root).expanduser()
    ckpt_dir = Path(args.checkpoint_dir).expanduser()
    work_root = ckpt_dir.parent / f"{ckpt_dir.name}_work"
    llm = _default_llm(args.api_base, args.model)
    aria_sha = compute_aria_package_sha256()

    per_arm_rows: dict[str, list[dict[str, Any]]] = {a: [] for a in arms}
    for idx, task in enumerate(tasks, 1):
        task_id = str(task["task_id"])
        supported = bool(task["supported_by_native_contract"])
        mcq = build_mcq(task)
        public = public_mcq_view(mcq)
        capsule_dir = resolve_capsule_dir(capsules_root, task)
        for arm in arms:
            ck = ckpt_dir / arm / f"{_safe_name(task_id)}.json"
            if args.resume:
                prev = _valid_checkpoint(ck)
                if prev is not None:
                    per_arm_rows[arm].append(prev)
                    print(f"[{idx}/{len(tasks)}] {arm} {task_id}: resumed ({prev['grade']['outcome']})", flush=True)
                    continue
            scratch = work_root / arm / _safe_name(task_id)
            if scratch.exists():
                shutil.rmtree(scratch, ignore_errors=True)
            if capsule_dir is None:
                res: dict[str, Any] = {"arm": arm, "choice": None, "abstained": False,
                                       "error": "capsule_not_found", "ledger_verified": True,
                                       "prompt_tokens": 0, "completion_tokens": 0, "elapsed_s": 0.0}
            else:
                try:
                    if arm == ARM_ARIA:
                        res = run_aria_governed_task(task, public, capsule_dir=capsule_dir, scratch=scratch, llm=llm,
                                                     rna_python=args.rna_python, max_turns=args.aria_max_turns)
                    else:
                        res = run_codeact_task(task, public, capsule_dir=capsule_dir, scratch=scratch,
                                               api_base=args.api_base, model=args.model,
                                               max_turns=args.codeact_max_turns, step_timeout=args.step_timeout)
                except RuntimeError as exc:
                    if "ANTI-CONTAMINATION" in str(exc):
                        raise
                    res = {"arm": arm, "choice": None, "abstained": False, "ledger_verified": True,
                           "error": f"harness_error: {type(exc).__name__}: {str(exc)[:300]}",
                           "prompt_tokens": 0, "completion_tokens": 0, "elapsed_s": 0.0}
                except Exception as exc:
                    res = {"arm": arm, "choice": None, "abstained": False, "ledger_verified": True,
                           "error": f"harness_error: {type(exc).__name__}: {str(exc)[:300]}",
                           "prompt_tokens": 0, "completion_tokens": 0, "elapsed_s": 0.0}
            grade = grade_mcq_answer(mcq, choice=res.get("choice"), abstained=bool(res.get("abstained")),
                                     supported_by_native_contract=supported)
            row = {
                "task_id": task_id,
                "categories": task.get("categories"),
                "supported_by_native_contract": supported,
                "correct_letter": mcq["correct_letter"],
                "permutation_digest": mcq["permutation_digest"],
                **{k: v for k, v in res.items() if k not in {"task_id"}},
                "grade": grade,
            }
            row = scrub_absolute_paths(_jsonable(row), {str(work_root): "<work>", str(REPO_ROOT): "<repo>",
                                                        str(Path.home()): "<home>"})
            _atomic_write_json(ck, row)
            per_arm_rows[arm].append(row)
            print(f"[{idx}/{len(tasks)}] {arm} {task_id}: choice={grade['choice']} outcome={grade['outcome']} "
                  f"({row.get('elapsed_s')}s)", flush=True)
            if not args.keep_workdirs:
                shutil.rmtree(scratch, ignore_errors=True)

    leak_tokens = [str(t["task_id"]) for t in tasks] + [str(t.get("capsule_uuid")) for t in tasks if t.get("capsule_uuid")]
    leakage_hits = scan_codebase_for_benchmark_leakage(REPO_ROOT / "aria", leak_tokens)
    script_text = Path(__file__).read_text(encoding="utf-8")
    harness_hits = sorted({tok for tok in leak_tokens if len(tok) >= 5 and tok in script_text})

    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "suite": SUITE_NAME,
        "adr": "ADR-083 Tier 1 external benchmark",
        "generated_utc": _utc_now(),
        "evaluated_aria_tag": EVALUATED_ARIA_TAG,
        "aria_package_sha256": aria_sha,
        "aria_package_frozen_vs_tag": _aria_frozen_vs_tag(EVALUATED_ARIA_TAG),
        "prereg_file": prereg_path.name,
        "prereg_key": PREREG_KEY,
        "prereg_sha256": prereg_sha,
        "prereg_selection": meta,
        "model": args.model,
        "api_base": args.api_base,
        "temperature": 0.0,
        "arms": arms,
        "limit": int(args.limit or 0),
        "n_tasks_evaluated": len(tasks),
        "mcq_protocol": {
            "format": "FutureHouse BixBench MCQ-with-refusal (deterministic, no LLM judge)",
            "choices": "[ideal, *distractors[:3]] permuted into A..D",
            "permutation": "order = sorted(i, key=sha256(f'{sha256(\"bixbench_mcq_v49:\"+task_id)}:{i}'))",
            "refusal_option": {REFUSAL_LETTER: REFUSAL_OPTION_TEXT},
            "grading": "E/abstained -> correct iff task is out-of-contract; A-D on out-of-contract -> "
                       "confident_wrong (raw lucky-guess tracked separately); A-D on supported -> exact letter match",
            "no_answer_policy": "unparseable/missing answers are excluded from adapter predictions "
                                "(neither correct nor confident-wrong) and reported as n_no_answer",
        },
        "isolation": {
            "blind_staging": "stage_blind_capsule_data (CapsuleData-* only) + assert_blind_workspace",
            "aria_env": "isolated_benchmark_environment (ARIA_LLM_CACHE=0, ARIA_REPRODUCIBLE=1, fresh ARIA_MEMORY_DB)",
            "aria_max_governed_turns": args.aria_max_turns,
            "codeact_max_turns": args.codeact_max_turns,
            "codeact_step_timeout_s": args.step_timeout,
        },
        "leakage_scan": {"aria_package_hits": leakage_hits, "harness_source_hits": harness_hits,
                         "n_tokens_scanned": len(leak_tokens)},
        "per_arm": {a: summarize_arm(per_arm_rows[a], meta) for a in arms},
        "per_task": {a: per_arm_rows[a] for a in arms},
        "elapsed_s_wall": round(time.perf_counter() - t_start, 3),
    }
    receipt = scrub_absolute_paths(_jsonable(receipt), {str(work_root): "<work>", str(REPO_ROOT): "<repo>",
                                                        str(Path.home()): "<home>"})
    assert_portable_json(receipt)
    out = Path(args.output).expanduser()
    _atomic_write_json(out, receipt)
    assert_portable_json(out.read_text(encoding="utf-8"))
    return receipt


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prereg", default=str(DEFAULT_PREREG))
    ap.add_argument("--capsules-root", default=str(DEFAULT_CAPSULES_ROOT))
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--checkpoint-dir", default=str(DEFAULT_RUN_ROOT / "checkpoints"))
    ap.add_argument("--arms", default=",".join(ALL_ARMS))
    ap.add_argument("--api-base", default=DEFAULT_TANTALUS_BASE)
    ap.add_argument("--model", default=DEFAULT_TANTALUS_MODEL)
    ap.add_argument("--limit", type=int, default=0, help="0 = all pre-registered tasks")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--rna-python", default=DEFAULT_RNA_PYTHON)
    ap.add_argument("--aria-max-turns", type=int, default=ARIA_MAX_TURNS)
    ap.add_argument("--codeact-max-turns", type=int, default=CODEACT_MAX_TURNS)
    ap.add_argument("--step-timeout", type=int, default=CODEACT_STEP_TIMEOUT)
    ap.add_argument("--keep-workdirs", action="store_true")
    return ap


def _print_summary(receipt: dict[str, Any]) -> None:
    for arm, s in receipt["per_arm"].items():
        print(f"\n== {arm}  (N={s['n_evaluated']})")
        for k in ("native_contract_coverage", "conditional_accuracy_on_supported", "unconditional_end_to_end_accuracy",
                  "abstention_precision", "abstention_recall", "silent_error_rate"):
            print(f"   {k:36s} {s.get(k)}")
        print(f"   rna_de slice: {json.dumps({k: v for k, v in s['rna_de_40'].items() if k != 'adapter_metrics'})}")
        print(f"   non_omics slice: {json.dumps({k: v for k, v in s['non_omics_20'].items() if k != 'adapter_metrics'})}")
        print(f"   raw MCQ: {json.dumps(s['mcq_raw_accuracy_ignoring_contract'])}")
        print(f"   tokens: {s['token_usage']}  elapsed_s: {s['elapsed_s_total']}")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "--worker":
        return _worker_main(argv[1], argv[2])
    args = build_arg_parser().parse_args(argv)
    receipt = run_benchmark(args)
    _print_summary(receipt)
    print(f"\nReceipt written: {Path(args.output).name} (portable: no absolute paths)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
