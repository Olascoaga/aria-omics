"""
ARIA RNA CellTypist Annotation Script
--------------------------------------
Database-backed cell type annotation. Implements the "code guarantees"
half of the LLM-proposes-code-guarantees pattern: CellTypist gives a
defensible label per cell from a pretrained model, the agent's LLM only
reinterprets / refines on top of that label.

Executed inside aria-rna-env by EnvironmentManager (celltypist 1.7.1+).

Input params:
    data_path:      str  — path to clustered .h5ad (log1p(CPM), sum=1e4)
    organism:       str  — "Homo sapiens" | "Mus musculus" (default: Homo sapiens)
    model:          str  (optional) — explicit CellTypist model name (.pkl)
                                       If omitted, picked from tissue_hint.
    tissue_hint:    str  (optional) — "pbmc" | "blood" | "brain" | "fetal"
                                       | "kidney" | "lung" | "intestine"
                                       Picks a default model if `model` is empty.
    allow_default_immune_model: bool (optional) — explicit ACK allowing the
                                       immune fallback when no model/hint matches.
    cluster_col:    str  (optional) — obs column for cluster IDs (default: "leiden")
    majority_voting: bool (optional) — collapse labels per cluster (default: True)
    over_clustering: str (optional) — column for majority voting groups
                                       (default: cluster_col)
    output_dir:     str  (optional)

Output:
    {
      "status":              "success" | "skipped" | "error",
      "model_used":          str,
      "predictions_path":    str   — CSV with per-cell labels
      "per_cluster": {
          cluster_id: {
              "label":        str,
              "frequency":    float,    # fraction of cluster with this label
              "n_cells":      int,
              "alt_labels":   [{label, frequency}, ...]   # runner-ups
          },
          ...
      },
      "output_path":         str   — annotated .h5ad (with cell_type_celltypist obs)
    }
"""

from __future__ import annotations
import sys, os
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


# Tissue hint → default CellTypist model.
# These models ship with celltypist's hub; the script downloads them on
# first use (~30MB each, cached under ~/.celltypist/data/models).
_DEFAULT_MODELS = {
    "pbmc":      "Immune_All_Low.pkl",     # 98 fine-grained immune types
    "blood":     "Immune_All_Low.pkl",
    "immune":    "Immune_All_Low.pkl",
    "broad":     "Immune_All_High.pkl",    # 32 broad immune types
    "fetal":     "Pan_Fetal_Human.pkl",
    "brain":     "Adult_Human_PrefrontalCortex.pkl",
    "cortex":    "Adult_Human_PrefrontalCortex.pkl",
    "kidney":    "Adult_Human_Kidney.pkl",
    "lung":      "Human_Lung_Atlas.pkl",
    "intestine": "Cells_Intestinal_Tract.pkl",
    "gut":       "Cells_Intestinal_Tract.pkl",
    "skin":      "Adult_Human_Skin.pkl",
}

_MOUSE_MODELS = {
    "brain":  "Mouse_Whole_Brain.pkl",
    "immune": "Mouse_Isocortex_Hippocampus.pkl",
}

# The catalog fallback when no explicit model and no usable tissue hint exist.
_IMMUNE_DEFAULT_MODEL = "Immune_All_Low.pkl"


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _resolve_celltypist_model(model_explicit, organism, tissue_hint):
    """N-ANNO3: resolve the CellTypist model and DISCLOSE a silent default.

    Returns ``(model_name, model_source, warning)`` where ``model_source`` is one
    of ``"explicit"``, ``"tissue_hint"``, ``"default_immune_fallback"``.

    Annotating an arbitrary tissue with the immune-only default
    (``Immune_All_Low.pkl``) produces biologically wrong labels that then define
    the per-cell-type DE groupings, so a fallback (no explicit model AND no
    matching tissue hint) is flagged with a ``warning`` string rather than chosen
    silently. The model is NOT forced — the caller decides; only the assumption
    is surfaced. Pure / dependency-light for unit testing.
    """
    if model_explicit:
        return str(model_explicit), "explicit", None

    org = (organism or "").lower()
    is_mouse = "musculus" in org or org == "mus musculus"
    catalog = _MOUSE_MODELS if is_mouse else _DEFAULT_MODELS
    hint = (tissue_hint or "").lower().strip()
    if hint and hint in catalog:
        return catalog[hint], "tissue_hint", None

    warning = (
        f"No CellTypist model or tissue hint was provided, so annotation fell "
        f"back to the immune default '{_IMMUNE_DEFAULT_MODEL}'. If this dataset "
        f"is not immune/PBMC tissue, the cell-type labels may be biologically "
        f"wrong, and they define the per-cell-type DE groupings. Provide an "
        f"explicit model or tissue hint to remove this assumption."
    )
    if is_mouse:
        warning += (
            " The dataset is mouse but the default is a HUMAN immune model — a "
            "species mismatch on top of the tissue mismatch."
        )
    return _IMMUNE_DEFAULT_MODEL, "default_immune_fallback", warning


def _celltypist_model_cached(model_name: str) -> bool:
    """Return True when the requested model is already available locally."""
    model_path = Path(str(model_name)).expanduser()
    if model_path.is_file():
        return True

    for cache_dir in (
        Path.home() / ".celltypist" / "data" / "models",
        Path.home() / ".celltypist" / "models",
    ):
        if (cache_dir / str(model_name)).is_file():
            return True
    return False


def _extract_assigned_confidences(prob_df, obs_names, assigned_labels):
    """T7: per-cell probability of each assigned label + a degradation flag.

    Returns ``(conf_scores, confidence_available, confidence_extraction_failed)``.

      confidence_available          ``prob_df`` was produced at all (the model
                                    emitted a probability matrix).
      confidence_extraction_failed  ``prob_df`` was produced but yielded ZERO
                                    finite per-cell scores — an exception, or
                                    label columns that do not align with the
                                    assigned labels. This is a SILENT degradation
                                    that previously collapsed to the same all-NaN
                                    result as "no probabilities available", so
                                    downstream treated it as full confidence and
                                    the N-ANNO1 caveat/cap never fired.

    Pure / dependency-light (pandas DataFrame in, lists out) so the misaligned-
    columns case is unit-testable without celltypist or a model download.
    """
    n = len(assigned_labels)
    conf_scores = [float("nan")] * n
    if prob_df is None:
        return conf_scores, False, False
    try:
        prob_aligned = prob_df.reindex(obs_names)
        col_index = {str(c): i for i, c in enumerate(prob_aligned.columns)}
        prob_mat = prob_aligned.to_numpy()
        conf_scores = [
            float(prob_mat[i, col_index[lab]])
            if lab in col_index and prob_mat[i, col_index[lab]] == prob_mat[i, col_index[lab]]
            else float("nan")
            for i, lab in enumerate(map(str, assigned_labels))
        ]
    except Exception:
        conf_scores = [float("nan")] * n
    # Degradation = probabilities WERE produced but no finite score survived.
    extraction_failed = n > 0 and not any(c == c for c in conf_scores)
    return conf_scores, True, extraction_failed


def _summarize_annotation_confidence(cluster_ids, raw_labels, assigned_labels,
                                     conf_scores, *, top_k: int = 3) -> dict:
    """N-ANNO2: build a GENUINE per-cluster annotation-confidence summary.

    The confidence proxy must not be structurally 1.0. With majority voting
    over the clustering partition, the assigned (collapsed) label is constant
    within each cluster, so a `frequency` computed over assigned labels is
    always 1.0. We instead measure:

      label           majority of the assigned (final) labels in the cluster.
      frequency       fraction of cells whose RAW per-cell CellTypist label
                      matches that majority label — the genuine within-cluster
                      agreement, NOT the collapsed majority.
      mean/median_confidence  central tendency of the model's per-cell
                      probability for the assigned label (None when no finite
                      probabilities are available — never faked as 0 or 1).
      alt_labels      runner-up RAW labels with their cluster fractions, so the
                      disagreement that majority voting hides stays visible.

    Pure / dependency-light (lists or arrays only) so it is unit-testable
    without celltypist or a model download.
    """
    from collections import Counter, defaultdict
    import math as _math

    by_cluster_raw: dict = defaultdict(list)
    by_cluster_assigned: dict = defaultdict(list)
    by_cluster_conf: dict = defaultdict(list)
    for cl, raw, assigned, conf in zip(
        cluster_ids, raw_labels, assigned_labels, conf_scores
    ):
        cl = str(cl)
        by_cluster_raw[cl].append(str(raw))
        by_cluster_assigned[cl].append(str(assigned))
        try:
            c = float(conf)
        except (TypeError, ValueError):
            c = _math.nan
        by_cluster_conf[cl].append(c)

    per_cluster: dict = {}
    for cl in sorted(by_cluster_raw, key=lambda x: (len(x), x)):
        raws = by_cluster_raw[cl]
        assigned = by_cluster_assigned[cl]
        confs = by_cluster_conf[cl]
        n = len(raws)
        if n == 0:
            continue

        # Final label = majority of the assigned (post-MV) labels.
        majority_label = Counter(assigned).most_common(1)[0][0]

        # Genuine agreement = fraction of RAW per-cell calls matching it.
        raw_counts = Counter(raws)
        frequency = round(raw_counts.get(majority_label, 0) / n, 3)

        finite = [c for c in confs if c == c and _math.isfinite(c)]
        if finite:
            finite_sorted = sorted(finite)
            mid = len(finite_sorted) // 2
            median_conf = (
                finite_sorted[mid] if len(finite_sorted) % 2
                else (finite_sorted[mid - 1] + finite_sorted[mid]) / 2.0
            )
            mean_confidence = round(sum(finite) / len(finite), 3)
            median_confidence = round(float(median_conf), 3)
        else:
            mean_confidence = None
            median_confidence = None

        alt_labels = [
            {"label": str(label), "frequency": round(c / n, 3)}
            for label, c in raw_counts.most_common()
            if str(label) != str(majority_label)
        ][: max(0, top_k - 1)]

        per_cluster[cl] = {
            "label":             majority_label,
            "frequency":         frequency,
            "mean_confidence":   mean_confidence,
            "median_confidence": median_confidence,
            "n_cells":           n,
            "alt_labels":        alt_labels,
        }
    return per_cluster


def rna_celltypist(params: dict) -> dict:
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/aria_numba_cache")

    data_path       = params["data_path"]
    organism        = params.get("organism", "Homo sapiens").lower()
    model_explicit  = params.get("model")
    tissue_hint     = (params.get("tissue_hint") or "").lower().strip()
    allow_default   = _as_bool(params.get("allow_default_immune_model"))
    cluster_col     = params.get("cluster_col", "leiden")
    majority_voting = bool(params.get("majority_voting", True))
    output_dir      = params.get("output_dir", str(Path(data_path).parent))
    out_dir         = Path(output_dir)

    # Resolve model (N-ANNO3: disclose a silent immune-default fallback).
    model_name, model_source, model_warning = _resolve_celltypist_model(
        model_explicit, organism, tissue_hint
    )

    if model_source == "default_immune_fallback" and not allow_default:
        return {
            "status":        "error",
            "error_type":    "DefaultImmuneModelRequiresAck",
            "details":       (
                f"CellTypist would use the immune default '{model_name}' because "
                "no explicit model or matching tissue hint was provided. This "
                "can create biologically wrong cell-type labels for non-PBMC "
                "tissue and must not define downstream cell-type DE groupings. "
                "Provide an explicit CellTypist model or matching tissue_hint, "
                "or set allow_default_immune_model=true to acknowledge this "
                "assumption explicitly."
            ),
            "model_used":    model_name,
            "model_source":  model_source,
            "model_warning": model_warning,
            "requires_ack":  True,
            "ack_param":     "allow_default_immune_model",
            "tissue_hint":   tissue_hint or None,
        }

    try:
        import celltypist
        from celltypist import models
    except ImportError:
        return {
            "status":     "error",
            "error_type": "CellTypistMissing",
            "details":    "celltypist not installed in aria-rna-env. "
                          "Run: pip install celltypist",
        }

    # Download only when the model is missing locally. Under W-PRIV,
    # ARIA_AIR_GAPPED governs this model-hub egress too, not only LLM calls.
    if not _celltypist_model_cached(model_name):
        try:
            from aria.utils.privacy import EgressBlocked, assert_egress_allowed
            assert_egress_allowed("celltypist_model_hub")
        except EgressBlocked as e:
            return {
                "status":        "error",
                "error_type":    "EgressBlocked",
                "details":       str(e),
                "model_used":    model_name,
                "model_source":  model_source,
                "model_warning": model_warning,
            }

        try:
            models.download_models(model=model_name, force_update=False)
        except Exception as e:
            return {
                "status":     "error",
                "error_type": "ModelDownloadFailed",
                "details":    f"Could not download CellTypist model "
                              f"'{model_name}': {e}",
            }

    import scanpy as sc
    import numpy as np
    from aria.utils.safe_h5ad import read_h5ad

    adata = read_h5ad(data_path)

    def _sample_matrix(a):
        n = min(int(a.n_obs), 2000)
        if n <= 0:
            return a.X
        idx = np.linspace(0, int(a.n_obs) - 1, n, dtype=int)
        return a.X[idx, :]

    def _looks_lognorm_10k(a) -> bool:
        try:
            x = _sample_matrix(a)
            if hasattr(x, "toarray"):
                x = x.toarray()
            if np.nanmin(x) < -1e-8 or np.nanmax(x) > 30:
                return False
            sums = np.expm1(x).sum(axis=1)
            sums = sums[np.isfinite(sums) & (sums > 0)]
            if len(sums) == 0:
                return False
            med = float(np.median(sums))
            return 7000 <= med <= 13000
        except Exception:
            return False

    # CellTypist validates both .X and .raw.X. Clustering leaves .X scaled
    # after PCA, so use a clean log1p-normalized AnnData copy and drop .raw
    # before annotation to avoid validating stale scaled/raw matrices.
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        raw_adata.obs = adata.obs.copy()
        if _looks_lognorm_10k(raw_adata):
            adata = raw_adata

    if not _looks_lognorm_10k(adata):
        if adata.X.min() < 0:
            return {
                "status": "error",
                "error_type": "InvalidExpressionMatrix",
                "details": (
                    "Expression matrix has negative/scaled values and no "
                    "valid log1p-normalized raw matrix for CellTypist."
                ),
            }
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    adata.raw = None

    # Predict
    pred = celltypist.annotate(
        adata,
        model=model_name,
        majority_voting=majority_voting,
        over_clustering=(cluster_col if cluster_col in adata.obs.columns else None),
    )
    pred_df = pred.predicted_labels.copy()  # cells × {predicted_labels, ...}

    # Write per-cell labels back into adata.obs
    prediction_label_col = (
        "majority_voting" if "majority_voting" in pred_df.columns
        else "predicted_labels"
    )
    assigned_labels = pred_df[prediction_label_col].astype(str).values
    adata.obs["cell_type_celltypist"] = assigned_labels

    # Raw per-cell predicted labels (PRE majority voting). When majority voting
    # is on, `cell_type_celltypist` is collapsed to one label per cluster, so a
    # frequency computed over it is structurally 1.0 (N-ANNO2). We keep the raw
    # calls to measure genuine within-cluster agreement.
    raw_labels = pred_df["predicted_labels"].astype(str).values

    # Per-cell probability of the ASSIGNED label, from CellTypist's probability
    # matrix (cells × labels). This is the model's own confidence, which was
    # previously discarded entirely (N-ANNO1 surfaces it downstream).
    # T7: distinguish "no probabilities available" from "probabilities existed
    # but extraction degraded" (exception / misaligned label columns). Both used
    # to collapse to silent all-NaN, hiding the degradation from N-ANNO1.
    conf_scores, confidence_available, confidence_extraction_failed = (
        _extract_assigned_confidences(
            getattr(pred, "probability_matrix", None),
            adata.obs_names,
            list(map(str, assigned_labels)),
        )
    )

    # Aggregate to a GENUINE per-cluster confidence summary (pure helper).
    per_cluster: dict = {}
    if cluster_col in adata.obs.columns:
        per_cluster = _summarize_annotation_confidence(
            adata.obs[cluster_col].astype(str).tolist(),
            list(map(str, raw_labels)),
            list(map(str, assigned_labels)),
            conf_scores,
        )

    # Persist outputs
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_csv = str(out_dir / "celltypist_predictions.csv")
    pred_df.to_csv(pred_csv, index=True)

    annotated_path = str(out_dir / "annotated.h5ad")
    adata.write_h5ad(annotated_path)

    return {
        "status":           "success",
        "model_used":       model_name,
        "model_source":     model_source,
        "model_warning":    model_warning,
        "label_col":        "cell_type_celltypist",
        "prediction_label_col": prediction_label_col,
        "n_cells":          int(adata.n_obs),
        "n_unique_labels":  int(adata.obs["cell_type_celltypist"].nunique()),
        "per_cluster":      per_cluster,
        "confidence_available":         confidence_available,
        "confidence_extraction_failed": confidence_extraction_failed,
        "predictions_path": pred_csv,
        "output_path":      annotated_path,
    }


if __name__ == "__main__":
    run_script(rna_celltypist)
