"""scATAC P3b — SnapATAC2 reference driver (PRODUCES the outputs P3a scores).

This is the infra-gated drop-in companion to the P3a concordance harness
(``benchmark_scatac_concordance.py``). P3a CONSUMES a reference tool's outputs;
this driver PRODUCES them by running SnapATAC2 inside ``aria-bench-env`` and
emitting the same JSON schema P3a reads:

    {"clusters": [<per-cell label>, ...],   # ALIGNED to ARIA's cell order
     "da_peaks": ["chr1:100-200", ...],     # reference DA-peak ids/coords
     "motifs":   ["MOTIF_A", ...]}          # reference motif ranking

P3b is gated on Samael: it requires installing SnapATAC2 (+ ArchR/Signac for the
other comparators) and staging public datasets. That install/download is NOT done
here. Until ``aria-bench-env`` has SnapATAC2, this driver returns a structured
honest ``not_run`` instead of fabricating reference outputs (ADR-002). The actual
SnapATAC2 calls live in ``_run_reference_pipeline`` — the single P3b fill-in slot;
its intended API is documented inline so wiring it does not touch the IPC/contract
boundary or the P3a scorer.

scATAC stays ``alpha`` + ``requires_ack``. This driver is benchmark scaffolding,
never a dispatch path for a user run and never a modality promotion.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from aria.scripts._base import run_script


def _run_reference_pipeline(
    input_path: Path, params: dict[str, Any]
) -> dict[str, Any]:
    """P3b fill-in slot: run SnapATAC2 and return the P3a output schema.

    Intended implementation (when ``aria-bench-env`` has SnapATAC2):

      import snapatac2 as snap
      adata = snap.read(str(input_path))                  # or read_h5ad
      snap.pp.select_features(adata)
      snap.tl.spectral(adata)                             # LSI / spectral embedding
      snap.pp.knn(adata); snap.tl.leiden(adata, resolution=params["resolution"])
      clusters = adata.obs["leiden"].astype(str).tolist() # ALIGN to ARIA's order
      da_peaks = _snapatac2_marker_peaks(adata)           # reference DA peaks
      motifs   = _snapatac2_motif_enrichment(adata)       # ranked reference motifs
      return {"clusters": clusters, "da_peaks": da_peaks, "motifs": motifs}

    The cell order MUST match the ARIA outputs (same input matrix, same cell
    filter) so cluster ARI/NMI compares the same cells. Until the tool is wired,
    raising keeps the no-fabrication contract; the caller turns it into a
    structured ``not_run``.
    """
    raise NotImplementedError(
        "P3b SnapATAC2 reference pipeline is not wired yet. Install SnapATAC2 in "
        "aria-bench-env and implement _run_reference_pipeline (see docstring); "
        "ARIA refuses to fabricate reference outputs."
    )


def benchmark_scatac_reference_snapatac2(params: dict) -> dict:
    tool = "SnapATAC2"
    input_path_str = params.get("input_path")
    if not input_path_str or not Path(input_path_str).exists():
        return {
            "status": "not_run",
            "error_type": "MissingDependency",
            "tool": tool,
            "reason": "input_path missing or unreadable (expected a .h5ad / "
                      ".h5mu peak matrix or fragments file matching the ARIA run).",
        }

    try:
        import snapatac2  # noqa: F401
    except ImportError:
        return {
            "status": "not_run",
            "error_type": "MissingDependency",
            "tool": tool,
            "reason": (
                "SnapATAC2 is not installed in this environment. The P3b reference "
                "driver runs in aria-bench-env, whose install (SnapATAC2/ArchR/"
                "Signac) + public datasets are infra-gated on Samael; ARIA does "
                "not fabricate reference outputs."),
        }

    try:
        outputs = _run_reference_pipeline(Path(input_path_str), params)
    except NotImplementedError as exc:
        return {
            "status": "not_run",
            "error_type": "NotImplemented",
            "tool": tool,
            "reason": str(exc),
        }

    result: dict[str, Any] = {
        "status": "success",
        "tool": tool,
        "clusters": outputs.get("clusters"),
        "da_peaks": outputs.get("da_peaks"),
        "motifs": outputs.get("motifs"),
    }
    output_json = params.get("output_json")
    if output_json:
        out = Path(output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
        result["artifacts"] = {"output_json": str(out)}
    return result


if __name__ == "__main__":
    run_script(benchmark_scatac_reference_snapatac2)
