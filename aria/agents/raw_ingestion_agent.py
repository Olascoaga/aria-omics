"""Raw input ingestion agent for canonical ARIA workspace artifacts."""

from __future__ import annotations

from pathlib import Path

from aria.agents.base_agent import BaseAgent
from aria.memory.memory import ARIAMemory
from aria.utils.provenance import hash_file
from aria.utils.raw_ingestion import (
    discover_10x_mtx_triplets,
    ingest_10x_mtx_triplet,
    scan_fastq_plan,
)


class RawIngestionAgent(BaseAgent):
    """Convert supported raw scRNA inputs into canonical workspace h5ads."""

    name = "raw_ingestion_agent"
    description = "Deterministically ingests raw scRNA matrices/FASTQs."

    def __init__(self, memory: ARIAMemory, llm=None, api_key: str = None):
        super().__init__(memory, llm=llm, api_key=api_key)

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx = context.get("exp_context", {}) or {}
        data_dir = Path(exp_ctx.get("data_dir") or context.get("data_dir") or ".")
        workspace = Path.home() / ".aria" / "workspace" / experiment_id / "ingested"
        workspace.mkdir(parents=True, exist_ok=True)

        self.publish_status(
            experiment_id,
            "RawIngestionAgent: scanning raw inputs...",
            0.0,
        )

        records = []
        generated_h5ads = []
        errors = []

        triplets = discover_10x_mtx_triplets(data_dir)
        for i, triplet in enumerate(triplets):
            self.publish_status(
                experiment_id,
                f"Converting 10X matrix {triplet.sample_id} to h5ad...",
                0.15 + (0.55 * (i / max(len(triplets), 1))),
            )
            try:
                rec = ingest_10x_mtx_triplet(
                    triplet.directory,
                    workspace,
                    sample_id=triplet.sample_id,
                )
                records.append(rec)
                generated_h5ads.append(rec["output_h5ad"])
            except Exception as exc:
                errors.append({
                    "mode": "10x_mtx",
                    "source_directory": str(triplet.directory),
                    "error_type": type(exc).__name__,
                    "details": str(exc),
                })

        fastq_plan = scan_fastq_plan(data_dir)
        if fastq_plan.get("fastq_count", 0):
            records.append(fastq_plan)

        if errors:
            return {
                "status": "error",
                "error_type": "RawIngestionFailed",
                "details": "; ".join(e["details"] for e in errors[:3]),
                "records": records,
                "errors": errors,
                "output_h5ads": generated_h5ads,
            }

        if not records:
            self.publish_status(
                experiment_id,
                "RawIngestionAgent: no raw-ingestion inputs detected.",
                1.0,
            )
            return {
                "status": "skipped",
                "reason": "No supported raw-ingestion inputs detected.",
                "records": [],
                "output_h5ads": [],
            }

        input_files = list(exp_ctx.get("input_files", []) or [])
        for rec in records:
            if rec.get("mode") != "10x_mtx":
                continue
            output = Path(rec["output_h5ad"])
            input_files.append({
                "modality": "scRNA_ingested_h5ad",
                "path": str(output),
                "size_bytes": output.stat().st_size,
                "sha256": hash_file(output),
                "source_mode": rec.get("mode"),
            })

        modalities = dict(exp_ctx.get("modalities", {}) or {})
        if generated_h5ads:
            modalities["scRNA"] = generated_h5ads

        exp_ctx["modalities"] = modalities
        exp_ctx["input_files"] = input_files
        exp_ctx["raw_ingestion"] = records

        self.publish_status(
            experiment_id,
            f"RawIngestionAgent: generated {len(generated_h5ads)} h5ad file(s).",
            1.0,
        )
        return {
            "status": "done",
            "records": records,
            "output_h5ads": generated_h5ads,
            "exp_context_updates": {
                "modalities": modalities,
                "input_files": input_files,
                "raw_ingestion": records,
            },
        }
