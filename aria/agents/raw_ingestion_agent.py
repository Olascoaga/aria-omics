"""Raw input ingestion agent for canonical ARIA workspace artifacts."""

from __future__ import annotations

import json
from pathlib import Path

from aria.agents.base_agent import BaseAgent
from aria.memory.memory import ARIAMemory
from aria.utils.provenance import hash_file
from aria.utils.raw_ingestion import (
    discover_10x_mtx_triplets,
    ingest_10x_mtx_triplet,
    scan_fastq_plan,
)


FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

KB_REQUIRED_FIELDS = (
    "fastq_files",
    "index_path",
    "index_sha256",
    "t2g_path",
    "t2g_sha256",
    "chemistry",
    "output_dir",
)


class RawIngestionAgent(BaseAgent):
    """Convert supported raw scRNA inputs into canonical workspace h5ads."""

    name = "raw_ingestion_agent"
    description = "Deterministically ingests raw scRNA matrices/FASTQs."

    def __init__(self, memory: ARIAMemory, llm=None, api_key: str = None):
        super().__init__(memory, llm=llm, api_key=api_key)
        from aria.utils.environment_manager import env_manager
        self.env = env_manager

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
        requested_scrna_fastqs = self._scrna_fastq_inputs(exp_ctx)

        def _progress(stage: str, report: dict) -> None:
            if stage == "scan_10x_mtx":
                self.publish_status(
                    experiment_id,
                    "RawIngestionAgent: scanning 10X MEX inputs "
                    f"(files={report.get('files_seen', 0)}, "
                    f"candidates={report.get('candidate_dirs', 0)})...",
                    0.05,
                )
            elif stage == "validate_10x_mtx":
                self.publish_status(
                    experiment_id,
                    "RawIngestionAgent: validating 10X MEX inputs "
                    f"({report.get('validated_triplets', 0)}/"
                    f"{report.get('candidate_dirs', 0)})...",
                    0.10,
                )
            elif stage == "scan_fastq":
                self.publish_status(
                    experiment_id,
                    "RawIngestionAgent: scanning scRNA FASTQ inputs "
                    f"(files={report.get('files_seen', 0)}, "
                    f"FASTQ={report.get('fastq_count', 0)})...",
                    0.72,
                )

        if self._should_scan_scrna_mtx(exp_ctx):
            self.publish_status(
                experiment_id,
                "RawIngestionAgent: scanning confirmed scRNA MEX inputs...",
                0.03,
            )
            triplets = self._discover_10x_mtx_triplets(data_dir, _progress)
        else:
            triplets = []
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

        fastq_plan = (
            self._scan_fastq_plan(data_dir, _progress)
            if self._should_scan_scrna_fastq(exp_ctx) else
            {"status": "no_fastq_found", "mode": "fastq_kb_plan", "fastq_count": 0}
        )
        if fastq_plan.get("fastq_count", 0):
            records.append(fastq_plan)
            kb_params = (
                exp_ctx.get("raw_ingestion_kb")
                or context.get("raw_ingestion_kb")
                or {}
            )
            if (not kb_params.get("execute")
                    or self._missing_kb_fields(kb_params)):
                checkpoint = self._resolve_fastq_kb_checkpoint(
                    experiment_id,
                    fastq_plan,
                    kb_params,
                )
                records.append(checkpoint)
                if checkpoint.get("status") == "error":
                    errors.append({
                        "mode": checkpoint.get("mode", "fastq_kb_checkpoint"),
                        "error_type": checkpoint.get(
                            "reason", "fastq_kb_checkpoint_error"),
                        "details": checkpoint.get(
                            "details", "FASTQ kb checkpoint was cancelled."),
                    })
                kb_params = checkpoint.get("raw_ingestion_kb", {})
            if kb_params.get("execute"):
                kb_result = self._run_kb_count(kb_params)
                records.append(kb_result)
                if kb_result.get("status") == "success":
                    generated_h5ads.append(kb_result["output_h5ad"])
                else:
                    errors.append(self._kb_error_record(kb_result))

        if (requested_scrna_fastqs and not generated_h5ads
                and not self._scrna_canonical_inputs(exp_ctx)):
            errors.append({
                "mode": "fastq_kb_count",
                "error_type": "CanonicalH5adMissing",
                "details": (
                    "scRNA FASTQ inputs were detected, but no canonical .h5ad "
                    "was generated. ARIA will not dispatch scRNAAgent on raw "
                    "FASTQ files; provide complete raw_ingestion_kb parameters "
                    "or a precomputed .h5ad/MEX input."
                ),
                "fastq_files": requested_scrna_fastqs,
            })

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
            if rec.get("status") != "success" or not rec.get("output_h5ad"):
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

    def _run_kb_count(self, params: dict) -> dict:
        env = getattr(self, "env", None)
        if env is None:
            from aria.utils.environment_manager import env_manager
            env = env_manager
            self.env = env
        return env.run_in_stack(
            stack="ingestion",
            script_path="aria/scripts/rna_kb_count.py",
            params=params,
            timeout=int(params.get("timeout_seconds") or 86400),
        )

    @staticmethod
    def _discover_10x_mtx_triplets(data_dir: Path, progress_callback):
        try:
            return discover_10x_mtx_triplets(
                data_dir,
                progress_callback=progress_callback,
            )
        except TypeError:
            # Tests and older embeddings may monkeypatch the helper with the
            # pre-progress one-argument signature.
            return discover_10x_mtx_triplets(data_dir)

    @staticmethod
    def _scan_fastq_plan(data_dir: Path, progress_callback) -> dict:
        try:
            return scan_fastq_plan(
                data_dir,
                progress_callback=progress_callback,
            )
        except TypeError:
            return scan_fastq_plan(data_dir)

    def _resolve_fastq_kb_checkpoint(
        self,
        experiment_id: str,
        fastq_plan: dict,
        kb_params: dict,
    ) -> dict:
        question = self._format_fastq_kb_checkpoint(fastq_plan, kb_params)
        _, decision = self.publish_blocking_escalation(
            experiment_id=experiment_id,
            checkpoint="raw_ingestion.fastq_kb",
            question=question,
            options=[
                "Skip FASTQ quantification for now (Recommended)",
                "Use supplied raw_ingestion_kb parameters",
                "Other: paste raw_ingestion_kb JSON",
                "Cancel analysis",
            ],
            context={
                "raw_ingestion_fastq_plan": fastq_plan,
                "raw_ingestion_kb": kb_params,
            },
        )
        choice = self.checkpoint_choice_text(decision)
        if "cancel" in choice.lower():
            return {
                "status": "error",
                "mode": "fastq_kb_checkpoint",
                "reason": "user_cancelled_fastq_ingestion",
            }
        if choice.strip().startswith("{"):
            parsed = self._parse_kb_json_choice(choice)
            return parsed
        if "use supplied" in choice.lower():
            prepared = self._prepare_kb_params(kb_params)
            missing = self._missing_kb_fields(prepared)
            if missing:
                return {
                    "status": "blocked",
                    "mode": "fastq_kb_checkpoint",
                    "reason": "raw_ingestion_kb_incomplete",
                    "missing_fields": missing,
                    "raw_ingestion_kb": {},
                }
            return {
                "status": "success",
                "mode": "fastq_kb_checkpoint",
                "decision": "use_supplied_parameters",
                "raw_ingestion_kb": prepared,
            }
        return {
            "status": "skipped",
            "mode": "fastq_kb_checkpoint",
            "decision": "skip_fastq_quantification",
            "reason": "user_skipped_fastq_quantification",
            "raw_ingestion_kb": {},
        }

    def _parse_kb_json_choice(self, choice: str) -> dict:
        try:
            parsed = json.loads(choice)
        except json.JSONDecodeError as exc:
            return {
                "status": "blocked",
                "mode": "fastq_kb_checkpoint",
                "reason": "invalid_raw_ingestion_kb_json",
                "details": str(exc),
                "raw_ingestion_kb": {},
            }
        if not isinstance(parsed, dict):
            return {
                "status": "blocked",
                "mode": "fastq_kb_checkpoint",
                "reason": "raw_ingestion_kb_json_must_be_object",
                "raw_ingestion_kb": {},
            }
        prepared = self._prepare_kb_params(parsed)
        missing = self._missing_kb_fields(prepared)
        if missing:
            return {
                "status": "blocked",
                "mode": "fastq_kb_checkpoint",
                "reason": "raw_ingestion_kb_incomplete",
                "missing_fields": missing,
                "raw_ingestion_kb": {},
            }
        return {
            "status": "success",
            "mode": "fastq_kb_checkpoint",
            "decision": "use_json_parameters",
            "raw_ingestion_kb": prepared,
        }

    def _prepare_kb_params(self, params: dict) -> dict:
        prepared = dict(params or {})
        prepared["execute"] = True
        return prepared

    @staticmethod
    def _is_fastq_path(path: str | Path) -> bool:
        lower = str(path).lower()
        return lower.endswith(FASTQ_SUFFIXES)

    @classmethod
    def _scrna_fastq_inputs(cls, exp_ctx: dict) -> list[str]:
        modalities = exp_ctx.get("modalities", {}) or {}
        files = modalities.get("scRNA", []) or []
        return [str(path) for path in files if cls._is_fastq_path(path)]

    @classmethod
    def _should_scan_scrna_fastq(cls, exp_ctx: dict) -> bool:
        modalities = exp_ctx.get("modalities", {}) or {}
        if not modalities:
            return True
        return bool(cls._scrna_fastq_inputs(exp_ctx))

    @classmethod
    def _should_scan_scrna_mtx(cls, exp_ctx: dict) -> bool:
        modalities = exp_ctx.get("modalities", {}) or {}
        if not modalities:
            return True
        files = modalities.get("scRNA", []) or []
        for path in files:
            lower = str(path).lower()
            if lower.endswith((".h5ad", ".h5")) or cls._is_fastq_path(lower):
                continue
            if any(
                token in lower
                for token in (
                    "matrix.mtx",
                    "barcodes.tsv",
                    "features.tsv",
                    "genes.tsv",
                    "filtered_feature_bc_matrix",
                    "raw_feature_bc_matrix",
                )
            ):
                return True
            p = Path(path)
            if p.is_dir() and (
                (p / "matrix.mtx.gz").exists()
                or (p / "matrix.mtx").exists()
            ):
                return True
        return False

    @classmethod
    def _scrna_canonical_inputs(cls, exp_ctx: dict) -> list[str]:
        modalities = exp_ctx.get("modalities", {}) or {}
        files = modalities.get("scRNA", []) or []
        canonical_suffixes = (".h5ad", ".h5")
        return [
            str(path) for path in files
            if str(path).lower().endswith(canonical_suffixes)
        ]

    @staticmethod
    def _kb_error_record(result: dict) -> dict:
        blockers = result.get("blockers") or []
        details = (
            result.get("details")
            or "; ".join(str(blocker) for blocker in blockers)
            or "kb count did not produce a canonical .h5ad output."
        )
        return {
            "mode": result.get("mode", "fastq_kb_count"),
            "error_type": result.get("error_type", "KbCountFailed"),
            "details": details,
        }

    @staticmethod
    def _missing_kb_fields(params: dict) -> list[str]:
        missing = []
        for field in KB_REQUIRED_FIELDS:
            value = (params or {}).get(field)
            if value in (None, "", []):
                missing.append(field)
        return missing

    def _format_fastq_kb_checkpoint(self, fastq_plan: dict,
                                    kb_params: dict) -> str:
        samples = fastq_plan.get("samples", []) or []
        sample_lines = []
        for sample in samples[:8]:
            files = sample.get("files", {}) or {}
            roles = ", ".join(sorted(files)) or "no parsed R1/R2 roles"
            lanes = ", ".join(sample.get("lanes", []) or []) or "lane not parsed"
            sample_lines.append(
                f"- {sample.get('sample_id', 'sample')}: {roles}; {lanes}"
            )
        if len(samples) > 8:
            sample_lines.append(f"- ... {len(samples) - 8} more sample groups")
        missing = self._missing_kb_fields(kb_params)
        blocker_lines = [f"- {b}" for b in (fastq_plan.get("blockers", []) or [])]
        required = "\n".join(f"- {field}" for field in KB_REQUIRED_FIELDS)
        supplied = ", ".join(sorted(kb_params)) if kb_params else "none"
        return (
            "ARIA detected scRNA FASTQ inputs, but FASTQ quantification cannot "
            "run until all chemistry/reference inputs are explicit.\n\n"
            "Detected sample groups:\n"
            f"{chr(10).join(sample_lines) if sample_lines else '- none parsed'}\n\n"
            "Current blockers:\n"
            f"{chr(10).join(blocker_lines) if blocker_lines else '- metadata incomplete'}\n\n"
            "Required raw_ingestion_kb JSON fields:\n"
            f"{required}\n\n"
            f"Currently supplied fields: {supplied}.\n"
            f"Missing fields: {', '.join(missing) if missing else 'none'}.\n\n"
            "Use 'Other' to paste a JSON object with those fields, or skip "
            "FASTQ quantification for this run. ARIA will not guess chemistry, "
            "download references, or fabricate a matrix."
        )
