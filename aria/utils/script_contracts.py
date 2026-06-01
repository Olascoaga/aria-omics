"""Typed IPC contracts for scripts executed by EnvironmentManager."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field


CONTRACT_VERSION = "1.0"

FieldType = Literal[
    "any",
    "str",
    "int",
    "float",
    "bool",
    "dict",
    "list",
    "path",
    "list[path]",
]


class ContractField(BaseModel):
    name: str
    type: FieldType = "any"
    required: bool = True
    allow_empty: bool = True
    path_must_exist: bool = False
    # Backward-compatible alternative keys. A payload satisfies the field if it
    # carries `name` OR any alias (the canonical `name` takes precedence). Used
    # to accept legacy param spellings without loosening the contract (P0-1).
    aliases: tuple[str, ...] = ()


class ContractIssue(BaseModel):
    field: str
    message: str
    severity: Literal["error", "warning"] = "error"


class ScriptContract(BaseModel):
    script_path: str
    contract_version: str = CONTRACT_VERSION
    validation_level: Literal["production", "beta", "scaffold"] = "production"
    inputs: tuple[ContractField, ...] = Field(default_factory=tuple)
    success_outputs: tuple[ContractField, ...] = Field(default_factory=tuple)

    def validate_params(self, params: dict[str, Any]) -> list[ContractIssue]:
        issues = self._version_issues(params)
        issues.extend(_validate_fields(params, self.inputs))
        return issues

    def validate_result(self, result: dict[str, Any]) -> list[ContractIssue]:
        issues = self._version_issues(result)
        if "status" not in result:
            issues.append(ContractIssue(
                field="status",
                message="Script output is missing required field 'status'.",
            ))
            return issues
        if result.get("status") == "error":
            issues.extend(_validate_fields(result, (
                ContractField(name="error_type", type="str", allow_empty=False),
                ContractField(name="details", type="str", allow_empty=False),
            )))
            return issues
        if result.get("status") != "success":
            issues.append(ContractIssue(
                field="status",
                message=(
                    "Script output status must be 'success' or 'error'; "
                    f"got {result.get('status')!r}."
                ),
            ))
            return issues
        issues.extend(_validate_fields(result, self.success_outputs))
        return issues

    def _version_issues(self, payload: dict[str, Any]) -> list[ContractIssue]:
        requested = (
            payload.get("_aria_contract_version")
            or payload.get("contract_version")
            or payload.get("ipc_contract_version")
        )
        if requested and str(requested) != self.contract_version:
            return [ContractIssue(
                field="_aria_contract_version",
                message=(
                    "Incompatible IPC contract version for "
                    f"{self.script_path}: expected {self.contract_version}, "
                    f"got {requested}."
                ),
            )]
        return []


def contract_for_script(script_path: str | Path) -> ScriptContract | None:
    key = _canonical_script_key(script_path)
    return SCRIPT_CONTRACTS.get(key)


def _canonical_script_key(script_path: str | Path) -> str:
    path = Path(script_path)
    parts = path.parts
    if "aria" in parts:
        idx = parts.index("aria")
        return str(Path(*parts[idx:]))
    return str(path).replace("\\", "/")


def _validate_fields(payload: dict[str, Any],
                     fields: tuple[ContractField, ...]) -> list[ContractIssue]:
    issues: list[ContractIssue] = []
    for field in fields:
        # Resolve the present key: canonical name first, then any alias.
        key = field.name if field.name in payload else next(
            (a for a in field.aliases if a in payload), None)
        if key is None:
            if field.required:
                issues.append(ContractIssue(
                    field=field.name,
                    message=f"Missing required field '{field.name}'.",
                ))
            continue
        value = payload.get(key)
        if _is_empty(value) and not field.allow_empty:
            issues.append(ContractIssue(
                field=field.name,
                message=f"Field '{field.name}' must not be empty.",
            ))
            continue
        if _is_empty(value) and field.allow_empty:
            continue
        type_issue = _type_issue(field, value)
        if type_issue:
            issues.append(type_issue)
            continue
        if field.path_must_exist:
            issues.extend(_path_issues(field, value))
    return issues


def _is_empty(value: Any) -> bool:
    return value is None or value == "" or value == []


def _type_issue(field: ContractField, value: Any) -> ContractIssue | None:
    ok = True
    if field.type == "str":
        ok = isinstance(value, str)
    elif field.type == "int":
        ok = isinstance(value, int) and not isinstance(value, bool)
    elif field.type == "float":
        ok = isinstance(value, (int, float)) and not isinstance(value, bool)
    elif field.type == "bool":
        ok = isinstance(value, bool)
    elif field.type == "dict":
        ok = isinstance(value, dict)
    elif field.type == "list":
        ok = isinstance(value, list)
    elif field.type == "path":
        ok = isinstance(value, str)
    elif field.type == "list[path]":
        ok = isinstance(value, list) and all(isinstance(v, str) for v in value)
    if ok:
        return None
    return ContractIssue(
        field=field.name,
        message=f"Field '{field.name}' must be {field.type}; got {type(value).__name__}.",
    )


def _path_issues(field: ContractField, value: Any) -> list[ContractIssue]:
    paths = value if field.type == "list[path]" else [value]
    issues = []
    for item in paths:
        if item in (None, ""):
            if not field.allow_empty:
                issues.append(ContractIssue(
                    field=field.name,
                    message=f"Field '{field.name}' contains an empty path.",
                ))
            continue
        if not Path(str(item)).expanduser().exists():
            issues.append(ContractIssue(
                field=field.name,
                message=f"Required path for '{field.name}' does not exist: {item}",
            ))
    return issues


def _f(name: str, type_: FieldType = "any", *, required: bool = True,
       allow_empty: bool = True, path_must_exist: bool = False,
       aliases: tuple[str, ...] = ()) -> ContractField:
    return ContractField(
        name=name,
        type=type_,
        required=required,
        allow_empty=allow_empty,
        path_must_exist=path_must_exist,
        aliases=aliases,
    )


SCRIPT_CONTRACTS: dict[str, ScriptContract] = {
    "aria/scripts/rna_qc.py": ScriptContract(
        script_path="aria/scripts/rna_qc.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("n_cells_before", "int"),
            _f("n_cells_after", "int"),
        ),
    ),
    "aria/scripts/rna_clustering.py": ScriptContract(
        script_path="aria/scripts/rna_clustering.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("output_path", "str", required=False),
        ),
    ),
    "aria/scripts/rna_bulk_de.py": ScriptContract(
        script_path="aria/scripts/rna_bulk_de.py",
        validation_level="production",
        inputs=(
            _f("files", "list[path]", allow_empty=False, path_must_exist=True),
            _f("design_factor", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("n_contrasts", "int"),
            _f("contrasts", "list"),
            _f("count_source", "str", required=False),
        ),
    ),
    "aria/scripts/rna_diff_abundance.py": ScriptContract(
        script_path="aria/scripts/rna_diff_abundance.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
            _f("groupby", "str", allow_empty=False),
            _f("condition_col", "str", allow_empty=False),
            _f("replicate_col", "str", allow_empty=False),
            _f("comparisons", "list", allow_empty=False),
        ),
        success_outputs=(
            _f("per_comparison", "dict"),
        ),
    ),
    "aria/scripts/rna_pseudobulk_de.py": ScriptContract(
        script_path="aria/scripts/rna_pseudobulk_de.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
            _f("groupby", "str", allow_empty=False),
            _f("condition_col", "str", allow_empty=False),
            _f("replicate_col", "str", allow_empty=False),
            _f("comparisons", "list", allow_empty=False),
        ),
        success_outputs=(
            _f("per_group", "dict"),
            _f("n_groups", "int"),
        ),
    ),
    "aria/scripts/rna_pathway_per_cluster.py": ScriptContract(
        script_path="aria/scripts/rna_pathway_per_cluster.py",
        validation_level="production",
        inputs=(
            _f("pseudobulk_result", "dict", required=False),
        ),
        success_outputs=(
            _f("per_cluster", "dict"),
        ),
    ),
    "aria/scripts/rna_cellcomm.py": ScriptContract(
        script_path="aria/scripts/rna_cellcomm.py",
        validation_level="beta",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
            _f("groupby", "str", allow_empty=False, aliases=("cell_type_col",)),
        ),
        success_outputs=(
            _f("n_interactions", "int", required=False),
        ),
    ),
    "aria/scripts/rna_trajectory.py": ScriptContract(
        script_path="aria/scripts/rna_trajectory.py",
        validation_level="beta",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
            _f("groupby", "str", required=False),
        ),
        success_outputs=(
            _f("paga", "dict", required=False),
        ),
    ),
    "aria/scripts/chromatin_qc.py": ScriptContract(
        script_path="aria/scripts/chromatin_qc.py",
        validation_level="scaffold",
        inputs=(
            _f("data_type", "str", allow_empty=False),
            _f("files", "list[path]", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("data_type", "str", required=False),
        ),
    ),
    "aria/scripts/chromatin_peaks.py": ScriptContract(
        script_path="aria/scripts/chromatin_peaks.py",
        validation_level="scaffold",
        inputs=(
            _f("data_type", "str", allow_empty=False),
            _f("files", "list[path]", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("n_peaks", "int", required=False),
        ),
    ),

    # ── P1-12: contracts for the remaining dispatchable scripts ──────────────
    # scRNA helper scripts (validated path).
    "aria/scripts/rna_apply_cluster_labels.py": ScriptContract(
        script_path="aria/scripts/rna_apply_cluster_labels.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("label_col", "str", required=False),
        ),
    ),
    "aria/scripts/rna_celltypist.py": ScriptContract(
        script_path="aria/scripts/rna_celltypist.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("output_path", "str", required=False),
        ),
    ),
    "aria/scripts/rna_concat.py": ScriptContract(
        script_path="aria/scripts/rna_concat.py",
        validation_level="production",
        inputs=(
            _f("samples", "list", required=False),
            _f("output_dir", "str", required=False),
        ),
        success_outputs=(
            _f("output_path", "str", required=False),
        ),
    ),
    "aria/scripts/rna_de_per_cluster.py": ScriptContract(
        script_path="aria/scripts/rna_de_per_cluster.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("de_genes_by_cluster", "dict", required=False),
        ),
    ),
    "aria/scripts/rna_inject_condition.py": ScriptContract(
        script_path="aria/scripts/rna_inject_condition.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
            _f("output_path", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("condition_col", "str", required=False),
        ),
    ),
    "aria/scripts/rna_integration.py": ScriptContract(
        script_path="aria/scripts/rna_integration.py",
        validation_level="production",
        inputs=(
            _f("data_path", "path", allow_empty=False, path_must_exist=True),
        ),
        success_outputs=(
            _f("output_path", "str", required=False),
        ),
    ),
    "aria/scripts/rna_figure_umap.py": ScriptContract(
        script_path="aria/scripts/rna_figure_umap.py",
        validation_level="production",
        inputs=(
            _f("output_dir", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("embedding_key", "str", required=False),
        ),
    ),
    "aria/scripts/rna_figure_paga.py": ScriptContract(
        script_path="aria/scripts/rna_figure_paga.py",
        validation_level="production",
        inputs=(
            _f("output_dir", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("figures", "dict", required=False),
        ),
    ),
    # FASTQ preprocessing flow (beta).
    "aria/scripts/rna_fastq_qc.py": ScriptContract(
        script_path="aria/scripts/rna_fastq_qc.py",
        validation_level="beta",
        inputs=(
            _f("fastq_dir", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("n_reads_raw", "int", required=False),
        ),
    ),
    "aria/scripts/rna_align.py": ScriptContract(
        script_path="aria/scripts/rna_align.py",
        validation_level="beta",
        inputs=(
            _f("samples", "list", allow_empty=False),
            _f("genome_dir", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("bam_files", "list", required=False),
        ),
    ),
    "aria/scripts/rna_quantify.py": ScriptContract(
        script_path="aria/scripts/rna_quantify.py",
        validation_level="beta",
        inputs=(
            _f("bam_files", "list", allow_empty=False),
            _f("gtf_file", "str", allow_empty=False),
        ),
        success_outputs=(
            _f("counts_matrix", "any", required=False),
        ),
    ),
    # Scaffold modalities (dispatch-gated; contract present for completeness).
    "aria/scripts/integration_wnn.py": ScriptContract(
        script_path="aria/scripts/integration_wnn.py",
        validation_level="scaffold",
        inputs=(
            _f("rna_files", "list", required=False),
            _f("atac_files", "list", required=False),
        ),
        success_outputs=(
            _f("n_cells", "int", required=False),
        ),
    ),
    "aria/scripts/integration_mofa.py": ScriptContract(
        script_path="aria/scripts/integration_mofa.py",
        validation_level="scaffold",
        inputs=(
            _f("output_dir", "str", required=False),
        ),
        success_outputs=(
            _f("factor_id", "any", required=False),
        ),
    ),
    "aria/scripts/integration_peak2gene.py": ScriptContract(
        script_path="aria/scripts/integration_peak2gene.py",
        validation_level="scaffold",
        inputs=(
            _f("rna_files", "list", required=False),
            _f("atac_files", "list", required=False),
        ),
        success_outputs=(
            _f("output_path", "str", required=False),
        ),
    ),
    "aria/scripts/hic_inspect.py": ScriptContract(
        script_path="aria/scripts/hic_inspect.py",
        validation_level="scaffold",
        inputs=(
            _f("files", "list", required=False),
        ),
        success_outputs=(
            _f("chromosomes", "list", required=False),
        ),
    ),
    "aria/scripts/hic_qc_and_balance.py": ScriptContract(
        script_path="aria/scripts/hic_qc_and_balance.py",
        validation_level="scaffold",
        inputs=(
            _f("files", "list", required=False),
        ),
        success_outputs=(
            _f("cis_trans_ratio", "float", required=False),
        ),
    ),
    "aria/scripts/hic_topology.py": ScriptContract(
        script_path="aria/scripts/hic_topology.py",
        validation_level="scaffold",
        inputs=(
            _f("files", "list", required=False),
        ),
        success_outputs=(
            _f("algorithm", "str", required=False),
        ),
    ),
}
