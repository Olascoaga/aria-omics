"""Tier 2.1 Typed Derived Query AST over Governed ARIA Output Tables (P2.1).

Provides deterministic, schema-checked relational queries (`filter`, `project`,
`rank`, `group_summary`, `validated_join`) over verified ARIA output tables for
benchmark adapters and post-run table inspection WITHOUT:
  - arbitrary `eval()` / `exec()` code execution,
  - new uncontracted statistical inference,
  - silent biological-unit conflation (e.g., joining cell-level rows onto donor-
    level contrasts without explicit projection), or
  - uncontrolled many-to-many (`m:m`) Cartesian row explosions.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

import pandas as pd

from aria.utils.provenance import hash_file, hash_params

FilterOperator = Literal["==", "!=", "<", "<=", ">", ">=", "in", "not_in", "is_not_null"]
JoinCardinality = Literal["1:1", "m:1", "1:m"]
AggFunc = Literal["count", "sum", "mean", "median", "min", "max"]

DERIVED_QUERY_SCHEMA = "aria.derived_query_ast.v1"


class DerivedQueryViolation(ValueError):
    """Raised when a DerivedQueryPlan violates schema, unit, or cardinality guards."""


@dataclass(frozen=True)
class TableContract:
    table_id: str
    biological_unit: str
    primary_key: tuple[str, ...]
    ledger_node_id: str
    source_path: str | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        if not self.table_id or not self.biological_unit or not self.ledger_node_id:
            raise DerivedQueryViolation("TableContract requires table_id, biological_unit, and ledger_node_id")
        if not self.ledger_node_id.startswith("ledger://"):
            raise DerivedQueryViolation(f"Invalid ledger_node_id {self.ledger_node_id!r}: must start with 'ledger://'")
        if not self.primary_key:
            raise DerivedQueryViolation("TableContract requires a non-empty primary_key tuple")


@dataclass(frozen=True)
class FilterOp:
    column: str
    op: FilterOperator
    value: Any = None


@dataclass(frozen=True)
class ProjectOp:
    columns: tuple[str, ...]


@dataclass(frozen=True)
class RankOp:
    by: tuple[str, ...]
    ascending: tuple[bool, ...]
    limit: int | None = None


@dataclass(frozen=True)
class GroupSummaryOp:
    by: tuple[str, ...]
    aggregations: dict[str, AggFunc]


@dataclass(frozen=True)
class ValidatedJoinOp:
    right_table_id: str
    on: tuple[str, ...]
    how: Literal["inner", "left"] = "inner"
    expected_cardinality: JoinCardinality = "1:1"
    allow_unit_projection: bool = False


QueryStep = FilterOp | ProjectOp | RankOp | GroupSummaryOp | ValidatedJoinOp


@dataclass(frozen=True)
class DerivedQueryPlan:
    base_table_id: str
    steps: tuple[QueryStep, ...] = field(default_factory=tuple)


def _validate_primary_key_unique(df: pd.DataFrame, pk: tuple[str, ...], table_id: str) -> None:
    missing = [c for c in pk if c not in df.columns]
    if missing:
        raise DerivedQueryViolation(f"Table '{table_id}' missing primary_key column(s): {missing}")
    if df.duplicated(subset=list(pk)).any():
        raise DerivedQueryViolation(f"Table '{table_id}' has duplicate keys on primary_key {pk}")


def _verify_join_cardinality(
    left: pd.DataFrame,
    right: pd.DataFrame,
    on: list[str],
    expected: JoinCardinality,
) -> None:
    if expected not in {"1:1", "m:1", "1:m"}:
        raise DerivedQueryViolation(f"Many-to-many or unknown join cardinality {expected!r} is forbidden")
    left_dup = bool(left.duplicated(subset=on).any())
    right_dup = bool(right.duplicated(subset=on).any())
    if expected == "1:1" and (left_dup or right_dup):
        raise DerivedQueryViolation(f"Join on {on} violated 1:1 cardinality (left_dup={left_dup}, right_dup={right_dup})")
    if expected == "m:1" and right_dup:
        raise DerivedQueryViolation(f"Join on {on} violated m:1 cardinality (right keys are not unique)")
    if expected == "1:m" and left_dup:
        raise DerivedQueryViolation(f"Join on {on} violated 1:m cardinality (left keys are not unique)")


def execute_derived_query(
    plan: DerivedQueryPlan,
    tables: dict[str, tuple[TableContract, pd.DataFrame]],
) -> dict[str, Any]:
    """Execute a validated `DerivedQueryPlan` and return a provenance-bound receipt."""
    if plan.base_table_id not in tables:
        raise DerivedQueryViolation(f"Unknown base_table_id: {plan.base_table_id!r}")

    base_contract, base_df = tables[plan.base_table_id]
    _validate_primary_key_unique(base_df, base_contract.primary_key, base_contract.table_id)

    cur = base_df.copy()
    input_denominator = int(len(cur))
    active_unit = base_contract.biological_unit
    active_pk = list(base_contract.primary_key)
    ledger_nodes = [base_contract.ledger_node_id]
    table_digests = {
        base_contract.table_id: (
            base_contract.source_sha256
            or (hash_file(base_contract.source_path) if base_contract.source_path else hash_params({"rows": len(base_df), "cols": list(base_df.columns)}))
        )
    }

    serialized_steps: list[dict[str, Any]] = []
    for step in plan.steps:
        if isinstance(step, FilterOp):
            if step.column not in cur.columns:
                raise DerivedQueryViolation(f"Filter column {step.column!r} not in {list(cur.columns)}")
            col = cur[step.column]
            if step.op == "==":
                mask = col == step.value
            elif step.op == "!=":
                mask = col != step.value
            elif step.op == "<":
                mask = col < step.value
            elif step.op == "<=":
                mask = col <= step.value
            elif step.op == ">":
                mask = col > step.value
            elif step.op == ">=":
                mask = col >= step.value
            elif step.op == "in":
                mask = col.isin(list(step.value or []))
            elif step.op == "not_in":
                mask = ~col.isin(list(step.value or []))
            elif step.op == "is_not_null":
                mask = col.notna()
            else:
                raise DerivedQueryViolation(f"Unsupported FilterOp operator: {step.op!r}")
            cur = cur.loc[mask].copy()
            serialized_steps.append({"op": "filter", **asdict(step)})

        elif isinstance(step, ValidatedJoinOp):
            if step.right_table_id not in tables:
                raise DerivedQueryViolation(f"Unknown right_table_id: {step.right_table_id!r}")
            r_contract, r_df = tables[step.right_table_id]
            _validate_primary_key_unique(r_df, r_contract.primary_key, r_contract.table_id)
            if (
                r_contract.biological_unit != active_unit
                and not step.allow_unit_projection
            ):
                raise DerivedQueryViolation(
                    f"Refusing silent biological-unit join ({active_unit!r} x {r_contract.biological_unit!r}); "
                    "set allow_unit_projection=True only when cross-unit key projection is valid."
                )
            join_cols = list(step.on)
            for c in join_cols:
                if c not in cur.columns or c not in r_df.columns:
                    raise DerivedQueryViolation(f"Join key {c!r} missing from left or right table")
            _verify_join_cardinality(cur, r_df, join_cols, step.expected_cardinality)
            cur = cur.merge(r_df, on=join_cols, how=step.how, suffixes=("", f"__{r_contract.table_id}"))
            if r_contract.biological_unit != active_unit:
                active_unit = f"{active_unit}_x_{r_contract.biological_unit}"
            if r_contract.ledger_node_id not in ledger_nodes:
                ledger_nodes.append(r_contract.ledger_node_id)
            table_digests[r_contract.table_id] = (
                r_contract.source_sha256
                or (hash_file(r_contract.source_path) if r_contract.source_path else hash_params({"rows": len(r_df), "cols": list(r_df.columns)}))
            )
            serialized_steps.append({"op": "validated_join", **asdict(step)})

        elif isinstance(step, GroupSummaryOp):
            by_cols = list(step.by)
            if not by_cols or any(c not in cur.columns for c in by_cols):
                raise DerivedQueryViolation(f"Invalid GroupSummaryOp columns: {by_cols}")
            allowed_aggs = {"count", "sum", "mean", "median", "min", "max"}
            for target_col, agg_fn in step.aggregations.items():
                if target_col not in cur.columns:
                    raise DerivedQueryViolation(f"Aggregation column {target_col!r} missing")
                if agg_fn not in allowed_aggs:
                    raise DerivedQueryViolation(f"Disallowed aggregation {agg_fn!r}")
            grouped = cur.groupby(by_cols, dropna=False, sort=True).agg(step.aggregations).reset_index()
            cur = grouped
            active_unit = f"grouped:{','.join(by_cols)}"
            active_pk = by_cols
            serialized_steps.append({"op": "group_summary", **asdict(step)})

        elif isinstance(step, RankOp):
            by_cols = list(step.by)
            asc = list(step.ascending)
            if not by_cols or len(by_cols) != len(asc):
                raise DerivedQueryViolation("RankOp requires matching `by` and `ascending` tuples")
            if any(c not in cur.columns for c in by_cols):
                raise DerivedQueryViolation(f"RankOp column missing in {list(cur.columns)}")
            tie_breakers = [c for c in active_pk if c in cur.columns and c not in by_cols]
            full_by = by_cols + tie_breakers
            full_asc = asc + [True] * len(tie_breakers)
            cur = cur.sort_values(full_by, ascending=full_asc, kind="mergesort")
            if step.limit is not None:
                if int(step.limit) < 0:
                    raise DerivedQueryViolation("RankOp limit must be non-negative")
                cur = cur.head(int(step.limit)).copy()
            serialized_steps.append({"op": "rank", **asdict(step)})

        elif isinstance(step, ProjectOp):
            proj_cols = list(step.columns)
            if not proj_cols or any(c not in cur.columns for c in proj_cols):
                raise DerivedQueryViolation(f"ProjectOp column missing from {list(cur.columns)}")
            cur = cur.loc[:, proj_cols].copy()
            active_pk = [c for c in active_pk if c in proj_cols]
            serialized_steps.append({"op": "project", **asdict(step)})

        else:
            raise DerivedQueryViolation(f"Unsupported AST step type: {type(step).__name__}")

    ast_sha256 = hash_params({
        "schema": DERIVED_QUERY_SCHEMA,
        "base_table_id": plan.base_table_id,
        "steps": serialized_steps,
    })
    return {
        "schema_version": DERIVED_QUERY_SCHEMA,
        "ast_sha256": ast_sha256,
        "source_tables_sha256": table_digests,
        "ledger_node_ids": ledger_nodes,
        "biological_unit": active_unit,
        "input_denominator": input_denominator,
        "output_row_count": int(len(cur)),
        "columns": list(cur.columns),
        "rows": cur.to_dict(orient="records"),
    }
