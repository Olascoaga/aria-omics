"""V48 Slice I1 — paired `.h5mu` contract, same-cell gate, joint QC, `/obs` join.

Layers (each gated on what the active environment can import):

1. Same-cell hard gate (`enforce_same_cell_integration_gate`) — pure Python,
   runs in any env. Agent / IPC-contract / registry checks additionally need
   the orchestration deps (`pydantic`, `litellm`) and skip without them
   (e.g. inside `aria-chromatin-env`).
2. HDF5 `/obs` metadata join (`read_h5ad_obs_columns`,
   `join_obs_labels_from_h5ad`) — needs `h5py` only; reference `.h5ad` files
   are hand-built with `h5py` so `anndata` is not required.
3. Synthetic paired `.h5mu` (`read_h5mu_paired`, `integration_paired_qc`) —
   needs `mudata` (chromatin stack).
4. Real-data validation on a paired multiome `.h5mu` — dataset-gated.

Synthetic fixtures use neutral names (`CellTypeA`, `GENE1`, `MT-CO1`) per
ADR-011.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from aria.utils.multiome_contracts import (
    enforce_same_cell_integration_gate,
    validate_multiome_contract,
)


GATE_REASON = "multiome_same_cell_gate_blocked"


def _checks(items: list[dict]) -> set[str]:
    return {str(i.get("check")) for i in items}


# ════════════════════════════════════════════════════════════════════════════
# 1. Same-cell hard gate (any env)
# ════════════════════════════════════════════════════════════════════════════

SPLIT_CTX = {
    "genome": "GRCh38",
    "modalities": {
        "scRNA": ["/data/sample1_rna.h5ad"],
        "scATAC": ["/data/sample1_atac.h5ad"],
    },
}
PAIRED_CTX = {"modalities": {"scATAC": ["/data/sample1_paired.h5mu"]}}


def test_gate_blocks_split_modalities_but_lane_validation_stays_yellow():
    validation = validate_multiome_contract(SPLIT_CTX)
    # Independent single-assay lanes must stay dispatchable (warning, not red).
    assert validation["status"] == "yellow"
    assert "multiome_contract_pairing_unconfirmed" in _checks(validation["findings"])

    gate = enforce_same_cell_integration_gate(SPLIT_CTX)
    assert gate["allowed"] is False
    assert gate["status"] == "blocked"
    assert gate["reason"] == GATE_REASON
    assert "multiome_contract_pairing_unconfirmed" in _checks(gate["blockers"])
    assert all(b["severity"] == "blocking" for b in gate["blockers"])


@pytest.mark.parametrize(
    "ctx",
    [
        {"modalities": {"scRNA": ["/data/sample1_rna.h5ad"]}},
        {"modalities": {"scATAC": ["/data/sample1_atac.h5ad"]}},
        {"modalities": {"bulk_RNA": ["/data/counts.tsv"],
                        "bulk_ATAC": ["/data/peaks.bed"]}},
        {"modalities": {}},
    ],
    ids=["scRNA_only", "scATAC_only", "bulk_rna_atac", "empty"],
)
def test_gate_blocks_not_applicable_contexts(ctx):
    assert validate_multiome_contract(ctx)["status"] == "not_applicable"
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert gate["reason"] == GATE_REASON
    assert gate["validation_status"] == "not_applicable"
    assert "multiome_same_cell_required" in _checks(gate["blockers"])


def test_gate_blocks_sample_mismatch():
    ctx = {
        "genome": "GRCh38",
        "multiome_contract": {
            "object_type": "paired_mudata",
            "same_cell": True,
            "cell_namespace": "shared_mudata_obs",
            "rna_sample_ids": ["sample1"],
            "atac_sample_ids": ["sample2"],
            "evidence_paths": ["/data/sample1_paired.h5mu"],
        },
    }
    assert validate_multiome_contract(ctx)["status"] == "red"
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert gate["reason"] == GATE_REASON
    assert "multiome_contract_sample_mismatch" in _checks(gate["blockers"])


@pytest.mark.parametrize("same_cell", [False, None, "true", 1])
def test_gate_requires_literal_same_cell_true(same_cell):
    ctx = {
        "genome": "GRCh38",
        "multiome_contract": {
            "object_type": "paired_mudata",
            "same_cell": same_cell,
            "cell_namespace": "shared_mudata_obs",
            "evidence_paths": ["/data/sample1_paired.h5mu"],
        },
    }
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert gate["reason"] == GATE_REASON


def test_gate_blocks_same_cell_without_barcode_namespace():
    ctx = {
        "genome": "GRCh38",
        "multiome_contract": {
            "object_type": "explicit_pairing",
            "same_cell": True,
        },
    }
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert "multiome_contract_barcode_namespace_missing" in _checks(gate["blockers"])


def test_gate_blocks_paired_mudata_without_backing_path():
    ctx = {
        "genome": "GRCh38",
        "multiome_contract": {
            "object_type": "paired_mudata",
            "same_cell": True,
            "cell_namespace": "shared_mudata_obs",
        },
    }
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is False
    assert "multiome_paired_path_missing" in _checks(gate["blockers"])


def test_gate_allows_paired_h5mu_genome_warning_by_default():
    gate = enforce_same_cell_integration_gate(PAIRED_CTX)
    assert gate["allowed"] is True
    assert gate["status"] == "allowed"
    assert gate["reason"] is None
    assert gate["blockers"] == []
    assert gate["contract"]["object_type"] == "paired_mudata"
    assert gate["contract"]["same_cell"] is True
    assert gate["evidence_paths"] == ["/data/sample1_paired.h5mu"]
    assert "multiome_contract_genome_feature_space_unconfirmed" in _checks(
        gate["warnings"]
    )


def test_gate_require_genome_promotes_missing_genome_to_blocker():
    gate = enforce_same_cell_integration_gate(PAIRED_CTX, require_genome=True)
    assert gate["allowed"] is False
    assert gate["reason"] == GATE_REASON
    assert "multiome_contract_genome_feature_space_unconfirmed" in _checks(
        gate["blockers"]
    )

    with_genome = {**PAIRED_CTX, "genome": "GRCh38"}
    gate_ok = enforce_same_cell_integration_gate(with_genome, require_genome=True)
    assert gate_ok["allowed"] is True
    assert gate_ok["warnings"] == []


def test_gate_allows_paired_h5mu_from_assay_detection():
    ctx = {
        "genome": "GRCh38",
        "modalities": {"scATAC": ["/data/sample1_paired.h5mu"]},
        "assay_detections": [{
            "path": "/data/sample1_paired.h5mu",
            "modality": "scATAC",
            "confidence": "high",
            "evidence": {"format": "h5mu", "modalities": ["rna", "atac"],
                         "paired_rna": True},
        }],
    }
    gate = enforce_same_cell_integration_gate(ctx, require_genome=True)
    assert gate["allowed"] is True
    assert gate["evidence_paths"] == ["/data/sample1_paired.h5mu"]


def test_gate_h5mu_listed_under_scrna_and_scatac_is_deduplicated():
    ctx = {"genome": "GRCh38", "modalities": {
        "scRNA": ["/data/sample1_paired.h5mu"],
        "scATAC": ["/data/sample1_paired.h5mu"],
    }}
    gate = enforce_same_cell_integration_gate(ctx)
    assert gate["allowed"] is True
    assert gate["evidence_paths"] == ["/data/sample1_paired.h5mu"]


# ── Agent / IPC contract / registry (need orchestration deps) ──────────────

def _require_orchestration_deps():
    pytest.importorskip("pydantic")
    pytest.importorskip("litellm")


class _RecordingEnv:
    def __init__(self, log: list, paired_result: dict):
        self.log = log
        self.calls: list[dict] = []
        self.paired_result = paired_result

    def run_in_stack(self, stack, script_path, params, timeout=None):
        self.calls.append({"stack": stack, "script_path": script_path,
                           "params": dict(params)})
        self.log.append(("env", Path(script_path).name))
        if Path(script_path).name == "integration_paired_qc.py":
            return dict(self.paired_result)
        return {"status": "success"}


def _make_agent(paired_result: dict | None = None):
    from aria.agents.integration_agent import IntegrationAgent

    log: list = []
    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.env = _RecordingEnv(log, paired_result or {
        "status": "success",
        "output_path": "/data/aria_integration/paired_qc_sample1.h5mu",
        "n_cells_before": 100, "n_cells_after": 90,
        "n_genes": 50, "n_peaks": 60, "label_join": {"ran": False},
    })
    agent._log = log
    agent._findings = []
    agent.publish_status = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None
    agent.publish_finding = lambda exp_id, payload, conf: agent._findings.append(payload)

    def _wnn(experiment_id, exp_ctx, intent, mods):
        log.append(("wnn", dict(mods)))
        return {"status": "success"}

    def _p2g(experiment_id, exp_ctx, intent, mods, rna_res, chrom_res):
        log.append(("peak_to_gene", dict(mods)))
        return {"status": "success"}

    def _mofa(experiment_id, exp_ctx, intent, mods):
        log.append(("mofa", dict(mods)))
        return {"status": "success"}

    agent._run_wnn = _wnn
    agent._run_peak_to_gene = _p2g
    agent._run_mofa = _mofa
    agent._resolve_conflicts = lambda *a, **k: []
    return agent


def test_agent_run_blocks_split_modalities_before_any_dispatch():
    _require_orchestration_deps()
    agent = _make_agent()
    out = agent.run("exp_split", {
        "exp_context": SPLIT_CTX,
        "biological_intent": {"user_question": "cell states"},
    })
    assert out["status"] == "blocked"
    assert out["reason"] == GATE_REASON
    assert out["multiome_gate"]["allowed"] is False
    assert agent.env.calls == []
    assert agent._log == []


def test_agent_run_dispatches_paired_qc_on_chromatin_stack_before_wnn():
    _require_orchestration_deps()
    agent = _make_agent()
    ctx = {**PAIRED_CTX, "genome": "GRCh38", "sample_id": "sample1"}
    out = agent.run("exp_paired", {
        "exp_context": ctx,
        "biological_intent": {"user_question": "cell states"},
    })
    assert out["status"] == "done"
    assert out["strategy"] == "wnn"
    assert out["findings"]["multiome_gate"]["allowed"] is True
    assert out["findings"]["paired_qc"]["status"] == "success"

    first = agent.env.calls[0]
    assert first["stack"] == "chromatin"
    assert first["script_path"] == "aria/scripts/integration_paired_qc.py"
    assert first["params"]["h5mu_path"] == "/data/sample1_paired.h5mu"
    assert first["params"]["genome"] == "GRCh38"
    assert first["params"]["sample_id"] == "sample1"

    order = [entry[0] if entry[0] != "env" else entry[1] for entry in agent._log]
    assert order.index("integration_paired_qc.py") < order.index("wnn")
    wnn_mods = next(e[1] for e in agent._log if e[0] == "wnn")
    validated = "/data/aria_integration/paired_qc_sample1.h5mu"
    assert wnn_mods["scRNA"] == [validated]
    assert wnn_mods["scATAC"] == [validated]


def test_agent_run_stops_when_paired_qc_fails():
    _require_orchestration_deps()
    agent = _make_agent({"status": "error", "error_type": "NoSharedBarcodes",
                         "details": "zero shared barcodes"})
    out = agent.run("exp_paired_fail", {
        "exp_context": {**PAIRED_CTX, "genome": "GRCh38"},
        "biological_intent": {"user_question": "cell states"},
    })
    assert out["status"] == "failed"
    assert out["reason"] == "paired_qc_failed"
    assert [e for e in agent._log if e[0] in ("wnn", "peak_to_gene", "mofa")] == []


def test_registry_integrity_and_integration_tier_still_scaffold():
    _require_orchestration_deps()
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION
    from aria.utils.registry_integrity import (
        AGENT_SCRIPT_CONTRACTS,
        check_registry_integrity,
    )

    issues = check_registry_integrity()
    assert [i for i in issues if i.severity == "error"] == [], issues
    assert "aria/scripts/integration_paired_qc.py" in (
        AGENT_SCRIPT_CONTRACTS["integration_agent"]["required"]
    )
    assert INTEGRATION_VALIDATION["level"] == "production"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is False
    assert INTEGRATION_VALIDATION["production_promoted"] is True


def test_ipc_contract_for_paired_qc(tmp_path):
    pytest.importorskip("pydantic")
    from aria.utils.script_contracts import contract_for_script

    contract = contract_for_script("aria/scripts/integration_paired_qc.py")
    assert contract is not None
    assert contract.validation_level == "production"

    h5mu = tmp_path / "sample1_paired.h5mu"
    h5mu.write_bytes(b"")
    assert contract.validate_params({"h5mu_path": str(h5mu)}) == []
    assert contract.validate_params({"data_path": str(h5mu)}) == []
    assert contract.validate_params({}) != []
    assert contract.validate_params({"h5mu_path": ""}) != []

    out = tmp_path / "paired_qc_sample1.h5mu"
    out.write_bytes(b"")
    success = {
        "status": "success", "output_path": str(out),
        "n_cells_before": 10, "n_cells_after": 8, "n_genes": 5, "n_peaks": 6,
        "label_join": {"ran": False}, "metrics_not_computed": ["frip"],
    }
    assert contract.validate_result(success) == []
    missing_peaks = {k: v for k, v in success.items() if k != "n_peaks"}
    assert contract.validate_result(missing_peaks) != []
    ghost = {**success, "output_path": str(tmp_path / "missing.h5mu")}
    assert contract.validate_result(ghost) != []
    assert contract.validate_result({"status": "error", "error_type": "NoSharedBarcodes",
                                     "details": "zero shared"}) == []


# ════════════════════════════════════════════════════════════════════════════
# 2. HDF5 `/obs` metadata join (h5py only)
# ════════════════════════════════════════════════════════════════════════════

def _str_ds(group, name, values):
    import h5py
    return group.create_dataset(
        name, data=[str(v).encode() for v in values],
        dtype=h5py.string_dtype(encoding="utf-8"),
    )


def _write_ref_h5ad(path: Path, barcodes: list[str], columns: dict, *,
                    encoding: str = "categorical") -> Path:
    """Write a minimal `.h5ad` holding only `/obs` (no `X`).

    `columns` maps name -> list of labels (None = missing). Encodings:
      categorical: modern `categories` + `codes` group (None -> code -1)
      legacy:      integer codes + `/obs/__categories/<col>`
      nullable:    `nullable-string-array` (`values` + `mask`)
      plain:       plain string dataset
    """
    import h5py
    import numpy as np

    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["column-order"] = np.array(list(columns), dtype=object).astype("S")
        _str_ds(obs, "_index", barcodes)
        for col, values in columns.items():
            cats = sorted({v for v in values if v is not None})
            codes = np.array([cats.index(v) if v is not None else -1 for v in values],
                             dtype=np.int8)
            if encoding == "categorical":
                grp = obs.create_group(col)
                grp.attrs["encoding-type"] = "categorical"
                grp.attrs["ordered"] = False
                _str_ds(grp, "categories", cats)
                grp.create_dataset("codes", data=codes)
            elif encoding == "legacy":
                legacy = obs.require_group("__categories")
                _str_ds(legacy, col, cats)
                obs.create_dataset(col, data=codes)
            elif encoding == "nullable":
                grp = obs.create_group(col)
                grp.attrs["encoding-type"] = "nullable-string-array"
                _str_ds(grp, "values", [v if v is not None else "" for v in values])
                grp.create_dataset("mask", data=np.array([v is None for v in values]))
            elif encoding == "plain":
                _str_ds(obs, col, [v if v is not None else "" for v in values])
            else:  # pragma: no cover - test helper misuse
                raise ValueError(encoding)
    return path


@pytest.fixture
def h5py_mod():
    return pytest.importorskip("h5py")


@pytest.mark.parametrize("encoding", ["categorical", "legacy", "nullable", "plain"])
def test_read_h5ad_obs_columns_encodings(tmp_path, h5py_mod, encoding):
    from aria.utils.mudata_io import read_h5ad_obs_columns

    ref = _write_ref_h5ad(
        tmp_path / f"ref_{encoding}.h5ad",
        ["AAAC-1", "AAAG-1", "AAAT-1"],
        {"subclass": ["CellTypeA", None, "CellTypeB"]},
        encoding=encoding,
    )
    out = read_h5ad_obs_columns(str(ref), ["subclass"])
    assert out["status"] == "success"
    assert out["obs_names"] == ["AAAC-1", "AAAG-1", "AAAT-1"]
    assert out["n_obs"] == 3
    assert out["columns"]["subclass"] == ["CellTypeA", None, "CellTypeB"]


def test_read_h5ad_obs_columns_errors(tmp_path, h5py_mod):
    from aria.utils.mudata_io import read_h5ad_obs_columns

    assert read_h5ad_obs_columns(str(tmp_path / "missing.h5ad"), ["subclass"])[
        "error_type"] == "ReferenceNotFound"
    wrong = tmp_path / "ref.txt"
    wrong.write_text("x")
    assert read_h5ad_obs_columns(str(wrong), ["subclass"])[
        "error_type"] == "InvalidReferenceFormat"
    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["AAAC-1"], {"subclass": ["CellTypeA"]})
    assert read_h5ad_obs_columns(str(ref), ["cell_type"])[
        "error_type"] == "ReferenceColumnMissing"


@pytest.mark.parametrize("encoding", ["categorical", "legacy", "nullable"])
def test_join_missing_reference_label_becomes_unassigned(tmp_path, h5py_mod, encoding):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["AAAC-1", "AAAG-1"],
                          {"subclass": ["CellTypeA", None]}, encoding=encoding)
    out = join_obs_labels_from_h5ad(["AAAC-1", "AAAG-1", "AAAT-1"], str(ref))
    assert out["status"] == "success"
    assert out["mechanism"] == "obs_barcode_join"
    assert out["match_strategy"] == "exact"
    assert out["labels"] == ["CellTypeA", "unassigned", "unassigned"]
    assert out["matched_mask"] == [True, False, False]
    assert out["n_matched"] == 1
    assert out["n_unmatched"] == 2


@pytest.mark.parametrize(
    "ref_barcodes, query, sample_id, strategy",
    [
        (["AAAC-1", "AAAG-1"], ["AAAC-1", "AAAG-1"], None, "exact"),
        (["sample1_AAAC-1", "sample1_AAAG-1"], ["AAAC-1", "AAAG-1"],
         "sample1", "sample_prefix_underscore"),
        (["SAMPLE1_AAAC-1", "SAMPLE1_AAAG-1"], ["AAAC-1", "AAAG-1"],
         "sample1", "sample_prefix_ci"),
        (["sample1#AAAC-1", "sample1#AAAG-1"], ["AAAC-1", "AAAG-1"],
         "sample1", "sample_prefix_hash"),
    ],
    ids=["exact", "prefix_underscore", "prefix_ci", "prefix_hash"],
)
def test_join_strategies(tmp_path, h5py_mod, ref_barcodes, query, sample_id, strategy):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ref_barcodes,
                          {"subclass": ["CellTypeA", "CellTypeB"]})
    out = join_obs_labels_from_h5ad(query, str(ref), sample_id=sample_id)
    assert out["status"] == "success"
    assert out["match_strategy"] == strategy
    assert out["labels"] == ["CellTypeA", "CellTypeB"]
    assert out["match_fraction"] == 1.0


def test_join_sample_prefix_does_not_match_longer_sample_token(tmp_path, h5py_mod):
    """`sample1` must never match `sample12_<bc>` (token boundary required)."""
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["sample12_AAAC-1"],
                          {"subclass": ["CellTypeA"]})
    out = join_obs_labels_from_h5ad(["AAAC-1"], str(ref), sample_id="sample1")
    assert out["status"] == "success"
    assert out["n_matched"] == 0
    assert out["labels"] == ["unassigned"]


def test_join_reference_sample_col_scopes_reference(tmp_path, h5py_mod):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    # Same barcode in two samples: unscoped exact join would be ambiguous.
    ref = _write_ref_h5ad(
        tmp_path / "ref.h5ad",
        ["AAAC-1", "AAAC-1", "AAAG-1"],
        {"subclass": ["CellTypeA", "CellTypeB", "CellTypeB"],
         "donor": ["sample1", "sample2", "sample1"]},
    )
    out = join_obs_labels_from_h5ad(["AAAC-1", "AAAG-1"], str(ref),
                                    sample_id="sample1", reference_sample_col="donor")
    assert out["status"] == "success"
    assert out["n_reference_in_scope"] == 2
    assert out["labels"] == ["CellTypeA", "CellTypeB"]

    ci = join_obs_labels_from_h5ad(["AAAC-1"], str(ref),
                                   sample_id="SAMPLE2", reference_sample_col="donor")
    assert ci["labels"] == ["CellTypeB"]

    unscoped = join_obs_labels_from_h5ad(["AAAC-1"], str(ref))
    assert unscoped["status"] == "error"
    assert unscoped["error_type"] == "AmbiguousBarcodeJoin"

    missing = join_obs_labels_from_h5ad(["AAAC-1"], str(ref),
                                        sample_id="sample9", reference_sample_col="donor")
    assert missing["error_type"] == "SampleIdNotFoundInReference"
    no_sid = join_obs_labels_from_h5ad(["AAAC-1"], str(ref),
                                       reference_sample_col="donor")
    assert no_sid["error_type"] == "InvalidSampleScope"


def test_join_trailing_barcode_within_sample_scope_and_collision_guard(tmp_path, h5py_mod):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    # Library-prefixed barcodes within an explicitly sample-scoped reference
    ref_ok = _write_ref_h5ad(
        tmp_path / "ref_lib.h5ad",
        ["donorA_lib1_AAACCC-1", "donorA_lib1_AAAGGG-1", "donorB_lib1_AAACCC-1"],
        {
            "subclass": ["CellTypeA", "CellTypeB", "CellTypeB"],
            "orig.ident": ["donorA", "donorA", "donorB"],
        },
    )
    ok_res = join_obs_labels_from_h5ad(
        ["AAACCC-1", "AAAGGG-1"],
        str(ref_ok),
        sample_id="donorA",
        reference_sample_col="orig.ident",
    )
    assert ok_res["status"] == "success"
    assert ok_res["match_strategy"] == "trailing_barcode_within_sample_scope"
    assert ok_res["labels"] == ["CellTypeA", "CellTypeB"]
    assert ok_res["match_fraction"] == 1.0

    # Two libraries under the same sample_id sharing a 16-bp barcode MUST fail closed
    ref_collide = _write_ref_h5ad(
        tmp_path / "ref_lib_collide.h5ad",
        ["donorA_lib1_AAACCC-1", "donorA_lib2_AAACCC-1"],
        {
            "subclass": ["CellTypeA", "CellTypeB"],
            "orig.ident": ["donorA", "donorA"],
        },
    )
    bad_res = join_obs_labels_from_h5ad(
        ["AAACCC-1"],
        str(ref_collide),
        sample_id="donorA",
        reference_sample_col="orig.ident",
    )
    assert bad_res["status"] == "error"
    assert bad_res["error_type"] == "AmbiguousBarcodeJoin"


def test_join_fails_closed_on_gem_well_collision(tmp_path, h5py_mod):
    """`AAAC-1` and `AAAC-2` collapse to `AAAC` after GEM stripping -> refuse."""
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["sample1_AAAC"],
                          {"subclass": ["CellTypeA"]})
    out = join_obs_labels_from_h5ad(["AAAC-1", "AAAC-2"], str(ref), sample_id="sample1")
    assert out["status"] == "error"
    assert out["error_type"] == "AmbiguousBarcodeJoin"


def test_join_fails_closed_on_duplicate_reference_barcodes(tmp_path, h5py_mod):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["AAAC-1", "AAAC-1"],
                          {"subclass": ["CellTypeA", "CellTypeB"]})
    out = join_obs_labels_from_h5ad(["AAAC-1"], str(ref))
    assert out["status"] == "error"
    assert out["error_type"] == "AmbiguousBarcodeJoin"


def test_join_rejects_duplicate_query_barcodes(tmp_path, h5py_mod):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["AAAC-1"], {"subclass": ["CellTypeA"]})
    out = join_obs_labels_from_h5ad(["AAAC-1", "AAAC-1"], str(ref))
    assert out["error_type"] == "DuplicateQueryBarcodes"


def test_join_no_match_returns_all_unassigned(tmp_path, h5py_mod):
    from aria.utils.mudata_io import join_obs_labels_from_h5ad

    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["TTTT-1"], {"subclass": ["CellTypeA"]})
    out = join_obs_labels_from_h5ad(["AAAC-1", "AAAG-1"], str(ref))
    assert out["status"] == "success"
    assert out["match_strategy"] is None
    assert out["n_matched"] == 0
    assert out["labels"] == ["unassigned", "unassigned"]


# ════════════════════════════════════════════════════════════════════════════
# 3. Synthetic paired `.h5mu` (mudata / chromatin stack)
# ════════════════════════════════════════════════════════════════════════════

LOW_QC = {
    "min_genes": 1, "min_rna_counts": 1, "max_pct_mt": 100.0,
    "min_peaks": 1, "min_atac_counts": 1,
}
RNA_GENES = [f"GENE{i}" for i in range(1, 9)] + ["MT-CO1", "MT-ND1"]
PEAKS = [f"chr1:{1000 * i}-{1000 * i + 500}" for i in range(1, 11)]


def _barcodes(n: int, prefix: str = "") -> list[str]:
    alphabet = "ACGT"
    out = []
    for i in range(n):
        core = "".join(alphabet[(i >> (2 * k)) & 3] for k in range(6))
        out.append(f"{prefix}{core}-1")
    return out


def _counts(n_rows: int, n_cols: int, seed: int):
    import numpy as np
    rng = np.random.default_rng(seed)
    return rng.poisson(3.0, size=(n_rows, n_cols)).astype(np.float32) + 1.0


def _write_paired_h5mu(
    path: Path,
    *,
    rna_obs: list[str],
    atac_obs: list[str] | None = None,
    rna_genes: list[str] = RNA_GENES,
    rna_x=None,
    rna_counts_layer=None,
    atac_x=None,
) -> Path:
    import anndata as ad
    import mudata as md
    import pandas as pd
    from scipy import sparse

    atac_obs = list(rna_obs) if atac_obs is None else list(atac_obs)
    if rna_x is None:
        rna_x = _counts(len(rna_obs), len(rna_genes), seed=7)
    if atac_x is None:
        atac_x = _counts(len(atac_obs), len(PEAKS), seed=11)
    rna = ad.AnnData(
        X=sparse.csr_matrix(rna_x),
        obs=pd.DataFrame(index=pd.Index(rna_obs)),
        var=pd.DataFrame(index=pd.Index(rna_genes)),
    )
    if rna_counts_layer is not None:
        rna.layers["counts"] = sparse.csr_matrix(rna_counts_layer)
    atac = ad.AnnData(
        X=sparse.csr_matrix(atac_x),
        obs=pd.DataFrame(index=pd.Index(atac_obs)),
        var=pd.DataFrame(index=pd.Index(PEAKS)),
    )
    md.MuData({"rna": rna, "atac": atac}).write_h5mu(str(path))
    return path


@pytest.fixture
def mudata_mod():
    return pytest.importorskip("mudata")


def test_read_h5mu_paired_preserves_rna_order_deterministically(tmp_path, mudata_mod):
    from aria.utils.mudata_io import read_h5mu_paired

    rna_obs = _barcodes(12)
    atac_obs = list(reversed(rna_obs[2:])) + ["TTTTTT-1"]  # permuted + extra
    path = _write_paired_h5mu(tmp_path / "perm.h5mu", rna_obs=rna_obs,
                              atac_obs=atac_obs)

    first = read_h5mu_paired(str(path))
    second = read_h5mu_paired(str(path))
    assert first["status"] == "success"
    assert first["rna_modality"] == "rna"
    assert first["atac_modality"] == "atac"
    assert first["exact_same_order"] is False
    expected = [bc for bc in rna_obs if bc in set(atac_obs)]
    assert first["shared_barcodes"] == expected
    assert second["shared_barcodes"] == first["shared_barcodes"]
    assert first["n_shared_cells"] == 10
    assert first["n_rna_only"] == 2
    assert first["n_atac_only"] == 1
    assert first["rna_counts_source"] == "X"
    assert first["rna_is_raw_counts"] is True


def test_read_h5mu_paired_rejects_fractional_counts(tmp_path, mudata_mod):
    from aria.utils.mudata_io import read_h5mu_paired

    rna_obs = _barcodes(6)
    frac = _counts(6, len(RNA_GENES), seed=3) + 0.25
    path = _write_paired_h5mu(tmp_path / "frac.h5mu", rna_obs=rna_obs, rna_x=frac)
    out = read_h5mu_paired(str(path))
    assert out["status"] == "error"
    assert out["error_type"] == "NonIntegerCounts"


def test_read_h5mu_paired_falls_back_to_counts_layer(tmp_path, mudata_mod):
    import numpy as np
    from aria.utils.mudata_io import read_h5mu_paired

    rna_obs = _barcodes(6)
    raw = _counts(6, len(RNA_GENES), seed=3)
    normalized = np.log1p(raw / raw.sum(axis=1, keepdims=True) * 1e4)
    path = _write_paired_h5mu(tmp_path / "layer.h5mu", rna_obs=rna_obs,
                              rna_x=normalized, rna_counts_layer=raw)
    out = read_h5mu_paired(str(path))
    assert out["status"] == "success"
    assert out["rna_counts_source"] == "layers['counts']"
    assert out["rna_is_raw_counts"] is True


def test_read_h5mu_paired_rejects_duplicate_barcodes(tmp_path, mudata_mod):
    import h5py
    from aria.utils.mudata_io import read_h5mu_paired

    rna_obs = _barcodes(6)
    path = _write_paired_h5mu(tmp_path / "dup.h5mu", rna_obs=rna_obs)
    # Inject a duplicate barcode at the storage level (MuData would refuse to
    # build it in memory), mimicking a malformed external `.h5mu`.
    with h5py.File(path, "r+") as f:
        idx_key = f["mod/atac/obs"].attrs["_index"]
        ds = f[f"mod/atac/obs/{idx_key}"]
        ds[1] = ds[0]
    out = read_h5mu_paired(str(path))
    assert out["status"] == "error"
    assert out["error_type"] == "DuplicateBarcodes"


def test_read_h5mu_paired_rejects_disjoint_barcodes(tmp_path, mudata_mod):
    from aria.utils.mudata_io import read_h5mu_paired

    rna_obs = _barcodes(6)
    atac_obs = [bc.replace("-1", "-2") for bc in rna_obs]  # no silent normalization
    path = _write_paired_h5mu(tmp_path / "disjoint.h5mu", rna_obs=rna_obs,
                              atac_obs=atac_obs)
    out = read_h5mu_paired(str(path))
    assert out["status"] == "error"
    assert out["error_type"] == "NoSharedBarcodes"


def test_read_h5mu_paired_missing_file(tmp_path, mudata_mod):
    from aria.utils.mudata_io import read_h5mu_paired

    out = read_h5mu_paired(str(tmp_path / "missing.h5mu"))
    assert out["status"] == "error"


def _run_qc(**params):
    from aria.scripts.integration_paired_qc import integration_paired_qc
    return integration_paired_qc({**LOW_QC, **params})


def test_paired_qc_ensembl_ids_never_fabricate_zero_mito(tmp_path, mudata_mod):
    rna_obs = _barcodes(8)
    genes = [f"ENSG{i:011d}" for i in range(1, 11)]
    path = _write_paired_h5mu(tmp_path / "ens.h5mu", rna_obs=rna_obs, rna_genes=genes)
    out = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"))
    assert out["status"] == "success"
    assert out["rna_qc_summary"]["n_mt_genes"] == 0
    assert out["rna_qc_summary"]["median_pct_mt_before"] is None
    assert out["rna_qc_summary"]["median_pct_mt_after"] is None
    assert out["thresholds_used"]["max_pct_mt"] is None
    assert any(m.startswith("pct_counts_mt") for m in out["metrics_not_computed"])
    assert out["atac_qc_summary"]["frip"] is None
    assert out["atac_qc_summary"]["tss_enrichment"] is None
    assert any(m.startswith("frip") for m in out["metrics_not_computed"])
    assert any(m.startswith("tss_enrichment") for m in out["metrics_not_computed"])


def test_paired_qc_measures_mito_from_mt_symbols(tmp_path, mudata_mod):
    rna_obs = _barcodes(8)
    path = _write_paired_h5mu(tmp_path / "mt.h5mu", rna_obs=rna_obs)
    out = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"))
    assert out["status"] == "success"
    assert out["rna_qc_summary"]["n_mt_genes"] == 2
    assert out["rna_qc_summary"]["mt_gene_source"] == "var_names"
    pct = out["rna_qc_summary"]["median_pct_mt_before"]
    assert pct is not None and 0.0 < pct < 100.0
    assert not any(m.startswith("pct_counts_mt") for m in out["metrics_not_computed"])
    assert out["atac_qc_summary"]["frip"] is None
    assert out["atac_qc_summary"]["tss_enrichment"] is None


def test_paired_qc_zero_cells_after_qc(tmp_path, mudata_mod):
    path = _write_paired_h5mu(tmp_path / "zero.h5mu", rna_obs=_barcodes(8))
    out = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"),
                  min_rna_counts=10 ** 9)
    assert out["status"] == "error"
    assert out["error_type"] == "ZeroCellsAfterQC"
    assert not (tmp_path / "out").exists() or not any((tmp_path / "out").iterdir())


def test_paired_qc_low_barcode_overlap(tmp_path, mudata_mod):
    rna_obs = _barcodes(10)
    atac_obs = rna_obs[:4] + _barcodes(6, prefix="GG")
    path = _write_paired_h5mu(tmp_path / "low.h5mu", rna_obs=rna_obs,
                              atac_obs=atac_obs)
    out = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"))
    assert out["status"] == "error"
    assert out["error_type"] == "LowBarcodeOverlap"

    partial = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"),
                      allow_partial_overlap=True)
    assert partial["status"] == "success"
    assert partial["n_cells_before"] == 4


def test_paired_qc_missing_input():
    out = _run_qc()
    assert out["status"] == "error"
    assert out["error_type"] == "MissingInput"


def test_paired_qc_joint_filter_atomic_write_and_label_join(tmp_path, mudata_mod, h5py_mod):
    import mudata as md
    import numpy as np
    from aria.utils.mudata_io import read_h5mu_paired

    rna_obs = _barcodes(12)
    atac_obs = list(reversed(rna_obs))  # permuted ATAC order
    rna_x = _counts(12, len(RNA_GENES), seed=21)
    rna_x[0, :] = 0.0  # fails RNA QC
    atac_x = _counts(12, len(PEAKS), seed=22)
    atac_x[0, :] = 0.0  # ATAC row 0 == rna_obs[-1] -> fails ATAC QC
    path = _write_paired_h5mu(tmp_path / "sample1_paired.h5mu", rna_obs=rna_obs,
                              atac_obs=atac_obs, rna_x=rna_x, atac_x=atac_x)

    labels = ["CellTypeA" if i % 2 == 0 else "CellTypeB" for i in range(1, 10)]
    ref = _write_ref_h5ad(
        tmp_path / "ref.h5ad",
        [f"sample1_{bc}" for bc in rna_obs[1:10]],  # rna_obs[10] has no reference
        {"subclass": labels},
    )

    out_dir = tmp_path / "out"
    params = dict(h5mu_path=str(path), output_dir=str(out_dir), sample_id="sample1",
                  reference_h5ad=str(ref), reference_label_col="subclass")
    out = _run_qc(**params)
    assert out["status"] == "success", out
    assert out["n_cells_before"] == 12
    assert out["n_pass_rna"] == 11
    assert out["n_pass_atac"] == 11
    assert out["n_cells_after"] == 10
    assert out["n_genes"] == len(RNA_GENES)
    assert out["n_peaks"] == len(PEAKS)
    assert out["exact_same_input_order"] is False

    out_path = Path(out["output_path"])
    assert out_path == out_dir / "paired_qc_sample1.h5mu"
    assert out_path.exists()
    leftovers = [p.name for p in out_dir.iterdir() if ".tmp." in p.name]
    assert leftovers == []

    join = out["label_join"]
    assert join["ran"] is True
    assert join["mechanism"] == "obs_barcode_join"
    assert join["match_strategy"] == "sample_prefix_underscore"
    assert join["n_query"] == 10
    assert join["n_matched"] == 9
    assert join["n_unmatched"] == 1

    mdata_out = md.read_h5mu(str(out_path))
    rna_names = list(mdata_out.mod["rna"].obs_names)
    assert rna_names == list(mdata_out.mod["atac"].obs_names)
    assert rna_names == rna_obs[1:11]
    rna_labels = list(mdata_out.mod["rna"].obs["subclass"].astype(str))
    expected = dict(zip(rna_obs[1:10], labels))
    assert rna_labels == [expected.get(bc, "unassigned") for bc in rna_names]
    assert list(mdata_out.mod["atac"].obs["subclass"].astype(str)) == rna_labels
    assert set(mdata_out.mod["rna"].obs["sample_id"].astype(str)) == {"sample1"}
    assert np.all(np.asarray(mdata_out.mod["atac"].obs["total_peak_counts"]) > 0)

    reread = read_h5mu_paired(str(out_path))
    assert reread["status"] == "success"
    assert reread["exact_same_order"] is True
    assert reread["n_shared_cells"] == 10

    # Reproducibility: an identical second run yields the identical result.
    out2 = _run_qc(**{**params, "output_dir": str(tmp_path / "out2")})
    strip = lambda d: {k: v for k, v in d.items() if k != "output_path"}  # noqa: E731
    assert strip(out2) == strip(out)
    second = md.read_h5mu(out2["output_path"])
    assert list(second.mod["rna"].obs_names) == rna_names


def test_paired_qc_require_reference_match_fails_closed(tmp_path, mudata_mod, h5py_mod):
    rna_obs = _barcodes(6)
    path = _write_paired_h5mu(tmp_path / "p.h5mu", rna_obs=rna_obs)
    ref = _write_ref_h5ad(tmp_path / "ref.h5ad", ["TTTTTT-9"], {"subclass": ["CellTypeA"]})
    out = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"),
                  reference_h5ad=str(ref))
    assert out["status"] == "error"
    assert out["error_type"] == "NoReferenceBarcodesMatched"

    lenient = _run_qc(h5mu_path=str(path), output_dir=str(tmp_path / "out"),
                      reference_h5ad=str(ref), require_reference_match=False)
    assert lenient["status"] == "success"
    assert lenient["label_join"]["n_matched"] == 0
    assert any("unassigned" in w for w in lenient["warnings"])


# ════════════════════════════════════════════════════════════════════════════
# 4. Real-data validation (dataset-gated)
# ════════════════════════════════════════════════════════════════════════════

_REAL_H5MU_CANDIDATES = (
    "/home/medusa/Samael/Erosion/data_inputs/muon_processed/hc11_paired.h5mu",
    "/home/medusa/Samael/Oligos/data/interim/paired/hc11_paired.h5mu",
)
_REAL_REF_CANDIDATES = (
    "/home/medusa/Samael/Erosion/raw/GSE278576_hippocampus_RNA.h5ad",
    "/home/medusa/Samael/Single/Samael_Final_Pipeline/data_inputs/reference/"
    "GSE278576_hippocampus_RNA.h5ad",
    "/home/medusa/Samael/Erosion/data_inputs/reference/GSE278576_hippocampus_RNA.h5ad",
)


def _first_existing(env_var: str, candidates: tuple[str, ...]) -> str | None:
    override = os.environ.get(env_var)
    if override:
        return override if Path(override).exists() else None
    return next((c for c in candidates if Path(c).exists()), None)


@pytest.fixture
def real_h5mu(mudata_mod):
    path = _first_existing("ARIA_V48_H5MU_VALIDATION_INPUT", _REAL_H5MU_CANDIDATES)
    if path is None:
        pytest.skip("real paired hc11 .h5mu not available")
    return path


def test_real_hc11_paired_qc(tmp_path, real_h5mu):
    from aria.utils.mudata_io import read_h5mu_paired

    params = dict(h5mu_path=real_h5mu, output_dir=str(tmp_path), sample_id="hc11",
                  organism="Homo sapiens", genome="GRCh38")
    ref = _first_existing("ARIA_V48_REF_H5AD_INPUT", _REAL_REF_CANDIDATES)
    if ref is not None:
        params.update(reference_h5ad=ref, reference_label_col="subclass",
                      reference_sample_col="orig.ident",
                      require_reference_match=True)

    from aria.scripts.integration_paired_qc import integration_paired_qc
    out = integration_paired_qc(params)
    print("\n[V48-I1 real hc11]", {k: out.get(k) for k in (
        "status", "n_cells_before", "n_pass_rna", "n_pass_atac", "n_cells_after",
        "pct_removed", "n_genes", "n_peaks", "rna_counts_source",
        "atac_counts_source", "atac_is_binary")})
    print("[V48-I1 real hc11] rna_qc", out.get("rna_qc_summary"))
    print("[V48-I1 real hc11] atac_qc", out.get("atac_qc_summary"))
    print("[V48-I1 real hc11] label_join", out.get("label_join"))

    assert out["status"] == "success", out
    assert out["n_cells_before"] == 3143
    assert out["n_genes"] == 36601
    assert out["n_peaks"] == 60990
    assert 0 < out["n_cells_after"] <= 3143
    assert out["atac_qc_summary"]["frip"] is None
    assert out["atac_qc_summary"]["tss_enrichment"] is None
    assert out["rna_qc_summary"]["n_mt_genes"] == 13
    assert out["rna_qc_summary"]["median_pct_mt_before"] is not None

    join = out["label_join"]
    if ref is not None:
        assert join["ran"] is True
        assert join["n_query"] == out["n_cells_after"]
        assert join["n_matched"] + join["n_unmatched"] == join["n_query"]
        assert sum(join["label_counts"].values()) == join["n_query"]

    reread = read_h5mu_paired(out["output_path"])
    assert reread["status"] == "success"
    assert reread["exact_same_order"] is True
    assert reread["n_shared_cells"] == out["n_cells_after"]


def test_real_hc11_paired_qc_via_environment_manager(tmp_path):
    """Production IPC E2E: `aria-env` dispatches `integration_paired_qc.py` into `aria-chromatin-env`."""
    _require_orchestration_deps()
    from aria.utils.environment_manager import EnvironmentManager, MissingEnvironment

    h5mu_path = _first_existing("ARIA_V48_H5MU_VALIDATION_INPUT", _REAL_H5MU_CANDIDATES)
    if h5mu_path is None:
        pytest.skip("real paired hc11 .h5mu not available")

    mgr = EnvironmentManager(workspace_dir=str(tmp_path / "aria_ws"))
    if not mgr.check_environments().get("chromatin", False):
        pytest.skip("aria-chromatin-env not installed")

    out_dir = tmp_path / "paired_out"
    params: dict = {
        "h5mu_path": h5mu_path,
        "output_dir": str(out_dir),
        "sample_id": "hc11",
        "organism": "Homo sapiens",
        "genome": "GRCh38",
    }
    ref = _first_existing("ARIA_V48_REF_H5AD_INPUT", _REAL_REF_CANDIDATES)
    if ref is not None:
        params.update(
            reference_h5ad=ref,
            reference_label_col="subclass",
            reference_sample_col="orig.ident",
            require_reference_match=True,
        )

    try:
        out = mgr.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/integration_paired_qc.py",
            params=params,
        )
    except MissingEnvironment:
        pytest.skip("aria-chromatin-env unavailable")

    print("\n[V48-I1 IPC E2E hc11]", {k: out.get(k) for k in (
        "status", "n_cells_before", "n_pass_rna", "n_pass_atac", "n_cells_after",
        "pct_removed", "n_genes", "n_peaks", "rna_counts_source",
        "atac_counts_source", "atac_is_binary")})
    print("[V48-I1 IPC E2E hc11] rna_qc", out.get("rna_qc_summary"))
    print("[V48-I1 IPC E2E hc11] atac_qc", out.get("atac_qc_summary"))
    print("[V48-I1 IPC E2E hc11] label_join", out.get("label_join"))

    assert out["status"] == "success", out
    assert out["n_cells_before"] == 3143
    assert out["n_genes"] == 36601
    assert out["n_peaks"] == 60990
    assert 0 < out["n_cells_after"] <= 3143
    assert out["atac_qc_summary"]["frip"] is None
    assert out["atac_qc_summary"]["tss_enrichment"] is None
    assert out["rna_qc_summary"]["n_mt_genes"] == 13
    assert Path(out["output_path"]).exists()

    if ref is not None:
        join = out["label_join"]
        assert join["ran"] is True
        assert join["n_matched"] > 0
        assert join["match_fraction"] > 0.50
