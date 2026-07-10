"""Preprint-readiness audit B1: technical replicates are not biological n.

The CP2.5 choice must change the executable design.  Technical libraries are
summed within a condition-scoped biological unit before inference; the original
library membership remains auditable.  An unresolved or under-replicated design
fails closed instead of treating libraries as independent observations.
"""
from __future__ import annotations

import json

import pytest


def _design_agent(groups):
    from aria.agents.design_agent import DesignAgent

    agent = object.__new__(DesignAgent)
    agent._confirmed_groups = groups
    agent._main_factor = "condition"
    agent._batch_covariate = None
    agent._organism = "synthetic"
    agent._genome = "test"
    agent._pseudobulk_design = {}
    agent._inferred_design = {}
    agent._replicate_handling = None
    agent._publish_confirm_checkpoint = lambda: None
    return agent


def _groups(n_biological=2):
    groups = {}
    for condition in ("ctrl", "stim"):
        members = []
        for biological in range(1, n_biological + 1):
            for technical in (1, 2):
                members.append(
                    f"{condition}_unit{biological}_rep{technical}"
                )
        groups[condition] = members
    return groups


def test_cp25_choice_changes_units_formula_declaration_and_df():
    biological = _design_agent(_groups())
    biological._handle_pseudorep_response(
        "Yes — they are independent biological replicates"
    )
    biological_design = biological._build_design()

    technical = _design_agent(_groups())
    technical._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )
    technical_design = technical._build_design()

    assert biological_design["replicates"] == {"ctrl": 4, "stim": 4}
    assert technical_design["replicates"] == {"ctrl": 2, "stim": 2}
    assert biological_design["n_total_samples"] == 8
    assert technical_design["n_total_samples"] == 4
    assert technical_design["n_input_libraries"] == 8
    assert biological_design["nominal_residual_df"] == 6
    assert technical_design["nominal_residual_df"] == 2
    assert (
        biological_design["analysis_design_formula"]
        != technical_design["analysis_design_formula"]
    )
    assert technical_design["design_formula"] == "~ condition"
    assert technical_design["pseudobulk"]["replicate_col"] == "biological_unit"
    assert technical_design["replicate_handling"]["mode"] == "technical_aggregate"
    assert set(technical_design["replicate_handling"]["sample_to_unit"]) == {
        sample for members in _groups().values() for sample in members
    }
    assert all(
        len(unit["members"]) == 2
        for unit in technical_design["replicate_handling"]["units"]
    )


def test_technical_suffix_detection_does_not_truncate_donor_ids():
    from aria.agents.design_agent import DesignAgent

    assert DesignAgent._technical_replicate_root("donor1") == "donor1"
    assert DesignAgent._technical_replicate_root("donor1_rep2") == "donor1"
    assert DesignAgent._technical_replicate_root("donor1-r2") == "donor1"
    assert DesignAgent._technical_replicate_root("donor1_2") == "donor1"


def test_cp25_uncertainty_fails_closed():
    agent = _design_agent(_groups())

    result = agent._handle_pseudorep_response("Not sure — proceed anyway")

    assert result["status"] == "cancelled"
    assert result["reason"] == "replicate_structure_unresolved"


def test_technical_merge_blocks_when_only_one_biological_unit_per_group():
    agent = _design_agent(_groups(n_biological=1))
    agent._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )

    design = agent._build_design()
    result = agent._handle_confirm_response("Yes — proceed")

    assert design["design_status"] == "blocking"
    assert design["replicates"] == {"ctrl": 1, "stim": 1}
    assert result["status"] == "cancelled"
    assert result["reason"] == "insufficient_biological_replicates"


def test_bulk_metadata_and_counts_use_biological_units(tmp_path):
    pd = pytest.importorskip("pandas")
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.scripts.rna_bulk_de import _aggregate_technical_replicates

    design_agent = _design_agent(_groups())
    design_agent._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )
    design = design_agent._build_design()
    group_labels = {
        sample: group
        for group, samples in design["groups"].items()
        for sample in samples
    }
    replicate_units = BulkRNAAgent._resolve_replicate_units(
        design, list(group_labels)
    )
    metadata_path = BulkRNAAgent._write_design_metadata(
        group_labels,
        "condition",
        tmp_path,
        replicate_units=replicate_units,
    )
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    assert list(metadata.columns) == ["condition", "biological_unit"]

    counts = pd.DataFrame(
        {
            sample: [idx + 1, (idx + 1) * 10]
            for idx, sample in enumerate(group_labels)
        },
        index=["feature1", "feature2"],
    )
    aggregated, unit_metadata, provenance = _aggregate_technical_replicates(
        counts,
        metadata,
        design_factor="condition",
        unit_col="biological_unit",
        covariates=[],
    )

    assert aggregated.shape == (2, 4)
    assert unit_metadata.shape == (4, 1)
    assert provenance["n_input_libraries"] == 8
    assert provenance["n_biological_units"] == 4
    assert provenance["replicates_per_condition"] == {"ctrl": 2, "stim": 2}
    assert provenance["residual_degrees_of_freedom"] == 2
    first_unit = metadata.iloc[0]["biological_unit"]
    first_members = metadata.index[metadata["biological_unit"] == first_unit]
    assert aggregated[first_unit].tolist() == counts[first_members].sum(axis=1).tolist()


def test_scrna_injection_writes_biological_unit_replicate_column(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    from aria.scripts.rna_inject_condition import inject
    from aria.utils.safe_h5ad import read_h5ad

    groups = _groups()
    samples = [sample for members in groups.values() for sample in members]
    obs_samples = [sample for sample in samples for _ in range(2)]
    adata = ad.AnnData(
        X=np.ones((len(obs_samples), 2)),
        obs=pd.DataFrame(
            {"sample_id": obs_samples},
            index=[f"cell{i}" for i in range(len(obs_samples))],
        ),
        var=pd.DataFrame(index=["feature1", "feature2"]),
    )
    source = tmp_path / "input.h5ad"
    output = tmp_path / "with_design.h5ad"
    adata.write_h5ad(source)

    agent = _design_agent(groups)
    agent._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )
    design = agent._build_design()
    result = inject(
        {
            "data_path": str(source),
            "groups": groups,
            "factor": "condition",
            "replicate_units": design["replicate_handling"]["sample_to_unit"],
            "output_path": str(output),
        }
    )

    assert result["status"] == "success"
    assert result["replicate_col"] == "biological_unit"
    assert result["n_biological_units"] == 4
    written = read_h5ad(output)
    assert "biological_unit" in written.obs
    assert written.obs["biological_unit"].nunique() == 4
    assert written.obs.groupby("biological_unit", observed=True)["condition"].nunique().max() == 1


def test_design_serializes_replicate_membership_for_audit():
    agent = _design_agent(_groups())
    agent._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )

    payload = json.loads(json.dumps(agent._build_design()))

    assert payload["experimental_unit"] == "biological_unit"
    assert payload["replicate_handling"]["source"] == "user_confirmed_cp2.5"


def test_bulk_methods_disclose_technical_aggregation():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    lines = BulkRnaNarrator().methods(
        "bulk_rna_agent",
        {
            "findings": {
                "design_used": "~condition",
                "contrasts": [],
                "technical_replicate_aggregation": {
                    "ran": True,
                    "method": "sum_raw_counts_by_biological_unit",
                    "n_input_libraries": 8,
                    "n_biological_units": 4,
                    "replicates_per_condition": {"ctrl": 2, "stim": 2},
                    "residual_degrees_of_freedom": 2,
                },
            }
        },
    )
    text = " ".join(lines)

    assert "8 technical libraries" in text
    assert "4 biological units" in text
    assert "residual degrees of freedom: 2" in text


def test_bulk_agent_forwards_cp25_contract_to_de_script(tmp_path):
    pd = pytest.importorskip("pandas")
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    groups = _groups()
    design_agent = _design_agent(groups)
    design_agent._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )
    design = design_agent._build_design()
    design["plan_contrasts"] = [
        {"numerator": "stim", "denominator": "ctrl"}
    ]
    samples = [sample for members in groups.values() for sample in members]
    counts_path = tmp_path / "counts.tsv"
    pd.DataFrame(
        {sample: [10 + idx, 20 + idx] for idx, sample in enumerate(samples)},
        index=["feature1", "feature2"],
    ).to_csv(counts_path, sep="\t")

    captured = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            captured.update({"stack": stack, "script_path": script_path, **params})
            return {"status": "success", "contrasts": []}

    class _Memory:
        def store_decision(self, **kwargs):
            return kwargs

    agent = object.__new__(BulkRNAAgent)
    agent.env = _Env()
    agent.memory = _Memory()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_finding = lambda *args, **kwargs: None
    agent._publish_findings = lambda *args, **kwargs: None
    agent._record_methodology_decisions = lambda *args, **kwargs: None

    result = agent.run(
        "experiment-b1",
        {
            "exp_context": {
                "modalities": {"bulk_RNA": [str(counts_path)]},
                "design": design,
            },
            "biological_intent": {},
        },
    )

    assert result["status"] == "done"
    assert captured["technical_replicate_col"] == "biological_unit"
    metadata = pd.read_csv(captured["metadata_file"], sep="\t")
    assert metadata["biological_unit"].nunique() == 4


def test_preflight_and_planning_count_biological_units_not_libraries():
    from aria.agents.design_intelligence import DesignIntelligence
    from aria.agents.modality_audit import ScRNAAuditAgent
    from aria.utils.design_power import assess_design_power

    agent = _design_agent(_groups())
    agent._pseudobulk_design = {
        "groupby_col": "cell_type",
        "comparisons": [["stim", "ctrl"]],
    }
    agent._handle_pseudorep_response(
        "No — technical replicates; merge by biological unit"
    )
    design = agent._build_design()
    design["plan_contrasts"] = [
        {"numerator": "stim", "denominator": "ctrl"}
    ]
    context = {
        "modalities": {"bulk_RNA": ["counts.tsv"], "scRNA": ["cells.h5ad"]},
        "design": design,
        "run_optional_supported": True,
    }

    profile = DesignIntelligence()._bulk_rna_profile("bulk_RNA", context, {})
    power = assess_design_power(context)
    scrna = ScRNAAuditAgent().audit(context)

    assert not any("with biological replicates" in item for item in profile["recommended"])
    assert any("low-power" in item for item in profile["optional"])
    assert power["checks"]["replicates_per_condition"] == {"ctrl": 2, "stim": 2}
    assert scrna["checks"]["pseudobulk_replicates"]["replicates_per_condition"] == {
        "ctrl": 2,
        "stim": 2,
    }
