"""Live DE/DA fits and narrative extraction across isolated scientific stacks.

PyDESeq2 lives in aria-rna-env (not every chromatin installation includes it).
Workers use that environment when necessary; no fitted statistics or LOBO
results are mocked. Narrative imports run separately in aria-env, keeping the
scientific processes free of agents, litellm, and pydantic.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
ENV_ROOT = Path(os.environ.get("CONDA_ENVS_PATH", Path.home() / "anaconda3" / "envs"))
SCRIPTS = (
    "rna_bulk_de", "rna_pseudobulk_de", "chromatin_diffacc",
    "chromatin_bulk_diffacc",
)


def _python_with(dependency, environment):
    if importlib.util.find_spec(dependency) is not None:
        return sys.executable
    python = ENV_ROOT / environment / "bin/python"
    assert python.exists(), f"Need {dependency} or the {environment} Python"
    return str(python)


def _worker_env():
    env = os.environ.copy()
    env.update({
        "PYTHONPATH": str(ROOT), "MPLBACKEND": "Agg",
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1",
    })
    return env


@pytest.fixture(scope="module")
def live_results(tmp_path_factory):
    work = tmp_path_factory.mktemp("live_lobo")
    proc = subprocess.run(
        [_python_with("pydeseq2", "aria-rna-env"), str(Path(__file__)), str(work)],
        cwd=ROOT, env=_worker_env(), capture_output=True, text=True, timeout=240,
    )
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-6000:]
    return json.loads((work / "results.json").read_text())


def _axis(claim, agent_results, exp_ctx=None):
    """Use the narrative environment without importing agents into a DE stack."""
    proc = subprocess.run(
        [_python_with("litellm", "aria-env"), "-c", """
import json, sys
from aria.agents.narrative.falsification_certificates import build_falsification_certificates
payload = json.load(sys.stdin)
bundle = build_falsification_certificates(**payload)
assert bundle['feeds_into_parameter_selection'] is False
assert bundle['causal_promotion_from_robustness'] is False
print(json.dumps(bundle['certificates'][0]['lobo_biological_unit_jackknife']))
"""],
        input=json.dumps({"claims": [claim], "agent_results": agent_results,
                          "exp_ctx": exp_ctx or {}}),
        cwd=ROOT, env=_worker_env(), capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def _claim(modality, analysis="differential_expression", tier="associative"):
    return {"claim_id": "live.de", "modality": modality, "analysis": analysis, "tier": tier}


def _assert_measured(summary, folds=8, units=4):
    assert summary["status"] == "measured", summary
    assert summary["n_jackknife_folds"] == folds
    assert summary["n_units_numerator"] == summary["n_units_denominator"] == units
    assert summary["biological_unit_enforced"] is True
    assert summary["n_features_universe"] == 120
    assert summary["n_features_evaluated"] > 0
    assert summary["n_features_lobo_stable"] + summary["n_single_unit_driven"] == summary["n_features_evaluated"]
    assert len(summary["per_feature"]) <= 200


def test_live_bulk_rna_biological_aggregation_and_outlier(live_results):
    result = live_results["bulk_rna"]
    assert result["technical_replicate_aggregation"]["n_biological_units"] == 8
    comp = result["contrasts"][0]
    assert comp["status"] == "success", comp
    lobo = comp["lobo_jackknife"]
    _assert_measured(lobo)
    assert "G_STABLE" in comp["all_sig_gene_ids"]
    assert lobo["per_feature"]["G_STABLE"]["single_unit_driven"] is False
    assert set(lobo["per_feature"]) == set(comp["all_sig_gene_ids"])
    # With no primary discoveries the live engine evaluates up to 50 genes,
    # including a planted outlier that DESeq2 may correctly filter via Cook's D.
    fallback = live_results["bulk_rna_no_hits"]["contrasts"][0]
    assert fallback["n_significant"] == 0
    assert fallback["lobo_jackknife"]["n_features_evaluated"] == 40
    assert fallback["lobo_jackknife"]["per_feature"]["G_OUTLIER"]["single_unit_driven"] is True


def test_live_pseudobulk_uses_each_finalized_block(live_results):
    pb = live_results["scrna"]
    assert pb["status"] == "success"
    assert pb["multiple_testing"]["fdr_strategy"] == "global"
    for group, units in (("A", 4), ("B", 3)):
        for name, comp in pb["per_group"][group]["per_comparison"].items():
            assert comp["status"] == "success", comp
            lobo = comp["lobo_jackknife"]
            _assert_measured(lobo, folds=units, units=units)
            assert lobo["design_mode"] == "paired_donor"
            assert name == f"{lobo['numerator']}_vs_{lobo['denominator']}"
            assert set(lobo["per_feature"]) == {r["gene"] for r in comp["all_sig"]}
            stable = lobo["per_feature"]["G_STABLE"]
            assert stable["single_unit_driven"] is False
            assert (stable["full_lfc"] > 0) == (lobo["numerator"] == "treated")


@pytest.mark.parametrize("key,folds,mode", [
    ("scatac_unpaired", 8, "unpaired_biological_replicate"),
    ("scatac_paired", 4, "paired_donor"),
])
def test_live_scatac_designs(live_results, key, folds, mode):
    result = live_results[key]
    assert result["status"] == "success"
    comp = result["pseudobulk"]["comparisons"][0]
    assert comp["status"] == "success", comp
    lobo = comp["lobo_jackknife"]
    _assert_measured(lobo, folds=folds)
    assert lobo["design_mode"] == mode
    assert comp["paired_block_modeled"] == (mode == "paired_donor")
    assert lobo["per_feature"]["chr1:0-100"]["single_unit_driven"] is False


def test_live_scatac_n2_is_not_testable(live_results):
    comp = live_results["scatac_n2"]["pseudobulk"]["comparisons"][0]
    assert comp["status"] == "success", comp
    lobo = comp["lobo_jackknife"]
    assert lobo["status"] == "not_testable_insufficient_replicates"
    assert lobo["minimum_required_per_arm"] == 3
    assert lobo["n_units_numerator"] == lobo["n_units_denominator"] == 2
    assert lobo["per_feature"] == {}


def test_live_bulk_chromatin_preserves_spikein_normalization(live_results):
    result = live_results["bulk_atac"]
    assert result["ran"] is True
    assert result["spikein_normalization"]["used"] is True
    lobo = result["comparisons"][0]["lobo_jackknife"]
    _assert_measured(lobo)
    assert lobo["normalization"] == "spikein_or_external_size_factor_log2"
    assert lobo["per_feature"]["chr1:0-100"]["single_unit_driven"] is False
    assert lobo == live_results["bulk_atac_expected_lobo"]


@pytest.mark.parametrize("engine", SCRIPTS)
def test_unexpected_lobo_error_does_not_abort_de_da(live_results, engine):
    comp = live_results["errors"][engine]
    assert comp["status"] == "success", comp
    assert comp["lobo_jackknife"] == {"status": "error", "details": "x" * 200}


@pytest.mark.parametrize("modality,source", [
    ("bulk_rna", "bulk_rna"), ("bulk RNA-seq", "bulk_rna"),
    ("scrna", "scrna"), ("scRNA-seq", "scrna"),
    ("chromatin", "scatac_paired"),
])
def test_certificates_extract_live_modality_matched_results(live_results, modality, source):
    results = {
        "bulk_rna_agent": {"findings": live_results["bulk_rna"]},
        "scrna_agent": {"findings": {"scRNA": {"findings": {"pseudobulk_de": live_results["scrna"]}}}},
        "chromatin_agent": {"findings": {"scATAC": {"findings": {"diffacc": live_results["scatac_paired"]}}}},
    }
    analysis = "differential_accessibility" if modality == "chromatin" else "differential_expression"
    axis = _axis(_claim(modality, analysis), results, {
        "design": {"donor_col": "donor", "replicates": {"ctrl": 4, "treated": 4}},
    })
    if source == "bulk_rna":
        summaries = [c["lobo_jackknife"] for c in live_results[source]["contrasts"]]
    elif source == "scrna":
        summaries = [c["lobo_jackknife"] for g in live_results[source]["per_group"].values()
                     for c in g["per_comparison"].values()]
    else:
        summaries = [c["lobo_jackknife"] for c in live_results[source]["pseudobulk"]["comparisons"]]
    assert axis["status"] == "measured"
    assert axis["biological_unit"] == "donor"
    assert axis["n_contrasts_measured"] == len(summaries)
    assert axis["n_jackknife_folds"] == max(s["n_jackknife_folds"] for s in summaries)
    for field in ("n_features_evaluated", "n_features_lobo_stable", "n_single_unit_driven"):
        assert axis[field] == sum(s[field] for s in summaries)
    assert axis["lobo_stable_fraction"] == round(axis["n_features_lobo_stable"] / axis["n_features_evaluated"], 4)


@pytest.mark.parametrize("path", [
    "diffacc.pseudobulk", "differential_accessibility.pseudobulk",
    "bulk_diffacc", "differential_accessibility",
])
def test_chromatin_extraction_paths_and_insufficient_replicates(live_results, path):
    container = live_results["scatac_n2"]["pseudobulk"]
    for key in reversed(path.split(".")):
        container = {key: container}
    axis = _axis(_claim("chromatin", "differential_accessibility"),
                 {"chromatin_agent": {"findings": container}})
    assert axis == {"status": "not_testable_insufficient_replicates",
                    "biological_unit": "biological_unit", "minimum_required_per_arm": 3}


def test_certificate_manual_override_and_non_de_claims(live_results):
    results = {"bulk_rna_agent": live_results["bulk_rna"]}
    manual = {"status": "measured", "n_jackknife_folds": 99,
              "n_features_evaluated": 7, "n_features_lobo_stable": 5}
    assert _axis(_claim("bulk_rna"), results, {"lobo_jackknife": manual}) == {
        **manual, "biological_unit": "biological_unit"}
    # Explicit non-measured status overrides live results; empty/None placeholders
    # must not shadow measured agent-result LOBO summaries.
    axis_err = _axis(_claim("bulk_rna"), results, {"lobo_jackknife": {"status": "error"}})
    assert axis_err["status"] == "error"
    assert axis_err.get("error_details") == ["LOBO jackknife failed"]
    for empty_placeholder in ({}, None):
        axis_live = _axis(_claim("bulk_rna"), results, {"lobo_jackknife": empty_placeholder})
        assert axis_live["status"] == "measured"
    assert _axis(_claim("chromatin", "differential_accessibility"), results)["status"] == "eligible_not_executed"
    assert _axis(_claim("bulk_rna", "qc", "descriptive"), results)["status"] == "not_applicable_descriptive"
    assert _axis(_claim("bulk_rna", "peak_to_gene"), results)["status"] == "eligible_not_executed"


def test_certificate_mixed_statuses_only_sum_measured_contrasts(live_results):
    measured = live_results["bulk_rna"]["contrasts"][0]
    insufficient = live_results["scatac_n2"]["pseudobulk"]["comparisons"][0]
    failed = {"lobo_jackknife": {"status": "error", "details": "diagnostic failed"}}
    results = {"bulk_rna_agent": {"contrasts": [measured, insufficient, failed]}}
    axis = _axis(_claim("bulk_rna"), results)
    assert axis["status"] == "measured"
    assert axis["n_contrasts_measured"] == 1
    assert axis["n_features_evaluated"] == measured["lobo_jackknife"]["n_features_evaluated"]
    results["bulk_rna_agent"]["contrasts"] = [insufficient, failed]
    axis_mixed_err = _axis(_claim("bulk_rna"), results)
    assert axis_mixed_err["status"] == "error"
    assert axis_mixed_err.get("error_details") == ["diagnostic failed"]


def test_script_import_isolation():
    for name in SCRIPTS:
        tree = ast.parse((ROOT / "aria/scripts" / f"{name}.py").read_text())
        for node in ast.walk(tree):
            modules = ([node.module or ""] if isinstance(node, ast.ImportFrom) else
                       [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            assert not any(m == p or m.startswith(p + ".") for m in modules
                           for p in ("aria.agents", "litellm", "pydantic"))
            if isinstance(node, ast.ImportFrom) and any(a.name == "compute_counts_lobo_jackknife" for a in node.names):
                assert node.module == "aria.utils.lobo_jackknife"


def _live_worker(work):
    """Fit actual PyDESeq2 models and save only JSON IPC results."""
    os.environ["NUMBA_CACHE_DIR"] = str(work / "numba_cache")
    os.environ["MPLCONFIGDIR"] = str(work / "mpl_config")
    import contextlib
    import importlib
    import multiprocessing
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scipy import sparse

    # Pseudobulk's older PyDESeq2 API has no worker parameter in its facade.
    # Bound the detected CPU count in this subprocess; inference remains real.
    multiprocessing.cpu_count = lambda: 1
    modules = {name: importlib.import_module(f"aria.scripts.{name}") for name in SCRIPTS}
    rng = np.random.default_rng(849)
    units = [f"{cond}_{i}" for cond in ("ctrl", "treated") for i in range(4)]
    means = rng.uniform(400, 1400, 120)
    raw = rng.negative_binomial(150, 150 / (150 + means[:, None]), size=(120, 8))
    raw[:10, 4:] *= 8
    raw[1] = [500, 510, 490, 505, 500, 510, 490, 30000]
    genes = ["G_STABLE", "G_OUTLIER"] + [f"G{i}" for i in range(2, 120)]
    counts = pd.DataFrame(raw, index=genes, columns=units)
    meta = pd.DataFrame({"condition": ["ctrl"] * 4 + ["treated"] * 4,
                         "donor": [f"d{i}" for i in range(4)] * 2}, index=units)
    out = {}

    def write_matrix(matrix, metadata, prefix):
        count_path, meta_path = work / f"{prefix}_counts.tsv", work / f"{prefix}_meta.tsv"
        matrix.to_csv(count_path, sep="\t", index_label="peak_id" if prefix == "atac" else "gene_id")
        metadata.to_csv(meta_path, sep="\t", index_label="sample_id" if prefix == "atac" else "sample")
        return str(count_path), str(meta_path)

    def write_cells(matrix, metadata, name, groups=("A",)):
        rows, obs = [], []
        for group in groups:
            for j, uid in enumerate(matrix.columns):
                if group == "B" and uid.endswith("_3"):
                    continue
                # Exact integer partition: aggregation reconstructs donor counts.
                cell_counts = rng.multinomial(1, [0.5, 0.5], size=(matrix.shape[0],))
                first = matrix.iloc[:, j].to_numpy() // 2
                # Both cells are raw integer counts; uneven partition avoids ties.
                first = np.maximum(0, first - cell_counts[:, 0])
                for k, row in enumerate((first, matrix.iloc[:, j].to_numpy() - first)):
                    rows.append(row)
                    obs.append({**metadata.loc[uid].to_dict(), "cell_type": group,
                                "leiden": "0", "cell": f"{group}_{uid}_{k}"})
        adata = ad.AnnData(sparse.csr_matrix(np.asarray(rows)),
                          obs=pd.DataFrame(obs).set_index("cell"),
                          var=pd.DataFrame(index=matrix.index))
        path = work / f"{name}.h5ad"
        adata.write_h5ad(path)
        return str(path)

    with (work / "fits.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        # Two technical libraries per biological sample must still yield 8 folds.
        libraries, library_meta = {}, []
        for uid in units:
            first = counts[uid] // 2
            for k, values in enumerate((first, counts[uid] - first)):
                sid = f"{uid}_lane{k}"
                libraries[sid] = values
                library_meta.append({"sample": sid, "condition": meta.at[uid, "condition"], "bio_unit": uid})
        count_path, meta_path = write_matrix(pd.DataFrame(libraries), pd.DataFrame(library_meta).set_index("sample"), "rna")
        bulk_params = {"files": [count_path], "metadata_file": meta_path,
                       "design_factor": "condition", "technical_replicate_col": "bio_unit",
                       "comparison": {"numerator": "treated", "denominator": "ctrl"},
                       "output_dir": str(work / "rna"), "run_pathways": False,
                       "lfc_threshold": 0.0, "lfc_shrink": False, "n_cpus": 1}
        bulk = modules["rna_bulk_de"].bulk_rna_de
        out["bulk_rna"] = bulk(bulk_params)
        fallback_counts = pd.DataFrame(libraries).iloc[:40].copy()
        treated_libraries = [s for s in fallback_counts.columns if s.startswith("treated_")]
        fallback_counts.loc[genes[2:10], treated_libraries] //= 8
        fallback_counts.loc["G_OUTLIER", ["treated_3_lane0", "treated_3_lane1"]] = 750
        fallback_path, _ = write_matrix(fallback_counts,
                                       pd.DataFrame(library_meta).set_index("sample"), "rna_small")
        out["bulk_rna_no_hits"] = bulk({**bulk_params, "files": [fallback_path],
                                       "padj_threshold": 0.0, "output_dir": str(work / "rna_no_hits")})

        pb_params = {"data_path": write_cells(counts, meta, "rna_cells", ("A", "B")),
                     "groupby": "cell_type", "condition_col": "condition", "replicate_col": "donor",
                     "comparisons": [["treated", "ctrl"], ["ctrl", "treated"]],
                     "min_cells_per_pseudosample": 2, "min_replicates_per_condition": 3,
                     "lfc_shrink": False, "lfc_min": 0.0, "fdr_strategy": "global",
                     "use_raw": False, "output_dir": str(work / "scrna")}
        out["scrna"] = modules["rna_pseudobulk_de"].rna_pseudobulk_de(pb_params)

        peaks = counts.copy()
        peaks.index = [f"chr1:{i * 200}-{i * 200 + 100}" for i in range(120)]
        sc_params = {}
        for key, paired, n in (("scatac_unpaired", False, 4), ("scatac_paired", True, 4), ("scatac_n2", False, 2)):
            selected = [u for u in units if int(u.rsplit("_", 1)[1]) < n]
            sc_meta = meta.loc[selected].copy()
            if not paired:
                sc_meta["donor"] = selected
            params = {"data_path": write_cells(peaks[selected], sc_meta, key),
                      "groupby": "leiden", "condition_col": "condition", "replicate_col": "donor",
                      "comparisons": [["treated", "ctrl"]], "skip_per_cluster": True,
                      "min_cells_per_pseudosample": 2, "min_replicates_per_condition": 2,
                      "lfc_min": 0.0, "n_cpus": 1, "output_dir": str(work / key)}
            out[key] = modules["chromatin_diffacc"].chromatin_diffacc(params)
            sc_params[key] = params

        peak_meta = meta.rename(columns={"donor": "replicate"})
        count_path, meta_path = write_matrix(peaks, peak_meta, "atac")
        # Heterologous exposure varies across libraries, forcing external SF use.
        exposures = dict(zip(units, [900, 1000, 1100, 1200, 950, 1050, 1150, 1250]))
        atac_params = {"counts_matrix_path": count_path, "sample_metadata_path": meta_path,
                       "comparisons": [["treated", "ctrl"]], "spikein_counts": exposures,
                       "lfc_min": 0.0, "n_cpus": 1, "top_n": 200, "output_dir": str(work / "bulk_atac")}
        out["bulk_atac"] = modules["chromatin_bulk_diffacc"].chromatin_bulk_diffacc(atac_params)
        from aria.utils.lobo_jackknife import compute_counts_lobo_jackknife
        rep_ids = [f"{peak_meta.at[u, 'condition']}::{peak_meta.at[u, 'replicate']}" for u in units]
        rep_counts, rep_meta = peaks.copy(), peak_meta.copy()
        rep_counts.columns, rep_meta.index = rep_ids, rep_ids
        sf = np.array(list(exposures.values()), dtype=float)
        sf /= np.exp(np.log(sf).mean())
        hits = [p["peak"] for p in out["bulk_atac"]["comparisons"][0]["top_peaks"]]
        out["bulk_atac_expected_lobo"] = compute_counts_lobo_jackknife(
            rep_counts, rep_meta, "condition", "treated", "ctrl", size_factors=sf,
            target_features=hits or list(peaks.index[:50]), primary_significant_features=hits,
        )

        def fail_lobo(*args, **kwargs):
            raise RuntimeError("x" * 250)

        out["errors"] = {}
        for name, params in (("rna_bulk_de", bulk_params), ("rna_pseudobulk_de", pb_params),
                             ("chromatin_diffacc", sc_params["scatac_unpaired"]),
                             ("chromatin_bulk_diffacc", atac_params)):
            modules[name].compute_counts_lobo_jackknife = fail_lobo
            func = bulk if name == "rna_bulk_de" else getattr(modules[name], name)
            result = func(params)
            if name == "rna_bulk_de":
                comp = result["contrasts"][0]
            elif name == "rna_pseudobulk_de":
                comp = result["per_group"]["A"]["per_comparison"]["treated_vs_ctrl"]
            elif name == "chromatin_diffacc":
                comp = result["pseudobulk"]["comparisons"][0]
            else:
                comp = result["comparisons"][0]
            out["errors"][name] = comp

    assert not any(n == "pydantic" or n == "litellm" or n.startswith("aria.agents") for n in sys.modules)
    (work / "results.json").write_text(json.dumps(out))


if __name__ == "__main__":
    _live_worker(Path(sys.argv[1]))
