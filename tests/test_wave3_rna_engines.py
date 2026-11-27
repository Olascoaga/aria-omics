"""Wave 3B scientific regressions; external model outputs are controlled fixtures."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

MOUSE_TOKENS = ['mouse', 'mm', 'mm10', 'mm39', 'musculus', 'Mus musculus']


def _adata(x, **obs):
    x = np.asarray(x, dtype=float)
    return ad.AnnData(x, obs=pd.DataFrame(obs, index=[f'c{i}' for i in range(len(x))]),
                      var=pd.DataFrame(index=[f'g{i}' for i in range(x.shape[1])]))


def _input(monkeypatch, tmp_path, a):
    import aria.utils.safe_h5ad as io
    path = tmp_path / 'input.h5ad'
    path.touch()
    monkeypatch.setattr(io, 'read_h5ad', lambda p: a.copy())
    return str(path)


def _strict(result):
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize('failure', ['nan_p', 'inf_p', 'exception_after_coef', 'nan_coef', 'none'])
def test_da_fisher_is_excluded_from_donor_fdr(monkeypatch, tmp_path, failure):
    import statsmodels.api as sm
    from aria.scripts.rna_diff_abundance import rna_diff_abundance
    a = _adata(np.ones((8, 2)), cell=['A', 'B'] * 4,
               condition=['test'] * 4 + ['ref'] * 4,
               donor=['t1'] * 2 + ['t2'] * 2 + ['r1'] * 2 + ['r2'] * 2)

    class Fit:
        # Negative adjusted effect despite equal raw proportions.
        params = {'is_test': np.nan if failure == 'nan_coef' else -2.0}

        @property
        def pvalues(self):
            if failure == 'exception_after_coef':
                raise ValueError('singular HC3')
            return {'is_test': {'nan_p': np.nan, 'inf_p': np.inf}.get(failure, 0.001)}

    monkeypatch.setattr(sm, 'OLS', lambda **kw: SimpleNamespace(fit=lambda **kw: Fit()))
    result = rna_diff_abundance(dict(data_path=_input(monkeypatch, tmp_path, a),
        groupby='cell', condition_col='condition', replicate_col='donor',
        comparisons=[['test', 'ref']]))
    comp = result['per_comparison']['test_vs_ref']
    for row in comp['per_cell_type']:
        assert row['dispersion'] is None
        if failure != 'none':
            assert row['model'] == 'fisher_exact_fallback'
            assert not row['fdr_included'] and row['padj'] is None
            assert row['clr_log_ratio'] is None and not row['significant']
        else:
            assert row['fdr_included'] and row['direction'] == 'down'
    assert comp['degraded'] is (failure != 'none')
    _strict(result)


@pytest.mark.parametrize('organism', MOUSE_TOKENS)
def test_celltypist_mouse_catalog_and_immune_rejection(organism):
    from aria.scripts.rna_celltypist import _resolve_celltypist_model, rna_celltypist
    assert _resolve_celltypist_model(None, organism, 'brain')[:2] == (
        'Mouse_Whole_Brain.pkl', 'tissue_hint')
    model, source, warning = _resolve_celltypist_model(None, organism, 'immune')
    assert source == 'default_immune_fallback' and 'HUMAN' in warning
    assert 'Isocortex' not in model
    result = rna_celltypist(dict(data_path='unused', organism=organism, tissue_hint='immune'))
    assert result['status'] == 'error' and result['requires_ack']


@pytest.mark.parametrize('raw_extra_genes', [False, True])
def test_celltypist_preserves_original_layers(monkeypatch, tmp_path, raw_extra_genes):
    pytest.importorskip('scanpy')
    import aria.scripts.rna_celltypist as engine
    full = _adata(np.log1p(np.full((4, 4), 2500.0)), leiden=['0', '0', '1', '1'])
    full.raw = full.copy()
    a = full[:, :2].copy() if raw_extra_genes else full.copy()
    a.X[:] = -1  # Annotation must use the full log-normalized raw matrix.
    for name in ['counts', 'spliced', 'unspliced']:
        a.layers[name] = sparse.csr_matrix(np.full(a.shape, 3))
    path = _input(monkeypatch, tmp_path, a)
    seen = {}

    def annotate(matrix, **kwargs):
        seen['n_vars'] = matrix.n_vars
        assert matrix.raw is None and matrix.X.min() >= 0
        return SimpleNamespace(predicted_labels=pd.DataFrame(
            {'predicted_labels': ['T'] * 4}, index=matrix.obs_names),
            probability_matrix=pd.DataFrame({'T': [0.9] * 4}, index=matrix.obs_names))

    monkeypatch.setitem(sys.modules, 'celltypist', SimpleNamespace(
        annotate=annotate, models=SimpleNamespace()))
    monkeypatch.setattr(engine, '_celltypist_model_cached', lambda m: True)
    result = engine.rna_celltypist(dict(data_path=path, model='explicit.pkl',
                                      output_dir=str(tmp_path)))
    assert result['status'] == 'success' and seen['n_vars'] == 4
    saved = ad.read_h5ad(result['output_path'])
    assert list(saved.var_names) == list(a.var_names)
    assert list(saved.obs['cell_type_celltypist']) == ['T'] * 4
    for name in a.layers:
        np.testing.assert_array_equal(saved.layers[name].toarray(), a.layers[name].toarray())


@pytest.mark.parametrize('organism', MOUSE_TOKENS + ['human'])
@pytest.mark.parametrize('has_raw', [True, False])
def test_cellcomm_species_raw_zero_pvalues_and_full_counts(monkeypatch, tmp_path, organism, has_raw):
    from aria.scripts.rna_cellcomm import rna_cellcomm
    a = _adata(np.ones((4, 2)), cell_type=['A', 'A', 'B', 'B'])
    if has_raw:
        a.raw = a.copy()
    n = 80
    frame = pd.DataFrame(dict(source=['A'] * 50 + ['B'] * 30,
        target=['B'] * 50 + ['A'] * 30, ligand_complex=[f'L{i}' for i in range(n)],
        receptor_complex=['R'] * n, specificity_rank=np.linspace(0.001, 0.1, n),
        cellphone_pvals=[0.0, np.nan, np.inf, True] + [0.2] * (n - 4)))
    frame.loc[n] = ['A', 'A', 'self', 'self', 0.0, 0.0]
    seen = {}

    def rank(matrix, **kwargs):
        seen.update(kwargs)
        matrix.uns['liana_res'] = frame

    monkeypatch.setitem(sys.modules, 'liana', SimpleNamespace(mt=SimpleNamespace(rank_aggregate=rank)))
    result = rna_cellcomm(dict(data_path=_input(monkeypatch, tmp_path, a),
                              organism=organism, output_dir=str(tmp_path)))
    assert seen['resource_name'] == ('consensus' if organism == 'human' else 'mouseconsensus')
    assert seen['use_raw'] is has_raw
    assert result['n_interactions'] == result['n_interactions_total'] == n
    assert result['n_significant'] == int((frame.iloc[:n].specificity_rank < 0.05).sum())
    assert result['n_autocrine_dropped'] == 1
    assert result['top_interactions'][0]['cellphone_pval'] == 0.0
    assert all(r['cellphone_pval'] is None for r in result['top_interactions'][1:4])
    assert len(result['top_interactions']) == 20
    assert result['top_pairs'] == ['A→B', 'B→A']
    assert len(pd.read_csv(result['output_path'])) == n
    _strict(result)


def test_trajectory_disconnects_and_explicit_paga_groups(monkeypatch, tmp_path):
    sc = pytest.importorskip('scanpy')
    from aria.scripts.rna_trajectory import rna_trajectory
    a = _adata(np.ones((5, 2)), leiden=pd.Categorical(['0', '0', '1', '1', '1']),
               annotation=pd.Categorical(['stem', 'stem', 'T', 'B', 'B']))
    a.uns['neighbors'] = {}
    a.obsm['X_umap'] = np.zeros((5, 2))
    seen = {}

    def paga(matrix, groups):
        seen['groups'] = groups
        matrix.uns['paga'] = {'connectivities': sparse.eye(3, format='csr')}

    def dpt(matrix):
        matrix.obs['dpt_pseudotime'] = [0.0, 0.4, np.inf, -np.inf, np.nan]

    monkeypatch.setattr(sc.tl, 'paga', paga)
    monkeypatch.setattr(sc.tl, 'diffmap', lambda *args, **kw: None)
    monkeypatch.setattr(sc.tl, 'dpt', dpt)
    result = rna_trajectory(dict(data_path=_input(monkeypatch, tmp_path, a),
        cell_type_col='annotation', root_cell_type='stem', output_dir=str(tmp_path)))
    assert seen['groups'] == result['groupby'] == 'annotation'
    pt = result['pseudotime']
    assert pt['computed'] and pt['n_disconnected_cells'] == 3
    assert pt['pseudotime_by_group'] == {'stem': 0.2, 'T': None, 'B': None}
    saved = ad.read_h5ad(result['output_path'])
    assert not np.isinf(saved.obs.dpt_pseudotime).any()
    _strict(result)


def _rgg(groupby='leiden', use_raw=False):
    def rec(values, dtype):
        return np.rec.fromarrays([values, values], names=['0', '1'], formats=[dtype, dtype])
    return dict(params=dict(groupby=groupby, use_raw=use_raw),
        names=rec(['g0', 'g1', 'g2', 'g3', 'g4', 'g5'], 'U2'),
        logfoldchanges=rec([2, np.inf, -np.inf, np.nan, 2, 2], 'f8'),
        pvals_adj=rec([0.001, 0.001, 0.001, 0.001, np.inf, np.nan], 'f8'))


@pytest.mark.parametrize('sparse_x', [False, True])
@pytest.mark.parametrize('cached_group,cached_raw', [('leiden', False), ('other', False), ('leiden', True)])
def test_cluster_de_finite_scaled_and_cache_compatibility(monkeypatch, tmp_path, sparse_x, cached_group, cached_raw):
    sc = pytest.importorskip('scanpy')
    from aria.scripts.rna_de_per_cluster import rna_de_per_cluster
    a = _adata(-np.ones((4, 6)), leiden=['0', '0', '1', '1'])
    if sparse_x:
        a.X = sparse.csr_matrix(a.X)
    a.uns['rank_genes_groups'] = _rgg(cached_group, cached_raw)
    calls = []

    def rank(matrix, **kwargs):
        calls.append(kwargs)
        matrix.uns['aria_de'] = _rgg()

    monkeypatch.setattr(sc.tl, 'rank_genes_groups', rank)
    result = rna_de_per_cluster(dict(data_path=_input(monkeypatch, tmp_path, a),
                                    min_pct_in=0.9, output_dir=str(tmp_path)))
    assert bool(calls) is (cached_group != 'leiden' or cached_raw)
    assert result['n_significant_total'] == 2
    for records in result['de_genes_by_cluster'].values():
        assert [r['gene'] for r in records] == ['g0']
        assert records[0]['pct_in'] is None and records[0]['pct_out'] is None
    assert result['warnings']
    _strict(result)


@pytest.mark.parametrize('values', [[np.nan], [np.inf], [-np.inf], [1e308, 1e308], []])
def test_qc_mad_nonfinite_never_rounds_nan_or_inf(values):
    from aria.scripts.rna_qc import _mad_bounds
    bounds = _mad_bounds(values, log_space=True)
    assert len(bounds) == 2
    assert all(np.isfinite(v) or np.isnan(v) for v in bounds)


def test_qc_none_overrides_and_cache_bounds(monkeypatch, tmp_path):
    from aria.scripts.rna_qc import _cache_params, _cache_matches, rna_qc
    nulls = {k: None for k in ['mt_threshold', 'min_genes', 'min_cells', 'min_counts',
                              'max_counts', 'max_genes', 'initial_min_genes', 'initial_min_cells']}
    assert _cache_params(nulls) == _cache_params({})
    for key in ['min_counts', 'max_counts', 'max_genes']:
        assert _cache_params({key: 500}) != _cache_params({})
    assert not _cache_matches({'cache_version': 6, 'cache_params': _cache_params({})}, _cache_params({}))
    pytest.importorskip('scanpy')
    a = _adata(np.ones((4, 3)), total_counts=[1000] * 4,
               n_genes_by_counts=[300] * 4, pct_counts_mt=[3] * 4)
    result = rna_qc(dict(data_path=_input(monkeypatch, tmp_path, a),
                        run_scrublet=False, output_dir=str(tmp_path), **nulls))
    assert result['status'] == 'success' and result['n_cells_after'] == 4
    assert result['mt_threshold_used'] is not None and result['min_genes_used'] >= 200


@pytest.mark.parametrize('organism', MOUSE_TOKENS)
@pytest.mark.parametrize('prefix', ['mt-', 'MT-'])
def test_qc_mouse_mitochondrial_prefix_and_case_fallback(monkeypatch, tmp_path, organism, prefix):
    pytest.importorskip('scanpy')
    from aria.scripts.rna_qc import rna_qc
    a = _adata(np.tile([1, 4, 5], (6, 1)))
    a.var_names = [prefix + 'Nd1', 'G1', 'G2']
    result = rna_qc(dict(data_path=_input(monkeypatch, tmp_path, a), organism=organism,
        initial_min_genes=1, initial_min_cells=1, min_genes=1, min_cells=1,
        mt_threshold=25, run_scrublet=False, output_dir=str(tmp_path)))
    assert result['status'] == 'success'
    assert result['mt_stats']['mean'] == 10.0
    if prefix == 'MT-':
        assert any('case-insensitive' in w for w in result['warnings'])


def _mock_scanpy_pipeline(monkeypatch):
    sc = pytest.importorskip('scanpy')
    seen = {}
    for name in ['normalize_total', 'log1p', 'highly_variable_genes', 'scale']:
        monkeypatch.setattr(sc.pp, name, lambda *a, **kw: None)

    def pca(matrix, **kwargs):
        seen['pca'] = kwargs['n_comps']
        assert 1 <= kwargs['n_comps'] < min(matrix.shape)
        matrix.obsm['X_pca'] = np.arange(matrix.n_obs * kwargs['n_comps']).reshape(matrix.n_obs, -1).astype(float)

    def neighbors(matrix, **kwargs):
        seen['neighbors'] = kwargs['n_pcs']
        assert 1 <= kwargs['n_pcs'] <= matrix.obsm[kwargs['use_rep']].shape[1]
        matrix.uns['neighbors'] = {}

    def leiden(matrix, **kw):
        matrix.obs['leiden'] = ['0', '0'] + ['1'] * (matrix.n_obs - 2)

    def rank(matrix, **kw):
        matrix.uns['rank_genes_groups'] = {'names': np.rec.fromarrays([['g0'], ['g1']], names=['0', '1'])}

    monkeypatch.setattr(sc.tl, 'pca', pca)
    monkeypatch.setattr(sc.pp, 'neighbors', neighbors)
    monkeypatch.setattr(sc.tl, 'umap', lambda matrix, **kw: matrix.obsm.__setitem__('X_umap', np.zeros((matrix.n_obs, 2))))
    monkeypatch.setattr(sc.tl, 'leiden', leiden)
    monkeypatch.setattr(sc.tl, 'rank_genes_groups', rank)
    monkeypatch.setitem(sys.modules, 'harmonypy', SimpleNamespace(
        run_harmony=lambda embedding, *a, **kw: SimpleNamespace(Z_corr=embedding)))
    return seen


@pytest.mark.parametrize('engine_name', ['rna_integration', 'rna_clustering'])
@pytest.mark.parametrize('shape', [(3, 5), (5, 3), (1, 4), (4, 1)])
def test_small_pca_dimensions(monkeypatch, tmp_path, engine_name, shape):
    import importlib
    engine = getattr(importlib.import_module('aria.scripts.' + engine_name), engine_name)
    a = _adata(np.ones(shape), batch=['a'] + ['b'] * (shape[0] - 1))
    seen = _mock_scanpy_pipeline(monkeypatch)
    result = engine(dict(data_path=_input(monkeypatch, tmp_path, a), resolution=0.5,
                         output_dir=str(tmp_path)))
    if min(shape) < 2:
        assert result['status'] == 'skipped'
    else:
        assert result['status'] == 'success'
        assert seen['pca'] == seen['neighbors'] == min(shape) - 1


@pytest.mark.parametrize('singleton,raise_score', [(True, False), (False, True), (False, False)])
def test_integration_precomputed_pcs_and_degenerate_silhouette(monkeypatch, tmp_path, singleton, raise_score):
    import sklearn.metrics
    from aria.scripts.rna_integration import rna_integration
    a = _adata(np.ones((6, 4)), batch=['a'] * (1 if singleton else 3) + ['b'] * (5 if singleton else 3))
    a.obsm['X_pca'] = np.arange(12).reshape(6, 2).astype(float)
    seen = _mock_scanpy_pipeline(monkeypatch)

    def score(*args, **kwargs):
        if singleton:
            pytest.fail('singleton labels must be guarded')
        if raise_score:
            raise ValueError('random sample contains only one label')
        return 0.2

    monkeypatch.setattr(sklearn.metrics, 'silhouette_score', score)
    result = rna_integration(dict(data_path=_input(monkeypatch, tmp_path, a), output_dir=str(tmp_path)))
    assert seen['neighbors'] == 2
    assert result['silhouette_before'] == (None if singleton or raise_score else 0.2)
    assert result['silhouette_after'] == result['silhouette_before']
    _strict(result)


@pytest.mark.parametrize('exc', [AttributeError, ValueError, RuntimeError])
@pytest.mark.parametrize('where', ['constructor', 'fit_size_factors', 'vst_fit', 'vst_transform'])
def test_vst_expected_failures_fall_back(monkeypatch, exc, where):
    from aria.scripts.rna_bulk.transforms import _run_vst

    class DDS:
        def __init__(self, **kw):
            if where == 'constructor':
                raise exc('small dataset')
            self.layers = {'normed_counts': kw['counts'].to_numpy()}

        def fit_size_factors(self):
            if where == 'fit_size_factors':
                raise exc('zero library')

        def vst_fit(self, **kw):
            if where == 'vst_fit':
                raise exc('no dispersion')

        def vst_transform(self):
            raise exc('unsupported transform')

        def deseq2(self):
            pass

    monkeypatch.setitem(sys.modules, 'pydeseq2.dds', SimpleNamespace(DeseqDataSet=DDS))
    counts = pd.DataFrame({'empty': [0, 0], 'sample': [10, 20]})
    result = _run_vst(counts, pd.DataFrame(index=counts.columns), [])
    assert result.shape == counts.shape and np.isfinite(result.to_numpy()).all()
    assert (result['empty'] == 0).all()


@pytest.mark.parametrize('all_zero', [False, True])
def test_vst_missing_dependency_zero_libraries(monkeypatch, all_zero):
    from aria.scripts.rna_bulk.transforms import _run_vst
    monkeypatch.setitem(sys.modules, 'pydeseq2.dds', None)
    counts = pd.DataFrame({'zero': [0, 0], 'other': [0, 0] if all_zero else [10, 20]})
    result = _run_vst(counts, pd.DataFrame(index=counts.columns), [])
    assert np.isfinite(result.to_numpy()).all() and (result['zero'] == 0).all()
    if not all_zero:
        np.testing.assert_allclose(result['other'], np.log2([11, 21]))


def test_tpm_invalid_lengths_empty_matches_and_zero_scaling():
    from aria.scripts.rna_bulk.transforms import _compute_tpm, _select_variable_genes
    counts = pd.DataFrame({'empty': [0] * 5, 'sample': [10] * 5}, index=list('abcde'))
    assert _compute_tpm(counts, {'other': 1000}, []) is None
    assert _compute_tpm(counts, {'a': 0, 'b': -1, 'c': np.inf}, []) is None
    result = _compute_tpm(counts, dict(a=1000, b=2000, c=0, d=-1, e=np.nan), [])
    assert list(result.index) == ['a', 'b']
    assert result['sample'].sum() == pytest.approx(1e6)
    assert (result['empty'] == 0).all() and np.isfinite(result.to_numpy()).all()
    _, n_pc, _ = _select_variable_genes(counts, biotype_map={'a': 'protein_coding'})
    assert n_pc == 0  # No protein-coding restriction actually applied.


@pytest.mark.parametrize('ext,sep', [('.csv', ','), ('.tsv', '\t'), ('.txt', '\t'),
    ('.counts', '\t'), ('.featureCounts', '\t'), ('.htseq', '\t'), ('.unknown', '\t')])
def test_bulk_count_delimiters_and_numeric_sample_ids(tmp_path, ext, sep):
    from aria.scripts.rna_bulk.counts import _load_counts, _load_or_infer_metadata
    frame = pd.DataFrame({'101': [10, 30, 40], '102': [20, 40, 50]}, index=['a', 'b', 'c'])
    path = tmp_path / ('counts' + ext)
    frame.to_csv(path, sep=sep)
    counts, warnings, provenance = _load_counts([str(path)])
    assert counts is not None, (warnings, provenance)
    pd.testing.assert_frame_equal(counts, frame)
    metadata_path = tmp_path / 'metadata.tsv'
    pd.DataFrame({'condition': ['A', 'B']}, index=[101, 102]).to_csv(metadata_path, sep='\t')
    meta, warnings = _load_or_infer_metadata(counts, str(metadata_path), 'condition')
    assert list(meta.index) == list(counts.columns) == ['101', '102']


@pytest.mark.parametrize('ext,sep', [('.csv', ','), ('.txt', '\t'), ('.counts', '\t')])
def test_bulk_multiple_per_sample_delimiters(tmp_path, ext, sep):
    from aria.scripts.rna_bulk.counts import _load_counts
    paths = []
    for sample in ['A', 'B']:
        path = tmp_path / (sample + ext)
        pd.Series([10, 30, 40], index=['a', 'b', 'c']).to_csv(path, sep=sep, header=False)
        paths.append(str(path))
    counts, warnings, provenance = _load_counts(paths)
    assert counts is not None, (warnings, provenance)
    assert counts.shape == (3, 2) and list(counts.columns) == ['A', 'B']


def test_bulk_qc_respects_design_factor_and_forwards_it(monkeypatch, tmp_path):
    import aria.scripts.rna_bulk.qc as engine
    counts = pd.DataFrame({f's{i}': [10 + i, 30 + i, 20 + i] for i in range(6)})
    metadata = pd.DataFrame({'donor_id': [f'd{i}' for i in range(6)],
        'condition': ['A'] * 3 + ['B'] * 3}, index=counts.columns)
    monkeypatch.setattr(engine, '_run_vst', lambda *args: np.log2(counts + 1))
    seen = {}

    def plot(matrix, meta, out, warnings, design_factor=None):
        seen['factor'] = design_factor
        return None, None, np.zeros((6, 2)), [0, 0]

    monkeypatch.setattr(engine, '_plot_pca_mds', plot)
    result = engine._sample_qc(counts, metadata, str(tmp_path), [], design_factor='condition')
    assert 'error' not in result
    assert seen['factor'] == 'condition'
    assert set(result['replicate_correlations']) == set(counts.columns)


def _fake_pseudobulk(monkeypatch):
    import aria.scripts.rna_pseudobulk_de as engine

    class DDS:
        def __init__(self, counts, metadata, **kwargs):
            self.counts = counts

        def deseq2(self):
            pass

    class Stats:
        def __init__(self, dds, **kw):
            self.results_df = pd.DataFrame({
                'pvalue': [0.001, 0.02, 0.003, 0.004, 0.005, 0.006],
                'padj': [np.nan, 0.03, 0.01, 0.01, -np.inf, 0.01],
                'log2FoldChange': [2, 0.1, np.inf, -np.inf, 3, -2],
                'baseMean': [np.inf, 10, 10, 10, 10, 10],
            }, index=dds.counts.columns)

        def summary(self):
            pass

    monkeypatch.setitem(sys.modules, 'pydeseq2.dds', SimpleNamespace(DeseqDataSet=DDS))
    monkeypatch.setitem(sys.modules, 'pydeseq2.ds', SimpleNamespace(DeseqStats=Stats))
    monkeypatch.setattr(engine, 'compute_counts_lobo_jackknife', lambda *a, **kw: {'status': 'skipped'})
    return engine


@pytest.mark.parametrize('strategy', ['global', 'per_cluster'])
def test_pseudobulk_block_composition_dispersion_finite_and_alpha(monkeypatch, tmp_path, strategy):
    engine = _fake_pseudobulk(monkeypatch)
    rows, obs = [], []
    # A varies across donors independently of condition; B remains constant.
    # Two nuisance types keep B's proportion fixed and A's proportion variable.
    for donor in range(6):
        n_a = [1, 2, 3, 1, 2, 3][donor]
        for cell_type, n in [('A', n_a), ('B', 2), ('Z', 8 - n_a)]:
            for _ in range(n):
                rows.append([10, 10, 10, 10, 10, [1, 20, 1, 20, 1, 20][donor]])
                obs.append((cell_type, 'test' if donor < 3 else 'ref', f'd{donor}'))
    a = _adata(rows)
    a.obs[['cell', 'condition', 'donor']] = pd.DataFrame(obs, index=a.obs_names)
    path = _input(monkeypatch, tmp_path, a)
    result = engine.rna_pseudobulk_de(dict(data_path=path, groupby='cell',
        condition_col='condition', replicate_col='donor', comparisons=[['test', 'ref']],
        min_cells_per_pseudosample=1, min_replicates_per_condition=3,
        composition_covariate=True, lfc_shrink=False, fdr_strategy=strategy,
        output_dir=str(tmp_path)))
    assert result['status'] == 'success'
    blocks = [g['per_comparison']['test_vs_ref'] for g in result['per_group'].values()]
    assert all(c['status'] == 'success' for c in blocks), blocks
    assert {c['corrected_for_composition'] for c in blocks} == {True, False}
    for group, group_result in result['per_group'].items():
        comp = group_result['per_comparison']['test_vs_ref']
        mask = a.obs.cell == group
        sums = pd.DataFrame(a.X[mask], index=a.obs.loc[mask, 'donor']).groupby(level=0).sum().T
        means = sums.mean(axis=1)
        dispersions = (sums.var(axis=1, ddof=1) - means) / means ** 2
        expected = dispersions[np.isfinite(dispersions) & (dispersions > 0)].median()
        axis = comp['robustness_multiverse']['axes']['composition_covariate']
        assert axis == ['included' if comp['corrected_for_composition'] else 'not_included']
        assert comp['dispersion_estimate'] == pytest.approx(expected)
        assert comp['effective_alpha_local'] == 0.02  # Small-LFC BH hit still sets alpha.
        assert comp['effective_alpha_primary'] == (comp['effective_alpha_global'] if strategy == 'global' else 0.02)
        assert all(np.isfinite(r['log2fc']) for r in comp['all_sig'])
        assert not {'g2', 'g3'} & {r['gene'] for r in comp['all_sig']}
        if strategy == 'global':
            record = next(r for r in comp['top_genes'] if r['gene'] == 'g0')
            assert record['padj_local'] is None and record['basemean'] is None
        else:
            assert {r['gene'] for r in comp['all_sig']} == {'g5'}
    _strict(result)
    full = pd.read_csv(result['output_csv'].replace('pseudobulk_de.csv', 'pseudobulk_de_all_tested.csv'))
    assert not np.isinf(full.select_dtypes(include='number').to_numpy()).any()


@pytest.mark.parametrize('bad_value', [np.nan, np.inf, 1e30])
def test_pseudobulk_recovery_rejects_unsampled_nonfinite_or_overflow(monkeypatch, tmp_path, bad_value):
    engine = _fake_pseudobulk(monkeypatch)
    a = _adata(np.log1p(np.full((12, 6), 10.0)), cell=['A'] * 12,
               condition=['test'] * 6 + ['ref'] * 6,
               donor=[f'd{i // 2}' for i in range(12)], total_counts=[10000.0] * 12)
    # Deliberately exclude the bad library from the classifier/recovery probe.
    a.obs.loc['c11', 'total_counts'] = bad_value
    monkeypatch.setattr(engine, 'sample_row_indices', lambda n: np.arange(10))
    monkeypatch.setattr(engine, 'classify_matrix', lambda *a, **kw: {
        'is_raw_counts': False, 'max': 3.0, 'kind': 'lognorm'})
    result = engine.rna_pseudobulk_de(dict(data_path=_input(monkeypatch, tmp_path, a),
        groupby='cell', condition_col='condition', replicate_col='donor',
        comparisons=[['test', 'ref']], min_cells_per_pseudosample=1,
        min_replicates_per_condition=3, lfc_shrink=False, output_dir=str(tmp_path)))
    assert result['status'] == 'error' and result['error_type'] == 'LognormRecoveryFailed'
    _strict(result)
