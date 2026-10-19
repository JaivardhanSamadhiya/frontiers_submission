"""Read-only replay of compact comparator metrics, not a fitted-model audit.

The full training/representation audit requires the locally retained models and
protein caches. This separate replay needs only the compact published tables.
It never fits models, alters predictions, or substitutes for that full audit.
"""
if not __debug__:
    raise RuntimeError('Run without Python optimization; verification assertions are mandatory')

from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score, matthews_corrcoef

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'data/results_v2/frontiers_strengthening_2026_10_16/comparator/official_full_budget_resumed'
sys.path.insert(0, str(ROOT / 'src'))
from precisionphage.benchmarks.screening import cluster_mean_ci, screening_curve


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return pd.read_csv(path, float_precision='round_trip')


def verify_development_selection(work):
    """Reconstruct official median/MCC cutoff choice using development only."""
    meta = json.loads((work / 'inputs.json').read_text())
    lock = json.loads((work / 'prediction_lock.json').read_text())
    predictions = read(work / 'development_predictions.csv')
    assert set(predictions.strain) <= set(meta['training_hosts'])
    assert not set(predictions.strain) & set(meta['test_hosts'])
    eligible = {int(path.stem.split('_')[-1]) for path in (work / 'tables').glob('*.csv')}
    observed_cutoffs = {int(str(cutoff).split('_')[-1]) for cutoff in predictions.cut_off.unique()}
    assert observed_cutoffs == eligible
    for _, subset in predictions.groupby('cut_off'):
        assert set(subset.run) == {f'run_{repeat}' for repeat in range(50)}
    medians = predictions.groupby(['cut_off', 'strain', 'phage', 'interaction']).Confidence.median().reset_index()
    scores = []
    for cutoff, subset in medians.groupby('cut_off'):
        scores.append((int(str(cutoff).split('_')[-1]),
                       matthews_corrcoef(subset.interaction, subset.Confidence.gt(.5))))
    chosen = sorted(scores, key=lambda item: (-item[1], item[0]))[0][0]
    assert chosen == int(lock['chosen_cutoff'])
    saved = read(work / 'development_summary/model_performance/model_performance_metrics.csv')
    saved = saved.assign(cutoff=saved.cut_off.str.split('_').str[-1].astype(int)).set_index('cutoff')
    assert set(saved.index) == {cutoff for cutoff, _ in scores}
    for cutoff, mcc in scores:
        np.testing.assert_allclose(mcc, saved.loc[cutoff, 'MCC'], rtol=0, atol=1e-12)
    assert int(saved.index[0]) == chosen
    expected_grid = {(iterations, rate, depth) for iterations in [500, 1000]
                     for rate in [.05, .1] for depth in [4, 6]}
    searches = 0
    for cutoff, _ in scores:
        for repeat in range(50):
            table = read(work / f'modeling/cutoff_{cutoff}/run_{repeat}/model_performance.csv')
            assert len(table) == 8 and np.isfinite(table.mcc).all()
            assert table.loss_function.eq('Logloss').all() and table.thread_count.eq(4).all()
            assert set(zip(table.iterations, table.learning_rate, table.depth)) == expected_grid
            searches += 1
    return {'fold': int(meta['fold']), 'chosen_cutoff_recomputed': chosen,
            'eligible_cutoffs': len(scores), 'eight_candidate_searches': searches}


def verify():
    audit = json.loads((RUN / 'audit.json').read_text())
    assert audit['passed'] and audit['folds'] == 5 and audit['selected_models_reloaded'] == 250
    for name, expected in audit['results_sha256'].items():
        assert Path(name).name == name and sha(RUN / name) == expected
    frame = read(RUN / 'verified_predictions.csv')
    assert len(frame) == 437 and not frame.duplicated(['strain', 'phage']).any()
    assert frame.strain.nunique() == 23 and frame.phage.nunique() == 19 and frame.group.nunique() == 17
    assert set(frame.fold) == set(range(5)) and frame.groupby('strain').phage.nunique().eq(19).all()
    np.testing.assert_array_equal(frame.interaction, frame.label)
    np.testing.assert_array_equal(frame.fold, frame.old_fold)
    assert set(frame.label) == {0, 1}
    development_checks = []
    for fold in range(5):
        work = RUN / f'0.001_fold{fold}'
        meta = json.loads((work / 'inputs.json').read_text())
        name = Path(meta['development_input'].replace('\\', '/'))
        assert not name.is_absolute() and '..' not in name.parts
        development = ROOT / name
        assert sha(development) == meta['development_input_sha256']
        assert sha(work / 'outer_query_features.csv') == meta['outer_query_sha256']
        assert set(read(development).strain) == set(meta['training_hosts'])
        assert not set(meta['training_hosts']) & set(meta['test_hosts'])
        observed = frame.loc[frame.fold.eq(fold)].set_index(['strain', 'phage']).sort_index()
        scores = read(work / 'outer_predictions.csv').set_index(['strain', 'phage']).sort_index()
        assert set(observed.index.get_level_values('strain')) == set(meta['test_hosts'])
        assert observed.index.equals(scores.index)
        np.testing.assert_allclose(observed.genophi, scores.Confidence, rtol=0, atol=1e-12)
        development_checks.append(verify_development_selection(work))
    methods = ['genophi', 'neighbour', 'prevalence']
    truth = frame.label.to_numpy()
    groups, inverse = np.unique(frame.group, return_inverse=True)
    point = {m: np.array([roc_auc_score(truth, frame[m]), average_precision_score(truth, frame[m])]) for m in methods}
    rng = np.random.default_rng(20261005)
    bootstrap = {m: [] for m in methods}
    skipped = 0
    for _ in range(2000):
        counts = np.bincount(rng.integers(len(groups), size=len(groups)), minlength=len(groups))
        weights = counts[inverse]
        if len(np.unique(truth[weights > 0])) < 2:
            skipped += 1
            continue
        for method in methods:
            bootstrap[method].append([roc_auc_score(truth, frame[method], sample_weight=weights),
                                      average_precision_score(truth, frame[method], sample_weight=weights)])
    rows = []
    for method in methods:
        bounds = np.quantile(bootstrap[method], [.025, .975], axis=0)
        for index, metric in enumerate(['auroc', 'average_precision']):
            rows.append((method, metric, point[method][index], bounds[0, index], bounds[1, index]))
    for method in methods[1:]:
        bounds = np.quantile(np.asarray(bootstrap['genophi']) - np.asarray(bootstrap[method]), [.025, .975], axis=0)
        for index, metric in enumerate(['auroc', 'average_precision']):
            rows.append(('genophi_minus_' + method, metric, point['genophi'][index] - point[method][index],
                         bounds[0, index], bounds[1, index]))
    host_rows = []
    for host, subset in frame.groupby('strain'):
        record = {'host': host, 'group': subset.group.iloc[0]}
        for method in methods:
            record[method] = screening_curve(subset.label, subset[method], subset.phage)['hit'][2]
        host_rows.append(record)
    hosts = pd.DataFrame(host_rows)
    saved_hosts = read(RUN / 'screening_hosts.csv')
    assert hosts.host.tolist() == saved_hosts.host.tolist()
    assert hosts.group.tolist() == saved_hosts.group.tolist()
    np.testing.assert_array_equal(hosts[methods], saved_hosts[methods])
    for method in methods:
        ci = cluster_mean_ci(hosts[method], hosts.group)
        rows.append((method, 'three_phage_union_hit', ci['mean'], ci['lo'], ci['hi']))
    for method in methods[1:]:
        ci = cluster_mean_ci(hosts.genophi - hosts[method], hosts.group)
        rows.append(('genophi_minus_' + method, 'three_phage_union_hit', ci['mean'], ci['lo'], ci['hi']))
    reproduced = pd.DataFrame(rows, columns=['comparison', 'metric', 'estimate', 'lo', 'hi']).set_index(['comparison', 'metric']).sort_index()
    saved = read(RUN / 'matched_metrics.csv').set_index(['comparison', 'metric']).sort_index()
    assert reproduced.index.equals(saved.index)
    np.testing.assert_allclose(reproduced, saved[['estimate', 'lo', 'hi']], rtol=0, atol=1e-10)
    assert skipped == audit['one_class_replicates_skipped']
    return {'passed': True, 'scope': 'Compact metric, paired interval, input-hash and alignment replay; no model reload or protein-representation reconstruction.',
            'pairs': len(frame), 'hosts': 23, 'phages': 19, 'host_groups': 17,
            'prevalence': float(frame.label.mean()), 'metrics_replayed': len(saved),
            'bootstrap_attempts': 2000, 'undefined_bootstrap_draws': skipped,
            'development_selection_checks': development_checks,
            'predictive_model_fits_started': 0, 'files_modified': 0}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
