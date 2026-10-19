"""Independent audit and matched comparison after all five full-budget folds."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import zipfile
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'data/results_v2/genophi_comparison'
RUN = BASE/'official_modeling_full_budget_v1'
sys.path.insert(0, str(ROOT/'src'))
from precisionphage.benchmarks.screening import cluster_mean_ci, screening_curve
from precisionphage.benchmarks.genophi_retention import validate_repair_manifest, verify_model_selection


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check_sources():
    repair = validate_repair_manifest(ROOT, RUN)
    manifest = json.loads((RUN/'run_manifest.json').read_text())
    for name, checksum in manifest['sha256'].items():
        actual = sha(ROOT/name)
        if actual != checksum:
            if not (name.endswith('47_official_genophi_outer_comparison.py') and
                    actual == repair['current_files_sha256'][name] and
                    sha(RUN/'driver_before_prediction_api_fix.py') == checksum):
                raise ValueError(f'Unexpected source change {name}')
    official = json.loads((BASE/'official_source_manifest.json').read_text())
    for f in official['files']:
        if sha(ROOT/'data/external_validation/genophi_comparison/official_source'/f['path']) != f['sha256']:
            raise ValueError('Official source changed')
    for marker in RUN.glob('0.001_fold*/stages/*.json'):
        row = json.loads(marker.read_text())
        if row['exit_code'] != 0: raise ValueError('Unsuccessful stage checkpoint')
        for path, value in row['sha256'].items():
            if sha(ROOT/path) != value: raise ValueError(f'Altered checkpoint {path}')


def verify_representation(work, meta):
    query = pd.read_csv(work/'outer_query_features.csv')
    full = pd.read_csv(ROOT/meta['development_input'])
    all_pairs = pd.concat([full.drop(columns='interaction'), query], ignore_index=True)
    for axis, cache_kind, names, fit in [('strain', 'host', sorted(set(meta['training_hosts']+meta['test_hosts'])), meta['training_hosts']),
                                        ('phage', 'phage', sorted(meta['phages']), meta['phages'])]:
        candidates = []
        for p in (ROOT/'data/external_validation/genophi_comparison/feature_cache_v2').glob(f'{cache_kind}_*/complete.json'):
            spec = json.loads(p.read_text())
            if spec['spec']['training'] == sorted(fit) and spec['spec']['all_names'] == names:
                candidates.append((p.parent, spec))
        if len(candidates) != 1: raise ValueError(f'Ambiguous reference cache {axis}')
        path, spec = candidates[0]
        if sha(path/'raw_features.npz') != spec['features_sha256']: raise ValueError('Reference cache corrupted')
        z = np.load(path/'raw_features.npz', allow_pickle=False)
        row_index = {n: i for i, n in enumerate(z['names'])}
        col_index = {n: i for i, n in enumerate(z['clusters'])}
        defs = pd.read_csv(work/f'{axis}_feature_definitions.csv')
        rows = np.array([row_index[n] for n in all_pairs[axis]])
        train_rows = np.array([row_index[n] for n in fit])
        for feature, members in defs.groupby('Feature'):
            columns = [col_index[c] for c in members.Cluster_Label]
            train_pattern = z['presence'][train_rows][:, columns]
            if not (train_pattern == train_pattern[:, :1]).all(): raise ValueError('Merged unequal development patterns')
            expected = z['presence'][rows][:, columns].max(axis=1)
            np.testing.assert_array_equal(expected, all_pairs[feature].to_numpy())
    return query


def main():
    for fold in range(5):
        if not (RUN/f'0.001_fold{fold}/prediction_lock.json').exists():
            raise RuntimeError(f'Incomplete full-budget comparison: missing outer fold {fold}')
    check_sources()
    folds = pd.read_csv(BASE/'matched_outer_folds.csv', dtype={'grouping': str})
    folds = folds[folds.grouping == '0.001']
    with zipfile.ZipFile(ROOT/'data/external_validation/genophi_2026/processed_data.zip') as z:
        raw = pd.read_csv(z.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'))
    truth = raw[raw.dataset == 'Pseudomonas'][['strain', 'phage', 'interaction']]
    rows, model_hashes, training_counts, retained_searches = [], {}, [], []
    for fold in range(5):
        work = RUN/f'0.001_fold{fold}'
        meta = json.loads((work/'inputs.json').read_text())
        expected_test = set(folds[folds.fold == fold].host)
        expected_train = set(folds[folds.fold != fold].host)
        if set(meta['test_hosts']) != expected_test or set(meta['training_hosts']) != expected_train: raise ValueError('Wrong outer fold')
        if set(folds[folds.fold == fold].group) & set(folds[folds.fold != fold].group): raise ValueError('Group leakage')
        dev = pd.read_csv(ROOT/meta['development_input'])
        if set(dev.strain) != expected_train: raise ValueError('Outer test entered development')
        dev_truth = dev.merge(truth, on=['strain', 'phage'], validate='one_to_one', suffixes=('', '_original'))
        if len(dev_truth) != len(expected_train)*19: raise ValueError('Incomplete development labels')
        np.testing.assert_array_equal(dev_truth.interaction, dev_truth.interaction_original)
        if len(list((work/'stages').glob('selection_*.json'))) != 25: raise ValueError('Reduced selection budget')
        for table in (work/'tables').glob('*.csv'):
            cutoff = int(table.stem.split('_')[-1])
            if len(list((work/'stages').glob(f'model_{cutoff}_*.json'))) != 50: raise ValueError('Reduced modeling budget')
            if not set(pd.read_csv(table).strain) <= expected_train: raise ValueError('External rows in feature selection')
            for repeat in range(50):
                directory = work/f'modeling/cutoff_{cutoff}/run_{repeat}'
                verified = verify_model_selection(directory)
                retained_searches.append({'fold': fold, 'cutoff': cutoff, 'repeat': repeat, **verified})
        query = verify_representation(work, meta)
        lock = json.loads((work/'prediction_lock.json').read_text())
        if sha(work/'outer_predictions.csv') != lock['prediction_sha256']: raise ValueError('Prediction changed')
        predictions = pd.read_csv(work/'outer_predictions.csv').set_index(['strain', 'phage']).sort_index()
        query = query.set_index(['strain', 'phage']).sort_index()
        if not predictions.index.equals(query.index): raise ValueError('Prediction alignment changed')
        scores = []
        for repeat in range(50):
            model_path = work/f'modeling/cutoff_{lock["chosen_cutoff"]}/run_{repeat}/best_model.pkl'
            model = joblib.load(model_path)
            scores.append(model.predict_proba(query[model.feature_names_], thread_count=1)[:, 1])
            model_hashes[str(model_path.relative_to(ROOT))] = sha(model_path)
            del model
        np.testing.assert_allclose(np.median(scores, axis=0), predictions.Confidence, rtol=0, atol=1e-12)
        frame = predictions.reset_index().rename(columns={'Confidence': 'genophi'})
        frame['fold'] = fold; frame['group'] = frame.strain.map(folds.set_index('host').group)
        rows.append(frame)
        training_counts.append({'fold': fold, 'selection_runs': 25,
                                'eligible_cutoffs': len(list((work/'tables').glob('*.csv'))),
                                'chosen_cutoff': lock['chosen_cutoff'], 'selected_models_reloaded': 50})
    frame = pd.concat(rows, ignore_index=True).merge(truth, on=['strain', 'phage'], how='outer', validate='one_to_one')
    if len(frame) != 437 or frame[['genophi', 'interaction']].isna().any().any(): raise ValueError('Cohort mismatch')
    baseline_path = ROOT/'data/results_v2/additional_banks/pair_predictions.csv'
    prep = json.loads((BASE/'preparation_checks.json').read_text())
    if sha(baseline_path) != prep['original_pair_predictions_sha256']: raise ValueError('Baseline changed')
    old = pd.read_csv(baseline_path, dtype={'grouping': str})
    old = old[(old.bank == 'Pseudomonas') & (old.grouping == '0.001')].rename(columns={'host': 'strain', 'score': 'neighbour', 'breadth_score': 'prevalence', 'fold': 'old_fold'})
    frame = frame.merge(old[['strain', 'phage', 'neighbour', 'prevalence', 'label', 'old_fold']], on=['strain', 'phage'], validate='one_to_one')
    np.testing.assert_array_equal(frame.interaction, frame.label)
    np.testing.assert_array_equal(frame.fold, frame.old_fold)
    frame.to_csv(RUN/'verified_predictions.csv', index=False)
    methods = ['genophi', 'neighbour', 'prevalence']
    y = frame.label.to_numpy(); groups, inverse = np.unique(frame.group, return_inverse=True)
    point = {m: np.array([roc_auc_score(y, frame[m]), average_precision_score(y, frame[m])]) for m in methods}
    rng = np.random.default_rng(20261005); boot = {m: [] for m in methods}; skipped = 0
    for _ in range(2000):
        counts = np.bincount(rng.integers(len(groups), size=len(groups)), minlength=len(groups))
        weights = counts[inverse]
        if len(np.unique(y[weights > 0])) < 2: skipped += 1; continue
        for m in methods:
            boot[m].append([roc_auc_score(y, frame[m], sample_weight=weights),
                            average_precision_score(y, frame[m], sample_weight=weights)])
    summaries = []
    for m in methods:
        bounds = np.quantile(boot[m], [.025, .975], axis=0)
        for j, metric in enumerate(['auroc', 'average_precision']):
            summaries.append({'comparison': m, 'metric': metric, 'estimate': point[m][j], 'lo': bounds[0, j], 'hi': bounds[1, j]})
    for m in methods[1:]:
        bounds = np.quantile(np.asarray(boot['genophi'])-np.asarray(boot[m]), [.025, .975], axis=0)
        for j, metric in enumerate(['auroc', 'average_precision']):
            summaries.append({'comparison': 'genophi_minus_'+m, 'metric': metric,
                              'estimate': point['genophi'][j]-point[m][j], 'lo': bounds[0, j], 'hi': bounds[1, j]})
    host_rows = []
    for host, z in frame.groupby('strain'):
        row = {'host': host, 'group': z.group.iloc[0]}
        for m in methods: row[m] = screening_curve(z.label, z[m], z.phage)['hit'][2]
        host_rows.append(row)
    hosts = pd.DataFrame(host_rows); hosts.to_csv(RUN/'screening_hosts.csv', index=False)
    for m in methods:
        ci = cluster_mean_ci(hosts[m], hosts.group)
        summaries.append({'comparison': m, 'metric': 'three_phage_union_hit', 'estimate': ci['mean'], 'lo': ci['lo'], 'hi': ci['hi']})
    for m in methods[1:]:
        ci = cluster_mean_ci(hosts.genophi-hosts[m], hosts.group)
        summaries.append({'comparison': 'genophi_minus_'+m, 'metric': 'three_phage_union_hit', 'estimate': ci['mean'], 'lo': ci['lo'], 'hi': ci['hi']})
    summary = pd.DataFrame(summaries); summary.to_csv(RUN/'matched_metrics.csv', index=False)
    audit = {'passed': True, 'audited_utc': datetime.now(timezone.utc).isoformat(), 'folds': 5,
             'selected_models_reloaded': 250, 'pairs': 437, 'hosts': 23, 'phages': 19, 'host_groups': len(groups),
             'requested_bootstrap': 2000, 'one_class_replicates_skipped': skipped,
             'training': training_counts, 'model_sha256': model_hashes,
             'all_model_searches_checked': len(retained_searches),
             'all_negative_model_searches_retained': sum(r['all_candidates_negative'] for r in retained_searches),
             'scope': 'Full-budget adapted GenoPHI modeling with negative-MCC retention repair and Windows reference adapter; retrospective matched fixed-bank comparison.'}
    (RUN/'model_search_audit.json').write_text(json.dumps(retained_searches, indent=2))
    report = '# Full-budget adapted GenoPHI modeling: matched results\n\n'
    report += 'Five outer folds and 250 selected ensemble models passed reload and representation checks. '
    report += 'The 437 pairs, 23 hosts, and 19 phages are matched across methods.\n\n'
    report += '| Comparison | Metric | Estimate | Lower 95% | Upper 95% |\n|---|---|---:|---:|---:|\n'
    for r in summary.itertuples(index=False):
        report += f'| {r.comparison} | {r.metric} | {r.estimate:.4f} | {r.lo:.4f} | {r.hi:.4f} |\n'
    report += '\n\nIntervals condition on saved fits and resample 17 host groups. They do not include all training uncertainty. '
    report += 'Three-phage union activity derives from individual assays, not mixture experiments. '
    report += 'This is not frozen cross-source transfer or an exact end-to-end published GenoPHI reproduction. '
    report += 'A disclosed initialization repair retains the best model even when all internal validation MCC scores are negative. '
    report += 'Read GENOPHI_OFFICIAL_FULL_BUDGET_PROTOCOL.md for all adaptations. '
    report += 'Intervals spanning zero do not establish an advantage. All results are retained regardless of direction.\n'
    (ROOT/'GENOPHI_OFFICIAL_FULL_BUDGET_RESULTS.md').write_text(report, encoding='utf-8')
    audit['auditor_sha256'] = sha(Path(__file__))
    audit['results_sha256'] = {name: sha(RUN/name) for name in ['matched_metrics.csv', 'verified_predictions.csv', 'screening_hosts.csv', 'model_search_audit.json']}
    (RUN/'audit.json').write_text(json.dumps(audit, indent=2))
    print(summary.to_string(index=False))


if __name__ == '__main__': main()
