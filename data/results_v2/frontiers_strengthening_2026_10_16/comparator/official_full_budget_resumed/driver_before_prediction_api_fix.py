"""Full repetition budget using unmodified official GenoPHI modeling routines.

See GENOPHI_OFFICIAL_FULL_BUDGET_PROTOCOL.md for outer holdouts and adaptations.
Subprocess isolation prevents long-lived native estimators accumulating state.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'data/external_validation/genophi_comparison/official_source'
BASE = ROOT/'data/results_v2/genophi_comparison'
RUN = BASE/'official_modeling_full_budget_v1'
CUTOFFS = [3,4,5,6,7,8,9,10,12,15,17,20,22,25,27,30,32,35,37,40,42,45,47,50]


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def validate_source():
    source = json.loads((BASE/'official_source_manifest.json').read_text())
    for row in source['files']:
        if sha(SOURCE/row['path']) != row['sha256']: raise ValueError('Official source modified')
    return source


def prepare(work, fold_number):
    import numpy as np
    import pandas as pd
    import zipfile
    from precisionphage.benchmarks.genophi import ProteinReferences
    from genophi.mmseqs2_clustering import feature_selection_optimized, get_genome_assignments_tables, feature_assignment, merge_feature_tables
    with zipfile.ZipFile(ROOT/'data/external_validation/genophi_2026/processed_data.zip') as z:
        all_data = pd.read_csv(z.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'))
    data = all_data[all_data.dataset == 'Pseudomonas'].copy()
    folds = pd.read_csv(BASE/'matched_outer_folds.csv', dtype={'grouping': str})
    split = folds[folds.grouping == '0.001']
    training = sorted(split[split.fold != fold_number].host)
    test = sorted(split[split.fold == fold_number].host)
    if set(split[split.fold != fold_number].group) & set(split[split.fold == fold_number].group): raise ValueError('Outer group leakage')
    hosts, phages = sorted(data.strain.unique()), sorted(data.phage.unique())
    reference = ProteinReferences(ROOT)
    tables, paths = {}, {}
    for axis, names, fit in [('strain', hosts, training), ('phage', phages, phages)]:
        raw, cluster_names, cache = reference.build(fit, names, 'host' if axis == 'strain' else 'phage')
        presence = pd.DataFrame(raw, columns=cluster_names)
        presence.insert(0, axis, names)
        fit_presence = presence[presence[axis].isin(fit)].copy()
        selected = feature_selection_optimized(fit_presence.copy(), axis, axis)
        # Official hash grouping must not merge unequal patterns by collision.
        for _, members in selected.groupby('Feature'):
            patterns = fit_presence[members.Cluster_Label.tolist()].to_numpy()
            if not (patterns == patterns[:, :1]).all(): raise ValueError('Official feature hash collision')
        assignments = get_genome_assignments_tables(presence.copy(), axis)
        _, wide = feature_assignment(assignments, selected, axis)
        wide = wide.set_index(axis).reindex(names).fillna(0).astype(int).reset_index()
        selected.to_csv(work/f'{axis}_feature_definitions.csv', index=False)
        wide.to_csv(work/f'{axis}_all_input_features.csv', index=False)
        train_table = work/f'{axis}_development_features.csv'
        wide[wide[axis].isin(fit)].to_csv(train_table, index=False)
        tables[axis], paths[axis] = wide, train_table
    phenotype = work/'development_labels.csv'
    data[data.strain.isin(training)][['strain', 'phage', 'interaction']].to_csv(phenotype, index=False)
    (work/'merged').mkdir(parents=True, exist_ok=True)
    merged = merge_feature_tables(str(paths['strain']), str(phenotype), str(work/'merged'),
                                 sample_column='strain', phage_features=str(paths['phage']),
                                 phenotype_column='interaction', remove_suffix=False,
                                 use_feature_clustering=False)
    # Input-only cartesian query table, no outer labels passed to official code.
    query = pd.MultiIndex.from_product([test, phages], names=['strain', 'phage']).to_frame(index=False)
    query = query.merge(tables['strain'], on='strain', validate='many_to_one').merge(tables['phage'], on='phage', validate='many_to_one')
    query.to_csv(work/'outer_query_features.csv', index=False)
    if len(pd.read_csv(merged)) != len(training)*len(phages): raise ValueError('Development rows lost')
    meta = {'fold': fold_number, 'training_hosts': training, 'test_hosts': test, 'phages': phages,
            'development_input': str(Path(merged).relative_to(ROOT)), 'num_features': 50,
            'development_input_sha256': sha(merged), 'outer_query_sha256': sha(work/'outer_query_features.csv')}
    (work/'inputs.json').write_text(json.dumps(meta, indent=2))


def worker(args):
    validate_source()
    sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(SOURCE))
    import pandas as pd
    work = RUN/f'0.001_fold{args.fold}'; work.mkdir(parents=True, exist_ok=True)
    if args.worker == 'prepare':
        prepare(work, args.fold); return
    meta = json.loads((work/'inputs.json').read_text())
    inp = ROOT/meta['development_input']
    if sha(inp) != meta['development_input_sha256']: raise ValueError('Development input changed')
    common = dict(threads=4, sample_column='strain', phenotype_column='interaction',
                  phage_column='phage', use_dynamic_weights=True, weights_method='inverse_frequency',
                  max_ram=4, use_clustering=True, cluster_method='hierarchical', n_clusters=20)
    if args.worker == 'selection':
        from genophi.feature_selection import run_feature_selection_iterations
        # Prior successful repetitions are skipped by the official function.
        run_feature_selection_iterations(str(inp), str(work/'selection'), num_features=50,
            filter_type='strain', num_runs=args.repeat+1, method='rfe', **common)
        path = work/f'selection/run_{args.repeat}/feature_importances.csv'
        if not path.exists() or pd.read_csv(path).empty: raise ValueError('Official selection did not complete')
    elif args.worker == 'tables':
        from genophi.feature_selection import generate_feature_tables
        generate_feature_tables(str(work/'selection'), str(inp), str(work/'tables'),
                                phenotype_column='interaction', sample_column='strain', cut_offs=CUTOFFS,
                                binary_data=True, max_features=50, min_features=5, filter_type='strain')
        if not list((work/'tables').glob('*.csv')): raise ValueError('No eligible official feature cutoff')
    elif args.worker == 'model':
        from genophi.select_feature_modeling import model_testing_select_MCC
        dest = work/f'modeling/cutoff_{args.cutoff}/run_{args.repeat}'
        dest.mkdir(parents=True, exist_ok=True)
        table = work/f'tables/select_feature_table_cutoff_{args.cutoff}.csv'
        model_testing_select_MCC(str(table), str(dest), random_state=args.repeat,
                                 set_filter='strain', binary_data=True, use_shap=False, **common)
        if not (dest/'best_model.pkl').exists(): raise ValueError('Official modeling returned no model')
    elif args.worker == 'predict':
        from genophi.select_feature_modeling import parse_model_predictions_and_performance, evaluate_model_performance
        from genophi.workflows.prediction_workflow import predict_with_models, calculate_median_predictions
        predictions, performance = parse_model_predictions_and_performance(str(work/'modeling'))
        path = work/'development_predictions.csv'; predictions.to_csv(path, index=False)
        performance.to_csv(work/'development_performance.csv', index=False)
        evaluate_model_performance(str(path), str(work/'development_summary'), sample_column='strain', phenotype_column='interaction')
        metrics = pd.read_csv(work/'development_summary/model_performance/model_performance_metrics.csv')
        cutoff = str(metrics.iloc[0].cut_off).split('_')[-1]
        model_dir = work/f'modeling/cutoff_{cutoff}'
        if len(list(model_dir.glob('run_*/best_model.pkl'))) != 50: raise ValueError('Incomplete ensemble')
        query_path = work/'outer_query_features.csv'
        if sha(query_path) != meta['outer_query_sha256']: raise ValueError('Outer input changed')
        query = pd.read_csv(query_path)
        pred = predict_with_models(str(model_dir), query, threads=4)
        if len(pred) != len(query)*50 or pred.run.nunique() != 50: raise ValueError('Missing official predictions')
        median = calculate_median_predictions(pred)
        pred.to_csv(work/'outer_model_predictions.csv', index=False)
        median.to_csv(work/'outer_predictions.csv', index=False)
        (work/'prediction_lock.json').write_text(json.dumps({'recorded_utc': datetime.now(timezone.utc).isoformat(),
            'fold': args.fold, 'chosen_cutoff': cutoff, 'selection_repetitions': 25, 'ensemble_models': 50,
            'prediction_sha256': sha(work/'outer_predictions.csv'), 'outer_outcomes_evaluated': False}, indent=2))


def run_stage(work, fold, stage, repeat=None, cutoff=None):
    key = '_'.join(str(v) for v in [stage, cutoff, repeat] if v is not None)
    marker = work/'stages'/f'{key}.json'
    if marker.exists():
        saved = json.loads(marker.read_text())
        for path, checksum in saved['sha256'].items():
            if sha(ROOT/path) != checksum: raise ValueError('Completed stage artifact changed')
        return
    args = [sys.executable, '-X', 'faulthandler', str(Path(__file__)), '--worker', stage, '--fold', str(fold)]
    if repeat is not None: args += ['--repeat', str(repeat)]
    if cutoff is not None: args += ['--cutoff', str(cutoff)]
    started = time.time()
    (work/'stages').mkdir(parents=True, exist_ok=True)
    state = {'state': 'running', 'fold': fold, 'stage': stage, 'repeat': repeat, 'cutoff': cutoff,
             'updated_utc': datetime.now(timezone.utc).isoformat(), 'pid': os.getpid()}
    (RUN/'status.json').write_text(json.dumps(state, indent=2)); print(json.dumps(state), flush=True)
    with (work/'stages'/f'{key}.log').open('w', encoding='utf-8') as log:
        subprocess.run(args, cwd=ROOT, stdout=log, stderr=log, check=True)
    targets = {'prepare': [work/'inputs.json', work/'outer_query_features.csv'],
               'selection': [work/f'selection/run_{repeat}/feature_importances.csv'],
               'tables': list((work/'tables').glob('*.csv')),
               'model': [work/f'modeling/cutoff_{cutoff}/run_{repeat}/best_model.pkl'],
               'predict': [work/'prediction_lock.json', work/'outer_predictions.csv']}[stage]
    marker.write_text(json.dumps({'exit_code': 0, 'seconds': time.time()-started,
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'sha256': {str(p.relative_to(ROOT)): sha(p) for p in targets}}, indent=2))


def main(args):
    validate_source(); RUN.mkdir(parents=True, exist_ok=True)
    files = [Path(__file__), ROOT/'GENOPHI_OFFICIAL_FULL_BUDGET_PROTOCOL.md',
             ROOT/'src/precisionphage/benchmarks/genophi.py', BASE/'matched_outer_folds.csv', BASE/'official_source_manifest.json']
    provenance = {'sha256': {str(p.relative_to(ROOT)): sha(p) for p in files},
                  'packages': {p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scikit-learn', 'catboost', 'plotnine', 'mizani', 'statsmodels', 'shap']}}
    manifest = RUN/'run_manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()) != provenance: raise ValueError('Run source/dependencies changed')
    if not manifest.exists(): manifest.write_text(json.dumps(provenance, indent=2))
    try:
        for fold in ([args.fold] if args.fold is not None else range(5)):
            work = RUN/f'0.001_fold{fold}'; work.mkdir(parents=True, exist_ok=True)
            run_stage(work, fold, 'prepare')
            for repeat in range(25): run_stage(work, fold, 'selection', repeat)
            run_stage(work, fold, 'tables')
            cutoffs = sorted(int(p.stem.split('_')[-1]) for p in (work/'tables').glob('*.csv'))
            for cutoff in cutoffs:
                for repeat in range(50): run_stage(work, fold, 'model', repeat, cutoff)
            run_stage(work, fold, 'predict')
        state = {'state': 'requested_folds_predicted_pending_independent_audit',
                 'folds': [args.fold] if args.fold is not None else list(range(5)),
                 'updated_utc': datetime.now(timezone.utc).isoformat()}
        (RUN/'status.json').write_text(json.dumps(state, indent=2)); print(json.dumps(state), flush=True)
    except Exception as exc:
        prior = json.loads((RUN/'status.json').read_text()) if (RUN/'status.json').exists() else {}
        prior.update(state='failed', error=repr(exc), updated_utc=datetime.now(timezone.utc).isoformat())
        (RUN/'status.json').write_text(json.dumps(prior, indent=2)); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--fold', type=int, choices=range(5))
    parser.add_argument('--worker', choices=['prepare', 'selection', 'tables', 'model', 'predict'])
    parser.add_argument('--repeat', type=int); parser.add_argument('--cutoff', type=int)
    args = parser.parse_args()
    worker(args) if args.worker else main(args)
