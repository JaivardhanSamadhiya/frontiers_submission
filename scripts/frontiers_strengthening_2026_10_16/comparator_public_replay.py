"""Read-only public replay of fitted GenoPHI models and cached representations.

This is deliberately narrower than the immutable historical local audit. It
does not replay private execution logs, the complete repair-registration
history, raw proteome downloads, or protein clustering. It performs no fitting
and writes no files. A companion archive may supply its enumerated dependencies.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError('Run without Python optimization; verification assertions are mandatory')

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import sys
import zipfile

sys.dont_write_bytecode = True

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUN_REL = Path('data/results_v2/frontiers_strengthening_2026_10_16/comparator/official_full_budget_resumed')
RUN = ROOT / RUN_REL
BASE = ROOT / 'data/results_v2/genophi_comparison'
HISTORICAL = Path('data/results_v2/genophi_comparison/official_modeling_full_budget_v1')
STEM = Path('scripts/frontiers_strengthening_2026_10_16')
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / STEM))
from precisionphage.benchmarks.genophi_retention import verify_model_selection
import verify_comparator_metrics as metrics

LIMITATIONS = [
    'No fitting, new predictions for an unobserved cohort, or outcome-blinded claim.',
    'Cached protein-family arrays and merged features are checked; raw-proteome annotation, alignment and clustering are not reconstructed.',
    'Private historical execution logs and the complete repair-registration/preserved-file history are intentionally not replayed or certified.',
    'This portable public replay is not a replacement for the immutable full local historical audit.',
    'Adapted official GenoPHI with Windows reference and negative-MCC retention adaptations; not an exact published end-to-end reproduction.',
]


def relpath(value):
    """Interpret recorded Windows separators on every platform, without escape."""
    text = str(value).replace('\\', '/')
    path = Path(text)
    if text.startswith('/') or path.is_absolute() or '..' in path.parts or ':' in text:
        raise ValueError(f'Non-relative scientific artifact path: {text}')
    return path


def recorded_file(value):
    """Alias old checkpoint forest paths to the byte-identical resumed forest."""
    path = relpath(value)
    if path.is_relative_to(HISTORICAL):
        tail = path.relative_to(HISTORICAL)
        if tail.parts and tail.parts[0].startswith('0.001_fold'):
            path = RUN_REL / tail
    return ROOT / path


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path, **kwargs):
    return pd.read_csv(path, float_precision='round_trip', **kwargs)


def selected_folds(folds):
    values = list(range(5)) if folds is None else list(folds)
    if not values or len(set(values)) != len(values) or not set(values) <= set(range(5)):
        raise ValueError('Provide distinct folds from 0 through 4')
    return sorted(values)


def bind_root(root):
    """Bind an extracted core-plus-companion tree without rewriting its files."""
    global ROOT, RUN, BASE
    ROOT = Path(root).resolve()
    RUN = ROOT / RUN_REL
    BASE = ROOT / 'data/results_v2/genophi_comparison'


def cached_references(meta):
    found = []
    for kind, names, training in [
        ('host', sorted(set(meta['training_hosts'] + meta['test_hosts'])), sorted(meta['training_hosts'])),
        ('phage', sorted(meta['phages']), sorted(meta['phages'])),
    ]:
        matches = []
        for path in (ROOT / 'data/external_validation/genophi_comparison/feature_cache_v2').glob(f'{kind}_*/complete.json'):
            record = read_json(path)
            if record['spec']['training'] == training and record['spec']['all_names'] == names:
                matches.append(path)
        if len(matches) != 1:
            raise ValueError(f'Expected one exact {kind} reference cache, found {len(matches)}')
        found.extend([matches[0], matches[0].parent / 'raw_features.npz'])
    return found


def source_files():
    """Scientific source identity checks, excluding preserved execution history."""
    paths = {BASE / 'official_source_manifest.json', RUN / 'run_manifest.json',
             RUN / 'execution_amendment.json', RUN / 'negative_mcc_repair_manifest.json',
             RUN / 'driver_before_prediction_api_fix.py',
             RUN / 'before_negative_mcc_repair/47_official_genophi_outer_comparison.py'}
    official = read_json(BASE / 'official_source_manifest.json')
    for record in official['files']:
        paths.add(ROOT / 'data/external_validation/genophi_comparison/official_source' / relpath(record['path']))
    for name in read_json(RUN / 'run_manifest.json')['sha256']:
        paths.add(ROOT / relpath(name))
    for name in read_json(RUN / 'negative_mcc_repair_manifest.json')['current_files_sha256']:
        paths.add(ROOT / relpath(name))
    return paths


def check_source_identity():
    official = read_json(BASE / 'official_source_manifest.json')
    for record in official['files']:
        path = ROOT / 'data/external_validation/genophi_comparison/official_source' / relpath(record['path'])
        assert sha(path) == record['sha256'], 'Pinned official source changed'
    repair = read_json(RUN / 'negative_mcc_repair_manifest.json')
    for name, checksum in repair['current_files_sha256'].items():
        assert sha(ROOT / relpath(name)) == checksum, 'Scientific repair source changed'
    manifest = read_json(RUN / 'run_manifest.json')
    for name, checksum in manifest['sha256'].items():
        path = ROOT / relpath(name)
        if sha(path) == checksum:
            continue
        assert relpath(name).as_posix() == 'experiments/47_official_genophi_outer_comparison.py'
        assert sha(RUN / 'driver_before_prediction_api_fix.py') == checksum
        assert sha(path) == repair['current_files_sha256'][name]
    amendment = read_json(RUN / 'execution_amendment.json')
    assert sha(RUN / 'driver_before_prediction_api_fix.py') == amendment['old_driver_sha256']
    assert sha(RUN / 'before_negative_mcc_repair/47_official_genophi_outer_comparison.py') == amendment['new_driver_sha256']
    return {'official_commit': official['commit'], 'official_source_files': len(official['files']),
            'current_scientific_source_hashes_checked': len(repair['current_files_sha256']),
            'original_driver_snapshots_checked': 2,
            'complete_private_repair_history_checked': False}


def enumerate_companion_files(folds=None, core_root=None):
    """Enumerate exact existing replay dependencies, optionally minus a core tree.

    Recorded merged development inputs retain their original repository paths.
    Old *forest* hash references instead resolve to the resumed copy: private
    historical forests are not duplicated. Hash equality is checked by replay.
    No archive or manifest is written by this function.
    """
    folds = selected_folds(folds)
    required = set(source_files())
    required.update(ROOT / name for name in [
        STEM / 'comparator_public_replay.py', STEM / 'verify_comparator_metrics.py',
        Path('experiments/49_audit_official_genophi_outer.py'),
        Path('src/precisionphage/__init__.py'), Path('src/precisionphage/benchmarks/__init__.py'),
        Path('src/precisionphage/benchmarks/screening.py'), Path('src/precisionphage/benchmarks/genophi_retention.py'),
        Path('requirements-genophi-comparator.txt'), Path('data/results_v2/genophi_comparison/matched_outer_folds.csv'),
        Path('data/results_v2/genophi_comparison/preparation_checks.json'),
        Path('data/results_v2/frontiers_strengthening_2026_10_16/comparator/source_data_license.json'),
        Path('data/external_validation/genophi_2026/processed_data.zip'),
        Path('data/results_v2/additional_banks/pair_predictions.csv')])
    for fold in folds:
        work = RUN / f'0.001_fold{fold}'
        meta = read_json(work / 'inputs.json')
        required.add(ROOT / relpath(meta['development_input']))
        required.update(cached_references(meta))
        for name in ['inputs.json', 'outer_query_features.csv', 'phage_feature_definitions.csv',
                     'strain_feature_definitions.csv', 'prediction_lock.json', 'outer_predictions.csv',
                     'development_predictions.csv', 'development_summary/model_performance/model_performance_metrics.csv']:
            required.add(work / name)
        for marker in (work / 'stages').glob('*.json'):
            required.add(marker)
            required.update(recorded_file(name) for name in read_json(marker)['sha256'])
        for table in (work / 'tables').glob('*.csv'):
            required.add(table)
            cutoff = int(table.stem.split('_')[-1])
            for repeat in range(50):
                for name in ['best_model.pkl', 'best_model_predictions.csv', 'model_performance.csv']:
                    required.add(work / f'modeling/cutoff_{cutoff}/run_{repeat}' / name)
    for name in ['audit.json', 'model_search_audit.json', 'verified_predictions.csv', 'matched_metrics.csv', 'screening_hosts.csv']:
        if (RUN / name).is_file():
            required.add(RUN / name)
    missing = [p.relative_to(ROOT).as_posix() for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError(f'Incomplete replay dependency closure: {missing[:12]}')
    if any(p.suffix.lower() == '.log' for p in required):
        raise ValueError('Execution logs must not enter the public replay closure')
    core_root = Path(core_root).resolve() if core_root is not None else None
    rows, core_count = [], 0
    for path in sorted(required):
        relative = path.relative_to(ROOT)
        checksum = sha(path)
        in_core = core_root is not None and (core_root / relative).is_file()
        if in_core:
            assert sha(core_root / relative) == checksum, f'Core content differs: {relative}'
            core_count += 1
            continue
        rows.append({'path': relative.as_posix(), 'bytes': path.stat().st_size, 'sha256': checksum})
    return {'scope': 'Exact existing dependencies of the narrower public fitted-model/cached-representation replay',
            'folds': folds, 'complete_five_fold_closure': folds == list(range(5)),
            'files_already_in_identical_core': core_count, 'extra_files': len(rows),
            'extra_bytes': sum(row['bytes'] for row in rows), 'files': rows,
            'limitations': LIMITATIONS, 'files_written': 0}


def dependency_paths(root, folds=range(5)):
    """Exact science-only closure; caller may subtract an identical core tree.

    Includes small scientific source snapshots and source-identity manifests,
    not logs, failed-run status history or duplicate historical model forests.
    """
    bind_root(root)
    inventory = enumerate_companion_files(folds)
    return [ROOT / relpath(row['path']) for row in inventory['files']]


def verify(root=ROOT, folds=range(5)):
    bind_root(root)
    folds = selected_folds(folds)
    full = folds == list(range(5))
    for fold in folds:
        if not (RUN / f'0.001_fold{fold}/prediction_lock.json').exists():
            raise RuntimeError(f'Fold {fold} has not completed locked predictions')
    if full:
        original_audit = read_json(RUN / 'audit.json')
        assert original_audit['passed'] and original_audit['folds'] == 5
        assert original_audit['selected_models_reloaded'] == 250
        assert sha(ROOT / 'experiments/49_audit_official_genophi_outer.py') == original_audit['auditor_sha256']
        for name, checksum in original_audit['results_sha256'].items():
            assert sha(RUN / relpath(name)) == checksum
    else:
        original_audit = None
    source_check = check_source_identity()
    spec = importlib.util.spec_from_file_location('public_replay_original_representation', ROOT / 'experiments/49_audit_official_genophi_outer.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    original.ROOT, original.RUN = ROOT, RUN
    metrics.ROOT, metrics.RUN = ROOT, RUN
    fold_table = read_csv(BASE / 'matched_outer_folds.csv', dtype={'grouping': str})
    fold_table = fold_table[fold_table.grouping.eq('0.001')]
    with zipfile.ZipFile(ROOT / 'data/external_validation/genophi_2026/processed_data.zip') as archive:
        truth = pd.read_csv(archive.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'))
    truth = truth.loc[truth.dataset.eq('Pseudomonas'), ['strain', 'phage', 'interaction']]
    assert len(truth) == 437 and truth.interaction.sum() == 158
    preparation = read_json(BASE / 'preparation_checks.json')
    assert sha(BASE / 'matched_outer_folds.csv') == preparation['fold_assignment_sha256']
    assert sha(ROOT / 'data/results_v2/additional_banks/pair_predictions.csv') == preparation['original_pair_predictions_sha256']
    expected_searches = {}
    if original_audit is not None:
        for row in read_json(RUN / 'model_search_audit.json'):
            key = (row['fold'], row['cutoff'], row['repeat'])
            assert key not in expected_searches
            expected_searches[key] = row
    checked, searches, selected_models, marker_count = [], 0, 0, 0
    for fold in folds:
        work = RUN / f'0.001_fold{fold}'
        meta = read_json(work / 'inputs.json')
        development = ROOT / relpath(meta['development_input'])
        assert sha(development) == meta['development_input_sha256']
        assert sha(work / 'outer_query_features.csv') == meta['outer_query_sha256']
        expected_test = set(fold_table.loc[fold_table.fold.eq(fold), 'host'])
        expected_train = set(fold_table.loc[~fold_table.fold.eq(fold), 'host'])
        assert set(meta['test_hosts']) == expected_test and set(meta['training_hosts']) == expected_train
        assert not set(fold_table.loc[fold_table.fold.eq(fold), 'group']) & set(fold_table.loc[~fold_table.fold.eq(fold), 'group'])
        dev = read_csv(development)
        assert set(dev.strain) == expected_train
        joined = dev.merge(truth, on=['strain', 'phage'], validate='one_to_one', suffixes=('', '_original'))
        assert len(joined) == len(expected_train) * 19
        np.testing.assert_array_equal(joined.interaction, joined.interaction_original)
        markers = list((work / 'stages').glob('*.json'))
        assert {p.stem for p in markers if p.stem.startswith('selection_')} == {f'selection_{i}' for i in range(25)}
        for marker in markers:
            record = read_json(marker)
            assert record['exit_code'] == 0
            for name, checksum in record['sha256'].items():
                assert sha(recorded_file(name)) == checksum, f'Checkpoint target changed in fold {fold}'
            marker_count += 1
        tables = list((work / 'tables').glob('*.csv'))
        for table in tables:
            cutoff = int(table.stem.split('_')[-1])
            assert {p.stem for p in markers if p.stem.startswith(f'model_{cutoff}_')} == {f'model_{cutoff}_{i}' for i in range(50)}
            assert set(read_csv(table).strain) <= expected_train
            for repeat in range(50):
                directory = work / f'modeling/cutoff_{cutoff}/run_{repeat}'
                result = verify_model_selection(directory)
                if expected_searches:
                    saved = expected_searches[(fold, cutoff, repeat)]
                    for key, value in result.items():
                        assert saved[key] == value, f'Original retained-search audit differs: fold {fold}'
                searches += 1
        portable_meta = dict(meta, development_input=relpath(meta['development_input']).as_posix())
        query = original.verify_representation(work, portable_meta).set_index(['strain', 'phage']).sort_index()
        lock = read_json(work / 'prediction_lock.json')
        assert sha(work / 'outer_predictions.csv') == lock['prediction_sha256']
        predictions = read_csv(work / 'outer_predictions.csv').set_index(['strain', 'phage']).sort_index()
        assert query.index.equals(predictions.index)
        ensemble = []
        for repeat in range(50):
            path = work / f'modeling/cutoff_{lock["chosen_cutoff"]}/run_{repeat}/best_model.pkl'
            model = joblib.load(path)
            ensemble.append(model.predict_proba(query[model.feature_names_], thread_count=1)[:, 1])
            if original_audit is not None:
                old_keys = {recorded_file(name): value for name, value in original_audit['model_sha256'].items()}
                assert sha(path) == old_keys[path]
            selected_models += 1
            del model
        np.testing.assert_allclose(np.median(ensemble, axis=0), predictions.Confidence, rtol=0, atol=1e-12)
        selected = metrics.verify_development_selection(work)
        assert selected['eligible_cutoffs'] == len(tables)
        checked.append(selected)
    if full:
        assert searches == original_audit['all_model_searches_checked'] == len(expected_searches)
        assert selected_models == len(original_audit['model_sha256']) == 250
        compact_metrics = metrics.verify()
    else:
        compact_metrics = {'performed': False, 'reason': 'Partial fold audit; aggregate five-fold metrics are not certified'}
    return {'passed': True, 'folds_complete': len(folds), 'complete_five_fold_public_replay': full,
            'scope': 'Fitted-model, retained-search, ensemble-prediction and cached-representation consistency replay; historical execution provenance excluded',
            'audited_utc': datetime.now(timezone.utc).isoformat(), 'folds': folds,
            'source_identity': source_check, 'fold_checks': checked,
            'retained_models_reloaded_for_search_checks': searches, 'selected_ensemble_models_reloaded': selected_models,
            'stage_checkpoints_replayed_with_historical_forest_path_aliases': marker_count,
            'original_final_search_and_selected_model_hashes_checked': full,
            'compact_metric_replay': compact_metrics,
            'packages': {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scikit-learn', 'catboost', 'joblib']},
            'limitations': LIMITATIONS, 'model_fits_started': 0, 'files_written': 0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folds', type=int, choices=range(5), nargs='+', help='Explicit partial replay; omission requires all five completed folds and final audit')
    parser.add_argument('--root', type=Path, default=ROOT, help='Extracted core-plus-companion repository root')
    parser.add_argument('--list-files', action='store_true', help='Print dependency closure instead of loading models')
    parser.add_argument('--core-root', type=Path, help='Exclude byte-identical files already in this core tree when listing dependencies')
    args = parser.parse_args()
    if args.core_root is not None and not args.list_files:
        parser.error('--core-root is only used with --list-files')
    bind_root(args.root)
    result = enumerate_companion_files(args.folds, args.core_root) if args.list_files else verify(args.root, args.folds)
    print(json.dumps(result, indent=2))
