"""Read-only inference portability diagnostics for the 174 saved outer models.

This is NOT an acceptance verifier and changes no numerical tolerance. It
retains the exact saved scores and fitted artifacts, never fits a model, and
reports numerical, ranking, metric and threshold differences without deciding
that any difference is scientifically acceptable. Run the independent strict
artifact/source/split checks as well. Output is JSON; no file is written.
"""
from __future__ import annotations

import argparse
from contextlib import nullcontext
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import platform
import stat
import sys
import zipfile

sys.dont_write_bytecode = True

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from threadpoolctl import threadpool_info, threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
BASE = 'data/results_v2/frontiers_strengthening_2026_10_16/models'
LIMIT = 1024 ** 3
HASH_KEYS = ('numpy', 'pandas', 'scipy', 'scikit-learn', 'xgboost', 'joblib', 'threadpoolctl')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def relative_name(value):
    value = str(value)
    path = PurePosixPath(value)
    require(bool(path.parts) and not path.is_absolute() and '..' not in path.parts
            and '\\' not in value and ':' not in value and '\x00' not in value,
            'Unsafe artifact-relative path')
    return path.as_posix()


class Inputs:
    """Bounded, byte-hash-checked read access to a known scientific snapshot."""
    def __init__(self, root, manifest=None):
        self.root = Path(root).resolve()
        self.read_bytes = 0
        self.checked = {}
        path = Path(manifest).resolve() if manifest else self.root / 'SCIENTIFIC_EXTENSION_MANIFEST.json'
        require(path.is_file() and path.stat().st_size <= 32 * 1024**2,
                'A bounded scientific extension manifest is required (use --artifact-manifest locally)')
        document = json.loads(path.read_text(encoding='utf-8-sig'))
        self.expected = {}
        for row in document['files']:
            name = relative_name(row['path'])
            require(name not in self.expected, 'Duplicate scientific manifest path')
            self.expected[name] = row

    def file(self, value):
        name = relative_name(value)
        path = self.root / name
        require(path.resolve().is_relative_to(self.root), 'Artifact resolves outside the scientific root')
        for item in [path, *path.parents]:
            require(not item.is_symlink(), 'Symlink artifact refused')
            if item.exists():
                attributes = getattr(item.lstat(), 'st_file_attributes', 0)
                require(not attributes & getattr(stat, 'FILE_ATTRIBUTE_REPARSE_POINT', 0), 'Reparse-point artifact refused')
            if item == self.root:
                break
        require(path.is_file() and path.stat().st_size <= LIMIT, 'Missing or oversized artifact')
        if name not in self.checked:
            require(name in self.expected, 'An artifact is not covered by the scientific extension manifest')
            record = self.expected[name]
            size = path.stat().st_size
            require(self.read_bytes + size <= LIMIT, 'The one-GiB distinct-artifact read budget was exceeded')
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(1024**2), b''):
                    digest.update(block)
            require(size == record['bytes'] and digest.hexdigest() == record['sha256'], 'Exact scientific artifact hash mismatch')
            self.checked[name] = {'bytes': size, 'sha256': digest.hexdigest()}
            self.read_bytes += size
        return path

    def json(self, name):
        return json.loads(self.file(name).read_text(encoding='utf-8-sig'))

    def csv(self, name):
        return pd.read_csv(self.file(name), float_precision='round_trip')

    def features(self, name, expected_shape):
        path = self.file(name)
        with zipfile.ZipFile(path) as archive:
            require(sum(row.file_size for row in archive.infolist()) <= LIMIT, 'Oversized decompressed feature archive')
        with np.load(path, allow_pickle=False) as archive:
            value = archive['X']
        require(value.shape == expected_shape and value.nbytes <= LIMIT and np.isfinite(value).all(),
                'Unexpected/non-finite feature matrix')
        return value


def metrics(y, p):
    y, p = np.asarray(y), np.asarray(p, dtype=np.float64)
    n = len(y)
    require(n > 0 and set(np.unique(y)) <= {0, 1} and np.isfinite(p).all()
            and ((p >= 0) & (p <= 1)).all(), 'Invalid diagnostic outcomes/probabilities')
    bins = np.minimum((p * 10).astype(np.int64), 9)
    ece = sum(np.mean(bins == b) * abs(p[bins == b].mean() - y[bins == b].mean())
              for b in range(10) if np.any(bins == b))
    return {'auroc': float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None,
            'average_precision': float(average_precision_score(y, p)),
            'brier': float(brier_score_loss(y, p)), 'ece_10_equal_width': float(ece)}


def pairwise_relations(y, original, replay):
    """Exact unordered-pair relations in bounded blocks; no tolerance or test."""
    y, original, replay = np.asarray(y), np.asarray(original), np.asarray(replay)
    n = len(y)
    counts = {key: 0 for key in ('order_reversals', 'original_ties_broken', 'new_ties_created',
                                 'unchanged_ties', 'positive_negative_order_reversals',
                                 'positive_negative_original_ties_broken', 'positive_negative_new_ties_created')}
    for start in range(0, n, 128):
        stop = min(n, start + 128)
        first = np.sign(original[start:stop, None] - original[None, :])
        second = np.sign(replay[start:stop, None] - replay[None, :])
        eligible = np.arange(start, stop)[:, None] < np.arange(n)[None, :]
        discordant = y[start:stop, None] != y[None, :]
        changes = {'order_reversals': first * second < 0,
                   'original_ties_broken': (first == 0) & (second != 0),
                   'new_ties_created': (first != 0) & (second == 0)}
        for key, mask in changes.items():
            counts[key] += int(np.count_nonzero(mask & eligible))
            counts['positive_negative_' + key] += int(np.count_nonzero(mask & eligible & discordant))
        counts['unchanged_ties'] += int(np.count_nonzero((first == 0) & (second == 0) & eligible))
    return {'unordered_pairs': n * (n - 1) // 2,
            'changed_pairwise_relations': counts['order_reversals'] + counts['original_ties_broken'] + counts['new_ties_created'],
            **counts}


def differences(y, original, replay, rows):
    original, replay = np.asarray(original, dtype=np.float64), np.asarray(replay, dtype=np.float64)
    require(original.shape == replay.shape == np.asarray(y).shape == np.asarray(rows).shape,
            'Diagnostic row alignment mismatch')
    old_metrics, new_metrics = metrics(y, original), metrics(y, replay)
    error = np.abs(original - replay)
    nonzero_original, nonzero_replay = original != 0, replay != 0
    original32, replay32 = original.astype(np.float32), replay.astype(np.float32)
    ulps = np.abs(original32.view(np.uint32).astype(np.int64) - replay32.view(np.uint32).astype(np.int64))
    top = []
    for k in sorted({min(len(original), k) for k in (1, 3, 5, 10, 20)}):
        # Explicit deterministic index tie-break, with boundary ties disclosed.
        old_order, new_order = np.lexsort((rows, -original)), np.lexsort((rows, -replay))
        old_set, new_set = set(np.asarray(rows)[old_order[:k]]), set(np.asarray(rows)[new_order[:k]])
        top.append({'k': k, 'symmetric_difference_count': len(old_set ^ new_set),
                    'original_boundary_tie': bool(k < len(original) and original[old_order[k-1]] == original[old_order[k]])})
    return {'n': len(original), 'positive': int(np.sum(y)),
            'max_absolute_error': float(error.max()), 'different_numeric_values': int(np.count_nonzero(error)),
            'max_relative_error_to_original_nonzero': float(np.max(error[nonzero_original] / abs(original[nonzero_original]))) if nonzero_original.any() else None,
            'max_relative_error_to_replay_nonzero': float(np.max(error[nonzero_replay] / abs(replay[nonzero_replay]))) if nonzero_replay.any() else None,
            'zero_original_nonzero_replay': int(np.count_nonzero(~nonzero_original & nonzero_replay)),
            'max_float32_ulp_difference': int(ulps.max()),
            'historical_rtol_1e_6_atol_1e_7_failing_rows': int(np.count_nonzero(~np.isclose(original, replay, rtol=1e-6, atol=1e-7))),
            'independent_rtol_1e_7_atol_1e_7_failing_rows': int(np.count_nonzero(~np.isclose(original, replay, rtol=1e-7, atol=1e-7))),
            'threshold_ge_0_5_decision_changes': int(np.count_nonzero((original >= .5) != (replay >= .5))),
            'threshold_gt_0_5_decision_changes': int(np.count_nonzero((original > .5) != (replay > .5))),
            'original_metrics': old_metrics, 'replayed_metrics': new_metrics,
            'metric_differences': {key: None if old_metrics[key] is None else new_metrics[key] - old_metrics[key] for key in old_metrics},
            'pairwise': pairwise_relations(y, original, replay), 'top_k_index_tiebroken': top}


def runtime():
    return {'python': platform.python_version(), 'system': platform.system(), 'machine': platform.machine(),
            'packages': {name: importlib.metadata.version(name) for name in HASH_KEYS},
            'threadpools': [{key: row.get(key) for key in ('user_api', 'internal_api', 'prefix', 'version',
                                                          'num_threads', 'threading_layer', 'architecture')}
                            for row in threadpool_info()]}


def diagnose(root=ROOT, manifest=None, threads=None):
    inputs = Inputs(root, manifest)
    reports, pooled = [], {}
    with threadpool_limits(limits=threads) if threads is not None else nullcontext():
        environment = runtime()
        sys.stderr.write(json.dumps({'diagnostic_only': True, 'runtime': environment}, sort_keys=True) + '\n')
        for analysis, expected_shape, families in [('secondsource', (1053, 448), ('phage_only', 'pair_only', 'full')),
                                                   ('nested_comparison', (1947, 448), ('GBM', 'MLP')),
                                                   ('mlp_budget_sensitivity', (1947, 448), ('MLP_extended',))]:
            source = 'nested_comparison' if analysis == 'mlp_budget_sensitivity' else analysis
            source_base, destination = BASE + '/' + source, BASE + '/' + analysis
            frame, X = inputs.csv(source_base + '/input_rows.csv'), inputs.features(source_base + '/input_features.npz', expected_shape)
            require(len(frame) == len(X), 'Input rows do not match feature matrix')
            splits = inputs.json(source_base + '/splits.json')
            require(len(splits) == (48 if analysis == 'secondsource' else 10), 'Unexpected split count')
            saved = inputs.csv(destination + '/predictions.csv')
            completion = inputs.json(destination + '/completion.json')
            require(completion.get('complete') is True, 'The scientific analysis is incomplete')
            if analysis != 'mlp_budget_sensitivity':
                protocol = inputs.file(destination + '/retrospective_protocol.json')
                require(hashlib.sha256(protocol.read_bytes()).hexdigest() == completion['protocol_sha256'], 'Protocol/completion hash mismatch')
            for split in splits:
                regime, fold, te = split['regime'], split['fold'], np.asarray(split['test'], dtype=np.int64)
                require(len(np.unique(te)) == len(te) and (te >= 0).all() and (te < len(X)).all(), 'Invalid held-out row indices')
                for family in families:
                    selected = saved.loc[saved.regime.eq(regime) & saved.fold.eq(fold) & saved.method.eq(family)]
                    seed = split.get('seed')
                    if seed is not None:
                        selected = selected.loc[selected.seed.eq(seed)]
                    require(np.array_equal(selected.row_index.to_numpy(), te), 'Saved and replayed held-out index order differs')
                    require(np.array_equal(selected.label.to_numpy(), frame.iloc[te].label.to_numpy()), 'Saved and replayed outcome labels differ')
                    if analysis == 'secondsource':
                        model_name = destination + '/fitted_models/{}_{}_{}_{}.joblib'.format(fold, seed, regime, family)
                        columns = {'phage_only': slice(0, 212), 'pair_only': slice(424, 448), 'full': slice(None)}[family]
                    elif analysis == 'nested_comparison':
                        model_name = destination + '/fitted_models/{}_{}_{}.joblib'.format(regime, fold, family)
                        columns = slice(None)
                    else:
                        model_name = destination + '/{}_{}/model.joblib'.format(regime, fold)
                        columns = slice(None)
                    model = joblib.load(inputs.file(model_name))
                    features = X[te, columns]
                    prediction = model.predict_proba(features)[:, 1]
                    require(len(prediction) == len(te), 'Model returned an unexpected score count')
                    raw_dtype = str(prediction.dtype)
                    detail = differences(selected.label.to_numpy(), selected.score.to_numpy(), prediction, te)
                    description = {'analysis': analysis, 'regime': regime, 'fold': fold, 'seed': seed, 'method': family,
                                   'artifact': model_name, 'artifact_sha256': inputs.checked[model_name]['sha256'],
                                   'input_dtype': str(features.dtype), 'prediction_dtype': raw_dtype}
                    if family.startswith('MLP'):
                        description.update(weight_dtypes=[str(value.dtype) for value in model['model'].coefs_],
                                           intercept_dtypes=[str(value.dtype) for value in model['model'].intercepts_],
                                           transformed_dtype=str(model['scale'].transform(features).dtype))
                    reports.append({**description, **detail})
                    key = (analysis, regime, family, seed)
                    pooled.setdefault(key, []).append((te, selected.label.to_numpy(), selected.score.to_numpy(), prediction))
        require(len(reports) == 174, 'Expected exactly 144 + 20 + 10 saved outer models')
        groups = []
        for key, parts in sorted(pooled.items(), key=lambda item: str(item[0])):
            rows, labels, old, new = (np.concatenate([part[column] for part in parts]) for column in range(4))
            require(len(np.unique(rows)) == len(rows), 'A pooled regime/model/seed repeats a held-out row')
            groups.append({'analysis': key[0], 'regime': key[1], 'method': key[2], 'seed': key[3],
                           **differences(labels, old, new, rows)})
    result = {'diagnostic_only': True, 'numerical_tolerance_not_changed': True, 'model_fits_started': 0, 'files_written': 0,
              'purpose': 'Report inference portability differences; no verification acceptance or scientific equivalence is asserted.',
              'runtime': environment, 'explicit_thread_limit': threads,
              'artifact_hashes_checked': len(inputs.checked), 'distinct_artifact_bytes': inputs.read_bytes,
              'completed_model_diagnostics': len(reports), 'completed_pooled_diagnostics': len(groups),
              'ece_definition': 'Ten equal-width [0,1] bins, with probability 1 included in the last bin; count-weighted absolute accuracy-confidence gap.',
              'pairwise_definition': 'Exact unordered score-order relations; changes are descriptive and not independent statistical tests.',
              'summary': {'max_absolute_score_error': max(row['max_absolute_error'] for row in reports),
                          'model_diagnostics_exceeding_historical_tolerance': sum(row['historical_rtol_1e_6_atol_1e_7_failing_rows'] > 0 for row in reports),
                          'historical_tolerance_failing_rows': sum(row['historical_rtol_1e_6_atol_1e_7_failing_rows'] for row in reports),
                          'threshold_ge_0_5_decision_changes': sum(row['threshold_ge_0_5_decision_changes'] for row in reports),
                          'fold_max_absolute_metric_differences': {metric: max(abs(row['metric_differences'][metric] or 0) for row in reports)
                                                                   for metric in reports[0]['metric_differences']},
                          'pooled_max_absolute_metric_differences': {metric: max(abs(row['metric_differences'][metric] or 0) for row in groups)
                                                                     for metric in reports[0]['metric_differences']}},
              'models': reports, 'pooled': groups}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--artifact-manifest', type=Path, help='Optional exact scientific extension inventory; default is the root manifest')
    parser.add_argument('--threads', type=int, choices=range(1, 65), help='Optional inference-only BLAS/OpenMP limit, reported explicitly')
    args = parser.parse_args()
    print(json.dumps(diagnose(args.root, args.artifact_manifest, args.threads), indent=2, allow_nan=False))
