"""Explicit, opt-in float32 MLP inference replay policy; no fitting or writes.

Default prediction comparisons retain their original tolerances. The optional
32-epsilon absolute cap applies only to actual fitted float32 MLPClassifier
pipeline outputs. An independently generated complete diagnostic must also
preserve discrimination, decision boundaries, opposite-label order relations
and inspected top-k sets. Saved scientific predictions are never replaced.
The cap is an empirical engineering policy selected after inspecting platform
diagnostics, not a formal roundoff bound or a statistical equivalence test.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import re

import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline

FLOAT32_ABSOLUTE_CAP = float(32 * np.finfo(np.float32).eps)
AUROC_CAP = 1e-12
AVERAGE_PRECISION_CAP = 1e-7
METRIC_CAPS = {'auroc': AUROC_CAP, 'average_precision': AVERAGE_PRECISION_CAP,
               'brier': FLOAT32_ABSOLUTE_CAP, 'ece_10_equal_width': FLOAT32_ABSOLUTE_CAP}


def require(value, message):
    if not value:
        raise AssertionError(message)


def policy():
    return {'name': 'explicit_float32_mlp_inference_replay_v1',
            'selection_basis': 'Post-diagnostic empirical engineering cap; not a formal floating-point error bound or statistical equivalence test',
            'default_strict_tolerances_unchanged': True,
            'portable_model': 'Pipeline named model step of exact MLPClassifier type; all coefficients/intercepts and observed scores float32',
            'portable_absolute_cap': FLOAT32_ABSOLUTE_CAP, 'portable_relative_tolerance': 0,
            'metric_difference_caps': dict(METRIC_CAPS),
            'thresholds': ['>=0.5', '>0.5'], 'opposite_label_order_and_tie_changes_allowed': 0,
            'top_k_symmetric_difference_allowed': 0,
            'same_label_pooled_tie_policy': 'Allowed only within score/metric caps; every change is counted and disclosed. Fold order/ties and pooled strict order reversals must remain unchanged.',
            'limitations': 'Inference reproducibility policy, not a new fit, recalibration, positive finding, metric replacement, or external validation.'}


def float32_mlp(model, observed):
    if not isinstance(model, Pipeline) or type(model.named_steps.get('model')) is not MLPClassifier:
        return False
    network = model.named_steps['model']
    weights = getattr(network, 'coefs_', None)
    intercepts = getattr(network, 'intercepts_', None)
    if not isinstance(weights, (list, tuple)) or not isinstance(intercepts, (list, tuple)):
        return False
    return (bool(weights) and len(weights) == len(intercepts)
            and all(isinstance(value, np.ndarray) and value.dtype == np.dtype('float32')
                    and value.size > 0 and np.isfinite(value).all() for value in [*weights, *intercepts])
            and observed.dtype == np.dtype('float32'))


def assert_saved_predictions(expected, observed, model, portable_float32=False,
                             strict_rtol=1e-6, strict_atol=1e-7):
    """Compare saved probabilities, keeping strict behavior unless opted in.

    Trees, non-pipelines, float64 networks and float64 replay outputs NEVER
    receive a relaxed comparison, even when the flag is true.
    """
    require(type(portable_float32) is bool, 'Portable inference must be explicitly boolean')
    expected, observed = np.asarray(expected), np.asarray(observed)
    require(expected.shape == observed.shape and expected.ndim == 1 and expected.size > 0,
            'Saved/replayed prediction vectors must have identical nonempty shapes')
    require(np.isfinite(expected).all() and np.isfinite(observed).all(), 'Non-finite prediction refused')
    require(np.isfinite(strict_rtol) and np.isfinite(strict_atol) and strict_rtol >= 0 and strict_atol >= 0,
            'Invalid original strict tolerances')
    apply_portable = portable_float32 and float32_mlp(model, observed)
    rtol, atol = (0, FLOAT32_ABSOLUTE_CAP) if apply_portable else (strict_rtol, strict_atol)
    np.testing.assert_allclose(expected, observed, rtol=rtol, atol=atol)
    return {'passed': True, 'portable_float32_mlp_applied': bool(apply_portable),
            'rtol': float(rtol), 'atol': float(atol),
            'max_absolute_error': float(np.max(np.abs(expected.astype(np.float64)-observed.astype(np.float64)))),
            'predictions_replaced': False}


def integer(value, name, minimum=0):
    require(type(value) is int and value >= minimum, 'Invalid diagnostic integer: '+name)
    return value


def finite(value, name):
    require(type(value) in (int, float) and np.isfinite(value), 'Invalid diagnostic number: '+name)
    return float(value)


def expected_keys(pooled=False):
    result = set()
    for regime, method, seed in itertools.product(
            ['shared_entities_allowed', 'host_cluster_excluded', 'phage_cluster_excluded', 'both_clusters_excluded'],
            ['phage_only', 'pair_only', 'full'], [20261016, 20261017, 20261018]):
        base = ('secondsource', regime, method, seed)
        result.update([base] if pooled else [base+(fold,) for fold in range(4)])
    for regime in ['host_cluster', 'dual_cluster']:
        for analysis, methods in [('nested_comparison', ['GBM', 'MLP']), ('mlp_budget_sensitivity', ['MLP_extended'])]:
            for method in methods:
                base = (analysis, regime, method, None)
                result.update([base] if pooled else [base+(fold,) for fold in range(5)])
    return result


def row_key(row, pooled=False):
    base = (row.get('analysis'), row.get('regime'), row.get('method'), row.get('seed'))
    return base if pooled else base+(row.get('fold'),)


def check_detail(row, pooled):
    n = integer(row.get('n'), 'n', 1)
    positive = integer(row.get('positive'), 'positive', 1)
    require(positive < n, 'Both classes must be represented in every diagnostic')
    error = finite(row.get('max_absolute_error'), 'max_absolute_error')
    require(error >= 0, 'Negative absolute error')
    mlp = row.get('method') in {'MLP', 'MLP_extended'}
    if mlp:
        require(error <= FLOAT32_ABSOLUTE_CAP, 'Float32 MLP absolute score cap exceeded')
        if not pooled:
            require(row.get('prediction_dtype') == 'float32', 'Portable MLP replay outputs must be float32')
            weights, intercepts = row.get('weight_dtypes'), row.get('intercept_dtypes')
            require(isinstance(weights, list) and isinstance(intercepts, list) and len(weights) > 0
                    and len(weights) == len(intercepts) and set(weights+intercepts) == {'float32'},
                    'Portable MLP weights and intercepts must all be float32')
    else:
        require(integer(row.get('historical_rtol_1e_6_atol_1e_7_failing_rows'), 'strict failing rows') == 0
                and integer(row.get('independent_rtol_1e_7_atol_1e_7_failing_rows'), 'independent failing rows') == 0,
                'Tree/non-MLP predictions must retain both original strict checks')
    for field in ['historical_rtol_1e_6_atol_1e_7_failing_rows', 'independent_rtol_1e_7_atol_1e_7_failing_rows', 'different_numeric_values']:
        require(integer(row.get(field), field) <= n, 'Diagnostic row count exceeds its cohort')
    for field in ['threshold_ge_0_5_decision_changes', 'threshold_gt_0_5_decision_changes']:
        require(integer(row.get(field), field) == 0, 'A 0.5 decision changed')
    deltas = {}
    for metric, cap in METRIC_CAPS.items():
        before = finite(row.get('original_metrics', {}).get(metric), 'original '+metric)
        after = finite(row.get('replayed_metrics', {}).get(metric), 'replayed '+metric)
        require(0 <= before <= 1 and 0 <= after <= 1, 'Metric outside [0,1]')
        difference = finite(row.get('metric_differences', {}).get(metric), 'difference '+metric)
        require(abs(difference-(after-before)) <= 1e-15, 'Reported metric difference was tampered with')
        require(abs(difference) <= cap, 'Diagnostic metric cap exceeded: '+metric)
        deltas[metric] = abs(difference)
    pairs = row.get('pairwise', {})
    require(integer(pairs.get('unordered_pairs'), 'unordered pairs') == n*(n-1)//2, 'Wrong unordered-pair population')
    changed = sum(integer(pairs.get(key), key) for key in ['order_reversals', 'original_ties_broken', 'new_ties_created'])
    require(integer(pairs.get('changed_pairwise_relations'), 'changed relations') == changed
            and changed <= pairs['unordered_pairs'], 'Pairwise change counts do not reconcile')
    for key in ['positive_negative_order_reversals', 'positive_negative_original_ties_broken', 'positive_negative_new_ties_created']:
        require(integer(pairs.get(key), key) == 0, 'An opposite-label order/tie relation changed')
    require(pairs['order_reversals'] == 0, 'A strict score-order relation reversed')
    if not pooled:
        require(changed == 0, 'Fold score-order/tie relations must remain unchanged')
    top = row.get('top_k_index_tiebroken')
    ks = sorted({min(n, k) for k in [1, 3, 5, 10, 20]})
    require(isinstance(top, list) and [entry.get('k') for entry in top] == ks, 'Missing or repeated top-k diagnostic')
    require(all(integer(entry.get('symmetric_difference_count'), 'top-k set changes') == 0 for entry in top),
            'A deterministic top-k set changed')
    return {'n': n, 'positive': positive, 'score_error': error, 'metric_errors': deltas,
            'same_label_ties_broken': pairs['original_ties_broken'], 'same_label_ties_created': pairs['new_ties_created']}


def validate_portable_diagnostic(report):
    """Validate the full 174-model/42-pool diagnostic, not its summary alone."""
    require(isinstance(report, dict) and report.get('diagnostic_only') is True
            and report.get('numerical_tolerance_not_changed') is True, 'Require the unchanged diagnostic-only report')
    require(integer(report.get('model_fits_started'), 'model fits') == integer(report.get('files_written'), 'files written') == 0,
            'Diagnostic must start no fits and write no files')
    require(integer(report.get('completed_model_diagnostics'), 'model diagnostics') == 174
            and integer(report.get('completed_pooled_diagnostics'), 'pooled diagnostics') == 42,
            'All 174 models and 42 pooled cohorts are required')
    require(integer(report.get('artifact_hashes_checked'), 'artifact hashes') >= 174, 'Insufficient hashed artifact checks')
    models, pools = report.get('models'), report.get('pooled')
    require(isinstance(models, list) and len(models) == 174 and isinstance(pools, list) and len(pools) == 42,
            'Full model and pooled diagnostic rows are required')
    require({row_key(row) for row in models} == expected_keys()
            and {row_key(row, True) for row in pools} == expected_keys(True), 'Diagnostic model/fold/regime/seed identities differ')
    artifacts = [row.get('artifact') for row in models]
    require(len(set(artifacts)) == 174 and all(isinstance(name, str) and name and not name.startswith('/')
            and '\\' not in name and ':' not in name and '..' not in name.split('/') for name in artifacts), 'Invalid/repeated artifact paths')
    require(all(re.fullmatch('[0-9a-f]{64}', row.get('artifact_sha256', '')) for row in models), 'Missing model artifact hashes')
    checked_models = [check_detail(row, False) for row in models]
    checked_pools = [check_detail(row, True) for row in pools]
    for pool, checked in zip(pools, checked_pools):
        members = [item for row, item in zip(models, checked_models) if row_key(row, True) == row_key(pool, True)]
        require(sum(item['n'] for item in members) == checked['n']
                and sum(item['positive'] for item in members) == checked['positive'], 'Pooled and fold cohort counts differ')
        require(abs(max(item['score_error'] for item in members)-checked['score_error']) <= 1e-15,
                'Pooled and fold maximum score errors differ')
    summary = report.get('summary', {})
    error = max(row['score_error'] for row in checked_models)
    require(abs(finite(summary.get('max_absolute_score_error'), 'summary score error')-error) <= 1e-15, 'Summary conceals a score difference')
    for label, rows in [('fold', checked_models), ('pooled', checked_pools)]:
        for metric in METRIC_CAPS:
            maximum = max(row['metric_errors'][metric] for row in rows)
            require(abs(finite(summary.get(label+'_max_absolute_metric_differences', {}).get(metric), 'summary '+metric)-maximum) <= 1e-15,
                    'Summary metric differences do not match the detailed rows')
    old_failing = sum(row['historical_rtol_1e_6_atol_1e_7_failing_rows'] for row in models)
    old_models = sum(row['historical_rtol_1e_6_atol_1e_7_failing_rows'] > 0 for row in models)
    require(integer(summary.get('historical_tolerance_failing_rows'), 'summary historical failing rows') == old_failing
            and integer(summary.get('model_diagnostics_exceeding_historical_tolerance'), 'summary historical failing models') == old_models
            and integer(summary.get('threshold_ge_0_5_decision_changes'), 'summary threshold changes') == 0,
            'Summary historical tolerance/decision counts differ from the detailed rows')
    same_label = [{'analysis': row['analysis'], 'regime': row['regime'], 'method': row['method'], 'seed': row['seed'],
                   'ties_broken': checked['same_label_ties_broken'], 'ties_created': checked['same_label_ties_created']}
                  for row, checked in zip(pools, checked_pools) if checked['same_label_ties_broken'] or checked['same_label_ties_created']]
    return {'passed': True, 'scope': 'Complete saved-model inference portability diagnostic only; original saved metrics and predictions retained',
            'models_checked': 174, 'pooled_cohorts_checked': 42, 'max_absolute_score_error': error,
            'model_diagnostics_exceeding_historical_tolerance': old_models,
            'historical_tolerance_failing_rows': old_failing,
            'fold_metric_max_absolute_differences': {metric: max(row['metric_errors'][metric] for row in checked_models) for metric in METRIC_CAPS},
            'pooled_metric_max_absolute_differences': {metric: max(row['metric_errors'][metric] for row in checked_pools) for metric in METRIC_CAPS},
            'same_label_pooled_tie_changes': same_label,
            'diagnostic_canonical_json_sha256': hashlib.sha256(json.dumps(report, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
            'model_fits_started': 0, 'files_written': 0, 'policy': policy()}
