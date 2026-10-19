"""Fit-free policy boundary/tamper tests; no artifacts or reports are written."""
import copy
import hashlib
import unittest

import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline

import prediction_portability_policy as p


def network(dtype=np.float32):
    estimator = MLPClassifier()
    estimator.coefs_ = [np.ones((1, 2), dtype=dtype), np.ones((2, 1), dtype=dtype)]
    estimator.intercepts_ = [np.zeros(2, dtype=dtype), np.zeros(1, dtype=dtype)]
    return Pipeline([('model', estimator)])


def detail(key, pooled=False, n=2):
    analysis, regime, method, seed = key[:4]
    result = {'analysis': analysis, 'regime': regime, 'method': method, 'seed': seed,
              'n': n, 'positive': n//2, 'max_absolute_error': 0., 'different_numeric_values': 0,
              'historical_rtol_1e_6_atol_1e_7_failing_rows': 0, 'independent_rtol_1e_7_atol_1e_7_failing_rows': 0,
              'threshold_ge_0_5_decision_changes': 0, 'threshold_gt_0_5_decision_changes': 0,
              'original_metrics': {'auroc': .8, 'average_precision': .9, 'brier': .1, 'ece_10_equal_width': .2},
              'replayed_metrics': {'auroc': .8, 'average_precision': .9, 'brier': .1, 'ece_10_equal_width': .2},
              'metric_differences': dict.fromkeys(p.METRIC_CAPS, 0.),
              'pairwise': {'unordered_pairs': n*(n-1)//2, 'changed_pairwise_relations': 0,
                           'order_reversals': 0, 'original_ties_broken': 0, 'new_ties_created': 0,
                           'positive_negative_order_reversals': 0, 'positive_negative_original_ties_broken': 0,
                           'positive_negative_new_ties_created': 0},
              'top_k_index_tiebroken': [{'k': k, 'symmetric_difference_count': 0}
                                      for k in sorted({min(n, k) for k in [1, 3, 5, 10, 20]})]}
    if not pooled:
        result.update(fold=key[4], artifact='models/'+hashlib.sha256(repr(key).encode()).hexdigest()+'.joblib',
                      artifact_sha256='a'*64, prediction_dtype='float32')
        if method.startswith('MLP'):
            result.update(weight_dtypes=['float32', 'float32'], intercept_dtypes=['float32', 'float32'])
    return result


def refresh_summary(report):
    report['summary'] = {'max_absolute_score_error': max(row['max_absolute_error'] for row in report['models']),
                         'historical_tolerance_failing_rows': sum(row['historical_rtol_1e_6_atol_1e_7_failing_rows'] for row in report['models']),
                         'model_diagnostics_exceeding_historical_tolerance': sum(row['historical_rtol_1e_6_atol_1e_7_failing_rows'] > 0 for row in report['models']),
                         'threshold_ge_0_5_decision_changes': 0}
    for label, field in [('fold', 'models'), ('pooled', 'pooled')]:
        report['summary'][label+'_max_absolute_metric_differences'] = {
            metric: max(abs(row['metric_differences'][metric]) for row in report[field]) for metric in p.METRIC_CAPS}


def diagnostic():
    report = {'diagnostic_only': True, 'numerical_tolerance_not_changed': True, 'model_fits_started': 0, 'files_written': 0,
              'artifact_hashes_checked': 188, 'completed_model_diagnostics': 174, 'completed_pooled_diagnostics': 42,
              'models': [detail(key) for key in sorted(p.expected_keys(), key=repr)],
              'pooled': [detail(key, pooled=True, n=8 if key[0]=='secondsource' else 10)
                         for key in sorted(p.expected_keys(True), key=repr)]}
    refresh_summary(report)
    return report


class PredictionPolicyTests(unittest.TestCase):
    def test_default_strict_unchanged_and_opt_in_boundary(self):
        observed = np.array([.25], dtype=np.float32)
        expected = observed.astype(np.float64)+p.FLOAT32_ABSOLUTE_CAP
        with self.assertRaises(AssertionError):
            p.assert_saved_predictions(expected, observed, network())
        checked = p.assert_saved_predictions(expected, observed, network(), portable_float32=True)
        self.assertTrue(checked['portable_float32_mlp_applied'])
        self.assertEqual(checked['atol'], 3.814697265625e-6)
        self.assertEqual(checked['rtol'], 0)
        with self.assertRaises(AssertionError):
            p.assert_saved_predictions(expected+np.finfo(np.float32).eps, observed, network(), portable_float32=True)

    def test_trees_and_nonpipeline_keep_strict(self):
        observed = np.array([.25], dtype=np.float32)
        expected = observed.astype(np.float64)+1e-6
        for estimator in [object(), Pipeline([('model', object())]), network()['model']]:
            with self.subTest(estimator=type(estimator)), self.assertRaises(AssertionError):
                p.assert_saved_predictions(expected, observed, estimator, portable_float32=True)

    def test_float64_weights_intercepts_or_predictions_not_relaxed(self):
        observed = np.array([.25], dtype=np.float32)
        expected = observed.astype(np.float64)+1e-6
        for estimator, replay in [(network(np.float64), observed), (network(), observed.astype(np.float64))]:
            with self.assertRaises(AssertionError):
                p.assert_saved_predictions(expected, replay, estimator, portable_float32=True)
        estimator = network()
        estimator['model'].intercepts_[0] = estimator['model'].intercepts_[0].astype(np.float64)
        with self.assertRaises(AssertionError):
            p.assert_saved_predictions(expected, observed, estimator, portable_float32=True)

    def test_shape_nonfinite_and_flag_guard(self):
        for expected, observed in [([.2], [[.2]]), ([np.nan], [np.nan]), ([], [])]:
            with self.assertRaises(AssertionError):
                p.assert_saved_predictions(expected, observed, network(), portable_float32=True)
        with self.assertRaises(AssertionError):
            p.assert_saved_predictions([.2], [.2], network(), portable_float32='true')

    def test_complete_diagnostic(self):
        checked = p.validate_portable_diagnostic(diagnostic())
        self.assertTrue(checked['passed'])
        self.assertEqual(checked['models_checked'], 174)
        self.assertEqual(checked['same_label_pooled_tie_changes'], [])

    def test_diagnostic_only_completeness_and_identity(self):
        for field, value in [('diagnostic_only', False), ('model_fits_started', 1), ('files_written', 1),
                             ('completed_model_diagnostics', 173), ('completed_pooled_diagnostics', 41)]:
            report = diagnostic()
            report[field] = value
            with self.subTest(field=field), self.assertRaises(AssertionError):
                p.validate_portable_diagnostic(report)
        report = diagnostic()
        report['models'][0] = copy.deepcopy(report['models'][1])
        with self.assertRaises(AssertionError):
            p.validate_portable_diagnostic(report)

    def test_diagnostic_score_cap_and_dtype_guard(self):
        for mutation in ['cap', 'weights', 'output']:
            report = diagnostic()
            row = next(row for row in report['models'] if row['method']=='MLP')
            if mutation == 'cap':
                row['max_absolute_error'] = p.FLOAT32_ABSOLUTE_CAP+1e-10
            elif mutation == 'weights':
                row['weight_dtypes'][0] = 'float64'
            else:
                row['prediction_dtype'] = 'float64'
            with self.subTest(mutation=mutation), self.assertRaises(AssertionError):
                p.validate_portable_diagnostic(report)

    def test_rank_tie_threshold_and_topk_failures(self):
        for mutation in ['reversal', 'opposite_tie', 'fold_tie', 'ge', 'gt', 'topk']:
            report = diagnostic()
            row = report['models'][0]
            if mutation == 'reversal':
                row['pairwise'].update(order_reversals=1, changed_pairwise_relations=1)
            elif mutation == 'opposite_tie':
                row['pairwise']['positive_negative_new_ties_created'] = 1
            elif mutation == 'fold_tie':
                row['pairwise'].update(new_ties_created=1, changed_pairwise_relations=1)
            elif mutation in ['ge', 'gt']:
                row['threshold_'+mutation+'_0_5_decision_changes'] = 1
            else:
                row['top_k_index_tiebroken'][0]['symmetric_difference_count'] = 2
            with self.subTest(mutation=mutation), self.assertRaises(AssertionError):
                p.validate_portable_diagnostic(report)

    def test_metric_caps_and_tampering(self):
        for metric, cap in p.METRIC_CAPS.items():
            report = diagnostic()
            row = report['models'][0]
            row['replayed_metrics'][metric] += cap*2
            row['metric_differences'][metric] = row['replayed_metrics'][metric]-row['original_metrics'][metric]
            with self.subTest(metric=metric), self.assertRaises(AssertionError):
                p.validate_portable_diagnostic(report)
        report = diagnostic()
        report['models'][0]['metric_differences']['auroc'] = 1e-13
        with self.assertRaisesRegex(AssertionError, 'tampered'):
            p.validate_portable_diagnostic(report)

    def test_same_label_pooled_tiny_tie_disclosed(self):
        report = diagnostic()
        pool = next(row for row in report['pooled'] if row['method']=='MLP_extended' and row['regime']=='host_cluster')
        model = next(row for row in report['models'] if row['method']=='MLP_extended' and row['regime']=='host_cluster')
        model['max_absolute_error'] = pool['max_absolute_error'] = 1e-6
        model['different_numeric_values'] = pool['different_numeric_values'] = 1
        pool['pairwise'].update(new_ties_created=1, changed_pairwise_relations=1)
        pool['replayed_metrics']['average_precision'] += 2e-8
        pool['metric_differences']['average_precision'] = pool['replayed_metrics']['average_precision']-pool['original_metrics']['average_precision']
        refresh_summary(report)
        checked = p.validate_portable_diagnostic(report)
        self.assertEqual(len(checked['same_label_pooled_tie_changes']), 1)
        self.assertEqual(checked['same_label_pooled_tie_changes'][0]['ties_created'], 1)

    def test_same_label_pooled_reversal_disclosed(self):
        report = diagnostic()
        pool = next(row for row in report['pooled'] if row['method']=='MLP_extended' and row['regime']=='host_cluster')
        model = next(row for row in report['models'] if row['method']=='MLP_extended' and row['regime']=='host_cluster')
        model['max_absolute_error'] = pool['max_absolute_error'] = 1e-6
        model['different_numeric_values'] = pool['different_numeric_values'] = 2
        pool['pairwise'].update(order_reversals=1, changed_pairwise_relations=1)
        refresh_summary(report)
        checked = p.validate_portable_diagnostic(report)
        self.assertEqual(checked['same_label_pooled_order_and_tie_changes'][0]['reversals'], 1)
        self.assertEqual(checked['same_label_pooled_tie_changes'], [])

    def test_opposite_label_pooled_reversal_rejected(self):
        report = diagnostic()
        pool = next(row for row in report['pooled'] if row['method']=='MLP_extended' and row['regime']=='host_cluster')
        pool['pairwise'].update(order_reversals=1, changed_pairwise_relations=1, positive_negative_order_reversals=1)
        with self.assertRaisesRegex(AssertionError, 'opposite-label'):
            p.validate_portable_diagnostic(report)

    def test_summary_and_pooled_count_tampering(self):
        report = diagnostic()
        report['summary']['max_absolute_score_error'] = 1e-7
        with self.assertRaises(AssertionError):
            p.validate_portable_diagnostic(report)
        report = diagnostic()
        report['pooled'][0]['positive'] += 1
        with self.assertRaisesRegex(AssertionError, 'cohort counts'):
            p.validate_portable_diagnostic(report)


if __name__ == '__main__':
    unittest.main()
