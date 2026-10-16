"""Generic saved-score metric/resampling tests; no biological model execution."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

path = Path(__file__).resolve().parents[1] / 'scripts/frontiers_revision/crossed_cocktail_uncertainty.py'
spec = importlib.util.spec_from_file_location('crossed_cocktail', path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class CrossedCocktailTests(unittest.TestCase):
    def test_canonical_groups(self):
        self.assertEqual(mod.canonical_composition('C|A|B'), 'A|B|C')
        with self.assertRaises(ValueError):
            mod.canonical_composition('A|A|B')

    def test_weighted_ties_against_sklearn(self):
        rng = np.random.default_rng(791)
        for _ in range(50):
            y = rng.integers(0, 2, 30)
            y[0], y[1] = 0, 1
            score = rng.integers(0, 6, 30)/5
            weight = rng.integers(0, 5, (4, 30))
            weight[:, :2] = 1
            actual = mod.weighted_metrics_many(y, score, weight)
            for index, w in enumerate(weight):
                ap = average_precision_score(y, score, sample_weight=w)
                prev = np.average(y, weights=w)
                expected = [roc_auc_score(y, score, sample_weight=w), ap, prev, ap-prev]
                np.testing.assert_allclose(actual[index], expected, rtol=0, atol=1e-13)

    def test_factor_product_and_paired_draws(self):
        y = np.array([0, 1, 1, 0, 1, 1])
        h = np.array([0, 0, 1, 2, 2, 3])
        c = np.array(['A', 'B', 'A', 'A', 'B', 'C'])
        a = mod.crossed_weights(y, h, c, replicates=200, seed=42)
        b = mod.crossed_weights(y, h, c, replicates=200, seed=42)
        w, hc, cc, hg, cg, hi, ci, attempts = a
        np.testing.assert_array_equal(w, hc[:, hi]*cc[:, ci])
        np.testing.assert_array_equal(w, b[0])
        self.assertTrue((hc.sum(axis=1) == len(hg)).all())
        self.assertTrue((cc.sum(axis=1) == len(cg)).all())
        self.assertEqual(w.shape, (200, 6))
        self.assertGreaterEqual(attempts, 200)
        self.assertTrue(((w*y).sum(axis=1) > 0).all())
        self.assertTrue(((w*(1-y)).sum(axis=1) > 0).all())

    def test_constant_score_ap_excess_is_zero(self):
        m = mod.weighted_metrics_many([0, 1, 1], [.4, .4, .4], [[1, 2, 3]])
        np.testing.assert_allclose(m[0], [.5, 5/6, 5/6, 0], rtol=0, atol=1e-14)

    def test_one_class_fails(self):
        with self.assertRaises(ValueError):
            mod.weighted_metrics_many([0, 1], [0, 1], [[1, 0]])
        with self.assertRaises(ValueError):
            mod.crossed_weights([1, 1], [0, 1], ['A', 'B'], replicates=1)


if __name__ == '__main__':
    unittest.main()
