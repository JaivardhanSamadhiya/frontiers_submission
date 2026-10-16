import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from precisionphage.cocktail.scoring import DevelopmentRankNormalizer, combined_ranking_score


class BatchInvariantScoringTests(unittest.TestCase):
    def test_single_host_batch_and_unrelated_queries_match(self):
        normalizer = DevelopmentRankNormalizer().fit([[.1, .2], [.5, .4], [.9, .8]])
        one = np.array([[.6, .3]])
        batch = np.vstack([one, [0, 1], [.3, .9], [.7, .1]])
        np.testing.assert_array_equal(normalizer.transform(one)[0], normalizer.transform(batch)[0])
        np.testing.assert_array_equal(normalizer.reference_, [[.1, .2], [.5, .4], [.9, .8]])

    def test_ties_monotonicity_and_extremes(self):
        normalizer = DevelopmentRankNormalizer().fit([[.1], [.1], [.8]])
        values = normalizer.transform([[0], [.1], [.2], [.8], [1]]).ravel()
        self.assertTrue((np.diff(values) >= 0).all())
        self.assertGreater(values.min(), 0)
        self.assertLess(values.max(), 1)

    def test_invalid_dimensions_and_nan_rejected(self):
        normalizer = DevelopmentRankNormalizer().fit([[.2, .3]])
        with self.assertRaises(ValueError):
            normalizer.transform([[.2]])
        with self.assertRaises(ValueError):
            normalizer.transform([[np.nan, .2]])
        with self.assertRaises(RuntimeError):
            DevelopmentRankNormalizer().transform([[.2]])

    def test_known_combination_and_invalid_inputs(self):
        np.testing.assert_allclose(combined_ranking_score([[.5,.5,.5], [0,0,0], [1,0,0]]), [.875, 0, 1])
        with self.assertRaises(ValueError):
            combined_ranking_score([[1.2]])


if __name__ == "__main__":
    unittest.main()
