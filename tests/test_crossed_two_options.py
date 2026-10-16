"""Tests for generic paired binary-policy metrics and crossed weights."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

p = Path(__file__).resolve().parents[1]/'scripts/frontiers_revision/crossed_two_options.py'
spec = importlib.util.spec_from_file_location('crossed_two_options', p)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class CrossedTwoOptionsTests(unittest.TestCase):
    def test_ties_use_actual_mean_outcome(self):
        # Both positive tie -> 1, both negative tie -> 0, discordant tie -> .5.
        vals, choices, observed = mod.policy_values([1, 0, 1, 0], [1, 0, 0, 1],
            [2, 0, 2, 0], [1, 0, 0, 1], [.5]*4, [.5]*4, np.ones(4))
        np.testing.assert_array_equal(choices, [.5]*4)
        np.testing.assert_array_equal(observed, [1, 0, .5, .5])
        self.assertEqual(vals[0, 0], .5)
        self.assertEqual(vals[0, 1], .5)
        self.assertEqual(vals[0, 8], 0)

    def test_perfect_choice_and_weighted_gains(self):
        vals, _, observed = mod.policy_values([1, 0, 1], [0, 1, 1],
            [2, 0, 1], [0, 2, 1], [1, 0, .5], [0, 1, .5], [[2, 1, 3]])
        np.testing.assert_array_equal(observed, [1, 1, 1])
        self.assertEqual(vals[0, 0], 1)
        self.assertEqual(vals[0, 1], 1)
        np.testing.assert_allclose(vals[0, 7], 1-4/6)

    def test_sampling_product_reproducibility(self):
        h, c = [0, 0, 1, 2], ['A', 'B', 'A', 'C']
        a = mod.sample_weights(h, c, [True, False, False, True], [True]*4, 100, 55)
        b = mod.sample_weights(h, c, [True, False, False, True], [True]*4, 100, 55)
        w, hc, cc, hg, cg, hi, ci, meta = a
        np.testing.assert_array_equal(w, hc[:, hi]*cc[:, ci])
        np.testing.assert_array_equal(w, b[0])
        self.assertTrue((hc.sum(axis=1) == len(hg)).all())
        self.assertTrue((cc.sum(axis=1) == len(cg)).all())
        self.assertTrue((w[:, [0, 3]].sum(axis=1) > 0).all())
        self.assertEqual(meta['valid'], 100)

    def test_weighted_equals_explicit_record_duplication(self):
        sy, by, sm, bm, ss, bs = map(np.asarray, [[1, 0, 1, 0], [0, 1, 1, 0],
            [2, 0, 3, 0], [0, 2, 2, 0], [1, 1, .5, .5], [0, 0, .5, .5]])
        weight = np.array([2, 3, 1, 2])
        a, _, _ = mod.policy_values(sy, by, sm, bm, ss, bs, weight)
        ix = np.repeat(np.arange(4), weight)
        b, _, _ = mod.policy_values(sy[ix], by[ix], sm[ix], bm[ix], ss[ix], bs[ix], np.ones(len(ix)))
        np.testing.assert_allclose(a, b, atol=1e-14, rtol=0)

    def test_missing_discordance_fails(self):
        with self.assertRaises(ValueError):
            mod.sample_weights([0, 1], ['A', 'B'], [False, False], [True, True], 10)


if __name__ == '__main__':
    unittest.main()
