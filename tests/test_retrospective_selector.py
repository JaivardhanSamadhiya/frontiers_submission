import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("selector_audit", ROOT / "experiments/36_nested_corrected_selector.py")
selector = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(selector)


class SelectionDiagnosticTests(unittest.TestCase):
    def test_unobserved_assay_does_not_become_a_positive(self):
        # Even a stored value must not be counted when its observed mask is false.
        lower, upper, fraction, complete = selector.outcome(0, np.array([0, 1, 2]),
                                                            np.array([[0, 0, 1]]), np.array([[True, True, False]]))
        self.assertEqual((lower, upper, fraction, complete), (0, 1, 0, False))

    def test_observed_positive_resolves_union_despite_missing_assay(self):
        lower, upper, fraction, complete = selector.outcome(0, np.array([0, 1, 2]),
                                                            np.array([[1, 0, 0]]), np.array([[True, True, False]]))
        self.assertEqual((lower, upper, complete), (1, 1, False))
        self.assertAlmostEqual(fraction, 1/3)

    def test_diversity_rule_uses_genus_and_isolation_host(self):
        selected = selector.selection(np.array([.9, .8, .7, .6]),
                                      np.array(["A", "A", "A", "B"]), np.array(["X", "X", "Y", "X"]))
        np.testing.assert_array_equal(selected, [0, 2, 3])

    def test_query_batch_does_not_change_raw_nearest_neighbour(self):
        distance = np.array([[0, .1, .3, .2], [.1, 0, .2, .4], [.3, .2, 0, .1], [.2, .4, .1, 0]])
        y = np.array([[1, 0], [0, 1], [1, 1], [0, 0]])
        observed = np.ones_like(y, dtype=bool)
        first = selector.nearest_scores(np.array([0, 1]), np.array([2]), y, observed, distance, 2)
        both = selector.nearest_scores(np.array([0, 1]), np.array([2, 3]), y, observed, distance, 2)
        np.testing.assert_array_equal(first[0], both[0])


if __name__ == "__main__":
    unittest.main()
