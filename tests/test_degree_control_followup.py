"""Training-only exact-entity control invariants."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

path = Path(__file__).resolve().parents[1] / "scripts/frontiers_revision/degree_control_followup.py"
spec = importlib.util.spec_from_file_location("degree_control_followup", path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class DegreeControlTests(unittest.TestCase):
    def setUp(self):
        self.train = pd.DataFrame({"phage": ["p1", "p1", "p2", "p3"],
                                   "host": ["h1", "h2", "h1", "h3"],
                                   "label": [1, 0, 1, 0]})
        self.test = pd.DataFrame({"phage": ["p1", "p3", "new"],
                                  "host": ["h1", "h3", "new"], "label": [1, 0, 1]})

    def test_fraction_denominator_uses_all_assayed_rows(self):
        scores, _ = mod.score_controls(self.train, self.test)
        np.testing.assert_array_equal(scores["phage_positive_fraction"], [.5, 0, .5])
        np.testing.assert_array_equal(scores["host_positive_fraction"], [1, 0, .5])
        np.testing.assert_array_equal(scores["additive_positive_fraction"], [.75, 0, .5])

    def test_unknown_is_prior_not_assayed_inactive(self):
        scores, details = mod.score_controls(self.train, self.test)
        self.assertEqual(scores["phage_positive_fraction"][1], 0)
        self.assertEqual(scores["phage_positive_fraction"][2], self.train.label.mean())
        np.testing.assert_array_equal(details.phage_seen, [True, True, False])

    def test_no_test_label_use(self):
        baseline, _ = mod.score_controls(self.train, self.test)
        other = self.test.drop(columns="label")
        repeated, _ = mod.score_controls(self.train, other)
        for key in baseline:
            np.testing.assert_array_equal(baseline[key], repeated[key])

    def test_legacy_positive_degree_reproduced(self):
        scores, _ = mod.score_controls(self.train, self.test)
        np.testing.assert_array_equal(scores["phage_positive_degree"], [1, 0, 0])
        np.testing.assert_array_equal(scores["host_positive_degree"], [2, 0, 0])
        np.testing.assert_array_equal(scores["additive_positive_degree"], [1, 0, 0])

    def test_shift_changes_no_within_fold_ranking(self):
        scores, _ = mod.score_controls(self.train, self.test)
        a = mod.metric_row(self.test.label, scores["additive_positive_fraction"])
        b = mod.metric_row(self.test.label, scores["centred_additive_positive_fraction_diagnostic"])
        self.assertEqual(a["auroc"], b["auroc"])
        self.assertEqual(a["average_precision"], b["average_precision"])

    def test_prior_constant_inside_one_fold(self):
        scores, _ = mod.score_controls(self.train, self.test)
        np.testing.assert_array_equal(scores["training_prior_diagnostic"], [.5, .5, .5])
        self.assertEqual(mod.metric_row(self.test.label, scores["training_prior_diagnostic"])["auroc"], .5)

    def test_empty_training_is_invalid(self):
        with self.assertRaises(ValueError):
            mod.score_controls(self.train.iloc[:0], self.test)


if __name__ == "__main__":
    unittest.main()
