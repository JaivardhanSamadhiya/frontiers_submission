"""Evaluation: metrics, calibration, bootstrap CIs, significance tests."""
from .metrics import (
    binary_metrics, aggregate_folds, bootstrap_ci, calibration_curve_ece,
)
from .significance import (
    benjamini_hochberg, bootstrap_auc_diff, delong_auc_ci, delong_test,
    mcnemar_test, permutation_auc_test,
)

__all__ = ["binary_metrics", "aggregate_folds", "bootstrap_ci",
           "calibration_curve_ece", "benjamini_hochberg", "bootstrap_auc_diff",
           "delong_auc_ci", "delong_test", "mcnemar_test", "permutation_auc_test"]
