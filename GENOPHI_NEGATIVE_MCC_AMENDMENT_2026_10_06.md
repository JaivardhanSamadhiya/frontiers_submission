# GenoPHI negative-MCC retention amendment — 6 October 2026

The full-budget run stopped in outer fold 0, cutoff 10, repetition 32. All eight
internal validation MCC values were negative (maximum -0.0254823596). The pinned
upstream classification grid search starts best_mcc at zero and accepts only
MCC >= best_mcc. Thus it returned no model. This was not a Windows access violation.
No outer test performance had been evaluated for the full-budget comparison.

Change exactly one classification-grid-search initialization from zero to negative
infinity. Apply this consistently to both feature-selection and modeling calls in
the remaining run. Preserve the upstream >= tie rule (last maximizing candidate),
all eight hyperparameters, seeds, splits, weights, early stopping, feature cutoffs,
25 feature-selection repetitions, and 50 modeling repetitions at every eligible
cutoff. Do not discard weak seeds, replace them with favorable seeds, reverse scores,
or use outer outcomes to guide the repair. Keep all finite negative-MCC candidates
eligible. Nonfinite/incomplete searches must still fail verification.

The original official source files stay byte-for-byte unchanged. A small, explicit
runtime adapter changes that function's initializer after checking its structure.
The result must be described as **full-budget adapted GenoPHI modeling with a
negative-MCC retention repair and Windows reference adapter**, not unmodified
official GenoPHI or an exact end-to-end reproduction.

Preserve the failed repetition, error log, old driver/auditor/continuation source,
and failure status before resuming. The 82 previously saved models are reusable
only after verifying complete grids, matching retained-model parameters and MCC,
and a nonnegative maximizing MCC. In such grids, initializing at negative infinity
rather than zero yields the same final maximizing candidate and tie decision.
The 25 accepted feature-selection runs likewise must show nonnegative selected
MCC and complete artifacts. A source/checkpoint manifest records this verification.

For retries, archive incomplete stage outputs and their logs before re-execution;
never treat an artifact without a successful checkpoint as a completed stage.
Preserve original run provenance and the earlier prediction-API amendment, and
record a separate hash-linked repair manifest. The independent final audit must
check every model repetition's full grid, selection rule, and matched predictions,
not just count model files. Poor validation scores remain in the record.

Completion still requires all five outer folds, all required repetitions, locked
outer predictions, and a passed independent matched-cohort audit. This repair
does not change the separate negative BASEL-ECOR result or the study's claim limits.
