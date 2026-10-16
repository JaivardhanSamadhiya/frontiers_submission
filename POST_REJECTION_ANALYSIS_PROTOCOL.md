# PrecisionPhage retrospective strengthening protocol

Date: 2026-10-05. Existing internal and external outcomes have already been
examined. All analyses in this protocol are retrospective and exploratory.
They cannot be represented as a new blinded validation or a preregistration.
Original frozen inputs, predictions, protocols, and submitted manuscripts are
preserved. No hyperparameter search is allowed on ColoColi or StaphStudy.

## External evidence audit

Compare the frozen 100-cocktail score with maximum constituent, logistic-only,
nearest-neighbour-only, training generalist, and mean predicted susceptibility
across all 96 phages. Report AUROC and average precision with their paired
host-cluster bootstrap differences using 2,000 resamples. Conditional inference
is for the fixed bank. The 82% outcome prevalence is the AP reference.
Composition shuffling is a diagnostic control, not an exchangeability test.
Within-host pair ranking is restricted to the three externally selected phages.

## Matched evaluation cohort

Reuse the existing five dual cold-start blocks and the fixed 400-tree XGBoost
configuration. Predict exactly the same held-out pairs in each comparison.
Training pools exclude the held-out rows and then either allow shared entities,
exclude host clusters, exclude phage clusters, or exclude both cluster axes.
Match training positive and negative counts within each block to the minimum
available counts across pools. Three deterministic seeds measure training-sample
sensitivity; they are not independent datasets. Report pooled and block-specific
AUROC and AP, actual overlap fractions, cohort counts, and input checksums.
This addresses the previous confounding by different test rows, sample sizes,
and training prevalence. It does not fully isolate a causal effect of leakage.

## Corrected selection evaluation

Use raw distance-weighted nearest-neighbour scores that remain comparable
across candidate phages. The same-genus AND same-isolation-host diversity rule
remains fixed. Select neighbour count from 5, 10, 20, 40 using inner grouped
validation for at-least-one observed constituent activity. Use constituent
positive fraction only as a secondary tie-breaker, then smaller neighbour count.
Outer validation uses five host-grouped folds, with three inner grouped folds.
Compare with a training-prevalence generalist triplet on the same hosts.
Report missing selected outcomes separately and paired cluster-bootstrap
differences. Examine core-genome grouping thresholds 1e-4, 1e-3, 1e-2 without
selecting a favourable threshold. Skip and document thresholds with fewer than
five connected components. Each rule is picked without outer test labels.
The nearest-neighbour family itself was chosen after past development audits,
so the analysis remains exploratory despite nested hyperparameter selection.
At-least-one constituent activity is a union-coverage proxy; joint cocktail
synergy, antagonism, and measured cocktail efficacy are not assessed here.

## Publication boundary

A new positive result will be retained alongside negative results. It will not
replace the inconclusive frozen cross-source result or imply bank/taxon transfer.
Current literature must be used to identify incremental contributions rather
than claim that leakage control or prediction-based cocktail design is new.
No cloud compute, payment, journal submission, or public repository push is
authorized by this protocol.
