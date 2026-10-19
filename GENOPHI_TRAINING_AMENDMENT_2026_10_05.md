# Compute-limited GenoPHI-derived comparator: training specification

Recorded before fitting the comparison. This adds implementation details to
GENOPHI_COMPARISON_PROTOCOL.md without replacing that record. All outcomes
are already public and have been analyzed before; this is retrospective.

## Name and budget

Use the label **GenoPHI-derived protein-family comparator (reduced budget)**.
This is not the authors' fitted model or an exact reproduction of their full
search. Use upstream GenoPHI's unmodified `perform_rfe`,
`compute_phage_weights`, and `build_row_weights` from commit
f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f. A separate driver handles exact folds,
Windows-compatible MMseqs2 calls, caching, and prediction. It must not import
or modify the unrelated yau_2026 checkout.

Primary grouping is patristic connected components at 0.001. Secondary
grouping is host identity. Use all five saved outer folds in each setting.
Within each outer training partition use three GroupKFold partitions for
hyperparameter selection. Rebuild protein references separately for each
inner training partition and the final outer training partition; validation
proteins are only queries, never reference construction inputs.

## Protein-family representation

Use MMseqs2 15-6f452, identity >=0.4, coverage >=0.8, sensitivity 7.5,
four CPU threads, and a 4 GB split-memory target. Give every protein a unique
entity-prefixed ID. Host reference clusters contain training hosts only;
known-bank phage references contain all 19 phages. No novel-phage claim is made.
Assign query proteins by the best qualifying alignment (highest bit score,
then lowest E-value, then stable target identifier) to training cluster
members. Retain host clusters present in at least two training host groups
and phage clusters present in at least two known-bank phages. Merge exactly
identical training presence patterns; an assigned member of any merged
cluster activates that feature. Remove features constant in fitting rows.
Do not use hash equality alone as proof of identical patterns.

## Supervised selection and fitting

For each fitting partition, perform three RFE repetitions on independently
seeded 80% samples of its host groups (seeds 20261005, 20261006, 20261007).
Use the published RFE base estimator: CatBoost 500 iterations, depth 4,
learning rate 0.1, approximately ten elimination steps. Select max(1,floor(P/20))
features per repetition, where P is the number of nonconstant input features.
Use training-only inverse-frequency positive class weights per phage, with
smoothing 1 and negative weight 1. Keep features having the maximum selection
frequency; cap their number at max(1,floor(P/20)), breaking ties by feature name.
If a subsample has only one class, record it and omit that repetition; if none
are usable, stop that fit rather than invent a trained classifier.

Tune the published eight-setting CatBoost grid (500/1000 iterations,
0.05/0.1 learning rate, depth 4/6) on mean inner-fold MCC at fixed threshold
0.5. Tie-break by fewer iterations, shallower depth, then smaller learning
rate. No outer-test early stopping or threshold tuning. Refit on all outer
training rows and take the median score of five seed fits (20261005–20261009).
Do not calibrate or interpret these weighted classifier scores as portable
probabilities. Preserve every grid result and selected feature list.

This reduces the published 25 feature-selection rounds and 50-model ensemble
to three and five, respectively. Grouped nested selection replaces the
published random/hierarchical repeated splitting, and host-group feature
presence replaces the reported 20-cluster Ward filter. Genome annotations
may differ from the original study. These are explicit adaptations; a loss
does not establish superiority over the full published GenoPHI pipeline.

## Evaluation and reporting

After prediction, evaluate all 437 held-out pairs per grouping, with paired
2,000 host-group bootstrap resamples (seed 20261005). Report AUROC and average
precision and their paired differences from saved nearest-neighbour and
prevalence scores. Rank phages within each host; take three with stable
identifier tie-breaking. Compute constituent-union coverage and paired
differences from saved neighbour and training-only greedy triplets. These
are single-phage assay unions, not measured mixture outcomes. Report exact
discordant host counts, all five fold outcomes and intervals regardless of
direction. Bootstrap intervals condition on fitted predictions and do not
include all training variation. With only 23 hosts/17 primary groups, avoid
strong claims from small or inconclusive differences.

No paid compute, external messages, public repository pushes, or submissions
are part of this run. If interrupted, retain checkpointed predictions but
do not summarize an incomplete subset as a completed comparison.
