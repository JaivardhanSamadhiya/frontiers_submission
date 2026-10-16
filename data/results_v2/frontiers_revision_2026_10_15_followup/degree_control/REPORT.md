# Retrospective degree-control audit, 15 October 2026

## Status and scope

This was a retrospective diagnostic, not a new blind validation. The analysis plan was saved before the new scores were computed (`retrospective_plan.json`); the historical outcomes were already known. No sequence model was refit, no held-out label was used to construct a score, and no best variant was selected on the test results. The new controls preserve the five original split designs and the original held-out row order. All **13,248 historical degree-prediction rows** reproduce (score tolerance 1e-15).

Seven unit tests verify the fraction denominator, distinction between unseen and observed inactive entities, lack of test-label use, historical score construction, within-fold shift invariance, constant training-prior control, and rejection of empty training sets.

## Definition of the added control

For each fold and exact phage/host identifier, the positive-fraction score is the number of positive *training* pairs divided by the number of *all observed training* pairs for that identifier. An identifier absent from training receives the training-set prevalence. The additive control is the equally weighted average of the two entity fractions. This denominator-aware control is reported beside, not in place of, the original count-degree control.

The training prior alone and additive fraction minus that fold's prior are **pooling diagnostics**. The latter is not a newly selected predictor. It does not change within-fold rankings, but avoids treating distinct fallback constants from different folds as evidence of within-fold discrimination.

## Principal results

Values are AUROC except the AP columns. All comparisons are descriptive; no independent-row inference is performed.

| Regime | Pairs; prevalence | Count additive: pooled / macro-fold / within-fold pair-weighted | Fraction additive: pooled / macro-fold / within-fold pair-weighted | Fraction AP | Fraction centred diagnostic, pooled |
| --- | --- | --- | --- | --- | --- |
| Species holdout | 1057; .760643 | .042824 / .088781 / .027967 | .889830 / .873009 / .937886 | .927460 | .896565 |
| Host-cluster holdout | 1082; .763401 | .047464 / .090545 / .049592 | .889866 / .871926 / .938034 | .928701 | .896681 |
| Phage-cluster holdout | 384; .364583 | .612339 / .612227 / .639151 | .587573 / .616860 / .570231 | .429866 | .603088 |
| Dual exclusion | 398; .738693 | .500000 / .500000 / .500000 | .422390 / .500000 / .500000 | .706239 | .500000 |
| Genus holdout (supplementary) | 1495; .730435 | .149049 / .199429 / .115563 | .866516 / .864991 / .893652 | .913678 | .889686 |

The within-fold pair-weighted AUROC pools concordance only among positive-negative pairs belonging to the **same** fold, with weights n-positive times n-negative. It is not the usual pooled out-of-fold AUROC. Macro-fold AUROC assigns equal weight to each eligible fold. Neither is reported as a confidence interval.

## Why the old count score was below chance

The failure is **not merely a cross-fold pooling artifact**: all 28 species folds and all 26 host-cluster folds have original additive-count AUROC below .5, including nine exactly-zero folds in each regime. Their within-fold pair-weighted AUROCs are .027967 and .049592, respectively.

Training exposure is strongly associated with the held-out label:

| Regime / held-out label | Training-seen phages | Mean observed training degree | Mean positive training degree | Mean phage fraction score |
| --- | --- | --- | --- | --- |
| Species / negative | 253/253 (100%) | 11.6206 | 2.1502 | .2331 |
| Species / positive | 79/804 (9.83%) | .7960 | .1617 | .7089 |
| Host-cluster / negative | 256/256 (100%) | 11.5313 | 2.1172 | .2267 |
| Host-cluster / positive | 82/826 (9.93%) | .7809 | .1550 | .7077 |

Most held-out positive pairs therefore lose all exact-phage training exposure, while the held-out negatives involve phages observed against other hosts. The raw positive-count score rewards that exposure and assigns zero to most positives. Moreover, the exposed phages have multiple assayed rows, many negative; a positive count is not the same as an estimated susceptibility fraction. Adding the all-assayed denominator corrects this particular mismatch and yields a strong simple baseline in the permissive regimes. This is evidence of **label-associated sparse ascertainment/training-entity structure**, not proof of the original source's selection intentions or a biological interaction mechanism.

The strong positive-fraction control and the phage-only descriptor ablation reinforce the need to separate entity/base-rate signal from transferable interaction signal. The original count control must not be presented as having ruled out entity-level shortcuts. Nor does the corrected fraction baseline establish cross-source or cold-start transfer.

## Fold prior and pooling

Under dual exclusion, both identifiers are absent from each fold's training data. Every fraction score in that fold is consequently the same training prior, so each fold's AUROC is exactly .5. Pooling these different constants produces AUROC .422390. The training-prior-only control gives exactly the same result; prior centring restores a single constant and pooled AUROC .5. This is an actual cross-fold pooling artifact, distinct from the count-degree anti-ranking observed within the permissive folds.

The host-only fraction is likewise a constant prior within every species/host-cluster fold; pooled values .292490/.307383 are not meaningful evidence of discrimination. Report the fold-constant property, not those pooled values as biological failure or success.

## Exact-zero GBM folds

The saved full GBM has one exact-zero species fold: *Leuconostoc pseudomesenteroides*, nine test pairs (five positive, four negative). The corresponding host-cluster fold is cluster 130 with the same nine pairs. In those folds all 20 positive-negative comparisons are reversed. Thus the reported 0–1 fold ranges reflect small eligible folds, not a calculation bug. The original eligibility criterion was at least three positive and three negative pairs per test fold.

## Suggested manuscript wording

Methods:

> Exact-identifier controls were computed from training data only. In addition to historical positive interaction counts, I evaluated the fraction of all assayed training pairs that were positive for each phage and host, averaging the two entity scores with equal weights. Unseen identifiers received the training prevalence. I report pooled and individual-fold metrics because different training-prior fallbacks can distort pooled rankings when a score is constant within folds.

Results / interpretation:

> The original count-degree control was anti-predictive under species and host-cluster holdout (pooled AUROC .043 and .047), and this reversal persisted within folds. All held-out negative pairs involved training-seen phages, compared with only 9.8–9.9% of positive pairs, showing that raw positive counts confounded exposure with susceptibility in this sparse matrix. Replacing counts by training positive fractions yielded additive AUROC .890 in both regimes (.588 under phage-cluster holdout). Under dual exclusion, the control was constant within each fold and had within-fold AUROC .500; its pooled value .422 arose solely from differing training-prior constants. These controls support, rather than exclude, a substantial contribution from entity-level and sampling structure.

No mechanism of the Gaborieau source's cocktail-selection process is inferred from this NCBI_HR audit.

## Artifacts

- `predictions.csv`: every new/legacy score and saved GBM reference on unchanged rows.
- `fold_metrics.csv`, `summary.csv`: all fold and pooled AUROC/AP/prevalence values.
- `row_training_degrees.csv`: training observed/positive counts, exact-ID presence and fractions for each held-out record.
- `entity_train_test_counts.csv`: direct per-fold entity exposure audit.
- `seen_status_by_outcome.csv`: label-stratified exposure summary (reporting diagnostic only).
- `fold_score_distributions.csv`: fold/method/outcome score ranges and moments.
- `manifest.json`: reproduction count and output checksums.

Run from the project root with `python scripts/frontiers_revision/degree_control_followup.py`; test with `python -m unittest discover -s tests -p test_degree_control_followup.py -v`.
