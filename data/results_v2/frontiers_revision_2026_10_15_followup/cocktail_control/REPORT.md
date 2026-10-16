# Training phage susceptibility control audit

Arithmetic mean over the three chosen phages of each phage's positive fraction among non-missing individual Picard assays; zero source cocktail outcomes or target-host outcomes are used to calculate this score.

The 100 records contain 82 active cocktails and 18 distinct compositions. The score is constant within each composition. AUROC is 0.326897, with crossed host-cluster by composition bootstrap interval 0.111649–0.570748. The interval includes 0.5.

## Observed score distribution

| Observed cocktail outcome | Records | Mean score | Median score |
| --- | --- | --- | --- |
| 0 | 18 | 0.576422 | 0.584660 |
| 1 | 82 | 0.519805 | 0.569697 |

In this selected cohort, inactive cocktails have higher mean training susceptibility than active cocktails. The control has only one value per composition and no host-specific term. Its crossed interval includes chance; a robust inverse relation is not established.

The upstream README section 2.6 states that the rescue step fills unfilled slots using generalist phages with high overall lytic coverage after earlier stages recommend fewer than three phages. This is documentary evidence of the selection recipe, not evidence that the target hosts were intrinsically harder or that the recipe caused the observed inverse score association.

The released records describe selection stages and repeated compositions, but do not identify the causal reason for the inverse association. They do not establish that broadly active phages were deliberately assigned to harder hosts.

No control direction is reversed or retuned using the target outcomes. Stage and menu summaries are descriptive, not adjustment models or causal analyses.

## Related reporting checks

Remove this row from main Table 1 to match Methods, Table 2 and Figure 1A; preserve count audit in Data S1.

The LOSO and host-cluster zero-AUROC folds each contain the same nine Leuconostoc pseudomesenteroides rows (5 positive, 4 negative); all 20 positive-negative score comparisons are reversed. Small eligible folds permit extreme AUROC values; eligibility requires only three observations of each class. This is not a missing-class calculation.

All input file SHA-256 hashes were unchanged. Exact training counts, menu summaries, released selection-stage patterns, and zero-AUROC folds accompany this report.
