# Measured doublets and composition-conditioned host discrimination

Retrospective prioritization among three measured two-phage subsets of externally selected triplets. No arbitrary cocktail design, unseen cohort, new organism or therapeutic efficacy validated.

There are 300 measured doublets across 100 hosts, 224 active results, and 26 hosts with mixed binary outcomes. The uniform-choice expected activity is 0.747; the observed oracle bound is 0.840. Only 3 parent compositions contain both active and inactive triplets. Five compositions contain negatives, but two of those have only one host and no positive observations.

## Primary doublet decision

| method | resampling | metric | estimate | ci_low | ci_high | valid_replicates | total_replicates |
| --- | --- | --- | --- | --- | --- | --- | --- |
| repaired_product | host_cluster | activity | 0.780000 | 0.680840 | 0.869585 | 2000 | 2000 |
| repaired_product | host_cluster | gain_vs_uniform | 0.033333 | -0.029725 | 0.085719 | 2000 | 2000 |
| repaired_product | host_cluster | gain_vs_training_mean | 0.030000 | -0.067961 | 0.117647 | 2000 | 2000 |

All methods and paired comparisons are in doublet_decision_summary.csv. Ties receive equal probability across tied maxima. A host-only phylogroup control is identical to uniform choice for this menu. Choosing among equal-size subsets avoids the structural superset advantage of product aggregation. Any favorable result remains conditional on the original study having chosen the three parent phages.

## Conditional host discrimination

| method | metric | resampling | estimate | ci_low | ci_high | valid_replicates | total_replicates |
| --- | --- | --- | --- | --- | --- | --- | --- |
| repaired_product | pair_weighted_conditional_auroc | host_cluster | 0.720930 | 0.522914 | 0.875000 | 2000 | 2000 |
| repaired_product | composition_mean_auroc | host_cluster | 0.688889 | 0.464076 | 0.851010 | 2000 | 2000 |
| repaired_product | pair_weighted_conditional_auroc | crossed | 0.720930 | 0.250000 | 0.940000 | 1863 | 2000 |
| repaired_product | composition_mean_auroc | crossed | 0.688889 | 0.263158 | 0.941176 | 1863 | 2000 |

The pair-weighted result compares active and inactive hosts only when they share exactly the same triplet. The composition-mean result weights each unique informative composition equally. Only three compositions identify these results, containing 72 hosts and 301 positive-negative host pairs; these 301 comparisons are not independent observations. Crossed intervals additionally multiply row weights by resampled parent-composition counts; the macro summary averages the unique informative compositions retained in each replicate. The repaired pair-weighted host-cluster interval excludes chance, but its crossed interval includes chance and is wide. Therefore, uncertainty across compositions remains substantial.

Options share constituents within hosts and phages recur across different compositions. Keeping each host menu intact and crossing parent compositions does not create independent phage banks or resample phage identities. Inference is conditional on this fixed phage bank; the intervals do not establish transfer to new phages or banks.

The primary doublet gain intervals include zero against both uniform choice and the training susceptibility control. The experiment does not establish better choice of a measured doublet. Repaired product and mean select the same highest-scoring doublet here: with three constituent scores and equal-size pairs, both ordinarily choose the top two constituents. Thus their agreement is algebraic, not independent support for the aggregation formula.

Leave-one-composition-out and largest-three exclusions, AN03 sensitivity, all raw option scores, observed outcomes, selection probabilities and resampling weights are included. No ordinal efficacy optimization was performed. All input hashes were unchanged.