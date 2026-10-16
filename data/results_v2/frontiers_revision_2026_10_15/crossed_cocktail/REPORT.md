# Retrospective crossed host-cluster × composition sensitivity

The historical host-only intervals remain unchanged. This additional analysis was planned after outcomes were known and is not preregistration.

Each draw independently resamples host-cluster and canonical-triplet labels; each original record receives the product of the two multiplicities. The same weights pair all methods. AP-minus-prevalence recomputes the prevalence in each draw. Exactly 2,000 two-class draws were retained per cohort.

| Cohort | Score | N | Compositions | AUROC [95% interval] | AP [95% interval] | AP minus prevalence [95% interval] |
| --- | --- | --- | --- | --- | --- | --- |
| historical_100 | pp_historical_product | 100 | 18 | 0.807 [0.543, 0.965] | 0.951 [0.852, 0.998] | 0.131 [0.024, 0.257] |
| historical_100 | pp_development_reference_product | 100 | 18 | 0.796 [0.528, 0.960] | 0.949 [0.844, 0.997] | 0.129 [0.023, 0.247] |
| historical_100 | pp_historical_max | 100 | 18 | 0.800 [0.522, 0.969] | 0.949 [0.845, 0.997] | 0.129 [0.020, 0.260] |
| historical_100 | pp_historical_mean | 100 | 18 | 0.797 [0.524, 0.965] | 0.947 [0.831, 0.998] | 0.127 [0.019, 0.255] |
| historical_100 | pp_logistic_raw_product | 100 | 18 | 0.788 [0.503, 0.974] | 0.944 [0.833, 0.998] | 0.124 [0.021, 0.254] |
| historical_100 | pp_knn_raw_product | 100 | 18 | 0.703 [0.441, 0.906] | 0.921 [0.764, 0.994] | 0.101 [0.007, 0.191] |
| historical_100 | bank_phage_prevalence_mean | 100 | 18 | 0.327 [0.112, 0.571] | 0.780 [0.628, 0.966] | -0.040 [-0.089, 0.048] |
| historical_100 | bank_phage_prevalence_product | 100 | 18 | 0.324 [0.110, 0.560] | 0.779 [0.627, 0.966] | -0.041 [-0.090, 0.047] |
| historical_100 | training_phylogroup_susceptibility | 100 | 18 | 0.705 [0.452, 0.944] | 0.893 [0.718, 0.996] | 0.073 [-0.008, 0.163] |
| exclude_AN03_99 | pp_historical_product | 99 | 18 | 0.805 [0.550, 0.973] | 0.950 [0.847, 0.999] | 0.131 [0.020, 0.247] |
| exclude_AN03_99 | pp_development_reference_product | 99 | 18 | 0.794 [0.540, 0.968] | 0.947 [0.840, 0.998] | 0.129 [0.020, 0.238] |
| exclude_AN03_99 | pp_historical_max | 99 | 18 | 0.798 [0.519, 0.975] | 0.948 [0.840, 0.998] | 0.129 [0.020, 0.252] |
| exclude_AN03_99 | pp_historical_mean | 99 | 18 | 0.795 [0.510, 0.971] | 0.945 [0.822, 0.999] | 0.127 [0.018, 0.243] |
| exclude_AN03_99 | pp_logistic_raw_product | 99 | 18 | 0.786 [0.497, 0.980] | 0.943 [0.836, 0.999] | 0.125 [0.019, 0.246] |
| exclude_AN03_99 | pp_knn_raw_product | 99 | 18 | 0.699 [0.444, 0.924] | 0.918 [0.765, 0.996] | 0.100 [0.006, 0.184] |
| exclude_AN03_99 | bank_phage_prevalence_mean | 99 | 18 | 0.322 [0.127, 0.596] | 0.777 [0.626, 0.972] | -0.041 [-0.093, 0.052] |
| exclude_AN03_99 | bank_phage_prevalence_product | 99 | 18 | 0.319 [0.123, 0.593] | 0.776 [0.625, 0.972] | -0.042 [-0.094, 0.052] |
| exclude_AN03_99 | training_phylogroup_susceptibility | 99 | 18 | 0.706 [0.462, 0.958] | 0.892 [0.712, 0.998] | 0.074 [-0.008, 0.166] |
| source_complete_64 | pp_historical_product | 64 | 14 | 0.795 [0.519, 1.000] | 0.967 [0.860, 1.000] | 0.092 [0.017, 0.250] |
| source_complete_64 | pp_development_reference_product | 64 | 14 | 0.797 [0.519, 1.000] | 0.968 [0.862, 1.000] | 0.093 [0.016, 0.253] |
| source_complete_64 | pp_historical_max | 64 | 14 | 0.756 [0.487, 1.000] | 0.961 [0.845, 1.000] | 0.086 [0.014, 0.236] |
| source_complete_64 | pp_historical_mean | 64 | 14 | 0.785 [0.494, 1.000] | 0.964 [0.839, 1.000] | 0.089 [0.015, 0.246] |
| source_complete_64 | pp_logistic_raw_product | 64 | 14 | 0.758 [0.475, 1.000] | 0.960 [0.845, 1.000] | 0.085 [0.016, 0.245] |
| source_complete_64 | pp_knn_raw_product | 64 | 14 | 0.718 [0.386, 1.000] | 0.948 [0.781, 1.000] | 0.073 [-0.004, 0.225] |
| source_complete_64 | bank_phage_prevalence_mean | 64 | 14 | 0.282 [0.038, 0.656] | 0.837 [0.658, 0.979] | -0.038 [-0.093, 0.072] |
| source_complete_64 | bank_phage_prevalence_product | 64 | 14 | 0.282 [0.038, 0.656] | 0.837 [0.658, 0.979] | -0.038 [-0.093, 0.072] |
| source_complete_64 | training_phylogroup_susceptibility | 64 | 14 | 0.802 [0.443, 1.000] | 0.956 [0.847, 1.000] | 0.081 [-0.003, 0.217] |
| source_complete_64 | source_released_product | 64 | 14 | 0.537 [0.030, 0.829] | 0.893 [0.658, 0.988] | 0.018 [-0.092, 0.120] |
| source_complete_exclude_AN03_63 | pp_historical_product | 63 | 14 | 0.791 [0.509, 1.000] | 0.966 [0.854, 1.000] | 0.093 [0.018, 0.240] |
| source_complete_exclude_AN03_63 | pp_development_reference_product | 63 | 14 | 0.793 [0.515, 1.000] | 0.966 [0.855, 1.000] | 0.093 [0.018, 0.243] |
| source_complete_exclude_AN03_63 | pp_historical_max | 63 | 14 | 0.751 [0.465, 1.000] | 0.960 [0.837, 1.000] | 0.087 [0.015, 0.222] |
| source_complete_exclude_AN03_63 | pp_historical_mean | 63 | 14 | 0.781 [0.480, 1.000] | 0.962 [0.835, 1.000] | 0.089 [0.016, 0.231] |
| source_complete_exclude_AN03_63 | pp_logistic_raw_product | 63 | 14 | 0.758 [0.473, 1.000] | 0.959 [0.848, 1.000] | 0.086 [0.017, 0.240] |
| source_complete_exclude_AN03_63 | pp_knn_raw_product | 63 | 14 | 0.713 [0.372, 1.000] | 0.945 [0.779, 1.000] | 0.072 [-0.012, 0.217] |
| source_complete_exclude_AN03_63 | bank_phage_prevalence_mean | 63 | 14 | 0.269 [0.009, 0.645] | 0.829 [0.643, 0.976] | -0.044 [-0.094, 0.051] |
| source_complete_exclude_AN03_63 | bank_phage_prevalence_product | 63 | 14 | 0.269 [0.009, 0.645] | 0.829 [0.643, 0.976] | -0.044 [-0.095, 0.051] |
| source_complete_exclude_AN03_63 | training_phylogroup_susceptibility | 63 | 14 | 0.801 [0.455, 1.000] | 0.955 [0.846, 1.000] | 0.082 [-0.003, 0.217] |
| source_complete_exclude_AN03_63 | source_released_product | 63 | 14 | 0.528 [0.000, 0.835] | 0.886 [0.667, 0.987] | 0.013 [-0.090, 0.122] |

## Audit of the existing composition shuffle

All 1,000 saved draws were independently reproduced. Each permutes 100 **whole triplet rows**, preserving every one of the 18 composition counts; no constituent is shuffled independently. The AUROC diagnostic central 95% range remains 0.603–0.749. This is not a formal permutation null or p-value.

## Limits

Additional post hoc dependence sensitivity, not a replacement for or original prespecified historical CI. Crossed resampling may be conservative/unstable with sparse, non-randomly selected cells and 18 compositions. It does not establish generalization to new banks, phages, selection procedures or composition populations. Shared phages across distinct compositions remain outside this grouping model. No multiplicity-adjusted confirmatory superiority test.

Only five of the 18 selected compositions have any negative outcomes, and three repeated compositions account for 72 of the 100 records. The non-random externally selected menu and sparse observed host/composition cells limit any composition-population interpretation. Distinct menus can share phages; grouping complete triplets does not eliminate all phage-level dependence. Host-only intervals are conditional on the saved menu, while these crossed intervals are an exploratory sensitivity to an additional grouping assumption. Neither is a guarantee of calibrated population coverage.

Exact source inputs, bootstrap row weights, factor multiplicities, per-method replicates, paired contrasts, full permutation indices and checksums are retained. No historical result has been overwritten.
