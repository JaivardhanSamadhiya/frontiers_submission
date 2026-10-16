# Two-measured-options crossed dependence sensitivity

This is a post hoc uncertainty amendment, not a new outcome-locked analysis. No saved source or historical interval was overwritten.

Each replicate independently resamples whole host-cluster labels and whole selected-menu labels. Product multiplicities reweight paired measured outcomes for each host. The single fixed baseline menu is conditioned upon. All metrics and methods share 2,000 valid draws. Preference accuracy uses weighted discordant hosts only; policy activity and gains use all reweighted hosts. Tied choices use the actual mean of the two observed outcomes.

| Cohort | Method | Preference on 17 discordant hosts [95% interval] | Chosen activity [95% interval] | Fixed-baseline activity [95% interval] | Change versus fixed baseline [95% interval] |
| --- | --- | --- | --- | --- | --- |
| historical_100 | pp_historical_product | 0.765 [0.286, 1.000] | 0.860 [0.682, 1.000] | 0.810 [0.606, 0.948] | 0.050 [-0.039, 0.194] |
| historical_100 | pp_development_reference_product | 0.647 [0.143, 1.000] | 0.840 [0.659, 0.970] | 0.810 [0.606, 0.948] | 0.030 [-0.056, 0.136] |
| historical_100 | pp_logistic_raw_product | 0.765 [0.316, 1.000] | 0.860 [0.700, 0.980] | 0.810 [0.606, 0.948] | 0.050 [-0.040, 0.182] |
| historical_100 | pp_knn_raw_product | 0.765 [0.333, 1.000] | 0.860 [0.701, 0.980] | 0.810 [0.606, 0.948] | 0.050 [-0.041, 0.182] |
| historical_100 | bank_phage_prevalence_product | 0.412 [0.000, 0.875] | 0.800 [0.628, 0.954] | 0.810 [0.606, 0.948] | -0.010 [-0.157, 0.177] |
| historical_100 | training_phylogroup_susceptibility | 0.500 [0.500, 0.500] | 0.815 [0.663, 0.945] | 0.810 [0.606, 0.948] | 0.005 [-0.071, 0.110] |
| exclude_AN03_99 | pp_historical_product | 0.765 [0.300, 1.000] | 0.859 [0.681, 1.000] | 0.808 [0.600, 0.947] | 0.051 [-0.043, 0.204] |
| exclude_AN03_99 | pp_development_reference_product | 0.647 [0.143, 1.000] | 0.838 [0.643, 0.979] | 0.808 [0.600, 0.947] | 0.030 [-0.056, 0.144] |
| exclude_AN03_99 | pp_logistic_raw_product | 0.765 [0.333, 1.000] | 0.859 [0.687, 0.987] | 0.808 [0.600, 0.947] | 0.051 [-0.041, 0.186] |
| exclude_AN03_99 | pp_knn_raw_product | 0.765 [0.333, 1.000] | 0.859 [0.688, 0.987] | 0.808 [0.600, 0.947] | 0.051 [-0.043, 0.194] |
| exclude_AN03_99 | bank_phage_prevalence_product | 0.412 [0.000, 0.903] | 0.798 [0.619, 0.955] | 0.808 [0.600, 0.947] | -0.010 [-0.161, 0.178] |
| exclude_AN03_99 | training_phylogroup_susceptibility | 0.500 [0.500, 0.500] | 0.813 [0.663, 0.950] | 0.808 [0.600, 0.947] | 0.005 [-0.074, 0.113] |

## Interpretation

Sensitivity to repeated selected menus, conditional on the existing fixed baseline menu, fixed predictions and observed paired assay records. Not a new frozen validation, prospective laboratory policy, selection-mechanism test or composition-population confidence guarantee. Only 17 hosts have discordant binary outcomes. Sparse crossed cells, 18 non-random selected menus and phages shared across distinct menus limit inference. No confirmatory p-value, multiplicity-adjusted superiority claim or therapeutic interpretation.

The historical and repaired point estimates are unchanged. These wider dependence sensitivities do not establish improved selection over the fixed baseline. The fixed baseline is one composition used across every host; a baseline-composition population effect cannot be estimated from it. Selected-menu resampling does not capture phage-level dependence across distinct selected menus. Rejected zero-discordance draws mean intervals condition on an estimable preference denominator; counts are retained in the manifest.
