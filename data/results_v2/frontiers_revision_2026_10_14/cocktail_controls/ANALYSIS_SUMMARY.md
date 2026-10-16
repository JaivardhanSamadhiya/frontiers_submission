# Retrospective cocktail controls for the Frontiers revision

This analysis reuses saved scores and measured outcomes. No predictor was refitted. Results and intervals condition on the fixed phage bank and the externally chosen composition menu.

## Released original-source scores on exactly common hosts

The source file contains 233 numeric decision scores and 67 missing rescue scores. Exactly 64 hosts (56 lytic and 8 non-lytic selected cocktails) have three numeric original scores; excluding ambiguous AN03 leaves 63. The same product formula is applied to the source decision scores, without pretending that the source publication reported a cocktail probability. All methods below use identical hosts.

| cohort | score | hosts | AUROC | AP | prevalence | AP minus prevalence |
| --- | --- | --- | --- | --- | --- | --- |
| source_complete_historical_100 | pp_historical_product | 64 | 0.795 [0.622, 0.940] | 0.967 | 0.875 | 0.092 [0.034, 0.168] |
| source_complete_historical_100 | pp_development_reference_product | 64 | 0.797 [0.625, 0.942] | 0.968 | 0.875 | 0.093 [0.035, 0.167] |
| source_complete_historical_100 | source_released_product | 64 | 0.537 [0.287, 0.801] | 0.893 | 0.875 | 0.018 [-0.049, 0.088] |
| source_complete_exclude_AN03 | pp_historical_product | 63 | 0.791 [0.619, 0.939] | 0.966 | 0.873 | 0.093 [0.034, 0.170] |
| source_complete_exclude_AN03 | pp_development_reference_product | 63 | 0.793 [0.617, 0.941] | 0.966 | 0.873 | 0.093 [0.035, 0.168] |
| source_complete_exclude_AN03 | source_released_product | 63 | 0.528 [0.277, 0.798] | 0.886 | 0.873 | 0.013 [-0.051, 0.086] |

On the 64 complete hosts, the historical score's paired AUROC difference from the released-source aggregate is 0.258 (95% host-group interval -0.042 to 0.563); the AP difference is 0.074 (interval -0.001 to 0.179). Both intervals include zero. The matched point estimates therefore do not establish superiority, and missing rescue scores prevent an equivalent point comparison over the whole cohort.

All-cohort identification bounds are also supplied: unknown rescue scores are allowed anywhere in [0,1]. The worst assignment takes positive cases to their lowest feasible cocktail score and negative cases to their highest; the best assignment reverses those choices. Because these assignments use known outcomes, they are bounds on an incompletely observed source comparator, not predicted scores or confidence intervals.

| cohort | case | hosts | AUROC | AP |
| --- | --- | --- | --- | --- |
| historical_100 | worst | 100 | 0.173 | 0.682 |
| historical_100 | best | 100 | 0.859 | 0.961 |
| exclude_AN03_99 | worst | 99 | 0.166 | 0.677 |
| exclude_AN03_99 | best | 99 | 0.858 | 0.959 |

## Choice between two experimentally measured triplets for a host

The source measured both the selected triplet and the fixed baseline LF82_P8/CLB_P2/LF73_P1 for all 100 hosts. Selected outcomes are 82/100 lytic and baseline outcomes 81/100: 73 both lytic, 10 neither, 9 selected-only and 8 baseline-only. Preference accuracy uses only the 17 binary-discordant hosts; score ties count one half. Chosen activity uses all hosts, with uniform random choice for score ties. These are retrospective comparisons within this two-option menu.

| cohort | score | discordant | preference accuracy | chosen activity | change vs baseline |
| --- | --- | --- | --- | --- | --- |
| historical_100 | pp_historical_product | 17 | 0.765 [0.556, 0.944] | 0.860 | 0.050 [0.000, 0.110] |
| historical_100 | pp_development_reference_product | 17 | 0.647 [0.412, 0.875] | 0.840 | 0.030 [-0.020, 0.083] |
| historical_100 | bank_phage_prevalence_product | 17 | 0.412 [0.176, 0.667] | 0.800 | -0.010 [-0.084, 0.070] |
| exclude_AN03_99 | pp_historical_product | 17 | 0.765 [0.556, 0.944] | 0.859 | 0.051 [0.000, 0.111] |
| exclude_AN03_99 | pp_development_reference_product | 17 | 0.647 [0.412, 0.875] | 0.838 | 0.030 [-0.020, 0.084] |
| exclude_AN03_99 | bank_phage_prevalence_product | 17 | 0.412 [0.176, 0.667] | 0.798 | -0.010 [-0.085, 0.071] |

There are 18 unique selected triplets and one unique baseline triplet. Canonical host/composition keys yield 200 unique measured triplet records. Repeated compositions on different hosts are retained; their outcomes are different assays. Single and double mixture outcomes are deliberately unused because the individual baseline b1/b2/b3 order has not been independently established.

## Interpretation and reproducibility

Host-group percentile intervals use 2,000 resamples. The AP-minus-prevalence reference is recomputed within each resample. Both historical test-batch percentiles and the later development-reference repair are reported, rather than replacing the historical result. The complete-source subset is selection-biased toward non-rescue recommendations; it cannot support an unqualified whole-cohort superiority claim.

Source: [Gaborieau et al. (2024)](https://doi.org/10.1038/s41564-024-01832-5), Methods ‘Assembling phages in cocktails’ and [public code](https://github.com/mdmparis/coli_phage_interactions_2023/tree/99656d8c020fcd1dbbc30a02ee10274a07655562). Endpoint: observed lytic activity (MLC > 0), which is not proof of productive replication at every titre.

- Every analysis is retrospective; external outcomes were known before these controls were defined.
- No predictor was refitted. None of these comparisons is a new independent frozen validation.
- The original source records 233 of 300 selected constituent decision scores; 67 rescue scores are missing, on 36 hosts.
- The 64 complete hosts are enriched for selections made before rescue. Complete-case source comparison may not represent all 100 hosts.
- Source y_pred_proba comes from several recommender stages. Its product aggregation is imposed here for comparison and is not a published end-to-end cocktail probability.
- Per-host lower/upper source-score columns bound missing constituent scores. Full-cohort identification bounds set unknown scores to outcome-dependent extremes solely to bound AUROC/AP. These bounds are not model predictions, confidence intervals, or substitutes for missing released scores.
- Historical PrecisionPhage percentiles depend on the external query cohort; the later development-reference repair is reported separately.
- Intervals resample entire historical host sequence groups, conditional on existing predictions and the fixed bank/composition menu. They do not generalize over new phages, banks, or selection procedures.
- AN03 labels refer to different development/external genome identifiers, while the published distance table has one AN03. The 99-host and 63-complete-host sensitivity excludes that ambiguity; it does not silently replace historical estimates.
- Only 17 hosts have discordant binary outcomes between the two measured triplets. Within-host binary preference accuracy is estimated on those hosts with ties worth one half.
- Repeated compositions across different hosts are legitimate distinct assays and are retained. Canonicalization removes ordering differences; only exact duplicate host/composition records are deduplicated.
- MLC > 0 is observed lytic activity; high-titre clearing does not establish productive phage replication, synergy, clinical benefit, or resistance suppression.
- Only the unordered fixed-baseline triplet is assigned. Published CLB_P2 maps to bank 55989_P2 by exact identity of the complete ENA OZ035781.1 sequence plus a 60-base terminal duplication of its prefix in the local assembly. No identity is inferred for b1/b2/b3 individual phages or two-phage baseline mixtures.
- AP-minus-prevalence intervals recompute the prevalence baseline in each identical bootstrap replicate. A finite-sample permutation AP reference can differ from prevalence.
- No permutation significance test, exchangeability claim, or confirmatory p-value is supplied.
- The ordinal MLC preference diagnostic treats the source scale as ordered, not equally spaced biological efficacy.

Run with: `python scripts/frontiers_revision/cocktail_source_comparison.py`. Exact source/script checksums, aliases, software versions, reconstructed historical groups, matched scores, all controls, and verification results are supplied alongside this report.
