# Full-budget adapted GenoPHI modeling: matched results

Five outer folds and 250 selected ensemble models passed reload and representation checks. The 437 pairs, 23 hosts, and 19 phages are matched across methods.

| Comparison | Metric | Estimate | Lower 95% | Upper 95% |
|---|---|---:|---:|---:|
| genophi | auroc | 0.6908 | 0.5960 | 0.7895 |
| genophi | average_precision | 0.5867 | 0.4264 | 0.7335 |
| neighbour | auroc | 0.5088 | 0.4359 | 0.5973 |
| neighbour | average_precision | 0.3712 | 0.2640 | 0.5063 |
| prevalence | auroc | 0.5765 | 0.4921 | 0.6653 |
| prevalence | average_precision | 0.4039 | 0.2756 | 0.5491 |
| genophi_minus_neighbour | auroc | 0.1821 | 0.0700 | 0.2834 |
| genophi_minus_neighbour | average_precision | 0.2155 | 0.0858 | 0.3467 |
| genophi_minus_prevalence | auroc | 0.1143 | 0.0184 | 0.2114 |
| genophi_minus_prevalence | average_precision | 0.1828 | 0.0530 | 0.3203 |
| genophi | three_phage_union_hit | 0.6522 | 0.4091 | 0.8462 |
| neighbour | three_phage_union_hit | 0.7391 | 0.5417 | 0.9130 |
| prevalence | three_phage_union_hit | 0.7391 | 0.5000 | 0.9200 |
| genophi_minus_neighbour | three_phage_union_hit | -0.0870 | -0.2857 | 0.0741 |
| genophi_minus_prevalence | three_phage_union_hit | -0.0870 | -0.2174 | 0.0000 |


Intervals condition on saved fits and resample 17 host groups. They do not include all training uncertainty. Three-phage union activity derives from individual assays, not mixture experiments. This is not frozen cross-source transfer or an exact end-to-end published GenoPHI reproduction. A disclosed initialization repair retains the best model even when all internal validation MCC scores are negative. Read GENOPHI_OFFICIAL_FULL_BUDGET_PROTOCOL.md for all adaptations. Intervals spanning zero do not establish an advantage. All results are retained regardless of direction.
