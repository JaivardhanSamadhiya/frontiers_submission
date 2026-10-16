# Historical cocktail reporting controls

Retrospective analysis of saved scores and already known assay outcomes. No model was fitted.

## Discrimination and prevalence

Intervals resample the same historical host clusters for all methods. AP minus prevalence is recalculated within each draw.

| Cohort | Method | AUROC (95% interval) | AP | Prevalence | AP minus prevalence (95% interval) |
| --- | --- | --- | --- | --- | --- |
| historical_100 | Historical product score | 0.807 (0.693–0.906) | 0.951 | 0.820 | 0.131 (0.072–0.195) |
| historical_100 | Maximum constituent | 0.800 (0.683–0.904) | 0.949 | 0.820 | 0.129 (0.071–0.195) |
| historical_100 | Mean constituent | 0.797 (0.678–0.901) | 0.947 | 0.820 | 0.127 (0.068–0.191) |
| historical_100 | Logistic component | 0.788 (0.673–0.892) | 0.944 | 0.820 | 0.124 (0.065–0.186) |
| historical_100 | Neighbour component | 0.703 (0.570–0.825) | 0.921 | 0.820 | 0.101 (0.046–0.159) |
| historical_100 | Training phage susceptibility | 0.327 (0.221–0.443) | 0.780 | 0.820 | -0.040 (-0.079–0.005) |
| historical_100 | Host phylogroup only | 0.705 (0.540–0.844) | 0.893 | 0.820 | 0.073 (0.016–0.139) |
| exclude_AN03_99 | Historical product score | 0.805 (0.689–0.905) | 0.950 | 0.818 | 0.131 (0.072–0.196) |
| exclude_AN03_99 | Maximum constituent | 0.798 (0.679–0.902) | 0.948 | 0.818 | 0.129 (0.071–0.196) |
| exclude_AN03_99 | Mean constituent | 0.795 (0.674–0.899) | 0.945 | 0.818 | 0.127 (0.068–0.192) |
| exclude_AN03_99 | Logistic component | 0.786 (0.672–0.890) | 0.943 | 0.818 | 0.125 (0.066–0.186) |
| exclude_AN03_99 | Neighbour component | 0.699 (0.568–0.822) | 0.918 | 0.818 | 0.100 (0.045–0.159) |
| exclude_AN03_99 | Training phage susceptibility | 0.322 (0.217–0.436) | 0.777 | 0.818 | -0.041 (-0.080–0.003) |
| exclude_AN03_99 | Host phylogroup only | 0.706 (0.540–0.841) | 0.892 | 0.818 | 0.074 (0.017–0.140) |

## Historical-score cohort triage

These are rank positions among already tested host–cocktail records; they are not validated laboratory cost savings.

| Cohort | Target recall | Positive records required | Score-ranked records screened | Expected random-order rank |
| --- | --- | --- | --- | --- |
| historical_100 | 50% | 41 | 43 | 49.89 |
| historical_100 | 80% | 66 | 72 | 80.31 |
| historical_100 | 90% | 74 | 85 | 90.05 |
| exclude_AN03_99 | 50% | 41 | 44 | 50.00 |
| exclude_AN03_99 | 80% | 65 | 71 | 79.27 |
| exclude_AN03_99 | 90% | 73 | 84 | 89.02 |

## Limits

- These are retrospective controls on already inspected external outcomes, not a new frozen validation or preregistration.
- Historical product, maximum and mean scores inherit normalization across the 100-host external cohort; subgroup sensitivity preserves those original scores rather than renormalizing.
- The outcome is the recorded positive assay of each externally selected historical triple. Naming the source column productive in a legacy CSV does not establish productive replication or clinical efficacy.
- One analyzed triple per host ranks records across the cohort; it does not evaluate alternative cocktails for the same host. Source measured pairs/control triples are handled in the separate source-comparison analysis.
- All uncertainty conditions on fixed scores, bank and historical host grouping; three missing-distance hosts receive singleton groups and the historical AN03 name is ambiguous.
- Seven methods and multiple reporting contrasts are exploratory; no multiplicity-adjusted superiority hypothesis or permutation p-value is asserted.
- Composition shuffle is reused as a diagnostic distribution; exchangeability sufficient for a confirmatory null test has not been established.
- Triage ranks use descending score and ascending exact host identifier to break ties. Results concern the already selected tested cohort, not prospective thresholds, laboratory cost savings or within-host selection.
- Random recovery expectation r*(N+1)/(P+1) is the expected rank of the rth positive record under a uniform permutation of this fixed labelled cohort; it is not a real-world cost estimate.

Exact per-method intervals, all 21 paired contrasts in each cohort, top-k results, rank-by-rank records, bootstrap weights, source hashes and training baseline counts are saved alongside this report.
