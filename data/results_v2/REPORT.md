# PrecisionPhage: results summary

Leakage-controlled phage–host interaction prediction and bounded fixed-bank cocktail prioritization. All values below are read from the frozen result artifacts in `data/results_v2`.

## 1. Within-source leakage hierarchy

Within the sequence-covered NCBI_HR subset, pooled AUROC declined as biological separation between training and test entities increased:

| Evaluation regime | Pooled GBM AUROC | Mean fold AUROC | Folds |
|---|---:|---:|---:|
| Eligible unseen host species | 0.959 | 0.903 | 28 |
| Unseen host sequence cluster | 0.950 | 0.899 | 26 |
| Unseen phage sequence cluster | 0.847 | 0.858 | 22 |
| Dual cold start | 0.785 | 0.799 | 5 |

Taxonomic and sequence-cluster holdouts measure different forms of novelty and should not be interpreted as one perfectly ordered scale. The pooled differences are descriptive. Because observations share phages and hosts, legacy row-wise DeLong, McNemar, bootstrap, and permutation p-values are not used for confirmatory claims. Fold-level variation and entity-aware uncertainty are emphasized in the manuscript.

## 2. Frozen cross-source sequence evaluation

The fixed 24-feature GBM was trained on 1,947 NCBI_HR pairs and applied without retraining or threshold optimization to all 1,053 sequence-covered StaphStudy pairs.

| Metric | Estimate |
|---|---:|
| AUROC | 0.637 |
| Two-way sequence-cluster bootstrap 95% CI | 0.498–0.739 |
| AUPRC | 0.391 |
| Brier score | 0.466 |
| Expected calibration error | 0.499 |

Because the primary AUROC interval includes 0.5 and calibration is poor, this result does not establish cross-source discrimination, portable calibrated probabilities, or transferable decision thresholds. The test is source-held-out but not prospectively blinded or independently collected because StaphStudy was already represented in VHRnet and the earlier four-feature audit.

## 3. Graph model comparison

Within NCBI_HR, the fixed feature-based GBM had higher pooled AUROC than the evaluated GraphSAGE model in every saved holdout regime. Message passing improved the otherwise identical no-graph neural model in descriptive row-pooled comparisons, including a +0.097 difference in dual cold start. These comparisons remain within-source and are not evidence of external model-class transfer.

## 4. Locked external fixed-bank cocktail validation

A separate protocol was locked before the ColoColi outcome and composition tables were examined. Models trained on the dense 402-host by 96-phage Picard matrix produced frozen scores for 100 unseen ColoColi hosts.

| Metric | Estimate |
|---|---:|
| Cocktail AUROC | 0.807 |
| Host-cluster bootstrap 95% CI | 0.693–0.906 |
| Average precision | 0.951 |
| Host-cluster bootstrap 95% CI | 0.904–0.983 |
| Brier score | 0.127 |
| Expected calibration error | 0.079 |

This result supports ranking externally selected, experimentally tested candidate cocktails for unseen hosts within one fixed *Escherichia* phage bank. Gaborieau et al., not PrecisionPhage, selected the tested compositions. It therefore does not validate de novo cocktail selection, unseen-phage transfer, another organism or bank, resistance suppression, or therapeutic efficacy.

## 5. Exploratory within-source set cover

The NCBI_HR set-cover analysis uses host-cluster-grouped out-of-fold scores. Its selected sets are transformations of sparse within-source predictions, not experimentally validated cocktail designs. Missing matrix cells are unobserved and are not treated as confirmed negatives.

## 6. Deterministic resistance sensitivity analysis

The temporal equations use assumed, unfitted mutation, fitness-cost, adsorption, and cross-resistance parameters. Across the sampled grid, modeled resistant rebound occurs under the evaluated phage strategies. This is conditional behavior of the imposed equations, not a fitted biological prediction, resistance-suppression validation, or efficacy result.

## Primary frozen artifacts

- `leakage_splits_results.json` — within-source leakage-controlled evaluations
- `coldstart_fold_aucs.json` — dual cold-start fold values
- `external_staph_validation.json` — frozen cross-source results and cluster-aware intervals
- `external_cocktail_validation/external_cocktail_validation_summary.json` — fixed-bank cocktail evaluation
- `release_manifest.json` — SHA-256 checksums for the frozen release artifacts
- `EXTERNAL_VALIDATION_PROTOCOL.md` and `EXTERNAL_COCKTAIL_VALIDATION_PROTOCOL.md` — dated analysis protocols

## Figures

- `figure_main.png` — within-source leakage hierarchy and model comparisons
- `fig_external_baselines.png` — external baseline results
- `cocktail_coverage.png` — exploratory within-source set-cover behavior
- `temporal_dynamics.png` — illustrative unfitted-parameter trajectories
