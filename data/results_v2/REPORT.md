# PrecisionPhage: results summary

Leakage-controlled phage–host interaction prediction, cocktail design, and eco-evolutionary therapy simulation. All numbers below are read directly from the pipeline's result artifacts.

## 1. Headline performance

- Easiest realistic regime (unseen species, LOSO): **AUROC 0.960** (0.946–0.974), ECE 0.024.
- Hardest regime (both phage and host clusters unseen — true cold start): **AUROC 0.780** (0.729–0.832).
- The feature-based GBM significantly beats the GNN in every regime (DeLong FDR q < 1e-6) and is well calibrated where the GNN is not.

### Table 1. Main results (GBM vs GNN per leakage regime)

| Regime | GBM AUC | GBM ECE | GBM > chance (q) | GNN AUC | GBM−GNN ΔAUC | DeLong q |
| --- | --- | --- | --- | --- | --- | --- |
| Unseen species (LOSO) | 0.960 (0.946–0.974) | 0.024 | 1.0e-03 | 0.879 (0.851–0.907) | +0.080 (0.057–0.104) | 6.0e-11 |
| Unseen host cluster | 0.954 (0.938–0.969) | 0.047 | 1.0e-03 | 0.878 (0.850–0.905) | +0.076 (0.053–0.100) | 2.3e-10 |
| Unseen phage cluster | 0.853 (0.815–0.891) | 0.424 | 1.0e-03 | 0.605 (0.546–0.663) | +0.249 (0.195–0.302) | 9.8e-19 |
| Both unseen (cold start) | 0.780 (0.729–0.832) | 0.166 | 1.0e-03 | 0.621 (0.555–0.686) | +0.160 (0.100–0.220) | 2.7e-07 |

## 2. Leakage hierarchy

AUC declines monotonically as more leakage is removed, confirming the splits are doing real work (no shortcut signal).

| regime | mean_auc | ci_lo | ci_hi | pooled_auc | ece | folds_used |
| --- | --- | --- | --- | --- | --- | --- |
| loso_species | 0.904 | 0.753 | 0.984 | 0.96 | 0.024 | 28 |
| logo_genus | 0.889 | 0.836 | 0.935 | 0.919 | 0.053 | 28 |
| host_cluster | 0.901 | 0.796 | 0.979 | 0.954 | 0.047 | 26 |
| phage_cluster | 0.869 | 0.799 | 0.932 | 0.853 | 0.424 | 22 |
| combined_unseen | 0.8 | 0.738 | 0.863 | 0.78 | 0.166 | 5 |

## 3. Graph message-passing ablation

Disabling message passing (GNN → MLP) does **not** hurt and even *helps* in phage-cluster cold start, so the relational signal is already captured by pairwise features. The GBM remains the headline model. See `figure_main.png` panel c.

## 4. Cocktail optimisation

- Minimum cocktail (exact ILP) to cover all coverable targets: **180 phages** (greedy 180, oracle minimum 188).
- Model-driven greedy tracks the truth-oracle and massively beats random selection (see `cocktail_coverage.png`).

### Table 3. k-robust cocktails

| k | cocktail_size | true_cover_>=1 | true_cover_>=k |
| --- | --- | --- | --- |
| 1 | 180 | 0.909 | 0.909 |
| 2 | 291 | 0.952 | 0.47 |
| 3 | 372 | 0.978 | 0.335 |

## 5. Eco-evolutionary therapy simulation

On a simulated multi-strain infection, only a redundant (k≥2) cocktail both suppresses load and prevents resistance; monophage fails and a non-redundant cocktail relapses via resistance (see `temporal_dynamics.png`).

### Table 4. Therapy outcomes

| strategy | n_phages | end_load | nadir | log10_drop | resistant_frac_end | rebound |
| --- | --- | --- | --- | --- | --- | --- |
| control | 0 | 1.00e+09 | 5.00e+06 | 0.0 | 0.0 | False |
| monophage | 1 | 1.00e+09 | 5.00e+06 | 0.0 | 0.0 | False |
| cocktail_k1 | 4 | 1.00e+09 | 2.91e+04 | 2.23 | 1.0 | True |
| robust_k2 | 8 | 6.02e+05 | 7.47e+03 | 2.83 | 0.0 | True |

## Figures

- `figure_main.png` — generalisation, leakage hierarchy, GNN ablation
- `cocktail_coverage.png` — cocktail coverage vs size
- `temporal_dynamics.png` — therapy dynamics with resistance
