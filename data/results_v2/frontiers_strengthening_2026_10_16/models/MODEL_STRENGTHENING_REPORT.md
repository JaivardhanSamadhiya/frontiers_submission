# Completed retrospective model strengthening

These are new retrospective analyses of previously accessed sources. They do not restore the original tuning history, create a new untouched test cohort, or change preserved historical predictions.

## StaphStudy second-source matched evaluation

The input contains 1053 pairs (333 positive), all 39 phages × 27 host reference identifiers. Canonical 21-mer sketches of size 2048 at distance 0.05 produce 8 phage and 26 host clusters. Four of five predetermined diagonal blocks pass the minimum three-positive/three-negative rule; they contain 200 test pairs (64 positive; prevalence 0.320). Every regime uses those exact test rows and fold-specific matched training class counts.

Assigned-block eligibility: [{"assigned_block": 0, "n": 70, "positive": 11, "negative": 59, "eligible": true}, {"assigned_block": 1, "n": 15, "positive": 3, "negative": 12, "eligible": true}, {"assigned_block": 2, "n": 100, "positive": 47, "negative": 53, "eligible": true}, {"assigned_block": 3, "n": 15, "positive": 3, "negative": 12, "eligible": true}, {"assigned_block": 4, "n": 15, "positive": 1, "negative": 14, "eligible": false}]

Three sampling seeds are sensitivity runs, not independent datasets. Values below average their pooled and within-fold metrics.

                    method                  regime  auroc_mean  auroc_min  auroc_max  AP_mean  macro_fold_auroc_mean  within_fold_weighted_auroc_mean
additive_positive_fraction  both_clusters_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
additive_positive_fraction   host_cluster_excluded    0.512753   0.505400   0.523208 0.348943               0.559260                         0.600664
additive_positive_fraction  phage_cluster_excluded    0.412473   0.398667   0.422909 0.317754               0.584728                         0.488117
additive_positive_fraction shared_entities_allowed    0.500766   0.479148   0.531824 0.341871               0.611794                         0.543120
                      full  both_clusters_excluded    0.376302   0.353286   0.399874 0.253380               0.558346                         0.342258
                      full   host_cluster_excluded    0.544347   0.505974   0.608513 0.363207               0.621358                         0.552719
                      full  phage_cluster_excluded    0.617456   0.573357   0.649874 0.437050               0.715928                         0.643939
                      full shared_entities_allowed    0.619926   0.551241   0.675379 0.437882               0.668245                         0.598848
    host_positive_fraction  both_clusters_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
    host_positive_fraction   host_cluster_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
    host_positive_fraction  phage_cluster_excluded    0.472024   0.445427   0.490636 0.336115               0.584728                         0.488117
    host_positive_fraction shared_entities_allowed    0.404661   0.387753   0.422679 0.304178               0.572721                         0.434672
                 pair_only  both_clusters_excluded    0.617819   0.605641   0.636604 0.372205               0.624141                         0.548101
                 pair_only   host_cluster_excluded    0.521618   0.352539   0.619773 0.375464               0.573293                         0.585876
                 pair_only  phage_cluster_excluded    0.669003   0.658146   0.686006 0.454202               0.714824                         0.658572
                 pair_only shared_entities_allowed    0.549517   0.525563   0.562155 0.343153               0.635855                         0.600560
                phage_only  both_clusters_excluded    0.638097   0.600299   0.675839 0.431774               0.453475                         0.550125
                phage_only   host_cluster_excluded    0.576631   0.565142   0.586052 0.397026               0.557692                         0.597862
                phage_only  phage_cluster_excluded    0.503332   0.446921   0.538028 0.329202               0.437974                         0.519925
                phage_only shared_entities_allowed    0.578048   0.541131   0.602022 0.386537               0.599704                         0.601079
   phage_positive_fraction  both_clusters_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
   phage_positive_fraction   host_cluster_excluded    0.579714   0.566406   0.589442 0.395571               0.559260                         0.600664
   phage_positive_fraction  phage_cluster_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
   phage_positive_fraction shared_entities_allowed    0.576325   0.539120   0.596680 0.385632               0.558975                         0.597551
            training_prior  both_clusters_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
            training_prior   host_cluster_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
            training_prior  phage_cluster_excluded    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000
            training_prior shared_entities_allowed    0.322495   0.322495   0.322495 0.277698               0.500000                         0.500000

Interpretation: the full model worsens from shared to dual exclusion, but the full four-step ordering is not monotonic and phage-only pooled AUROC does not consistently fall. The training-prior control is 0.322 pooled and exactly 0.500 within every fold, demonstrating a cross-fold prevalence artifact. Phage-only pooled and macro-fold trends disagree. Therefore, report pooled, macro-fold and positive-negative-pair-weighted within-fold metrics together; this is evidence of source-dependent evaluation behavior, not a universal hierarchy. Pair-only features outperform the full model under dual exclusion descriptively. Only four eligible blocks limit inference.

Runtime: 180.7 seconds, 144 CPU fits; 144 fitted tree artifacts retained. Input reconstruction reproduces the historical frozen StaphStudy predictions to maximum absolute error 2.98e-08.

## Fresh compatible-input nested GBM versus MLP

Both families use the same 448 columns, identical outer/inner rows, train-only scaling and three predeclared candidates. Host evaluation uses five group-size-balanced folds; dual evaluation uses the five original dual-exclusion blocks. Inner folds follow the same exclusion axes. Candidate selection uses only mean inner AUROC; there is no random-row fallback, random-row early stopping, calibration or outer-test-based selection. MLP uses 150 fixed epochs; GBM uses 200 or 400 trees. Equal candidate and split counts are not equal wall-clock, capacity, or exhaustive tuning budgets.

      regime method    n  positive  prevalence    auroc  average_precision  macro_fold_auroc  within_fold_pair_weighted_auroc
dual_cluster    GBM  398       294    0.738693 0.763409           0.901207          0.792539                         0.738857
dual_cluster    MLP  398       294    0.738693 0.723116           0.875983          0.770788                         0.664684
host_cluster    GBM 1947      1488    0.764253 0.935105           0.974189          0.941165                         0.941281
host_cluster    MLP 1947      1488    0.764253 0.891437           0.960826          0.889353                         0.888941

Per-fold results:

      regime  fold family  chosen_candidate   n  positive    auroc  average_precision
host_cluster     0    GBM                 2 390       302 0.909241           0.961037
host_cluster     0    MLP                 0 390       302 0.855132           0.935903
host_cluster     1    GBM                 2 390       296 0.907292           0.962702
host_cluster     1    MLP                 0 390       296 0.839815           0.935369
host_cluster     2    GBM                 2 389       292 0.954314           0.984329
host_cluster     2    MLP                 2 389       292 0.871593           0.958542
host_cluster     3    GBM                 0 389       305 0.959758           0.987701
host_cluster     3    MLP                 0 389       305 0.937568           0.981921
host_cluster     4    GBM                 1 389       293 0.975220           0.990967
host_cluster     4    MLP                 1 389       293 0.942655           0.973601
dual_cluster     0    GBM                 0  74        49 0.915102           0.956853
dual_cluster     0    MLP                 1  74        49 0.823265           0.860891
dual_cluster     1    GBM                 0 158       122 0.665984           0.862218
dual_cluster     1    MLP                 0 158       122 0.564435           0.813009
dual_cluster     2    GBM                 2  50        35 0.860952           0.951277
dual_cluster     2    MLP                 1  50        35 0.855238           0.946046
dual_cluster     3    GBM                 1  57        49 0.721939           0.948653
dual_cluster     3    MLP                 2  57        49 0.854592           0.974620
dual_cluster     4    GBM                 1  59        39 0.798718           0.899288
dual_cluster     4    MLP                 2  59        39 0.756410           0.877595

Runtime: 1031.8 seconds; 180 inner fits plus 20 outer refits. 100 recorded warnings reflect fixed-epoch MLP runs reaching their budget; convergence is not claimed. All 20 outer fitted pipelines are retained, with inner scores, exact configurations and all split indices. Inner fitted objects themselves are not retained.

This answers a new, bounded model-family comparison. It does not prove that the historical GBM was selected by a documented comparison, nor does it establish one family as universally superior. Host and dual outer cohorts differ, so compare families within each regime.

## Verification

{
  "unit_tests": 4,
  "artifact_reload_predictions": 164,
  "split_checks": 88,
  "input_hash_checks": 8,
  "staph_assembly_matches_frozen_prediction_max_error": 2.975769042290466e-08,
  "secondsource": {
    "prediction_rows": 16800,
    "split_records": 48,
    "completion": {
      "complete": true,
      "elapsed_seconds": 180.73500000000058,
      "protocol_sha256": "a3f42760999b1f5563da6de60307492cc68b9743b3ee5049da100e8b1d928a50",
      "fits": 144
    }
  },
  "nested_comparison": {
    "prediction_rows": 4690,
    "split_records": 10,
    "completion": {
      "complete": true,
      "elapsed_seconds": 1031.7660000000033,
      "protocol_sha256": "9e2d6b1b0df21bdda4a084008fbc875f13234fff9f04082a139d1d92de6c04c2",
      "inner_fits": 180,
      "outer_refits": 20,
      "convergence_warning_count": 100
    }
  },
  "elapsed_seconds": 4.39100000000326,
  "passed": true
}

No GPU or paid compute was used. No manuscript, portal declaration, public release, or historical result was overwritten by these scripts.