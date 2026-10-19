# Scientific snapshot v2026.10.18.1

This release preserves the preceding tags and frozen results while adding
retrospective scientific analyses. Its exact commit is recorded in the release
metadata and verified archive inventory, rather than self-embedded in this file.

## Measured decisions and source dependence

Three measured doublets per host provide a same-size choice among 300 assays.
Repaired scoring selects an active option for 78 of 100 hosts, compared with
74.67 expected under uniform choice; the paired gain intervals include zero.
Only three triplet compositions contain both outcome classes. Conditional
pair-weighted concordance is 0.721 with crossed interval 0.250–0.940.

The matched StaphStudy analysis uses the same 200 eligible test pairs and
training class counts across exclusion conditions. Full-model mean pooled
AUROC changes from 0.620 to 0.376 under shared-to-dual exclusion, while the
macro-fold values are 0.668 and 0.558. Intermediate regimes and phage-only
trends are not consistently ordered. A constant fold prior can have pooled
AUROC below chance despite exactly 0.500 within every fold.

Fresh compatible-input grouped GBM/MLP fits retain all inner scores, outer
pipelines and split indices. Host/dual AUROC is 0.935/0.763 for GBM and
0.891/0.723 for the initial fixed-epoch MLP. One disclosed training-loss stopping
sensitivity gives 0.882/0.713; all 100 later fits stop before the 1,000-epoch cap.
This does not reconstruct historical selection or establish universal superiority.

## Full-budget adapted GenoPHI comparison

Five held-out-host folds cover 437 measured Pseudomonas pairs, 23 hosts, 19
fixed-bank phages and 17 host groups; the positive prevalence is 158/437 (0.362).
The unchanged budget has 25 selection repetitions per fold, 50 modeling
partition repetitions per eligible cutoff and the eight-setting official grid.
There are 2,200 retained search models and 250 selected ensemble models.
All passed the local search, prediction, source and cached-representation replay.

Adapted GenoPHI AUROC/AP is 0.691/0.587, versus nearest-neighbour 0.509/0.371
and training phage breadth 0.577/0.404. Its paired AUROC differences are
0.182 (95% interval 0.070–0.283) and 0.114 (0.018–0.211), respectively.
These host-group intervals condition on saved fits. Top-three single-assay
union hit is 0.652 versus 0.739 for each baseline; these comparisons do not
establish a screening advantage and are not measured mixture efficacy.
The Windows reference adaptation and negative-MCC retention repair are disclosed.
This is not an exact end-to-end published reproduction or frozen cross-source test.

## Reproduction and limitations

The core includes compact scientific inputs, predictions, selections, manifests,
source terms and 174 fitted outer pipelines. Read SCIENTIFIC_EXTENSION_README.md
for exact commands and the separate comparator environment. The read-only
verify_extension_analysis.py checks the saved analyses without fitting or report
regeneration; verify_comparator_metrics.py checks compact comparator evidence.
A separately identified fitted-model/cached-representation companion enables
comparator_public_replay.py when extracted into the same directory. The companion
is not claimed available until its archival file inventory identifies it.
Raw proteome preprocessing and private execution history are outside that replay.

The authentic all-100 Gaborieau comparison remains partial: 67 rescue selections
have no released numeric constituent score. No probabilities were invented.
No new untouched compatible Picard-bank outcome cohort has been established.
Calibration transport, arbitrary cocktail design and clinical efficacy remain
unvalidated. Historical scientific logs are retained as provenance; their legacy
non-independent p-values are not confirmatory evidence.

The earlier v2026.10.15.2 archive remains at https://doi.org/10.5281/zenodo.23426457
and does not identify the new extension. A new archive will be linked only after
publication and independent byte verification. No manuscript, private review,
correspondence, writing script, credential or private runtime transcript is included.
