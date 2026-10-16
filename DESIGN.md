# PrecisionPhage v2: implemented analysis design

The current dated scientific snapshot is 2026.10.15. Historical outputs remain
preserved, but retained code does not make all older analyses current evidence.
Set-cover/resistance illustrations and the RaFAH-inspired proxy are excluded
from current comparative claims. New controls and repairs are retrospective.

This document describes the code that exists in this repository. It does not
list planned components as though they were completed.

## Data layers

The full VHRnet table distributed with the VHIP publication contains 8,849
assayed pairs from three study labels. A
four-feature baseline evaluates LOSO, LOGO, and leave-one-study-out transfer on
that table. The training sequence subset contains 1,947 NCBI_HR pairs. A locked
test then applies the unchanged full sequence GBM to all 1,053 StaphStudy pairs
without retraining. A separate fixed-bank analysis trains on the dense Picard
Escherichia matrix and scores published, experimentally tested cocktails for
100 unseen ColoColi hosts. GNN, external-baseline, NCBI_HR set-cover, and
temporal analyses remain NCBI_HR-only.

Genome resolution is exact after normalization, with deterministic aliases for
NCBI accession version suffixes. Unresolved pairs are excluded from sequence
models rather than zero-imputed.

## Features

Node features comprise canonical 4-mer composition, dinucleotide composition,
codon composition, GC fraction, and genome length. PCA and scaling are fitted
inside each outer training fold.

The 24 saved edge features combine:

- four pair features recomputed from the two frozen genomes;
- four columns supplied by the VHIP input (`k3dist`, `k6dist`, `GCdiff`, and
  `Homology`); and
- sixteen cached sequence-pair features, including multiscale exact-word,
  CRISPR-like, and translated-protein proxies.

Because four inputs are source-supplied rather than regenerated here, the
headline feature vector is not a fully portable predictor for arbitrary new
pairs. The frozen StaphStudy test therefore includes a prespecified
sequence-only sensitivity model excluding those columns; its AUROC (0.639) was
similar to the full model (0.637).

## Splits

Taxonomic LOSO/LOGO folds require at least three positives and three negatives,
so the reported taxonomic pooled estimates cover eligible groups rather than all
323 sequence-covered host taxa.

Sequence-aware grouping uses an in-house forward-strand bottom-k MinHash sketch,
a Mash-form distance transformation, and single-linkage clustering at 0.05. It
must not be described as the Mash program or as a direct 95% ANI calculation.

The dual cold-start routine assigns phage and host clusters to five seeded bins.
Fold `i` tests pairs whose phage and host bins both equal `i`, trains pairs whose
two bins both differ from `i`, and discards one-axis-overlap pairs as a leakage
buffer. Assertions verify that strict test entity clusters do not occur in the
corresponding training set.

## Models

The reference flat baseline is a fixed XGBoost classifier with 400 trees,
depth 5, learning rate 0.05, row/column sampling 0.8 and seed 42. The inherited
settings predate the staged source test, but the original selection rationale
cannot be fully reconstructed. The two June tuning logs use 20 pair inputs,
not the final 24, and lack exact selected per-fold configurations; they do not
justify superiority of this final fixed reference. Better four-feature MLP
source-audit performance is retained. A full-input frozen MLP/GNN source test
is not claimed.

The neural model is an inductive GraphSAGE encoder plus an edge-feature decoder.
PCA, scaling, and the message-passing graph are fitted on the inner-training
partition. The graph contains inner-training positives only; inner-validation
rows do not enter preprocessing or message passing. Isotonic calibration is
fitted on that inner validation split.

Saved `significance_*` and `gnn_ablation_*` artifacts came from separate neural
runs and must not be silently combined. New manuscript tables should name the
artifact family used or rerun both after code changes.

## Evaluation and uncertainty

Reported discrimination metrics include AUROC, AUPRC, and ECE. Saved DeLong,
McNemar, row-permutation, and row-bootstrap tests assume independent rows, an
assumption violated because pairs share phages and hosts. They are retained as
exploratory diagnostics, not definitive inferential evidence. Fold-level
variation is emphasized. The frozen StaphStudy test uses 2,000 valid two-way
phage- and host-sequence-cluster bootstrap replicates.

## Frozen cross-source sequence test

The protocol was locked before held-out prediction. The StandardScaler and
fixed 400-tree XGBoost model were fitted only on NCBI_HR, then applied once to
StaphStudy. Neither labels nor test-derived thresholds entered fitting,
selection, calibration, or exclusion. The primary AUROC was 0.637 (95% CI
0.498-0.739), so cluster-aware uncertainty includes chance and cross-source
discrimination remains inconclusive. Probability transport was poor (ECE
0.499), and exact sequence overlap was zero. This is a source-held-out
evaluation, not a prospectively blinded or independently collected cohort.

## External comparisons

PHIST is the genuine published software applied to staged genomes. Only 35.3%
of evaluated pairs receive a nonzero shared-25-mer score.

The retained `RaFAH_style` method is an in-house random-forest proxy built from
six-frame translation and feature-hashed amino-acid 6-mers. It does not use
published RaFAH weights or HMM protein clusters and is excluded from the current
comparative evidence. No released-tool superiority to unrun RaFAH or CHERRY is
inferred.

## Cocktail analysis

Host-cluster-grouped outer-OOF GBM probabilities are converted to binary
decisions using fold-specific thresholds selected on group-aware inner-OOF
predictions from the corresponding outer training partition. Outer test labels
do not select their thresholds. Set cover is optimized on the predicted binary
matrix and evaluated on observed labels. Cells without an assay are unavailable,
not validated negatives. The exact optimizer is SciPy `milp`/HiGHS, not PuLP/CBC.

The verified nested-threshold k=1 solution contains 176 phages and covers 89.6%
of eligible taxa on observed labels. k=2 and k=3 solutions improve at-least-one
coverage but achieve only 45.2% and 31.3% observed k-fold coverage. The analysis
is exploratory and does not propose a clinically practical 176-phage formulation.

## Locked fixed-bank cocktail scoring

The Picard-to-ColoColi protocol is separate from the VHIP/NCBI_HR pipeline.
For each of 96 known phages, it combines a regularized categorical host model
with a core-genome-distance nearest-neighbour model. Hyperparameters are chosen
globally by host-grouped Picard cross-validation. Per-phage scores are
rank-normalized across the 100 test hosts and averaged. A cocktail receives the
prespecified ranking score `1 - product(1 - p_i)`; this formula is a scoring
rule and not a claim that resistance or infection events are biologically
independent.

Predictions and checksums are written before the held-out outcome files are
opened. The primary result is AUROC 0.807 (95% host-cluster bootstrap CI
0.693-0.906) for observed lytic activity of 100 tested host-triplet records,
against 0.820 prevalence. Because the published study selected the compositions,
this is bounded fixed-bank scoring, not de novo external cocktail design. The
historical normalization needs a cohort: singleton ranks are all 0.5, producing
the same triplet score 0.875. The later development-reference repair is
batch-invariant, but its external evaluation is retrospective.

The 18 repeated selected compositions receive an additional crossed
host-cluster × composition bootstrap sensitivity (2,000 valid draws), which
gives AUROC 0.807 with interval 0.543–0.965. Shared phages across distinct menus
remain outside that grouping model. The 1,000 whole-triplet shuffles preserve
the frequency of all 18 menus but are only a diagnostic, not a formal
exchangeability-based permutation significance test.

The measured paired comparison evaluates two observed triplets per host,
conditional on one fixed baseline composition. Only 17 hosts have discordant
binary outcomes. Historical policy gain 0.050 has crossed interval
−0.039–0.194; meaningful selection improvement is unresolved. The proposed
laboratory use is to order already specified host-cocktail assays under limited
screening capacity, not choose new constituents or prescribe treatment. No
prospective assay-allocation trial or validated laboratory-cost saving occurred.

Experiment 31 documents why the same rank-normalized scores cannot be reused
unchanged to select among phages: they discard between-phage breadth and
underperform a generalist baseline. Experiment 32 is explicitly
development-only and shows that raw cross-phage-comparable nearest-neighbour or
logistic estimates repair this selector error in grouped Picard folds. These
development results are not promoted to locked external efficacy evidence.

## Temporal sensitivity model

The deterministic ODE couples sensitive and resistant host populations to free
phage. The five taxa are selected because they have the most predicted candidate
phages among eligible taxa. A prespecified grid spans complete cross-resistance,
an intermediate dependence, and independent resistance using
`mu ** (1 + alpha * (n_targeting - 1))` for `alpha` in 0, 0.5, and 1. These
structures and all kinetic parameters are assumed rather than fitted.

All temporal parameters now come from `configs/default.yaml`. The analysis is a
mechanistic sensitivity illustration, not biological validation.

## Reproducibility invariants

1. Test labels never enter training tensors or training-positive graphs.
2. Strict split training and test sets share no held-out entity clusters.
3. Missing genomes exclude a pair; they do not become zero features.
4. Every result claim names its data layer and saved artifact family.
5. A clean release passes `experiments/00_validate_release.py` and unit tests.
