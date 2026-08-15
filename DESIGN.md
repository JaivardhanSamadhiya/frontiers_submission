# PrecisionPhage v2: implemented analysis design

This document describes the code that exists in this repository. It does not
list planned components as though they were completed.

## Data layers

The full VHIP table contains 8,849 assayed pairs from three study labels. A
four-feature baseline evaluates LOSO, LOGO, and leave-one-study-out transfer on
that table. The sequence-covered subset contains 1,947 NCBI_HR pairs only; it is
used by all genome, GNN, external-baseline, cocktail, and temporal analyses.

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
headline feature vector is not yet a fully portable predictor for arbitrary new
pairs. A future release should reproduce or remove those four columns.

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

The headline flat baseline is a fixed XGBoost classifier. Nested tuning and
isotonic calibration exist as an ablation but were not used for saved headline
GBM estimates.

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
variation and future two-way entity-aware resampling are more defensible.

## External comparisons

PHIST is the genuine published software applied to staged genomes. Only 35.3%
of evaluated pairs receive a nonzero shared-25-mer score.

The `RaFAH_style` method is an in-house random-forest proxy built from six-frame
translation and feature-hashed amino-acid 6-mers. It does not use the published
RaFAH pretrained weights or HMM protein clusters and is labeled
“RaFAH-inspired proxy” in submission text.

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

## Temporal sensitivity model

The deterministic ODE couples sensitive and resistant host populations to free
phage. The five taxa are selected because they have the most predicted candidate
phages among eligible taxa. The model assumes independent resistance across
targeting phages, making the effective resistance rate `mu ** n_targeting`.
This assumption structurally favors redundant cocktails and is not learned by
the GNN or supported by cross-resistance measurements. Despite that favorable
assumption, neither verified cocktail strategy prevented rebound or modeled
resistant takeover.

All temporal parameters now come from `configs/default.yaml`. The analysis is a
mechanistic sensitivity illustration, not biological validation.

## Reproducibility invariants

1. Test labels never enter training tensors or training-positive graphs.
2. Strict split training and test sets share no held-out entity clusters.
3. Missing genomes exclude a pair; they do not become zero features.
4. Every result claim names its data layer and saved artifact family.
5. A clean release passes `experiments/00_validate_release.py` and unit tests.
