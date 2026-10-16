# PrecisionPhage: five-priority strengthening stage

Recorded 5 October 2026 after previous study outcomes were inspected. Existing
results are retrospective. The BASEL–ECOR outcome table has been located in the
official GenoPHI repository but its pair-level labels have not been inspected
in this stage. This document does not claim preregistration or independently
certified historical blindness. Preserve original predictions and failed runs.

## 1. Candidate frozen external evaluation

Pin official GenoPHI commit f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f.
Inventory BASEL_ECOR_interaction_matrix.csv by identifiers/schema only through
a script that quarantines the original and emits no labels or performance.
Check host/phage overlap with training sources, assay definition, accession
provenance, genome availability and feature compatibility. Public aggregate
paper findings are already known; they are not this project's performance.

The historical full model requires four VHIP-supplied columns unavailable for
arbitrary new pairs. The already saved sequence-only NCBI_HR model is the
candidate frozen model if its precise input schema can be reproduced. Do not
impute unavailable columns or relabel a refit as frozen. The fixed Picard model
cannot score phages outside its fitted bank. Exact genome and identifier overlap
must qualify any unseen-host/unseen-phage description.

Before outcome evaluation, fix eligible rows based only on input resolution,
record model/input/code checksums and write all predictions. AUROC, average
precision, Brier, ten-bin ECE, two-way host/phage cluster uncertainty, all
coverage/exclusion counts and distance strata are retained regardless of
direction. Do not tune against BASEL–ECOR outcomes. If provenance or model
compatibility blocks a defensible test, report the specific blocker.

## 2. Official full-budget comparator

Inspect the pinned official workflow and paper settings before selecting a
command. The requested budget is 25 feature-selection repetitions and 50
modeling runs. These numbers alone do not make the prior adapter an official
reproduction. Inventory clustering, split, weighting, feature-filter and
prediction/ensemble semantics. Use official unmodified routines and distinguish
the packaged workflow from the earlier strictly nested reduced-budget adapter.
Do not compare published scores on different cohorts as paired results.

Diagnose the Windows Python access violation with cached training inputs and
isolated native-fitting processes. Preserve the failed run. Changing native
process lifetime is an execution repair; record it before using checkpoints.
No paid service is launched. Check local feasibility and actual timings before
committing to expensive full-bank repetitions. Official native Windows support,
full-pipeline completion and biological performance must be verified separately.

## 3. Applicability domain

Use existing out-of-fold predictions and saved training distances, and report
each bank/grouping separately. Use fixed distance bins <=1e-4, (1e-4,1e-3],
(1e-3,1e-2], >1e-2 in each source's own units, plus missing-distance counts.
No favorable cutoff selection. Report host counts, class counts, AUROC/AP where
estimable, coverage and uncertainty. Report descriptive distance-error
association and selective screening curves; do not claim a validated deployment
threshold from outcome-informed analysis. Host exclusions already investigated
in ColoColi remain retrospective sensitivities.

## 4. Screening efficiency

On complete assay banks, rank candidate phages for each held-out host using
saved OOF model scores and training-only prevalence, breaking score ties by
stable phage identifier. Evaluate fixed budgets 1, 3, 5, 10 (where available),
hit rate, fraction of active phages recovered, precision, first-active rank and
normalized area under the cumulative hit curve. Include hosts with no active
phage in overall hit rates; condition first-hit rank on actionable hosts and
label that restriction. Compare with exact random-screening expectations from
the observed positive count, and report paired host-group bootstrap differences.

For incomplete Picard labels, do not score unassayed pairs as negative. Report
selection-outcome bounds or an explicitly assay-restricted diagnostic. Do not
claim actual experiments or dollars saved without an experimental cost model.

## 5. Biological error analysis

Use available source metadata such as host phylogroup, sequence type, LPS/O/H
types and phage genus/isolation host; retain missing and small groups. For
externally selected constituent pairs, calculate descriptive errors at the
fixed 0.5 score threshold, and explicitly state that rank scores are not
calibrated and the evaluated pairs were externally selected. Add composition,
distance and prevalence context. Use group-aware uncertainty where estimable.

Do not infer a receptor, defense-system or resistance mechanism from an error
pattern alone. Functional gene annotations or experimental perturbation data
can support hypotheses only after exact entity/feature linkage. Avoid searching
many subgroups for a favorable p-value. New analyses remain exploratory and
must retain failures and limitations in the manuscript/report.
