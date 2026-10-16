# Official Gaborieau comparator recovery protocol

Date: 15 October 2026 (America/Los_Angeles). This protocol is written after the historical selected-cocktail outcomes and the earlier 64-host comparison were already known. It does not establish a new blind validation.

## Question and allowable evidence

Determine whether the authors' released predictor or released numeric outputs can provide an authentic comparator for all 100 historically selected triplets (300 constituent records), without using external cocktail outcomes for fitting, choosing a predictor, filling missing outputs, or tuning an aggregation rule.

The current public upstream tree will be checked against the pinned source revision `99656d8c020fcd1dbbc30a02ee10274a07655562`. Inventory code, fitted-model files, feature inputs, original step-score output, and compressed cross-validation predictions. Inspect an unmodified-script dry run with output-directory creation intercepted and source files unchanged. No downloaded pickle is to be executed.

## Decision gate before any new outcome evaluation

An official rerun requires either (a) the authors' fitted predictor, preprocessing, and complete external features, or (b) executable released training/inference code, sufficient inputs and a documented way to reproduce the published predictor without selecting it on the 100 cocktail outcomes. A reconstruction that chooses unspecified settings, substitutes features, refits UMAP on the query cohort, or supplies population coverage as a fabricated missing probability does not pass this gate.

The source four-stage recommender and the source pointwise cross-validation classifier are not interchangeable. Rescue selections may have no per-host predicted probability by design. A new reconstruction could be separately labelled, but must not be presented as recovery of the original missing scores or used to claim the required exact comparison is complete.

## If the gate passes

Save a locked all-300 table with exact host, genome identifier, phage, model/stage, score interpretation, source/input/script hashes, software and seed before loading cocktail outcomes. Use the already documented product score only for a common, explicitly imposed aggregation, alongside a source-native cocktail score if one exists. Evaluate all 100 hosts and the 99-host AN03-excluded sensitivity on identical rows as PrecisionPhage, with host-cluster and composition-dependence sensitivities coordinated with the main revision. Do not optimize on target outcomes.

## If the gate fails

Save the exact missing resources, verified audit outputs and errors, a missing-score request table, public contact provenance, and an unsent author-request draft. Preserve the prior 233 released values, 64 complete-host analysis and identification bounds. Do not claim that the full 100-host required comparison has been completed or that the authors have been contacted.

## Restrictions

CPU-only, no paid compute, no changes to upstream code/data, no sending email/issues, and no publishing. The audit does not calculate new performance metrics or read cocktail outcome files unless the gate passes.
