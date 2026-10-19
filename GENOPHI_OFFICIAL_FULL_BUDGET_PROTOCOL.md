# Official GenoPHI modeling routines at full repetition budget

This evaluation is distinct from the failed reduced-budget adapter. It invokes
unmodified official feature-assignment, repeated selection, cutoff generation,
modeling, cutoff assessment, and median-prediction routines at commit
f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f. There are 25 feature-selection
repetitions and 50 modeling repetitions for EVERY eligible feature-frequency
cutoff, not just 50 random seeds at one chosen parameter setting.

Use the previously fixed five Pseudomonas outer folds at patristic grouping
0.001. Each outer test contains only host groups absent from development.
The official package selects hyperparameters and the feature-frequency cutoff
using its own internal evaluation subsets. All such subsets must lie within
outer development; its internally selected scores must not be reported as
unbiased final test performance. The outer test labels cannot enter any
official selection or fitting command. Final test predictions are locked before
the separate audit joins the already-public outcomes.

Protein references use the existing training-only Windows MMseqs2 adapter,
with 0.4 identity, 0.8 coverage and sensitivity 7.5. Exact-query deduplication,
stable best-hit tie resolution, and native scratch handling are documented
execution/representation adaptations. Consequently this is an **outer-held-out
evaluation of official GenoPHI modeling routines**, not a byte-identical
end-to-end reproduction of the published study. Do not label it an exact
official reproduction or claim that every comparator-fairness question is gone.

Unlike the older reduced-budget adapter, representation merging and assignment
use the official feature_selection_optimized and feature_assignment functions.
They see only development hosts when deriving merged feature definitions;
those definitions are then applied to queries. All 19 known-bank phages can
enter their reference construction. The official function uses hash grouping;
explicitly verify that no unequal training patterns were merged by collision.

For the small development matrices (under 500 pairs), use the protein workflow's
automatic 50-feature RFE target, 4 CPU threads, 4 GB CatBoost memory request,
host filtering, hierarchical clustering with requested 20 clusters (the official
routine reduces this for smaller cohorts), dynamic inverse-frequency modeling
weights, and the official eight-setting grid. Preserve the official behavior
that its RFE wrapper does not forward weights_method and therefore uses the
RFE default log10 weights. This behavior is a source-code observation, not an
undocumented repair. SHAP plots are disabled; this does not change fitting.

Frequency cutoffs and feature-count eligibility follow the protein workflow:
3,4,5,6,7,8,9,10,12,15,17,20,22,25,27,30,32,35,37,40,42,45,47,50;
at least five and fewer than 50 retained features. The chosen cutoff maximizes
official development MCC, breaking ties by the smaller cutoff. All 50 models
at that cutoff contribute their median score to the untouched outer hosts.

Each selection/modeling repetition runs in a fresh subprocess. A stage is
accepted only after normal exit and artifact validation. Preserve failed
outputs and stop, rather than silently accepting a partially written model or
silently reducing the budget. A one-fold completion is only a pilot; final
comparative claims require all five outer folds and an independent audit of
model reloads, feature provenance, rows, labels, and matched baselines.

No paid GPU is needed for this CPU configuration. Runtime depends primarily on
protein preprocessing and the number of eligible cutoffs. Report observed
timings instead of promising completion before measuring the full workload.
