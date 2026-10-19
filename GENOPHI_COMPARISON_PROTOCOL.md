# GenoPHI comparison: input preparation and evaluation boundaries

Recorded 2026-10-05, before fitting a new predictor. This is a retrospective
extension: the Pseudomonas interaction outcomes have already been analyzed in
experiment 41. It is not preregistered, prospective, or outcome-blinded.

## First stage: resolve inputs without fitting

Start with the complete published Pseudomonas bank (23 hosts, 19 phages,
437 interactions). Choose this bank because every matrix identifier is a
public assembly or nucleotide accession and its size permits a CPU-only
implementation. Do not select a bank according to new predictor performance.
Retain the published matrix and existing experiment-41 predictions unchanged.

Resolve each host through NCBI Datasets. Use the submitted GenBank assembly's
protein annotation when supplied. If it has no protein annotation, use only
the explicitly paired RefSeq assembly reported by NCBI, retaining both
accessions and the full source report. This is an annotation-source adaptation,
not proof that the proteomes are identical to those used by GenoPHI's authors.
Resolve phage coding-sequence translations from their exact GenBank accessions.
Preserve raw downloads, source URLs, versioned identifiers, retrieval times,
protein counts and SHA-256 checksums. Do not fuzzy-match names or silently
drop entities. Missing proteins block a complete-bank comparison.

Inventory the official GenoPHI release and published data archive. Summary
performance tables are not held-out prediction files and must not be used as
a paired comparison. A generic CatBoost model using PrecisionPhage features
must not be labeled GenoPHI. Preserve upstream software and license separately
from the unrelated yau_2026 project. MMseqs2 15-6f452 is the version reported
in the GenoPHI paper and has a Windows CPU binary; no paid GPU is authorized
or required for input preparation.

## Requirements before the prediction stage

Use the same existing experiment-41 five outer folds, primarily host
patristic connected components at 0.001; host-identity folds are a secondary
sensitivity. Do not change folds to improve either method. Rebuild host
protein clusters and all supervised feature selection using outer-training
hosts only. Assign held-out proteins to training references without outcomes.
Known-bank phages may all be represented, but phage novelty is not tested.
Hyperparameter selection must occur within training-only grouped folds.
Record the implemented search and ensemble budget in an amendment before
fitting, including every departure from the published GenoPHI workflow.

Evaluate paired pair-level AUROC/AP and top-three constituent-union coverage
against the saved neighbour, prevalence and greedy-coverage controls. Use
2,000 paired host-group bootstrap replicates, conditional on fitted models;
small numbers of host groups limit inference. These are measured individual
phage activities, not experimentally measured mixture outcomes. Report all
results, including losses, ties, skipped fits and computational limits.

Completion of sequence preparation or a software smoke test does not complete
the comparison. Do not revise manuscript performance claims until verified
held-out predictions exist. This comparison cannot establish cross-source
frozen-weight transfer or clinical readiness.
